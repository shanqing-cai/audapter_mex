#include "mex.h"
#include "Audapter.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include "audioIO.h"

#include <windows.h>
#pragma comment(lib, "Winmm.lib")

//#include "procWarpAtom.h"

using namespace std;

#define N_ACTIONS 16
char *actionNames[N_ACTIONS] = {"info", "start", "stop", "setParam", "getParam", 
								"getData", "runFrame", "reset", "outFrame",
								"ost", "pcf", "playTone", "playWave",
								"playToneSeq", "writeToneSeq", "deviceName"};
int actionCode[N_ACTIONS] = {0, 1, 2, 3, 15,
							 4, 5, 6, 7,
							 8, 9, 11, 12, 
							 13, 14, 100};

int getActionNum(int nActions, char **actionNames, int *actionCode, char *actionName) {
	int act = -1;

	for (int i0 = 0; i0 < N_ACTIONS; i0++) {
		if (!strcmp(actionNames[i0], actionName)) {
			act = actionCode[i0];
			break;
		}
	}

	return act;
}

void printHelp() {
	mexPrintf("Audapter - Online digital speech signal perturbation\n");
	mexPrintf("\t[Author: Shanqing Cai (shanqing.cai@gmail.com)]\n\n");
	mexPrintf("Usage: Audapter(command, arg1, arg2, ...)\n");
	mexPrintf("\tCommand list:\n");
	mexPrintf("\t\t0 / info:		Print audio device info\n");
	mexPrintf("\t\t1 / start:		Start audio I/O\n");
	mexPrintf("\t\t2 / stop:		Stop audio I/O\n");
	mexPrintf("\t\t3 / setParam:	Set parameter value\n");
	mexPrintf("\t\t4 / getData:	Get data from the last run\n");
	mexPrintf("\t\t5 / runFrame:	Supply a frame of data for offline processing\n");
	mexPrintf("\t\t6 / reset:		Software state reset for new run\n");
	mexPrintf("\t\t7 / outFrame:	Read outFrameBufPS\n");
	mexPrintf("\t\t8 / ost:		Read online sentence tracking (ost) configuration file\n");
	mexPrintf("\t\t9 / pcf:		Read perturbation configuration (pcf) file\n");
	mexPrintf("\t\t\n");
	mexPrintf("\t\t11 / playTone:	Play tone\n");
	mexPrintf("\t\t12 / playWav:	Play pre-supplied sound waveform (see 3 - datapb)\n");
	mexPrintf("\t\t13 / playToneSeq:	Play pre-configured tone sequence\n");
	mexPrintf("\t\t14 / writeToneSeq:	Write the waveform of the last tone sequence to wav file\n");	
	mexPrintf("\t\t\n");
	mexPrintf("\t\t100 / deviceName:Set audio device name\n");
	mexPrintf("\t\t\n");
	mexPrintf("\n");
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	unsigned int i0;
	int action;	
	dtype  *signal_ptr,*data_ptr;
	const dtype  *algosignal_ptr, *algodata_ptr, *algobuf_ptr;
	//Parameter::paramType retval;
	dtype   value= 0;
	int size=1;
	int vecsize=1;
	int toPrompt=0;

	static Audapter audapter;		//SC audapter is a Audapter class object. Notice it is a static pointer.
	static audioIO audio_obj;	//SC Refer to the audioIO project
	static	DeviceParams devpar;//SC DeviceParams is a struct, fields: num, chans, fs, set
	static bool started = 0;

	int frame_len;
	int t_downFact;
	int activeDeviceNum;
	int bVerbose;
	bool bVerb = true;
	char actionStr[512];
	char tmpMsg[512];

	mwSize inParNDims;
	const mwSize *inParSize;
	size_t nPars;
	bool isActionString;

	char *paramName = NULL;

	if (nrhs == 0) {
		printHelp();
		return;
	}
	else if (nrhs > 4) {
		mexErrMsgTxt("Invalid syntax");
	}

	/* Determine if the first input argument is a string */
	inParNDims = mxGetNumberOfDimensions(prhs[0]);
	if (inParNDims == 2) {
		inParSize = mxGetDimensions(prhs[0]);

		if (inParSize[0] == 1) {
			if (!mxGetString(prhs[0], actionStr, inParSize[1] + 1)) {
				isActionString = true;
			}
			else {
				isActionString = false;
			}
		}
	}
	else {
		mexErrMsgTxt("Erroneous number of dimensions in the first input argument");
	}
	
	if (isActionString) { /* Character-string action */
		action = getActionNum(N_ACTIONS, actionNames, actionCode, actionStr);
		if (action == -1) {
			sprintf_s(tmpMsg, sizeof(tmpMsg), "Unrecognized command: %s", actionStr);
			mexErrMsgTxt(tmpMsg);
		}

		if (!strcmp(actionStr, "deviceName")) {
			
		}
	}
	else {
		action = (int) floor(*(double*) mxGetPr(prhs[0]));	//SC Figure out the input argument
	}
		/* if (~mxGetString(prhs[0], actionStr, )) {
			action = (int)floor(*(double*)mxGetPr(prhs[0]));	//SC Figure out the input argument
		} */

		//LPCWSTR in_wav_fn = L"E:\\speechres\\blueshift\\mcode\\test1.wav";

	if (action == 1 || action == 11 || action == 12 || action == 13) { /* Audio I and/or O required: determine the device number */
		if (strlen(audapter.deviceName) == 0) { /* In which case we select the first suitable device */
			for(i0 = 0; i0 < audio_obj.devices.size(); i0++) {
				if ((audio_obj.devices[i0].inputChannels > 0) && (audio_obj.devices[i0].outputChannels > 0)) { // Find the first active output device
					activeDeviceNum = i0 + 1;
					break;
				}
			}

			if (i0 == audio_obj.devices.size()) {
				mexErrMsgTxt("Cannot find valid audio device");
			}
		}
		else {
			bool devFound = false;
			for(i0 = 0; i0 < audio_obj.devices.size(); i0++) {
				if (strstr(audio_obj.devices[i0].name.c_str(), audapter.deviceName) != NULL
					&& audio_obj.devices[i0].inputChannels > 0
					&& audio_obj.devices[i0].outputChannels > 0) { // Find the first active output device
					activeDeviceNum = i0 + 1;
					devFound = true;
					break;
				}
			}

			if (!devFound) {
				sprintf_s(tmpMsg, sizeof(tmpMsg), "Cannot find device with specified name: \"%s\"", audapter.deviceName);
				mexErrMsgTxt(tmpMsg);
			}
		}
	}
		
	switch (action){
		case 0:		//SC enumerate all the audio devices
			mexPrintf("Information: \n");
			mexPrintf("User-selected deviceName = \"%s\"\n", audapter.deviceName);

			for(i0 = 0; i0 < audio_obj.devices.size(); i0++)
			{
				if (audio_obj.devices[i0].inputChannels > 0)	//SC input devices
					mexPrintf("[I]: %d : %s; ", i0, audio_obj.devices[i0].name.c_str());
				if (audio_obj.devices[i0].outputChannels > 0)	//SC input devices
					mexPrintf("[O]: %d : %s", i0, audio_obj.devices[i0].name.c_str());

				if (audio_obj.devices[i0].inputChannels > 0 && audio_obj.devices[i0].outputChannels > 0) {
					if (strlen(audapter.deviceName) == 0) {
						mexPrintf("\n");
					}
					else {
						if (strstr(audio_obj.devices[i0].name.c_str(), audapter.deviceName) != NULL) {
							mexPrintf(" (Match)\n");
						}
						else {
							mexPrintf(" (Non-match)\n");
						}
					}
				}
			}

			frame_len = *((int *) audapter.getParam("framelen"));
			t_downFact = *((int *) audapter.getParam("downfact"));

			audio_obj.setcallbackparams(frame_len * t_downFact, (void *)&algoCallbackFunc, (void *)&audapter);			
			break;

		case 1:			//SC set audio device parameters and start the action
			frame_len = *((int *) audapter.getParam("framelen"));
			t_downFact = *((int *) audapter.getParam("downfact"));

			audio_obj.setcallbackparams(frame_len * t_downFact, (void *)&algoCallbackFunc, (void *)&audapter);
			
			devpar.fs = *((int *) audapter.getParam("srate")) * t_downFact;	//SC device sampling rate

			devpar.num = activeDeviceNum;
			devpar.chans = 1;
			audio_obj.setdevparams(&devpar, 1);
			devpar.num = activeDeviceNum;
			devpar.chans = 2;
			audio_obj.setdevparams(&devpar, 2);

			if (!started)
			{
				printf("Action Start:  %d\n", action);
				audapter.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;

		case 2:			//SC Stop the audio device
			if (started){
				printf("Action End:  %d\n", action);
				if (audio_obj.started)
					audio_obj.stopdev();
			}
			else {
				printf("Not started\n");
			}
			started = 0;
			break;
			
		case 3:			//SC Set paramters of the audapter object
			inParNDims = mxGetNumberOfDimensions(prhs[2]);
			inParSize = mxGetDimensions(prhs[2]);

			if (inParNDims > 2) {
				mexErrMsgTxt("ERROR: Input parameter has too many dimensions (> 2)");
				return;
			}

			if (inParNDims == 2) {
				if (inParSize[0] == 1)
					nPars = inParSize[1];				
				else if (inParSize[1] == 1)
					nPars = inParSize[0];
				else
					mexErrMsgTxt("ERROR: Input parameter is not a scalar or a row or column vector");
			}

			if (nrhs >=4)
				bVerb = ((int) floor(*(double*) mxGetPr(prhs[3])) == 1);

			if (nrhs >= 3)
				audapter.setParam(mxArrayToString(prhs[1]), (void *)mxGetPr(prhs[2]), nPars, bVerb);

			break;

		case 15:			//SC Get paramters of the audapter object
			if (nrhs < 2)
				mexErrMsgTxt("getParam: name not supplied");
			
			inParNDims = mxGetNumberOfDimensions(prhs[1]);
			if (inParNDims != 2) {
				mexErrMsgTxt("Unexpected format in parameter name");
			}

			inParSize = mxGetDimensions(prhs[1]);
			paramName = new char[inParSize[1] + 1];
			mxGetString(prhs[1], paramName, inParSize[1] + 1);
			
			if (nlhs == 1) {
				audapter.queryParam(paramName, &(plhs[0]));
			}
			else if (nlhs > 1) {
				mexErrMsgTxt("Too many output arguments");
			}

			if (paramName) {
				delete [] paramName; paramName = NULL;
			}

			break;

		case 4:			//SC Manually get three outputs: 1
			// signal
			algosignal_ptr = audapter.getSignal(size);	//SC size is updated in this process: size = frame_counter*p.frameLen
			if (size > 0)
			{
				plhs[0] = mxCreateDoubleMatrix(size, 2, mxREAL);
				signal_ptr = mxGetPr(plhs[0]);
				for(int ii = 0; ii < size; ii++)
				{
					signal_ptr[ii] = algosignal_ptr[ii];
					signal_ptr[size + ii] = algosignal_ptr[audapter.getMaxRecSize() + ii];
				}

				// data
				algodata_ptr = audapter.getData(size,vecsize);
				plhs[1] = mxCreateDoubleMatrix(size,vecsize, mxREAL);
				data_ptr = mxGetPr(plhs[1]);

				for(int ii = 0; ii < size; ii++)
				{
					for(int jj = 0; jj < vecsize; jj++)
						data_ptr[jj * size + ii] = algodata_ptr[jj * audapter.getMaxDataSize() + ii];
				}				
			}
			else
			{
				plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
				plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
				plhs[2] = mxCreateDoubleMatrix(1,1, mxREAL);
			}

			break;

		case 5:				//SC manually supply the buffer and run the process
			frame_len = *((int *) audapter.getParam("framelen"));
			t_downFact = *((int *) audapter.getParam("downfact"));
			data_ptr = (double*)mxGetPr(prhs[1]);	//SC pointer to buffer
			algoCallbackFuncMono((char *)data_ptr, frame_len * t_downFact, (void *)&audapter);
			break;

		case 6:				//SC manually reset audapter object
			audapter.reset();
			break;

		case 7:				// Read out the outFrameBufPS
			t_downFact = *((int *) audapter.getParam("downfact"));
			algobuf_ptr = audapter.getOutFrameBufPS();
			plhs[0] = mxCreateDoubleMatrix(audapter.getMaxFrameLen() * t_downFact * audapter.getMaxDelayFrames(), 1, mxREAL);
			signal_ptr = mxGetPr(plhs[0]);
			for (int i = 0; i < (audapter.getMaxFrameLen() * t_downFact * audapter.getMaxDelayFrames()); i++) {
				signal_ptr[i] = algobuf_ptr[i];
			}
			break;

		case 8:				// Set OST file name and read the file (10/18/2012)			
			if (nrhs == 2 || nrhs == 3) {
				if (nrhs == 2)
					bVerbose = 1;
				else
					bVerbose = (int)(*(double *) mxGetPr(prhs[2]));

				strcpy_s(audapter.ostfn, sizeof(audapter.ostfn), mxArrayToString(prhs[1]));
				if (bVerbose)
					printf("ostfn = %s\n", audapter.ostfn);
				//audapter.readOSTTab(bVerbose);//Marked
				audapter.readOSTTab(bVerbose);
			}
			else {
				printf("Syntax error.\n");
			}

			break;
				
		case 9:				// Set PIP cfg file name and read the file (10/19/2012)
			if (nrhs == 2 || nrhs == 3) {
				if (nrhs == 2)
					bVerbose = 1;
				else
					bVerbose = (int)(*(double *) mxGetPr(prhs[2]));

				strcpy_s(audapter.pipcfgfn, sizeof(audapter.pipcfgfn), mxArrayToString(prhs[1]));				
				if (bVerbose)
					printf("pipcfgfn = %s\n", audapter.pipcfgfn);
				audapter.readPIPCfg(bVerbose);
			}
			else {
				printf("Syntax error.\n");
			}

			break;

		case 11:			//SC Sine wave generator
			frame_len = *((int *) audapter.getParam("framelen"));
			t_downFact = *((int *) audapter.getParam("downfact"));

			audio_obj.setcallbackparams(frame_len * t_downFact, (void *)&algoCallbackFuncSineGen, (void *)&audapter);		

			devpar.num = activeDeviceNum;
			devpar.chans = 2;
			devpar.fs = *((int *) audapter.getParam("srate")) * t_downFact;
			audio_obj.setdevparams(&devpar, 1);
			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar, 2);


			if (!started)
			{
				printf("Action Start:  %d\n", action);
				audapter.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;

		case 12:			//SC wave playback
			frame_len = *((int *) audapter.getParam("framelen"));
			t_downFact = *((int *) audapter.getParam("downfact"));

			audio_obj.setcallbackparams(frame_len * t_downFact, (void *)&algoCallbackFuncWavePB, (void *)&audapter);
					
			devpar.fs = *((int *) audapter.getParam("srate")) * t_downFact;	//SC device sampling rate
			devpar.num = activeDeviceNum;
			devpar.chans = 1;
			audio_obj.setdevparams(&devpar,1);
			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar,2);

			if (!started)
			{
				printf("Action Start:  %d\n", action);
				audapter.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;

		case 13:		//SC Tone sequence generator
			frame_len = *((int *) audapter.getParam("framelen"));
			t_downFact = *((int *) audapter.getParam("downfact"));

			audio_obj.setcallbackparams(frame_len * t_downFact, (void *)&algoCallbackFuncToneSeq, (void *)&audapter);
				
			devpar.fs = *((int *) audapter.getParam("srate")) * t_downFact;	//SC device sampling rate

			audapter.tsgRecCounter = 0;

			devpar.num = activeDeviceNum;
			devpar.chans = 2;
			audio_obj.setdevparams(&devpar, 1);
			devpar.num = activeDeviceNum;
			audio_obj.setdevparams(&devpar, 2);

			if (!started)
			{
				printf("Action Start:  %d\n", action);
				audapter.reset();			//SC Reset audio device params
				audio_obj.startdev();	//SC Start audio device. Callback function will be automatically called
			}
			else
				printf("Already started\n");
			started = 1;
			break;

		case 14: /* Write the waveform of the result of last tone sequence generation to text file */
			/* TODO */
			char *tsgWavFN; /* ASCII format */
			FILE *tsgWavF;
			int nwrit, n;

			if (nrhs != 2) {
				mexPrintf("Usage: TransShiftMex(14, tsgWavASCIIFN);\n");
				return;
			}
			tsgWavFN = mxArrayToString(prhs[1]);

			/* if ((tsgWavF = fopen(tsgWavFN, "wt")) == NULL) { */
			if (fopen_s(&tsgWavF, tsgWavFN, "wt")) {
				mexPrintf("ERROR: Failed to open file %s for writing\n", tsgWavFN);
				return;
			}

			for (n = 0; n < audapter.tsgRecCounter; n++) {
				if ((nwrit = fprintf(tsgWavF, "%.6f\n", audapter.tsg_wf[n])) < 0) {
					mexPrintf("ERROR: Error occurred during writing to file: %s\n", tsgWavFN);
					fclose(tsgWavF);
					return;
				}
			}

			fclose(tsgWavF);
			mexPrintf("Sucessfully wrote %d samples of tone-sequence waveform to ASCII file: %s\n", n, tsgWavFN);

			break;

		case 21:	//SC(2012/02/29) Write data to binary file
			algosignal_ptr = audapter.getSignal(size);	//SC size is updated in this process: size = frame_counter*p.frameLen
			if (size > 0){
				//printf("Writing data to file data0.bin: size = %d (%d bytes)\n", size, size * (sizeof dtype));
				printf("Writing data to file data0.bin: size = %d (%d bytes)\n", audapter.getMaxRecSize(), audapter.getMaxRecSize() * (sizeof dtype));
				ofstream dataFileCl ("data0.bin");
				dataFileCl.close();

				fstream dataFile("data0.bin", ios::binary | ios::out);
				if (!dataFile){
					printf("WARNING: Cannot open file data0.bin\n");
				}

				dataFile.write((char *) algosignal_ptr, audapter.getMaxRecSize() * (sizeof dtype));
				dataFile.close();
			}
			break;

		case 22:			//SC Manually get the entire recording buffer
			// signal
			//algosignal_ptr = audapter.getSignal(size);	//SC size is updated in this process: size = frame_counter*p.frameLen
			plhs[0] = mxCreateDoubleMatrix(audapter.getMaxRecSize(), 1, mxREAL);
			signal_ptr = mxGetPr(plhs[0]);
			algosignal_ptr = audapter.getSignal(size);
			for(int i = 0; i < audapter.getMaxRecSize(); i++) {
				signal_ptr[i] = algosignal_ptr[i];
			}
				
			break;

		case 23:			// SC Manually get the entire playback buffer
			plhs[0] = mxCreateDoubleMatrix(audapter.getMaxRecSize(), 1, mxREAL);
			signal_ptr = mxGetPr(plhs[0]);
			algosignal_ptr = audapter.getSignal(size);
			for(int i = 0; i < audapter.getMaxRecSize(); i++)
			{
				signal_ptr[i] = algosignal_ptr[audapter.getMaxRecSize() + i];
			}

				
			break;


		case 30: // Set the output base;
			if (nrhs != 2) {
				printf("Usage: TransShiftMex(30, wavOutBase);\n");
				return;
			}			

			sprintf_s(audapter.wavFileBase, "%s", mxArrayToString(prhs[1]));
			printf("audapter.wavFileBase = %s\n", audapter.wavFileBase);

			break;

		case 31:
			if (nrhs != 1) {
				printf("Usage: TransShiftMex(31);\n");
				return;
			}
			audapter.writeSignalsToWavFile();
			break;

		case 100:
			if (nrhs != 2) {
				mexErrMsgTxt("deviceName: Invalid syntax");
			}
			else {
				inParNDims = mxGetNumberOfDimensions(prhs[1]);
				if (inParNDims != 2) {
					mexErrMsgTxt("Input deviceName is not a char string");
				}
				else {
					inParSize = mxGetDimensions(prhs[1]);
					if (inParSize[0] != 1) {
						mexErrMsgTxt("Input deviceName is not a char string");
					}
					else {
						if (mxGetString(prhs[1], audapter.deviceName, inParSize[1] + 1)) {
							mexErrMsgTxt("Input deviceName is not a char string");
						}
						else {
							sprintf_s(tmpMsg, sizeof(tmpMsg), "Set deviceName to: %s\n", audapter.deviceName);
							mexPrintf(tmpMsg);
						}
					}
				}
			}
			break;

		default:
			printf("Action unknown:  %d\n", action);
				
	}


}
