//#include "noisecancel.h"
#include "audioIO.h"
#include "RTAudio.h"
#include <math.h>
#include <stdio.h>

#define FORMAT RTAUDIO_FLOAT64
RtAudio *audio = NULL;

audioIO::audioIO()
{
	input_dev = new DeviceParams;
	output_dev = new DeviceParams;
	input_dev->set = 0;
	output_dev->set = 0;
	params_set = 0;
	//int a = 1;
	//a = 2;

	RtAudioDeviceInfo info;
	try {
		audio = new RtAudio();
	}  
	catch (RtError &error) {
		error.printMessage();
	}
	int num_devices = audio->getDeviceCount();
	devices.resize(num_devices);
	for (int i=0; i<num_devices; i++) {
		try {
			info = audio->getDeviceInfo(i+1);
		}
		catch (RtError &error) {
			error.printMessage();
			break;
		}
		devices[i].name = info.name;
		devices[i].outputChannels = info.outputChannels;
		devices[i].inputChannels	= info.inputChannels;
		devices[i].duplexChannels = info.duplexChannels;
		devices[i].isDefault		= info.isDefault;
		devices[i].sampleRates	= info.sampleRates;
      }
	delete audio;
	audio = NULL;
	started = 0;
}

audioIO::~audioIO()
{
	delete input_dev;
	delete output_dev;
	if (audio != NULL)
		delete audio;
}

int audioIO::setcallbackparams(int size, void* func, void* data)
{
	algosize = size;
	algofunc = func;
	algodata = data;
	return 0;
}


int audioIO::setdevparams(DeviceParams* devpar, int dev_code)
{
	if ((dev_code == 0) || (dev_code == 1))
	{
		input_dev->num = devpar->num;
		input_dev->fs  = devpar->fs;
		input_dev->chans = devpar->chans;
		input_dev->set = 1;
	}
	if ((dev_code == 0) || (dev_code == 2))
	{
		output_dev->num = devpar->num;
		output_dev->fs  = devpar->fs;
		output_dev->chans = devpar->chans;
		output_dev->set = 1;
	}
	if (input_dev->set && output_dev->set)
		params_set = 1;
	return 0;
}

int audioIO::startdev()
{
	buffer_size = algosize; //int(floor(input_dev->fs*0.008));
	//printf("startdev: buffer_size = %d\n", buffer_size);
	try {
		if (audio == NULL)
			//printf("startdev: creating new audio object with buffer_size = %d\n", buffer_size);
			audio = new RtAudio(output_dev->num, output_dev->chans, input_dev->num, input_dev->chans,FORMAT, input_dev->fs, &buffer_size, 1);
	}
	catch (RtError &error) {
		error.printMessage();
		//	  exit(EXIT_FAILURE);
		return 1;
	}
	try {
		audio->setStreamCallback((RtAudioCallback)algofunc, algodata);
		audio->startStream();
	}
	catch (RtError &error) {
		error.printMessage();
		stopdev(1);
		return 1;
	}
	started = 1;
	return 0;
}

int audioIO::stopdev(int close)
{
	if (close == 0)
	{
		try {
			audio->stopStream();
		}
		catch (RtError &error) {
			error.printMessage();
			return 1;
		}
	}
	audio->closeStream();
	delete audio; audio = NULL;
	started = 0;
	return 0;
}