/*
pcf.cpp
	Perturbation configuration 
	Part of Audapter

	Shanqing Cai, 2013
*/

#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>

#include "mex.h"
#include "pcf.h"
#include "utils.h"

using namespace std;

/* pvocWarpAtom: Default constructor */
pvocWarpAtom::pvocWarpAtom() {
	tBegin = 0.25;
	dur1 = 0.25; 
	rate1 = 0.5;
	durHold = 1;
	rate2 = 2;

	/* ostInitState = -1; */
	ostInitState = -1; /* Use a negative number to effectively disable time warping by default */
			
	dur2 = (1 - rate1) / (rate2 - 1) * dur1;

	tEnd = tBegin + dur1 + durHold + dur2;
}

/* pvocWarpAtom: with input arguments: Version 1 */
pvocWarpAtom::pvocWarpAtom(double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2) {
	tBegin = t_tBegin;
	rate1 = t_rate1;
	dur1 = t_dur1;
	durHold = t_durHold;
	rate2 = t_rate2;

	if ( !(rate1 > 0.0 && rate1 < 1.0) || !(rate2 > 1.0) )
		throw rateError();

	ostInitState = 0;

	dur2 = (1 - rate1) / (rate2 - 1) * dur1;

	tEnd = tBegin + dur1 + durHold + dur2;
}

/* pvocWarpAtom: with input arguments: Version 2 */
pvocWarpAtom::pvocWarpAtom(int t_ostInitState, double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2) {
	tBegin = t_tBegin;
	rate1 = t_rate1;
	dur1 = t_dur1;
	durHold = t_durHold;
	rate2 = t_rate2;

	if ( !(rate1 > 0.0 && rate1 < 1.0) || !(rate2 > 1.0) )
		throw rateError();

	ostInitState = t_ostInitState;

	dur2 = (1 - rate1) / (rate2 - 1) * dur1;

	tEnd = tBegin + dur1 + durHold + dur2;
}

/* pvocWarpAtom: Test if the input time t is within the time-shift period: Variant 2: with initial state number */
const bool pvocWarpAtom::procTimeWarp(const int stat, const int statOnsetIndex, 
										  const int nDelay, const double frameDur, 
										  const double & t, double & wt) const {
/* Input arguments: 
		stat:			current OST status number
		statOnsetIndex: frame index at which the current status is first entered
		nDelay:			Audapter nDelay
		frameDur:		Audapter frame duration (in s, = frameLen / sr)
		t:				current time (s)
		wt:				warped time (s)

	Return value:
		true if the current time (t) is during this time-warping event
		false otherwise
*/

	double tt = t;
	if ((ostInitState < 0) || (stat < ostInitState)) {
		return false;
	}
	else {
		double t01;
		bool duringTimeWarp;

		t01 = tt - static_cast<double>(statOnsetIndex) * frameDur;
		//t01 = tt - static_cast<double>(statOnsetIndex - (nDelay - 1)) * frameDur;
		duringTimeWarp = (t01 >= tBegin) && (t01 < tEnd);

		if (duringTimeWarp)
			tt = t01;			

		/* Determine warped time */
		if (tt < tBegin + dur1){ /* Time dilation (deceleration) */
			wt = (tt - tBegin) * rate1 + tBegin;
		}
		else if (tt < tBegin + dur1 + durHold){ /* Time shifting (no compression or dilation) */
			wt = rate1 * dur1 - dur1 + tt;
		}
		else if (tt < tBegin + dur1 + durHold + dur2){ /* Time compression (acceleration) at the end of the warp interval */			
			wt = tBegin + dur1 + durHold + dur2;
			wt -= (tBegin + dur1 + durHold + dur2 - tt) * rate2;
		}
			
		wt += static_cast<double>(statOnsetIndex) * frameDur;
		//wt += static_cast<double>(statOnsetIndex - (nDelay - 1)) * frameDur;
		/* ~Determine warped time */

		return duringTimeWarp;
	}

}

/* PERT_CFG: Constructor */
PERT_CFG::PERT_CFG() {
	n = 0;

	pitchShift = NULL;
	intShift = NULL;
	fmtPertAmp = NULL;
	fmtPertPhi = NULL;

	/*warpCfg = new pvocWarpAtom();*/
}

/* PERT_CFG: Destructor */
PERT_CFG::~PERT_CFG() {
	if (pitchShift)
		free(pitchShift);

	if (intShift)
		free(intShift);

	if (fmtPertAmp)
		free(fmtPertAmp);

	if (fmtPertPhi)
		free(fmtPertPhi);

	/*if (warpCfg)
		delete warpCfg;*/
}

/* Subroutine: readline */
int readline(FILE *fp, char *line) {
	int lineLen = 0;
	char c;

	c = fgetc(fp);
	while (c != '\n' && c != '\r' && c != EOF) {		
		line[lineLen ++] = c;
		c = fgetc(fp);
	}

	line[lineLen] = '\0';

	return lineLen;
}

/* Subroutine: deblank a char string */
int deblank(char * line) {
	int lineLen = strlen(line);
	int newLineLen, i;
	int ine0 = 0;
	int ine1 = lineLen - 1;

	while (line[ine0] == ' ' || line[ine0] == '\t')
		ine0++;

	while (line[ine1] == ' ' || line[ine1] == '\t') {
		ine1--;
		if (ine1 == ine0)
			break;
	}

	newLineLen = ine1 - ine0 + 1;
	for (i = 0; i < newLineLen; i++)
		line[i] = line[i + ine0];
	line[newLineLen] = '\0';
	
	return newLineLen;
}


/* Subroutine: string_count_char */
int string_count_char(char *str, char c) {
	unsigned int i0;
	int cnt = 0;

	for (i0 = 0; i0 < strlen(str); i0++)
		if (str[i0] == c)
			cnt++;

	return cnt;
}

/* Subroutine: sscanf_floatArray */
int sscanf_floatArray(char *str, double *xs, int nx) {
	int xc = 0, j;
	unsigned int i;
	bool tBreak = false;
	char tmpstr[256];

	i = 0;
	while (xc < nx) {
		j = 0;

		while (str[i] == ',' || str[i] == ' ' || str[i] == '\t' || str[i] == ';')
			i++;

		tmpstr[j++] = str[i++];
		while (tmpstr[j - 1] != ',' && tmpstr[j - 1] != ' ' && tmpstr[j - 1] != ';' && tmpstr[j - 1] != '\t' 
			   && tmpstr[j - 1] != EOF && tmpstr[j - 1] != '\0')
			tmpstr[j++] = str[i++];
	
		if (i >= strlen(str))
			tBreak = true;

		tmpstr[j] = '\0';
		xs[xc ++] = atof(tmpstr);

		if (tBreak)
			break;
	}

	return xc;
}

/* Check whether the input warp event has any time overlap with exisiting warp events */
const bool PERT_CFG::checkWarpIntervalsOverlap(const pvocWarpAtom t_warpCfg) {
/* Return value: 
	true if there is one or more overlaps 
	false if there is none 

	Limited to checking within the same ostInitStates
*/
	for (list<pvocWarpAtom>::const_iterator w_it = warpCfg.cbegin();
		 w_it != warpCfg.cend();
		 ++w_it) {
		if (t_warpCfg.ostInitState != w_it->ostInitState) {
			continue;
		}
		else {
			if ( !((t_warpCfg.tBegin > w_it->tEnd) || 
				   (t_warpCfg.tEnd < w_it->tBegin)) )
				return true;
		}
	}

	return false;
}

/* Add time-warping event */
void PERT_CFG::addWarpCfg(double t_tBegin, double t_rate1, 
						  double t_dur1, double t_durHold, double t_rate2)
	throw(overlappingWarpIntervalsError) {

	try {
		pvocWarpAtom t_warpCfg(t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);

		/* Check to make sure that there is no overlapping time warp periods */
		if ( !checkWarpIntervalsOverlap(t_warpCfg) )
			warpCfg.push_back(t_warpCfg);
		else
			throw overlappingWarpIntervalsError();
	}
	catch (pvocWarpAtom::rateError) {
		throw warpCfgInitError();
	}

}

/* Add time-warping event */
void PERT_CFG::addWarpCfg(int t_ostInitState, double t_tBegin, double t_rate1, 
						  double t_dur1, double t_durHold, double t_rate2)
	throw(overlappingWarpIntervalsError) {

	try {
		pvocWarpAtom t_warpCfg(t_ostInitState, t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);

		/* Check to make sure that there is no overlapping time warp periods */
		if ( !checkWarpIntervalsOverlap(t_warpCfg) )
			warpCfg.push_back(t_warpCfg);
		else
			throw overlappingWarpIntervalsError();
	}
	catch (pvocWarpAtom::rateError) {
		throw warpCfgInitError();
	}

}


/* Read from file */
void PERT_CFG::readFromFile(const string pertCfgFN, const int bVerbose) 
	throw(pcfFileReadingError, pcfFileSyntaxError)
{
	int nTimeWarpAtoms = -1; /* SCai: currently a dummy variable. TODO: implement multiple time warping atoms */
	int t_ostInitState;
	double t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2;

	list<string> lines_0 = readLinesFromFile(pertCfgFN);
	if (lines_0.empty())
		throw pcfFileReadingError();

	/* Trim lines; remove empty lines; remove commented lines */
	list<string> lines_1;
	for (list<string>::const_iterator lit = lines_0.begin(); 
		 lit != lines_0.end(); ++lit) {
		string t_str = trimString(*lit);

		if (t_str.size() == 0) /* Skip empty lines */
			continue;

		if ( (t_str.size() > 0) && (t_str[0] == commentChar) ) /* Skip commented lines */
			continue;

		lines_1.push_back(t_str);
	}

	// Free previously existing fields of ostTab
	if (pitchShift) {
		free(pitchShift);
		pitchShift = NULL;
	}
	if (intShift) {
		free(intShift);
		intShift = NULL;
	}

	/* 1. Read the time warping section */
	list<string>::const_iterator lit = lines_1.begin();
	vector<string> items = removeComments(splitStringToVector(*lit), commentChar);

	if ( items.size() != 1 )
		throw pcfFileSyntaxError(*lit);
	nTimeWarpAtoms = atoi(items[0].c_str());

	/* Clean up existing time-warp events */
	warpCfg.clear();

	/* Read the time warp details, one by one {tBegin, rate1, dur1, durHold, rate2} */
	/*for (i0 = 0; i0 < nTimeWarpAtoms; i0 ++) {*/ /* TODO */
	int nReadAtoms = 0;
	while (nReadAtoms < nTimeWarpAtoms) {
		items = removeComments(splitStringToVector(*(++lit)), commentChar);

		if (items.size() == 5) {
			t_tBegin = atof(items[0].c_str());
			t_rate1 = atof(items[1].c_str());
			t_dur1 = atof(items[2].c_str());
			t_durHold = atof(items[3].c_str());
			t_rate2 = atof(items[4].c_str());

			try {
				addWarpCfg(t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);
			}
			catch (overlappingWarpIntervalsError) {
				ostringstream oss;
				oss << "ERROR: Overlapping time intervals between time-warping events. " 
					<< "Time-warp event #" << nReadAtoms + 1 << " cannot be loaded.";
				
				mexErrMsgTxt(oss.str().c_str());
			}
			catch (warpCfgInitError) {
				ostringstream oss;
				oss << "ERROR: Failed to initialize time-warping event # " 
					<< nReadAtoms + 1;

				mexErrMsgTxt(oss.str().c_str());
			}
		}
		else if (items.size() == 6) {
			t_ostInitState = static_cast<int>(atof(items[0].c_str()));
			t_tBegin = atof(items[1].c_str());
			t_rate1 = atof(items[2].c_str());
			t_dur1 = atof(items[3].c_str());
			t_durHold = atof(items[4].c_str());
			t_rate2 = atof(items[5].c_str());

			try {
				addWarpCfg(t_ostInitState, t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);
			}
			catch (overlappingWarpIntervalsError) {
				ostringstream oss;
				oss << "ERROR: Overlapping time intervals between time-warping events. " 
					<< "Time-warp event #" << nReadAtoms + 1 << " cannot be loaded.";
				
				mexErrMsgTxt(oss.str().c_str());
			}
			catch (warpCfgInitError) {
				ostringstream oss;
				oss << "ERROR: Failed to initialize time-warping event # " 
					<< nReadAtoms + 1;

				mexErrMsgTxt(oss.str().c_str());
			}
		}
		else {
			throw pcfFileSyntaxError(*lit);
		}
		
		nReadAtoms++;
	}
	
	/* 2. Read the formant/pitch/intensity shift section */
	items = removeComments(splitStringToVector(*(++lit)), commentChar);
	if ( items.size() != 1 )
		throw pcfFileSyntaxError(*lit);
	n = atoi(items[0].c_str());

	if ((pitchShift = (float *)calloc(n, sizeof(float))) == NULL) {
		printf("ERROR: failed to allocate memor for pertCfg.pitchShift\n");
		return;
	}
	if ((intShift = (float *)calloc(n, sizeof(float))) == NULL) {
		printf("ERROR: failed to allocate memor for pertCfg.intShift\n");
		return;
	}
	if ((fmtPertAmp = (float *)calloc(n, sizeof(float))) == NULL) {
		printf("ERROR: failed to allocate memor for pertCfg.fmtPertAmp\n");
		return;
	}
	if ((fmtPertPhi = (float *)calloc(n, sizeof(float))) == NULL) {
		printf("ERROR: failed to allocate memor for pertCfg.fmtPertPhi\n");
		return;
	}

	for (int i0 = 0; i0 < n; i0++) {
		items = removeComments(splitStringToVector(*(++lit)), commentChar);
		if ( items.size() != 5 )
			throw pcfFileSyntaxError(*lit);

		if ( atoi(items[0].c_str()) != i0 )
			throw pcfFileSyntaxError(*lit);

		pitchShift[i0] = static_cast<float>(atof(items[1].c_str()));
		intShift[i0] = static_cast<float>(atof(items[2].c_str()));
		fmtPertAmp[i0] = static_cast<float>(atof(items[3].c_str()));
		fmtPertPhi[i0] = static_cast<float>(atof(items[4].c_str()));
		
		if (bVerbose) {
			printf("\tstat=%d: pitchShift = %f s.t.; intShift = %f dB;\n", 
				   i0, pitchShift[i0], intShift[i0]);
			printf("\t\tfmtPertAmp = %f; fmtPertPhi = %f rad;\n", 
				   fmtPertAmp[i0], fmtPertPhi[i0]);
		}
	}

}


void PERT_CFG::nullify() {
	n = 0;

	if (pitchShift)	{ free(pitchShift);	pitchShift = NULL; }
	if (intShift)	{ free(intShift);	intShift = NULL; }
	if (fmtPertAmp)	{ free(fmtPertAmp);	fmtPertAmp = NULL; }
	if (fmtPertPhi)	{ free(fmtPertPhi);	fmtPertPhi = NULL; }

	warpCfg.clear();
}

const bool PERT_CFG::procTimeWarp(const int stat, const int * statOnsetIndices, 
							      const int nDelay, const double frameDur, 
								  double & t, double & wt) const {
	if (warpCfg.empty())
		return false;

	for (list<pvocWarpAtom>::const_iterator w_it = warpCfg.cbegin(); 
		 w_it != warpCfg.cend(); 
		 ++w_it) {
		int statOnsetIndex = statOnsetIndices[w_it->ostInitState];
		if ( w_it->procTimeWarp(stat, statOnsetIndex, nDelay, frameDur, t, wt) )
			return true;
	}

	return false;
}