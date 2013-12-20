/* pcf.cpp
	Perturbation configuration 
	Part of Audapter

	Shanqing Cai, 2013
*/

#include <cstdlib>
#include <cstdio>
#include <string>

#include "pcf.h"

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
}

/* pvocWarpAtom: with input arguments: Version 1 */
pvocWarpAtom::pvocWarpAtom(double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2) {
	tBegin = t_tBegin;
	rate1 = t_rate1;
	dur1 = t_dur1;
	durHold = t_durHold;
	rate2 = t_rate2;

	ostInitState = 0;

	dur2 = (1 - rate1) / (rate2 - 1) * dur1;
}

/* pvocWarpAtom: with input arguments: Version 2 */
pvocWarpAtom::pvocWarpAtom(int t_ostInitState, double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2) {
	tBegin = t_tBegin;
	rate1 = t_rate1;
	dur1 = t_dur1;
	durHold = t_durHold;
	rate2 = t_rate2;

	ostInitState = t_ostInitState;

	dur2 = (1 - rate1) / (rate2 - 1) * dur1;
}

/* pvocWarpAtom: Test if the input time t is within the time-shift period: Variant 2: with initial state number */
const bool pvocWarpAtom::isDuringTimeWarp(const int stat, const int statOnsetIndex, 
										  const int nDelay, const double frameDur, 
										  double & t, double & wt) const {
/* Input arguments: 
		stat:			current OST status number
		statOnsetIndex: frame index at which the current status is first entered
		nDelay:			Audapter nDelay
		frameDur:		Audapter frame duration (in s, = frameLen / sr)
		t:				current time (s)
		wt:				warped time (s)
*/

	if ((ostInitState < 0) || (stat < ostInitState)) {
		return false;
	}
	else {
		double t01;
		bool duringTimeWarp;

		t01 = t - static_cast<double>(statOnsetIndex) * frameDur;
		//t01 = t - static_cast<double>(statOnsetIndex - (nDelay - 1)) * frameDur;
		duringTimeWarp = (t01 >= tBegin) && (t01 < tBegin + dur1 + durHold + dur2);

		if (duringTimeWarp)
			t = t01;
			//t = t;

		/* Determine warped time */
		if (t < tBegin + dur1){ /* Time dilation (deceleration) */
			wt = (t - tBegin) * rate1 + tBegin;
		}
		else if (t < tBegin + dur1 + durHold){ /* Time shifting (no compression or dilation) */
			wt = rate1 * dur1 - dur1 + t;
		}
		else if (t < tBegin + dur1 + durHold + dur2){ /* Time compression (acceleration) at the end of the warp interval */
			//t1 = (t0 - (warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold)) * warpCfg->rate2 + warpCfg->tBegin + warpCfg->dur1 + warpCfg->durHold;
			wt = tBegin + dur1 + durHold + dur2;
			wt -= (tBegin + dur1 + durHold + dur2 - t) * rate2;
		}
			
		//if (pertCfg.warpCfg->ostInitState >= 0)
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

	warpCfg = new pvocWarpAtom();
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

	if (warpCfg)
		delete warpCfg;
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
int deblank(char *line) {
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

/* Set time-warping configuration */
void PERT_CFG::setWarpCfg(double t_tBegin, double t_rate1, 
						  double t_dur1, double t_durHold, double t_rate2) {
	if (warpCfg) {
		delete warpCfg;
		warpCfg = NULL;
	}

	warpCfg = new pvocWarpAtom(t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);
}


/* Read from file */
void PERT_CFG::readFromFile(const string pertCfgFN, const int bVerbose) {
	FILE *fp;
	int i0;
	int lineWidth;
	int nTimeWarpAtoms = -1; /* SCai: currently a dummy variable. TODO: implement multiple time warping atoms */
	char line[512];
	double tmpx[16];
	int t_ostInitState;
	double t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2;

	// Free previously existing fields of ostTab
	if (pitchShift) {
		free(pitchShift);
		pitchShift = NULL;
	}
	if (intShift) {
		free(intShift);
		intShift = NULL;
	}

	/*fp = fopen(pertCfgFN, "r");*/
	if (fopen_s(&fp, pertCfgFN.c_str(), "r")) {
		printf("ERROR: Unable to open ost file: %s\n", pertCfgFN.c_str());
		return;
	}

	/* 1. Read the time warping section */
	while (nTimeWarpAtoms == -1) {
		lineWidth = readline(fp, line);
		lineWidth = deblank(line);

		if (lineWidth == 0) 
			continue;
		if (line[0] == '#') 
			continue;

		nTimeWarpAtoms = atoi(line);
	}

	/* Read the time warp details, one by one {tBegin, rate1, dur1, durHold, rate2} */
	/*for (i0 = 0; i0 < nTimeWarpAtoms; i0 ++) {*/ /* TODO */
	for (i0 = 0; i0 < 1 && i0 < nTimeWarpAtoms; i0 ++) {
		lineWidth = readline(fp, line);
		lineWidth = deblank(line);

		if (lineWidth == 0) 
			continue;
		if (line[0] == '#') 
			continue;

		/*sscanf(line, "%f, %f, %f, %f, %f", 
			   &t_tBegin, &t_rate1, &t_dur1, &t_durHold, &t_rate2);*/

		if (string_count_char(line, ',') == 4) {
			sscanf_floatArray(line, tmpx, 5);

			t_tBegin = tmpx[0];
			t_rate1 = tmpx[1];
			t_dur1 = tmpx[2];
			t_durHold = tmpx[3];
			t_rate2 = tmpx[4];

			/* TODO: Implement multiple time warping events */
			if (warpCfg)
				delete warpCfg;
			warpCfg = new pvocWarpAtom(t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);		
			
		}
		else {
			sscanf_floatArray(line, tmpx, 6);

			t_ostInitState = (int) tmpx[0];
			t_tBegin = tmpx[1];
			t_rate1 = tmpx[2];
			t_dur1 = tmpx[3];
			t_durHold = tmpx[4];
			t_rate2 = tmpx[5];

			warpCfg = new pvocWarpAtom(t_ostInitState, t_tBegin, t_rate1, t_dur1, t_durHold, t_rate2);		
		}
		
	}
	
	/* 2. Read the formant/pitch/intensity shift section */
	n = -1;

	while (n == -1) {
		lineWidth = readline(fp, line);
		lineWidth = deblank(line);

		if (lineWidth == 0) 
			continue;
		if (line[0] == '#') 
			continue;

		n = atoi(line);
		
		if (bVerbose)
			printf("pertCfg.n = %d\n", n);
	}

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

	for (i0 = 0; i0 < n; i0++) {
		lineWidth = readline(fp, line);
		lineWidth = deblank(line);

		if (lineWidth == 0) 
			continue;
		if (line[0] == '#') 
			continue;

		sscanf_floatArray(line, tmpx, 5);
		pitchShift[i0] = (float) tmpx[1];
		intShift[i0] = (float) tmpx[2];
		fmtPertAmp[i0] = (float) tmpx[3];
		fmtPertPhi[i0] = (float) tmpx[4];
		
		if (bVerbose) {
			printf("\tstat=%d: pitchShift = %f s.t.; intShift = %f dB;\n", 
				   i0, pitchShift[i0], intShift[i0]);
			printf("\t\tfmtPertAmp = %f; fmtPertPhi = %f rad;\n", 
				   fmtPertAmp[i0], fmtPertPhi[i0]);
		}
	}

	fclose(fp);
}