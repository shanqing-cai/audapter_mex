/* ost.cpp
	Online sentence tracking 
	Part of Audapter
	
	Shanqing Cai, 2013
*/


#include "ost.h"

#include <cstdlib>
#include <cstdio>
#include <string>

using namespace std;

/* Constructor */
OST_TAB::OST_TAB() {
	n = 0;
	mode = NULL;
	stat0 = NULL;
	prm1 = NULL;
	prm2 = NULL;
	prm3 = NULL;

	statOnsetIndices = NULL;	

	maxIOICfg.n = 0;
	maxIOICfg.stat0 = NULL;
	maxIOICfg.maxInterval = NULL;
	maxIOICfg.stat1 = NULL;

	stretchCnt = 0;
}

/* Destructor */
OST_TAB::~OST_TAB() {
	if (mode) 
		free(mode);

	if (stat0)
		free(stat0);

	if (prm1)
		free(prm1);

	if (prm2)
		free(prm2);

	if (prm3)
		free(prm3);

	if (statOnsetIndices)
		free(statOnsetIndices);
}

/* Constructor for maxInterOnsetIntervalCfg */
OST_TAB::maxInterOnsetIntervalCfg::maxInterOnsetIntervalCfg() {
	n = 0;

	stat0 = NULL;
	maxInterval = NULL;
	stat1 = NULL;
}

/* Destructor for maxInterOnsetIntervalCfg */
OST_TAB::maxInterOnsetIntervalCfg::~maxInterOnsetIntervalCfg() {
	if (stat0)
		free(stat0);

	if (maxInterval)
		free(maxInterval);

	if (stat1)
		free(stat1);
}

void OST_TAB::reset() {
	stretchCnt = 0;
	lastStatEnd = 0;
	stretchSpanAccum = 0.0;
}

void OST_TAB::readFromFile(const string ostFN, const int bVerbose) {
	FILE *fp;
	int i0, t_n;
	static int maxStatesPerLine = 4;
	char c0[128], c1;

	// Free previously existing fields of ostTab
	if (stat0) {
		free(stat0);
		stat0 = NULL;
	}
	if (mode) {
		free(mode);
		mode = NULL;
	}
	if (prm1) {		
		free(prm1);
		prm1 = NULL;
	}
	if (prm2) {
		free(prm2);
		prm2 = NULL;
	}
	
	if (maxIOICfg.stat0) {
		free(maxIOICfg.stat0);
		maxIOICfg.stat0 = NULL;
	}
	if (maxIOICfg.maxInterval) {
		free(maxIOICfg.maxInterval);
		maxIOICfg.maxInterval = NULL;
	}
	if (maxIOICfg.stat1) {
		free(maxIOICfg.stat1);
		maxIOICfg.stat1 = NULL;
	}

	if (statOnsetIndices) {
		free(statOnsetIndices);
		statOnsetIndices = NULL;
	}

	if (bVerbose)
		printf("ostFN = %s\n", ostFN.c_str());

	/* fp = fopen(this->ostFN, "r"); */
	if (fopen_s(&fp, ostFN.c_str(), "r")) {
		printf("ERROR: Unable to open ost file: %s\n", ostFN.c_str());
		return;
	}

	for (i0 = 0; i0 < 3; i0++)
		/* fscanf(fp, "%s", c0); */
		fscanf_s(fp, "%s", c0, sizeof(c0));
	rmsSlopeWin = atof(c0);
	
	if (bVerbose)
		printf("rmsSlopeWin = %f\n", rmsSlopeWin);

	for (i0 = 0; i0 < 3; i0++)
		/* fscanf(fp, "%s", c0); */
		fscanf_s(fp, "%s", c0, sizeof(c0));
	n = atoi(c0);

	if (bVerbose)
		printf("ostTab.n = %d\n", n);

	if ((stat0 = (int *)calloc(n, sizeof(int))) == NULL) {
		printf("ERROR: failed to allocate memor for ostTab.stat0");
		return;
	}
	if ((mode = (int *)calloc(n, sizeof(int))) == NULL) {
		printf("ERROR: failed to allocate memor for ostTab.mode");
		return;
	}
	if ((prm1 = (double *)calloc(n, sizeof(double))) == NULL) {
		printf("ERROR: failed to allocate memor for ostTab.prm1");
		return;
	}
	if ((prm2 = (double *)calloc(n, sizeof(double))) == NULL) {
		printf("ERROR: failed to allocate memor for ostTab.prm2");
		return;
	}

	for (i0 = 0; i0 < n; i0++) {
		/* fscanf(fp, "%s", c0); */
		fscanf_s(fp, "%s", c0, sizeof(c0));
		stat0[i0] = atoi(c0);
		//printf("\tstat0[%d] = %d\n", i0, stat0[i0]);

		/*fscanf(fp, "%s", c0);*/
		fscanf_s(fp, "%s", c0, sizeof(c0));
		mode[i0] = atoi(c0);
		//printf("\tmode[%d] = %d\n", i0, mode[i0]);

		/*fscanf(fp, "%s", c0);*/
		fscanf_s(fp, "%s", c0, sizeof(c0));
		prm1[i0] = atof(c0);
		//printf("\tprm1[%d] = %f\n", i0, prm1[i0]);

		/* fscanf(fp, "%s", c0); */
		fscanf_s(fp, "%s", c0, sizeof(c0));
		prm2[i0] = atof(c0);
		//printf("\tprm2[%d] = %f\n", i0, prm2[i0]);

		c1 = '\0';
		while (!(c1 == '\n' || c1 == '\r' || c1 == EOF)) {
			c1 = fgetc(fp);
		}
	
		if (bVerbose)
			printf("\tSeg %d: stat0=%d; mode=%d; prm1=%f; prm2=%f\n", 
				   i0, stat0[i0], mode[i0], prm1[i0], prm2[i0]);
	}

	if ((statOnsetIndices = (int *) calloc(n * maxStatesPerLine, sizeof(int))) == NULL) {
		printf("ERROR: failed to allocate memor for statOnsetIndices");
		return;
	}

	if (c1 == EOF) {
		fclose(fp);
		return;
	}
			
	while (c1 != EOF) {
		c1 = fgetc(fp);

		if (c1 == EOF) {
			fclose(fp);
			return;
		}

		if (c1 == 'n')
			break;
	}

	/* Process maxInterOnsetIntervalCfg (maxIOICfg)  */
	for (i0 = 0; i0 < 2; i0++)
		/* n = fscanf(fp, "%s", c0); */
		t_n = fscanf_s(fp, "%s", c0, sizeof(c0));

	if (atoi(c0) > 0) {
		maxIOICfg.n = atoi(c0);
		if (bVerbose)
			printf("ostTab.maxIOICfg.n = %d\n", maxIOICfg.n);

		if ((maxIOICfg.stat0 = (int *)calloc(maxIOICfg.n, sizeof(int))) == NULL) {
			printf("ERROR: failed to allocate memor for ostTab.maxIOICfg.stat0");
			fclose(fp);
			return;
		}
		if ((maxIOICfg.maxInterval = (double *)calloc(maxIOICfg.n, sizeof(double))) == NULL) {
			printf("ERROR: failed to allocate memor for ostTab.maxIOICfg.maxInterval");
			fclose(fp);
			return;
		}
		if ((maxIOICfg.stat1 = (int *)calloc(maxIOICfg.n, sizeof(int))) == NULL) {
			printf("ERROR: failed to allocate memor for ostTab.maxIOICfg.stat1");
			fclose(fp);
			return;
		}

		for (i0 = 0; i0 < maxIOICfg.n; i0++) {
			/* fscanf(fp, "%s", c0); */
			fscanf_s(fp, "%s", c0, sizeof(c0));
			maxIOICfg.stat0[i0] = atoi(c0);

			/* fscanf(fp, "%s", c0); */
			fscanf_s(fp, "%s", c0, sizeof(c0));
			maxIOICfg.maxInterval[i0] = atof(c0);

			/* fscanf(fp, "%s", c0); */
			fscanf_s(fp, "%s", c0, sizeof(c0));
			maxIOICfg.stat1[i0] = atoi(c0);

			c1 = '\0';
			while (!(c1 == '\n' || c1 == '\r' || c1 == EOF)) {
				c1 = fgetc(fp);
			}

			if (bVerbose)
				printf("maxIOICfg %d: stat0=%d; maxInterval=%d; stat1=%d\n", 
					   i0, maxIOICfg.stat0[i0], maxIOICfg.maxInterval[i0], maxIOICfg.stat1[i0]);
		}
	}
	else {
		fclose(fp);
		return;
	}

	fclose(fp);
}

int OST_TAB::osTrack(const int stat, const int data_counter, const int frame_counter, 
					 const double rms_o, const double rms_o_slp, const double rms_ratio, const double *rms_rec, 
					 const double frameDur) {
/* Input: stat: current status number */
	int k, i, j;
	int t_stat0, t_mode;
	int bIsGOET;
	int nLBDelay;
	int minDurN;

	int stat_out = stat;
	/*double elapTime;*/

	const double rmsLBDelay = 5 * 0.002;

	nLBDelay = (int)(floor(rmsLBDelay / frameDur + 0.5));

	if (n == 0) {
		return 0;
	}

	// Determine which segment of the tab (ostCfg) we are currently in
	k = -1;
	for (i = 0; i < n; i++) {
		if (stat >= stat0[i] && stat < stat0[i + 1]) {
			k = i;
			break;
		}
	}

	if (k >= 0) {
		t_stat0 = stat0[k];
		t_mode = mode[k];

		if (t_mode == 1) { // (+1) Elapsed time from previous state. prm1: duration (s)
			if ((data_counter - statOnsetIndices[stat]) * frameDur > prm1[k]) {
				stat_out = stat + 1;
				statOnsetIndices[stat] = frame_counter;
			}
		}
		else if (t_mode == 5) { // (+2) Crossing an rmsThresh (from below) and hold. prm1: rmsThresh; prm2: minHoldDur (s)			
			if (stat == t_stat0) {
				if (rms_o > prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					stretchCnt = 1;
				}
			}
			else {
				minDurN = (int) floor(prm2[k] / frameDur + 0.5);

				if (rms_o > prm1[k]) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}

		}
		else if (t_mode == 6) { // (+2) Crossing an rmsThresh (from below) and hold, during positive RMS slopes. prm1: rmsThresh; prm2: minHoldDur (s)			
			if (stat == t_stat0) {
				if (rms_o > prm1[k] && 
					rms_o_slp > 0) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					stretchCnt = 1;
				}
			}
			else {
				minDurN = (int) floor(prm2[k] / frameDur + 0.5);

				if (rms_o > prm1[k] &&
					rms_o_slp > 0) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}

		}
		else if (t_mode == 10) { // (+2) Stretch of of positive rms slope, with only a stretch count threshold
			if (stat == t_stat0) {
				if (rms_o_slp > 0) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					stretchCnt = 1;
				}
			}
			else {
				if (rms_o_slp > 0) {
					stretchCnt++;
					if (stretchCnt > prm1[k]) {
						stat_out = stat + 1;
						statOnsetIndices[stat] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}

		}
		else if (t_mode == 11) { // Stretch of negative rms slope, with a stretch count thresold and a stretch span threshold
			if (stat == t_stat0) {
				if (rms_o_slp < 0) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					stretchCnt = 1;
					stretchSpanAccum = rms_o_slp;
				}
			}
			else {
				if (rms_o_slp < 0) {
					stretchCnt++;
					stretchSpanAccum += rms_o_slp;

					if ((stretchCnt > prm1[k]) && (stretchSpanAccum < prm2[k])) {
						stat_out = stat + 1;
						statOnsetIndices[stat] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat + 1;
				}
			}

		}
		else if (t_mode == 20) { // (+1) Fall from a certain RMS threshold
			minDurN = static_cast<int>(floor(prm2[k] / frameDur + 0.5));

			bIsGOET = 0;
			for (j = 0; j < nLBDelay; j++) {
				if (data_counter - j < 0) {
					bIsGOET = 1;
					break;
				}
				if (rms_rec[data_counter - j] >= prm1[k]) {
					bIsGOET = 1;
					break;
				}
			}

			if ((bIsGOET == 0) && ((data_counter - lastStatEnd) > minDurN)) {
				stat_out = stat + 1;
				statOnsetIndices[stat] = frame_counter;
				lastStatEnd = data_counter;
			}

		}
		else if (t_mode == 30) { // (+2) RMS ratio cross, hold and fall. prm1: rms_ratio threshold; prm2: minDurN
			if (stat == t_stat0) {
				if (1. / rms_ratio > prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					stretchCnt = 0;
				}
			}
			else if (stat - t_stat0 == 1) {
				minDurN = (int)floor(prm2[k] / frameDur);

				if (1. / rms_ratio > prm1[k]) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat] = frame_counter;
					}
				}
				else {
					stat_out = stat + 1;
				}
			}
			else {
				if (1. / rms_ratio < prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					lastStatEnd = data_counter;
				}
			}			
		}
		else if (t_mode == 31) { // (+2) RMS ratio fall from a threshold, hold and fall. prm1: rms_ratio threshold; prm2: minDurN
			if (stat == t_stat0) {
				if (1. / rms_ratio < prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					stretchCnt = 0;
				}
			}
			else if (stat - t_stat0 == 1) {
				minDurN = (int)floor(prm2[k] / frameDur);

				if (1. / rms_ratio < prm1[k]) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat] = frame_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}
			else {
				if (1. / rms_ratio < prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat] = frame_counter;
					lastStatEnd = data_counter;
				}
			}
		}
			 
	}

	/* Process maxIOICfg */
	if (maxIOICfg.n > 0) {
		for (i = 0; i < maxIOICfg.n; i++) {
			if ((stat >= maxIOICfg.stat0[i]) && (stat < maxIOICfg.stat1[i])) {
				if ((frame_counter - statOnsetIndices[maxIOICfg.stat0[i]]) * frameDur > maxIOICfg.maxInterval[i]) {
					for (j = stat + 1; j <= maxIOICfg.stat1[i]; j++) {
						statOnsetIndices[stat] = frame_counter;
					}

					stat_out = maxIOICfg.stat1[i];
					
				}
			}
		}
	}

	return stat_out;
}