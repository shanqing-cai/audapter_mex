/* ost.cpp
	Online sentence tracking 
	Part of Audapter
	
	Shanqing Cai, 2013
*/

#include <cstdio>
#include <string>
#include <list>
#include <vector>

#include "ost.h"
#include "utils.h"

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

	/* Create map for OST modes */
	ostModeMap[string("OST_END")] = OST_END;
	ostModeMap[string("ELAPSED_TIME")] = ELAPSED_TIME;
	ostModeMap[string("INTENSITY_RISE_HOLD")] = INTENSITY_RISE_HOLD;
	ostModeMap[string("INTENSITY_RISE_HOLD_POS_SLOPE")] = INTENSITY_RISE_HOLD_POS_SLOPE;
	ostModeMap[string("POS_INTENSITY_SLOPE_STRETCH")] = POS_INTENSITY_SLOPE_STRETCH;
	ostModeMap[string("NEG_INTENSITY_SLOPE_STRETCH_SPAN")] = NEG_INTENSITY_SLOPE_STRETCH_SPAN;
	ostModeMap[string("INTENSITY_FALL")] = INTENSITY_FALL;
	ostModeMap[string("INTENSITY_RATIO_RISE")] = INTENSITY_RATIO_RISE;
	ostModeMap[string("INTENSITY_RATIO_FALL_HOLD")] = INTENSITY_RATIO_FALL_HOLD;
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



void OST_TAB::readFromFile(const string ostFN, const int bVerbose) 
	throw(ostFileReadingError, 
	      unrecognizedOSTModeError, 
		  ostFileSyntaxError) 
{
	//FILE *fp;
	int i0;
	const int maxStatesPerLine = 4;

	list<string> ostLines_0 = readLinesFromFile(ostFN);
	if (ostLines_0.empty()) {
		throw ostFileReadingError();
	}
		

	/* Trim lines; remove empty lines; remove commented lines */
	list<string> ostLines_1;
	for (list<string>::const_iterator lit = ostLines_0.begin(); 
		 lit != ostLines_0.end(); ++lit) {
		string t_str = trimString(*lit);

		if (t_str.size() == 0) /* Skip empty lines */
			continue;

		if ( (t_str.size() > 0) && (t_str[0] == '#') ) /* Skip commented lines */
			continue;

		ostLines_1.push_back(t_str);
	}

	/* Iterate through the lines */
	/*for (list<string>::const_iterator lit = ostLines_1.begin(); 
		 lit != ostLines_1.end(); ++lit) {
		list<string> items = splitStringToList(*lit);
	}*/

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

	list<string>::const_iterator lit = ostLines_1.begin();
	vector<string> items = splitStringToVector(*lit);
	if ( (items.size() != 3) || 
		 (items[0] != string("rmsSlopeWin")) || 
		 (items[1] != string("=")) )
		throw ostFileSyntaxError(*lit);

	rmsSlopeWin = atof(items[2].c_str());

	if (bVerbose)
		printf("rmsSlopeWin = %f\n", rmsSlopeWin);

	items = splitStringToVector(*(++lit));
	if ( (items.size() != 3) || 
		 (items[0] != string("n")) || 
		 (items[1] != string("=")) )
		throw ostFileSyntaxError(*lit);

	n = atoi(items[2].c_str());

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
		items = splitStringToVector(*(++lit));
		if ( (items.size() != 5) || 
			 (items[4] != string("{}")) )
			throw ostFileSyntaxError(*lit);
		stat0[i0] = atoi(items[0].c_str());


		mode[i0] = atoi(items[1].c_str());
		if ( (mode[i0] == 0) && (items[1] != string("0")) && (items[1] != string("OST_END")) ) {
			/* Test if the string is in the ostModeMap */
			try {
				mode[i0] = ostModeMap.at(items[1]); /* WARNING: Assume C++11 is available */
			}
			catch (out_of_range) {
				//fclose(fp);
				throw unrecognizedOSTModeError(items[1]);
			}
		}

		prm1[i0] = atof(items[2].c_str());
		prm2[i0] = atof(items[3].c_str());
	
		if (bVerbose)
			printf("\tSeg %d: stat0=%d; mode=%d; prm1=%f; prm2=%f\n", 
				   i0, stat0[i0], mode[i0], prm1[i0], prm2[i0]);
	}

	if ((statOnsetIndices = (int *) calloc(n * maxStatesPerLine, sizeof(int))) == NULL) {
		printf("ERROR: failed to allocate memor for statOnsetIndices");
		return;
	}

	if (++lit == ostLines_1.end())
		return;

	items = splitStringToVector(*lit);
	if ( (items.size() != 3) || 
		 (items[0] != string("n")) || 
		 (items[1] != string("=")) )
		throw ostFileSyntaxError(*lit);

	/* Process maxInterOnsetIntervalCfg (maxIOICfg)  */
	maxIOICfg.n = atoi(items[2].c_str());

	if (maxIOICfg.n > 0) {
		if (bVerbose)
			printf("ostTab.maxIOICfg.n = %d\n", maxIOICfg.n);

		if ((maxIOICfg.stat0 = (int *)calloc(maxIOICfg.n, sizeof(int))) == NULL) {
			printf("ERROR: failed to allocate memor for ostTab.maxIOICfg.stat0");
			return;
		}
		if ((maxIOICfg.maxInterval = (double *)calloc(maxIOICfg.n, sizeof(double))) == NULL) {
			printf("ERROR: failed to allocate memor for ostTab.maxIOICfg.maxInterval");
			return;
		}
		if ((maxIOICfg.stat1 = (int *)calloc(maxIOICfg.n, sizeof(int))) == NULL) {
			printf("ERROR: failed to allocate memor for ostTab.maxIOICfg.stat1");
			return;
		}

		for (i0 = 0; i0 < maxIOICfg.n; i0++) {
			items = splitStringToVector(*(++lit));

			if ( items.size() != 3 )
				throw ostFileSyntaxError(*lit);

			maxIOICfg.stat0[i0] = atoi(items[0].c_str());
			maxIOICfg.maxInterval[i0] = atof(items[1].c_str());
			maxIOICfg.stat1[i0] = atoi(items[2].c_str());

			if (bVerbose)
				printf("maxIOICfg %d: stat0=%d; maxInterval=%d; stat1=%d\n", 
					   i0, maxIOICfg.stat0[i0], maxIOICfg.maxInterval[i0], maxIOICfg.stat1[i0]);
		}
	}
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

		if (t_mode == ELAPSED_TIME) { // (+1) Elapsed time from previous state. prm1: duration (s)
			if ((data_counter - statOnsetIndices[stat]) * frameDur > prm1[k]) {
				stat_out = stat + 1;
				statOnsetIndices[stat_out] = frame_counter;
			}
		}
		else if (t_mode == INTENSITY_RISE_HOLD) { // (+2) Crossing an rmsThresh (from below) and hold. prm1: rmsThresh; prm2: minHoldDur (s)			
			if (stat == t_stat0) {
				if (rms_o > prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
					stretchCnt = 1;
				}
			}
			else {
				minDurN = (int) floor(prm2[k] / frameDur + 0.5);

				if (rms_o > prm1[k]) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat_out] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}

		}
		else if (t_mode == INTENSITY_RISE_HOLD_POS_SLOPE) { // (+2) Crossing an rmsThresh (from below) and hold, during positive RMS slopes. prm1: rmsThresh; prm2: minHoldDur (s)			
			if (stat == t_stat0) {
				if (rms_o > prm1[k] && 
					rms_o_slp > 0) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
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
						statOnsetIndices[stat_out] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}

		}
		else if (t_mode == POS_INTENSITY_SLOPE_STRETCH) { // (+2) Stretch of of positive rms slope, with only a stretch count threshold
			if (stat == t_stat0) {
				if (rms_o_slp > 0) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
					stretchCnt = 1;
				}
			}
			else {
				if (rms_o_slp > 0) {
					stretchCnt++;
					if (stretchCnt > prm1[k]) {
						stat_out = stat + 1;
						statOnsetIndices[stat_out] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}

		}
		else if (t_mode == NEG_INTENSITY_SLOPE_STRETCH_SPAN) { // Stretch of negative rms slope, with a stretch count thresold and a stretch span threshold
			if (stat == t_stat0) {
				if (rms_o_slp < 0) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
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
						statOnsetIndices[stat_out] = frame_counter;
						lastStatEnd = data_counter;
					}
				}
				else {
					stat_out = stat + 1;
				}
			}

		}
		else if (t_mode == INTENSITY_FALL) { // (+1) Fall from a certain RMS threshold
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
				statOnsetIndices[stat_out] = frame_counter;
				lastStatEnd = data_counter;
			}

		}
		else if (t_mode == INTENSITY_RATIO_RISE) { // (+2) RMS ratio cross, hold and fall. prm1: rms_ratio threshold; prm2: minDurN
			if (stat == t_stat0) {
				if (1. / rms_ratio > prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
					stretchCnt = 0;
				}
			}
			else if (stat - t_stat0 == 1) {
				minDurN = (int)floor(prm2[k] / frameDur);

				if (1. / rms_ratio > prm1[k]) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat_out] = frame_counter;
					}
				}
				else {
					stat_out = stat + 1;
				}
			}
			else {
				if (1. / rms_ratio < prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
					lastStatEnd = data_counter;
				}
			}			
		}
		else if (t_mode == INTENSITY_RATIO_FALL_HOLD) { // (+2) RMS ratio fall from a threshold, hold and fall. prm1: rms_ratio threshold; prm2: minDurN
			if (stat == t_stat0) {
				if (1. / rms_ratio < prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
					stretchCnt = 0;
				}
			}
			else if (stat - t_stat0 == 1) {
				minDurN = (int)floor(prm2[k] / frameDur);

				if (1. / rms_ratio < prm1[k]) {
					stretchCnt++;
					if (stretchCnt > minDurN) {
						stat_out = stat + 1;
						statOnsetIndices[stat_out] = frame_counter;
					}
				}
				else {
					stat_out = stat - 1;
				}
			}
			else {
				if (1. / rms_ratio < prm1[k]) {
					stat_out = stat + 1;
					statOnsetIndices[stat_out] = frame_counter;
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