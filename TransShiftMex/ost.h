/* ost.h
	Online sentence tracking 
	Part of Audapter
	
	Shanqing Cai, 2013
*/

#ifndef OST_H
#define OST_H

#include <string>

class OST_TAB {
private:
	int stretchCnt;
	int lastStatEnd;
	double stretchSpanAccum;

	int n;	// Number of segments

	int *stat0; // Initial stat number
	int *mode; // Mode number
	double *prm1; // First parameter
	double *prm2; // Second parameter
	double *prm3; // Third parameter

public:
	double rmsSlopeWin;
	int *statOnsetIndices; /* Onset index of all state numbers. By definition, statOnsets[0] = 0 */

	class maxInterOnsetIntervalCfg {
	public:
		int n;
	
		int *stat0;
		double *maxInterval; /* Unit: s */
		int *stat1;

		/* Constructor */
		maxInterOnsetIntervalCfg();

		/* Destructor */
		~maxInterOnsetIntervalCfg();
	};

	maxInterOnsetIntervalCfg maxIOICfg;

	/* Member functions */
	/* Constructor */
	OST_TAB();

	/* Reading from a file */
	void readFromFile(const std::string ostFN, const int bVerbose);

	/* Main function: online status tracking */
	int osTrack(const int stat, const int data_counter, const int frame_counter, 
		     	const double rms_o, const double rms_o_slp, const double rms_ratio, const double *rms_rec, 
				const double frameDur);

	/* Reset status */
	void reset();

	/* Destructor */
	~OST_TAB();
};

#endif