/* ost.h
	Online sentence tracking 
	Part of Audapter
	
	Shanqing Cai, 2013
*/

#ifndef OST_H
#define OST_H

#include <string>
#include <map>

class OST_TAB {
private:
	static const char commentChar = '#';

	int stretchCnt;
	int lastStatEnd;
	double stretchSpanAccum;

	int n;	// Number of segments

	int *stat0; // Initial stat number
	int *mode; // Mode number
	double *prm1; // First parameter
	double *prm2; // Second parameter
	double *prm3; // Third parameter

	/* OST heuristic modes */
	typedef enum {
		OST_END = 0, 
		ELAPSED_TIME = 1, 
		INTENSITY_RISE_HOLD = 5, 
		INTENSITY_RISE_HOLD_POS_SLOPE = 6, 
		POS_INTENSITY_SLOPE_STRETCH = 10, 
		NEG_INTENSITY_SLOPE_STRETCH_SPAN = 11, 
		INTENSITY_FALL = 20, 
		INTENSITY_RATIO_RISE = 30, 
		INTENSITY_RATIO_FALL_HOLD = 31
	} OST_MODE_NAME;

	std::map<std::string, int> ostModeMap;

public:
	/* Error classes */
	class ostFileReadingError {};

	class unrecognizedOSTModeError {
	public:
		std::string modeStr;

		unrecognizedOSTModeError(std::string t_modeStr) 
			: modeStr(t_modeStr) {};
	};

	class ostFileSyntaxError {
	public:
		std::string errLine;
		
		ostFileSyntaxError(std::string t_line)
			: errLine(t_line) {};
	};

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
	void readFromFile(const std::string ostFN, const int bVerbose) 
		throw(ostFileReadingError, 
			  unrecognizedOSTModeError, 
			  ostFileSyntaxError);

	/* Nullify */
	void nullify();

	/* Main function: online status tracking */
	int osTrack(const int stat, const int data_counter, const int frame_counter, 
		     	const double rms_s, const double rms_o_slp, const double rms_ratio, const double *rms_rec, 
				const double frameDur);

	/* Reset status */
	void reset();

	/* Destructor */
	~OST_TAB();
};

#endif