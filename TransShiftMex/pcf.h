/* pcf.h
	Perturbation configuration 
	Part of Audapter

	Shanqing Cai, 2013
*/

#ifndef PIP_H
#define PIP_H

#include <string>

class pvocWarpAtom {
public:
	double tBegin;
	double rate1;	// Should aways be in the interval of (0, 1)
	double dur1;
	double durHold;
	double rate2;	// Should aways be > 1
	double dur2;
	int ostInitState; // The ost (online speech tracking) state number at which the clock for pvoc starts. -1: use the same clock as TransShiftMex (default)

	/* Default constructor */
	pvocWarpAtom();
	
	/* Constructor: with input arguments: version 1 */
	pvocWarpAtom(double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2);

	/* Constructor: with input arguments: version 2 */
	pvocWarpAtom(int t_ostInitState, double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2);

	///* Test if the input time t is within the time-shift period: Variant 1: without initial state number */
	//const bool isDuringTimeWarp(const double t) const;

	/* Test if the input time t is within the time-shift period: Variant 2: with initial state number */
	const bool isDuringTimeWarp(const int stat, const int statOnsetIndex, 
								const int nDelay, const double frameDur, 
								double & t, double & wt) const;
};

class PERT_CFG { // Pitch and intensity perturbation configurtion
public:
	int n; // Number of stats
	
	float *pitchShift;	// Unit: st
	float *intShift;	// Unit: dB

	float *fmtPertAmp;
	float *fmtPertPhi;	// Unit: phi

	pvocWarpAtom *warpCfg; // Time warping events

	/* Member functions */
	/* Constructor */
	PERT_CFG();

	/* Destructor */
	~PERT_CFG();

	/* Read PERT_CFG from file */
	void readFromFile(const std::string pertCfgFN, int bVerbose);

	/* Set time-warping configuration */
	void setWarpCfg(double t_tBegin, double t_rate1, double t_dur1, double t_durHold, double t_rate2);


};

#endif