/* 
	lpc_formant.h
	Linear prediction and formant tracking 

	Shanqing Cai 12/2013

*/

#ifndef LPC_FORMANT_H
#define LPC_FORMANT_H

typedef double dtype;

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

/* Right shift */
#define RSL(INTEGER,SHIFT) (int)( ( (unsigned)INTEGER ) >> SHIFT )

struct CepstralPitchTrackerConfig {
	bool activated;
	dtype pitchLowerBoundHz;
	dtype pitchUpperBoundHz;

	CepstralPitchTrackerConfig() :
		activated(false), pitchLowerBoundHz(0.0), pitchUpperBoundHz(0.0) {}

	CepstralPitchTrackerConfig(
		const bool activated,
		const dtype pitchLowerBoundHz,
		const dtype pitchUpperBoundHz) :
		activated(activated), pitchLowerBoundHz(pitchLowerBoundHz),
		pitchUpperBoundHz(pitchUpperBoundHz) {}
};

class LPFormantTracker {
	//TODO:
	//		resetting properly
	//		confirm temp_fame and R size issue
	//		incorporate trackPhi()
	//		below-RMS-thresh reset for realRoots, etc.
	//		Constructor t_nFFT power of 2 check
	//		Check that avgLen is shorter than maxAvgLen (i.e., maxPitchLen)

private:
	/* Constants */
	static const int maxNLPC = 20;
	static const int maxNPoles = maxNLPC / 2 + 2;
	static const int maxFmtTrackJump = maxNPoles;
	static const int maxNTracks = 5;	
	static const int maxAvgLen = 100;
	
	/* Linear prediction parameters */
	int nLPC;		/* Order of LP */
	int sr;			/* Sampling rate */
	int bufferSize; /* Input buffer size (# of samples) */
	//int winSize;	/* Window length */
	int nFFT;		/* Size of FFT frame */

	/* Parameters related to cepstral liftering */
	bool bCepsLift;
	int cepsWinWidth;	

	/* DP algorithm parameters */
	int nTracks;
	int nCands;
	dtype aFact, bFact, gFact;
	dtype fn1, fn2;

	/* Pitch tracking */
	int bTrackPitch;
	dtype pitchLowerBoundHz;
	dtype pitchUpperBoundHz;

	/* Moving weighted average (MWA) parameters */
	bool bMWA;
	int avgLen;	

	dtype trackFF;

	int nLPC_SQR;

	dtype * winFunc; /* Window function */

	/* Working date fields for FFT (cepstral operations) */
	dtype * ftBuf1;
	dtype * ftBuf2;
	dtype * fftc;

	/* Working data fields for LP */
	dtype * temp_frame;
	dtype * R;

	dtype * Acompanion;			// companion matrix for (eigenvalue ---> )roots calculation //Marked
	dtype * AHess;

	/* Intermediate data fields */
	dtype * realRoots;
	dtype * imagRoots;

	dtype * cumMat;
	dtype * costMat;

	dtype * weiMatPhi;
	dtype * weiMatBw;
	dtype * weiVec;
	dtype * sumWeiPhi;
	dtype * sumWeiBw;
	dtype sumWei;

	dtype * radius_us;		/* Unsmoothed pole radii */
	dtype * phi_us;			/* Unsmoothed pole angles */
	dtype * bandwidth_us;	/* Unsmoothed pole bandwidths */
	//dtype * phi_s;			/* Smoothed pole angles */

	dtype latestPitchHz;  /* Latest-tracked pitch value, in Hz. */

	int mwaCtr;		/* Moving average counter */
	int mwaCircCtr; /* Moving average circular counter */
	
	/* Private member functions */
	/* LPC and formant tracking subroutines */
	/* Levinson-Durbin recursion */
	void levinson(dtype * R, dtype * aa, const int size);

	/* Solve for the roots of a polynomial */
	int hqr_roots(dtype * c, dtype * wr, dtype * wi);

	void getAi(dtype * xx, dtype * aa);

	/* Get the angle (Phi) and magnitude (Bw) of the roots  */
	void getRPhiBw(dtype * wr, dtype * wi, 
				   dtype * radius, dtype * phi, dtype * bandwith);

	/* DP formant tracking */
	void trackPhi(dtype * r_ptr, dtype * phi_ptr);

	/* Weighted moving average (mwa) of phi and R */
	int mwa(dtype * phi_ptr, dtype * bw_ptr , dtype * wmaPhi_ptr);

	void genHanningWindow();

public:
	dtype * lpcAi;

	/* Constructor */
	LPFormantTracker(const int t_nLPC, const int t_sr, const int t_bufferSize, 					 
					 const int t_nFFT, const int t_cepsWinWidth, 
					 const int t_nTracks, 
					 const dtype t_aFact, const dtype t_bFact, const dtype t_gFact, 
					 const dtype t_fn1, const dtype t_fn2, 
					 const bool t_bMWA, const int t_avgLen,
		             const CepstralPitchTrackerConfig& pitchTrackerConfig);

	/* Destructor */
	~LPFormantTracker();

	/* Error classes */
	class initializationError {};
	class nLPCTooLargeError {};

	/* Member functions */
	
	/* Reset after a supra-threhsold interval */
	void postSupraThreshReset();

	/* Reset status */
	void reset();	

	/* LPC */
	//void procFrame(dtype * xx, dtype st_rms, dtype * radius, dtype * phi, dtype * bandwidth);
	void procFrame(dtype * xx, dtype st_rms, 
				   dtype * radius, dtype * phi, dtype * bandwidth, 
				   dtype * fmts);

	/* Setters and getters (inline) */
	void setNLPC(const int t_nLPC);
	void setBufferSize(const int t_bufferSize);
	void setNFFT(const int t_nFFT);	

	void setBCepsLift(const bool t_bCepsLift);
	void setCepsWinWidth(const int t_cepsWinWidth);

	const int getNLPC() const;
	
	const dtype getLatestPitchHz() const; 
};

#endif