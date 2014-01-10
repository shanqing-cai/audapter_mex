#ifndef PHASE_VOCODER_H
#define PHASE_VOCODER_H

#define M_PI       3.14159265358979323846

typedef double dtype;

class PhaseVocoder {
	/* TODOs: */

public:
	typedef enum {
		PITCH_SHIFT_ONLY,
		TIME_WARP_ONLY, 
		TIME_WARP_WITH_FIXED_PITCH_SHIFT, 
	} operMode; /* Mode of operation */

	operMode mode;

	/* Error classes */
	class initializationError {};
	class timeWarpFuturePredError {};
	class fixedPitchShiftNotSpecifiedErr {};

	/* Default constructor */
	PhaseVocoder();

	/* Constructor */
	PhaseVocoder(operMode t_operMode, 
				 const int t_nDelay, 
				 const dtype t_sr, 
				 const int t_frameLen, 
				 const int t_pvocFrameLen, 
				 const int t_pvocHop);	

	/* Destructor */
	~PhaseVocoder();

	/* Memory clean up */
	void cleanup();

	/* Configure parameters */
	void config(operMode t_operMode, 
				const int t_nDelay, 
				const dtype t_sr, 
				const int t_frameLen, 
				const int t_pvocFrameLen, 
				const int t_pvocHop);

	/* Main interface function: process a frame */
	void procFrame(const dtype * inBuf, const dtype shift);

	/* Status reset */
	void reset();

	/* Getters */
	const operMode getMode() const;

	/* Setters */
	void setFixedPitchShiftST(const dtype t_fixedPitchShiftST);

private:
	/* Parameters */
	int nDelay;			/* Amount of global delay (# of frames) */

	dtype sr;			  /* Sampling rate */
	int frameLen;	  /* Length of individual input frames */
	int pvocFrameLen; /* Phase vocoder frame length: must be power of 2 and be multiples of frameLen) */
	int	pvocHop;	  /* Step size of analysis window (must be multiples of frameLen) */

	int internalBufLen;
	int maxNDelayFrames;

	int outFrameBuf_circPtr;

	/* Fixed pitch shift (semitones) 
		For operMode TIME_WARP_WITH_FIXED_PITCH_SHIFT */
	dtype fixedPitchShiftST; 

	/* Internal variables */
	int pvCtr;	/* Phase vocoder call counter */ 

	dtype expct;
	dtype osamp;
	dtype freqPerBin;

	dtype * hWin;			/* Windowing function, length = pvocFrameLen */
	dtype * xFrameW;		/* Windowed input frame, length = pvocFrameLen */	
	
	dtype * fftc;			/* Coefficients for Fourier transform */

	/* Fourier transform buffers */
	dtype * ftBuf1;			/* Length = pvocFrameLen * 2 */

public: /* DEBUG */
	dtype * ftBuf2;

private: /* DEBUG */
	dtype * lastPhase;
	dtype * lastPhase_nps;	/* No pitch shifting version of lastPhase */
	dtype * lastPhase_ntw;	/* No time warp version of lastPhase */
	dtype * sumPhase;

	bool lastPhasePrimed;

	dtype * outFrameBufPV;	/* Phase vocoder output buffer */

	dtype * X_magn;
	dtype * X_phase;
	dtype * anaMagn;
	dtype * anaFreq;
	dtype * synMagn;
	dtype * synFreq;

	/* For time-warping only */
	dtype * warpCacheMagn;		/* Cache for time warping */
	dtype * warpCachePhase;	/* Cache for time warping */

	dtype * outFrameBuf;

	/* Status variables */
	bool bPitchShift_prev;
	bool bPitchShift;

	bool bWarp_prev;
	bool bWarp;


};

/* Non-member utility functions */
const bool check_power2(const int & x);

#endif