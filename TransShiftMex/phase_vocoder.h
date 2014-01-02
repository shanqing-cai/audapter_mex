#ifndef PHASE_VOCODER_H
#define PHASE_VOCODER_H

#define M_PI       3.14159265358979323846

typedef double dtype;

class PhaseVocoder {
	/* TODOs: */

public:
	/* Error classes */
	class initializationError {};

	/* Constructor */
	PhaseVocoder(const dtype t_sr, 
				 const int t_frameLen, 
				 const int t_pvocFrameLen, 
				 const int t_pvocHop);

	/* Destructor */
	~PhaseVocoder();

	/* Main interface function: process a frame */
	void procFrame(const dtype * inBuf, const dtype shift);

	/* Status reset */
	void reset();

private:
	/* Parameters */
	dtype sr;			  /* Sampling rate */
	int frameLen;	  /* Length of individual input frames */
	int pvocFrameLen; /* Phase vocoder frame length: must be power of 2 and be multiples of frameLen) */
	int	pvocHop;	  /* Step size of analysis window (must be multiples of frameLen) */

	int internalBufLen;
	int outFrameBuf_circPtr;

	/* Internal variables */
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
	dtype * pvocWarpCache;	/* Cache for time warping */

	dtype * lastPhase;
	dtype * lastPhase_ntw;	/* No time warp version of lastPhase */
	dtype * sumPhase;

	dtype * outFrameBufPV;	/* Phase vocoder output buffer */

	dtype * X_magn;
	dtype * X_phase;
	dtype * anaMagn;
	dtype * anaFreq;
	dtype * synMagn;
	dtype * synFreq;

	dtype * outFrameBuf;

};

/* Non-member utility functions */
const bool check_power2(const int & x);

#endif