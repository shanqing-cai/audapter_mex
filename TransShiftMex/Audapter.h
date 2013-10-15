/* 
Audapter.h

Defines the interface for the [a]-->[i] transition shifting algorithm

(c) 2007 Marc Boucek
(c) 2013 Shanqing Cai
speech communication group, RLE @MIT
*/

#pragma once

//#include "pvocWarpAtom.h"
#include <windows.h>
#include <process.h>
#include "mex.h"

typedef double dtype;

//#define N_COEFFS_SRFILT 21
//#define DOWNSAMP_FACT_DEFAULT 3
//#define maxFrameLen (768 / DOWNSAMP_FACT)
//#define MAX_FRAMELEN (1152 / DOWNSAMP_FACT_DEFAULT)
//#define MAX_BUFLEN (15 * MAX_FRAMELEN)  // max ndelay< 8
//#define MAX_NWIN   16
//#define M_1_SQRTPI 0.564189583547756
#define M_PI       3.14159265358979323846
//#define MAX_NFMTS 5
//#define MAX_NLPC	20
//#define MAX_NPOLES (MAX_NLPC / 2 + 2)
//#define MAX_NLPC_SQR (MAX_NLPC * MAX_NLPC)
//#define MAX_FMT_TRACK_JUMP MAX_NPOLES
//#define MAX_NTRACKS 5
//#define MAX_PITCHLEN 100	//100
//#define MAX_NTRANS	10		//SC How many dipthong transitions are allowed at most? In practice, only one is used.
//#define PF_NPOINTS	257
//#define PF_NBIT	8
//
//#define MAX_NTONES 64
//#define MAX_TONESEQ_REC_LEN 230400
//
//#define	MAX_ITER	10		//SC Maximum number of iteractions for the polynomial perpendicular foot finding
//
//#define NFFT 1024			//SC(2008/05/06)
//#define MAX_NFFT 4096		//SC(2012/03/05): for frequency/pitch shifting
//
//#define MAX_DATA_VEC (2 * MAX_NTRACKS + 10 + MAX_NLPC + 5)
//#define MAX_REC_SIZE  230400 // 10 sec recording at p.sr = 48kHz/4 
//#define MAX_DATA_SIZE 230400 // 10 sec data p.sr = 48kHz/4 and frameshift>=1
//
//#define MAX_PB_SIZE   230400 // 2.5 sec, playback rate = 48 kHz.
//
//#define MAX_DELAY_FRAMES  600 // The maximum number of added frame delays
//
//#define MAX_N_VOICES 4
//
//#define MAX_RMS_SLOPE_N 100

// subfunction roots specific definitions
#define SIGN(a,b) ((b)>=0.0 ? fabs((a)) : -fabs((a)))
#define IMAX(k,j) ((k)<=(j) ? (j) : (k))
#define aMat(k,j) AHess[((j)-1)*nLPC+(k)-1]
#define ISABOVE(a,b)((a)>=(b) ? true : false)
#define ZEROIFABOVE(a,b) ((a)>=(b) ? 0 : (a))

int algoCallbackFunc(char *buffer, int buffer_size, void * data); // algorithm callback function: stereo: for online runs
int algoCallbackFuncMono(char *buffer, int buffer_size, void * data); // algorithm callback function: mono: for offline simulations only
int algoCallbackFuncSineGen(char *buffer, int buffer_size, void * data); //SC algorithm sine wave generator
int algoCallbackFuncWavePB(char *buffer, int buffer_size, void * data); //SC algorithm sine wave generator

int algoCallbackFuncToneSeq(char *buffer, int buffer_size, void * data); //SC(2009/12/01) algorithm tone sequence generation

// DSPLib routines - in case we have to go back to the TI board
void	DSPF_dp_blk_move(const dtype * x, dtype * r, const int nx);
dtype	DSPF_dp_vecsum_sq(const dtype *x,int n);                           
void	DSPF_dp_biquad(dtype * x, dtype * b, dtype * a, dtype * delay, dtype * r, int nx);
void	DSPF_dp_vecmul(const dtype * x, const dtype * y, dtype * r, int n);
void	DSPF_dp_autocor(dtype * r, dtype * x, int nx, int nr);

typedef struct tag_thrWriteWavStruct {
	void *pThis;
	char *wavfn_in;
	char *wavfn_out;
} thrWriteWavStruct;

typedef struct tag_maxInterOnsetIntervalCfg {
	int n;

	int *stat0;
	dtype *maxInterval; /* Unit: s */
	int *stat1;
} maxInterOnsetIntervalCfg;

typedef struct tag_ostTab {
	int n;	// Number of segments

	int *stat0; // Initial stat number
	int *mode; // Mode number
	dtype *prm1; // First parameter
	dtype *prm2; // Second parameter
	dtype *prm3; // Third parameter

	maxInterOnsetIntervalCfg maxIOICfg;
} OST_TAB;

typedef struct tag_pipCfg { // Pitch and intensity perturbation configurtion
	int n; // Number of stats
	
	float *pitchShift;	// Unit: st
	float *intShift;	// Unit: dB

	float *fmtPertAmp;
	float *fmtPertPhi;	// Unit: phi

} PIP_CFG;

class pvocWarpAtom {
public:
	dtype tBegin;
	dtype rate1;	// Should aways be in the interval of (0, 1)
	dtype dur1;
	dtype durHold;
	dtype rate2;	// Should aways be > 1
	dtype dur2;
	int ostInitState; // The ost (online speech tracking) state number at which the clock for pvoc starts. -1: use the same clock as TransShiftMex (default)

	pvocWarpAtom(){
		tBegin = 0.25;
		dur1 = 0.25; 
		rate1 = 0.5;
		durHold = 1;
		rate2 = 2;
		/* ostInitState = -1; */
		ostInitState = 9999; /* Use a large number to effectively disable time warping by default */
			
		dur2 = (1 - rate1) / (rate2 - 1) * dur1;
	}
	
	pvocWarpAtom(dtype t_tBegin, dtype t_rate1, dtype t_dur1, dtype t_durHold, dtype t_rate2){
		tBegin = t_tBegin;
		rate1 = t_rate1;
		dur1 = t_dur1;
		durHold = t_durHold;
		rate2 = t_rate2;

		dur2 = (1 - rate1) / (rate2 - 1) * dur1;
	}

	pvocWarpAtom(int t_ostInitState, dtype t_tBegin, dtype t_rate1, dtype t_dur1, dtype t_durHold, dtype t_rate2){
		tBegin = t_tBegin;
		rate1 = t_rate1;
		dur1 = t_dur1;
		durHold = t_durHold;
		rate2 = t_rate2;

		ostInitState = t_ostInitState;

		dur2 = (1 - rate1) / (rate2 - 1) * dur1;
	}
};

// Class Audapter
// Implements the audio input - process - output queue
class Audapter {
private:
	/* Static constants */
	/* Core configuration */
	static const int nCoeffsSRFilt = 21;
	static const int downSampFact_default = 3;
	static const int maxFrameLen = 1152 / downSampFact_default;
	static const int maxBufLen = maxFrameLen * 15;
	static const int maxNWin = 16;
	static const int maxNFmts = 5;
	static const int maxNLPC = 20;
	static const int maxNPoles = maxNLPC / 2 + 2;
	static const int maxNLPC_squared = maxNLPC * maxNLPC;
	static const int maxFmtTrackJump = maxNPoles;
	static const int maxNTracks = 5;
	static const int maxPitchLen = 100;
	static const int nFFT = 1024;
	static const int max_nFFT = 4096;

	/* Data recorder */
	static const int maxDataVec = 2 * maxNTracks + 10 + maxNLPC + 5;
	static const int maxRecSize = 230400;
	static const int maxDataSize = 230400;	
	
	/* Perturbation field */
	static const int pfNPoints = 257;	/* Number of points in the perturbation field */
	static const int pfNBit = 8;

	/* Waveform player */
	static const int maxPBSize = 230400;			/* Maximum length (samples) for waveform playback */

	/* Tone sequence player */
	static const int maxNTones = 64;			/* Maximum number of tones in the sequence */
	static const int maxToneSeqRecLen = 230400;	/* Maximum length (samples) of the tone sequence waveform */

	/* Global delay */
	static const int maxDelayFrames = 600;	/* Maximum global delay (in frames) */
	static const int maxNVoices = 4;			/* Maximum number of blended voices */
	
	static const int internalBufLen = maxFrameLen * downSampFact_default * maxDelayFrames;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VARIABLES & COUNTERS  *****************************************************%%%%%%%		
	// Variables
	thrWriteWavStruct thrwws;

	// (10/19/2012): OST sentence status
	int stat;
	int lastStatEnd;
	int stretchCnt;
	int *statOnsetIndices; /* Onset index of all state numbers. By definition, statOnsets[0] = 0 */
	dtype stretchSpanAccum;

	int doEdge;
	bool   bTransReset;			// set true when resetted, set to false after first process

	bool	bFrameLenShown;

	OST_TAB ostTab;
	dtype rmsSlopeWin;
	int rmsSlopeN;

	PIP_CFG pipCfg;	

	// counters
	int frame_counter;		// frame counter : increments at every new frame from soundcard	(frame rate)						
	int data_counter;		// data counter  : increment rate is process rate ,i.e. frame rate * nWIn
	int circ_counter;		// circular data counter , loops over MAX_PITCHLEN

	/* unsigned long int frame_counter_nowarp; */
	long int frame_counter_nowarp;

	// 
	dtype time_step;       // process time unit 
	dtype ma_rms1;			// moving average rms 
	dtype ma_rms2;	
	dtype ma_rms_fb;		// moving average rms for feedback mode 4

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BUFFERS  *****************************************************%%%%%%%%%%%	

	// buffers stores samples at original sample rate
	dtype inFrameBuf[maxFrameLen * downSampFact_default];      // stores incoming  frame from soundcard @p.sr*downsamp_fact
	//dtype outFrameBuf[maxFrameLen*DOWNSAMP_FACT];     // stroes outgoing frame to soundcard @p.sr*downsamp_fact
	//SC(2012/02/28) DAF: expanded outFrameBuf for DAF
	dtype outFrameBuf[internalBufLen]; // Shifting output buffer after formant shifting, before pitch shifting
	dtype outFrameBufPS[maxNVoices][internalBufLen]; // Shifting output buffer after pitch shifting

	int outFrameBuf_circPtr;

	//SCai (2012/09/08) BlueShift
	dtype outFrameBufSum[maxFrameLen * downSampFact_default];
	dtype outFrameBufSum2[maxFrameLen * downSampFact_default];

	dtype srfilt_buf[maxFrameLen * downSampFact_default];      // multiple purpose before for down / upsampling

	// buffers stores samples at downsampled rate
	dtype filtbuf[maxFrameLen];                       // filter buffer (formant shift)
	dtype oBuf[maxBufLen];							// Buffer stores original input samples
	dtype pBuf[maxBufLen];							//SC-Mod(2008/05/15) Buffer stores preemphasized input samples
	dtype inBuf[maxFrameLen];							// Buffer stores input samples after downsampling
	dtype outBuf[maxFrameLen];							// Buffer stores processed output
	dtype fakeBuf[maxFrameLen];							// trash buffer to flush filter states
	dtype zeros[maxFrameLen];						    // utility buffer for filtering


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VARIOUS EXTRACTED DATA  *****************************************************%%%%%%%%%%%	


	// other buffers : process rate = downsampled rate * nWin
	dtype lpcAi[maxNLPC+1];						// lpc coeeficients	
	dtype hwin[maxBufLen];						// Hanning window
	dtype hwin2[maxBufLen];						// Hanning window for frequency/pitch shifting
	dtype realRoots[maxNLPC];						// real part of roots
	dtype imagRoots[maxNLPC];						// imag part of roots
	dtype Acompanion[maxNLPC_squared];				// companion matrix for (eigenvalue ---> )roots calculation
	dtype amps[maxNPoles];						// radius of poles
	dtype orgPhis[maxNPoles];						// angle of original poles
	dtype bw[maxNPoles];							// bandwith of poles

	dtype newPhis[maxNPoles];						// new (shifted) pole angles 
	dtype gtot[maxNWin];							// gain factor (for gain adaption)

	dtype weiMatPhi[maxNTracks][maxPitchLen];	// weithed matrix of past pole angles (for moving average formant smoothing)
	dtype weiMatBw[maxNTracks][maxPitchLen];		// weithed matrix of past bandwith (for moving average formant smoothing)
	dtype weiVec[maxPitchLen];					// weigthing vector of past weigths (default weigth : short time rms ... )
	dtype sumWeiPhi[maxNTracks];					// sum of past weights (rms) *  angles (over avglen )
	dtype sumWeiBw[maxNTracks];                   // sum of past weights (rms) *  angles (over avglen )
	dtype sumWei;									// sum of weigths (over pitchlen)

	dtype wmaPhis[maxNTracks];					// weighted moving average pole angles
	dtype wmaR[maxNTracks];						// weighted moving average pole radius
	dtype deltaPhi[maxNFmts];						// derivative of angle
	dtype fmts[maxNTracks];						// estimated formant tracks
	dtype dFmts[maxNFmts];						// formant derivatives
	dtype sFmts[maxNFmts];						// shifted formants (=0 when no shift)

	// getDFmt
	dtype deltaFmt[2];								// formant derivatives
	dtype deltaMaFmt[2];							// moving average deltafmt
	dtype lastFmt[2];								// last formant estimate

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FILTER COEEFS AND DELAYS  *****************************************************	

	// preemphasis 
	dtype a_preemp[2];								// denominator coefficients 
	dtype b_preemp[2];								// numerator coefficients
	dtype preemp_delay[1];							// filter delay

	// deemphasis
	dtype a_deemp[2];								// denominator coefficients 
	dtype b_deemp[2];								// numerator coefficients
	dtype deemp_delay[1];								// filter delay


	// first filter (f1 shift)
	dtype a_filt1[2];								// denominator coefficients 
	dtype b_filt1[3];								// numerator coefficients
	dtype filt_delay1[2];							// filter delays


	// second filter (f2 shift)
	dtype a_filt2[2];								// denominator coefficients 
	dtype b_filt2[3];								// numerator coefficients
	dtype filt_delay2[2];							// filter delays

	// sample rate conversion filter 
	dtype srfilt_a[nCoeffsSRFilt];				// denominator coefficients
	dtype srfilt_b[nCoeffsSRFilt];				// numerator coefficients
	dtype srfilt_delay_up[nCoeffsSRFilt-1];		// filter delays for upsampling
	dtype srfilt_delay_down[nCoeffsSRFilt-1];    // filter delays for downsampling

	//SC(2008/05/07)
	dtype ftBuf1[nFFT*2];
	dtype ftBuf2[nFFT*2];
	dtype ftBuf1ps[max_nFFT * 2]; // For frequency/pitch shifting
	dtype ftBuf2ps[2][max_nFFT * 2]; // For frequency/pitch shifting: (Normal / nps)	

	//SC(2012/03/05): Pitch Shifting
	dtype lastPhase[2][max_nFFT]; // (Normal | nps)
	dtype lastPhase_ntw[max_nFFT]; // ntw: no time warping	
	dtype sumPhase[2][maxNVoices][max_nFFT]; // (Normal | nps)	
	dtype outputAccum[max_nFFT];

	dtype phase0, phase1;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DATA   RECORDING  *****************************************************%%%%%%%%%%%	


	//recorded data is stored here
	dtype signal_recorder[2][maxRecSize];				// downsampled input signal and downsampled output signal
	dtype data_recorder[maxDataVec][maxDataSize];		// stores other data
	dtype transdata[2];								// only transition data

	dtype a_rms_o[maxDataSize];
	dtype a_rms_o_slp[maxDataSize];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DATA PLAYBACK     **************************************************************
	//SC-Mod(2007/12/28)
	dtype data_pb[maxPBSize];

	//SC-Mod(2008/04/07) Vowel related
	dtype f1m;     // F1 (mel)
	dtype f2m;     // F2 (mel)
	dtype f2mp;	// F2 (mel) of the foot of perpendicular to the fit polynomial.

	//SC-Mod(2008/05/15) FFT related
	dtype fftc[nFFT*2];
	dtype fftc_ps0[max_nFFT * 2];
	dtype fftc_ps[max_nFFT * 2];	//SC(2012/03/05) For frequency/pitch shifting

	//SC(2009/02/02)
	int bLastFrameAboveRMS;

	//SC(2012/02/29) For data writing:
	int dataFileCnt;

	//SC(2012/03/13) For PVOC time warping
	dtype pvocWarpCache[internalBufLen / 64][1024 * 2];

	dtype rms_ratio;

	//SCai(2012/10/19) PIP intensity shift
	dtype intShiftRatio;

	//SC(2009/12/01) Tone sequence generator
	dtype tsgToneOnsets[maxNTones];		// sec

	//SC(2013/07/19) Feedback mode 4 (modulated noise)
	dtype rmsFF_fb_now;
	int fb4_status;			
	// 0 - Before voicing onset; 1 - During onset ramping (Kernel size ramping up); 
	// 2 - During utterance; 3 - Tentative voice ending; 4 - Voice ended; 5 - Kernel size ramping down.

	int fb4_counter;

	dtype amp_ratio;
	dtype amp_ratio_prev;

	//pvocWarpAtom warpCfg;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS  *****************************************************%%%%%%%%%%%

	// Externally adjustable parameters
	struct{
		// framing and processing
		int	   frameLen;				// length of one frame ( framelen = nWin * frameshift)
		int	   nWin;					// number of processes per frame (of frameLen samples)
		int	   frameShift;				// number of samples shift between two processes ( = size of processed samples in 1 process)	
		int	   bufLen;					// main buffer length : buflen stores (2*nDelay -1)*frameLen samples
		int    nDelay;					// number of delayed framas (of size framelen) before an incoming frame is sent back
		// thus the overall process latency (without souncard) is :Tproc=nDelay*frameLen/sr
		int	   anaLen;					// size of lpc analysis (symmetric around window to be processed)
		int	   pvocFrameLen;			// size of analysis window for pitch shifting (aim for power of 2) 
		int	   pvocHop;					// step size of analysis window (must be multiples of frameLen)

		int    avgLen;				// length of smoothing ( should be approx one pitch period, 
		// can be greater /shorter if you want more / lesss moothing)
		// avgLen = 1 ---> no smoothing ( i.e. smoothing over one value)

		int		downFact;				// Down-sampling factor (Default: downSampFact_default)
		int    sr;						// internal samplerate

		int    nLPC;					// lpc order ... number of lpc coeefs= nLPC +1
		int    nFmts;					// originally the number of formants you want shift ( hardseted = 2)

		dtype dRMSThresh;				// RMS threshhold for voicing detection
		dtype dRMSRatioThresh;			// preemp / original RMS ratio threshhold, for fricative detection 
		dtype rmsFF;                   // rms fogetting factor for long time rms 
		dtype rmsFF_fb[4];             // rms fogetting factor for fb mode 4 {ff_onset, ff_final, trans_time (s), ending time (s)} 

		dtype fb3Gain;					// gain (scaling factor) of the blended noise under fb mode 3

		dtype fb4GainDB;				// gain (in dB) of the speech-modulated noise under fb mode 4
		dtype fb4Gain;					// gain (scaling factor) of the speech-modulated noise under fb mode 4

		dtype dPreemp;					// preemphasis factor
		dtype dScale;					// scaling the output (when upsampling) (does not affect internal signal

		// for transition detection
		dtype dFmtsFF;					// formant forgeeting factor for s
		dtype maxDelta;				// maximal allowed formant derivative 		
		dtype fmtDetectStart[2];		// formant frequencies at start of transition (goal region),i.e. [a]		
		int    minDetected;				// min transition detection in row before starting transition		


		// for formant tracking algorithm
		dtype trackFF;
		int	   nTracks;					// number of tracked formants 
		int    nCands;					// number of possible formant candiates  ( > ntracks     but  < p.nLPC/2!!!! (choose carefully : not idiot proofed!)

		dtype trackIntroTime;          // intro tracking time 

		//dtype aFact1;				// 
		dtype bFact;
		//dtype gFact1;
		dtype aFact;
		//dtype bFact2;
		dtype gFact;

		dtype fn1;		// Neutral F1: used in formant tracking
		dtype fn2;		// Neutral F2: used in formant tracking			

		// for shifting
		dtype dDev;					// deviation vector (in angle euclidien distance) @ max of hanning window

		// booleans						
		int    bRecord;					// record signal
		int    bTrack;					// use formant tracking algorithm
		int    bShift;					// do shifting
		int    bGainAdapt;				// use gain adaption
		int    bDetect;					// detect transition	
		int    bRelative;				// shift relative to actual formant point, (otherwise absolute coordinate)			
		int    bWeight;					// do weighted moving average formant smoothing (over pitchlen) , otherwise not weigthed (= simple moving average)				
		int	   bCepsLift;				//SC-Mod(2008/05/15) Whether the cepstral lifting is done before the autocorrelation

		int	   bRatioShift;				//SC(2009/01/20). 
		int	   bMelShift;				//SC(2009/01/20). 

		//SC(2012/03/05) Frequency/pitch shifting
		int		bPitchShift;
		dtype	pitchShiftRatio[maxNVoices];

		
		//SC-Mod(2008/01/05). Arrays: for intensity correction during formant shifting (mainly for the use of parallel shifts)		
		
		dtype wgFreq;										//SC Wave generator frequency (Hz)
		dtype wgAmp;										//SC Wave generator amplitude (digitized peak amplitude)
		dtype wgTime;										//SC Wave generator current phase (rad)


		//SC(2008/04/03). Perturbation-related variables		
		int transCounter;
		dtype F2Min;
		dtype F2Max;
		//dtype triF2MinMin;
		//dtype triF2MaxMax;
		dtype F1Min;
		dtype F1Max;
		dtype LBk;
		dtype LBb;
		dtype pertF2[pfNPoints];		
		dtype pertPhi[pfNPoints];
		dtype pertAmp[pfNPoints];
		dtype minVowelLen;
		
		bool transDone;

		//SC(2008/06/20) Feedback type: fb=0: silent; fb=1: voice only (pert. or unpert.); fb=2: noise only; fb=3: voice+noise;
		int fb;		
		
		//SC(2008/06/22)
		int bRepData;		// Replicate data or replay, in mode 1. 0: no; 1: yes.

		//SC(2008/05/15) Cepstral lifting related
		int	cepsWinWidth;

		//SC(2008/06/22)
		dtype trialLen;
		dtype rampLen;

		//SC(2009/02/06)
		int bRMSClip;
		dtype rmsClipThresh;

		//SC(2012/02/28)
		int delayFrames[maxNVoices];	// For DAF: the number of frames delayed. 0 correponds to no _added_ delay (~11 ms).
		dtype gain[maxNVoices];
		int mute[maxNVoices];

		int bDownSampFilt;

		//SC(2012/09/08) BlueShift
		int nFB;	// Number of feedback voices

		//SC(2012/09/24) For saving to wav files

		//SC(2009/12/01) Parameters related to the tone sequence generator
		int tsgNTones;
		dtype tsgToneDur[maxNTones];	// sec
		dtype tsgToneFreq[maxNTones];	// Hz
		dtype tsgToneAmp[maxNTones];	// Peak value
		dtype tsgToneRamp[maxNTones];	// sec
		dtype tsgInt[maxNTones];		// sec The interval between the onsets of tones (NOT between the offset of a tone and the onset of the next one)

		/*SC(2013/04/07) Options to bypass the formant tracker (useful for situations in which lower latency under pitch shifting or time warping is required */
		int bBypassFmt;

		/* SC (2013-08-06) stereoMode */
		int stereoMode;		/* 0 - left only; 1 - left-right identical; 2 - left audio + right simulate TTL */

		int bPvocAmpNorm;	/* Pitch vocoder amplitude normalization */		
		int pvocAmpNormTrans; /* Length of the amplitude normalization factor transitional period */
	} p;


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  *****************************************************%%%%%%%%%%%

	void     getAi(dtype* xx, dtype* aa, const int size, const int nlpc);
	dtype  getGain(dtype * r, dtype * ophi,dtype * sphi, int nfmts);
	int     hqr_roots (	dtype *c, 	dtype *wr, dtype *wi,	dtype *Acompanion, const int nLPC);
	void    getRPhiBw (dtype *wr,  dtype *wi, dtype *radius,  dtype *phi ,dtype *bandwith);
	void trackPhi(dtype *r_ptr,dtype *phi_ptr,dtype time);
	void myFilt (dtype *xin_ptr, dtype* xout_ptr,dtype *oldPhi_ptr,dtype *newPhi_ptr,dtype *r_ptr,const int size);
	dtype calcRMS1(const dtype *xin_ptr, int size);	
	dtype calcRMS2(const dtype *xin_ptr, int size);
	dtype calcRMS_fb(const dtype *xin_ptr, int size, bool above_rms);
	int getWma(dtype *phi_ptr, dtype *bw_ptr , dtype * wmaPhi_ptr, dtype * wmaR_ptr);
	void    levinson(dtype *R, dtype* aa, int size);
	int gainAdapt(dtype *buffer,dtype *gtot_ptr,int framelen, int frameshift);
	int gainPerturb(dtype *buffer,dtype *gtot_ptr,int framelen, int frameshift);

	int sign(dtype x);
	bool  detectTrans(dtype *fmt_ptr, dtype *dFmt_ptr,int datcnt, dtype time);
	int getDFmt(dtype *fmt_ptr,dtype *dFmt_ptr, dtype time);	
	void upSampSig (dtype *b, dtype *a, dtype *x, dtype *buffer, dtype *r,dtype  *d,const int nr, const int n_coeffs, const int upfact,const dtype scalefact);
	void downSampSig(dtype *b, dtype *a, dtype *x, dtype *buffer, dtype *r,dtype  *d,const int nr, const int n_coeffs, const int downfact);
	void downSampSig_noFilt(dtype *b, dtype *a, dtype *x, dtype *buffer, dtype *r,dtype  *d,const int nr, const int n_coeffs, const int downfact);
	void iir_filt (dtype *b, dtype *a,  dtype *x, dtype *r,dtype  *d,const int nr, const int n_coeffs,  dtype g);
			
	dtype Audapter::hz2mel(dtype hz);
	dtype Audapter::mel2hz(dtype hz);
	dtype Audapter::locateF2(dtype f2);

	void	DSPF_dp_cfftr2(int n, dtype * x, dtype * w, int n_min);
	void	DSPF_dp_icfftr2(int n, double * x, double * w, int n_min);
	void	gen_w_r2(double* w, int n);
	void	bit_rev(double* x, int n);

	void	smbFft(dtype *fftBuffer, double fftFrame_Size, int sign);
	void	Audapter::calcRMSSlope();
	void	Audapter::osTrack();


public:
	// Constructor initializes all variables
	Audapter();

	// Destructor 
	~Audapter();

	// The Reset function reinitializes all internal buffers
	void reset();

	// The main function that controls the processing
	int handleBuffer(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size, bool bSingleOutputBuffer);
	int handleBufferSineGen(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size);	//SC Since wave generator
	int handleBufferWavePB(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size);	//SC Wave playback

	int handleBufferToneSeq(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size);	//SC(2009/12/01) Tone sequence generator

	// Allows external entities to set / get various parameters of the process
	int setparams(void * name, void * value, int nPars);
	int getparams(void * name);
	const dtype* getsignal(int & size);
	const dtype* getdata(int & size,int & vecsize);	
	const dtype* getOutFrameBufPS();

	int pbCounter;	//SC The integer counter used in wave playback 

	pvocWarpAtom *warpCfg;

	char deviceName[256];

	char wavFileBase[256];

	//void writeSignalsToWavFile(char *wavfn_input, char *wavfn_output);
	void writeSignalsToWavFile();

	static unsigned __stdcall Audapter::thrStatEntPnt(void *pThrStruct) { // For threading: writing to wav files
		thrWriteWavStruct *thrStruct = (thrWriteWavStruct *)pThrStruct;
		Audapter *p_this = (Audapter *)(thrStruct->pThis);
			
		//mexPrintf("Thread executing!\n");
		//p_this->writeSignalsToWavFile(thrStruct->wavfn_in, thrStruct->wavfn_out);
		p_this->writeSignalsToWavFile();
		fflush(stdout);

		return 0;
	}
	
	// OST related (10/18/2012)
	char ostfn[512];
	char pipcfgfn[512];

	bool duringTimeWarp, duringTimeWarp_prev;
	bool duringPitchShift, duringPitchShift_prev;

	void Audapter::readOSTTab(int bVerbose);
	void Audapter::readPIPCfg(int bVerbose);

	dtype tsg_wf[maxToneSeqRecLen];
	int tsgRecCounter;

	/* Inline constant getters */
	const int getMaxRecSize() const { return maxRecSize; };
	const int getMaxDataSize() const { return maxDataSize; };
	const int getMaxFrameLen() const { return maxFrameLen; };
	const int getMaxDelayFrames() const {return maxDelayFrames; };
};

void init_ostTab(OST_TAB *ostTab);
void init_pipCfg(PIP_CFG *pipCfg);