/* 
Audaper.cpp

Speech auditory feedback manipulation program capable of the following types 
of perturbations:
	1) Formants frequencies (F1 and F2)
	2) Pitch
	3) Fine-scale timing (time warping)
	4) Global time delay (DAF)
	5) Local intensity
	6) Global intensity

Heuristics-based online utterance status tracking (OST) is incorporated

Requires ASIO compatible sound cards

Also incorporated 
	1) Cepstral liftering (optional)
	2) Dynamic programming style formant tracking 
		(Xia and Espy-Wilson, 2000, ICSLP)
	3) Gain adaptation after shifting (optional)
	4) RMS-based formant smoothing (to reduce the influence of subglottal 
		coupling on formant estimates) (optional)

Authors:
2007 Marc Boucek, Satrajit Ghosh
2008-2013 Shanqing Cai (shanqing.cai@gmail.com)

Developed at:
Speech Communication Group, RLE, MIT
Speech Laboratory, Boston University
*/

#pragma once

#include <windows.h>
#include <process.h>
#include <math.h>
#include <memory>
#include <vector>

#include "mex.h"

#include "DSPF.h"
#include "filter.h"
#include "lpc_formant.h"
#include "ost.h"
#include "pcf.h"
#include "phase_vocoder.h"
#include "time_domain_shifter.h"

typedef double dtype;

#define M_PI       3.14159265358979323846

/* Utility inline functions */
inline dtype mul_sign(const dtype &a, const dtype &b) {
	return (b >= 0.0 ? fabs(a) : -fabs(a));
}

inline int sign(const dtype &x) {
	if (x>0)
		return 1;
	else if (x<0)
		return -1;
	else
		return 0;
}

inline int imax(const int &k, const int &j) {
	return (k <= j ? j : k);
}

inline bool isabove(const dtype &a, const dtype &b) {
	return (a >= b);
}

/* Callback function declarations */
int audapterCallback(char *buffer, int buffer_size, void * data); // algorithm callback function: stereo: for online runs


typedef struct tag_thrWriteWavStruct {
	void *pThis;
	char *wavfn_in;
	char *wavfn_out;
} thrWriteWavStruct;


class Parameter {
public:
	typedef enum {
		TYPE_NULL = 0,
		TYPE_BOOL, 
		TYPE_INT, 
		TYPE_DOUBLE, 
		TYPE_BOOL_ARRAY,
		TYPE_INT_ARRAY, 
		TYPE_DOUBLE_ARRAY, 
		TYPE_DOUBLE_2DARRAY,
		TYPE_PVOC_WARP, 
		TYPE_SMN_RMS_FF,
        TYPE_TIME_DOMAIN_PITCH_SHIFT_SCHEDULE,
        TYPE_TIME_DOMAIN_PITCH_SHIFT_ALGORITHM,
	} paramType;

private:
	std::vector<const char *> names;
	std::vector<const char *> helpMsgs;	
	std::vector<paramType> types;
	int nParams;

public:
	Parameter() { nParams = 0; };

    void addParam(const char *name, const char * helpMsg, const paramType type) {
        names.push_back(name);
        helpMsgs.push_back(helpMsg);
        types.push_back(type);

        nParams = names.size();
    };
    void addBoolParam(const char* name, const char* helpMsg) {
        addParam(name, helpMsg, TYPE_BOOL);
    };
    void addIntParam(const char* name, const char* helpMsg) {
        addParam(name, helpMsg, TYPE_INT);
    }
    void addIntArrayParam(const char* name, const char* helpMsg) {
        addParam(name, helpMsg, TYPE_INT_ARRAY);
    }
    void addDoubleParam(const char* name, const char* helpMsg) {
        addParam(name, helpMsg, TYPE_DOUBLE);
    }
    void addDoubleArrayParam(const char* name, const char* helpMsg) {
        addParam(name, helpMsg, TYPE_DOUBLE_ARRAY);
    }
	void addDouble2DArrayParam(const char* name, const char* helpMsg) {
		addParam(name, helpMsg, TYPE_DOUBLE_2DARRAY);
	}

	paramType checkParam(const char *name);
};

// Class Audapter
// Implements the audio input - process - output queue
class Audapter {
private:
	/* Static constants */
	/* Core configuration */
	static const int nCoeffsSRFilt = 21;
	static const int downSampFact_default = 3;
	static const int maxFrameLen = 2880 / downSampFact_default;
	static const int maxBufLen = maxFrameLen * 15;
	static const int maxNWin = 16;
	static const int maxNFmts = 5;
	static const int maxNLPC = 20;
	static const int maxNPoles = maxNLPC / 2 + 2;
	static const int maxNTracks = 5;
	static const int maxPitchLen = 100;
	static const int nFFT = 1024;  // TODO(cais): Make into adjustable parameter. Original: 1024.
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

	/* Parameter names and help info */
	Parameter params;

	/* Display related */
	static const int maxDisplayElem = 6;

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VARIABLES & COUNTERS  *****************************************************%%%%%%%		
	// Variables
	thrWriteWavStruct thrwws;

	// (10/19/2012): OST sentence status
	int stat;

	int doEdge;
	bool   bTransReset;			// set true when resetted, set to false after first process

	bool	bFrameLenShown;

	OST_TAB ostTab;
	dtype rmsSlopeWin;
	int rmsSlopeN;

	PERT_CFG pertCfg;

	// counters
	int frame_counter;		// frame counter : increments at every new frame from soundcard	(frame rate)						
	int data_counter;		// data counter  : increment rate is process rate ,i.e. frame rate * nWIn
	int circ_counter;		// circular data counter , loops over MAX_PITCHLEN

	/* unsigned long int frame_counter_nowarp; */
	long int frame_counter_nowarp;

	/* LP formant tracker object */
	std::unique_ptr<LPFormantTracker> fmtTracker;

	/* Phase vocoder object */
    std::vector<std::unique_ptr<PhaseVocoder>> pVocs;

    /* Time-domain pitch shifter */
    std::unique_ptr<audapter::TimeDomainShifter> timeDomainShifter;

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
    dtype f0Buf[maxFrameLen];  // Buffer for F0-related operations.
	dtype zeros[maxFrameLen];						    // utility buffer for filtering

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VARIOUS EXTRACTED DATA  *****************************************************%%%%%%%%%%%	
	// other buffers : process rate = downsampled rate * nWin

	dtype amps[maxNPoles];						// radius of poles
	dtype orgPhis[maxNPoles];						// angle of original poles
	dtype bw[maxNPoles];							// bandwith of poles

	dtype newPhis[maxNPoles];						// new (shifted) pole angles 
	dtype gtot[maxNWin];							// gain factor (for gain adaption)

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

	/* Filters */
	IIR_Filter<dtype> preEmpFilter; /* Pre-emphasis filter */
	IIR_Filter<dtype> deEmpFilter; /* De-emphasis filter */
	IIR_Filter<dtype> downSampFilter;
	IIR_Filter<dtype> upSampFilter;
	dtype downSampBuffer[maxFrameLen * downSampFact_default];
	dtype upSampBuffer[maxFrameLen * downSampFact_default];

	IIR_Filter<dtype> shiftF1Filter;	/* Biquad filter for shifting F1 */
	IIR_Filter<dtype> shiftF2Filter;	/* Biquad filter for shifting F2 */
	IIR_Filter<dtype> f0Filter;  // Filter for pitch (e.g., bandpass filtering around pitch).

	// first filter (f1 shift)
	dtype a_filt1[3];								// denominator coefficients 
	dtype b_filt1[3];								// numerator coefficients

	// second filter (f2 shift)
	dtype a_filt2[3];								// denominator coefficients 
	dtype b_filt2[3];								// numerator coefficients

	// sample rate conversion filter 

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

	//SC(2009/02/02)
	int bLastFrameAboveRMS;

	//SC(2012/02/29) For data writing:
	int dataFileCnt;

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

		dtype fb2Gain;					// gain (scaling factor) of the noise under fb mode 2

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

		// Parameters related to the pitch tracker.
		int    bTimeDomainShift;
		dtype  pitchLowerBoundHz;
		dtype  pitchUpperBoundHz;

		int	   bRatioShift;				//SC(2009/01/20). 
		int	   bMelShift;				//SC(2009/01/20). 

		//SC(2012/03/05) Frequency/pitch shifting
		int		bPitchShift;
		dtype	pitchShiftRatio[maxNVoices];
        audapter::TimeDomainShifter::PitchShiftSchedule timeDomainPitchShiftSchedule;
        audapter::TimeDomainShifterAlgorithm timeDomainPitchShiftAlgorithm;

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
		dtype pertF1[pfNPoints];
		dtype pertF2[pfNPoints];		
		dtype pertPhi[pfNPoints];
		dtype pertAmp[pfNPoints];
		dtype pertPhi2D[pfNPoints][pfNPoints];
		dtype pertAmp2D[pfNPoints][pfNPoints];
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

		// Switch for using F1 and F2 formant perturbation, instead of just F2.
		int bShift2D;

		/* SC (2013-08-06) stereoMode */
		int stereoMode;		/* 0 - left only; 1 - left-right identical; 2 - left audio + right simulate TTL */

		int bPvocAmpNorm;	/* Pitch vocoder amplitude normalization */		
		int pvocAmpNormTrans; /* Length of the amplitude normalization factor transitional period */
	} p;


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  *****************************************************%%%%%%%%%%%	
	dtype  getGain(dtype * r, dtype * ophi,dtype * sphi, int nfmts);
		
	void formantShiftFilter(dtype* xIn, dtype* xOut,
		                    dtype oldPhis[2], dtype newPhis[2], dtype mags[2], const int size);
	void f0BandpassFilter(dtype* xIn, dtype* xOut, dtype f0, dtype sr, const int size);
	dtype calcRMS1(const dtype *xin_ptr, int size);	
	dtype calcRMS2(const dtype *xin_ptr, int size);
	dtype calcRMS_fb(const dtype *xin_ptr, int size, bool above_rms);	
	
	int gainAdapt(dtype *buffer,dtype *gtot_ptr,int framelen, int frameshift);
	int gainPerturb(dtype *buffer,dtype *gtot_ptr,int framelen, int frameshift);

	bool  detectTrans(dtype *fmt_ptr, dtype *dFmt_ptr,int datcnt, dtype time);
	int getDFmt(dtype *fmt_ptr,dtype *dFmt_ptr, dtype time);	
	/*void upSampSig (dtype *b, dtype *a, dtype *x, dtype *buffer, dtype *r,dtype  *d,const int nr, const int n_coeffs, const int upfact,const dtype scalefact);*/
	void upSampSig(IIR_Filter<dtype> & usFilt, dtype *x, dtype *r, const int nr, const int upfact, const dtype scalefact);
	/*void downSampSig(dtype *b, dtype *a, dtype *x, dtype *buffer, dtype *r,dtype  *d,const int nr, const int n_coeffs, const int downfact);*/
	void downSampSig(dtype *x, dtype *r, const int nr, const int downfact, const bool bFilt);
	/*void downSampSig_noFilt(dtype *b, dtype *a, dtype *x, dtype *buffer, dtype *r,dtype  *d,const int nr, const int n_coeffs, const int downfact);*/
	
	dtype Audapter::hz2mel(dtype hz);
	dtype Audapter::mel2hz(dtype hz);
	dtype Audapter::locateF1(dtype f1);
	dtype Audapter::locateF2(dtype f2);

	void	DSPF_dp_cfftr2(int n, dtype * x, dtype * w, int n_min);
	void	DSPF_dp_icfftr2(int n, double * x, double * w, int n_min);
	void	bit_rev(double* x, int n);

	void	smbFft(dtype *fftBuffer, double fftFrame_Size, int sign);
	void	Audapter::calcRMSSlope();
	void	Audapter::osTrack();

	void *Audapter::setGetParam(bool bSet, const char *name, void * value, int nPars, bool bVerbose, int *length);

	void initializePreEmpFilter();

    void checkParameters() const;

public:
	/* Action modes */
	enum {
		PROC_AUDIO_INPUT_OFFLINE, 
		PROC_AUDIO_INPUT_ONLINE, 
		GEN_SINE_WAVE, 
		WAV_PLAYBACK, 
		GEN_TONE_SEQ
	} actionMode;

	// Constructor initializes all variables
	Audapter();

	// Destructor 
	~Audapter();

	/* Functions for parameter query */
	void queryParam(char const *name, mxArray **output);

	// The Reset function reinitializes all internal buffers
	void reset();

	// The main function that controls the processing
	int handleBuffer(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size, bool bSingleOutputBuffer);
	int handleBufferSineGen(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size);	//SC Since wave generator
	int handleBufferWavePB(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size);	//SC Wave playback

	int handleBufferToneSeq(dtype *inFrame_ptr, dtype *outFrame_ptr, int frame_size);	//SC(2009/12/01) Tone sequence generator

	// Allows external entities to set / get various parameters of the process
	void setParam(const char *name, void * value, int nPars, bool bVerbose=false);
	void *getParam(const char *name);

	const dtype* getSignal(int & size) const;
	const dtype* getData(int & size,int & vecsize) const;
	const dtype* getOutFrameBufPS() const;

	int pbCounter;	//SC The integer counter used in wave playback 

	char deviceName[256];

	char wavFileBase[256];

	//void writeSignalsToWavFile(char *wavfn_input, char *wavfn_output);
	void writeSignalsToWavFile();

	/* For threading: writing to wav files */
	static unsigned __stdcall Audapter::thrStatEntPnt(void *pThrStruct) { 
		thrWriteWavStruct *thrStruct = (thrWriteWavStruct *)pThrStruct;
		Audapter *p_this = (Audapter *)(thrStruct->pThis);
			
		//mexPrintf("Thread executing!\n");
		//p_this->writeSignalsToWavFile(thrStruct->wavfn_in, thrStruct->wavfn_out);
		p_this->writeSignalsToWavFile();
		fflush(stdout);

		return 0;
	}
	
	// OST related (10/18/2012)
	char ostFN[512];
	char pertCfgFN[512];

	bool duringTimeWarp, duringTimeWarp_prev;
	bool duringPitchShift, duringPitchShift_prev;

	void Audapter::readOSTTab(int bVerbose);
	void Audapter::readPertCfg(int bVerbose);

	dtype tsg_wf[maxToneSeqRecLen];
	int tsgRecCounter;

	/* Inline constant getters */
	const int getMaxRecSize() const { return maxRecSize; };
	const int getMaxDataSize() const { return maxDataSize; };
	const int getMaxFrameLen() const { return maxFrameLen; };
	const int getMaxDelayFrames() const {return maxDelayFrames; };

	const int getMaxPBSize() const { return maxPBSize; };
};
