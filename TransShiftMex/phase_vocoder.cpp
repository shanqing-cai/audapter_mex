#include <cmath>

#include "phase_vocoder.h"
#include "DSPF.h"

/* Default constructor */
PhaseVocoder :: PhaseVocoder() :
	mode(operMode::PITCH_SHIFT_ONLY), nDelay(0), 
	sr(16000), frameLen(32), 
	pvocFrameLen(1024), pvocHop(256)
{
	internalBufLen = pvocFrameLen * 2;

	/* Windowing function */
	hWin = new dtype[pvocFrameLen];
	for(int i = 0; i < pvocFrameLen; ++i){
		hWin[i] = 0.5 * cos((2 * M_PI * static_cast<dtype>(i)) / static_cast<dtype>(pvocFrameLen)); 
		hWin[i] = 0.5 - hWin[i];
	}

	xFrameW = new dtype[pvocFrameLen];

	/* FFT coefficients */
	fftc = new dtype[pvocFrameLen * 2];
	gen_w_r2(fftc, pvocFrameLen);

	/* FFT intermediate data fields */
	ftBuf1 = new dtype[pvocFrameLen * 2];
	ftBuf2 = new dtype[pvocFrameLen * 2];

	lastPhase = new dtype[pvocFrameLen]; /* TODO: Check size bound is tight */
	lastPhase_nps = new dtype[pvocFrameLen];
	lastPhase_ntw = new dtype[pvocFrameLen];
	sumPhase = new dtype[pvocFrameLen];

	outFrameBufPV = new dtype[internalBufLen];

	X_magn = new dtype[pvocFrameLen];
	X_phase = new dtype[pvocFrameLen];
	anaMagn = new dtype[pvocFrameLen];
	anaFreq = new dtype[pvocFrameLen];
	synMagn = new dtype[pvocFrameLen];
	synFreq = new dtype[pvocFrameLen];

	outFrameBuf = new dtype[internalBufLen];

	/* Internal variables */
	expct = 2.* M_PI * pvocHop / pvocFrameLen; /* Expected phase change */
	osamp = pvocFrameLen / pvocHop; /* Oversampling factor */
	freqPerBin = sr / pvocFrameLen; /* Frequency per bin */

	if (mode == TIME_WARP_ONLY || 
		mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		maxNDelayFrames = pvocFrameLen * 8;	/* TODO: Check that the delay doesn't exceed this limit */

		warpCacheMagn = new dtype[maxNDelayFrames * pvocFrameLen];
		warpCachePhase = new dtype[maxNDelayFrames * pvocFrameLen];
	}
	else {
		maxNDelayFrames = 0;

		warpCacheMagn = 0;
		warpCachePhase = 0;
	}

	if (mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		fixedPitchShiftST = static_cast<dtype>(0xffffffff); /* NAN */
	}

	reset();
}

/* Constructor 
	Input argument:
		t_operMode:		operation mode of the phase vocoder
		nDelay:			amount of global delay (# of frames)
		t_sr:			sampling rate (Hz) 
		t_frameLen:		frame length (# of samples) 
		t_pvocFrameLen: phase vocoder frame length (# of samples) 
		t_pvocHop:		phase vocoder frame hop (# of samples, < t_pvocFrameLen 
							and be one of its divisors) 
*/
PhaseVocoder :: PhaseVocoder(const operMode t_operMode, 
							 const int t_nDelay, 
							 const dtype t_sr, 
							 const int t_frameLen, 
							 const int t_pvocFrameLen, 
							 const int t_pvocHop) :
	mode(t_operMode), nDelay(t_nDelay), 
	sr(t_sr), frameLen(t_frameLen), 
	pvocFrameLen(t_pvocFrameLen), pvocHop(t_pvocHop)
{
	/* Input sanity checks */
	if ( sr <= 0 )
		throw initializationError();

	if ( (t_pvocFrameLen <= 0) || (t_frameLen <= 0) || (t_pvocHop <= 0) )
		throw initializationError();

	/* Make sure that pvocFrameLen is power of 2 */
	if ( !check_power2(t_pvocFrameLen) )
		throw initializationError();

	/* Make sure that pvocFrameLen is a multiple of pvocHop */
	if ( t_pvocFrameLen & t_pvocHop != 0 )
		throw initializationError();
	
	internalBufLen = pvocFrameLen * 2;

	/* Windowing function */
	hWin = new dtype[pvocFrameLen];
	for(int i = 0; i < pvocFrameLen; ++i){
		hWin[i] = 0.5 * cos((2 * M_PI * static_cast<dtype>(i)) / static_cast<dtype>(pvocFrameLen)); 
		hWin[i] = 0.5 - hWin[i];
	}

	xFrameW = new dtype[pvocFrameLen];

	/* FFT coefficients */
	fftc = new dtype[pvocFrameLen * 2];
	gen_w_r2(fftc, pvocFrameLen);

	/* FFT intermediate data fields */
	ftBuf1 = new dtype[pvocFrameLen * 2];
	ftBuf2 = new dtype[pvocFrameLen * 2];

	lastPhase = new dtype[pvocFrameLen]; /* TODO: Check size bound is tight */
	lastPhase_nps = new dtype[pvocFrameLen];
	lastPhase_ntw = new dtype[pvocFrameLen];
	sumPhase = new dtype[pvocFrameLen];

	outFrameBufPV = new dtype[internalBufLen];

	X_magn = new dtype[pvocFrameLen];
	X_phase = new dtype[pvocFrameLen];
	anaMagn = new dtype[pvocFrameLen];
	anaFreq = new dtype[pvocFrameLen];
	synMagn = new dtype[pvocFrameLen];
	synFreq = new dtype[pvocFrameLen];

	outFrameBuf = new dtype[internalBufLen];

	/* Internal variables */
	expct = 2.* M_PI * pvocHop / pvocFrameLen; /* Expected phase change */
	osamp = pvocFrameLen / pvocHop; /* Oversampling factor */
	freqPerBin = sr / pvocFrameLen; /* Frequency per bin */

	if (mode == TIME_WARP_ONLY || 
		mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		maxNDelayFrames = pvocFrameLen * 8;	/* TODO: Check that the delay doesn't exceed this limit */

		warpCacheMagn = new dtype[maxNDelayFrames * pvocFrameLen];
		warpCachePhase = new dtype[maxNDelayFrames * pvocFrameLen];
	}
	else {
		maxNDelayFrames = 0;

		warpCacheMagn = 0;
		warpCachePhase = 0;
	}

	if (mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		fixedPitchShiftST = static_cast<dtype>(0xffffffff); /* NAN */
	}

	reset();
}

/* Configure parameters */
void PhaseVocoder :: config(operMode t_operMode, 
							const int t_nDelay, 
							const dtype t_sr, 
							const int t_frameLen, 
							const int t_pvocFrameLen, 
							const int t_pvocHop) 
{
	/* Clean up */
	cleanup();

	mode = t_operMode;
	nDelay = t_nDelay; 
	sr = t_sr;
	frameLen = t_frameLen;
	pvocFrameLen = t_pvocFrameLen;
	pvocHop = t_pvocHop;

	/* Input sanity checks */
	if ( sr <= 0 )
		throw initializationError();

	if ( (t_pvocFrameLen <= 0) || (t_frameLen <= 0) || (t_pvocHop <= 0) )
		throw initializationError();

	/* Make sure that pvocFrameLen is power of 2 */
	if ( !check_power2(t_pvocFrameLen) )
		throw initializationError();

	/* Make sure that pvocFrameLen is a multiple of pvocHop */
	if ( t_pvocFrameLen & t_pvocHop != 0 )
		throw initializationError();
	
	internalBufLen = pvocFrameLen * 2;

	/* Windowing function */
	hWin = new dtype[pvocFrameLen];
	for(int i = 0; i < pvocFrameLen; ++i){
		hWin[i] = 0.5 * cos((2 * M_PI * static_cast<dtype>(i)) / static_cast<dtype>(pvocFrameLen)); 
		hWin[i] = 0.5 - hWin[i];
	}

	xFrameW = new dtype[pvocFrameLen];

	/* FFT coefficients */
	fftc = new dtype[pvocFrameLen * 2];
	gen_w_r2(fftc, pvocFrameLen);

	/* FFT intermediate data fields */
	ftBuf1 = new dtype[pvocFrameLen * 2];
	ftBuf2 = new dtype[pvocFrameLen * 2];

	lastPhase = new dtype[pvocFrameLen]; /* TODO: Check size bound is tight */
	lastPhase_nps = new dtype[pvocFrameLen];
	lastPhase_ntw = new dtype[pvocFrameLen];
	sumPhase = new dtype[pvocFrameLen];

	outFrameBufPV = new dtype[internalBufLen];

	X_magn = new dtype[pvocFrameLen];
	X_phase = new dtype[pvocFrameLen];
	anaMagn = new dtype[pvocFrameLen];
	anaFreq = new dtype[pvocFrameLen];
	synMagn = new dtype[pvocFrameLen];
	synFreq = new dtype[pvocFrameLen];

	outFrameBuf = new dtype[internalBufLen];

	/* Internal variables */
	expct = 2.* M_PI * pvocHop / pvocFrameLen; /* Expected phase change */
	osamp = pvocFrameLen / pvocHop; /* Oversampling factor */
	freqPerBin = sr / pvocFrameLen; /* Frequency per bin */

	if (mode == TIME_WARP_ONLY || 
		mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		maxNDelayFrames = pvocFrameLen * 8;	/* TODO: Check that the delay doesn't exceed this limit */

		warpCacheMagn = new dtype[maxNDelayFrames * pvocFrameLen];
		warpCachePhase = new dtype[maxNDelayFrames * pvocFrameLen];
	}
	else {
		maxNDelayFrames = 0;

		warpCacheMagn = 0;
		warpCachePhase = 0;
	}

	if (mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		fixedPitchShiftST = static_cast<dtype>(0xffffffff); /* NAN */
	}

	reset();
}

void PhaseVocoder :: cleanup() {
	if (hWin)	{ delete [] hWin; hWin = 0; };
	if (xFrameW){ delete [] xFrameW; xFrameW = 0; };

	if (fftc)	{ delete [] fftc; fftc = 0; }
	if (ftBuf1)	{ delete [] ftBuf1; ftBuf1 = 0; }
	if (ftBuf2)	{ delete [] ftBuf2; ftBuf2 = 0; }

	if (lastPhase)		{ delete [] lastPhase; lastPhase = 0; }
	if (lastPhase_nps)	{ delete [] lastPhase_nps; lastPhase_nps = 0; }
	if (lastPhase_ntw)	{ delete [] lastPhase_ntw; lastPhase_ntw = 0; }
	if (sumPhase)		{ delete [] sumPhase; sumPhase = 0; }
	if (outFrameBufPV)	{ delete [] outFrameBufPV; outFrameBufPV = 0; }

	if (X_magn)			{ delete [] X_magn; X_magn = 0; }
	if (X_phase)		{ delete [] X_phase; X_phase = 0; }
	if (anaMagn)		{ delete [] anaMagn; anaMagn = 0; }
	if (anaFreq)		{ delete [] anaFreq; anaFreq = 0; }
	if (synMagn)		{ delete [] synMagn; synMagn = 0; }
	if (synFreq)		{ delete [] synFreq; synFreq = 0; }

	if (warpCacheMagn)	{ delete [] warpCacheMagn; warpCacheMagn = 0; }
	if (warpCachePhase)	{ delete [] warpCachePhase; warpCachePhase = 0; }

	if (outFrameBuf)	{ delete [] outFrameBuf; outFrameBuf = 0; }
}

/* Destructor */
PhaseVocoder :: ~PhaseVocoder() {
	cleanup();
}

/* Main interface function: process a frame 
	Input arguments: 
		inBuf: Input buffer (assume length == pvocFrameLen) 
		shift: When shiftMode = PITCH_SHIFT: the amount of pitch shift, in semitones;
			   When shiftMode = TIME_WARP:   the amount of timing delay, in s.

*/
void PhaseVocoder :: procFrame(const dtype * inBuf, const dtype shift) {
	/* Declaration of local variables */
	dtype ms_in = 0.0;
	dtype p_tmp, magn, phase, dp;
	int qpd, index;

	/* Multiply the input frame with the window */
	DSPF_dp_vecmul(inBuf, hWin, xFrameW, pvocFrameLen);

	for (int i = 0; i < pvocFrameLen; ++i){
		ftBuf1[i * 2] = xFrameW[i];
		ftBuf1[i * 2 + 1] = 0.0;
	}	

	/* Amplitude calculation: mean square */
	for (int i = 0; i < pvocFrameLen; ++i) {
		ms_in += xFrameW[i] * xFrameW[i];
	}
	ms_in /= pvocFrameLen;

	/* Perform FFT */
	DSPF_dp_cfftr2(pvocFrameLen, ftBuf1, fftc, 1);
	bit_rev(ftBuf1, pvocFrameLen);

	for (int i = 0; i < pvocFrameLen; ++i) {
		X_magn[i] = 2.0 * sqrt(ftBuf1[i * 2] * ftBuf1[i * 2] + 
							   ftBuf1[i * 2 + 1] * ftBuf1[i * 2 + 1]);
		X_phase[i] = atan2(ftBuf1[i * 2 + 1], ftBuf1[i * 2]);
	}

	/* Record data into warpCache */
	int cidx0 = 0;
	if (mode == TIME_WARP_ONLY) {
		cidx0 = pvCtr;
		//cidx0 = pvCtr / (pvocHop / frameLen);

		for (int i = 0; i < pvocFrameLen; ++i){
			warpCacheMagn[(cidx0 % maxNDelayFrames) + i * maxNDelayFrames] = X_magn[i];
			warpCachePhase[(cidx0 % maxNDelayFrames) + i * maxNDelayFrames] = X_phase[i];
		}
	}

	if (mode == PITCH_SHIFT_ONLY || mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {  /* TODO: combined tWarp and pShift */
		/*------------------------*/
		/* Perform pitch shifting */
		/*------------------------*/

		/* Calculate the pitch shift ratio 
			Convert semitones to pitch ratio */
		dtype pitchShiftRatio;
		if (mode == PITCH_SHIFT_ONLY) {
			pitchShiftRatio = pow(2.0, shift / 12.0);
		}
		else {
			if (fixedPitchShiftST == static_cast<dtype>(0xffffffff)) { /* Test for NaN */
				throw fixedPitchShiftNotSpecifiedErr();
			}
			pitchShiftRatio = pow(2.0, fixedPitchShiftST / 12.0);
		}

		bPitchShift = (pitchShiftRatio != 1.0);

		if ( !bPitchShift && bPitchShift_prev ) {
			/* Recover from the last interval of pitch shifting */ 
			for (int i = 0; i <= pvocFrameLen / 2; ++i) {
				sumPhase[i] = lastPhase[i];
			}
		}

		for (int i0 = 0; i0 <= pvocFrameLen / 2; ++i0){
			/*for (int i1 = 0; i1 < 2; ++i1) {*/
			p_tmp = X_phase[i0] - lastPhase[i0];

			p_tmp -= (dtype)i0 * expct;

			qpd = static_cast<int>(p_tmp / M_PI);

			if (qpd >= 0) 
				qpd += qpd & 1;
			else 
				qpd -= qpd & 1;

			p_tmp -= M_PI * static_cast<dtype>(qpd);

			p_tmp = osamp * p_tmp / (2. * M_PI);
			p_tmp = (dtype)i0 * freqPerBin + p_tmp * freqPerBin;

			anaMagn[i0] = X_magn[i0];
			anaFreq[i0] = p_tmp;

			//if (i1 == 0) { // DEBUG: amp
			//	ss_anaMagn += anaMagn[i1][i0] * anaMagn[i1][i0];
			//} 

			lastPhase[i0] = X_phase[i0];
			//}						
		}

		/* Pitch shifting synthesis */
		for (int i0 = 0; i0 < pvocFrameLen; ++i0){
			synMagn[i0] = 0.0;
			synFreq[i0] = 0.0;
		}

		for (int i0 = 0; i0 <= pvocFrameLen / 2; ++i0){
			//for (int i1 = 0; i1 < 2; i1++) {
			//if (i1 == 0)
				index = static_cast<int>(i0 / pitchShiftRatio);
			//else
				//index[i1] = i0;
							
			if (index <= pvocFrameLen / 2) {
				synMagn[i0] += anaMagn[index];

				//if (i1 == 0)
				synFreq[i0] = anaFreq[index] * pitchShiftRatio;
				//else
					//synFreq[i1][i0] = anaFreq[i1][index[i1]];
			}
			//}
		}

		for (int i0 = 0; i0 <= pvocFrameLen / 2; ++i0) {
			//for (int i1 = 0; i1 < 2; i1++) {
			magn = synMagn[i0]; // get magnitude and true frequency from synthesis arrays						

			p_tmp = synFreq[i0];
							
			p_tmp -= static_cast<dtype>(i0) * freqPerBin;	// subtract bin mid frequency		
			p_tmp /= freqPerBin;	// get bin deviation from freq deviation
			p_tmp = 2. * M_PI * p_tmp / osamp;	// take osamp into account
			p_tmp += static_cast<double>(i0) * expct;		// add the overlap phase advance back in
							
			sumPhase[i0] += p_tmp;		// accumulate delta phase to get bin phase
			phase = sumPhase[i0];
			//sumPhase[i0] = X_phase[i0]; /* DEBUG */

			/* get real and imag part and re-interleave */
			ftBuf2[2 * i0] = magn * cos(phase);
			ftBuf2[2 * i0 + 1] = magn * sin(phase);			// What causes the sign reversal here?
			//}
		}

		/* Time warping on top of a fixed pitch shift 
		      (specified in fixedPitchShiftST) */
		if (mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
			cidx0 = pvCtr;
			//cidx0 = pvCtr / (pvocHop / frameLen);

			for (int i = 0; i < pvocFrameLen; ++i){
				warpCacheMagn[(cidx0 % maxNDelayFrames) + i * maxNDelayFrames] = synMagn[i];
				warpCachePhase[(cidx0 % maxNDelayFrames) + i * maxNDelayFrames] = sumPhase[i];
			}
		} /* TODO: Finish and test it */

		bPitchShift_prev = bPitchShift;

	}
	
	if (mode == TIME_WARP_ONLY || mode == TIME_WARP_WITH_FIXED_PITCH_SHIFT) {
		/*----------------------*/
		/* Perform time warping */
		/*----------------------*/

		dtype warp_t = shift;
		if (warp_t > 0.0)
			throw timeWarpFuturePredError();
		
		int cidx1;
		dtype cidx1_f;
		if (warp_t != 0.0) {
			/* Each PVOC hop takes a period of (pvocHop / sr) */
			dtype warp_i = 
				warp_t * sr / static_cast<dtype>(pvocHop);
			dtype cidx1_d = static_cast<dtype>(cidx0) + warp_i;
			cidx1 = static_cast<int>(floor(cidx1_d));
			cidx1_f = cidx1_d - static_cast<dtype>(cidx1);

			bWarp = true;
		}
		else {
			if (bWarp_prev) {
				 /* Recover from the last interval of time warping */ 
				for (int i = 0; i <= pvocFrameLen / 2; ++i) {
					lastPhase[i] = lastPhase_ntw[i];
				}
			}

			cidx1 = cidx0;
			cidx1_f = 0.0;

			bWarp = false;
		}

		while (cidx1 < 0)
			cidx1 += maxNDelayFrames;

		int k0, k1;
		k0 = cidx1 % maxNDelayFrames;
		k1 = (cidx1 + 1) % maxNDelayFrames;

		for (int i0 = 0; i0 <= pvocFrameLen / 2; ++i0) {
			magn = warpCacheMagn[k0 + i0 * maxNDelayFrames] * (1 - cidx1_f) + 
				   warpCacheMagn[k1 + i0 * maxNDelayFrames] * cidx1_f;
			/* No time warp */
			/*magn = pvocWarpCache[cidx0 % (internalBufLen / 64)][i0];*/

			if (!lastPhasePrimed) {
				lastPhase[i0] = warpCachePhase[k0 + i0 * maxNDelayFrames];
			}

			ftBuf2[2 * i0] = 0.5 * magn * cos(lastPhase[i0]);
			ftBuf2[2 * i0 + 1] = 0.5 * magn * sin(lastPhase[i0]);

			if (cidx1_f == 0.0) {
				dtype interpPhase 
					= warpCachePhase[k0 + i0 * maxNDelayFrames] * (1 - cidx1_f) + 
   					  warpCachePhase[k1 + i0 * maxNDelayFrames] * cidx1_f;
				dp = interpPhase - lastPhase[i0];
			}
			else {
				dp = warpCachePhase[k1 + i0 * maxNDelayFrames] -
					 warpCachePhase[k0 + i0 * maxNDelayFrames];
			}
						
			dp -= static_cast<dtype>(i0) * expct;
						
			qpd = (int)(dp / M_PI);
			if (qpd >= 0) 
				qpd += qpd & 1;
			else 
				qpd -= qpd & 1;
			dp -= M_PI * static_cast<dtype>(qpd);

			/* DEBUG */
			/*ftBuf2[2 * i0] = ftBuf1[2 * i0];
			ftBuf2[2 * i0 + 1] = ftBuf1[2 * i0 + 1];*/
			/* ~DEBUG */

			/* No time warp */
			/*ftBuf2ps[2 * i0] = magn * cos(X_phase[i0]);
			ftBuf2ps[2 * i0 + 1] = magn * sin(X_phase[i0]);*/

			lastPhase[i0] += static_cast<dtype>(i0) * expct + dp;

			/* No time warp */
			/*lastPhase[i0] = X_phase[i0];*/

			/* No-time-warp backup */
			lastPhase_ntw[i0] = X_phase[i0];
			sumPhase[i0] = X_phase[i0];
		}

		if (!lastPhasePrimed)
			lastPhasePrimed = true;

		/* Record time warp status during the current frame */
		bWarp_prev = bWarp;
	}

	for (int i0 = pvocFrameLen + 1; i0 < pvocFrameLen * 2; ++i0)
		/*for (int i1 = 0; i1 < 2; i1++) {*/
		ftBuf2[i0] = 0.0;
		//}

	/* Inverse Fourier transform */				
	//for (int i1 = 0; i1 < 2; i1 ++ ) {
	DSPF_dp_icfftr2(pvocFrameLen, ftBuf2, fftc, 1);
	bit_rev(ftBuf2, pvocFrameLen);
	for (int i0 = 0; i0 < pvocFrameLen; ++i0){
		ftBuf2[i0 * 2] /= pvocFrameLen;
		ftBuf2[i0 * 2 + 1] /= pvocFrameLen;
	}
	//}

	/* Scaling */
	for (int i = 0; i < pvocFrameLen; ++i)
		ftBuf2[2 * i] = 2 * ftBuf2[2 * i] * hWin[i] / (osamp / 2.0);

	///* Accumulate to output buffer */
	//for (int i = 0; i < pvocFrameLen; ++i){
	//	outFrameBuf[(outFrameBuf_circPtr + i) % internalBufLen] += 
	//		ftBuf2[2 * i];
	//}			

	//// Front zeroing
	//for (int i = 0; i < pvocHop; ++i)
	//	outFrameBuf[(outFrameBuf_circPtr + pvocFrameLen + i) % internalBufLen] = 0.0;
	//			
	//// Back Zeroing
	//for (int i0 = 1; i0 <= p.pvocHop; i0++){ // Ad hoc alert!
	//	outFrameBuf[(outFrameBuf_circPtr - i0) % (internalBufLen)] = 0.;						
	//}

	pvCtr++;
}

/* Status reset */
void PhaseVocoder :: reset() {
	/* Related to pitch shifting */
	for (int i0 = 0; i0 < pvocFrameLen * 2; i0++) {
		ftBuf1[i0] = 0.0;
		ftBuf2[i0] = 0.0;
	}

	for (int i0 = 0; i0 < pvocFrameLen; i0++) {
		lastPhase[i0] = 0.0;
		lastPhase_nps[i0] = 0.0;
		lastPhase_ntw[i0] = 0.0;

		sumPhase[i0] = 0.0;
	}

	for (int i = 0; i < internalBufLen; ++i) {
		outFrameBuf[i] = 0.0;
	}
	outFrameBuf_circPtr = 0;

	if (maxNDelayFrames > 0) {
		for (int i = 0; i < maxNDelayFrames * pvocFrameLen; ++i) {
			warpCacheMagn[i] = 0.0;
			warpCachePhase[i] = 0.0;
		}
	}

	bPitchShift = false;
	bPitchShift_prev = false;

	bWarp = false;
	bWarp_prev = false;	
	lastPhasePrimed = false;

	pvCtr = 0;
}

const bool check_power2(const int & x) {
  /* Number of bits */
  int nBits = sizeof(int) * 8;

  int nnzb = 0; /* Number of non-zero bits */
  unsigned int y = x;

  for (int i = 0; i < nBits; ++i) {      
    nnzb += (y & 1);
    y >>= 1;
  }

  return (nnzb == 1);
}


const PhaseVocoder::operMode PhaseVocoder::getMode() const {
	return mode;
}

void PhaseVocoder::setFixedPitchShiftST(const dtype t_fixedPitchShiftST) {
	fixedPitchShiftST = t_fixedPitchShiftST;
}
