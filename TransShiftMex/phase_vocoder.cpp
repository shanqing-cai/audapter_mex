#include <cmath>

#include "phase_vocoder.h"
#include "DSPF.h"

/* Constructor */
PhaseVocoder :: PhaseVocoder(const dtype t_sr, 
							 const int t_frameLen, 
							 const int t_pvocFrameLen, 
							 const int t_pvocHop) 
	//: sr(t_sr), frameLen(t_frameLen), 
	//  pvocFrameLen(t_pvocFrameLen), pvocHop(t_pvocHop)
{
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



	reset();
}

/* Destructor */
PhaseVocoder :: ~PhaseVocoder() {
	if (hWin)	delete [] hWin;
	if (xFrameW)	delete [] xFrameW;

	if (fftc)	delete [] fftc;
	if (ftBuf1)	delete [] ftBuf1;
	if (ftBuf2)	delete [] ftBuf2;

	if (lastPhase)		delete [] lastPhase;
	if (lastPhase_ntw)	delete [] lastPhase_ntw;
	if (sumPhase)		delete [] sumPhase;
	if (outFrameBufPV)	delete [] outFrameBufPV;

	if (X_magn)			delete [] X_magn;
	if (X_phase)		delete [] X_phase;
	if (anaMagn)		delete [] anaMagn;
	if (anaFreq)		delete [] anaFreq;
	if (synMagn)		delete [] synMagn;
	if (synFreq)		delete [] synFreq;

	if (outFrameBuf)	delete [] outFrameBuf;
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
	dtype p_tmp, magn, phase;
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

	/* Calculate the pitch shift ratio 
	      Convert semitones to pitch ratio */
	dtype pitchShiftRatio = pow(2.0, shift / 12.0);

	for (int i0 = 0; i0 <= pvocFrameLen / 2; ++i0){
		/*for (int i1 = 0; i1 < 2; ++i1) {*/
		p_tmp = X_phase[i0] - lastPhase[i0];

		p_tmp -= (dtype)i0 * expct;

		qpd = static_cast<int>(p_tmp / M_PI);

		if (qpd >= 0) 
			qpd += qpd&1;
		else 
			qpd -= qpd&1;

		p_tmp -= M_PI * (dtype)qpd;

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

		/* get real and imag part and re-interleave */
		ftBuf2[2 * i0] = magn * cos(phase);
		ftBuf2[2 * i0 + 1] = magn * sin(phase);			// What causes the sign reversal here?
		//}
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
		lastPhase_ntw[i0] = 0.0;

		sumPhase[i0] = 0.0;
	}

	for (int i = 0; i < internalBufLen; ++i) {
		outFrameBuf[i] = 0.0;
	}
	outFrameBuf_circPtr = 0;
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
