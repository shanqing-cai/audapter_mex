/* 
	lpc_formant.cpp
	Linear prediction and formant tracking 

	Shanqing Cai 12/2013

*/

#include <cmath>

#include "lpc_formant.h"
#include "DSPF.h"

/* Utility inline functions */
inline dtype mul_sign(const dtype &a, const dtype &b) {
	return (b >= 0.0 ? fabs(a) : -fabs(a));
}

inline int imax(const int &k, const int &j) {
	return (k <= j ? j : k);
}

/* Constructor
	
	Input arguments:
		t_nLPC: LPC order
		t_sr: sampling rate of the input signal (Hz)
		t_bufferSize: input buffer size (# of samples)
		t_nFFT: FFT length (must be power of 2)
		t_cepsWinWidth: cepstral liftering window width (dimensionless)
		t_nTracks: Number of formants to be tracked
		t_aFact: alpha parameter in the DP formant tracking algorithm (Xia and Espy-Wilson, 2000)
		t_bFact: beta
		t_gFact: gamma
		t_fn1: prior value of F1 (Hz)
		t_fn2: prior value of F2 (Hz)

		If cepsWinWidth is <= 0, cepstral liftering will be disabled.
*/
LPFormantTracker::LPFormantTracker(const int t_nLPC, const int t_sr, const int t_bufferSize, 					 
								   const int t_nFFT, const int t_cepsWinWidth, 
								   const int t_nTracks, 
								   const dtype t_aFact, const dtype t_bFact, const dtype t_gFact, 
			   					   const dtype t_fn1, const dtype t_fn2, 
								   const bool t_bMWA, const int t_avgLen) :
	nLPC(t_nLPC), sr(t_sr), 
	bufferSize(t_bufferSize), 
	nFFT(t_nFFT), cepsWinWidth(t_cepsWinWidth), 
	nTracks(t_nTracks), 
	aFact(t_aFact), bFact(t_bFact), gFact(t_gFact), 
	fn1(t_fn1), fn2(t_fn2), 
	bMWA(t_bMWA), avgLen(t_avgLen)	
{
	/* Input sanity checks */
	if ( (nLPC <= 0) || (bufferSize <= 0) || (nFFT <= 0) ||
	     (nTracks <= 0) )
		 throw initializationError();

	if (nLPC > maxNLPC)
		throw nLPCTooLargeError();

	bCepsLift = (cepsWinWidth > 0);
	
	winFunc = new dtype[bufferSize];

	/* Initialize window */
	genHanningWindow();

	nLPC_SQR = nLPC * nLPC;

	Acompanion = new dtype[nLPC_SQR];
	AHess = new dtype[nLPC_SQR];

	temp_frame = new dtype[bufferSize + nLPC + 2]; // +2 for playing it safe! TODO: improve
	R = new dtype[nLPC * 2]; // *2 for playing it safe! TODO: improve

	/* Initalize FFT working date fields */
	ftBuf1 = new dtype[nFFT * 2];
	ftBuf2 = new dtype[nFFT * 2];
	fftc = new dtype[nFFT * 2];

	gen_w_r2(fftc, nFFT);

	/* Intermediate data fields */
	lpcAi = new dtype[maxNLPC + 1];
	realRoots = new dtype[maxNLPC];
	imagRoots = new dtype[maxNLPC];

	cumMat = new dtype[maxFmtTrackJump * maxNTracks];
	costMat = new dtype[maxFmtTrackJump * maxNTracks];

	nCands = 6; 
	/* number of possible formant candiates  
	   (should be > ntracks but  < p.nLPC/2!!!! (choose carefully : not fool-proof) 
	   TODO: Implement automatic checks */

	weiMatPhi = new dtype[nTracks * maxAvgLen];
	weiMatBw = new dtype[nTracks * maxAvgLen];
	weiVec = new dtype[maxAvgLen];
	sumWeiPhi = new dtype[nTracks];
	sumWeiBw = new dtype[nTracks];

	trackFF = 0.95;

	radius_us = new dtype[maxNLPC]; /* TODO: Tighten the size */
	phi_us = new dtype[maxNLPC];
	bandwidth_us = new dtype[maxNLPC];
	//phi_s = new dtype[maxNLPC];

	/* Call reset */
	reset();
}

/* Destructor */
LPFormantTracker::~LPFormantTracker() 
{
	if (winFunc)
		delete [] winFunc;

	if (Acompanion)
		delete [] Acompanion;
	if (AHess)
		delete [] AHess;

	if (temp_frame)
		delete [] temp_frame;
	if (R)
		delete [] R;

	if (ftBuf1)
		delete [] ftBuf1;
	if (ftBuf2)
		delete [] ftBuf2;
	if (fftc)
		delete [] fftc;

	if (lpcAi)
		delete [] lpcAi;
	if (realRoots)
		delete [] realRoots;
	if (imagRoots)
		delete [] imagRoots;

	if (cumMat)
		delete [] cumMat;
	if (costMat)
		delete [] costMat;

	if (weiMatPhi)
		delete [] weiMatPhi;
	if (weiMatBw)
		delete [] weiMatBw;
	if (weiVec)
		delete [] weiVec;
	if (sumWeiPhi)
		delete [] sumWeiPhi;
	if (sumWeiBw)
		delete [] sumWeiBw;

	if (radius_us)
		delete [] radius_us;
	if (phi_us)
		delete [] phi_us;
	if (bandwidth_us)
		delete [] bandwidth_us;
	/*if (phi_s)
		delete [] phi_s;*/
}

/* Reset after a supra-threshold interval */
void LPFormantTracker::postSupraThreshReset() {
	for (int i = 0; i < maxNLPC; ++i) {
		realRoots[i] = 0.0;
		imagRoots[i] = 0.0;
	}

	for(int j0 = 0; j0 < nTracks; ++j0)
	{
		sumWeiPhi[j0] = 0.0;
		sumWeiBw[j0] = 0.0;
		for(int i0 = 0; i0 < maxAvgLen; ++i0)
		{
			weiMatPhi[j0 + nTracks * i0] = 0.0;
			weiMatBw[j0 + nTracks * i0] = 0.0;
		}
	}

	for(int i0 = 0; i0 < maxAvgLen; ++i0)
		weiVec[i0]=0.0;
	//weiVec[mwaCircCtr] = 0.0;
	sumWei = 0.0;

	mwaCtr = 0;
	mwaCircCtr = 0;
}

/* Resetting */
void LPFormantTracker::reset() {
	for (int i = 0; i < nFFT * 2; ++i) {
		ftBuf1[i] = 0.0;
		ftBuf2[i] = 0.0;
	}

	for (int i = 0; i < nLPC_SQR; ++i) {
		Acompanion[i] = 0.0;
		AHess[i] = 0.0;
	}

	for (int i = nLPC - 2; i >= 0; --i) {
		Acompanion[(nLPC + 1) * i + 1] = 1.0;
		AHess[(nLPC + 1) * i + 1] = 1.0;
	}


	for(int j0 = 0; j0 < nTracks; ++j0)
	{
		sumWeiPhi[j0] = 0.0;
		sumWeiBw[j0] = 0.0;
		for(int i0 = 0; i0 < maxAvgLen; ++i0)
		{
			
			weiMatPhi[j0 + nTracks * i0] = 0.0;
			weiMatBw[j0 + nTracks * i0] = 0.0;
		}
	}

	for(int i0 = 0; i0 < maxAvgLen; ++i0)
		weiVec[i0]=0.0;

	sumWei = 0.0;

	postSupraThreshReset();
}

/* Levinson recursion for linear prediction (LP) */
void LPFormantTracker::levinson(dtype * R, dtype * aa, const int size) {
	dtype ki, t;
    dtype E = R[0];
    int   i0, j0;

	if (R[0] == 0.0) {
		for(i0=1; i0<size; i0++)
			aa[i0] = 0.0;

		aa[0] = 1;
		return;
	}

    for(i0=1; i0<size; i0++) { 
        ki = R[i0]; 
      
        // Update reflection coefficient: 
        for (j0=1; j0<i0; j0++) 
			ki += aa[j0] * R[i0-j0]; 
      
        ki   /= -E; 
        E    *= (1 - ki*ki); 
        
        // Update polynomial: 
        for (j0 = 1; j0 <= RSL(i0 - 1, 1); j0++) {
            t = aa[j0];
            aa[j0] += ki * aa[i0 - j0];
            aa[i0 - j0] += ki * t; 
        } 
  
        if (i0%2 == 0) aa[RSL(i0, 1)] *= 1+ki; 
  
        // Record reflection coefficient
        aa[i0] = ki; 

    } // end of for loop

    aa[0] = 1.0;
}


/* The following takes in a polynomial stored in *c, and yields the roots of
   this polynomial (*wr stores the real comp, and *wi stores the imag comp)
   It forms a companion matrix, then uses the hqr algorithm 
   VV 19 June 2003 
	
	Input arguments:
		c: coefficients of the polynomial
		wr: real parts of the roots (output)
		wi: imaginary parts of the roots (output)

*/
int LPFormantTracker::hqr_roots(dtype * c, dtype * wr, dtype * wi) {
#ifndef aMat
	#define aMat(k, j) AHess[((j) - 1) * nLPC + (k) - 1]
#endif

	int nn, m, l, k, j, i, its, mmin, nLPC_SQR = nLPC * nLPC;

    /* dtype AHess[maxNLPC_squared]; */
    dtype z, y, x, w, v, u, t, s, r, q, p, anorm = 0.0F;
        
	/* generate companion matrix, starting off with an intialized version */
	DSPF_dp_blk_move(Acompanion, AHess, nLPC_SQR);

	for (i = 0; i < nLPC; i++)
		AHess[nLPC * i] = -c[i + 1];

    /* end of companion matrix generation  */

    /* the following performs the hqr algoritm  */
    /* NOTE:  This was taken from Numerical Recipes in C, with modification
       Specifically, the wr and wi arrays were assumed to number from 1..n
       in the book.  I modified calls to these arrays so that they number 0..n-1
       Additionally, n (the order of the polynomial) is hardset to be 8 
       VV 19 June 2003 */

    for (i = 1; i < nLPC + 1; i++)
        for (j = imax(i - 1, 1); j < nLPC + 1; j++)
            anorm += fabs(aMat(i, j)); 
			/*anorm += fabs(AHess[(j - 1) * nLPC + i - 1]);*/

    nn = nLPC;
    t = 0.0;
    while (nn >= 1) {
		its=0;
        do {
           for (l = nn; l >= 2; l--) {
               s = fabs(aMat(l - 1, l - 1)) + fabs(aMat(l, l));
			   /*s = fabs(AHess[(l - 2) * nLPC + l - 2]) + fabs(AHess[(l - 1) * nLPC + l - 1]);*/

               if (s == 0.0) s = anorm;

               if ((dtype) (fabs(aMat(l, l - 1)) + s) == s) 
			   /*if ((dtype) (fabs(AHess[(l - 2) * nLPC + l - 1]) + s) == s)*/
				   break;
           }

         x=aMat(nn, nn);
		 /*x = AHess[(nn - 1) * nLPC + nn - 1];*/

         if (l == nn) {
			wr[(-1) + nn] = x + t;
			wi[(-1) + nn--] = 0.0;
		 }
		 else {
			y = aMat(nn - 1 ,nn - 1);
			/*y = AHess[(nn - 2) * nLPC + nn - 2];*/

            w = aMat(nn, nn - 1) * aMat(nn - 1, nn);
			/*w = AHess[(nn - 2) * nLPC + nn - 1] * AHess[(nn - 2) * nLPC + nn - 1];*/
                    
			if (l == (nn-1)) {
				p = 0.5 * (y - x);
                q = p * p + w;
                z=sqrt(fabs(q));
                x += t;
                
				if (q >= 0.0) {
					z = p + mul_sign(z, p);
					wr[(-1) + nn - 1] = wr[(-1) + nn] = x + z;
                    if (z) wr[(-1) + nn] = x - w / z;
                    wi[(-1) + nn - 1] = wi[(-1) + nn] = 0.0;
                } 
				else {
					wr[(-1) + nn - 1] = wr[(-1) + nn] = x + p;
                    wi[(-1) + nn - 1] = -(wi[(-1) + nn] = z);
                }
				nn -= 2;
			} 
			else {
				if (its == 10 || its == 20) {
					t += x;
                    for (i = 1; i <= nn; i++) 
						aMat(i, i) -= x;
						/*AHess[(i - 1) * nLPC + i - 1] -= x;*/

                    s = fabs(aMat(nn, nn - 1)) + fabs(aMat(nn - 1, nn - 2));
					/*s = fabs(AHess[(nn - 2) * nLPC + nn - 1]) + fabs(AHess[(nn - 3) * nLPC + nn - 2]);*/

					y = x = 0.75 * s;
                    w = -0.4375 * s * s;
                }
				
				++its;
                for (m = nn - 2; m >= l; m--) {
					z = aMat(m, m);
					/*z = AHess[(m - 1) * nLPC + m - 1];*/

					r = x - z;
					s = y - z;
					
					p = (r * s - w) / aMat(m + 1, m) + aMat(m, m + 1);
					/*p = (r * s - w) /  AHess[(m - 1) * nLPC + m] + AHess[m * nLPC + m - 1];*/

					q = aMat(m + 1, m + 1) - z - r - s;
					/*q = AHess[m * nLPC + m] - z - r - s;*/

					r = aMat(m + 2, m + 1);
					/*r = AHess[m * nLPC + m + 1];*/

					s = fabs(p) + fabs(q) + fabs(r);
					p /= s;
					q /= s;
					r /= s;
					if (m == l) break;

					u = fabs(aMat(m, m - 1)) * (fabs(q) + fabs(r));
					/*u = fabs(AHess[(m - 2) * nLPC + m - 1]) * (fabs(q) + fabs(r));*/

					v = fabs(p) * (fabs(aMat(m - 1, m - 1)) + fabs(z) + fabs(aMat(m + 1, m + 1)));
					/*v = fabs(p) * (fabs(AHess[(m - 2) * nLPC + m - 2]) + fabs(z) + fabs(AHess[m * nLPC + m]));*/

					if ((dtype) (u+v) == v) break;
                }
                
				for (i = m + 2; i <= nn; i++) {
					aMat(i, i - 2) = 0.0F;
					//AHess[(i - 3) * nLPC + i - 1] = 0.0F;

                    if (i != (m + 2)) 
						aMat(i, i - 3) = 0.0F;
						/*AHess[(i - 4) * nLPC + i - 1] = 0.0F;*/
                }

				for (k = m; k <= nn - 1; k++) {
					if (k != m) {
						p = aMat(k, k - 1);
						/*p = AHess[(k - 2) * nLPC + k - 1];*/

                        q = aMat(k + 1,k - 1);
						/*p = AHess[(k - 2) * nLPC + k];*/

                        r = 0.0F;
                        if (k != (nn - 1)) 
							r = aMat(k + 2, k - 1);
							/*r = AHess[(k - 2) * nLPC + k + 1];*/

                        if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
                            p /= x;
                            q /= x;
                            r /= x;
                        }
					}
                    
					if ((s = mul_sign(sqrt(p * p + q * q + r * r), p)) != 0.0) {
						if (k == m) {
                            if (l != m)
                            aMat(k,k-1) = -aMat(k,k-1);
							/*AHess[(k - 2) * nLPC + k - 1] *= -1.0;*/
                        } 
						else
                            aMat(k, k - 1) = -s * x;
							/*AHess[(k - 2) * nLPC + k - 1] = -s * x;*/

                        p += s;
                        x=p/s;
                        y=q/s;
                        z=r/s;
                        q /= p;
                        r /= p;

                        for (j=k;j<=nn;j++) {
                            p=aMat(k,j)+q*aMat(k+1,j);
							/*p = AHess[(j - 1) * nLPC + k - 1] + q * AHess[(j - 1) * nLPC + k];*/

                            if (k != (nn-1)) {
                                p += r * aMat(k + 2, j);
								/*p += r * AHess[(j - 1) * nLPC + k + 1];*/

                                aMat(k + 2, j) -= p * z;
							/*	AHess[(j - 1) * nLPC + k + 1] -= p * z;*/
                            }
                            aMat(k+1,j) -= p*y;
							/*AHess[(j - 1) * nLPC + k] -= p * y;*/

                            aMat(k,j) -= p*x;
							/*AHess[(j - 1) * nLPC + k - 1] -= p * x;*/
                        }
                        mmin = nn<k+3 ? nn : k+3;
                        for (i=l;i<=mmin;i++) {
                            p = x * aMat(i, k) + y * aMat(i, k + 1);
							/*p = x * AHess[(k - 1) * nLPC + i - 1] + y * AHess[k * nLPC + i - 1];*/

                            if (k != (nn-1)) {
                                p += z*aMat(i,k+2);
								/*p += z * AHess[(k + 1) * nLPC + i - 1];*/

                                aMat(i,k+2) -= p*r;
								/*AHess[(k + 1) * nLPC + i - 1] -= p * r;*/
                            }

                            aMat(i,k+1) -= p*q;
							/*AHess[k * nLPC + i - 1] -= p * q;*/

                            aMat(i,k) -= p;
							/*AHess[(k - 1) * nLPC + i - 1] -= p;*/
                        }
                    }
                }
			}
		}
		
		} 
	while (l < nn - 1);
	}

	if (nn == 0) 
		return 1;
	else
		return 0;
}


/* Autocorrelation-based LPC analysis
	Performs the LPC analysis on a given frame... returns the lpc coefficients
	Input arguments
	xx: input buffer pointer
	aa: LPC coefficients pointer (output)
	size: size of the input buffer (xx)
	nlpc: LPC order, i.e., number of LPC coefficients
	SC (2008/05/06): Incoroporating cepstral lifting */
void LPFormantTracker::getAi(dtype* xx, dtype* aa) {
	int i0;	

	//dtype temp_frame[maxBufLen + maxNLPC];	// Utility buffer for various filtering operations
	//dtype R[maxNLPC];					// lpc fit Autocorrelation estimate
	
	int nlpcplus1 = nLPC + 1;

	for(i0 = 0; i0 < nlpcplus1; i0++)
		temp_frame[i0] = 0;
	
	// Window input
	//SC vecmu1: vector multiply
	//SC	hwin: a Hanning window
	DSPF_dp_vecmul(xx, winFunc, &temp_frame[nlpcplus1], bufferSize);	// Apply a Hanning window
	// TODO: Check temp_frame size

	////////////////////////////////////////////////////
	//SC(2008/05/07) ----- Cepstral lifting -----
	if (bCepsLift){
		for (i0=0;i0<nFFT;i0++){
			if (i0 < bufferSize){
				ftBuf1[i0 * 2] = temp_frame[nlpcplus1 + i0];
				ftBuf1[i0 * 2 + 1]=0;
			}
			else{
				ftBuf1[i0 * 2] = 0;
				ftBuf1[i0 * 2 + 1] = 0;
			}
		}
		DSPF_dp_cfftr2(nFFT, ftBuf1, fftc, 1);	
		bit_rev(ftBuf1, nFFT);
		// Now ftBuf1 is X
		for (i0 = 0; i0 < nFFT; i0++){
			if (i0 <= nFFT / 2){
				ftBuf2[i0 * 2] = log(sqrt(ftBuf1[i0*2]*ftBuf1[i0*2]+ftBuf1[i0*2+1]*ftBuf1[i0*2+1]));	// Optimize
				ftBuf2[i0 * 2 + 1]=0;
			}
			else{
				ftBuf2[i0 * 2] = ftBuf2[(nFFT - i0) * 2];
				ftBuf2[i0 * 2 + 1]=0;
			}
		}
		DSPF_dp_icfftr2(nFFT, ftBuf2, fftc, 1);
		bit_rev(ftBuf2, nFFT);

		// Now ftBuf2 is Xceps: the cepstrum
		for (i0 = 0; i0 < nFFT; i0++){
			if (i0 < cepsWinWidth || i0 > nFFT - cepsWinWidth){			// Adjust! 
				ftBuf1[i0 * 2] = ftBuf2[i0 * 2] / nFFT;		// Normlize the result of the previous IFFT
				ftBuf1[i0 * 2 + 1] = 0;
			}
			else{
				ftBuf1[i0 * 2] = 0;
				ftBuf1[i0 * 2 + 1] = 0;
			}
		}
		// Now ftBuf1 is Xcepw: the windowed cepstrum
		DSPF_dp_cfftr2(nFFT,ftBuf1,fftc,1);
		bit_rev(ftBuf1,nFFT);
		for (i0=0;i0<nFFT;i0++){
			if (i0<=nFFT/2){
				ftBuf2[i0*2]=exp(ftBuf1[i0*2]);
				ftBuf2[i0*2+1]=0;
			}
			else{
				ftBuf2[i0*2]=ftBuf2[(nFFT-i0)*2];
				ftBuf2[i0*2+1]=0;
			}
		}
		DSPF_dp_icfftr2(nFFT,ftBuf2,fftc,1);	// Need normalization
		bit_rev(ftBuf2,nFFT);
		
		for (i0 = 0; i0 < bufferSize / 2; i0++){
			temp_frame[nlpcplus1 + bufferSize / 2 + i0] = ftBuf2[i0 * 2] / nFFT;
		}
		for (i0 = 1; i0 < bufferSize / 2; i0++){
			temp_frame[nlpcplus1 + bufferSize / 2 - i0] = ftBuf2[i0 * 2] / nFFT;
		}
		temp_frame[nlpcplus1] = 0.0;
	}
	//~SC(2008/05/07) ----- Cepstral lifting -----
	///////////////////////////////////////////////////////////

	// Find autocorrelation values
	DSPF_dp_autocor(R, temp_frame, bufferSize, nlpcplus1);		//SC Get LPC coefficients by autocorrelation

	// Get unbiased autocorrelation
	for(i0 = 0; i0 < nlpcplus1; i0++)
		R[i0] /= bufferSize;

	// levinson recursion
	levinson(R, aa, nlpcplus1);
}

/* Get the angle (Phi) and magnitude (Bw) of the roots 
	Input argments:
	wr: real part of roots
	wi: imag part of roots
	radius: root radii (output)
	phi: root angle (output)
	bandwidth: root bandwidth (output) */
void LPFormantTracker::getRPhiBw(dtype * wr, dtype * wi, 
								 dtype * radius, dtype * phi, dtype * bandwidth) {
  /* The following sorts the roots in wr and wi.  It is adapted from a matlab script that Reiner wrote */
	/*const int maxNLPC = 64;*/

	//if (nLPC > maxNLPC) 
	//	mexErrMsgTxt("getRPhiBw: nLPC too large");
	dtype arc[maxNLPC], arc2[maxNLPC], wr2[maxNLPC], wi2[maxNLPC];
	dtype wreal, wimag, warc, wmag, mag[maxNLPC], mag2[maxNLPC];
	int numroots, i0, j0, nmark;

	/* calculate angles for all defined roots */
	numroots = 0;

	for (i0=0; i0 < nLPC; i0++) {
		arc[i0] = atan2(wi[i0],wr[i0]);
		mag[i0] = sqrt(wi[i0]*wi[i0] + wr[i0]*wr[i0]);
	    if  ( /*(arc[i0] > F1_min) && */ (wi[i0]>0) /*&& 
    	 (mag[i0] > 0.9) && (mag[i0] < 1.0) */ )
        /* only store positive arc root of conjugate pairs */
        {
            mag2[numroots]	 = mag[i0];
            arc2[numroots]   = arc[i0];  
            wr2[numroots]    = wr[i0];
            wi2[numroots++]  = wi[i0];
		}
	}

	/* sort according to arc using a stupid sort algorithm. */
	for (i0=0; i0<numroots; i0++)  /* look for minimal first */
	{
		nmark = i0;
		for (j0=i0+1; j0<numroots; j0++)  /* find smallest arc (frequency) */
			if (arc2[j0] < arc2[nmark]) nmark = j0;
		if (nmark != i0) /* switch places if smaller arc */
        {
			wreal = wr2[i0];
            wimag = wi2[i0];
            warc  = arc2[i0];
            wmag  = mag2[i0];
            wr2[i0] = wr2[nmark];
            wi2[i0] = wi2[nmark];
            arc2[i0] = arc2[nmark];
            mag2[i0] = mag2[nmark];
            wr2[nmark] = wreal;
            wi2[nmark] = wimag;
            arc2[nmark] = warc;
            mag2[nmark] = wmag;
        }
    }

	for (i0=0; i0<numroots; i0++) {
		radius[i0]=mag2[i0];
		bandwidth[i0] = -log(mag2[i0]) * dtype(sr) / M_PI;
		phi[i0]=arc2[i0];
	}
}


/* Dynamic programming based formant tracking (Xia and Espy-Wilson, 2000, ICSLP) 
	Input: r_ptr: array of amplitudes of the roots
		   phi_ptr: array of phase angles of the roots
		   
		   In-place operations are done on r_ptr and phi_ptr.
*/
void LPFormantTracker::trackPhi(dtype *r_ptr, dtype *phi_ptr)
{
	//dtype cumMat[maxFmtTrackJump][maxNTracks];// cumulative cost matrix
	//dtype costMat[maxFmtTrackJump][maxNTracks];// local cost Mat // just for debugging
	const dtype fmts_min[maxNTracks] = {0, 350, 1200, 2000, 3000}; // defines minimal value for each formant
	const dtype fmts_max[maxNTracks] = {1500, 3500, 4500, 5000, 7000};// defines maximal value for each formant
	dtype fn[maxNTracks] = {500, 1500, 2500, 3500, 4500};// neutral formant values (start formants : here vowel [a])
	static dtype last_f[maxNTracks]={500,1500,2500,3500,4500};// last moving average estimated formants 
	int f_list[maxNTracks] = {0, 0, 0, 0, 0};
	const int tri[maxNTracks] = {0, 1, 2, 3, 4};
	dtype this_fmt, this_bw, this_cost, this_cum_cost, low_cost, min_cost, inf_cost=10000000;
	bool in_range = false;
	int k = 0, i0 = 0, j0;
	int bound, indx = 0, new_indx;

	k = 0;
	i0 = 0;
	j0 = 0;

	fn[0] = fn1;
	fn[1] = fn2;

	// loop builds the cumulative cost Matrix cumMat
	// each column represents the cumulative cost for each node which will be the cost entry value for the next column
	for (k = 0; k < nTracks; k++)
	{
		low_cost=inf_cost;
			for(i0=k;i0<(nCands-nTracks+k+2);i0++)
			{
				//cumMat[i0-k][k]=inf_cost;
				this_cum_cost=inf_cost;
				this_fmt=phi_ptr[i0]*sr/(2*M_PI);
				this_bw=-log(r_ptr[i0])*sr/M_PI;
				if((this_fmt>fmts_min[k]) && (this_fmt<fmts_max[k]))// check if actual formant is in range
				{
					in_range=true;
					this_cost=aFact*this_bw+bFact*fabs(this_fmt-fn[k])+gFact*fabs(this_fmt-last_f[k]);//calc local cost
					costMat[(i0 - k) + k * maxFmtTrackJump]=this_cost;
					if (k==0)// build first column: cumulative cost = local cost
						this_cum_cost=this_cost;
					else// build all other columns: cumulative cost(this column) = cumulative cost(previous column)+local cost
						this_cum_cost = cumMat[(i0 - k) + (k - 1) * maxFmtTrackJump] + this_cost;

					if (this_cum_cost<low_cost)
						low_cost=this_cum_cost;// low_cost is the lowest cumulative cost of all elements in this column until element [i0] (included)
											  // therefore, for each column :i0=0 low_cost=cumulative cost
											  //							:i0=n low_cost=min(cumulative_cost(from 0 to [n]))                           
				
				}
				if (k<nTracks-1)// for all columns except last
					cumMat[(i0 - k) + k * maxFmtTrackJump] = low_cost;//ATTENTION: represents the minimal cost that serves as entry cost for the next node (same row element, but next column)
				else// last column  
					cumMat[(i0 - k) + k * maxFmtTrackJump] = this_cum_cost;//shows the overall accumulated cost... from here will start the viterbi traceback
			}
	}
	
	bound=nCands-nTracks+2;// VERY IMPORTANT!!! because values of cumMat beyond this point are not referenced !!
	indx=0;
	// viterbi traceback updates index vector f_list
	// ATTENTION!!! f_list is not the definitive index vector.. has to be diagonalised
	for(k=0;k<nTracks;k++)
	{	
		min_cost=inf_cost;
		for(i0=0;i0<bound;i0++)
		{
			if(cumMat[i0 + (nTracks - k - 1) * maxFmtTrackJump] < min_cost)
			{
				min_cost=cumMat[i0 + (nTracks - k - 1) * maxFmtTrackJump];
				indx=i0;
			}
		}
		if(indx==0)
			break;
		else
		{
			bound=indx+1;
			f_list[nTracks-k-1]=indx;
		}
	}

	// update r, phi and last_f
	for(k=0;k<nTracks;k++)
	{
		new_indx = f_list[k]+k;// rediagonalize index vector
		r_ptr[k] = r_ptr[new_indx];
		phi_ptr[k] = phi_ptr[new_indx];
		last_f[k] = (1 - trackFF) * phi_ptr[k] * sr / (2 * M_PI) + trackFF * last_f[k];
	}
}

/* Computational efficient weithed moving average 
	Input arguments: 
		phi_ptr: unaveraged phi
		bw_ptr: unaveraged bw
		wmaPhi_ptr: moving-averaged phi

	Return value:
		
*/
int LPFormantTracker::mwa(dtype * phi_ptr, dtype * bw_ptr, dtype * wmaPhi_ptr)
{
	int circ_indx_sub = 0;

	circ_indx_sub = (mwaCtr++ - avgLen) % maxAvgLen; // points to the data to withdraw from sum 
	if (circ_indx_sub < 0)
		circ_indx_sub += maxAvgLen;

	sumWei += weiVec[mwaCircCtr] - weiVec[circ_indx_sub];  // update weighting sum

	for (int i0 = 0; i0 < nTracks; ++i0) {
		weiMatPhi[i0 + nTracks * mwaCircCtr] = weiVec[mwaCircCtr] * phi_ptr[i0];
		weiMatBw[i0 + nTracks * mwaCircCtr] = weiVec[mwaCircCtr] * bw_ptr[i0];
			
		sumWeiPhi[i0] += weiMatPhi[i0 + nTracks * mwaCircCtr] - weiMatPhi[i0 + nTracks * circ_indx_sub];

		sumWeiBw[i0] += weiMatBw[i0 + nTracks * mwaCircCtr] - weiMatBw[i0 + nTracks * circ_indx_sub];

		if (sumWei > 0.0000001) {
			wmaPhi_ptr[i0] = sumWeiPhi[i0] / sumWei;
		}
		else {
			mwaCircCtr = mwaCtr % maxAvgLen;
			return 1;
		}
	}

	mwaCircCtr = mwaCtr % maxAvgLen;
	return 0;
}

/* Public interface: process incoming frame 
	Input arguments: 
		xx: input signal frame (assume length == bufferSize)
		st_rms: short-time RMS intensity, 
		        used by the moving weighted average algorithm if bMWA == true

		TODO: radius, phi, bandwidth

		fmts: (output) formant frequencies
*/
void LPFormantTracker::procFrame(dtype * xx, dtype st_rms, 
							     dtype * radius, dtype * phi, dtype * bandwidth, 
								 dtype * fmts) { //Marked
	for (int i0 = 0; i0 < nLPC + 1; i0++) {		//SC Initialize the order LPC coefficients
		lpcAi[i0] = 0.0; // TODO: Incorporate into fmtTracker
	}

	getAi(xx, lpcAi);

	int quit_hqr = hqr_roots(lpcAi, realRoots, imagRoots);

	/*getRPhiBw(realRoots, imagRoots, radius, phi, bandwidth);*/
	getRPhiBw(realRoots, imagRoots, radius, phi_us, bandwidth);

	/*trackPhi(radius, phi);*/
	trackPhi(radius, phi_us);

	/* Moving window average (MWA) (Optional, but highly recommended) */
	weiVec[mwaCircCtr] = st_rms;
	if (bMWA) {
		mwa(phi_us, bandwidth, phi);
	}
	else {
		/* Copy from phi_us to phi */
		for (int i = 0; i < nTracks; ++i) {
			phi[i] = phi_us[i];
		}
	}

	for (int i = 0; i < nTracks; ++i) //Marked
		fmts[i] = phi[i] * sr /(2 * M_PI);
	
	/* DEBUG */
	//dtype * n_radius = new dtype[maxNPoles];
	//dtype * n_phi = new dtype[maxNPoles];

	//for (int i = 0; i < maxNPoles; ++i) {
	//	n_radius[i] = radius[i];
	//	n_phi[i] = phi[i];
	//}
	//
	//trackPhi(n_radius, n_phi); /* DEBUG */
	//
	///* DEBUG */
	//delete [] n_radius;
	//delete [] n_phi;
	
}


/* Setters and getters (inline) */
void LPFormantTracker::setBCepsLift(const bool t_bCepsLift) {
	bCepsLift = t_bCepsLift;
}

void LPFormantTracker::setNLPC(const int t_nLPC) {
	nLPC = t_nLPC;
}

void LPFormantTracker::setBufferSize(const int t_bufferSize) {
	bufferSize = t_bufferSize;
}

void LPFormantTracker::setNFFT(const int t_nFFT) {
	nFFT = t_nFFT;

	/* TODO: re-initiliaze windows */
}

void LPFormantTracker::setCepsWinWidth(const int t_cepsWinWidth) {
	cepsWinWidth = t_cepsWinWidth;

	/* TODO: re-initialize FFT data */
}

void LPFormantTracker::genHanningWindow() {
	for(int i0 = 0; i0 < bufferSize; i0++){ /* Total length: p.anaLen */
		winFunc[i0] = 0.5*cos(dtype(2 * M_PI * (i0 + 1)) / dtype(bufferSize + 1));  //Marked
		winFunc[i0] = 0.5 - winFunc[i0]; //Marked
		winFunc[bufferSize - i0 - 1] = winFunc[i0];
	}
}

const int LPFormantTracker::getNLPC() const {
	return nLPC;
}