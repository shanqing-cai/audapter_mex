/* DSPF.cpp 

DSP functions

*/

#include <cmath>
#include <memory>

#include "DSPF.h"

/* copy a block of data from one address x to r */
void DSPF_dp_blk_move(const double * x, double * r, const int nx)	
{
    memcpy(r, x, nx * sizeof(double));
    /*for (int i = 0; i < nx; i++) {
        r[i] = x[i];
    }*/
}
double DSPF_dp_vecsum_sq(const double *x,int n)    //SC get norm (power) of a vector                       
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}              
void DSPF_dp_fir_gen(const double *  x, const double *  h, double *  r, int nh, int nr)	//SC convolution or correlation
{
int i, j;
double sum;


 for(i=0; i < nr; i++)
  {
    sum = 0;/*Clear the accumulator*/
    for(j=0; j < nh; j++)
    {
     sum += x[i+j] * h[j];
    }
   r[i] = sum; /*Storing the Sum*/
  }
}

double DSPF_dp_maxval(const double* x, const int &nx)		//SC maximum of a vector
{
// x: input array
// nx: size of input array

int i;
double max;

*((int *)&max) = 0x00000000;
*((int *)&max+1) = 0xfff00000;

for (i = 0; i < nx; i++)
if (x[i] > max)
{
max = x[i];
}

return max;
}

void DSPF_dp_vecmul(const double * x, const double * y, double * r, const int &n)		//SC vector multiplication
//SC Input arguments
//SC	x: vector 1 pointer
//SC	y: vecotr 2 poniter
//SC	r: output pointer
//SC	n: length of vectors 1 and 2
{
    int i;

    for(i = 0; i < n; i++)
        r[i] = x[i] * y[i];
}

void DSPF_dp_autocor(double * r, double * x, const int & nx, const int & nr)			//SC autocorrelation
//SC Input arguments
//SC	r: output pointer
//SC	x: input signal pointer
//SC	nx: length of input signal (x)
//SC	nr: number of delays
{
	int i, k;
	double sum;

	for (i = 0; i < nr; i++) {
	    sum = 0.0;
		for (k = nr; k < nx + nr; k++)
			sum += x[k] * x[k - i];

		r[i] = sum;
	}
}

void gen_w_r2(double* w, int n)		//SC An FFT subroutine
{
	int i, j=1;
	double pi = 4.0*atan(1.0);
	double e = pi*2.0/n;
	for(j=1; j < n; j <<= 1)
	{
		for(i=0; i < ( n>>1 ); i += j)
		{
			*w++ = cos(i*e);
			*w++ = -sin(i*e);
		}
	}
}

void bit_rev(double* x, int n)	//SC Bit reversal: an FFT subroutine
{
	int i, j, k;
	double rtemp, itemp;

	j = 0;
	for(i=1; i < (n-1); i++)
	{
		k = n >> 1;
		while(k <= j)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
		if(i < j)
		{
			rtemp = x[j*2];
			x[j*2] = x[i*2];
			x[i*2] = rtemp;
			itemp = x[j*2+1];
			x[j*2+1] = x[i*2+1];
			x[i*2+1] = itemp;
		}
	}
}


void DSPF_dp_cfftr2(const int n, double * x, double * w, const int n_min)	//SC Fast Fourier transform on complex numbers
{
	int n2, ie, ia, i, j, k, m;
	double rtemp, itemp, c, s;
	n2 = n;

	ie = 1;

	for(k = n; k > n_min; k >>= 1)
	{
		n2 >>= 1;
		ia = 0;
		for(j=0; j < ie; j++)
		{
			for(i=0; i < n2; i++)
			{
				c = w[2*i];
				s = w[2*i+1];
				m = ia + n2;
				rtemp = x[2*ia] - x[2*m];
				x[2*ia] = x[2*ia] + x[2*m];
				itemp = x[2*ia+1] - x[2*m+1];
				x[2*ia+1] = x[2*ia+1] + x[2*m+1];
				x[2*m] = c*rtemp - s*itemp;
				x[2*m+1] = c*itemp + s*rtemp;
				ia++;
			}
			ia += n2;
		}	
		ie <<= 1;
		w = w + k;
	}
}

void DSPF_dp_icfftr2(const int n, double * x, double * w, const int n_min)	//SC Inverse Fast Fourier transform on complex numbers
{
	int n2, ie, ia, i, j, k, m;
	double rtemp, itemp, c, s;
	n2 = n;
	ie = 1;
	for(k = n; k > n_min; k >>= 1)
	{
		n2 >>= 1;
		ia = 0;
		for(j=0; j < ie; j++)
		{
			for(i=0; i < n2; i++)
			{
				c = w[2*i];
				s = w[2*i+1];
				m = ia + n2;
				rtemp = x[2*ia] - x[2*m];
				x[2*ia] = x[2*ia] + x[2*m];
				itemp = x[2*ia+1] - x[2*m+1];
				x[2*ia+1] = x[2*ia+1] + x[2*m+1];
				x[2*m] = c*rtemp + s*itemp;
				x[2*m+1] = c*itemp - s*rtemp;
				ia++;
			}
			ia += n2;
		}
		ie <<= 1;
		w = w + k;
	}
}