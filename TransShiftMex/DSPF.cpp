/* DSPF.cpp 

DSP functions

*/

#include <cmath>

void DSPF_dp_blk_move(const double * x, double * r, const int nx)	//SC copy a block of data from one address x to r
{
	//memmove((void*)r,(void*)x,nx*sizeof(double));
	//memcpy((void*)r,(void*)x,nx*sizeof(double));
	int i;
	for (i = 0 ; i < nx; i++)
		r[i] = x[i];
}
double DSPF_dp_vecsum_sq(const double *x,int n)    //SC get norm (power) of a vector                       
{
	int i;
    double sum=0;
    for(i = 0;  i < n; i++ )
		sum += x[i]*x[i];
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

