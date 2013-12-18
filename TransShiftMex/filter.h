#ifndef FILTER_H
#define FILTER_H

template <class D_TYPE>
class IIR_Filter {
public:
	/* Member variables */
	int filtLen;  /* Filter length: Number of coefficients */
	int buffLen;  /* Buffer length */ 
	D_TYPE * a;   /* Denominator coefficients */
	D_TYPE * b;   /* Numerator coefficients */	
	D_TYPE * buff; /* Buffer */
	D_TYPE * delay; /* Filter delay */

	/* Error classes */
	class filterLengthError {};
	class bufferLengthError {};

	/* Member functions */
	/* Constructor */
	IIR_Filter(int t_filtLen, int t_buffLen) : filtLen(t_filtLen), buffLen(t_buffLen) {
		if (t_filtLen > 1) {
			a = new D_TYPE[filtLen];
			b = new D_TYPE[filtLen];
			delay = new D_TYPE[filtLen - 1];
		}
		else {
			throw filterLengthError();
		}

		if (t_buffLen > 0) {
			buff = new D_TYPE[buffLen];
		}
		else {
			throw bufferLengthError();
		}
	}

	/* Set filter coefficients */
	void setCoeff(int aLen, const D_TYPE ta[], int bLen, const D_TYPE tb[]) 
		throw(filterLengthError) {
		if ( (aLen > filtLen) || (bLen > filtLen) )
			throw filterLengthError();

		for (int i = 0; i < aLen; ++i)
			a[i] = ta[i];
		for (int i = 0; i < bLen; ++i)
			b[i] = tb[i];
	}	

	/* Reset buffer */
	void resetBuffer() {
		for (int i = 0; i < buffLen; ++i) {
			buff[i] = 0.0;
		}
		for (int i = 0; i < filtLen - 1; ++i) {
			delay[i] = 0.0;
		}
	}

	/* Filter signal frame */
	void filter (D_TYPE *x, const int nr, D_TYPE g=1.0) {// IIR Direct II transposed form
		// Input: 
		//		x : input frame
		//		nr: length of output frame
		//		g : additional gain
		// Output: 
		//		buff (public member)
		int m, k;

		for(m = 0; m < nr; m++) {
			buff[m] = g * b[0] * x[m] + delay[0];
			for(k = 0; k < filtLen - 2; k++) {// start delay recursion
				delay[k] = g * b[k + 1] * x[m] + delay[k + 1] - a[k + 1] * buff[m];
			}	
			delay[filtLen - 2] = g * b[filtLen - 1] * x[m] - a[filtLen - 1] * buff[m]; 
		}
	}

	/* Destructor */
	~IIR_Filter() {
		if (a)
			delete [] a;
		if (b) 
			delete [] b;
		if (delay) 
			delete [] delay;
		if (buff) 
			delete [] buff;
	}

};

#endif