#ifndef DSPF_H
#define DSPF_H

typedef double dtype;

void	DSPF_dp_blk_move(const dtype * x, dtype * r, const int nx);
dtype	DSPF_dp_vecsum_sq(const dtype *x,int n);        
void	DSPF_dp_vecmul(const dtype * x, const dtype * y, dtype * r, const int &n);
void	DSPF_dp_autocor(dtype * r, dtype * x, const int & nx, const int & nr);

void gen_w_r2(double* w, int n);
void bit_rev(double* x, int n);

void DSPF_dp_cfftr2(const int n, double * x, double * w, const int n_min);
void DSPF_dp_icfftr2(const int n, double * x, double * w, const int n_min);

#endif 