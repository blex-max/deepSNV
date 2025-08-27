#ifdef __cplusplus
extern "C" {
#endif

void bam2R(char **bamfile, char **ref,
          int *beg, int *end, int *counts,
          int *q, int *mq, int *s,
          int *head_clip, int *maxdepth, int *verbose,
          int *mask, int *keepflag, int *maxmismatches);

void dbetabinom(double* p, int* lp, int *x, int* lx, int *n, int* ln, double *mu, int* lmu, double *disp, int* ldisp, int *logp);

void pbetabinom(double* p, int*lp, int *x, int* lx, int *n, int* ln, double *mu, int* lmu, double *disp, int* ldisp, int *logp);

#ifdef __cplusplus
}
#endif
