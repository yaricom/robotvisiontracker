//
//  onlineSVM.h
//  RobotVisionTracker
//
//  Created by Iaroslav Omelianenko on 4/3/15.
//  Copyright (c) 2015 Nologin. All rights reserved.
//

#ifndef RobotVisionTracker_onlineSVM_h
#define RobotVisionTracker_onlineSVM_h

//
//-----------------------------Start TRON------------------------------------------------------
typedef struct { float r, i; } fcomplex;
typedef struct { double r, i; } dcomplex;
typedef int blasbool;

#ifdef F2C_COMPAT

void cdotc_(fcomplex *dotval, int *n, fcomplex *cx, int *incx,
            fcomplex *cy, int *incy);

void cdotu_(fcomplex *dotval, int *n, fcomplex *cx, int *incx,
            fcomplex *cy, int *incy);

double sasum_(int *n, float *sx, int *incx);

double scasum_(int *n, fcomplex *cx, int *incx);

double scnrm2_(int *n, fcomplex *x, int *incx);

double sdot_(int *n, float *sx, int *incx, float *sy, int *incy);

double snrm2_(int *n, float *x, int *incx);

void zdotc_(dcomplex *dotval, int *n, dcomplex *cx, int *incx,
            dcomplex *cy, int *incy);

void zdotu_(dcomplex *dotval, int *n, dcomplex *cx, int *incx,
            dcomplex *cy, int *incy);

#else

fcomplex cdotc_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

fcomplex cdotu_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

float sasum_(int *n, float *sx, int *incx);

float scasum_(int *n, fcomplex *cx, int *incx);

float scnrm2_(int *n, fcomplex *x, int *incx);

float sdot_(int *n, float *sx, int *incx, float *sy, int *incy);

float snrm2_(int *n, float *x, int *incx);

dcomplex zdotc_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

dcomplex zdotu_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

#endif

/* Remaining functions listed in alphabetical order */

int caxpy_(int *n, fcomplex *ca, fcomplex *cx, int *incx, fcomplex *cy,
           int *incy);

int ccopy_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

int cgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           fcomplex *alpha, fcomplex *a, int *lda, fcomplex *x, int *incx,
           fcomplex *beta, fcomplex *y, int *incy);

int cgemm_(char *transa, char *transb, int *m, int *n, int *k,
           fcomplex *alpha, fcomplex *a, int *lda, fcomplex *b, int *ldb,
           fcomplex *beta, fcomplex *c, int *ldc);

int cgemv_(char *trans, int *m, int *n, fcomplex *alpha, fcomplex *a,
           int *lda, fcomplex *x, int *incx, fcomplex *beta, fcomplex *y,
           int *incy);

int cgerc_(int *m, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *a, int *lda);

int cgeru_(int *m, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *a, int *lda);

int chbmv_(char *uplo, int *n, int *k, fcomplex *alpha, fcomplex *a,
           int *lda, fcomplex *x, int *incx, fcomplex *beta, fcomplex *y,
           int *incy);

int chemm_(char *side, char *uplo, int *m, int *n, fcomplex *alpha,
           fcomplex *a, int *lda, fcomplex *b, int *ldb, fcomplex *beta,
           fcomplex *c, int *ldc);

int chemv_(char *uplo, int *n, fcomplex *alpha, fcomplex *a, int *lda,
           fcomplex *x, int *incx, fcomplex *beta, fcomplex *y, int *incy);

int cher_(char *uplo, int *n, float *alpha, fcomplex *x, int *incx,
          fcomplex *a, int *lda);

int cher2_(char *uplo, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *a, int *lda);

int cher2k_(char *uplo, char *trans, int *n, int *k, fcomplex *alpha,
            fcomplex *a, int *lda, fcomplex *b, int *ldb, float *beta,
            fcomplex *c, int *ldc);

int cherk_(char *uplo, char *trans, int *n, int *k, float *alpha,
           fcomplex *a, int *lda, float *beta, fcomplex *c, int *ldc);

int chpmv_(char *uplo, int *n, fcomplex *alpha, fcomplex *ap, fcomplex *x,
           int *incx, fcomplex *beta, fcomplex *y, int *incy);

int chpr_(char *uplo, int *n, float *alpha, fcomplex *x, int *incx,
          fcomplex *ap);

int chpr2_(char *uplo, int *n, fcomplex *alpha, fcomplex *x, int *incx,
           fcomplex *y, int *incy, fcomplex *ap);

int crotg_(fcomplex *ca, fcomplex *cb, float *c, fcomplex *s);

int cscal_(int *n, fcomplex *ca, fcomplex *cx, int *incx);

int csscal_(int *n, float *sa, fcomplex *cx, int *incx);

int cswap_(int *n, fcomplex *cx, int *incx, fcomplex *cy, int *incy);

int csymm_(char *side, char *uplo, int *m, int *n, fcomplex *alpha,
           fcomplex *a, int *lda, fcomplex *b, int *ldb, fcomplex *beta,
           fcomplex *c, int *ldc);

int csyr2k_(char *uplo, char *trans, int *n, int *k, fcomplex *alpha,
            fcomplex *a, int *lda, fcomplex *b, int *ldb, fcomplex *beta,
            fcomplex *c, int *ldc);

int csyrk_(char *uplo, char *trans, int *n, int *k, fcomplex *alpha,
           fcomplex *a, int *lda, fcomplex *beta, fcomplex *c, int *ldc);

int ctbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           fcomplex *a, int *lda, fcomplex *x, int *incx);

int ctbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           fcomplex *a, int *lda, fcomplex *x, int *incx);

int ctpmv_(char *uplo, char *trans, char *diag, int *n, fcomplex *ap,
           fcomplex *x, int *incx);

int ctpsv_(char *uplo, char *trans, char *diag, int *n, fcomplex *ap,
           fcomplex *x, int *incx);

int ctrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, fcomplex *alpha, fcomplex *a, int *lda, fcomplex *b,
           int *ldb);

int ctrmv_(char *uplo, char *trans, char *diag, int *n, fcomplex *a,
           int *lda, fcomplex *x, int *incx);

int ctrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, fcomplex *alpha, fcomplex *a, int *lda, fcomplex *b,
           int *ldb);

int ctrsv_(char *uplo, char *trans, char *diag, int *n, fcomplex *a,
           int *lda, fcomplex *x, int *incx);

int daxpy_(int *n, double *sa, double *sx, int *incx, double *sy,
           int *incy);

int dcopy_(int *n, double *sx, int *incx, double *sy, int *incy);

int dgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           double *alpha, double *a, int *lda, double *x, int *incx,
           double *beta, double *y, int *incy);

int dgemm_(char *transa, char *transb, int *m, int *n, int *k,
           double *alpha, double *a, int *lda, double *b, int *ldb,
           double *beta, double *c, int *ldc);

int dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
           int *lda, double *x, int *incx, double *beta, double *y,
           int *incy);

int dger_(int *m, int *n, double *alpha, double *x, int *incx,
          double *y, int *incy, double *a, int *lda);

int drot_(int *n, double *sx, int *incx, double *sy, int *incy,
          double *c, double *s);

int drotg_(double *sa, double *sb, double *c, double *s);

int dsbmv_(char *uplo, int *n, int *k, double *alpha, double *a,
           int *lda, double *x, int *incx, double *beta, double *y,
           int *incy);

int dscal_(int *n, double *sa, double *sx, int *incx);

int dspmv_(char *uplo, int *n, double *alpha, double *ap, double *x,
           int *incx, double *beta, double *y, int *incy);

int dspr_(char *uplo, int *n, double *alpha, double *x, int *incx,
          double *ap);

int dspr2_(char *uplo, int *n, double *alpha, double *x, int *incx,
           double *y, int *incy, double *ap);

int dswap_(int *n, double *sx, int *incx, double *sy, int *incy);

int dsymm_(char *side, char *uplo, int *m, int *n, double *alpha,
           double *a, int *lda, double *b, int *ldb, double *beta,
           double *c, int *ldc);

int dsymv_(char *uplo, int *n, double *alpha, double *a, int *lda,
           double *x, int *incx, double *beta, double *y, int *incy);

int dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx,
          double *a, int *lda);

int dsyr2_(char *uplo, int *n, double *alpha, double *x, int *incx,
           double *y, int *incy, double *a, int *lda);

int dsyr2k_(char *uplo, char *trans, int *n, int *k, double *alpha,
            double *a, int *lda, double *b, int *ldb, double *beta,
            double *c, int *ldc);

int dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha,
           double *a, int *lda, double *beta, double *c, int *ldc);

int dtbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           double *a, int *lda, double *x, int *incx);

int dtbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           double *a, int *lda, double *x, int *incx);

int dtpmv_(char *uplo, char *trans, char *diag, int *n, double *ap,
           double *x, int *incx);

int dtpsv_(char *uplo, char *trans, char *diag, int *n, double *ap,
           double *x, int *incx);

int dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, double *alpha, double *a, int *lda, double *b,
           int *ldb);

int dtrmv_(char *uplo, char *trans, char *diag, int *n, double *a,
           int *lda, double *x, int *incx);

int dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, double *alpha, double *a, int *lda, double *b,
           int *ldb);

int dtrsv_(char *uplo, char *trans, char *diag, int *n, double *a,
           int *lda, double *x, int *incx);


int saxpy_(int *n, float *sa, float *sx, int *incx, float *sy, int *incy);

int scopy_(int *n, float *sx, int *incx, float *sy, int *incy);

int sgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           float *alpha, float *a, int *lda, float *x, int *incx,
           float *beta, float *y, int *incy);

int sgemm_(char *transa, char *transb, int *m, int *n, int *k,
           float *alpha, float *a, int *lda, float *b, int *ldb,
           float *beta, float *c, int *ldc);

int sgemv_(char *trans, int *m, int *n, float *alpha, float *a,
           int *lda, float *x, int *incx, float *beta, float *y,
           int *incy);

int sger_(int *m, int *n, float *alpha, float *x, int *incx,
          float *y, int *incy, float *a, int *lda);

int srot_(int *n, float *sx, int *incx, float *sy, int *incy,
          float *c, float *s);

int srotg_(float *sa, float *sb, float *c, float *s);

int ssbmv_(char *uplo, int *n, int *k, float *alpha, float *a,
           int *lda, float *x, int *incx, float *beta, float *y,
           int *incy);

int sscal_(int *n, float *sa, float *sx, int *incx);

int sspmv_(char *uplo, int *n, float *alpha, float *ap, float *x,
           int *incx, float *beta, float *y, int *incy);

int sspr_(char *uplo, int *n, float *alpha, float *x, int *incx,
          float *ap);

int sspr2_(char *uplo, int *n, float *alpha, float *x, int *incx,
           float *y, int *incy, float *ap);

int sswap_(int *n, float *sx, int *incx, float *sy, int *incy);

int ssymm_(char *side, char *uplo, int *m, int *n, float *alpha,
           float *a, int *lda, float *b, int *ldb, float *beta,
           float *c, int *ldc);

int ssymv_(char *uplo, int *n, float *alpha, float *a, int *lda,
           float *x, int *incx, float *beta, float *y, int *incy);

int ssyr_(char *uplo, int *n, float *alpha, float *x, int *incx,
          float *a, int *lda);

int ssyr2_(char *uplo, int *n, float *alpha, float *x, int *incx,
           float *y, int *incy, float *a, int *lda);

int ssyr2k_(char *uplo, char *trans, int *n, int *k, float *alpha,
            float *a, int *lda, float *b, int *ldb, float *beta,
            float *c, int *ldc);

int ssyrk_(char *uplo, char *trans, int *n, int *k, float *alpha,
           float *a, int *lda, float *beta, float *c, int *ldc);

int stbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           float *a, int *lda, float *x, int *incx);

int stbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           float *a, int *lda, float *x, int *incx);

int stpmv_(char *uplo, char *trans, char *diag, int *n, float *ap,
           float *x, int *incx);

int stpsv_(char *uplo, char *trans, char *diag, int *n, float *ap,
           float *x, int *incx);

int strmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, float *alpha, float *a, int *lda, float *b,
           int *ldb);

int strmv_(char *uplo, char *trans, char *diag, int *n, float *a,
           int *lda, float *x, int *incx);

int strsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, float *alpha, float *a, int *lda, float *b,
           int *ldb);

int strsv_(char *uplo, char *trans, char *diag, int *n, float *a,
           int *lda, float *x, int *incx);

int zaxpy_(int *n, dcomplex *ca, dcomplex *cx, int *incx, dcomplex *cy,
           int *incy);

int zcopy_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

int zdscal_(int *n, double *sa, dcomplex *cx, int *incx);

int zgbmv_(char *trans, int *m, int *n, int *kl, int *ku,
           dcomplex *alpha, dcomplex *a, int *lda, dcomplex *x, int *incx,
           dcomplex *beta, dcomplex *y, int *incy);

int zgemm_(char *transa, char *transb, int *m, int *n, int *k,
           dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b, int *ldb,
           dcomplex *beta, dcomplex *c, int *ldc);

int zgemv_(char *trans, int *m, int *n, dcomplex *alpha, dcomplex *a,
           int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y,
           int *incy);

int zgerc_(int *m, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *a, int *lda);

int zgeru_(int *m, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *a, int *lda);

int zhbmv_(char *uplo, int *n, int *k, dcomplex *alpha, dcomplex *a,
           int *lda, dcomplex *x, int *incx, dcomplex *beta, dcomplex *y,
           int *incy);

int zhemm_(char *side, char *uplo, int *m, int *n, dcomplex *alpha,
           dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta,
           dcomplex *c, int *ldc);

int zhemv_(char *uplo, int *n, dcomplex *alpha, dcomplex *a, int *lda,
           dcomplex *x, int *incx, dcomplex *beta, dcomplex *y, int *incy);

int zher_(char *uplo, int *n, double *alpha, dcomplex *x, int *incx,
          dcomplex *a, int *lda);

int zher2_(char *uplo, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *a, int *lda);

int zher2k_(char *uplo, char *trans, int *n, int *k, dcomplex *alpha,
            dcomplex *a, int *lda, dcomplex *b, int *ldb, double *beta,
            dcomplex *c, int *ldc);

int zherk_(char *uplo, char *trans, int *n, int *k, double *alpha,
           dcomplex *a, int *lda, double *beta, dcomplex *c, int *ldc);

int zhpmv_(char *uplo, int *n, dcomplex *alpha, dcomplex *ap, dcomplex *x,
           int *incx, dcomplex *beta, dcomplex *y, int *incy);

int zhpr_(char *uplo, int *n, double *alpha, dcomplex *x, int *incx,
          dcomplex *ap);

int zhpr2_(char *uplo, int *n, dcomplex *alpha, dcomplex *x, int *incx,
           dcomplex *y, int *incy, dcomplex *ap);

int zrotg_(dcomplex *ca, dcomplex *cb, double *c, dcomplex *s);

int zscal_(int *n, dcomplex *ca, dcomplex *cx, int *incx);

int zswap_(int *n, dcomplex *cx, int *incx, dcomplex *cy, int *incy);

int zsymm_(char *side, char *uplo, int *m, int *n, dcomplex *alpha,
           dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta,
           dcomplex *c, int *ldc);

int zsyr2k_(char *uplo, char *trans, int *n, int *k, dcomplex *alpha,
            dcomplex *a, int *lda, dcomplex *b, int *ldb, dcomplex *beta,
            dcomplex *c, int *ldc);

int zsyrk_(char *uplo, char *trans, int *n, int *k, dcomplex *alpha,
           dcomplex *a, int *lda, dcomplex *beta, dcomplex *c, int *ldc);

int ztbmv_(char *uplo, char *trans, char *diag, int *n, int *k,
           dcomplex *a, int *lda, dcomplex *x, int *incx);

int ztbsv_(char *uplo, char *trans, char *diag, int *n, int *k,
           dcomplex *a, int *lda, dcomplex *x, int *incx);

int ztpmv_(char *uplo, char *trans, char *diag, int *n, dcomplex *ap,
           dcomplex *x, int *incx);

int ztpsv_(char *uplo, char *trans, char *diag, int *n, dcomplex *ap,
           dcomplex *x, int *incx);

int ztrmm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b,
           int *ldb);

int ztrmv_(char *uplo, char *trans, char *diag, int *n, dcomplex *a,
           int *lda, dcomplex *x, int *incx);

int ztrsm_(char *side, char *uplo, char *transa, char *diag, int *m,
           int *n, dcomplex *alpha, dcomplex *a, int *lda, dcomplex *b,
           int *ldb);

int ztrsv_(char *uplo, char *trans, char *diag, int *n, dcomplex *a,
           int *lda, dcomplex *x, int *incx);



int dscal_(int *n, double *sa, double *sx, int *incx) {
    long int i, m, nincx, nn, iincx;
    double ssa;
    
    /* scales a vector by a constant.
     uses unrolled loops for increment equal to 1.
     jack dongarra, linpack, 3/11/78.
     modified 3/93 to return if incx .le. 0.
     modified 12/3/93, array(1) declarations changed to array(*) */
    
    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    ssa = *sa;
    
    if (nn > 0 && iincx > 0)
    {
        if (iincx == 1) /* code for increment equal to 1 */
        {
            m = nn-4;
            for (i = 0; i < m; i += 5)
            {
                sx[i] = ssa * sx[i];
                sx[i+1] = ssa * sx[i+1];
                sx[i+2] = ssa * sx[i+2];
                sx[i+3] = ssa * sx[i+3];
                sx[i+4] = ssa * sx[i+4];
            }
            for ( ; i < nn; ++i) /* clean-up loop */
                sx[i] = ssa * sx[i];
        }
        else /* code for increment not equal to 1 */
        {
            nincx = nn * iincx;
            for (i = 0; i < nincx; i += iincx)
                sx[i] = ssa * sx[i];
        }
    }
    
    return 0;
}

double dnrm2_(int *n, double *x, int *incx) {
    long int ix, nn, iincx;
    double norm, scale, absxi, ssq, temp;
    
    /*  DNRM2 returns the euclidean norm of a vector via the function
     name, so that
     
     DNRM2 := sqrt( x'*x )
     
     -- This version written on 25-October-1982.
     Modified on 14-October-1993 to inline the call to SLASSQ.
     Sven Hammarling, Nag Ltd.   */
    
    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    
    if( nn > 0 && iincx > 0 )
    {
        if (nn == 1)
        {
            norm = fabs(x[0]);
        }
        else
        {
            scale = 0.0;
            ssq = 1.0;
            
            /* The following loop is equivalent to this call to the LAPACK
             auxiliary routine:   CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */
            
            for (ix=(nn-1)*iincx; ix>=0; ix-=iincx) {
                if (x[ix] != 0.0 && !isnan(x[ix]))
                {
                    absxi = fabs(x[ix]);
                    if (scale < absxi)
                    {
                        temp = scale / absxi;
                        ssq = ssq * (temp * temp) + 1.0;
                        scale = absxi;
                    }
                    else if (scale != 0.0)
                    {
                        temp = absxi / scale;
                        ssq += temp * temp;
                    }
                }
                if (isnan(ssq)) {
                    Printf("Epic fail. ssq is NaN");
                }
            }
            norm = scale * sqrt(ssq);
        }
    }
    else
        norm = 0.0;
    
    return norm;
    
}

double ddot_(int *n, double *sx, int *incx, double *sy, int *incy) {
    long int i, m, nn, iincx, iincy;
    double stemp;
    long int ix, iy;
    
    /* forms the dot product of two vectors.
     uses unrolled loops for increments equal to one.
     jack dongarra, linpack, 3/11/78.
     modified 12/3/93, array(1) declarations changed to array(*) */
    
    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    iincy = *incy;
    
    stemp = 0.0;
    if (nn > 0)
    {
        if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
        {
            m = nn-4;
            for (i = 0; i < m; i += 5)
                stemp += sx[i] * sy[i] + sx[i+1] * sy[i+1] + sx[i+2] * sy[i+2] +
                sx[i+3] * sy[i+3] + sx[i+4] * sy[i+4];
            
            for ( ; i < nn; i++)        /* clean-up loop */
                stemp += sx[i] * sy[i];
        }
        else /* code for unequal increments or equal increments not equal to 1 */
        {
            ix = 0;
            iy = 0;
            if (iincx < 0)
                ix = (1 - nn) * iincx;
            if (iincy < 0)
                iy = (1 - nn) * iincy;
            for (i = 0; i < nn; i++)
            {
                stemp += sx[ix] * sy[iy];
                ix += iincx;
                iy += iincy;
            }
        }
    }
    
    return stemp;
}

int daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy) {
    long int i, m, ix, iy, nn, iincx, iincy;
    register double ssa;
    
    /* constant times a vector plus a vector.
     uses unrolled loop for increments equal to one.
     jack dongarra, linpack, 3/11/78.
     modified 12/3/93, array(1) declarations changed to array(*) */
    
    /* Dereference inputs */
    nn = *n;
    ssa = *sa;
    iincx = *incx;
    iincy = *incy;
    
    if( nn > 0 && ssa != 0.0 )
    {
        if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
        {
            m = nn-3;
            for (i = 0; i < m; i += 4)
            {
                sy[i] += ssa * sx[i];
                sy[i+1] += ssa * sx[i+1];
                sy[i+2] += ssa * sx[i+2];
                sy[i+3] += ssa * sx[i+3];
            }
            for ( ; i < nn; ++i) /* clean-up loop */
                sy[i] += ssa * sx[i];
        }
        else /* code for unequal increments or equal increments not equal to 1 */
        {
            ix = iincx >= 0 ? 0 : (1 - nn) * iincx;
            iy = iincy >= 0 ? 0 : (1 - nn) * iincy;
            for (i = 0; i < nn; i++)
            {
                sy[iy] += ssa * sx[ix];
                ix += iincx;
                iy += iincy;
            }
        }
    }
    
    return 0;
}

static void default_print(const char *buf)
{
    fputs(buf,stdout);
    fflush(stdout);
}

class tronFunction
{
public:
    virtual double fun(double *w) = 0 ;
    virtual void grad(double *w, double *g) = 0 ;
    virtual void Hv(double *s, double *Hs) = 0 ;
    
    virtual int get_nr_variable(void) = 0 ;
    virtual ~tronFunction(void){}
};

class TRON
{
public:
    TRON(const tronFunction *fun_obj, double eps = 0.1, int max_iter = 1000);
    
    void tron(double *w);
    void set_print_string(void (*i_print) (const char *buf));
    
private:
    int trcg(double delta, double *g, double *s, double *r);
    double norm_inf(int n, double *x);
    
    double eps;
    int max_iter;
    tronFunction *fun_obj;
    void info(const char *fmt,...);
    void (*tron_print_string)(const char *buf);
};

void TRON::info(const char *fmt,...)
{
    char buf[BUFSIZ];
    va_list ap;
    va_start(ap,fmt);
    vsprintf(buf,fmt,ap);
    va_end(ap);
    (*tron_print_string)(buf);
}

TRON::TRON(const tronFunction *fun_obj, double eps, int max_iter)
{
    this->fun_obj=const_cast<tronFunction *>(fun_obj);
    this->eps=eps;
    this->max_iter=max_iter;
    tron_print_string = default_print;
}

void TRON::tron(double *w)
{
    // Parameters for updating the iterates.
    double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;
    
    // Parameters for updating the trust region size delta.
    double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;
    
    int n = fun_obj->get_nr_variable();
    int i, cg_iter;
    double delta, snorm, one=1.0;
    double alpha, f, fnew, prered, actred, gs;
    int search = 1, iter = 1, inc = 1;
    double *s = new double[n];
    double *r = new double[n];
    double *w_new = new double[n];
    double *g = new double[n];
    
    double *w0 = new double[n];
    for (i=0; i<n; i++)
        w0[i] = 0;
    fun_obj->fun(w0);
    fun_obj->grad(w0, g);
    double gnorm1 = dnrm2_(&n, g, &inc);
    delete [] w0;
    
    f = fun_obj->fun(w);
    fun_obj->grad(w, g);
    delta = dnrm2_(&n, g, &inc);
    double gnorm = delta;
    
    if (gnorm <= eps*gnorm1)
        search = 0;
    
    iter = 1;
    
    while (iter <= max_iter && search)
    {
        cg_iter = trcg(delta, g, s, r);
        
        memcpy(w_new, w, sizeof(double)*n);
        daxpy_(&n, &one, s, &inc, w_new, &inc);
        
        gs = ddot_(&n, g, &inc, s, &inc);
        prered = -0.5*(gs-ddot_(&n, s, &inc, r, &inc));
        fnew = fun_obj->fun(w_new);
        
        // Compute the actual reduction.
        actred = f - fnew;
        
        // On the first iteration, adjust the initial step bound.
        snorm = dnrm2_(&n, s, &inc);
        if (iter == 1)
            delta = min(delta, snorm);
        
        // Compute prediction alpha*snorm of the step.
        if (fnew - f - gs <= 0)
            alpha = sigma3;
        else
            alpha = max(sigma1, -0.5*(gs/(fnew - f - gs)));
        
        // Update the trust region bound according to the ratio of actual to predicted reduction.
        if (actred < eta0*prered)
            delta = min(max(alpha, sigma1)*snorm, sigma2*delta);
        else if (actred < eta1*prered)
            delta = max(sigma1*delta, min(alpha*snorm, sigma2*delta));
        else if (actred < eta2*prered)
            delta = max(sigma1*delta, min(alpha*snorm, sigma3*delta));
        else
            delta = max(delta, min(alpha*snorm, sigma3*delta));
        
        info("iter %2d act %5.3e pre %5.3e delta %5.3e f %5.3e |g| %5.3e CG %3d\n", iter, actred, prered, delta, f, gnorm, cg_iter);
        
        if (isnan(actred)) {
            Printf("Epic failure in TRON::tron. Skipping.\n");
            break;
        }
        
        if (actred > eta0*prered)
        {
            iter++;
            memcpy(w, w_new, sizeof(double)*n);
            f = fnew;
            fun_obj->grad(w, g);
            
            gnorm = dnrm2_(&n, g, &inc);
            if (gnorm <= eps*gnorm1)
                break;
        }
        if (f < -1.0e+32)
        {
            info("WARNING: f < -1.0e+32\n");
            break;
        }
        if (fabs(actred) <= 0 && prered <= 0)
        {
            info("WARNING: actred and prered <= 0\n");
            break;
        }
        if (fabs(actred) <= 1.0e-12*fabs(f) &&
            fabs(prered) <= 1.0e-12*fabs(f))
        {
            info("WARNING: actred and prered too small\n");
            break;
        }
    }
    
    delete[] g;
    delete[] r;
    delete[] w_new;
    delete[] s;
}

int TRON::trcg(double delta, double *g, double *s, double *r)
{
    int i, inc = 1;
    int n = fun_obj->get_nr_variable();
    double one = 1;
    double *d = new double[n];
    double *Hd = new double[n];
    double rTr, rnewTrnew, alpha, beta, cgtol;
    
    for (i=0; i<n; i++)
    {
        s[i] = 0;
        r[i] = -g[i];
        d[i] = r[i];
    }
    cgtol = 0.1 * dnrm2_(&n, g, &inc);
    if (isnan(cgtol)) {
        Printf("Epic fail. cgtol is NaN");
    }
    
    int cg_iter = 0;
    rTr = ddot_(&n, r, &inc, r, &inc);
    while (1)
    {
        if (dnrm2_(&n, r, &inc) <= cgtol)
            break;
        cg_iter++;
        fun_obj->Hv(d, Hd);
        
        alpha = rTr/ddot_(&n, d, &inc, Hd, &inc);
        daxpy_(&n, &alpha, d, &inc, s, &inc);
        if (dnrm2_(&n, s, &inc) > delta)
        {
            info("cg reaches trust region boundary\n");
            alpha = -alpha;
            daxpy_(&n, &alpha, d, &inc, s, &inc);
            
            double std = ddot_(&n, s, &inc, d, &inc);
            double sts = ddot_(&n, s, &inc, s, &inc);
            double dtd = ddot_(&n, d, &inc, d, &inc);
            double dsq = delta*delta;
            double rad = sqrt(std*std + dtd*(dsq-sts));
            if (std >= 0)
                alpha = (dsq - sts)/(std + rad);
            else
                alpha = (rad - std)/dtd;
            daxpy_(&n, &alpha, d, &inc, s, &inc);
            alpha = -alpha;
            daxpy_(&n, &alpha, Hd, &inc, r, &inc);
            break;
        }
        alpha = -alpha;
        daxpy_(&n, &alpha, Hd, &inc, r, &inc);
        rnewTrnew = ddot_(&n, r, &inc, r, &inc);
        beta = rnewTrnew/rTr;
        dscal_(&n, &beta, d, &inc);
        daxpy_(&n, &one, r, &inc, d, &inc);
        rTr = rnewTrnew;
    }
    
    delete[] d;
    delete[] Hd;
    
    return(cg_iter);
}

double TRON::norm_inf(int n, double *x)
{
    double dmax = fabs(x[0]);
    for (int i=1; i<n; i++)
        if (fabs(x[i]) >= dmax)
            dmax = fabs(x[i]);
    return(dmax);
}

void TRON::set_print_string(void (*print_string) (const char *buf))
{
    tron_print_string = print_string;
}
//-----------------------------End TRON-------------------------------------------------------
//-----------------------------Start SVM------------------------------------------------------
struct feature_node {
    int index;
    double value;
};

struct problem {
    int l, n;
    double *y;
    struct feature_node **x;
    double bias;            /* < 0 if no bias term */
};

enum { L2R_LR, L2R_L2LOSS_SVC_DUAL, L2R_L2LOSS_SVC, L2R_L1LOSS_SVC_DUAL, MCSVM_CS, L1R_L2LOSS_SVC, L1R_LR, L2R_LR_DUAL, L2R_L2LOSS_SVR = 11, L2R_L2LOSS_SVR_DUAL, L2R_L1LOSS_SVR_DUAL }; /* solver_type */

struct parameter {
    int solver_type;
    
    /* these are for training only */
    double eps;	        /* stopping criteria */
    double C;
    int nr_weight;
    int *weight_label;
    double* weight;
    double p;
};

struct model {
    struct parameter param;
    int nr_class;		/* number of classes */
    int nr_feature;
    double *w;
    int *label;		/* label of each class */
    double bias;
};
//---------------------------------------------------------------------------------------------
struct model* warm_start_train(const struct problem *prob, const struct parameter *param, const struct model *wsmodel);
int check_regression_model(const struct model *model);
int check_probability_model(const struct model *model);
double predict(const struct model *model_, const struct feature_node *x);

typedef signed char schar;
template <class S, class T> static inline void clone(T*& dst, S* src, int n)
{
    dst = new T[n];
    memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define INF HUGE_VAL

static void print_string_stdout(const char *s)
{
    fputs(s,stdout);
    fflush(stdout);
}

static void (*liblinear_print_string) (const char *) = &print_string_stdout;

#if 1
static void info(const char *fmt,...)
{
    char buf[BUFSIZ];
    va_list ap;
    va_start(ap,fmt);
    vsprintf(buf,fmt,ap);
    va_end(ap);
    (*liblinear_print_string)(buf);
}
#else
static void info(const char *fmt,...) {}
#endif

class l2r_lr_fun: public tronFunction
{
public:
    l2r_lr_fun(const problem *prob, double *C);
    ~l2r_lr_fun();
    
    double fun(double *w);
    void grad(double *w, double *g);
    void Hv(double *s, double *Hs);
    
    int get_nr_variable(void);
    
private:
    void Xv(double *v, double *Xv);
    void XTv(double *v, double *XTv);
    
    double *C;
    double *z;
    double *D;
    const problem *prob;
};

l2r_lr_fun::l2r_lr_fun(const problem *prob, double *C)
{
    int l=prob->l;
    
    this->prob = prob;
    
    z = new double[l];
    D = new double[l];
    this->C = C;
}

l2r_lr_fun::~l2r_lr_fun()
{
    delete[] z;
    delete[] D;
}


double l2r_lr_fun::fun(double *w)
{
    int i;
    double f=0;
    double *y=prob->y;
    int l=prob->l;
    int w_size=get_nr_variable();
    
    Xv(w, z);
    
    for(i=0;i<w_size;i++)
        f += w[i]*w[i];
    f /= 2.0;
    for(i=0;i<l;i++)
    {
        double yz = y[i]*z[i];
        if (yz >= 0)
            f += C[i]*log(1 + exp(-yz));
        else
            f += C[i]*(-yz+log(1 + exp(yz)));
    }
    
    return(f);
}

void l2r_lr_fun::grad(double *w, double *g)
{
    int i;
    double *y=prob->y;
    int l=prob->l;
    int w_size=get_nr_variable();
    
    for(i=0;i<l;i++)
    {
        z[i] = 1/(1 + exp(-y[i]*z[i]));
        D[i] = z[i]*(1-z[i]);
        z[i] = C[i]*(z[i]-1)*y[i];
    }
    XTv(z, g);
    
    for(i=0;i<w_size;i++)
        g[i] = w[i] + g[i];
}

int l2r_lr_fun::get_nr_variable(void)
{
    return prob->n;
}

void l2r_lr_fun::Hv(double *s, double *Hs)
{
    int i;
    int l=prob->l;
    int w_size=get_nr_variable();
    double *wa = new double[l];
    
    Xv(s, wa);
    for(i=0;i<l;i++)
        wa[i] = C[i]*D[i]*wa[i];
    
    XTv(wa, Hs);
    for(i=0;i<w_size;i++)
        Hs[i] = s[i] + Hs[i];
    delete[] wa;
}

void l2r_lr_fun::Xv(double *v, double *Xv)
{
    int i;
    int l=prob->l;
    feature_node **x=prob->x;
    
    for(i=0;i<l;i++)
    {
        feature_node *s=x[i];
        Xv[i]=0;
        while(s->index!=-1)
        {
            Xv[i]+=v[s->index-1]*s->value;
            s++;
        }
    }
}

void l2r_lr_fun::XTv(double *v, double *XTv)
{
    int i;
    int l=prob->l;
    int w_size=get_nr_variable();
    feature_node **x=prob->x;
    
    for(i=0;i<w_size;i++)
        XTv[i]=0;
    for(i=0;i<l;i++)
    {
        feature_node *s=x[i];
        while(s->index!=-1)
        {
            XTv[s->index-1]+=v[i]*s->value;
            s++;
        }
    }
}

class l2r_l2_svc_fun: public tronFunction
{
public:
    l2r_l2_svc_fun(const problem *prob, double *C);
    ~l2r_l2_svc_fun();
    
    double fun(double *w);
    void grad(double *w, double *g);
    void Hv(double *s, double *Hs);
    
    int get_nr_variable(void);
    
protected:
    void Xv(double *v, double *Xv);
    void subXv(double *v, double *Xv);
    void subXTv(double *v, double *XTv);
    
    double *C;
    double *z;
    double *D;
    int *I;
    int sizeI;
    const problem *prob;
};

l2r_l2_svc_fun::l2r_l2_svc_fun(const problem *prob, double *C)
{
    int l=prob->l;
    
    this->prob = prob;
    
    z = new double[l];
    D = new double[l];
    I = new int[l];
    this->C = C;
}

l2r_l2_svc_fun::~l2r_l2_svc_fun()
{
    delete[] z;
    delete[] D;
    delete[] I;
}

double l2r_l2_svc_fun::fun(double *w)
{
    int i;
    double f=0;
    double *y=prob->y;
    int l=prob->l;
    int w_size=get_nr_variable();
    
    Xv(w, z);
    
    for(i=0;i<w_size;i++)
        f += w[i]*w[i];
    f /= 2.0;
    for(i=0;i<l;i++)
    {
        z[i] = y[i]*z[i];
        double d = 1-z[i];
        if (d > 0)
            f += C[i]*d*d;
    }
    
    return(f);
}

void l2r_l2_svc_fun::grad(double *w, double *g)
{
    int i;
    double *y=prob->y;
    int l=prob->l;
    int w_size=get_nr_variable();
    
    sizeI = 0;
    for (i=0;i<l;i++)
        if (z[i] < 1)
        {
            z[sizeI] = C[i]*y[i]*(z[i]-1);
            I[sizeI] = i;
            sizeI++;
        }
    subXTv(z, g);
    
    for(i=0;i<w_size;i++) {
        g[i] = w[i] + 2*g[i];
        
        if (isnan(g[i])) {
            Printf("Epic fail at: %i, w: %.2f, g: %.2f", i, w[i], g[i]);
            g[i] = 0;
        }
    }
}

int l2r_l2_svc_fun::get_nr_variable(void)
{
    return prob->n;
}

void l2r_l2_svc_fun::Hv(double *s, double *Hs)
{
    int i;
    int w_size=get_nr_variable();
    double *wa = new double[sizeI];
    
    subXv(s, wa);
    for(i=0;i<sizeI;i++)
        wa[i] = C[I[i]]*wa[i];
    
    subXTv(wa, Hs);
    for(i=0;i<w_size;i++)
        Hs[i] = s[i] + 2*Hs[i];
    delete[] wa;
}

void l2r_l2_svc_fun::Xv(double *v, double *Xv)
{
    int i;
    int l=prob->l;
    feature_node **x=prob->x;
    
    for(i=0;i<l;i++)
    {
        feature_node *s=x[i];
        Xv[i]=0;
        while(s->index!=-1)
        {
            Xv[i]+=v[s->index-1]*s->value;
            s++;
        }
    }
}

void l2r_l2_svc_fun::subXv(double *v, double *Xv)
{
    int i;
    feature_node **x=prob->x;
    
    for(i=0;i<sizeI;i++)
    {
        feature_node *s=x[I[i]];
        Xv[i]=0;
        while(s->index!=-1)
        {
            Xv[i]+=v[s->index-1]*s->value;
            s++;
        }
    }
}

void l2r_l2_svc_fun::subXTv(double *v, double *XTv)
{
    int i;
    int w_size=get_nr_variable();
    feature_node **x=prob->x;
    
    for(i=0;i<w_size;i++)
        XTv[i]=0;
    for(i=0;i<sizeI;i++)
    {
        feature_node *s=x[I[i]];
        while(s->index!=-1)
        {
            XTv[s->index-1]+=v[i]*s->value;
            s++;
        }
    }
}

class l2r_l2_svr_fun: public l2r_l2_svc_fun
{
public:
    l2r_l2_svr_fun(const problem *prob, double *C, double p);
    
    double fun(double *w);
    void grad(double *w, double *g);
    
private:
    double p;
};

l2r_l2_svr_fun::l2r_l2_svr_fun(const problem *prob, double *C, double p):
l2r_l2_svc_fun(prob, C)
{
    this->p = p;
}

double l2r_l2_svr_fun::fun(double *w)
{
    int i;
    double f=0;
    double *y=prob->y;
    int l=prob->l;
    int w_size=get_nr_variable();
    double d;
    
    Xv(w, z);
    
    for(i=0;i<w_size;i++)
        f += w[i]*w[i];
    f /= 2;
    for(i=0;i<l;i++)
    {
        d = z[i] - y[i];
        if(d < -p)
            f += C[i]*(d+p)*(d+p);
        else if(d > p)
            f += C[i]*(d-p)*(d-p);
    }
    
    return(f);
}

void l2r_l2_svr_fun::grad(double *w, double *g)
{
    int i;
    double *y=prob->y;
    int l=prob->l;
    int w_size=get_nr_variable();
    double d;
    
    sizeI = 0;
    for(i=0;i<l;i++)
    {
        d = z[i] - y[i];
        
        // generate index set I
        if(d < -p)
        {
            z[sizeI] = C[i]*(d+p);
            I[sizeI] = i;
            sizeI++;
        }
        else if(d > p)
        {
            z[sizeI] = C[i]*(d-p);
            I[sizeI] = i;
            sizeI++;
        }
        
    }
    subXTv(z, g);
    
    for(i=0;i<w_size;i++)
        g[i] = w[i] + 2*g[i];
}

// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
static void group_classes(const problem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm) {
    int l = prob->l;
    int max_nr_class = 16;
    int nr_class = 0;
    int *label = Malloc(int,max_nr_class);
    int *count = Malloc(int,max_nr_class);
    int *data_label = Malloc(int,l);
    int i;
    
    for(i=0;i<l;i++)
    {
        int this_label = (int)prob->y[i];
        int j;
        for(j=0;j<nr_class;j++)
        {
            if(this_label == label[j])
            {
                ++count[j];
                break;
            }
        }
        data_label[i] = j;
        if(j == nr_class)
        {
            if(nr_class == max_nr_class)
            {
                max_nr_class *= 2;
                label = (int *)realloc(label,max_nr_class*sizeof(int));
                count = (int *)realloc(count,max_nr_class*sizeof(int));
            }
            label[nr_class] = this_label;
            count[nr_class] = 1;
            ++nr_class;
        }
    }
    
    //
    // Labels are ordered by their first occurrence in the training set.
    // However, for two-class sets with -1/+1 labels and -1 appears first,
    // we swap labels to ensure that internally the binary SVM has positive data corresponding to the +1 instances.
    //
    if (nr_class == 2 && label[0] == -1 && label[1] == 1)
    {
        swap(label[0],label[1]);
        swap(count[0],count[1]);
        for(i=0;i<l;i++)
        {
            if(data_label[i] == 0)
                data_label[i] = 1;
            else
                data_label[i] = 0;
        }
    }
    
    int *start = Malloc(int,nr_class);
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+count[i-1];
    for(i=0;i<l;i++)
    {
        perm[start[data_label[i]]] = i;
        ++start[data_label[i]];
    }
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+count[i-1];
    
    *nr_class_ret = nr_class;
    *label_ret = label;
    *start_ret = start;
    *count_ret = count;
    free(data_label);
}

static void train_one(const problem *prob, const parameter *param, double *w, double Cp, double Cn)
{
    double eps=param->eps;
    int pos = 0;
    int neg = 0;
    for(int i=0;i<prob->l;i++)
        if(prob->y[i] > 0)
            pos++;
    neg = prob->l - pos;
    
    double primal_solver_tol = eps*max(min(pos,neg), 1)/prob->l;
    
    tronFunction *fun_obj=NULL;
    switch(param->solver_type)
    {
        case L2R_LR:
        {
            double *C = new double[prob->l];
            for(int i = 0; i < prob->l; i++)
            {
                if(prob->y[i] > 0)
                    C[i] = Cp;
                else
                    C[i] = Cn;
            }
            fun_obj=new l2r_lr_fun(prob, C);
            TRON tron_obj(fun_obj, primal_solver_tol);
            tron_obj.set_print_string(liblinear_print_string);
            tron_obj.tron(w);
            delete fun_obj;
            delete[] C;
            break;
        }
        case L2R_L2LOSS_SVC:
        {
            double *C = new double[prob->l];
            for(int i = 0; i < prob->l; i++)
            {
                if(prob->y[i] > 0)
                    C[i] = Cp;
                else
                    C[i] = Cn;
            }
            fun_obj=new l2r_l2_svc_fun(prob, C);
            TRON tron_obj(fun_obj, primal_solver_tol);
            tron_obj.set_print_string(liblinear_print_string);
            tron_obj.tron(w);
            delete fun_obj;
            delete[] C;
            break;
        }
        default:
            fprintf(stderr, "ERROR: unknown solver_type\n");
            break;
    }
}

struct label_index
{
    int label;
    int index;
};
int label_compare(const void *a, const void *b)
{
    int a_label = ((struct label_index*)a)->label;
    int b_label = ((struct label_index*)b)->label;
    if(a_label > b_label)
        return 1;
    else if(a_label == b_label)
        return 0;
    else
        return -1;
}

//
// Interface functions
//
model* train(const problem *prob, const parameter *param)
{
    return warm_start_train(prob, param, NULL);
}

model* warm_start_train(const problem *prob, const parameter *param, const model *wsmodel)
{
    int i,j;
    int l = prob->l;
    int n = prob->n;
    int w_size = prob->n;
    model *model_ = Malloc(model,1);
    
    if(prob->bias>=0)
        model_->nr_feature=n-1;
    else
        model_->nr_feature=n;
    model_->param = *param;
    model_->bias = prob->bias;
    
    if(check_regression_model(model_))
    {
        model_->w = Malloc(double, w_size);
        model_->nr_class = 2;
        model_->label = NULL;
        train_one(prob, param, &model_->w[0], 0, 0);
    }
    else
    {
        int nr_class;
        int *label = NULL;
        int *start = NULL;
        int *count = NULL;
        int *perm = Malloc(int,l);
        
        // group training data of the same class
        group_classes(prob,&nr_class,&label,&start,&count,perm);
        
        model_->nr_class=nr_class;
        model_->label = Malloc(int,nr_class);
        for(i=0;i<nr_class;i++)
            model_->label[i] = label[i];
        
        // calculate weighted C
        double *weighted_C = Malloc(double, nr_class);
        for(i=0;i<nr_class;i++)
            weighted_C[i] = param->C;
        for(i=0;i<param->nr_weight;i++)
        {
            for(j=0;j<nr_class;j++)
                if(param->weight_label[i] == label[j])
                    break;
            if(j == nr_class)
                fprintf(stderr,"WARNING: class label %d specified in weight is not found\n", param->weight_label[i]);
            else
                weighted_C[j] *= param->weight[i];
        }
        
        // constructing the subproblem
        feature_node **x = Malloc(feature_node *,l);
        for(i=0;i<l;i++)
            x[i] = prob->x[perm[i]];
        
        int k;
        problem sub_prob;
        sub_prob.l = l;
        sub_prob.n = n;
        sub_prob.x = Malloc(feature_node *,sub_prob.l);
        sub_prob.y = Malloc(double,sub_prob.l);
        
        for(k=0; k<sub_prob.l; k++)
            sub_prob.x[k] = x[k];
        
        
        if(nr_class == 2) {
            model_->w=Malloc(double, w_size);
            
            int e0 = start[0]+count[0];
            k=0;
            for(; k<e0; k++)
                sub_prob.y[k] = +1;
            for(; k<sub_prob.l; k++)
                sub_prob.y[k] = -1;
            
            if(wsmodel != NULL)
            {
                int min_nr_feature = min(w_size, wsmodel->nr_feature);
                if(wsmodel->label[0] == model_->label[0])
                    for(i=0;i<min_nr_feature;i++)
                        model_->w[i] = wsmodel->w[i];
                else
                    for(i=0;i<min_nr_feature;i++)
                        model_->w[i] = -wsmodel->w[i];
                for(i=min_nr_feature;i<w_size;i++)
                    model_->w[i] = 0;
            }
            else
                for(i=0;i<w_size;i++)
                    model_->w[i] = 0;
            
            train_one(&sub_prob, param, &model_->w[0], weighted_C[0], weighted_C[1]);
        }
        else
        {
            model_->w=Malloc(double, w_size*nr_class);
            double *w=Malloc(double, w_size);
            
            int min_nr_feature = w_size;
            int nr_matched_label = 0;
            struct label_index *label_map = NULL;
            if(wsmodel != NULL)
            {
                min_nr_feature = min(w_size, wsmodel->nr_feature);
                label_map = Malloc(struct label_index, wsmodel->nr_class);
                for(i=0;i<wsmodel->nr_class;i++)
                {
                    label_map[i].label = wsmodel->label[i];
                    label_map[i].index = i;
                }
                qsort(label_map, wsmodel->nr_class, sizeof(struct label_index), label_compare);
            }
            
            for(i=0;i<nr_class;i++)
            {
                int si = start[i];
                int ei = si+count[i];
                
                k=0;
                for(; k<si; k++)
                    sub_prob.y[k] = -1;
                for(; k<ei; k++)
                    sub_prob.y[k] = +1;
                for(; k<sub_prob.l; k++)
                    sub_prob.y[k] = -1;
                
                int index = -1;
                if(wsmodel != NULL)
                {
                    struct label_index key;
                    struct label_index *found;
                    key.label = model_->label[i];
                    found = (struct label_index*)bsearch(&key, label_map, wsmodel->nr_class, sizeof(struct label_index), label_compare);
                    if(found != NULL)
                        index = found->index;
                }
                if(index >= 0)
                {
                    for(j=0;j<min_nr_feature;j++)
                        w[j] = wsmodel->w[j*wsmodel->nr_class+index];
                    for(j=min_nr_feature;j<w_size;j++)
                        w[j] = 0;
                    nr_matched_label++;
                }
                else
                    for(j=0;j<w_size;j++)
                        w[j] = 0;
                
                train_one(&sub_prob, param, w, weighted_C[i], param->C);
                
                for(int j=0;j<w_size;j++)
                    model_->w[j*nr_class+i] = w[j];
            }
            free(w);
            if(wsmodel != NULL)
            {
                if(nr_matched_label != nr_class || nr_class != wsmodel->nr_class)
                    fprintf(stderr,"WARNING: class labels in training data do not match those in the initial model.\n");
                free(label_map);
            }
        }
        
        free(x);
        free(label);
        free(start);
        free(count);
        free(perm);
        free(sub_prob.x);
        free(sub_prob.y);
        free(weighted_C);
    }
    return model_;
}

double predict_values(const struct model *model_, const struct feature_node *x, double *dec_values)
{
    int idx;
    int n;
    if(model_->bias>=0)
        n=model_->nr_feature+1;
    else
        n=model_->nr_feature;
    double *w=model_->w;
    int nr_class=model_->nr_class;
    int i;
    int nr_w;
    if(nr_class==2 && model_->param.solver_type != MCSVM_CS)
        nr_w = 1;
    else
        nr_w = nr_class;
    
    const feature_node *lx=x;
    for(i=0;i<nr_w;i++)
        dec_values[i] = 0;
    for(; (idx=lx->index)!=-1; lx++)
    {
        // the dimension of testing data may exceed that of training
        if(idx<=n)
            for(i=0;i<nr_w;i++)
                dec_values[i] += w[(idx-1)*nr_w+i]*lx->value;
    }
    
    if(nr_class==2)
    {
        if(check_regression_model(model_))
            return dec_values[0];
        else
            return (dec_values[0]>0)?model_->label[0]:model_->label[1];
    }
    else
    {
        int dec_max_idx = 0;
        for(i=1;i<nr_class;i++)
        {
            if(dec_values[i] > dec_values[dec_max_idx])
                dec_max_idx = i;
        }
        return model_->label[dec_max_idx];
    }
}

double predict(const model *model_, const feature_node *x)
{
    double *dec_values = Malloc(double, model_->nr_class);
    double label=predict_values(model_, x, dec_values);
    free(dec_values);
    return label;
}

double predict_probability(const struct model *model_, const struct feature_node *x, double* prob_estimates)
{
    if(check_probability_model(model_))
    {
        int i;
        int nr_class=model_->nr_class;
        int nr_w;
        if(nr_class==2)
            nr_w = 1;
        else
            nr_w = nr_class;
        
        double label=predict_values(model_, x, prob_estimates);
        for(i=0;i<nr_w;i++)
            prob_estimates[i]=1/(1+exp(-prob_estimates[i]));
        
        if(nr_class==2) // for binary classification
            prob_estimates[1]=1.-prob_estimates[0];
        else
        {
            double sum=0;
            for(i=0; i<nr_class; i++)
                sum+=prob_estimates[i];
            
            for(i=0; i<nr_class; i++)
                prob_estimates[i]=prob_estimates[i]/sum;
        }
        
        return label;
    }
    else
        return 0;
}

void free_model_content(struct model *model_ptr)
{
    if(model_ptr->w != NULL)
        free(model_ptr->w);
    if(model_ptr->label != NULL)
        free(model_ptr->label);
}

void free_and_destroy_model(struct model **model_ptr_ptr)
{
    struct model *model_ptr = *model_ptr_ptr;
    if(model_ptr != NULL)
    {
        free_model_content(model_ptr);
        free(model_ptr);
    }
}

void destroy_param(parameter* param)
{
    if(param->weight_label != NULL)
        free(param->weight_label);
    if(param->weight != NULL)
        free(param->weight);
}

const char *check_parameter(const problem *prob, const parameter *param) {
    if(param->eps <= 0)
        return "eps <= 0";
    
    if(param->C <= 0)
        return "C <= 0";
    
    if(param->p < 0)
        return "p < 0";
    
    if(param->solver_type != L2R_LR
       && param->solver_type != L2R_L2LOSS_SVC)
        return "unknown solver type";
    
    return NULL;
}

int check_probability_model(const struct model *model_)
{
    return (model_->param.solver_type==L2R_LR ||
            model_->param.solver_type==L2R_LR_DUAL ||
            model_->param.solver_type==L1R_LR);
}

int check_regression_model(const struct model *model_)
{
    return (model_->param.solver_type==L2R_L2LOSS_SVR ||
            model_->param.solver_type==L2R_L1LOSS_SVR_DUAL ||
            model_->param.solver_type==L2R_L2LOSS_SVR_DUAL);
}

void set_print_string_function(void (*print_func)(const char*))
{
    if (print_func == NULL)
        liblinear_print_string = &print_string_stdout;
    else
        liblinear_print_string = print_func;
}

//-----------------------------End SVM---------------------------------------------------------

void readProblem(const VVD &samples, const VD &dv, struct problem *prob, struct feature_node *x_space) {
    prob->l = (int)samples.size();
    prob->n = (int)samples[0].size();
    prob->y = Malloc(double, prob->l);
    size_t elements = prob->l * prob->n;
    prob->x = Malloc(struct feature_node *, prob->l);
    x_space = Malloc(struct feature_node, elements + prob->l);
    
    // fill with values
    int j = 0;
    for (int i = 0; i < prob->l; i++) {
        // set Y
        prob->y[i] = dv[i];
        
        // set X
        prob->x[i] = &x_space[j];
        for (int k = 0; k < prob->n; k++) {
            x_space[j].index = k + 1;// starts from 1
            x_space[j].value = samples[i][k];
            
            //            Printf("%i : %.2f\n", x_space[j].index, x_space[j].value);
            
            ++j;
        }
        // mark end of row
        x_space[j++].index = -1;
    }
    Printf("Problem parsed, samples: %i, features: %i\n", prob->l, prob->n);
}

#endif
