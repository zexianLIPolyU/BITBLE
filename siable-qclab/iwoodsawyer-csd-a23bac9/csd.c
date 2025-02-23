/*
 * Cosine Sine Decomposition (CSD)
 *
 * [C,S,U1,V1,U2,V2] = csd(X11,X12,X21,X22)
 *
 * compile command:
 * mex -O csd.c libmwblas.lib libmwlapack.lib (>= R2012A)
 *
 * calls the SORCSD/DORCSD/CUNCSD/ZUNCSD named LAPACK function
 *
 * Ivo Houtzager
 */

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "csd.h"

void csd_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *X11p, *X11pr, *X12p, *X12pr, *X21p, *X21pr, *X22p, *X22pr;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    double *X11pi, *X12pi, *X21pi, *X22pi;
    #endif
    double *U1p, *U1pr, *U2p, *U2pr, *Cpr, *Spr, *V1Tp, *V1pr, *V2Tp, *V2pr;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    double *U1pi, *U2pi, *V1pi, *V2pi;
    #endif
    double *theta, *work, *rwork, size, rsize;
    mwSignedIndex lwork, lrwork, *iwork, info = 1, cplx = 0, dc = 1;
    mwSignedIndex i, j, p, q, m, r;
    mwSignedIndex ldu1, ldu2, ldv1t, ldv2t, ldx11, ldx12, ldx21, ldx22;
    mwSize mx11, nx11, mx12, nx12, mx21, nx21, mx22, nx22;
    char jobu1='Y', jobu2='Y', jobv1t='Y', jobv2t='Y', trans='O', signs='D';
    mxClassID classid = mxDOUBLE_CLASS;
    mxComplexity cplxflag = mxREAL;
    size_t element_size = sizeof(double);
    
    /* check complex */
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) ||
            mxIsComplex(prhs[2]) || mxIsComplex(prhs[3])) {
        if (!mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) ||
                !mxIsComplex(prhs[2]) || !mxIsComplex(prhs[3])) {
            mexErrMsgTxt("All inputs must be real or all must be complex.");
        }
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }
    
    /* check for proper dimensions*/
    mx11 = mxGetM(prhs[0]);
    nx11 = mxGetN(prhs[0]);
    mx12 = mxGetM(prhs[1]);
    nx12 = mxGetN(prhs[1]);
    mx21 = mxGetM(prhs[2]);
    nx21 = mxGetN(prhs[2]);
    mx22 = mxGetM(prhs[3]);
    nx22 = mxGetN(prhs[3]);
    if (mx11 != mx12) {
        mexErrMsgTxt("Number of rows of X11 and X12 must be equal." );
    }
    if (mx21 != mx22) {
        mexErrMsgTxt("Number of rows of X21 and X22 must be equal." );
    }
    if (nx11 != nx21) {
        mexErrMsgTxt("Number of columns of X11 and X21 must be equal." );
    }
    if (nx12 != nx22) {
        mexErrMsgTxt("Number of columns of X12 and X22 must be equal." );
    }
    if (mx11+mx21 != nx11+nx12) {
        mexErrMsgTxt("The matrix [X11 X12; X21 X22] must be square." );
    }

    /* leading dimensions */
    p = mx11;
    q = nx11;
    m = mx11+mx21;
    
    /* allocate output matrices */
    if (nlhs < 3) {
        jobu1 = 'O';
        ldu1 = 0;
        U1p = NULL;
    }
    else {
        ldu1  = max(1,p);
        U1p = mxCalloc(dc*ldu1*ldu1,element_size);
    }
    if (nlhs < 5) {
        jobu2 = 'O';
        ldu2 = 0;
        U2p = NULL;
    }
    else {
        ldu2  = max(1,m-p);
        U2p = mxCalloc(dc*ldu2*ldu2,element_size);
    }
    if (nlhs < 4) {
        jobv1t = 'O';
        ldv1t = 0;
        V1Tp = NULL;
    }
    else {
        ldv1t = max(1,q);
        V1Tp = mxCalloc(dc*ldv1t*ldv1t,element_size);
    }
    if (nlhs < 6) {
        jobv2t = 'O';
        ldv2t = 0;
        V2Tp = NULL;
    }
    else {
        ldv2t = max(1,m-q);
        V2Tp = mxCalloc(dc*ldv2t*ldv2t,element_size);
    }   
    
    /* deep copy input matrices */
    ldx11 = max(1,mx11);
    ldx12 = max(1,mx12);
    ldx21 = max(1,mx21);
    ldx22 = max(1,mx22);
    X11p = mxMalloc(dc*ldx11*nx11*element_size);
    X12p = mxMalloc(dc*ldx12*nx12*element_size);
    X21p = mxMalloc(dc*ldx21*nx21*element_size);
    X22p = mxMalloc(dc*ldx22*nx22*element_size);
    X11pr = mxGetData(prhs[0]);
    X12pr = mxGetData(prhs[1]);
    X21pr = mxGetData(prhs[2]);
    X22pr = mxGetData(prhs[3]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        X11pi = mxGetImagData(prhs[0]);
        X12pi = mxGetImagData(prhs[1]);
        X21pi = mxGetImagData(prhs[2]);
        X22pi = mxGetImagData(prhs[3]);
        for (j=0; j<nx11; j++) {
            for (i=0; i<mx11; i++) {
                X11p[j*2*mx11+2*i]   = X11pr[j*mx11+i];
                X11p[j*2*mx11+2*i+1] = X11pi[j*mx11+i];
            }
        }
        for (j=0; j<nx12; j++) {
            for (i=0; i<mx12; i++) {
                X12p[j*2*mx12+2*i]   = X12pr[j*mx12+i];
                X12p[j*2*mx12+2*i+1] = X12pi[j*mx12+i];
            }
        }
        for (j=0; j<nx21; j++) {
            for (i=0; i<mx21; i++) {
                X21p[j*2*mx21+2*i]   = X21pr[j*mx21+i];
                X21p[j*2*mx21+2*i+1] = X21pi[j*mx21+i];
            }
        }
        for (j=0; j<nx22; j++) {
            for (i=0; i<mx22; i++) {
                X22p[j*2*mx22+2*i]   = X22pr[j*mx22+i];
                X22p[j*2*mx22+2*i+1] = X22pi[j*mx22+i];
            }
        }
    }
    else {
    #endif
        memcpy(X11p,X11pr,dc*mx11*nx11*element_size);
        memcpy(X12p,X12pr,dc*mx12*nx12*element_size);
        memcpy(X21p,X21pr,dc*mx21*nx21*element_size);
        memcpy(X22p,X22pr,dc*mx22*nx22*element_size);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif
    
    /* allocate iwork matrix */
    r = max(1,min(min(p,m-p),min(q,m-q)));
    theta = mxMalloc(r*element_size);
    iwork = mxCalloc(max(1,m-r),sizeof(mwSignedIndex));
    
    /* determine blocksize */
    lwork = -1;
    if (cplx) {
        lrwork = -1;
        zuncsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                &size, &lwork, &rsize, &lrwork, iwork, &info);
    }
    else {
        dorcsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                &size, &lwork, iwork, &info);
    }
    if (info != 0) {
        mxFree(theta);
        mxFree(iwork);
        mxFree(X11p);
        mxFree(X12p);
        mxFree(X21p);
        mxFree(X22p);
        mxFree(U1p);
        mxFree(U2p);
        mxFree(V1Tp);
        mxFree(V2Tp);
        if (cplx) {
            mexPrintf("ZUNCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("ZUNCSD not successful.");
        }
        else {
            mexPrintf("DORCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("DORCSD not successful.");
        }
    }
    
    /* workaround with issue MKL and visual studio */
    size = size + 1;
       
    /* allocate work matrix */
    lwork = max(1,(mwSignedIndex)size);
    work = mxCalloc(dc*lwork,element_size);
    if (cplx) {
        if (rsize > 0)
            lrwork = max(1,(mwSignedIndex)rsize);
        else {
            /* We compute the size of lrwork because LAPACK does not return it on a query call. */
            lrwork = 2 + 9*ldv1t + 9*ldv2t + max(1,8*p) + max(1,8*q);
        }
        rwork = mxCalloc(lrwork,element_size);
    }
    
    /* calls the DORCSD/ZUNCSD function */
    if (cplx) {
        zuncsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                work, &lwork, rwork, &lrwork, iwork, &info);
        mxFree(rwork);
    }
    else {
        dorcsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                work, &lwork, iwork, &info);
    }
    mxFree(work);
    mxFree(iwork);
    mxFree(X11p);
    mxFree(X12p);
    mxFree(X21p);
    mxFree(X22p);
    if (info < 0) {
        mxFree(theta);
        mxFree(U1p);
        mxFree(U2p);
        mxFree(V1Tp);
        mxFree(V2Tp);
        if (cplx) {
            mexPrintf("ZUNCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("ZUNCSD not successful.");
        }
        else {
            mexPrintf("DORCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("DORCSD not successful.");
        }
    }
    else if (info > 0) {
        if (cplx) {
            mexPrintf("ZUNCSD returned INFO=%d.\n",info);
            mexWarnMsgTxt("ZUNCSD did not converge.");
        }
        else {
            mexPrintf("DORCSD returned INFO=%d.\n",info);
            mexWarnMsgTxt("DORCSD did not converge.");
        }
    }
    
    /* create C and S matrix */
    plhs[0] = mxCreateNumericMatrix(r,r,classid,mxREAL);
    plhs[1] = mxCreateNumericMatrix(r,r,classid,mxREAL);
    Cpr = mxGetData(plhs[0]);
    Spr = mxGetData(plhs[1]);
    for (i=0; i<r; i++) {
        Cpr[i*r+i] = cos(theta[i]);
        Spr[i*r+i] = sin(theta[i]);
    }
    mxFree(theta);
    
    /* create U1 matrix */
    if (nlhs > 2)
    {
        plhs[2] = mxCreateNumericMatrix(p,p,classid,cplxflag);
        U1pr = mxGetData(plhs[2]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            U1pi = mxGetImagData(plhs[2]);
            for (j=0; j<p; j++) {
                for (i=0; i<p; i++) {
                    U1pr[j*p+i] = U1p[j*2*ldu1+2*i];
                    U1pi[j*p+i] = U1p[j*2*ldu1+2*i+1];
                }
            }
        }
        else {
        #endif
            if (ldu1 == p) {
                memcpy(U1pr,U1p,dc*p*p*element_size);
            }
            else {
                for (j=0; j<p; j++) {
                    for (i=0; i<dc*p; i++) {
                        U1pr[j*dc*p+i] = U1p[j*dc*ldu1+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
        mxFree(U1p);
    }
    
    /* create V1 matrix */
    if (nlhs > 3)
    {
        plhs[3] = mxCreateNumericMatrix(q,q,classid,cplxflag);
        V1pr = mxGetData(plhs[3]);
        if (cplx) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            for (j=0; j<q; j++) {
                for (i=0; i<q; i++) {
                    V1pr[i*2*q+2*j]   =  V1Tp[j*2*ldv1t+2*i];
                    V1pr[i*2*q+2*j+1] = -V1Tp[j*2*ldv1t+2*i+1];
                }
            }
            #else
            V1pi = mxGetImagData(plhs[3]);
            for (j=0; j<q; j++) {
                for (i=0; i<q; i++) {
                    V1pr[i*q+j] =  V1Tp[j*2*ldv1t+2*i];
                    V1pi[i*q+j] = -V1Tp[j*2*ldv1t+2*i+1];
                }
            }
            #endif
        }
        else {
            for (j=0; j<q; j++) {
                for (i=0; i<q; i++) {
                    V1pr[i*q+j] = V1Tp[j*ldv1t+i];
                }
            }
        }
        mxFree(V1Tp);
    }
    
    /* create U2 matrix */
    if (nlhs > 4)
    {
        mwSignedIndex l = m-p;
        plhs[4] = mxCreateNumericMatrix(l,l,classid,cplxflag);
        U2pr = mxGetData(plhs[4]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            U2pi = mxGetImagData(plhs[4]);
            for (j=0; j<l; j++) {
                for (i=0; i<l; i++) {
                    U2pr[j*l+i] = U2p[j*2*ldu2+2*i];
                    U2pi[j*l+i] = U2p[j*2*ldu2+2*i+1];
                }
            }
        }
        else {
        #endif
            if (ldu2 == l) {
                memcpy(U2pr,U2p,dc*l*l*element_size);
            }
            else {
                for (j=0; j<l; j++) {
                    for (i=0; i<dc*l; i++) {
                        U2pr[j*dc*l+i] = U2p[j*dc*ldu2+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
        mxFree(U2p);
    }
    
    /* create V2 matrix */
    if (nlhs > 5)
    {
        mwSignedIndex k = m-q;
        plhs[5] = mxCreateNumericMatrix(k,k,classid,cplxflag);
        V2pr = mxGetData(plhs[5]);
        if (cplx) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            for (j=0; j<k; j++) {
                for (i=0; i<k; i++) {
                    V2pr[i*2*k+2*j]   =  V2Tp[j*2*ldv2t+2*i];
                    V2pr[i*2*k+2*j+1] = -V2Tp[j*2*ldv2t+2*i+1];
                }
            }
            #else
            V2pi = mxGetImagData(plhs[5]);
            for (j=0; j<k; j++) {
                for (i=0; i<k; i++) {
                    V2pr[i*k+j] =  V2Tp[j*2*ldv2t+2*i];
                    V2pi[i*k+j] = -V2Tp[j*2*ldv2t+2*i+1];
                }
            }
            #endif
        }
        else {
            for (j=0; j<k; j++) {
                for (i=0; i<k; i++) {
                    V2pr[i*k+j] = V2Tp[j*ldv2t+i];
                }
            }
        }
        mxFree(V2Tp);
    }
}

void csd_single(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float *X11p, *X11pr, *X12p, *X12pr, *X21p, *X21pr, *X22p, *X22pr;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    float *X11pi, *X12pi, *X21pi, *X22pi;
    #endif
    float *U1p, *U1pr, *U2p, *U2pr, *Cpr, *Spr, *V1Tp, *V1pr, *V2Tp, *V2pr;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    float *U1pi, *U2pi, *V1pi, *V2pi;
    #endif
    float *theta, *work, *rwork, size, rsize;
    mwSignedIndex lwork, lrwork, *iwork, info = 1, cplx = 0, dc = 1;
    mwSignedIndex i, j, p, q, m, r;
    mwSignedIndex ldu1, ldu2, ldv1t, ldv2t, ldx11, ldx12, ldx21, ldx22; 
    mwSize mx11, nx11, mx12, nx12, mx21, nx21, mx22, nx22;
    char jobu1='Y', jobu2='Y', jobv1t='Y', jobv2t='Y', trans='O', signs='D';
    mxClassID classid = mxSINGLE_CLASS;
    mxComplexity cplxflag = mxREAL;
    size_t element_size = sizeof(float);
    
    /* check complex */
    if (mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) ||
            mxIsComplex(prhs[2]) || mxIsComplex(prhs[3])) {
        if (!mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) ||
                !mxIsComplex(prhs[2]) || !mxIsComplex(prhs[3])) {
            mexErrMsgTxt("All inputs must be real or all must be complex.");
        }
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }
    
    /* check for proper dimensions*/
    mx11 = mxGetM(prhs[0]);
    nx11 = mxGetN(prhs[0]);
    mx12 = mxGetM(prhs[1]);
    nx12 = mxGetN(prhs[1]);
    mx21 = mxGetM(prhs[2]);
    nx21 = mxGetN(prhs[2]);
    mx22 = mxGetM(prhs[3]);
    nx22 = mxGetN(prhs[3]);
    if (mx11 != mx12) {
        mexErrMsgTxt("Number of rows of X11 and X12 must be equal." );
    }
    if (mx21 != mx22) {
        mexErrMsgTxt("Number of rows of X21 and X22 must be equal." );
    }
    if (nx11 != nx21) {
        mexErrMsgTxt("Number of columns of X11 and X21 must be equal." );
    }
    if (nx12 != nx22) {
        mexErrMsgTxt("Number of columns of X12 and X22 must be equal." );
    }
    if (mx11+mx21 != nx11+nx12) {
        mexErrMsgTxt("The matrix [X11 X12; X21 X22] must be square." );
    }

    /* leading dimensions */
    p = mx11;
    q = nx11;
    m = mx11+mx21;
    
    /* allocate output matrices */
    if (nlhs < 3) {
        jobu1 = 'O';
        ldu1 = 0;
        U1p = NULL;
    }
    else {
        ldu1  = max(1,p);
        U1p = mxCalloc(dc*ldu1*ldu1,element_size);
    }
    if (nlhs < 5) {
        jobu2 = 'O';
        ldu2 = 0;
        U2p = NULL;
    }
    else {
        ldu2  = max(1,m-p);
        U2p = mxCalloc(dc*ldu2*ldu2,element_size);
    }
    if (nlhs < 4) {
        jobv1t = 'O';
        ldv1t = 0;
        V1Tp = NULL;
    }
    else {
        ldv1t = max(1,q);
        V1Tp = mxCalloc(dc*ldv1t*ldv1t,element_size);
    }
    if (nlhs < 6) {
        jobv2t = 'O';
        ldv2t = 0;
        V2Tp = NULL;
    }
    else {
        ldv2t = max(1,m-q);
        V2Tp = mxCalloc(dc*ldv2t*ldv2t,element_size);
    }   
    
    /* deep copy input matrices */
    ldx11 = max(1,mx11);
    ldx12 = max(1,mx12);
    ldx21 = max(1,mx21);
    ldx22 = max(1,mx22);
    X11p = mxMalloc(dc*ldx11*nx11*element_size);
    X12p = mxMalloc(dc*ldx12*nx12*element_size);
    X21p = mxMalloc(dc*ldx21*nx21*element_size);
    X22p = mxMalloc(dc*ldx22*nx22*element_size);
    X11pr = mxGetData(prhs[0]);
    X12pr = mxGetData(prhs[1]);
    X21pr = mxGetData(prhs[2]);
    X22pr = mxGetData(prhs[3]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        X11pi = mxGetImagData(prhs[0]);
        X12pi = mxGetImagData(prhs[1]);
        X21pi = mxGetImagData(prhs[2]);
        X22pi = mxGetImagData(prhs[3]);
        for (j=0; j<nx11; j++) {
            for (i=0; i<mx11; i++) {
                X11p[j*2*mx11+2*i]   = X11pr[j*mx11+i];
                X11p[j*2*mx11+2*i+1] = X11pi[j*mx11+i];
            }
        }
        for (j=0; j<nx12; j++) {
            for (i=0; i<mx12; i++) {
                X12p[j*2*mx12+2*i]   = X12pr[j*mx12+i];
                X12p[j*2*mx12+2*i+1] = X12pi[j*mx12+i];
            }
        }
        for (j=0; j<nx21; j++) {
            for (i=0; i<mx21; i++) {
                X21p[j*2*mx21+2*i]   = X21pr[j*mx21+i];
                X21p[j*2*mx21+2*i+1] = X21pi[j*mx21+i];
            }
        }
        for (j=0; j<nx22; j++) {
            for (i=0; i<mx22; i++) {
                X22p[j*2*mx22+2*i]   = X22pr[j*mx22+i];
                X22p[j*2*mx22+2*i+1] = X22pi[j*mx22+i];
            }
        }
    }
    else {
    #endif
        memcpy(X11p,X11pr,dc*mx11*nx11*element_size);
        memcpy(X12p,X12pr,dc*mx12*nx12*element_size);
        memcpy(X21p,X21pr,dc*mx21*nx21*element_size);
        memcpy(X22p,X22pr,dc*mx22*nx22*element_size);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif
    
    /* allocate iwork matrix */
    r = max(1,min(min(p,m-p),min(q,m-q)));
    theta = mxMalloc(r*element_size);
    iwork = mxCalloc(max(1,m-r),sizeof(mwSignedIndex));
    
    /* determine blocksize */
    lwork = -1;
    if (cplx) {
        lrwork = -1;
        cuncsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                &size, &lwork, &rsize, &lrwork, iwork, &info);
    }
    else {
        sorcsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                &size, &lwork, iwork, &info);
    }
    if (info != 0) {
        mxFree(theta);
        mxFree(iwork);
        mxFree(X11p);
        mxFree(X12p);
        mxFree(X21p);
        mxFree(X22p);
        mxFree(U1p);
        mxFree(U2p);
        mxFree(V1Tp);
        mxFree(V2Tp);
        if (cplx) {
            mexPrintf("CUNCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("CUNCSD not successful.");
        }
        else {
            mexPrintf("SORCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("SORCSD not successful.");
        }
    }
    
    /* workaround with issue MKL and visual studio */
    size = size + 1;
    
    /* allocate work matrix */
    lwork = max(1,(mwSignedIndex)size);
    work = mxCalloc(dc*lwork,element_size);
    if (cplx) {
        if (rsize > 0)
            lrwork = max(1,(mwSignedIndex)rsize);
        else {
            /* We compute the size of lrwork because LAPACK does not return it on a query call. */
            lrwork = 2 + 9*ldv1t + 9*ldv2t + max(1,8*p) + max(1,8*q);
        }
        rwork = mxCalloc(lrwork,element_size);
    }
    
    /* calls the SORCSD/CUNCSD function */
    if (cplx) {
        cuncsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                work, &lwork, rwork, &lrwork, iwork, &info);
        mxFree(rwork);
    }
    else {
        sorcsd( &jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q,
                X11p, &ldx11, X12p, &ldx12, X21p, &ldx21, X22p, &ldx22,
                theta, U1p, &ldu1, U2p, &ldu2, V1Tp, &ldv1t, V2Tp, &ldv2t,
                work, &lwork, iwork, &info);
    }
    mxFree(work);
    mxFree(iwork);
    mxFree(X11p);
    mxFree(X12p);
    mxFree(X21p);
    mxFree(X22p);
    if (info < 0) {
        mxFree(theta);
        mxFree(U1p);
        mxFree(U2p);
        mxFree(V1Tp);
        mxFree(V2Tp);
        if (cplx) {
            mexPrintf("CUNCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("CUNCSD not successful.");
        }
        else {
            mexPrintf("SORCSD returned INFO=%d.\n",info);
            mexErrMsgTxt("SORCSD not successful.");
        }
    }
    else if (info > 0) {
        if (cplx) {
            mexPrintf("CUNCSD returned INFO=%d.\n",info);
            mexWarnMsgTxt("CUNCSD did not converge.");
        }
        else {
            mexPrintf("SORCSD returned INFO=%d.\n",info);
            mexWarnMsgTxt("SORCSD did not converge.");
        }
    }
    
    /* create C and S matrix */
    plhs[0] = mxCreateNumericMatrix(r,r,classid,mxREAL);
    plhs[1] = mxCreateNumericMatrix(r,r,classid,mxREAL);
    Cpr = mxGetData(plhs[0]);
    Spr = mxGetData(plhs[1]);
    for (i=0; i<r; i++) {
        Cpr[i*r+i] = cosf(theta[i]);
        Spr[i*r+i] = sinf(theta[i]);
    }
    mxFree(theta);
    
    /* create U1 matrix */
    if (nlhs > 2)
    {
        plhs[2] = mxCreateNumericMatrix(p,p,classid,cplxflag);
        U1pr = mxGetData(plhs[2]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            U1pi = mxGetImagData(plhs[2]);
            for (j=0; j<p; j++) {
                for (i=0; i<p; i++) {
                    U1pr[j*p+i] = U1p[j*2*ldu1+2*i];
                    U1pi[j*p+i] = U1p[j*2*ldu1+2*i+1];
                }
            }
        }
        else {
        #endif
            if (ldu1 == p) {
                memcpy(U1pr,U1p,dc*p*p*element_size);
            }
            else {
                for (j=0; j<p; j++) {
                    for (i=0; i<dc*p; i++) {
                        U1pr[j*dc*p+i] = U1p[j*dc*ldu1+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
        mxFree(U1p);
    }
    
    /* create V1 matrix */
    if (nlhs > 3)
    {
        plhs[3] = mxCreateNumericMatrix(q,q,classid,cplxflag);
        V1pr = mxGetData(plhs[3]);
        if (cplx) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            for (j=0; j<q; j++) {
                for (i=0; i<q; i++) {
                    V1pr[i*2*q+2*j]   =  V1Tp[j*2*ldv1t+2*i];
                    V1pr[i*2*q+2*j+1] = -V1Tp[j*2*ldv1t+2*i+1];
                }
            }
            #else
            V1pi = mxGetImagData(plhs[3]);
            for (j=0; j<q; j++) {
                for (i=0; i<q; i++) {
                    V1pr[i*q+j] =  V1Tp[j*2*ldv1t+2*i];
                    V1pi[i*q+j] = -V1Tp[j*2*ldv1t+2*i+1];
                }
            }
            #endif
        }
        else {
            for (j=0; j<q; j++) {
                for (i=0; i<q; i++) {
                    V1pr[i*q+j] = V1Tp[j*ldv1t+i];
                }
            }
        }
        mxFree(V1Tp);
    }
    
    /* create U2 matrix */
    if (nlhs > 4)
    {
        mwSignedIndex l = m-p;
        plhs[4] = mxCreateNumericMatrix(l,l,classid,cplxflag);
        U2pr = mxGetData(plhs[4]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            U2pi = mxGetImagData(plhs[4]);
            for (j=0; j<l; j++) {
                for (i=0; i<l; i++) {
                    U2pr[j*l+i] = U2p[j*2*ldu2+2*i];
                    U2pi[j*l+i] = U2p[j*2*ldu2+2*i+1];
                }
            }
        }
        else {
        #endif
            if (ldu2 == l) {
                memcpy(U2pr,U2p,dc*l*l*element_size);
            }
            else {
                for (j=0; j<l; j++) {
                    for (i=0; i<dc*l; i++) {
                        U2pr[j*dc*l+i] = U2p[j*dc*ldu2+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
        mxFree(U2p);
    }
    
    /* create V2 matrix */
    if (nlhs > 5)
    {
        mwSignedIndex k = m-q;
        plhs[5] = mxCreateNumericMatrix(k,k,classid,cplxflag);
        V2pr = mxGetData(plhs[5]);
        if (cplx) {
            #if MX_HAS_INTERLEAVED_COMPLEX
            for (j=0; j<k; j++) {
                for (i=0; i<k; i++) {
                    V2pr[i*2*k+2*j]   =  V2Tp[j*2*ldv2t+2*i];
                    V2pr[i*2*k+2*j+1] = -V2Tp[j*2*ldv2t+2*i+1];
                }
            }
            #else
            V2pi = mxGetImagData(plhs[5]);
            for (j=0; j<k; j++) {
                for (i=0; i<k; i++) {
                    V2pr[i*k+j] =  V2Tp[j*2*ldv2t+2*i];
                    V2pi[i*k+j] = -V2Tp[j*2*ldv2t+2*i+1];
                }
            }
            #endif
        }
        else {
            for (j=0; j<k; j++) {
                for (i=0; i<k; i++) {
                    V2pr[i*k+j] = V2Tp[j*ldv2t+i];
                }
            }
        }
        mxFree(V2Tp);
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs < 4) {
        mexErrMsgTxt("CSD requires four input arguments.");
    }
    if (nlhs < 2) {
        mexErrMsgTxt("CSD requires at least two output arguments.");
    }
    if (nlhs > 6) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input X11 must be a full matrix." );
    }
    if (!mxIsNumeric(prhs[1]) || mxIsSparse(prhs[1])) {
        mexErrMsgTxt( "Input X12 must be a full matrix." );
    }
    if (!mxIsNumeric(prhs[2]) || mxIsSparse(prhs[2])) {
        mexErrMsgTxt( "Input X21 must be a full matrix." );
    }
    if (!mxIsNumeric(prhs[3]) || mxIsSparse(prhs[3])) {
        mexErrMsgTxt( "Input X22 must be a full matrix." );
    }
    if (mxIsDouble(prhs[0]) && mxIsDouble(prhs[1]) &&
            mxIsDouble(prhs[2]) && mxIsDouble(prhs[3])) {
        csd_double(nlhs, plhs, nrhs, prhs);
    }
    else if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1]) &&
            mxIsSingle(prhs[2]) && mxIsSingle(prhs[3])) {
        csd_single(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt("All inputs must be single or all must be double.");
    }
}
