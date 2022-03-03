/******************************************
Fix arrythmias using the inverse gaussian likelihood
Feb 10 2011
Luca Citi

% Copyright (C) Luca Citi and Riccardo Barbieri, 2010-2011.
% All Rights Reserved. See LICENSE.TXT for license details.
% {lciti,barbieri}@neurostat.mit.edu
% http://users.neurostat.mit.edu/barbieri/pphrv

******************************************/


/*
* To compile this mex file manually:
*
* Matlab for Linux, Octave for Linux, and Octave for Windows (with gcc/mingw compiler):
* mex arrythmia_invnorm_mex.c
*
* Matlab for Windows (with the default LCC compiler):
* mex arrythmia_invnorm_mex.c
*/


#include "mex.h"
#include <math.h>

typedef double real;


#if defined(DEBUG)
#define dbg(ARGS) mexPrintf ARGS ;
#else
#define dbg(ARGS) ;
#endif

#if !defined(M_PI)
#define M_PI (3.14159265358979323846)
#endif

#if !defined(MAXLFOAT)
#define MAXFLOAT (1e38)
#endif


real loglikel_invnorm(real w, real mu, real k)
{
    real dd = w / mu - 1.;
    real sqd = dd * dd / w;
    return .5 * log(k / (2.*M_PI)) - k/2. * sqd - 1.5 * log(w);
}


void dloglikel_invnorm(real w, real mu, real k, real *dll_dmu, real *dll_dw, real *d2ll_dmu2, real *d2ll_dw2, real *d2ll_dmudw)
{
    real dd = w / mu - 1.;
    real mu2 = mu * mu;
    real w2 = w * w;
    *dll_dmu = k * dd / mu2;
    *d2ll_dmu2 = k * (2*mu - 3*w) / mu2 / mu2;
    *dll_dw = .5 * k * (1/w2 - 1/mu2) - 1.5 / w;
    *d2ll_dw2 = 1.5 / w2 - k / w / w2;
    *d2ll_dmudw = k / mu2 / mu;
}


#define min_mu (1e-5)
real regr_mu(int nregr, real *thetap, real *x)
{
    int i;
    real mu;
    for (i=0, mu=0; i<nregr; ++i)
        mu += thetap[i] * x[i];
    return (mu > min_mu ? mu : min_mu);
}


void shift_x(int nregr, int hasTheta0, real last, real *x, real *out)
{
    int i;
    out[0] = x[0];
    for (i=nregr-1; i>hasTheta0; --i)
        out[i] = x[i-1];
    out[hasTheta0] = last;
}


real loglikel_step(int nregr, int hasTheta0, real w, real *thetap, real *x, real k)
{
    real mu = regr_mu(nregr, thetap, x);
    real ll = loglikel_invnorm(mu, w, k);
dbg(("ll(%f %f %f) = %f\n", mu, w, k, ll));
    shift_x(nregr, hasTheta0, w, x, x);
    return ll;
}


int maximize_loglikel_wt(int nregr, real *xt, real *xt1, real w_sum, int hasTheta0, real *thetap, real k, int steps, real *pwt)
{
    int s;
    real wt = *pwt;
    real a1 = thetap[hasTheta0];

/*    for (wt=0; wt<w_sum; wt+=.01) {*/
/*        int i;*/
/*        real mut, mut1, llt, llt1;*/
/*        real wt1 = w_sum - wt;*/
/*        xt1[hasTheta0] = wt;*/
/*        for (i=0, mut=0, mut1=0; i<nregr; ++i) {*/
/*            mut += thetap[i] * xt[i];*/
/*            mut1 += thetap[i] * xt1[i];*/
/*        }*/
/*        llt = loglikel_invnorm(wt, mut, k);*/
/*        llt1 = loglikel_invnorm(wt1, mut1, k);*/
/*        dbg(("%f %f %f\n", wt, llt, llt1));*/
/*    }*/
/*    wt = w_sum / 2;*/

    for (s=0; s<steps; ++s) {
        real dll_dwt, d2ll_dwt2;
        real dllt_dmu, dllt_dw, d2llt_dmu2, d2llt_dw2, d2llt_dmudw;
        real dllt1_dmu, dllt1_dw, d2llt1_dmu2, d2llt1_dw2, d2llt1_dmudw;
        real mut, mut1, dx;
        real wt1 = w_sum - wt;
        xt1[hasTheta0] = wt;
        mut = regr_mu(nregr, thetap, xt);
        mut1 = regr_mu(nregr, thetap, xt1);
/*mexPrintf("%g %g %g\n", loglikel_invnorm(wt, mut, k), loglikel_invnorm(wt, mut, k), loglikel_invnorm(wt, mut, k) + loglikel_invnorm(wt, mut, k));*/
        dloglikel_invnorm(wt, mut, k, &dllt_dmu, &dllt_dw, &d2llt_dmu2, &d2llt_dw2, &d2llt_dmudw);
        dloglikel_invnorm(wt1, mut1, k, &dllt1_dmu, &dllt1_dw, &d2llt1_dmu2, &d2llt1_dw2, &d2llt1_dmudw);
        dll_dwt = dllt_dw - dllt1_dw + a1 * dllt1_dmu;
        d2ll_dwt2 = d2llt_dw2 + d2llt1_dw2 - 2*a1 * d2llt1_dmudw + a1*a1 * d2llt1_dmu2;
        dx = - dll_dwt / d2ll_dwt2;
        if ((fabs(dx) < 1e-6) || (fabs(dll_dwt) < 1e-6))
            break;
        while (dx < -wt || dx > wt1)
            dx = dx / 2;
        wt = wt + dx;
    }
    *pwt = wt;
    return s;
}


int arrythmia_invnorm(int nregr, real *xt, int nR, real *R, int hasTheta0, real *thetap, real k, real *thresholds, int steps, real *tmp, real *loglikel, real *Rt1)
{
    int i, s, action = 0;
    real mut, mut1, mut2, mu_sum, mu_sum2, a1, a2, w_sum, w_sum2, k_sum, k_sum2, wt, wt_opt, loglikel_max, loglikel_best_reg, loglikel_reg, loglikel_reg0, max;
    real eta = 1;
    real *xt1 = tmp; /* tmp should be at least 2 * nregr long */
    real *xt2 = tmp + nregr;
    mut = regr_mu(nregr, thetap, xt);

    shift_x(nregr, hasTheta0, mut, xt, xt1);
    mut1 = regr_mu(nregr, thetap, xt1);

    shift_x(nregr, hasTheta0, mut1, xt1, xt2);
    mut2 = regr_mu(nregr, thetap, xt2);

    a1 = thetap[hasTheta0];
    a2 = thetap[hasTheta0+1];
    mu_sum = mut + mut1;
    mu_sum2 = mut + mut1 + mut2;
    k_sum = k * (mu_sum*mu_sum*mu_sum) / ((1+a1)*(1+a1)*mut*mut*mut + mut1*mut1*mut1);
    k_sum2 = k * (mu_sum2*mu_sum2*mu_sum2) / ((1+a1+a2)*(1+a1+a2)*mut*mut*mut + (1+a1)*(1+a1)*mut1*mut1*mut1 + mut2*mut2*mut2);

    wt = R[1] - R[0];
    w_sum = R[2] - R[0];
    w_sum2 = R[3] - R[0];

dbg(("wt, w_sum, R[2]-R[1], R[3]-R[4], k\n"));
dbg(("%f %f %f %f %f\n", wt, w_sum, R[2]-R[1], R[3]-R[2], k));

    loglikel[0] = loglikel_invnorm(mut, wt, k);
    loglikel[1] = loglikel_invnorm(mu_sum, wt, k_sum);
    loglikel[2] = loglikel_invnorm(mut, w_sum, k);
    loglikel[3] = loglikel_invnorm(mu_sum, w_sum, k_sum);
    loglikel[4] = loglikel_invnorm(mu_sum2, w_sum2, k_sum2);
    loglikel[5] = loglikel_invnorm(mut, w_sum - wt, k);
dbg(("loglikel = %f %f %f %f %f %f\n", loglikel[0], loglikel[1], loglikel[2], loglikel[3], loglikel[4], loglikel[5]));
    
dbg(("ll(%f %f %f) = %f\n", mut, wt, k, loglikel[0]));
    loglikel_reg0 = loglikel[0];
    xt1[hasTheta0] = wt;
    for (i=1; i<nR-2; ++i)
        loglikel_reg0 += loglikel_step(nregr, hasTheta0, R[i+1]-R[i], thetap, xt1, k);
dbg(("loglikel_reg0 %f /(nR-2) %f thresholds[end] %f\n", loglikel_reg0, loglikel_reg0/(nR-2), thresholds[10]));
    loglikel_best_reg = loglikel_reg0;

    if (loglikel[3] - loglikel[0] > thresholds[2]) { /* move event */
        real w_min;
        shift_x(nregr, hasTheta0, mut, xt, xt1);
        wt_opt = w_sum / 2;
        s = maximize_loglikel_wt(nregr, xt, xt1, w_sum, hasTheta0, thetap, k, steps, &wt_opt);
        w_min = wt < w_sum-wt ? wt : w_sum-wt;
        if (s < steps && wt_opt > w_min && wt_opt < w_sum-w_min) {
            loglikel_reg = loglikel_invnorm(mut, wt_opt, k);
dbg(("ll(%f %f %f) = %f\n", mut, wt_opt, k, loglikel_reg));
/*!!            xt1[hasTheta0] = wt_opt; */
            loglikel_reg += loglikel_step(nregr, hasTheta0, w_sum-wt_opt, thetap, xt1, k);
            for (i=2; i<nR-2; ++i)
                loglikel_reg += loglikel_step(nregr, hasTheta0, R[i+1]-R[i], thetap, xt1, k);
dbg(("loglikel_reg move %f /(nR-2) %f\n", loglikel_reg, loglikel_reg/(nR-2)));
            if ((loglikel_reg > loglikel_best_reg) && (loglikel_reg > loglikel_reg0 + (nR-2) * thresholds[7])) {
                loglikel_best_reg = loglikel_reg;
                *Rt1 = R[0] + wt_opt;
                action = 3;
            }
        }
    }
    if (loglikel[4] - loglikel[0] > thresholds[2]+thresholds[3] && loglikel[4] - loglikel[3] > thresholds[3]) { /* move 2 events */
        real w_min, wt_opt_min, wt2_opt = w_sum2/3, w_sum_opt = 2*w_sum2/3;
        int ds;
        wt_opt = wt2_opt;
        shift_x(nregr, hasTheta0, mut, xt, xt1);
        s = 0;
/*mexPrintf("\nmove2 %d   %g %g %g\n", ds, wt, w_sum-wt, w_sum2-w_sum);*/
        do {
            ds = maximize_loglikel_wt(nregr, xt, xt1, w_sum_opt, hasTheta0, thetap, k, steps, &wt_opt);
            shift_x(nregr, hasTheta0, wt_opt, xt1, xt2);
/*mexPrintf("move2a %d   %g %g %g\n", ds, wt_opt, w_sum_opt-wt_opt, w_sum2-w_sum_opt);*/
            ds += maximize_loglikel_wt(nregr, xt1, xt2, w_sum2-wt_opt, hasTheta0, thetap, k, steps, &wt2_opt);
            w_sum_opt = wt_opt + wt2_opt;
            s += ds;
/*mexPrintf("move2b %d   %g %g %g\n", ds, wt_opt, wt2_opt, w_sum2-wt_opt-wt2_opt);*/
        } while (ds > 0 && s < steps);
        w_min = wt < w_sum-wt ? wt : w_sum-wt;
        w_min = w_min < w_sum2-w_sum ? w_min : w_sum2-w_sum;
        wt_opt_min = wt_opt < wt2_opt ? wt_opt : wt2_opt;
        wt_opt_min = wt_opt_min < w_sum2-w_sum_opt ? wt_opt_min : w_sum2-w_sum_opt;
        if (s < steps && wt_opt_min > w_min) {
            loglikel_reg = loglikel_invnorm(mut, wt_opt, k);
dbg(("ll(%f %f %f) = %f\n", mut, wt_opt, k, loglikel_reg));
/*!!            xt1[hasTheta0] = wt_opt; */
/*!!            xt2[hasTheta0] = wt2opt; */
            loglikel_reg += loglikel_step(nregr, hasTheta0, wt2_opt, thetap, xt1, k);
            loglikel_reg += loglikel_step(nregr, hasTheta0, w_sum2-w_sum_opt, thetap, xt2, k);
            for (i=3; i<nR-2; ++i)
                loglikel_reg += loglikel_step(nregr, hasTheta0, R[i+1]-R[i], thetap, xt1, k);
dbg(("loglikel_reg move %f /(nR-2) %f\n", loglikel_reg, loglikel_reg/(nR-2)));
            if (loglikel_reg > loglikel_best_reg + (nR-2) * thresholds[8]) { /* for this to work "move 1 ev" should be first and this second, basically we require a further improvement over "move 1 ev" */
                loglikel_best_reg = loglikel_reg;
                Rt1[0] = R[0] + wt_opt;
                Rt1[1] = R[0] + w_sum_opt;
                action = 4;
            }
        }
    }
    if (loglikel[2] - loglikel[0] > thresholds[1]) { /* remove event */
        shift_x(nregr, hasTheta0, mut, xt, xt1);
        loglikel_reg = loglikel_invnorm(mut, w_sum, k);
dbg(("ll(%f %f %f) = %f\n", mut, w_sum, k, loglikel_reg));
        xt1[hasTheta0] = w_sum;
        for (i=2; i<nR-1; ++i)
            loglikel_reg += loglikel_step(nregr, hasTheta0, R[i+1]-R[i], thetap, xt1, k);
dbg(("loglikel_reg remove %f /(nR-2) %f\n", loglikel_reg, loglikel_reg/(nR-2)));
        if ((loglikel_reg > loglikel_best_reg) && (loglikel_reg > loglikel_reg0 + (nR-2) * thresholds[6])) {
            loglikel_best_reg = loglikel_reg;
            action = 2;
        }
    }
    if (loglikel[1] - loglikel[0] > thresholds[0]) { /* insert event */
        shift_x(nregr, hasTheta0, mut, xt, xt1);
        wt_opt = wt / 2.;
        s = maximize_loglikel_wt(nregr, xt, xt1, wt, hasTheta0, thetap, k, steps, &wt_opt);
        if (s < steps) {
            loglikel_reg = loglikel_invnorm(mut, wt_opt, k);
dbg(("ll(%f %f %f) = %f\n", mut, wt_opt, k, loglikel_reg));
            xt1[hasTheta0] = wt_opt;
            loglikel_reg += loglikel_step(nregr, hasTheta0, wt-wt_opt, thetap, xt1, k);
            if (nR > 4)
                loglikel_reg += loglikel_step(nregr, hasTheta0, w_sum-wt, thetap, xt1, k);
            for (i=2; i<nR-3; ++i)
                loglikel_reg += loglikel_step(nregr, hasTheta0, R[i+1]-R[i], thetap, xt1, k);
dbg(("loglikel_reg insert %f /(nR-2) %f\n", loglikel_reg, loglikel_reg/(nR-2)));
/*dbg(("xt[i], xt1[i], thetap[i]\n"));*/
/*for (i=0; i<nregr; ++i)*/
/*dbg(("%f %f %f\n", xt[i], xt1[i], thetap[i]));*/
/*dbg(("mut, mut1\n"));*/
/*dbg(("%f %f\n", mut, mut1));*/
/*dbg(("loglikel[0], loglikel[1], loglikel[2], loglikel_best_reg, loglikel_reg\n"));*/
/*dbg(("%f %f %f %f %f\n", loglikel[0], loglikel[1], loglikel[2], loglikel_best_reg, loglikel_reg));*/
            if ((loglikel_reg > loglikel_best_reg) && (loglikel_reg > loglikel_reg0 + (nR-2) * thresholds[5])) {
                loglikel_best_reg = loglikel_reg;
                *Rt1 = R[0] + wt_opt;
                action = 1;
dbg(("insert action\n"));
            }
        }
    }
    loglikel_max = (loglikel[0] > loglikel[1]) ? loglikel[0] : loglikel[1];
    loglikel_max = (loglikel_max > loglikel[2]) ? loglikel_max : loglikel[2];
    loglikel_max = (loglikel_max > loglikel[3]) ? loglikel_max : loglikel[3];
    if (loglikel[5] - loglikel_max > thresholds[4]) { /* ignore event */
        shift_x(nregr, hasTheta0, mut, xt, xt1);
        loglikel_reg = loglikel_invnorm(mut, w_sum-wt, k);
dbg(("ll(%f %f %f) = %f\n", mut, w_sum-wt, k, loglikel_reg));
        xt1[hasTheta0] = w_sum - wt;
        for (i=2; i<nR-1; ++i)
            loglikel_reg += loglikel_step(nregr, hasTheta0, R[i+1]-R[i], thetap, xt1, k);
dbg(("loglikel_reg ignore %f /(nR-2) %f\n", loglikel_reg, loglikel_reg/(nR-2)));
        if ((loglikel_reg > loglikel_best_reg) && (loglikel_reg > loglikel_reg0 + (nR-2) * thresholds[9])) {
            loglikel_best_reg = loglikel_reg;
            action += 16;
        }
    }

    if (loglikel_best_reg / (nR-2) > thresholds[10]) /* cannot trust the regressors enough to improve the current event */
        return 0;

if (action) dbg(("ACTION %d\n", action));
    return action;
}




char syntax[] = "Syntax:\n[action, Rt1, loglikel] = arrythmia_invnorm_mex(xt, R(j-1:j+P+1), hasTheta0, thetap, k, thresholds);";

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i, nregr, nR, hasTheta0, action, steps=1000;
    double *xt, *R, *thetap, k, *tmp;
    double thresholds[] = {3., 3., 3., 3., 3., 10., 10., 10., 10., 10., MAXFLOAT};
    mxArray *m_loglikel, *m_Rt;

    if (((nrhs < 5) || (nrhs > 7)) || (nlhs < 1) || (nlhs > 3))
        mexErrMsgIdAndTxt("arrythmia_invnorm_mex:syntax", syntax);
    for (i=0; i<nrhs; ++i)
        if (!mxIsDouble(prhs[i]))
            mexErrMsgIdAndTxt("arrythmia_invnorm_mex:double", "All inputs should be of class double.");

    nregr = mxGetNumberOfElements(prhs[0]);
    nR = mxGetNumberOfElements(prhs[1]);
    if (nR < 4 || mxGetNumberOfElements(prhs[3]) != nregr)
        mexErrMsgIdAndTxt("arrythmia_invnorm_mex:nregr", "The argument 'xn' should be N elements, 'R' should be at least 4 elements, 'thetap' should be N elements.");

    if (nrhs >= 6) {
        double *athresholds = mxGetPr(prhs[5]);
        int nthresholds = mxGetNumberOfElements(prhs[5]);
        if (nthresholds > 10)
            mexErrMsgIdAndTxt("arrythmia_invnorm_mex:nregr", "The argument 'thresholds' should be [DL_add, DR_add; DL_remove, DR_remove; DL_move, DR_move; DL_move2, DR_move2; DL_ignore, DR_ignore].");
        for (i=0; i<nthresholds; ++i)
            thresholds[i] = athresholds[i];
        if (nrhs >= 7)
            thresholds[10] = mxGetScalar(prhs[6]);
    }

    xt = mxGetPr(prhs[0]);
    R = mxGetPr(prhs[1]);
    hasTheta0 = mxGetScalar(prhs[2]);
    thetap = mxGetPr(prhs[3]);
    k = mxGetScalar(prhs[4]);
    m_loglikel = mxCreateDoubleMatrix(6, 1, mxREAL);
    m_Rt = mxCreateDoubleMatrix(2, 1, mxREAL);
    tmp = mxMalloc(2*nregr*sizeof(double));

    action = arrythmia_invnorm(nregr, xt, nR, R, hasTheta0, thetap, k, thresholds, steps, tmp, mxGetPr(m_loglikel), mxGetPr(m_Rt));

    mxFree(tmp);
    plhs[0] = mxCreateDoubleScalar((double)action);
    if (nlhs >= 2) {
        if (action == 1 || action == 3 || action == 17 || action == 19)
            mxSetM(m_Rt, 1);
        else if(action != 4 && action != 20) 
            mxSetM(m_Rt, 0);
        plhs[1] = m_Rt;
    } else
        mxDestroyArray(m_Rt);
    if (nlhs >= 3)
        plhs[2] = m_loglikel;
    else
        mxDestroyArray(m_loglikel);
}
