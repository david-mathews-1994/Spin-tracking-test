/* DOP853

This is my adaptation of the DOP853 code from J .Colinge (COLINGE@DIVSUN.UNIGE.CH) which in turn is an adaptation of E. Hairer's
original FORTRAN code. The main changes are the removal of pointer arrays for passing between functions
and using std::vector instead. The other main change is that the intergrator is now a class instead of a
function which relied on delcaring a bunch of static variables. I did this because it should make it cleaner
to use and expand in the future, and becuase the static variables will create huge problems when trying to
parallelize the spin simulation this is intended for. In addition to modernizing the C code I took the liberty
to remove certain parts that we don't really need for our purposes. Most of the error checking is still present
but more could probably be added.

Author: Matt Morano
Date: 1/4/2023

Credit for the original implementation goes to:

    Authors : E. Hairer & G. Wanner
      Universite de Geneve, dept. de Mathematiques
      CH-1211 GENEVE 4, SWITZERLAND
      e-mail : HAIRER@DIVSUN.UNIGE.CH, WANNER@DIVSUN.UNIGE.CH

The code is described in : E. Hairer, S.P. Norsett and G. Wanner, Solving
ordinary differential equations I, nonstiff problems, 2nd edition,
Springer Series in Computational Mathematics, Springer-Verlag (1993).

*/

#pragma once
#include "../include/coeff.h"
#include "../include/options.h"
#include <openacc.h>
#include <functional>
#include <iostream>

typedef void (*OBS)(long nr, double xold, double x, double* y, int* irtrn);

class DOP853
{
public:

    DOP853() {}

    DOP853(std::function<void(const double x, const double* y, double* f)> fcn, double* y0, OBS obs, options OPT)
    : n(3), y(y0), fcn(fcn), rtol(OPT.rtol), atol(OPT.atol), obs(obs), iout(OPT.iout), uround(OPT.uround), safe(OPT.safe), 
      fac1(OPT.fac1), fac2(OPT.fac2), beta(OPT.beta), hmax(OPT.hmax), h(OPT.h), nmax(OPT.nmax)
      //, yy1(), k1(), k2(), k3(), k4(), k5(), k6(), k7(), k8(), 
      // k9(), k10(), rcont1(), rcont2(), rcont3(), rcont4(), rcont5(), rcont6(), rcont7(), rcont8() 
      {}

    int integrate(double x0, double xf);

    double sign(double a, double b);

    double min_d(double a, double b);

    double max_d(double a, double b);

    double hinit(std::function<void(const double x, const double* y, double* f)> fcn, double x, double* y,
      double posneg, double* f0, double* f1, double* yy1, int iord, double hmax, double atoler, double rtoler);

private:
    unsigned int n;
    std::function<void(const double x, const double* y, double* f)> fcn;
    double* y;
    double x;
    double xend;
    double rtol;
    double atol;
    OBS obs;
    int iout;
    double uround;
    double safe;
    double fac1;
    double fac2;
    double beta;
    double hmax;
    double h;
    long nmax;
    double facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
    double atoli, rtoli, hlamb, err, sk, hnew, yd0, ydiff, bspl;
    double stnum, stden, sqr, err2, erri, deno;
    int iasti, iord, irtrn, reject, last, nonsti;
    unsigned int i, j;
    unsigned int nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
    unsigned int *indir;
    double yy1[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], k8[3], k9[3], k10[3];
    double rcont1[3], rcont2[3], rcont3[3], rcont4[3], rcont5[3], rcont6[3], rcont7[3], rcont8[3];
    double hout, xold, xout;
    int arret, idid;
};