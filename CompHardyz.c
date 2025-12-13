// -------------------------------------------------------------------
// Program last modified November 30, 2025. 
// Copyright (c) 2024-2025 Terrence P. Murphy
// MIT License -- see hardyz.c for details.
// -------------------------------------------------------------------

#include <quadmath.h>
#define MPFR_WANT_FLOAT128 1

#include <time.h>
#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <pthread.h>
#include <mpfr.h>

#include "hgt.h"
#include "hardyz.h"

struct HZ_RPT 	hzReport;


// *******************************************************************
// We compute the Hardy Z values here.  The loop is over the count of
// different 't' values to compute (based on hz.iCount). Inside the 
// loop, we call ComputeSingleHardyZ. We then printf the result,
// and then increase 't' by hz.dIncr and repeat hz.iCount times.
// *******************************************************************
int ComputeAllHardyZ(struct HZ hz)
{
mpfr_t				t, Incr;

// -------------------------------------------------------------------
// Initialize the MPFR system, and the variables that will hold the 
// Riemann-Siegel coefficients.
// -------------------------------------------------------------------
InitMPFR(hz.DefaultBits, hz.Threads, hz.DebugFlags, true);

// -------------------------------------------------------------------
// Initialize the 't' and 'Incr' MPFR variables.  The 't' value 
// is obtained from the hz.tBuf string, and the amount to increment 
// 't' is obtained from the hz.incrBuf string.  
// We use the mpfr_set_str function because the MPFR FAQ says the
// other variations such as mpfr_set_ld are much less accurate.
// -------------------------------------------------------------------
mpfr_inits2 (hz.DefaultBits, t, Incr, (mpfr_ptr) 0);
mpfr_set_str (t, hz.tBuf, 10, MPFR_RNDN);
mpfr_set_str (Incr, hz.incrBuf, 10, MPFR_RNDN);

hzReport.OutputDPt = hz.OutputDPt;
hzReport.OutputDPz = hz.OutputDPz;
hzReport.Verbose   = hz.Verbose;

HardyZWithCount(t, Incr, hz.Count, 0, HardyZCallback);

CloseCoeffMPFR();
CloseMPFR();
return(1);
}

// We avoid the unused parameter warning with these pragma 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

// *******************************************************************
// This is the callback function called by a libhgt.a function
// computing Hardy Z values. This function is called, in order, for each
// of the values computed.  The passed 'i' is the index (zero-based)
// of the current incremented 't' being processed. 
// *******************************************************************
int HardyZCallback(mpfr_t t, mpfr_t HardyZ, int i, int CallerID) 
{
if(hzReport.Verbose == true) { // true = Verbose report
	mpfr_printf("Item = %d, t = %.*Rf, Z(t) = %.*Rf \n", 
		++i, hzReport.OutputDPt, t, hzReport.OutputDPz, HardyZ);	
	}
else {
	mpfr_printf("%.*Rf, %.*Rf \n", 
		hzReport.OutputDPt, t, hzReport.OutputDPz, HardyZ);	
	}
return(1);
}
#pragma GCC diagnostic pop
