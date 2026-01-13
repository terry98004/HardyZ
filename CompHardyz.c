// -------------------------------------------------------------------
// Program last modified January 13, 2026. 
// Copyright (c) 2024-2026 Terrence P. Murphy
// MIT License -- see hardyz.c for details.
// -------------------------------------------------------------------

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
// We have gathered the user's command line input and placed most of 
// those values (plus various default values the user did not override)
// into the passed struct HZ. We first do some further light processing 
// of those values.  Then, we call the libHGT library function 
// "HardyZWithCount" to compute the "Count" number of Hardy Z values.
// As part of the "HardyZWithCount" processing, after each Hardy Z value 
// is calculated, there is a call to the "HardyZCallback" function 
// (below) to printf a report of that Hardy Z value.
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
