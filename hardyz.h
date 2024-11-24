// -------------------------------------------------------------------
// Program last modified November 23, 2024. 
// Copyright (c) 2024 Terrence P. Murphy
// MIT License -- see hardyZ.c for details.
// -------------------------------------------------------------------

#define _USE_MATH_DEFINES	// This includes the define M_PI for pi -- needed??
//              3.14159265358979323846264338327950288419716939937510
#define	M_PI_X	3.141592653589793238462643383279502884L // for 80-bit long double
#define	M_2PI_X	6.283185307179586476925286766559005768L // for 80-bit long double

#include <quadmath.h>
#define MPFR_WANT_FLOAT128 1

#include <time.h>
#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>			// for command line argument processing
#include <mpfr.h>

#define		MY_DEFAULT_PRECISION	256
#define		MAXT_POWER3_TERM		31000000000L	// long double  = 31,000,000,000
#define		T_BUF_SIZE				100

struct HZ {
	char		tBuf[T_BUF_SIZE];
	long double	t;
	int 		iCount;     // number of t values to check
	double		dIncr;		// amount to increment t
	bool		bVerbose;	// verbose report? true or false
	bool		bSeconds;	// report compute second taken T/F
	int			iDebug;		// used for debugging
	int			iActualDPt;	// digits after '.' in t
	int			iWholeDt;	// digits before '.' in t
	int			iActualDPi;	// digits after '.' in dIncr
	int			iWholeDi;	// digits before '.' in dIncr
	int			iWholeD;	// max(iWholeDt, iWholeDi)
	int			iOutputDPt;	// digits after '.' in reported t value 
	int			iOutputDPz;	// digits after '.' in reported HardyZ value 
	int			iFloatBits;	// Bits for MPFR floating point 
}; 

int			strCheckAndCount(const char *string);
int 		GetFloatBits(const char *str);
int 		ComputeHardyZ(struct HZ hz);
int			ComputeMain(mpfr_t *Result, mpfr_t t, unsigned int N, struct HZ hz);
int 		InitMPFR(struct HZ hz);
int 		CloseMPFR(void);
long double ComputeRemainder(unsigned int uiN, long double tFraction, long double dP);
int 		ComputeRemainder128(mpfr_t *Result, mpfr_t P, mpfr_t tFraction, unsigned int N, struct HZ hz);
void 		Show128(__float128 x);