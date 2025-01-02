// -------------------------------------------------------------------
// Program last modified January 1, 2025. 
// Copyright (c) 2024-2025 Terrence P. Murphy
// MIT License -- see hardyZ.c for details.
// -------------------------------------------------------------------

#define _USE_MATH_DEFINES	// This includes the define M_PI for pi -- needed??

#include <quadmath.h>
#define MPFR_WANT_FLOAT128 1

#include <time.h>
#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <unistd.h>			// for command line argument processing
#include <mpfr.h>

#define		MY_DEFAULT_PRECISION	256
#define		MAX_T_POWER3_STRING		"1.1e12"
#define		MAX_T_STRING			"1.15e20"	// (otherwise N > (4-byte) unsigned int)
												// would need (8-byte) unsigned long long int
#define		T_BUF_SIZE				100
#define		INCR_BUF_SIZE			80
#define		NUM_POWERS_P			40

#define		NUM_POWERS_P_GABCKE		88
#define		NUM_Cj_TERMS			5
#define		COEFF_PER_Cj			44
#define		DECIMAL_PLACES_GABCKE	50

	// debug flags 
#define	COEFF_128					2
#define	PRINT_COEFF					3
#define COS_ARG_2PI					5
#define	COS_ARG_NOT_SAVED			7
#define	PRINT_REMAINDER				11
#define	DEBUG_MAX_VALUE				2310


struct HZ {
	char		tBuf[T_BUF_SIZE];		// holds entered '-t' value
	char		incrBuf[INCR_BUF_SIZE];	// holds entered '-d' value
	int 		iCount;     // number of t values to check
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
int 		ComputeTheta(mpfr_t *Result, mpfr_t t, struct HZ hz);
int 		InitMPFR(struct HZ hz);
int 		CloseMPFR(struct HZ hz);
int 		ComputeRemainder128(mpfr_t *Result, mpfr_t P, mpfr_t tFraction, unsigned int N, struct HZ hz);
bool 		DebugMode(struct HZ hz, int DebugNum);
int 		ThetaAxB(mpfr_t *Result, mpfr_t Big, mpfr_t Small, struct HZ hz);
void 		Show128(__float128 x);

int 		BuildCoefficientsMPFR(struct HZ hz);
int 		CoeffStrToMPFR(mpfr_t *Result, const char *strCoeff);
int 		InitCoeffMPFR(int FloatBits);
int 		CloseCoeffMPFR(void);
int 		ComputeRemainder256(mpfr_t *Result, mpfr_t P, mpfr_t tFraction, unsigned int N, struct HZ hz);