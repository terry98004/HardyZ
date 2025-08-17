// -------------------------------------------------------------------
// Program last modified August 9, 2025. 
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
#include <pthread.h>
#include <mpfr.h>

#define		MY_DEFAULT_PRECISION	256
#define		MAX_T_POWER3_STRING		"1.1e14"
#define		MAX_T					1.15e20		// (otherwise N > (4-byte) unsigned int)
												// would need (8-byte) unsigned long long int
#define		MAX_INCREMENT			1.15e10		// (so incremented 't' is not too large)

#define		T_BUF_SIZE				100
#define		INCR_BUF_SIZE			80
#define		NUM_POWERS_P			40
#define		NUM_POWERS_P_GABCKE		88
#define		NUM_Cj_TERMS			5
#define		COEFF_PER_Cj			44
#define		DECIMAL_PLACES_GABCKE	50
#define		MAX_THREADS				8

	// debug flags 
#define	DEBUG_RESERVED1				2
#define	PRINT_COEFF					3
#define COS_ARG_2PI					5
#define	COS_ARG_NOT_SAVED			7
#define	PRINT_REMAINDER				11
#define	USE_THREADS					13
#define	DEBUG_MAX_VALUE				30030


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
	int			iThreads;	// Number of threads to use 
}; 

struct computeHZ {
	mpfr_t		t; 				// 't' value to compute
	mpfr_t		Result; 		// To hold mpfr computed value
	mpfr_t		*ptrResult;		// To hold pointer to mpfr computed value
	int			iFloatBits;		// Bits for MPFR floating point 
	int			DebugFlagsSet;	// Based on user input
}; 

int			strCheckAndCount(const char *string, double dMax);
int 		GetFloatBitsMPFR(const char *str);
int 		ComputeAllHardyZ(struct HZ hz);
int 		ComputeSingleHardyZ(struct computeHZ * comphz);
void * 		ComputeSingleHardyZThreaded(void * comphz);
int 		ComputeMain(mpfr_t *Result, mpfr_t t, unsigned int N, int iFloatBits, int DebugFlagsSet);
int 		ComputeTheta(mpfr_t *Theta, mpfr_t t, int iFloatBits);
int 		InitMPFR(int iFloatBits, int DebugFlagsSet);
int 		CloseMPFR(void);
bool 		DebugMode(int FlagsSet, int FlagToTest);
int 		ThetaAxB(mpfr_t *Result, mpfr_t Big, mpfr_t Small, int iFloatBits);
void 		Show128(__float128 x);

int 		BuildCoefficientsMPFR(bool Debug);
int 		CoeffStrToMPFR(mpfr_t *Result, const char *strCoeff);
int 		InitCoeffMPFR(int iFloatBits);
int 		CloseCoeffMPFR(void);
int 		ComputeRemainderMPFR(mpfr_t *Result, mpfr_t tOver2Pi, unsigned int N, mpfr_t P, int iFloatBits);