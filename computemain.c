// -------------------------------------------------------------------
// Program last modified August 11, 2025. 
// Copyright (c) 2024-2025 Terrence P. Murphy
// MIT License -- see hardyz.c for details.
// -------------------------------------------------------------------
#include "hardyz.h"

mpfr_t	myPi, my2Pi, myLog2;

// *******************************************************************
// We call this function before using any MPFR functions.  We set the
// default MPFR precision and create global variables holding the
// values of Pi 2Pi and Log(2), 
// *******************************************************************
int InitMPFR(int iFloatBits, int DebugFlagsSet)
{
// -------------------------------------------------------------------
// Set default precision for MPFR
// -------------------------------------------------------------------
mpfr_set_default_prec (iFloatBits);

// -------------------------------------------------------------------
// Initialize the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_inits2 (iFloatBits, myPi, my2Pi, myLog2, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Initialize the MPFR variables that will hold the coefficients.
// -------------------------------------------------------------------
InitCoeffMPFR(iFloatBits);	
BuildCoefficientsMPFR(DebugMode(DebugFlagsSet, PRINT_COEFF));

// -------------------------------------------------------------------
// Set the value of the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_const_pi (myPi, MPFR_RNDN); 
mpfr_mul_2ui (my2Pi, myPi, 1, MPFR_RNDN); /* 2Pi */
mpfr_const_log2 (myLog2, MPFR_RNDN);
return(0);
}

// *******************************************************************
// We call this function after we are done using all MPFR functions.  
// We clear the remaining memory and any cache used by MPFR. 
// *******************************************************************
int CloseMPFR(void)
{
// -------------------------------------------------------------------
// Free the space used by the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_clears (myPi, my2Pi, myLog2, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Free the MPFR variables that held the coefficients.
// -------------------------------------------------------------------
CloseCoeffMPFR();

// -------------------------------------------------------------------
// Clear the cache used by MPFR.
// -------------------------------------------------------------------
mpfr_free_cache ();
return(0);
}

// *******************************************************************
// We compute the Hardy Z values here.  The loop is over the count of
// different 't' values to compute (based on hz.iCount). Inside the 
// loop, we call ComputeSingleHardyZ. We then printf the result,
// and then increase 't' by hz.dIncr and repeat hz.iCount times.
// *******************************************************************
int ComputeAllHardyZ(struct HZ hz)
{
mpfr_t				t, Incr;
struct computeHZ 	comphz[MAX_THREADS];
int					i;

// -------------------------------------------------------------------
// Initiate use of the MPFR system.  
// -------------------------------------------------------------------
InitMPFR(hz.iFloatBits, hz.iDebug);

// -------------------------------------------------------------------
// Initialize several MPFR variables.  The 't' value is obtained
// from the hz.tBuf string, and the amount to increment 't' is obtained
// from the hz.incrBuf string.  
// We use the mpfr_set_str function because the MPFR FAQ says that the
// other variations such as mpfr_set_ld are much less accurate.
// -------------------------------------------------------------------
mpfr_inits2 (hz.iFloatBits, t, Incr, (mpfr_ptr) 0);

for(i=0; i< hz.iThreads; i++) { 
	mpfr_inits2 (hz.iFloatBits, comphz[i].t, comphz[i].Result, (mpfr_ptr) 0);
	comphz[i].ptrResult = &comphz[i].Result;
	comphz[i].iFloatBits = hz.iFloatBits;		
	comphz[i].DebugFlagsSet = hz.iDebug;		
	}

mpfr_set_str (t, hz.tBuf, 10, MPFR_RNDN);
mpfr_set_str (Incr, hz.incrBuf, 10, MPFR_RNDN);

// -------------------------------------------------------------------
// We loop hz.iCount times, processing 't' in the first (i = 1) loop.
// If hz.iCount is greater than one, we increment 't' by the Incr
// amount and loop again.  
//
// In the loop, we call ComputeSingleHardyZ and print the result.
// -------------------------------------------------------------------
int			j, m, iRemaining;
pthread_t	thread_id[MAX_THREADS];

i = 0;
while(i < hz.iCount)
	{
	iRemaining = hz.iCount - i;
	m = hz.iThreads > iRemaining ? iRemaining : hz.iThreads;
	
	for (j = 0; j < m; j++)
		{
		// ---------------------------------------------------------------
		// Update comphz values for the current 't'
		// ---------------------------------------------------------------	
		mpfr_set (comphz[j].t, t, MPFR_RNDN);
		mpfr_set_ui (comphz[j].Result, 0, MPFR_RNDN);
		
		if(hz.iThreads > 1) {
			pthread_create( &thread_id[j], NULL, ComputeSingleHardyZThreaded, &comphz[j] );
			}
		else {
			ComputeSingleHardyZ(&comphz[0]);
			}
		mpfr_add (t, t, Incr, MPFR_RNDN);
		}
	for (j = 0; j < m; j++)
		{
		if(hz.iThreads > 1) {
			pthread_join( thread_id[j], NULL); 
			}
		}		
	for (j = 0; j < m; j++)
		{
		hz.bVerbose ? mpfr_printf("For t = %.*Rf, Z(t) = %.*Rf \n", 
					hz.iOutputDPt, comphz[j].t, hz.iOutputDPz, comphz[j].Result)
			: mpfr_printf("%.*Rf, %.*Rf \n", hz.iOutputDPt, comphz[j].t, 
					hz.iOutputDPz, comphz[j].Result);
		i++;	// we need 'i' here so it is incremented 'm' times 
		}			
	}	// end of outer for loop

// -------------------------------------------------------------------
// We are done with the MPFR system.  Clear our local MPFR variables
// and then close down the MPFR system.
// -------------------------------------------------------------------
mpfr_clears (t, Incr, (mpfr_ptr) 0);

for(i=0; i< hz.iThreads; i++) {
	mpfr_clears (comphz[i].t, comphz[i].Result, (mpfr_ptr) 0);
	}

CloseMPFR();
return(1);
}

// *******************************************************************
// For correct typing on the pthread_cerate call, we need this
// pass-through function.
// *******************************************************************
void * ComputeSingleHardyZThreaded(void * comphz)
{
ComputeSingleHardyZ((struct computeHZ *) comphz);
return(NULL);
}

// *******************************************************************
// We compute a single Hardy Z values here.  
// *******************************************************************
int ComputeSingleHardyZ(struct computeHZ * comphz)
{
mpfr_t			tOver2Pi, T, N, P, Main, Remainder;
unsigned int	uiN;

mpfr_inits2 (comphz->iFloatBits, 
				tOver2Pi, T, N, P, Main, Remainder, (mpfr_ptr) 0);

// ---------------------------------------------------------------
// Compute N and P for the given 't'. 
// NOTE: Because N is (currently) an unsigned int (we assume 32-bit),
// then 0 <= N <= 4,294,967,295.  
// Thus, t cannot exceed about 1.15 * 10^{20} in the calcs below.
// (We assume here that the value of 't' was previously checked).
// ---------------------------------------------------------------	
mpfr_div (tOver2Pi, comphz->t, my2Pi, MPFR_RNDN);
mpfr_sqrt (T, tOver2Pi, MPFR_RNDN);
mpfr_modf (N, P, T, MPFR_RNDN);
uiN = (unsigned int) mpfr_get_ui (N, MPFR_RNDN);
	
// ---------------------------------------------------------------
// Compute the remainder term.
// ---------------------------------------------------------------		
ComputeRemainderMPFR(&Remainder, tOver2Pi, uiN, P, comphz->iFloatBits);

if(DebugMode(comphz->DebugFlagsSet, PRINT_REMAINDER)) {
	mpfr_printf("Remainder R(4): %.50Rf \n", Remainder);
	}
	
// ---------------------------------------------------------------
// Now compute the Main term and add to Remainder to get HardyZ.
// ---------------------------------------------------------------	
ComputeMain(&Main, comphz->t, uiN, comphz->iFloatBits, comphz->DebugFlagsSet);
mpfr_add (*comphz->ptrResult, Main, Remainder, MPFR_RNDN);

// -------------------------------------------------------------------
// Clear our local MPFR variables.
// -------------------------------------------------------------------
mpfr_clears (tOver2Pi, T, N, P, Main, Remainder, (mpfr_ptr) 0);
return(1);
}

// *******************************************************************
// We compute the main term of the Riemann-Siegel formula.  
//
// First Step: compute theta(t).  
//
// Second Step: a loop (from n = 1 to n = N), where we calculate and add 
// together the results from each n:
//    
//		sqrt(1/n) * cos[theta(t) - t log n]
//
// Third Step: return 2 times the sum of those N terms.
// *******************************************************************
int ComputeMain(mpfr_t *Result, mpfr_t t, unsigned int N, int iFloatBits, int DebugFlagsSet)
{
mpfr_t	Theta;

// -------------------------------------------------------------------
// Step 0: Check the case N < 1 (nothing to do so return 0 in Result).
// -------------------------------------------------------------------
if(N < 1)
	{
	mpfr_set_ui (*Result, 0, MPFR_RNDN);
	return(1);
	}

// -------------------------------------------------------------------
// Step 1: Compute Theta.
// -------------------------------------------------------------------
mpfr_inits2 (iFloatBits, Theta, (mpfr_ptr) 0);
ComputeTheta(&Theta, t, iFloatBits);
//mpfr_printf("Theta = %.40Rf\n", Theta);

// -------------------------------------------------------------------
// Step 2: loop n = 1 to N.  We will handle the n = 1 and
// n = 2 cases separately before entering the loop.
// ------------------------------------------------------------------
mpfr_t	Temp1, Temp2;
mpfr_t	Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne;

mpfr_inits2 (iFloatBits, Temp1, Temp2, 
	Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// For the n = 1 term, we set the initial value of Main to cos(theta).
// -------------------------------------------------------------------
mpfr_cos (Main, Theta, MPFR_RNDN); 

// -------------------------------------------------------------------
// Next, the n = 2 term.  Why outside the loop? Depending on which 
// algorithm we use, we MAY need the value of the n = 2 
// CosArg = [theta - t log 2] and of Log 2 --> LogMinusOne BEFORE we 
// enter the main loop.  
// For the n = 2 term, we calculate sqrt(1/2) * cos(Theta - t log 2)
// and add that value to Main.
// -------------------------------------------------------------------
if(N >= 2)  
	{
	mpfr_set_ui (Temp1, 2, MPFR_RNDN);
	mpfr_rec_sqrt (RecipSqrtn, Temp1, MPFR_RNDN);	// sqrt(1/2)
	mpfr_mul (Temp1, t, myLog2, MPFR_RNDN);			// t * log 2
	mpfr_sub (CosArg, Theta, Temp1, MPFR_RNDN);		// theta - [t * log 2]
	mpfr_fmod (CosArg, CosArg, my2Pi, MPFR_RNDN);	// CosArg % 2Pi = OPTIONAL
	mpfr_cos (Temp1, CosArg, MPFR_RNDN);			// Temp1 = cos(CosArg)
	mpfr_mul (FullTerm, RecipSqrtn, Temp1, MPFR_RNDN);
	mpfr_add (Main, Main, FullTerm, MPFR_RNDN);
	mpfr_set (LognMinusOne, myLog2, MPFR_RNDN); 	// for Version 2, below
	}
	
// -------------------------------------------------------------------
// Now process the n = 3 through n = N terms
// -------------------------------------------------------------------
for (unsigned int n = 3; n <= N; ++n) { // need C99 ro define 'n' here 
	// ---------------------------------------------------------------
	// We need an mpfr_t version of n, so save in Temp1
	// ---------------------------------------------------------------	
	mpfr_set_ui (Temp1, n, MPFR_RNDN);
	
	// ---------------------------------------------------------------
	// First, compute the square root of 1/n
	// ---------------------------------------------------------------	
	mpfr_rec_sqrt (RecipSqrtn, Temp1, MPFR_RNDN);	

	if(DebugMode(DebugFlagsSet, COS_ARG_NOT_SAVED)) { // do NOT do Version 2 (saving CosArg, etc.)
		//#####################################################################################	
		// ---------------------------------------------------------------
		// Second, compute the argument to the cosine term. (Version 1)
		// That is, CosArg = [theta(t) - t log n].  
		// Then (further below) compute cos(CosArg).
		// ---------------------------------------------------------------	
		mpfr_log (Temp2, Temp1, MPFR_RNDN);			// log n
		mpfr_mul (Temp2, t, Temp2, MPFR_RNDN); 		// t * log n
		mpfr_sub (CosArg, Theta, Temp2, MPFR_RNDN); // theta(t) - [t * log n]	
		//#####################################################################################
		}
	else { // this is the default case:
		//#####################################################################################	
		// ---------------------------------------------------------------
		// Second, compute the argument to the cosine term. (Version 2)
		// We have saved CosArg(of n-1) and log(n-1), allowing:  
		// CosArg(of n) = CosArg(of n-1) + (t * [log(n-1) - log(n)])
		// ---------------------------------------------------------------	
		mpfr_log (Temp1, Temp1, MPFR_RNDN); 				// this sets Temp1 to log(n)
		mpfr_sub (Temp2, LognMinusOne, Temp1, MPFR_RNDN); 	// Temp2 = log(n-1) - log(n)
		mpfr_mul (Temp2, t, Temp2, MPFR_RNDN); 				// Temp2 = t * [log(n-1) - log(n)]
		mpfr_add (CosArg, CosArg, Temp2, MPFR_RNDN); 		// CosArg = Temp2 + CosArg(of n-1)
		mpfr_set (LognMinusOne, Temp1, MPFR_RNDN); 			// update for the next 'n' - use "swap"?
		//#####################################################################################
		}

	if(!DebugMode(DebugFlagsSet, COS_ARG_2PI)) { // skip this if indicated debug flag is set
		//#####################################################################################	
		//----------------------------------------------------------------
		// divide CosArg by 2 pi and keep the remainder
		// Or, should we just let mpfr_cos do the same thing?
		//----------------------------------------------------------------
		mpfr_fmod (CosArg, CosArg, my2Pi, MPFR_RNDN);
		//#####################################################################################	
		}

	//----------------------------------------------------------------
	// We are now ready to compute the cosine value = CosCalc.
	//----------------------------------------------------------------
	mpfr_cos (CosCalc, CosArg, MPFR_RNDN);
	//----------------------------------------------------------------
	// For the full term, multiply CosCalc by RecipSqrtn, then
	// add to Main.
	//----------------------------------------------------------------
	mpfr_mul (FullTerm, RecipSqrtn, CosCalc, MPFR_RNDN);
	mpfr_add (Main, Main, FullTerm, MPFR_RNDN);
	} // end of for loop

// -------------------------------------------------------------------
// We have calculated Main.  Now, multiply by 2 and return result.
// -------------------------------------------------------------------
mpfr_mul_2ui (*Result, Main, 1, MPFR_RNDN);

// -------------------------------------------------------------------
// Free the space used by the local mpfr (constant) variables
// -------------------------------------------------------------------	
mpfr_clears (Theta, Temp1, Temp2,  
  Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne, (mpfr_ptr) 0);

return(1);
}
