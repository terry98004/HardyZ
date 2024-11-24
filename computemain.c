// -------------------------------------------------------------------
// Program last modified November 23, 2024. 
// Copyright (c) 2024 Terrence P. Murphy
// MIT License -- see hardyz.c for details.
// -------------------------------------------------------------------
#include "hardyz.h"

mpfr_t	myPi, my2Pi, myLog2;

// -------------------------------------------------------------------
// We call this function before using any MPFR functions.  We set the
// default MPFR precision and create global variables holding the
// values of Pi 2Pi and Log(2), 
// -------------------------------------------------------------------
int InitMPFR(struct HZ hz)
{
// -------------------------------------------------------------------
// Set default precision for MPFR
// -------------------------------------------------------------------
mpfr_set_default_prec (hz.iFloatBits);

// -------------------------------------------------------------------
// Initialize the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_inits2 (hz.iFloatBits, myPi, my2Pi, myLog2, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Set the value of the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_const_pi (myPi, MPFR_RNDN); 
mpfr_mul_2ui (my2Pi, myPi, 1, MPFR_RNDN); /* 2Pi */
mpfr_const_log2 (myLog2, MPFR_RNDN);
return(0);
}


// -------------------------------------------------------------------
// We call this function after we are done using all MPFR functions.  
// We clear the remaining memory and any cache used by MPFR. 
// -------------------------------------------------------------------
int CloseMPFR(void)
{
// -------------------------------------------------------------------
// Free the space used by the global mpfr (constant) variables
// -------------------------------------------------------------------
mpfr_clears (myPi, my2Pi, myLog2, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// Clear the cache used by MPFR.
// -------------------------------------------------------------------
mpfr_free_cache ();
return(0);
}

// -------------------------------------------------------------------
// We compute the Hardy Z values here.  The loop is over the count of
// different 't' values to compute (based on hz.iCount). Inside the 
// loop, we call ComputeMain and ComputeRemainder and add those 
// together to get the Hardy Z value. We then printf the result,
// and then increase 't' by hz.dIncr and repeat hz.iCount times.
// -------------------------------------------------------------------
int ComputeHardyZ(struct HZ hz)
{
mpfr_t			t, tOver2Pi, T, Incr, N, P, Temp1, Main, Remainder, HardyZ;
int				i;
unsigned int	uiN;
long double		tFraction, dP, ldRemainder;

// -------------------------------------------------------------------
// Initiate use of the MPFR system.  
// -------------------------------------------------------------------
InitMPFR(hz);

// -------------------------------------------------------------------
// Initialize several MPFR variables.  The 't' value is obtained
// from the hz.tBuf string, and the amount to increment 't' is obtained
// from hz.dIncr.  Both are then stored in MPFR vaiables.  
// -------------------------------------------------------------------
mpfr_inits2 (hz.iFloatBits, t, tOver2Pi, T, Incr, N, P, Temp1, 
		Main, Remainder, HardyZ, (mpfr_ptr) 0);
mpfr_set_str (t, hz.tBuf, 10, MPFR_RNDN);
mpfr_set_d (Incr, hz.dIncr, MPFR_RNDN);

// -------------------------------------------------------------------
// We loop hz.iCount times, processing 't' in the first (i = 1) loop.
// If hz.iCount is greater than one, we increment 't' by the Incr
// amount and loop again.  
//
// In the loop, we compute N and P for the given 't'.  We then call
// ComputeRemainder and ComputeMain to do the "heavy lifting" in
// computing the Hardy Z value.
// -------------------------------------------------------------------
for (i = 1; i <= hz.iCount; i++ ) {	
	// ---------------------------------------------------------------
	// Compute t/(2 pi) - needed in multiple places
	// ---------------------------------------------------------------	
	mpfr_div (tOver2Pi, t, my2Pi, MPFR_RNDN);

	// ---------------------------------------------------------------
	// Compute N and P in-line 
	// ---------------------------------------------------------------	
	mpfr_sqrt (T, tOver2Pi, MPFR_RNDN);
	mpfr_modf (N, P, T, MPFR_RNDN);
	uiN = (unsigned int) mpfr_get_ui (N, MPFR_RNDN);
	dP  = mpfr_get_ld (P, MPFR_RNDN);
	
	// ---------------------------------------------------------------
	// Compute the Remainder term.  The first step is to compute
	// t/(2* pi)]^{-1/4}.  We use the "square root of the reciprocal
	// then square root" method (unless iDebug % 2 == 0)
	// ---------------------------------------------------------------	
	if( !(hz.iDebug % 2 )) { // do this if iDebug is an even multiple of 2
	mpfr_set_ld (Temp1, -0.25, MPFR_RNDN);
	mpfr_pow (Temp1, tOver2Pi, Temp1, MPFR_RNDN);
	}
	else { // the default case
	mpfr_rec_sqrt (Temp1, tOver2Pi, MPFR_RNDN);
	mpfr_sqrt (Temp1, Temp1, MPFR_RNDN);
	}	
	// ---------------------------------------------------------------
	// We now call ComputeRemainder128, which uses a combination of
	// quadmath (128-bit floats) plus MPFR. But, (if iDebug % 3 == 0), 
	// we call ComputeRemainder, which uses only 80-bit long doubles.
	// ---------------------------------------------------------------		
	if( !(hz.iDebug % 3 )) { // do this if iDebug is an even multiple of 3
	tFraction = mpfr_get_ld (Temp1, MPFR_RNDN);
	ldRemainder = ComputeRemainder(uiN, tFraction, dP);	
	mpfr_set_ld (Remainder, ldRemainder,  MPFR_RNDN);
	}
	else { // the default case
	ComputeRemainder128(&Remainder, P, Temp1, uiN, hz);
	}	

	// ---------------------------------------------------------------
	// Now compute the Main term and add to Remainder to get HardyZ.
	// ---------------------------------------------------------------	
	ComputeMain(&Main, t, uiN, hz);
	mpfr_add (HardyZ, Main, Remainder, MPFR_RNDN);

	// ---------------------------------------------------------------
	// Print the results for our current t 
	// ---------------------------------------------------------------		
	hz.bVerbose ? mpfr_printf("For t = %.*Rf, Z(t) = %.*Rf \n", 
					hz.iOutputDPt, t, hz.iOutputDPz, HardyZ)
		: mpfr_printf("%.*Rf, %.*Rf \n", hz.iOutputDPt, t, hz.iOutputDPz, HardyZ);

	// ---------------------------------------------------------------
	// Increment 't' by Incr for the next pass (if any) through the loop.
	// ---------------------------------------------------------------	
	mpfr_add (t, t, Incr, MPFR_RNDN);
 	} // end of for loop

// -------------------------------------------------------------------
// We are done with the MPFR system.  Clear our local MPFR variables
// and then close down the MPFR system.
// -------------------------------------------------------------------
mpfr_clears (t, tOver2Pi, T, Incr, N, P, Temp1, 
		Main, Remainder, HardyZ, (mpfr_ptr) 0);
CloseMPFR();
return(1);
}


// -------------------------------------------------------------------
// We compute the main term of the Riemann-Siegel formula.  
// NOTE: 
//
// We pass in the Result variable (to hold the result of our calculation),
// the 't' to calculate, the previously calculated 'N' variable, 
// plus the HZ structore.
//
// First Step: compute theta(t).  
//
// Second Step: a loop (from n = 1 to n = N), where we calculate and add 
// together the results from each n:
//    
//		sqrt(1/n) * cos[theta(t) - t log n]
//
// Third Step: return 2 times the sum of those N terms.
// -------------------------------------------------------------------
int ComputeMain(mpfr_t *Result, mpfr_t t, unsigned int N, struct HZ hz)
{
mpfr_t		tOver2, PiOver8, LogOftOver2Pi;
mpfr_t		Recip48t, Power3Term, Temp1, Temp2, Theta;

// -------------------------------------------------------------------
// Step 0: Check the case N < 1 (nothing to do so return 0 in Result).
// -------------------------------------------------------------------
if(N < 1)
	{
	mpfr_set_ui (*Result, 0, MPFR_RNDN);
	return(1);
	}

// -------------------------------------------------------------------
// First step: compute theta(t).  The formula from our book is:
//	
//	dTheta = ((t/2) * (log(t / (2 * M_PI_X)))) - M_PI_X/8 - t/2
//		+ 1/(48 * t) + 7/(5760 * powl(t, 3));
// ------------------------------------------------------------------

// -------------------------------------------------------------------
// initialize all mpfr_t variables used in this first step
// -------------------------------------------------------------------
mpfr_inits2 (hz.iFloatBits, tOver2, PiOver8, LogOftOver2Pi, 
	Recip48t, Power3Term, Temp1, Temp2, Theta, (mpfr_ptr) 0);

// set tOver2
mpfr_div_ui (tOver2, t, 2, MPFR_RNDN);

// set LogOftOver2Pi
mpfr_div (Temp1, tOver2, myPi, MPFR_RNDN);
mpfr_log (LogOftOver2Pi, Temp1, MPFR_RNDN);

// set PiOver8
mpfr_div_ui (PiOver8, myPi, 8,  MPFR_RNDN);

// set Recip48t
mpfr_mul_ui (Temp1, t, 48, MPFR_RNDN);
mpfr_ui_div (Recip48t, 1, Temp1, MPFR_RNDN);

// -------------------------------------------------------------------
// We now have the pre-calculations necessary to compute all terms of 
// theta(t) (except for the Power3Term)
// -------------------------------------------------------------------
mpfr_mul (Theta, tOver2, LogOftOver2Pi, MPFR_RNDN);
mpfr_sub (Theta, Theta, PiOver8, MPFR_RNDN);
mpfr_sub (Theta, Theta, tOver2, MPFR_RNDN);
mpfr_add (Theta, Theta, Recip48t, MPFR_RNDN);

// -------------------------------------------------------------------
// We calculate and add the Powers3Term ONLY IF if the computed value 
// is large enough to matter (that is, only if t is not too large).
// -------------------------------------------------------------------
if(mpfr_cmp_ld (t, MAXT_POWER3_TERM) < 0)
	{
	mpfr_pow_ui (Temp1, t, 3, MPFR_RNDN);	
	mpfr_mul_ui (Temp1, Temp1, 5760, MPFR_RNDN);	
	mpfr_ui_div (Temp1, 1, Temp1, MPFR_RNDN);	
	mpfr_mul_ui (Power3Term, Temp1, (unsigned long int) 7, MPFR_RNDN);
	mpfr_add (Theta, Theta, Power3Term, MPFR_RNDN);
//	mpfr_printf("Power3Term = %.40Rf\n", Power3Term);
	}

// -------------------------------------------------------------------
// Second step: loop n = 1 to N.  We will handle the n = 1 and
// n = 2 cases separately before entering the loop.
// ------------------------------------------------------------------
mpfr_t	Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne;

mpfr_inits2 (hz.iFloatBits, 
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
	mpfr_sqrt_ui (Temp1, 2, MPFR_RNDN);
	mpfr_ui_div (RecipSqrtn, 1, Temp1, MPFR_RNDN);
	mpfr_mul (Temp1, t, myLog2, MPFR_RNDN);
	mpfr_sub (CosArg, Theta, Temp1, MPFR_RNDN);

if( hz.iDebug % 5) { // do this unless iDebug is an even multiple of 5
	//#########################################################################################	
	//----------------------------------------------------------------
	// divide CosArg by 2 pi and keep the remainder
	// Or, should we just let mpfr_cos do the same thing?
	//----------------------------------------------------------------
	mpfr_fmod (CosArg, CosArg, my2Pi, MPFR_RNDN);
	//#########################################################################################	
	}

	mpfr_cos (Temp1, CosArg, MPFR_RNDN);
	mpfr_mul (FullTerm, RecipSqrtn, Temp1, MPFR_RNDN);
	mpfr_add (Main, Main, FullTerm, MPFR_RNDN);
	mpfr_set (LognMinusOne, myLog2, MPFR_RNDN); // for Version 3, below
	}
	
// -------------------------------------------------------------------
// Now process the n = 3 through n = N terms
// -------------------------------------------------------------------
for (unsigned int n = 3; n <= N; ++n) {
	// ---------------------------------------------------------------
	// We will need an mpfr_t version of n in multiple places, so
	// save in Temp1
	// ---------------------------------------------------------------	
	mpfr_set_ui (Temp1, n, MPFR_RNDN);
	
	// ---------------------------------------------------------------
	// First, compute the square root of 1/n
	// ---------------------------------------------------------------	
	mpfr_rec_sqrt (RecipSqrtn, Temp1, MPFR_RNDN);	

if( !(hz.iDebug % 7 )) { // do this if iDebug is an even multiple of 7
	//#########################################################################################	
	// ---------------------------------------------------------------
	// Second, compute the argument to the cosine term. (Version 1)
	// That is, CosArg = [theta(t) - t log n].  
	// Then (further below) compute cos(CosArg).
	// ---------------------------------------------------------------	
	mpfr_log (Temp2, Temp1, MPFR_RNDN);
	mpfr_mul (Temp2, t, Temp2, MPFR_RNDN); 
	mpfr_sub (CosArg, Theta, Temp2, MPFR_RNDN); 	
	//#########################################################################################
	}
	else { // this is the default case:
	//#########################################################################################	
	// ---------------------------------------------------------------
	// Second, compute the argument to the cosine term. (Version 2)
	// We have saved CosArg(of n-1) and log(n-1), allowing:  
	// CosArg(of n) = CosArg(of n-1) + (t * [log(n-1) - log(n)])
	// ---------------------------------------------------------------	
	mpfr_log (Temp1, Temp1, MPFR_RNDN); 				// this sets Temp1 to log(n)
	mpfr_sub (Temp2, LognMinusOne, Temp1, MPFR_RNDN); 	// Temp2 = log(n-1) - log(n)
	mpfr_mul (Temp2, t, Temp2, MPFR_RNDN); 				// t * [log(n-1) - log(n)]
	mpfr_add (CosArg, CosArg, Temp2, MPFR_RNDN); 		// Add Temp 2 to CosArg(of n-1)
	mpfr_set (LognMinusOne, Temp1, MPFR_RNDN); 			// update for the next 'n' - use "swap"?
	//#########################################################################################
	}

if( hz.iDebug % 5) { // do this unless iDebug is an even multiple of 5
	//#########################################################################################	
	//----------------------------------------------------------------
	// divide CosArg by 2 pi and keep the remainder
	// Or, should we just let mpfr_cos do the same thing?
	//----------------------------------------------------------------
	mpfr_fmod (CosArg, CosArg, my2Pi, MPFR_RNDN);
	//#########################################################################################	
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
// We have calculated Main.  Now, multiply by 2 and then convert
// to a long double (which will be the return value)
// -------------------------------------------------------------------
mpfr_mul_2ui (Main, Main, 1, MPFR_RNDN);
mpfr_set (*Result, Main, MPFR_RNDN);

// -------------------------------------------------------------------
// Free the space used by the local mpfr (constant) variables
// -------------------------------------------------------------------	
mpfr_clears ( tOver2, PiOver8, LogOftOver2Pi, 
	Recip48t, Power3Term, Temp1, Temp2, Theta, (mpfr_ptr) 0);
mpfr_clears (Main, RecipSqrtn, CosArg, CosCalc, FullTerm, LognMinusOne, (mpfr_ptr) 0);

return(1);
}

