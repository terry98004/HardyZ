// -------------------------------------------------------------------
// Program last modified August 9, 2025. 
// Copyright (c) 2024-2025 Terrence P. Murphy
// MIT License -- see hardyz.c for details.
// -------------------------------------------------------------------
#include "hardyz.h"

extern mpfr_t	myPi, my2Pi;

// -------------------------------------------------------------------
// We compute \theta(t) as used in the Riemann-Siegel formula.  
// The formula from our book is:
//	
// Theta = ((t/2) * (log(t / (2 * PI)))) - PI/8 - t/2
//		+ 1/(48 * t) + 7/(5760 * powl(t, 3));
//
// Using the variable names below, the formula is:
//    Theta = tOver2 * [LogOftOver2Pi] - PiOver8 - tOver2
//         + Recip48t + Power3Term
// 
// Or, equally (as implemented below):
//    Theta = tOver2 * [LogOftOver2Pi - 1] 
//         + Recip48t - PiOver8 + Power3Term
// -------------------------------------------------------------------
int ComputeTheta(mpfr_t *Theta, mpfr_t t, int iFloatBits)	
{
mpfr_t		tOver2, PiOver8, LogOftOver2Pi;
mpfr_t		Recip48t, Power3Term, Temp1, MinorTerms;

// -------------------------------------------------------------------
// initialize all mpfr_t variables used in computing Theta
// -------------------------------------------------------------------
mpfr_inits2 (iFloatBits, tOver2, PiOver8, LogOftOver2Pi, 
	Recip48t, Power3Term, Temp1, MinorTerms, (mpfr_ptr) 0);

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
// Compute the minor terms in the \theta(t) formula:
//    MinorTerms = Recip48t - PiOver8 + Power3Term
// -------------------------------------------------------------------
mpfr_set (MinorTerms, Recip48t, MPFR_RNDN);
mpfr_sub (MinorTerms, MinorTerms, PiOver8, MPFR_RNDN);

// -------------------------------------------------------------------
// Calculate and add the Powers3Term UNLESS t is so large that the 
// computed value of this term will be too small to matter.
// We set Temp1 to the maximum allowed value and then do the
// comparison.
// -------------------------------------------------------------------
mpfr_set_str (Temp1, MAX_T_POWER3_STRING, 10, MPFR_RNDN);
if(mpfr_cmp (t, Temp1) < 0)
	{
	mpfr_pow_si (Temp1, t, -3, MPFR_RNDN);
	mpfr_mul_ui (Temp1, Temp1, 7, MPFR_RNDN);
	mpfr_div_ui (Power3Term, Temp1, 5760, MPFR_RNDN);
	mpfr_add (MinorTerms, MinorTerms, Power3Term, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// Now calculate the major term = tOver2 * [LogOftOver2Pi - 1]
// -------------------------------------------------------------------
mpfr_sub_ui (Temp1, LogOftOver2Pi, 1, MPFR_RNDN);
ThetaAxB (&Temp1, tOver2, Temp1, iFloatBits);
mpfr_add (*Theta, Temp1, MinorTerms, MPFR_RNDN);	

// -------------------------------------------------------------------
// Free the space used by the local mpfr (constant) variables
// -------------------------------------------------------------------	
mpfr_clears ( tOver2, PiOver8, LogOftOver2Pi, 
	Recip48t, Power3Term, Temp1, MinorTerms, (mpfr_ptr) 0);

return(1);
}

// -------------------------------------------------------------------
// In computing theta, there is a point where we multiply a large
// amount (t/2) by a smaller amount.  We do that calculation here.
// At the same time, our use of theta in all later computations 
// needs only the remainder after dividing by 2 * \pi.  
// We use that property here to keep the maximum value of this 
// multiplication well under the value of 't'.  This code is 
// in place of:
//				mpfr_mul (*Result, Big, Small, MPFR_RNDN);
// -------------------------------------------------------------------
int ThetaAxB(mpfr_t *Result, mpfr_t Big, mpfr_t Small, int iFloatBits)
{											
mpfr_t	FracBig2Pi, IntSmall, FracSmall, Temp1, Temp2;

mpfr_inits2 (iFloatBits, FracBig2Pi, IntSmall, FracSmall, 
   Temp1, Temp2, (mpfr_ptr) 0);

//----------------------------------------------------------------
// Divide Big by 2 pi and keep the fractional part.
// Obtain the integer and fractional part of Small.
//----------------------------------------------------------------
mpfr_fmod (FracBig2Pi, Big, my2Pi, MPFR_RNDN);
mpfr_modf (IntSmall, FracSmall, Small, MPFR_RNDN);

//----------------------------------------------------------------
// Multiply FracBig2Pi * IntSmall and then Big * FracSmall.
// Add the two results together.
//----------------------------------------------------------------
mpfr_mul (Temp1, FracBig2Pi, IntSmall, MPFR_RNDN);
mpfr_mul (Temp2, Big, FracSmall, MPFR_RNDN);
mpfr_add (*Result, Temp1, Temp2, MPFR_RNDN);

mpfr_clears ( FracBig2Pi, IntSmall, FracSmall, 
   Temp1, Temp2, (mpfr_ptr) 0);
return(1);
}
