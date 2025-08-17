// -------------------------------------------------------------------
// Program last modified August 9, 2025. 
// Copyright (c) 2024-2025 Terrence P. Murphy
// MIT License -- see hardyZ.c for details.
// -------------------------------------------------------------------
#include "hardyz.h"

// -------------------------------------------------------------------
// In support of computing the remainder term of the Riemann-Siegel 
// formula, we use a table from Gabcke (here "coeffMPFR"), where he 
// computed power series coefficients of the first five Cj terms. 
// As a reminder, each Cj is an entire function, with the C0, C2, 
// and C4 terms even functions (only even coefficients are non-zero) 
// and with C1 and odd.  
// We use the first 44 non-zero coefficients for all five Cj terms,
// with each term having coefficients of 50 decimal places.
//
// With a leading (whole number) term of 0 or -0, the MPFR (256 bit) 
// data type should be accurate to more than 70 decimal places.
// -------------------------------------------------------------------
extern mpfr_t	coeffMPFR[5][44];

// -------------------------------------------------------------------
// The function computes a remainder FACTOR and a remainder 
// SUM, and then returns the product of FACTOR * SUM.
// -------------------------------------------------------------------
int ComputeRemainderMPFR(mpfr_t *Result, mpfr_t tOver2Pi, 
		unsigned int N, mpfr_t P, int iFloatBits)
{
mpfr_t		tFraction, Factor, Total, Temp1, Temp2, AdjP, Cj;
mpfr_t		PowersOfP[NUM_POWERS_P_GABCKE];
unsigned int	j;

// -------------------------------------------------------------------
// initialize all mpfr_t variables
// -------------------------------------------------------------------
mpfr_inits2 (iFloatBits, tFraction, Factor, Total, Temp1, Temp2, AdjP, Cj, (mpfr_ptr) 0);	

for(j = 0; j < NUM_POWERS_P_GABCKE; j++)	// init all P powers slots
	{
	mpfr_init2 (PowersOfP[j], iFloatBits);	// note: singular "init"
	}	

// ---------------------------------------------------------------
// Compute // t/(2* pi)]^{-1/4}.  We use the "square root of the 
// reciprocal and then the square root" method.
// ---------------------------------------------------------------	
mpfr_rec_sqrt (tFraction, tOver2Pi, MPFR_RNDN);
mpfr_sqrt (tFraction, tFraction, MPFR_RNDN);

// -------------------------------------------------------------------
// Compute FACTOR = tFraction * (-1)^{N - 1}
// -------------------------------------------------------------------
mpfr_set (Factor, tFraction, MPFR_RNDN);
if(N % 2 == 0) {
	mpfr_neg (Factor, Factor, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// We next compute the SUM. Before we start, we must deal with P.
//
// Due to the construction of the table of coefficients, we need
// to adjust P, with AdjP = 1 - (2 * P). Also, in the SUM formula, 
// we will need to compute AdjP to powers 0 through 87, most of them 
// multiple times.  For efficiency, we will do the computations 
// once and store the results for later use.
// -------------------------------------------------------------------
	// Compute AdjP = 1 - (2 * P).
mpfr_mul_ui (AdjP, P, (unsigned long int) 2, MPFR_RNDN);
mpfr_ui_sub (AdjP, (unsigned long int) 1, AdjP, MPFR_RNDN);

	// initialize the PowersOfP[0] slot; AdjP^{0} = 1
mpfr_init_set_ui (PowersOfP[0], (unsigned long int) 1, MPFR_RNDN);
	
	// Compute and save AdjP^{1} through AdjP^{87}
for(unsigned int k=1; k < NUM_POWERS_P_GABCKE; k++) {
	mpfr_mul (PowersOfP[k], PowersOfP[k-1], AdjP, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// We are now ready to compute the SUM by adding the individual Cj 
// terms. The j-loop is over the 5 individual Cj terms. Inside the
// j-loop, the i-loop is over the Gabcke coefficients for the
// given Cj.  If j is even, Cj is an even function and the
// coefficients go in power series slots 0, 2,...; if j is odd, 
// they go in power series slots 1, 3,...
// After computing each Cj, we multiply that Cj by tFraction^{2j}
// and add that value to the Total value.
// -------------------------------------------------------------------
unsigned int	i, uiEvenOdd, uiPowersJ;

mpfr_set_zero (Total, 1);	// Start with a Total of zero

	////////////////////////////////////
	// Outer loop calculates a Cj term 
	////////////////////////////////////
for(j=0, uiPowersJ = 0; j < NUM_Cj_TERMS; j++, uiPowersJ += 2) {
	// we are beginning a new Cj term, so set Cj = 0
	mpfr_set_zero (Cj, 1);
		//////////////////////////////////////////////////////////////////
		// Inner loop is through the non-zero coefficients of this Cj term
		//////////////////////////////////////////////////////////////////
	for(i=0; i < COEFF_PER_Cj; i++) {
		// Determine if non-zero coefficients are in even or odd slots 
		uiEvenOdd = (j % 2);
		// Multiply the correct AdjP power by the correct coefficient; save in Temp1
		mpfr_mul (Temp1, PowersOfP[(2 * i) + uiEvenOdd], coeffMPFR[j][i], MPFR_RNDN);
		// Add the saved Temp1 valur to Cj 
		mpfr_add (Cj, Cj, Temp1, MPFR_RNDN);
		} // end of inner loop 
	// back to the outer loop; set Temp1 to tFraction^{uiPowersJ} 
	mpfr_pow_ui (Temp1, tFraction, uiPowersJ, MPFR_RNDN);
	// Multiply Cj by Temp1 and again save in Temp1 
	mpfr_mul (Temp1, Temp1, Cj, MPFR_RNDN);
	// Add Temp1 to Total and go to the top of the outer loop 
	mpfr_add (Total, Total, Temp1, MPFR_RNDN);
	} // end of outer loop 

	// After exiting the outer loop, save Total in the passed Result variable
mpfr_mul (*Result, Total, Factor, MPFR_RNDN);
	// Clear the MPFR variables 
mpfr_clears (tFraction, Factor, Total, Temp1, Temp2, AdjP, Cj, (mpfr_ptr) 0);	
	
for(j = 0; j < NUM_POWERS_P_GABCKE; j++)	// clear all P powers slots
	{
	mpfr_clear (PowersOfP[j]);				// Note singular "clear"
	}		
return(1);
}
