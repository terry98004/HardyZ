// -------------------------------------------------------------------
// Program last modified December 13, 2024. 
// Copyright (c) 2024 Terrence P. Murphy
// MIT License -- see hardyZ.c for details.
// -------------------------------------------------------------------
#include "hardyz.h"

// -------------------------------------------------------------------
// In support of computing the remainder term of the Riemann-Siegel 
// formula, we use a table from Haselgrove (here "CjCoeff128"), where he 
// computed power series coefficients of the Cj terms. As a reminder, 
// each Cj is an entire function, with the C0, C2, and C4 terms even 
// functions (only even coefficients are non-zero) and with C1 and 
// C3 odd.  
// We use the first 20 non-zero coefficients for C0 and C1 the first
// 19 for C2, the first 17 for C3 and the first 16 for C4. 
// (see "CjCoeffCount128").
//
// With a leading term of 0, the __float128 data type should be 
// accurate to about 36 decimal places.
// -------------------------------------------------------------------

const unsigned int	CjCoeffCount128[5] = {20, 20, 19, 17, 16};
const __float128	CjCoeff128[5][20] = {
{0.38268343236508977173Q,
 0.43724046807752044936Q,
 0.13237657548034352333Q,
-0.01360502604767418865Q,
-0.01356762197010358088Q,
-0.00162372532314446528Q,
 0.00029705353733379691Q,
 0.00007943300879521469Q,
 0.00000046556124614504Q,
-0.00000143272516309551Q,
-0.00000010354847112314Q,
 0.00000001235792708384Q,
 0.00000000178810838577Q,
-0.00000000003391414393Q,
-0.00000000001632663392Q,
-0.00000000000037851094Q,
 0.00000000000009327423Q,
 0.00000000000000522184Q,
-0.00000000000000033506Q,
-0.00000000000000003412Q},

{0.02682510262837535Q,
-0.01378477342635185Q,
-0.03849125048223508Q,
-0.00987106629906208Q,
 0.00331075976085840Q,
 0.00146478085779542Q,
 0.00001320794062488Q,
-0.00005922748701847Q,
-0.00000598024258537Q,
 0.00000096413224562Q,
 0.00000018334733722Q,
-0.00000000446708757Q,
-0.00000000270963509Q,
-0.00000000007785289Q,
 0.00000000002343763Q,
 0.00000000000158302Q,
-0.00000000000012120Q,
-0.00000000000001458Q,
 0.00000000000000029Q,
 0.00000000000000009Q},
 
{ 0.005188542830293Q,
  0.000309465838807Q,
 -0.011335941078229Q,
  0.002233045741958Q,
  0.005196637408862Q,
  0.000343991440762Q,
 -0.000591064842747Q,
 -0.000102299725479Q,
  0.000020888392217Q,
  0.000005927665493Q,
 -0.000000164238384Q,
 -0.000000151611998Q,
 -0.000000005907803Q,
  0.000000002091151Q,
  0.000000000178157Q,
 -0.000000000016164Q,
 -0.000000000002380Q,
  0.000000000000054Q,
  0.000000000000020Q,
  0.000000000000000Q},
  
 { 0.0013397160907Q,
  -0.0037442151364Q,
   0.0013303178920Q,
   0.0022654660765Q,
  -0.0009548499998Q,
  -0.0006010038459Q,
   0.0001012885828Q,
   0.0000686573345Q,
  -0.0000005985366Q,
  -0.0000033316599Q,
  -0.0000002191929Q,
   0.0000000789089Q,
   0.0000000094147Q,
  -0.0000000009570Q,
  -0.0000000001876Q,
   0.0000000000045Q,
   0.0000000000022Q,
   0.0000000000000Q,
   0.0000000000000Q,
   0.0000000000000Q},
   
{  0.00046483389Q, 
  -0.00100566074Q,
   0.00024044856Q,
   0.00102830861Q,
  -0.00076578609Q, 
  -0.00020365286Q, 
   0.00023212290Q, 
   0.00003260215Q, 
  -0.00002557905Q, 
  -0.00000410746Q, 
   0.00000117812Q, 
   0.00000024456Q, 
  -0.00000002392Q, 
  -0.00000000750Q, 
   0.00000000013Q, 
   0.00000000014Q, 
   0.00000000000Q, 
   0.00000000000Q, 
   0.00000000000Q, 
   0.00000000000Q }};


// -------------------------------------------------------------------
// The function computes a remainder FACTOR and a remainder 
// SUM, and then returns the product of FACTOR * SUM.
// -------------------------------------------------------------------
int ComputeRemainder128(mpfr_t *Result, mpfr_t P, mpfr_t tFraction, unsigned int N, struct HZ hz)
{
mpfr_t		Factor, Total, Temp1, Temp2;
__float128	P128, AdjP128;
__float128	PowersP[NUM_POWERS_P];

// -------------------------------------------------------------------
// initialize all mpfr_t variables
// -------------------------------------------------------------------
mpfr_inits2 (hz.iFloatBits, Factor, Total, Temp1, Temp2, (mpfr_ptr) 0);	

// -------------------------------------------------------------------
// Compute FACTOR = tFraction * (-1)^{N - 1}
// -------------------------------------------------------------------
mpfr_set (Factor, tFraction, MPFR_RNDN);
if(N % 2 == 0) {
	mpfr_neg (Factor, Factor, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// We next compute the SUM. Before we start, we must deal with P.
// For now, we are using the gcc "quadmath" 128-bit float for P
// instead of a mpfr_t variable, so we convert P to P128.
//
// Because we are using Haselgrove's table of coefficients, we need
// to adjust P128, with AdjP128 = 1 - (2 * P128). Also, in the formula, 
// we will need to compute AdjP128 to powers 0 through 39, most of them 
// multiple times.  For efficiency, we will do the computations 
// once and store the results for later use.
// -------------------------------------------------------------------
P128 = mpfr_get_float128 (P, MPFR_RNDN);
AdjP128 	= (1 - (2 * P128));
PowersP[0] 	= 1;		// initialize the '0' slot; AdjP128^{0} = 1

for(unsigned int k=1; k < NUM_POWERS_P; k++) {
	PowersP[k] = PowersP[k-1] * AdjP128;
	}

// -------------------------------------------------------------------
// We are now ready to compute the SUM by adding the individual Cj 
// terms. The j-loop is over the 5 individual Cj terms. Inside the
// j-loop, the i-loop is over the Haselgrove coefficients for the
// given Cj.  If j is even, Cj is an even function and the
// coefficients go in power series slots 0, 2,...; if j is odd, 
// they go in power series slots 1, 3,...
// After computing each Cj, we multiply that Cj by tFraction^{2j}
// and adding that value to the Total value.
// -------------------------------------------------------------------
__float128	 	Cj128;
unsigned int	i, j, uiEvenOdd, uiPowersJ = 0;

mpfr_set_zero (Total, 1);

for(j=0; j <= 4; j++, uiPowersJ += 2) {
	for(i=0, Cj128 = 0; i < CjCoeffCount128[j]; i++) {
		uiEvenOdd = (j % 2);
		Cj128 += (CjCoeff128[j][i] * PowersP[(2 * i) + uiEvenOdd]);
		}
	mpfr_pow_ui (Temp1, tFraction, uiPowersJ, MPFR_RNDN);
	mpfr_set_float128 (Temp2, Cj128, MPFR_RNDN);
	mpfr_mul (Temp1, Temp1, Temp2, MPFR_RNDN);
	mpfr_add (Total, Total, Temp1, MPFR_RNDN);
	}

mpfr_mul (*Result, Total, Factor, MPFR_RNDN);
mpfr_clears (Factor, Total, Temp1, Temp2, (mpfr_ptr) 0);	
return(1);
}

// -------------------------------------------------------------------
// For testing, this does a printf of a __float128 value.
// -------------------------------------------------------------------
void Show128(__float128 x)
{
char buf[128];

quadmath_snprintf (buf, sizeof buf, "%Qf", x);
printf("%s \n", buf);	
}
