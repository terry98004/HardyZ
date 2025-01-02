// -------------------------------------------------------------------
// Program last modified November 23, 2024. 
// Copyright (c) 2024 Terrence P. Murphy
// MIT License -- see hardyZ.c for details.
// -------------------------------------------------------------------
#include "hardyz.h"

// -------------------------------------------------------------------
// In support of computing the remainder term of the Riemann-Siegel 
// formula, we use a table from Haselgrove (here "CjCoeff"), where he 
// computed power series coefficients of the Cj terms. As a reminder, 
// each Cj is an entire function, with the C0, C2, and C4 terms even 
// functions (only even coefficients are non-zero) and with C1 and 
// C3 odd.  
// We use the first 20 non-zero coefficients for C0 and C1 the first
// 19 for C2, the first 17 for C3 and the first 16 for C4. 
// (see "CjCoeffCount").
// -------------------------------------------------------------------
const unsigned int	CjCoeffCount[5] = {20, 20, 19, 17, 16};
const long double	CjCoeff[5][20] = {
{0.38268343236508977173,
 0.43724046807752044936,
 0.13237657548034352333,
-0.01360502604767418865,
-0.01356762197010358088,
-0.00162372532314446528,
 0.00029705353733379691,
 0.00007943300879521469,
 0.00000046556124614504,
-0.00000143272516309551,
-0.00000010354847112314,
 0.00000001235792708384,
 0.00000000178810838577,
-0.00000000003391414393,
-0.00000000001632663392,
-0.00000000000037851094,
 0.00000000000009327423,
 0.00000000000000522184,
-0.00000000000000033506,
-0.00000000000000003412},

{0.02682510262837535,
-0.01378477342635185,
-0.03849125048223508,
-0.00987106629906208,
 0.00331075976085840,
 0.00146478085779542,
 0.00001320794062488,
-0.00005922748701847,
-0.00000598024258537,
 0.00000096413224562,
 0.00000018334733722,
-0.00000000446708757,
-0.00000000270963509,
-0.00000000007785289,
 0.00000000002343763,
 0.00000000000158302,
-0.00000000000012120,
-0.00000000000001458,
 0.00000000000000029,
 0.00000000000000009},	
 
{ 0.005188542830293,
  0.000309465838807,
 -0.011335941078229,
  0.002233045741958,
  0.005196637408862,
  0.000343991440762,
 -0.000591064842747,
 -0.000102299725479,
  0.000020888392217,
  0.000005927665493,
 -0.000000164238384,
 -0.000000151611998,
 -0.000000005907803,
  0.000000002091151,
  0.000000000178157,
 -0.000000000016164,
 -0.000000000002380,
  0.000000000000054,
  0.000000000000020,
  0.000000000000000},
  
 { 0.0013397160907,
  -0.0037442151364,
   0.0013303178920,
   0.0022654660765,
  -0.0009548499998,
  -0.0006010038459,
   0.0001012885828,
   0.0000686573345,
  -0.0000005985366,
  -0.0000033316599,
  -0.0000002191929,
   0.0000000789089,
   0.0000000094147,
  -0.0000000009570,
  -0.0000000001876,
   0.0000000000045,
   0.0000000000022,
   0.0000000000000,
   0.0000000000000,
   0.0000000000000},
   
{  0.00046483389, 
  -0.00100566074,
   0.00024044856,
   0.00102830861,
  -0.00076578609, 
  -0.00020365286, 
   0.00023212290, 
   0.00003260215, 
  -0.00002557905, 
  -0.00000410746, 
   0.00000117812, 
   0.00000024456, 
  -0.00000002392, 
  -0.00000000750, 
   0.00000000013, 
   0.00000000014, 
   0.00000000000, 
   0.00000000000, 
   0.00000000000, 
   0.00000000000 }};


// -------------------------------------------------------------------
// We have precomputed [t/(2* pi)]^{-1/4} and passed that value to
// the function as tFraction.  The function computes a FACTOR and 
// a SUM, and then returns the product of FACTOR * SUM.
// -------------------------------------------------------------------
long double ComputeRemainder(unsigned int N, long double tFraction, long double P)
{
// -------------------------------------------------------------------
// We will compute the FACTOR first. 
// -------------------------------------------------------------------
int			iPM     = N % 2 == 0 ? -1 : 1;  // iPM = (-1)^{N - 1}
long double	dFactor = iPM * tFraction; 		   

// -------------------------------------------------------------------
// We next compute the SUM. Before we start, we must deal with P.
// Because we are using Haselgrove's table of coefficients, we need
// to adjust P, with dAdjP = 1 - (2 * P). Also, in the formula, we 
// will need to compute dAdjP to powers 0 through 39, most of them 
// multiple times.  For efficiency, we will do the computations 
// once and store the results for later use.
// -------------------------------------------------------------------
long double	dAdjP        = (1 - (2 * P));
long double	dPowersP[40] = {1};	// initialize the '0' slot; dAdjP^{0} = 1

for(unsigned int k=1; k <= 39; k++) {
	dPowersP[k] = dPowersP[k-1] * dAdjP;
	}
// -------------------------------------------------------------------
// We are now ready to compute the SUM by adding the individual Cj 
// terms. The j-loop is over the 5 individual Cj terms. Inside the
// j-loop, the i-loop is over the Haselgrove coefficients for the
// given Cj.  If j is even, Cj is an even function and the
// coefficients go in power series slots 0, 2,...; if j is odd, 
// they go in power series slots 1, 3,...
// After computing each Cj, we multiply Cj by tFraction^{2j}
// and add that to the Total value.
// -------------------------------------------------------------------
long double 	dTotal = 0, dPowerJ = 0, dCj;
unsigned int	i, j, uiEvenOdd;

for(j=0; j <= 4; j++, dPowerJ += 2) {
	for(i=0, dCj = 0; i < CjCoeffCount[j]; i++) {
		uiEvenOdd = (j % 2);
		dCj += (CjCoeff[j][i] * dPowersP[(2 * i) + uiEvenOdd]);
		}
	dTotal += powl(tFraction, dPowerJ) * dCj;
	}
return(dFactor * dTotal);
}


//printf("Size of unsigned int: %zu bytes\n", sizeof(unsigned int));						4 bytes
//printf("Size of unsigned long int: %zu bytes\n", sizeof(unsigned long int));				4 bytes
//printf("Size of unsigned long long int: %zu bytes\n", sizeof(unsigned long long int));	8 bytes

