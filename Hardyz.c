// -------------------------------------------------------------------
// Program last modified January 9, 2026. 
// -------------------------------------------------------------------

/*
MIT License

Copyright (c) 2025-2026 Terrence P. Murphy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <time.h>
#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>			// for command line argument processing
#include <mpfr.h>

#include "hgt.h"
#include "hardyz.h"


const char sUsage[] = "Command Line Parameters\n" \
 "-t [positive number]	The t value for Z(t) - this parameter is required. (Digits and '.' only).\n" \
 "-i [positive number]	Amount to increment t (if checking multiple t values) - defaults to 1.\n" \
 "-c [positive integer]	Count of the number of t values to check - defaults to 1.\n" \
 "-p [positive integer]	Decimal point digits of Z(t) to show in report - defaults to 6.\n" \
 "-b [positive integer]	Floating point bits: 128 <= b <= 1024 - defaults to 256.\n" \
 "-d [positive integer]	Used for debugging only.  Please disregard.\n" \
 "-k [positive integer]	Number of threads to use - defaults to 1, maximum of 8.\n" \
 "-h			Show command line parameters.  All other parameters will be ignored.\n" \
 "-s			Report the total seconds taken to compute the Hardy Z values.\n"\
 "-v			Verbose report (otherwise CSV only)."; 

const char sCopyright[] = "Copyright 2025-2026 by Terrence P. Murphy." \
" Licensed under MIT License.\n\n"; 


int main( int argc, char *argv[] )  
{
int			c;
struct HZ 	hz;
int			tDecimalDigits = -1;
int			iDecimalDigits = 0;

hz.Verbose 		= false;
hz.ShowSeconds	= false;
hz.Count 		= 1;
hz.DebugFlags	= 2311; // (2 * 3 * 5 * 7 * 11) + 1
hz.OutputDPt 	= 0;
hz.OutputDPz 	= 6;
hz.DefaultBits	= HGT_PRECISION_DEFAULT;
hz.Threads		= 1;

strcpy(hz.incrBuf, "1");

fprintf(stderr, "%s", sCopyright);
if(argc == 1) {
	printf("%s\n", sUsage);
	exit(EXIT_FAILURE);
	}
opterr = 0; // To prevent _getopt from printing an error message on standard error

while ((c = getopt (argc, argv, "t:i:c:k:p:b:d:hvs")) != -1)
	switch (c)
		{
		case 'h':
			printf("%s\n", sUsage);		
//			intmax_t	iMax;
//			int64_t		i64;
//			printf("Sizeof intmax_t: %zu, Sizeof int64_t: %zu", sizeof(iMax), sizeof(i64));
			exit(EXIT_SUCCESS);
		case 'v':
			hz.Verbose = true;
			break;
		case 's':
			hz.ShowSeconds = true;
			break;			
		case 't':
			if (ValidateHardyT (optarg) < 1) {
				printf("Invalid argument to -t \n");
				return(EXIT_FAILURE);
				}
			strcpy(hz.tBuf, optarg);
			tDecimalDigits =  GetDecimalDigits(hz.tBuf);
			break;
		case 'i':
		if (ValidateIncr (optarg) < 1) {
				printf("Invalid argument to -i \n");
				return(EXIT_FAILURE);
				}
			strcpy(hz.incrBuf, optarg);
			iDecimalDigits =  GetDecimalDigits(hz.incrBuf);
			break;
		case 'c':
			hz.Count = ValidateCount(optarg);	
			if(hz.Count < 1){
				printf("Invalid argument to -c \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'k':
			hz.Threads = ValidateThreads(optarg);	
			if(hz.Threads < 1){
				printf("Invalid argument to -k \n");
				return(EXIT_FAILURE);
				}
			break;	
		case 'd':
			hz.DebugFlags = ValidateDebugFlags(optarg);	
			if( hz.DebugFlags < 1){
				printf("Invalid argument to -d \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'p':
			hz.OutputDPz = ValidateReportDecimalPlaces(optarg);	
			if( hz.OutputDPz < 1){
				printf("Invalid argument to -z \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'b':
			hz.DefaultBits = ValidatePrecisionMPFR(optarg);	
			if(hz.DefaultBits < 1) {
				printf("Invalid argument to -b \n");
				return(EXIT_FAILURE);
				}
			break;
		case '?':
		case ':':
		default:
			printf("Option -%c is either unknown or missing its argument\n", optopt);
			return (EXIT_FAILURE);
		}
if(tDecimalDigits == -1) {
	printf("The t parameter is required.\n");
	return(EXIT_FAILURE);
	}
hz.OutputDPt 	= tDecimalDigits > iDecimalDigits ? tDecimalDigits : iDecimalDigits;
hz.Threads		= hz.Threads > hz.Count ? hz.Count : hz.Threads; 

// -------------------------------------------------------------------
// We have finished validating the command line parameters.  Now compute 
// and printf our HardyZ results.  Report the time it takes to do the
// computations if the user enters the -s command line parameter.
// -------------------------------------------------------------------
clock_t	t;
double 	time_taken;

t = clock(); 
ComputeAllHardyZ(hz); 
t = clock() - t; 
time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 

if(hz.ShowSeconds){
	printf("Compute took %f seconds to execute \n", time_taken); 
	}
return(EXIT_SUCCESS);
}
