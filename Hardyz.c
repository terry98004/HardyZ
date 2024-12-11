// -------------------------------------------------------------------
// Program last modified December 10, 2024. 
// -------------------------------------------------------------------

/*
MIT License

Copyright (c) 2024 Terrence P. Murphy

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
#include "hardyz.h"

const char sUsage[] = "Command Line Parameters\n" \
 "-t [positive number]	The t value for Z(t) - this parameter is required. (Digits and '.' only).\n" \
 "-i [positive number]	Amount to increment t (if checking multiple t values) - defaults to 1.\n" \
 "-c [positive integer]	Count of the number of t values to check - defaults to 1.\n" \
 "-z [positive integer]	Decimal points of Z(t) to show in report - defaults to 6.\n" \
 "-b [positive integer]	Floating point bits: 128 <= b <= 1024, divisible by 64 - defaults to 256.\n" \
 "-d [positive integer]	Used for debugging only.  Please disregard.\n" \
 "-h			Show command line parameters.  All other parameters will be ignored.\n" \
 "-s			Report the total seconds taken to compute the Hardy Z values.\n"\
 "-v			Verbose report (otherwise CSV only)."; 

const char sCopyright[] = "Copyright 2024 by Terrence P. Murphy." \
" Licensed under MIT License.\n\n"; 


int main( int argc, char *argv[] )  
{
int			c, i;
const char  sDot[] = "."; 
//char *		end;
struct HZ 	hz;

hz.bVerbose 	= false;
hz.bSeconds		= false;
hz.iCount 		= 1;
hz.iDebug		= 2311; // (2 * 3 * 5 * 7 * 11) + 1
hz.iActualDPt 	= -1;
hz.iActualDPi 	= 0;
hz.iOutputDPt 	= 0;
hz.iOutputDPz 	= 6;
hz.iFloatBits	= MY_DEFAULT_PRECISION;
hz.iWholeDt 	= 0;
hz.iWholeDi 	= 0;
hz.iWholeD 		= 0;

strcpy(hz.incrBuf, "1");

fprintf(stderr, "%s", sCopyright);
if(argc == 1) {
	printf("%s\n", sUsage);
	exit(EXIT_FAILURE);
	}
opterr = 0; // To prevent _getopt from printing an error message on standard error

while ((c = getopt (argc, argv, "t:i:c:z:b:d:hvs")) != -1)
	switch (c)
		{
		case 'h':
			printf("%s\n", sUsage);
			exit(EXIT_SUCCESS);
		case 'v':
			hz.bVerbose = true;
			break;
		case 's':
			hz.bSeconds = true;
			break;			
		case 't':
			hz.iActualDPt = strCheckAndCount(optarg);
			if(hz.iActualDPt == -1 || strlen(optarg) > (T_BUF_SIZE -2)){
				printf("Invalid argument to -t \n");
				return(EXIT_FAILURE);
				}
			strcpy(hz.tBuf, optarg);
			hz.iWholeDt = strcspn(optarg, sDot);
			break;
		case 'i':
			hz.iActualDPi = strCheckAndCount(optarg);
			if(hz.iActualDPi == -1 || strlen(optarg) > (INCR_BUF_SIZE -2)){
				printf("Invalid argument to -i \n");
				return(EXIT_FAILURE);
				}
			strcpy(hz.incrBuf, optarg);
			hz.iWholeDi = strcspn(optarg, sDot);
			break;
		case 'c':
			hz.iCount = atoi(optarg);
			if(!(hz.iCount >= 1)){
				printf("Invalid argument to -c \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'd':
			i = atoi(optarg);
			if( i < 2 || i > DEBUG_MAX_VALUE ){
				printf("Invalid argument to -d \n");
				return(EXIT_FAILURE);
				}
			hz.iDebug = i;
			break;
		case 'z':
			hz.iOutputDPz = atoi(optarg);
			if(!(hz.iOutputDPz > 0)){
				printf("Invalid argument to -z \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'b':
			hz.iFloatBits = GetFloatBits(optarg);
			if(!(hz.iFloatBits > 0)){
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
if(hz.iActualDPt == -1) {
	printf("The t parameter is required.\n");
	return(EXIT_FAILURE);
	}
hz.iOutputDPt 	= hz.iActualDPt > hz.iActualDPi ? hz.iActualDPt : hz.iActualDPi;
hz.iWholeD 		= hz.iWholeDt > hz.iWholeDi ? hz.iWholeDt : hz.iWholeDi;

// -------------------------------------------------------------------
// We have finished validating the command line parameters.  Now compute 
// and printf our HardyZ results.  Time the seconds it takes to do the
// computations if the user enters the -s command line parameter.
// -------------------------------------------------------------------
if(hz.bSeconds){
	clock_t	t;
	double 	time_taken;

	t = clock(); 
	ComputeHardyZ(hz); 
	t = clock() - t; 
	time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
	printf("Compute took %f seconds to execute \n", time_taken); 
	}
else{
	ComputeHardyZ(hz); 
	}
return(EXIT_SUCCESS);
}

// -------------------------------------------------------------------
// We are given a string that must be validated as a number (i.e,
// only digits and at most one decimal point allowed -- no commas).  
// If validated, we return the number of decimal places in the number.
// -------------------------------------------------------------------
int strCheckAndCount(const char *str)
{
const char * ptr       = strchr(str, '.'); // point to decimmal point, if any
const char   sNumDot[] = "0123456789."; 
const char   sNumNZ[]  = "123456789"; 
int  		NonZero    = 0;
int	i;

for (i=0; !NonZero && i < 9; i++) {
	if(strchr(str, sNumNZ[i]) != NULL) {
		NonZero = 1;
		}
	}

if(!NonZero || strspn(str, sNumDot) != strlen(str) 
	   || strchr(str, '.') != strrchr(str, '.')) {
	return(-1); 		// either invalid char or too may '.'
	}
if(!ptr) return(0); 	// no decimal point so no decimal digits
return(strlen(++ptr));	// found decimal point, count decimal digits
}


// -------------------------------------------------------------------
// We are given a string that must be validated as a number. In this
// case, we are looking for a whole number (digits only, no decimal 
// point).  We are looking for a number between 128 and 1024, so the
// string length must be between 3 and 4 bytes.  Finally, the number
// must be divisible by 64; that is, 128, 192, 256, ..., 1024.
// If validated, we return that value.
// -------------------------------------------------------------------
int GetFloatBits(const char *str)
{
char * 		ptr;
size_t		Len;
int			iValue;
const char  sNum[] = "0123456789"; 

Len = strlen(str);

if(Len < 3 || Len > 4 || strspn(str, sNum) != Len){
	return(-1); 
	}

iValue = (int) strtol(str, &ptr, 10);
if(iValue % 64 != 0 || iValue < 128 || iValue > 1024){
	return(-1); 
	}
return(iValue);
}

// -------------------------------------------------------------------
// We check whether the value passed on the command line (using the
// -d debug parameter, and saved in hz.iDebug) "matches" the DebugNum 
// parameter passed in here.  There is a "match" if there is no 
// remainder when you divide hz.iDebug by Debug Num.
// -------------------------------------------------------------------
bool DebugMode(struct HZ hz, int DebugNum)
{
return(hz.iDebug % DebugNum == 0 ? true : false);
}
