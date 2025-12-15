# HardyZ
Using the libHGT static library, calculate values of the Hardy Z function

## Overview

We provide here a C program to calculate *Z(t)*, the Hardy Z function, for 't' a positive real number.  The actual calculation work is done by calls to functions in the **libHGT** static library.  The code here is a front end to that calculating code that: (1) gathers the user-requested parameters for the calculations via the command line, (2) validates those command line parameters (by calls to the library), (3) passes those parameters to the calculating code in the library), and (4) prints (via stdout) a report of those calculations.

## Building the Executable

For Windows 11 users, an executable is included with any release posted on GitHub.

For other operating systems, you will need to build the executable, as follows.

*  You need the [**gcc**][gcc-gnu-link] C compiler installed on your system.  That installation must include the **libmpfr.a** and **libgmp.a** (floating point) static libraries.

*  From [**libHGT**][libhgt-link], you need to: (1) create the **libhgt.a** static library file, and (2) make that library file plus **hgt.h** visible to the **gcc** compiler.

*  Following the build logic in the **makehardyz.bat** file, you need to create the necessary 'makefile', in the form that applies to your operating system and the **gcc** compiler.

You can then build the HardyZ executable from the provided source files.

## Files

This distribution consists of the following files:

  * [README.md][readme-link]. The file you are currently reading.
  
  * [Hardyz.c][hardyz-c-link]. The entry pont for our program.  This source code file is 
  quite straightforward. We validate the user’s command line input, save their
  choices and then call the 'ComputeAllHardyZ' function, which in turn calls the **libHGT** functions that do the actual calculations.
 
  * [CompHardyz.c][CompHardyz-c-link]. This source code file includes: (1) the 'ComputeAllHardyZ' function, which in calls the **libHGT** functions that do the actual calculations, and (2) the 'HardyZCallback' function which is called from the **libHGT** functions and is used to print the output reports.
  
  * [hardyz.h][hardyz-h-link]. The is the only (local) include file for the program.  
  
  * [makehardyz.bat][makehardyz-bat-link]. The is the "makefile" for the program.  Currently,
  this file is a Windows batch file (**not** an actual makefile), but can be easily converted to 
  a standard makefile.

## Command Line Parameters

 *  -t [positive number]	The t value for Z(t) - this parameter is required. (Digits and '.' only).
 *  -i [positive number]	Amount to increment t (if checking multiple t values) - defaults to 1.
 *  -c [positive integer]	Count of the number of t values to check - defaults to 1.
 *  -p [positive integer]	Decimal point digits of Z(t) to show in report - defaults to 6.
 *  -b [positive integer]	Floating point bits: 128 <= b <= 1024 - defaults to 256.
 *  -d [positive integer]	Used for debugging only.  Please disregard.
 *  -k [positive integer]	Number of threads to use - defaults to 1, maximum of 8.
 *  -h			Show command line parameters.  All other parameters will be ignored.
 *  -s			Report the total seconds taken to compute the Hardy Z values.
 *  -v			Verbose report (otherwise CSV only).

## Terms of use

This **Hardy Z Function Calculator** is free and distributed under the
**MIT License** (MIT). 

We used the [**gcc**][gcc-gnu-link] compiler and the [**MPFR**][mpfr-link] floating point library.
We also used the [**msys2**][msys2-link] software distribution and building platform for windows.
See their respective links for theirs terms of license.  

[website-link]:			https://riemann1859.com
[readme-link]:			https://github.com/terry98004/HardyZ/blob/master/README.md
[hardyz-c-link]:		https://github.com/terry98004/HardyZ/blob/master/hardyz.c
[CompHardyz-c-link]:	https://github.com/terry98004/HardyZ/blob/master/CompHardyz.c
[hardyz-h-link]:		https://github.com/terry98004/HardyZ/blob/master/hardyz.h
[makehardyz-bat-link]:	https://github.com/terry98004/HardyZ/blob/master/makehardyz.bat
[mpfr-link]:			https://www.mpfr.org/
[gcc-gnu-link]:			https://gcc.gnu.org/
[msys2-link]:			https://www.msys2.org/
[libhgt-link]:			https://github.com/terry98004/libHGT/
