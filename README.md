# HardyZ
Calculate Values of the Hardy Z Function

## Overview

We provide here a C program to calculate *Z(t)*, the Hardy Z function, for 't' a positive real number.  Our program uses 
the Riemann-Siegel Formula.  

The included PDF file briefly describes the mathematics behind the calculations.  Here we will just mention
that, for any given 't', the calculations compute a "main" term *M* and an approximate "remainder" term *R(4)*,
with *R(4)* a very close approximation of the actual remainder term *R*. We then have *Z(t) = M + R(4)*.

Much greater detail on the mathematics is given in my book: 
*A Study of Riemann's Zeta Function* by Terrence P. Murphy.  The book is available on Amazon.  Go [here][website-link] to 
view the Table of Contents and Preface for all of my books.

The source code is intended for use with the **gcc** compiler and the **MPFR** floating point library.  Although only tested on 
64-bit Windows 11, I believe the source code is portable to other operating systems.

## Files

This distribution consists of the following files:

  * [README.md][readme-link]. The file you are currently reading.
  
  * [license.txt][license-link] A text file containing the MIT License which applies to this work.
	
  * [Hardyz.c][hardyz-c-link]. The `main` function entry pont for our program.  This source code file is 
  quite straightforward. We validate the user’s command line input, save their
  choices and then call the `ComputeHardyZ` function.
 
  * [remainder128.c][remainder128-c-link]. This source code file provides the 
  `ComputeRemainder128` function that computes the remainder term *R(4)*.
  We use the **gcc** quadmath 128-bit floating point library for the Haselgrove
  coefficients, and the **MPFR** library (set to 256-bit by default and changeable by the `-b` command line
  parameter) for all other floating point calculations.
 
  * [remainder.c][remainder-c-link]. This source code file provides the 
  `ComputeRemainder` function, an alternate function for computing the remainder term *R(4)*.
  The difference: here we use the **gcc** compiler’s built-in long double floating
  point type (with 80-bit precision on our computer system).
  Our testing showed no meaningful increase in accuracy when using `ComputeRemainder128`.
  However, the time cost of using `ComputeRemainder128` was very minimal. For that reason,
  we default to using `ComputeRemainder128`. To test using `ComputeRemainder` instead, use the `-d`
  command line parameter, setting the value to any positive integer that is both divisible by 3 and less than 211.
 
  * [computemain.c][computemain-c-link]. This source code file is “central control” for computing 
  *Z(t)* values. The `ComputeHardyZ` function is the entry point. In that function, we
  do needed initialization of the **MPFR** floating point system, call either `ComputeRemainder128`
  or `ComputeRemainder`, call the `ComputeMain` function (found in this source code file), and 
  print the results to `stdout`.
  
  * [hardyz.h][hardyz-h-link]. The is the only (local) include file for the program.  
  
  * [makehardyz.bat][makehardyz-bat-link]. The is the "makefile" for the program.  Currently,
  this file is a Windows batch file (**not** an actual makefile), but can be easily converted to 
  a standard makefile.
  
  * [hardyz.pdf][hardyz-pdf-link]. A PDF file with further discussion of our software, the building
  of our software program, and the mathematics behind our software. The PDF is created using LaTex, 
  needed to allow proper display / layout of the mathematics.

## Terms of use

This **Hardy Z Function Calculator** is free and distributed under the
[**MIT License**][license-link] (MIT). 

We used the [**gcc**][gcc-gnu-link] compiler and the [**MPFR**][mpfr-link] floating point library.
We also used the [**msys2**][msys2-link] software distribution and building platform for windows.
See their respective links for theirs terms of license.  

[website-link]:			https://riemann1859.com
[license-link]:			https://github.com/terry98004/HardyZ/blob/master/license.txt
[readme-link]:			https://github.com/terry98004/HardyZ/blob/master/README.md
[hardyz-c-link]:		https://github.com/terry98004/HardyZ/blob/master/hardyz.c
[remainder128-c-link]:	https://github.com/terry98004/HardyZ/blob/master/remainder128.c
[remainder-c-link]:		https://github.com/terry98004/HardyZ/blob/master/remainder.c
[computemain-c-link]:	https://github.com/terry98004/HardyZ/blob/master/computemain.c
[hardyz-h-link]:		https://github.com/terry98004/HardyZ/blob/master/hardyz.h
[makehardyz-bat-link]:	https://github.com/terry98004/HardyZ/blob/master/makehardyz.bat
[hardyz-pdf-link]:		https://github.com/terry98004/HardyZ/blob/master/hardyz.pdf
[mpfr-link]:			https://www.mpfr.org/
[gcc-gnu-link]:			https://gcc.gnu.org/
[msys2-link]:			https://www.msys2.org/
