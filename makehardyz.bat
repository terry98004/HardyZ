gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c hardyz.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c comphardyz.c 
gcc --static -o hardyz hardyz.o comphardyz.o  -L. -l:libhgt.a -l:libquadmath.a -l:libmpfr.a -l:libgmp.a
