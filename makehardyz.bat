gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c hardyz.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c computeMain.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c remainder.c 
gcc -Wall -Wextra -std=c99 -posix -c remainder128.c 
gcc --static -o hardyz hardyz.o computeMain.o remainder.o remainder128.o -l:libquadmath.a -l:libmpfr.a -l:libgmp.a
