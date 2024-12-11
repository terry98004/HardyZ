gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c hardyz.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c computemain.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c computetheta.c 
gcc -Wall -Wextra -std=c99 -posix -c remainder128.c 
gcc --static -o hardyz hardyz.o computemain.o computetheta.o remainder128.o -l:libquadmath.a -l:libmpfr.a -l:libgmp.a
