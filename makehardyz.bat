gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c hardyz.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c computemain.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c computetheta.c 
gcc -Wall -Wextra -std=c99 -posix -c remainderMPFR.c 
gcc -Wall -Wextra -std=c99 -posix -c buildcoeff.c 
gcc --static -o hardyz hardyz.o computemain.o computetheta.o remainderMPFR.o buildcoeff.o -l:libmpfr.a -l:libgmp.a
