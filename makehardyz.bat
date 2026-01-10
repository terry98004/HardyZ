gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c hardyz.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c comphardyz.c 
gcc -static -pthread -o hardyz hardyz.o comphardyz.o  -L. -l:libhgt.a -l:libmpfr.a -l:libgmp.a
