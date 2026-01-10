CC = gcc
CFLAGS = -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c
LFLAGS = -static -pthread -L. -o 
SRCS = Hardyz.c CompHardyz.c 
LIBS = -l:libhgt.a -l:libmpfr.a -l:libgmp.a
OBJS = $(SRCS:.c=.o)
DEPS = hgt.h hardyz.h
TARGET = hardyz

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(TARGET) $(OBJS) $(LIBS)

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)
