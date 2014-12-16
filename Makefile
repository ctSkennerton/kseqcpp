CC = g++
CFLAGS = -ggdb -Wall -O2
LIBS = -lbz2 -lz

test: test.cpp
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
