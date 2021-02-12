CC = gcc
CFLAGS = -O -std=gnu11 -Wall
LDFLAGS = -I.
LDLIBS = -lm

out.txt: main vector
	./$< > $@
	./vector >> out.txt

main: main.c

vector: vector.c

clean:
	$(RM) main vector out.txt


