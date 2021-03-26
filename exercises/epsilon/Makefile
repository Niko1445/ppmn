CC = gcc
CFLAGS = -O -std=gnu11 -fwrapv
LDLIBS = -lm

out.txt: main
	./main > out.txt
	cat out.txt

main: main.c

.PHNOEY: clean
clean:
	$(RM) main out.txt
