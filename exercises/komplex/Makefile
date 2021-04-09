CC = gcc
CFLAGS = -O -std=gnu99
LDLIBS = -lm

out.txt: main
	./main > out.txt
	cat out.txt

main: main.o komplex.o

main.o komplex.o: komplex.h

.PHNOEY: clean
clean:
	$(RM) main out.txt *.o
