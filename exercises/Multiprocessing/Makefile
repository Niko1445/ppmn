CC = gcc
CFLAGS = -O -std=gnu99
CFLAGS+= -pthread
LDLIBS = -lm

out.txt: main
	./main > out.txt

main: main.c


.PHONEY: clean
clean:
	$(RM) main
