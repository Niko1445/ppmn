CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

out.txt: main
	./main > out.txt

main: main.c

.PHONEY: clean
clean:
	$(RM) main
