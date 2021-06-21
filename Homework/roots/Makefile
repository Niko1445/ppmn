CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


out.txt: main
	./main > out.txt
	cat out.txt

main: main.o vec_calc.o

.PHONEY: clean
clean:
	$(RM) main *.txt
