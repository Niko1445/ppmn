CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


out.txt: main
	./main > out.txt
	cat out.txt

.PHONEY: clean
clean:
	$(RM) main *.txt

ODE.png: out.txt Makefile
	echo '\
set terminal png;\
set output "$@";\
set key outside;\
set tics out;\
set xlabel "x";\
set ylabel "u";\
set title "u'' = -u";\
plot \
   "data_A3.txt" using 1:4 with line linewidth 1 title "rk32" \
' | pyxplot
