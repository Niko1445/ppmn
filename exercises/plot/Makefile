CC = gcc
#CFLAGS = -O -std=gnu99
CFLAGS = $(shell gsl-config --cflags)
#LDLIBS = -lm
LDLIBS = $(shell gsl-config --libs)

all: data.txt erf.pyxplot.png gamma.pyxplot.png

erf.pyxplot.png: erf.data.txt erf.tab.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"erf.tab.data.txt" using 1:2 with points pointtype 3 title "tabulated data"\
	' | pyxplot

gamma.pyxplot.png: gamma.data.txt gamma.tab.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set yrange [-5:5];\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "gamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "homemade Gamma"\
	,"gamma.tab.data.txt" using 1:2 with points pointtype 3 title "tabulated data"\
	' | pyxplot


data.txt: main
	./main

main: main.o erf.o gamma.o


.PHNOEY: clean
clean:
	$(RM) main erf.data.txt gamma.data.txt *.o
