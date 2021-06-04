CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

exp.pdf: exp.tex exp.pyxplot.png
	pdflatex $<
	pdflatex $<

exp.pyxplot.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Exponential function";\
	plot \
	 "$<" using 1:3 with line title "exp from math.h"\
	,"$<" using 1:2 with points pointtype 3 pointsize 0.5 title "Approximated homemade exp"\
	' | pyxplot

data.txt: main
	./main > data.txt

main: main.c

.PHONEY: clean
clean:
	main *.txt
