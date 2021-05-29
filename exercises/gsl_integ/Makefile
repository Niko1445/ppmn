CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)



out.txt: main
	./main > out.txt
	cat out.txt


main: main.c

erf.pyxplot.png: out.txt
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
	' | pyxplot

.PHNOEY: clean
clean:
	$(RM) main out.txt
