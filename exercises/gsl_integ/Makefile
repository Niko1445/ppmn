CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


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
	 "erf_data.txt" using 1:2 with line title "integral erf"\
	' | pyxplot


out.txt: main
	./main > out.txt
	cat out.txt


main: main.c



.PHNOEY: clean
clean:
	$(RM) main out.txt
