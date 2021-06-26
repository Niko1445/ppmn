CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


annInterpolation.png: out.txt Makefile
	echo '\
set terminal png;\
set output "$@";\
set key outside;\
set tics out;\
set xlabel "x";\
set ylabel "y";\
set title "Neural Network Interpolation";\
plot \
   "data_1.txt" using 1:2 with line title "Random data" \
	 ,"points.txt" using 1:2 with points pointtype 2 title "interpolated point" \
 \
' | pyxplot


out.txt: main
	./main > out.txt
	cat out.txt

main: main.o vec_calc.o minimization.o

.PHONEY: clean
clean:
	$(RM) main *.txt
