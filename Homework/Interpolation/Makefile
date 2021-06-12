CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)



all: spline.png integ.png

spline.png: out.txt Makefile
	echo '\
set terminal png;\
set output "$@";\
set key top left;\
set tics out;\
set xlabel "x";\
set ylabel "y";\
set title "Linear spline";\
plot \
  "data.txt" using 1:2 with line title "data" \
 ,"data.txt" using 1:2 with points pointsize 0.5 title "data points" \
 ,"data.txt" using 3:4 with points pointtype 3 pointsize 0.5 title "interpolated data" \
' | pyxplot

integ.png: out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Linear integral";\
	plot \
	"data.txt" using 1:2 with line title "data" \
	,"gsl_lin_integ.txt" using 1:3 with line title "my integrated data" \
	,"gsl_lin_integ.txt" using 1:2 with line title "gsl integrated data" \
' | pyxplot
	cat out.txt
out.txt: main
	./main > out.txt




.PHONEY: clean
clean:
	$(RM) main *.txt
