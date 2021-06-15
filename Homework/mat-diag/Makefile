CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

particleInBox.png: out.txt Makefile
	echo '\
set terminal png;\
set output "$@";\
set key outside;\
set tics out;\
set xlabel "x";\
set ylabel "Wavefunction";\
set title "Particle in a box";\
plot \
   "data_B.txt" using 1:2 with line linewidth 3 title "Numeric n=0" \
  ,"data_B.txt" using 1:3 with line linewidth 3 title "Numeric n=1" \
  ,"data_B.txt" using 1:4 with line linewidth 3 title "Numeric n=2" \
	,"data_B.txt" using 1:5 with line title "Analytical n=0" \
  ,"data_B.txt" using 1:6 with line title "Analytical n=1" \
  ,"data_B.txt" using 1:7 with line title "Analytical n=2" \
' | pyxplot

out.txt: main
	./main > out.txt
	cat out.txt

.PHONEY: clean
clean:
	$(RM) main *.txt
