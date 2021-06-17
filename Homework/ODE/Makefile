CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: ODE.png SIR.png

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

SIR.png: out.txt Makefile
	echo '\
set terminal png;\
set output "$@";\
set key outside;\
set tics out;\
set yrange [0:5.8e6];\
set xlabel "time (days)";\
set ylabel "Infectious";\
set title "SIR";\
plot \
  "data_SIR1.txt" using 1:5 with line linewidth 1 title "rk32 Tc = 6" \
	,"data_SIR2.txt" using 1:5 with line linewidth 1 title "rk32 Tc = 5" \
	,"data_SIR3.txt" using 1:5 with line linewidth 1 title "rk32 Tc = 4" \
	,"data_SIR4.txt" using 1:5 with line linewidth 1 title "rk32 Tc = 3" \
	,"data_SIR5.txt" using 1:5 with line linewidth 1 title "rk32 Tc = 2" \
	,"data_SIR6.txt" using 1:5 with line linewidth 1 title "rk32 Tc = 1" \
' | pyxplot

out.txt: main
	./main > out.txt
	cat out.txt

.PHONEY: clean
clean:
	$(RM) main *.txt
