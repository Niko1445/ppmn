CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

leastsq.png: out.txt Makefile
	echo '\
set terminal png;\
set output "$@";\
set key top right;\
set tics out;\
set xlabel "time (days)";\
set ylabel "ln(activity)";\
set title "Least squares fit";\
plot \
  "fit_data.txt" using 1:2 with line title "fit" \
  ,"fit_data.txt" using 1:3 with line title "fit - delta c" \
  ,"fit_data.txt" using 1:4 with line title "fit + delta c" \
 ,"tab_data.txt" using 1:2:3 with yerrorbars pointsize 0.5 title "data points" \
' | pyxplot

out.txt: main
	./main > out.txt
	cat out.txt

.PHONEY: clean
clean:
	$(RM) main *.txt
