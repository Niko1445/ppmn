CC = gcc
CFLAGS = -O -std=gnu11 -fwrapv
LDLIBS = -lm

out.txt: main
	./main > out.txt
	cat out.txt

main: equal.o epsilon.o
	$(CC) $^ -o $@

equal.o: equal.c
	$(CC) -c $^ -o $@ $(CFLAGS)

epsilon.o: epsilon.c
	$(CC) -c $^ -o $@ $(CFLAGS)

.PHNOEY: clean
clean:
	$(RM) main out.txt *.o
