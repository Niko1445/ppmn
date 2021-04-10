CC = gcc
CFLAGS = -O -std=gnu11
LDLIBS = -lm

all: read_cmdline read_stdin read_file
	./read_cmdline 1 2 3 4 5 > cmd.out.txt
	echo 1 2 3 4 5 | ./read_stdin > stdin.out.txt
	./read_file input.txt file.out.txt




.PHNOEY: clean
clean:
	$(RM) read_file read_stdin read_cmdline *out.txt
