CC=gcc


out.txt: hello
	./hello > out.txt

hello: hello.o
	$(CC) -o hello hello.o

hello.o: hello.c
	$(CC) -c hello.c

clean:
	$(RM) hello.o out.txt
