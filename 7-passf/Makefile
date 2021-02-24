CFLAGS = -Wall -Werror -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY: default
default: out.main.txt out.params.txt out.mygamma.png
out.main.txt:main
	./$< > $@
out.params.txt:params
	./$< > $@

out.mygamma.png: out.mygamma.txt Makefile
	graph \
	--output-format png --bitmap-size 1024x768 \
	--x-label "x" --y-label "y" \
	--top-label "Gamma function" \
	--line-mode 1 $< \
	> $@

out.mygamma.txt:mygamma
	./$< > $@

.PHONEY: clean
clean:
	$(RM) out* main params mygamma
