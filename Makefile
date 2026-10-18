CC ?= cc
CFLAGS ?= -O2 -Wall -Wextra

all: lzsa2d

lzsa2d: lzsa2.c lzsa2d.c lzsa2.h
	$(CC) $(CFLAGS) -o $@ lzsa2.c lzsa2d.c

# vectors made with upstream `lzsa -f 2`; each is also truncated and
# bit-flipped 2000 times under ASan/UBSan
test: lzsa2.c test.c lzsa2.h
	$(CC) -O1 -g -fsanitize=address,undefined -o $@ lzsa2.c test.c
	./test $(foreach f,$(wildcard tests/*.bin),$(f:.bin=.lzsa2) $(f))

clean:
	rm -f lzsa2d test

.PHONY: all clean
