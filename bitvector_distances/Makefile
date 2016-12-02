CFLAGS = -O3 -Wall -Wextra -Werror
CC     = gcc

incFiles = bitvector_subsets.h utilities.h

default: bitvector_distances

bitvector_distances: bitvector_distances.o bitvector_subsets.o utilities.o

%.o: %.c Makefile ${incFiles}
	${CC} -c ${CFLAGS} $< -o $@

clean: cleano
	rm -f bitvector_distances

cleano:
	rm -f *.o

