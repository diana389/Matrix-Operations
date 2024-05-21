# !!! Do not modify this file !!!
CC=gcc
CFLAGS=-Wall -Werror -O0
LIBS=-lsatlas
LIBDIRS=-L/usr/lib64/atlas

all: compare tema3_blas tema3_neopt tema3_opt_m

tema3_blas: solver_blas.c main.c utils.h
	$(CC) $(CFLAGS) $^ $(LIBDIRS) $(LIBS) -o $@

tema3_neopt: solver_neopt.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

tema3_opt_m: solver_opt.c main.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

compare: compare.c utils.h
	$(CC) $(OPT_CFLAGS) $^ -lm -o $@

neopt_test: solver_neopt.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

blas_test: solver_blas.c utils.h
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -rf compare tema3_blas tema3_neopt tema3_opt_m
