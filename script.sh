#!/bin/bash

# ./src/tema3_neopt input/input > output.txt
# ./src/tema3_opt_m input/input >> output.txt
# ./src/tema3_blas input/input >> output.txt

valgrind --tool=memcheck --leak-check=full src/tema3_blas input/input_valgrind 2> memory/blas.memory
valgrind --tool=memcheck --leak-check=full src/tema3_neopt input/input_valgrind 2> memory/neopt.memory
valgrind --tool=memcheck --leak-check=full src/tema3_opt_m input/input_valgrind 2> memory/opt_m.memory


# valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes src/tema3_blas input/input_valgrind 2> cache/blas.cache
# valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes src/tema3_neopt input/input_valgrind 2> cache/neopt.cache
# valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes src/tema3_opt_m input/input_valgrind 2> cache/opt_m.cache