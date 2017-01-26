CPPFLAGS=-march=native -std=c++11 -m64 -O2 -g -Wextra -Wall -Wshadow
LDFLAGS=-lstdc++ -m64 -g -march=native -flto
LDFLAGS+=-L/usr/lib -lgsl -lgslcblas -lm
CC=g++

main: main.o newfunc.o oldfunc.o

profile: main
	valgrind --tool=callgrind --callgrind-out-file=profile.out -v --dump-instr=yes --trace-jump=yes --smc-check=all-non-file ./main

cache: main
	valgrind --tool=cachegrind --cachegrind-out-file=cache.out --cache-sim=yes --branch-sim=no ./main
	valgrind --tool=cachegrind --cachegrind-out-file=branch.out --cache-sim=no --branch-sim=yes ./main
