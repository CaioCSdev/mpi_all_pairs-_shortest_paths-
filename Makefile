.PHONY: run

build: main.c
	mpicc -Wall --pedantic -g -o run main.c

run: build
	mpirun -np 4 run

