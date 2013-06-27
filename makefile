CC=g++
all:
	$(CC) GeneticTSP.cpp -o GeneticTSP -g
clear:
	rm -f output.txt route.txt GeneticTSP
