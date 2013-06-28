CC=g++
all:
	$(CC) GeneticTSP.cpp -o GeneticTSP -g -std=c++0x
clear_no:
	rm -f output_no_opt*.txt route_no_opt*.txt GeneticTSP
clear_opt:
	rm -f output_opt*.txt route_opt*.txt GeneticTSP
