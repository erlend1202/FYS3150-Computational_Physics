cppFiles: compile link run


compile:
	c++ -c main.cpp penningTrap.cpp utilities.cpp -std=c++11

link:
	c++ -o main.out main.o penningTrap.o utilities.o -larmadillo

run:
	./main.out

omp: omp_compile omp_link omp_run

omp_compile:
	c++ -c main.cpp penningTrap.cpp utilities.cpp -std=c++11 -O3 -Xpreprocessor -fopenmp
omp_link:
	c++ -o main.out main.o penningTrap.o utilities.o -larmadillo -lomp

omp_run:
	./main.out



ubuntu:	Ecompile Elink Erun

Ecompile:
	g++ -c penningTrap.cpp -std=c++11

Elink:
	g++ penningTrap.o -o penningTrap.exe -larmadillo

Erun:
	./penningTrap.exe
