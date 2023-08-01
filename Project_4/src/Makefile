cppFiles: compile link run


compile:
	c++ -c main.cpp IsingModel.cpp utilities.cpp -std=c++11 

link:
	c++ -o main.exe main.o IsingModel.o utilities.o -larmadillo  

run:
	./main.exe

omp: omp_compile omp_link omp_run

omp_compile:
	c++ -c main.cpp IsingModel.cpp utilities.cpp omp_rng.cpp -std=c++11 -O3 -Xpreprocessor -fopenmp

omp_link:
	c++ -o main.exe main.o IsingModel.o utilities.o omp_rng.o -larmadillo -lomp

omp_run:
	./main.exe -v


OMP: OMPcompile OMPlink OMPrun

OMPcompile:
	c++ -fopenmp -c main.cpp IsingModel.cpp utilities.cpp omp_rng.cpp -std=c++11 

OMPlink:
	c++ -fopenmp -o main.exe main.o IsingModel.o utilities.o omp_rng.o -larmadillo  

OMPrun:
	./main.exe



