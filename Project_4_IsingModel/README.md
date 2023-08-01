# Project_4

- Repo URL : https://github.uio.no/jonathel/Project_4

## Authors

- Tov Uberg Tyvold (tovut@math.uio.no)
- Sophus B Gullbekk (sophusbg@math.uio.no)
- Jonathan Larsen (jonathel@math.uio.no) 
- Erlend Kristensen (erlek@math.uio.no)

## Compile and run

For Linux users, in the file <code>src</code> write <code>$make OMP</code>
For Mac users, in the file <code>src</code> write <code>$make omp</code>

To recreate the results we found in our report the <code>laticeSizeComp</code> function in <code>main.cpp</code> has to be run with 8 threads.

## File structure

All the code for the project is located in the <code>src</code> folder.

The code for implementing the Ising model is located in the <code>IsingModel.cpp</code> and <code>IsingModel.hpp</code> files.

Other helpful functions for the Ising model is located in the <code>utilities.cpp</code> and <code>utilities.hpp</code> files.

The random number generator is located in the <code>omp_rng.cpp</code> and <code>omp_rng.hpp</code> files.

The code for generating results is located in the <code>main.cpp</code> file.

The code for plotting the figures used in the report is located in <code>plot.py</code>.

The data from our simulations is stored in the <code>textfiles</code> folder.
