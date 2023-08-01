# Project 3: The Penning Trap
Third project in FYS3150 

- Repo url : https://github.uio.no/jonathel/Project_3

## Authors

- Tov Uberg Tyvold (tovut@math.uio.no)
- Sophus B Gullbekk (sophusbg@math.uio.no)
- Jonathan Larsen (jonathel@math.uio.no) 
- Erlend Kristensen (erlek@math.uio.no)

## Main concept
The purpose of this project is to simulate the phenomena of trapping particles using a Penning Trap with numerical methods. Reviewing the errors and the speed of our computation, we are able to make an estimate of how accurate our method really is. In the figures we attain models for the particles movement and the containment. 

## Technology
For this project we were required to run simulations that required a lot of evaluations. The approximations had to be precise, thus we needed a tool to go through several calculations. We also wanted to have parameters that were unscathed during the numerous evaluations and manageable across classes, so a <code>class</code> with <code>public</code> and <code>private</code> parameters were quite useful. For these reasons we obtained all the data with <code>c++</code>. We used version <code>c++11</code>

For all the data we gathered, we took use of Python as a tool for visualizing. We used version <code>Python 3.7.3</code>

## Code 
In this project we aimed to have structure within our code and a nice formatting for our readers. We included comments in our code which are easy to interpret and gives a nice idea of what is going on. 

_________________________

The main part of the code lies within folder src in <code>penningTrap.cpp</code> and <code>main.cpp</code>. Here we do all the necessary calculations regarding different scenarios such as interaction and a time-dependent potential field. We use the <code>class</code> <code>PenningTrap</code> to represent the Penning trap and 
the <code>class</code <code>Particle</code> to represent the particles inside the trap. 

_________________________

Having this in mind, we made a <code>Makefile</code> for linking and compiling all the <code>.cpp</code> files. It makes it easy to compile and link the full project, and hopefully motivates the reader to run our code. For compiling in the terminal one can navigate in the Terminal to the <code>src</code> folder and write <code>$ make</code>. This will write data to the <code>.txt</code> files and can then later be plotted. If you want to compile using OpenMp you can write 
<code>$ make omp</code>. 

_________________________

When the <code>Makefile</code> is ran, the file <code>plot.py</code> can plot the results written to text files.

_________________________

We utilized <code>utilities.cpp</code> for writing data regarding the 3D-plots with a header file <code>utilities.hpp</code> for linking.

_________________________

We also made a <code>test_analytical.py</code> file for testing our methods and implementation.

## Figures
All the figures that we have produced and put in our paper are located here.


## Latex
This folder contains the <code>.tex</code> and <code>.pdf</code> file along with the auxiliary files that latex needs.

 
