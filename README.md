# Calculation of the Hubbard Model

### Warning: This code is untested for Windows and MacOS. Works on linux Ubuntu.

## Installation
This is a Python code. Use the latest version of language as of 21/12/2023 to avoid issues.
To run this code, you will need to install the following package(s) for your Python. The links to their documentation is also provided :
- numpy--> (https://numpy.org/install/)
- matplotlib--> (https://matplotlib.org/stable/)

## Description
This code originated in another project of mine. I decided to redo the entire thing in order to make the calculation faster, and also make the code way more readable. This code uses OOP in order to seperate distinct calculations in the grand scheme of things, and then uses the object has a holder for the results of the calculation.

The hubbard class calculates all the blocks associated with a specific system. Said blocks are obtained anatically by diagonalizing per matrix block the hamiltonian of the system. This class uses another class called "block". Objects instantiated by this class hold one of the block matrices of the hamiltonian, in addition to its lowest eigen value and its associated eigen vector. It also holds other useful attributes like the number of electrons, the spin, etc.

All blocks can be accessed using the getters, or by accessing the list "blocks" directly from the hubbard object.

The hubbard class needs to be fed a hopping matrix. Said matrix can be custom made, but it can also be automatically calculated by the function included in this code called "hopping_matrix", located in the `hop_matrix.py` file. Other parameters of the constructor of the hubbard class are : N - the number of sites, t - the hopping energy, U - the potential energy and mu - the chemical potential energy.

In the case that the calculation yields more than one block with the ground state energy, all are accessible through by methods, or by accessing the list "gs_blocks" directly. Two methods of the class can get you the block with the biggest or smallest spin.

A basic example of how to use the class is presented in `run.py`.
