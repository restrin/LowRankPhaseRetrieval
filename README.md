# LowRankPhaseRetrieval

A collection of scripts for solving phase-retrieval problems using Gauge-Duality based on [Insert Link].

We want to recover a vectorized image ![equation](https://latex.codecogs.com/gif.latex?x%20%5Cin%20%5Cmathbb%7BR%7D%5En)
, given measurements of the form ![equation](https://latex.codecogs.com/gif.latex?b_i%20%3D%20%7Ca_i%5ET%20x%7C%5E2) for ![equation](https://latex.codecogs.com/gif.latex?i%3D1%2C%5Cdots%2Cm).

Given the primal PhaseLift problem:

![equation](https://latex.codecogs.com/gif.latex?%5Cmin_%7BX%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bn%5Ctimes%20n%7D%7D%20%5Cmbox%7Btr%7D%28X%29%20%5Censpace%5Cmbox%7Bsubject%20to%7D%5Censpace%20a_i%20X%20a_i%20%3D%20b_i)

our approach focuses on solving the gauge dual:

![equation](https://latex.codecogs.com/gif.latex?%5Cmin_%7By%20%5Cin%20%5Cmathbb%7BR%7D%5Em%7D%20%5Clambda_%7B%5Cmax%7D%5Cleft%28%5Csum_%7Bi%3D1%7D%5Em%20y_i%20a_i%20a_i%5ET%20%5Cright%29%20%5Censpace%5Cmbox%7Bsubject%20to%7D%5Censpace%20b%5ET%20y%20%3D%201)

# Repository Structure

## data

Contains various test problems and images used to construct them.

## methods

Contains several methods for solving phase retrieval problems:
- ```wirtinger_flow``` : Simple implementation of Wirtinger flow [(paper)](https://arxiv.org/abs/1407.1065)
- ```coorddescent```: Coordinate descent on gauge dual (supports uniform gradient descent, and an approximate Gauss-Southwell rule)
- ```reducedgrad```: Reduced gradient method on gauge dual
- ```projgrad```: Projected gradient method on gauge dual

## scripts

Contains scripts that run phase retrieval experiments:
- ```run_wf```: Runs Wirtinger Flow on given problem instance
- ```run_phase_retrieval```: Runs any of the above methods on a given phase retrieval instance
- ```initialize_wf```: Initializes Wirtinger Flow by running one of the above methods on the gauge dual first

## utils
- ```generate_problem```: Create a phase retrieval problem instance when given an image (using Hadamard measurements)
- ```opA```: Function used to represent measurement operator
- ```primal_from_dual```: Recovers primal solution ```x``` given an (approximate) dual solution ```y```
