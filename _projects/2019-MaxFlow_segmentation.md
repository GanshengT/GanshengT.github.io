---
title: "MaxFlow_segmentation"
collection: projects
type: "innovation project"
permalink: /projects/2019-MaxFlow_segmentation
venue: "CentraleSupélec, MICS"
date: 2019-05-28
location: "Paris, France"
---




This project is served as an example about representing energy function as multi-layer graph.

MaxFlow_segmentation
============================

Project Data
-------------------------------------------------------
images available on request

Implementation
--------

|  | Entry                                                        | Description                                                 |
|--| --------                                                     |------------------------------------------------------------ |
|01| [Image restoration and stereo](/files/01_intro_optim_en.pdf) | Introduction to optimisation                                |



Tutorials
---------

|  | Entry                                                  | Description                                                 |
|--| --------                                               |------------------------------------------------------------ |
|01| [Tutorial 1 text](/files/TD1-algo_en.pdf)              | Simplex algorithm, examples, formulations                   |
|02| [Tutorial 2 text](/files/TD2_optim_en.pdf)             | Solving LP problems with spreadsheets. Duality              |
|03| [Tutorial 3 text](/files/TD3-algo_en.pdf)              | This tutorial is on integer programming                     |
|04| [Tutorial 4 text](/files/TD4_cs_en.pdf)                | This tutorial is *assignment 1*                             |
|05| [solving the TSP](/files/TD5-tsp.pdf)                  | Solving the TSP. This is *assignment 2*                     |


Solutions and code
---------

|  | Entry                                                  | Description                                                 |
|--| --------                                               |------------------------------------------------------------ |
|01| [Tutorial 1 solution](/files/TD1-solution.pdf)         | Solution to the first tutorial |
|02| [a Python Simplex solver](/files/simplexe.py)          | Basic, commented Simplex solver |
|03| [Sudoku solver](/files/Sudoku_ilp.ipynb)               | This code requires [cvxopt](http://cvxopt.org/install/index.html). |


Code 1
------

Here you will find verbose, straightforward, numpy-based code for the
simplex.

Here is the basic code, [a Python Simplex solver](/files/simplexe.py),
with no claim with respect to efficiency. Here is a 
[Python Notebook](/files/Simplexe.ipynb), with worked out examples.

I recommend you try the Python Notebook version. Here is the
[online rendering](https://nbviewer.jupyter.org/urls/hugues-talbot.github.io/files/Simplexe.ipynb)
of this notebook. 



Thanks
------

Special thank to Hugues Talbot





Challenges 
---------
* almost-continue discretising step for stereo image to find visual correspodance

Other relevant challenges will be posted here.

Alternative project
-------------------

The Kaggle challenge is obviously very challenging. An alternative project is to implement
a peg-solitaire solver.

Here are a couple of articles on how this might be done [Article 1](/files/Peg_Solitaire_1.pdf) ; [Article 2](/files/Peg_Solitaire_2.pdf).

