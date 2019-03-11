[![Build Status](https://travis-ci.org/hachmannlab/chemlg.svg?branch=master)](https://travis-ci.org/hachmannlab/chemlg)
[![codecov](https://codecov.io/gh/hachmannlab/chemlg/branch/master/graph/badge.svg)](https://codecov.io/gh/hachmannlab/chemlg)
# ChemLG – A Molecular and Materials Library Generator for the Enumeration and Exploration of Chemical Space
ChemLG is a smart and massive parallel molecular library generator for chemical and materials sciences.

Program Version: 0.2

Release Date: Feb 20, 2019

With contributions by:
Janhavi Abhay Dudwadkar (UB): Jupyter GUI

## Code Design:
ChemLG is developed in the Python 3 programming language and uses OpenBabel and its Python extension, Pybel for handling molecules. The development follows a strictly modular and object-oriented design to make the overall code as flexible and versatile as possible.

## Documentation:
ChemLG documentation can be found here https://hachmannlab.github.io/chemlg

## Installation and Dependencies:
The dependencies for ChemLG are OpenBabel and MPI4Py. It is recommended that these two dependencies are installed in a virtual environment prior to installing ChemLG. We suggest the following conda installations:


    conda create --name my_chemlg_env python=3.6
    source activate my_chemlg_env
    conda install -c openbabel openbabel
    conda install -c anaconda mpi4py
    
You can download ChemLG from Python Package Index (PyPI) via pip. 

    pip install chemlg


You can test the installation with:

    pytest -v



## Citation:
Please cite the use of ChemLG as:


    (1) M.A.F. Afzal, J.A. Dudwadkar, J. Hachmann, ChemLG – A Program Suite for the Generation of Compound Libraries and the Survey of Chemical Space, in preparation.
    (2) M.A.F. Afzal, J. Hachmann, ChemLG – A Molecular and Materials Library Generator for the Enumeration and Exploration of Chemical Space. Available from: https://hachmannlab.github.io/chemlg. 
    (3) J. Hachmann, M.A.F. Afzal, M. Haghighatlari, Y. Pal, Building and Deploying a Cyberinfrastructure for the Data-Driven Design of Chemical Systems and the Exploration of Chemical Space, Mol. Simul. 44 (2018), 921-929. DOI: 10.1080/08927022.2018.1471692

## Acknowledgement
ChemLG is based upon work supported by the U.S. National Science Foundation under grant #OAC-1751161. It was also supported by start-up funds provided by UB's School of Engineering and Applied Science and UB's Department of Chemical and Biological Engineering, the New York State Center of Excellence in Materials Informatics through seed grant #1140384-8-75163, and the U.S. Department of Energy under grant #DE-SC0017193.

## License and Copyright:
ChemLG is copyright (C) 2015-2018 Johannes Hachmann and Mohammad Atif Faiz Afzal, all rights reserved. 
ChemLG is distributed under 3-Clause BSD License (https://opensource.org/licenses/BSD-3-Clause).

(C) 2015-2018 Johannes Hachmann, Mohammad Atif Faiz Afzal
University at Buffalo - The State University of New York (UB)
Contact: hachmann@buffalo.edu, m27@buffalo.edu
http://hachmannlab.cbe.buffalo.edu
