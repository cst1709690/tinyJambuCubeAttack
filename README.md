# Background

This repository stores the collection of C++ source codes for the implementation of cube attack on the first version of TinyJAMBU (Wu & Huang, 2019), a finalist cipher in the Lightweight Cryptography (LWC) Standardisation Project organised by the National Institute of Standards and Technology (NIST). Additional results that are not presented in the paper are also included in the repository.

# General Description of the Implementation for Each Attack

The implementation for each attack is separated into two phases: `Pre-processing Phase` and `Online Phase`.

- `Pre-Processing Phase`: The pre-processing phase computes the successful cubes with the parameters listed in the output file `Parameters.txt`. The successful cubes are cubes that passed the designated number of linearity tests and they are recorded in `successfulCubes.txt` with their corresponding coefficients recorded in `LHS.txt`. Each line in `successfulCubes.txt` corresponds to 32 lines of coefficients in `LHS.txt`. If a line in `LHS.txt` is -1, this indicates that the superpoly does not pass the designated number of linearity tests in that particular bit of the 32-bit keystream. Otherwise, each line in `LHS.txt` contains a constant, followed by the coefficients of the key bits, if any, that are present in the superpoly for that cube.

- `Online Phase`: The online phase verifies the cubes by checking the values of the distinguishers or the superpolies obtained in the pre-processing phase with a random key. The output files are divided into two groups of three files:
  - Filtered Cubes: If the values of the superpolies (or distinguishers) computed by the cube summation in the online phase are consistent with the values of the random key (or distinguishers), then the corresponding cubes are stored in `filteredCubes.txt`. The superpolies (and distinguisher values) are stored in `filteredSuperpolies.txt` whereas the values of the superpolies (and the distinguishers) computed in the online phase are recorded in `filteredRHS.txt`.
  - False Cubes: In contrast, If the values of the superpolies (or distinguishers) computed by the cube summation in the online phase are not consistent with the values of the random key (or distinguishers), the corresponding cubes are then stored in `falseCubes.txt`. Their respective superpolies and distinguishers are stored in `falseSuperpolies.txt`, whereas the values computed in the online phase are recorded in `falseRHS.txt`.

## Further Information on the Implementation of `DA2 and KRA2 (Reduced Cube Space)`

- `Pre-Processing Phase`: The pre-processing phase for this technique in the implementation for DA2 and KRA2 requires an input file called `cubeSpace.txt` that lists the cube spaces where the  new cubes are generated from. The formatting of `cubeSpace.txt` follows the formatting in the `successfulCubes.txt` file of the other implementation. In fact, the input file used in the experiments of this technique is actually the renamed `successfulCubes.txt` file from the other DA2 and KRA2 experiments. Additionally, the implementation of this technique outputs one extra file named `successfulCubeSpace.txt` to keep track of the cube spaces where the new cubes that passed the linearity tests originate from.

- `Online Phase`: Similarly, the online phase for the implementation of this technique outputs an additional pair of files, `falseCubeSpace.txt` and `filteredCubeSpace.txt`, that keeps track of the cube spaces where the cubes in the two groups ("filtered" and "false"; see `Online Phase` above) originate from.

## Source Codes for Each Attack

The source code for each phase contains three files according to the respective attack: `main.cpp`, `Cube.h` and `components.h`.

- `components.h`: Implementation of the functions neccessary for the execution of the cipher's algorithm (adapted from the original implementation created by the authors of the first version of TinyJAMBU, Wu and Huang (2019)).
- `Cube.h`: Implementation of cube attack for the cipher using the functions in `components.h`.
- `main.cpp`: Driver code using the implementations in `Cube.h` and `components.h` where the parameters of the experiments are set in this file.

# Additional Results for DA1, DA2, KRA1, and KRA2

Additional results for the attacks of DA1 and DA2 or KRA2 (both types of techniques) are provided in the Excel files under the folder `Additional Results`. Each worksheet in the files contains the results of an experiment with the provided parameters shown in the first few rows of the cells. For the parameters, note that NROUND1 and NROUND2 actually represent the values of r<sub>2</sub> and r<sub>1</sub> respectively, which are the notations used in the paper. Correspondingly, NROUND3 represents the value of r<sub>3</sub>. It is advisable to view the results without clicking the buttons in the worksheets. The results can be viewed even if the macros are disabled.

- `[DA1].xlsm`: The results in this file show all the cubes and cube testers with their corresponding non-constant superpolies and distinguisher values that are obtained for the experiments in DA1. The cube sizes vary from 3 to 20.
- `[DA2] Full Cube Space.xlsm`: The results in this file contain all of the 652 cube testers of size 15 and one cube tester of size 25.
- `[DA2] Reduced Cube Space (Technique 1).xlsm` and `[DA2] Reduced Cube Space (Technique 2).xlsm`: The results in these two files contain the cubes and cube testers with their corresponding superpolies and distinguisher values according to the two techniques discussed in the paper. Cubes that are highlighted in red indicate that the superpolies of the cubes fail the linearity tests when they are tested for a higher number of linearity tests than the number designated by the parameters.

## Discussion of Each Attack

For more information on each of the attacks, please refer to the paper.

# Reference

[First Version of TinyJAMBU (Wu & Huang, 2019).](https://csrc.nist.gov/CSRC/media/Projects/Lightweight-Cryptography/documents/round-1/spec-doc/TinyJAMBU-spec.pdf)

# Feedback

Any feedbacks or suggestions towards this work are welcomed. They can be sent to the email addresses of cst1709690@xmu.edu.my or iftekhar.salam@xmu.edu.my.