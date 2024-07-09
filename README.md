# AbInitioPython
This is a simple Python program for basic *ab initio* calculations. It was written just as a proof of concept and a test of the author's understanding of Hartree-Fock and post-Hartree-Fock methods. <ins>**AbInitioPython** is slow and unintuitive.</ins> For a real Python implementation of computational chemistry methods, the author recommends [PySCF](https://github.com/pyscf/pyscf).

## Capabilities

The current version of **AbInitioPython** supports the following methods:
* RHF (restricted Hartree-Fock)
* MP2 (Moller-Plesset 2<sup>nd</sup> order)
* MP3 (Moller-Plesset 3<sup>rd</sup> order)
* LCCD (Linear Coupled Cluster with doubles)
* CID (Configuration Interaction with doubles)
* Basic calculations of electric properties using the finite-field approach

The program does not work on its own! It needs the integrals generated from the `ACES2` software package.

## Usage

Execute the command to start the program:
```sh
python main.py
```
