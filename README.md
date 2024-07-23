# AbInitioPython
This is a simple Python program for basic *ab initio* calculations. It was written just as a proof of concept and a test of the author's understanding of Hartree-Fock and post-Hartree-Fock methods. <ins>**AbInitioPython** is slow, unoptimized, and unintuitive.</ins> For a real Python implementation of computational chemistry methods, the author recommends [PySCF](https://github.com/pyscf/pyscf).

AbInitioPython is unable to prepare integrals needed in calculations on its own -- it uses integrals generated using the [ACES2](https://github.com/ajithperera/ACES-II) software.

## Capabilities

The current version of **AbInitioPython** supports the following methods:
* RHF (restricted Hartree-Fock)
* MP2 (Moller-Plesset 2<sup>nd</sup> order)
* MP3 (Moller-Plesset 3<sup>rd</sup> order)
* LCCD (Linear Coupled Cluster with doubles)
* CID (Configuration Interaction with doubles)
* Basic calculations of RHF, MP2 and MP3 electric properties using the finite-field approach

The program does not work on its own! It needs the integrals generated from the `ACES2` software package.

## Requirements

The program works only on the Linux 64bit machine. It also needs Python 3+ and the `numpy` package, which can be installed using `pip`:

```sh
pip install numpy
```

## Compatibility issues

If you ever encounter the error related to the `libgfortran.so.3` library, e.g.:

```sh
error while loading shared libraries: libgfortran.so.3: cannot open shared object file: No such file or directory
```
Then, you should execute the `compatibility.sh` script from the `compatibility` directory and provide the sudo password.

## Usage

Execute the command to start the program:
```sh
./abipy.sh
```

You can modify the `ZMAT` file to change the studied molecule and/or the basis set used in calculations.

## Disclaimer

The `ACES2` Program is a Product of the Quantum Theory Project; University of Florida: Gainesville, FL, USA, 2005. For more, see the LICENSE file in the `ACES2` directory and their [GitHub repository](https://github.com/ajithperera/ACES-II).
