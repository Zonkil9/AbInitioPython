# AbInitioPython
This is a simple Python program for basic *ab initio* calculations. It was written as a proof of concept and a test of the author's understanding of Hartree-Fock and post-Hartree-Fock methods. For a <ins>real</ins> Python implementation of computational chemistry methods, the author recommends [PySCF](https://github.com/pyscf/pyscf) and [PyBEST](https://fizyka.umk.pl/~pybest/).

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

You may encounter the error related to the `libgfortran.so.3` library, e.g.:

```sh
error while loading shared libraries: libgfortran.so.3: cannot open shared object file: No such file or directory
```

It is due to `libgfortran.so.3` is not distributed with recent Debian and Ubuntu distributions. If you are running the Debian machine,
you should execute the `compatibility.sh` script from the `compatibility` directory and provide the sudo password -- the script installs
the required library.

## Usage

Execute the command to start the program:
```sh
./abipy.sh
```

You can modify the `ZMAT` file to change the studied molecule and/or the basis set used in calculations.

## Acknowledgements

The author wishes to gratefully acknowledge the support and helpful advice of Prof. Stanisław A. Kucharski during the development of this code.

The `ACES2` Program is a Product of the Quantum Theory Project; University of Florida: Gainesville, FL, USA, 2005. For more, see their [GitHub repository](https://github.com/ajithperera/ACES-II).
