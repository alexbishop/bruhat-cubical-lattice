# Spanning cubical lattice

[![DOI](https://zenodo.org/badge/903949776.svg)](https://doi.org/10.5281/zenodo.14498455)

This repository contains a SageMath script which generates the bruhat graph for Coxeter groups, and attempts to find spanning cubical lattices.

## Usage

```
Usage: sage generate-lattice.sage [OPTIONS] TYPE NUM

Attempts to find a spanning cubical lattice in a given group

Arguments:

    TYPE    one of A, B, C, D, E, F, G, or H
    NUM     this is the rank of the finite Coxeter (sub)group

Options:

    --help              print these usage instructions
    --draw-graph        creates a graph of the embedding if one is found
    --affine=UPPER      use the affine Coxeter group. In this case, you
                         need to give an upper bound for the Bruhat graph

Example:
    
    sage generate-lattice.sage D 5

      verifies that is a cubical lattice for D5

    sage generate-lattice.sage --draw-graph --affine=01020102 A 2

      shows that there is no cubical lattice in the given interval

    sage generate-lattice.sage --draw-graph --affine=1021020 A 2

      shows that there is a cubical lattice in the given interval
```

## Requirements

In order to run this script, you will require the following.

 - SageMath
 - [Coxeter](http://math.univ-lyon1.fr/~ducloux/coxeter/coxeter3/english/coxeter3_e.html) by Fokko du Cloux
 - `coxeter3-sage` (see [PyPi](https://pypi.org/project/coxeter3-sage/))
