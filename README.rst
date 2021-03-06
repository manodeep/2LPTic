**Warning** I generated 2048^3 ICs (subsampled from a 3072^3) with this IC code and those simulations were wrong. However, 1024^3 simulations generated with this were fine. Presumably, something is wrong with *(at least)* the subsample fraction logic - **do not use the subsampling features in this code** -- i.e., set both ``GlassTileFacSampleNumerator`` and ``GlassTileFacSampleNumerator = 1``.

Description
===========

This repo contains 2LPT IC generation code from Roman Scoccimarro. 
Original code can be found `here <http://cosmo.nyu.edu/roman/2LPT/>`__

Improvements
============

1. The code has been updated with lots of checks about INT overflow and 64 bit ID overflows. 
2. Dependency on Numerical Recipes in C has been replaced with ``gsl`` 
3. This code can accept P(k) generated by `CAMB <http://camb.info/>`__ (contribution from Greg Poole) 
4. Can generate consistent particle IDs across simulations of different resolutions (contribution from Greg Poole)

**Note** this version produces a header compatible with standard public Gadget2. In particular, 
the number of particles in excess of INT_MAX (~2e9) is stored in ``npartTotalHighWord`` in the header, rather than ``npartTotal[2]``.

Installation
============

Pre-requisites
--------------

1. An MPI capable compiler ``mpicc``
2. Double precision MPI FFTW2 libraries (``fftw2.1.5``)
3. ``gsl`` library

Compilation
-----------

1. Setup the ``Makefile`` options according to the simulation IC. Default options are matched to running a single-species dark matter only cosmological simulation. 
2. Check the paths to ``MPI``, ``FFTW`` and ``GSL`` libraries
3. Typing ``make`` should generate the executable 


Running 
=======

1. Edit the parameter file. See example parameter file provided in `example.params <example.params>`__
2. Launch the MPI process ``mpirun -np Ncpus ./2LPTic simulation.param`` (refer to your computing cluster documentation for submitting and running jobs)


Author
======

Author is Roman Scoccimarro (I think based on N-GenIC developed by
Volker Springel). Current repo contains contribution from 
`Greg Poole <https://github.com/gbpoole/>`__ for the CAMB P(k)
handling and consistent particle IDs generation. 

Maintained by `Manodeep Sinha <mailto:manodeep@gmail.com>`__. Any new bugs
are probably my fault. 





