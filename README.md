Personal code for calculating the solvent-accessible surface area of a molecule
from GROMACS simulation output using a Monte Carlo technique. Random points are
generated around each site located at distance equal to the summation of the
solvent probe radius and VDW radius for that site. A generated point is rejected
if it is closer to another site than the site it was generated around. The ratio
of points accepted to the number of tries gives the fraction of the surface area
of a sphere for each atomic site, the summation of which is the
solvent-accessible surface area of the molecule.

## Requirements

Requires [libgmxfort](https://github.com/wesbarnett/libgmxfort) and
[fortran-json](https://github.com/jacobwilliams/json-fortran).

## Compilation

    $ make

## Installation

    # make install

## Usage

    $ sasa config.json

## Configuration file

Configuration files use the JSON format. The following options are available:

* *input.xtcfile* - The trajectory file from your simulation.
* *input.ndxfile* - The index file corresponding with the *xtcfile*, containing
  *ndxgroup*.
* *input.ndxgroup* - The index group indicating the solute molecule.
* *config.r* - Radius (nm) of solvent probe. The default is 0.14 nm.
* *config.rvdw* - Array of VDW radii for the solute, which will be added to the
  solvent probe radius in the calculation. These need to be in the exact order
of the sites listed in *ndxgrp*. Alternatively, only specify one VDW radius and
it will be used for all sites. If no VDW radii are specified, 0.2 nm is used for
each site.
* *config.nrand* - Number of attempts for inserting a point around each site.
  The default is 1000.
* *config.nblocks* - Number of blocks used in finding the blocked standard
  deviation, which is output as the uncertainty. The default is 5.
* *output.file* - Where the output is going to be saved. The default is
  "sasa.txt". Results are in nm².

## Example

An example configuration file and other input files is in the `example` folder.

To run the example after compilation do:

    $ ./bin/sasa example/config.json

The output file will be located at `example/sasa.txt`.
