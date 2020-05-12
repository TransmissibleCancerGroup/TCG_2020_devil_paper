## Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils - supporting code

This repository contains code required to reproduce figures and analysis in the paper "Evolution 
and lineage dynamics of a transmissible cancer in Tasmanian devils".

### Installation
The most direct way to get this repository onto your computer is using git:

    git clone --recursive https://github.com/TransmissibleCancerGroup/TCG_2020_devil_paper.git

The `--recursive` flag is there so that the submodules `dftdLowCov` and `additional_scripts` will
be downloaded all in one go. If `--recursive` is left out, it is still possible to get the
submodules by using the command `git submodule update --init --recursive`.

### R package
Some of the code relies on `dftdLowCov` being installed as an R package. This can be done from
an R terminal using the package `devtools`.

#### Install from the remote repository

    # Inside an R terminal...
    require(devtools)
    devtools::install_git(
        "https://github.com/TransmissibleCancerGroup/TCG_2020_devil_paper.git",
        subdir = "dftdLowCov",
        args = "--recursive")

#### Install from a local copy of the repository

    # Inside an R terminal...
    require(devtools)
    devtools::install_local("/path/to/dftdLowCov")  # replace with actual path


