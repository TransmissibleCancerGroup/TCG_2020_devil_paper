#/bin/bash
##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# This bash script launches the phylogenetic analysis. Requires IQTree (>=1.6.12). Required data
# files are provided in the repository that accompanies the paper.

iqtree -spp partition.best.models.nex -g constraint.tree -pre constrained -bb 1000 --alrt 1000 -wbtl
iqtree -spp partition.best.models.nex -pre unconstrained -bb 1000 --alrt 1000 -wbtl
