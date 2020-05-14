#/bin/bash
##################################################################################################
# Supporting code for Evolution and lineage dynamics of a transmissible cancer in Tasmanian devils
# Author: Kevin Gori
# Date: May 2020
#
# WHAT THIS FILE DOES:
# This bash script uses IQTree (>=1.6.12) to run the approximately unbiased test (AU test) to
# reject any trees that are significantly worse than the maximum likelihood tree. Required data
# files are provided in the repository that accompanies the paper.

iqtree -spp partition.best.models.nex -z mltrees.txt -n 0 -zb 10000 -au -pre AUtest
