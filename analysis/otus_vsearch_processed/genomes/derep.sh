#!/bin/bash

#===============================================#
# Derep reads
#
# March 25 2016
#===============================================#

vsearch --derep_fulllength HL_16S_lib.fna \
--output HL_16S_lib.derep.fna --uc HL_16S_lib.derep.uc \
--sizeout