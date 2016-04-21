#!/bin/bash

module load R

R CMD BATCH --vanilla "--args $SLURM_NTASKS" reuse-biowulf.R

