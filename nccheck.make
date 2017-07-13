#!/bin/bash

gfortran -fbounds-check -g -Wextra -Wall -pedantic -o ../nccheck0820 netcdf.check.v0.8.20.f90 -I/home/marco/prog/include -L/home/marco/prog/lib -lnetcdff

