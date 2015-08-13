#!/bin/bash

cd ../tau
make clean
cd ../ellip
make clean
cd ../example
rm *.hed *.tbl ttim1.lis
