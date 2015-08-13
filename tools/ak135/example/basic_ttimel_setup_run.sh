#!/bin/bash

# Build stuff
echo "Building stuff"
sleep 3
cd ../tau
make
cd ../ellip
make
cd ../example

# Setup tau tables
echo "Creating tau tables"
sleep 3
../tau/remodlv ak135 
../tau/setbrn

# Setup elcordir.tbl
echo "Creating elcordir tables (don't know what these do)"
sleep 3
../ellip/direct < ELCOR.dat

# Now run ttimel
echo "Running ttimel - we would need to do feed this with each event"
sleep 3
../ellip/ttimel




