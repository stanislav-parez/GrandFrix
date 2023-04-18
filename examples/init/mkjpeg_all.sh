#!/bin/bash
for f in bpage*.ps; do
    csh -f ./pstojpeg.dws.sh -portrait -xsize 1800 ${f%.*}
done
