#! /bin/bash
for f in bpage*; do
	./pstojpeg.dws.sh -portrait -xsize 450 $f
done
