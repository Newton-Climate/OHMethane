#!/bin/csh -f
# Run this script to get the AGAGE MCF data
rm ./*-gcmd.mon
wget -r --no-parent -A '*-gcmd.mon' https://agage2.eas.gatech.edu/data_archive/agage/gc-md/monthly/
mv agage2.eas.gatech.edu/data_archive/agage/gc-md/monthly/*-gcmd.mon ./.
rm -rf agage2.eas.gatech.edu
exit(0)
