#!/bin/csh -f
# Run this script to get the GAGE MCF data
rm ./*-gage.mon
wget -r --no-parent -A '*-gage.mon' https://agage2.eas.gatech.edu/data_archive/gage/monthly/
mv agage2.eas.gatech.edu/data_archive/gage/monthly/*-gage.mon ./.
rm -rf agage2.eas.gatech.edu
exit(0)
