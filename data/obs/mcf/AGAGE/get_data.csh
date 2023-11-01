#!/bin/csh -f
# Download only CH3CCl3 data from the AGAGE data server

# Remove any existing files
rm ./*ch3ccl3_mon.txt

# Recursively download files matching the specific pattern for CH3CCl3
wget -r --no-parent -A '*ch3ccl3_mon.txt' https://agage2.eas.gatech.edu/data_archive/agage/gc-md/monthly/

# Move the downloaded files to the current directory and clean up
find agage2.eas.gatech.edu/data_archive/agage/gc-md/monthly/ -name '*ch3ccl3_mon.txt' -exec mv {} ./ \;
rm -rf agage2.eas.gatech.edu

exit 0

