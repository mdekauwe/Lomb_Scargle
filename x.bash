	#!/bin/bash

X=155
Y=2
nice -19 lombscargle -dir "/users/eow/mgdk/DATA/REGRIDDED_LST/gridded_03/" -fname "lsta_daily_" -fext ".gra" -segs 3 -scargle -seed -31123 -of 4.0 -b -r $Y -c $X -nrows 332 	-ncols 667 -fp -my < input_files >x
cat x

