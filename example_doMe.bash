#!/bin/bash
: '
	Example of the code working...
'
LOMB="./lombscargle"
INFILE=sunspot.dat
IMG_FILE=test.eps
OUTFILE=x
$LOMB < $INFILE > $OUTFILE

echo "generating plot"
WHERE=$(pwd)
IMG_FILE="$WHERE/test1.eps"
rm -f $IMG_FILE
gnuplot <<END
	set terminal postscript eps enhanced colour "Times" 16
	set output "$IMG_FILE"
	set size 1.5,1.5
	set origin 0.0, 0.0
	set multiplot
	set size 1.5, 0.75
	set origin 0.0, 0.75
	set yrange[*:*]
	set xrange[*:*]
	set ylabel "Sunspot Number"
	set xlabel "Year"
	plot \
		"sunspot.dat" u 1:2 t "" w l
	
	set size 1.5, 0.75
	set origin 0.0, 0.0
	set yrange[0:*]
	set xrange[0:0.5]
	set ylabel "Power"
	set xlabel "Frequency (cycles per year)"
	plot\
		"$OUTFILE" u 1:2 t "Bias corrected spectrum" w l,\
		"$OUTFILE" u 1:3 t "90 % significance level" w l 1 -1
END

display $IMG_FILE
#ggv $IMG_FILE
exit
