#!/bin/bash

# script to look at spectra of the central england temperature record...
# sanity check on lombscargle code... comparable to Weedon 2003, pg. 169.

# data set from http://www.metoffice.gov.uk/research/hadleycentre/CR_data/Daily/HadCET_act.txt
# in file CET.txt


# edited CET.txt by hand to remove commentry and final year total info -> central_england_temp.txt

# summed data up so it's quaterly temp as in Weedon example...

flatten.awk central_england_temp.txt | gawk '{

												temp[NR] = $1
												quat     = 0.25
		
											} END {

												for( i=1; i<=NR; i+=3 ) {
	
													print quat, temp[i] + temp[i+1] + tem[i+2]
													quat += 0.25
	
												}
											}'  > tmp_quaterly

lombscargle -fconf99 < tmp_quaterly > tmp

# rescale spectra so it's relative power...
gawk '{
		
		day[NR] = $1
		pow[NR] = $2
		sig[NR] = $3
		a      += $2
		b      += $3
	
	} END {
		
		for ( i=1; i<=NR; i++ ) 
			print day[i], pow[i] / a, sig[i] / b
	}' tmp > tmp2

# rescale so plot can be made of log relative power vs. freq
gawk '{print $1, log($2) / log(10), log($3) / log(10)}' tmp2 > tmp

# make plot - can clearly see 1.0 and 0.5 year cycles...or annual and second harmonic ;P

echo "generating plot"
WHERE=$(pwd)
IMG_FILE="$WHERE/central_england_temp_spectra.eps"
rm -f $IMG_FILE
~/sun4u/bin/gnuplot <<END
	set terminal postscript eps enhanced colour "Times" 14
	set output "$IMG_FILE"
	set yrange[*:*]
	set xrange[*:*]
	set size 1.5,1.5
	set origin 0.0, 0.0
	set multiplot
	set size 1.5, 0.75
	set origin 0.0, 0.75
	set ylabel "Temperature (deg C)"
	set xlabel "Time 1559-2007 (quaterly month totals)"
	plot \
		"tmp_quaterly" u 1:2 t "" w l lw 3 lt 1 lc rgb "blue"
		
	set size 1.5, 0.75
	set origin 0.0, 0.0
	set yrange[*:*]
	set xrange[0:2.5]
	set ylabel "Log Relative Power"
	set xlabel "Frequency (Cycles per year)"
	plot \
		"tmp" u 1:2 t "" w l lw 3 lt 1 lc rgb "blue",\
		"tmp" u 1:3 t "99% confidence interval" w l lw 3 lt 1 lc rgb "red"
		
END

# analysis of sunspots data...
lombscargle -fconf99 < sunspot.dat > tmp

# rescale spectra so it's relative power...
gawk '{
		
		day[NR] = $1
		pow[NR] = $2
		sig[NR] = $3
		a      += $2
		b      += $3
	
	} END {
		
		for ( i=1; i<=NR; i++ ) 
			print day[i], pow[i] / a, sig[i] / b
	}' tmp > tmp2

# make plot - can clearly see 1.0 and 0.5 year cycles...or annual and second harmonic ;P

echo "generating plot"
WHERE=$(pwd)
IMG_FILE="$WHERE/sunspots_spectra.eps"
rm -f $IMG_FILE
~/sun4u/bin/gnuplot <<END
	set terminal postscript eps enhanced colour "Times" 14
	set output "$IMG_FILE"
	set yrange[*:*]
	set xrange[*:*]
	set size 1.5,1.5
	set origin 0.0, 0.0
	set multiplot
	set size 1.5, 0.75
	set origin 0.0, 0.75
	set ylabel "Sunspots Number"
	set xlabel "Year"
	plot \
		"sunspot.dat" u 1:2 t "" w l lw 3 lt 1 lc rgb "blue"
		
	set size 1.5, 0.75
	set origin 0.0, 0.0
	set yrange[*:*]
	set xrange[0:0.5]
	set ylabel "Relative Power"
	set xlabel "Frequency (Cycles per year)"
	plot \
		"tmp2" u 1:2 t "" w l lw 3 lt 1 lc rgb "blue",\
		"tmp2" u 1:3 t "99% confidence interval" w l lw 3 lt 1 lc rgb "red"
		
END
