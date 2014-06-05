#!/usr/bin/env python
"""
Test out interpolation?	

"""

import numpy as np
import scipy, scipy.interpolate 
import matplotlib.pyplot as plt
import sys, random
from scipy import fftpack
import subprocess as sp

class SortInputFiles:
	
	def __init__(self):
		pass
	
	# run bash cmds
	def runBash(self, cmd):
		p   = sp.Popen(cmd, shell = True, stdout = sp.PIPE)
		out = p.stdout.read().strip()
	
		return out
	
	def fillUpArrays(self, filename):
	
		f = open(filename, 'r')
		data = []
		num_elements = 0
	
		for line in f:
			if line[0] == '#': 
				continue
			line = line[:-1]
			elems = map(float, line.split())
			if len(elems) == 2: 
				data.append(elems[0])
				data.append(elems[1])
				data.append(np.nan)
				data.append(np.nan)
			elif len(elems) == 4: 
				data.append(elems[0])
				data.append(elems[1])
				data.append(elems[2])
				data.append(elems[3])
			elif len(elems) == 0:
				data.append(num_elements)
				data.append(np.nan)
				data.append(np.nan)
				data.append(np.nan)
		
			num_elements += 1
		
		f.close()
		data = np.asarray(data).reshape(num_elements, 4)
		
		return (data)


if __name__ == "__main__":
	
	S = SortInputFiles()
	lst = S.fillUpArrays('agoufou_LST.dat')
	lst_lomb = S.fillUpArrays('agoufou_LST_SPEC.dat')
	tot_power_lst = lst_lomb[:,1].sum()
	x = lst[:,0]
	y = lst[:,1]
	
	""" interpolate btw data gaps """
	tck = scipy.interpolate.splrep(x, y) 
	x2 = np.arange(1,143, dtype=float)
	y2 = scipy.interpolate.splev(x2, tck) 
	
	""" write interpolated LST data to a file """
	f = open('agoufou_LST_interpolated.dat', 'w')
	for i in xrange(len(y2)):
		print >>f, x2[i], y2[i]
	f.flush()
	f.close
	
	""" run lombscargle code over interpolated data set """
	seed = random.randint(1, 332 * 667) * -1
	cmd = 'nice -19 lombscargle -seed ' + str(seed) + ' -of 4.0 < agoufou_LST_interpolated.dat > agoufou_LST_INTERPOLATED_SPEC.dat'
	S.runBash(cmd)
	
	""" read the data back in """
	lst_lomb_interpolated = S.fillUpArrays('agoufou_LST_INTERPOLATED_SPEC.dat')
	tot_power_lst_interpolated = lst_lomb_interpolated[:,1].sum()	
	""" change defaults dimensions, fonts etc """
	params = {'axes.labelsize' : 14,
			  'text.fontsize'  : 10,
			  'legend.fontsize': 8,
			  'xtick.labelsize': 12,
			  'ytick.labelsize': 12}
	plt.rcParams.update(params)
	
	
	plt.subplot(211)
	plt.title('Agoufou')
	plt.plot(x, y, 'ro', x2, y2, 'b-')
	plt.ylabel('LSTA (deg C)')
	plt.xlabel("Day of year since May 22nd 2006")
	plt.legend(('LSTA - with gaps', 'LSTA - interpolated'), 'best', numpoints = 1).draw_frame(True)
	plt.xlim(1, 142)
	
	plt.subplot(212)
	plt.plot(lst_lomb[:,0], lst_lomb[:,1] / tot_power_lst, 'g-',
			 lst_lomb[:,0], lst_lomb[:,2] / tot_power_lst, 'r--', 
			 lst_lomb_interpolated[:,0], lst_lomb_interpolated[:,1] / tot_power_lst_interpolated, 'b-',
			 lst_lomb_interpolated[:,0], lst_lomb_interpolated[:,2] / tot_power_lst_interpolated, 'r--') 
	plt.ylabel("Relative Power")
	plt.xlabel("Frequency (cycles per day)")
	plt.xlim(0.1, 0.5)
	plt.ylim(0, 0.01)
	plt.legend(('Spectrum - gaps', '90% confidence level', 'Spectrum - no gaps'), 'best').draw_frame(True)
	plt.savefig('Effect_of_interpolation.png', dpi = 150, figsize = (8,6))
	#plt.show()
	
	plt.clf()
	
	
	#Y = np.fft.fft(y2)
	#n = len(Y)
	#power = abs(Y[1:(n / 2)])**2
	#tot_power = power.sum()
	#nyquist = 1. / 2.
	#freq = np.array(xrange(0, (n / 2) - 1, )) / (n / 2.) * nyquist
	
	
	pxx, freqs = np.asarray(plt.psd(y2, NFFT=256, Fs=1, detrend=plt.mlab.detrend_linear, 
									window=plt.mlab.window_hanning, pad_to=1024,scale_by_freq = True))
	tot_pxx = pxx.sum()
	
	plt.subplot(211)
	plt.title('Agoufou')
	plt.plot(x, y, 'ro', x2, y2, 'b-')
	plt.ylabel('LSTA (deg C)')
	plt.xlabel("Day of year since May 22nd 2006")
	plt.legend(('LSTA - with gaps', 'LSTA - interpolated'), 'best', numpoints = 1).draw_frame(True)
	plt.xlim(1, 142)
	
	plt.subplot(212)
	plt.plot(lst_lomb[:,0], lst_lomb[:,1] / tot_power_lst, 'g-',
			 lst_lomb[:,0], lst_lomb[:,2] / tot_power_lst, 'r--')
	plt.plot(freqs, pxx / tot_pxx, 'b-')
	plt.legend(('Spectrum - lombscargle', '90% confidence level', 'Spectrum - Python FFT'), 'best').draw_frame(True)
	plt.xlim(0.1, 0.5)
	plt.ylim(0, 0.008)
	plt.savefig('FFT_vs_lombscargle.png', dpi = 150, figsize = (8,6))
	plt.show()
