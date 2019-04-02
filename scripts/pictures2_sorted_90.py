#!/usr/bin/env python
# -*- coding: utf-8 -*-

# pictures.py
# code to draw graphs of output from clustering efforts
# based on efforts in ipython notebook, but adapted for
# use on the server.  Should generate a series of pdfs...
# 21/2/17 added code so that clusters are sorted and drawn
# based on fasta file size


import statistics
import matplotlib
matplotlib.use('pdf')
#%matplotlib inline
import numpy as np
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import glob, os

wd = '/users/jpjc1/biol-chong-2019/2019-02-05-CRISPRAD/james_madness/results90/'	# working directory
filenames = sorted(glob.glob(wd+'Cluster_*.fasta'), key=os.path.getsize, reverse=True) # sort based on fasta file size
file_number = len(glob.glob1(wd, 'Cluster_*.fasta'))


for tracker in range(0, file_number, 30):
    if (tracker/30) == int(file_number/30):
        bottom = tracker
        top = file_number
    else:
        bottom = tracker
        top = tracker + 30
    gs = gsp.GridSpec(5,6)
    gsplace = 0
    for e in filenames[bottom:top]:	# try first 30, but this then needs to loop
        f = str(e[:-5])+'csv' # e has .fasta extension, replace this with .csv
        data = np.genfromtxt(f, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80), delimiter=',')
        info_cov = np.genfromtxt(f, usecols=81, delimiter=',')
        info_len = np.genfromtxt(f, usecols=82, delimiter=',')
        info_gc = np.genfromtxt(f, usecols=83, delimiter=',')
        x_range = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80]  #data[0]
        y_range = data[1:]
        axes1 = plt.subplot(gs[gsplace])
        count = 0
        for lines in y_range:
            axes1.plot(x_range,lines)
            count += 1
        na = str(f[86:-4])  # affects name label[39:-4] [33:-4] for other samples
        nu = str(count)
        len = str('{0:.1f}'.format((sum(info_len[1:]))/1000))
        cov = str(int(statistics.mean(info_cov[1:])))
        covsd = str(int(statistics.stdev(info_cov[1:])))
        gc = str('{0:.1f}'.format(statistics.mean(info_gc[1:])))
        gcsd = str('{0:.1f}'.format(statistics.stdev(info_gc[1:])))
#    plt.axis('off')
        plt.semilogy(linewidth=0.5,mew=0.5)
        plt.axhline(y=0.01, ls='--', lw = 0.5, c = 'black')
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,x2+10,0.00001,1))
        plt.text(5, 0.4, na, fontsize = 4)
        plt.text(5, 0.0001, nu+' cov:'+cov+'+/-'+covsd, fontsize=4)
        plt.text(5, 0.00045, len+'kb', fontsize=4)
        plt.text(5, 0.00002, 'GC% '+gc+'+/-'+gcsd, fontsize=4)
        plt.tick_params(axis='both', labelsize=3, pad=-1, direction='out', length=2, width=0.5)
        plt.tick_params(right='off', top='off')
#    plt.tick_params(axis='x', which='minor', direction='out', length=1, width=0.5)
        plt.tight_layout()
        gsplace += 1
    strange = str(bottom)+'_'+str(top)
    if float(len) >= 1000:
        plt.savefig('pix_sorted_90/large/CRISPR_'+strange+'.pdf', type='pdf')
        print('Generated large '+strange)
        plt.close('all')
    elif float(len) >= 100 and float(len) < 1000:
        plt.savefig('pix_sorted_90/medium/CRISPR_'+strange+'.pdf', type='pdf')
        print('Generated medium '+strange)
        plt.close('all')
    else:
        plt.savefig('pix_sorted_90/small/CRISPR_'+strange+'.pdf', type='pdf')
        print('Generated small '+strange)
        plt.close('all')
print('Rendering complete!')
