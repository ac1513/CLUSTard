import statistics
import matplotlib
matplotlib.use('pdf')
#%matplotlib inline
import numpy as np
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import glob, os
import pandas as pd

in_file = "test.txt"

with open(in_file, 'r') as text_file:
    files = text_file.read().strip().split()
    
if len(files) <= 30:
    gs = gsp.GridSpec(5,6)
    gsplace = 0
    for file in files:
        with open(file, 'r') as f:
            df = pd.read_csv(f, index_col='contig')
            av_cov = str('{0:.1f}'.format(statistics.mean(df['cover'].tolist())))
            sd_cov = str('{0:.1f}'.format(statistics.stdev(df['cover'].tolist())))
            tot_len = str('{0:.1f}'.format(sum(df['length'].tolist())/1000))
            av_gc =  str('{0:.1f}'.format(statistics.mean(df['GC'].tolist())))
            sd_gc = str('{0:.1f}'.format(statistics.stdev(df['GC'].tolist())))
            na = str(file.split('.')[0][8:])
            nu = str(len(df))
            axes1 = plt.subplot(gs[gsplace])
            x_range = [x for x in range(len(df.keys()[:-3]))]
            for index, row in df.iterrows():
                y_range = row[:-3]
                axes1.plot(x_range, y_range)
            plt.semilogy(linewidth=0.25,mew=0.25)
            x1,x2,y1,y2 = plt.axis()
            plt.axis((x1,x2+10,0.00001,1))
            plt.axhline(y=0.01, ls='--', lw = 0.5, c = 'black')
            plt.text(5, 0.4, na, fontsize = 3)
            plt.text(5, 0.0001, nu+' cov:'+av_cov+'±'+sd_cov, fontsize=3)
            plt.text(5, 0.00045, tot_len +'kb', fontsize=3)
            plt.text(5, 0.00002, 'GC% '+ av_gc +'±'+ sd_gc, fontsize=3)
            plt.tick_params(axis='both', labelsize=3, pad=-1, direction='out', length=2, width=0.5)
            plt.tick_params(right='off', top='off')
#            plt.tight_layout()
            gsplace += 1
    plt.savefig('picture.pdf', type='pdf')
    print('Generated plot')
    plt.close('all')
    
    
#    
#else:
#    for item in range(0, len(files), 30):
    