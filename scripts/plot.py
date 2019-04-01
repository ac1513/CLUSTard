import statistics
import matplotlib
matplotlib.use('pdf')
#%matplotlib inline
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='usage = python entrez_down.py file_list_of_queries')
parser.add_argument('in_file', help='the name of the file containing a list of csv files', type=str)
parser.add_argument('prefix', help='prefix of the jobs', type=str)
args = parser.parse_args()
in_file = args.in_file
prefix = args.prefix

with open(in_file, 'r') as text_file:
    files = text_file.read().strip().split()

counter = 0
for i in range(0, len(files), 30):
    gs = gsp.GridSpec(5,6)
    gsplace = 0
    sub_files = files[i:i+30]
    counter += 1
    for file in sub_files:
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
    plt.savefig('plots/' + str(counter) + '_' + prefix + '_plot.pdf', type='pdf')
    print('Generated plot number ' + str(counter) + ' -> ' + str(counter) + '-' + prefix + '_plot.pdf')
    plt.close('all')
