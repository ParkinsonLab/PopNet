'''
Created on Sep 18, 2013

@author: javi
'''
import os
from os import listdir
from os.path import isfile, join
import re
import numpy as np
from subprocess import call
import time
import random
import csv
import SimilarityContigs as simcon
import Proteins.HMMParser as hmp
import Proteins.DomainFamily as dmf
import Proteins.SequenceSelect as ss
import Proteins.SRSCytoscape as srsce
import multiprocessing as mp


def list_maker():
    return [1,2,3,4,5], [6,7,8,9,0]

def c(w, x, y, z):
    print((w, x, y, z))
    return x
    
def b(var_tuple):
    
    var_a = var_tuple[0]
    var_b = var_tuple[1]
    var_c = var_tuple[2]

#     var_a = 1
#     var_b = 2
#     var_c = 3

    var_d = 99
    
    return c(var_a, var_b, var_c, var_d)

def a(i):
    
    def list_maker2(var_c, list_1, list_2):
#         return [(var_c, x, list_2) for x in list_1]
        return [(var_c, i, list_2) for i in list_1]
    
    n = 3
    m = 4
    o = 5
    list_a, list_b = list_maker()

    pool = mp.Pool()
    param_list = list_maker2(n, list_a, list_b)
    
    results = pool.map(b, param_list)
    print(results)

if __name__ == '__main__':
    
    output = '/data/new/javi/test_graph.pdf'
    
    import matplotlib as mpl
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    
    bar_width = 0.2
    step = 2000
    inds = np.arange(5)
    
    title = 'test'
    
    fig = mpl.figure.Figure(figsize=(15,10))
    canvas = FigureCanvas(fig)
    fig.suptitle(title)
    fig.subplots_adjust(wspace = 0.5, hspace = 0.5)

    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('xaxis')
    ax.set_ylabel('yaxis')

#     bars = []
#     bottom = np.zeros((len(coords_list[0][1]),))
#     for i, (coords, stack_label) in enumerate(zip(coords_list, stack_label_list)):
# #         if i == 0:
# #             bar = ax.bar(inds, coords[1], color=colorTable[colorRange.index(coords[0])], align='center', linewidth=0)
# #         else:
# #             bar = ax.bar(inds, coords[1], color=colorTable[colorRange.index(coords[0])], align='center', bottom = bottom, linewidth=0)
#         
#         if i == 0:
#             bar = ax.bar(inds, coords[1], color=colorTable[i], align='center', linewidth=0)
#         else:
#             bar = ax.bar(inds, coords[1], color=colorTable[i], align='center', bottom = bottom, linewidth=0)
#             
#         bars.append(bar)
#         bottom = bottom + np.array(coords[1])    
    
    ax.set_xticks(inds)       
    ax.set_xticklabels([1,5,8,13,15], rotation=45, ha='right', fontsize=24)
#     ax.legend(bars, stack_label_list, loc=1, bbox_to_anchor=(1.0, 1), ncol=1, markerscale=0.2, fontsize=6)
#     fig.colorbar(sm)
    pdf = PdfPages(output)
    pdf.savefig(fig)
    pdf.close()



