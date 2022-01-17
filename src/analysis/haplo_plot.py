from os import replace
import numpy as np
from numpy.core.fromnumeric import size
import pandas as pd
import glob
import getopt, sys
import re
import matplotlib.pyplot as plt

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]
short_options = "hp:"
long_options = ["help", "path="]
try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    print (str(err))
    sys.exit(2)
for current_argument, current_value in arguments:
    if current_argument in ("-h", "--help"):
        print("-p/--path specifies the path for the .csv-files (add slash or backslash at the end)")
    elif current_argument in ("-p", "--path"):
        path = str(current_value)

def haploPlot(input: str):

#Importing and simplifying tree sequence

    df = pd.read_csv(input)

    def drawPieMarker(xs, ys, ratios, sizes, colors, ax=None):

        if ax is None:
            fig, ax = plt.subplots(figsize=(10,8))

        if(sum(ratios) <= 1): pass

        markers = []
        previous = 0
        # calculate the points of the pie pieces
        for color, ratio in zip(colors, ratios):
            this = 2 * np.pi * ratio + previous
            x  = [0] + np.cos(np.linspace(previous, this, 10)).tolist() + [0]
            y  = [0] + np.sin(np.linspace(previous, this, 10)).tolist() + [0]
            xy = np.column_stack([x, y])
            previous = this
            markers.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(sizes),'alpha': 0.6, 'facecolor':color, 'edgecolor': 'none'})

        # scatter each of the pie pieces to create pies
        for marker in markers:
            ax.scatter(xs, ys, **marker)
        
        return ax

    fig, ax = plt.subplots()

    for i in range(len(df.index)):
        if df.at[i,'Group'] == 123:
            drawPieMarker(xs=1+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p1Leaves)']/df.at[i,'n(Leaves)'], df.at[i,'n(p2Leaves)']/df.at[i,'n(Leaves)'], df.at[i,'n(p3Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['cyan', 'orange', 'maroon'],
            ax=ax)
        elif df.at[i,'Group'] == 12:
            drawPieMarker(xs=2+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p1Leaves)']/df.at[i,'n(Leaves)'], df.at[i,'n(p2Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['cyan', 'orange'],
            ax=ax)
        elif df.at[i,'Group'] == 13:
            drawPieMarker(xs=3+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p1Leaves)']/df.at[i,'n(Leaves)'], df.at[i,'n(p3Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['cyan', 'maroon'],
            ax=ax)
        elif df.at[i,'Group'] == 23:
            drawPieMarker(xs=4+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p2Leaves)']/df.at[i,'n(Leaves)'], df.at[i,'n(p3Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['orange', 'maroon'],
            ax=ax)
        elif df.at[i,'Group'] == 1:
            drawPieMarker(xs=5+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p1Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['cyan'],
            ax=ax)
        elif df.at[i,'Group'] == 2:
            drawPieMarker(xs=6+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p2Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['orange'],
            ax=ax)
        elif df.at[i,'Group'] == 3:
            drawPieMarker(xs=7+np.random.uniform(-0.25,0.25), ys=df.at[i,'Time'], 
            ratios=[df.at[i,'n(p3Leaves)']/df.at[i,'n(Leaves)']],
            sizes=df.at[i,'n(Leaves)'],
            colors=['maroon'],
            ax=ax)
        else:
            raise ValueError("wtf is this node?")

    parameters = re.findall("\d+", input)
    plt.axhline(y=int(parameters[-1]), color='green')
    plt.axhline(y=int(parameters[-3])+int(parameters[-1]), color='red')
    plt.xlim(0,8)
    ax.set_xticklabels(['','123','12','13','23','1', '2', '3'])
    plt.xlabel("Group")
    plt.ylabel("T[generations ago]")
    plt.savefig(inputFile.replace('.csv','.png'))
    plt.close()

for filepath in glob.iglob(r"{p}*.csv".format(p = path)):
    inputFile = filepath
    haploPlot(inputFile)