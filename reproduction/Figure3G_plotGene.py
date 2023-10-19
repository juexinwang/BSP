#############################################################
# Plot Figure 3G of the manuscript
#############################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

coordFilename = '../SlideSeq2/spa.csv'
expFilename = '../SlideSeq2/count.csv'

coord=pd.read_csv(coordFilename)
t1= pd.read_csv(expFilename,index_col=0)

## Debug info
print("Loading ready.")

def visualGene(geneName):    
    exp = t1[geneName].tolist()
    # exp = exp/np.max(exp)
    df1 = pd.DataFrame({'x': coord['x'],
                    'y': coord['y'],
                    'z': exp})
    plt.figure(dpi=300)
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='gist_yarg', s=1)
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='YlGn', s=1)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.show()

def saveGenePlot(geneName):    
    exp = t1[geneName].tolist()
    # exp = exp/np.max(exp)
    df1 = pd.DataFrame({'x': coord['x'],
                    'y': coord['y'],
                    'z': exp})
    plt.figure(dpi=300)
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='gist_yarg', s=1)
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='bwr', s=1)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.savefig(geneName+'.jpg')
    plt.close()

def saveGenePlotLog(geneName):    
    exp = t1[geneName].tolist()
    expo = []
    for elem in exp:
        expo.append(np.log(float(elem+1)))
    df1 = pd.DataFrame({'x': coord['x'],
                    'y': coord['y'],
                    'z': expo})
    plt.figure(dpi=300)
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='gist_yarg', s=1)
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='YlGn', s=1)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.savefig(geneName+'.jpg')
    plt.close()

# Select part genes 
with open('selectgenelist.txt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        # saveGenePlot(line)
        saveGenePlotLog(line)
