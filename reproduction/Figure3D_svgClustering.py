###############################################################################
# Perform clustering the identified SVGs from the AKI kidney and 3dst RA 
# Figure 3D, Supplementary 25
###############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster import hierarchy

######################################################################
# Kidney
######################################################################
df1 = pd.read_csv('/Users/wangjuex/projects/svgGranularity/kidneysvglist.txt',header=None)
genelist= df1[0].tolist()

expFilename = '/Users/wangjuex/kidneycounts/kidney085_XY01_20-0038/count.csv'
coordFilename = '/Users/wangjuex/kidneycounts/kidney085_XY01_20-0038/spa.csv'
df = pd.read_csv(expFilename)

df_svg = df[[c for c in df.columns if c in genelist]]
df_svg = df_svg.T
df_svg = df_svg.transform(lambda x: np.log(x+1))

## Plot Figure 3D
clusters = hierarchy.linkage(df_svg, method="ward")
plt.figure(figsize=(20, 6))
dendrogram = hierarchy.dendrogram(clusters, labels=df_svg.index, leaf_font_size=4)
## Plotting a horizontal line based on the first biggest distance between clusters 
# plt.axhline(90, color='red', linestyle='--'); 
## Plotting a horizontal line based on the second biggest distance between clusters 
# plt.axhline(40, color='crimson'); 
plt.savefig('kidneysvg.png',dpi=600)
plt.close()

clustering_model = AgglomerativeClustering(n_clusters=2, linkage="ward")
# clustering_model = AgglomerativeClustering(n_clusters=4, linkage="ward")
clustering_model.fit(df_svg)
labels = clustering_model.labels_

glist = df_svg.iloc[labels==0].index.tolist()
geneName=df_svg.iloc[labels==3].index[0]

with open('kidney_0_genes.txt', 'w') as f:
    f.write('\n'.join(glist))

glist = df_svg.iloc[labels==1].index.tolist()

with open('kidney_1_genes.txt', 'w') as f:
    f.write('\n'.join(glist))

######################################################################
# plot
######################################################################

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
    plt.autoscale()
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='gist_yarg', s=1)
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='YlGn', s=10)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.show()

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
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='YlGn', s=10)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.savefig(geneName+'.jpg')
    plt.close()

# Select part genes 
with open('/Users/wangjuex/projects/svgGranularity/kidneysvglist.txt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        # saveGenePlot(line)
        saveGenePlotLog(line)


######################################################################
# RA 3dst
######################################################################
df1 = pd.read_csv('/Users/wangjuex/projects/svgGranularity/RAsvglist.txt',header=None)
genelist= df1[0].tolist()

expFilename = '/Users/wangjuex/projects/svgGranularity/3dst/3dst_RA1_exp.csv'
coordFilename = '/Users/wangjuex/projects/svgGranularity/3dst/3dst_RA1_loc.csv'
df = pd.read_csv(expFilename)

df = df.T
df_svg = df[[c for c in df.columns if c in genelist]]
df_svg = df_svg.T
# df_svg = df_svg.transform(lambda x: np.log(x+1))

## Supplementary 25
clusters = hierarchy.linkage(df_svg, method="ward")
plt.figure(figsize=(20, 6))
dendrogram = hierarchy.dendrogram(clusters, labels=df_svg.index, leaf_font_size=1)
## Plotting a horizontal line based on the first biggest distance between clusters 
# plt.axhline(90, color='red', linestyle='--'); 
## Plotting a horizontal line based on the second biggest distance between clusters 
# plt.axhline(40, color='crimson'); 
plt.savefig('RAsvg.png',dpi=600)
plt.close()

# clustering_model = AgglomerativeClustering(n_clusters=2, linkage="ward")
clustering_model = AgglomerativeClustering(n_clusters=3, linkage="ward")
clustering_model.fit(df_svg)
labels = clustering_model.labels_

glist = df_svg.iloc[labels==0].index.tolist()
# geneName=df_svg.iloc[labels==3].index[0]

with open('RA_0_genes.txt', 'w') as f:
    f.write('\n'.join(glist))

glist = df_svg.iloc[labels==1].index.tolist()

with open('RA_1_genes.txt', 'w') as f:
    f.write('\n'.join(glist))

glist = df_svg.iloc[labels==2].index.tolist()

with open('RA_2_genes.txt', 'w') as f:
    f.write('\n'.join(glist))

glist = df_svg.iloc[labels==3].index.tolist()

with open('RA_3_genes.txt', 'w') as f:
    f.write('\n'.join(glist))

#########
# plot
#########

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
    plt.autoscale()
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='gist_yarg', s=1)
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='YlGn', s=10)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.show()

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
    plt.scatter(df1.x, df1.y, c=df1.z, cmap='YlGn', s=10)
    # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    # plt.scatter(df1.x, df1.y, c=df1.z, cmap='RdYlGn')
    plt.colorbar()
    plt.savefig(geneName+'.jpg')
    plt.close()

# Select part genes 
with open('/Users/wangjuex/projects/svgGranularity/kidneysvglist.txt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip()
        # saveGenePlot(line)
        saveGenePlotLog(line)