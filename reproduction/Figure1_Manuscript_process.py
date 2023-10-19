# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 23:42:59 2022

@author: lijinp
"""
import numpy as np
import pandas as pd
from scipy.spatial import distance



# load in data
InputPath = "./../../Data/Main_Figure/"
InputFiles = ["SVG", "NULL"]

for InputFile in InputFiles:
    InputExp =  pd.read_csv(InputPath + InputFile + ".csv")
    
    spatialMatrix = InputExp[['x','y']]
    expMatrix = InputExp.drop(['x','y'],  axis = 1)   
    expMatrix = expMatrix.transpose()
    expMatrix = expMatrix.to_numpy()
    spatialMatrix = spatialMatrix.to_numpy()
    OutputExp = InputExp
    for DK in [3, 10]:
        def patch(CenterPoint):
            CenterCoord = spatialMatrix[CenterPoint,:].reshape(1,-1)
            DistList = distance.cdist(CenterCoord, spatialMatrix, 'euclidean')
            PatchCells =[i for i,v in enumerate(DistList[0]) if (v <= DK and v>0.0)]
            if len(PatchCells)==0:               
                PatchCells.append(CenterPoint)
            return(PatchCells)    
        def patchmeans(CenterPoint):
            X_kpj = np.mean(expMatrix[:,PatchesCells[CenterPoint]],axis=1)
            return(X_kpj)
    
        PatchesCells = [patch(p) for p in np.arange(spatialMatrix.shape[0])] 
    
        X_kj = np.zeros([spatialMatrix.shape[0], expMatrix.shape[0]])
        for p in np.arange(spatialMatrix.shape[0]):
            X_kj[p,:] = patchmeans(p)  
        OutputExp['PatchMean_' + str(DK)] = X_kj
    OutputExp.to_csv("./../../Data/Main_Figure/"+ InputFile + "_Mean.csv")


 



