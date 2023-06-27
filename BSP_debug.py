import numpy as np
import pandas as pd
import sys, time, datetime, os
from scipy import stats
from scipy.spatial import distance
from sklearn.preprocessing import minmax_scale
from scipy.stats import gmean
import glob
import argparse

parser = argparse.ArgumentParser(description='Detecting Spatial Variable Gene using Big-Small Patch algorithm')
# Input and output
parser.add_argument('--datasetName', type=str, default='MOB', help='dataset: subfolder under data/. MOB as the example')
parser.add_argument('--inputDir', type=str, default='data/', help='Directory of input: default(data/)')
parser.add_argument('--spaLocFilename', type=str, default='Rep11_MOB_spa.csv', help='spatial location inputfile: default(%(default)s). Format: Loc*(x,y) with header')
parser.add_argument('--expFilename', type=str, default='Rep11_MOB_count.csv', help='expression inputfile: default(%(default)s). Format: Loc*Gene with header')
parser.add_argument('--outputDir', type=str, default='./', help='Directory of output: default(./)')
parser.add_argument('--outputFilename', type=str, default='P_values.csv', help='Filename of output: default(%(default)s)')

# Parameters
parser.add_argument('--D1', type=float, default=1.0,  help='radius of small patch  (default: %(default)s)')
parser.add_argument('--D2', type=float, default=3.0, help='radius of big patch (default: %(default)s)')
parser.add_argument('--distType', type=str, default='euclidean',  help='euclidean/cityblock (default: %(default)s)')
parser.add_argument('--scaleFactor', type=float, default=1.0,  help='alpha>0 (default: %(default)s)')
parser.add_argument('--edgeFillingTag', action='store_true', default=False, help='Test: fill nodes on the edge')

# For 3D transcriptomics
parser.add_argument('--for3DTag', action='store_true', default=False, help='For 3D')

# IO for large scale analysis 
parser.add_argument('--oneFileTag', action='store_true', default=False, help='Use One file')
parser.add_argument('--inputFile', type=str, default='Rep11_MOB_merge.csv', help='Input file combining X Y and expression')
parser.add_argument('--useDirTag', action='store_true', default=False, help='Use all csv in the dir')

# Normalization
parser.add_argument('--normalize', type=str, default='minmax',  help='scale (default: %(default)s)')
parser.add_argument('--noinputCellTag', action='store_true', default=False,  help='Whether input file has no cell column (default: %(default)s)')
parser.add_argument('--notransTag', action='store_true', default=False,  help='Whether input file has transposed (default: %(default)s)')
parser.add_argument('--logTransform', action='store_true', default=False,  help='Whether logtransform (default: %(default)s)')

# Select Fitting Distribution
parser.add_argument('--fitDist', type=str, default='lognormal',  help='Select fitting distribution lognormal/beta  (default: %(default)s)')
parser.add_argument('--adjustP', action='store_true', default=False,  help='Whether adjust Pvalue: No for lognormal distribution, Yes for beta distribution')

# Adjust for low throughput tech, like starmap
parser.add_argument('--noExtraNullGenes', action='store_true', default=False, help='Default: add extra 1000* null genes if the input genes are fewer than extraNullGeneNumberThres')
parser.add_argument('--extraNullGeneNumberThres', type=int, default=50, help='if less than extra Null Gene Number, add null genes (default: %(default)s)')
parser.add_argument('--nullGeneNumber', type=int, default=1000, help='extra Null Gene Number (default: %(default)s)')
parser.add_argument('--ManyNullGenes', action='store_true', default=False, help='1000 or 1000* null genes (default: %(default)s)')
parser.add_argument('--nullDebug', action='store_true', default=False, help='Debug information for output null genes (default: %(default)s)')

# If the data has many isolated patches
parser.add_argument('--isolated', action='store_true', default=False, help='Default: if the input is very sparse and has many isolated patches')

# Scaling 
parser.add_argument('--noScaling', action='store_true', default=False, help='Default: scaling')

# debug
parser.add_argument('--memory', action='store_true', default=False, help='Output perform info in mem usage')
parser.add_argument('--debugperform', action='store_true', default=False, help='Output perform info in mem usage and computation time')

args = parser.parse_args()
scaleFactor = args.scaleFactor

if not os.path.exists(args.outputDir):
   os.makedirs(args.outputDir)

def check_param():
    '''Check valid parameters'''
    if args.D1>=args.D2:
        print("Radius of small patch D1 should be smaller than big patch D2, please check")
        sys.exit(0)

def debuginfoStr(info):
    '''
    Default: output memory and computational time
    '''
    print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time)))+'---'+info)
    if args.memory:
        if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        # linux or OSX
            import resource
            mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print('Mem consumption: '+str(mem))
        else:
            print('Non-linux system does not support memory consumption statistics')

def debugperformStr():
    '''
    For performance comparison usage, output memory and computational time
    '''
    if args.debugperform:
        print('Running time: '+str(time.time()-alg_time))
        if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        # linux or OSX
            import resource
            mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            print('Mem consumption: '+str(mem))
        else:
            print('Non-linux system does not support memory consumption statistics') 

def checkNullGenes(expMatrix):
    '''
    Add Null genes if the input gene number is so small, especially for low-throughput techs as FISH
    '''
    # print(args.noExtraNullGenes)
    # print(args.extraNullGeneNumberThres)
    addNullNum = 0
    if not args.noExtraNullGenes:
        if expMatrix.shape[0] < args.extraNullGeneNumberThres:
            print('Add 1000* null genes')
            nullMatrix = expMatrix.copy()
            # print(nullMatrix.shape)            
            # Option 1: (Not used) 1000 null genes
            if not args.ManyNullGenes:
                nullMatrix = np.tile(nullMatrix,(1 + (1000 // nullMatrix.shape[0]), 1))[:1000, :] 
            # Option 2: 1000*genes null genes
            else:
                nullMatrix = np.repeat(nullMatrix,args.nullGeneNumber,axis=0)
            addNullNum = nullMatrix.shape[0]
            # print(nullMatrix.shape)
            # print('Start Permutation')
            for i in range(len(nullMatrix)):
                nullMatrix[i] = np.random.permutation(nullMatrix[i])                    
            expMatrix = np.concatenate((expMatrix,nullMatrix),axis=0)
            # print(expMatrix.shape)
            # expMatrix: (originalGeneNumber+1000,locNumber)
    return expMatrix,addNullNum

def granp(spatialMatrix, expMatrix, D1 = 1.0, D2 = 2.5 ):
    '''calculate the p-values for a given array'''
    # calculate the statistics for a given radius D
    def targetstatisticsbyD(radiusD):
        # define the patch for each centerpoint
        def patch(CenterPoint):
           # calculate the distance between cells for each center cell
            CenterCoord = spatialMatrix[CenterPoint,:].reshape(1,-1)
            DistList = distance.cdist(CenterCoord, spatialMatrix, args.distType)
            ## Option 2: select kth nearest cells
            # DistListOrder = DistList.argsort()[0]
            # PatchCells = [DistListOrder[i] for i in np.arange(0,K+1) if DistList[0,DistListOrder[i]] <= radiusD]            
            ## Distance
            PatchCells =[i for i,v in enumerate(DistList[0]) if (v <= radiusD and v>0.0)]
            # Fix a bug, if the patchCells has 0 cells, then just include itself, to aviod NaN error.
            if len(PatchCells)==0:               
                PatchCells.append(CenterPoint)
            return(PatchCells)    
        
        # define the function to calculate variance of X_kpj for each gene j
        def patchvar():
            # define the function to calcualte X_kpj for each gene j and each cell p
            def patchmeans(CenterPoint):
                # Test for edge effects: May be delete
                if args.edgeFillingTag:
                    if radiusD == 1.0:
                        edgeFixNum = 4
                    elif radiusD == 4.3:
                        edgeFixNum = 48
                    if len(PatchesCells[CenterPoint])>=edgeFixNum:
                        X_kpj = np.mean(Norm_Exp[:,PatchesCells[CenterPoint]],axis=1)
                    else:
                        padding_num = edgeFixNum-len(PatchesCells[CenterPoint])
                        tmp_mat = np.mean(Norm_Exp,axis=1)
                        tmp_mat = tmp_mat.reshape(tmp_mat.shape[0],1)
                        ori_mat = Norm_Exp[:,PatchesCells[CenterPoint]]
                        for i in range(padding_num):                        
                            ori_mat = np.concatenate((ori_mat,tmp_mat),axis=1)
                        X_kpj = np.mean(ori_mat,axis=1)
                else:
                    X_kpj = np.mean(Norm_Exp[:,PatchesCells[CenterPoint]],axis=1)

                return(X_kpj)
                
            X_kj = np.zeros([spatialMatrix.shape[0], expMatrix.shape[0]])
            # calculate the variance
            for p in np.arange(spatialMatrix.shape[0]):
                X_kj[p,:] = patchmeans(p)  
            return(np.var(X_kj,axis=0))   
        debuginfoStr('Start Define Patches')
        # define the cell of patches
        PatchesCells = [patch(p) for p in np.arange(spatialMatrix.shape[0])] 
        debuginfoStr('Calculating Variances of Each Patch')
        
        # min max normalization
        if args.normalize == 'minmax':
            Norm_Exp = minmax_scale(expMatrix,axis=1)
        else:
            Norm_Exp = expMatrix
        var_X_K = patchvar()

        return(var_X_K)
    
    # p-value corrections
    def adjpvalue(p):
        p = np.asfarray(p)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / np.arange(float(len(p)), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]
    
    debuginfoStr('Start BSP Calculating')
       
    # calculate the summary statistics
    Var_X = [targetstatisticsbyD(K) for K in [D1] + [D2]]
    Var_X_0_Add = [np.var(expMatrix[i,:]) for i in range(expMatrix.shape[0])]
    Var_X_0_Add = (Var_X_0_Add/max(Var_X_0_Add))**scaleFactor
    Var_X = np.asarray(Var_X)
    # Calculate the test statistics:
    ## Debug for isolated data, replace the 0 value, but may takes some time
    if args.isolated:
        np.where(Var_X[0]==0,Var_X[0]+0.0001*np.random.rand(),Var_X[0])
    T_matrix_sum = Var_X[1]/Var_X[0]*Var_X_0_Add
    ## Debug for further exploration, some data return error
    if args.isolated:
        np.where(T_matrix_sum==0,T_matrix_sum+0.0001*np.random.rand(),T_matrix_sum)
 
    debuginfoStr('Start calculate p-value')
    # print(np.max(T_matrix_sum))
    # print(np.min(T_matrix_sum))
    # calculate p-values
    T_matrix_sum_upper90 = np.quantile(T_matrix_sum, 0.90)
    T_matrix_sum_mid = [i for i in T_matrix_sum if i < T_matrix_sum_upper90]
    # print(np.max(T_matrix_sum_mid))
    # print(np.min(T_matrix_sum_mid))
    
    ## Select the fitting distribution 
    if args.fitDist == 'lognormal':
        debuginfoStr('Start lognormal fit')
        LogNormPar = [np.mean(np.log(T_matrix_sum_mid)), np.std(np.log(T_matrix_sum_mid))]
        debuginfoStr('Start calculate pvalue')
        pvalues = [1 - stats.lognorm.cdf(i, scale=np.exp(LogNormPar[0]), s = LogNormPar[1]) for i in T_matrix_sum]    
    elif args.fitDist == 'beta':
        debuginfoStr('Start beta fit')
        BetaPar = stats.beta.fit(T_matrix_sum_mid, floc=0, fscale=1)
        a0 = BetaPar[0]
        b0 = BetaPar[1]
        debuginfoStr('Start calculate pvalue')
        pvalues = [1-stats.beta.cdf(i, a0, b0) for i in T_matrix_sum]    
    else:
        print('Please select appropriate fitting distribution: lognormal/beta')
        sys.exit()
    ## whether adjust pvalue, or not
    if args.adjustP:
        adjpvalues = adjpvalue(pvalues)
        debuginfoStr('End adjust pvalue')
    else:
        adjpvalues = pvalues
    return(adjpvalues)

if __name__ == "__main__":
    start_time = time.time()

    debuginfoStr('Start')
    check_param()
    scalFactor = 1.0 # Function f

    # if 3D: 
    # The example data comes from Vickovic, S., Schapiro, D., Carlberg, K. et al. Three-dimensional spatial transcriptomics uncovers cell type localizations in the human rheumatoid arthritis synovium. Commun Biol 5, 129 (2022). https://doi.org/10.1038/s42003-022-03050-3
    # https://github.com/mssanjavickovic/3dst
    if args.for3DTag:
        if not args.useDirTag:
            spaLocFile = args.inputDir + args.datasetName + '/' + args.spaLocFilename
            expFile = args.inputDir + args.datasetName + '/' + args.expFilename
            outputFile = args.outputDir + args.datasetName + '_' +args.outputFilename
            # data pre-processing
            debuginfoStr('Loading data and preprocessing')

            spatialMatrix = pd.read_csv(spaLocFile)
            spatialMatrix = spatialMatrix[['x','y','z']]
            spatialMatrix = spatialMatrix.to_numpy() #For MOB: (262,2)

            if args.noinputCellTag:
                expMatrix = pd.read_csv(expFile)    
            else:
                expMatrix = pd.read_csv(expFile,index_col=0)
                            
            print('Input dim: '+str(expMatrix.shape))
            # drop genes with all zeros
            expMatrix = expMatrix.loc[~expMatrix.apply(lambda row: (row==0.0).all(), axis=1)]
            # expMatrix = expMatrix.loc[:, expMatrix.any()]
            print('Nonzero dim: '+str(expMatrix.shape))
            if not args.notransTag:          
                expMatrix = expMatrix.transpose()

            # print(expMatrix)
            geneIndex = expMatrix.index
            expMatrix = expMatrix.to_numpy() 

            # for data with very few genes, add null genes from permutations
            expMatrix, addNullNum=checkNullGenes(expMatrix)

            # logTransform:
            if args.logTransform:
                expMatrix = np.log(expMatrix+1)

            # Scaling to normalize the spots in different scales: e.g.100 spots in 100*100: sqrt(100*100)/sqrt(100)=10
            if not args.noScaling:
                scalFactor = gmean(np.max(spatialMatrix,axis=0)-np.min(spatialMatrix,axis=0))/(spatialMatrix.shape[0])**(1/3)
                D1 = args.D1*scalFactor
                D2 = args.D2*scalFactor
                print("Scaled small patch radius:"+str(D1)+"\tScaled big patch radius:"+str(D2))

            alg_time = time.time()
            #============================ calculate granularity ==========================#
            P_values,Var_X, T_matrix_sum = granp(spatialMatrix, expMatrix, D1 = D1, D2 = D2)
            np.save(outputFile+'_score.npy',T_matrix_sum)
            np.save(outputFile+'_variance.npy',Var_X)
            debugperformStr()

            if not args.nullDebug:
                if addNullNum !=0 :
                    # 1000* null genes
                    if args.ManyNullGenes:
                        P_values=P_values[:-addNullNum]
                    # 1000 null genes
                    else:
                        P_values=P_values[:-1000]

            #================================= generate outputs ==========================#
            debuginfoStr('Post processing')
            if args.nullDebug:
                    outputData = pd.DataFrame(P_values, columns = ["p_values"])
            else:
                outputData = pd.DataFrame(P_values, 
                                        columns = ["p_values"],
                                        index = geneIndex)
            outputData.to_csv(outputFile)
            debuginfoStr('BSP Finished')        
        else:
            path = args.inputDir
            files = glob.glob(path+"*.csv")
            for filename in files:
                InputData =  pd.read_csv(filename)
                debuginfoStr('Loading data and preprocessing:'+filename)
                spatialMatrix = InputData[['x','y','z']]
                spatialMatrix = spatialMatrix.to_numpy()
                expMatrix = InputData.drop(['x','y','z'],  axis = 1)
                print('Input dim: '+str(expMatrix.shape))
                # drop genes with all zeros
                expMatrix = expMatrix.loc[~expMatrix.apply(lambda row: (row==0.0).all(), axis=1)]
                # expMatrix = expMatrix.loc[:, expMatrix.any()]
                print('Nonzero dim: '+str(expMatrix.shape))
                expMatrix = expMatrix.transpose()
                geneIndex = expMatrix.index
                expMatrix = expMatrix.to_numpy()

                # for data with very few genes, add null genes from permutations
                expMatrix, addNullNum=checkNullGenes(expMatrix)

                # logTransform:
                if args.logTransform:
                    expMatrix = np.log(expMatrix+1)
                
                # Scaling to normalize the spots in different scales: e.g.100 spots in 100*100: sqrt(100*100)/sqrt(100)=10
                if not args.noScaling:
                    scalFactor = gmean(np.max(spatialMatrix,axis=0)-np.min(spatialMatrix,axis=0))/(spatialMatrix.shape[0])**(1/3)
                    D1 = args.D1*scalFactor
                    D2 = args.D2*scalFactor
                    print("Scaled small patch radius:"+str(D1)+"\tScaled big patch radius:"+str(D2))

                #============================ calculate granularity ==========================#
                P_values,Var_X,T_matrix_sum = granp(spatialMatrix, expMatrix, D1 = D1, D2 = D2)

                if not args.nullDebug:
                    if addNullNum !=0 :
                        # 1000* null genes
                        if args.ManyNullGenes:
                            P_values=P_values[:-addNullNum]
                        # 1000 null genes
                        else:
                            P_values=P_values[:-1000]

                #================================= generate outputs ==========================#
                debuginfoStr('Post processing')
                if args.nullDebug:
                    outputData = pd.DataFrame(P_values, columns = ["p_values"])
                else:
                    outputData = pd.DataFrame(P_values, 
                                        columns = ["p_values"],
                                        index = geneIndex)
                outputData.to_csv(args.outputDir+filename.split('/')[-1].split('.csv')[0]+'_P_values.csv')
                np.save(args.outputDir+filename.split('/')[-1].split('.csv')[0]+'_P_values.csv'+'_score.npy',T_matrix_sum)
                np.save(args.outputDir+filename.split('/')[-1].split('.csv')[0]+'_P_values.csv'+'_variance.npy',Var_X)
            debuginfoStr('BSP Finished')

    # 2D
    else:
        if not args.useDirTag:
            spaLocFile = args.inputDir + args.datasetName + '/' + args.spaLocFilename
            expFile = args.inputDir + args.datasetName + '/' + args.expFilename
            outputFile = args.outputDir + args.datasetName + '_' +args.outputFilename

            # data pre-processing
            debuginfoStr('Loading data and preprocessing')

            # Simulation data only in one file
            if args.oneFileTag:
                InputData =  pd.read_csv(args.inputDir + args.datasetName + '/' + args.inputFile)
                spatialMatrix = InputData[['x','y']]
                spatialMatrix = spatialMatrix.to_numpy()
                expMatrix = InputData.drop(['x','y'],  axis = 1)
                expMatrix = expMatrix.transpose() 
            else:
                spatialMatrix = pd.read_csv(spaLocFile)
                spatialMatrix = spatialMatrix[['x','y']]
                spatialMatrix = spatialMatrix.to_numpy()

                if args.noinputCellTag:
                    expMatrix = pd.read_csv(expFile)    
                else:
                    expMatrix = pd.read_csv(expFile,index_col=0)    
                
                if not args.notransTag:
                    # drop genes with all zeros
                    expMatrix = expMatrix.loc[:, expMatrix.any()]
                    expMatrix = expMatrix.transpose()
            
            print('Input dim:'+str(expMatrix.shape)) #(Genes,locus) (12602,262) for MOB

            # filtering
            mat = expMatrix.to_numpy()
            pp=np.sum(mat,axis=1)
            ppuse = [i for i in range(mat.shape[0]) if pp[i] >= 15]
            expMatrix=expMatrix.iloc[ppuse]

            geneIndex = expMatrix.index
            expMatrix = expMatrix.to_numpy()
            print('Input dim after filter: '+str(expMatrix.shape)) #(Genes,locus) (12602,262) for MOB

            # for data with very few genes, add null genes from permutations
            expMatrix, addNullNum=checkNullGenes(expMatrix)

            # logTransform:
            if args.logTransform:
                expMatrix = np.log(expMatrix+1)
            
            # Scaling to normalize the spots in different scales: e.g.100 spots in 100*100: sqrt(100*100)/sqrt(100)=10
            if not args.noScaling:
                scalFactor = gmean(np.max(spatialMatrix,axis=0)-np.min(spatialMatrix,axis=0))/np.sqrt(spatialMatrix.shape[0])
                D1 = args.D1*scalFactor
                D2 = args.D2*scalFactor
                print("Scaled small patch radius:"+str(D1)+"\tScaled big patch radius:"+str(D2))

            alg_time = time.time()
            #============================ calculate granularity ==========================#
            P_values,Var_X,T_matrix_sum = granp(spatialMatrix, expMatrix, D1 = D1, D2 = D2)
            np.save(outputFile+'_score.npy',T_matrix_sum)
            np.save(outputFile+'_variance.npy',Var_X)
            debugperformStr()
                
            if not args.nullDebug:
                if addNullNum !=0 :
                    # print(len(P_values))
                    # 1000* null genes
                    if args.ManyNullGenes:
                        P_values=P_values[:-addNullNum]
                    # 1000 null genes
                    else:
                        P_values=P_values[:-1000]
                    # print(len(P_values))

            #================================= generate outputs ==========================#
            debuginfoStr('Post processing')
            if args.nullDebug:
                outputData = pd.DataFrame(P_values, columns = ["p_values"])
            else:
                outputData = pd.DataFrame(P_values, 
                                        columns = ["p_values"],
                                        index = geneIndex)
            outputData.to_csv(outputFile)
            debuginfoStr('BSP Finished')
        else:
            # Use files in the dir, for the simulation
            # path = args.inputDir + args.datasetName + '/'
            path = args.inputDir
            files = glob.glob(path+"*.csv")
            for filename in files:
                InputData =  pd.read_csv(filename)
                spatialMatrix = InputData[['x','y']]
                spatialMatrix = spatialMatrix.to_numpy()
                expMatrix = InputData.drop(['x','y'],  axis = 1)
                expMatrix = expMatrix.transpose()
                geneIndex = expMatrix.index
                expMatrix = expMatrix.to_numpy()

                # for data with very few genes, add null genes from permutations
                expMatrix, addNullNum=checkNullGenes(expMatrix)

                # logTransform:
                if args.logTransform:
                    expMatrix = np.log(expMatrix+1)

                # Scaling to normalize the spots in different scales: e.g.100 spots in 100*100: sqrt(100*100)/sqrt(100)=10
                if not args.noScaling:
                    scalFactor = gmean(np.max(spatialMatrix,axis=0)-np.min(spatialMatrix,axis=0))/np.sqrt(spatialMatrix.shape[0])
                    D1 = args.D1*scalFactor
                    D2 = args.D2*scalFactor
                    print("Scaled small patch radius:"+str(D1)+"\tScaled big patch radius:"+str(D2))

                #============================ calculate granularity ==========================#
                P_values,Var_X,T_matrix_sum = granp(spatialMatrix, expMatrix, D1 = D1, D2 = D2)

                if not args.nullDebug:
                    if addNullNum !=0 :
                        # 1000* null genes
                        if args.ManyNullGenes:
                            P_values=P_values[:-addNullNum]
                        # 1000 null genes
                        else:
                            P_values=P_values[:-1000]

                #================================= generate outputs ==========================#
                if args.nullDebug:
                    outputData = pd.DataFrame(P_values, columns = ["p_values"])
                else:
                    outputData = pd.DataFrame(P_values, 
                                            columns = ["p_values"],
                                            index = geneIndex)
                outputData.to_csv(args.outputDir+filename.split('/')[-1].split('.csv')[0]+'_P_values.csv')
                np.save(args.outputDir+filename.split('/')[-1].split('.csv')[0]+'_P_values.csv'+'_score.npy',T_matrix_sum)
                np.save(args.outputDir+filename.split('/')[-1].split('.csv')[0]+'_P_values.csv'+'_variance.npy',Var_X)
