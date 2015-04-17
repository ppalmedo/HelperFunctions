'''
helper Func: Created July 2014
    helper Func contains helper functions for Python.  Developed primarily for bioinformatic purposes.

@author: ppalmedo
'''

import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class myData:
    'Common class for RNAseq Data'
    
    def __init__(self,GeneList,CleanSampIDs,Data,Type):
        self.geneList = GeneList
        self.cleanSampIDs = CleanSampIDs
        self.data = Data
        self.type = Type
        
class combinedData:
    
    def __init__(self,GeneList,CleanSampIDs,rnaData,cnData):
        self.geneList = GeneList
        self.cleanSampIDs = CleanSampIDs
        self.rnaData = rnaData
        self.cnData = cnData

# function to take two lists (i.e. sample or genes) and return indices for them
def findOLIndices(listA,listB):
    indexList = [(listA.index(g), listB.index(g)) for g in listA if g in listB]
    if len(indexList) == 0:
        return [],[]
    iA, iB = map(list, zip(*indexList))
    return iA,iB


# function to combine tumor types
def combineTT(data1,data2):
    iA,iB = findOLIndices(data1.geneList,data2.geneList)
    newGeneList = [data1.geneList[i] for i in iA]
    newSampleList = data1.cleanSampIDs + data2.cleanSampIDs
    newCnData = np.column_stack((data1.cnData[iA,:],data2.cnData[iB,:]))
    newRnaData = np.column_stack((ss.zscore(data1.rnaData[iA,:],1),ss.zscore(data2.rnaData[iB,:],1)))
    newData = combinedData(newGeneList,newSampleList,newRnaData,newCnData)
    return newData


# function to get the data for rna and cn data sets
def getOLdata(data1,data2):
    if data1.type == 'rna' and data2.type == 'cn':
        #find overlapping samples
        numUnqRnaSamp = len(set(data1.cleanSampIDs))
        numUnqCnSamp = len(set(data2.cleanSampIDs))
        print "Finding overlapping samples..."
        rnaSampI,cnSampI = findOLIndices(data1.cleanSampIDs,data2.cleanSampIDs)
        if len(rnaSampI)>0:
            print "There are %s overlapping samples of %s for RNAseq and %s for CN..." % (len(rnaSampI), numUnqRnaSamp, len(data2.cleanSampIDs))
        else:
            print "ERROR: No overlapping samples.  Check Sample IDs!"
            return
        rnaSampI = np.array(rnaSampI)
        cnSampI = np.array(cnSampI)
        
        #find overlapping genes
        print "Finding overlapping genes..."
        rnaGeneI,cnGeneI = findOLIndices(data1.geneList,data2.geneList)
        if len(rnaGeneI)>0:
            print "There are %s overlapping genes of %s for RNAseq and %s for CN..." % (len(rnaGeneI), len(data1.geneList),len(data2.geneList))
        else:
            print "ERROR: No overlapping genes.  Check Gene Name Format!"
            return
        rnaGeneI = np.transpose(np.array(rnaGeneI))
        cnGeneI = np.transpose(np.array(cnGeneI))
        
        #select relevant data given overlap (won't reshape both at once?)
        OLrnaData = np.array(data1.data[rnaGeneI,:])
        OLcnData = data2.data[cnGeneI,:]
        OLrnaData = OLrnaData[:,rnaSampI]
        OLcnData = OLcnData[:,cnSampI]
        OLgeneList = [data1.geneList[rnaID] for rnaID in rnaGeneI]
        OLsampList = [data1.cleanSampIDs[rnaID] for rnaID in rnaSampI]
        
        return combinedData(OLgeneList,OLsampList,OLrnaData,OLcnData)
    else:
        print 'Error: RNA data must be first and CN data second'
        return


def getGeneOfInterest(nameGOI,geneList):
    print 'Finding gene of interest...'
    inData = 1
    if nameGOI in geneList:
        Goi = geneList.index(nameGOI)
    else:
        print "ERROR: Gene of Interest not in Overlapping Data"

    return Goi,inData


def loadRNAseqv2(GEfileName):
    # load gene expression (expression data from TCGA in \t delimited form)
    print 'Reading expression data from file...'
    rnaData = []
    rnaGeneList = []
    with open(GEfileName,'r') as theFile:
        rnaSampIDs = theFile.readline()
        rnaSampIDs = rnaSampIDs.split('\t')
        rnaSampIDs = rnaSampIDs[1:]
        rnaSignals = theFile.readline()
        rnaSignals = rnaSignals.split('\t')
        rnaSignals = rnaSignals[1:]
        for line in theFile:
            thisLine = line.split('\t')
            rnaGeneList.append(thisLine[0])
            rnaData.append([float(s) for s in thisLine[1:]])
        rnaGeneList = [gID.split('|')[0] for gID in rnaGeneList]
    
    rnaData = np.array(rnaData)
    rnaCleanSampIDs = [sampId[0:16] for sampId in rnaSampIDs]
    
    # get rid of genes with std < 1, because these will return nan for correlation
    dataGene = [(x,y) for x,y in zip(rnaGeneList,rnaData) if np.std(y)>1]
    rnaGeneList, rnaData = map(list, zip(*dataGene))
    rnaData = np.log(np.array(rnaData)+1)
    #rnaData = ss.zscore(rnaData,1)
    
    # print warning if non-unique samples
    numUnqRnaSamp = len(set(rnaCleanSampIDs))
    if numUnqRnaSamp < len(rnaCleanSampIDs):
        print 'WARNING: Replicates in expression data may inflate correlation!'
    
    print 'Finished reading expression data from file...'
    
    return rnaData,rnaGeneList,rnaCleanSampIDs


def createStratifiedFigure(gCNdata,tCNdata,gRnaData,tRnaData,nameGOI,newGeneName):
    # create necessary variables
    tCnDist = ss.itemfreq(tCNdata)
    tCnLevels = tCnDist[:,0]
    tCnCounts = tCnDist[:,1]
    tMeanCn = np.mean(tCNdata)
    tStdCn = np.std(tCNdata)
    tMeanRnaStrat = [np.mean(tRnaData[tCNdata == i]) for i in tCnLevels]
    tStdRnaStrat = [np.std(tRnaData[tCNdata == tCnLevels[i]])/np.sqrt(tCnCounts[i]) for i in range(len(tCnLevels))]
    
    gCnDist = ss.itemfreq(gCNdata)
    gCnLevels = gCnDist[:,0]
    gCnCounts = gCnDist[:,1]
    gMeanCn = np.mean(gCNdata)
    gStdCn = np.std(gCNdata)
    gMeanRnaStrat = [np.mean(tRnaData[gCNdata == i]) for i in gCnLevels]
    gStdRnaStrat = [np.std(tRnaData[gCNdata == gCnLevels[i]])/np.sqrt(gCnCounts[i]) for i in range(len(gCnLevels))]
    
    # create figure with 4 subplots
    fig = plt.figure(figsize=(13,7),dpi=100)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    
    # plot tRnaData vs tCnData
    ax1.scatter(tCNdata, tRnaData, s = 100, alpha = 0.3)
    ax1.errorbar(tCnLevels,tMeanRnaStrat,tStdRnaStrat,marker='_',markersize=40,markeredgewidth=2,color='r',elinewidth=3,capsize=10)
    ax1.errorbar(tMeanCn,min(tRnaData),xerr=tStdCn,marker='^',markersize=10,elinewidth=3,color='k')
    ax1.set_xlabel('%s CN' %newGeneName)
    ax1.set_ylabel('%s Expression' %newGeneName)
    
    # plot tRnaData vs gCnData
    ax2.scatter(gCNdata, tRnaData, s=100, alpha = 0.3)
    ax2.errorbar(gCnLevels,gMeanRnaStrat,gStdRnaStrat,marker='_',markersize=40,markeredgewidth=2,color='r',elinewidth=3,capsize=10)
    ax2.errorbar(gMeanCn,min(tRnaData),xerr=gStdCn,marker='^',markersize=10,elinewidth=3,color='k')
    ax2.set_xlabel('%s CN' %nameGOI)
    ax2.set_ylabel('%s Expression' %newGeneName)
    
    # plot the two stratifications
    plotStratification(tCNdata,gCNdata,tRnaData,newGeneName,nameGOI,ax3,1)
    plotStratification(gCNdata,tCNdata,tRnaData,newGeneName,nameGOI,ax4,2)
    ax3.set_xlabel('Big Separation: %s CN' %nameGOI)
    ax4.set_xlabel('Big Separation: %s CN' %newGeneName)
    mark1 = plt.Line2D((0,0),(0,0),color = 'k',marker = '<',linestyle = '')
    mark2 = plt.Line2D((0,0),(0,0),color = 'k',marker = 'v',linestyle = '')
    mark3 = plt.Line2D((0,0),(0,0),color = 'k',marker = 'o',linestyle = '')
    mark4 = plt.Line2D((0,0),(0,0),color = 'k',marker = '^',linestyle = '')
    mark5 = plt.Line2D((0,0),(0,0),color = 'k',marker = '>',linestyle = '')
    mark6 = mpatches.Patch(color = 'c')
    mark7 = mpatches.Patch(color = 'b')
    mark8 = mpatches.Patch(color = 'g')
    mark9 = mpatches.Patch(color = 'r')
    mark10 = mpatches.Patch(color = 'm')
    
    legLab1 = [nameGOI + ' -2+',nameGOI + ' -1',nameGOI +' diploid',nameGOI + ' +1',nameGOI + ' +2+']
    legLab2 = [newGeneName + ' -2+',newGeneName + ' -1',newGeneName +' diploid',newGeneName + ' +1',newGeneName + ' +2+']
    legLab3 = [newGeneName + ' -2+',newGeneName + ' -1',newGeneName +' diploid',newGeneName + ' +1',newGeneName + ' +2+']
    legLab4 = [nameGOI + ' -2+',nameGOI + ' -1',nameGOI +' diploid',nameGOI + ' +1',nameGOI + ' +2+']
    legend3 = ax3.legend([mark1,mark2,mark3,mark4,mark5,mark6,mark7,mark8,mark9,mark10],legLab1+legLab3,numpoints=1,loc='upper left',prop={'size':6})
    legend4 = ax4.legend([mark1,mark2,mark3,mark4,mark5,mark6,mark7,mark8,mark9,mark10],legLab1+legLab3,numpoints=1,loc='upper left',prop={'size':6})
    
def plotStratifiedScatter(gCNdata,tCNdata,gRnaData,tRnaData,nameGOI,newGeneName):  
    # color, shape, and alpha schemes for the stratification
    myColorScheme = ['c','b','g','r','m']*5
    myShapeScheme = ['<']*5+['v']*5+['o']*5+['^']*5+['>']*5
    myAlphaScheme = [0.01]*5+[0.01]*5+[0.01]*5+[0.3]*5+[0.01]*5
    
    fig = plt.figure(figsize=(13,6),dpi=100)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    theCorr = ss.pearsonr(tRnaData,gRnaData)[0]
    theCorr = "{:.2f}".format(theCorr)
    
    sumCNdata = 10*tCNdata + 2*gCNdata
    theColorCn = sumCNdata + 24
    colorDist = ss.itemfreq(theColorCn)
    colorDist = colorDist[:,0]
    for level in colorDist:
        thisIndex = level/2
        ax1.scatter(gRnaData[theColorCn == level],tRnaData[theColorCn == level],s=100,alpha=myAlphaScheme[int(thisIndex)],color=myColorScheme[int(thisIndex)],marker=myShapeScheme[int(thisIndex)])
    ax1.set_xlabel(nameGOI + ' RNA Expression')
    ax1.set_ylabel(newGeneName + ' RNA Expression')
    ax1.set_title('Expression Correlation: ' + str(theCorr))
    
    cnPair = zip(gCNdata,tCNdata)
    listVals = [-2,-1,0,1,2]
    potentMatch = [(i,j) for i in listVals for j in listVals]
    matchCount = [cnPair.count(r) for r in potentMatch]
    
    theGcn, theTcn = map(list, zip(*potentMatch))
    for i in range(len(theGcn)):
        ax2.scatter(theGcn[i],theTcn[i], s = 10*matchCount[i],alpha=0.5)
        if matchCount[i] > 0:
            ax2.text(theGcn[i],theTcn[i], '%i' %matchCount[i])
    ax2.set_xticks(np.arange(-2.5,3,0.5))
    ax2.set_yticks(np.arange(-2.5,3,0.5))
    ax2.set_xlabel(nameGOI + ' CN Level')
    ax2.set_ylabel(newGeneName+ ' CN Level')
    
    
    
def plotStratification(gCNdata,tCNdata,tRnaData,newGeneName,nameGOI,theAxes,strat):
    # color, shape, and alpha schemes for the stratification
    if strat == 1:
        myColorScheme = ['c','b','g','r','m']*5
        myShapeScheme = ['<']*5+['v']*5+['o']*5+['^']*5+['>']*5
        myAlphaScheme = [0.1]*5+[0.1]*5+[0.1]*5+[0.1]*5+[0.1]*5
    elif strat == 2:
        myColorScheme = ['c']*5 + ['b']*5 + ['g']*5 + ['r']*5 + ['m']*5
        myShapeScheme = ['<','v','o','^','>']*5
        myAlphaScheme = [0.1]*5+[0.1]*5+[0.1]*5+[0.1]*5+[0.1]*5
    
    sumCNdata = 10*tCNdata + 2*gCNdata
    theColorCn = sumCNdata + 24
    colorDist = ss.itemfreq(theColorCn)
    colorDist = colorDist[:,0]
    for level in colorDist:
        thisIndex = level/2
        theAxes.scatter(sumCNdata[theColorCn == level],tRnaData[theColorCn == level],s=100,alpha=0.3,color=myColorScheme[int(thisIndex)],marker=myShapeScheme[int(thisIndex)])
    sumCnDist = ss.itemfreq(sumCNdata)
    sumCnLevels = sumCnDist[:,0]
    sumCnCounts = sumCnDist[:,1]
    meanCN = np.mean(gCNdata)
    stdCN = np.std(gCNdata)
    tMeanCN = np.mean(tCNdata)
    tStdCN = np.std(tCNdata)
    sumMeanExpT = [np.mean(tRnaData[sumCNdata == i]) for i in sumCnLevels]
    sumStdExpT = [np.std(tRnaData[sumCNdata == sumCnLevels[i]])/np.sqrt(sumCnCounts[i]) for i in range(len(sumCnLevels))]
    theAxes.errorbar(sumCnLevels,sumMeanExpT,sumStdExpT,marker='_',markersize=15,markeredgewidth=2,color='k',elinewidth=3,capsize=4)
    theAxes.errorbar(10*tMeanCN,min(tRnaData),xerr=10*tStdCN,marker='^',markersize=10,elinewidth=3,color='k')
    theAxes.errorbar(10*meanCN,min(tRnaData)-0.1,xerr=10*stdCN,marker='^',markersize=10,elinewidth=3,color='k')
    theAxes.set_xticks(np.arange(-25,25,5))
    theAxes.grid()
    theAxes.set_ylabel('RNA Expression of %s' %newGeneName)
    

def loadCNdata(CNfileName):
    # load copy number data (output of GISTIC in gene form)
    print 'Reading copy-number data from file...'
    cnGeneList = []
    cnGeneProp = []
    cnData = []
    with open(CNfileName, 'r') as theFile:
        cnLabels = theFile.readline()
        cnLabels = cnLabels.split('\t')
        cnGeneLabels = cnLabels[0:3]
        cnSampIDs = cnLabels[3:]
        for line in theFile:
            thisLine = line.split('\t')
            cnGeneList.append(thisLine[0])
            cnGeneProp.append(thisLine[1:3])
            cnData.append([float(s) for s in thisLine[3:]])
    cnData = np.array(cnData)
    cnCleanSampIDs = [sampId[0:16] for sampId in cnSampIDs]
    
    numUnqCnSamp = len(set(cnCleanSampIDs))
    if numUnqCnSamp < len(cnCleanSampIDs):
        print 'WARNING: Replicates in CN data may inflate correlation!'
        
    print 'Finished reading copy-number data from file...'
    
    return cnData,cnCleanSampIDs,cnGeneList

def returnSigGenes(gID,cnSigVal,rnaSigVal,data):
    geneList = []
    geneListNeg = []
    
    for cType in data:
        regCn,regRna = runRegressionModel(gID,cType.cnData,cType.rnaData)
        
        betaHatCn, t_statCn = map(list, zip(*regCn))
        tfTstatCn = [w[1] for w in t_statCn]
        
        betaHatRna, t_statRna = map(list, zip(*regRna))
        tfTstatRna = [w[1] for w in t_statRna]
        
        theList = [(geneID,tValCn,tValRna) for geneID,tValCn,tValRna in zip(cType.geneList,tfTstatCn,tfTstatRna) if tValCn > cnSigVal and tValRna > rnaSigVal]
        theListNeg = [(geneID,tValCn,tValRna) for geneID,tValCn,tValRna in zip(cType.geneList,tfTstatCn,tfTstatRna) if tValCn < -cnSigVal and tValRna < -rnaSigVal]
        
        geneList.append(set(theList))
        geneListNeg.append(set(theListNeg))
        print theList
        
    print list(set(geneList[1]) & set(geneList[2]))
        

def runRegressionModel(gID,cnData,rnaData):
    print 'Running regressions...'
    numGenes = len(cnData[:,0])
    cnReg = [algebraRegress(cnData[gID,:],cnData[i,:],rnaData[i,:]) for i in range(numGenes)]
    rnaReg = [algebraRegress(rnaData[gID,:],cnData[i,:],rnaData[i,:]) for i in range(numGenes)]
    return cnReg, rnaReg

def runReverseRegressionModel(gID,cnData,rnaData):
    print 'Running regressions...'
    numGenes = len(cnData[:,0])
    cnReg = [algebraRegress(cnData[i,:],cnData[gID,:],rnaData[gID,:]) for i in range(numGenes)]
    rnaReg = [algebraRegress(rnaData[i,:],cnData[gID,:],rnaData[gID,:]) for i in range(numGenes)]
    return cnReg, rnaReg

# function to calculate the multivariate regression for two variables
def algebraRegress(v1,v2,responseV):
    # test to make sure not totally correlated - return 0 if so
    theCorr = ss.spearmanr(v1,v2)
    if np.abs(theCorr)[0] == 1:
        betaHat = np.zeros(3)
        t_stat =  np.zeros(3)
        return (betaHat,t_stat)
    
    # calculate linear regression
    n =responseV.size
    
    # zscore the response data
    #v1 = ss.zscore(v1)
    #v2 = ss.zscore(v2)
    #responseV = ss.zscore(responseV)
    
    # add column of ones for intercept
    inputV = np.column_stack((v1,v2))
    inputV = np.column_stack((np.ones(shape=(n,1)),inputV))
    
    # perform the regression
    iVprime = np.transpose(inputV)
    iVpiV = iVprime.dot(inputV)
    iVpiVT = np.linalg.inv(iVpiV)
    iVpiVTiV = iVpiVT.dot(iVprime)
    betaHat = iVpiVTiV.dot(responseV)
    
    # calculate standard error
    residV = (responseV - inputV.dot(betaHat))
    sigHatSq = residV.dot(residV) / (n - 3)
    varMat = sigHatSq * (iVpiVT)
    t_stat = np.array([betaHat[i]/np.sqrt(varMat[i,i]) for i in range(3)])
    
    return (betaHat,t_stat)