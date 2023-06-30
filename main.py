import argparse
import glob
import sys
sys.path.append("../DataSet")
from ReadGraph import _Graph_
from ReadPCs import _pCohesion_
from ReadELV import _ELV_

from Triangle import *
########################################################################################################################
#################### global variables ####################
dataResultFolder = "../DataSet/"
global resultFolder
global figureFolder
global subgraphDataFolder
global subgraphFigureFolder
global noiseStr
global randomTime
global ResultsFolder
########################################################################################################################
def __DDPCountSubgraphsSub__(lock, checkset, **kwargs):
    varEpsilon1 = checkset[0]
    varEpsilon2 = checkset[1]
    hp = checkset[2]
    rt = checkset[3]

    PCSetTag = kwargs['PCSetTag']
    ELVSetTag = kwargs['ELVSetTag']
    subgraphType = kwargs['subgraphType']

    dictrl1 = {'varEpsilon1': varEpsilon1, 'varEpsilon2': varEpsilon2,
            'hPrime': hp, 'rt': rt, 'ResultsFolder': ResultsFolder}
    dictrl1.update(kwargs)

    paraStr = '_var1=' + str(varEpsilon1) + '_var2=' + str(varEpsilon2)
    if subgraphType in ['Triangle', 'kClique']:
        paraStr = '_var1=' + str(varEpsilon1) + '_var2=' + str(varEpsilon2) + '_hprime=' + str(hp)
    ######################################## For Triangle Counting ########################################
    if subgraphType == 'Triangle':
        ######################################## Files ########################################
        triFileCSV = subgraphDataFolder + 'Triangle' + paraStr + '.csv'
        triInfoFile = resultFolder + 'Triangle_Info.pickle'
        lambdaELVFile = ResultsFolder + 'Triangle_ELV_Lambda'
        if PCSetTag:
            triFileCSV = subgraphDataFolder + 'Triangle' + paraStr + '_PC-ELV.csv'
            triInfoFile = resultFolder + 'Triangle_PC-ELV_Info.pickle'
        ######################################## run ########################################
        runTag = checkRunTag(lock, triFileCSV, triInfoFile, rt)
        if runTag:
            ctrl = {'randomState': np.random.RandomState(), 'dataInfoFile': triInfoFile, 'lambdaELVFile':lambdaELVFile}
            ctrl.update(dictrl1)
            triangleObj = Triangle(**ctrl)
            MREValPC, upperBoundPC, MREValELV, upperBoundELV, lambdaPC, lambdaELV = triangleObj.twoPhaseApproachTriangle(
                                triFile = triFileCSV, hInput = hp, PCSetTag = PCSetTag, localViewSetTag = ELVSetTag,
                                dataFolder = resultFolder, parallelTag = parallelTag)
            del triangleObj
            dataInfo = [[varEpsilon1, varEpsilon2, hp, rt, round(kwargs['deltaVal'],2), MREValPC,
                         MREValELV, lambdaPC, lambdaELV, upperBoundPC, upperBoundELV]]
            writeCSV(lock, rt, dataSet=dataInfo, filename=triFileCSV, randomTime = randomTime)

    elif subgraphType == 'kClique':
        ######################################## Files ########################################
        triFolder = resultFolder + 'Triangle' + '_' + noiseStr + "_" + dpmechanism + "/";  mymkdir(triFolder)
        triFileCSV = triFolder + 'Triangle' + paraStr + '.csv'
        kCliqueCSV = subgraphDataFolder + 'kClique' + paraStr + '.csv'

        triInfoFile = resultFolder + 'Triangle_Info.pickle'
        kCliqueInfoFile = resultFolder + 'kClique_Info.pickle'
        lambdaELVFileTri = ResultsFolder + 'Triangle_ELV_Lambda'

        if PCSetTag:
            triFileCSV = triFolder + 'Triangle' + paraStr + '_PC-ELV.csv'
            kCliqueCSV = subgraphDataFolder + 'kClique' + paraStr + '_PC-ELV.csv'

            triInfoFile = resultFolder + 'Triangle_PC-ELV_Info.pickle'
            kCliqueInfoFile = resultFolder + 'kClique_PC-ELV_Info.pickle'
        lambdaELVFile = ResultsFolder + 'kClique_ELV_Lambda'

        ######################################## run ########################################
        runClique = checkRunTag(lock, kCliqueCSV, kCliqueInfoFile, rt)
        lambdaPC, lambdaELV = 1, 1
        randomState = np.random.RandomState()
        ########## triangle ##########
        if runClique:
            if not check_file_exists(triInfoFile, rmTag=False):
                runTriangle = True
            else:
                if not check_file_exists(triFileCSV, rmTag=False):
                    runTriangle = True
                else:
                    dataInfo = getDataInfo(lock, triFileCSV)
                    if rt not in dataInfo[:, 3]:
                        runTriangle = True
                    else:
                        runTriangle = False
                        for i in range(len(dataInfo)):
                            dVal = dataInfo[i, :]
                            if dVal[0] == varEpsilon2 and dVal[1] == varEpsilon2 and dVal[2] == hp and dVal[3] == rt:
                                lambdaPC = dVal[7]
                                lambdaELV = dVal[8]
                                break
        else:
            runTriangle = False
        ########## k clique ##########

        #################### Triangle ####################
        if runTriangle and runClique:
            ctrl = {'randomState': np.random.RandomState(), 'dataInfoFile': triInfoFile, 'lambdaELVFile':lambdaELVFileTri}
            ctrl.update(dictrl1)
            triangleObj = Triangle(**ctrl)

            MREValPCTri, UpperPCTri, MREValELVTri, UpperELVTri, lambdaPC, lambdaELV = triangleObj.twoPhaseApproachTriangle(
                triFile=triFileCSV, PCSetTag=PCSetTag, localViewSetTag=ELVSetTag, dataFolder=resultFolder,
                parallelTag=parallelTag)
            dataInfoTri = [[varEpsilon1, varEpsilon2, hp, rt, round(kwargs['deltaVal'], 2), MREValPCTri,
                            MREValELVTri, lambdaPC, lambdaELV, UpperPCTri,
                            UpperELVTri]]
            writeCSV(lock, rt, dataSet=dataInfoTri, filename=triFileCSV, randomTime=randomTime)

            # lambdaPC, lambdaELV = triangleObj.Phase1Triangle(hInput=hp, PCSetTag = PCSetTag, localViewSetTag = ELVSetTag,
            #                                                  dataFolder = resultFolder, parallelTag = parallelTag)
            del triangleObj
        maxPC = math.ceil(lambdaPC * varEpsilon2 / 3)
        maxELV = math.ceil(lambdaELV * varEpsilon2 /3)

        UpperPC = k * ncr(maxPC, k - 2)
        UpperELV = k * ncr(maxELV, k-2)

        lambdaPC = UpperPC / varEpsilon2
        lambdaELV = UpperELV / varEpsilon2
        #################### k-clique ####################
        if runClique:
            ctrl = {'randomState': np.random.RandomState(), 'dataInfoFile': kCliqueInfoFile, 'lambdaELVFile':lambdaELVFile}
            ctrl.update(dictrl1)
            kcliqueObj = kClique(**ctrl)

            MREValPC, upperBoundPC, MREValELV, upperBoundELV, lambdaPC, lambdaELV = kcliqueObj.twoPhaseApproachKCliquePath(
                                        kCliqueFile=kCliqueCSV, PCSetTag = PCSetTag, ELVSetTag = ELVSetTag,
                                        dataFolder = resultFolder, lambdaPC=lambdaPC, lambdaELV=lambdaELV, UpperPC=UpperPC, UpperELV=UpperELV)
            del kcliqueObj
            dataInfo = [[varEpsilon1, varEpsilon2, hp, rt, round(kwargs['deltaVal'], 2), MREValPC,
                         MREValELV, lambdaPC, lambdaELV, upperBoundPC, upperBoundELV]]
            writeCSV(lock, rt, dataSet=dataInfo, filename=kCliqueCSV, randomTime=randomTime)
    else:
        sys.exit('!!!!!!!!!!!!!!!!!!!! please give the correct subgraph type !!!!!!!!!!!!!!!!!!!!')

def __DPCountSubgraphsSubPre__(**kwargs):
    subgraphType = kwargs['subgraphType']
    varSepSetTime1 = kwargs['varSepSetTime1']
    hPrimeSet = kwargs['hPrimeSet']
    parallelTag = kwargs['parallelTag']
    ########################################
    checkSet = []
    checkFile = resultFolder + 'checkSets.pickle'
    if not check_file_exists(checkFile, rmTag = True):
        if subgraphType in ['Triangle', 'kClique']:
            for varSep in varSepSetTime1:
                var1, var2 = varSep[0], varSep[1]
                for hP in hPrimeSet:
                    for rt in range(randomTime):
                        checkSet.append([round(var1,2), round(var2,2), round(hP), rt])
        else:
            for varSep in varSepSetTime1:
                var1, var2 = varSep[0], varSep[1]
                for rt in range(randomTime):
                    checkSet.append([round(var1,2), round(var2,2), -1, rt])
        dataInfo = {'checkSet': checkSet}; savePickle(dataSet = dataInfo, filename = checkFile)
    else:
        with open(checkFile, 'rb') as fval:
            fdata = pickle.load(fval)
            checkSet = fdata['checkSet']
            del fdata
    ######################################## Parallel ########################################
    tryResult=True
    if parallelTag:
        try:
            lock = Manager().Lock()
            pool = Pool(processes = cpu_count())
            for checkset in checkSet:
                pool.apply_async(__DDPCountSubgraphsSub__, (lock, checkset,), kwargs)
            pool.close(); pool.join()
        except:
            tryResult = False
    if not parallelTag or not tryResult:
        lock = Manager().Lock()
        for checkset in checkSet:
            __DDPCountSubgraphsSub__(lock, checkset, **kwargs)
    ######################################## Plot Results ########################################


def __plotMultiVar1__(p, dataSetName, ReadPara, ObjParSub, subgraphFigureFolder, varTimesSet1, folderName = 'multiTimes0.1', picklefile = None):
    subfigureFolderMulti1 = subgraphFigureFolder + folderName + '/'; mymkdir(subfigureFolderMulti1)
    print("========== Run For {} ==========".format(folderName))
    ObjPar = {'figureFolder': subfigureFolderMulti1}
    ObjPar.update(ObjParSub)

    plotObj = DuplexPlot(varTimesSet1, **ObjPar)
    if not check_file_exists(picklefile):
        plotObj.ReadFiles(**ReadPara)
        dataInfo = {'plotMREPCSet': plotObj.plotMREPCSet, 'plotMREPCSDSet': plotObj.plotMREPCSDSet,
                    'plotMREELVSet': plotObj.plotMREELVSet, 'plotMREELVSDSet': plotObj.plotMREELVSDSet,
                    'plotUpperPCSet': plotObj.plotUpperPCSet, 'plotUpperPCSDSet': plotObj.plotUpperPCSDSet,
                    'plotUpperELVSet': plotObj.plotUpperELVSet, 'plotUpperELVSDSet': plotObj.plotUpperELVSDSet,
                    'plotLambdaPCSet': plotObj.plotLambdaPCSet, 'plotLambdaPCSDSet': plotObj.plotLambdaPCSDSet,
                    'plotLambdaELVSet': plotObj.plotLambdaELVSet, 'plotLambdaELVSDSet': plotObj.plotLambdaELVSDSet,
                    'var1SetTimes': plotObj.var1SetTimes, 'hPrimeSet': plotObj.hPrimeSet}
        savePickle(dataSet = dataInfo, filename = picklefile)
    else:
        with open(picklefile, 'rb') as fval:
            fdata = pickle.load(fval)
            plotObj.plotMREPCSet = fdata['plotMREPCSet']
            plotObj.plotMREPCSDSet = fdata['plotMREPCSDSet']
            plotObj.plotMREELVSet = fdata['plotMREELVSet']
            plotObj.plotMREELVSDSet = fdata['plotMREELVSDSet']
            plotObj.plotUpperPCSet = fdata['plotUpperPCSet']
            plotObj.plotUpperPCSDSet = fdata['plotUpperPCSDSet']
            plotObj.plotUpperELVSet = fdata['plotUpperELVSet']
            plotObj.plotUpperELVSDSet = fdata['plotUpperELVSDSet']

            plotObj.plotLambdaPCSet = fdata['plotLambdaPCSet']
            plotObj.plotLambdaPCSDSet = fdata['plotLambdaPCSDSet']
            plotObj.plotLambdaELVSet = fdata['plotLambdaELVSet']
            plotObj.plotLambdaELVSDSet = fdata['plotLambdaELVSDSet']

    plotObj.DuplexPlots(p, dataSetName)

def __DDPCountSubgraphs__(p, dataSetName, verDeg, neiSet, PCSetTag = False, ELVSetTag = False,
                          deltaVal = 1e-10, subgraphType = "Triangle", k = 4,
                          noiseCal = True, dpmechanism = 'Laplace', openCheck = False, parallelTag = False,
                          InitCandidates = 'Partial'):
    # Optimized Two-Phase Approach for Triangle counting
    print("========================================================================================================================")
    print("==================== Running for {} ====================".format(subgraphType))
    ######################################## variables ########################################
    global randomTime
    if platform.system() in ['Windows', 'Darwin']:
        marVar = 3; maxhPrime = 5; randomTime = 10
        varTimesSet1, varTimesSet2, varTimesSet3, varTimesSet4, varTimesSet5, varTimesSet6 = [],[],[],[],[],[]
    else:
        marVar = 20; maxhPrime = 30; randomTime = 3000
        varTimesSet1 = [[round(i, 1), round(1.0 - i, 1)] for i in np.arange(0.1, 1.0, 0.1)]
        varTimesSet2 = [[round(i, 1), round(2.0 - i, 1)] for i in np.arange(0.1, 2.0, 0.1)]
        varTimesSet3 = [[round(i, 1), round(3.0 - i, 1)] for i in np.arange(0.1, 3.0, 0.1)]
        varTimesSet4 = [[round(i, 1), round(4.0 - i, 1)] for i in np.arange(0.1, 4.0, 0.1)]
        varTimesSet5 = [[round(i, 1), round(5.0 - i, 1)] for i in np.arange(0.1, 5.0, 0.1)]
        varTimesSet6 = [[round(i, 1), round(10.0 - i, 1)] for i in np.arange(0.1, 10.0, 0.1)]

    if subgraphType in ['Triangle', 'kClique']: hPrimeSet = list(np.arange(1, maxhPrime+1, 1))
    else: hPrimeSet = [0]
    ########################################
    varSet = [i for i in np.arange(1, marVar+1, 1)]

    varSepSetTime1 = [[0.1*var, 0.9*var] for var in varSet]
    varSepSetTime1.extend(varTimesSet1)
    varSepSetTime1.extend(varTimesSet2)
    varSepSetTime1.extend(varTimesSet3)
    varSepSetTime1.extend(varTimesSet4)
    varSepSetTime1.extend(varTimesSet5)
    varSepSetTime1.extend(varTimesSet6)

    ctrlInfo = {'verDeg': verDeg, 'neiSet': neiSet, 'PCSetTag': PCSetTag, 'ELVSetTag': ELVSetTag, 'deltaVal': deltaVal,
                'subgraphType': subgraphType, 'k': k, 'noiseCal': noiseCal, 'dpmechanism': dpmechanism,
                'openCheck': openCheck, 'parallelTag': parallelTag, 'InitCandidates': InitCandidates,
                'varSepSetTime1': varSepSetTime1, 'hPrimeSet': hPrimeSet}
    __DPCountSubgraphsSubPre__(**ctrlInfo)

######################################## For 3-Path Counting ########################################
########################################################################################################################
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    ########################################################################################################################
    parser = argparse.ArgumentParser(description='Process some parameters.')
    parser.add_argument('--dataSetName', type = str, default = "HepPh", help='The data to process')
    ########## for p-cohesion ##########
    parser.add_argument('--p', type=float, default = 0.5, help ='The fraction p of p-Cohesion')
    parser.add_argument('--initCand', type=str, default = "Partial", choices=['All', 'ALL', 'all', 'Partial', 'partial', 'PARTIAL'], help='Choose all/Partial neighours of the query vertex')
    ########## for DDP ##########
    parser.add_argument('--var1', type=float, default=1, help='VarEpsilon 1')
    parser.add_argument('--var2', type=float, default=9, help='VarEpsilon 2')
    parser.add_argument('--delt', type=float, default=0, help='Delta')
    parser.add_argument('--hprime', type=int, default=100, help='h Prime')
    parser.add_argument('--k', type=int, default=4, help='kClique')
    parser.add_argument('--subgraph', type=str, default='Triangle', choices = ['Triangle', 'kClique'], help='The case study of subgraph counting')
    parser.add_argument('--noiseCal', type=bool, default=False, help='calibrate noise or not')
    parser.add_argument('--dpm', type=str, default='Laplace', help='DP-Mechanisms')
    ########## Public use ##########
    parser.add_argument('--openCheck', type=bool, default=False, help='Check the medium results')
    parser.add_argument('--parTag', type=bool, default=True, help='Parallel')
    ########################################################################################################################
    ########## Start to read the input parameters  ##########
    args = parser.parse_args()
    ########## Graph Info  ##########
    dataSetName = args.dataSetName
    ########## P-Cohesion  ##########
    p = args.p
    InitCandMethod = args.initCand
    ########## Two Phase DDP  ##########
    var1 = args.var1
    var2 = args.var2
    delt = args.delt
    hprime = args.hprime
    k = args.k
    subgraphType = args.subgraph
    noiseCal = args.noiseCal
    dpmechanism = args.dpm
    ########## Public Use ##########
    openCheck = args.openCheck
    parallelTag = args.parTag

    ########################################################################################################################
    warnings.simplefilter('error', UserWarning)
    warnings.simplefilter("ignore")
    global noiseStr
    if noiseCal: noiseStr = 'InterNoiseCalibrated'
    else: noiseStr = ''
    global ResultsFolder
    ResultsFolder = dataResultFolder + dataSetName + "/Result/"; mymkdir(ResultsFolder)
    global resultFolder
    resultFolder = ResultsFolder + "p" + str(p) + "/"; mymkdir(resultFolder)
    global subgraphDataFolder
    subgraphDataFolder = resultFolder + subgraphType + '_' + noiseStr + "_" + dpmechanism + "/"; mymkdir(subgraphDataFolder)
    global figureFolder
    figureFolder = ResultsFolder + "figures_p" + str(p) + "/"; mymkdir(figureFolder)
    global subgraphFigureFolder
    subgraphFigureFolder = figureFolder + subgraphType + '_' + noiseStr + "_" + dpmechanism + '_3000/'; mymkdir(subgraphFigureFolder)

    if platform.system() == 'Darwin': import pickle5 as pickle
    else: import pickle
    ######################################## Read graph information ########################################
    verDeg, neiSet, verOriTags = _Graph_(dataSetName, dataResultFolder)
    ########################################################################################################################
    ######################################## p-cohesion computation ########################################
    PCSet = _pCohesion_(p, verDeg, neiSet, verOriTags = verOriTags, InitCandidates = InitCandMethod, openCheck = openCheck, parallelTag = parallelTag, resultFolder = resultFolder)
    localViewSet = _ELV_(verDeg, neiSet, verOriTags = verOriTags, neiHop = 2, openCheck = openCheck, parallelTag = parallelTag, ResultsFolder = ResultsFolder)

    ######################################################su##################################################################
    ######################################## Optimized Two-Phase Approach ########################################
    if delt == 0: delt = 1 / len(verDeg)
    paraInfo = {'PCSetTag': True, 'ELVSetTag': True, 'deltaVal':delt,
                'subgraphType': subgraphType, 'k': k, 'noiseCal': noiseCal, 'dpmechanism': dpmechanism,
                'openCheck': openCheck, 'parallelTag': parallelTag, 'InitCandidates': InitCandMethod}
    __DDPCountSubgraphs__(p, dataSetName, verDeg, neiSet, **paraInfo)
    ########################################################################################################################
