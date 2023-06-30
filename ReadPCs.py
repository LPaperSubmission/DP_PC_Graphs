from pCohesion_Old import *
from pCohesion import *
########################################################################################################################
def _pCohesion_(p, verDeg, neiSet, verOriTags = None, InitCandidates = 'All', openCheck = False, parallelTag = False, pcChoose = 'new', resultFolder = None):
    print("========================================================================================================================")
    #################### Start to generate p-cohesion for all node ####################
    if pcChoose in ['old', 'OLD', 'Old']:
        PCFile = resultFolder + 'PC_' + InitCandidates + '_Nei_p' + str(p) + 'Old.pickle'
        pcFolder = resultFolder + 'PCsOld/'; mymkdir(pcFolder)
    else:
        PCFile = resultFolder + 'PC_' + InitCandidates + '_Nei_p' + str(p) + '.pickle'
        pcFolder = resultFolder + 'PCs/'; mymkdir(pcFolder)

    pcSets = [[] for i in range(len(verDeg))]
    maxPCSize = 0
    minPCSize = len(verDeg)
    meanPCSize = 0.0

    if not check_file_exists(PCFile, rmTag=False):
        ########################################  generate p-cohesion for each node ########################################
        pcObj = pCohesion(verDeg, neiSet, p, openCheck)
        if pcChoose in ['old', 'OLD', 'Old']:
            pcObj = pCohesionOld(verDeg, neiSet, p, openCheck)
        tryResult = True
        if parallelTag:
            try:
                results = []
                pool = Pool(processes=cpu_count())
                for queryID in range(len(verDeg)):
                    queryPCFile = pcFolder + 'PC_' + str(queryID) + '_' + InitCandidates + '.pickle'
                    if not check_file_exists(queryPCFile):
                        results.append(pool.apply_async(pcObj.compute_pCohesion, (queryID,), {'candMethod': InitCandidates, 'queryPCFile': queryPCFile}))
                pool.close()
                pool.join()
                computedTag = [0] * len(verDeg)
                for res in results:
                    queryID, pcSet = res.get()
                    pcSets[int(queryID)] = pcSet; computedTag[int(queryID)] = 1
                    if len(pcSet) > maxPCSize: maxPCSize = len(pcSet)
                    if len(pcSet) < minPCSize: minPCSize = len(pcSet)
                    meanPCSize += len(pcSet)
                for i in range(len(verDeg)):
                    queryPCFile = pcFolder + 'PC_' + str(i) + '_' + InitCandidates + '.pickle'
                    if computedTag[i] == 0:
                        if check_file_exists(queryPCFile):
                            with open(queryPCFile, 'rb') as fval:
                                fdata = pickle.load(fval)
                                pcSet = fdata['pcSet']
                                queryID = fdata['queryID']
                                del fdata
                        else:
                            queryID, pcSet = pcObj.compute_pCohesion(i, candMethod=InitCandidates, queryPCFile=queryPCFile)
                        if queryID != i:
                            sys.exit("!!!!!!!!!!!!!!!!!!!! The computation/save is wrong for {} (real: {}) !!!!!!!!!!!!!!!!!!!!".format(i, queryID))
                        pcSets[queryID] = pcSet
                        if len(pcSet) > maxPCSize: maxPCSize = len(pcSet)
                        if len(pcSet) < minPCSize: minPCSize = len(pcSet)
                        meanPCSize += len(pcSet)
                meanPCSize = meanPCSize/len(verDeg)
                dataInfo = {'pcSets': pcSets, 'maxPCSize': maxPCSize, 'minPCSize': minPCSize, 'meanPCSize': meanPCSize}; savePickle(dataSet = dataInfo, filename = PCFile)
            except:
                tryResult = False
        if not parallelTag or not tryResult:
            pcSets, maxPCSize, minPCSize, meanPCSize = pcObj.pCohesionAll(PCFile, InitCandidates, verOriTags = verOriTags)
        del pcObj
    else:
        with open(PCFile, 'rb') as fval:
            pdata = pickle.load(fval)
            pcSets = pdata["pcSets"]
            minPCSize = pdata["minPCSize"]
            maxPCSize = pdata["maxPCSize"]
            meanPCSize = pdata["meanPCSize"]
            del fval

    #################### p-cohesion computation finished ####################
    print('==================== The p({:.2f})-Cohesion({}) size is as follows ==================== \n max = {:d} \n min = {:d} \n meanSize = {:.2f}'.format(p, pcChoose, maxPCSize, minPCSize, meanPCSize))
    return pcSets
