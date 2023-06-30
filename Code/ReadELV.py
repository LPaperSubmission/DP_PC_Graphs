from ELV import *
########################################################################################################################
def _ELV_(verDeg, neiSet, verOriTags = None, neiHop = 2, openCheck = False, parallelTag = False, ResultsFolder = None):
    print("========================================================================================================================")
    #################### Start to generate p-cohesion for all node ####################
    ELVFile = ResultsFolder + 'ELV_Nei_' + str(neiHop) + '-hop.pickle'
    ELVFolder = ResultsFolder + 'ELVs/'; mymkdir(ELVFolder)
    if not check_file_exists(ELVFile, rmTag=False):
        ########################################  generate p-cohesion for each node ########################################
        elvObj = ELV(verDeg, neiSet, neiHop, openCheck)
        tryResult = True
        if parallelTag:
            try:
                ELVsets = [[] for i in range(len(verDeg))]
                maxELVSize = 0
                minELVSize = len(verDeg)
                meanELVSize = 0.0
                results = []
                pool = Pool(processes=cpu_count())
                for queryID in range(len(verDeg)):
                    queryELVFile = ELVFolder + 'ELV_' + str(queryID) + '.pickle';
                    if not check_file_exists(queryELVFile):
                        results.append(pool.apply_async(elvObj.localExpand, (),{"_queryId": queryID, 'curHopNum': neiHop, 'filename': queryELVFile}))
                pool.close()
                pool.join()
                computedTag = [0] * len(verDeg)
                for res in results:
                    queryID, elvSet = res.get()
                    ELVsets[int(queryID)] = elvSet; computedTag[int(queryID)] = 1
                    if len(ELVsets) > maxELVSize: maxELVSize = len(elvSet)
                    if len(elvSet) < minELVSize: minELVSize = len(elvSet)
                    meanELVSize += len(elvSet)
                for i in range(len(verDeg)):
                    queryELVFile = ELVFolder + 'ELV_' + str(i) + '.pickle'
                    if computedTag[i] == 0:
                        if check_file_exists(queryELVFile):
                            with open(queryELVFile, 'rb') as fval:
                                fdata = pickle.load(fval)
                                elvSet = fdata['elvSet']
                                queryID = fdata['queryID']
                                del fdata
                        else:
                            queryID, elvSet = elvObj.localExpand(_queryId = i, curHopNum=neiHop, filename=queryELVFile)
                        if queryID != i:
                            sys.exit("!!!!!!!!!!!!!!!!!!!! ELV:The computation/save is wrong for {} (real: {}) !!!!!!!!!!!!!!!!!!!!".format(i, queryID))
                        ELVsets[int(queryID)] = elvSet
                        if len(elvSet) > maxELVSize: maxELVSize = len(elvSet)
                        if len(elvSet) < minELVSize: minELVSize = len(elvSet)
                        meanELVSize += len(elvSet)
                meanELVSize = meanELVSize/len(verDeg)
                dataInfo = {'ELVsets': ELVsets, 'maxELVSize': maxELVSize, 'minELVSize': minELVSize, 'meanELVSize': meanELVSize}; savePickle(dataSet = dataInfo, filename = ELVFile)
            except:
                tryResult = False
        if not parallelTag or not tryResult:
            ELVsets, maxELVSize, minELVSize, meanELVSize = elvObj.ELVAll(ELVFile, verOriTags = verOriTags)
        del elvObj
    else:
        with open(ELVFile, 'rb') as fval:
            pdata = pickle.load(fval)
            ELVsets = pdata["ELVsets"]
            maxELVSize = pdata["maxELVSize"]
            minELVSize = pdata["minELVSize"]
            meanELVSize = pdata["meanELVSize"]
            del fval
    #################### ELV computation finished ####################
    print('==================== The ELV size is as follows ==================== \n max = {:d} \n min = {:d} \n meanSize = {:.2f}'.format(maxELVSize, minELVSize, meanELVSize))
    return ELVsets
