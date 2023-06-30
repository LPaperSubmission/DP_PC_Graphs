from publicFus import *
########################################################################################################################
######################################## For Extend Local View (ELV) computation ########################################
class ELV(object):
    def __init__(self, verDeg, neiSet, neiHop, openCheck):
        self.verDeg = verDeg
        self.neiSet = neiSet
        self.neiHop = neiHop
        self.openCheck = openCheck

    #################### Local expand ####################
    def localExpand(self, _queryId = 0, curHopNum = 2, filename = None):
        elvSet = []; elvSet.append(_queryId)
        while curHopNum >= 1:
            lenElvSet = len(elvSet)
            for vIdx in range(lenElvSet):
                vid = elvSet[vIdx]
                for nid in self.neiSet[vid]:
                    if nid not in elvSet:
                        elvSet.append(nid)
            curHopNum -= 1

        if self.openCheck:
            if self.neiHop == 2:
                localELVSet = []; localELVSet.append(_queryId)
                for nid in self.neiSet[_queryId]:
                    if nid not in localELVSet:
                        localELVSet.append(nid)
                lelv = len(localELVSet)
                for vIdx in range(lelv):
                    vid = localELVSet[vIdx]
                    for nid in self.neiSet[vid]:
                        if nid not in localELVSet:
                            localELVSet.append(nid)
                if len(elvSet) != len(localELVSet):
                    sys.exit("!!!!!!!!!!!!!!!!!!!! The extend is wrong ({} != {}) !!!!!!!!!!!!!!!!!!!!".format(len(elvSet), len(localELVSet)))
                else:
                    for vid in localELVSet:
                        if vid not in elvSet:
                            sys.exit("!!!!!!!!!!!!!!!!!!!! Node are not same ({} should in elvSet) !!!!!!!!!!!!!!!!!!!!".format(vid))

        if filename is not None:
            dataInfo = {'elvSet': elvSet.copy(), 'queryID': _queryId}; savePickle(dataSet=dataInfo, filename=filename)
        return _queryId, elvSet.copy()

    #################### generate p-cohesion ####################
    def ELVAll(self, ELVFile, verOriTags = None):
        ELVsets = [[] for i in range(len(self.verDeg))]
        maxELVSize = 0
        minELVSize = len(self.verDeg)
        meanELVSize = 0.0
        for queryID in range(len(self.verDeg)):
            if self.openCheck and verOriTags is not None: print("==================== The query ID is {:d} (OID: {:d}) ====================".format(queryID, verOriTags.index(queryID)))
            _, ELVsets[queryID] = self.localExpand(_queryId = queryID, curHopNum = self.neiHop)
            if len(ELVsets[queryID]) > maxELVSize: maxELVSize = len(ELVsets[queryID])
            if len(ELVsets[queryID]) < minELVSize: minELVSize = len(ELVsets[queryID])
            meanELVSize += len(ELVsets[queryID])
            if self.openCheck and verOriTags is not None: print('The size of ELV for {:d} (OID: {:d}) is {:d} ({})'.format(queryID, verOriTags.index(queryID), len(ELVsets[queryID]), list([verOriTags.index(qid) for qid in ELVsets[queryID]])))
        meanELVSize /= len(self.verDeg)
        ########## save ELV for future use ##########
        dataInfo = {'ELVsets': ELVsets, 'maxELVSize': maxELVSize, 'minELVSize': minELVSize, 'meanELVSize': meanELVSize}; savePickle(dataSet = dataInfo, filename = ELVFile)
        return ELVsets, maxELVSize, minELVSize, meanELVSize