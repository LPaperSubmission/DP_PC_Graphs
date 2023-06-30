import multiprocessing
from queue import PriorityQueue
from publicFus import *
########################################################################################################################
#################### For p-cohesion computation ####################
class pCohesionOld(object):
    def __init__(self, verDeg, neiSet, p, openCheck):
        self.verDeg = verDeg
        self.neiSet = neiSet
        self.p = p
        self.verPDeg = verDeg
        self.openCheck = openCheck
    ########## sorted neighbours based on candidate degree ##########
    def compCtri(self, vid, candDeg, candTag):
        if candTag[vid] > 0: return 0
        contriCount = 0
        for nid in self.neiSet[vid]:
            if candTag[nid] > 0 and candDeg[nid] < self.verPDeg[nid]:
                contriCount += 1
        return contriCount

    def compCommonNei(self, queryID, vid):
        return len([val for val in self.neiSet[queryID] if val in self.neiSet[vid]])

    def sortedNeiBasedContriDeg(self, vid, candDeg, candTag):
        nSetLocal = self.neiSet[vid]
        nContriCount = [0] * len(nSetLocal)
        cDeg = [candDeg[nid] for nid in nSetLocal]
        Idx = 0
        for nid in nSetLocal:
            nContriCount[Idx] = self.compCtri(nid, candDeg, candTag)
            Idx += 1
        # lzipped = zip(nContriCount, nSetLocal)

        lzipped = zip(cDeg, nSetLocal)
        res = zip(*sorted(lzipped, reverse= True)) # contribution: 9, 8, ...
        _, nSet_ = [list(x) for x in res]
        return nSet_

    def sortNeiBasedCommonNeiCandDegOriDeg(self, queryID, vid, candDeg):
        cnDeg = [0] * len(self.neiSet[vid])
        pDeg = [0] * len(self.neiSet[vid])
        cNeiSets = [0] * len(self.neiSet[vid])
        for i in range(len(self.neiSet[vid])):
            nid = self.neiSet[vid][i]
            cNeiSets[i] = len([val for val in self.neiSet[queryID] if val in self.neiSet[nid]]) # has more common neighbours with query vertex
            cnDeg[i] = candDeg[nid]
            pDeg[i] = self.verPDeg[nid]
        lzipped = zip(cNeiSets, cnDeg, pDeg, self.neiSet[vid])
        res = zip(*sorted(lzipped, reverse=True))  # commonneighbours: 9, 8, 7, candidate degree: 9, 8, 7, ...; original degree: 9, 8, 7, ....; ID: 9, 8, 7, ...
        _, _, _, nSet_ = [list(x) for x in res]
        return list(nSet_)

    def sortNeiBasedCandDegOriDeg(self, vid, candDeg):
        cnDeg = [0]*len(self.neiSet[vid])
        pDeg = [0] * len(self.neiSet[vid])
        for i in range(len(self.neiSet[vid])):
            nid = self.neiSet[vid][i]
            cnDeg[i] = candDeg[nid]
            pDeg[i] = self.verPDeg[nid]
        lzipped = zip(cnDeg, pDeg, self.neiSet[vid])
        res = zip(*sorted(lzipped, reverse = True)) # candidate degree: 9, 8, 7, ...; original degree: 9, 8, 7, ....; ID: 9, 8, 7, ...
        _, _, nSet_ = [list(x) for x in res]
        return list(nSet_)

    ########## update merit and penalty ##########
    def updateMeritPenalty(self, queryID, vid, candTag, candDeg):
        ########## When added vid's neighbours into candidate set, vid needs less number of neighbours to meet the p threshold
        ########## Merit computation ##########
        contriCounts = self.compCtri(vid, candDeg, candTag)
        commonNeiWithQuery = self.compCommonNei(queryID, vid)
        if commonNeiWithQuery == 0: commonNeiWithQuery = 1
        meritVid = contriCounts
        ########## Penalty computation ##########
        nSet = self.sortNeiBasedCommonNeiCandDegOriDeg(queryID, vid, candDeg)
        penaltyCount = 0; idx = 0; threshold = self.verPDeg[vid] - candDeg[vid]
        for nid in nSet:
            if idx < threshold:
                if candTag[nid] < 1:
                    idx += 1
                    penaltyCount += candDeg[nid]/len(self.neiSet[nid])
            else:
                break
        penaltyVid = penaltyCount  # update penalty
        if penaltyVid < 0: penaltyVid = 0
        ########## return result ##########
        return meritVid, penaltyVid

    ########## update candidate set ##########
    def updateCandSet(self, queryID, vid, candSet, candDeg, candTag, Merit, Penalty):
        # when add vid to candidate set, sets need to be updated
        candSet.append(vid)
        candTag[vid] = 1
        # nSetLocal = self.sortedNeiBasedContriDeg(vid, candDeg, candTag)
        nSetLocal = self.neiSet[vid]
        for nid in nSetLocal:
            candDeg[nid] += 1
            nDeg = self.verPDeg[nid] - candDeg[nid] # the number of neighbours need to be added to the candidate set
            if nDeg < 0 and candTag[nid] > 0 and candTag[nid] != 2:
                nDeg = 0
                candTag[nid] = 2
            #################### merit and penalty ####################
            Merit[nid], Penalty[nid] = self.updateMeritPenalty(queryID, nid, candTag, candDeg)

        Merit[vid], Penalty[vid] = self.updateMeritPenalty(queryID, vid, candTag, candDeg)
        # #candDeg[vid] / math.pow(len(self.neiSet[vid]), 2) * meritCount
        # Penalty[vid] = 0 #(max(self.verPDeg[vid] - candDeg[vid], 0)) / len(self.neiSet[vid]) / penaltyCount

    def expandCheck(self, candDeg, candTag):
        for v in range(len(candDeg)):
            vd = 0
            for u in self.neiSet[v]:
                if candTag[u] > 0:
                    vd += 1
            if vd != candDeg[v]:
                sys.exit("!!!!!!!!!!!!!!!!!!!! The degree of node {:d} is wrong ({:d} != {:d}), please check !!!!!!!!!!!!!!!!!!!!".format(v, vd, candDeg[v]))

    ########## local expend ##########
    def localExpand(self, queryID, candSet, candDeg, candTag, merit, penalty):
        pq = PriorityQueue()

        for cid in candSet:
            pq.put((-1 * candDeg[cid], -1 * self.verPDeg[cid], -1 * cid)) # candidate degree: 9, 8, 7, ...; original degree: 9, 8, 7, ....; ID: 9, 8, 7, ...
        ## First, expand the vertex in candidate set with largest candidate degree
        while not pq.empty():
            vid = -1 * (list(pq.get())[-1])
            if candDeg[vid] >= self.verPDeg[vid]: continue

            npq = PriorityQueue()
            ########## generate candidate sets ##########
            for vnid in self.neiSet[vid]:
                npq.put((-1 * merit[vnid], penalty[vnid], (-1) * vnid)) # merit: 9, 8, 7, ...; penalty: 1, 2, 3, ....; ID: 9, 8, 7, ...
            ## the neighbours with the largest merit and less penalty will be added to the candidate set
            ########## start to expand ##########
            while not npq.empty():
                nid = -1 * (list(npq.get())[-1])

                if candTag[nid] < 1: # the neighbour is not in the candidate set
                    deg = 0 # used to check if the candidate degree is correct
                    # pq.put((-1*merit[nid], penalty[nid], -1 * nid)) # merit the most and penalty the less, nid the large
                    pq.put((-1 * candDeg[nid], -1*self.verPDeg[nid], -1 * nid))  # merit the most and penalty the less, nid the large
                    candSet.append(nid) # add the 'nid' to the candidate set
                    candTag[nid] = 1 # notice the candTag

                    for nnid in self.neiSet[nid]: # update the new added vertex's neighbours
                        candDeg[nnid] += 1
                        if candTag[nnid] > 0: # the neighbour already in the candidate set
                            deg += 1
                            if candDeg[nnid] == self.verPDeg[nnid]: # helps it's neighbour meet the p threshold
                                for nnnid in self.neiSet[nnid]: # update the neighbours merit and penalty
                                    merit[nnnid] -= candDeg[nnnid]
                                    if merit[nnnid] <= 0:
                                        merit[nnnid] = candDeg[nnnid]
                        merit[nnid], penalty[nnid] = self.updateMeritPenalty(queryID, nnid, candTag, candDeg)
                    if candDeg[nid] != deg:
                        print("Wrong add id degree, please check, should be {:d} (but {:d})".format(deg, candDeg[nid]))
                if candDeg[vid] == self.verPDeg[vid]:
                    candTag[vid] = 1
                    break
            if self.openCheck:
                self.expandCheck(candDeg, candTag)
        return candSet, candTag, candDeg

    ########## check minimal ##########
    def minimalCheckSub(self, vid, candTag, candDeg):
        nodeNum = len(self.neiSet)
        stopTag = False; recoverTag = False
        stopID = 0; stopaid = 0
        if self.openCheck:
            bacDeg = [deg for deg in candDeg]
        ########## can delete vertex ##########
        canDelSet = []
        canDelTag = [0] * nodeNum # 0: not delete; 1:? 2: to be deleted; 3: deleted
        ########## the first is the input vertex ID ##########
        canDelSet.append(vid)
        canDelTag[vid] = 2
        ########## trying to remove the input vertex and others ##########
        for delID in canDelSet:
            candTag[delID] = 0 # removed from candidates
            canDelTag[delID] = 3 # removed
            for nid in self.neiSet[delID]: # check its neighbours
                candDeg[nid] -= 1
                if candTag[nid] > 0: # is already in the candidate set
                    if candDeg[nid] < self.verPDeg[nid]:
                        if candTag[nid] == 1 and canDelTag[nid] == 0:
                            canDelSet.append(nid)
                            canDelTag[nid] = 2
                        elif candTag[nid] == 3:
                            stopTag = True
                            stopID = nid
                            stopaid = delID
                            break
            if stopTag:
                recoverTag = True
                break
        if recoverTag:
            for cdID in canDelSet:
                candTag[cdID] = 1
                if canDelTag[cdID] != 3:
                    break
                for ncdID in self.neiSet[cdID]:
                    candDeg[ncdID] += 1
                    if ncdID == stopID and cdID == stopaid:
                        break
            candTag[vid] = 3
            if candDeg[vid] == self.verPDeg[vid]:
                tmpSet = []
                tmpSet.append(vid)
                while len(tmpSet) > 0:
                    tid = tmpSet.pop(-1)
                    for ntid in self.neiSet[tid]:
                        if candTag[ntid] == 1:
                            candTag[ntid] = 3
                            if candDeg[ntid] == self.verPDeg[ntid]:
                                tmpSet.append(ntid)
            if self.openCheck:
                for cIdx in range(nodeNum):
                    if candDeg[cIdx] != bacDeg[cIdx]:
                        print("minimalCheckSub: wrong recover, please check!")
        return candTag, candDeg

    def minimalCheck(self, candSet, candTag, candDeg, merit, penalty, mSet = None):
        mpq = PriorityQueue()
        cTag = [x for x in candTag]
        for cid in candSet:
            mpq.put((merit[cid], -1* penalty[cid], candDeg[cid], cid))
            if mSet is not None and cid in mSet:
                cTag[cid] = 3
        while not mpq.empty():
            vid = list(mpq.get())[-1]

            if cTag[vid] == 1:
                candTag, candDeg = self.minimalCheckSub(vid, cTag, candDeg)
        return candTag, candDeg

    ########## generated p-cohesion for a queryID ##########
    def compute_pCohesion(self, queryID, candMethod = "All", queryPCFile = None):
        nodeNum = len(self.verDeg)
        #################### generate p-cohesion degrees for each vertex ####################
        self.verPDeg = [math.ceil(deg * self.p) for deg in self.verDeg]
        #################### Initialize the penalty and merit sets based on the density (\frac{|E|}{|V|}) ####################
        Penalty = [0.0] * nodeNum
        Merit = [0.0] * nodeNum
        #################### Generate candidates and candidate degrees ####################
        CandDeg = [0] * nodeNum
        CandTag = [0] * nodeNum # 0: not in the candidate set; 1: in the candidate set; 3: must include
        CandSet = []
        mustSet = [queryID]
        #################### add one node to the candidate set ####################
        self.updateCandSet(queryID, queryID, CandSet, CandDeg, CandTag, Merit, Penalty)
        #################### Start to add must include nodes ####################
        if candMethod in ['All', 'all', 'ALL']:
            mustSet.extend(self.neiSet[queryID])
            for nid in self.neiSet[queryID]:
                self.updateCandSet(queryID, nid, CandSet, CandDeg, CandTag, Merit, Penalty)

        CandSet, CandTag, CandDeg = self.localExpand(queryID, CandSet, CandDeg, CandTag, Merit, Penalty)
        CandTag, _ = self.minimalCheck(CandSet, CandTag, CandDeg, Merit, Penalty, mSet= mustSet)
        ########## p-Cohesion ##########
        pcSet = [tID for tID in range(len(CandTag)) if CandTag[tID] > 0]
        if queryPCFile is not None:
            dataInfo = {'pcSet': pcSet.copy(), 'queryID': queryID}; savePickle(dataSet=dataInfo, filename=queryPCFile)
        return queryID, pcSet.copy()

    ########## generate p-cohesion ##########
    def pCohesionAll(self, PCFile, InitCandidates, verOriTags = None):
        pcSets = [[] for i in range(len(self.verDeg))]
        maxPCSize = 0
        minPCSize = len(self.verDeg)
        meanPCSize = 0
        for queryID in range(len(self.verDeg)):
            now = datetime.datetime.now()
            if self.openCheck:
                if verOriTags is not None:
                    print("==================== The query ID is {:d} (OID: {:d}) ====================".format(queryID, verOriTags.index(queryID)))
                else:
                    print("==================== The query ID is {:d} ====================".format(queryID))
            _, pcSets[queryID] = self.compute_pCohesion(queryID, candMethod = InitCandidates)
            end = datetime.datetime.now()
            print('==================== {} used {:.2f}s to finish ({}) ===================='.format(queryID, end-now, len(pcSets[queryID])))
            if len(pcSets[queryID]) > maxPCSize: maxPCSize = len(pcSets[queryID])
            if len(pcSets[queryID]) < minPCSize: minPCSize = len(pcSets[queryID])
            meanPCSize += len(pcSets[queryID])

            if self.openCheck:
                if verOriTags is not None:
                    print('The size of p-cohesion for {:d} (OID: {:d}), with p = {:.2f} is {:d} ({})'.format(queryID, verOriTags.index(queryID), self.p, len(pcSets[queryID]), list([verOriTags.index(qid) for qid in pcSets[queryID]])))
                else:
                    print('The size of p-cohesion for {:d}, with p = {:.2f} is {:d}'.format(queryID, self.p, len(pcSets[queryID])))
        meanPCSize /= len(self.verDeg)
        ########## save p-cohesions for future use ##########
        dataInfo = {'pcSets': pcSets, 'maxPCSize': maxPCSize, 'minPCSize': minPCSize, 'meanPCSize': meanPCSize}; savePickle(dataSet = dataInfo, filename = PCFile)
        return pcSets, maxPCSize, minPCSize, meanPCSize