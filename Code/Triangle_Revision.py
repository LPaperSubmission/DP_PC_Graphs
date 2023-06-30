from publicFus import *
########################################################################################################################
######################################## For triangle counting ########################################
class Triangle(object):
    def __init__(self, **kwargs):
        ########## Public data info ##########
        self.verDeg = getVals('verDeg', **kwargs)
        self.NodeNum = len(self.verDeg)
        self.neiSet = getVals('neiSet', **kwargs)
        self.ResultsFolder = getVals('ResultsFolder', **kwargs)
        ########## Parameters for edge-rldp ##########
        self.randomstate = getVals('randomState', **kwargs)
        ########## Phase 1 ##########
        self.hPrime = getVals('hPrime', **kwargs)
        self.varEpsilon1 = getVals('varEpsilon1', **kwargs)
        self.deltaVal1 = getVals('deltaVal1', **kwargs)
        self.alphaVal = getVals('alphaVal', **kwargs)
        self.h = 0
        ########## Phase 2 ##########
        self.varEpsilon2 = getVals('varEpsilon2', **kwargs)
        self.deltaVal2 = getVals('deltaVal2', **kwargs)
        self.spVal = getVals("spVal", **kwargs)
        self.r4 = getVals('r4', **kwargs)
        self.prob = getVals('prob', **kwargs)
        ########## Others ##########
        self.openCheck = getVals('openCheck', **kwargs)
        self.noiseCal = getVals('noiseCal', **kwargs)

        self.rt = getVals('rt', **kwargs)
        self.dpmechanism = getVals('dpmechanism', **kwargs)
        ########## triangle & wedge Info ##########
        self.triInfoFile = getVals('dataInfoFile', **kwargs)
        self.lambdaELVFile = getVals('lambdaELVFile', **kwargs)
        ########## Wedge ##########
        self.cSets = None # all nodes in the local view
        self.cSetsPC = None
        self.PCcMaxVid = None
        self.cpSetsBar = None
        self.cSetsELV = None
        self.ceSetsBar = None
        self.ELVcMaxVid = None
        ########## Triangle
        self.triSets = None
        self.triSetsPC = None
        self.triSetsELV = None
        ########## Differential privacy ##########
        self.upperBoundPC = 0.0
        self.upperBoundELV = 0.0
        self.h = 0
        ########## pcohesion & ELV Related ##########
        self.PCSetsTag = False
        self.ELVSetsTag = False
        self.parallelTag = False
        self.resultFolder = None
        self.InitCandidates = getVals('InitCandidates', **kwargs); #kwargs['InitCandidates']

    ######################################## count common neighbours ########################################
    #################### triangle counting ####################
    def toExpandCandGen(self, curSet = None, queryID = -1):
        toExpandAll = [0] * self.NodeNum
        toExpandPC = [0] * self.NodeNum
        toExpandELV = [0] * self.NodeNum
        ########## CHECK CORRECTNESS ##########
        if queryID != curSet[0]: sys.exit('!!!!!!!!!!!!!!!!!!!! Wrong curSet !!!!!!!!!!!!!!!!!!!!')
        ########## READ PCOHESION AND ELV ##########
        pcSetTmp, ELVSetTmp = [], []
        if self.PCSetsTag:
            _, pcSetTmp, _ = ReadPCELV(curSet[0], self.resultFolder, PCSetsTag=self.PCSetsTag,
                                       ReadTag='PC', InitCandidates = self.InitCandidates, ELVFolder=self.ResultsFolder)
        if self.ELVSetsTag:
            _, _, ELVSetTmp = ReadPCELV(curSet[0], self.resultFolder, ELVSetsTag=self.ELVSetsTag, ReadTag='ELV', ELVFolder=self.ResultsFolder)
        ########## CHECK NODES CAN BE ADDED TO CURRENT SET ##########
        for vid in curSet:
            for nid in self.neiSet[vid]:
                if nid in curSet: continue
                if self.PCSetsTag and nid in pcSetTmp:
                    toExpandPC[nid] += 1
                if self.ELVSetsTag and nid in ELVSetTmp:
                    toExpandELV[nid] += 1
                toExpandAll[nid] += 1

        toExpandSet = [vid for vid in range(self.NodeNum) if toExpandAll[vid] == len(curSet)]
        toExpandSetPC = [vid for vid in range(self.NodeNum) if toExpandPC[vid] == len(curSet)]
        toExpandSetELV = [vid for vid in range(self.NodeNum) if toExpandELV[vid] == len(curSet)]

        del pcSetTmp
        del ELVSetTmp

        return toExpandSet, toExpandSetPC, toExpandSetELV

    def expandTokClique(self, kSet, toExpandSet, tri = None, triPC = None, triELV = None, processTag = 'tri', queryID = -1):
        ########## check correctness ##########
        if kSet[0] != queryID: sys.exit('!!!!!!!!!!!!!!!!!!!! Wrong k Set !!!!!!!!!!!!!!!!!!!!')
        ########## backup ##########
        kSetTmp = kSet.copy()
        ########## to expand more
        for nid in toExpandSet:
            kSet.append(nid)
            if len(kSet) < 3:
                toExpandSetTmp, toExpandSetTmpPC, toExpandSetTmpELV = self.toExpandCandGen(curSet=kSet, queryID = queryID)
                if processTag == 'PC' and self.PCSetsTag:
                    self.expandTokClique(kSet, toExpandSetTmpPC, tri = tri, triPC = triPC, triELV = triELV, processTag = processTag, queryID = queryID)
                elif processTag == 'ELV' and self.ELVSetsTag:
                    self.expandTokClique(kSet, toExpandSetTmpELV, tri = tri, triPC = triPC, triELV = triELV, processTag = processTag, queryID = queryID)
                else:
                    self.expandTokClique(kSet, toExpandSetTmp, tri = tri, triPC = triPC, triELV = triELV, processTag = processTag, queryID = queryID)
            else:
                if processTag == 'PC' and self.PCSetsTag: triPC[kSet[0]] += 1
                elif processTag == 'ELV' and self.ELVSetsTag: triELV[kSet[0]] += 1
                else: tri[kSet[0]] += 1
            kSet = kSetTmp.copy()

    def kCliqueGen(self, queryNode, tri = None, triPC = None, triELV = None):
        toexpand, toexpandPC, toexpandELV = self.toExpandCandGen(curSet=[queryNode], queryID = queryNode)
        self.expandTokClique([queryNode], toexpand, tri = tri, triPC = triPC, triELV = triELV, processTag = 'Normal', queryID = queryNode)
        self.expandTokClique([queryNode], toexpandPC, tri = tri, triPC = triPC, triELV = triELV, processTag = 'PC', queryID = queryNode)
        self.expandTokClique([queryNode], toexpandELV, tri = tri, triPC = triPC, triELV = triELV, processTag = 'ELV', queryID = queryNode)
    #################### wedge counting ####################
    def countWedgeOutSideLV(self):
        # this function is suitable for pCohesion and ELV
        pcSetTmp, ELVSetTmp = [], []
        for vid in range(self.NodeNum):
            cCountsPC, cCountsELV = 0, 0
            ########## p Cohesion ##########
            nidPC = self.PCcMaxVid[vid]
            for kid in self.neiSet[nidPC]:
                _, pcSetTmp, _ = ReadPCELV(vid, self.resultFolder, PCSetsTag = self.PCSetsTag, ReadTag='PC',
                                           InitCandidates=self.InitCandidates, ELVFolder=self.ResultsFolder)
                if kid == vid or kid in pcSetTmp: continue
                if kid in self.neiSet[vid]: cCountsPC += 1
            ########## ELV ##########
            nidELV = self.ELVcMaxVid[vid]
            for kid in self.neiSet[nidELV]:
                _, _, ELVSetTmp = ReadPCELV(vid, self.resultFolder, PCSetsTag = self.PCSetsTag, ELVSetsTag = self.ELVSetsTag,
                                            ReadTag='ELV', InitCandidates = self.InitCandidates, ELVFolder=self.ResultsFolder)
                if kid == vid or kid in ELVSetTmp: continue
                if kid in self.neiSet[vid]: cCountsELV += 1

            self.cpSetsBar[vid] = cCountsPC
            self.ceSetsBar[vid] = cCountsELV

            if cCountsELV > 0:  sys.exit("!!!!!!!!!!!!!!!!!!!! The cCountsELV  should be 0 !!!!!!!!!!!!!!!!!!!!")
            if self.cpSetsBar[vid] + self.cSetsPC[vid] > self.verDeg[vid]:
                sys.exit("!!!!!!!!!!!!!!!!!!!! Wrong pc-wedge computation !!!!!!!!!!!!!!!!!!!!")
            if self.ceSetsBar[vid] + self.cSetsELV[vid] > self.verDeg[vid]:
                sys.exit("!!!!!!!!!!!!!!!!!!!! Wrong pc-wedge computation !!!!!!!!!!!!!!!!!!!!")

            del pcSetTmp
            del ELVSetTmp

    def countWedgeParSub(self, vid):
        cCounts, cCountsPC, cCountsELV = 0, 0, 0
        pcNidMax, elvNidMax = -1, -1
        ELVSetTmp, pcSetTmp = [], []
        ########## Normal wedge Number ##########
        for nid in range(self.NodeNum):
            if nid == vid: continue
            csn = len([val for val in self.neiSet[vid] if val in self.neiSet[nid]])
            if csn > cCounts: cCounts = csn
        ########## Wedges in ELV ##########
        if self.ELVSetsTag:
            _, _, ELVSetTmp = ReadPCELV(vid, self.resultFolder, ELVSetsTag = self.ELVSetsTag, ReadTag='ELV', ELVFolder=self.ResultsFolder)
            for nid in ELVSetTmp:
                if nid == vid: continue
                csnELV = len([val for val in self.neiSet[vid] if val in self.neiSet[nid] if val in ELVSetTmp])
                if csnELV > cCountsELV:
                    cCountsELV = csnELV
                    elvNidMax = nid

        ########## Wedges in pCohesion ##########
        if self.PCSetsTag:
            _, pcSetTmp, _ = ReadPCELV(vid, self.resultFolder, PCSetsTag = self.PCSetsTag, ReadTag='PC',
                                       InitCandidates = self.InitCandidates, ELVFolder=self.ResultsFolder)
            for nid in pcSetTmp:
                if nid == vid: continue
                csnPC = len([val for val in self.neiSet[vid] if val in self.neiSet[nid] if val in pcSetTmp])
                if csnPC > cCountsPC:
                    cCountsPC = csnPC
                    pcNidMax = nid

        ########## CORRECTNESS CHECK ##########
        if cCounts > self.verDeg[vid]:
            sys.exit(
                "!!!!!!!!!!!!!!!!!!!! {} wedge number is larger than its degree: wedge = {}, deg = {}".format(
                    vid, cCounts, self.verDeg[vid]))
        if self.PCSetsTag:
            pcDeg = len([nid for nid in self.neiSet[vid] if nid in pcSetTmp])
            if cCountsPC > pcDeg:
                sys.exit(
                    "!!!!!!!!!!!!!!!!!!!! PC: {} wedge number is larger than its degree: wedge = {}, deg = {}".format(
                        vid, cCountsPC, pcDeg))
        if self.ELVSetsTag:
            ELVDeg = len([nid for nid in self.neiSet[vid] if nid in ELVSetTmp])
            if cCountsELV > ELVDeg:
                sys.exit(
                "!!!!!!!!!!!!!!!!!!!! ELV: {} wedge number is larger than its degree: wedge = {}, deg = {}".format(
                    vid, cCountsELV, ELVDeg))

        del pcSetTmp
        del ELVSetTmp

        return vid, cCounts, cCountsPC, cCountsELV, pcNidMax, elvNidMax

    def countWedgePar(self):
        if self.parallelTag:
            results = []
            try:
                pool = Pool(processes=cpu_count())
                for vid in range(self.NodeNum):
                    results.append(pool.apply_async(self.countWedgeParSub, (vid,), {}))
                pool.close(); pool.join()
                for res in results:
                    vid, cCounts, cCountsPC, cCountsELV, pcNidMax, elvNidMax = res.get()
                    self.cSets[vid] = cCounts
                    self.cSetsPC[vid] = cCountsPC
                    self.cSetsELV[vid] = cCountsELV
                    self.PCcMaxVid[vid] = pcNidMax
                    self.ELVcMaxVid[vid] = elvNidMax
            except:
                for vid in range(self.NodeNum):
                    vid, cCounts, cCountsPC, cCountsELV, pcNidMax, elvNidMax = self.countWedgeParSub(vid)
                    self.cSets[vid] = cCounts
                    self.cSetsPC[vid] = cCountsPC
                    self.cSetsELV[vid] = cCountsELV
                    self.PCcMaxVid[vid] = pcNidMax
                    self.ELVcMaxVid[vid] = elvNidMax
        else:
            for vid in range(self.NodeNum):
                vid, cCounts, cCountsPC, cCountsELV, pcNidMax, elvNidMax = self.countWedgeParSub(vid)
                self.cSets[vid] = cCounts
                self.cSetsPC[vid] = cCountsPC
                self.cSetsELV[vid] = cCountsELV
                self.PCcMaxVid[vid] = pcNidMax
                self.ELVcMaxVid[vid] = elvNidMax

    def CountWeges(self):
        self.cSets = [0] * self.NodeNum
        self.cSetsPC = [0] * self.NodeNum
        self.cSetsELV = [0] * self.NodeNum
        self.PCcMaxVid = [-1] * self.NodeNum
        self.ELVcMaxVid = [-1] * self.NodeNum
        #################### Call functions ====================
        self.countWedgePar()
        #################### Count wedges outsite the lv ####################
        self.cpSetsBar = [0] * self.NodeNum
        self.ceSetsBar = [0] * self.NodeNum
        #################### Call Functions ====================
        if self.PCSetsTag and self.ELVSetsTag: self.countWedgeOutSideLV()
    #################### Triangle computation ####################
    def countTrianglesOutSideLV(self):
        triSetsPCBar = [0.0] * self.NodeNum
        triSetsELVBar = [0.0] * self.NodeNum
        # only for p-cohesion and ELV
        for vid in range(self.NodeNum):
            cCountsPC, cCountsELV = 0, 0
            for nid in self.neiSet[vid]:
                for kid in self.neiSet[nid]:
                    if kid in self.neiSet[vid]:
                        ########## read p-cohesion and ELV ##########
                        rID, pcSetTmp, elvSetTmp = ReadPCELV(vid, self.resultFolder, PCSetsTag = self.PCSetsTag, ELVSetsTag = self.ELVSetsTag,
                                                             ReadTag='All', InitCandidates = self.InitCandidates, ELVFolder=self.ResultsFolder)
                        if rID != vid: sys.exit('!!!!!!!!!!!!!!!!!!!! Wrong pcohesion elv reading !!!!!!!!!!!!!!!!!!!!')
                        ########## pcohesion ##########
                        if nid in pcSetTmp and kid in pcSetTmp: pass
                        else: cCountsPC += 1
                        ########## ELV ##########
                        if nid in elvSetTmp and kid in elvSetTmp: pass
                        else: cCountsELV += 1
                        del elvSetTmp
                        del pcSetTmp
            triSetsPCBar[vid] = cCountsPC/2
            triSetsELVBar[vid] = cCountsELV/2
            if cCountsELV > 0: sys.exit("!!!!!!!!!!!!!!!!!!!! The cCountsELV  should be 0 !!!!!!!!!!!!!!!!!!!!")
            if self.triSetsPC[vid] + triSetsPCBar[vid] != self.triSets[vid]:
                sys.exit("!!!!!!!!!!!!!!!!!!!! Triangle Bar computation wrong {} ({} != {}) !!!!!!!!!!!!!!!!!!!!".format(
                    vid, self.triSetsPC[vid] + triSetsPCBar[vid], self.triSets[vid]))
            if self.triSetsELV[vid] + triSetsELVBar[vid] != self.triSets[vid]:
                sys.exit("!!!!!!!!!!!!!!!!!!!! Triangle Bar computation wrong {} ({} != {}) !!!!!!!!!!!!!!!!!!!!".format(
                    vid, self.triSetsELV[vid] + triSetsELVBar[vid], self.triSets[vid]))

    def CountTriangles(self):
        self.triSets = [0] * self.NodeNum
        self.triSetsPC = [0] * self.NodeNum
        self.triSetsELV = [0] * self.NodeNum
        for vid in range(self.NodeNum):
            pcSetTmp, ELVSetTmp = [], []
            if self.PCSetsTag:
                _, pcSetTmp, _ = ReadPCELV(vid, self.resultFolder, PCSetsTag=self.PCSetsTag,
                                           ELVSetsTag=self.ELVSetsTag, ReadTag='PC', InitCandidates = self.InitCandidates,
                                           ELVFolder=self.ResultsFolder)
            if self.ELVSetsTag:
                _, _, ELVSetTmp = ReadPCELV(vid, self.resultFolder, PCSetsTag=self.PCSetsTag,
                                            ELVSetsTag=self.ELVSetsTag, ReadTag='ELV', InitCandidates = self.InitCandidates,
                                            ELVFolder=self.ResultsFolder)
            for nid in self.neiSet[vid]:
                self.triSets[vid] += len([val for val in self.neiSet[vid] if val in self.neiSet[nid]])
            self.triSets[vid] /= 2

            if self.PCSetsTag:
                for nid in pcSetTmp:
                    if nid == vid: continue
                    self.triSetsPC[vid] += len([val for val in self.neiSet[vid] if nid in self.neiSet[vid] if val in self.neiSet[nid] if val in pcSetTmp]) # request common neighbours in the PCohesion
                self.triSetsPC[vid] /= 2

            if self.ELVSetsTag:
                for nid in ELVSetTmp:
                    if nid == vid: continue
                    self.triSetsELV[vid] += len([val for val in self.neiSet[vid] if nid in self.neiSet[vid] if val in self.neiSet[nid] if val in ELVSetTmp]) # request common neighbours in the PCohesion
                self.triSetsELV[vid] /= 2

            del pcSetTmp
            del ELVSetTmp

        if self.openCheck:
            triSetTmp = [0] * self.NodeNum
            triSetPCTmp = [0] * self.NodeNum
            triSetELVTmp = [0] * self.NodeNum
            for vid in range(self.NodeNum):
                self.kCliqueGen(vid, tri=triSetTmp, triPC=triSetPCTmp, triELV=triSetELVTmp)
                triSetTmp[vid] /= 2
                triSetPCTmp[vid] /= 2
                triSetELVTmp[vid] /= 2

                if triSetTmp[vid] != self.triSets[vid]:
                    sys.exit('!!!!!!!!!!!!!!!!!!!! Triangles computation wrong {} ({} != {})!!!!!!!!!!!!!!!!!!!!'.format(vid, triSetTmp[vid], self.triSets[vid]))

        if self.openCheck:
            if self.PCSetsTag and self.ELVSetsTag:
                self.countTrianglesOutSideLV()

    def readPCELV(self, vid, ReadTag = 'PC'):
        lcSet = None
        if ReadTag == 'PC':
                _, lcSet, _ = ReadPCELV(vid, self.resultFolder, PCSetsTag=self.PCSetsTag,
                                           ELVSetsTag=self.ELVSetsTag, ReadTag= ReadTag, InitCandidates = self.InitCandidates,
                                           ELVFolder=self.ResultsFolder)
        elif ReadTag == 'ELV':
            _, _, lcSet = ReadPCELV(vid, self.resultFolder, PCSetsTag=self.PCSetsTag,
                                        ELVSetsTag=self.ELVSetsTag, ReadTag=ReadTag, InitCandidates = self.InitCandidates,
                                        ELVFolder=self.ResultsFolder)
        else:
            sys.exit('!!!!!!!!!!!!!!!!!!!! Please give the correct read tag (PC, ELV) !!!!!!!!!!!!!!!!!!!!')
        deg = len([nid if nid in lcSet for nid in self.neiSet[vid]])
        return deg

    ######################################## two-phase Approach ########################################

    def readTriangleInfo(self):
        with open(self.triInfoFile, 'rb') as fval:
            fdata = pickle.load(fval)
            ########## Wedges ##########
            self.cSets = fdata['cSets']
            self.cSetsPC = fdata['cSetsPC']
            self.cpSetsBar = fdata['cpSetsBar']
            self.cSetsELV = fdata['cSetsELV']
            self.ceSetsBar = fdata['ceSetsBar']
            ########## triangles ##########
            self.triSets = fdata['triSets']
            self.triSetsPC = fdata['triSetsPC']
            self.triSetsELV = fdata['triSetsELV']
    
    def GlobalDataCorrelationMeasuring(self, hInput = None):
        ########## PC ELV ##########
        pcDegSet, elvDegSet = None, None
        if self.PCSetsTag:
            pcDegSet = [self.readPCELV(i, ReadTag = 'PC') for i in range(self.NodeNum)]
            elvDegSet = [self.readPCELV(i, ReadTag = 'ELV') for i in range(self.NodeNum)]
            assert np.array(elvDegSet).all() == np.array(self.verDeg).all()
        ########## Public ##########
        noiseScale = 2 / self.alphaVal / self.varEpsilon1
        logDelta = np.log((self.hprime + 1) / self.deltaVal1)
        LapLambdaD = np.round(self.randomstate.laplace(0, noiseScale, self.NodeNum), 3)
        ########## Upper bound ##########
        verDegPCUpper, verDegELVUpper = None, None
        ########## normal
        verDegPrime = [self.verDeg[i] + LapLambdaD[i] for i in range(self.NodeNum)]
        verDegUpper = [verDegPrime[i] + noiseScale * logDelta for i in range(self.NodeNum)]
        ########## pc 
        verPCDegPrime = [pcDegSet[i] + LapLambdaD[i] for i in range(self.NodeNum)]
        verDegPCUpper = [verPCDegPrime[i] + noiseScale * logDelta for i in range(self.NodeNum)]
        ########## elv
        verELVDegPrime = [elvDegSet[i] + LapLambdaD[i] for i in range(self.NodeNum)]
        verDegELVUpper = [verELVDegPrime[i] + noiseScale * logDelta for i in range(self.NodeNum)]
        #==================== call function ====================
        if self.noiseCal:
            verDegUpper = noiseCalibrate(verDegUpper.copy(), self.verDeg.copy())
            if self.PCSetsTag:
                verDegPCUpper = noiseCalibrate(verDegPCUpper.copy(), pcDegSet.copy())
                verDegELVUpper = noiseCalibrate(verDegELVUpper.copy(), elvDegSet.copy())
        ########## sort S{v_i} ##########
        ########## normal
        verDegUpper = np.array(verDegUpper)
        Ind = verDegUpper.argsort()[::-1] # descending order
        sortedDegIDSet = np.array([[verDegUpper[idx], idx] for idx in Ind])
        sortedDegIDSetPC, sortedDegIDSetELV = None, None
        if self.PCSetsTag:
            ########## PC
            verDegPCUpper = np.array(verDegPCUpper)
            Ind = verDegPCUpper.argsort()[::-1] # descending order
            sortedDegIDSetPC = np.array([[verDegPCUpper[idx], idx] for idx in Ind])
            ########## ELV
            verDegELVUpper = np.array(verDegELVUpper)
            Ind = verDegELVUpper.argsort()[::-1] # descending order
            sortedDegIDSetELV = np.array([[verDegELVUpper[idx], idx] for idx in Ind])
        ########## compute h ##########
        self.h = hInput
        if self.openCheck: print('==================== h = {:d} ===================='.format(self.h))
        ########## get S ##########
        S, SPC, SELV = [], [], []
        for i in range(self.h+1):
            S.append(int(sortedDegIDSet[i, 1]))
            if self.PCSetsTag:
                SPC.append(int(sortedDegIDSetPC[i, 1]))
                SELV.append(int(sortedDegIDSetELV[i, 1]))

        ##########
        lambdaC = self.h / ((1-self.alphaVal) * self.varEpsilon1)
        cUpperSet = [0] *self.NodeNum
        cUpperSetPC = [0] *self.NodeNum
        cUpperSetELV = [0] * self.NodeNum
        kStarSetPC = [0] * self.NodeNum
        kStarSetELV = [0] * self.NodeNum
        ########## laplace noise generation ##########
        LapLambdaC = np.round(self.randomstate.laplace(0, lambdaC, len(S)), 3)
        ########## 
        if not self.PCSetsTag:
            cUpperSetPCPart = [self.cSets[S[sIdx]] + LapLambdaC[sIdx] + lambdaC * logDelta for sIdx in range(len(S))]
            truth = [self.cSets[vid] for vid in S]
            # ==================== call function ====================
            if self.noiseCal: cUpperSetPCPart = noiseCalibrate(cUpperSetPCPart.copy(), truth.copy())
            kStarSetPC = [min(verDegUpper[i], cUpperSetPC[i]) + 2 for i in rang(self.NodeNum)]
            kStarSetELV = kStarSetPC.copy()
        else:
            ########## p cohesion ##########
            cUpperSetPCPart = [self.cSetsPC[SPC[sIdx]] + LapLambdaC[sIdx] + lambdaC * logDelta for sIdx in range(len(SPC))]
            truth = [self.cSetsPC[vid] for vid in SPC]
            # ==================== call function ====================
            if self.noiseCal: cUpperSetPCPart = noiseCalibrate(cUpperSetPCPart.copy(), truth.copy())
            ########## ELV ##########
            cUpperSetELVPart = [self.cSetsELV[SELV[sIdx]] + LapLambdaC[sIdx] + lambdaC * logDelta for sIdx in range(len(SELV))]
            truth = [self.cSetsELV[vid] for vid in SELV]
            # ==================== call function ====================
            if self.noiseCal: cUpperSetELVPart = noiseCalibrate(cUpperSetELVPart.copy(), truth.copy())
            ########## k* ##########
            kStarSetPC = [cUpperSetPC[i] + 2 if i not in SPC else min(verDegPCUpper[i], cUpperSetPC[i]) + 2 for i in rang(self.NodeNum)]
            kStarSetELV = [cUpperSetELV[i] + 2 if i not in SELV else min(verDegELVUpper[i], cUpperSetELV[i]) + 2 for i in rang(self.NodeNum)]
        ####################
        self.lambdaELVFile = self.lambdaELVFile + '_h' + str(self.h) + '_var1' + str(round(self.varEpsilon1, 2)) + '_var2' + str(round(self.varEpsilon2,2)) + '.pickle'

        if check_file_exists(self.lambdaELVFile):
            with open(self.lambdaELVFile, 'rb') as fval:
                fdata = pickle.load(fval)
                kStarELV = fdata['kStarELV']
        else:
            dataInfo = {'kStarELV': np.max(kStarSetELV)}
            savePickle(dataInfo, self.lambdaELVFile)
            
        return np.max(kStarSetPC), np.max(kStarSetELV), verDegPrime.copy(), verPCDegPrime.copy(), verELVDegPrime.copy()


    def sampleTriNo(self, vid, ci):
        ciSamp = 0
        for cidx in range(ci):
            if np.random.rand() > self.prob:
                ciSamp += 1
        return ciSamp


    def GlobalTriangleCountingCollection(self, hInput = None, kStarPC = 0, kStarELV = 0, verDegPrime = None, verPCDegPrime = None, verELVDegPrime = None):
        ######################################## Phase 1 ########################################
        kStarPC, kStarELV, verDegPrime, verPCDegPrime, verELVDegPrime= self.GlobalDataCorrelationMeasuring(hInput = hInput)

        ######################################## Phase 2 ########################################
        #################### compute e_2' ####################
        deltaVal3 = self.deltaVal1

        km2PC = kStarPC - 2
        km2ELV = kStarELV -2

        sqrtPC = np.sqrt(km2PC * np.log(1/deltaVal3))
        sqrtELV = np.sqrt(2 * km2ELV * np.log(1/deltaVal3))

        e2pPC = sy.symbols('e2pPC', real=True)
        e2pELV = sy.symbols('e2pELV', real=True)
        exprPC = 2 * e2pPC + sy.log(1 + self.prob * (sy.exp(e2pPC) - 1)) * (sqrtPC + km2PC * self.prob * (sy.exp(e2pPC - 1)) / (self.prob * (sy.exp(e2pPC + 1)) + 2)) - varEpsilon2
        exprELV = 2 * e2pELV + sy.log(1 + self.prob * (sy.exp(e2pELV) - 1)) * (sqrtPC + km2PC * self.prob * (sy.exp(e2pELV - 1)) / (self.prob * (sy.exp(e2pELV + 1)) + 2)) - varEpsilon2

        varEpsilon2PrimePC = sy.solve(exprPC, e2pPC)
        varEpsilon2PrimeELV = sy.solve(exprELV, e2pELV)
        #################### upper bound of local sensitivity ####################
        LSUppers = [0] * self.NodeNum
        LSUpperPC = [0] * self.NodeNum
        LSUpperELV = [0] * self.NodeNum

        triNoSamp = [0] * self.NodeNum
        triNoSampPC = [0] * self.NodeNum
        triNoSampELV = [0] * self.NodeNum

        cPrime, cPCPrime, cELVPrime = 0, 0, 0
        for vid in range(self.NodeNum):
            LSUpper[vid] = verDegPrime[vid] + 2 / self.alphaVal / self.varEpsilon1 * np.log((self.h + 1) / self.deltaVal1 + self.r4 * kStarPC * self.prob * (1 - self.prob))
            LSUpperPC[vid] = verPCDegPrime[vid] + 2 / self.alphaVal / self.varEpsilon1 * np.log((self.h + 1) / self.deltaVal1 + self.r4 * kStarPC * self.prob * (1 - self.prob))
            LSUpperELV[vid] = verELVDegPrime[vid] + 2 / self.alphaVal / self.varEpsilon1 * np.log((self.h + 1) / self.deltaVal1 + self.r4 * kStarELV * self.prob * (1 - self.prob))

            triNoSamp[vid] = self.sampleTriNo(vid, self.triSets[vid]) + np.round(np.random.laplace(0, LSUpper[vid]), 3); cPrime += triNoSamp[vid]
            triNoSampPC[vid] = self.sampleTriNo[vid, self.triSetsPC[vid]] + np.round(np.random.laplace(0, LSUpperPC[vid]), 3) + (self.triSets[vid] - self.triSetsPC[vid]); cPCPrime += triNoSampPC[vid]
            triNoSampELV[vid] = self.sampleTriNo[vid, self.triSetsELV[vid]] + np.round(np.random.laplace(0, LSUpperELV[vid]), 3) + (self.triSets[vid] - self.triSetsELV[vid]); cELVPrime += triNoSampELV[vid]


        cT = np.sum(self.triSets) / 3

        cPrime = cPrime / 3 / self.prob
        cPCPrime = cPCPrime / 3 / self.prob
        cELVPrime = cELVPrime / 3 / self.prob


        MREValPC = round(abs(cPCPrime - cT) / cT, 3)
        MREValELV = round(abs(cELVPrime - cT) / cT, 3)

        return MREValPC, MREValELV
