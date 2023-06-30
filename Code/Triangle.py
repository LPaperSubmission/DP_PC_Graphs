from publicFus import *
########################################################################################################################
######################################## For triangle counting ########################################
class Triangle(object):
    def __init__(self, **kwargs):
        # randomState, verDeg, neiSet, varEpsilon1, varEpsilon2, deltaVal, hPrime, openCheck, noiseCal,
        # triInfoFile, InitCandidates, rt, dpmechanism, ResultsFolder):
        ########## Public data info ##########
        self.randomstate = getVals('randomState', **kwargs); #kwargs['randomState']
        self.verDeg = getVals('verDeg', **kwargs); #kwargs['verDeg']
        self.NodeNum = len(self.verDeg)
        self.neiSet = getVals('neiSet', **kwargs); #kwargs['neiSet']
        self.varEpsilon1 = getVals('varEpsilon1', **kwargs); #kwargs['varEpsilon1']
        self.varEpsilon2 = getVals('varEpsilon2', **kwargs); #kwargs['varEpsilon2']
        self.deltaVal = getVals('deltaVal', **kwargs); #kwargs['deltaVal']
        self.hPrime = getVals('hPrime', **kwargs); #kwargs['hPrime']
        self.openCheck = getVals('openCheck', **kwargs); #kwargs['openCheck']
        self.noiseCal = getVals('noiseCal', **kwargs); #kwargs['noiseCal']
        self.rt = getVals('rt', **kwargs); #kwargs['rt']
        self.dpmechanism = getVals('dpmechanism', **kwargs); #kwargs['dpmechanism']
        self.ResultsFolder = getVals('ResultsFolder', **kwargs); #kwargs['ResultsFolder']
        ########## triangle & wedge Info ##########
        self.triInfoFile = getVals('dataInfoFile', **kwargs); #kwargs['dataInfoFile']
        self.lambdaELVFile = getVals('lambdaELVFile', **kwargs)
        ########## Wedge
        self.cSets = None
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
        if queryID != curSet[0]:
            sys.exit('!!!!!!!!!!!!!!!!!!!! Wrong curSet !!!!!!!!!!!!!!!!!!!!')
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

    def Phase1Triangle(self, hInput = None, PCSetTag = False, localViewSetTag = False, parallelTag = False, dataFolder = None):
        self.PCSetsTag = PCSetTag
        self.ELVSetsTag = localViewSetTag
        self.resultFolder = dataFolder
        self.parallelTag = parallelTag
        ######################################## generate max common neighbors for all nodes ########################################
        if check_file_exists(self.triInfoFile): self.readTriangleInfo()
        else:
            self.CountTriangles()
            self.CountWeges()
            dataInfo = {'cSets': self.cSets, 'cSetsPC': self.cSetsPC, 'cpSetsBar': self.cpSetsBar, 'cSetsELV': self.cSetsELV,
                        'ceSetsBar': self.ceSetsBar, 'triSets': self.triSets, 'triSetsPC': self.triSetsPC,
                        'triSetsELV': self.triSetsELV}
            savePickle(dataInfo, self.triInfoFile)
        ######################################## Phase 1 ########################################
        lambdaD = 2 / (0.5 * self.varEpsilon1)
        deltaPrime = self.deltaVal / (2 * self.hPrime + 2)
        logDelta = np.log(1 / (2 * deltaPrime))
        ########## id and degrees ##########
        LapLambdaD = np.round(self.randomstate.laplace(0, lambdaD, self.NodeNum), 3)
        ########## get sorted d_upperband ##########
        verDegUpper = [self.verDeg[i] + LapLambdaD[i] + lambdaD * logDelta for i in range(self.NodeNum)]
        pcDegSet, elvDegSet, verDegPCUpper, verDegELVUpper = None,None,None,None
        if self.PCSetsTag:
            pcDegSet = [self.readPCELV(i, ReadTag = 'PC') for i in range(self.NodeNum)]
            elvDegSet = [self.readPCELV(i, ReadTag = 'ELV') for i in range(self.NodeNum)]

            assert np.array(elvDegSet).all() == np.array(self.verDeg).all()

            verDegPCUpper = [pcDegSet[i] + LapLambdaD[i] + lambdaD * logDelta for i in range(self.NodeNum)]
            verDegELVUpper = [elvDegSet[i] + LapLambdaD[i] + lambdaD * logDelta for i in range(self.NodeNum)]
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
        self.lambdaELVFile = self.lambdaELVFile + '_h' + str(self.h) + '_var1' + str(round(self.varEpsilon1, 2)) + '_var2' + str(round(self.varEpsilon2,2)) + '.pickle'
        if self.openCheck: print('==================== h = {:d} ===================='.format(self.h))
        ########## get S ##########
        S, SPC, SELV = [], [], []
        for i in range(self.h+1):
            S.append(int(sortedDegIDSet[i, 1]))
            if self.PCSetsTag:
                SPC.append(int(sortedDegIDSetPC[i, 1]))
                SELV.append(int(sortedDegIDSetELV[i, 1]))
        ##########
        lambdaC = self.h / (0.5 * self.varEpsilon1)
        cUpperSetPC = [0] *self.NodeNum
        cUpperSetELV = [0] * self.NodeNum
        ########## laplace noise generation ##########
        LapLambdaC = np.round(self.randomstate.laplace(0, lambdaC, len(S)), 3)

        if not self.PCSetsTag:
            SPC = S.copy()
            SELV = S.copy()
            cUpperSetPCPart = [self.cSets[S[sIdx]] + LapLambdaC[sIdx] + lambdaC * logDelta for sIdx in range(len(S))]
            truth = [self.cSets[vid] for vid in S]
            # ==================== call function ====================
            if self.noiseCal: cUpperSetPCPart = noiseCalibrate(cUpperSetPCPart.copy(), truth.copy())
            for sIdx in range(len(S)):
                cUpperSetPC[S[sIdx]] = cUpperSetPCPart[sIdx]
            cUpperSetELV = cUpperSetPC.copy()
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
            ########## noised wedge + left wedges ##########
            for sIdx in range(len(S)):
                cUpperSetPC[SPC[sIdx]] = cUpperSetPCPart[sIdx] + self.cpSetsBar[SPC[sIdx]]
                cUpperSetELV[SELV[sIdx]] = cUpperSetELVPart[sIdx] + self.ceSetsBar[SELV[sIdx]]
        ####################
        cUpperSetFinalPC = [verDegPCUpper[vid] if vid not in SPC else min(cUpperSetPC[vid], verDegPCUpper[vid]) for vid in range(self.NodeNum)]
        cUpperSetFinalELV = [verDegELVUpper[vid] if vid not in SELV else min(cUpperSetELV[vid], verDegELVUpper[vid]) for vid in range(self.NodeNum)]

        self.upperBoundPC = max(cUpperSetFinalPC)
        self.upperBoundELV = max(cUpperSetFinalELV)

        if self.upperBoundPC == 0 or self.upperBoundELV == 0:
            sys.exit("!!!!!!!!!!!!!!!!!!!! The common neighbours of nodes in S should be larger than 0 (|S| = {:d}) !!!!!!!!!!!!!!!!!!!!".format(len(S)))
            # NOTE: The reason maybe the chosen of parameter: hprime, varEpsilon1, varEpslon2; re-run may help

        lambda_PC = (3 / self.varEpsilon2) * max(sortedDegIDSetPC[self.h + 1, 0], self.upperBoundPC)
        lambda_ELV = (3 / self.varEpsilon2) * max(sortedDegIDSetELV[self.h + 1, 0], self.upperBoundELV)
        if check_file_exists(self.lambdaELVFile):
            with open(self.lambdaELVFile, 'rb') as fval:
                fdata = pickle.load(fval)
                lambda_ELV = fdata['lambda_ELV']
        else:
            dataInfo = {'lambda_ELV': lambda_ELV}
            savePickle(dataInfo, self.lambdaELVFile)

        return lambda_PC, lambda_ELV

    def twoPhaseApproachTriangle(self, triFile = './triangle.pickle', hInput = None, PCSetTag = False, localViewSetTag = False, dataFolder = None, parallelTag = False):
        ######################################## phase 1 start ########################################
        lambdaPC, lambdaELV= self.Phase1Triangle(hInput = hInput, PCSetTag = PCSetTag, localViewSetTag = localViewSetTag, parallelTag = parallelTag, dataFolder = dataFolder)
        if self.openCheck: print("==================== The lambda returned by Phase one is (PC: {:.5f}, ELV: {:.5f}) ====================".format(lambdaPC, lambdaELV))

        ######################################## phase 2 start ########################################
        ########## add laplace noise to triangle number ##########
        lapNoisePC = np.round(self.randomstate.laplace(0, lambdaPC, self.NodeNum), 3)
        if lambdaPC == lambdaELV: lapNoiseELV = lapNoisePC.copy()
        else: lapNoiseELV = np.round(self.randomstate.laplace(0, lambdaELV, self.NodeNum), 3)

        if not PCSetTag:
            trisetsNoised = np.array(self.triSets) + lapNoisePC
            #===================== call function ====================
            trisetsNoised = noiseCalibrate(trisetsNoised.copy(), self.triSets.copy())
            MREValSet = [abs(trisetsNoised[i] - self.triSets[i]) / self.triSets[i] for i in range(self.NodeNum)]
            ##########
            MREValELV = sum(MREValSet) / self.NodeNum
            MREValPC = MREValELV

            triPCNoised = trisetsNoised.copy()
            triELVNoised = trisetsNoised.copy()
        else:
            ########## p cohesion ##########
            triPCNoised = self.triSetsPC + lapNoisePC
            triPCNoised = noiseCalibrate(triPCNoised.copy(), self.triSetsPC.copy())
            triPCNoised = np.array(triPCNoised) + np.array(self.triSets)-np.array(self.triSetsPC)

            ########## ELV ##########
            triELVNoised = self.triSetsELV + lapNoiseELV
            triELVNoised = noiseCalibrate(triELVNoised.copy(), self.triSetsELV.copy())
            triELVNoised = np.array(triELVNoised) + np.array(self.triSets)-np.array(self.triSetsELV)

            MREValSetPC = [abs(triPCNoised[i] - self.triSets[i]) / self.triSets[i] if self.verDeg[i] > 2 and self.triSets[i] > 0
                           else 0 for i in range(self.NodeNum)]
            MREValSetELV = [abs(triELVNoised[i] - self.triSets[i]) / self.triSets[i] if self.verDeg[i] > 2 and self.triSets[i] > 0
                           else 0 for i in range(self.NodeNum)]

            MREValELV = sum(MREValSetELV) / self.NodeNum
            MREValPC = sum(MREValSetPC) / self.NodeNum

        if self.openCheck:
            print("==================== MREPC (Mean Relative Error) = {:.2f} ({}, {}) ====================".format(MREValPC, self.varEpsilon1, self.h))
            print("==================== MREELV (Mean Relative Error) = {:.2f} ({}, {}) ====================".format(MREValELV, self.varEpsilon1, self.h))
        ######################################## save pickle file ########################################
        MREValPC = round(MREValPC, 3)
        MREValELV = round(MREValELV, 3)
        lambdaPC = round(lambdaPC, 3)
        lambdaELV = round(lambdaELV, 3)
        self.upperBoundPC = round(self.upperBoundPC, 3)
        self.upperBoundELV = round(self.upperBoundELV, 3)
        ######################################## return ########################################
        return MREValPC, self.upperBoundPC, MREValELV, self.upperBoundELV, lambdaPC, lambdaELV
