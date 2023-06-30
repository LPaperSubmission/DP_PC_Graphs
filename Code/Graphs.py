from publicFus import *
########################################################################################################################
#################### For graph reading ####################
class Graphs(object):
    def __init__(self, dataFileName, pickleFilename):
        self.dataFileName = dataFileName
        self.pickleFilename = pickleFilename

    def get_vertex_id(self, vid, num, verOriTags):
        if verOriTags[vid] < 0:
            num = num + 1
            verOriTags[vid] = num
        return verOriTags[vid], num

    def readGraphData(self):
        check_file_exists(self.pickleFilename, rmTag=True)  # if file exists, remove
        #################### read graph edges ####################
        # print(self.dataFileName)
        if self.dataFileName[-3:] == 'csv':
            Edges = pd.read_csv(self.dataFileName)
        else:
            Edges = pd.read_table(self.dataFileName, sep="	", header=None)
        Edges = Edges.values  # print(type(edges), edges.shape)
        #################### re-generate neighbor set ####################
        maxID = 1 + Edges.max() # get the max id of the original graph

        num = -1  # the start id of all vertex
        edgeNum = 0  # the real edge number
        verOriTags = [-1] * maxID  ## OriIDs --> newIDs
        neiOriSet = [[] for i in range(maxID)]

        for edge in Edges:
            v_ = edge[0]
            u_ = edge[1]

            if v_ == u_: continue  # remove duplicate edges
            v, num = self.get_vertex_id(v_, num, verOriTags)
            u, num = self.get_vertex_id(u_, num, verOriTags)
            neiOriSet[v].append(u)
            neiOriSet[u].append(v)

            edgeNum = edgeNum + 1

        #################### re generate degree file and neighbour set file ####################
        num = num + 1
        nT = 0
        eNum = 0
        neiSet = [[] for i in range(num)]
        nTag = [-1] * num
        verDeg = [0] * num

        minDeg, maxDeg, sumDeg = len(neiSet), 0, 0
        for i in range(num):
            deg = 0
            for j in range(len(neiOriSet[i])):
                nid = neiOriSet[i][j]
                if nTag[nid] != nT and nid != i:
                    neiSet[i].append(nid)
                    nTag[nid] = nT
                    verDeg[i] = verDeg[i] + 1
                    deg += 1
                    eNum = eNum + 1
            sumDeg += deg
            if deg > maxDeg: maxDeg = deg
            if deg < minDeg: minDeg = deg
            nT = nT + 1
        ########## save files ##########
        dataInfo = {'verDeg': verDeg, 'neiSet': neiSet, 'verOriTags': verOriTags, 'edgeNum': int(eNum/2), 'nodeNum': num, 'minDeg': minDeg, 'maxDeg': maxDeg, 'avgDeg': sumDeg/len(verDeg)}
        savePickle(dataSet = dataInfo, filename = self.pickleFilename)