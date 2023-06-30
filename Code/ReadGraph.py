from Graphs import *
########################################################################################################################
def _Graph_(dataSetName, dataResultFolder):
    print("========================================================================================================================")
    #################### graph files ####################
    dataInfoFile = dataResultFolder + dataSetName + "/dataInfo.pickle"
    #################### get graph information ####################
    if not check_file_exists(dataInfoFile):
        ######################################## Generate neighbour set and degree set ########################################
        graphName = '/Graph.txt'
        if dataSetName in ['WIKIChameleon', 'WIKICrocodile', 'WIKISquirrel', 'Dolphins', 'WIKIVote', 'USAirport', 'Bitcoin', 'SisterCity', 'Yeast', 'Celegans']:
            graphName = '/Graph.csv'
        filename = dataResultFolder + dataSetName + graphName
        graphObj = Graphs(filename, dataInfoFile)
        graphObj.readGraphData()
        del graphObj
    #################### reading graph data information ####################
    with open(dataInfoFile, 'rb') as fval:
        fdata = pickle.load(fval)
        verDeg = fdata['verDeg']
        neiSet = fdata['neiSet']
        edgeNum = fdata['edgeNum']
        verOriTags = fdata['verOriTags']
        avgDeg = fdata['avgDeg']
        del fdata
    print("==================== The processing graph is as follows ====================")
    print("Details of {} is as follows: \n Nodes \t {} \n Edges \t {} (or {}) \n average deg \t {:.2f} \n max deg \t {}".format(dataSetName, len(verDeg), int(np.sum(verDeg) / 2), edgeNum, avgDeg, max(verDeg)))

    return verDeg, neiSet, verOriTags
