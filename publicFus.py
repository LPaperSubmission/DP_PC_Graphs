import csv

import numpy as np
np.finfo(np.dtype("float32"))
np.finfo(np.dtype("float64"))
import pickle
import os
import pandas as pd
import datetime
import math
import sys
import scipy as sp
import operator as op
from functools import reduce
from multiprocessing import Pool, cpu_count, Lock, Manager
import platform
import warnings
import sympy as sy

import matplotlib
import shutil
import matplotlib.font_manager
matplotlib.use("Agg")
matplotlib.use("PS")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
#rcParams['font.size'] = 80
rcParams['font.weight'] = 'normal'
rcParams['lines.linewidth'] = 5
rcParams['legend.columnspacing'] = 0.4
rcParams['legend.labelspacing'] = 0.4
plt.rcParams['figure.dpi'] = 600
#plt.rcParams['figure.figsize'] = [23.5,20]
rcParams['font.size'] = 76
matplotlib.rcParams.update({'font.size': 76})
plt.rcParams['figure.figsize'] = [20,18]
########################################################################################################################
#################### For public using ####################

def mymkdir(folderName):
    if not os.path.exists(folderName):
        os.makedirs(folderName)

def check_file_exists(filename, rmTag = False):
    if os.path.exists(filename):
        file_stat = os.stat(filename)
        if file_stat.st_size == 0:
            os.remove(filename)
            return False
        if rmTag:
            os.remove(filename)
            return False
        return True
    else:
        return False

def savePickle(dataSet = None, filename = None):
    with open(filename, 'wb') as fval:
        pickle.dump(dataSet, fval, pickle.HIGHEST_PROTOCOL)
######################################## Write TXT ########################################
def writeTXT(argsIn, TXTfile):
    rt = argsIn['rt']
    vare1 = argsIn['varEpsilon1']
    vare2 = argsIn['varEpsilon2']
    hprime = argsIn['hPrime']
    MREValPC = argsIn['MREValPC']
    MREValELV = argsIn['MREValELV']
    UpperPC = argsIn['UpperPC']
    UpperELV = argsIn['UpperELV']
    lambdaPC = argsIn['lambdaPC']
    lambdaELV = argsIn['lambdaELV']
    ############### save file ###############
    lineData = str(rt) + "\t" + str(round(vare1,3)) + "\t" + str(round(vare2,3)) + "\t" + str(hprime) + "\t" + str(round(MREValPC, 3)) + "\t" + str(round(MREValELV, 3)) + "\t" + str(round(UpperPC, 3)) + "\t" + str(round(UpperELV, 3)) + "\t" + str(round(lambdaPC, 3)) + "\t" + str(round(lambdaELV, 3)) + "\n"
    with open(TXTfile, 'a+') as f:
        f.write(lineData)
######################################## Calibrate noise: change negative value ########################################
def noiseCalibrate(noised, truth):
    fun = lambda n: lambda x: sp.linalg.norm(x-n, ord = 1)
    bnds = lambda v: np.full(list(v.shape) + [2], fill_value=(0, None))
    res = sp.optimize.minimize(fun(noised), truth, bounds = bnds(np.array(truth)))

    resX = np.round(res.x, 3)
    if min(resX) < 0: sys.exit('!!!!!!!!!!!!!!!!!!!! Wrong noise calibration !!!!!!!!!!!!!!!!!!!!')

    return resX
######################################## combination numer ########################################
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom  # or / in Python 2
######################################## Read p-cohesion and ELV ########################################
def ReadPCELV(vid, resultFolder, PCSetsTag = True, ELVSetsTag = True, ReadTag = 'All', InitCandidates = 'Partial', ELVFolder = None):
    pcSet, elvSet = [], []
    if PCSetsTag and ReadTag in ['All', 'pc', 'PC']:
        pcFolder = resultFolder + 'PCs/'
        PCfilename = pcFolder + 'PC_' + str(vid) + '_' + InitCandidates + '.pickle'
        with open(PCfilename, 'rb') as fval:
            fdata = pickle.load(fval)
            pcSet = fdata['pcSet']
            queryID = fdata['queryID']
            del fdata
    if ELVSetsTag and ReadTag in ['All', 'ELV']:
        ELVFolders = ELVFolder + 'ELVs/'
        ELVfilename = ELVFolders + 'ELV_' + str(vid) + '.pickle'
        with open(ELVfilename, 'rb') as fval:
            fdata = pickle.load(fval)
            elvSet = fdata['elvSet']
            queryID = fdata['queryID']
            del fdata
    return queryID, pcSet.copy(), elvSet.copy()
######################################## check if need to run the class or not ########################################
def getVals(key, **kwargs):
    if key in kwargs.keys():
        return kwargs[key]
    else:
        sys.exit('!!!!!!!!!!!!!!!!!!!! you should give {} from the input !!!!!!!!!!!!!!!!!!!!'.format(key))

def getDataInfo(lock, filename):
    headCheck = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    dataInfoDF = pd.read_csv(filename, header=None)

    dataInfoNodup = dataInfoDF.drop(dataInfoDF[(dataInfoDF.loc[:,0] == headCheck[0]) &
                                          (dataInfoDF.loc[:,1] == headCheck[1]) &
                                          (dataInfoDF.loc[:,2] == headCheck[2]) &
                                          (dataInfoDF.loc[:,3] == headCheck[3]) &
                                          (dataInfoDF.loc[:,4] == headCheck[4])].index)

    df = dataInfoNodup.round(3)
    if len(dataInfoNodup) != len(dataInfoDF):
        lock.acquire()
        df.to_csv(filename, index=False, header=None)
        lock.release()

    return df.values


def dropDFIndx(lock, rt, fileName):
    dataInfoDF = pd.read_csv(fileName, header=None)
    dataInfoNodup = dataInfoDF.drop_duplicates(keep=False)

    df = dataInfoNodup.drop(dataInfoNodup[dataInfoNodup.loc[:,3] == rt].index)
    df = df.round(3)

    lock.acquire()
    df.to_csv(fileName, index=False, header=None)
    lock.release()

def writeCSV(lock, rt, dataSet = None, columLabel = None, filename = None, randomTime = 200):
    df = pd.DataFrame(data=dataSet, columns=columLabel)
    if check_file_exists(filename) and os.path.getsize(filename) > 0:
        dataInfo = getDataInfo(lock, filename)
        if rt in dataInfo[:, 3]: return

    lock.acquire()
    df.to_csv(filename, index=False, header=None, mode='a')
    lock.release()

def checkRunTag(lock, dataCSV, dataInfoFile, rt):
    runTag = False
    if not check_file_exists(dataInfoFile, rmTag=False):
        runTag = True
    else:
        if not check_file_exists(dataCSV, rmTag=False):
            runTag = True
        else:
            dataInfo = getDataInfo(lock, dataCSV)
            if rt not in dataInfo[:, 3]:
                runTag = True
    return runTag
