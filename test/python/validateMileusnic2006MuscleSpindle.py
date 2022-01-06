import matplotlib.pyplot as plt
import pandas as pd
import sys

def formatAxis(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


def rampTest(inputDataDir, modelOutputDir):
    nrows = 5
    ncols = 3
    fig, ax = plt.subplots(nrows=nrows,ncols=ncols, sharex=True, sharey=False)
    figure1Dir = inputDataDir + "/mileusnic2006/Figure1/"
    dataList = pd.read_csv(figure1Dir+"description.txt")
    dataList.sort_values("ramping_velocity", inplace=True)
    print(dataList)
    enum = 0
    vels = [0.11, 0.66, 1.55]
    fusimotors = [(0.0, 0.0), (70.0, 0.0), (0.0, 70.)]
    for index, row in dataList.iterrows():
        data = pd.read_csv(figure1Dir+row['filename'], header=None)
        if(row['secondary_afferent']):
            r = nrows-2
        else:
            r = fusimotors.index((row['fusimotor_dynamic'],row['fusimotor_static'] ))
            
        axx= ax[r][vels.index(row['ramping_velocity'])]
        data.plot(x=0, y=1, kind='scatter', ax=axx)
        formatAxis(axx)
        enum = enum + 1 

    outDataList = pd.read_csv(modelOutputDir + "/Mileusnic2006Figure1_description.txt")
    outDataList.sort_values("ramping_velocity", inplace=True)
    enum = 0
    for index, row in outDataList.iterrows():
        data = pd.read_csv(modelOutputDir+'/'+row['filename'], header=1,delimiter='\t')
        r = fusimotors.index((row['fusimotor_dynamic'],row['fusimotor_static'] ))
        axx= ax[r][vels.index(row['ramping_velocity'])]
        yy = 'primary_afferent'
        data.plot(x='Time', y=yy, kind='line', ax=axx)
        formatAxis(axx)
        if row['fusimotor_dynamic'] + row['fusimotor_static'] ==0.:
            yy = 'secondary_afferent' 
            c = vels.index(row['ramping_velocity'])
            axx= ax[-2][c]
            data.plot(x='Time', y=yy, kind='line', ax=axx)
            formatAxis(axx)
            yy = 'length'
            axx= ax[-1][c]
            data.plot(x='Time', y=yy, kind='line', ax=axx)
            formatAxis(axx)
        enum = enum + 1 
    
    plt.show()

def __init__():
    rampTest('C:/Users/s2849511/coding/versioning/ceinms2/data', 'C:/Users/s2849511/coding/versioning/ceinms2/build/test')


def __init1__():
    if len(sys.argv) < 3:
        print("Usage: python validateMileusnic2006MuscleSpindle.py inputDataDir modelOutputDir")
    else:
        inputDataDir = sys.argv[1] 
        modelOutputDir = sys.argv[2]

__init__()