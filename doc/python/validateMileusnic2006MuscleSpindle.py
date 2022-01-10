import matplotlib.pyplot as plt
import pandas as pd

def formatAxis(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

def rampTest(inputDataDir, modelOutputDir, figuresOutputDir):
    nrows = 5
    ncols = 3
    fig, ax = plt.subplots(nrows=nrows,ncols=ncols, sharex=True, sharey=False, figsize=(8, 6), dpi=80)
    rampDir = inputDataDir + "/mileusnic2006/ramp/"
    dataList = pd.read_csv(rampDir+"description.txt")
    dataList.sort_values("ramping_velocity", inplace=True)
    enum = 0
    vels = [0.11, 0.66, 1.55]
    fusimotors = [(0.0, 0.0), (70.0, 0.0), (0.0, 70.)]
    for index, row in dataList.iterrows():
        data = pd.read_csv(rampDir+row['filename'], header=None)
        if(row['secondary_afferent']):
            r = nrows-2
        else:
            r = fusimotors.index((row['fusimotor_dynamic'],row['fusimotor_static'] ))
            
        axx= ax[r][vels.index(row['ramping_velocity'])]
        data.plot(x=0, y=1, kind='scatter', ax=axx, legend=False, ylabel='', s=5)
        formatAxis(axx)
        enum = enum + 1 

    outDataList = pd.read_csv(modelOutputDir + "/Mileusnic2006MuscleSpindle_rampStretches_description.txt")
    outDataList.sort_values("ramping_velocity", inplace=True)
    enum = 0
    for index, row in outDataList.iterrows():
        data = pd.read_csv(modelOutputDir+'/'+row['filename'], header=1,delimiter='\t')
        r = fusimotors.index((row['fusimotor_dynamic'],row['fusimotor_static'] ))
        axx= ax[r][vels.index(row['ramping_velocity'])]
        yy = 'primary_afferent'
        data.plot(x='Time', y=yy, kind='line', ax=axx, legend=False, ylabel='')
        formatAxis(axx)
        if row['fusimotor_dynamic'] + row['fusimotor_static'] ==0.:
            yy = 'secondary_afferent' 
            c = vels.index(row['ramping_velocity'])
            axx= ax[-2][c]
            data.plot(x='Time', y=yy, kind='line', ax=axx, legend=False, ylabel='')
            formatAxis(axx)
            yy = 'length'
            axx= ax[-1][c]
            data.plot(x='Time', y=yy, kind='line', ax=axx, legend=False, ylabel='')
            formatAxis(axx)
        enum = enum + 1 

    rows_text = [r'$\lambda_d$=0'+'\n'+r'$\lambda_s$=0', r'$\lambda_d$=70'+'\n'+r'$\lambda_s$=0', r'$\lambda_d$=0'+'\n'+r'$\lambda_s$=70', r'$\lambda_d$=0'+'\n'+r'$\lambda_s$=0']
    pad = 15 # in points
    for axx, row in zip(ax[:,0], rows_text):
        axx.annotate(row, xy=(0, 0.5), xytext=(-axx.yaxis.labelpad - pad, 0),
                xycoords=axx.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center')
        axx.set_ylabel('pps')
        axx.get_yaxis().set_label_coords(-0.2,0.5)
    ax[-1][0].set_ylabel('fascicle\nlength')
    ax[-1][0].get_yaxis().set_label_coords(-0.2,0.5)
    ax[0][1].set_title(f'Primary afferents', fontweight='semibold')
    ax[3][1].set_title(f'Secondary afferents', fontweight='semibold')

  #  ax[0][0].set_ylabel(,rotation=0, labelpad=20)
  #  ax[1][0].set_ylabel(,rotation=0, labelpad=20)
  #  ax[2][0].set_ylabel(,rotation=0, labelpad=20)
  #  ax[3][0].set_ylabel(,rotation=0, labelpad=20)
    for axx in ax[4]:
        axx.set_xlabel("time (s)")
        axx.set_xlim([0,3.3])
        axx.set_xticks([0., 3.3])
        axx.set_yticks([])
    
    ax[-1][0].set_yticks([0.95, 1.08])
    fig.subplots_adjust(bottom=0.091,
                        left=0.146,
                        right=0.985,
                        hspace=0.248,
                        wspace=0.221)
    fig.tight_layout()
    fig.savefig(figuresOutputDir+'/Mileusnic2006MuscleSpindle_rampStretches.png')
    #plt.show()


def __init__():
    rampTest('../../data/', '../../build/test/','../fig/')

__init__()