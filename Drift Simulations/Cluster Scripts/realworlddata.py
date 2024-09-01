import simpledrift as sd
import matrixmaker as mm
import pandas as pd
import numpy as np
import imp
import matplotlib.pyplot as plt
import os
import math
import datatoEnv as dtv

i = os.getenv('SLURM_ARRAY_TASK_ID')
# if i!>0:
# i=1
print('Checking array variable')
print(i)
print(type(i))
i=int(i)
print(i)


# If not running as part of a batch, then just set i=1
# if not isinstance(i, int):
#    i=1 

print(i)

#Parameters for driftvsbethedging matrix
bhlower=0
bhupper=.5
bhinterval=10
driftlower=0
driftupper=.1
driftinterval=10

#Parameters for each trial
numberofbins=100
numberofdays=100

# Parameters for generating stimuli
# Currently set to 



# lowerbound=.1
# upperbound=.5
power=0
# fsp='$SCRATCH/debivort_lab/Ryan/figs'
fsp='../Results/figs'
csvpath="AllChunks.csv"
batchsize=100
for j in range(batchsize):
  index=i*batchsize+j
  dtv.runSimulationOnline('AllChunks.csv', fsp, index)
# envimeanvariance=.3
# envivariance=.1 

# for j in range(10):
#   k=int(i*10+j)
#   readcsv=pd.read_csv(csvpath, usecols=[k+1], skiprows=[1,2])
#   f=readcsv.columns[0]

#   csvnp=np.array(readcsv)
#   env=dtv.convertTimeSeriestoEnv(csvnp, envimeanvariance=envimeanvariance,   envivariance=envivariance, length=100)


#   mm.matrixmaker(env, bhlower, bhupper, bhinterval, driftlower, driftupper,  driftinterval, runindex=k, fband=f, envimeanvariance=envimeanvariance, envivariance=envivariance)

