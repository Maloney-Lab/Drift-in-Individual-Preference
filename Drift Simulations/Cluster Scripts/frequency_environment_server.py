
import sys
import os
sys.path.insert(0, '../Scripts/')

import simpledrift as sd
import matrixmaker as mm
import pandas as pd
import numpy as np

i = os.getenv('SLURM_ARRAY_TASK_ID')
# if i!>0:
# i=1
print(i)


# fbands=1/2**np.arange(1,11, 1)[:0:-1]
# per=1/fbands
# print(per)
# matureagematch=np.int64(per[-2:0:-1])


# Version for September 2022
# fbands=1/2**np.arange(1,10.5, .25)[:0:-1]
# linear version for testing on 12-8
fbands=1/np.arange(0,72, 2)[:0:-1]

per=1/fbands
print(per)
matureagematch=np.int64(per[-2:0:-1])
# print(matureagematch)
phasearray_size=fbands.shape[0]
envimeanvariance=np.linspace(0, .3, phasearray_size)
# #simple test
# fbands=1/2**np.arange(1,10, 2)[:0:-1]
# per=1/fbands
# print(per)
# matureagematch=np.int64(per[-2:0:-1])
# print(matureagematch)

savepath="../Results/FreqxEnv_Run7_lin/"
if not os.path.exists(savepath):
  os.makedirs(savepath)

mm.frequency_phaseplane(matureage=[10],
fbands=fbands,
numberofdays=1001, bhmax=.15, driftmax=.05,
strategy_resolution=11,
numberofbins=300,
envivariance=[.2],
filename_prefix="FreqvsEnvironmentMean_1001days_x3ba_.2ev_R"+str(i),
envimeanvariance=envimeanvariance,
savepath=savepath,
savealldays=False,
birthrate=[20, 40, 60],
i=i
)