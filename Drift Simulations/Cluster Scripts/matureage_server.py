
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

mode='log'
mode='lin'

if mode=='log':
  fbands=1/2**np.arange(1,10.5, .25)[:0:-1]
elif mode=='lin':
  fbands=1/np.arange(0,72, 2)[:0:-1]

per=1/fbands
print(per)
matureagematch=np.unique(np.array(np.round(per[-1:0:-1]),dtype=int))
# print(matureagematch)

# #simple test
# fbands=1/2**np.arange(1,10, 2)[:0:-1]
# per=1/fbands
# print(per)
# matureagematch=np.int64(per[-2:0:-1])
# print(matureagematch)

savepath="../Results/MatureAge_Run7_"+mode+"_matchinterestingemv/"
if not os.path.exists(savepath):
  os.makedirs(savepath)


# mm.frequency_phaseplane(matureage=matureagematch,
# fbands=fbands,
# numberofdays=1001, bhmax=.3, driftmax=.1,
# strategy_resolution=41,
# numberofbins=100,
# envivariance=[.2],
# filename_prefix="FreqvsMatureAge_1001days_.2ev_.1emv_R"+str(i),
# envimeanvariance=[.1],
# savepath=savepath,
# savealldays=True,
# birthrate=[40],
# i=i
# )
envivariance=.125
numberofdays=1001
emv=.1

mm.frequency_phaseplane(matureage=matureagematch,
fbands=fbands,
numberofdays=numberofdays, bhmax=.15, driftmax=.05,
strategy_resolution=11,
numberofbins=300,
envivariance=[.125],
filename_prefix="FreqvsMatureAge_"+str(numberofdays)+"days_"+str(envivariance)+"ev_"+str(emv)+"emv_R",
envimeanvariance=[emv],
savepath=savepath,
savealldays=False,
birthrate=[40],
i=i
)

print('Done')