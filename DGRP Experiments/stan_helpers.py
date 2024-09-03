import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy as sci
import seaborn.objects as so

def plotResults(fit, poplabels="DGRP line", poplabel="Populations",Des2="SD of Drift Gaussian", Des1="SD of initial preference",
                 check=False, var1="BH", var2="D", var3="", Des3="Phi", showallbivariateplots=False, addvar1tovar2=False, 
                 var1_limits=[0,1.6], var2_limits=[.2,.6], var3_limits=[0,1],
                 dropcolumns=[]):
  # label=skt["LineNumber"]
  # numfits=fit.BH.shape[1]



  numfits=poplabels.shape[0]
  if check:
    print("Checking Diagnose")
    print(fit.diagnose())

# Check if Var3 isn't an empty string (i.e. to show Phi)
  if var3=="":
    fig, axd=plt.subplot_mosaic([['upper left', 'right'],
                                ['left', 'right']],
                                figsize=(12, 6), constrained_layout=True)
    

  else:
    if showallbivariateplots:
      fig, axd=plt.subplot_mosaic([['upper left', 'center left', 'center right', 'right'],
                                  ['left','center left', 'center right', 'right'],
                                  ['lower left','center left', 'center right',  'right']],
                                  figsize=(24, 6), constrained_layout=True)
  
    else:  
      fig, axd=plt.subplot_mosaic([['upper left', 'right'],
                                  ['left', 'right'],
                                  ['lower left', 'right']],
                                  figsize=(12, 6), constrained_layout=True)

  fd=fit.draws_pd()

  start1=var1+"["+"1"+"]"
  end1=var1+"["+str(numfits)+"]"
  D_parameters=fd.loc[:,start1:end1]

  # print(D_parameters)
  D_parameters.columns=poplabels
  for i in dropcolumns:
    D_parameters.drop(columns=D_parameters.columns[i], axis=1, inplace=True)




  dpm=D_parameters.melt()
  # print(dpm)
  dpm=dpm.rename(columns={"variable":poplabel, "value":Des1})
  # return dpm
  # print(dpm)
  d_violin=sns.violinplot(data=dpm, x=poplabel, y=Des1, ax=axd["upper left"])
  d_violin.set_ylim(var1_limits)


  start2=var2+"["+"1"+"]"
  end2=var2+"["+str(numfits)+"]"

  BH_parameters=fd.loc[:,start2:end2]
  BH_parameters.columns=poplabels
  for i in dropcolumns:
    print("Droping column "+str(i)+" from BH_parameters")
    BH_parameters.drop(columns=D_parameters.columns[i], axis=1, inplace=True)

  if addvar1tovar2:
    BH_parameters=BH_parameters+D_parameters

  bhpm=BH_parameters.melt()
  bhpm=bhpm.rename(columns={"variable":poplabel, "value":Des2})
  
    

  
  bh_violin=sns.violinplot(data=bhpm, x=poplabel, y=Des2, ax=axd["left"])
  bh_violin.set_ylim(var2_limits)

  if var3!="":
    start3=var3+"["+"1"+"]"
    end3=var3+"["+str(numfits)+"]"
    Phi_parameters=fd.loc[:,start3:end3]
    Phi_parameters.columns=poplabels
    for i in dropcolumns:
      Phi_parameters.drop(columns=D_parameters.columns[i], axis=1, inplace=True)

    phipm=Phi_parameters.melt()
    phipm=phipm.rename(columns={"variable":poplabel, "value":Des3})
    phi_violin=sns.violinplot(data=phipm, x=poplabel, y=Des3, ax=axd["lower left"])
    phi_violin.set_ylim(var3_limits)



  # dpm
  # bhpm[Des2]
  if var3=="":
    bhd=pd.concat([dpm, bhpm[Des2]], axis=1)
  else:
    bhd=pd.concat([dpm, bhpm[Des2], phipm[Des3]], axis=1)
  bhd=bhd.astype({poplabel: "category"})
  bhd_r=bhd.sample(frac=1)


  if showallbivariateplots:
    sp=sns.kdeplot(data=bhd_r, x=Des2, y=Des1, hue=poplabel,
    ax=axd["center left"])
    sp.set_xlim(var2_limits)
    sp.set_ylim(var1_limits)
    
    sp1=sns.kdeplot(data=bhd_r, x=Des2, y=Des3, hue=poplabel,
    ax=axd["center right"])
    sp1.set_xlim(var2_limits)
    sp1.set_ylim(var3_limits)


    sp2=sns.kdeplot(data=bhd_r, x=Des1, y=Des3, hue=poplabel,
    ax=axd["right"])
    sp2.set_xlim(var1_limits)
    sp2.set_ylim(var3_limits)

  else:
    sp=sns.kdeplot(data=bhd_r, x=Des1, y=Des2, hue=poplabel,
    ax=axd["right"],
  )


  return fig, fd, sp

def makedictfrompandastable(table):
  return
  #Expects longform data with Index, BatchNumber, Tray Number, FlyID, Sex, uID, Day, n, r, x

# def flattenfits

def toydatagenerator(popid='Test', days=[1,2,3, 8,9,10, 15, 16,17], nflies=100, bh=.1, drift=.01, bound=.0001):
  index=["uID", "PopID", "Day", "x", "nTurns", "rbias"]
  # uid=0
  # r=.5
  firstrun=True

  nres=10000
  px=np.linspace(0,1,nres)

  for uid in range(nflies):
    alive=True
    r=min(max(np.random.normal(.5, bh), 0),1)

    for d in range(max(days)):
      
      if d>1:
        # print(px)
        jd=sci.stats.norm.pdf(px,.5, bound)*sci.stats.norm.pdf(px,r, drift)
        # print(a)
        jd=jd/sum(jd)
        jdc=jd.cumsum()
        r=np.searchsorted(jdc, np.random.random())/nres
        # print(r[0])
        # r=.5
        if alive:
          if np.random.rand()>.95:
            alive=False
      
      n=max(0, int(sci.stats.beta.rvs( 2, 3)*1200)-100)
      
      if d<2:
        n=max(50,n)



      if not alive:
        n=0

      if n>0:
        x=np.random.binomial(n,r)
      else:
        x=np.nan

      # print(x,n)
      if (np.array(days)==d).any():
        if firstrun:
          a=pd.Series([uid, popid, d, x, n, r], index=index)
          firstrun=False
        else:
          b=pd.Series([uid, popid, d, x, n, r], index=index)
          a=pd.concat([a,b], axis=1)


  testtable=a.T.reset_index(drop=True)
  testtable["uID"]=pd.Categorical(testtable["uID"])
  testtable["Day"]=np.array(testtable["Day"], dtype=int)
  testtable["rbias"]=np.array(testtable["rbias"], dtype=float)
  # return
  plot=so.Plot(data=testtable, x="Day", y="rbias", group="uID").add(so.Line()).limit(y=(0,1))
  plot.show()
  vplot=sns.violinplot(data=testtable, x="Day", y="rbias")
  vplot

  return testtable

# def randomnselector(nmean)

def populateStanDictionary(datatable, imputex=False, poplabel="Drug Treatment", nturnadjustment=1, timeiadjustment=1):
  datatable=datatable.sort_values(["uID", "Day"])
  if imputex:
    datatable["x"]=np.array(((datatable["nTurns"]+1)*datatable["rBias"]), dtype=int)

  datatable["uID"]=np.array(datatable["uID"], dtype=int)

  # datatable
  datatable_mis=datatable[datatable.isna().any(axis=1)]
  datatable_obs=datatable.dropna()
  
  dummyrow=datatable_mis.iloc[0,:]
  
  actualdays=np.unique(datatable_mis["Day"])
  fullrange=np.arange(0,np.max(actualdays), dtype=float)+1

  missingdays=np.array(pd.Series(fullrange).dropna(), dtype=int)

  # for ii, i in enumerate(actualdays):
  #   for jj, j in enumerate(fullrange):
  #     if j==i:
  #       # print(j)
  #       # print(i)
  #       fullrange[jj]=np.nan

  # np.unique(datatable_obs["Day"])
  for i in range(datatable["uID"].unique().shape[0]):
    dummyrownew=dummyrow.copy()
    dummyrownew.loc['uID']=int(i)
    for j in missingdays:
      dummyrownew.loc['Day']=int(j)
      datatable_mis=pd.concat([datatable_mis, pd.DataFrame(dummyrownew).T])
      # skt_unmelted_mis_b.append(dummyrow)
  
  datatable_u=datatable.drop_duplicates(subset="uID").sort_values("uID")

  datadict={
    "N": datatable["uID"].unique().shape[0], # skt.shape[0], #Total number of flies
    "S": datatable[poplabel].unique().shape[0], #Total number of Lines

    "T": int(datatable["Day"].max()+timeiadjustment), #Total number of Days. Changed this to 2 as hack, but I forget why (1 is because stan is 1 indexing). Maybe unmelted mis is higher?
    "s": pd.Categorical(datatable_u[poplabel]).codes+1, # list of line number per fly
    "fly_i_obs": np.array(datatable_obs["uID"]+1), # indexes of fly number for flattened observations
    "time_i_obs": np.array(datatable_obs["Day"], dtype=int)+timeiadjustment, # indexes of fly number for flattened observations
    "t_obs": datatable_obs['nTurns'].shape[0], # total number of observed trials
    "x_obs": np.array(datatable_obs['x'], dtype=int), # flattened array with number of right turns for each obs
    "n_obs": np.array(datatable_obs['nTurns']+nturnadjustment, dtype=int), # flattened array with number of total turns for each obs
    "fly_i_mis": np.array(datatable_mis["uID"], dtype=int)+1, # indexes of fly number for flattened observations
    "time_i_mis": np.array(datatable_mis["Day"], dtype=int)+timeiadjustment, # indexes of day number for flattened observations
    "t_mis": datatable_mis['nTurns'].shape[0] # total number of observed trials
    # "x_mis": np.array(skt_unmelted['x'], dtype=int), # flattened array with number of right turns for each obs
    # "n_mis": np.array(skt_unmelted['n'], dtype=int) # flattened array with number of total turns for each obs
  }
  if datadict["time_i_mis"].min()<1:
    print("error, time_i_mis<1, will cause issues")

  return datadict, datatable_obs

def teststandata(datadict):
  R=np.zeros([datadict["N"], datadict["T"]])
  R=R*np.nan
  R
  for i in np.arange(0,datadict['t_obs']):
      # R[datadict['fly_i_obs'][i]-1, datadict["time_i_obs"][i]-1]=True
      R[datadict['fly_i_obs'][i]-1, datadict["time_i_obs"][i]-1]=(datadict['x_obs'][i])/datadict['n_obs'][i]

  if R[:,0].all():
    print("All Day 1 are fixed")
  else:
    print("Nan in Day 1")
  return R

