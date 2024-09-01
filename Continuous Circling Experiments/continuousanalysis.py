from scipy.io import loadmat
import pandas as pd
import glob
import numpy as np
import xarray as xr
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import lombscargle
from astropy.timeseries import LombScargle
from joblib import Parallel, delayed
import astropy.units as u
import re
import os
from joblib import Parallel, delayed


def test():
  print("test")

def showpsd(angle_dataset_all, batch=0, fly=0, freq=False, freqmode="Standard",
 fres=10000, normalization='psd', 
 a_bool=True, x_bool=False, y_bool=False,
 overwrite=False,
 showGraph=True, saveResult=False):
  # t=angle_dataset_all[batch]['timestamps']*u.second
  # t=t-t[0]
  if saveResult:
    filename="PowerSpectrum/Power_b"+str(batch)+"_f"+str(fly)+".nc"

    if glob.glob(filename) and overwrite:
      print("Skipping, file already exists")
      return
    else:
      print("Working on "+filename)

  ang=angle_dataset_all[batch].Angles[fly].dropna(dim='t')

  if ang.shape[0]<100:
    print("Aborting: Insufficient data for "+filename)
    return

  t=angle_dataset_all[batch].timestamps[(~np.isnan(angle_dataset_all[batch]['Angles']))[fly]]*u.second
  # t=t[ang['t']]

  numpowers=sum([a_bool, x_bool, y_bool])
  names=["Foo", "Foo", "Foo"]

  if a_bool:
    ls1=LombScargle(t, ang.values)
    names[0]="Angle"
  elif y_bool:
    y=np.cos(ang)
    ls1=LombScargle(t, y.values)
    names[0]='y'
  elif x_bool:
    x=np.sin(ang)
    ls1=LombScargle(t, x.values)
    names[0]='x'
  else:
    print("No measure selected to calculate power")
    return

  if freqmode=="Auto":
    freq, power=ls1.autopower(normalization=normalization)
  elif (freqmode=="Standard"):
    print("Standard Mode")
    freq=calcfreq(fres=fres)
    power=ls1.power(freq, normalization=normalization)
  elif (freqmode=="freq"):
    power=ls1.power(freq, normalization=normalization)

  power=power[:,np.newaxis].repeat(numpowers, axis=1)


  if numpowers>1: 
    if a_bool and y_bool:
      ls2=LombScargle(t, np.cos(ang).values)
      names[1]='y'
    else:
      ls2=LombScargle(t, np.sin(ang).values)
      names[1]='x'
    # else:
    #   ls2=LombScargle(t, np.sin(ang).values)
    #   names[1]='x'
    power[:,1]=ls2.power(freq,normalization=normalization)
    
  if numpowers==3:
    ls3=LombScargle(t, np.sin(ang).values)
    names[2]='x'
    power[:,3]=ls3.power(freq,normalization=normalization)  


  power_df=pd.DataFrame(power, columns=names[0:numpowers], index=freq)
  # powerall=pd.DataFrame(power1, index=freq, columns=[name[0]]])

  # for i in range(numpowers-1):
  #   powerx=pd.DataFrame(power1, index=freq, columns=[name[0]]])
  #   pd.concat()



  # print(power.shape, powerx.shape, powery.shape)
  # print([power[:, np.newaxis], powerx[:, np.newaxis], powery[:, np.newaxis]])
  # powerall=np.concatenate([power[:, np.newaxis], powerx[:, np.newaxis], powery[:, np.newaxis]], axis=1)

  if showGraph:
    sns.lineplot(data=power_df)

  if saveResult:
    power_df.to_parquet(filename)
  else: 
    return power_df
  # print("Done")
  return 


# def makeparallelindex(angle_dataset_all):

def calcfreq(fmin=(1/60/60/60/24), fmax=6, fres=10000, log=False):
  if log:
    return np.e**np.linspace(np.log(fmin), np.log(fmax), fres)
  else:
    return np.linspace(fmin, fmax, fres)

  

def circulardistance(angle_vector, t):
  dist = np.pi - abs(np.pi -abs(angle_vector[1:]-angle_vector[0:-1]))
  tmid = (t[1:]+t[0:-1])/2
  return dist, tmid


def windowaverage(anglevector, t, window_size=10):
  #The goal here is to either lowpass or highpass with a sliding window
  lowpass=np.zeros(anglevector.shape)
  highpass=np.zeros(anglevector.shape)
  x=np.cos(anglevector)
  y=np.sin(anglevector)

  for i in range(anglevector.shape[0]):
    t_bool=abs(t-t[i])>window_size
    lowpass[i]=np.arctan(np.nanmean(y[t_bool])/np.nanmean(x[t_bool]))
    # highpass[i]=anglevector[i]-np.arctan(np.nanmean(y[~t_bool])/np.nanmean(x[~t_bool]))
    highpass[i]=anglevector[i]-lowpass[i]
  
  return lowpass, highpass


def ywindowaverage(y, t, startindex=0, endindex=-1, window_size=10):

  if endindex==-1 or endindex >y.shape[0]:
    endindex=y.shape[0]
  # x=np.cos(anglevector)
  # y=np.sin(anglevector)



  y=y[startindex:endindex]
  t=t[startindex:endindex]
  lowpass=np.zeros(y.shape)
  highpass=np.zeros(y.shape)


  for i in range(y.shape[0]):
    t_bool=abs(t-t[i])>window_size*u.second
    lowpass[i]=np.nanmean(y[t_bool])
    # highpass[i]=anglevector[i]-np.arctan(np.nanmean(y[~t_bool])/np.nanmean(x[~t_bool]))
    highpass[i]=y[i]-lowpass[i]
  df=pd.DataFrame({"Lowpass":lowpass, "Highpass":highpass}, index=pd.Index(t, name="time (s)"))

  sns.lineplot(data=df)
  return df
  
# def movingaveragefilter(i, t, anglevector):
#   return

def loadall_saved(freq=False, SearchPath="PowerSpectrum/Power*", drop_a=True, drop_x=True):
  if type(freq)==bool:
    freq=calcfreq()
    freq_series=pd.Series(freq, name="Freq")

  a=glob.glob(SearchPath)
  for file_index, file in enumerate(a):
    # print(a[file_index])
    b=pd.read_parquet(a[file_index])
    # b=pd.concat(b, freqpd)
    b=pd.concat([b,freq_series],axis=1)
    # b=b[1:]
    if drop_a:
      b.drop("A", axis=1, inplace=True)
    if drop_x:
      b.drop("x", axis=1, inplace=True)
    
    c=b.melt(id_vars="Freq")
    c["batch"]=int(re.search(pattern="_b([0-9])", string=a[file_index]).group(1))
    c["fly"]=int(re.search(pattern="_f([0-9]*)", string=a[file_index]).group(1))
    if file_index==0:
      df_all=c
    else:
      df_all=pd.concat([df_all, c])

  df_all=df_all.reset_index().drop("index", axis=1)
  return df_all

from random import choices
def rolling_bootstrap(window):
  return np.mean(choices(window.values, k=100))

def make_rollingaverages(angle_dataset_all, batch=0, fly=0,subsample='1min', window='1D', bootstrap_window='1H'):
  ang=angle_dataset_all[batch]['Angles'].sel(fly=fly)
  y=np.sin(ang)
  t=angle_dataset_all[batch]["timestamps"]
  t=pd.to_datetime(t, origin="unix", unit='s')
  t-=t[0]

  # y_randomization=

  df=pd.DataFrame(data=y, columns=pd.Index(['y'], name="Measure"), index=pd.Index(t, name="Time(s)")).sort_index()
  
  smoothed=df.resample("1min").mean()
  # scrambled_t=smoothed.sample(frac=1)
  

  print(smoothed)

  rolling=smoothed.rolling(window=bootstrap_window, center=True, min_periods=10**3).mean()

  nsample=100
  ndp=smoothed.dropna().values.shape[0]
  bootstrap= np.zeros([nsample, ndp])

  index=(~np.isnan(smoothed.values)).nonzero()
  print(index)
  print(smoothed)
  for i in range(nsample):
    smoothed["Bootstrap"]=np.empty(smoothed.shape)
    smoothed["Bootstrap"].iloc[index]=smoothed.dropna().sample(frac=1).values

    rolling["Bootstrap_Smoothed"+str(i)]=df["Bootstrap"].rolling(window=bootstrap_window, center=True, min_periods=10**1).mean()

  df["Error"]=df["y"].rolling(window=bootstrap_window, center=True, min_periods=10**1).std()

  df["Smoothed"]=smoothed

  df=df.resample("1min").mean()
  x=xr.DataArray(rolling)
  # x=xr.DataArray(df)
  x["FB_key"]="B"+str(batch)+"F"+str(fly)
  # x["Fly"]=fly
  # x["Batch"]=batch
  # x=x.expand_dims("Fly")
  # x=x.expand_dims("Batch")
  x=x.expand_dims("FB_key")
  return rolling

def multiple_rollingaverages(angle_dataset_all, batch=-1, fly=-1, subsample='1min', window='1D'):
  if type(batch)== int and batch==-1:
    batch=np.arange(len(angle_dataset_all))
  else:
    batch=np.array(batch, ndmin=1)

  if type(fly) == int and fly == -1:
    fly= np.arange(angle_dataset_all[0].Angles['fly'].shape[0])
  else:
    fly=np.array(fly, ndmin=1)
  print(batch)
  print(fly)
  index=pd.MultiIndex.from_product([batch, fly])

  df_all=[]
  for i,ii in index:
    print(f"B{i}F{ii}")
    # print("B%dF%d", i, ii)
    df_all.append(make_rollingaverages(angle_dataset_all, i, ii,subsample=subsample, window=window))

  # df_all=Parallel(n_jobs=-1, verbose=5)(delayed(make_rollingaverages)(angle_dataset_all, i[0], i[1]) for i in index)
  # return df_all
  return xr.combine_by_coords(df_all)

# def
def mapfunc(a):
  a=a.dropna(dim="t")
  t=a["Timestamps"].squeeze()
  
  y=a["Angles"].squeeze()
  # print(t)
  # print(y)
  ls=LombScargle(t,y)
  # f,p=ls.autopower()
  f=ca.calcfreq()
  p=ls.power(f)

  freq=xr.DataArray(p.reshape([p.shape[0],1,1,1]),  dims=["Freq", "Fly", "Batch", "Trial"], coords={"Fly":a["Fly"],"Batch":a["Batch"],"Trial":a["Trial"], "Freq":f}, name="Power")
  # freq=xr.DataArray(p, dims=["Freq"], coords={"Fly":a["Fly"],"Batch":a["Batch"],"Trial":a["Trial"], "Freq":f}, name="Power")
  # freq=freq.expand_dims(["Fly", "Batch", "Trial"])

  # print(freq)
  # a["freq"]=f
  a=a.merge(freq, compat='override')
  a=a.chunk()
  # a["power"]=p

  return a

# def 