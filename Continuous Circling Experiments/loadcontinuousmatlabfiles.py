import numpy as np
import glob
import h5py
import re
import xarray as xr
import continuousanalysis as ca
from astropy.timeseries import LombScargle
from joblib import Parallel, delayed
import astropy.units as u
import pandas as pd
import time
import scipy as sci
from pandarallel import pandarallel
import seaborn.objects as so
from seaborn import axes_style



def Test():
  print('test')

def load_raw_data(path):
  centroid_path=glob.glob(path+'*centroid.bin')
  centroid=np.fromfile(centroid_path[0], dtype='single')
  c_length=centroid.shape[0]
  c_length_t=int(c_length/2/42)
  centroid=centroid.reshape([c_length_t, 2, 42])
  time_path=glob.glob(path+'*time.bin')
  time=np.fromfile(time_path[0], dtype='single')
  # l=glob.glob(path)
  # print(l)
  return centroid, time

def load_24hrmatfiles(path, num=2):
  tpath=glob.glob(path+str(num)+'TimeArray24h_u*.mat')
  apath=glob.glob(path+str(num)+'*AngleArray*.mat')
  tfile=h5py.File(tpath[0], 'r')
  # t_varname=tfile.keys()
  # print(t_varname)
  tmatrix=tfile.get('timearray')
  afile=h5py.File(apath[0], 'r')
  # a_varname=afile.keys()
  amatrix=afile.get('angle')
  return tmatrix, amatrix

def convertmattonc(path, duration=24, savepath='../Data', calcpower=True, vector=False, circumvert=False):
  print("Getting all the stuff")
  tpaths=np.sort(glob.glob(path+'*TimeArray'+str(duration)+'h_u.mat'))
  if vector:
    apaths=np.sort(glob.glob(path+'*VectorArray'+str(duration)+'h.mat'))
  else:
    apaths=np.sort(glob.glob(path+'*AngleArray'+str(duration)+'h.mat'))
  
  cpaths=np.sort(glob.glob(path+'*CentroidArray'+str(duration)+'h.mat'))
  if len(tpaths)==len(apaths) & len(tpaths)==len(cpaths):
    for i, tpath in enumerate(tpaths):
      print(i)
      batchi=int(re.search('([0-9]*)Time',tpath).group(1))
      th5=h5py.File(tpath)
      tmat=np.array(th5['timearray']).squeeze()
      print(th5['timearray'].shape)

      print(tmat.shape)
      # tmat
      if len(tmat.shape)>1:
        numtrials=tmat.shape[0]
        numt=tmat.shape[1]
      else:
        numtrials=1
        numt=tmat.shape[0]
      print("numtrials: ", numtrials)
      print("numt: ", numt)   
      tmat=np.reshape(tmat,[numtrials, numt])

      tar=xr.DataArray(tmat, 
      dims=['Trial', 't'],
      coords={"Trial":np.arange(0,numtrials), "t":np.arange(0,numt)},
      name="Timestamps")
      tar.coords["Batch"]=batchi
      tar=tar.expand_dims("Batch")
      # tsave=str(batchi)+"_Timestamps"+str(duration)+"h.nc"
      # tar.to_netcdf(tsave)
      # if i==0:
      #   tar_all=tar
      # else:
      #   tar_all=xr.concat([tar_all,tar], dim="Batch")
      
      
      ah5=h5py.File(apaths[i])
      print(apaths[i])
      if vector:
        amat=np.array(ah5['vector']).squeeze()
      else:
        amat=np.array(ah5['angle']).squeeze()

      print(amat.shape)
      if len(amat.shape)==3:
        numtrials=amat.shape[0]
        numt=amat.shape[2]
        numflies=amat.shape[1]
      elif len(amat.shape)==2:
        numtrials=1
        numt=amat.shape[1]
        numflies=amat.shape[0]
      amat=np.reshape(amat,[numtrials, numflies, numt])
      aar=xr.DataArray(amat, 
        dims=['Trial', 'Fly', 't'],
        coords={"Trial":np.arange(0,numtrials), "t":np.arange(0,numt)},
        name="Angles")
      aar.coords["Batch"]=batchi
      aar=aar.expand_dims("Batch")
      asave=str(batchi)+"_Angles_"+str(duration)+"h.nc"
      # print(aar)
      Warning("Error")
      # print(asave)
      # print(aar.shape)
      # aar.to_netcdf(asave)
      bar=aar.to_dataset().merge(tar)
      # Version without Trial
      # print("crashes right after here?")
      # print(asave)
      # bar.to_netcdf(asave)
      # return bar
    
      anglesxtimestamp=xr.concat(Parallel(n_jobs=-1, verbose=4)(delayed(removetrials)(bar, fi) for fi, f in enumerate(bar["Fly"])), dim="Fly")
      # anglesxtimestamp=bar.to_dataframe().dropna(subset=["Angles", "Timestamps"]).reset_index().drop(columns=["t", "Trial"]).set_index(["Timestamps", "Fly", "Batch"]).to_xarray()
      print(asave)
      # return anglesxtimestamp
      anglesxtimestamp.to_netcdf(asave)
      # downsamplesave=str(batchi)+"_DownsampleMeans"+str(duration)+"h.csv"
      # angledf_resampled=dailyaverage(anglesxtimestamp)
      # angledf_resampled.to_csv(downsamplesave)
      # lsarray[:,:,batchi-1]=np.array(calculatepowerall(anglesxtimestamp, batchi=1)).T
      # return anglesxtimestamp
      # return bar
      if calcpower:
        # f=ca.calcfreq(fres=1000, log=True, fmax=5, fmin=1/60/60/24/30)
        
        lsarray=calculatepowerall(anglesxtimestamp, batchi=batchi, computeshuffled=True, circumvert=circumvert)

        if vector:
          lssave=str(batchi)+"_VectorPower"+str(duration)+"h.nc"
        else:
          lssave=str(batchi)+"_Power"+str(duration)+"h.nc"
        lsarray.to_netcdf(lssave)
      

  return 


def removetrials(bar, fly):
  return bar.sel(Fly=fly).to_dataframe().dropna(subset=["Angles", "Timestamps"]).reset_index().drop(columns=["t", "Trial"]).set_index(["Timestamps", "Batch"]).to_xarray()

# def calcpower

def dailyaverage(anglexarray, circling=False):
  tic=time.perf_counter()
  anglexarray["Timestamps"]=pd.to_datetime(anglexarray["Timestamps"], origin="unix", unit='s')

  if circling:
    anglexarray["Angles"]=np.sin(anglexarray["Angles"])

  angledf=anglexarray["Angles"].to_dataframe().reset_index().set_index("Timestamps")
  angledf["Fly"]=pd.Categorical(angledf["Fly"])

  toc=time.perf_counter()
  print(toc-tic) 
  # return angledf
  angledf_pivot=angledf.pivot(columns=["Batch", "Fly"], values="Angles")
  angledf_pivot_resampled=angledf_pivot.resample('1D').mean()
  toc=time.perf_counter()
  print(toc-tic) 
  angledf_pivot_resampled_error=angledf_pivot.resample('1D').sem()
  toc=time.perf_counter()
  print(toc-tic) 
  angledf_resampled=angledf_pivot_resampled.melt(ignore_index=False, value_name="Angles").reset_index()
  angledf_resampled_error=angledf_pivot_resampled_error.melt(ignore_index=False, value_name="Angles").reset_index()
  angledf_resampled["Fly"]=pd.Categorical(angledf_resampled["Fly"])
  angledf_resampled["Angle SEM"]=angledf_resampled_error["Angles"]

  angledf_resampled["A_lb"]=angledf_resampled["Angles"]-angledf_resampled["Angle SEM"]
  angledf_resampled["A_ub"]=angledf_resampled["Angles"]+angledf_resampled["Angle SEM"]
  # return angledf_resampled
  angledf_resampled["Day"]=np.array(angledf_resampled["Timestamps"].values, dtype=int)
  angledf_resampled["Day"]=angledf_resampled["Day"]-angledf_resampled["Day"][0]
  angledf_resampled["Day"]=angledf_resampled["Day"]//86400000000000
  angledf_resampled=angledf_resampled[~np.isnan(angledf_resampled["Angles"])]
  toc=time.perf_counter()
  print(toc-tic) 
  return angledf_resampled

def calculatepowerall(anglesxtimestamp, batchi=1, computeshuffled=False, circumvert=True):

  f=ca.calcfreq(fres=100, log=True, fmax=5, fmin=1/60/60/24/30)
  # lsarray=xr.DataArray(np.zeros([f.shape[0], bar["Fly"].values.shape[0], bar["Trial"].values.shape[0],1]),
  #  coords={"Freq":f, "Fly":bar["Fly"].values, "Trial":bar["Trial"].values, "Batch":[batchi]},
  #  dims=["Freq", "Fly", "Trial", "Batch"],
  #  name="Power")

  if computeshuffled:
    lsarray=xr.DataArray(np.zeros([f.shape[0], anglesxtimestamp["Fly"].values.shape[0],1, 2]),
      coords={"Freq":f, "Fly":anglesxtimestamp["Fly"].values, "Batch":[batchi], "Shuffled":["Unshuffled", "Shuffled"]},
      dims=["Freq", "Fly","Batch", "Shuffled"],
      name="Power")
  else:
    lsarray=xr.DataArray(np.zeros([f.shape[0], anglesxtimestamp["Fly"].values.shape[0],1]),
      coords={"Freq":f, "Fly":anglesxtimestamp["Fly"].values, "Batch":[batchi]},
      dims=["Freq", "Fly","Batch"],
      name="Power")
  # print(lsarray)
  # return
  


  # for trial in bar["Trial"].values:
  output=Parallel(n_jobs=1, verbose=4)(delayed(calculatepowerforonefly)(f, anglesxtimestamp, fly, circumvert=circumvert) for fly, flynum in enumerate(anglesxtimestamp["Fly"]))
  
  if computeshuffled:
    lsarray[:,:,0, 0]=np.array(output).T
    # print("Computing Shuffled Version Now")
    output=Parallel(n_jobs=-1, verbose=4)(delayed(calculatepowerforonefly)(f, anglesxtimestamp, fly, shuffle=True, circumvert=circumvert) for fly, flynum in enumerate(anglesxtimestamp["Fly"]))
    lsarray[:,:,0, 1]=np.array(output).T

  else:
    lsarray[:,:,0]=np.array(output).T

  # lsarray["Fly"]=anglesxtimestamp["Fly"]
  return lsarray


def calculatepowerforonefly(f, anglesxtimestamp, fly, shuffle=False, circumvert=True):
  # print(anglesxtimestamp)
  y=anglesxtimestamp["angle"].squeeze()[fly,:]
  if circumvert:
    y=np.sin(y)

  # print("Power Shapes")
  # print(y.shape)
  # print(fly)
  yi=np.isnan(y) #Problem is here, wrong dimensions (42 instead of)
  # print(yi.shape)

  t=anglesxtimestamp["timestamps"][~yi]*u.second

  y=anglesxtimestamp["angle"].squeeze()[fly,~yi]
  if shuffle:
    rng = np.random.default_rng()
    rng.shuffle(y.values)

  ls=LombScargle(t=t, y=y)
  # slicesize=lsarray[:,fly,:].shape
  try:
    p=ls.power(f)
  except:
    # print(trial, fly)
    # when there isn't enough non nan flies
    p=np.empty(f.shape)
    p.fill(np.NaN)
    # return([t,y])
  # print(p.shape)
  # print(slicesize)
  # print("Calculated Power for Fly "+str(int(fly)))
  return p
  

def lombscarglelowpass(dataframe, maxsamples=-1, plotres=1000, lp=60*60*24):
  x=dataframe.dropna().index.to_numpy()
  y=dataframe.dropna()["Angles"]
  # min(x)
  yi=~np.isnan(y)
  yi
  if maxsamples<1:
    imax=x.shape[0]-1
  else:
    imax=min(x.shape[0]-1, maxsamples)
  # yi=np.arange(0,imax)
  x=np.array(x,dtype=float)

  t=x*10**-9
  t-=t[0]
  t=t*u.second
  print("Total time is", t[-1].to(u.day))
  print("Processing until", t[imax].to(u.day))

  ls=LombScargle(t[0:imax], y[0:imax], nterms=1)
  f, p=ls.autopower()
  print(imax)
  tmax=t[imax]
  print(tmax)
  t_fit=np.linspace(0, np.array(tmax), plotres)*u.second
  if lp>0:
    numlowpass=sum(f<(1/(lp*u.second)))
  else:
    numlowpass=len(f)-1

  fsum=np.zeros([t_fit.shape[0], numlowpass])
  fweighted=np.zeros([t_fit.shape[0], numlowpass])

  for i, fi in enumerate(f[0:numlowpass]):
    # print(i)
    fsum[:,i]=np.array(ls.model(t_fit, f[i]))
    fweighted[:,i]=fsum[:,i]*p[i]
  
  fweighted/=np.sum(np.array(p[0:numlowpass]))
  # fweighted=fsum*p
  pdfsum=pd.DataFrame(fsum[:,0:numlowpass])
  pdfsum.index=pd.TimedeltaIndex(np.array(t_fit), unit='s', name="Time")
  pdfsum.columns=pd.Index(f[0:numlowpass], name="Frequency (1/s)")
  pdfsum=pd.DataFrame(pdfsum.mean(axis=1), columns=["Averaged"])
  pdfsum["Weighted"]=np.sum(fweighted, axis=1)
  return pdfsum

def exponentialnan(x, sigma=24):
  mp=int(np.floor(x.shape[0]/2))
  # sigma=60*24
  win=1/(sigma*np.sqrt(2)*np.pi)*np.exp(-.5*np.abs(np.arange(x.shape[0])-mp)/sigma)**2
  win*=~np.isnan(x)
  return np.nansum(x*win/np.nansum(win))

def nandownsample(x):
  return np.nanmean(x)

def ntermslombscargle(dataframe, nterms=30, maxsamples=-1, basefreq=1/60/60/24/30):
  x=dataframe.dropna().index.to_numpy()
  y=dataframe.dropna()["Angles"]
  # min(x)
  yi=~np.isnan(y)
  yi
  if maxsamples<1:
    imax=x.shape[0]-1
  else:
    imax=min(x.shape[0]-1, maxsamples)
  # yi=np.arange(0,imax)
  x=np.array(x,dtype=float)

  t=x*10**-9
  t-=t[0]
  t=t*u.second
  print("Total time is", t[-1].to(u.day))
  print("Processing until", t[imax].to(u.day))

  ls=LombScargle(t[0:imax], y[0:imax], nterms=nterms)

  t_fit=np.linspace(0,14,1000)*u.day
  baseperiod=(1/(basefreq/u.second).to(1/u.day)).to_string()
  fit=ls.model(t_fit, basefreq/u.second)
  method_str="Fit of "+str(nterms)+" terms, Bp="+baseperiod  
  fitdf=pd.DataFrame({"Turning":fit, "Method":method_str}, index=pd.TimedeltaIndex(t_fit.value, unit='D', name="Time")).reset_index()
  # plt.plot(fit)
  return fitdf

def loadncdata(batch, hour):

  a=xr.load_dataset("1_Angles24h.nc")
  a["Fly"]
  b=xr.load_dataset("1_Timestamps24h.nc")
  c=b.merge(a)
  c["Timestamps"]
  d=c.to_dataframe().dropna(subset=["Angles", "Timestamps"]).reset_index().drop(columns=["t", "Trial"]).set_index(["Timestamps", "Fly", "Batch"]).to_xarray()

def loadallpowerdata(path="./", hour='2'):
  # powerpath=glob.glob(path+'*.Power'+hour+'h.nc')
  allpower=xr.open_mfdataset(path+'*Power'+hour+'h.nc')
  return allpower

def plotallpowers(allpower):
  # Expects an xarray with Power x hours x fly x batch
  return

def blacknan(x, threshold=.5):
  window=sci.signal.windows.blackman(x.shape[0])
  if np.nansum(window[~np.isnan(x).squeeze()]/np.sum(window))>(threshold):
    return np.nanmean(np.multiply(window,x.squeeze()/np.mean(window[~np.isnan(x).squeeze()])))
  else:
    return np.nan
  
def blacknanerror(x):
  if np.nansum(sci.signal.windows.blackman(x.shape[0])[~np.isnan(x).squeeze()])>(threshold):
    weights=sci.signal.windows.blackman(x.shape[0])*x
    rng=rng
    return np.nanmean(sci.signal.windows.blackman(x.shape[0])*x)
  else:
    return np.nan

def blacknanscramble(x,nbootstrap=1000):
  #Scramble the data to make a smoothed white noise version to show what's going on
  dfx=x.resample('60s').mean().rolling(window=3001, min_periods=1, center=True).apply(blacknan)
  rng = np.random.default_rng()
  dfx["Method"]=Blacknan
  for i in range(nbootstrap):
    rng.shuffle(x.values)

  return

def getpercentilesfromshuffleddata(shuffleddata, raw, threshold=.5):
  numshuffles=shuffleddata.shape[1]
  p50t=int(numshuffles*.5)
  p5t=int(numshuffles*.05)
  p95t=int(numshuffles*.95)
  shuffleddatasort=np.sort(shuffleddata, axis=1)
  mean=shuffleddatasort.mean(axis=1)
  p50=shuffleddatasort[:,p50t]
  p5=shuffleddatasort[:,p5t]
  p95=shuffleddatasort[:,p95t]
  resampled1h=raw.resample('1H').mean()
  percentiledf=pd.DataFrame({"Shuffled":mean, "Median":p50, "lb":p5, "ub":p95}, index=pd.Index(data=resampled1h.index, name="Hours"))
  percentiledf=percentiledf.resample('1H').mean()
  percentiledf.reset_index(inplace=True)
  

  percentiledf["Raw"]=resampled1h.values
  percentiledf["Low Pass Filtered"]=resampled1h.rolling(window=51, min_periods=10, center=True).apply(blacknan, kwargs={"threshold":threshold}).values
  
  percentiledf["Chunk"]=0
  percentiledf.loc[1:,"Chunk"]=np.cumsum(np.diff(np.isnan(percentiledf["Low Pass Filtered"])))
  percentiledf=percentiledf.drop(columns=["Median",]).melt(id_vars=["Hours", "lb", "ub", "Chunk", "Raw"])
  percentiledf.Hours=np.array(np.array(percentiledf.Hours, dtype='timedelta64[h]'), dtype="int")
  # percentiledf=percentiledf[[]]
  return percentiledf
  
def parallelshuffle(x):
  rng = np.random.default_rng()
  return rng.permutation(x).squeeze()

def parallelchunkandshuffle(x, chunksize=1000):
  x=x.squeeze()
  rng = np.random.default_rng()
  chunks=np.array_split(x, chunksize)


def parallelshuffleandfilter(x_df, threshold=.5):
  # rng = np.random.default_rng()
  x=x_df.sample(frac=1).reset_index(drop=True)
  x.index=x_df.index
  x=x.resample('1H').mean()
  return x.rolling(window=51, min_periods=10, center=True).apply(blacknan, kwargs={"threshold":threshold})

def shuffleandfilter(rawdata, numshuffles=1000, threshold=.5, parallel=True):
  time_len=rawdata.shape[0]
  shuffleddata=np.zeros([numshuffles, time_len])
  if parallel:
    njobs=-1
  else:
    njobs=1
  shuffleddata[:,:]=Parallel(n_jobs=1)(delayed(parallelshuffle)(rawdata) for i in range(numshuffles))
  
  # return shuffleddata
  shuffleddata=pd.DataFrame(shuffleddata.T, index=rawdata.index) #s index=nonanindex
  shuffleddata=shuffleddata.resample('1H').mean()

  if parallel:
    pandarallel.initialize(progress_bar=False)
    shuffleddata=shuffleddata.rolling(window=51, min_periods=10, center=True, ).parallel_apply(blacknan, kwargs={"threshold":threshold})
  else:
    shuffleddata=shuffleddata.rolling(window=51, min_periods=10, center=True, ).apply(blacknan, kwargs={"threshold":threshold})
  
  # shuffledata.

  percentiledf=getpercentilesfromshuffleddata(shuffleddata, rawdata, threshold)
  return percentiledf
  
def getallflies(xarray, numshuffles=1000,  maxflies=50, save=False):
  #get the xarray with xr.load_dataset
  flyindex=pd.MultiIndex.from_product([xarray["Batch"].values, xarray["Fly"].values])

  if save:
    Parallel(n_jobs=1, verbose=5)(delayed(getonefly)(xarray=xarray, f=f, b=b, numshuffles=numshuffles, save=save) for b, f in flyindex)
    return
  else:
    pdall=Parallel(n_jobs=4, verbose=5)(delayed(getonefly)(xarray=xarray, f=f, b=b, numshuffles=numshuffles, save=save) for b, f in flyindex)
    pdall=pd.concat(pdall)
    return pdall.reset_index()

def getonefly(xarray, f, b, numshuffles, save=False):
  # print("F",f,"B",b, end="")
  name=str('F'+str(f)+'B'+str(b)+'.csv')
  if save & len(glob.glob(name))>0:
    return

  ddf=xarray.sel(Fly=f, Batch=b).squeeze().to_dataframe().drop(columns="Batch")
  ddf.index=pd.to_datetime(ddf.index, origin="unix", unit='s')
  ddfs=np.sin(ddf)
  ddfs0=ddfs
  ddfs0.index=ddfs.index-ddfs.index[0]
  nonanraw=ddfs0.dropna()

  pdsubset=shuffleandfilter(nonanraw, numshuffles=numshuffles, parallel=False)
  pdsubset["Fly"]=f
  pdsubset["Batch"]=b
  if save:
    pdsubset.to_csv(name)
    return
  else:
    return pdsubset.dropna()

def getallfliesbatch():
  pdall=Parallel(n_jobs=-1, verbose=5)(delayed(getonefly)(xarray=xarray, f=f, numshuffles=numshuffles) for f in xarray["Fly"].values[0:maxflies])
  pdall=pd.concat(pdall)

  return

def nonchunkedmemoryhelper(list):
  for l in list:
    print(l)
    xarray=xr.load_dataset(l)
    flyindex=pd.MultiIndex.from_product([xarray["Batch"].values, xarray["Fly"].values])
    print(flyindex)
    Parallel(n_jobs=1, verbose=5, temp_folder='partemp')(delayed(getonefly)(xarray=xarray, f=f, b=b, numshuffles=1000, save=True) for b, f in flyindex)

def daskcompute(xarraychunk):
  xarraychunk=xarraychunk.dropna(dim="Timestamps").to_dataframe().drop(columns="Batch")
  xarraychunk.index=pd.to_datetime(xarraychunk.index, origin="unix", unit='s')
  xarraychunk=np.sin(xarraychunk)
  if len(xarraychunk)>0:
    xarraychunk.index=xarraychunk.index-xarraychunk.index[0]
  list=Parallel(n_jobs=-1, verbose=5)(delayed(parallelshuffleandfilter)(xarraychunk) for i in range(1000))
  percentiledata=getpercentilesfromshuffleddata(listruns, xarraychunk)

def getmeanandstdofalltimestamps(xarraydataset):
  xarraydataset=xarraydataset.to_dataframe().reset_index().drop(columns=["Trial","t"])
  xarraydataset["Timestamps"]=pd.to_datetime(xarraydataset["Timestamps"], origin="unix", unit='s')
  xarraydataset["batchxfly"]="f"+xarraydataset["Fly"].astype(str)+"b"+xarraydataset["Batch"].astype(str)
  xarraydataset=xarraydataset.drop(columns=["Fly","Batch"])
  xarraydataset["Timestamps"]=xarraydataset.groupby("batchxfly")["Timestamps"].transform(lambda x: x-x.iloc[0]).dropna()
  return xarraydataset, xarraydataset.groupby("batchxfly").agg(["mean","std", "count"])
                                                           
def import_matlab_netcdf_files(file_path):
  imported_xarray=xr.load_dataset(file_path, engine="netcdf4")
  imported_xarray['timestamp'].values=pd.to_datetime(imported_xarray['timestamp'].values.squeeze(),unit='s', origin='unix').values.reshape(1,1,max(imported_xarray['timestamp'].values.shape))
  return imported_xarray

def import_dask_matlab_files(file_path="", duration=[2, 24], batch=np.arange(1,8), fly=np.arange(1,43), override=False):
  # Loop through all the files and calculate the summary statistics. Parallize the loop

  iterables=[duration, batch, fly]
  index=pd.MultiIndex.from_product(iterables, names=["Duration", "Batch", "Fly"])
  # index

  Parallel(n_jobs=8, verbose=5)(delayed(importonefile)(file_path, i[0],i[1],i[2], override=override) for i in index)
      # return

  
# def nestedload(path):
#   xr.open_mfdataset
#   returnb

def importonefile(file_path, d, b, f, save_timeseries=False, save_summary=True, bins_and_power=True, override=False):
  #load the raw data and create a summary file that has the relevant data we're interested in
  oneflylist=glob.glob(file_path+"CirclingData_"+str(d)+"h_B"+str(b)+"_F"+str(f)+"_T*.nc")
  print(d,b, f)
  savename="Circling_Summary_"+str(d)+"h_B"+str(b)+"_F"+str(f)+".nc"
  if len(glob.glob("./Summary Data/"+savename))>0:
    if not override:
      # print(override)
      print(savename+" already exists")
      return
  # print(file_path+"CirclingData_"+str(d)+"h_B"+str(b)+"_F"+str(f)+"_T*.nc")
  # xarraysubset=[
  if len(oneflylist)>0:
    if len(oneflylist)==1:
      print("Only One trial detected for ", savename)
      xarraysubset=xr.open_mfdataset(oneflylist, combine="nested", concat_dim="t").squeeze().compute().stack({"thack":["t"]}).swap_dims({"thack":"timestamps"})
    elif len(oneflylist)>0:
    # xarraysubset=xr.open_mfdataset(oneflylist, combine="nested", concat_dim="t", chunks=-1).squeeze().compute().stack({"t_hack":["t", "Trial"]}).swap_dims({"t_hack":"timestamps"})
    # return xarraysubset
    # xarraysubset=xr.open_mfdataset(oneflylist, combine="nested", concat_dim="t", chunks=-1).squeeze()
      # return xarraysubset
    # print("Trial Testing:", d,b, f)
      xarraysubset=xr.open_mfdataset(oneflylist, combine="nested", concat_dim="t").squeeze().compute().stack({"thack":["Trial", "t"]}).swap_dims({"thack":"timestamps"})
    # return xarraysubset
    # xarraysubset=xr.open_mfdataset(oneflylist, combine="nested", concat_dim="t").compute().squeeze().swap_dims({"t":"timestamps"})
    
    xarraysubset["timestamps"]=pd.to_datetime(xarraysubset["timestamps"], unit='s', origin='unix')
    # print("Calculated datetime")
    xarraysubset=xarraysubset.sortby("timestamps")
    # print("Sortedtimestamp")

    xarraysubset["timestamps"]=xarraysubset["timestamps"]-xarraysubset["timestamps"][0]
    # print("Zeroed Timestamps")

    # xarraysubset["timestamps"].expand_dims(["Fly", "Batch"])
    # xarraysubset=xarraysubset.swap_dims({"timestamps":"t_hack"})
    # print("reswapped dims")

    xarraysubset=xarraysubset.unstack().expand_dims(["Fly", "Batch"]).drop("thack")
    # return xarraysubset
    if save_timeseries:
      # print(d, b, f)
      # print("CirclingData_"+str(d)+"h_B"+str(b)+"_F"+str(f)+".nc")
      pathname="CirclingData_"+str(d)+"h_B"+str(b)+"_F"+str(f)+".nc"
      # print(type(pathname))
      xarraysubset.to_netcdf(pathname)
      print("CirclingData_"+str(d)+"h_B"+str(b)+"_F"+str(f)+".nc")

    if save_summary:
      xarraysubset=xarrayone_processing(xarraysubset)
    if bins_and_power:
      savename="Circling_Summary_"+str(d)+"h_B"+str(b)+"_F"+str(f)+".nc"
      print(savename)
      xarraysubset=get_bins(xarraysubset, save_name=savename)
  else:
    print("Couldn't Find any files here: "+file_path+"CirclingData_"+str(d)+"h_B"+str(b)+"_F"+str(f)+"_T*.nc")
    # xarraysubset=[]
  print("End", d, b, f)
  
  # return xarraysubset #for debugging,

def xarrayone_processing(xarraysubset):
  varnames=list(xarraysubset)
  # print(varnames)
  while 'inx' in varnames: varnames.remove('inx')
  while 'iny' in varnames: varnames.remove('iny')
  if 'speed' in varnames: 
    varnames.append('speed_logged')
    xarraysubset["speed_logged"]=np.log(xarraysubset["speed"])
    # print("Speed Logged")
  
  for foobar in varnames:
    # print(foobar)
    xarraysubset[foobar+"_summary"]=xr.concat([xarraysubset[foobar].mean(dim="timestamps", skipna=True).compute().expand_dims(dim={"Measure":["Mean"]}),
    xarraysubset[foobar].std(dim="timestamps").compute().expand_dims(dim={"Measure":["Std"]}),
    xarraysubset[foobar].min(dim="timestamps").compute().expand_dims(dim={"Measure":["Min"]}),
    xarraysubset[foobar].max(dim="timestamps").compute().expand_dims(dim={"Measure":["Max"]}),
    ], dim="Measure")

    # if foobar="speed_logged":

  # for bin

  #set to time delta
  # xarraysubset["timestamps"]=xarraysubset["timestamps"]-xarraysubset["timestamps"][0]
  return xarraysubset


def get_bins(x_array, numbins=201, save_name="", savepath="./Summary Data/"):
  
  f=ca.calcfreq(fres=1000, log=True, fmax=5, fmin=1/60/60/24/30)/u.second
  t=x_array["timestamps"]

  bin_ds=xr.Dataset()
  a=x_array["direction"].plot.hist(bins=np.linspace(-np.pi, np.pi, numbins))
  # print("Got here")
  bin_ds["direction_bins"]=xr.DataArray(np.expand_dims(a[0], axis=[1,2]), coords={"Bins_direction":np.mean([a[1][0:-1], a[1][1:]], axis=0),"Fly":x_array.Fly, "Batch":x_array.Batch})
  # x_array["direction_bins"].expand_dims(dim={"Fly":x_array.Fly})
  p=calcpowerforfly(x_array["direction"],x_array["timestamps"])
  ps=calcpowerforfly(x_array["direction"],x_array["timestamps"],shuffle=True)
  bin_ds["direction_psd"]=xr.concat([xr.DataArray(np.expand_dims(p, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":["Non-shuffled Data"]}), 
  xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})], dim="Shuffled")
  
  a=x_array["angle"].plot.hist(bins=np.linspace(-1, 1, numbins))
  bin_ds["angle_bins"]=xr.DataArray(np.expand_dims(a[0], axis=[1,2]), {"Bins_angle":np.mean([a[1][0:-1], a[1][1:]], axis=0),"Fly":x_array.Fly, "Batch":x_array.Batch})

  p=calcpowerforfly(x_array["angle"],x_array["timestamps"])
  ps=calcpowerforfly(x_array["angle"],x_array["timestamps"],shuffle=True)
  bin_ds["angle_psd"]=xr.concat([xr.DataArray(np.expand_dims(p, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":["Non-shuffled Data"]}), 
  xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})], dim="Shuffled")
  # bin_ds["angle_psd"]=xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})

  a=x_array["theta"].plot.hist(bins=np.linspace(-np.pi, np.pi, numbins))
  bin_ds["theta_bins"]=xr.DataArray(np.expand_dims(a[0], axis=[1,2]), {"Bins_theta":np.mean([a[1][0:-1], a[1][1:]], axis=0),"Fly":x_array.Fly, "Batch":x_array.Batch})
  bin_ds["theta_summary"]=x_array["theta_summary"]

  p=calcpowerforfly(x_array["theta"],x_array["timestamps"])
  ps=calcpowerforfly(x_array["theta"],x_array["timestamps"],shuffle=True)
  bin_ds["theta_psd"]=xr.concat([xr.DataArray(np.expand_dims(p, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":["Non-shuffled Data"]}), 
  xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})], dim="Shuffled")

  a=x_array["turning"].plot.hist(bins=np.linspace(-2*np.pi, 2*np.pi, numbins))
  bin_ds["turning_bins"]=xr.DataArray(np.expand_dims(a[0], axis=[1,2]), {"Bins_turning":np.mean([a[1][0:-1], a[1][1:]], axis=0),"Fly":x_array.Fly, "Batch":x_array.Batch})
  bin_ds["turning_summary"]=x_array["turning_summary"]

  p=calcpowerforfly(x_array["turning"],x_array["timestamps"])
  ps=calcpowerforfly(x_array["turning"],x_array["timestamps"],shuffle=True)
  bin_ds["turning_psd"]=xr.concat([xr.DataArray(np.expand_dims(p, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":["Non-shuffled Data"]}), 
  xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})], dim="Shuffled")

  a=x_array["speed_logged"].plot.hist(bins=np.linspace(-7, 7, numbins))
  bin_ds["speed_logged_bins"]=xr.DataArray(np.expand_dims(a[0], axis=[1,2]), {"Bins_speed_logged":np.mean([a[1][0:-1], a[1][1:]], axis=0),"Fly":x_array.Fly, "Batch":x_array.Batch})
  bin_ds["speed_logged_summary"]=x_array["speed_logged_summary"]

  p=calcpowerforfly(x_array["speed"],x_array["timestamps"])
  ps=calcpowerforfly(x_array["speed"],x_array["timestamps"],shuffle=True)
  bin_ds["speed_psd"]=xr.concat([xr.DataArray(np.expand_dims(p, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":["Non-shuffled Data"]}), 
  xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})], dim="Shuffled")

  # print(bin_ds["speed_psd"])
  # a=x_array["turning"].plot.hist(bins=np.linspace(-1, 1, 201))
  # x_array["turning_bins"]=xr.DataArray(a[0], {"Bins_turning":np.mean([a[1][0:-1], a[1][1:]], axis=0)})

  a=x_array["r"].plot.hist(bins=np.linspace(0, 1.2, numbins))
  bin_ds["r_bins"]=xr.DataArray(np.expand_dims(a[0], axis=[1,2]), {"Bins_r":np.mean([a[1][0:-1], a[1][1:]], axis=0),"Fly":x_array.Fly, "Batch":x_array.Batch})

  p=calcpowerforfly(x_array["r"],x_array["timestamps"])
  ps=calcpowerforfly(x_array["r"],x_array["timestamps"],shuffle=True)
  bin_ds["r_psd"]=xr.concat([xr.DataArray(np.expand_dims(p, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":["Non-shuffled Data"]}), 
  xr.DataArray(np.expand_dims(ps, axis=[1,2,3]), coords={"Freq":f,"Fly":x_array.Fly, "Batch":x_array.Batch, "Shuffled":['Shuffled']})], dim="Shuffled")

  for varitem in list(x_array.variables):
    if "summary" in varitem:
        bin_ds[varitem]=x_array[varitem]
  # return bin_ds
  if len(save_name):
    bin_ds.to_netcdf(savepath+save_name)


  return bin_ds


# def makehistogram()

def calcpowerforfly(y, t, f=ca.calcfreq(fres=1000, log=True, fmax=5, fmin=1/60/60/24/30)/u.second, shuffle=False):
  y=y.squeeze().values
  yi=np.isnan(y)

  # f=ca.calcfreq(fres=1000, log=True, fmax=5, fmin=1/60/60/24/30)/u.second
  # print(y.shape)
  t=np.array(t.squeeze()[~yi], dtype=float)*u.second*10**-9

  # return yi
  t=t.squeeze()
  y=y[~yi]
  
  # return y,t

  if shuffle:
    rng = np.random.default_rng()
    rng.shuffle(y)

  ls=LombScargle(t=t, y=y)
  try:
    p=ls.power(f)
  except:
    # print(trial, fly)
    # when there isn't enough non nan flies
    p=np.empty(f.shape)
    p.fill(np.NaN)
    # return([t,y])
  # print(p.shape)
  # print(slicesize)
  # print("Calculated Power for Fly "+str(int(fly)))
  return p
  
def loadallsummarydata(path="", duration=2):
  # loads all summary data, as produced by get_bins
  # nestedlist_batch.append(glob.glob("CirclingData_2h_B"+str(batch)+"_F*.nc"))
  nestedlist_batch=[] 
  batch=1
  all_xarray=xr.open_mfdataset(glob.glob(path+"Circling_Summary_"+str(duration)+"h_B"+str(batch)+"_F*.nc"), combine="nested", concat_dim=["Fly"]).stack(flyid=["Batch", "Fly"])
  # sum1=test1.stack(flykey=["Fly", "Batch"])
  print(all_xarray["flyid"].shape)
  for batch in np.arange(2,7):
    # nestedlist_batch=[]
    print(batch)
    # nestedlist_batch.append(glob.glob("CirclingData_2h_B"+str(batch)+"_F*.nc"))
    if len(glob.glob(path+"Circling_Summary_"+str(duration)+"h_B"+str(batch)+"_F*.nc"))>0:
      new_xarray=xr.open_mfdataset(glob.glob(path+"Circling_Summary_"+str(duration)+"h_B"+str(batch)+"_F*.nc"), combine="nested", concat_dim=["Fly"]).stack(flyid=["Batch", "Fly"])
      all_xarray=xr.concat([all_xarray,new_xarray], dim="flyid")
    print(all_xarray["flyid"].shape)
  return all_xarray

def makehistogram(summary_xarray, data_category, xlabel="", ylabel="Counts", title=""):
  if "_bins" in data_category and len(xlabel)<1:
    # print(data_label)
    data_category=data_category[:-5]
  if len(xlabel)<1:
    xlabel=data_category  .capitalize()



  if len(title)<1:
    title="Distribution of "+data_category.capitalize()+" across all Fly Behavior"
  
  summary_df=summary_xarray[data_category+"_bins"].sum(dim="flyid").to_dataframe().reset_index()
  summary_df["recording_length"]=pd.Categorical(summary_df["recording_length"])
  # return summary_df

  summary_df.sort_values(by="recording_length", ascending=False, inplace=True)
  # return summary_df
  histogram=(
    so.Plot(data=summary_df, x="Bins_"+data_category, y=data_category+"_bins", color="recording_length").add(so.Bars())
    .label(x=xlabel, y=ylabel, title=title, color="Recording Length")
    .theme({**axes_style("ticks")})
    )
  
  return histogram

def plotsummary(summary_xarray):
  return

def calculate_bouts(xarray, target_var="speed"):
  speed_series=xarray['speed'].squeeze().to_series()
  target_series=xarray[target_var].squeeze().to_series()

  bouts=(((np.log(speed_series)>.5).rolling(10, center=True).sum())>1)
  boutnumber=((np.round(bouts).diff())>0).cumsum()
  groups=boutnumber*bouts
  # return [target_series, groups]
  chunkmeans=target_series.reset_index().groupby(np.array(groups)).mean(numeric_only=False).set_index('timestamps')
  # return chunkmeans
  chunkmeans.dropna(inplace=True)
  # valuename=target_var
  valuename="value"
  chunkmeans.columns=["Chunk "+target_var]
  chunkDB=chunkmeans.melt(ignore_index=False, value_name=valuename)
  # return chunkDB
  # # turning_all=pd.concat([chunkDB.melt(),allDB.melt(), dayDB.melt(), hourDB.melt(), minDB.melt()], axis=0)

  byframe=pd.DataFrame({"Frame "+target_var:target_series[bouts]}).melt(ignore_index=False, value_name=valuename)
  rollingaverage_hour=(target_series[bouts].resample('1H').mean()).dropna()
  hourDB=pd.DataFrame({"Hour "+target_var:rollingaverage_hour}).melt(ignore_index=False, value_name=valuename)
  rollingaverage_hour=(target_series[bouts].resample('1D').mean()).dropna()
  dayDB=pd.DataFrame({"Day "+target_var:rollingaverage_hour}).melt(ignore_index=False, value_name=valuename)
# dayDB
  rollingaverage_min=(target_series[bouts].resample('60s').mean()).dropna()
  minDB=pd.DataFrame({"Minute "+target_var:rollingaverage_min}).melt(ignore_index=False, value_name=valuename)
  dataframe=pd.concat([chunkDB, byframe, dayDB, hourDB, minDB], axis=0)
  
  dataframe["Fly"]=int(xarray["Fly"])
  dataframe["Batch"]=int(xarray["Batch"])

  dataframe["Fly"]=pd.Categorical(dataframe["Fly"])
  dataframe["Batch"]=pd.Categorical(dataframe["Batch"])

  dataframe["Recording Length"]=int(xarray["recording_length"])
  dataframe["Recording Length"]=pd.Categorical(dataframe["Recording Length"])

  return dataframe

def bootstrapintervals(xarray, target_var="speed"):
  speed_series=xarray['speed'].squeeze().to_series()
  

  if target_var=="log_speed":
    target_series=np.log(speed_series)
  else:
    target_series=xarray[target_var].squeeze().to_series()
    

  flynum=xarray["Fly"]
  batchnum=xarray["Batch"]
  recordinglength=int(xarray["recording_length"])

  bouts=(((np.log(speed_series)>.5).rolling(10, center=True).sum())>1)
  boutnumber=((np.round(bouts).diff())>0).cumsum()
  groups=boutnumber*bouts
  # return [target_series, groups]
  chunkmeans=target_series.reset_index().groupby(np.array(groups)).mean(numeric_only=False).set_index('timestamps')
  # return chunkmeans
  chunkmeans.dropna(inplace=True)
  # valuename=target_var
  valuename="value"
  chunkmeans.columns=["Chunk "+target_var]
  # return chunkmeans
  # return target_series[bouts]
  chunkpd=bootstrap_helper(chunkmeans, "Activity Bout", target_var, fly=flynum, batch=batchnum, recording_length=recordinglength)
  # return chunkpd
  # ci=sci.stats.bootstrap((chunkmeans,), np.mean, method='basic').confidence_interval
  # return pd.DataFrame({"Variable":target_var, "Method":"Chunked", "Fly":"Fly1", "Batch":"BatchX", "CI_low":ci[0], "Mean":np.mean(chunkmeans), "CI_high":ci[1]})

  # chunkDB=chunkmeans.melt(ignore_index=False, value_name=valuename)
  # return chunkDB
  # turning_all=pd.concat([chunkDB.melt(),allDB.melt(), dayDB.melt(), hourDB.melt(), minDB.melt()], axis=0)
  # return target_series[bouts].resample('1H').mean().dropna()
  hourpd=bootstrap_helper((target_series[bouts].resample('1H').mean()).dropna(), "Hour", target_var, fly=flynum, batch=batchnum, recording_length=recordinglength)
  framepd=bootstrap_helper(target_series[bouts].dropna()[1:], "Frame", target_var, fly=flynum, batch=batchnum, recording_length=recordinglength)
  minutepd=bootstrap_helper((target_series[bouts].resample('60s').mean()).dropna(), "Minute", target_var, fly=flynum, batch=batchnum, recording_length=recordinglength)
  daypd=bootstrap_helper((target_series[bouts].resample('1D').mean()).dropna(), "Day", target_var, fly=flynum, batch=batchnum, recording_length=recordinglength)

  dfall=pd.concat([framepd, chunkpd, minutepd, hourpd, daypd], ignore_index=True)
  dfall["FlyID"]=pd.Categorical("B"+dfall["Batch"].astype(str)+"-F"+dfall["Fly"].astype(str))
  return dfall



  byframe=pd.DataFrame({"Frame "+target_var:target_series[bouts]}).melt(ignore_index=False, value_name=valuename)
  rollingaverage_hour=(target_series[bouts].resample('1H').mean()).dropna()
  hourDB=pd.DataFrame({"Hour "+target_var:rollingaverage_hour}).melt(ignore_index=False, value_name=valuename)
  rollingaverage_hour=(target_series[bouts].resample('1D').mean()).dropna()
  dayDB=pd.DataFrame({"Day "+target_var:rollingaverage_hour}).melt(ignore_index=False, value_name=valuename)
# dayDB
  rollingaverage_min=(target_series[bouts].resample('60s').mean()).dropna()
  minDB=pd.DataFrame({"Minute "+target_var:rollingaverage_min}).melt(ignore_index=False, value_name=valuename)
  dataframe=pd.concat([chunkDB, byframe, dayDB, hourDB, minDB], axis=0)
  
  dataframe["Fly"]=int(xarray["Fly"])
  dataframe["Batch"]=int(xarray["Batch"])

  dataframe["Fly"]=pd.Categorical(dataframe["Fly"])
  dataframe["Batch"]=pd.Categorical(dataframe["Batch"])

  dataframe["Recording Length"]=int(xarray["recording_length"])
  dataframe["Recording Length"]=pd.Categorical(dataframe["Recording Length"])

  return dataframe
  

def bootstrap_helper(series, method_name, variable_name, fly=0, batch=0, recording_length=9):
  ci=sci.stats.bootstrap((series.values,), np.mean, method='basic', batch=100).confidence_interval

  return pd.DataFrame({"Variable":variable_name, "Block Method":method_name, "Fly":fly, "Batch":batch, "Recording Length":recording_length, "CI_low":ci[0], "Mean":np.mean(series.values), "CI_high":ci[1]}, index=[0])


def compareflies(filelist, range='all', variables=['direction']):
  if range=='all':
    range=np.arange(len(filelist))
  

  # needs_initialization=True

  # for fly in range:
  b=Parallel(n_jobs=1, verbose=10)(delayed(compareflies_parallel)(file, variables) for file in filelist)
  return pd.concat(b, ignore_index=True )

  # return all

def compareflies_parallel(filepath, variables=['direction']):
  # print(variables)
  xarray=xr.open_dataset(filepath)
  a=Parallel(n_jobs=1, verbose=0)(delayed(bootstrapintervals)(xarray, v) for v in variables)
  return pd.concat(a, ignore_index=True)