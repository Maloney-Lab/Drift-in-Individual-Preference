# make the matrix
# call simpledrift to fill in the values for different drift and bh
# pass on to heatmap
# pay attention to axes so they give the correct values
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.npyio import save
import simpledrift as sd
import math
from joblib import Parallel, delayed
import os
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import xarray as xr
from datetime import date
from sys import getsizeof
import sys
import gc
import seaborn as sns

# import netCDF4

def matrixmaker(envi, bhstrats=np.linspace(0, 1, 4), driftstrats=np.linspace(0, 1, 4), showgraphs=False, figuresavepath='../Results/figs',
     runindex=0, environtype=-1, freqmin=-1, freqmax=-1, power=-1, envimeanvariance=[.3], envivariance=[.1], birthrate=[40], phi=[1], matureage=[10],
     numberofbins=100,
     driftmaxdistribution=[.3],
     nameprefix="", saveallprefs=False, savealldays=True, savedata=True, saveenv=True):

    prefmean=0
    adaptivetracking=0


    iterables=[matureage, np.round(birthrate,3), np.round(envivariance,3), np.round(envimeanvariance,3), np.round(bhstrats,3), np.round(driftstrats,3), np.round(phi,3), np.round(driftmaxdistribution,3)]
#     interables=np.round(iterables,3)
#     print(iterables)
    index=pd.MultiIndex.from_product(iterables, names=["matureage", "birthrate", "envvar", "envmeanvar", "bet-hedging", "drift", "phi", "boundinggaussian"])
    # print(index)
    print("Number of combinations: ", len(index))

    allseries=Parallel(n_jobs=-1, verbose=5)(delayed(sd.driftmodeling)(envi, prefmean,
        prefvariance=i[4],
        phi=i[6],

        driftvariance=i[5], 
        adaptivetracking=adaptivetracking, 
        birthrate=i[1],
        matureage=i[0],
        driftmaxdistribution=i[7],
        # driftmaxdistribution=0,

        showgraphs=showgraphs,
        figuresavepath=figuresavepath,
        envimeanvariance=i[3],
        envivariance=i[2],
        numberofbins=numberofbins,
        saveallprefs=saveallprefs,
        # envivariance=.001,
        savealldays=savealldays) for i in index)
            # print(["Drift is: ", driftvariance[y], "Bet-hedging is: ", driftvariance[x]] )
            # print(matrix)
    
    allseriesshape=np.shape(allseries)

    if saveallprefs and savealldays:
        dataframe=xr.DataArray(data=allseries,
        coords={"Parameters":index, "Preferences":np.arange(allseriesshape[1]), "day":np.arange(allseriesshape[2])},
        dims=["Parameters", "Preferences", "day"])

        # np.save('index_findme', index)
        # dataframe.reset_index("Parameters").to_netcdf('find_me_please.nc')
    elif saveallprefs:
        allseriesshape=np.shape(allseries)
        dataframe=xr.DataArray(data=allseries,
        coords={"Parameters":index, "Preferences":np.arange(allseriesshape[1])},
        dims=["Parameters", "Preferences"])
        # np.save('index_findme', index)
        # dataframe.reset_index("Parameters").to_netcdf('find_me_please.nc')
        dataframe=dataframe.reset_index("Parameters")
    elif savealldays:
        # dataframe=pd.DataFrame(allseries, index=index, dtype='float')
        dataframe=xr.DataArray(data=allseries,
        coords={"Parameters":index, "day":np.arange(allseriesshape[1])},
        dims=["Parameters", "day"])

    ## This version is for day names, rather than indexes
    # if savealldays:
    #     columnnames=np.empty(dataframe.shape[1], dtype='object')
    #     for i in range(dataframe.shape[1]):
    #         columnnames[i]='Day '+str(i)
    #     dataframe.columns=columnnames
    # else:
    #     columnnames=np.empty(dataframe.shape[1], dtype='object')
    #     for i in range(dataframe.shape[1]):
    #         columnnames[i]='Day '+str(envi.shape[0])
    #     dataframe.columns=columnnames 
    else:
        # dataframe=pd.DataFrame(allseries, index=in    dex)
        # print(np.array(allseries).squeeze())
        # return allseries
        dataframe=xr.DataArray(data=np.array(allseries).squeeze(),
            coords={"Parameters":index},
            dims=["Parameters"])

    # if savealldays:
    #     columnnames=np.empty(dataframe.shape[1], dtype='object')
    #     for i in range(dataframe.shape[1]):
    #         columnnames[i]=i
    #     dataframe.columns=columnnames
    # else:
    #     columnnames=np.empty(dataframe.shape[1], dtype='object')
    #     for i in range(dataframe.shape[1]):
    #         columnnames[i]=envi.shape[0]-1
    #         print(i)
    #     dataframe.columns=columnnames 

    if os.path.exists(figuresavepath):

        if nameprefix=="":
          filename='R'+str(runindex)
        else:
          filename=nameprefix+'_R'+str(runindex)
        
        if environtype!=-1:
          filename+='_T'+str(environtype)+'Populations.parquet'

        ncname=filename+"_Populations"
        # dataframe=xr.DataArray(dataframe)
        # dataframe=dataframe.unstack().rename({"dim_1":"day"})

        if savedata:

            dataframe.unstack().to_netcdf(ncname+".nc")
        dataframe=dataframe.unstack()
        
        enviname=filename+'_env'
        if saveenv:
            np.save(os.path.join(figuresavepath,enviname), envi)

    else:
        print('not saving, no valid path')
    # print(getsizeof(dataframe))
    return dataframe

def frequency_phaseplane(i=1, fbands=1/(np.arange(1,21)[:0:-1]*2), strategy_resolution=51, numberofdays=101,
 matureage=[10], envimeanvariance=[.1], birthrate=[20], envivariance=[.1], bhmax=.5, driftmax=.1, numberofbins=100,
 savealldays=False, saveallprefs=False,
 savepath="../Results/", filename_prefix="", savechunks=True):

    control_res=3

    # savepath="../Results/figs_1-7-21/"

    # fbands=np.arange(frequency_resolution)

    # fbands=fbands[:0:-1]
    # # fbands=np.array([4.0])
    # fbands=2*fbands
    # bands=1/fbands
    fbands.round(3)

    # Each run should run all fbands, each iteration is one run 

    # output3=mm.matrixmaker(sa, bhstrats=np.linspace(0, .5, 51), driftstrats=np.linspace(0, .1, 51),
    #                       envimeanvariance=np.linspace(.06, .3, 5), envivariance=np.linspace(.06, .3, 5), birthrate=np.linspace(20,60,3, dtype=int), matureage=np.linspace(5,20,4, dtype=int), nameprefix="test3")

    first_iteration=True
    #Step 1: Make the core signals
    for fi in np.arange(1,len(fbands)-1):
        f=fbands[fi]
        print(str("F"+str(f)+"_R"+str(i)))
        envi, us=sd.makefilterednoise(lowerbound=fbands[fi-1], upperbound=fbands[fi+1], numberofdays=numberofdays)

        dataframe=matrixmaker(us, bhstrats=np.linspace(0, bhmax, strategy_resolution),
        driftstrats=np.linspace(0, driftmax, strategy_resolution),
        numberofbins=numberofbins,
        figuresavepath=savepath,
        # envivariance=np.linspace(.1, .3, 3),
        envimeanvariance=envimeanvariance,
        envivariance=envivariance,
        # birthrate=np.linspace(20,60,3, dtype=int),
        birthrate=birthrate,
        savealldays=savealldays,
        saveallprefs=saveallprefs,
        matureage=matureage,
        # matureage=[10],
        nameprefix=str("F"+str(np.round(f, 3))),
        runindex=i,
        savedata=False
        )
        # print('dataframe')
        # print(dataframe)


        #make environment into array that can be  merged into xarray total
        
        envda=xr.DataArray(np.expand_dims(us, 1), coords={"day": np.arange(numberofdays), "freq":[f]}, dims=["day", "freq"], name="envi")
        envda.swap_dims()
        # print(envda)
        if first_iteration:
            # if not savealldays:
            #     print('Is here the problem?')
            #     print(savealldays)
            #     x_all=dataframe.isel(day=0)
            # else:
            x_all=dataframe


            # x_all=dataframe.to_xarray()
            x_all=x_all.expand_dims({"freq":[f]})
            first_iteration=False

            x_all=xr.Dataset({"environment_mean":envda, "sim_results":x_all})
            if savechunks:
                x_all.to_netcdf(savepath+filename_prefix+str(i)+"_F"+str(f)+"_R"+str(i)+".nc")
        else:
            # if not savealldays:
            #     x=dataframe.isel(day=0)
            # else:
            #     x=dataframe
            x=dataframe
        # x_all=x_all+x #But actual xarray concatenation
            # x.merge(x_all, join=outer)
            x=x.expand_dims({"freq":[f]})
            first_iteration=False
            # print(envda)
            x=xr.Dataset({"environment_mean":envda, "sim_results":x})
            if savechunks:
                x.to_netcdf(savepath+filename_prefix+str(i)+"_F"+str(f)+"_R"+str(i)+".nc")
            else:
                x_all=xr.merge([x_all, x])
                        # x_all=xarray.merge([x_all, x], dim="freq")

            # print('merged')conda
            # x_all=xarray.merge([x_all, x])

    # Step 2, collate each run into one xarray
    # print(x_all)
    # savepath="../Results/figs/"
    today = date.today()
    if not savechunks:
        filename=savepath+filename_prefix+str(today)+"_R"+str(i)+".nc"
   
        filename=savepath+filename_prefix+str(today)+"_R"+str(i)+".nc"
        print(f"Run saved at {filename}")
        x_all.to_netcdf(filename)
        print(actualsize(x_all))
    return

def getallruns(listofinputs, day=1000):
    for i, input in enumerate(listofinputs):
        print(i,input)
        sin=xr.load_dataset(input)
        sin.coords["Run"]=i
        sin=sin.expand_dims("Run")
        if type(day) == int:
            sin=sin.isel(day=day)
        if i==0:
            allruns=sin
        else:
            allruns=allruns.merge(sin)
    return allruns

def calcgeomean(allruns):
    geomean=np.log10(allruns['sim_results'].squeeze()).mean(dim="Run")
    return geomean  


def makeMatureAgexFreqFigure2(phaseplaneoutput):

    for j, freq in enumerate(np.array(phaseplaneoutput['freq'])):
        envi=sd.meantofullenvi(envimeans=phaseplaneoutput["environment_mean"][:, j], envimeanvariance=float(phaseplaneoutput["envmeanvar"]), envivariance=float(phaseplaneoutput["envvar"]), numberofbins=100) 
        fig=px.imshow(envi, aspect='auto')
        fig.show()
    for i in range(0,len(phaseplaneoutput['matureage'])):

        fig=px.imshow(np.log(phaseplaneoutput["sim_results"].isel({"day":100, "matureage":i}).squeeze()), facet_col="freq")

        fig.layout.title.text="Mature Age: "+str(int(phaseplaneoutput["matureage"][i]))
        fig.layout.coloraxis.colorbar.title.text="Final Population (Log)"
        for j, freq in enumerate(np.array(phaseplaneoutput['freq'])):
            # print(freq)
            fig.layout.annotations[j]['text']= 'Period = %d' %int(1/freq)
        fig.show()

def makeMatureAgexFreqFigure(phaseplaneoutput):

    gap=.005

    numrows=len(phaseplaneoutput['freq'])+1
    print(numrows)
    # numrows=1
    numrows=20
    numcols=len(phaseplaneoutput['matureage'])
    print(numcols)
    numcols=25
    # numcols=1
    print(numcols, numrows)
    
    fig = make_subplots(rows=numrows, cols=numcols, shared_yaxes=False, vertical_spacing=gap, horizontal_spacing=gap)
    # fig.update_layout(height=1200, width=1200)
    # fig.show()
    fig.update_layout(height=2400, width=2400)

    for j, freq in enumerate(np.array(phaseplaneoutput['freq'])):
        if j< numrows:
            envi=sd.meantofullenvi(envimeans=phaseplaneoutput["environment_mean"][:, j], envimeanvariance=float(phaseplaneoutput["envmeanvar"]), envivariance=float(phaseplaneoutput["envvar"]), numberofbins=100) 
            # envi=np.zeros([4,4])
            print(j)
            fig.add_trace(go.Heatmap(z=envi), row=j+1, col=1)
            # fig.show()
            # print(j)
            for i in range(0,len(phaseplaneoutput['matureage'])):
                if i<numcols-1:

                    if "day" in phaseplaneoutput["sim_results"].dims:
                        fig.add_trace(go.Heatmap(z=np.log(phaseplaneoutput["sim_results"].isel({"day":100, "matureage":i, "freq":j}).squeeze())), col=i+2, row=j+1)
                    else:
                        fig.add_trace(go.Heatmap(z=np.log(phaseplaneoutput["sim_results"].isel({"matureage":i, "freq":j}).squeeze())), col=i+1, row=j+1)
                    fig.layout.coloraxis.colorbar.showticklabels=False
                    # fig.layout.title.text="Mature Age: "+str(int(phaseplaneoutput["matureage"][i]))
                    # fig.layout.coloraxis.colorbar.title.text="Final Population (Log)"
                    # for j, freq in enumerate(np.array(phaseplaneoutput['freq'])):
                    #     # print(fre
                    #     fig.layout.annotations[j]['text']= 'Period = %d' %int(1/freq)
    print("Done with heatmaps")
    # fig.update_layout(height=1200, width=1200)

    for anno in fig['layout']['annotations']:
        anno['text']=""
    for axis in fig.layout:
        
        if type(fig.layout[axis]) == go.layout.YAxis:
            fig.layout[axis].title.text = ''
            fig.layout[axis].tickfont = dict(color = 'rgba(0,0,0,0)')
            
        if type(fig.layout[axis]) == go.layout.XAxis:
            fig.layout[axis].title.text = ''
            fig.layout[axis].tickfont = dict(color = 'rgba(0,0,0,0)')
        if type(fig.layout[axis]) == go.layout.coloraxis:
            fig.layout[axis].title.text = ''
            fig.layout[axis].tickfont = dict(color = 'rgba(0,0,0,0)')
    # fig.show()
    # fig.write_image('../Results/figs/MatureAgexFreq-test.pdf')
    return fig

def dailychangehistogram(phaseplaneoutput):
    return phaseplaneoutput


def realworlddataPhaseSpace(realworlddata, i=1, strategy_resolution=51, numberofdays=101,
 matureage=[10], envimeanvariance=[.1], birthrate=[20], envivariance=[.1], bhmax=.5, driftmax=.1, numberofbins=100,
 savealldays=False, saveallprefs=False,
 savepath="../Results/", filename_prefix=""):
    f=0
    first_iteration=True
    #Step 1: Make the core signals
    try:
        realworlddata=realworlddata.expand_dims("siteID")
    except:
        print("Processing multiple sites")
    try:
        realworlddata=realworlddata.expand_dims("variable")
    except:
        print("Processing Multiple Variables")
    try:
        realworlddata=realworlddata.expand_dims("chunk")
    except:
        print("Processing multiple Chunks")


    for ii, i in enumerate(realworlddata['siteID']):
        # print(i)
        for ji, j in enumerate(realworlddata['variable']):
            for ki, k in enumerate(realworlddata['chunk']):
                # print(k)  
                # 
                if ~np.isnan(realworlddata[ii,ji,ki]).all():
                    l=realworlddata[ii,ji,ki].copy()

                    siteID=i["siteID"].values[()]
                    variable=j["variable"].values[()]
                    chunk=k["chunk"].values[()]
                    # print(siteID, variable, chunk)
                    # print(l)
                    # print(i)
                    # print(j)
                    # print(k)
                # f=fbands[fi]
        # print(str("F"+str(f)+"_R"+str(i)))
        # envi, us=sd.makefilterednoise(lowerbound=fbands[fi-1], upperbound=fbands[fi+1], numberofdays=numberofdays)
                    # return l.values
                    # l[0]=0
                    #Run Simulations
                    dataframe=matrixmaker(l.values, bhstrats=np.linspace(0, bhmax, strategy_resolution),
                    driftstrats=np.linspace(0, driftmax, strategy_resolution),
                    numberofbins=numberofbins,
                    figuresavepath=savepath,
                    # envivariance=np.linspace(.1, .3, 3),
                    envimeanvariance=envimeanvariance,
                    envivariance=envivariance,
                    # birthrate=np.linspace(20,60,3, dtype=int),
                    birthrate=birthrate,
                    savealldays=savealldays,
                    matureage=matureage,
                    saveenv=False,
                    # matureage=[10],
                    # nameprefix=str("F"+str(np.round(f, 3))),
                    runindex=0,
                    savedata=False,
                    saveallprefs=saveallprefs,
                    )
                    print(dataframe)
                    dataframe=dataframe.squeeze()

                    #make environment into array that can be merged into xarray total
                    # print(l)
                    # envda=xr.DataArray(np.expand_dims(l, 1), coords={"day": np.arange(numberofdays), "siteID":siteID, "chunk":chunk, "variable":variable}, dims=["day", "siteID", "chunk", "variable"], name="envi")
                    envda=l
                    envda.swap_dims()
                    if first_iteration:
                        # x_all=xr.DataArray(dataframe, dims=("MultiIndex", "day")).unstack("MultiIndex")
                        x_all=dataframe

                        # x_all=dataframe.to_xarray()
                        # print(x_all)
                        # x_all=x_all.expand_dims(coords={"day": np.arange(numberofdays),"siteID":siteID, "chunk":chunk, "variable":variable}, dims=["day", "siteID", "chunk", "variable"])
                        first_iteration=False

                        x_all=xr.Dataset({"environment_mean":envda, "sim_results":x_all})
                        x_all=x_all.expand_dims(dim=["siteID", "variable", "chunk"])
                    else:
                        # x=xr.DataArray(dataframe, dims=("MultiIndex", "day")).unstack("MultiIndex")
                        x=dataframe
                    # x_all=x_all+x #But actual xarray concatenation
                        # x.merge(x_all, join=outer)
                        # x=x.expand_dims({"siteID":siteID, "chunk":chunk, "variable":variable})
                        first_iteration=False
                        x=xr.Dataset({"environment_mean":envda, "sim_results":x})
                        x=x.expand_dims(dim=["siteID", "variable", "chunk"])
                        # print(x_all)
                        # print(x)
                        try:
                            x_all=xr.combine_by_coords([x_all, x])
                            f+=1
                            # print(f)
                        except:
                            try:
                                x_all=xr.concat([x_all,x], dim="siteID")
                            except:
                                x_all=xr.concat([x_all,x], dim="variable")
                            # print("bug trying to combine runs")
                            # return x_all, x

                        # x_all=xr.merge([x_all, x])
                                    # x_all=xarray.merge([x_all, x], dim="freq")
        #                 break
        #         else:
        #             continue
        #         break
        #     else:
        #         continue
        #     break
        # else:
        #     continue
        # break
               
                        


                        # print('merged')conda
                        # x_all=xarray.merge([x_all, x])

    # Step 2, collate each run into one xarray
    # print(x_all)
    # savepath="../Results/figs/"
    today = date.today()
    filename=savepath+filename_prefix+str(today)+"_R"+str(i)+".nc"
    # x_all.to_netcdf(filename)
    return x_all

def actualsize(input_obj):
    memory_size = 0
    ids = set()
    objects = [input_obj]
    while objects:
        new = []
        for obj in objects:
            if id(obj) not in ids:
                ids.add(id(obj))
                memory_size += sys.getsizeof(obj)
                new.append(obj)
        objects = gc.get_referents(*new)
    return memory_size

def mapwideformforseabornheatmap(data, color=0):
  wideform= data.pivot(columns='drift', index='bet-hedging', values='sim_results')
  # print(wideform)
  # return wideform
  return sns.heatmap(wideform.iloc[-1::-1], cmap='RdBu')

# def plotidealstrategy(data, dim="drift"):
#     argmaxes=data["sim_results"].squeeze().argmax(dim=["drift", "bet-hedging"])
#     a=

def plotstrategyplots(data, sumstat="none", param1="period", param2="matureage", sharey=True):
    argmaxes=data.squeeze().argmax(dim=["drift", "bet-hedging"])
    f, axs = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw=dict(width_ratios=[4, 4]), sharey=sharey)

    # argmaxes=data.squeeze().argmax(dim=["drift", "bet-hedging"])
    if sumstat=="none":
        argmax_bh_df=argmaxes["bet-hedging"].to_dataframe()["sim_results"].reset_index()
    elif sumstat=="mean":
        argmax_bh_df=argmaxes["bet-hedging"].mean(dim="Run").to_dataframe()["sim_results"].reset_index()
    elif sumstat=="std":
        argmax_bh_df=argmaxes["bet-hedging"].std(dim="Run").to_dataframe()["sim_results"].reset_index()
    argmax_bh_df["period"]=1/argmax_bh_df["freq"]
    argmax_bh_df["period"]=np.round(argmax_bh_df["period"], 1)

    argmax_bh_df["sim_results"]=argmax_bh_df["sim_results"]*data["bet-hedging"][-1].values/len(data["bet-hedging"])
    argmax_bh_df.drop(columns=["freq"], inplace=True)
    # return argmax_bh_df
    sns.heatmap(argmax_bh_df.pivot(index=param1, columns=param2, values="sim_results"), ax=axs[0], vmax=.1, vmin=0, cmap="rocket")

    if sumstat=="none":
        argmax_d_df=argmaxes["drift"].to_dataframe()["sim_results"].reset_index()
    elif sumstat=="mean":
        argmax_d_df=argmaxes["drift"].mean(dim="Run").to_dataframe()["sim_results"].reset_index()
    elif sumstat=="std":
        argmax_d_df=argmaxes["drift"].std(dim="Run").to_dataframe()["sim_results"].reset_index()
    # argmax_d_df=argmaxes["drift"].to_dataframe()["sim_results"].reset_index()
    argmax_d_df["period"]=1/argmax_d_df["freq"]
    argmax_d_df["period"]=np.round(argmax_d_df["period"], 1)

    argmax_d_df["sim_results"]=argmax_d_df["sim_results"]*data["drift"][-1].values/len(data["drift"])
    argmax_d_df.drop(columns=["freq"], inplace=True)
    sns.heatmap(argmax_d_df.pivot(index=param1, columns=param2, values="sim_results"), ax=axs[1], vmax=.05, vmin=0, cmap="flare_r")


    # sns.heatmap(argmaxes["bet-hedging"]*data["bet-hedging"][-1], ax=axs[0])
    # sns.heatmap(argmaxes["drift"]*data["drift"][-1], ax=axs[1])

    # axs[0].xlabel("Mature Age", )
    # plt.ylabel("Period", ax=axs[0])
    # plt.title("Ideal amount of bet-hedging for each period and mature age", ax=axs[0])
    # label x axis of first subplot 
    axs[0].set_xlabel("Age of Reproductive Maturity (days)")

    if param1=="period":
        axs[0].set_ylabel("Environmental Fluctuation Period (days)")
        axs[1].set_ylabel("Environmental Fluctuation Period (days)")
    elif param1=="matureage":
        axs[0].set_ylabel("Age of Reproductive Maturity (days)")
        axs[1].set_ylabel("Age of Reproductive Maturity (days)")
    
    if param2=="matureage":
        axs[1].set_xlabel("Age of Reproductive Maturity (days)")
        axs[0].set_xlabel("Age of Reproductive Maturity (days)")
    elif param2=="envmeanvar":
        axs[1].set_xlabel("Environmental Fluctuation Strength (AU)")
        axs[0].set_xlabel("Environmental Fluctuation Strength (AU)")


    if sumstat=='std':
        axs[0].set_title("Standard Deviation of Ideal amount of bet-hedging")
        axs[1].set_title("Standard Deviation of Ideal amount of drift")
    else:
        axs[0].set_title("Ideal amount of bet-hedging")
        axs[1].set_title("Ideal amount of drift")

    # axs[0].yaxis.set_major_formatter('{x:0<2.1f}')
    # axs[1].yaxis.set_major_formatter('{x:0<2.1f}')
    axs[0].invert_yaxis()
    if not sharey:
        axs[1].invert_yaxis()
    # axs[1].invert_yaxis()

    return f, axs

# def plotsamplestrategyplots(data, param1="period", param2="matureage"):
#     strats_subset=data.isel({param1:})