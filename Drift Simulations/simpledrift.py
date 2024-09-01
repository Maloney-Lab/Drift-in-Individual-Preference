import numpy as np
import scipy.stats as sci
import matplotlib.pyplot as plt
import time
import math
import os
import scipy.stats as stat
import pandas as pd
import plotly.express as px
import colorednoise as cn
import xarray as xr


def sinwaveinput(numberofbins, numberofdays, envimean, envivariance, maxsurvivalrate, gain, per):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
    envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate

    for t in range(1,numberofdays):
        envi[:,t]=sci.norm.pdf(x,(envimean+gain*np.sin(t*np.pi*2/per)),envivariance) # Making envi a sin wave that changes over time
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def diffusion(numberofbins, numberofdays, envimean,envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
    envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
    blur=np.zeros((numberofbins,numberofbins)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

    for b in range(numberofbins):
        blur[b,:]=sci.norm.pdf(x,x[b],dailydrift)

    for t in range(1,numberofdays):
        for b in range (numberofbins):
            # blur=sci.norm.pdf(x,x[b],dailydrift)
            envi[:,t]+=envi[b,t-1]*blur[b,:]/np.sum(blur[b,:])
        # envi[:,t]=sci.norm.pdf(x,(envimean,envivariance) # Making envi a sin wave that changes over time
        # envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def randomwalk(numberofbins, numberofdays, envimean,envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envi[:,0]=sci.norm.pdf(x,envimean,envivariance) # A gaussian of environment with center around 0
    envi=envi/(np.max(envi))*maxsurvivalrate # Normalizing the maximum envi value and factoring in deathrate
    blur=np.zeros((numberofbins,numberofbins)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]

    for t in range(1,numberofdays):
        envimean+=np.random.normal(0,dailydrift)
        envi[:,t]=sci.norm.pdf(x,envimean,envivariance) # Making envi a sin wave that changes over time
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def metropolishastingsdrift(numberofbins, numberofdays, envimeanvariance, envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    envimean=np.random.normal(0,envimeanvariance)
    envi[:,0]=sci.norm.pdf(x, envimean, envivariance) # A gaussian of environment with center around 0

    for t in range(1, numberofdays):
        pcurrentvalue=stat.norm.pdf(envimean,0,envimeanvariance)
        # proposedvalue=np.random.normal(0,variability)
        proposedvalue=np.random.normal(0,envimeanvariance)*dailydrift+envimean*(1-dailydrift)

        pproposedvalue=stat.norm.pdf(proposedvalue,0,envimeanvariance)
        if pproposedvalue/pcurrentvalue>(1-np.random.rand()):
            # return proposedvalue*percentdrift+currentvalue*(1-percentdrift)
            envimean=proposedvalue
        envi[:,t]=sci.norm.pdf(x, envimean, envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def whitedrift(numberofbins, numberofdays, envimeanvariance, envivariance, maxsurvivalrate, dailydrift):
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    # envimean=np.random.normal(0,envimeanvariance)

    for t in range(1, numberofdays):
        envimean=np.random.normal(0,envimeanvariance)
        envi[:,t]=sci.norm.pdf(x, envimean, envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate

    return(envi)

def powerdrift(numberofbins=100, numberofdays=10, envimeanvariance=1, envivariance=.3, maxsurvivalrate=1, power=1):
    # beta = 1 # the exponent
    samples = 100 # number of samples to generate
    y = cn.powerlaw_psd_gaussian(power, numberofdays)*envimeanvariance

    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)
    # envimean=np.random.normal(0,envimeanvariance)
    # stat.norm.ppf()
    for t in range(numberofdays):
        envi[:,t]=sci.norm.pdf(x, y[t], envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    return(envi)

def makefilterednoise(numberofbins=100, numberofdays=50, envimeanvariance=.1, envivariance=.3, maxsurvivalrate=1, lowerbound=-1, upperbound=-1,  filtertype='bandpass', oversamplerate=4, lengthbuffer=2, power=0, normalizevariance='true'):
    totallength=numberofdays*oversamplerate*lengthbuffer
    s=cn.powerlaw_psd_gaussian(power, totallength)*envimeanvariance

    frequencies=np.fft.fftfreq(len(s))


    if upperbound>.5:
        upperbound=.5
        print('anchoring to nyquist')

    fs=np.zeros(len(s), dtype=complex)
    if (upperbound<=0) and (lowerbound<=0):
        filtertype='none'
    elif upperbound<=0:
        filtertype='highpass'
    elif lowerbound<=0:
        filtertype='lowpass'
    # upperbound=.2
    # lowerbound=.1
    lowerbound/=oversamplerate
    upperbound/=oversamplerate
    
    if lowerbound<frequencies[1]:
        lowerbound=frequencies[1]

    freqrange=upperbound-lowerbound
    # print(frequencies)
    
    if filtertype !='none':
        # reallowerneg=np.min(np.abs(-frequencies+lowerbound))
        # lowerindexneg=int(np.where(np.abs(frequencies+lowerbound)==reallowerneg)[0][0])
        # reallowerpos=np.min(np.abs(frequencies-lowerbound))
        # lowerindexpos=int(np.where(frequencies==-frequencies[lowerindexneg])[0][0])
        # realupper=np.min(np.abs(frequencies-upperbound))
        # upperindexpos=int(np.where(np.abs(frequencies-upperbound)==realupper)[0][0])
        # upperindexneg=int(np.where(frequencies==-frequencies[upperindexpos])[0][0])
        # if upperindexneg-upperindexpos<=2:
        #     upperindexpos+=2

        # if lowerbound<=0:
        #     lowerindexneg=0
        # else:
        #     lowerindexneg = (np.abs(frequencies - lowerbound)).argmin()
        # reallower=frequencies[lowerindexneg]
        # upperindexpos= (np.abs(frequencies - upperbound)).argmin()
        # realupper=frequencies[upperindexpos]
        # print(realupper)
        # # print('LI+', lowerindexpos, frequencies[lowerindexpos])
        # # print('LI-', lowerindexneg, frequencies[lowerindexneg])
        # # print('UI+',upperindexpos, frequencies[upperindexpos])
        # # print('UI-', upperindexneg, frequencies[upperindexneg])
        fs=np.fft.fft(s)

        lowerindexpos=-1
        upperindexpos=-1
        lowerindexneg=-1
        upperindexneg=-1

        if upperbound>max(np.abs(frequencies)):
            upperbound=max(np.abs(frequencies)) #Set to nyquist frequency if above

        for i in range(len(frequencies)):
            if frequencies[i]>=lowerbound and lowerindexpos==-1:
                reallower=frequencies[i]
                lowerindexpos=i
            if frequencies[i]>=upperbound and upperindexpos==-1:
                realupper=frequencies[i]
                upperindexpos=i
            if -frequencies[-i]>=lowerbound and lowerindexneg==-1:
                # reallower=frequencies[-i]
                lowerindexneg=len(frequencies)-i
            if -frequencies[-i]>=upperbound and upperindexneg==-1:
                # realupper=frequencies[-i]
                upperindexneg=len(frequencies)-i
        if upperindexpos==-1:
                upperindexpos=upperindexneg-1

        # print('LI+', lowerindexpos, frequencies[lowerindexpos])
        # print('LI-', lowerindexneg, frequencies[lowerindexneg])
        # print('UI+',upperindexpos, frequencies[upperindexpos])
        # print('UI-', upperindexneg, frequencies[upperindexneg])
        if filtertype=='notch':
            fs[lowerindexpos:upperindexpos].real=0
            fs[upperindexneg:lowerindexneg].real=0
        # fs[upperindexneg:lowerindexneg].real=np.random.uniform(0, 2*np.pi, (lowerindexneg-upperindexneg,))
            fs[lowerindexpos:upperindexpos].imag=0
            fs[upperindexneg:lowerindexneg].imag=0
        elif filtertype=='bandpass':
            fs[0:lowerindexpos]=0
            fs[upperindexpos:upperindexneg]=0
            fs[lowerindexneg:]=0
        elif filtertype=='lowpass':
            fs[upperindexpos:upperindexneg]=0
        elif filtertype=='highpass':
            fs[0:lowerindexpos]=0
            fs[lowerindexneg:]=0
        # fs[lowerindexpos:upperindexpos].imag=np.random.uniform(0, 2*np.pi, (upperindexpos-lowerindexpos,))
        # fs[upperindexneg:lowerindexneg].imag=np.random.uniform(0, 2*np.pi, (lowerindexneg-upperindexneg,))
        s=np.fft.ifft(fs)
    us= np.real(s[int((lengthbuffer-1)*numberofdays/2*oversamplerate): int((lengthbuffer+1)*numberofdays/2*oversamplerate): oversamplerate])

    if normalizevariance:
        us-=np.mean(us)
        rms=np.mean(us**2)**.5
        if rms>0:
            print(f"rms was {rms}")
            us/=rms
        us*=envimeanvariance

    # return s, us, fs
    envi=np.zeros((numberofbins,numberofdays))
    x=np.linspace(-1,1,numberofbins)

    for t in range(numberofdays):
        envi[:,t]=sci.norm.pdf(x, us[t], envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate

    return(envi, us)
    
def meantofullenvi(envimeans, envimeanvariance, envivariance, numberofbins=100, maxsurvivalrate=1, normalizevariance=True):
    numberofdays=envimeans.shape[0]
    envi=np.zeros([numberofbins, numberofdays])
    x=np.linspace(-1,1,numberofbins)

    if normalizevariance:
        envimeans-=np.mean(envimeans)
        rms=np.mean(envimeans**2)**.5
        if rms>0:
            envimeans/=rms
        envimeans*=envimeanvariance

    for t in range(numberofdays):
        envi[:,t]=sci.norm.pdf(x, envimeans[t], envivariance) # A gaussian of environment with center around 0
        envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate
    
    return envi

# def driftmodeling(flynum, numberofbins, numberofdays, prefmean, prefvariance, envimean, envivariance, driftvariance, adaptivetracking, gain, per, maxsurvivalrate, birthrate, matureage, percentbh, showgraphs, figuresavepath):
    # adaptivetracking=0

def driftmodeling(envi, prefmean=0, prefvariance=0, driftvariance=0, adaptivetracking=0, birthrate=40,
                   matureage=10, percentbh=.1, showgraphs=False, figuresavepath='', driftmaxdistribution=.3,
                     envimeanvariance=1, envivariance=1, numberofbins=1000, savealldays=True, saveallprefs=False, phi=1, deathscale=1):

    if len(envi.squeeze().shape)==1:

        envi=meantofullenvi(envimeans=envi, envimeanvariance=envimeanvariance, envivariance=envivariance, numberofbins=numberofbins) 
        # print(envi.shape)
        # fig.

    
    flynum=1
    numberofbins=envi.shape[0]
    numberofdays=envi.shape[1]

    x=np.linspace(-1,1,numberofbins) # number of bins between -1 and 1
    maxage=2*matureage+10 #maximum age
    prefvariance=np.array([prefvariance])
    prefmean=np.array([prefmean])
    driftvariance=np.array([driftvariance])
    adaptivetracking=np.array([adaptivetracking])
    numconditions=max(prefvariance.shape) # Number of conditions based on the total conditions we're running
    finalpop=np.zeros((numconditions))
    #print(np.floor(numberofbins/2))

    for q in range(numconditions): # for loop for each condition\
        # print(numconditions)
        # print(prefvariance[q])
        # print(driftvariance[q])
        pref=np.zeros((numberofbins,numberofdays,maxage)) # Matrix, [bins, days, maxage]
        #reducebethedge=np.zeros((numberofbins,numberofdays,maxage,2))
        # Set pref[:,0,0,0], which is the "reduced bet hedge version"
        # if prefvariance[q]>=0.015:  #Check if variance is so small to just eliminate bet-hedging
        if prefvariance[q]>=.05/numberofbins:  #Check if variance is so small to just eliminate bet-hedging

            pref[:,0,0]=sci.norm.pdf(x,prefmean[q],prefvariance[q]) # A fly's first day preference gaussian of preference with center around 0
        else: # Make the bin in the middle have all the flies
            #print('Zero bet-hedging')
            #pref[50,0,0,0]=flynum
            pref[math.floor(numberofbins/2),0,0]=flynum
            #print(pref[:,0,0,0])
#         print(np.sum(pref[:,0,0]))
        pref[:,0,0]=pref[:,0,0]/np.sum(pref[:,0,0])*flynum # total # of flies=flynum
        
        #This line is to tweak the initial preferences to spread out the first generation to avoid generation aliasing
        for i in range(1, matureage):
          pref[:,0,i]=pref[:,0,0]/matureage

        pref[:,0,0]=pref[:,0,0]/matureage

        blur=np.zeros((numberofbins,numberofbins)) # [which bin profile it's for, what the distribution is between that bin and all other bins, 2]
        # print(blur.shape)
        for b in range(numberofbins):
            if driftvariance[q]>0:
                blur[b,:]=sci.norm.pdf(x,x[b],driftvariance[q])
            else:
                blur[b, b]=1
            # print(blur)
            blur[b,:]/=np.sum(blur[b,:])
        if driftmaxdistribution!=0:
            for b in range(numberofbins):
                # boundingdistribution=sci.norm.pdf(x,0,driftmaxdistribution)
                blur[b,:]=blur[b,:]*sci.norm.pdf(x,0,driftmaxdistribution)
                if np.sum(blur[b,:]*sci.norm.pdf(x,0,driftmaxdistribution))>0:
                    # blur[b,:]=blur[b,:]*sci.norm.pdf(x,0,driftmaxdistribution)
                    blur[b,:]/=np.sum(blur[b,:])
        
        phishift=np.zeros([numberofbins, numberofbins])
        # print(numberofbins)
        for b in range(numberofbins):
            # print(b)
            newmean=x[b]*phi
            # print(newmean)
            closestbin=int(np.searchsorted(x,newmean, 'left'))
            # print("closestbin%f "%closestbin)
            closestbinright=int(min(numberofbins-1, np.searchsorted(x,newmean, 'right')))
            # print("closestbin%f "%closestbinright)

            lefterror=x[closestbin]-newmean
            righterror=x[closestbinright]-newmean
            # print(lefterror)
            # print(righterror)
            if righterror>0:
                phishift[b,closestbin]=1*righterror/(lefterror+righterror)
            if lefterror>0:
                phishift[b,min(numberofbins, closestbinright)]=1*lefterror/(lefterror+righterror)
        
        # fig, (ax1, ax2, ax3)=plt.subplots(1,3)
        # ax1.pcolormesh(phishift)
        # ax2.pcolormesh(blur)
        # blur=(phishift@blur)*2
        # ax3.pcolormesh(blur)
        # fig
        # return fig
        # fig.show()


        # return blur

        for t in range(1,numberofdays):

            # envi[:,t]=sci.norm.pdf(x,(envimean+gain*np.sin(t*np.pi*2/per)),envivariance) # Making envi a sin wave that changes over time
            # envi[:,t]=envi[:,t]/np.max(envi[:,t])*maxsurvivalrate # Normalizing the envi and multiplying by maxsurvival rate
            # # print('t is: '+str(t))
            # for w in range(2):
#             print(pref[:,t-1,0:-1])
            if driftvariance > .05/numberofbins/10:
              pref[:,t,1:]=blur[:,:].T @ pref[:,t-1,0:-1] #drift
            elif driftvariance > 0:
              pref[:,t,1:]=pref[:,t-1,0:-1]
              print('rounded to zero')
            else:
              pref[:,t,1:]=pref[:,t-1,0:-1]
#             pref[:,t,:]=np.multiply(pref[:,t,:], envi[:,t]) # Multiplying the preference to the environment
#             deadflies=np.sum(pref[:,t,:])
            pref[:,t,:]=envi[:,t:t+1] * pref[:,t,:]
#             print("Dead flies")
#             print(deadflies-np.sum(pref[:,t,:]))
#             print(envi[:,t])
            pref[:,t,0]=pref[:,0,0]*birthrate/flynum*np.sum(pref[:,t-1,matureage:]) # Calculate newborn flies
            if adaptivetracking[q]>0:
                pref[:,t,0]=pref[:,t,0]*(1-adaptivetracking[q])+adaptivetracking[q]*birthrate*np.sum(pref[:,t-1,matureage:],1)


        if showgraphs:
            #before = time.perf_counter()
            fig, (ax0, ax1,  ax1d, ax2) = plt.subplots(4, 1)
#             fig, (ax0, ax1,  ax1d, ax2, ax3, ax4) = plt.subplots(6, 1)

            fig.set_figwidth(10)
            fig.set_figheight(12)
            fig.tight_layout()
            plt.subplots_adjust(hspace=.6)
            c=ax0.pcolormesev=h(envi)
            fig.colorbar(c,ax=ax0)
            ax0.set_title('Environment (color is fraction of flies of given pref survive)')
            ax0.set_ylabel('Preference')
            ax0.set_xlabel('Day')

            c=ax1.pcolormesh(np.sum(pref[:,:,:],axis=2))
            fig.colorbar(c,ax=ax1)
            ax1.set_title('Fly Preference (color is log(num) flies each day)')
            ax1.set_ylabel('Preference')
            ax1.set_xlabel('Day')

            c=ax1d.pcolormesh(np.sum(pref[:,:,:],axis=2)/np.max(np.sum(pref[:,:,:],axis=2),axis=0))
            fig.colorbar(c,ax=ax1d)
            ax1d.set_title('Fly Preference Distribution (color is %daily flies)')
            ax1d.set_ylabel('Preference')
            ax1d.set_xlabel('Day')

            ax2.plot(np.log(np.sum(pref[:,:,:],axis=(0,2))))
            ax2.set_title('total log(num) flies)') # lowest value is 0.0001 (prefvariance = 0.01 with percent bh=0.01)
            ax2.set_ylabel('log(num) flies)')
            ax2.set_xlabel('Day')
            ax2.set_xlim(0,numberofdays)

            # ax3.plot(driftadvantage) #Commented out because not used for simulations
#             ax3.plot(driftadvantage/np.sum(pref[:,:,:],axis=(0,2)))
#             ax3.set_title('Change in death rate due to last day\'s drift ')
#             ax3.set_ylabel('∆surviving flies/total flies')
#             ax3.set_xlabel('Day')
#             ax3.set_xlim(0,numberofdays)
# 
#             ax4.plot(betadvantage/np.sum(pref[:,:,:],axis=(0,2)))
#             ax4.set_title('Change in death rate due to last day\'s bethedging ')
#             ax4.set_ylabel('∆surviving flies/total flies')
#             ax4.set_xlabel('Day')
#             ax4.set_xlim(0,numberofdays)

            fig.colorbar(c,ax=ax2)
#             fig.colorbar(c,ax=ax3)
#             fig.colorbar(c,ax=ax4)

            fig.suptitle('Bet-hedge variance: '+str(prefvariance[q])+', Drift variance: '+str(driftvariance[q])+', Adaptive Tracking: '+str(adaptivetracking[q]), y=-.05, fontsize=16)

            plt.show()

            if os.path.exists(figuresavepath):
                fig.savefig(os.path.join(figuresavepath,'bh'+str(prefvariance[q])+'dv'+str(driftvariance[q])+'.png'),bbox_inches='tight', pad_inches=.3)
                print(figuresavepath)
            else:
                print('not saving, no valid path')
        #after = time.perf_counter()
        #print(after-before)
        if savealldays & ~saveallprefs:
            finalpop=pd.Series(np.sum(pref[:,:,:], axis=(0,2)), name='Populations')
        elif savealldays & saveallprefs:
            # finalpop=xr.DataArray(pref, dims=["", "Day", ""])
            finalpop=pd.DataFrame(np.sum(pref[:,:,:], axis=2))
        else:
            finalpop=pd.Series(np.sum(pref[:,-1,:], axis=(0,1)), name='Populations')

    return finalpop