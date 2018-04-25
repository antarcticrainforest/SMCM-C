#!/usr/bin/env python2.7


import pandas as pd
import sys,os,datetime
import numpy as np
from Driver import Reader
from netCDF4 import num2date
from Master.configdir import Config
from run_cloudmodel import GetNames



def get_df(box):
    start = datetime.datetime(1998,1,15,0,0)
    end = datetime.datetime(2015,8,15,18,0)
    dt = (end - start).days
    R=Reader(start,dt,loc=box.lower())
    
    times,tr = R.boxdata('trigger',freq='6H')
    times = pd.DatetimeIndex(num2date(times,'Minutes since 1998-01-01 00:00'))
    qi = pd.DataFrame({'qi':2*R.boxdata('qi',freq='6H')[1]},index=times).resample('D').mean()
    ki = pd.DataFrame({'ki':R.boxdata('kindex',freq='6H')[1]/12.5},index=times).resample('D').mean()
    precip = R.boxdata('precip',freq='3H')
    p_times=pd.DatetimeIndex(num2date(precip[0],'Minutes since 1998-01-01 00:00'))
    pi = pd.DataFrame({'pi':precip[1]},index=p_times).resample('D').sum()
    tr=pd.DataFrame({'tr':np.fabs(tr)},index=times).resample('D').mean()
    return pd.concat([qi,ki,pi,tr],axis=1).dropna()
    
def dayli_mean(df):

    qi=df['qi'].quantile(.3)
    ki=df['ki'].quantile(.3)
    pi=df['pi'].quantile(.7)
    data=df.resample('M').mean()

    


    return qi,ki,pi,data#df.resample('M').mean()

def plot(df,region,loc,name):

    if len(df.index[:]) == 0:
        return
    elif df['tr'].mean() == 0:
        return
    elif len(df.index[:]) > 3:
        df = df.loc[df.index[0:3]]
    
    import matplotlib
    from matplotlib import pyplot as plt,dates
    from calendar import monthrange
    from matplotlib import pyplot as plt
    
    #Plot the candidates
    num=len(df.index)
    fig = plt.figure()
    T=[]
    Data=[]
    for date  in df.index[:]:
        try:
            year,month = date.year,date.month
            ndays = monthrange(date.year,date.month)[-1]
            start = datetime.datetime(date.year,date.month,1,0,0)
            R = Reader(start,ndays,loc=loc)
            qi = R.boxdata('qi',freq='6H')[1]
            times_i,ki = R.boxdata('kindex',freq='6H')
            times_r,rain = R.boxdata('precip',freq='3H')
            times_tr,trigger = R.trigger()
            Data.append((qi,ki/22.5,rain,trigger))
            T.append((times_i,times_r,times_tr))
        except (IndexError,OSError):
            pass
    N = len(Data)
    nn = 1
    for i in xrange(len(Data)):
        qi,ki,rain,trigger=Data[i]
        times_i,times_r,times_tr = T[i]


        tday_tr = num2date(times_tr,'Minutes since 1998-01-01 00:00:00')
        tday_r = num2date(times_r,'Minutes since 1998-01-01 00:00:00')
        tday_i = num2date(times_i,'Minutes since 1998-01-01 00:00:00')
        time_tr = dates.date2num(tday_tr)
        time_r = dates.date2num(tday_r)
        time_i = dates.date2num(tday_i)
        
        
        ax = fig.add_subplot(1,N,nn)
        axes = [ax,ax.twinx(),ax.twinx()]
        axes[-1].spines['right'].set_position(('axes', -0.1))
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)
        axes[-1].tick_params(axis='y', colors='g')
        
        a1 = axes[0].plot(time_r,rain,'k-',lw=3.5,label='Rain')
        a2 = axes[1].plot(time_i,qi,'b-',lw=3.5,label='Moisture (right)')
        a3 = axes[1].plot(time_i,ki,'r-',lw=3.5,label='Instability (right)')
        a4 = axes[2].plot(time_tr,trigger,'g-',lw=3.5,label='Trigger function')
        
        axes[2].set_ylabel('Coastal effects []',color='g',labelpad=-80)
        axes[2].tick_params(axis='y', colors='g')
        axes[1].set_ylabel('Moisture/Instability []')
        axes[0].set_ylabel('Rain-rate [mm/3h]')
        
        axes[1].set_ylim(0,2)
        axes[0].set_ylim(0,rain.max()+0.05*rain.max())
        #axes[2].set_ylim(trigger.min()-0.05*trigger.min(),trigger.max()*0.05*trigger.max())
        ax.set_xlim(time_tr.min(),time_tr.max())
        hfmt = dates.DateFormatter('%d/%m/%y %H LT')
        ax.xaxis.set_major_formatter(hfmt)
        font = {'family' : 'normal','weight' : 'normal','size'   : 20}
        ax.set_title('Conditions in %s'%name)
        lns=a1+a2+a3+a4
        labs = [l.get_label() for l in lns]
        plt.legend(lns,labs,loc=2)
        nn += 1
        matplotlib.rc('font', **font)
    plt.show()
    #exit()
def main(boxes):
    Cfg=Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'boxes.txt'))
    for b in boxes:
        name=GetNames(Cfg[b.lower()])
        #name=b.lower()
        sys.stdout.flush()
        sys.stdout.write("Working on %s (%s) ... "%(name,b.lower()))
        sys.stdout.flush()
        
        p_q,p_k,p_p,df = dayli_mean(get_df(b))
        

        sys.stdout.write('ok\n')

        ps = df.loc[(df['qi'] <= p_q)  & (df['ki'] <= p_k) &(df['pi'] >= p_p)]
        plot(ps,name,b,name)
        if len(ps.index[:]) >= 1:
            print ps.sort_values(['tr','pi'],ascending=False)
    return



if __name__ == "__main__":

    import sys

    try:
        boxes=sys.argv[1].split(',')
    except IndexError:
        boxes=['coast_%02i'%i for i in xrange(1,12)]
    

    main(boxes)
