#!/usr/bin/env python2.7

import os,sys,datetime
import pandas as pd
import numpy as np
from netCDF4 import num2date,date2num,Dataset as nc
from Driver import check
import Plot
from timezonefinder import TimezoneFinder

def __writePrecip(*args):
    """
    Helper function to write cmoroph data to netcdf-file
    """
    database = args[0]
    f = args[1]
    stype=args[2]
    
    lookup = dict(sea=' over ocean',land=' over land',all='')
    precip_kw={'precip_%s'%stype:'total','lsp_%s'%stype:'coastally affected'}
    try:
        f.createDimension('precip_time',None)
    except (ValueError,RuntimeError):
        pass
    for i in database.keys():
        try:
            f.createVariable(i,'f',('precip_time',))
            f.variables[i].units='mm/3h'
            f.variables[i].standard_name=precip_kw[i]+' precipitation'+lookup[stype]
            f.variables[i][:]=database[i].values
        except (ValueError,RuntimeError):
            f.variables[i].units='mm/3h'
            f.variables[i].standard_name=precip_kw[i]+' precipitation'+lookup[stype]
            f.variables[i][:]=database[i].values
    try:
        f.createVariable('precip_time','f',('precip_time',))
    except (ValueError,RuntimeError):
        pass
    f.variables['precip_time'].units='Minutes since 1998-01-01 00:00'
    f.variables['precip_time'].standard_name='time'
    f.variables['precip_time'].axis='T'
    f.variables['precip_time'][:]=date2num(database.index.to_pydatetime(),
            f.variables['precip_time'].units)



##############################################################################
##############################################################################

def __readPrecip(*args,**kwargs):
    """
    Helper function to read the cmorph data
    """
    box = args[0]
    f = args[1]
    if 'slmf' in kwargs.keys():
        from Master.configdir import Config
        from Master.Useful import LonLatBox
        clat,clon = Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),'boxes.txt'))[box]
        from LargeScale.Distribute import GetRect
        rect=GetRect((clat,clon,box))
        lats,lons,slat,slon,elat,elon=LonLatBox(kwargs['slmf'],box=rect)
        mask = nc(kwargs['slmf']).variables['slm'][slat:elat,slon:elon]
    else:
        stype='all'
        mask= 1
    if 'stype' in kwargs.keys():
        stype=kwargs['stype']
        if stype == 'total':
            stype='all'
        elif stype == 'ocean':
            stype = 'sea'
        tmp=dict(land=1,sea=0,all=None)
        if tmp[stype] != type(None):
            mask = np.ma.masked_not_equal(mask,tmp[stype])
            num=mask.mean()
            if num == 0:
                mask+=1
        else:
            mask = np.ones_like(mask)

    else:
        mask = 1
        stype='all'

    mask = 1
    #Open the precip data for correlcatoin
    precipdir = os.path.join(os.environ['HOME'],'PhD','Data','LargeScale')
    precipfile = os.path.join(precipdir,'Database_fm_%s.pkl'%box.lower())
    #Get the start time
    try:
        times = f.variables['time'][:].filled(-1)
        times = times[times>0]
    except AttributeError:
        times = f.variables['time'][:]
    times = num2date(times,f.variables['time'].units)
    if 'slmf' in kwargs.keys():
        #Read the data from cmorph
        fnames=[]
        precip=[]
        lsp=[]
        ptimes=[]
        datadir=os.path.join(os.environ['HOME'],'Data','CMORPH')
        ldatadir=os.path.join(os.environ['HOME'],'Data','PatternDetect','CMORPH')
        for t in times:
            fn=os.path.join(datadir,'netcdf',str(t.year),'%s'%t.strftime('Cmorph_%Y_%m_%d.nc'))
            fl=fn.replace('Cmorph_','CmorphPattern_').replace(datadir,ldatadir)
            if os.path.isfile(fn) and os.path.isfile(fl) and fn not in fnames:
                fnames.append(fn)
                with nc(fn) as fp:
                    with nc(fl) as ffl:
                        units=fp.variables['time'].units
                        for tt in xrange(len(fp.variables['precip'][:])):
                            ptimes.append(num2date(fp.variables['time'][tt],units))
                            precip.append(\
                                    (fp.variables['precip'][tt,slat:elat,slon:elon]\
                                    *mask).mean())
                            lsp.append(\
                                    (ffl.variables['lsp'][tt,slat:elat,slon:elon]\
                                    *mask).mean())

        ptimes=pd.DatetimeIndex(ptimes)
        precip=np.array(precip)
        lsp=np.array(lsp)

        data=pd.DataFrame({'lsp_%s'%stype:lsp,'precip_%s'%stype:precip},index=ptimes)

    else:


        start,end = pd.DatetimeIndex([times[0],times[-1]])
        database = pd.read_pickle(precipfile)
        si = np.fabs((database.index - start).total_seconds())
        si = np.argmin(si)

        ei = np.fabs((database.index - end).total_seconds())
        ei = np.argmin(ei)+1
    
        data=database[['precip','lsp']].loc[database.index[si:ei]]
        data.columns=['precip_all','lsp_all']
    
    __writePrecip(data,f,stype)

##############################################################################
##############################################################################

##############################################################################
##############################################################################
def __getTimeS(f,stype):
    """
    Helpber function to create a timeseries from the setups

    Variables:
        f (netcdf-ojbect) the netcdf file-object where the data is stored
    """

    g = f['interact']
    g2 = f['nointeract']
    smcm={}
    for i in ('con','deep','strat'):
        smcm[i]=g2.variables[i][:].mean(axis=(1,2))
    
    T = pd.DatetimeIndex(num2date(f.variables['time'][:],f.variables['time'].units))
    shape=g2.variables['deep'][0].shape
    mask_l=np.ones(shape)
    mask_o=np.ones_like(mask_l)
    mask_o[:,:shape[0]/2]=0
    mask_l[:,shape[0]/2:]=0
    mask_l = np.ma.masked_equal(mask_l,0)
    mask_o = np.ma.masked_equal(mask_o,0)
    if str(stype) == 'land':
        mask = mask_l
    elif str(stype) == 'ocean':
        mask = mask_o
    else:
        mask = np.ones(shape)
        
    smcm = pd.DataFrame(smcm,index=T).tz_localize(f.tz)
    try:
        Time = f.variables['time'][:].filled(-1)
        time = Time[Time>0]
    except AttributeError:
        Time = f.variables['time'][:]
        time =f.variables['time'][:]
    time = pd.DatetimeIndex(num2date(time,f.variables['time'].units))

    try:
        new_time = f.variables['precip_time'][:]
        unit = f.variables['precip_time'].units
    except KeyError:
        new_time = f.variables['ptime'][:]
        unit = f.variables['ptime'].units
    new_time = pd.DatetimeIndex(num2date(new_time,unit))

    setup = g.groups
    con,deep,strat={},{},{}
    nlon = int(len(f.variables['lon'][:]))/2 + 2
    for s in xrange(len(setup)):
        #try:
        if True:
            G = g["%05i"%s]
            sys.stdout.write('Adding %4i / %4i (%3.2f %%)\n'%(s+1,len(setup),
                100 * float(s)/(len(setup)-1)))
            for ary,form in ((con,'con'),(deep,'deep'),(strat,'strat')):
                d = (G.variables[form][:]*mask).mean(axis=2).mean(axis=1)
                
                ary[s] = pd.Series(d,index=time).tz_localize(f.tz)

    
    precip = pd.DataFrame(f.variables['precip_%s'%stype][:],index=new_time,\
            columns=['precip']).tz_localize('UTC').tz_convert(f.tz)
    data = {'con':pd.DataFrame(con),'deep':pd.DataFrame(deep),\
            'strat':pd.DataFrame(strat),'precip':precip,'smcm':smcm}
    params = pd.DataFrame(f.variables['params'][:],columns=['mul','add_c'])
    
    df = pd.concat(data.values(),axis=1,keys=data.keys())
    
    return df,params

def __corr(df,params,lag=-3,mins=90):

    """
    Method to calculate the correlation of lag lag

    """
    new_time = df['precip'].dropna().index
    start=new_time-pd.Timedelta('%i min'%mins)
    end=new_time+pd.Timedelta('%i min'%mins)
    
    start=new_time-pd.Timedelta('%i min'%mins)
    end=new_time+pd.Timedelta('%i min'%mins)
    new_index=[]
    data={'con':[],'deep':[],'strat':[]}
    smcm={'con':[],'deep':[],'strat':[]}
    
    for i in xrange(len(new_time)):
        st,en=start[i],end[i]
        for cl in data.keys():
            data[cl].append(df[cl].loc[st:en].values.mean(axis=0))
            smcm[cl].append(df.smcm[cl].loc[st:en].values.mean(axis=0))
    for i in data.keys():
        
        data[i] = pd.DataFrame(data[i],index=df.precip.dropna().index).shift(lag)
        smcm[i] = pd.DataFrame(smcm[i],index=df.precip.dropna().index).shift(lag)
    
    data['rain']=df.precip.dropna()
    data=pd.concat(data.values(),axis=1,keys=data.keys()).dropna()
    corr={}
    mul={'con':1,'strat':1,'deep':1}
    for cl in ('con','deep','strat'):
        corr['%s_corr'%cl]=[mul[cl]*data[cl][i].corr(data.rain.precip) for i in data[cl].keys()]
        corr['%s_mean'%cl]=df[cl].mean(axis=0).values
    corr['mul']=params['mul'].values
    corr['add_c']=params['add_c'].values
    corr = pd.DataFrame(corr,index=data['con'].keys())
    t = 'deep'
    R = (params['mul'] >= 0.5) & (params['add_c'] >= 0.3)
    coastal = corr['%s_corr'%t].loc[R].abs().argmax()
    no = corr['%s_corr'%t].abs().argmin()
    no = 284
    noC = smcm[t][0]#.shift(-10)
    
    print(coastal,corr['%s_corr'%t][coastal],params['mul'][coastal],params['add_c'][coastal])
    rain = pd.concat([data[t][coastal],noC,data['rain']['precip']],axis=1,\
            keys=('Coastal','Original','Rain'))
    #S = pd.concat([data[t],data['rain']['precip']],axis=1,keys=('SMCM','Rain'))
    #S.to_pickle('/home/wilfred/test_sim.pkl')
    #exit()

    return corr,rain

def __avg(d1,d2,la,avg,cltype,Dt):
    """
        Method to create an average of avg hours of the cloud data set and
        combine it with the largescale atmosphere data
    """

    df1 = pd.read_pickle(d1)[[cltype,'thc']]
    df1.columns = ['%s_coastal'%cltype,'thc']

    df2 = pd.read_pickle(d2)[cltype]
    df2=pd.DataFrame(df2.values,index=df2.index,columns=['%s_original'%cltype])
    df = pd.concat([df1,df2],axis=1)
    
    dfla = pd.read_pickle(la)
    
    dt = int((df.index[1]-df.index[0]).total_seconds()/60.)
    itt = range(0,avg+dt,dt)
    end = len(itt)
    T = np.arange(0,24,3)
    a = T[np.argmin(np.fabs(T-df.index[end].hour))]
    b = T[np.argmin(np.fabs(T-df.index[-1].hour))]
    s = '%04i-%02i-%02i %02i:00' %(df.index[end].year,df.index[end].month,
            df.index[end].day,a)
    e = '%04i-%02i-%02i %02i:00' %(df.index[-1].year,df.index[-1].month,
            df.index[-1].day,a)

    times = pd.date_range(s,e,freq='3 H')
    df['thc']=df['thc'].abs()
    df_1=df.rolling(len(itt)).mean().loc[times]
    
    df=(pd.concat([df_1,dfla[['kindex','wwind_600','qi']]],axis=1)).dropna()
    df.index=np.arange(Dt,Dt+len(df.index))
    Dt+=len(df.index)
    return df,Dt

################################################################################
def variability(**kwargs):
    """
        Method to calculte the spectrum of the simulated clouds and compar+e it
        to obs
    """

    from Master.configdir import Config
    smcmdir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    armdir = os.path.join(os.environ['HOME'],'Data','Darwin','ARM','179567')
    head = "twp15swfanalskyrad1longC3.c1*.cdf"
    target = os.path.join(smcmdir,'arm_data.pkl')
    if not os.path.isfile(target):
        import glob
        fnames = glob.glob(os.path.join(armdir,head))
        caf = None
        for nn,fn in enumerate(fnames):
            with nc(fn,'a') as f:
                time=f.variables['base_time']
                T = pd.DatetimeIndex(num2date(time[:]+f.variables['time_offset']\
                        ,time.units),tz='UTC')
                f.variables['cloudfraction'].missing_value=-999.
                data = pd.DataFrame({'caf':f.variables['cloudfraction'][:]},index=T)
            if type(caf) == type(None):
                caf = data
            else:
                caf=pd.concat([caf,data],axis=0).sort_index()
        

        caf.to_pickle(target)
    else:
        caf = pd.read_pickle(target).tz_convert('Australia/Darwin')
    print(caf)

    box = 'coast_06'
    interact ="Areamean-Interact"
    nointeract = "Areamean-Nointeract"
    d_i = os.path.join(smcmdir,"%s_%s.pkl"%(interact,box))
    d_ni = os.path.join(smcmdir,"%s_%s.pkl"%(nointeract,box))
    fnames={'Coastal':d_i,'Original':d_ni}
    intp=True
    for typ,fn in fnames.iteritems():
        smcm=pd.read_pickle(fn).tz_localize('Australia/Darwin')
        df = smcm['deep']+smcm['strat']+smcm['con']
        if intp :
            times=pd.date_range(df.index[0],df.index[-1],freq='5min')
            intp=False
        
        df = df.reindex(times).interpolate(method='linear')
        
        fnames[typ]=df.rolling(3).mean().loc[caf.index]

    fnames['Obs']=caf.caf
    de = 3 * 60 
    T = pd.date_range(caf.index[0],caf.index[-1],freq='%imin'%de)
    data = pd.DataFrame(fnames).reindex(T).interpolate('linear')
    from scipy import stats,signal
    import matplotlib
    from matplotlib import pyplot as plt
    
    #ax = plt.subplot2grid((1,1),(1,1),colspan=1)
    plt.gca().invert_xaxis()
    font = {'family' : 'sans','size'   : 24}
    matplotlib.rc('font', **font)
    co={'Original':('k',0.001),'Coastal':('g',0.05),"Obs":('b',2)}
    for i in ['Obs','Coastal','Original']:
        f,psd = signal.welch(data[i].values,nperseg=data[i].values.size)
        S = psd/psd.sum() * 100
        f = (de/f)/(60*24)
        """
        if i == 'Coastal':

            idx=np.where(f>50)
            Min=np.argmin(psd[idx])
            psd[Min]*=250
            idx=(f>6)
            Min=np.argmin(psd[idx])
            psd[Min]*=250
        """
        plt.plot(f,S,\
                '-',color=co[i][0],label=i)
        #plt.xscale('log')
        #plt.yscale('log')
        plt.ylabel('Spectral variance [1/day]',size=32)
        plt.xlabel('1/f [days]',size=32)
        plt.xlim(700,0.3)
    plt.legend(loc=0)
    plt.show()

################################################################################
def large_scale(**kwargs):
    """
        Method to create a plot of the large-scale to cloudiness relationship
    """

    from Master.configdir import Config

    datadir = os.path.join(os.environ['HOME'],'PhD','Data','LargeScale')
    smcmdir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    interact ="Areamean-Interact"
    nointeract = "Areamean-Nointeract"
    C = Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'boxes.txt'))
    avg = int(check(kwargs,'avgtime',120))
    cltype = check(kwargs,'cltype','deep')
    Data = None
    dt = 1
    for box in ('coast_%02i'%i for i in xrange(1,12)):
        d_i = os.path.join(smcmdir,"%s_%s.pkl"%(interact,box))
        d_ni = os.path.join(smcmdir,"%s_%s.pkl"%(nointeract,box))
        la = os.path.join(datadir,'Database_fm_%s.pkl'%box)
        if os.path.isfile(d_i) and os.path.isfile(d_ni) :
            lat,lon = C[box]
            if type(Data) == type(None) :
                Data,dt = __avg(d_i,d_ni,la,avg,cltype,dt)

            else:
                data,dt = __avg(d_i,d_ni,la,avg,cltype,dt)
                Data=pd.concat([Data, data],axis=0)
            return Plot.large_scale(Data)
    


        

    
################################################################################
################################################################################
def stochasticity(**kwargs):
    """
        Method to plot the influence of stochasticity
    """

    from Master.configdir import Config
    from matplotlib import pyplot as plt,dates
    from netCDF4 import Dataset as nc, num2date

    C = Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'constants.config'))
    box = C.obs.upper()
    cltype=C.animate
    if str(cltype).lower() == 'none':
        cltype='deep'


    dirn=os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    if not os.path.isdir(dirn):
        dirn=os.path.join(os.environ['HOME'],'Data','SMCM')

    fn = os.path.join(dirn,'Diurnalcycle_%s.nc'%box.upper())
    
    F = nc(fn)
    time=pd.DatetimeIndex(num2date(F.variables['time'][:],F.variables['time'].units))
    out={}
    for g in F.groups:
        if g.lower() != 'none':
            c=F[g].variables[cltype][:].mean(axis=(1,2))
            mc=1-c.mean()
            c=((c+mc)**3.5-mc - 0.03 ) * 0.85/2
            c[c<0]=0
            out[g]=c
        else:
            c=F[g].variables[cltype][:].mean(axis=(1,2))
            mc=1-c.mean()
            c=((c+mc)**3.5-mc - 0.03 ) * 0.85/2
            c[c<0]=0
            norm = c

    F.close()

    df=pd.DataFrame(out,index=time,columns=out.keys())
    N = pd.Series(norm,index=time)

    Min,Max=[],[]
    for t in time:
        Min.append(df.loc[t].min())
        Max.append(df.loc[t].max())

    hfmt = dates.DateFormatter('%d/%m/%y %H LT')
    fig,ax = plt.subplots()
    ax.xaxis.set_major_formatter(hfmt)
    ax.grid()
    ax.plot(N.index, N.values, 'k')
    ax.fill_between(N.index, Min, Max, color='b', alpha=0.2)
    ax.set_ylabel('CAF of %s clouds'%cltype)



    #ax.set_xlim(tday[0],tday[idx][-1])
    #ax3[0].set_ylim(0,3)
    #ax2.set_ylim(0.06,0.14)
    #ax3[-1].set_ylim(-1.8,1.8)
    font = {'family' : 'sans','weight' : 'normal','size'   : 34}
    import matplotlib
    fig.autofmt_xdate()
    ax.set_xlim(time[0],time[-1])
    #plt.legend(lns,labs,bbox_to_anchor=(0,-.02,1,-0.102),mode='expand',ncol=4, loc=2, borderaxespad=0.,prop={'size':22})
    #plt.legend(lns, labs, loc=2)
    matplotlib.rc('font', **font)
    plt.show()



################################################################################
################################################################################
def diurnal_cycle(**kwargs):
    """
        Method to plot the diurnal cycle in a map
    """
    from Master.configdir import Config
    from matplotlib import pyplot as plt,dates
    import glob
    datadir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    C = Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'constants.config'))
    box = C.obs.upper()
    
    smcmdir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    armdir = os.path.join(os.environ['HOME'],'Data','Darwin','ARM','179567')
    precipdir=os.path.join(os.environ['HOME'],'PhD','Data')
    head = "twp15swfanalskyrad1longC3.c1*.cdf"
    target = os.path.join(precipdir,'Darwin','cpol_03-2005_cp.pkl')
    if not os.path.isfile(target):
        fnames = glob.glob(os.path.join(armdir,head))
        obs = None
        for nn,fn in enumerate(fnames):
            with nc(fn,'a') as f:
                time=f.variables['base_time']
                T = pd.DatetimeIndex(num2date(time[:]+f.variables['time_offset']\
                        ,time.units),tz='UTC')
                f.variables['cloudfraction'].missing_value=-999.
                data = pd.DataFrame({'caf':f.variables['cloudfraction'][:]},index=T)
            if type(obs) == type(None):
                obs = data
            else:
                obs=pd.concat([obs,data],axis=0).sort_index()
        

        obs.to_pickle(target)
    else:
        obs = pd.read_pickle(target)#.tz_localize('UTC').tz_convert('Australia/Darwin')
    #obs=obs.shift(11*6)
    pdata=pd.read_csv(os.path.join(precipdir,'Database_%s.csv'%C.obs.upper()))
    precip=pdata['precip'].shift(6).dropna()
    precip.index=num2date(precip.index*3,'Hours since 1998-01-01 00:00:00')

    
    
    start=datetime.datetime.strptime(C.start.strip(),'%Y-%m-%d_%H:%M')
    end=start+datetime.timedelta(hours=int(C.tend))
    obs.index=pd.DatetimeIndex([i.replace(tzinfo=None) for i in obs.index])
    obs = obs[start:end]
    fname = 'Diurnalcycle_%s_%s-*00.nc'%(box,start.strftime('%Y.%m.%d_%H_%M'))
    fname = glob.glob(os.path.join(datadir,fname))[0]

    if not os.path.isfile(fname):
        sys.stderr.write('File %s does not exsist....exit\n'%fname)
    
    hfmt = dates.DateFormatter('%d/%m/%y %H LT')
    from Driver import Reader

    R=Reader(C.start.replace('_',' '),C.tend/24.,loc=C.obs,method='linear')
    
    group='%1.2f'%C.J00
    interact='nointeract'
    print(fname)
    with nc(fname,'r') as nf :
        time = num2date(nf.variables['time'][:]+3*60,nf.variables['time'].units)
        tday= dates.date2num(time)
        g = nf['interact'][group]
        lon =  nf.variables['lon'][:]
        deep1 = g.variables['deep'][:].mean(axis=-2)
        strat1 = g.variables['deep'][:].mean(axis=-2)
        deep2 = nf['nointeract'][group].variables['deep'][:].mean(axis=-2)
        strat2 = nf['nointeract'][group].variables['deep'][:].mean(axis=-2)
        caf=((deep1+strat1)/2,(deep2+strat2)/2)
    
    x2,tr = R.trigger()
    ki = R.boxdata('kindex')[-1]/12.5
    x3,qi = R.boxdata('qi',shift=0)
    qi *= 2
    #x3+=120
    tday2 = num2date(x2,'Minutes since 1998-01-01 00:00:00')
    tday3 = num2date(x3,'Minutes since 1998-01-01 00:00:00')
    #qi = obs.values
    tday31 = (obs.index+pd.Timedelta(hours=2)).to_pydatetime()
    tday1 = dates.date2num(obs.index.to_pydatetime())
    tday2 = dates.date2num(tday2)
    tday3 = dates.date2num(tday3)
    tday31 = dates.date2num(tday31)
    plotend = dates.date2num([datetime.datetime(2005,3,27,12)])[0]
    plotstart = dates.date2num([datetime.datetime(2005,3,22,0)])[0]
    idx = np.where(tday<plotend)
    idx1 = np.where(tday1<plotend)
    idx2 = np.where(tday2<plotend)
    idx3 = np.where(tday3<plotend)
    idx31 = np.where(tday31<plotend)
    fig = plt.figure()
    ax1 = plt.subplot2grid((3,4), (0,0), colspan=2,rowspan=2)
    #X,Y=np.meshgrid(lon,tday)
    v=np.linspace(0,0.2,10)
    pcolor = ax1.contourf(lon,tday[idx],caf[0][idx],v,extend='both',cmap='Blues')
    m=int(len(lon))
    dx=np.fabs((lon[1]-lon[1])/2.)
    ax1.set_xticks([])
    ax1_2 = plt.subplot2grid((3,4), (0,2), colspan=2,rowspan=2)
    pcolor_2 = ax1_2.contourf(lon,tday[idx],caf[1][idx],v,extend='both',cmap='Blues')
    if C.form == 'i':
        ax1.text(0.025,1.035,'Ocean',verticalalignment='top', horizontalalignment='left',\
                transform=ax1.transAxes,fontsize=24)
        ax1.text(0.5,1.035,'Land',verticalalignment='top', horizontalalignment='center',\
                transform=ax1.transAxes,fontsize=24)
        ax1.text(1-0.025,1.035,'Ocean',verticalalignment='top', horizontalalignment='right',\
                transform=ax1.transAxes,fontsize=24)
        ax1.axvline(x=lon[int(m/6)],linewidth=3,color='r')
        ax1.axvline(x=lon[5*int(m/6)],linewidth=3,color='r')
        ax1.yaxis.set_major_formatter(hfmt)

        ax1_2.text(0.025,1.035,'Ocean',verticalalignment='top', horizontalalignment='left',\
                transform=ax1_2.transAxes,fontsize=24)
        ax1_2.text(0.5,1.035,'Land',verticalalignment='top', horizontalalignment='center',\
                transform=ax1_2.transAxes,fontsize=24)
        ax1_2.text(1-0.025,1.035,'Ocean',verticalalignment='top', horizontalalignment='right',\
                transform=ax1_2.transAxes,fontsize=24)
        ax1_2.axvline(x=lon[int(m/6)],linewidth=3,color='r')
        ax1_2.axvline(x=lon[5*int(m/6)],linewidth=3,color='r')
    elif C.form == 'v':
        mm = (lon[-1]-lon[0])/2. + lon[0]
        ax1.text(0.25,1.035,'Ocean',verticalalignment='top', horizontalalignment='left',\
                transform=ax1.transAxes,fontsize=24)
        ax1.text(0.75,1.035,'Land',verticalalignment='top', horizontalalignment='center',\
                transform=ax1.transAxes,fontsize=24)
        ax1.axvline(x=mm,linewidth=3,color='r')
        ax1.yaxis.set_major_formatter(hfmt)

        ax1_2.text(0.25,1.035,'Ocean',verticalalignment='top', horizontalalignment='left',\
                transform=ax1_2.transAxes,fontsize=24)
        ax1_2.text(0.75,1.035,'Land',verticalalignment='top', horizontalalignment='center',\
                transform=ax1_2.transAxes,fontsize=24)
        ax1_2.axvline(x=mm,linewidth=3,color='r')

    ax1_2.set_xticks([])
    ax1_2.set_yticks([])

    cbar_ax = fig.add_axes([0.11, 1-0.64, 0.8, 0.02])
    cbar = fig.colorbar(pcolor,orientation='horizontal',\
            format='%.2f',cax=cbar_ax)
    
    shift=pd.DatetimeIndex(['2005-03-23 18:15'])
    s_idx=np.fabs((obs.index - shift[0]).astype(np.int64)).argmin()
    from matplotlib import patches
    import matplotlib.patches as mpatches
    ax2 = plt.subplot2grid((3,4),(2,0),colspan=4)
    ax3=[ax2.twinx(),ax2.twinx()]
    patch = [patches.Rectangle((tday[0],0),tday[s_idx],.09,alpha=0.2,facecolor='y',label='Deep West Regime'),
            patches.Rectangle((tday[s_idx],0),tday[-1],.09,alpha=0.2,facecolor='#d01c8b',label='East Regime')]
    l1=ax2.add_patch(patch[0])
    l2=ax2.add_patch(patch[1])

    #l1=mpatches.Patch(color='r',alpha=0.2,label='Deep West Regime')
    #l2=mpatches.Patch(color='y',alpha=0.2,label='East Regime')

    ax3[-1].spines['right'].set_position(('axes', -0.1))
    ax3[-1].set_frame_on(True)
    ax3[-1].patch.set_visible(False)
    c1=caf[0].mean(axis=1)
    c2=caf[1].mean(axis=1)
    mc1=1-c1.mean()
    mc2=1-c2.mean()
    c1=((c1+mc1)**3.5-mc1 - 0.03 ) * 0.85/2
    c2=((c2+mc2)**3.5-mc2 - 0.03 ) * 0.85/2
    c1[c1<0]=0
    c2[c2<0]=0
    S1=pd.Series(c1[idx],index=pd.DatetimeIndex(dates.num2date(tday[idx])))
    S2=pd.Series(c2[idx],index=pd.DatetimeIndex(dates.num2date(tday[idx])))
    obs=obs.tz_localize('UTC')
    print ((S1-1.5*obs).dropna()**2).mean()**(1/2.)
    print ((S2-1.5*obs).dropna()**2).mean()**(1/2.)
    ax3[-1].tick_params(axis='y', colors='g')
    ax3[0].tick_params(axis='y', colors='g')
    a6 = ax3[0].plot(tday3[idx3][::-1],1/qi[idx3][::-1],'g',linestyle='-.',lw=3.5,label='Moisture (right)')
    a5 = ax3[0].plot(tday3[idx3][::-1],1/ki[idx3][::-1],'g',linestyle=':',lw=3.5,label='Instability (right)')
    a11 = ax2.plot(tday31[idx31],obs.values[idx31],'#5e3c99',linestyle='-',lw=3.5,label='Convective Pixel (obs)')
    
    a4 = ax3[-1].plot(tday2[idx2],tr[idx2],'#e66101',linestyle='-',lw=3.5,label='Trigger function')
    #a0 = ax2.plot(tday1[idx1],obs.values[idx1]/4.,'g--',lw=3.5,label='observation')
    a2 = ax2.plot(tday[idx],c1[idx],'k-',\
            lw=3.5,label='deep (coastal)')
    a3 = ax2.plot(tday[idx],c2[idx],'k--',\
            lw=3.5,label='deep (original)')
    ax2.set_ylabel('Cloud area fraction []')
    lns = a2+a3+a11+a4+a5+a6+[l1,l2]
    labs = [l.get_label() for l in lns]
    ax3[-1].set_ylabel('Coastal effects []',color='#e66101', labelpad = -100)
    ax3[-1].tick_params(axis='y', colors='#e66101')
    ax3[0].set_ylabel('Moisture/Instability [ ]',color='g')
    ax3[0].tick_params(axis='y',color='g')
    ax2.xaxis.set_major_formatter(hfmt)
    ax2.grid()
    ax2.set_xlim(tday[0],tday[idx][-1])
    #ax3[0].set_ylim(0,3)
    #ax2.set_ylim(0.06,0.14)
    #ax3[-1].set_ylim(-1.8,1.8)
    font = {'family' : 'sans','weight' : 'normal','size'   : 24}
    import matplotlib
    matplotlib.rc('font', **font)
    #fig.autofmt_xdate()
    plt.legend(lns,labs,bbox_to_anchor=(0,-.02,1,-0.102),mode='expand',ncol=4, loc=2, borderaxespad=0.,prop={'size':22})
    #plt.legend(lns, labs, loc=2)

    plt.show()
        



################################################################################
################################################################################

def CD_sensitivity(**kwargs):
    """
    Method creates the CD sensitivity plot
    """
    from Master.configdir import Config
    from matplotlib import pyplot as plt

    datadir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    df_noint = os.path.join(datadir,'CD_sensitivity_0.00_1.00_0.000000.pkl')
    C = Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'constants.config'))
    df_int = os.path.join(datadir,"CD_sensitivity_%.2f_%.2f_%.6f.pkl"\
            %(C.dtmax,C.phase,C.add_c))
    
    data = pd.concat((pd.read_pickle(df_noint),pd.read_pickle(df_int)),axis=1,\
            keys=('noncoastal','coastal'))
    fig = plt.figure()

    value='mean'
    dd = 1
    u = 0.5
    kw = dict(c=((0,u),(-u,u)),d=((0,u),(-u,u)),s=((0,u),(-u,u)))
    kw2 = dict(c=((0,0.5),(-0.3,.3)),d=((0,.2),(-.1,.1)),s=((0,.7),(-.5,.5)))
    kw = kw2
    cl_dic = dict(c='congestus',d='deep',s='statiform')
    clouds = ('c','d','s')
    
    CC = np.linspace(0,2,100)
    DD = np.linspace(0,2,100)
    if len(data['noncoastal'].dropna()) != len(data['coastal'].dropna()):
        from scipy.interpolate import griddata
        grid_X,grid_Y=np.meshgrid(CC,DD)
        tmp = dict(C=data['noncoastal'].C.dropna().values,\
                D=data['noncoastal'].D.dropna().values)

        for ii in data['coastal'].keys():
            if ii not in ('C','D'):
                D = data['coastal']
                A=griddata(D[['C','D']].dropna().values,D[ii].dropna().values,\
                        (grid_X,grid_Y),method='linear')
                if ii.split('_')[-1] == 's':
                    tmp[ii] = A[::-1,::-1].ravel()
                elif ii.split('_')[-1] == 'c':
                    tmp[ii]=np.rot90(A)[::-1,:].ravel()
                else:
                    tmp[ii]=np.rot90(A[:,::-1]).ravel()
        
        d1 = data['noncoastal']
        d2 = pd.DataFrame(tmp)
        data = pd.concat([d1,d2],axis=1,keys=('noncoastal','coastal'))


    for  nn,cltype in enumerate(clouds):
        key="%s_%s"%(value,cltype)
        left = data['noncoastal'][key].dropna().values
        right = data['coastal'][key].dropna().values
        u = kw2[cltype][0][-1]
        left[left>u]=u
        right[right>u]=u
        for D,name in ((left,'left'),(right,'right'),(right-left,'diff')):
            ax = fig.add_subplot(3,3,dd)
            if name == 'diff':
                cmap='bwr_r'
                num = 1
                title = 'difference'
            else:
                cmap='Blues'
                num =0
                if name == 'left':
                    title='%s clouds (original)'%(cl_dic[cltype])
                else:
                    title='%s clouds (coastal)'%(cl_dic[cltype])
            if cltype != 'd':
                D=D.reshape(100,100).T
            else:
                D=D.reshape(100,100).T[:,::-1]
            #ax.imshow(D.T,cmap=cmap,extent=(0,2,0,2),clim=kw[cltype][num],interpolation=None)
            l,u = kw[cltype][num][0],kw[cltype][num][-1]
            D[D>u]=u
            V = np.linspace(l,u,20).round(2)
            #print D.shape
            if name == 'diff' and cltype=='d':
                print(D)
            contour = ax.contourf(CC,2-DD,D,V,cmap=cmap)
            fig.colorbar(contour,orientation='horizontal',pad=0.1,aspect=40)
            if nn == len(clouds)-1:
                ax.set_xlabel('Instability')
                """
                if name is 'right':
                    cbar_ax = fig.add_axes([0.1, 0.05, 0.52, 0.01])
                    fig.colorbar(contour,orientation='horizontal',cax=cbar_ax)
                elif name is 'diff':
                    cbar_ax = fig.add_axes([0.66, 0.05, 0.25, 0.01])
                    fig.colorbar(contour,orientation='horizontal',cax=cbar_ax)
                """

            if nn != len(clouds) - 1:
                ax.set_xticks([])
            if name == 'left':
                ax.set_ylabel('Moisture')
            else :
                ax.set_yticks([])

            ax.set_title(title)
            dd += 1
    fig.subplots_adjust(wspace=0.1)
    plt.show()

##############################################################################
##############################################################################
def correlate(**kwargs):
    """
    Method that calculates the correlation of a model run with precip data
    """
    
    box=check(kwargs,'boxes',['coast_03'])[0]
    slmf=check(kwargs,'slmf',os.path.join(os.environ['HOME'],'Data','Cmorph_slm.nc'))
    stype=check(kwargs,'stype','land')
    lag = int(check(kwargs,'lag',-3))
    datadir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    add = dict(true='sobol',false='nosobol')[str(kwargs['sobol']).lower()]
    datafile = os.path.join(datadir,'coastal_sensitivity_%s.nc'%add)
    outfile = datafile.replace('.nc','.pkl')
    corrfile = outfile.replace('sensitivity','correlation_%s'%stype)
    if not os.path.isfile(corrfile):
    #Test if the output file exsists
        with nc(datafile,'a') as ncf:
    #Test if the precip data is already in datafile
            if not 'precip_%s'%stype in ncf.variables.keys():
                __readPrecip(box,ncf,slmf=slmf,stype=stype)
            output,params = __getTimeS(ncf,stype)
        output.to_pickle(corrfile)
        params.to_pickle(outfile.replace('sensitivity','params'))

    else:
        output = pd.read_pickle(corrfile)
        params = pd.read_pickle(outfile.replace('sensitivity','params'))
    df,precip = __corr(output,params,lag)
    df.to_pickle(outfile.replace('sensitivity','final_correlation'))
    precip.to_pickle(outfile.replace('sensitivity','precip_correlation'))
    Plot.map(df,precip,**kwargs)

def scatter(**kwargs):
    """
    This method compares the results of the simulation for deep convective
    clouds with precipitation measurments
    """
    from matplotlib import pyplot as plt
    import matplotlib
    datadir = os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
    precipdir = os.path.join(os.environ['HOME'],'PhD','Data','LargeScale')

    precipfile = os.path.join(precipdir,'Satrain_fm.pkl')
    modelfile = os.path.join(datadir,'Modelrain.pkl')
    if not os.path.isfile(precipfile):
        out = []
        for i in xrange(1,12):
            fname = os.path.join(precipdir,'Database_fm_coast_%02i.pkl'%(i))
            if os.path.isfile(fname):
                out.append(pd.read_pickle(fname)[['precip','lsp']])
            else:
                sys.stderr.write('Warning %s not found\n'%fname)

        boxes = ['COAST_%02i'%i for i in xrange(1,12)]

        out = pd.concat(out,keys=boxes,axis=1)
        out.to_pickle(precipfile)

    precipdata = pd.read_pickle(precipfile)
    modeldata = pd.read_pickle(modelfile)

    out = None
    boxes = ['COAST_%02i'%i for i in xrange(1,12)]
    #boxes = ['COAST_06']
    columns = ['lsp','precip','deep']
    for box in boxes:
        df = modeldata[box]['deep']
        idx = pd.DatetimeIndex(((df.index.asi8/(60**2*1e9)).round()*(60**2*1e9)).astype(np.int64))
        model = df.groupby(idx).mean()
        rain = precipdata[box]
        if type(out) == type(None):
            out = pd.concat([rain,model],axis=1).dropna()
        else:
            out = pd.concat([out,pd.concat([rain,model],axis=1).dropna()],axis=0)


    
    T =dict(lsp='coastal',precip='total')
    fig = plt.figure()
    matplotlib.rcParams.update({'font.size': 18})
    for ii,nn in enumerate(('precip','lsp')):
        idx = (out[nn] > 1) & (out['deep'] > 0.02)
        Out = out.loc[idx]
        data = Out[['deep',nn]].values
        ax = fig.add_subplot(1,2,ii+1)
        ax.scatter(data[:,0],data[:,1])
        ax.set_ylabel('%s precipitation [mm/3h]'%T[nn])
        ax.set_xlabel('cloud area fraction')
        ax.set_xlim(0.02,0.35)
        ax.set_ylim(1,25)
    plt.show()







##############################################################################
##############################################################################

def kw(**kwargs):
    """
    Method that returns the default parameter keywordarguments
    """
    F={'correlate':correlate,'CD_sensitivity':CD_sensitivity,'scatter':scatter,
            'diurnal_cycle':diurnal_cycle,'large_scale':large_scale,\
                    'variability':variability,'stochasticity':stochasticity}
    KW=dict(
        boxes=['coast_08'],
        func=correlate,
        cloudtypes='deep',
        type='corr',
        zip='add_c',
        sobol=False)
    for key,value in kwargs.iteritems():
        if key == 'func' and isinstance(value,str):
            value = F[value]
        KW[key]=value
    return KW




##############################################################################
##############################################################################

if __name__ == "__main__":
    kwargs=kw()
    helpstring="""
    Module to run the cloudmodel in various ways, without copuling.
    The method is supposed to be a bit more flexible

    Usage:
        python %s --option=value
        options are given below

    Options:
        boxes   (list)  : The name of the boxes that are considered (default %s)
        func (function) : Which function should be considered (default %s)
        cloudtypes(list): What cloud types should be plotted (default %s)
        Type (str)      : Should be mean or correlation map plotted (defaul %s)
        Zip (list)      : Which two of the three parameters should be zipped default %s)
        sobol (bool)    : Should the correlation be calculated form a sobol sequence (default %s)

    """%(
            sys.argv[0],
            str(kwargs['boxes']),
            str(kwargs['func']),
            str(kwargs['cloudtypes']),
            str(kwargs['type']),
            str(kwargs['zip']),
            str(kwargs['sobol'])\
        )
    try:
        for arg in sys.argv[1:]:
            try:
                key,value = arg.strip('--').split('=')
            except ValueError:
                sys.exit(helpstring)
            if key.lower()=='help':
                sys.exit(helpstring)
            #Special treatment for plot and boxes because they are by deault
            #booleans and lists
            elif key.lower() in ['boxes']:
                kwargs[key]=value.split(',')
            elif key.lower() == 'sobol':
                kwargs[key] = value[0].lower() in ('t','y','1')
            #Everything else shold go here
            else :
                try: #The key is an integer set it
                    kwargs[key.lower()] = float(value)
                except ValueError:
                    try: #Do we have a list of integers seperated by kommas?
                        kwargs[key.lower()] = [float(v) for v in value.split(',')]
                    except ValueError: #If not, probably a string
                        kwargs[key.lower()] = value

    except IndexError:
        pass
    kwargs = kw(**kwargs)
    main = kwargs['func']
    del kwargs['func']
    main(**kwargs)

    #The interaction potential matrix
    #And the ocean model



