#!/usr/bin/env python2.7

"""
Module to make experitmens with the real cloudmodel
"""
import sys,os,datetime
from coarsgraining import coarsgraining
from model import SMCM
import pandas as pd
from configdir import Config
from Driver import check
import geocoder
import numpy as np

def integrate(cl,**kwargs):
    """
    This method should integrate the cloud model only once

    Variables:
        cl (coarsegraining object) = the coarsgraining object
    """

    #Time-vector
    t = [] #The times 
    tt = int(0)
    while tt <= int(cl.tend * 60):
        t.append(float(tt)/60.)
        tt += int(cl.dt * 60)
    t = np.array(t)/24
    #What should be printed
    Print=check(kwargs,'Print',None)
    printkw = check(kwargs,'printkw',{})
    #Should the cloudmap be saved to a file
    fname=check(kwargs,'fname',None)
    
    cl._update(0*24,m=None,n=None)
    ary_con=np.empty([len(t),cl.m,cl.m])
    ary_deep=np.empty_like(ary_con)
    ary_strat=np.empty_like(ary_con)
    Boxconfig = Config(os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'boxes.txt'))

    try:
        RegionName='in %s' %GetNames(Boxconfig[cl.obs.lower()])
    except (KeyError,AttributeError):
        RegionName=''

    if isinstance(cl.animate,str):
        import matplotlib
        matplotlib.use('TkAgg')
        from matplotlib import pyplot as plt
        ani=plt.figure()
        if isinstance(cl.form,tuple) and not cl.plotmask:
            from mpl_toolkits.basemap import Basemap
            l_lon=cl.form[0]
            r_lon=cl.form[1]
            l_lat=cl.form[2]
            u_lat=cl.form[3]
            m=Basemap(llcrnrlon=l_lon,llcrnrlat=l_lat,urcrnrlon=r_lon,urcrnrlat=u_lat,\
                resolution='f',area_thresh=10.,\
                projection='cyl',lon_0=0,suppress_ticks=True)
            #m.drawcoastlines(ax=ax,color='b')

            ax=ani.add_axes([0.02,0.03,0.87,0.96])
        else:
            plt.axis([0,cl.m,0,cl.m])
            plt.text(cl.m/4.,-0.25,'Ocean',fontsize='14')
            plt.text(3*cl.m/4.,-0.25,'Land',fontsize='14')
            plt.plot([cl.m/2.,cl.m/2.],[0.,cl.m],'r-',lw=1)
        plt.xticks([])
        plt.yticks([])
        plt.ion()
        plt.show()
    pcolor = None
    #The time list
    tday, thc = [], []
    #A little helper for getting the right appendix of numbers
    App = {1:'st',2:'nd',3:'rd'}
    y = cl.start
    cbar=None
    for ii,tt in enumerate(t):
        #Call the coarsegraining procedure
        cl.birthdeath(cl.dt,24*tt)
        #And integrate the convection-scheme
        cl._update(tt*24,m=None,n=None)
        tdiff=cl.RD.thc(cl.dtmax,24*tt)
        if type(Print) != type(None):
            sys.stdout.flush()
            sys.stdout.write(Print(cl,tt,RegionName,tdiff,**printkw))
            sys.stdout.flush()
        
        #Output the results
        tday.append(cl.start + datetime.timedelta(hours=24*tt))
        thc.append(tdiff)
        cafd=(cl.Ndcg)/cl.q**2
        cafc=(cl.Nccg)/cl.q**2
        cafs=(cl.Nscg)/cl.q**2
        
        ary_strat[ii]=cafs
        ary_deep[ii]=cafd
        ary_con[ii]=cafc
        #Create an animation
        if isinstance(cl.animate,str):
            
            Ani = dict(deep=cafd,strat=cafs,con=cafd,mask=cl.c_f(cl.lsm*tdiff))[cl.animate]
            vmin,vmax=dict(deep=(0,0.25),strat=(0,1),con=(0,0.75),mask=(0,1.5))[cl.animate]
            if pcolor:
                pcolor.remove()
            elif isinstance(cl.form,tuple) and not cl.plotmask:
                m.drawcoastlines(ax=ax,color='b')
            if isinstance(cl.form,tuple) and not cl.plotmask:
                pcolor=m.imshow(Ani[::-1,:],cmap=matplotlib.cm.gray,
                        vmin=vmin,vmax=vmax,interpolation='bicubic',ax=ax)
                st = y.strftime('%b. %Y %H:%M LT')
                day = y.day
                try:
                    app=App[day]
                except KeyError:
                    app='th'
                st='%2i%s %s' %(day,app,st)


                fname='/home/wilfred/test'
                plt.title('d: %02.2f c: %02.2f - %2.3f %% at %s %s'\
                        %(cl.D[0,0],cl.C[0,0],Ani.mean()*100))#st,#RegionName))

            else:
                if not cl.obs:
                    day = (y - datetime.datetime(\
                        cl.start.year,cl.start.month,cl.start.day,0,0)).days
                    ext=''
                else:
                    day=y.day
                    ext=y.strftime(' %b %Y')
                try:
                    app=App[day]
                except KeyError:
                    app='th'
                day = "%3i%s%s"%(day,app,ext)
                hour = y.hour
                Min = y.minute
                pcolor=plt.pcolor(Ani,cmap=matplotlib.cm.gray,vmin=vmin,vmax=vmax)
                if cl.animate == 'mask':
                    plt.title("Local time %02i:%02i "%(hour,Min),size=28)
                    #plt.vline(0.5,0,1,r,lw=3)
                    if type(cbar) == type(None):
                        cbar=plt.colorbar(pcolor)
                        cbar.ax.tick_params(labelsize=28)
                    plt.plot(Ani.shape[0]/2.*np.ones(100),np.linspace(0,Ani.shape[1],100),'r-',lw=3)
                    plt.text(Ani.shape[0]/4,-0.5,'Land',size=28,verticalalignment='bottom', horizontalalignment='center')
                    plt.text(3*Ani.shape[0]/4,-0.5,'Ocean',size=28,verticalalignment='bottom', horizontalalignment='center')
                    fname='/home/wilfred/test'
                else:
                    fname='/home/wilfred/test'
                    plt.title('d: %02.2f $\\omega_c$: %02.2f %2.3f %% at %s %02i:%02i'\
                            %(cl.D[0,0],cl.C[0,0],Ani.mean()*100,day,hour,Min))#,RegionName))
            if fname:
                fout=os.path.join(fname,'fig_%04i_%i_%s.png'\
                        %(ii,int(cl.interact),cl.type))
                ani.savefig(fout,format='png',dpi=72,bbox_inches='tight')

            
            ani.canvas.draw()
        y += datetime.timedelta(seconds=cl.dt*60**2)
            
    if isinstance(cl.animate,str):
        plt.close(),ani.clf()
    if not cl.plotmask:
        return ary_con,ary_deep,ary_strat,np.array(thc),pd.DatetimeIndex(tday)
    else:
       return ary_deep,pd.DatetimeIndex(tday)

##############################################################################
##############################################################################
def init(**kwargs):
    """
    Method that creates and returns the coarsgraining object
    """
    #Scaling values cape and moisture
    M0 = check(kwargs,'M0', 20)
    C0 = check(kwargs,'C0', 1500)
    conf = check(kwargs,'conf', 'constants.config')
    multicol = check(kwargs, 'multicol',True)
    J0 = check(kwargs,'J0', np.array([[0.75, 0.25, 0.],
                                     [0.4,0.85,0.],
                                     [0,0.2,1]]))
    #Initialize the cloud model
    seed=check(kwargs, 'seed', None)
    return coarsgraining(conf, C0 ,M0, multicol=multicol, J0=J0, seed=seed)

##############################################################################
##############################################################################
def __getdates(**kwargs):
    """
    Helper method to distribute dates

    """
    start=kwargs['start']
    end=kwargs['end']
    del kwargs['start'],kwargs['end']
    Cfg = Config('constants.config')
    if start and end:
        #Take the date information from the arguments
        dates = pd.date_range(start,end,freq='H')
        nn = np.unique(dates.year).shape[0]
        AA = np.array_split(dates,nn)
        add = 3
        Cfg.start = (dates[0]).strftime('%Y-%m-%d_%H:%M')

        delta = dates[-1] - dates[0]
        hours = delta.days*24 + delta.seconds/60**2 + add
        Cfg.tend=hours
    else:
        #Take the information from the config file
        add = 0
    end = (datetime.datetime.strptime(Cfg.start.strip(),'%Y-%m-%d_%H:%M')+\
            datetime.timedelta(hours=Cfg.tend-add))\
            .strftime('%Y.%m.%d_%H:%M')

    #Ok, create the array containing all dates
    dates = pd.date_range(Cfg.start.replace('_','-'),end.replace('_','-'),
            freq='D')
    return Cfg,dates

##############################################################################
##############################################################################


def run_once(**kwargs):
    """
    This method only runs the cloudmodel once and saves the output to a file

    """
    def Print(*arg):
        RegionName=arg[-2]

        thc=arg[-1]
        time=arg[0].start+datetime.timedelta(hours=24*arg[1])
        time=time.strftime('%d %b %Y %H:%M')
        return 'Iterating %.2f (tdiff) %.2f (phase) + %.2f (const) on %s %s\n'\
                %(thc,arg[0].phase,arg[0].add_c,time,RegionName)
    from mpi4py import MPI
    #Process based mpi
    root = 0
    if type(kwargs['rank']) == type(None):
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        name = comm.name
        scatter = True
    else:
        rank = int(kwargs['rank'])
        size = int(kwargs['size'])
        scatter = False
        if size <= rank:
            size = rank+1



    Cfg,dates  = __getdates(**kwargs)
    #Split the days for the MPI
    samples = np.array_split(dates,size)
    if scatter :
        if rank == root:
            sendbuf=samples
        else:
            sendbuf = []
    #Scatter the buffer to all processes, if empty it is dropped
        X = comm.scatter(sendbuf,root)
    else:
        X = samples[rank]

    if len(X) > 0:
        Cfg.start = (X[0]).strftime('%Y-%m-%d_%H:%M')
        end = (X[-1]).strftime('%Y-%m-%d_%H:%M')
        delta = X[-1] - X[0]
        hours = max(delta.days*24 + delta.seconds/60.**2,24)
        Cfg.tend=hours
        for box in kwargs['boxes']:
            if box:
                Cfg.obs = box
            if type(Cfg.obs) != type(None):
                Str=os.path.join(os.environ['HOME'],'PhD','Data','SMCM',
                        'Modelrain_%s_%s-%s.pkl'%(Cfg.obs.upper(),
                            Cfg.start.replace('-','.'),end.replace('-','.')))
            #"""
            else:
                Str=os.path.join(os.environ['HOME'],'PhD','Data','SMCM',
                        'Modelrain_None_%s-%s.pkl'%(
                            Cfg.start.replace('-','.'),end.replace('-','.')))
            if os.path.isfile(Str):
                os.remove(Str)
            
            ######################################################################3y
            ######################################################################3y
            data=None

            cl = init(conf=Cfg,**kwargs)
            con,deep,strat,thc,time = integrate(cl,Print=Print)
            startidx = time[0]+ pd.Timedelta(hours=3)
            idx = np.where(time >= startidx)[0]
            deep = deep.mean(axis=(1,2))[idx]
            con = con.mean(axis=(1,2))[idx]
            strat = strat.mean(axis=(1,2))[idx]
            time = time[idx]
            thc=thc[idx]
            data=pd.DataFrame({'con':con,'deep':deep,'strat':strat,
                'thc':thc},index=time)
            if not os.path.isfile(Str):
                data.to_pickle(Str)
            else:
                old=pd.read_pickle(Str)
                os.remove(Str)
                (old.append(data)).to_pickle(Str)

            del deep,con,strat,idx,cl,data
        #"""
        out=(X[0],X[-1])
    else:
        out=(None,None)
    if scatter:
        recbuf = comm.gather(out,root)
        if rank == root:
            sendbuf ="done"
            recvbuf = comm.bcast(sendbuf,root)
    else:
        return
    if type(recbuf) != type(None):
        if len(samples) <= 1:
            #Don't worry only on process here
            return
        return
        first = min(i for i,j in recbuf if type(i) is not type(None))
        last =  max(j for i,j in recbuf if type(j) is not type(None))

        sys.stdout.write('Gathering all data .... ')

        for box in kwargs['boxes']:
            first = first.strftime('%Y.%m.%d_%H:%M')
            last = last.strftime('%Y.%m.%d_%H:%M')
            Write = os.path.join(os.environ['HOME'],'PhD','Data','SMCM',
                    'Modelrain_%s_%s-%s.pkl'%(box.upper(),first,last))
            if os.path.isfile(Write):
                os.remove(Write)
            data = None
            for X in samples:
                if len(X) > 0:
                    start = (X[0]).strftime('%Y.%m.%d_%H:%M')
                    end = (X[-1]).strftime('%Y.%m.%d_%H:%M')

                    Read = os.path.join(os.environ['HOME'],'PhD','Data','SMCM',
                            'Modelrain_%s_%s-%s.pkl'%(box.upper(),start,end))

                    if not os.path.isfile(Read):
                        sys.stderr.write('\n File %s does not exsist, skipping\n'%Read)
                    else:
                        if type(data) == type(None):
                            data = pd.read_pickle(Read)
                        else:
                            data = pd.concat((data,pd.read_pickle(Read)),axis=0)

            data.index =pd.DatetimeIndex((((data.index.asi8+30*1e9)/\
                    (1e9*60)).round()*(1e9*60)).astype(np.int64))
            data.to_pickle(Write)
            sys.stdout.write(' ok\n')

##############################################################################
##############################################################################
def diurnal_cylce(**kwargs):

    """
    Methd to calculate the dirunal cycle and save to 

    """
    #Process based mpi
    from mpi4py import MPI
    from netCDF4 import Dataset as nc,date2num
    from netCDF4 import MFDataset as mnc

    #Process based mpi
    root = 0
    if type(kwargs['rank']) == type(None):
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        name = comm.name
        scatter = True
    else:
        rank = int(kwargs['rank'])
        size = int(kwargs['size'])
        scatter = False
        if size <= rank:
            size = rank+1



    def Print(*arg,**kw):
        RegionName=arg[-2]

        thc=arg[-1]
        time=arg[0].start+datetime.timedelta(hours=24*arg[1])
        time=time.strftime('%d %b %Y %H:%M')
        return 'Iterating %.2f (tdiff) %.2f (phase) + %.2f (const) on %s %s\n'\
                %(thc,arg[0].phase,arg[0].add_c,time,RegionName)

    Cfg,dates = __getdates(**kwargs)
    #Split the days for the MPI
    boxes=kwargs['boxes']
    del kwargs['start'],kwargs['end'],kwargs['boxes']
    samples = np.array_split(dates,size)
    if scatter:
        if rank == root:
            sendbuf=samples
        else:
            sendbuf = []
        #Scatter the buffer to all processes, if empty it is dropped
        X = comm.scatter(sendbuf,root)
    else:
        X = samples[rank]
    try:
        Cfg.seed=kwargs['seed']
        del kwargs['seed']
    except KeyError:
        pass
    if len(X) > 0:
        Cfg.start = (X[0]).strftime('%Y-%m-%d_%H:%M')
        end = (X[-1]).strftime('%Y-%m-%d_%H:%M')
        delta = X[-1] - X[0]
        hours = max(delta.days*24 + delta.seconds/60.**2,24)
        Cfg.tend=hours
        datadir=os.path.join(os.environ['HOME'],'PhD','Data','SMCM')
        if not os.path.isdir(datadir):
            datadir=os.path.join(os.environ['HOME'],'Data','SMCM')

        for box in boxes:
            if box:
                Cfg.obs = box
            Str=os.path.join(datadir,
                    'Diurnalcycle_%s_%s-%s_%s_tmp.nc'%(Cfg.obs.upper(),
                        Cfg.start.replace(':','_'),end.replace(':','_'),str(Cfg.seed)))
            #"""
            if os.path.isfile(Str):
                os.remove(Str)

            ######################################################################3y
            ######################################################################3y
            data=None
            
            cl = init(conf=Cfg,**kwargs)
            con,deep,strat,thc,time = integrate(cl,Print=Print)
            startidx = time[0]+ pd.Timedelta(hours=3)
            idx = np.where(time >= startidx)[0]
            deep = deep[idx]
            con = con[idx]
            strat = strat[idx]
            time = time[idx]
            thc=thc[idx]
            geoinfo = os.path.join(os.path.dirname(os.path.abspath(__file__)),'boxes.txt')
            lat_c,lon_c = Config(geoinfo)[Cfg.obs.lower()]
            lon = np.linspace(lon_c-1.5,lon_c+1.5,deep.shape[-1])
            lat = np.linspace(lat_c-1.5,lat_c+1.5,deep.shape[1])
            attrs =dict(\
                    time=dict(\
                        units='Minutes since 1998-01-01 00:00:00',
                        axis='T',
                        long_name='Time'),
                    lon=dict(\
                        units='degrees_east',
                        axis='X',
                        long_name='lngitude'),
                    lat=dict(\
                        units='degrees_north',
                        axis='Y',
                        long_name='latitude'))
            with nc(Str,'w',format="NETCDF4_CLASSIC") as h5:
                for i,j in (('lon',lon),('lat',lat),('time',None)):
                    if i == 'time':
                        h5.createDimension(i,None)
                        h5.createVariable(i,'f',(i,))
                        h5.variables['time'][:]=date2num(time.to_pydatetime(),\
                                attrs['time']['units']).round(0)
                    else:
                        h5.createDimension(i,len(j))
                        h5.createVariable(i,'f',(i,))
                        h5.variables[i][:] = j
                    for key,value in attrs[i].items():
                        setattr(h5.variables[i],key,value)

                for i,j in (('deep',deep),('con',con),('strat',strat)):
                    h5.createVariable(i,'f',('time','lat','lon'))
                    h5.variables[i][:]=j
                    h5.variables[i].units='[]'
                    h5.variables[i].long_name='CAF of %s clouds'%i
                    h5.variables[i].grid='latlon'
            
            del deep,con,strat,idx,cl,data
        #"""
        out=(X[0],X[-1])
    else:
        out=(None,None)
    if scatter:
        recbuf = comm.gather(out,root)
        if comm.rank == root:
            sendbuf ="done"
            recvbuf = comm.bcast(sendbuf,root)
    else:
        return
        
    if type(recbuf) != type(None):
        
        first = min(i for i,j in recbuf if type(i) is not type(None))
        last =  max(j for i,j in recbuf if type(j) is not type(None))

        sys.stdout.write('Gathering all data .... ')
        if Cfg.interact:
            interact='interact'
        else:
            interact='nointeract'
        for box in boxes:
            first = first.strftime('%Y.%m.%d_%H_%M')
            last = last.strftime('%Y.%m.%d_%H_%M')
            Write = os.path.join(datadir,
                    'Diurnalcycle_%s_%s-%s_%s.nc'%(Cfg.obs.upper(),first,last,str(Cfg.seed)))
            if os.path.isfile(Write):
                mode='a'
            else:
                mode = 'w'

            with nc(Write,mode) as h5:
                for i,j in (('lon',lon),('lat',lat),('time',None)):
                    if i == 'time':
                        ll = None
                    else:
                        ll = len(j)
                    try:
                        h5.createDimension(i,ll)
                    except (RuntimeError,ValueError,IOError,OSError):
                        pass
                    try:
                        h5.createVariable(i,'f',(i,))
                    except (RuntimeError,ValueError,IOError,OSError):
                        pass
                try:
                    g = h5.createGroup(interact)
                except (RuntimeError,ValueError,IOError,OSError):
                    g = h5.groups[interact]

                try:
                    w = g.createGroup('%1.2f'%Cfg.J00)
                except (RuntimeError,ValueError,IOError,OSError):
                    w = g.groups('%1.2f'%Cfg.J00)

                for i in ('deep','con','strat'):
                    try:
                        w.createVariable(i,'f',('time','lat','lon'))
                    except (RuntimeError,ValueError,IOError,OSError):
                        pass
                Files=[]
                for s, e in recbuf:
                    start = s.strftime('%Y-%m-%d_%H_%M')
                    end = e.strftime('%Y-%m-%d_%H_%M')
                    file = os.path.join(datadir,
                    'Diurnalcycle_%s_%s-%s_%s_tmp.nc'%(Cfg.obs.upper(),start,end,str(Cfg.seed)))
                    
                    if os.path.isfile(file):
                        Files.append(file)

                with mnc(Files,'r') as f:
                    for v in ('deep','con','strat','time','lat','lon'):
                        if v in ('deep','con','strat'):
                            var = w
                        else:
                            var = h5
                        var.variables[v][:]=f.variables[v][:]
                        for attr in ('units','axis','long_name','grid'):
                            try:
                                value=getattr(f.variables[v],attr)
                                setattr(var.variables[v],attr,value)
                            except AttributeError:
                                pass

            sys.stdout.write(' ok\n')


##############################################################################
##############################################################################

def coastal_sensitivity(**kwargs):

    """
    This method creates a variance based sensitivity analysis
    based on a various values for the values of dtmax, pahse and add_c
    the model should be driven by observational data from:
        2004-01-14 - 2004-01-16
        in the Bight of Panama
    
    """
    from mpi4py import MPI
    from netCDF4 import date2num,Dataset as nc
    
    #Process based mpi
    root = 0
    if type(kwargs['rank']) == type(None):
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        name = comm.name
        scatter = True
    else:
        rank = int(kwargs['rank'])
        size = int(kwargs['size'])
        scatter = False
        if size <= rank:
            size = rank+1

    def Print(*arg,**kw):
        RegionName=arg[-2]

        thc=arg[-1]
        time=arg[0].start+datetime.timedelta(hours=24*arg[1])
        time=time.strftime('%d. %b %Y %H:%M')
        rank=kw['rank']
        perc=kw['perc']
        return 'Iterating C: %.2f D: %.2f'%(arg[0].C[0,0],arg[0].D[0,0])\
                +' with tdiff: %.2f with c: %.2f at %s (%2i |%3.2f)\n'\
                %(thc*arg[0].mul,arg[0].add_c,time,rank,perc)
    #Read the configfile
    Cfg = Config('constants.config')
    Cfg.dt=10./60.
    #Now setup the model
    start=kwargs['start']
    end=kwargs['end']
    
    #Cfg.start='2000-09-08_18:00'
    Cfg.start='2000-08-31_18:00'
    #Six hours are considered as spin off
    #Calculate how long the integration should be
    
    NN=int(1/Cfg.dt*6)
    start=datetime.datetime.strptime(Cfg.start,'%Y-%m-%d_%H:%M')
    #end=datetime.datetime(2000,9,10,0,0)
    end=datetime.datetime(2000,10,1,0,0) - datetime.timedelta(hours=start.hour)
    Cfg.tend = (end-start).total_seconds()/60.**2
    #Set the observation box to Borneo
    Cfg.interact = True
    Cfg.obs = 'coast_03'
    lat,lon=Config('boxes.txt')[Cfg.obs]
    #set the filename for the output
    fname = 'coastal_sensitivity_%s.nc'%(\
            dict(true='sobol',false='nosobol')[str(kwargs['sobol']).lower()])
    
    fname=os.path.join(os.environ['HOME'],'PhD','Data','SMCM',fname)
    attrs=dict(
            lat=dict(long_name='latitude',units='degrees_east',axis='X'),
            lon=dict(long_name='longitude',units='degrees_nort',axis='Y'),
            params=dict(long_name='sobol parameter',units='dtmax add_c'),
            time=dict(long_name='time',units='Minutes since 1998-01-01 00:00:00',axis='T'))
    #Check for preconditioned dtmax
    #if so check if it is 0
    if 'dtmax' in kwargs.keys():
        if kwargs['dtmax'] == 0:
            sys.stdout.write("Running only once\n")
            Cfg.dtmax = 0
            Cfg.mul=0
            Cfg.add_c = 0
            Cfg.interact = False
            del kwargs['dtmax']
            #We need acually only modle run
            if rank == root:
                printkw=dict(rank=rank,perc=100)
                cl = init(conf=Cfg,**kwargs)
                con,deep,strat,thc,time_m = integrate(cl,Print=Print,printkw=printkw)
                if os.path.isfile(fname):
                    mode='a'
                else:
                    mode='w'
                with nc(fname,mode) as h5 :
                    try:
                        h5.createGroup('nointeract')
                    except (ValueError,RuntimeError):
                        pass
                    f = h5['nointeract']
                    lon=np.linspace(lon-1.5,lon+1.5,cl.m)
                    lat=np.linspace(lat-1.5,lat+1.5,cl.m)
                    for i,j in (('lat',lat),('lon',lon),('time',time_m[NN:])):
                        if i == 'time':
                            s = None
                            Values=date2num(j.to_pydatetime(),attrs[i]['units'])
                        else:
                            s=len(j)
                            Values=j
                        try:
                            h5.createDimension(i,s)
                        except (ValueError,RuntimeError):
                            pass

                        try:
                            h5.createVariable(i,'f',(i,))
                        except (ValueError,RuntimeError):
                            pass
                        for key,value in attrs[i].items():
                            setattr(h5[i],key,value)
                        h5[i][:]=Values
                    for i,j in (('con',con[NN:]),('deep',deep[NN:]),('strat',strat[NN:])):
                        try:
                            f.createVariable(i,'f',('time','lat','lon'))
                        except (ValueError,RuntimeError):
                            pass
                        
                        f[i][:] = j
                        f[i].units=' '
                        f[i].standard_name='CAF of %s. clouds'%i
                return



    #Create the sobol sequence
    problem={
            'num_vars':2,
            'names':['mul','add_c'],
            'bounds':[[0.,2],[0,2]]}
    bounds=[]
    num_b = 25
    if num_b % 2 == 0:
        n2 = int(num_b/2)
    else:
        n2 = (num_b+1)/2
    bounds.append(np.r_[np.linspace(problem['bounds'][0][0],1,n2)[:-1]\
            ,np.linspace(1,problem['bounds'][0][1],n2)])
    bounds.append(np.linspace(problem['bounds'][1][0],\
            problem['bounds'][1][1],num_b))
    if kwargs['sobol']:
    #Create the random input data and split it into n samples
        from SALib.sample import saltelli
        Global_param_values = saltelli.sample(problem,500)
    else:
        Global_param_values = \
                np.array([(ii,jj) for ii in bounds[0] for jj in bounds[1]])
    samples = np.array_split(Global_param_values,size)
    if scatter:
        if rank == root:
            sendbuf=samples
        else:
            sendbuf = []
        #Scatter the buffer to all processes, if empty it is dropped
        X = comm.scatter(sendbuf,root)
    else:
        X = samples[rank]
    
    #Now instanciate the model setup
    try:
        cl = init(conf=Cfg,**kwargs)
        #Create tmporary file
        with nc(fname.replace('.nc','_%02i.nc'%rank),'w') as h5:
            lon=np.linspace(lon-1.5,lon+1.5,cl.m)
            lat=np.linspace(lat-1.5,lat+1.5,cl.m)
            h5.createDimension('lon',len(lon))
            h5.createDimension('lat',len(lat))
            h5.createDimension('time',None)
            h5.createVariable('lon','f',('lon',))
            h5.createVariable('time','f',('time',))
            h5.createVariable('lat','f',('lat',))
            h5.variables['lon'][:]=lon
            h5.variables['lat'][:]=lat
        #Count the lengths of the prevous arrays:
        tt = 0
        for nn in xrange(rank):
            tt += len(samples[nn])
        out = []
        with nc(fname.replace('.nc','_%02i.nc'%rank),'a') as h5:
            for ii,setup in enumerate(X):
                printkw=dict(rank=rank,perc=float(ii)/len(X)*100)
                multi,add_c = tuple(setup)
                cl.mul=multi
                cl.add_c=add_c
                con,deep,strat,thc,time_m = integrate(cl,Print=Print,printkw=printkw)
                h5.createGroup(str(tt+ii))
                f=h5.groups[str(tt+ii)]
                f.createVariable('con','f',('time','lat','lon'))
                f.createVariable('deep','f',('time','lat','lon'))
                f.createVariable('strat','f',('time','lat','lon'))
                f.variables['con'][:]=con[NN:]
                f.variables['deep'][:]=deep[NN:]
                f.variables['strat'][:]=strat[NN:]
                out.append(tt+ii)
                del f,con,deep,strat,thc
            h5.variables['time'][:]=date2num(time_m[NN:].to_pydatetime(),
                    'Minutes since 1998-01-01 00:00:00')
        #Send the output back to the root communicator
        out={rank:out}
    except IndexError:
        if os.path.isfile(fname.replace('.nc','_%02i.nc'%rank)):
            os.remove(fname.replace('.nc','_%02i.nc'%rank))
        out={rank:[]}
    if scatter:
        recbuf = comm.gather(out,root)
        if comm.rank == root:
            sendbuf ="done"
            recvbuf = comm.bcast(sendbuf,root)
    else:
        return
    if type(recbuf) != type(None):
        sys.stdout.write("Gathering all information into one file.....\n")
        if os.path.isfile(fname):
            mode = 'a'
        else:
            mode = 'w'
        with nc(fname,mode) as h5 :
            try:
                h5.createGroup('interact')
            except (ValueError,RuntimeError):
                pass
            
            for i,j in (('lat',lat),('lon',lon),('time',time_m[NN:]),
                    ('num',Global_param_values[0]),
                    ('params',Global_param_values[:,0])):
                if i == 'time':
                    s = None
                    Values=date2num(j.to_pydatetime(),attrs[i]['units'])
                elif i == 'params' or i == 'num':
                    s=len(j)
                    Values=None
                else:
                    s=len(j)
                    Values=j
                try:
                    h5.createDimension(i,s)
                except (ValueError,RuntimeError):
                    pass
                if type(Values) != type(None):
                    try:
                        h5.createVariable(i,'f',(i,))
                    except (ValueError,RuntimeError):
                        pass
                    for key,value in attrs[i].items():
                        setattr(h5[i],key,value)
                    h5[i][:]=Values
            try:
                h5.createVariable('params','f',('params','num'))
            except (ValueError,RuntimeError):
                pass
            h5['params'][:]=Global_param_values
            [setattr(h5['params'],k,v) for k,v in attrs['params'].items()]
            
            g = h5['interact']
        
            for ii in xrange(len(recbuf)):
                try:
                    with nc(fname.replace('.nc','_%02i.nc'%ii),'r') as sourcef:
                        for tt in recbuf[ii][ii]:
                            In = Global_param_values[tt]
                            sourced=sourcef[str(tt)]
                            sys.stdout.flush()
                            sys.stdout.write('%02i - %i ... ' %(ii,int(tt)))
                            sys.stdout.flush()
                            con = sourced['con'][:]
                            deep = sourced['deep'][:]
                            strat = sourced['strat'][:]
                            try:
                                g.createGroup('%05i'%tt)
                            except (ValueError,RuntimeError):
                                pass
                            f = g['%05i'%tt]

                            for i,j in (('con',con),('deep',deep),('strat',strat)):
                                try:
                                    f.createVariable(i,'f',('time','lat','lon'))
                                except (ValueError,RuntimeError):
                                    pass
                                
                                f[i][:] = j
                                f[i].units=' '
                                f[i].standard_name='CAF of %s. clouds'%i
                                f[i].mul=In[0]
                                f[i].add_c=In[1]
                            del sourced,con,deep,strat
                            sys.stdout.write('ok\n')
                            sys.stdout.flush()

                    #if os.path.isfile(fname.replace('.nc','_%02i.nc'%ii)):
                    #    os.remove(fname.replace('.nc','_%02i.nc'%ii))
                except (IOError,OSError):
                    pass


##############################################################################
##############################################################################

def CD_sensitivity(**kwargs):

    """
    This method creates a variance based sensitivity analysis
    based on a various C,D values as a sobol sequence
    """
    from mpi4py import MPI
    root = 0

    if type(kwargs['rank']) == type(None):
        #Process based mpi
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        name = comm.name
        scatter = True
    else:
        rank = int(kwargs['rank'])
        size = int(kwargs['size'])
        scatter = False
        if size <= rank:
            size = rank+1

    def Print(*arg,**kw):
        RegionName=arg[-2]

        thc=arg[-1]
        time=arg[0].start+datetime.timedelta(hours=24*arg[1])
        rank = kw['rank']
        perc = kw['perc']
        time=time.strftime('%j %H:%M')
        return 'Iterating C: %.2f D: %.2f'%(arg[0].C[0,0],arg[0].D[0,0])\
                +' with tdiff: %.2f in  phase %.2f and + %02f at %s (%02i | %3.2f)\n'\
                %(thc,arg[0].phase,arg[0].add_c,time,rank,perc)

    problem={
            'num_vars':2,
            'names':['C','D'],
            'bounds':[[0,2],[0,2]]}
    #Create the random input data and split it into n samples
    if kwargs['sobol']:
        #Create the sobol sequence
        if not os.path.isfile('saltelli_CD.txt'):
            from SALib.sample import saltelli
            Global_param_values = saltelli.sample(problem,500)
            np.savetxt('saltelli_CD.txt',Global_param_values)
        else:
            Global_param_values=np.loadtxt('saltelli_CD.txt')
    else:
        CC = np.linspace(0,2,100)
        DD = np.linspace(0,2,100)
        Global_param_values = np.array([np.array([i,j]) for i in CC for j in DD])
    samples = np.array_split(Global_param_values,size)

    if scatter:
        if rank == root:
            sendbuf=samples
        else:
            sendbuf = []
        #Scatter the buffer to all processes, if empty it is dropped
        X = comm.scatter(sendbuf,root)
    else:
        X = samples[rank]
    
    #Now setup the model
    start=kwargs['start']
    end=kwargs['end']
    boxes=kwargs['boxes']
    del kwargs['boxes'],kwargs['start'],kwargs['end']
    Cfg = Config('constants.config')
    Cfg.start='1999-01-01_00:00'
    #One day should be considered as spin-off how many indices is one day
    nn=int(1/Cfg.dt*6)
    Cfg.tend += 6
    #Check for preconditioned dtmax, phase and add_c
    for checkv  in ('dtmax','phase','add_c'):
        if checkv in kwargs.keys():
            Cfg[checkv]=kwargs[checkv]
            del kwargs[checkv]
    #Now instanciate the model setup
    try:
        C0,D0=X[0][0],X[0][-1]
        cl = init(conf=Cfg,C0=C0,D0=D0,**kwargs)
        mean,std = None,None
        for ItNum,setup in enumerate(X):
            printkw=dict(rank=rank,perc=float(ItNum)/len(X)*100)
            C,D = tuple(setup)
            cl.D=D*np.ones_like(cl.D)
            cl.C=C*np.ones_like(cl.C)

            con,deep,strat,thc,time = integrate(cl,Print=Print,printkw=printkw)
            con=con[nn:].mean(axis=(1,2))
            deep=deep[nn:].mean(axis=(1,2))
            strat=strat[nn:].mean(axis=(1,2))
            if type(mean) == type(None) and type(std) == type(None):
                mean=[np.array([con.mean(),deep.mean(),strat.mean()])]
                std=[np.array([con.std(),deep.std(),strat.std()])]
            else:
                mean=np.append(mean,[np.array([con.mean(),deep.mean(),strat.mean()])],axis=0)
                std=np.append(std,[np.array([con.std(),deep.std(),strat.std()])],axis=0)
        out = [mean,std]
        #Send the output back to the root communicator
    except IndexError:
        out=[]

    if scatter:
        recbuf = comm.gather(out,root)
        if comm.rank == root:
            sendbuf ="done"
            recvbuf = comm.bcast(sendbuf,root)
    else:
        return

    if type(recbuf) != type(None):
        nn = 0
        while True:
            try:
                mean,std=recbuf[nn]
                In = samples[nn]
                break
            except ValueError:
                nn+=1
        if len(recbuf) > 1 and len(samples) >0 :
            for ii in xrange(1,len(recbuf)):
                if len(recbuf[ii]) > 0:
                    mean=np.append(mean,recbuf[ii][0],axis=0)
                    std=np.append(std,recbuf[ii][1],axis=0)
                    In=np.append(In,samples[ii],axis=0)
        frame = np.c_[In,mean,std]
            
        out_names = ['%s_%s'%(i,j) for i in ('mean','std') for j in ('c','d','s')]
        out = pd.DataFrame(frame,columns=problem['names']+out_names)
        fname = 'CD_sensitivity_%0.2f_%0.2f_%02f.pkl'%(Cfg.dtmax,Cfg.phase,Cfg.add_c)
        
        fname=os.path.join(os.environ['HOME'],'PhD','Data','SMCM',fname)
        out.to_pickle(fname)
        

##############################################################################
##############################################################################


def GetNames(ko,fallback='Unkown'):
    """
        Method that returns information about the location
        of given kooardinates.

        Variabls:
            ko (tuple) : lat,lon of the location
        Keywords:
            fallback   : if the google api doesn't find information
                        which string should be returned
        Returns str
    """
    
    lat,lon=float(ko[0]),float(ko[1])
    request=r=geocoder.reverse([lat,lon],short_name=False)
    name = r.country
    if type(r.state) != type(None):
        name = '%s, %s' %(r.state,r.country)

    if type(r.state) is type(None) and type(r.country) is type(None):
        ele=geocoder.elevation([lat,lon],short_name=False).meters
        if ele < 0:
            if lon > 48 and lon < 76 and lat > 3.9 and lat < 27.9:
                name='Arabian Sea'
            elif lon>79 and lon < 98 and lat > 3.9 and lat < 27.9:
                name='Bay of Bengal'
            elif lon > 40 and lon < 122:
                name='Indian Ocean'
            elif lon > 145 and lon < 163 and lat < -11 and lat > -30:
                name = 'Coral Sea'
            elif lon > 121 and lon < 180.1:
                if lat > 0 :
                    name = 'North West Pacific'
                else:
                    name = 'South West Pacific'
            elif lon > -180.1 and lon < -70:
                if lat > 0:
                    name = 'North East Pacific'
                else:
                    name = "South East Pacific"
            elif lon > -97 and lon < -90 and lat < 30 and lat > 18:
                name = 'Gulf of Mexico'
            elif lon >-89 and lon < -80 and lat <30 and lat > 21:
                name = "Gulf of Mexico"
            elif lon >-87 and lon < -73 and lat < 21 and lat > 9:
                name = 'Caribbean Sea'
            elif lon > -80 and lon < 13:
                if lat > 0:
                    name = 'North Atlantic'
                else:
                    name = 'South Atlantic'
            else:
                name = '%02f %02f' %(float(lat),float(lon))
        else:
                name = '%02f %02f' %(float(lat),float(lon))

    return name



##############################################################################
##############################################################################

def kw(**kwargs):
    import multiprocessing
    """
    Method that returns the default parameter keywordarguments
    """
    F={'run_once':run_once,'CD_sensitivity':CD_sensitivity,\
            'coastal_sensitivity':coastal_sensitivity,\
            'diurnal_cycle':diurnal_cylce}
    KW=dict(
        boxes=[None],
        start=None,
        end=None,
        seed=None,
        J0 = np.array([\
                [0.45,0.65,0.00,],
                [0.65,0.35,0.20,],
                [0.00,0.10,1.00,]\
                        ]),
        func=run_once,
        sobol=False,
        size=(multiprocessing.cpu_count() * 2) -1,
        rank=None)
    for key,value in kwargs.items():
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
    start (str)     : The start of the integration (dafualt %s)
    end   (str)     : The end of the integration (default %s)
    boxes   (list)  : The name of the boxes that are considered (default %s)
    J0  (2d-array)  : The interaction potential 
                      (default %s
                               %s
                               %s)
    func (function) : Which experiment should be done (default %s)
    sobol (bool)    : Should a sobol sample be created (dfault %s)
    rank (int)      : The number of process (default %s)
    size (int)      : The number of total proceses (default %s)
    seed (int)      : The random seed number for rand int generation (default %s)
    """%(
            sys.argv[0],
            str(kwargs['end']),
            str(kwargs['end']),
            str(kwargs['boxes']),
            np.array2string(kwargs['J0'][0], precision=4, formatter={'float_kind':lambda x: "%.2f" % x}),
            np.array2string(kwargs['J0'][1], precision=4, formatter={'float_kind':lambda x: "%.2f" % x}),
            np.array2string(kwargs['J0'][2], precision=4, formatter={'float_kind':lambda x: "%.2f" % x}),
            str(kwargs['func'].__name__),
            str(kwargs['sobol']),
            str(kwargs['rank']),
            str(kwargs['size']),
            str(kwargs['seed'])
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
            elif key.lower().startswith('box'):
                kwargs['boxes']=[i.upper() for i in value.split(',')]
            elif key.lower() == 'j0':
                for i,j in (('(','{'),(')','}')):
                    value=value.replace(i,']').replace(j,'[')
                c=value.replace('[[','').replace(']]','').split('],[')
                kwargs['J0']=np.array([i.split(',') for i in c],dtype=np.float32)
            elif key.lower() == 'sobol':
                kwargs['sobol'] = value[0].lower() in ('t','y','1',1)
            elif key.lower() == 'seed':
                if value.lower() == 'none':
                    kwargs['seed']=None
                else:
                    kwargs['seed'] = int(value)
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



