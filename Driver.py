import sys,os,datetime,pytz,warnings
from timezonefinder import TimezoneFinder
import numpy as np
import pandas as pd
from netCDF4 import Dataset as nc
from netCDF4 import date2num,num2date
from configdir import Config

def open_file(file,function,**kwargs):
    """
    This method tries to open a file and check if there is any I/O Error that
    prevents the file to be openend
    
    Varialbes:
        file (str) the filename
        function (function) the function that is used to open the file
    """
    if not os.path.isfile(file):
        raise OSError("File %s does not exsist"%file)
    while True:
        try:
            f = function(file,'r',**kwargs)
            return f
        except TypeError:
            f = function(file)
            return f
        except (IOError,OSError) as e :
            n = 10*np.random.rand(10)[np.random.randint(10)]
            warnings.warn('%s -- Sleeping for %f secs'%(repr(e), n))
            time.sleep(n)


def check(kwarg,key,default):
    """
        This function gets a dictionary and returns the value
        or if not present the default value

        Variables:
            kwarg (dict) : the key words in a dictionary
            key : the key whos value should be returned
            default : the value returne if key is not present
    """
    try:
        return kwarg[key]
    except KeyError:
        return default

def close(loc):
    """
        This function returns the name of a box that is closest to given
        coordinates

        Variables:
            loc (dict) = the locations

    """
    x,y = loc['lon'],loc['lat']
    X = np.array([x,y])

    Cfg = Config('boxes.txt')
    Y = np.array(Cfg.values())[:,::-1]
    d = np.linalg.norm(X - Y,axis=1)
    i = np.argmin(d)

    return Cfg.keys()[i]


class Reader(object):

    def __init__(self,start,titt,loc={'lon':131.63,'lat':-12.31},method='linear'):

        """
        This class gets provides some data to be fed into the SMCM
        
        Variables:
            start = the starting date and time of the integration
            titt = how many days are integrated
        Keywords:
            loc = at which geographical location should the integration take place
        Instances:
            tzone = time difference between localtime and UTC
            loc = at which geographical location should the integration take place
            start = the start as a datetime-object
            end = the end as a datetime-object


        """
        if type(loc) == type({}):
            clon,clat = loc['lon'],loc['lat']
            self.box = close(loc).upper()
        else:
            Cfg = Config('boxes.txt') #The configuration of the boxes
            clat,clon=Cfg[loc.lower()]
            self.box = loc.upper()

        tzone = TimezoneFinder().timezone_at(clon,clat)
        self.tzone=pytz.timezone(tzone)
        self.loc=dict(lon=clon,lat=clat)
        self.times=pd.date_range(start,periods=titt*24*60,freq='min',tz=tzone)
        self.start, self.end = self.times[0],self.times[-1]

        self.method = method

    @staticmethod
    def __gettdiff(t,slm,loc):
        """
        
        Static method to calculte the land-sea thermal heating contrast
        Variables:
            t (netcdf4-object): the temperature data in netcdf4 format
            slm (netcdf4-object): the land-sea mask in netcdf4 format
            loc (tuple) : indices for reading the values
        
        Returns (nd-array) 
        """

        idx,idy=loc
        #Get the land-sea mask
        mask=slm[...,idy-2:idy+2,idx-2:idx+2]
        t=t[:,idy-2:idy+2,idx-2:idx+2]
        #Get only the water and land temps
        tsea = t * (np.ma.masked_equal(mask,1)+1)
        tland = t *  np.ma.masked_equal(mask,0)
        #Average over the domain.
        tsea=np.mean(tsea,axis=1).mean(axis=1)
        tland=np.mean(tland,axis=1).mean(axis=1)
        #Return the temp. difference
        return tland - tsea
        
    def clim(self,dt=6):
        """
        This method gets a climatology of the thermal heating contrast drive the
        SMCM

        Keywords:
            dt = the time resolution of the data
        """
        head = 'Erai_t2'
        #The leading string of the data
        head = 'Erai_t2_'
        #Where is data located
        folder=os.path.join(os.environ['HOME'],'Data','ERAI','netcdf')
        
        
        fnames=[]
        slmfile = nc(os.path.join(os.environ['HOME'],'Data',
            'Erai_slm_tropics_0.75.nc'))
        
        
        idx=np.argmin(np.fabs(slmfile.variables['lon'][:]-self.loc['lon']))
        idy=np.argmin(np.fabs(slmfile.variables['lat'][:]-self.loc['lat']))
        mask=slmfile.variables['lsm']
        #The months
        months=dict([(i,np.zeros(24/int(dt))) for i in range(1,13)])
        nmonths=dict([(i,0) for i in range(1,13)])
        while self.start<=self.end:
            f=os.path.join(folder,'%i'%(self.start.year),\
                    self.start.strftime(head+'%Y_%m_%d.nc'))
            if os.path.isfile(f):
                nf = nc(f)
                t=nf.variables['t2']
                months[self.start.month]+=self.__gettdiff(0,t,mask,(idx,idy))
                nmonths[self.start.month]+=1
            self.start += datetime.timedelta(days=1)
        R.x=np.arange(0,24+dt,dt)*60
        xnew=np.arange(0,24*60,1)
        ary=np.zeros(len(R.x))
        R.dtemp={}
        out={}
        for mon in months.keys():
            dtemp=months[mon]/nmonths[mon]
            ary[:-1]=dtemp
            ary[-1]=dtemp[0]
            R.dtemp[mon]=ary
            f=interp1d(R.x,ary,kind=self.method)
            out[mon]=f(xnew)
        slmfile.close()
        return xnew,out
    
    def get_time(self,time,tseries):
        """
            This method get's the information in a time series for certain
            times

            Variabels:
                time (int) : the times that are considered
                tseries (pandas time series): the time series
        """
        return (tseries.index.hour == time) & (tseries.index.minute == 0)
    def winddirchange(self,uwind,vwind):
        """

        This method calculates the wind-direction change

        Variables:
            uwind (pandas time series): the uwind series
            vwind (pandas time series): the vwind series
        """
        #For simblicity convert the Dir result from radians to deg
        DperR = 180/np.pi
        #Calculate Direction
        Dir = np.arctan2(-uwind,-vwind) * DperR
        #Calculate the change in winddirection
        return (Dir.diff()+180) % 360 -180
    def windspeedchange(self,uwind,vwind):
        """

        This method calculates the wind-speed change

        Variables:
            uwind (pandas time series): the uwind series
            vwind (pandas time series): the vwind series
        """
        #Calculate Windspeed
        WS = (uwind**2+vwind**2).apply(np.sqrt)
        #Calculate the change in windspeed
        return WS

    def __time_adopt(self,**kwargs):
        """
        This private method increses or decreases the timperiod that is
        considered by a certain magnitude

        Keywords:
            the amounts of the shift
        """
        
        self.start-=datetime.timedelta(**kwargs)
        self.end+=datetime.timedelta(**kwargs)
        self.times=pd.date_range(self.start,self.end,freq='min',tz=self.times.tz)
    
    @staticmethod
    def __winddir(dataframe):
        """
            Private static method to calculate the change in winddirection

            Variables:
                dataframe the dataframe that contains the data
        """
        #For simblicity convert the Dir result from radians to deg
        DperR = 180/np.pi
        #Calculate Direction
        Dir = np.arctan2(-dataframe['uwind'],-dataframe['vwind']) * DperR
        #Calculate the change in winddirection
        return (Dir.diff()+180) % 360 -180

    @staticmethod
    def __windspeed(dataframe):
        """
        Method to calculate the windspeed

        Variables:
            dataframe : the dataframe that contains the data

        """

        return (dataframe['uwind']**2 + dataframe['vwind']**2).apply(np.sqrt)
    
    @classmethod
    def __filter(cls,dataframe,old_data,times,windspeed=8,windchange=3,thc=1.5,\
            dirchange=90.):
        """
        This private class method applies the actual filter from Borne et. al

        Variables:
            dataframe = the dataframe the contains the data where the filter
                        is applied
            old_data   = after filter the data is rolled back to the old index
                        what does the old index look like
            times      = the times (self.times)
        Keywords:
            windsppeed  = the max windspeed for the filter [8 m/s]
            windchange  = the max change of windspeed [3 m/s]
            thc       = the min thermal heating contrast [1.5 K]
            dirchange   = the change in wind directioin [90 deg]

        """

        #Calculate the change in wind direction and the speed
        Dir_c = cls.__winddir(dataframe)
        WS = cls.__windspeed(old_data)
        WS_c = WS[dataframe.index].diff()
        #Apply the boolean comparsion
        Dir_c_b = (Dir_c.abs() <= dirchange ).astype(np.int32) #Change in Winddir
        WS_b    = (WS[dataframe.index].abs() <= windspeed).astype(np.int32) #Only calm winds
        WS_c_b  = (WS_c.abs() <= windchange).astype(np.int32) #small wind change
        Td_b    = (dataframe['tdiff'].abs() >= thc).astype(np.int32)#high THC

        #Now evaluate the booleasn
        idx = (Dir_c_b ==1 ) & (WS_b == 1) & (WS_c_b == 1) & (Td_b == 1)
        Bo = pd.Series(np.zeros(len(dataframe)),index=dataframe.index)
        #Now scale the strength of the sea breeze
        Bo[idx] = 1





        #Now reindex everything
        tindex = pd.date_range(times[0],times[-1],freq='min',
                tz=dataframe.index.tz)
        Bo = Bo.reindex(WS.index).interpolate(method='nearest')
        old_data = pd.concat([old_data['tdiff'],WS,WS_c,Dir_c],axis=1,\
                keys=('tdiff','ws','wc','dc'))
               # .reindex(tindex).interpolate(method='from_derivatives')
        sign = np.sign(old_data.tdiff.values)
        data = (windspeed - old_data.ws.abs()).abs()/windspeed * \
                old_data.tdiff/thc
        return (Bo * data).reindex(tindex).interpolate(method='linear')
    @staticmethod
    def __getIndex(dataframe,hours):
        """
        This method returns only certain time steps according to given hours

        Variabes:
            dataframe = the pandas dataframe or timeseries that is treated
            hours     = the full hours that should be returned from the input
        """
        
        #we need only full hours, get them
        idx=np.where((pd.Index(hours).get_indexer(dataframe.index.hour) >= 0) &\
        (dataframe.index.minute == 0))[0]
        if not len(idx) :
            #We have hourly data but probably a timeseries within an half hour
            #round timezone to full hours
            idx=np.where((pd.Index(hours).get_indexer(dataframe.index.hour) >= 0) &\
                    (dataframe.index.minute == dataframe.index[0].minute))[0]
        idx2=np.r_[idx,[len(dataframe.index)]]
        mean =np.array([dataframe.loc[dataframe.index[idx2[i:i+1]]].values.mean(axis=0)
            for i in xrange(0,len(idx2)-1)])

        return pd.DataFrame(mean,index=dataframe.index[idx],\
                columns=dataframe.keys())


    @classmethod
    def __gettimes(cls,dataframe,time,period):
        """
        This methods evaluates only timesteps that match given local times

        Variables 
            time (int) = the local times that are evaluated
            period (int) = the duration of a period

        """
        
        #We need to get the times that are of interest this can be calculated
        #from the period and the start of a period time
        S,E=time,(time - 1)%24

        #Make the periods an even cycle for 24 h
        p = 24/(24/int(period))
        if int(p) != int(period):
            sys.stdout.write('Waring changed period from %2.2f to %i\n' 
                    %(p,period))

        hours = np.unique(pd.date_range('2000-01-01 %02i'%S,'2000-01-02 %02i'%E,
                freq='%iH'%p).hour)
        return cls.__getIndex(dataframe,hours)

    def trigger(self,**kwargs):
        """
        This method calculates the trigger of the sea breezel

        Keywords:
            level (int) : the wind level in hPa (default 850)
            box (str) : the box in the database where the info is extracted
            time (in) : the time since when the change occured

        """
        box=check(kwargs,'box','COAST_04')
        time=check(kwargs,'time',13)
        level=check(kwargs,'level',700)
        dirchange=check(kwargs,'dirchange',90.)
        period=check(kwargs,'perdiod',24.)
        tdiff=check(kwargs,'tdiff',.85)
        windchange=check(kwargs,'windchange',6.)
        windspeed=check(kwargs,'windspeed',11.)
        dt = check(kwargs,'freq','1H')
        #Increase the timeperiod by two days because the trigger functions
        #truncates the data
        self.__time_adopt(days=2)
        #Get u and v winds in the from the level we are interested in
        vwind=self.boxdata('vwind_%i'%(level),times=False,freq=dt)
        uwind=self.boxdata('uwind_%i'%(level),times=False,freq=dt)
        #Get the thermal heating contrast
        thc = self.tdiff(times=False,freq=dt)
        #And create a dataframe
        data=pd.concat([uwind,vwind,thc],axis=1)
        data.columns=['uwind','vwind','tdiff']
        #Now get only the times of interest from the dataframe
        new_data = self.__gettimes(data,time,period)

        filtered = self.__filter(new_data,data,(self.times[0],self.times[-1]),
                    windspeed=windspeed,thc=tdiff,
                    dirchange=dirchange,windchange=windchange)
        self.__time_adopt(days=-2)
        return self.pddate2num(filtered[self.times].index)\
                ,filtered[self.times].values

    def boxdata(self,varname,times=True,freq='1min',shift=0):
        """
        This method should read a timeseries of a variables from a pandas
        dataframe

        Variables:
            varname (str) : the name of the variable to be considered
        Keywords:
            times (bool) : should the time vector be returned or only the data
        """
        home=os.environ['HOME']
        dbfile=os.path.join(home,'PhD','Data','LargeScale',
                'Database_fm_%s.pkl'%(self.box.lower()))
        dbfile_a = os.path.join(home,'Data','LargeScale',
                'Database_fm_%s.pkl'%(self.box.lower()))
        if not os.path.isfile(dbfile):
            dbfile=dbfile_a

        if os.path.isfile(dbfile):
            try:
                data = open_file(dbfile,pd.read_pickle)[varname].shift(shift)
            except:
                print(open_file(dbfile,pd.read_pickle).keys())
                raise(KeyError,'key not found')

        else:
            

            sys.stderr.write('File %s does not excist\n'%dbfile)
            raise(OSError,'File %s does not excist'%dbfile)

        start = (self.times[0].to_pydatetime()) - datetime.timedelta(days=1)
        end = (self.times[-1].to_pydatetime()) + datetime.timedelta(days=1)
        
        tindex = pd.date_range(
                start.strftime("%Y-%m-%d 00:00"),end.strftime("%Y-%m-%d 00:00"),
                freq = '6H')

        data = data[tindex].tz_localize('UTC')
        data = data.reindex(
                pd.date_range(data.index[0],data.index[-1],freq=freq,tz='UTC')
                )
        data = data.tz_convert(self.times.tz)
        end_idx = np.argmin(np.fabs(data.index.asi8 - self.times.asi8[-1]))
        start_idx = np.argmin(np.fabs(data.index.asi8 - self.times.asi8[0]))
        
        data = data.loc[data.index[start_idx:end_idx+1]].interpolate(method=self.method)
        data=data.dropna()
        if not times:
            return data
        return self.pddate2num(data.index),data.values
    


    
    def tdiff(self,dt=6,times=True,freq='1min'):

        """
        This method should get the thermal heating contrast from ERAI reanlysis
        data and returns a time series of the thermal heating contrast
        
        
            
        Keywords
            dt = the time resolution of the data
            times (bool) : should the time vector be returned or only the data
            frqe str : the interpolation period
        """
        
        #The leading string of the data
        head = 'Erai_t2_'
        #Where is data located
        folder=os.path.join(os.environ['HOME'],'Data','ERAI','netcdf')
        
        
        fnames=[]
        slmfile = os.path.join(os.environ['HOME'],'Data',
            'Erai_slm_tropics_0.75.nc')
       
        slm = open_file(slmfile,nc)
        #Get the indices for the locations fo avoid passing the whole
        #array
        idx=np.argmin(np.fabs(slm.variables['lon'][:]-self.loc['lon']))
        idy=np.argmin(np.fabs(slm.variables['lat'][:]-self.loc['lat']))
        mask=slm.variables['lsm']

        #convert the times back to utc because eri is in utc
        times_utc=self.times.tz_convert('UTC')
        x1 = datetime.datetime(times_utc[0].year,
                times_utc[0].month,times_utc[0].day,0,0)
        x2 = datetime.datetime(times_utc[-1].year,
                times_utc[-1].month,times_utc[-1].day,23,59)
        x1 -= datetime.timedelta(days=1)
        x2 += datetime.timedelta(days=1)
        start='%04i-%02i-%02i 00:00' %(x1.year,x1.month,x1.day)
        while x1 <= x2:
            f=os.path.join(folder,'%i'%(x1.year),\
                    x1.strftime(head+'%Y_%m_%d.nc'))
            fn=open_file(f,nc)
                
            t=fn.variables['t2']
            delta= self.__gettdiff(t,mask,(idx,idy))
            try:
                dtemp=np.append(dtemp,delta)
            except UnboundLocalError:
                dtemp=delta
            fn.close()
            x1 += datetime.timedelta(days=1)

        slm.close()
        del slm
        #Create a pandas time series because interpolation is fast
        tindex=pd.date_range(start,x2,freq='%iH'%dt,tz='UTC')
        tindex_min=pd.date_range(start,x2,freq=freq,tz='UTC')
        #Tseries should be in localtime again
        tdiff_loc = pd.Series(dtemp,index=tindex).tz_convert('UTC')
        #Reindex to the minutely output
        
        tdiff_new = tdiff_loc.reindex(tindex_min).tz_convert(self.times.tz)
        end_idx = np.argmin(np.fabs(tdiff_new.index.asi8 - self.times.asi8[-1]))
        start_idx = np.argmin(np.fabs(tdiff_new.index.asi8 - self.times.asi8[0]))
        
        tdiff_new = tdiff_new.loc[tdiff_new.index[start_idx:end_idx+1]].\
                interpolate(method=self.method).dropna()
        if not times:
            return tdiff_new
        return self.pddate2num(tdiff_new.index),tdiff_new.values
    
    @staticmethod
    def pddate2num(index):
        """
        This method converts a pandas datetime array to a num array
        
        Variable:
            index (pd index) = pandas datetime index
        """
        idx = [i.replace(tzinfo=None).to_pydatetime() for i in index]
        return date2num(idx,'Minutes since 1998-01-01 00:00:00')

    
if __name__ == "__main__":
    import matplotlib,datetime
    from matplotlib import pyplot as plt,dates
    matplotlib.rcParams.update({'font.size':18})
    form='instant'
    varname=('trigger','tdiff')
    C=Config('constants.config')
    C.start='2005-03-22_00:00'
    C.obs = 'coast_06'
    method ='linear'
    end = 365
    s = datetime.datetime.strptime(C.start,'%Y-%m-%d_%H:%M')
    #e = datetime.datetime.strptime(end,'%Y-%m-%d_%H:%M')
    #dt = (e - s).days
    dt = end
    try:
        key,value = sys.argv[1].replace('--','').split('=')
        if key.lower() == 'method':
            method = value.lower()
    except IndexError:
        pass
    R=Reader(C.start.replace('_',' '),dt,loc=C.obs,method=method)

    if form == 'instant':
        if varname[0] == 'tdiff':
            x1,var=R.tdiff()

        elif varname[0].startswith('trigger'):
            x1,var=R.trigger()
            x2,var2=R.boxdata('uwind_700')
            x4,var4=R.boxdata('vwind_700')
            var2=np.sqrt(var2**2+var4**2)
            x3,var3=R.tdiff()
            #x3,var3=R.boxdata('qi')
            #var2 /= 12.5
            #var3 = 2 * (var3)
            var2 = [var2,var3]
            co = ['r','b']

            label=('Coastal Trigger (f) []','$|V|$ [$\\frac{m}{s}$]','$\\Delta T$ [$K$]')
            
            #label=('Sea-breeze strenght[]','Thermal heating contrast [$^\\circ C$]')
        else:
            x1,var=R.boxdata(varname[0],freq='6H')
            x2,var2=R.boxdata(varname[1],freq='6H')
            var = 2- (2 * ( 1 - var))
            var2 /= 12.5
            var2[var2>2]=2
            plt.scatter(var,var2)
            plt.xlabel('Moisture')
            plt.ylabel('Instability')
            plt.show()
            label=('Dryness []','Instability []')

    else:
        x1,var=R.clim()
        var=var[1]
        
    tday=num2date(x1,'Minutes since 1998-01-01 00:00:00')
    tday2=num2date(x2,'Minutes since 1998-01-01 00:00:00')
    tday = dates.date2num(tday)
    tday2 = dates.date2num(tday2)
    fig, ax = plt.subplots()

    #v=var/-0.00604959292619*0.0692523393478
    ax.plot(tday,var,'g-',label='Trigger')
    ax2 = [ax.twinx() for i in xrange(len(label)-1)]
    for e,v in enumerate(var2):
        ax2[e].plot(tday2,v,co[e],label=label[e+1])
        if e >= 1:
            ax2[e].spines['right'].set_position(('axes', 1.05))
            ax2[e].set_ylabel(label[e+1],color=co[e],labelpad=-10)
            #ax2[e].set_frame(True)
            #ax2[e].patch.set_visible(False)
            ax2[e].tick_params(axis='y', color=co[e])
        else:
            ax2[e].set_ylabel(label[e+1])

    if varname[0].startswith('wind'):
        ax.fill_between(tday, 0, 1, where=var == 1, facecolor='blue', alpha=0.5) 
    hfmt = dates.DateFormatter('%d.%m.%Y %H LT')
    ax.xaxis.set_major_formatter(hfmt)
    
    ax.set_ylabel(label[0])
    #[ax2[e].set_ylabel(label[e]) for e in xrange(len(label)-1)]
    fig.autofmt_xdate()
    plt.legend(loc=0)
    plt.show()
    



