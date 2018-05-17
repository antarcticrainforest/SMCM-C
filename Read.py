"""
 Module to read data
"""
import numpy as np
import datetime,os
from configdir import Config
from netCDF4 import date2num,num2date,Dataset as nc

class Data(object):

    def __init__(self,config):

        """
        Class to read external data that is needed to run the SMCM
            Variables:
                config (config-object): the config object that contains all
                                        vital information
        """

        self.config = config
        self.startdate = config.start
        self.idx=0

    def thc(self,tmax,time,k=1):
        """
            Method to calculate the thermal heating contrast
                Variables:
                    tmax = the maximal strength of the sea-breeze or 
                            thermal heating contrast tmax
                    
                    time (int) = the integration step in hours
                Keywords:
                    k = the period for calculating the cosine
                Returns:
                    float = the thermal heating contrast
        """
        #If observations then read them first
        if not isinstance(tmax,(float,int)):
            tmax = self.atmos(tmax,time)
        else:
            #Calculate the time of the day:
            t = (time + self.startdate.hour+self.startdate.minute/60.) % 24
            tmax *= np.cos(k*np.pi/12.*t - 3/12. * np.pi)

        return tmax #* np.cos(k*np.pi/12.*t - 3/12. * np.pi)
    
    def atmos(self,var,time):
        """
            This method reads atmospheric variables like qi and omega
            Variables:
                var (tuple) = the time and variable saved in a vector
                time (int) = the integration step in hours
            Returns:
                float = the variables that is closest to the according timestep
        """
        #time and variable are saved in a tuple
        obs_time,obs_var = var

        #Create a time number that is in the same format as the time vector 
        date = self.startdate+datetime.timedelta(seconds=int(time*60**2))
        num = date2num(date, "Minutes since 1998-01-01 00:00:00")
        #Closest index an return
        idx = np.argmin(np.fabs(obs_time-num))
        
        return obs_var[idx]
            

    def lsm(self,n,form='v'):
        """
        Creates an coastline can be either artificial or from real data
        Variables :
            n (int) the number grid cells in one direction (n x n)
        Keywords:
            frome (str) the form of the coastline:
                v vertical line seperating the grid in left and right
                h horizontal line seperating the grid in top and bottome
                i a square island in the middle of the grid
                if form is a tuple get the real world data, tuples
                should contain (lon_left,lon_right,lat_lower,lat_lupper)
        """

        #if n % 2 == 0:
        #    n+=1
        n=int(n)
        if form == 'v' or form =='i' or form =='h':
            lsm = np.ones([n,n])
            if form == 'v':
                lsm[:,:int(n/2)] = -1
                lsm[:,int(n/2):] = 1
            elif form == 'h':
                lsm[:int(n/2),:] = -1
                lsm[int(n/2):,:] = 1
            else:
                lsm *= -1
                a = 1/np.sqrt(2) * lsm.shape[0]
                b = float(lsm.shape[0])
                s = int(round((b - a)/2,0))
                e = s+int(round(a,0))

                lsm[s:e,s:e] = 1
        else:
            lsmf = nc(os.path.expanduser(landmask))
            lsm = lsmf.variables['slm'][:]
            lon = lsmf.variables['lon'][:]
            lat = lsmf.variables['lat'][:]
            if not isinstance(form,tuple):
                c=Config('boxes.txt')
                clat,clon=c[self.config.form]
                form=(clon-1.5,clon+1.5,clat-1.5,clat+1.5)
            slon=np.argmin(np.fabs(lon-form[0]))
            elon=np.argmin(np.fabs(lon-form[1]))
            slat=np.argmin(np.fabs(lat-form[2]))
            elat=np.argmin(np.fabs(lat-form[3]))
            
            self.lons=lon[slon:elon+1]
            self.lats=lat[slat:elat+1]
            lsm=lsm[slat:elat+1,slon:elon+1]
            lsm=lsm[::-1]*2 - 1

        '''
        coast = canny(lsm,sigma=1).astype(np.int8)
        coast[-1,n/2] = 1
        coast[0,n/2] = 1
        points=np.where(coast == 1)
        points=np.array([points[0],points[1]]).T
        d=np.zeros_like(coast)
        for i in xrange(d.shape[0]):
            for j in xrange(d.shape[1]):
                d[i,j]=self.dist(i,j,points)
        '''
        #return np.ones_like(lsm),np.ones_like(lsm)
        return lsm, lsm 


    def dist(self,i,j,points):

        d=np.array([i,j])
        sqrt = np.sum((d-points)**2,axis=1)
        return np.sqrt(np.min(sqrt))
