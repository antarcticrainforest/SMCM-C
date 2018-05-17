#!/usr/bin/env python2.7

import numpy as np
from configdir import Config
from Read import Data
from Driver import Reader
import datetime
"""
    This python module provides the smcm class to run the
    smcm and gather all information
"""

def _kwargs(kwargs,key,alt):
    """
    Description:
        Helper function the checks for the existance of a keyword argument and 
        returns either the keyword argument or an given alternative

    Variables:
        kwarg (dict-type) : dict of type **kwargs
        key : key of the kwarg
        alt: value that should be returned if keword argument is not set
    """

    try:
        return kwargs[key]
    except KeyError:
        return alt

class SMCM(object):
    @staticmethod
    def __gamma(x,y=0,c=0):
        """
        Helper-function to define the gamma function:
        """
        x+=c*y
        return max(1 - np.exp(-1*x),0.001)
    @staticmethod
    def __c(thc):
        """
        Static mehtod to calculate a skewed dryness parameter
        """
        d = np.arctan(thc+1.5574077246549023)+np.pi/2.
        #d2 = 2*np.exp(thc)/(1+np.exp(thc))
        return (d * 11./(9.*np.pi))**2

    @staticmethod
    def __rho(R01,R02,R12,R23,R30,R20,R10,num=0):
        """
        This static method should calculate the equilibrium distribution rho
        """
        
        #the equilibrium dist rho
        #rho = np.zeros(4)
        if (R10+R12) == 0:
            return  0
        else:
            rho = R01/(R10+R12)
        if num == 1 or num == 2:
            if R12*R01 == 0:
                rho = 1./(R20+R23)*(R02)
            else:
                rho = (R02 + R12*rho)/(R20+R23)
            if  num == 2:
                return  R23*rho / R30
        elif num == -1 or num == 3:
            rho = 1
        return rho
    @staticmethod
    def _makeMatrix(array):
        """
            This function takes an array and enlarges it by 2 fields
        """
        mask=np.zeros([array.shape[0]+2,array.shape[1]+2])
        mask[1:-1,1:-1]=array
        mask[-1] = mask[-2]
        mask[0] = mask[1]
        mask[:,-1] = mask[:,-2]
        mask[:,0] = mask[:,1]

        return mask


    def __init__(self,configfile='constants.config',C=0.5,D=0.55,scaling=False,
            **kwargs):
        """
        The smcm calls provides is a class calculates stand alone transition
        potentials without the interaction of neighbors or coars graining.
        The numbers refer to the equations in the article by Khouider 2014.
        (DOI: http://dx.doi.org/10.4310/CMS.2014.v12.n8.a1) 

        Variables:
            configfile: the filename where all constants are stored
            D = the dryness (1 - humidity)
            C = the instability constant

        Instances:
            R10,R20,R30,R12,R23,R01,R02 = transition rates without respect to
            interaction with other cells see eqn 2.17

            rho = the equilibrium interaction distribution, see eqn 2.16
            N squared number of the microscopic sites
            time,nt,dtau,timescg = intial values of the time iteration
            J0 = interaction-matix
            RD = read-data object to read all external data
        """
        if isinstance(configfile,str):
            conf=Config(configfile)
        else:
            conf=configfile
        #Convert the startdate to a datetime object
        conf.start = datetime.datetime.strptime(conf.start,'%Y-%m-%d_%H:%M')
        for i,j in conf.items():
            setattr(self,i,j)
       
        self.J0 = _kwargs(kwargs,'J0',np.array([[1.,0.,0.],[0.,.5,.2],\
                [0,0.2,0.5]]))
        
        #Calculate the number of days
        self.ndays = self.tend/24.
        if int(self.ndays) < self.ndays:
            self.ndays = int(self.ndays)+1
        else:
            self.ndays = int(self.ndays)
        
        #Get the land-sea-mask if one is given
        #Create the Read data object
        self.RD=Data(conf)
        self.m = self.n/self.q
        self.lsm,self.dist = self.RD.lsm(self.m,form=self.form)
        self.mask = SMCM._makeMatrix(self.lsm)
        #check if we prescirbing with an real or an artificial thermal heating
        #contrast
        #If the model is driven with observastions we need to read the data
        if self.obs:
            self.Reader = Reader(self.start,self.ndays,loc=self.obs.upper())
            self.qi= self.Reader.boxdata('qi')
            self.ki = self.Reader.boxdata('kindex')
            c=Config('boxes.txt')
            lats,lons=c[conf.obs.lower()]
            D = 2*(1 - np.ones_like(self.mask)*self.RD.atmos(self.qi,0))
            C = min(self.RD.atmos(self.ki,0)/12.5,2)*np.ones_like(self.mask)

            #There are two possibilities run the model with a real-world coast
            #or with an idialized coast but real world data (faster)
            #This is set in form
            if self.form.lower().startswith('coast_'):
                #Update the coast form to the geographical components
                self.form=(lons-1.5,lons+1.5,lats-1.5,lats+1.5)
        else:
            #No observations highly idealized case
            C = np.ones_like(self.mask)*C
            D = np.ones_like(self.mask)*D
            if scaling:
                C /= c**2
                D[self.mask==1] = D[self.mask==1] / self.tebmtembar_l
                D[self.mask==-1] = D[self.mask==-1] / self.tebmtembar_o

            self.qi,self.ki = None,None
        #If dtmax is not set this means we are having either observations or 
        #no thermal heating contrast at all
        if  type(self.dtmax) == type(None) and self.interact:
            self.dtmax=self.Reader.trigger()
        elif not self.interact:
            self.dtmax = 0

        #Should some output be plotted:
        self.plot=_kwargs(kwargs,'plot',True)
        #What plottype is prefered:
        self.plottype = _kwargs(kwargs,'plottype','ts')
        #Background rates calculate via equ. 2.17:
        #total number of micorscopic sites
        #Calculate the normalized values of D and C
        c = np.sqrt(self.N2)*self.ZT/np.pi

        #Calculate the number of coarse-grain-cells
        if isinstance(self.form,tuple) or self.form not in ('v','h','i'):
            self.m=self.lsm.shape[0]
            self.n=self.m*self.q

        thc = self.mask*self.RD.thc(self.dtmax,0,k=self.phase)
        alpha_bar=self.Hm*self.N2*self.theta0/self.g

        self.N = self.n**2

        #Define all timesteps needed for the iteration of the process
        self.time,self.nt,self.dtau=self.start.hour,1,0
        self.tend+=self.time

        #Calculate the equilibrium potential matrix
        self.J0=self.J00*self.J0



        self.gamma_f = np.vectorize(SMCM.__gamma)
        self.rho_f = np.vectorize(SMCM.__rho)
        self.c_f = np.vectorize(SMCM.__c)
        thc *= self.mul
        y=self.c_f(thc)
        c=self.add_c
        if self.interact:
            #birth of congestus:
            self.R01 = y*self.gamma_f(C,thc,c)*self.gamma_f(D,thc,-c)/self.tau01
            #birth of deep:
            self.R02 = y*self.gamma_f(C,thc,c)*(1-self.gamma_f(D,thc,-c))/self.tau02
            #conversion of congestus to deep:
            self.R12 = y*self.gamma_f(C,thc,c)*(1-self.gamma_f(D,thc,-c))/self.tau12
        else:
            #birth of congestus
            self.R01 = self.gamma_f(C,y=0)*self.gamma_f(D,y=0)/self.tau01
            #birth of deep
            self.R02 = self.gamma_f(C,y=0)*(1-self.gamma_f(D,y=0))/self.tau02
            #conversion of congestus to deep:
            self.R12 = self.gamma_f(C,y=0)*(1-self.gamma_f(D,y=0))/self.tau12

        #conversion from deep to stratiform
        self.R23 = np.ones_like(C)/self.tau23
        #decay of congestus:
        self.R10 = self.gamma_f(D) / self.tau10
        #decay of deep:
        self.R20 = (1-self.gamma_f(C)) / self.tau20
        #decay of stratiform:
        self.R30 = np.ones_like(C)/self.tau30

        self.D,self.C = D,C
        
        self.rho = np.array([self.rho_f(self.R01,self.R02,self.R12,self.R23,self.R30,\
                self.R20,self.R10,num=x) for x in xrange(4)])
        self.rho = self.rho.T.swapaxes(0,1)
        
        self.a01 = (self.rho[...,1]*self.R20-self.rho[...,0]*self.R12)
        self.a02 = self.rho[...,2] * self.R30
        
        self.norm_rho = self.rho/self.rho.sum(axis=-1)[...,np.newaxis]
    def _update(self,time,m=0,n=0):
        """
        Method that is called when the initialized D or some entries are set

        Arguments:
            n,m = the indices of the arrays
        """
        #Get the instance that has to be changed
        if self.interact:
            thc = self.mask*self.RD.thc(self.dtmax,time,k=self.phase)
        else:
            thc = 0
        if self.obs:
            D = 2*(1 - np.ones_like(self.mask)*self.RD.atmos(self.qi,time))
            C = min(self.RD.atmos(self.ki,time)/12.5,2)*np.ones_like(self.mask)
        else:
            D,C=self.D,self.C
            C[C<0]=0
        mulclol=True
        c=self.add_c
        thc *= self.mul
        y=self.c_f(thc)

        #self.__setattr__(key,change) #Overwirte the instance
        if self.interact:
            #birth of congestus:
            R01 = y*self.gamma_f(C,thc,c)*self.gamma_f(D,thc,-c)/self.tau01
            #birth of deep:
            R02 = y*self.gamma_f(C,thc,c)*(1-self.gamma_f(D,thc,-c))/self.tau02
            #conversion of congestus to deep:
            R12 = y*self.gamma_f(C,thc,c)*(1-self.gamma_f(D,thc,-c))/self.tau12
        else:
            #birth of congestus
            R01 = self.gamma_f(C,y=0)*self.gamma_f(D,y=0)/self.tau01
            #birth of deep
            R02 = self.gamma_f(C,y=0)*(1-self.gamma_f(D,y=0))/self.tau02
            #conversion of congestus to deep:
            R12 = self.gamma_f(C,y=0)*(1-self.gamma_f(D,y=0))/self.tau12
        
        #conversion from deep to stratiform
        R23 = np.ones_like(C)/self.tau23
        #decay of congestus:
        R10 = self.gamma_f(D) / self.tau10
        #decay of deep:
        R20 = (1-self.gamma_f(C)) / self.tau20
        #decay of stratiform:
        R30 = np.ones_like(C)/self.tau30

        rho = np.array([self.rho_f(R01,R02,R12,R23,R30,R20,R10,num=x)\
                for x in xrange(4)])

        self.rho[:] = rho.T.swapaxes(0,1)
        R30 = R30 * np.ones_like(D)
        self.norm_rho[:] = self.rho/self.rho.sum(axis=-1)[...,np.newaxis]
        self.a01 = (self.rho[...,1]*R20-self.rho[...,0]*R12)
        self.a02 = self.rho[...,2] * R30

        self.R01 = R01
        self.R02 = R02
        self.R12 = R12
        self.R23 = R23
        self.R30 = R30
        self.R20 = R20
        self.R10 = R10
        
        setattr(self,'D',D)
        setattr(self,'C',C)

if __name__ == "__main__":
    import sys
    from coarsgraining import Coarsgraining
    D = 0.2 # Test-run, the dryness parameter
    C = 0.6 # Test-run, the instability parameter
    CG = Coarsgraining('constants.config', C, D) #The model object
    #Time-vector
    tt = int(0)
    t = []
    #Construct the time vector
    while tt <= int(CG.tend * 60):
        t.append(float(tt)/60.)
        tt += int(CG.dt * 60)
    t = np.array(t)
    caf = np.zeros([3,len(t)])


    for ii,tt in enumerate(t):
        sys.stdout.flush()
        sys.stdout.write('\rRunning model %03i/%03i '%(tt,len(t)))
        sys.stdout.flush()
        CG.birthdeath(CG.dt, tt) #Call the birth-death process
        CG._update(tt, m=None, n=None) #And update all parameters
        caf[0,ii] = np.mean(CG.Ndcg)/CG.q**2 #CAF of congestus clouds
        caf[1,ii] = np.mean(CG.Nccg)/CG.q**2 #CAF of deep clouds
        caf[2,ii] = np.mean(CG.Nscg)/CG.q**2 #CAF of stratiform clouds

    #Plot the results
    import matplotlib
    from matplotlib import pyplot as plt
    font = {'family' : 'normal', 'weight' : 'normal',  'size'   : 22}
    matplotlib.rc('font', **font)
    plt.plot(t, caf[0], label='Congestus', lw=2)
    plt.plot(t, caf[1], label='Deep', lw=2)
    plt.plot(t, caf[2], label='Stratiform',lw=2)
    plt.legend(loc=0)
    plt.xlabel('Time [hours]')
    plt.ylabel('Cloud Area Fraction []')
    plt.show()

