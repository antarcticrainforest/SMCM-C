#!/usr/bin/env python2.7

"""
This module describes the coarsgraining porcesses for the SMCM

"""

import numpy as np
from model import SMCM
from interaction import Thermal
#from land import Land
#from ocean import Ocean

class Coarsgraining(SMCM):
    def __init__(self,configfile,C,D,multicol=False,**kwargs):
        """
        The coarsgraining provides is a class to calculate interaction
        potentials by the coars graining methods developed by Khouider 2014.
        (DOI: http://dx.doi.org/10.4310/CMS.2014.v12.n8.a1) the numbers given
        in this describtion refer to the mentioned math article.


        Variables:
            configfile: the filename where all constants are stored
            D = the dryness (1 - humidity)
            C = the instability constant
        Keywordes:
            **kwargs : additional keyword arguments passed to the model

        Class methods:
            birthdeath:
                method to calculate the transitions and interaction between 
                possible states in the coars grain lettice
            coarseratefthct:
                method to calculate the Hemiltonian to determine the interaction
                between the coarse grain particles
            gctranstion:
                mthod that updates the number of cloud types in a cell
        
        Instances:
            m = number of coars grain cells 
            M = m^2
            Q = squared number of the micro lettice cells
            nbp,n1,n2,n3 = constants for the calculation of the conditional
                            expectation of the Hemiltonian of coars grain 
                            interaction (3.21/3.22)
            nn0,nn1,nn2 = parameters for tha calculation of the Hemiltoninan
                         (3.27)
            Nccg,Ndcg,Nscg = mxm-Arrays representing the number of congestus
                            deep and stratiform clouds in each coars cell
            timeseries = dictionary that contains the evolution of the cloud
            area fraction of the different cloudtypes
            ntime number of total timesteps
            X  = nxn micro-lettice describing the cloudtypes
        """
        self.multicol = multicol
        try:
            self.seed = kwargs['seed']
            del kwargs['seed']
        except KeyError:
            self.seed = None
        #Get all information from the SMCM class
        SMCM.__init__(self,configfile,C,D,**kwargs)
        
        #Calculate the total number of coarse graining cells
        self.M = (self.m)**2
        self.Q = self.q**2
        #Get the thermal heating contrast:



        #Define the constants for equation 3.21 (nearest neighbor interact. dep)
        if self.nb == 8:
            self.nbp,self.n1,self.n2,self.n3=3,3,2,1
        elif self.nb == 4:
            self.nbp,self.n1,self.n2,self.n3=1,2,1,0
        self.nbv = np.arange(1,4)
        
        #Calculate the parameters to determine the coarse grain interaction Potential
        self.nn0=self.nb*(self.q-2)**2+4*(self.nb-self.nbp)*(self.q-2)+4*self.n1
        self.nn1=self.nbp*(self.q-2)+2*self.n2
        self.nn2=self.n3
    
        #Initialize the microscopic field containing the cloudtypes
        self.X = np.zeros([self.n,self.n])
        #np.random.seed(42)
        self.X = np.random.randint(4,size=(self.n,self.n))
        
        #Initialize the matrices of the coarse grain area fractions
        self.Nccg = np.zeros([self.m,self.m]) #for congestus clouds
        self.Ndcg = np.zeros_like(self.Nccg) #for deep clouds
        self.Nscg = np.zeros_like(self.Nccg) #for stratiform clouds
        
        #Initialize the dynamic core for the coasrse graining cells
        self.core=np.empty_like(self.Nccg,dtype=object)

        for I in range(self.m):
            for J in range(self.m):#Set the cloud area frac. from the micro-state
                Xq=self.X[I*self.q:(I+1)*self.q,J*self.q:(J+1)*self.q]
                self.Nccg[I,J] = np.sum(np.sum(Xq*(2 - Xq)*(3-Xq)/2.)) #???????
                self.Ndcg[I,J] = np.sum(np.sum(Xq*(Xq-1)*(3-Xq)/2.))
                self.Nscg[I,J] = np.sum(np.sum(Xq*(Xq-1)*(Xq-2)/6.))
                #Assing each cell it's one single column
                #if self.multicol:
                #    if self.lsm[I,J] == 1: #We are over land
                #        self.core[I,J] = Land(self,C0=self.C,D0=self.D,I=I,J=J)
                #    else: #We are over water
                #        self.core[I,J] = Ocean(self,C0=self.C,D0=self.D,I=I,J=J)


        ntime = np.round(self.tend/self.dt,0).astype(np.int32)
        
        #Initialize a dictionary containing the time evolution of the cloud area frac.
        self.timeseries={}
        for typ in ('congestus','deep','stratiform'):
            self.timeseries[typ]=[]
        self.timeseries['congestus'].append(self.Nccg.sum())
        self.timeseries['deep'].append(self.Ndcg.sum())
        self.timeseries['stratiform'].append(self.Nscg.sum())

    def birthdeath(self,DT,hours):
        """
        Description:
            This function calculates the birth and death rates of clouds as a 
            stochastic process 
        Variables:
            DT the lenght of a timestep in the cloud scheme (in hours)
            hours the number of hours the model has been run previously
        """
        
        self.time = 0
        #Define a macro-lattice field that is two rows and cols. bigger than
        #the original macro-lattice (working lattice)
        WK=np.zeros([self.m+2,self.m+2,3])
        WK[1:self.m+1,1:self.m+1,0]=self.Nccg
        WK[1:self.m+1,1:self.m+1,1]=self.Ndcg
        WK[1:self.m+1,1:self.m+1,2]=self.Nscg

        #set boundaries of the working lattice
        for i,ary in enumerate([self.Nccg,self.Ndcg,self.Nscg]):
            WK[0,1:self.m+1,i]        = ary[self.m-1,:]
            WK[self.m+1,1:self.m+1,i] = ary[0,:]
            WK[1:self.m+1,0,i]        = ary[:,self.m-1]
            WK[1:self.m+1,self.m+1,i] = ary[:,0]
            WK[0,0,i]                 = ary[self.m-1,self.m-1]
            WK[self.m+1,self.m+1,i]   = ary[0,0]
        WK[0,self.m+1,:]=WK[0,1,:]
        WK[self.m+1,0]=WK[self.m-1,0,:]
        
        #Get the thermal-heating contrast form the data
        tdiff=self.RD.thc(self.dtmax,self.time+hours,k=self.phase)*np.ones_like(WK)
        #There are seven possible transition rates in each macro-lattice
        #site (3.5)
        Rates = np.zeros([self.m,self.m,7])
        for I in range(self.m):
            for J in range(self.m):
                Xn = WK[I:I+3,J:J+3,:]
                #Get the neighborhood in the  land-sea-coast-mask
                mask = self.mask[I:I+3,J:J+3]
                thc = tdiff[I:I+3,J:J+3]
                #Calculate the rates according to state and thier neighbors
                Rates[I,J,:]= self.coarseratesfct(Xn,mask,thc,I,J)
        
        #Now the lattice is initializede, now iteration is done
        while self.time < DT:
            #Get the thermal-heating contrast form the data
            tdiff=self.RD.thc(self.dtmax,self.time+hours,k=self.phase)*\
                    np.ones_like(WK)

            #Calculate the propability dist. of cloudchanges on the lattice
            RminIJ = np.sum(Rates,axis=2)
            Rmin = np.sum(RminIJ)
            RminIJ = np.reshape(RminIJ, self.m**2,1)
            #Calculate a random number to wait for a transition
            np.random.seed(self.seed)
            r = np.random.rand(1)[0]
            np.seterr(all='warn')
            self.dtau = -np.log(r)/Rmin
            if np.isinf(self.dtau):
                self.dtau=2
            #Get teh CCF of the propability dist.
            probdist = np.cumsum(RminIJ/Rmin,0)
            
            #Calculate a new random number to determine where the transition 
            #occurs
            np.random.seed(self.seed)
            r = np.random.rand(1)[0]
            k = 1
            while r > probdist[k-1]:
                    k+=1
            n1 = int(k % self.m)
            if n1 == 0:
                n1 = self.m
            n2 = int((k - n1)/self.m + 1)
            
            #Calculate the CCF of the phase transitions at this specific macro
            #site
            np.random.seed(self.seed)
            #np.random.seed(None)
            r = np.random.rand(1)[0]
            RatesN1N2=Rates[n1-1,n2-1,:]
            RminN1N2=np.sum(RatesN1N2 ,axis=0)
            probdist = np.cumsum(RatesN1N2/RminN1N2)
            
            k=1
            while(r>probdist[k-1]):
                k+=1
            
            #Where the transition occurs and what type of transition is now
            #determined, let's do the transition at this site and update
            #the number of cloud types at this site
            WK[n1,n2,:] = self.gctranstion(WK[n1,n2,:],k)
            
            #Update the cloud area fractions for the cloudtypes
            self.Nccg[n1-1,n2-1] = WK[n1,n2,0]
            self.Ndcg[n1-1,n2-1] = WK[n1,n2,1]
            self.Nscg[n1-1,n2-1] = WK[n1,n2,2]
            #And update the timestep
            self.time += self.dtau
            
            # Update the boundary conditions
            if n1 == 1:
                WK[self.m+1,n2,:]=WK[1,n2,:]
            elif n1 ==self.m:
                WK[0,n2,:]=WK[self.m,n2,:]
            if n2 == 1:
                WK[n1,self.m+1,:] = WK[n1,1,:]
            elif  n2 == self.m:
                WK[n1,0,:] = WK[n1,self.m,:]
            if n1 == 1 and n2==1:
                WK[self.m+1,self.m+1,:]=WK[1,1,:]
            elif  n1 == 1 and n2 == self.m:
                WK[self.m+1,0,:] = WK[1,self.m,:]
            elif n1 == self.m and n2 ==1:
                WK[0,self.m+1,:] = WK[self.m,1,:]
            elif n1 == self.m  and n2== self.m:
                WK[0,0,:] = WK[self.m,self.m,:]

            #Now update the transtion rates according to the previous state
            #and the neighbors again
            for I in range(max(n1-1,1),min(n1+1,self.m)+1):
                for J in range(max(n2-1,1),min(n2+1,self.m)+1):
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:] = self.coarseratesfct(Xn,mask,thc,I,J)
                if n2 == 1:
                    J = self.m
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:] = self.coarseratesfct(Xn,mask,thc,I,J)
                elif n2 == self.m: 
                    J = 1
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:] = self.coarseratesfct(Xn,mask,thc,I,J)
            
            #SPECIAL TREATMENT FOR BOUNDARY SIDES
            #If the transitions should occur in a boundary lattice side
            #it needs special treatment
            if n1 == 1 :
                I = self.m
                for J in range(max(n2-1,1),min(n2+1,self.m)+1):
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:]= self.coarseratesfct(Xn,mask,thc,I,J)

                if  n2 == 1:
                    J = self.m
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:] = self.coarseratesfct(Xn,mask,thc,I,J)
                elif n2 == self.m:
                    J = 1
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:] = self.coarseratesfct(Xn,mask,thc,I,J)
            elif  I == self.m:
                I = 1
                for J in range(max(n2-1,1),min(n2+1,self.m)+1):
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:]= self.coarseratesfct(Xn,mask,thc,I,J)
                if n2  == 1:
                    J = self.m
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:]= self.coarseratesfct(Xn,mask,thc,I,J)
                elif n2 == self.m :
                    J = 1
                    Xn = WK[I-1:I+2,J-1:J+2,:]
                    #Get the potential for the coast
                    mask = self.mask[I-1:I+2,J-1:J+2]
                    thc = tdiff[I-1:I+2,J-1:J+2]
                    Rates[I-1,J-1,:] = self.coarseratesfct(Xn,mask,thc,I,J)
        
    def coarseratesfct(self,Xn,mask,thc,i,j):
        """
            Description:
                This method calculates the interaction Hamiltonian of the coarse
                grain cells and calculates from there the transition rates

            Variables:
                Xn = The coarse grained cell of clouds with its neigbourhoods
                mask = The neighborhood of the point in the land-sea-coast-mask
                thc = The thermal heating-contrast of the point
        """
        #Calculate the coars grain interaction Hamiltonian (eq. 2.27)
        Hb0 = np.zeros(3)
        #Update the interaction potential according to the sorrounding of the 
        #point
        J=Thermal(mask,thc,self.J00,self.J0,self.time,self.tend,self.ndays,\
                interact=self.interact)
        for k in range(3):
            Hb0[k] = (self.nn0*(\
                    J.j0((1,1),(k,0))*Xn[1,1,0]+J.j0((1,1),(k,1))*Xn[1,1,1]+\
                    J.j0((1,1),(k,2))*Xn[1,1,2])\
                    +self.nn1*(\
                    J.j0((0,1),(k,0))*Xn[0,1,0]+J.j0((2,1),(k,0))*Xn[2,1,0]+\
                    J.j0((1,0),(k,0))*Xn[1,0,0]+J.j0((1,2),(k,0))*Xn[1,2,0]+\
                    J.j0((0,1),(k,1))*Xn[0,1,1]+J.j0((2,1),(k,1))*Xn[2,1,1]+\
                    J.j0((1,0),(k,1))*Xn[1,0,1]+J.j0((1,2),(k,1))*Xn[1,2,1]+\
                    J.j0((0,1),(k,2))*Xn[0,1,2]+J.j0((2,0),(k,2))*Xn[2,0,2]+\
                    J.j0((1,0),(k,2))*Xn[1,0,2]+J.j0((1,2),(k,2))*Xn[1,2,2])\
                    +self.nn2*(\
                    J.j0((0,0),(k,0))*Xn[0,0,0]+J.j0((2,0),(k,0))*Xn[2,0,0]+\
                    J.j0((0,2),(k,0))*Xn[0,2,0]+J.j0((2,2),(k,0))*Xn[2,2,0])+\
                    J.j0((0,0),(k,1))*Xn[0,0,1]+J.j0((2,0),(k,1))*Xn[2,0,1]+\
                    J.j0((0,2),(k,1))*Xn[0,2,1]+J.j0((2,2),(k,1))*Xn[2,2,1]+\
                    J.j0((0,0),(k,2))*Xn[0,0,2]+J.j0((2,0),(k,2))*Xn[2,0,2]+\
                    J.j0((0,2),(k,2))*Xn[0,2,2]+J.j0((2,2),(k,2))*Xn[2,2,2])
            Hb0[k]/=self.Q**2
        
        N0 = self.Q-np.sum(Xn[1,1,:]) #number of clear sky sites
        #Death rates (2.17)
        mR10=self.R10[i,j]*Xn[1,1,0]
        mR20=self.R20[i,j]*Xn[1,1,1]
        mR30=self.R30[i,j]*Xn[1,1,2]
        
        #Birth rates (2.17)
        mR01=self.R01[i,j]*np.exp(Hb0[0])*N0
        mR02 = (self.a01[i,j]*np.exp(Hb0[1])+ self.a02[i,j]*np.exp(Hb0[2]) )* N0

        #conversion from deep to stratiform (Hamiltonian from 2.22)
        Hb23p = Hb0[2] - Hb0[1] + self.nn0*(self.J0[1,1] - self.J0[2,1])/self.Q**2
        mR23=Xn[1,1,1]*self.R23[i,j]*np.exp(Hb23p) #2.17

        ##conversion from congestus to deep (Hamiltonian from 2.22)
        Hb12p = Hb0[1] - Hb0[0] + self.nn0*(self.J0[0,0] - self.J0[1,0])/self.Q**2
        mR12=self.R12[i,j]*Xn[1,1,0]*np.exp(Hb12p) #2.17

        return np.array([mR01,mR02,mR12,mR23,mR10,mR20,mR30])

    def gctranstion(self,Xn,k):
        """
            Description:
                this method updates the number of clouds in a coars grain cell
                according to the transition that should be made
            Variables:
                Xn : array containing the cloud numbers
                k : number representing the transition that occurs
        """
        switch = {1:[1.,0,0], #birth of congestus
                2:[0,1.,0], #birth of deep
                3:[-1.,1.,0], # congestus -> deep
                4:[0,-1.,1.], # deep -> stratiform
                5:[-1.,0,0], # death of congestus
                6:[0,-1.,0], #death of deep
                7:[0,0,-1.]} #death of stratiform
        return Xn+np.array(switch[k])
