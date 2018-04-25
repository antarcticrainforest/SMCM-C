#!/usr/bin/env python2.7


"""
Module to calculate any kind of new interaction potentials for coastal
points
"""

import numpy as np
import sys

class Thermal(object):
    """
    Class to calculate interaction by the thermal heating contrast in
    the coastal vicinity. 

    """

    def __init__(self,mask,deltaT,J00,J0,time,end,ndays,interact=True):

        """
        Vriables:
            mask (2d-array): the point and its neighborhood as a 
                                land/ocean/coast-mask type
            deltaT (2d-array): the thermal heating contrast of the neighborhood
            J00 (float): the background interaction potential
            J0 (nd-array): the interaction potential
        Methods:
            j0 calculates the interaction potential form the thermal heating 
                contrast
        Instances:
            p = the background interaction potential when the thermal heating
                contrast is considered
            J0 = the interaction potential as given by the SMCM

        """

        #land/ocean = 0, coastal land = 1 , coastal sea = -1
        
        #If were are in a pure continental or oceanic region we don't need
        #to do anyting here 
        """
        if len(np.unique(mask)) == 1 and mask.mean() == 0. :
            self.p=np.ones_like(mask)*J00

        elif 0. in mask: #Point has at least one ocean/land point
            olmask = (np.ma.masked_not_equal(mask,0)+1).filled(0)
            cmask = np.ma.masked_equal(mask,0)*deltaT
            jj0 = 1/np.pi*(np.pi + np.arctan(cmask*np.pi/2.))
            self.p=(jj0*J00).filled(0) + olmask*J00
        else: #A pure coastal point
            try:
                cmask = np.ma.masked_equal(mask,0)*deltaT
            except ValueError:
                print mask.shape,deltaT.shape
                exit()
            jj0 = 1/np.pi*(np.pi + np.arctan(cmask*np.pi/2.))
        """
        #jj0 = 1/np.pi*(np.pi/2 + np.arctan(mask*deltaT*np.pi/14.))
        #jj0 = ((deltaT/5.)+1) * 0.5

        j00=np.ones_like(mask)
        """
        if not J00 == 0 and interact:
            for i in xrange(mask.shape[0]):
                for j in xrange(mask.shape[1]):
                    tdiff = -deltaT[i,j,0]*mask[i,j]
                    if tdiff >= 0:
                        j00[i,j]= 1
                    else:
                        j00[i,j] = 0
        """
        self.p=j00
        self.J0 = J0
    
    def j0(self,yx,kl):
        """
            This method returns one value for the interaction potential
            depending on the location and the type transion

            Variables:
                xy (tuple) : the location in the land-sea-coast-mask
                kl (tuple) : the indices for the interaction potential
            Returns:
                float : the value for the interaction potential tacking the
                thermal heating contrast into account
        """

        y,x=yx
        k,l=kl
        #Get the transition times the thermal heating contrast at that spot
        return self.J0[k,l] * self.p[y,x]




