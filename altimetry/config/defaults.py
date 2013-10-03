# -*- coding: utf-8 -*-
'''
Created on 7 dec 2012

@author: rdussurg
'''
import os
from altimetry import __file__ as root_file

#Path storage class
class subclass(object):
    def __init__(self,path):
        self.path=path
        self.base=os.path.basename(self.path)
        self.dir=os.path.dirname(self.path)
        self.ext=os.path.splitext(self.base)[1]
        self.set=self.set()
    def set(self):
        return os.path.exists(self.path)

#Main defaults class
class defaults(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.etopo=subclass(os.path.join('externals%sbathy' %os.path.Sep,'ETOPO2v2g_f4.nc')
        self.menor=subclass(os.path.join('externals%sbathy'  %os.path.Sep,'bathy_menor.mat')
        self.rootDir=os.path.dirname(root_file)
        self.cptDir=os.path.join(self.rootDir,'externals{0}cpt-city{0}'.format(os.path.sep))