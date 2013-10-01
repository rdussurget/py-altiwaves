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
        self.etopo=subclass(os.path.join('externels%sbathy%sETOPO2v2g_f4.nc' % os.path.sep))
        self.menor=subclass(os.path.join('externels%sbathy%sbathy_menor.mat' % os.path.sep))
        self.rootDir=os.path.dirname(root_file)
        self.cptDir=os.path.join(self.rootDir,'externals%scpt-city%s' %os.path.sep)