'''
Created on Apr 6, 2016

@author: javi
'''

class ParamWrapper:
    '''
    classdocs
    '''
    

    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def setIMax(self, imax):
        self.imax = imax
        
    def getIMax(self):
        return self.imax
    
    def setIMin(self, imin):
        self.imin = imin
    
    def getIMin(self):
        return self.imin
    
    def setIStep(self, istep):
        self.istep = istep
    
    def getIStep(self):
        return self.istep
    
    def setPiMax(self, pimax):
        self.pimax = pimax
    
    def getPiMax(self):
        return self.pimax
    
    def setPiMin(self, pimin):
        self.pimin = pimin
        
    def getPiMin(self):
        return self.pimin
    
    def setPiStep(self, pistep):
        self.pistep = pistep
    
    def getPiStep(self):
        return self.pistep
    
    def setOutputFolder(self, outfolder):
        self.outfolder = outfolder
    
    def getOutputFolder(self):
        return self.outfolder
    
    def setIVal(self, iVal):
        self.iVal = iVal
        
    def getIVal(self):
        return self.iVal
    
    def setPiVal(self, piVal):
        self.piVal = piVal
    
    def getPiVal(self):
        return self.piVal
    
    def setSectionLength(self, sect_length):
        self.sect_length = sect_length
    
    def getSectionLength(self):
        return self.sect_length
    