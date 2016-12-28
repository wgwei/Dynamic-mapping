# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 11:35:13 2013
The backgroun noise model
To calculate the background noise level, especially positions inside a enclosed
yard. The algorithm is based on the reflection model and diffraction model

new features:
    1. loading line_x_poly model to calculate the distance from source to 
        building or from receiver to building
@author: W.Wei
"""

import os
import shutil
import numpy as np
import time
import dbflib, shapelib
import math
import pickle
import logging
from BGM2_utils import dist2Dpyx, packBuildingToPKL, get_pti_pt2_ext_pt_by_dist, \
        distance_to_building, fitModel_refl, fitModel_scatter,\
        get_blocked_building_num

class Receivers():    
    def __init__(self, vertix, identify, absoluteHeight=0.0):
        self.vertix = vertix
        self.absoluteHeight = absoluteHeight
        self.identify = identify  # unique id to identify the receiver
        
class Sources():    
    def __init__(self, vertix, immissonSpectrum, absoluteHeight=0.0):
        self.vertix = vertix
        self.immissonSpectrum = immissonSpectrum
        self.absoluteHeight = absoluteHeight        
        
            
def packSourceToPKL(sourceShapeFile):
    shp = shapelib.open(sourceShapeFile)
    dbf = dbflib.open(sourceShapeFile)
    specField = ['L_63', 'L_125', 'L_250', 'L_500', 'L_1000', 'L_2000', 'L_4000', 'L_8000']
    sourceObjects = []    
    for r in xrange(shp.info()[0]):
        shpObj = shp.read_object(r)
        p_source = shpObj.vertices()[0]  # (x, y)  
        rec = dbf.read_record(r)
        specD = [rec[f] for f in specField]
        spec = [specD]
        sObj = Sources(p_source, spec)
        sourceObjects.append(sObj)
    shp.close()
    dbf.close()    
    return sourceObjects  # list of soruce object. values defined in "Sources" can be queried

def packReceiverToPKL(receiverShapeFile, IDFieldName):
    shp = shapelib.open(receiverShapeFile)
    dbf = dbflib.open(receiverShapeFile)    
    receiverObjects = []    
    for r in xrange(shp.info()[0]):
        shpObj = shp.read_object(r)
        vtx = shpObj.vertices()[0]
        identify = dbf.read_record(r)[IDFieldName]
        rObj = Receivers(vtx, identify)
        receiverObjects.append(rObj)  
    shp.close()
    dbf.close()
    return receiverObjects
        
class Model():
    def __init__(self,Cvsq, Ctsq, buildingFile, receiverFile, sourceFile, resultOutFile,\
        rZoneSize=2000.0, r2sDist=1500.0, flags=['D', 'E', 'N'], modelType='scattering'):
        ''' buildingFile, receiverFile, sourceFile and resultOutFile 
            are shape file names
            resultOutFile is also the new folder
            receiverZone is the vertix of the smaller receiver region
            SBzone is the corresponding zone of the receivers
            flags -> 'D', 'E', 'N' represent for day, evening and Night
        '''
        print "preparing to write to file: ", resultOutFile
        if not os.path.exists(resultOutFile):
            os.mkdir(resultOutFile)
        logging.basicConfig(filename=os.path.join(resultOutFile, 'log.log'),level=logging.DEBUG)
        self.sourceFile = sourceFile
        self.receiverFile = receiverFile
        self.buildingFile = buildingFile
        self.Cvsq = Cvsq
        self.Ctsq = Ctsq
        self.i = 0   # acounter to write results out
        print 'initialing...'   
        # common constants
        self.fr = np.array([63, 125, 250, 500, 1000, 2000, 4000, 8000])
        self.waveLength = 340.0/np.array(self.fr)
        self.Aweight =  np.array([-26.2, -16.1,-8.6, -3.2, 0, 1.2, 1, -1.1])   # {63:-26.2, 125:-16.1, 250:-8.6, 500:-3.2, 1000:0, 2000:1.2, 4000:1, 8000:-1.1}
        self.sHi = 0.0   # source Height♠
        self.rHi = 4.5   # receiver height
        self.modelType = modelType        
        
        # preparing to write results out        
        shutil.copy(receiverFile+'.shp', resultOutFile)
        shutil.copy(receiverFile+'.shx', resultOutFile)
        self.outw = open(os.path.join(resultOutFile, 'results.txt'), 'wb')
        
        logging.info('stat time in second')
        logging.info(time.time())
        
    def prepare(self):
        print 'Prepare source, receiver and buildings'  
        if not os.path.exists('pkData'):
            os.mkdir('pkData')        
        if not os.path.exists(os.path.join('pkData',self.sourceFile+'.pkl')):
            sourceObjects = packSourceToPKL(self.sourceFile)
            sw = open(os.path.join('pkData',self.sourceFile+'.pkl'), 'wb')
            pickle.dump(sourceObjects, sw, 2)  # protocal 2 for verion 2.x, 3for 3.x
            sw.close()
        else:
            sr = open(os.path.join('pkData', self.sourceFile+'.pkl'), 'rb')
            sourceObjects = pickle.load(sr)
            sr.close()
        if not os.path.exists(os.path.join('pkData', self.receiverFile+'.pkl')):
            receiverObjects = packReceiverToPKL(self.receiverFile, 'ID')
            rw = open(os.path.join('pkData',self.receiverFile+'.pkl'), 'wb')
            pickle.dump(receiverObjects, rw, 2)
            rw.close()
        else:
            rr = open(os.path.join('pkData',self.receiverFile+'.pkl'), 'rb')
            receiverObjects = pickle.load(rr)
            rr.close()
        if not os.path.exists(os.path.join('pkData', self.buildingFile+'.pkl')):
            buildingObjects = packBuildingToPKL(self.buildingFile, 'REL_HEIGHT', 'ID')
            bw = open(os.path.join('pkData', self.buildingFile+'.pkl'), 'wb')
            pickle.dump(buildingObjects, bw, 2)
            bw.close()
        else:
            br = open(os.path.join('pkData', self.buildingFile+'.pkl'), 'rb')
            buildingObjects = pickle.load(br)
            br.close()     
        
        print 'calculating...'       
        self.sourceObjects = sourceObjects
        self.receiverObjects = receiverObjects
        self.buildingObjects = buildingObjects
        logging.info('preparing time in second')
        logging.info(time.time())
        
    def runModel(self, distLimit):   
        logging.info('start calculation time in sedond')
        logging.info(time.time())
        if len(self.receiverObjects)>0:
            for i, receiverObj in enumerate(self.receiverObjects):
                # write arguments out. Modified on April 3rd, 2013
                print i, ' receiver'
                p_receiver = receiverObj.vertix   # (x, y) 
                IDay = np.array([0.0] * 8)
                for j, sourceObj in enumerate(self.sourceObjects):   # loop for all sources
                    # contribution of a source
                    p_source = sourceObj.vertix
                    distSR = dist2Dpyx(p_receiver[0], p_receiver[1], p_source[0], p_source[1]) 
                    if distSR<=distLimit:
                        bnums = get_blocked_building_num(p_source[0], p_source[1], p_receiver[0], p_receiver[1], self.buildingObjects, distLim=50.0)
                        ds, dr, recs, recr = distance_to_building(p_source[0], p_source[1], p_receiver[0], p_receiver[1], self.buildingObjects, bnums)
                        if recs!=-1 and recr!=-1:
                            if (ds>distSR or dr>distSR):
                                print 'ds, dr, distSR', ds, dr, distSR
                            N1point = [ds, self.buildingObjects[recs].relativeHeight]
                            N2point = [distSR-ds, self.buildingObjects[recr].relativeHeight]
                            barrierWidth = distSR-ds-dr
                            spos = [0.0,self.sHi] 
                            rpos = [distSR, self.rHi]
                            if barrierWidth>2.5:  # 29June, 2nd oct
                                rs = np.sqrt((N1point[1]-self.sHi)**2+ds**2.) # modified Jan31 weigang  
                                rr = np.sqrt((N2point[1]-self.rHi)**2+dr**2.)
                                h1 = max(N1point[1]-self.sHi, 0.01) # avoid ZeroDivisionError
                                h2 = max(N2point[1]-self.rHi, 0.01) # avoid ZeroDivisionError
                                phis = math.asin(ds/rs)
                                phir = math.asin(dr/rr)
                                
                                xrs, yrs = get_pti_pt2_ext_pt_by_dist(p_receiver[0], p_receiver[1], p_source[0], p_source[1], 1000.0)
                                bnums = get_blocked_building_num(p_source[0], p_source[1], xrs, yrs, self.buildingObjects, distLim=50.0)
                                ds3, dr3, recs3, recr3 = distance_to_building(p_source[0], p_source[1],  xrs, yrs, self.buildingObjects, bnums)
                                if ds3>200 or recs3==-1:
                                    Hs = 0.0
                                    Ws = 1000.0
                                else:
                                    Ws = max(ds3+ds, 3.)
                                    Hs = self.buildingObjects[recs3].relativeHeight
                                    
                                xsr, ysr = get_pti_pt2_ext_pt_by_dist(p_source[0], p_source[1], p_receiver[0], p_receiver[1], 1000.0)
                                bnums = get_blocked_building_num(p_source[0], p_source[1], xsr, ysr, self.buildingObjects, distLim=50.0)
                                ds4, dr4, recs4, recr4 = distance_to_building(p_receiver[0], p_receiver[1],  xsr, ysr, self.buildingObjects, bnums)
                                if ds4>200 or recr4==-1:
                                    Hr = 0.0
                                    Wr = 1000.0
                                else:
                                    Wr = max(ds4+dr, 3.)
                                    Hr = self.buildingObjects[recs4].relativeHeight   
                                                                                                 
                                if self.modelType=='scattering':
                                    sposVertical = (0.0, 0.0)
                                    rposVertical = (np.sqrt(distSR**2+(self.rHi-self.sHi)**2), self.rHi-self.sHi)
                                    mixAttenuation = fitModel_scatter(self.Cvsq, self.Ctsq, self.sHi, \
                                                self.rHi, self.fr,distSR,\
                                                sposVertical[0], sposVertical[1], N1point[0], N1point[1], rposVertical[0], rposVertical[1], N2point[0], N2point[1],\
                                                phis, phir, Ws, Wr, (self.buildingObjects[recs].relativeHeight+self.buildingObjects[recr].relativeHeight)/2.)
                                    IL = -mixAttenuation + 20.0*np.log10(distSR) + 11   # insertion loss                                                                                       
                                elif self.modelType=='FDTDfitting':
                                    [Abar, AcanE, Ainter, AroofNocan, AroofCan] = fitModel_refl(h1, h2, Hs, Hr, barrierWidth,\
                                                self.sHi, self.rHi, distSR, \
                                                N1point[0], N1point[1], N2point[0], N2point[1],\
                                                Ws, Wr, rs, rr, spos[0], spos[1],  rpos[0], rpos[1],\
                                                phis, phir,  self.waveLength)
                                    IL = -10.*np.log10(10.**(0.1*(Abar+AroofNocan))+10.**(0.1*(AcanE+AroofCan))) + Ainter + 20.0*np.log10(distSR)+11   # insertion loss
                                #calculating source contribution to immission level
                                spower = np.loadtxt(self.sourceFile+'.txt')
                                spowerDEN = spower[0, :] # use the 0h data of every category
                                IDay += 10**((spowerDEN -3. + self.Aweight - IL)/10.0)   # -3.0 is to remove the ground effect. Day level                                        
                try:
                    if IDay.all() <=0:
                        LD = np.zeros(8)
                        LweightedD = 0.0
                    else:
                        LD = 10*np.log10(IDay)                    
                        LweightedD = 10*np.log10(sum(10**(0.1*(LD))))
                except ValueError:
                    LD = np.zeros(8)
                    LweightedD = 0.0
                finally:
                    print 'invalid value detected'
                print ("Noise level: ", LD)
                print ("total A weighted: ", LweightedD)     
                self.i += 1     
                string = str(receiverObj.identify) + ' ' + ' '.join('%0.1f' %v for v in LD) + ' ' + str(LweightedD) + '\r\n'
                self.outw.write(string)
        else:
            print "No receiver in this zone!"
            print "Program will continue without writing anything out!"


def call_Model(sourceShape, receiverShape, buildingShape, resultOutFile, modelType, Ctsq=0.002, Cvsq=0.01):
    ''' modelType = 'scattering' or 'FDTDfitting' '''
    mObj = Model(Cvsq, Ctsq, buildingShape, receiverShape, sourceShape, resultOutFile, flags=['D'], modelType=modelType)
    print 'write log file...'
    argins = vars(mObj)
    for item in argins.items():
        logging.info(item)
    mObj.prepare()
    mObj.runModel(1500.0)
    mObj.outw.close()     
    
    logging.info('finish calculation time in sedond')
    logging.info(time.time())
    logging.disable(logging.CRITICAL) # stop logging        
    
if __name__=="__main__":      
    sourceShape = 'industry'
    receiverShape = 'Pmeas'
    buildingShape = 'SSM_buildings'
    resultOutFile = 'resultOutFile'
    modelType = 'FDTDfitting'
    call_Model(sourceShape, receiverShape, buildingShape, resultOutFile, modelType, Ctsq=0.002, Cvsq=1.01)
    
        

            

