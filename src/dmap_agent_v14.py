# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 17:42:05 2013
Using LMS interpolation to smooth the noise map
@author: W.Wei

Map2Web - Grid Interpollation by Luc Dekoninck
"""
import sys, os
import numpy

#from PIL import Image
#from PIL import ImageOps
#from PIL import ImageFont
#from PIL import ImageDraw

# append the core path
#sys.path.append(r'U:\SSM-Project\core\knowledge-extraction\python')
#here = os.path.dirname(os.path.realpath(__file__))
#sys.path.append(here + "\\..\\..\\..\\core\\knowledge-extraction\\python")
#sys.path.append(here + "\\..\\..\\..\\visualization\\Map2Web")
#import connect
#import database
#import agents
import datetime, time
from dnmap_v14 import MapLMS, Add_Extra_Source, write_log, calc_squared_err, exp_Lx_Lmeas, timehour_to_den
#import Map2Web as mwebsrcCateg_contribution(foiIDs, epsilon, delta, timeHour,  Cn=4, prop_j=3, freqn=8)
    
#========================================================
#==================  run agent  =========================
class LocaMapLMSagent():
#class MapLMSagent(agents.AcousticsAgent):
    """ to calculate the dynamic map every 15min by LMS interpolation method
    """
    agentID = 201406 # change to correct procedure_id in warehousedb
    server = 'ideadb' # therion4 for testing, ideadb for roll-out
    foiIDs = [202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 220] # 206 and 220 for cheking, so remove them from the list
    PROCID  = 3#201302 #used to querry data from ideadb server, warehousedb database   
    
    folderName = input('file name: ')
    
    yn = raw_input('input betaEps, muEps manually? y-->yes    n-->no: ')
    if yn.lower()=='y':
        REF = input('reference case 1 -- yes; 0 -- no: ')
        if REF==1:
            betaEps =  1
            muEps = 0
            betaDelta= 1
            muDelta = 0
            muEpsExt = [0]
            muDeltaExt = [0]
            betaEpsExt = [1]
            betaDeltaExt = [1]
        else:
            REF = 0
            betaEps =  input('betaEps  : ')# 0.25
            muEps =    input('muEps    : ')#0.05
            betaDelta= input('betaDelta: ')#0.004
            muDelta=   input('muDelta  : ')#15
            
            betaEpsExt = [0.0]#[0.001] #this is used in version 13
            muEpsExt = [5e-4]#[8e-3]
            betaDeltaExt = [0.0]#[2e-6] # this is used in version 13
            muDeltaExt = [5e-4]#[8e-3]
    else:
        REF = 0
        betaEps =  0.25 #0.25 # this is used in version 13
        muEps = 0.001 # 0.05  #this is used in version 13
        betaDelta= 0.004
        muDelta=   0.0015
        
        betaEpsExt = [0.001]#[0.001]
        muEpsExt = [2e-4] #[2e-6]
        betaDeltaExt = [2e-6]#[2e-6] # this is used in version 13
        muDeltaExt = [8e-4]
    
    extObj = Add_Extra_Source(nameOfMeasFoiPKL='foi_ATTMdict_industry.pkl', \
        nameOfAffectRcvPKL='receiver_categATTMdict_IDs_Leq_industry.pkl',\
        sourceSpecFile='extraSourceSpec.txt', sourceNum=4)
    extObjs = [extObj]
    mapObj = MapLMS(folderName, extObjs,  betaEps, muEps, betaDelta, muDelta,\
                    betaEpsExt, muEpsExt,betaDeltaExt, muDeltaExt)

    firstRun = 'yes'
    if firstRun=='yes':   
        epsilonDay = numpy.zeros(mapObj.Cnx) # 4 is soruce category, 8 is freq. [[63,125,250,...], [], [], ...]
        deltaDay = numpy.zeros((mapObj.Cnx, len(mapObj.propa_j)))
        epsilonEve = numpy.zeros(mapObj.Cnx) # 4 is soruce category, 8 is freq. [[63,125,250,...], [], [], ...]
        deltaEve = numpy.zeros((mapObj.Cnx, len(mapObj.propa_j)))
        epsilonNig = numpy.zeros(mapObj.Cnx) # 4 is soruce category, 8 is freq. [[63,125,250,...], [], [], ...]
        deltaNig = numpy.zeros((mapObj.Cnx, len(mapObj.propa_j)))
        
        epsilonSpecDay = numpy.zeros([mapObj.Cnx, 8])
        epsilonSpecEve = numpy.zeros([mapObj.Cnx, 8])
        epsilonSpecNig = numpy.zeros([mapObj.Cnx, 8])
        
                
    fileEps, fileDelta, fileEpsSpec = [], [], []
    for den in ['Day', 'Eve', 'Nig']:
        feps = open(mapObj.fmap+'\\'+'epsilon'+den+'.txt', 'wb')
        fdelta = open(mapObj.fmap+'\\'+'delta'+den+'.txt', 'wb')
        fepsSpec = open(mapObj.fmap+'\\'+'epsilonSpec'+den+'.txt', 'wb')
        fileEps.append(feps)
        fileDelta.append(fdelta)
        fileEpsSpec.append(fepsSpec)
        
    
    fileObjs = []
    for foiID in foiIDs:
        fobj = open(mapObj.fmap+'\\'+'Lx-Lm-foi'+str(foiID)+'.txt', 'wb')
        fileObjs.append(fobj)
    fileObj10 = []
    for foiID in foiIDs:
        fobj10 = open(mapObj.fmap+'\\'+'Lx-Lm-foi'+str(foiID)+'_L10.txt', 'wb')
        fileObj10.append(fobj10)
    fileObj90 = []
    for foiID in foiIDs:
        fobj90 = open(mapObj.fmap+'\\'+'Lx-Lm-foi'+str(foiID)+'_L90.txt', 'wb')
        fileObj90.append(fobj90)
    fileLeqx, fileL10x, fileL90x= [], [], []
    for foiID in foiIDs:
        fLeqx = open(mapObj.fmap+'\\'+'Leqx-foi'+str(foiID)+'.txt', 'wb')
        fL10x = open(mapObj.fmap+'\\'+'L10x-foi'+str(foiID)+'.txt', 'wb')
        fL90x = open(mapObj.fmap+'\\'+'L90x-foi'+str(foiID)+'.txt', 'wb') 
        fileLeqx.append(fLeqx)
        fileL10x.append(fL10x)
        fileL90x.append(fL90x)
    fileLeqm, fileL10m, fileL90m, fileSigma = [], [], [], []
    for foiID in foiIDs:
        fLeqm = open(mapObj.fmap+'\\'+'Leqm-foi'+str(foiID)+'.txt', 'wb')
        fL10m = open(mapObj.fmap+'\\'+'L10m-foi'+str(foiID)+'.txt', 'wb')
        fL90m = open(mapObj.fmap+'\\'+'L90m-foi'+str(foiID)+'.txt', 'wb') 
        fsigma = open(mapObj.fmap+'\\'+'sigma-foi'+str(foiID)+'.txt', 'wb') 
        fileLeqm.append(fLeqm)
        fileL10m.append(fL10m)
        fileL90m.append(fL90m)
        fileSigma.append(fsigma)
    
    # load local data    
    import pickle
    flvl = open('levels_saved.pkl', 'r')
    ftstamp = open('timestamp_saved.pkl', 'r')
    ffoi = open('validFoiIDs_saved.pkl', 'r')
    fncn = open('ncn_saved.pkl', 'r')
    
    levels_saved = pickle.load(flvl)
    timestamp_saved = pickle.load(ftstamp)
    validFoiIDs_saved = pickle.load(ffoi)
    ncn_saved = pickle.load(fncn)
    
    flvl.close()
    ftstamp.close()
    ffoi.close()
    fncn.close()   
    
    
    def run_local(self):
        
#    def performTask(self, task):
        """ calculate the coefficiets for dynamic map """         
#        obsID = None
        # load measurement
#        wdb = database.WarehouseDB(connect.connStr.conStrings['dhw'])
##        wdb2 = database.WarehouseDB(connect.connStr.conStrings['dhw_therion'])
#        timestamp = task.measurement_time
#        levelsEq, levels10, levels90, ncn = {}, {}, {}, {}
#        for attempt in range(3):
#            print 'getting data from warehouse, attempt:', attempt 
#            # three attempts to get valid data, each waiting 30 seconds
#            validFoiIDs = []
#            for foiID in self.foiIDs:
#                dataLeq = wdb.observations(foiID = foiID, 
#                                          phenID = 31, #Leq, 
#                                          offID = task.offering_id, 
#                                          t1 = task.measurement_time, 
#                                          t2 = task.measurement_time+datetime.timedelta(seconds = task.periodicity), 
#                                          procID = self.PROCID, 
#                                          field = 'text',
#                                          quality = None)
#                dataL10 = wdb.observations(foiID = foiID, 
#                                          phenID = 76, #task.phenomenon_id, 
#                                          offID = task.offering_id, 
#                                          t1 = task.measurement_time, 
#                                          t2 = task.measurement_time+datetime.timedelta(seconds = task.periodicity), 
#                                          procID = 2, 
#                                          field = 'numeric',
#                                          quality = None)
#                dataL90 = wdb.observations(foiID = foiID, 
#                                          phenID = 72, #task.phenomenon_id, 
#                                          offID = task.offering_id, 
#                                          t1 = task.measurement_time, 
#                                          t2 = task.measurement_time+datetime.timedelta(seconds = task.periodicity), 
#                                          procID = 2, 
#                                          field = 'numeric',
#                                          quality = None)
#                                          
#                datancn = wdb.observations(foiID = foiID, 
#                                          phenID = 350, #task.phenomenon_id, 
#                                          offID = task.offering_id, 
#                                          t1 = task.measurement_time, 
#                                          t2 = task.measurement_time+datetime.timedelta(seconds = task.periodicity), 
#                                          procID = 201, 
#                                          field = 'numeric',
#                                          quality = None)
#                    
#                if dataLeq and dataL10 and dataL90 and datancn: # if data are not empty
#                    validFoiIDs.append(foiID)
#                    timestamp = dataLeq[0][1]                
#                    level = map(float, dataLeq[0][2].split(','))
#                    octaveLevel = []
#                    for n in xrange(8): # 8 central frequencies
#                        esum = sum(10.**(0.1*numpy.asarray([level[5+n*3-1], level[5+n*3], level[5+n*3+1]])))
#                        octaveLevel.append(10.*numpy.log10(esum))
#                    levelsEq.update({foiID:numpy.array(octaveLevel)})
#                    levels10.update({foiID:float(dataL10[0][2])})
#                    levels90.update({foiID:float(dataL90[0][2])})
#                    ncn.update({foiID:int(datancn[0][2])})
#                        
#            if len(validFoiIDs)> 0:
#                print 'valid data after attempt', attempt
#                break
#            else:
#                time.sleep(5)               
#       
#        levels = [levelsEq, levels10, levels90]     
        
        for levels, validFoiIDs, ncn, timestamp in zip(self.levels_saved, self.validFoiIDs_saved, self.ncn_saved, self.timestamp_saved): 
            
#            print #->oooooooooooooooooooooooooooooooooooooooooooooooooo
#            # 206 and 220 are used for checking ONLY. not used for fitting            
#            print timestamp
#            if 206 in validFoiIDs:
#                validFoiIDs.remove(206)
#                foi206 = [206]
#            else:
#                foi206 = []                
#            if 220 in validFoiIDs:
#                validFoiIDs.remove(220)                
#                foi220= [220]
#            else:
#                foi220 = []
#                
##            if 202 in validFoiIDs:
##                validFoiIDs.remove(202)
#                
#            validFoisIDs_err2 = validFoiIDs+foi206+foi220
            
            validFoisIDs_err2 = validFoiIDs
                
            print 'validFoiIDs -> ', validFoiIDs 
#            
            if timehour_to_den(timestamp.hour)==0:
                epsilon = self.epsilonDay
                delta = self.deltaDay
                epsilonSpec = self.epsilonSpecDay
            if timehour_to_den(timestamp.hour)==1:
                epsilon = self.epsilonEve
                delta = self.deltaEve
                epsilonSpec = self.epsilonSpecEve
            if timehour_to_den(timestamp.hour)==2:
                epsilon = self.epsilonNig
                delta = self.deltaNig
                epsilonSpec = self.epsilonSpecNig

                
            # add temporary sources and correct the measurement                 
    #            print 'levelsEq[foi] -> ', levelsEq[validFoiIDs[0]]            
            
            
            [epsilonUpdate, deltaUpdate, epsilonSpecUpdate] = \
            self.mapObj._runLMS(validFoiIDs, levels, epsilon, delta, epsilonSpec, timestamp.hour, self.firstRun, timestamp, ncn)
    #        [epsilonUpdate, deltaUpdate, epsilonSpecUpdate] = limit_over_weight(epsilonUpdate, deltaUpdate, epsilonSpecUpdate, self.mapObj.Cn, len(self.mapObj.propa_j))
            
    #        print '10log industryContr -> ', 10*numpy.log10(industrContrFoi[validFoiIDs[0]])
            
            if timehour_to_den(timestamp.hour)==0:
                self.epsilonDay = epsilonUpdate
                self.deltaDay = deltaUpdate
                self.epsilonSpecDay = epsilonSpecUpdate
            if timehour_to_den(timestamp.hour)==1:
                self.epsilonEve = epsilonUpdate
                self.deltaEve = deltaUpdate
                self.epsilonSpecEve = epsilonSpecUpdate
            if timehour_to_den(timestamp.hour)==2:
                self.epsilonNig = epsilonUpdate
                self.deltaNig = deltaUpdate
                self.epsilonSpecNig = epsilonSpecUpdate
                
                
    
            if validFoiIDs:
                self.firstRun = 'no'
            
            # 206 and 220 are included when calculate  the squared errors
            [squaredErrLeq, squaredErrL10, squaredErrL90] = calc_squared_err(self.extObjs, epsilonUpdate, deltaUpdate, epsilonSpecUpdate,\
                        validFoisIDs_err2, timestamp.hour, levels, timestamp, ncn, \
                         self.mapObj.Lx_Laeq_srcCateg, self.REF)
            [Leqx, L10x, L90x , Leqm, L10m, L90m, sigmas] = exp_Lx_Lmeas(self.extObjs, epsilonUpdate, deltaUpdate, epsilonSpecUpdate,\
                        validFoisIDs_err2, timestamp.hour, levels, timestamp,ncn,\
                        self.mapObj.Lx_Laeq_srcCateg)
    #            print 'squaredErrLeq -> ', squaredErrLeq[validFoiIDs[0]]
             
            # store result in warehousedb
            strheadEps =str(timestamp)+' '+ str(self.mapObj.Cnx) +' '
            strheadDelta =str(timestamp)+' '+ str(self.mapObj.Cnx)+' '+str(len(self.mapObj.propa_j))+' '
            streps = ' '.join(str(v) for v in epsilonUpdate.reshape([1, -1])[0])
            strdelta = ' '.join(str(v) for v in numpy.asarray(deltaUpdate).reshape([1,-1])[0])            
            strEpsilonSpec = ' '.join(str(v) for v in numpy.asarray(epsilonSpecUpdate).reshape([1,-1])[0])     
    #            for p, phenID in enumerate([201311, 201312, 201313]): # 201301 is epsilon, 201302 is delta, 201303 is mu and beta     
    #                obsID = wdb2.addObservation(timeStamp = task.measurement_time, 
    #                                           foiID = task.location_id,
    #                                           procID = self.agentID, 
    #                                           phenID = phenID, 
    #                                           offID = task.offering_id,
    #                                           value = strValues[p])
            
            # write to local file
            for foi in validFoisIDs_err2:
                self.fileObjs[self.foiIDs.index(foi)].write(str(timestamp)+' '+' '.join(str(v) for v in squaredErrLeq[foi])+'\r\n')
                self.fileObj10[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(squaredErrL10[foi])+'\r\n')
                self.fileObj90[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(squaredErrL90[foi])+'\r\n')
                
                self.fileLeqx[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(Leqx[foi]) + '\r\n')
                self.fileL10x[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(L10x[foi]) + '\r\n')
                self.fileL90x[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(L90x[foi]) + '\r\n')
                self.fileLeqm[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(Leqm[foi]) + '\r\n')
                self.fileL10m[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(L10m[foi]) + '\r\n')
                self.fileL90m[self.foiIDs.index(foi)].write(str(timestamp)+' '+str(L90m[foi]) + '\r\n')
                self.fileSigma[self.foiIDs.index(foi)].write(str(sigmas[foi]) + '\r\n')
                
    #            print 'epsilon[3] -> ', epsilonUpdate[3]
    #            print 'delta[3][0] -> ', deltaUpdate[3][0]
            t = timehour_to_den(timestamp.hour)
            self.fileEps[t].write(strheadEps+streps+'\r\n')
            self.fileDelta[t].write(strheadDelta+strdelta+'\r\n')
            self.fileEpsSpec[t].write(strheadEps+strEpsilonSpec+'\r\n')
            
            
            print ('\n')
            print # ->oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

        
#        if obsID != None:
#            return 'complete'
#        else:
#            return 'complete'
#    
if __name__ == '__main__': 
#    agents.runAgent(MapLMSagent) # change to name of agent
    localObj = LocaMapLMSagent()
    localObj.run_local()
