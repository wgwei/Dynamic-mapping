# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 17:42:05 2013
Using LMS interpolation to smooth the noise map
@author: W.Wei
"""
import sys, os
import numpy
import pickle
import dbflib
import logging
# append the core path
sys.path.append(r'U:\SSM-Project\core\knowledge-extraction\python')
here = os.path.dirname(os.path.realpath(__file__))
sys.path.append(here + "\\..\\..\\..\\..\\core\\knowledge-extraction\\python")
import time, datetime

def timehour_to_den(timeHour):
    if timeHour in  [7,8,9,10,11,12,13,14,15,16,17, 18]:
        t = 0
    if timeHour in [19, 20,21, 22]:
        t = 1
    if timeHour in  [23, 0, 1, 2, 3, 4, 5, 6]:
        t = 2
    return t
        
def calc_L10(Leqa, Leqb, sigmaa, sigmab):
    Lx10 = 10*numpy.log10(10**(0.1*Leqa)+10**(0.1*Leqb)+1.28*numpy.sqrt((sigmaa**2+sigmab**2)))
    return Lx10
   
def calc_L90(Leqa, Leqb, sigmaa, sigmab):
    Lx90 = 10*numpy.log10(abs(10**(0.1*Leqa)+10**(0.1*Leqb)-1.28*numpy.sqrt((sigmaa**2+sigmab**2))))
    return Lx90
    
def write_log(logFileName, objectName):
    ''' only the self.arg will be logged.
        objectName is a class object
    '''
    logging.basicConfig(filename=logFileName, level=logging.DEBUG)    
    argins = vars(objectName)
    for item in argins.items():
        if 'ATTM' not in item[0] and 'ID' not in item[0]: 
            logging.info(item)
    logging.disable(logging.CRITICAL)

def get_Lx_Leq(ncn, dist2road, c, Lx_Laeq_srcCateg):
    ''' ncn is noticed noise events. Usually notices passing by
        dist2road is the distance from receiver to the closest road
        c is the source category from 0 to 3
        Lx_Laeq_srcCateg is the index for the second column L10-Leq and the third
        column L90-Leq. The first colum is use to find the index, which is calcuated
        by ncn+dist2road+1000
    '''
    indexValue = Lx_Laeq_srcCateg[c][:,0]
    qry = ncn+dist2road+1000.
    try: 
        idx = numpy.where(indexValue>qry)[0][0] # numpy.where return (array([2, 3, 10..], dtype=..)
        L10_Leq = Lx_Laeq_srcCateg[c][idx, 1]
        L90_Leq = Lx_Laeq_srcCateg[c][idx, 2] 
    except: # if the 
        L10_Leq = 0.
        L90_Leq = 0.
    if L10_Leq>5:
        L10_Leq = 5.
    if L10_Leq<=0:
        L10_Leq = 0.
    if L90_Leq>=0:
        L90_Leq = 0.
    if L90_Leq<-5:
        L90_Leq = -5.
    return [L10_Leq, L90_Leq]

def sigma_traffic(L10_Leq, L90_Leq, Leq):
    if L10_Leq<=0:
        L10_Leq = 0.1
    sigmai = abs(10**(0.1*(Leq+L10_Leq)) - 10**(0.1*Leq))/1.28
    return sigmai

def sigma_traffic90(L10_Leq, L90_Leq, Leq):
    if L90_Leq>=0:
        L90_Leq = -0.1
    sigma2 = abs( 10**(0.1*Leq) - 10**(0.1*(Leq+L90_Leq)))/1.28
    return sigma2

def calc_levelCategDict(foiIDs, extObjs, epsilonSpec, crct, epsilon, delta, t, C_f_c_Num, foi_categATTMdict,  Cn=4, propa_j=['direct', 'refl', 'scat']):
    ''' calculate the level caused by diffrenret source category separately
        return levelCategDict[sourceCateg][foiID]    
    '''
    
    levelCategDict = []
    for c in xrange(Cn):
        level = []  # 10log sum_{i}^{N_i} sum_{j}^{N_j}10^{0.1[L^{'}_{w,f,i}(t)-A^{'}_{f,i,j}(p_{meas})]}
        for foiID in foiIDs:   
            ATTMcj = numpy.zeros(8)            
            for j in xrange(len(propa_j)):
                ATTMcj += 10.**(0.1*(C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c])) * 10.**(-0.1*delta[c][j]) * foi_categATTMdict[foiID][c][j][0]
            level.append(10*numpy.log10(sum(ATTMcj)))        
        levelDict = {foiIDs[m]:level[m] for m in xrange(len(foiIDs))} # make a dict: foi ->  LxMLmOverDenom
        levelCategDict.append(levelDict)
    if len(extObjs)>0:
        level = []
        for cx, exts in enumerate(extObjs):
            for fioID in foiIDs:
                ATTMcj = numpy.zeros(8)  
                for j in xrange(len(propa_j)):
                    ATTMcj += 10.**(0.1*(exts.extraSourceSpec+epsilon[Cn+cx])) * 10.**(-0.1*delta[Cn+cx][j]) * exts.foi_ATTMdict_extraSource[foiID][j][0]
                level.append(10*numpy.log10(sum(ATTMcj)))
            levelDict = {foiIDs[m]:level[m] for m in xrange(len(foiIDs))}
            levelCategDict.append(levelDict)
    return levelCategDict
    
def calc_denominatorsDict(foiIDs, Cn, extObjs, epsilonSpec, crct, epsilon, delta, t, C_f_c_Num, foi_categATTMdict, propa_j):
    denominators = []  # sum_{i}^{N_i} sum_{j}^{N_j}10^{0.1[L^{'}_{w,f,i}(t)-A^{'}_{f,i,j}(p_{meas})]}
    for foiID in foiIDs:   
        ATTMcj = numpy.zeros(8)
        for c in xrange(Cn):
            for j in xrange(len(propa_j)):
                ATTMcj += 10.**(0.1*(C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c])) * 10.**(-0.1*delta[c][j]) * foi_categATTMdict[foiID][c][j][0]
        if len(extObjs)>0:
            for cx, exts in enumerate(extObjs):
                for j in xrange(len(propa_j)):
                    ATTMcj += 10.**(0.1*(exts.extraSourceSpec+epsilon[Cn+cx])) * 10.**(-0.1*delta[Cn+cx][j]) * exts.foi_ATTMdict_extraSource[foiID][j][0]
        denominators.append(ATTMcj)     
    denominatorsDict = {foiIDs[m]:denominators[m] for m in xrange(len(foiIDs))} # make a dict: foi ->  LxMLmOverDenom
    return denominatorsDict
    
def calc_level_categPropa_Dict(foiIDs, extObjs, epsilonSpec, crct, epsilon, delta, t, C_f_c_Num, foi_categATTMdict,  Cn=4, propa_j=['direct', 'refl', 'scat']):
    ''' calculate the level of every category for every propagation path '''
    levelCPdict = {}   
    for foiID in foiIDs:
        level = [[] for _ in range(Cn+len(extObjs))]
        for c in xrange(Cn):           
            for j in xrange(len(propa_j)):
                ATTMcj = 10.**(0.1*(C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c])) * 10.**(-0.1*delta[c][j]) * foi_categATTMdict[foiID][c][j][0]
                level[c].append(10*numpy.log10(sum(ATTMcj)))        
        if len(extObjs)>0:
            for cx, exts in enumerate(extObjs):
                for j in xrange(len(propa_j)):
                    ATTMcj = 10.**(0.1*((exts.extraSourceSpec+epsilon[Cn+cx]))) * 10.**(-0.1*delta[Cn+cx][j]) * exts.foi_ATTMdict_extraSource[foiID][j][0]
                    level[Cn+cx].append(10*numpy.log10(sum(ATTMcj)))
        levelCPdict.update({foiID:level})
    return levelCPdict

def spec_correction(Cn, extObjs, C_f_c_Num, epsilonSpec):
    cspec = numpy.zeros([Cn+len(extObjs), 3]) # category 1 day, evening and night
    for c in range(Cn):
        for tt in range(3):
            newlvl = 10.*numpy.log10(sum(10.**(0.1*(C_f_c_Num[c][0][tt]+epsilonSpec[c]))))
            cspec[c,tt] = newlvl - C_f_c_Num[c][1][tt]
    if len(extObjs)>0:
        for cx, exts in enumerate(extObjs):
            for tt in range(3):
                newlvl = 10.*numpy.log10(sum(10.**(0.1*(exts.extraSourceSpec+epsilonSpec[c]))))
                cspec[Cn+cx,tt] = newlvl - exts.extraSourceTotl
    return cspec
        
def srcCateg_contribution(foiIDs, extObjs, C_f_c_Num, foi_categATTMdict, epsilon,\
                         delta,epsilonSpec, timeHour, timestamp,  Cn=4, propa_j=['direct', 'refl', 'scat'], freqn=8):
    mobj = Map_Feature()
    corrections = mobj.corrections    
    wd =  datetime.date(timestamp.year, timestamp.month, timestamp.day).weekday()
    crct = corrections[wd][timeHour] # correction for a specific day and hour
    t = timehour_to_den(timeHour)                    

    crct = 0.00    
    
    denominatorsDict = calc_denominatorsDict(foiIDs, Cn, extObjs, epsilonSpec, crct, epsilon, delta, t, C_f_c_Num, foi_categATTMdict, propa_j)
    
    # calculate the level by seperate source category
    levelCPdict = calc_level_categPropa_Dict(foiIDs, extObjs, epsilonSpec, crct, epsilon, delta, t, C_f_c_Num, foi_categATTMdict,  Cn, propa_j)    
            
    return [denominatorsDict, levelCPdict]

def exp_Lx_Lmeas(extObjs, epsilon, delta, epsilonSpec, validFois, timeHour, Lxfmeaspt, timestamp, ncn, Lx_Laeq_srcCateg):
    ''' status: approved''' 
    LeqMeasPt, L10MeasPt, L90MeasPt = Lxfmeaspt[0], Lxfmeaspt[1], Lxfmeaspt[2]
    # load source 
    f1 = open('foi_categATTMdict.pkl', 'r') # precalculated attenuation
    f2 = open('C_f_c_Num.pkl', 'r') # source info         
    [foi_categATTMdict, foiIDs] = pickle.load(f1)
    C_f_c_Num = pickle.load(f2)        
    f1.close()
    f2.close()
            
    [denominatorsDict, levelCPdict] =\
    srcCateg_contribution(foiIDs, extObjs, C_f_c_Num, foi_categATTMdict, epsilon, delta, epsilonSpec, timeHour, timestamp)
    squaredErrLeq, squaredErrL10, squaredErrL90 = {}, {}, {}
    Leqx, L10x, L90x = {}, {}, {} # calculated Leq, L10 and L90
    Leqm, L10m, L90m = {}, {}, {} # measured Leq, L10 and L90 
    sigmas = {}
    for foi in validFois:       
        # squared error of Leq
#        print 'foi: ', foi
        Leqcalc = 10*numpy.log10(sum(denominatorsDict[foi]))
        Leqmeas = 10*numpy.log10(sum(10.**(0.1*LeqMeasPt[foi])))
        Leqx.update({foi : Leqcalc})
        Leqm.update({foi : Leqmeas})        
#        print 'Leq, calc -> ', Leqcalc
    
        mobj = Map_Feature()        
        [L10_Leq, L90_Leq] = get_Lx_Leq(ncn[foi], mobj.dist2roadFoi[foi], 0, Lx_Laeq_srcCateg)       
        
        # calculate the total standard deviation (sigma in normal distribution)        
#        sigmaa = sigma_traffic(L10_Leq, L90_Leq, levelCPdict[foi][0][0]) # here actually ca is forced to 0. only category 1 affect these fois                
        sigmaa = sigma_traffic(L10_Leq, L90_Leq ,Leqcalc)
        sigmas.update({foi:sigmaa})
        Lx10 = calc_L10(10*numpy.log10(sum(denominatorsDict[foi])), 0.0, sigmaa, sigmab=1)        
        L10x.update({foi : Lx10})
        L10m.update({foi : L10MeasPt[foi]})
#        print 'L10 calc -> ', Lx10
        
        # squared error of L90
        Lx90 = calc_L90(10*numpy.log10(sum(denominatorsDict[foi])), 0.0, sigmaa, sigmab=1)
        L90x.update({foi : Lx90})
        L90m.update({foi : L90MeasPt[foi]})
#        print 'L90 calc -> ', Lx90
    return [ Leqx, L10x, L90x , Leqm, L10m, L90m, sigmas]
    
    
#====================================================================
def calc_squared_err(extObjs, epsilon, delta, epsilonSpec, validFois, timeHour, Lxfmeaspt, timestamp, ncn, Lx_Laeq_srcCateg, REF):
    ''' status: approved'''
    LeqMeasPt, L10MeasPt, L90MeasPt = Lxfmeaspt[0], Lxfmeaspt[1], Lxfmeaspt[2]
    # load source 
    f1 = open('foi_categATTMdict.pkl', 'r') # precalculated attenuation
    f2 = open('C_f_c_Num.pkl', 'r') # source info         
    [foi_categATTMdict, foiIDs] = pickle.load(f1)
    C_f_c_Num = pickle.load(f2)        
    f1.close()
    f2.close()
    DIFF = {202:[1.5378, -3.2952], 203:[1.0887, -2.1591], 206:[1.4987, -3.5869], 208:[1.5639, -2.7063],\
            209:[1.4537, -2.2128], 211:[1.5884, -5.6491], 212:[1.3568, -1.9183], 220:[2.1461, -5.203]}
    t = timehour_to_den(timeHour)
    
    [denominatorsDict, levelCPdict] = \
    srcCateg_contribution(foiIDs, extObjs, C_f_c_Num, foi_categATTMdict, epsilon, delta, epsilonSpec, timeHour, timestamp)
    
    squaredErrLeq, squaredErrL10, squaredErrL90 = {}, {}, {}
    for foi in validFois:       
        # squared error of Leq
        print 'foi: ', foi
        Laeq_err2 = (10*numpy.log10(sum(denominatorsDict[foi])) - 10*numpy.log10(sum(10.**(0.1*LeqMeasPt[foi]))))**2.
        freq_err2 = (10*numpy.log10(denominatorsDict[foi]) - LeqMeasPt[foi])**2.
        err2 = list(freq_err2) + [Laeq_err2]
        squaredErrLeq.update( {foi : err2})
        print 'Leq, Laeq_err2 -> ', 10*numpy.log10(sum(denominatorsDict[foi])), Laeq_err2
        
        mobj = Map_Feature()
        # since in the Katendrach district, the affected source category to the fois
        # are only category 1, so the category is forced to 1. i.e. index 0
        [L10_Leq, L90_Leq] = get_Lx_Leq(ncn[foi], mobj.dist2roadFoi[foi], 0, Lx_Laeq_srcCateg)       
        
        # calculate the total standard deviation (sigma in normal distribution)        
        sigmaa = sigma_traffic(L10_Leq, L90_Leq, levelCPdict[foi][0][0]) # here actually ca is forced to 0. only category 1 affect these fois
#        sigmaa = sigma_traffic(L10_Leq, L90_Leq ,10*numpy.log10(sum(denominatorsDict[foi])))
        Lx10 = calc_L10(10*numpy.log10(sum(denominatorsDict[foi])), 0.0, sigmaa, sigmab=1)
        print Lx10
        L10_err2 = (Lx10 - L10MeasPt[foi])**2.
        squaredErrL10.update( {foi : L10_err2})
        print 'L10_err2 -> ', L10_err2
        if Lx10>100:
            assert(False)
        # squared error of L90
        Lx90 = calc_L90(10*numpy.log10(sum(denominatorsDict[foi])), 0.0, sigmaa, sigmab=1)
        L90_err2 = (Lx90 - L90MeasPt[foi])**2.
        squaredErrL90.update({foi : L90_err2})
        print 'L90_err2 -> ', L90_err2
    return [squaredErrLeq, squaredErrL10, squaredErrL90]
    
#====================================================================
#==================  preparation  ===================================
class FOIfeature():
    def __init__(self):
        ''' will load foi_categATTMdict.pkl and C_f_c_Num.pkl
        
            foi_categATTMdict.pkl contains the ATTM of all the source category 
            and all the propagation manner for spectrum and total
            # data structrure for foi_categATTMdict:
            # foi_categATTMdict->[{foi:[A], foi:[],...}, [foi1, foi2, ...]];
            # A-> [B],[],[],[]; B->  [C],[],[]]; C-> [D],E
            # B is source categories, C is propagation manner, D is ATTM of spectrum
            # E is ATTM of total
            # self.foi_categATTMdict[0][foiID][c][j][0] foiID-> source caegory c
              -> propagation manner  j -> spectrum            

            C_f_c_Num contains the averaged source power for every category, 
            the total La and the number of sources in each category
            # data structure of C_f_c_Num
            C_f_c_Num-> {0:[A], 1:[], ...} 0-> source category; A->[array_spec(a list for 24h), total_La(a list for 24h), num_sources]
            C_f_c_Num[c][0][0] will be: for category c -> array spec -> the spectrum of day
            C_f_c_Num[c][1][1] will be: for category c -> total Leq -> the level of evening
            C_f_c_Num[c][2] -> the category c -> the source amount           
            
            status: approved
        '''
        f1 = open('foi_categATTMdict.pkl', 'r') # precalculated attenuation
        f2 = open('C_f_c_Num.pkl', 'r') # source info         
        [self.foi_categATTMdict, foiIDs] = pickle.load(f1)
        self.C_f_c_Num = pickle.load(f2) 
        f1.close()
        f2.close()

class Map_Feature():
    def __init__(self):
        self.Cn = 4
        self.propa_j = ['direct', 'refl', 'scat']
        self.hoursDay = [7,8,9,10,11,12,13,14,15,16,17, 18]
        self.hoursEvening = [19,20,21, 22]
        self.hoursNight = [23, 0, 1, 2, 3, 4, 5, 6] 
        # correction for every hour from Monday to Sunday
        # this patten is extracted from measurement
        self.corrections = [[-1.98, -3.40, -1.91, -2.54, -1.41, 2.68, 4.03, 0.79, 1.40, 0.55, -0.50, 2.54, -1.34, -1.92, -1.24, 0.46, -0.70, -0.57, -2.01, 0.76, 1.22, -0.90, -1.72, -2.17],\
                            [-1.53, -2.00, -1.98, -2.26, -0.43, 2.06, 3.37, 0.21, 2.38, -0.31, -0.26, -0.76, -0.56, -0.40, 0.55, -0.41, 0.79, -0.70, -2.08, 1.38, 0.09, -0.89, -1.03, -1.24],\
                            [-3.00, 0.50, 1.60, -1.45, -3.80, 0.11, 1.45, 0.10, 1.05, 0.99, 0.01, -0.69, -0.51, -0.17, -1.29, -0.09, 0.07, -0.45, 0.42, 0.56, -0.71, 0.36, -0.33, 1.39],\
                            [-0.09, -4.48, -0.53, -1.56, -0.51, 2.59, 2.94, -1.30, -0.33, -0.02, 1.29, 2.84, 2.23, 1.44, 1.36, -2.46, -3.21, -5.03, -5.55, -0.23, -1.42, 2.21, -1.74, -3.92],\
                            [-1.64, -3.37, -4.19, -4.48, -2.69, 1.54, 3.41, 0.63, 1.19, -0.45, -1.80, 1.08, -1.47, -0.40, -0.12, -0.57, 0.20, 0.95, -0.38, -0.86, -0.29, 1.36, -0.57, 2.92],\
                            [0.86, -0.82, 3.03, -1.28, -2.06, -1.48, -0.27, -1.48, 0.09, 1.30, 2.08, 1.61, 2.38, -0.21, -2.42, -1.91, -1.15, -1.83, -2.73, -1.03, -1.90, 2.70, -1.52, -0.35],\
                            [0.68, -1.51, -2.03, -3.45, -2.75, -1.95, -1.57, -5.03, -2.47, -2.84, -2.87, -1.81, -2.14, 0.04, 1.57, 2.20, 2.47, 1.92, 1.42, 0.36, 0.25, -0.56, -0.11, 5.02 ]]
#        self.corrections = numpy.zeros([7, 24])
        self.Nsecall = [{202:228.6, 203:503.7,  204:228.6,  205:368.9,  206:266.8, 208:454.5,  209:454.5,  210:368.9,  211:427.6,  212:454.5,  213:266.8, 220:526.3},\
                        {202:171.5, 203:402.1,  204:171.5,  205:284.9,  206:201.7, 208:358.2,  209:358.2,  210:284.9,  211:334.8,  212:358.2,  213:201.7, 220:422.8},\
                        {202:38.3,  203:103.1,  204:38.3,   205:67.7,   206:45.7,  208:89.2,   209:89.2,   210:67.7,   211:82.1,   212:89.2,   213:45.7,  220:110.0}]
           
        self.segms = {202:5, 203:14.,  204:5,  205:9,  206:6, 208:12.,  209:12.,  210:9,  211:11,  212:12.,  213:6, 220:15}
        
        self.flowPer15min = [[27.5, 20.7, 3.2], [146.3, 111.1, 21.3], [326.3, 219.1, 50.9], [603.1, 382.2, 94.3]] # categ1, 2, 3, 4 for day evening and night
        self.dist2roadFoi = {202:18.6, 203:9999.99, 205:18.8, 206:18.5, 208:308, 209:308, 210:23.7, 211:32.0, 212:38.0, 213:18.0, 220:9999.999}



##================== add extra source ==========================
##==============================================================
class Add_Extra_Source():
    def __init__(self, nameOfMeasFoiPKL='foi_ATTMdict_industry.pkl',\
                    nameOfAffectRcvPKL='receiver_categATTMdict_IDs_Leq_industry.pkl', sourceSpecFile='extraSourceSpec.txt', sourceNum=4):
        f3 = open(nameOfMeasFoiPKL, 'r')
        [self.foi_ATTMdict_extraSource, self.fois] = pickle.load(f3)
        f4 = open(nameOfAffectRcvPKL, 'r')
        [self.receiver_categATTMdict_extraSource, self.affectedIDs] = pickle.load(f4)             
        f3.close()
        f4.close()
        self.extraSourceSpec = numpy.loadtxt(sourceSpecFile)
        self.sourceNum = sourceNum
        self.extraSourceTotl = 10*numpy.log10(sum(10**(0.1*self.extraSourceSpec)))
           
#====================================================================
#==================  calculate map  =================================
class MapLMS(FOIfeature, Map_Feature):
    def __init__(self, fmap,  extObjs, betaEps = 0.00001, muEps = 1, betaDelta=0.00001, muDelta=1, \
                betaEpsExt=0.01, muEpsExt=2e-7, betaDeltaExt=0.0001, muDeltaExt=8e-5):
        """ test only for LAeq """
        FOIfeature.__init__(self) # load the pre-calculated data
        Map_Feature.__init__(self)
        self.Cnx = self.Cn + len(extObjs)
        self.extObjs = extObjs
        self.fmap = fmap
        self.freq = numpy.array([63, 125, 250, 500, 1000, 2000, 4000, 8000])
        self.epsilon0 = numpy.zeros(self.Cnx) 
        self.delta0 =  numpy.zeros([self.Cnx, len(self.propa_j)])
        self.epsilonSpec = numpy.zeros(8)
        self.betaEps = betaEps
        self.muEps = muEps
        self.betaDelta = betaDelta
        self.muDelta = muDelta       
        self.betaEpsExt = betaEpsExt
        self.muEpsExt = muEpsExt
        self.betaDeltaExt = betaDeltaExt
        self.muDeltaExt = muDeltaExt
        fmapf = open('receiver_categATTMdict_IDs_Leq.pkl', 'r')
        [self.receiver_categATTMdict, self.IDs] = pickle.load(fmapf)
        fmapf.close()
        
        os.mkdir(fmap)
        self.sumderiv_eps, self.sumderiv_delta = [], []
        for c in range(self.Cnx):
            self.sumderiv_eps.append(open(os.path.join(fmap, 'sumderiv_eps_c'+str(c+1)+'.txt'), 'wb'))
            self.sumderiv_delta.append(open(os.path.join(fmap, 'sumderiv_delta_c'+str(c+1)+'.txt'), 'wb'))
        
        # from Bert's simulation
        Lx_Laeq_srcCateg1 = numpy.loadtxt('Lx_Laeq_srcCateg1.txt')
        Lx_Laeq_srcCateg2 = numpy.loadtxt('Lx_Laeq_srcCateg2.txt')
        Lx_Laeq_srcCateg3 = numpy.loadtxt('Lx_Laeq_srcCateg3.txt')
        Lx_Laeq_srcCateg4 = numpy.loadtxt('Lx_Laeq_srcCateg4.txt')
        self.Lx_Laeq_srcCateg = [Lx_Laeq_srcCateg1, Lx_Laeq_srcCateg2, Lx_Laeq_srcCateg3, Lx_Laeq_srcCateg4]
    
    def calc_denominatorsDict(self, foiIDs, epsilonSpec, crct, epsilon, delta, t, ncn):
        denominators = []  # sum_{i}^{N_i} sum_{j}^{N_j}10^{0.1[L^{'}_{w,f,i}(t)-A^{'}_{f,i,j}(p_{meas})]}
        for foiID in foiIDs:   
            ATTMcj = numpy.zeros(8)
            for c in xrange(self.Cn):
                for j in xrange(len(self.propa_j)):
                    ATTMcj += 10.**(0.1*(self.C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c])) * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][0]
            if len(self.extObjs)>0:
                for cx, exts in enumerate(self.extObjs):
                    for j in xrange(len(self.propa_j)):
                        ATTMcj += 10.**(0.1*(exts.extraSourceSpec+epsilon[self.Cn+cx])) * 10.**(-0.1*delta[self.Cn+cx][j]) * exts.foi_ATTMdict_extraSource[foiID][j][0]
            denominators.append(ATTMcj)     
#            print 'double sum -> ', 10*numpy.log10(ATTMcj)
        denominatorsDict = {foiIDs[m]:denominators[m] for m in xrange(len(foiIDs))} # make a dict: foi ->  LxMLmOverDenom
        return denominatorsDict
              
    def sumDeriv_of_Lfeq_Lmeas_sq_2eps(self, foiID, c, t, epsilonSpec, denominatorsDict, LeqMeasPt, crct, epsilon, delta):
        ''' foiID is a interger number
            c -> source category from 0 to 3
            t -> is integer number 0 to 2 to clarify day everning and night
        '''
        # for Leq
        numeratorEps = numpy.asarray([0.0]*8) # sum_{j}^{N_j}\frac{1}{C_c}sum_{n=1}^{C_c}10^{0.1[L^{'}_{w,f,C_c}(t)-A^{'}_{f,n,j}(p_{meas})]} # 24 hours
        if len(self.extObjs)>0 and c>=self.Cn:
            for j in xrange(len(self.propa_j)):                
                numeratorEps += 10.**(0.1*(self.extObjs[c-self.Cn].extraSourceSpec+epsilonSpec[c]+epsilon[c])) \
                * 10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][0]
        else:
            for j in xrange(len(self.propa_j)):
                numeratorEps +=  10.**(0.1*(self.C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c]))\
                * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][0] # bug detected cspec[c,t] -> epsilonSpec[c]
                
        Lxfpt = 10.*numpy.log10(denominatorsDict[foiID])          
        Lxpt = 10.*numpy.log10(sum(denominatorsDict[foiID])) # total value used for L10 and L90 calculation
        sumDerivLeq = 1./8.*sum((Lxfpt-LeqMeasPt[foiID])*numeratorEps/denominatorsDict[foiID])
        sumDerivLeq2nd = (Lxfpt-LeqMeasPt[foiID])*numeratorEps/denominatorsDict[foiID]
        return [Lxpt, sumDerivLeq, sumDerivLeq2nd]
    
    def sumDeriv_of_Lfeq_Lmeas_sq_2delta(self, foiID, c, j, t, epsilonSpec, denominatorsDict, LeqMeasPt, crct, epsilon, delta):
        if len(self.extObjs)>0 and c>=self.Cn:
            numeratorDelta = -10.**(0.1*(self.extObjs[c-self.Cn].extraSourceSpec+epsilonSpec[c]+epsilon[c]))\
             * 10.**(-0.1*delta[c][j]) *self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][0]
        else:
            numeratorDelta = -10.**(0.1*(self.C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c]))\
            * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][0]
        Lxfpt = 10.*numpy.log10(denominatorsDict[foiID])
        sumDerivLeq = 1./8.*sum((Lxfpt-LeqMeasPt[foiID])*numeratorDelta/denominatorsDict[foiID])
        return sumDerivLeq
    
    def sumDeriv_of_Leq_Lmeas_sq_2eps(self, foiID, c, t, cspec, denominatorsDict, LeqMeasPt, crct, epsilon, delta):
        # for Leq
        numeratorEps = 0.0 # sum_{j}^{N_j}\frac{1}{C_c}sum_{n=1}^{C_c}10^{0.1[L^{'}_{w,f,C_c}(t)-A^{'}_{f,n,j}(p_{meas})]} # 24 hours
        if len(self.extObjs)>0 and c>=self.Cn:
            for j in xrange(len(self.propa_j)):
                numeratorEps +=  (self.extObjs[c-self.Cn].extraSourceTotl+cspec[c,t]) \
                *10.**(0.1*((self.extObjs[c-self.Cn].extraSourceTotl*+cspec[c,t])*(1.+epsilon[c]))) \
                * 10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][1] 
        else:
            for j in xrange(len(self.propa_j)):
                numeratorEps +=  (self.C_f_c_Num[c][1][t]+cspec[c,t]+crct) \
                *10.**(0.1*(self.C_f_c_Num[c][1][t]+cspec[c,t]+crct+epsilon[c])) \
                * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
        Lxpt = 10.*numpy.log10(sum(denominatorsDict[foiID])) # total value used for L10 and L90 calculation
        LeqMeasPtTotal = 10.*numpy.log10(sum(10**(0.1*LeqMeasPt[foiID])))
        sumDerivLeq = (Lxpt-LeqMeasPtTotal)*numeratorEps/(10**(0.1*Lxpt))
        return [Lxpt, sumDerivLeq]
    
    def deriv_of_Leq_2eps(self, foiID, c, t, cspec, denominatorsDict, crct, epsilon, delta):
        # for Leq
        numeratorEps = 0.0 # sum_{j}^{N_j}\frac{1}{C_c}sum_{n=1}^{C_c}10^{0.1[L^{'}_{w,f,C_c}(t)-A^{'}_{f,n,j}(p_{meas})]} # 24 hours
        if len(self.extObjs)>0 and c>=self.Cn:
            for j in xrange(len(self.propa_j)):
                numeratorEps += 10.**(0.1*(self.extObjs[c-self.Cn].extraSourceTotl+cspec[c,t]+epsilon[c])) \
                * 10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][1] 
        else:
            for j in xrange(len(self.propa_j)):
                numeratorEps += 10.**(0.1*(self.C_f_c_Num[c][1][t]+cspec[c,t]+crct+epsilon[c])) * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
        Lxpt = 10.*numpy.log10(sum(denominatorsDict[foiID])) # total value used for L10 and L90 calculation
        derivLeq = numeratorEps/(10**(0.1*Lxpt))
        return [Lxpt, derivLeq]
    
    def deriv_of_Leq_c_2eps(self, foiID, c, t, cspec, levelCategDict, crct, epsilon, delta):
        # for Leq
        numeratorEps = 0.0 # sum_{j}^{N_j}\frac{1}{C_c}sum_{n=1}^{C_c}10^{0.1[L^{'}_{w,f,C_c}(t)-A^{'}_{f,n,j}(p_{meas})]} # 24 hours
        if len(self.extObjs)>0 and c>=self.Cn:
            for j in xrange(len(self.propa_j)):
                numeratorEps +=  (self.extObjs[c-self.Cn].extraSourceTotl+cspec[c,t])\
                *10.**(0.1*(self.extObjs[c-self.Cn].extraSourceTotl*+cspec[c,t]+epsilon[c])) \
                * 10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][1] 
        else:
            for j in xrange(len(self.propa_j)):
                numeratorEps +=  (self.C_f_c_Num[c][1][t]+cspec[c,t]+crct) * \
                10.**(0.1*(self.C_f_c_Num[c][1][t]+cspec[c,t]+crct+epsilon[c]))\
                * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
        derivLeqc = numeratorEps/(10**(0.1*levelCategDict[c][foiID]))
        return derivLeqc
        
    def deriv_of_Leq_c_j_2eps(self, foiID, c, j, t, cspec, levelCPdict, crct, epsilon, delta):
        ''' this derivative comes from the partial derivative of sigma when calculating L10
            this sigma of j>0 is constant. so when j>0, the derivative is zero
        '''
        # for Leq
        if j>0:
            derivLeqcj = 0.
        else:
            if len(self.extObjs)>0 and c>=self.Cn:
                numeratorEps =  (self.extObjs[c-self.Cn].extraSourceTotl+cspec[c,t])\
                *10.**(0.1*(self.extObjs[c-self.Cn].extraSourceTotl*+cspec[c,t]+epsilon[c])) \
                * 10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][1] 
            else:
                numeratorEps =  (self.C_f_c_Num[c][1][t]+cspec[c,t]+crct) * \
                10.**(0.1*(self.C_f_c_Num[c][1][t]+cspec[c,t]+crct+epsilon[c])) \
                * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
        derivLeqcj = numeratorEps/(10**(0.1*levelCPdict[foiID][c][j]))
        return derivLeqcj
    
    def deriv_of_Leq_2delta(self, foiID, c,j, t, cspec, denominatorsDict, crct, epsilon, delta):
        if len(self.extObjs)>0 and c>=self.Cn:
            numeratorDelta = -10.**(0.1*(self.extObjs[c-self.Cn].extraSourceTotl+cspec[c,t]+epsilon[c]))\
            * 10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][1] 
        else:
            numeratorDelta = -10.**(0.1*(self.C_f_c_Num[c][1][t]+cspec[c, t]+crct+epsilon[c]))\
            * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
        derivLeq = numeratorDelta/sum(denominatorsDict[foiID])  # sum all frequency
        return derivLeq
        
    def deriv_of_Leq_c_2delta(self, foiID, c,j, t, cspec, levelCategDict, crct, epsilon, delta):
        if len(self.extObjs)>0 and c>=self.Cn:
            numeratorDelta = -10.**(0.1*(self.extObjs[c-self.Cn].extraSourceTotl+cspec[c,t]+epsilon[c]))\
            *10.**(-0.1*delta[c][j]) * self.extObjs[c-self.Cn].foi_ATTMdict_extraSource[foiID][j][1] 
        else:
            numeratorDelta = -10.**(0.1*(self.C_f_c_Num[c][1][t]+cspec[c, t]+crct+epsilon[c])) \
            * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
        derivLeqc = numeratorDelta/(10**(0.1*levelCategDict[c][foiID]))  # sum all frequency
        return derivLeqc
    
    def deriv_L10_Lmeas_sq_of_allsrc(self, sourceNum, L10, L10Meas, Lw, Leq, Leqc, Leqcj, sigmaa, sigmab, derivLeq, asswitch, delta10, delta90):
        ''' L10 is the calculated total L10
            L10Meas is the measured L10
            sigmaa is the standard deviation of all traffic sources on condition of L10
            sigmab is the standard deviation of the industry source
            derivLeq is the derivative of Leq of all traffic sources
            asswitch is a swith here. In Katendrech, only the first category is used. if category is 1, its on otherwise its off
            delta10 is the level difference between L10simu - Leq . where L10simu is simulated by Bert
            delta90 is the level difference between L90simu - Leq
        '''
        numeratorP1 = 10**(0.1*Leqc)
        numeratorp2 = sigmaa/numpy.sqrt(sigmaa**2.+sigmab**2.) * (10**(0.1*delta10)-1) * 10**(0.1*Leqcj) * asswitch
        denominator = 10**(0.1*Leq) + 1.28*numpy.sqrt(sigmaa**2.+sigmab**2.)
        return (L10-L10Meas)*(numeratorP1+numeratorp2)/denominator
    
    def deriv_L90_Lmeas_sq_of_allsrc(self, sourceNum, L90, L90Meas, Lw, Leq, Leqc, Leqcj, sigmaa, sigmab, derivLeq, asswitch, delta10, delta90):
        numeratorP1 = 10**(0.1*Leqc)
        numeratorp2 = sigmaa/numpy.sqrt(sigmaa**2.+sigmab**2.) * (1 - 10**(0.1*delta90)) * 10**(0.1*Leqcj) * asswitch
        denominator = 10**(0.1*Leq) - 1.28*numpy.sqrt(sigmaa**2.+sigmab**2.)
        return (L90-L90Meas)*(numeratorP1+numeratorp2)/denominator
    
    def deriv_L10_Lmeas_sq_of_allsrc_delta(self, sourceNum, L10, L10Meas, Lw, Leq, Leqc, Leqcj, sigmaa, sigmab, derivLeq, asswitch, delta10, delta90):
        numeratorP1 = -10**(0.1*Leqc)
        numeratorp2 = -sigmaa/numpy.sqrt(sigmaa**2.+sigmab**2.) * (10**(0.1*delta10)-1) * 10**(0.1*Leqcj) * asswitch
        denominator = 10**(0.1*Leq) + 1.28*numpy.sqrt(sigmaa**2.+sigmab**2.)
        return (L10-L10Meas)*(numeratorP1+numeratorp2)/denominator
    
    def deriv_L90_Lmeas_sq_of_allsrc_delta(self, sourceNum, L90, L90Meas, Lw, Leq, Leqc, Leqcj, sigmaa, sigmab, derivLeq, asswitch, delta10, delta90):
        numeratorP1 = -10**(0.1*Leqc)
        numeratorp2 = -sigmaa/numpy.sqrt(sigmaa**2.+sigmab**2.) * (1 - 10**(0.1*delta90)) * 10**(0.1*Leqcj) * asswitch
        denominator = 10**(0.1*Leq) - 1.28*numpy.sqrt(sigmaa**2.+sigmab**2.)
        return (L90-L90Meas)*(numeratorP1+numeratorp2)/denominator
    
    
    def update_eps(self, foiIDs, epsilonSpec, crct, cspec, epsilon, delta, t, \
            denominatorsDict, LeqMeasPt, L10MeasPt, L90MeasPt, betaEps, muEps, betaEpsExt, muEpsExt,Ti, ncn):   
        # epsilon updating
        epsilonUpdate = numpy.zeros(self.Cnx)
        epsilonSpecUpdate = numpy.zeros([self.Cnx, 8]) 
        levelCPdict = calc_level_categPropa_Dict(foiIDs, self.extObjs, epsilonSpec, crct, epsilon, delta, t, self.C_f_c_Num, self.foi_categATTMdict,  self.Cn, self.propa_j)
        levelCategDict = calc_levelCategDict(foiIDs, self.extObjs, epsilonSpec, crct, epsilon, delta, t, self.C_f_c_Num, self.foi_categATTMdict,  Cn=4, propa_j=['direct', 'refl', 'scat'])
        for c in xrange(self.Cnx): 
            if c>=self.Cn:
                muEps = muEpsExt[c-self.Cn]
                betaEps = betaEpsExt[c-self.Cn]
                sourceNum = self.extObjs[c-self.Cn].sourceNum
            else:
                sourceNum = self.C_f_c_Num[c][2]
            sumDerivLeq = 0.0
            sumDeriv = 0.0
            sumDerivL10 = 0.0
            sumDerivL90 = 0.0
            sumDerivLeq2nd = numpy.zeros(8) # 8 central frequency                 
            for foiID in foiIDs:
                # for Leq
                [LxptTot, sumDerivLeqin, sumDerivLeq2ndin] = self.sumDeriv_of_Lfeq_Lmeas_sq_2eps(foiID, c, t, epsilonSpec, denominatorsDict, LeqMeasPt, crct, epsilon, delta)
                sumDerivLeq = sumDerivLeq+sumDerivLeqin
                sumDerivLeq2nd = sumDerivLeq2nd+sumDerivLeq2ndin
#                
                # check whether all sources are constant
                # in this katendlercht case, only c=0 can contribute to L10
                const = 0
                for ct in xrange(self.Cnx-1):
                    if levelCategDict[ct][foiID]>levelCategDict[0][foiID]:
                        const = 1
                        break                       
                    
                if const==1:
                    sumDerivL10in = 0.0
                    sumDerivL90in = 0.0
                else:
                    # # calculate contribution of L10 and L90
                    if c>0:
                        L10_Leq = 0.0
                        L90_Leq = 0.0
                    else:                    
                        [L10_Leq, L90_Leq] = get_Lx_Leq(ncn[foiID], self.dist2roadFoi[foiID], 0,  self.Lx_Laeq_srcCateg) # for the updating fois, only category 1 affect L10 and L90. so here c is fixed with 0                
                    sigmaa = sigma_traffic(L10_Leq, L90_Leq, levelCPdict[foiID][0][0]) 
    #                sigmaa = sigma_traffic(L10_Leq, L90_Leq, LxptTot) #levelCategDict is  forced to category 1 and the propagation is forced to direct
                    sigmab = 1. # force it to 1 a small number for the industry source
    
                    # for L10 
                    [LxptTot, derivLeq] = self.deriv_of_Leq_2eps(foiID, c, t, cspec, denominatorsDict, crct, epsilon, delta)
    #                print 'derivLeq: ', derivLeq
                    L10 = calc_L10(LxptTot, 0, sigmaa, sigmab)
    #                print 'L10 ', L10
                    if c==0:
                        asswitch = 1
                    else:
                        asswitch = 0
                    if c>=self.Cn:
                        Lw = self.extObjs[c-self.Cn].extraSourceTotl
                    else:
                        Lw = self.C_f_c_Num[c][1][t]       
                    
                    
                    if (L10_Leq>1. and L90_Leq<-1):
                        sumDerivL10in = self.deriv_L10_Lmeas_sq_of_allsrc(sourceNum, L10, L10MeasPt[foiID], Lw, \
                                LxptTot, levelCategDict[c][foiID], levelCPdict[foiID][c][0], sigmaa, sigmab, derivLeq,\
                                asswitch, L10_Leq, L90_Leq)
                    else:
                        sumDerivL10in = 0.0                    
                        
                    # for L90
                    L90 = calc_L90(LxptTot, 0, sigmaa, sigmab)
                    if (L10_Leq>1. and L90_Leq<-1):
                        sumDerivL90in = self.deriv_L90_Lmeas_sq_of_allsrc(sourceNum, L90, L90MeasPt[foiID], Lw, \
                                LxptTot, levelCategDict[c][foiID], levelCPdict[foiID][c][0], sigmaa, sigmab, derivLeq, \
                                asswitch, L10_Leq, L90_Leq)    
                    else:
                        sumDerivL90in = 0.0
                        
                sumDerivL10 = sumDerivL10 + sumDerivL10in
                sumDerivL90 = sumDerivL90 + sumDerivL90in
            
            sumDeriv = sumDerivLeq + sumDerivL10 + sumDerivL90                
            self.sumderiv_eps[c].write(str(sumDeriv)+'\r\n')
            epsilonUpdate[c] = epsilon[c]*(1.-betaEps) - muEps * sumDeriv            
            epsilonSpecUpdate[c, :] = epsilon[c]*(1.-betaEps) - 0.5*muEps*sumDerivLeq2nd
        return [epsilonUpdate, epsilonSpecUpdate]
       
        
    def update_delta(self, foiIDs, epsilonSpec, crct, cspec, epsilon, delta, t,\
            denominatorsDict, LeqMeasPt, L10MeasPt, L90MeasPt, betaDelta, muDelta, betaDeltaExt, muDeltaExt, Ti, ncn):        
        # delta updating
        levelCPdict = calc_level_categPropa_Dict(foiIDs, self.extObjs, epsilonSpec, crct, epsilon, delta, t, self.C_f_c_Num, self.foi_categATTMdict,  self.Cn, self.propa_j)
        levelCategDict = calc_levelCategDict(foiIDs, self.extObjs, epsilonSpec, crct, epsilon, delta, t, self.C_f_c_Num, self.foi_categATTMdict,  Cn=4, propa_j=['direct', 'refl', 'scat'])
        deltaUpdate = numpy.zeros([self.Cnx, len(self.propa_j)])        
        for c in xrange(self.Cnx):        
            if c>=self.Cn:
                muDelta = muDeltaExt[c-self.Cn]
                betaDelta = betaDeltaExt[c-self.Cn]
                sourceNum = self.extObjs[c-self.Cn].sourceNum
            else:
                sourceNum = self.C_f_c_Num[c][2]
            for j in xrange(len(self.propa_j)):
                sumDeriv = 0.0                
                
                # for Leq 
                sumDerivLeq =  0.0
                for foiID in foiIDs:
                    sumDerivLeqin = self.sumDeriv_of_Lfeq_Lmeas_sq_2delta(foiID, c, j, t, epsilonSpec, denominatorsDict, LeqMeasPt, crct, epsilon, delta)
                    sumDerivLeq = sumDerivLeq+sumDerivLeqin               
                
                # for L10 and L90
                sumDerivL10 = 0.0
                sumDerivL90 = 0.0
                for foiID in foiIDs:
                    
                    # check whether all sources are constant
                    # in this katendlercht case, only c=0 can contribute to L10
                    const = 0
                    for jj in xrange(len(self.propa_j)):
                        if levelCPdict[foiID][c][jj]>levelCPdict[foiID][0][0]:
                            const = 1
                            break                            
                        
                    if const==1:
                        sumDerivL10in = 0.0
                        sumDerivL90in = 0.0
                    else:
                        LxptTot = 10.*numpy.log10(sum(denominatorsDict[foiID]))
                        if c>0:
                            L10_Leq = 0.0
                            L90_Leq = 0.0
                        else:
                            [L10_Leq, L90_Leq] = get_Lx_Leq(ncn[foiID], self.dist2roadFoi[foiID], 0, self.Lx_Laeq_srcCateg)
                        sigmaa = sigma_traffic(L10_Leq, L90_Leq, levelCPdict[foiID][0][0])
    #                    sigmaa = sigma_traffic(L10_Leq, L90_Leq, LxptTot)
                        sigmab = 1. # force it to 1
                        derivLeq = self.deriv_of_Leq_2delta(foiID, c,j, t, cspec, denominatorsDict, crct, epsilon, delta)
                        
                        # L10
                        L10 = calc_L10(LxptTot, 0, sigmaa, sigmab)
                        if c==0 and j==0:
                            asswitch = 1
                        else:
                            asswitch = 0
                        if c>=self.Cn:
                            Lw = self.extObjs[c-self.Cn].extraSourceTotl
                        else:
                            Lw = self.C_f_c_Num[c][1][t]
                        sumDerivL10in = self.deriv_L10_Lmeas_sq_of_allsrc_delta(sourceNum, L10, L10MeasPt[foiID], Lw,\
                                LxptTot, levelCategDict[c][foiID], levelCPdict[foiID][c][0], sigmaa, sigmab, derivLeq, \
                                asswitch, L10_Leq, L90_Leq)                   
                    
                        # L90
                        L90 = calc_L90(LxptTot, 0, sigmaa, sigmab)
                        sumDerivL90in = self.deriv_L90_Lmeas_sq_of_allsrc_delta(sourceNum, L90, L90MeasPt[foiID], Lw,\
                                LxptTot, levelCategDict[c][foiID], levelCPdict[foiID][c][0], sigmaa, sigmab, derivLeq,\
                                asswitch, L10_Leq, L90_Leq)  
                            
                    sumDerivL10 = sumDerivL10+sumDerivL10in                                                                               
                    sumDerivL90 = sumDerivL90+sumDerivL90in
                             
                sumDeriv = sumDerivLeq +sumDerivL10+sumDerivL90
                self.sumderiv_delta[c].write(str(sumDeriv)+' ')
                deltaUpdate[c,j] = delta[c,j]*(1.-betaDelta) - muDelta * sumDeriv
            self.sumderiv_delta[c].write('\r\n')
        return deltaUpdate
    
    def calc_nsec(self, ncn, segms, speed, Ti=900.,segLen=10.):
        ''' ncn is number of event. number of cars passing-by
            segms: number of segments involded to calculated L10
            speed: the speed limit of the roald km/h
            segLen: segment length when split line source to point sources
        '''
        secpcar = segLen/(numpy.asarray(speed)/3.6)
        nsecpseg = ncn*secpcar
        if nsecpseg >=Ti:
            nsecall=Ti
        else:
            nsecall = (1-(1-nsecpseg/Ti)**segms)*Ti 
        
        return nsecall
        
    def _getEpsilonDelta(self, foiIDs, Lxfmeaspt, epsilon, delta, epsilonSpec, betaEps, muEps, betaDelta, muDelta, betaEpsExt, muEpsExt,betaDeltaExt, muDeltaExt, timeHour, timestamp, ncn):
        """ calculate the epsilon by LMS method. 
            foiIDs is the feature_of_interest ID number in list [foi1, foi2, foi3,...]
            Lxfmeaspt is the measured sound power level at time stamp t, in 8 central frequencies
            which is a dictionary respect to the foiID. i.e. Lxfmeaspt = {foiID: [L_f63, L_f125, ...], foiID2:[], ...}
            epsilon is the epsion needed optimized in a list or array [e1, e2, ] where e1 = [e_f63, e_f125,...]
                 every element is the epsilon corresponding for every source category
            delta is the delta needed to be optimized in a list or array
                 [[delta_c1j1, delta_c1j2,...],[delta_c2j1, delta_c2j2,...],[]], where delta_cnjm = [delta_f63, delta_f125,...]
            betaEps is a small number which is  constrainer of epsilon
            muEps is a small number which is the step of the derivative of epsilon
            betaDelta is a small number which is the constrainer of delta
            muEps is a small number which is the step of the derivative of delta
            timeHour is the time stamp in hour to determine the traffic flow
            
            status: approved
        """ 
        if foiIDs:
            Ti = 900.0 # time interval
            LeqMeasPt, L10MeasPt, L90MeasPt = Lxfmeaspt[0], Lxfmeaspt[1], Lxfmeaspt[2]
            wd =  datetime.date(timestamp.year, timestamp.month, timestamp.day).weekday()
            print ('weekday: %d'  %wd)
            crct = self.corrections[wd][timeHour] # correction for a specific day and hour           
            
            crct = 0.00
                        
            t =     timehour_to_den(timeHour)  
            
            print ('corrections %0.1f' %crct)
            
            cspec = spec_correction(self.Cn, self.extObjs, self.C_f_c_Num, epsilonSpec)

            # calculate the denominator     status: approved
            denominatorsDict = calc_denominatorsDict(foiIDs, self.Cn, self.extObjs, epsilonSpec, crct, epsilon, delta, t, self.C_f_c_Num, self.foi_categATTMdict, self.propa_j)
      
            # epsilon updating
            [epsilonUpdate, epsilonSpecUpdate] = self.update_eps(foiIDs, epsilonSpec, crct, cspec, epsilon, delta, t, denominatorsDict, LeqMeasPt, L10MeasPt, L90MeasPt, betaEps, muEps, betaEpsExt, muEpsExt,Ti, ncn)
                
            # delta updatingfoi_categATTMdict
            deltaUpdate = self.update_delta(foiIDs, epsilonSpec, crct, cspec, epsilon, delta, t, denominatorsDict, LeqMeasPt, L10MeasPt, L90MeasPt, betaDelta, muDelta, betaDeltaExt, muDeltaExt, Ti, ncn)
            
        else:
            epsilonUpdate = epsilon    
            deltaUpdate = delta
            epsilonSpecUpdate = epsilonSpec
        return [epsilonUpdate, deltaUpdate, epsilonSpecUpdate]
        
    
    def transfer_str_eps_to_value_extr(self, strEps):
        ''' strEps = "2013-09-01 03:15:00+02:00 4 6.63171868177e-08 0.178829209559 0.0188434443955 -0.0210514145411, ....." 
            return an array of epsilon
        '''
        strs = strEps.strip('\n').split(' ')
        epsValue = numpy.array(map(float, strs[3::])).reshape(int(strs[2]), 8) # 63, 125, 250, 500, 1000, 2000, 4000, 8000
        return epsValue 
        
    def transfer_str_delta_to_value_extr(self, strDelta):
        ''' strdelta = '2013-09-01 01:15:00+02:00 4 3 -1.08793626823e-09 -6.07292988079e-08 -6.22135299037e-08 -0.0332977496145 -0.0333839441106 -0.0333686794129 -0.00717665740818 -0.00682591812668 -0.00614515962004 0.000589467395165 -0.00146700345391 0.0016736650452'
            return an array of delta
        '''
        strs = strDelta.strip('\n').split(' ')
        values = numpy.array(map(float, strs[4::])).reshape([int(strs[2]), -1])
        deltaValue = [values[n].reshape(int(strs[3]), 8) for n in range(int(strs[2]))]
        return deltaValue
        
    def transfer_str_eps_to_value(self, strEps):
        ''' strEps = "2013-09-01 03:15:00+02:00 4 6.63171868177e-08 0.178829209559 0.0188434443955 -0.0210514145411" 
            return an array of epsilon
        '''
        strs = strEps.strip('\n').split(' ')[3::]
        epsValue = numpy.array(map(float, strs))
        return epsValue 
        
    def transfer_str_delta_to_value(self, strDelta):
        ''' strdelta = '2013-09-01 01:15:00+02:00 4 3 -1.08793626823e-09 -6.07292988079e-08 -6.22135299037e-08 -0.0332977496145 -0.0333839441106 -0.0333686794129 -0.00717665740818 -0.00682591812668 -0.00614515962004 0.000589467395165 -0.00146700345391 0.0016736650452'
            return an array of delta
        '''
        strs = strDelta.strip('\n').split(' ')
        deltaValue = numpy.array(map(float, strs[4::])).reshape(int(strs[2]), int(strs[3]))
        return deltaValue
        
    def _runLMS(self, foiIDs, LxfmeasptByFoi, epsilon, delta, epsilonSpec, timeHour ,firstRun, timestamp, ncn):
        if firstRun=='yes':
            # writing log 
            flog = open(self.fmap+'\\'+'log' +'.txt', 'wb')
            flog.write(time.ctime()+'\r\n')
            flog.write('Source category number -> ' + str(self.Cn)+'\r\n')
            flog.write('propagation manner -> ' + str(self.propa_j)+'\r\n')
            flog.write('betaEps -> '+str(self.betaEps)+'\r\n')
            flog.write('muEps -> '+ str(self.muEps)+'\r\n')
            flog.write('betaDelta -> '+str(self.betaDelta)+'\r\n')
            flog.write('muDleta -> '+str(self.muDelta)+'\r\n')
            flog.close()
            [epsilonUpdate, deltaUpdate, epsilonSpecUpdate] = self._getEpsilonDelta(foiIDs, LxfmeasptByFoi, self.epsilon0, self.delta0, self.epsilonSpec, self.betaEps, self.muEps, self.betaDelta, self.muDelta, self.betaEpsExt, self.muEpsExt, self.betaDeltaExt, self.muDeltaExt, timeHour, timestamp, ncn)
        else:  
            [epsilonUpdate, deltaUpdate, epsilonSpecUpdate] = self._getEpsilonDelta(foiIDs, LxfmeasptByFoi, epsilon, delta, epsilonSpec, self.betaEps, self.muEps, self.betaDelta, self.muDelta, self.betaEpsExt, self.muEpsExt, self.betaDeltaExt, self.muDeltaExt, timeHour, timestamp, ncn)        
        return [epsilonUpdate, deltaUpdate, epsilonSpecUpdate] 

    def _updateMap(self, epsilon, delta, epsilonSpec, timeHour, timestamp):
        wd =  datetime.date(timestamp.year, timestamp.month, timestamp.day).weekday()
        crct = self.corrections[wd][timeHour] # correction for a specific day and hour
        
        crct = 0.00
        
        
        t =     timehour_to_den(timeHour)  
        receiverLevels = []  # sum_{i}^{N_i} sum_{j}^{N_j}10^{0.1[L^{'}_{w,f,i}(t)-A^{'}_{f,i,j}(p_{meas})]}
        for ID in self.IDs:   
            powerSum = 0.
            for c in xrange(self.Cn):
                for j in xrange(len(self.propa_j)):
                    powerSum += 10.**(0.1*(self.C_f_c_Num[c][0][t]+epsilonSpec[c]+crct+epsilon[c])) * 10.**(-0.1*delta[c][j]) * self.receiver_categATTMdict[ID][c][j][0]
            if len(self.extObjs)>0:
                for cx, exts in enumerate(self.extObjs):
                    for j in xrange(len(self.propa_j)):
                        powerSum += 10.**(0.1*(exts.extraSourceSpec+epsilonSpec[self.Cn+cx]+epsilon[self.Cn+cx])) * 10.**(-0.1*delta[self.Cn+cx][j]) * self.receiver_categATTMdict_extraSource[ID][j][0]
                    
            receiverLevels.append(10.*numpy.log10(powerSum))
        return receiverLevels    
        
    def _redrawMap(self, receiverFileName, epsilon, delta, timeHour, minute, timestamp):        
        receiverLevels = self._updateMap(epsilon, delta, timeHour, timestamp)
        totalAlvl = []
        for n, rlvl in enumerate(receiverLevels):
            receiverLevels[n] = 10.*numpy.log10(10.**(0.1*numpy.array(rlvl)))
            totalAlvl.append(10.*numpy.log10(sum(10.**(0.1*receiverLevels[n]))))
        DBF = dbflib.open(receiverFileName)
        RECN = DBF.record_count()
        newDBF = dbflib.create(receiverFileName+'_update_'+str(timeHour)+'u') #+'_'+str(minute))
        # copy the old field and record
        for f in xrange(DBF.field_count()):
            fi = DBF.field_info(f)
            newDBF.add_field(fi[1], fi[0], fi[2], fi[3])
        for r in xrange(RECN):
            rec = DBF.read_record(r)
            newDBF.write_record(r, rec)
        # add new field and attribute
        moreFd = ['L_63', 'L_125', 'L_250', 'L_500', 'L_1000', 'L_2000', 'L_4000', 'L_8000']
        [newDBF.add_field(fd, dbflib.FTDouble, 16, 1) for fd in moreFd]
        newDBF.add_field('LA', dbflib.FTDouble, 16, 1)
        for r in xrange(RECN):
            [newDBF.write_attribute(r, f+1+m, receiverLevels[r][m]) for m in range(len(moreFd))]
            newDBF.write_attribute(r, f+1+m+1, totalAlvl[r])
        DBF.close()
        newDBF.close()

        
    def _redrawMap2spatialwdb(self, receiverFileName, epsilon, delta, timeHour, extraSrcCorrect, immissionGrid_info):
        receiverLevels = self._updateMap(epsilon, delta, timeHour)
        receiverLevels = 10.*numpy.log10(10.**(0.1*numpy.array(receiverLevels))) # + extraSrcCorrect)
        
        print '------------ report receiver Levels -------------------------'
        print min(receiverLevels), sum([x for x in receiverLevels if x > 0])/len([x for x in receiverLevels if x > 0]), max(receiverLevels)
        
        """        
        DBF = dbflib.open(receiverFileName)
        RECN = DBF.record_count()
        newDBF = dbflib.create(receiverFileName+'_update_'+str(timeHour)+'u')
        # copy the old field and record
        for f in xrange(DBF.field_count()):
            fi = DBF.field_info(f)
            newDBF.add_field(fi[1], fi[0], fi[2], fi[3])
        for r in xrange(RECN):
            rec = DBF.read_record(r)
            newDBF.write_record(r, rec)
        # add new field and attribute
        newDBF.add_field('DBA_D', dbflib.FTDouble, 9, 1)
        for r in xrange(RECN):
            newDBF.write_attribute(r, f+1, receiverLevels[r])
        DBF.close()
        newDBF.close()
        """
        
                          
        table = immissionGrid_info['table2update']
        field = immissionGrid_info['field2update']
        IDfield = immissionGrid_info['ID_field']
        
        # this is dangerous
        # the ID field in the imission grid is a range() 
        # it would be nicer if the receiverLevels would be a set of (id, value) tuples !!!
        sql_store_list = []
        for id, res in enumerate(receiverLevels):
            if res> 0:
                query = 'UPDATE "%s" SET %s = %f WHERE %s = %i' % (table, field, res, IDfield, id)
                sql_store_list.append(query)
            
            
        print sql_store_list[:10]
          
        immissionGrid_info['mapdhw'].customSQL_list(sql_store_list, echo=True)
        

    
