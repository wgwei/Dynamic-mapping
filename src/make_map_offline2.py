# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 17:42:05 2013
Using LMS interpolation to smooth the noise map
offline test version
necessary to online:
    1 -> stote the procedure and phenominen ot database
    2 -> change status of first run or not
@author: W.Wei
"""
import sys, os
import numpy
import pickle
import dbflib
# append the core path
sys.path.append(r'U:\SSM-Project\core\knowledge-extraction\python')
here = os.path.dirname(os.path.realpath(__file__))
sys.path.append(here + "\\..\\..\\..\\..\\core\\knowledge-extraction\\python")
import time

sys.path.append(r"U:\svn\core\knowledge-extraction\python")

#====================================================================
#==================  preparation  ===================================
class FOIfeature():
    def __init__(self):
        ''' will load foi_categATTMdict.pkl and C_f_c_Num.pkl
        
            foi_categATTMdict.pkl contains the ATTM of all the source category 
            and all the propagation manner for spectrum and total
            # data structrure for foi_categATTMdict:
            # foi_categATTMdict->{foi:[A], foi:[],...};
            # [A]-> [[B],[],[],[]]; [B]->  [[C],[],[]]; [C]-> [[D],[E]]
            # B is source categories, C is propagation manner, D is ATTM of spectrum
            # E is ATTM of total
            # self.foi_categATTMdict[foiID][c][j][0] foiID-> source caegory c
              -> propagation manner  j -> spectrum            

            C_f_c_Num contains the averaged source power for every category, 
            the total La and the number of sources in each category
            # data structure of C_f_c_Num
            C_f_c_Num-> {0:[A], 1:[], ...} 0-> source category; [A]->[array_spec(a list for 24h), total_La(a list for 24h), num_sources]
            C_f_c_Num[c][0][0] will be: for category c -> array spec -> the spectrum of 0:00h
            C_f_c_Num[c][1][0] will be: for category c -> total Leq -> the level of 0:00h
            C_f_c_Num[c][2] -> the category c -> the source amount
            
            status: approved
        '''
        f1 = open('foi_categATTMdict.pkl', 'r') # precalculated attenuation
        f2 = open('C_f_c_Num.pkl', 'r') # source info 
        [self.foi_categATTMdict, foiIDs] = pickle.load(f1)
        self.C_f_c_Num = pickle.load(f2)           
        f1.close()
        f2.close()
    
#====================================================================
#==================  calculate map  =================================
class MapLMS(FOIfeature):
    def __init__(self, fmap,  Cn=4, propa_j = ['direct', 'refl', 'scat'], betaEps = 0.00001, muEps = 1, betaDelta=0.00001, muDelta=1):
        """ test only for LAeq """
        FOIfeature.__init__(self) # load the pre-calculated data
        self.freq = numpy.array([63., 125., 250., 500., 1000., 2000., 4000., 8000.])
        self.Cn = Cn
        self.propa_j = propa_j 
        self.epsilon0 = numpy.zeros(Cn) # 4 is soruce category, 8 is freq. [[63,125,250,...], [], [], ...]
        self.delta0 = numpy.zeros((Cn, len(propa_j)))#[[[63,125,250,500,...], [], []], [[], [], []],... ]
        self.betaEps = betaEps
        self.muEps = muEps
        self.betaDelta = betaDelta
        self.muDelta = muDelta         
        fmapf = open('receiver_categATTMdict_IDs_Leq.pkl', 'r')
        [self.receiver_categATTMdict, self.IDs] = pickle.load(fmapf)
        fmapf.close()
        self.fmap = fmap
        
    def _getEpsilonDelta(self, foiIDs, Lxfmeaspt, epsilon, delta, betaEps, muEps, betaDelta, muDelta, timeHour):
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
        """ 
        
        # calculate the denominator     status: approved
        denominators = []  # sum_{i}^{N_i} sum_{j}^{N_j}10^{0.1[L^{'}_{w,f,i}(t)-A^{'}_{f,i,j}(p_{meas})]}
        for foiID in foiIDs:   
            ATTMcj = 0.
            for c in xrange(self.Cn):
                for j in xrange(len(self.propa_j)):
                    ATTMcj += 10.**(0.1*(self.C_f_c_Num[c][1][timeHour]*(1.+epsilon[c]))) * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1] # 24 hours
            denominators.append(ATTMcj)                
        denominatorsDict = {foiIDs[m]:denominators[m] for m in xrange(len(foiIDs))} # make a dict: foi ->  LxMLmOverDenom
        
        # calculate the difference between calculated and simulated results
        # write out the error for every foiIDs
        LxfMinusLxfmeas =[]
        for foiID in foiIDs:
            LxfMinusLxfmeas.append([foiID, (10.*numpy.log10(denominatorsDict[foiID])-Lxfmeaspt[foiID])**2.])
        
        # epsilon updating
        epsilonUpdate = numpy.zeros(self.Cn)        
        for c in xrange(self.Cn): 
            sumDeriv = 0.0
            for foiID in foiIDs:
                numeratorEps = 0.0 # sum_{j}^{N_j}\frac{1}{C_c}sum_{n=1}^{C_c}10^{0.1[L^{'}_{w,f,C_c}(t)-A^{'}_{f,n,j}(p_{meas})]} # 24 hours
                for j in xrange(len(self.propa_j)):
                    numeratorEps += 1./self.C_f_c_Num[c][2]*10.**(0.1*(self.C_f_c_Num[c][1][timeHour]*(1.+epsilon[c]))) * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
                Lxfpt = 10.*numpy.log10(denominatorsDict[foiID])                
                sumDeriv += (Lxfpt-Lxfmeaspt[foiID])*numeratorEps/denominatorsDict[foiID]
#                print 'logup vs logdw ', 10.*numpy.log10(numeratorEps), ' <-> ', 10.*numpy.log10(denominatorsDict[foiID])
#                print 'ratio -> ', numeratorEps/denominatorsDict[foiID]
#                print 'Lx-Lm -> ', Lxfpt-Lxfmeaspt[foiID]
#                print 'eps sumDeriv -> ', sumDeriv
            epsilonUpdate[c] = epsilon[c]*(1.-betaEps) - muEps*(self.C_f_c_Num[c][1][timeHour]) * sumDeriv
            
        # delta updating
        deltaUpdate = numpy.zeros([self.Cn, len(self.propa_j)])
        for c in xrange(self.Cn):
            sumDeriv =0.0
            for j in xrange(len(self.propa_j)):
                for foiID in foiIDs:
                    numeratorDelta = -1./self.C_f_c_Num[c][2]*10.**(0.1*(self.C_f_c_Num[c][1][timeHour]*(1.+epsilon[c]))) * 10.**(-0.1*delta[c][j]) * self.foi_categATTMdict[foiID][c][j][1]
                    Lxfpt = 10.*numpy.log10(denominatorsDict[foiID])
                    sumDeriv += (Lxfpt-Lxfmeaspt[foiID])*numeratorDelta/denominatorsDict[foiID] 
#                    print 'delta sum deriv -> ', sumDeriv
                delta_c_j = delta[c][j]*(1.-betaDelta) - muDelta * sumDeriv
                deltaUpdate[c][j] = delta_c_j
                
        return [epsilonUpdate, deltaUpdate, LxfMinusLxfmeas]
        
    def _runLMS(self, foiIDs, LxfmeasptByFoi, epsilon, delta, timeHour ,firstRun):
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
            
            [epsilonUpdate, deltaUpdate, LxfMinusLxfmeas] = self._getEpsilonDelta(foiIDs, LxfmeasptByFoi, self.epsilon0, self.delta0, self.betaEps, self.muEps, self.betaDelta, self.muDelta, timeHour)
        else:  
            [epsilonUpdate, deltaUpdate, LxfMinusLxfmeas] = self._getEpsilonDelta(foiIDs, LxfmeasptByFoi, epsilon, delta, self.betaEps, self.muEps, self.betaDelta, self.muDelta, timeHour)        
        return [epsilonUpdate, deltaUpdate, LxfMinusLxfmeas] 

    def _updateMap(self, epsilon, delta, timeHour):
        receiverLevels = []  # sum_{i}^{N_i} sum_{j}^{N_j}10^{0.1[L^{'}_{w,f,i}(t)-A^{'}_{f,i,j}(p_{meas})]}
        for ID in self.IDs:   
            powerSum = 0.
            for c in xrange(self.Cn):
                for j in xrange(len(self.propa_j)):
                    powerSum += 10.**(0.1*(self.C_f_c_Num[c][1][timeHour]*(1.+epsilon[c]))) * 10.**(-0.1*delta[c][j]) * self.receiver_categATTMdict[ID][c][j][1]
            receiverLevels.append(10.*numpy.log10(powerSum))
        return receiverLevels
        
    def _redrawMap(self, receiverFileName, epsilon, delta, timeHour, day):
        receiverLevels = self._updateMap(epsilon, delta, timeHour)
        DBF = dbflib.open(receiverFileName)
        RECN = DBF.record_count()
        if not os.path.exists(str(day)):
            os.mkdir(str(day))
        newDBF = dbflib.create(str(day)+'\\'+receiverFileName+'_update_'+str(timeHour)+'u')
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

def get_24h_params(folder, days):
    ''' @params: days -> how many days we want to calculate
        @return: return epsilon and delta for every day. [day1, day2,...]
                day1 = [epsilon_by_hour, delta_by_hour]
    '''    
    feps = open(os.path.join(folder, 'epsilon.txt'))
    ss = feps.readlines()
    fdelta = open(os.path.join(folder, 'delta.txt'))
    sdt = fdelta.readlines()
    epsilon, delta = [], []
    for s in ss:
        sd = s.strip('\n').split(' ')
        thour = int(sd[1][0:2])
        epsilon.append([thour, map(float, sd[3::])])
    for s in sdt:        
        sd = s.strip('\n').split(' ')
        thour = int(sd[1][0:2])
        hDelt = [numpy.reshape(map(float, sd[4::]), [int(sd[2]), int(sd[3])])]
        delta.append([thour, hDelt])
    del sdt, ss # releas memory    
    n = 0
    dayEps = []
    dayDelta = []
    for day in range(days):    
        print 'day -> ', day
        hourEps = {n:[] for n in range(24)}
        hourDelta = {n:[] for n in range(24)}
        for h in range(24):
            print 'hour -> ', h
            while (23-h)!=epsilon[-1-n][0]:
                n += 1
                print 'increased n -> ', n
            hourEps[23-h] += epsilon[-1-n][1]
            hourDelta[23-h] += delta[-1-n][1]
        dayEps.append(hourEps)
        dayDelta.append(hourDelta)
#    import matplotlib.pylab as plt
#    [plt.plot(range(24), [hourEps[h][r] for h in range(24)]) for r in range(4)]
#    plt.legend(['C1', 'C2', 'C3', 'C4'], loc='best')
#    plt.grid()
#    for j in range(3):
#        plt.figure()
#        [plt.plot(range(24), [hourDelta[h][0][r,j] for h in range(24)]) for r in range(4)]
#        plt.legend(['C1_j'+str(j), 'C2_j'+str(j), 'C3_j'+str(j), 'C4_j'+str(j)], loc='best')
#        plt.grid()
#    plt.show()
    return [dayEps, dayDelta]
        
if __name__=='__main__':    
    [dayEps, dayDelta] = get_24h_params('Run36#', 2)
    mapObj = MapLMS('Map2',  Cn=4, propa_j = ['direct', 'refl', 'scat'], betaEps = 0.001, muEps = 0.15, betaDelta=0.001, muDelta=15)
    for day in [0, 1]:
        hourEps = dayEps[day]
        hourDelta = dayDelta[day]
        for timeHour in xrange(24):
            print 'dawing map of hour -> ', timeHour
            mapObj._redrawMap('receiverMapGrid5', hourEps[timeHour], hourDelta[timeHour][0], timeHour, day)






