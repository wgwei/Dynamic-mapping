# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 18:55:44 2013
Main file to prepare the necessary file to calculate the dynamic noise map

modified from prepare_dmap_main.py in version1
@author: weigang

"""

import os, logging, sys         
from SSM_utils_1_2_1 import Road_Sources, Source_Assign_From_Buffer, Categorize_Source,\
         Reciever_Grid, CalculateATTM, Calculate_ATTM_temp_srcCategory, \
         split_polyline_to_points, assign_lvl_from_polyLine_to_points
#from calc_direct_copy_LUC import runBass
from BGM2_v1 import call_Model
#from calc_direct_copy_LUC import runBass

sys.path.append(r'C:\Python27\Lib\site-packages')
class DMap_Main():
    def __init__(self):
        print 'initialize basic input info\n'         
        self.lineShapeFileName = os.path.join(r'calculate', 'tram02_copy')
        self.distance = 10.0
        
        self.trafficLinkFile = os.path.join(r'calculate', 'link_test')        
        self.emissionModel = 'CNOSSOS'
        self.linkPowerOutput = os.path.join(r'calculate', 'lvl_octave_link_CNOSSOS')        
        
        self.linkBufferFile = os.path.join(r'calculate', 'buffer_test2')
        self.prefixEmissionPointFile = os.path.join(r'calculate',  'SSM_sources')
        self.linkBufferID = 'LINKNR'        
    
        self.cagegorizeIntervals = [[55.8, 75.7], [75.8, 88.4], [88.5, 95.4], [95.5, 114.6]]
        self.categorizeByHour = 8
        self.categorizeByHz = 1000.        
        
        self.upLeftPt =[91945.,  435514.] # the receiver zone
        self.lowRightPt = [93531.,  434461.]
        self.refineZones = [[(92223., 435082.), (92466., 435082.), (92466., 434903.), (92223., 434903.), (92223.,435082)], \
                            [(92534., 434898.), (92755., 434898.), (92755., 434788.), (92534., 434788.), (92534., 434898.)], \
                            [(92724., 435205.), (92878., 435205.), (92878., 435095.), (92724., 435095.), (92724., 435205.)], \
                            [(93147., 435091.), (93306., 435091.), (93306., 435020.), (93147., 435020.), (93147., 435091.)]]
        self.circularCenter = [92794., 435073.]
        self.refineCenterList = [[92294., 435009.], [92587., 434836.], [92730., 435194.], [93231., 435056.]]
        self.refineGridSize = 5.
        self.avoidInsideOf = os.path.join('calculate','SSM_buildings')
        self.receiverFileName = 'SSM_receivers'
        self.gridSize = 20.
        
        self.box_min = (103250,  192500)   #lower left corner of the calculation region of direct sound
        self.box_max = (105650, 195200)   # upper right corner of the calculation region
        self.pathServer = 'U:\PhD\projects\SSM\src\Rotterdam'
        self.buildingName4direct = 'SSM_buildings'
        self.buildingName4refl = 'SSM_buildings'        
        
        self.pMeasFoiIDs = [202,203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 220]
        self.Cn = 4 # numer of categorized sources
        self.propa_j = ['direct', 'refl', 'scat']
        self.pMeasFile = 'calculate\\Pmeas' # which include the fois for every measurement nodes        
        
        self.saveTag = 'industry'
        self.tempPMeasFile = 'calculate\\Pmeas'
        self.tmepPMeasFoiIDs = [202,203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 220]
        
    
    def call(self):         
        print ('\t\t1 Split polyline to point sources')
        print ('\t\t2 Categorize sources')
        print ('\t\t3 Generate receivers')
        print ('\t\t4 call background noise model -- reflections')
        print ('\t\t5 call background noise model -- scattering')
        print ('\t\t100 call direct noise model')
        print ('\t\t300 calculate the temporary sources ATTM')
        print '\n\nCall processes...\n\n'                  
        print 'Step 1: prepare sources'        
        print '    create source emission points from polyline'
        rk = input('input your option: ')
        
        if rk==1:
            assign_lvl_from_polyLine_to_points( os.path.join(r'calculate', 'tram02_copy2'), self.distance, os.path.join(r'calculate', 'tram02_copy2_pt'))   
            assign_lvl_from_polyLine_to_points( os.path.join(r'calculate', 'lokalewegen02_copy2'), self.distance, os.path.join(r'calculate', 'lokalewegen02_copy2_pt'))    
            # merge the sources to one shape file in DAY, EVE and NGT separately. The categorizing code will deal with it.
        if rk==2:
            print 'step 2: categorize the sources'
            cs = Categorize_Source(self.prefixEmissionPointFile, self.cagegorizeIntervals,  'LDAY_1K')
            cs.write_medium_for_category()
        
        if rk==3:
            print 'step 3: generate the receiver grid'
            rvr = Reciever_Grid(self.upLeftPt, self.lowRightPt, self.avoidInsideOf, self.gridSize, self.refineGridSize, self.receiverFileName)
            rvr.generate_circular_grid([92775, 435068], 600., self.refineCenterList, [50, 50, 50, 50])
        if rk==100:
            print 'step 8: call direct, reflect and scat model to calculat the sound and copy results to calculate folder'
            # call direct sound calculation. This is only available on THEREON
            runBass(self.box_min,self.box_max, self.pathServer, self.sourceD_shapeFile, self.sourceE_shapeFile, self.sourceN_shapeFile, \
                    self.immissionShapeFile, self.buildingShapeFile, outputFile = 'Rotterdam_direct',\
                    model='iso9613', misc=2, nReflections=1, nDiffractions=1, report_fileName='sim_report.txt', maxDistance=500., \
                    ground_impedantie = 0.5, meteoOrFlat = 'flat', \
                    meteo_distribution = True, meteo_propagation = 'mdist', writeSpectral=True, begin=0, end=-1)
        if rk==4:
            rcvFile = input('\t receiver file name: ')
            srcFile = input('\t source file name: ')
            outFile = input('\t output file name: ')
            print ('\t\t\t loading ', srcFile)
            print ('\t\t\t saving to ', 'calculate'+'\\'+outFile)
            call_Model(srcFile, rcvFile, self.buildingName4refl, 'calculate'+'\\'+outFile, modelType='FDTDfitting')
        if rk==5:
            rcvFile = input('\t input receiver file name: ')
            srcFile = input('\t source file name: ')
            outFile = input('\t output file name: ')
            print ('\t\t\t loading ', srcFile)
            print ('\t\t\t saving to ', 'calculate'+'\\'+outFile)
            call_Model(srcFile, rcvFile, self.buildingName4refl, 'calculate'+'\\'+outFile, modelType='scattering')
        if rk==200:
            print 'step 9: calcluate the attenuation with munus sign -> ATTM'
            attm = CalculateATTM(self.pMeasFoiIDs, self.Cn, self.propa_j, self.pMeasFile)
            attm.callDBFandSaveAsPKL()
            attm.callMapDBFandSaveAsPKLLeq()
        if rk==300:
            print 'step 10: calculate the temporary sources ATTM'
            attmT = Calculate_ATTM_temp_srcCategory(self.saveTag, self.tempPMeasFile, self.tmepPMeasFoiIDs)
#            attmT.callDBFandSaveAsPKL(['calculate\\industry_direct_Pmeas.txt', 
#                                    'calculate\\industry_refl_Pmeas.txt', 
#                                    'calculate\\industry_scat_Pmeas.txt'])
            attmT.calc_affect_receiv_save_pkl(['calculate\\industry_direct.txt', 
                                            'calculate\\industry_refl.txt', 
                                            'calculate\\industry_scat.txt'])

def main(logName):
    logging.basicConfig(filename=logName, level=logging.DEBUG)    
    if not os.path.exists('for_map'):
        os.mkdir('for_map')
        logging.info('Create the for_map folder')
    if not os.path.exists('calculate'):
        os.mkdir('calculate')  
        logging.info('Create the calculate folder')
        logging.warning('Necessary files are missing! Program fails!\n\n')
        print 'prepare the link file, the emission points and the link buffer first'
        return 0
    else:        
        dmp = DMap_Main()
        args = vars(dmp)
        for arg in args:
            logging.info(arg+' : '+str(args[arg]))
        dmp.call()

if __name__=='__main__':
    main('run_log.log')
