# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 09:49:34 2013
modified from proj_util2.py in version1
@author: weigang
"""
import os
import dbflib
import shapelib
import lass
import numpy as np
import numpy
import ast
import pickle
import shutil
import math
#from harmonoiseRoadEmission import HREbaseTotal, from13ToISOOctave
#from cnossosRoadEmission import CREbaseTotal, aWeightISOOctave

def dist2D(pti, pt2):
    return np.sqrt((pti[0]-pt2[0])**2.+(pti[1]-pt2[1])**2.)
    
class CopyDBF():    
    def __init__(self, srcFielname, aimFilename):       
        """ copy DBF file and generate a new file
            1, generate the field by self._copyDBFfield()
            2, copy the selected record by recID by 
            self._copyRecord(recID, dbfNum)
        """
        self.nDBF = dbflib.create(aimFilename)
        self.DBF = dbflib.open(srcFielname)        
    def _copyDBFfield(self):        
        # creat field
        for f in xrange(self.DBF.field_count()):
            fi = self.DBF.field_info(f)
            self.nDBF.add_field(fi[1], fi[0], fi[2], fi[3])
    def _copyRecord(self,recID, dbfNum):        
        dbfRec = self.DBF.read_record(recID)
        self.nDBF.write_record(dbfNum, dbfRec)
    def _writeRecord(self, dbfNum, dbfRec):
        self.nDBF.write_record(dbfNum, dbfRec)
    def _close(self):
        self.nDBF.close()
        self.DBF.close()

class CreateDBF():
    def __init__(self, filename):
        self.ndbf = dbflib.create(filename)
        # add a default field name as ID        

class CreateNewMap():
    def __init__(self, fileName, tp, defaultDBFcreate=1, basedOnSHP=None):
        ''' tp = shpTypes
            shpTypes = {'1--Point':       shapelib.SHPT_POINT,
                        '8--PointList':   shapelib.SHPT_MULTIPOINT,
                        '3--LineSegment': shapelib.SHPT_ARC,
                        '3--Line':        shapelib.SHPT_ARC,
                        '3--LineList':    shapelib.SHPT_ARC,
                        '5--Polygon':     shapelib.SHPT_POLYGON}
            defaultDBFcreate == 1->  it will create a default DBF file with one field named 'ID'
            basedOnShp -> if the new shape is created based on a old shape file, 
            the name should be mentioned here.
        '''
        self.shpType = {1:shapelib.SHPT_POINT,3:shapelib.SHPT_ARC,5:shapelib.SHPT_POLYGON}
        self.basedOnSHP = basedOnSHP
        self.w2shp = shapelib.create(fileName, tp)        
        self.tp = tp 
        self.defaultDBFcreate = defaultDBFcreate        
        if defaultDBFcreate==1 and basedOnSHP==None:            
            self.dbfobj = dbflib.create(fileName)
            self.dbfobj.add_field("ID", dbflib.FTInteger, 9, 0)       
        else:
            if basedOnSHP!=None:
                self.dbfobj = dbflib.create(fileName)
                self.DBF0 = dbflib.open(basedOnSHP)
                for n in xrange(self.DBF0.field_count()):
                    fi = self.DBF0.field_info(n)
                    self.dbfobj.add_field(fi[1], fi[0], fi[2], fi[3])
            else:
                print ' .dbf will not be generated automatically!'
                    
                
    def write_shp(self, idn, vertix):
        if self.tp==1: 
            vertix = [vertix]
        shpObj = shapelib.SHPObject(self.shpType[self.tp], self.tp, vertix)
        self.w2shp.write_object(idn, shpObj) 
        if self.defaultDBFcreate==1:
            self.dbfobj.write_record(idn, {"ID":idn})
        else:
            if self.basedOnSHP!=None:
                record = self.DBF0.read_record(idn)
                self.dbfobj.write_record(idn, record)                        
        
    def _close(self):
        self.w2shp.close() 
        if self.defaultDBFcreate==1 or self.basedOnSHP!=None:
            self.dbfobj.close()



def readBuffer(iFilename, ID_field):
        print 'reading simple polygon file', iFilename
        shp = shapelib.ShapeFile(iFilename)
        dbf = dbflib.open(iFilename)
        polygons = []
        for i in xrange(shp.info()[0]):
            if len(shp.read_object(i).vertices())>0:
                vs = shp.read_object(i).vertices()[0]
                vs = [(v[0] , v[1]) for v in vs[1:]]
                vs.reverse()
                polygon = SimplePolygon2D(vs) # polygon.vertices return the vertices of the  
                polygon.pid = dbf.read_record(i)[ID_field]
                polygons.append(polygon)
        shp.close()
        dbf.close()
        return polygons
        
def readPolygonFile(iFilename, ID_field):
    print 'reading simple polygon file', iFilename
    relativeHeihgtField = 'MEANREL_H'
    shp = shapelib.ShapeFile(iFilename)
    dbf = dbflib.open(iFilename)
    polygons = []
    for i in xrange(shp.info()[0]):
        if len(shp.read_object(i).vertices())>0:
            vs = shp.read_object(i).vertices()[0]
            vs = [(v[0] , v[1]) for v in vs[1:]] # [(x,y), (x,y), ...]
            vs.reverse()
            polygon = SimplePolygon2D(vs) # polygon.vertices return the vertices of the  
            polygon.relativeHeight = dbf.read_record(i)[relativeHeihgtField]
            polygon.pid = dbf.read_record(i)[ID_field]
            polygons.append(polygon)
    shp.close()
    dbf.close()
    return polygons 
        
def readPolygonByVertices(verticesList):
    ''' an rectangle example for verticesList is:
        = [[(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)] ,
        [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]]
    '''
    polygons = []
    for vs in verticesList:
        polygons.append(SimplePolygon2D(vs))
    return polygons    

class SimplePolygon2D:
    def __init__(self, iVertices):
        self.vertices = iVertices
        if not self._isCounterClockwise():
            self.vertices.reverse()
        self.dim2 = range(2)
    def aabb(self):
        aabbMin = tuple([min([p[k] for p in self.vertices]) for k in self.dim2])
        aabbMax = tuple([max([p[k] for p in self.vertices]) for k in self.dim2])
        return (aabbMin, aabbMax)   
    def contains(self, point):
        c = False
        x = point[0]
        y = point[1]
        for i in range(len(self.vertices)):
            a = self.vertices[i]
            b = self.vertices[i - 1]
            if ((a[1] <= y and y < b[1]) or (b[1] <= y and y < a[1])) and x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0]:
                c = not c
        return c
    def cross2D(self, a, b):
        return a[0] * b[1] - a[1] * b[0]  
    def _isCounterClockwise(self):
        return sum([self.cross2D(self.vertices[k - 1], self.vertices[k]) for k in xrange(len(self.vertices))]) >= 0


class KDTreeManualPolygon:
    def __init__(self, listPolygonObject):
        '''listPolygonObject: list of polygon objects generated by SimplePolygon2D'''
        self.polygonTree= lass.AabbTree2D(listPolygonObject)
#        print len(listPolygonObject), 'added to KDTree'

    def findPolygonsInside(self, point):
        ''' return polygonObj
            polygonObj can be used to query info defined in 'readPolygonFile'
            if len(polyonObj)==0, indicates the point is not inside the polygon
        '''
        objectSequence = self.polygonTree.find(point)
        return objectSequence   # objectSequence = [polygonObject1, polygonObject2, ...] polygonObj can be used to query info defined in 'readPolygonFile'

    def pointInside(self, point):
        polygonObj = self.polygonTree.find(point)
        if len(polygonObj)>0:
            return True
        else:
            return False        

class KDTreePoints:
    def __init__(self, iFileName, ID_field= 'ID'):        
        verts, pointIDs = self.shapeLoadpoints(iFileName, ID_field)
        KDTReePrep = []
        for  i, v in  zip(pointIDs, verts):
            KDTReePrep.append(lass.KdTree2D.Vertex2D(v, i))
        self.KDTRee = lass.KdTree2D(KDTReePrep)
        print len(verts), 'points added to KDTree'

    def get_closest(self, p):
        psFound = self.KDTRee.nearestNeighbour(p)
        return psFound #psFound.value is the ID value of the nearest point, psFound.position  is the vertice of the nearest point

    def shapeLoadpoints(self, iFilename, fieldID='ID', subset_fields = None):
        print 'reading file', iFilename
        shp = shapelib.ShapeFile(iFilename)
        dbf = dbflib.open(iFilename,'r')
        verts = []
        pointIDs =[]
        for i in range(shp.info()[0]): # loop all shape objects
            pointIDs.append(dbf.read_record(i)[fieldID])
            verts.append(shp.read_object(i).vertices()[0]) # the append format should be (x, y) NOT [(x1, y1), (x2, y2), ...]
        shp.close()
        dbf.close()
        return verts, pointIDs


class KDTreePoints2():
    def __init__(self, pointVerticeList, corrspIDList):
        ''' find the closest point in the input pointVerticeList and return its 
            vertice by .position and ID by .value
            pointVerticeList is [(x1, y1), (x2, y2), ...]
            corrspIDList is [ID1, ID2, ...]            
        '''
        KDTReePrep = []
        for pt, ID in zip(pointVerticeList, corrspIDList):
            KDTReePrep.append(lass.KdTree2D.Vertex2D(pt, ID))
        self.KDTRee = lass.KdTree2D(KDTReePrep)
        print len(pointVerticeList), ' points added to KDTree'
    def get_closestNear(self, p):
        ''' find the closest point to p in piontVerticeList 
            where p is (x, y)
        '''
        psFound = self.KDTRee.nearestNeighbour(p)
        return psFound #psFound.value is the ID value of the nearest point, psFound.position  is the vertice of the nearest point
         
def get_centriod(polygon):
    """ polygon = [[x1, y1], [x2, y2], [x3, y3],...]
    """
    xsum = 0.0
    ysum = 0.0
    c = 0.0
    for p in polygon:
        xsum += p[0]
        ysum += p[1]
        c += 1
    return [xsum/c, ysum/c]           
    

class Sources_Clean():
    def __init__(self, originSourceShape, avoidInsideOfBuilding, cleanedSources):
        ''' to kick out the emissions which are inside buildings 
            the originSourceShape may include some emission points which are inside
            the building and are needed to remove from the emission positions
        '''
        self.SHP0 = shapelib.ShapeFile(originSourceShape)
        self.DBF0 = dbflib.open(originSourceShape)
        self.SHPcleaned = shapelib.create(cleanedSources, 1) # 1 is shapelib.SHPT_POINT
        self.DBFcleaned = dbflib.create(cleanedSources)
        listPolygonObject = readPolygonFile(avoidInsideOfBuilding, "ID")
        self.inBuilding = KDTreeManualPolygon(listPolygonObject)   
        
        self.add_field_to_newSHP()
        
    def add_field_to_newSHP(self):
        for n in range(self.DBF0.field_count()):
            fi = self.DBF0.field_info(n)
            self.DBFcleaned.add_field(fi[1], fi[0], fi[2], fi[3])
                
    def clean_sources_in_building(self):        
        cnt = 0
        for b in xrange(self.SHP0.info()[0]):
            shpobj = self.SHP0.read_object(b)
            dbfrec = self.DBF0.read_record(b)
            vts = shpobj.vertices()
            if not self.inBuilding.pointInside(vts[0]):
                self.SHPcleaned.write_object(cnt, shpobj)
                self.DBFcleaned.write_record(cnt, dbfrec)
                cnt += 1
        print '%d sources have been removed. ' % (self.SHP0.info()[0]-cnt)
        self.close_file()
                
    def close_file(self):
        self.SHP0.close()
        self.SHPcleaned.close()
        self.DBF0.close()
        self.DBFcleaned.close()        


#==============================================================================
#============= parent class used to initialized all the classes
#================================================================================  
class Sources():
    def __init__(self):
        self.trafficLinkFile = os.path.join(r'S:\projecten\2013\dynkmap\src\loaded\traffic_link', 
                                    'links_vkc_antw_twin_link_traffic_ssgent')        
        self.bufferFile = os.path.join(r'S:\projecten\2013\dynkmap\src\loaded\traffic_link', 'link_buffer_20m_copy')        
        self.linkEmissionOutput = os.path.join(r'S:\projecten\2013\dynkmap\src\loaded\traffic_link',
                                        'lvl_octave_link_CNOSSOS')
        self.emissionFile = os.path.join(r'loaded\emission_road', 'navstr_gent_eqsources_25m_cleaned')
        
        self.octaveFreq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
        
#================================================================================
#============  calculate the total sound power of the links======================
#================================================================================  
class Road_Sources():
    def __init__(self, trafficLinkFile='links_vkc_antw_twin_link_traffic_ssgent', model='CNOSSOS'):
        ''' Created on Tue Oct 15 18:35:42 2013
            to calculate the emision of every link of road(including double directions)
            the Harmnoise module is found in the server: 
                \\acoustserv\acoust\literatuur\calculations\harmonoise\bronmodel
            According to the comparison of different harmonoise versioins, and CNOSSOS, 
            the CNOSSOS is the latest emmision model which should be recommended.
            @author: Wei
            calculate the total power level of link(or road) 
        '''
        self.model = model
        self.vehicleTypeCar = 0  # correspoding to the category 1 in harmonise
        self.vehicleTypeHeaveyTruck = 2 # correspoding to the category 3 in harmonise
        self.prefix = {1:'au', 3:'vr'} # the prefix of the field name to querry the info
        self.prefixSpeed = 'sn'        
        self.polylineVerticeFieldName = 'SE_ORIG' 
        self.DBF0 = dbflib.open(trafficLinkFile)
        self.octaveFreq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
        
    def length_of_polyLine(self, recordID):
        # vertices is stored as string. Convert them to a value
        strVts = self.get_polyLine_vertice(recordID) 
        vertices = ast.literal_eval(strVts)
        pt1 = vertices[0]
        pt2 = vertices[1]
        return numpy.sqrt(sum([(v1-v2)**2. for v1, v2 in zip(pt1, pt2)]))    
        
    def get_polyLine_vertice(self, recordID):
        return self.DBF0.read_record(recordID)[self.polylineVerticeFieldName]    
    
    def get_trafficFlow_of_type(self, vehicleType, recordID, hour):
        ''' vehicleType: 1-> lieght vehicle
                         2 -> medium heavy vehicle
                         3-> heavy vehicle
        '''
        prefix = self.prefix[vehicleType]
        return self.DBF0.read_record(recordID)[prefix+str(hour)]   
        
    def get_speed_byTimeStamp(self, recordID, hour):
        return self.DBF0.read_record(recordID)[self.prefixSpeed+str(hour)]
        
    def source_emission_single_vehicle_harmonoise(self, vehicleType, speed):
        ''' to calculate the emission of a single vehicle by vehicleType and speed
            vehicle type : 1 -> for car
                           2 -> light truck
                           3 -> for heavy truck
            speed -> km/h
            return power levels in Octave bands from 63Hz to 8000Hz
        '''
        levelTerts = HREbaseTotal(vehicleType, speed)
        levelOctave = from13ToISOOctave(levelTerts)
        return levelOctave     
    
    def source_emission_single_vehicle_cnossos(self, vehicleType, speed):
        ''' to calculate the emission of a single vehicle by vehicleType and speed
            vehicle type : 1 -> for car
                           2-> for light truck
                           3 -> for heavy truck
            speed -> km/h
            return power levels in Octave bands from 63Hz to 8000Hz
        '''
        levelOctave = CREbaseTotal(vehicleType, speed)
        return levelOctave 
        
    def intergrate_sources_of_type(self, vehicleType, recordID, hour):
        Qm = self.get_trafficFlow_of_type(vehicleType, recordID, hour)
        if Qm==0:
            Qm = 0.01
        veq = self.get_speed_byTimeStamp(recordID, hour)
        if veq>130:
            veq = 130.0
        if self.model=='CNOSSOS':
            splevel = self.source_emission_single_vehicle_cnossos(vehicleType-1, veq)  # vehicleType-1 is used to match the list format which starting from 0
        if self.model=='HARMONOISE':
            splevel = self.source_emission_single_vehicle_harmonoise(vehicleType-1, veq)  # vehicleType-1 is used to match the list format which starting from 0
        return splevel + 10.*numpy.log10(Qm/1000./veq)
        
    def intergrate_sources_allTypes(self, recordID, hour):        
        lvlcar = self.intergrate_sources_of_type(1, recordID, hour) # cars
        lvlvr = self.intergrate_sources_of_type(3, recordID, hour)  # heavy trucks
        return 10.*numpy.log10(10.**(0.1*lvlcar)+10.**(0.1*lvlvr))     
    
    def calculate_total_source_power_level(self, recordID, hour): 
        ''' return source power level of a road   '''
        lvlperm = self.intergrate_sources_allTypes(recordID, hour)
        lenPL = self.length_of_polyLine(recordID)
        return 10.*numpy.log10(10.**(0.1*lvlperm)*lenPL)   
            
    def define_field_info_link(self):
        ''' used to create the field for the new dbf file'''
        fieldInfo = [(1,'ID', 16, 0), # 1 -> integer
                     (1, 'LINKNR', 16, 0), 
                     (0, self.polylineVerticeFieldName, 40, 0)] # 0 -> string
        fieldInfo2 = [(2, 'L_'+str(fr), 16, 1) for fr in self.octaveFreq]  # 2 -> double
        return fieldInfo+fieldInfo2
        
    def sumspec_to_Leq(self, levelSpectrum):
        ''' sum up all the central frequencies to the total Leq '''
        return 10.*numpy.log10(sum([10.**(0.1*v) for v in levelSpectrum]))
        
    def sumspec_to_Laeq(self, levelSpectrum):
        ''' sum up all the central frequencies to the total LAeq '''
        lvlspecAweighted = aWeightISOOctave(levelSpectrum)
        return self.sumspec_to_Leq(lvlspecAweighted)
        
    def form_record(self, recordID, hour):
        ''' 
        '''
        LINKNR = self.DBF0.read_record(recordID)['LINKNR']
        VTSstr = self.DBF0.read_record(recordID)[self.polylineVerticeFieldName]
        lvlspec = self.calculate_total_source_power_level(recordID, hour)
        fieldValue = [recordID, LINKNR, VTSstr]+list(lvlspec)
        fieldName = ['ID', 'LINKNR', self.polylineVerticeFieldName] + ['L_'+str(fr) for fr in self.octaveFreq]  
        return {fieldName[m]:fieldValue[m] for m in range(len(fieldName))}        
        
    def close_file(self):
        self.DBF0.close()
        for n in xrange(len(self.fwobj)):
            self.fwobj[n].close()
            
    def save_emission_hour(self, outputFile):
        # add field to 24h
        self.fwobj = []
        fieldInfo = self.define_field_info_link()
        for h in xrange(24):
            fw = dbflib.create(outputFile+'_'+str(h)+'h')
            for fi in fieldInfo:
                fw.add_field(fi[1], fi[0], fi[2], fi[3])
            self.fwobj.append(fw)
        # write record out
        for h in xrange(24):
            print 'writing -> %s hour' %str(h)
            for r in xrange(self.DBF0.record_count()):            
                record = self.form_record(r, h)
                self.fwobj[h].write_record(r, record)
            self.fwobj[h].close()

#================================================================================
#============  assign the source power to emission points which are in the buffer
#================================================================================        
class Source_Assign_From_Buffer():
    def __init__(self, bufferFile, emissionFile, linkBufferID='LINKNR'):
        ''' link the emission position with the LINKNR to assign source power to the emission 
            positoin. 

            if a position is linked to more than one LINKNR(more than one links), the total
            source power should be added together by energy. 
        '''
        self.bufferFile = bufferFile
        self.emissionFile = emissionFile
        self.DBFbf = dbflib.open(bufferFile)
        self.SHPbf = shapelib.ShapeFile(bufferFile)
        self.DBFem = dbflib.open(emissionFile)
        self.SHPem = shapelib.ShapeFile(emissionFile)
        self.bfID = linkBufferID        
        self.IDlabel = 'ID'    
        self.octaveFreq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
        self.HOUR = 8 # the timestamp used to order results
        self.FREQ = 1000 # the central frequency used for ordering
        
    def define_field_info(self):
        ''' used to create the field for the new dbf file'''
        fieldInfo = [(1,self.IDlabel, 16, 0)] # 1 -> integer
        fieldInfo2 = [(2, 'L_'+str(v), 16, 1) for v in self.octaveFreq]  # 2 -> double
        return fieldInfo+fieldInfo2
        
    def form_record(self, ID, L8octave):
        fieldnames = [self.IDlabel] + ['L_'+str(v) for v in self.octaveFreq]
        fieldvalues = [ID] + list(L8octave)
        record = {fieldnames[n]:fieldvalues[n] for n in xrange(len(fieldnames))}
        return record   
    
    def put_all_buffers_to_KDtree(self):
        listPolygonObject = readBuffer(self.bufferFile, self.bfID)
        inBuffers = KDTreeManualPolygon(listPolygonObject) 
        return inBuffers
        
    def put_buffer_to_KDTree(self, bufferShapeID):
        ''' return the KDTree which is used to judge whether a emission point 
            is in the buffer zone or not
        '''
        vs = self.SHPbf.read_object(bufferShapeID).vertices()[0]
        vs = [(v[0] , v[1] ) for v in vs[1:]]
        vs.reverse()
        polygon = SimplePolygon2D(vs) # polygon.vertices return the vertices of the  
        polygon.pid = self.DBFbf.read_record(bufferShapeID)[self.bfID]
        inBuffer = KDTreeManualPolygon([polygon] )  
        return inBuffer
    
    def calc_linknr_division_portion(self):
        ''' loop all buffers. For every buffer, loop all the emission points 
            count the emissions which are in the buffer and the counted number 
            is the divisioin portion of the link specified by the linknr
            return a dictionary
        '''
        linknrBufferDivision = {}
        for b in xrange(self.SHPbf.info()[0]):
            print '     Analyzing buffer -> ', b
            inBuffer = self.put_buffer_to_KDTree(b)
            insideEmissionCount = 0
            for bem in xrange(self.SHPem.info()[0]):
                shpObjEms = self.SHPem.read_object(bem) 
                if inBuffer.pointInside(shpObjEms.vertices()[0]):
                    insideEmissionCount += 1
            LINKNR = self.DBFbf.read_record(b)[self.bfID]
            linknrBufferDivision.update({LINKNR:insideEmissionCount})
        fp = open('linknrBufferDivision.pkl', 'wb')
        pickle.dump(linknrBufferDivision, fp)
        fp.close()
        return linknrBufferDivision       
        
    def dict_emission_linknr_list(self):
        ''' loop the emission points and find the points inside
            the specifided buffer and put hang them to the emission_linknr dictionary
            return the dictionary for every 
        '''
        emissionLinknr = {r:[] for r in range(self.DBFem.record_count())}  
        inBuffers = self.put_all_buffers_to_KDtree()
        for bem in xrange(self.SHPem.info()[0]):
            shpObjEms = self.SHPem.read_object(bem) 
            linknrsSequ = inBuffers.findPolygonsInside(shpObjEms.vertices()[0]) # return the sequence of polygons if the point is inside it
            if linknrsSequ: # if linknrSequ is not empty
                for obj in linknrsSequ:
                    emissionLinknr[bem].append(obj.pid)   
        fp = open('emissionLinknr.pkl', 'wb')
        pickle.dump(emissionLinknr, fp)
        fp.close()
        return emissionLinknr
    
    def load_link_power_as_dict(self, linkWithPowerFile):
        ''' load the link dbf file and dict them as LINKNR:numpy.asarray([L_63, L_125, ..., L_8000])
        '''        
        DBFlinkpower = dbflib.open(linkWithPowerFile)
        linknrPower8octave = {}
        fieldNames = ['L_'+str(v) for v in self.octaveFreq]
        for r in xrange(DBFlinkpower.record_count()): 
            linknrPower8octave.update({DBFlinkpower.read_record(r)[self.bfID] : np.asarray([DBFlinkpower.read_record(r)[v] for v in fieldNames])})
        DBFlinkpower.close()
        return linknrPower8octave
        
    def close_file(self):
        self.DBFbf.close()
        self.SHPbf.close()
        self.DBFem.close()
        self.SHPem.close()
        
    def assign_sourcePower_to_emission_of(self, linkWithPowerFile, hour):
        ''' The original link file includes only the traffic info of 24 hours
            According to the CNOSSOS or HARMONIOSE model, the sound power per link
            is calculated and the spectrum is stored for all the 24 hours in different
            file name.
            linkWithPowerFile is the link file with spectrum of sound power level of an hour
            hour is the time stamp which is used to identify the output file names
        '''
        print 'processing -> emission_recordID:[LINKNR, LINKNR2, ...]'
        if not os.path.exists('emissionLinknr.pkl'):
            print '     calculating...' 
            emissionLinknr = self.dict_emission_linknr_list() # emission_recordID:[LINKNR, LINKNR2, ...]
        else:
            print '     loading...'  
            fp = open('emissionLinknr.pkl', 'r')
            emissionLinknr = pickle.load(fp)   
            fp.close()
            
        print 'processing -> LINKNR:buffer division'
        if not os.path.exists('linknrBufferDivision.pkl'):
            print '     calculating...' 
            linknrBufferDivision = self.calc_linknr_division_portion() # LINKNR:buffer division
        else:
            print '     loading...'  
            fp2 = open('linknrBufferDivision.pkl', 'r')
            linknrBufferDivision = pickle.load(fp2)        
            fp2.close()
            
        print 'processing -> LINKNR:[L_63, L_125,..., L_8000] in the power file'
        linknrPower8octave = self.load_link_power_as_dict(linkWithPowerFile) # LINKNR:[L-63, L_125] in the power file
        
        print 'processing -> assign emissioin points...'
        newDBFem = dbflib.create(self.emissionFile+'_'+str(hour)+'h')
        for fi in self.define_field_info():         
            newDBFem.add_field(fi[1], fi[0], fi[2], fi[3])
        for r in xrange(self.DBFem.record_count()):
            lvl8oct = np.asarray([0.0]*8)
            linknrs = emissionLinknr[r]
            if linknrs: # if not empty
                for linknr in linknrs:
                    lvl8oct += 10.**(0.1*linknrPower8octave[linknr]) / linknrBufferDivision[linknr]
                record = self.form_record(r, 10.*np.log10(lvl8oct))
            else:
                record = self.form_record(r, lvl8oct)
            newDBFem.write_record(r, record)
        newDBFem.close()
        
    def assign_sourcePower_to_emission(self, prefix):
        for h in xrange(24):
            print '\nassign %s:00' %str(h)
            linkWithPowerFileName = prefix + '_' + str(h)+'h'
            self.assign_sourcePower_to_emission_of(linkWithPowerFileName, h)
    
    def select_min_power_level(self, hour):
        ''' sort the source power at emission points for every hour  and select the min record
        '''
        self.DBF8hOctave = dbflib.open(self.emissionFile+'_'+str(hour)+'h')
        lvlHour = []
        for r in xrange(self.DBF8hOctave.record_count()):
            lvl8octave = [self.DBF8hOctave.read_record(r)[v] for v in ['L_'+str(fr) for fr in self.octaveFreq]]
            if not np.asarray(lvl8octave).any()==0:
                lvlHour.append(lvl8octave)
        minLvl = sorted(lvlHour, key=lambda x: x[self.octaveFreq.index(self.FREQ)])[0] # sort the list by 1000Hz and select the min record
        return minLvl
    
    def replace_0_byMin_of(self, hour):
        ''' Most of the emission points will be assigned by the closest link power, however, 
            some of them are not inside the link buffer and they will keep zeros as the power
            level. This function replaces these zeros by the min power level of the assigned 
            emission points.
        '''
        print 'replaceing 0 of -> %s hour' %str(hour)
        minLevel = self.select_min_power_level(hour)
        self.DBFatH = dbflib.open(self.emissionFile+'_'+str(hour)+'h')
        DBFnew = dbflib.create(self.emissionFile+'_'+str(hour)+'h'+'_replace_0_byMin')
        for fi in self.define_field_info():            
            DBFnew.add_field(fi[1], fi[0], fi[2], fi[3])
        for r in xrange(self.DBFatH.record_count()):
            lvl8oct = [self.DBFatH.read_record(r)[v] for v in ['L_'+str(fr) for fr in self.octaveFreq]]
            if np.asarray(lvl8oct).any()==0:
                record = self.form_record(r, minLevel)
            else:
                record = self.form_record(r, lvl8oct) 
            DBFnew.write_record(r, record)
        DBFnew.close()        
        self.DBFatH.close()
        
    def replace_0_byMin(self):
        ''' loop 24 hours to find the min level for every hour '''
        for h in xrange(24):
            self.replace_0_byMin_of(h)
        # copy other necessary shape file for source categorization
        shutil.copyfile(self.emissionFile+'.shp', self.emissionFile+'_'+str(self.HOUR)+'h'+'_replace_0_byMin.shp' )
        shutil.copyfile(self.emissionFile+'.shx', self.emissionFile+'_'+str(self.HOUR)+'h'+'_replace_0_byMin.shx' )
        
    def put_min_level_to_unassigned_emission_pt(self, prefix):
        if os.path.exists(self.emissionFile+'_0h_replace_0_byMin.dbf'):
            print 'found ', self.emissionFile+'_xxh_replace_0_byMin.dbf', ' start to replace'
            self.replace_0_byMin()
        else:
            print 'Basice emission is unavailable, calculting it first..'
            self.assign_sourcePower_to_emission(prefix)
            self.replace_0_byMin()
        self.close_file()
        
    
#==============================================================================
#============================ categorize sources ==============================
#==============================================================================    
def load_config(fn='init.txt'):
    ''' return : values = [parentSourceFileName, outputReceiverGrid, buildingsShape, [uppper left receiver zone],
                           [lower left receiver zone], [[source category intervals], [], ...]]
    '''
    f = open(fn, 'r')
    linestrs = f.readlines()
    values = []
    for n, s in enumerate(linestrs):
        if  len(s.strip('\n'))>0: # avoid including the enter at the end
            stringValue = s.strip('\n').split(':')[1] # get the string after :
            if n>2:
                values.append(ast.literal_eval(stringValue.strip(' ')))
            else:
                values.append(stringValue.strip(' '))
    f.close()
    return values 
    
class Categorize_Source():
    def __init__(self,prefixEmissionFile, categorizeIntervals, categorizeFieldKey):
        ''' put all the initialized file into init folder
            put the loaded file in to loaded folder.
            the categorization is based on 8:00 and 1000 Hz in default. 
            The new generated category shapes are empty. The shape file is only 
            showing the emission posisions of every category. The txt file will
            store the spectrum for every hour, which is in style 24rows * 8 columns. 
            Every row is the spectrum of a specified hour of a specified category. 
        '''       
        self.sourceMap = prefixEmissionFile
        self.attribute = categorizeFieldKey
        self.intervals = categorizeIntervals
        self.octaveFreq = [63, 125, 250, 500, 1000, 2000, 4000, 8000]
        
    def inside(self, value, interval):
        if value>=interval[0] and value<interval[1]:
            return True
        else:
            return False 
        
    def copy_shp_exclude(self, sourceMap, aimPath, exceptID):
        SHP0 = shapelib.ShapeFile(sourceMap)
        SHPnew = shapelib.create(aimPath, 1) # 1 for point
        DBFnew = dbflib.create(aimPath)
        otherFields = ['ID', 'POINTID', 'POLYID', 'L_63', 'L_125', 'L_250', 'L_500', 'L_1000', 'L_2000', 'L_4000', 'L_8000']
        for n, ff in enumerate(otherFields):
            if n<3:
                DBFnew.add_field(ff, 1, 9, 0) # 1 is dbflib.FTInteger
            else:
                DBFnew.add_field(ff, 2, 16, 1) # 2 is dbflib.FTDouble
        cnt = 0
        for b in xrange(SHP0.info()[0]):
            if b not in exceptID:
                shpobj = SHP0.read_object(b)
                SHPnew.write_object(cnt, shpobj)
                attributes = [cnt, cnt, cnt]+[100.0]*8
                record = {otherFields[n]:attributes[n] for n in range(len(otherFields))}                    
                DBFnew.write_record(cnt, record)
                cnt += 1
        SHP0.close()
        SHPnew.close()
        DBFnew.close()
        
    def categrize_sources_by(self, intervals, sourceMap, attributeName):
        """ the intervals is extracted from ArcGIS and categrorized by Natrual break jenks method
            attributeName is the attribute used in ArcGIS to categorize. 
        """
        print "Categorization based on: ", sourceMap
        print 'intervals -> ', intervals
        exceptIDs = []
        categNUM = len(intervals)
        for n in xrange(categNUM):
            exceptIDs.append([]) # exceptIDs = [[],[],[]...] generate empty list. every one is corresponding to a interval
        srcDBF = dbflib.open(sourceMap) 
        print "reading source map.."
        # scan every record in the shape file, categorize them and put the except index in the list
        for n in xrange(srcDBF.record_count()):
            v = srcDBF.read_record(n)[attributeName]
            for m,intv in enumerate(intervals): # m is coresponding the intervals index. intv = intervals[m]
                if not self.inside(v, intv):
                    exceptIDs[m].append(n)
        print len(exceptIDs[0]), len(exceptIDs[1]), len(exceptIDs[2]), len(exceptIDs[3])
        print "Categarizing source map.."
        for idx, exceptList in enumerate(exceptIDs):
            print "    source categrory ", idx+1
            self.copy_shp_exclude(sourceMap, "calculate\\srcCateg"+str(idx+1), exceptList)
        return exceptIDs
    
        
    def calc_medium_for_categrory(self, exceptID):
        ''' calculate the medium value of a specified category looping all the hours.
            The methodology is to exclude the IDs in exceptID and select 
            the medium value of the remaining IDs. 
            Here supposes the source categories share the same positions in different hours
            but with different spectrum
        '''
        frStr = ['63', '125', '250', '500', '1K', '2K', '4K', '8K']
        denStr = ['DAY', 'EVE', 'NGT']
        mediumValueSpec = []
        for den in denStr:
            print '         procecing file ', self.sourceMap+'_'+den
            csrcDBF = dbflib.open(self.sourceMap+'_'+den)
            lvlspec = []
            for n in xrange(csrcDBF.record_count()):
                if n not in exceptID:
                    lvlspec.append([csrcDBF.read_record(n)['L'+den+'_'+fr] for fr in frStr])           
            lvlspecDEN = []
            for p, fr in enumerate(frStr):
                sortByfr = sorted(lvlspec, key=lambda x: x[p])
                lvlspecDEN.append(sortByfr[int((csrcDBF.record_count()-len(exceptID))/2.)][p]) # p is the index of fr. the lambda is used to sort by taking position p
            mediumValueSpec.append(lvlspecDEN)
            csrcDBF.close()
        return mediumValueSpec
            
    def write_medium_for_category(self):     
        ''' prefixOfSourceFile used to specify which file is used to calculate 
            the medium value of every category
        '''
        exceptIDs = self.categrize_sources_by(self.intervals, self.sourceMap+'_DAY', self.attribute)
        print "Writing out.."
        for c, exceptID in enumerate(exceptIDs):
            print '     prepareing category -> ', c+1
            mediumValueSpec = self.calc_medium_for_categrory(exceptID) 
            fc = open("calculate\\srcCateg"+str(c+1)+".txt", 'wb')
            for medValhour in mediumValueSpec:
                fc.write(" ".join("%0.1f" % v for v in medValhour)+'\r\n')
            fc.close()
        print "Done"        
                    
#==============================================================================
#============================ simplify polygons ==============================
#==============================================================================   
def distance(p1, p2):
    return np.sqrt(sum((p[0]-p[1])**2. for p in zip(p1, p2)))

def remove_small_segments(pointList, tolerance): # may be not used
    IDs = [0]
    i, j = 0, 0
    while j< len(pointList)-2:
        d12 = distance(pointList[j], pointList[j+1])
        d23 = distance(pointList[j+1], pointList[j+2])
        if d12>tolerance and d23>tolerance:
            IDs.append(j+1)
            j += 1     
        else:
            j += 2       
            IDs.append(j)
    IDs.append(0) # to make the polygon points style
    newPoints = [pointList[n] for n in IDs]
    return newPoints
    
def order_weighting(pointList, percentage, verticeIndex):
    IDsWeightsIndex = []
    pointList = pointList[0:-1]
    for j in xrange(len(pointList)):            
        if j<len(pointList)-1:
            s1 = distance(pointList[j-1], pointList[j])
            s2 = distance(pointList[j], pointList[j+1])
            a = np.asarray(pointList[j-1]) - np.asarray(pointList[j])
            b = np.asarray(pointList[j+1]) - np.asarray(pointList[j])                
           
        else:
            s1 = distance(pointList[j-1], pointList[j])
            s2 = distance(pointList[j], pointList[0])
            a = np.asarray(pointList[j-1]) - np.asarray(pointList[j])
            b = np.asarray(pointList[0]) - np.asarray(pointList[j])
        try:
            beta = math.acos((a[0]*b[0]+a[1]*b[1])/(np.sqrt(a[0]**2.+a[1]**2.)*np.sqrt(b[0]**2.+b[1]**2.)))/2./np.pi*360
        except ValueError:
            beta = 90.0
        IDsWeightsIndex.append([j, beta*s1*s2/(s1+s2), verticeIndex]) #
    return IDsWeightsIndex
        
def simplify_polygon_by_evaluate_weighting(verticesList, percentage):
    ''' pointList is the list of polygons which in format like [[(x0, y0), (x1, y1), ...
        (xn, yn), (x0, y0)], [], []]
        percentage is the percentage of ignorance
        return a same style but simplified polygon
    '''
    if percentage==0:
        return verticesList
    else:
        IDsWeights = []
        for idx, vertices in enumerate(verticesList):
            IDsWeightsIndex = order_weighting(vertices, percentage, idx)
            IDsWeights += IDsWeightsIndex 
        IDsWeights = sorted(IDsWeights, key=lambda x: x[1])
        # collect IDs for all sub polygons
        try:
            excludingList = [[IDsWeights[jj][0], IDsWeights[jj][2]] for jj in xrange(int(len(IDsWeights)*percentage))]
            excludingIDs = [[] for _ in xrange(len(verticesList))] # every elements in the list is corresponding to a sub-polygon
            for x in excludingList:
                excludingIDs[x[1]].append(x[0])
            newVerticesList = []
            for idx, vertices in enumerate(verticesList):
                newVertices = []
                if len(vertices)-int(len(vertices)*percentage)>7: # only simlifiy the polygon when the edge number is greater than 5 after simplification
                    for j in xrange(len(vertices)-1):
                        if j not in excludingIDs[idx]:
                            newVertices.append(vertices[j])
                    newVertices += [newVertices[0]] # to keep the start and end point are the same one
                else:
                    newVertices = vertices
                newVerticesList.append(newVertices)
        except:
            newVerticesList = verticesList                            
    return newVerticesList

    
class Polygons():
    def __init__(self, buildingFileName):
        self.polygons = buildingFileName
        self.SHP0 = shapelib.ShapeFile(buildingFileName)
        
    def run_simplify(self, verticesList, percentage, run=2):
        print ('\tbefore simplify %d nodes' %(len(verticesList[0])))
        verticesList = simplify_polygon_by_evaluate_weighting(verticesList, percentage)
       
        if run>1:
            for r in xrange(run-1):
                verticesList = simplify_polygon_by_evaluate_weighting(verticesList, percentage)
        print ('\tafter simplify %d nodes' %(len(verticesList[0])))
        return verticesList
            
    def simplify_polygons(self, percentage, run=2):
        cnt = 0
        newPolyObj = CreateNewMap(self.polygons+'_simplified', 5, defaultDBFcreate=0, basedOnSHP=self.polygons)        
        for b in xrange(self.SHP0.info()[0]):
            print 'simplify -> ', b
            shpobj = self.SHP0.read_object(b)            
            smplPolygon = self.run_simplify(shpobj.vertices(), percentage, run=run)
            newPolyObj.write_shp(cnt, smplPolygon)
            cnt += 1    
        newPolyObj._close()

#==============================================================================
#============================ receiver ==============================
#==============================================================================   

class Reciever_Grid():
    def __init__(self, upLeft, lowRight, avoidInsideOf, gridSize, refineGridSize, newFileName):        
        """ the upLeft = [x1, y1] is the vertix of the start point which should be the 
            upper left of a region; 
            lowRight = [x2, y2] is the lower right of the region
            avoidInsideOf = 'string': is the file name of the buildings to avoid puting
            receivers inside the buildings.
            gridSize = 10. length of the grid square
            refineGridSize = 5. is the length of the refining grid square
            newFileName = 'string' shape file name
            
        """
        self.startPt = upLeft
        self.endPt = lowRight
        self.avoidInsideOf = avoidInsideOf
        self.gridSize = gridSize
        self.refineGridSize = refineGridSize        
        self.newFileName = newFileName
        print '\nupLeft -> ', self.startPt
        print 'lowRight ->   ', self.endPt
        
    def generate_rectangular_grid(self, refineZones):
        ''' frefineZones: [[(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)] ,
                            [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]]
        '''
        cobj = CreateNewMap(self.newFileName, 1, 0)
        DBFOut = dbflib.create(self.newFileName)
        DBFOut.add_field("ID", dbflib.FTInteger, 16, 0)  # add the idential ID
        DBFOut.add_field("X", dbflib.FTDouble, 16, 2)   
        DBFOut.add_field("Y", dbflib.FTDouble, 16, 2)
        DBFOut.add_field("Z", dbflib.FTDouble, 16, 2)
        DBFOut.add_field('MEANREL_H', dbflib.FTDouble, 16, 2)   
        if self.avoidInsideOf != 0:
            listPolygonObject = readPolygonFile(self.avoidInsideOf, "ID")
            inBuilding = KDTreeManualPolygon(listPolygonObject)            
        listPolygonOjbect2 = readPolygonByVertices(refineZones)
        inRefine = KDTreeManualPolygon(listPolygonOjbect2)
        idn = 0
        # generate the refine zones
        for rfZone in refineZones:
            xmin = np.min(np.asarray(rfZone)[:,0])
            xmax = np.max(np.asarray(rfZone)[:,0])
            ymin = np.min(np.asarray(rfZone)[:,1])
            ymax = np.max(np.asarray(rfZone)[:,1])
            for x in np.linspace(xmin, xmax, (xmax-xmin)/self.refineGridSize):
                for y in np.linspace(ymin, ymax, (ymax-ymin)/self.refineGridSize):
                    if not inBuilding.pointInside((x,y)):
                        cobj.write_shp(idn, [(x,y)])
                        xyzh = {'ID':idn, 'X':x, 'Y':y, 'Z':4., 'MEANREL_H':4.}
                        DBFOut.write_record(idn, xyzh)                    
                        idn += 1
        
        # generate normal grid
        x = self.startPt[0]        
        while x<self.endPt[0]:
            y = self.startPt[1]
            while y>self.endPt[1]:
                vt = (x, y)
                if (not inBuilding.pointInside(vt)) and (not inRefine.pointInside(vt)):
                    cobj.write_shp(idn, [vt])
                    xyzh = {'ID':idn, 'X':vt[0], 'Y':vt[1], 'Z':4., 'MEANREL_H':4.}
                    DBFOut.write_record(idn, xyzh)                    
                    idn += 1
                y = y-self.gridSize
            x = x+self.gridSize            
        cobj._close()
        DBFOut.close()
        
    def generate_circular_grid(self,circularCenter, radius, refineCenterList, refineRadiusList):
        print ('circular center is:', circularCenter)
        print ('radius is: ', radius)
        cobj = CreateNewMap(self.newFileName, 1, 0)
        DBFOut = dbflib.create(self.newFileName)
        DBFOut.add_field("ID", dbflib.FTInteger, 16, 0)  # add the idential ID
        DBFOut.add_field("X", dbflib.FTDouble, 16, 2)   
        DBFOut.add_field("Y", dbflib.FTDouble, 16, 2)
        DBFOut.add_field("Z", dbflib.FTDouble, 16, 2)
        DBFOut.add_field('MEANREL_H', dbflib.FTDouble, 16, 2)   
        if self.avoidInsideOf != 0:
            listPolygonObject = readPolygonFile(self.avoidInsideOf, "ID")
            inBuilding = KDTreeManualPolygon(listPolygonObject)            
        idn = 0
        # generate the refine zones
        for cen, rad in zip(refineCenterList, refineRadiusList):
            xmin = cen[0]-rad
            xmax = cen[0]+rad
            ymin = cen[1]-rad
            ymax = cen[1]+rad
            for x in np.linspace(xmin, xmax, (xmax-xmin)/self.refineGridSize):
                for y in np.linspace(ymin, ymax, (ymax-ymin)/self.refineGridSize):
                    r = dist2D(cen, (x, y))
                    if (not inBuilding.pointInside((x,y))) and (r<rad):
                        cobj.write_shp(idn, [(x,y)])
                        xyzh = {'ID':idn, 'X':x, 'Y':y, 'Z':4., 'MEANREL_H':4.}
                        DBFOut.write_record(idn, xyzh)                    
                        idn += 1
        
        # generate normal grid
        x = circularCenter[0]-radius       
        while x<circularCenter[0]+radius:
            y = circularCenter[1]+radius
            while y>circularCenter[1]-radius:
                vt = (x, y)
                r = dist2D(circularCenter, vt)
                rfr = np.asarray([dist2D(rfcen, vt) for rfcen in refineCenterList])
                if (not inBuilding.pointInside(vt)) and r<radius :
                    cobj.write_shp(idn, [vt])
                    xyzh = {'ID':idn, 'X':vt[0], 'Y':vt[1], 'Z':4., 'MEANREL_H':4.}
                    DBFOut.write_record(idn, xyzh)                    
                    idn += 1
                y = y-self.gridSize
            x = x+self.gridSize            
        cobj._close()
        DBFOut.close()
        
    
#==============================================================================
#============================ calculate ATTM = sum_n^{C_c} 10.^(-0.1A_{f,n,j}(p))
#==============================================================================                   
class CalculateATTM():
    def __init__(self, PmeasFOIs=[202,203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 220], Cn=4, propa_j = ['direct', 'refl', 'scat'], PmeasFile = 'Pmeas'):
        """ This class is based on BGL_16 and direct sounf field calculation
        
            calculate the sum of the attenuation. In latex format, it is:
            sum_{i}^{N_i} 10^{-0.1A_{f,i,j}(p)}
            ATT is attenuation; M is to specify the minus symbol in the power
        """        
        self.fieldSpec = ['L_63', 'L_125', 'L_250', 'L_500', 'L_1000', 'L_2000', 'L_4000', 'L_8000', 'LAD']
        self.fieldSpec2 = ['L_63', 'L_125', 'L_250', 'L_500', 'L_1000', 'L_2000', 'L_4000', 'L_8000', 'LADAY']
        self.AWeighting = np.array([-26.2, -16.1, -8.6, -3.2, 0.0, 1.2, 1.0, 1.1]) # 63 125 250 500 1000 2000 4000 8000Hz  
        self.PmeasFOIs = PmeasFOIs
        self.Cn = Cn
        self.propa_j = propa_j
        self.crspfoi_simID = self.link_foi_simID(PmeasFile)
        self.sourceNumInCategory = self.get_source_num_by_category() # source category1, 2, 3, and 4    
        
    def get_source_num_by_category(self):
        ''' calclualte how many sources in one category '''
        snum = []
        for n in xrange(self.Cn):
            DBF = dbflib.open('srcCateg'+str(n+1))
            snum.append(DBF.record_count())
        return snum
        
    def link_foi_simID(self, DBFpmeas):
        print DBFpmeas
        DBF = dbflib.open(DBFpmeas)
        crspfoi_simID = {v:[] for v in self.PmeasFOIs}
        for foi in self.PmeasFOIs:
            [crspfoi_simID[foi].append(r) for r in range(DBF.record_count()) if DBF.read_record(r)['foi']==foi]
        DBF.close()
        return crspfoi_simID    
        
    def _loadsrc(self, srcCategory=4):
        """ calculate the source from sepctrum to total LA
        """
        Cc = [[] for _ in xrange(srcCategory)]
        Cc_f = [[] for _ in xrange(srcCategory)]
        CfcNumList = []
        for n in xrange(srcCategory):
            for h in xrange(3):
                LCn = np.loadtxt('srcCateg'+str(n+1)+'.txt')[h, :]# Fro rotterdam they are probably weighed +self.AWeighting # h row is for day, evening and night
                v = 10.*np.log10(sum(10**(0.1*LCn)))
                Cc[n].append(v)
                Cc_f[n].append(LCn)
            CfcNumList.append([Cc_f[n], Cc[n], self.sourceNumInCategory[n]])
        if not os.path.exists('for_map'):
            os.mkdir('for_map')
        f = open('for_map\\C_f_c_Num.pkl', 'wb')
        pickle.dump(CfcNumList, f)
        f.close()
        return [Cc, Cc_f]
        
    def callDBFandSaveAsPKL(self):
        ''' status: approved
            put the calcualted results in the calculate folder!!
        '''
        CCF = numpy.asarray([100.0]*8)  # when calculating the direct field, the spectrum is supposed to be 100.0dB
        CCT = 10.0*numpy.log10(sum(10.**(0.1*CCF)))
        [Cc, Cc_f] = self._loadsrc()
        foiCategATTM = []
        for foi in self.PmeasFOIs: # loop foi positions     
            print 'processing foi -> ', foi
            categATTM = [[] for _ in xrange(self.Cn)]
            recIDs = self.crspfoi_simID[foi]
            for c in xrange(self.Cn): # loop source categories
                for j in xrange(len(self.propa_j)): # loop propagation manner
                    fn  = 'calculate'+'\\'+'C'+str(c+1)+'_'+self.propa_j[j]+'_Pmeas' # 'C1_direct_Pmeas'
                    dbf = dbflib.open(fn)
                    if 'direct' in fn:
                        toatalL, sepctrumL = [], []
                        LcalDirectSpectrum = np.loadtxt(fn+'.txt')                        
                        for r in recIDs:
                            rec = dbf.read_record(r)
                            toatalL.append(rec['DBA_D'])
                            sepctrumL.append(LcalDirectSpectrum[r, 2:10])  # 1st cloumn is the record ID, 2nd to 8th colum are sepctrum
                        ATTMspectrum = 10.**(0.1*(np.mean(sepctrumL, 0)-CCF))    # category c, 0 -> spectrum, 0 -> hour      
                        ATTMtotal = 10.**(0.1*(np.mean(toatalL)-CCT))           
                    else:
                        attributes = [[] for _ in xrange(len(recIDs))]
                        for numb, r in enumerate(recIDs):
                            rec = dbf.read_record(r)                            
                            for fd in self.fieldSpec:
                                attributes[numb].append(rec[fd])       
                        meanAttributes = np.mean(attributes, 0)                        
                        ATTMspectrum = 10**(0.1*(meanAttributes[0:8] - Cc_f[c][0]))
                        ATTMtotal = 10.**(0.1*(meanAttributes[8] - Cc[c][0]))
                    categATTM[c].append([ATTMspectrum, ATTMtotal])
            foiCategATTM.append(categATTM)
        foiCategATTMdict = {self.PmeasFOIs[n]:foiCategATTM[n] for n in xrange(len(self.PmeasFOIs))}
        f = open('for_map\\foi_categATTMdict.pkl', 'wb')
        pickle.dump([foiCategATTMdict, self.PmeasFOIs], f)
        f.close()
        return foiCategATTMdict 
        
    def callMapDBFandSaveAsPKLLeq(self):
        ''' the functioin is the same as _callMapDBFandSaveAsPKL excpet changing 
            the spectrum to [0] to decrease the size of the saving file
        '''
        CCF = numpy.asarray([100.0]*8)  # when calculating the direct field, the spectrum is supposed to be 100.0dB
        CCT = 10.0*numpy.log10(sum(10.**(0.1*CCF)))
        [Cc, Cc_f] = self._loadsrc()
        # load direct field. This is big file so load it first and use it later
        LcalDirectSpectrum =[np.loadtxt('calculate\\C1_direct.txt'), np.loadtxt('calculate\\C2_direct.txt'), np.loadtxt('calculate\\C3_direct.txt'), np.loadtxt('calculate\\C4_direct.txt')]
        receiverCategATTM = []
        dbf0 = dbflib.open('calculate\\C1_direct')
        RECNUM = dbf0.record_count()
        IDs = []
        for r in xrange(RECNUM):
            IDs.append(dbf0.read_record(r)['ID'])
            print 'processing receiv ', r
            categATTM = [[] for _ in xrange(self.Cn)]
            for c in xrange(self.Cn): # loop source categories
                for j in xrange(len(self.propa_j)): # loop propagation manner
                    fn  = 'calculate'+'\\'+'C'+str(c+1)+'_'+self.propa_j[j] # 'C1_direct'
                    dbf = dbflib.open(fn)
                    if 'direct' in fn:
                        rec = dbf.read_record(r)
                        toatalL = rec['DBA_D']
                        sepctrumL = LcalDirectSpectrum[c][r, 2:10]  # 1st cloumn is the record ID, 2nd to 8th colum are sepctrum
                        ATTMspectrum = 10.**(0.1*(sepctrumL-CCF)) #the comment part is used for optimizing the spectrum
                        ATTMtotal = 10.**(0.1*(toatalL-CCT))      
                        categATTM[c].append([ATTMspectrum, ATTMtotal])
                    else:
                        attributes = []
                        rec = dbf.read_record(r)                            
                        for fd in self.fieldSpec2:
                            attributes.append(rec[fd])       
                        ATTMspectrum = 10**(0.1*(np.asarray(attributes[0:8]) - Cc_f[c][0])) # the comment part is used for optimizing the spectrum
                        ATTMtotal = 10.**(0.1*(attributes[8] - Cc[c][0]))
                        categATTM[c].append([ATTMspectrum, ATTMtotal])
            receiverCategATTM.append(categATTM)        
        print 'arrangin and writing out'
        receiverCategATTMdict = {ID:receiverCategATTM[n] for n, ID in enumerate(IDs)}
        f = open('for_map\\receiver_categATTMdict_IDs_Leq.pkl', 'wb')
        pickle.dump([receiverCategATTMdict, IDs], f)
        f.close()
        return receiverCategATTMdict
    
def calc_ATTM(calculatedSpectrum, sourceSpectrum):
    ATTMspectrum = 10**(0.1*(calculatedSpectrum - sourceSpectrum))
    ATTMtotal = 10.**(0.1*(10.*numpy.log10(numpy.sum(10.**(0.1*calculatedSpectrum))) - \
                           10.*numpy.log10(numpy.sum(10.**(0.1*sourceSpectrum)))))
    return [ATTMspectrum, ATTMtotal]
    
def get_mean(total, rowsOfMean):
    gm = []
    for r in rowsOfMean:
        gm.append(total[r])
    return np.mean(gm, 0)
    

class Extra_sources():
    def __init__(self, sourceDBFFile):
        dbf = dbflib.open(sourceDBFFile)
        self.recNum = dbf.record_count()

    
class Calculate_ATTM_temp_srcCategory():
    def __init__(self, saveTag, receiverDBFfile, foiIDs):
        ''' calculate the ATTM for the temporary souce category and for the 
            specialized receiver
            the extra sources are assigned 100dB for all frequencies at first
        '''
        self.foiIDs = foiIDs
        self.saveTag = saveTag        
        self.fieldSpec = ['L_63D', 'L_125D', 'L_250D', 'L_500D', 'L_1000D', 'L_2000D', 'L_4000D', 'L_8000D', 'LAD']
        self.crspfoi_simID = self.link_foi_simID_tempsrc(foiIDs, receiverDBFfile)
    
    def link_foi_simID_tempsrc(self, foiIDs, DBFpmeas):
        DBF = dbflib.open(DBFpmeas)
        crspfoi_simID = {v:[] for v in foiIDs}
        for foi in foiIDs:
            [crspfoi_simID[foi].append(r) for r in range(DBF.record_count()) if DBF.read_record(r)['foi']==foi]
        DBF.close()
        return crspfoi_simID             
        
    def callDBFandSaveAsPKL(self, propaFileName):
        ''' status: approved
            put the calcualted results in the calculate folder!!
            @ propaFileName is file stored the results of every propagation. 
                 it should be in a txt file format   
                 The file is in the order direct, refl, scat
            @return the dictionary and the foi list [foiCategATTMdict, foiIDs]
                    write the file out [foiCategATTMdict, foiIDs]
                    foiCategATTMdict = {foi:[[A], [], []]}
                    A -> [[array(spectrum)], total]
                    array(spectrum) -> array([level1, level2, ...])
        '''
        CCF = np.asarray([123.5,	117.6,	114.6,	112.4,	108.2,	99.9,	93.7,	52.4])
        foiCategATTM = []
        for foi in self.foiIDs: # loop foi positions     
            print 'processing foi -> ', foi
            categATTM = []
            recIDs = self.crspfoi_simID[foi]
            for fn in propaFileName: # loop propagation manner
                calc = np.loadtxt(fn)
                meanspec = get_mean(calc, recIDs)
                if ('refl' in fn) or ('scat' in fn): 
                    [ATTMspectrum, ATTMtotal] = calc_ATTM(meanspec[2:10], CCF)
                else:
                    [ATTMspectrum, ATTMtotal] = calc_ATTM(meanspec[1:9], CCF)
                categATTM.append([ATTMspectrum, ATTMtotal])
            foiCategATTM.append(categATTM)
        foiCategATTMdict = {self.foiIDs[n]:foiCategATTM[n] for n in xrange(len(self.foiIDs))}
        f = open('for_map\\foi_ATTMdict_'+self.saveTag+'.pkl', 'wb')
        pickle.dump([foiCategATTMdict, self.foiIDs], f)
        f.close()
        return [foiCategATTMdict, self.foiIDs]
    
    def get_affected_receiver_id(self, reflResultFile):
        data = np.loadtxt(reflResultFile)
        IDs = []
        for n, d in enumerate(data):
            if not d[1::].any()==0:
                IDs.append(n)
        del data
        return IDs
            
    def calc_affect_receiv_save_pkl(self, propaFileName):
        ''' calculate the ATTM of the affected receiver by the temporary sources
            save out the results out to a pickle
            @propaFileName -> results file by all propagation manner
        '''
#        IDs = self.get_affected_receiver_id(propaFileName[1])
        IDs = map(int, np.loadtxt(propaFileName[1])[:,1]) # consider all the IDs as affected
        Mij = [np.loadtxt(propaFileName[0]), np.loadtxt(propaFileName[1]), np.loadtxt(propaFileName[2])]   
        print len(Mij[0]), len(Mij[1]), len(Mij[2])
        CCF = np.asarray([123.5,	117.6,	114.6,	112.4,	108.2,	99.9,	93.7,	52.4])
        # load direct field. This is big file so load it first and use it later
        affectedRcvATTM = []
        for r in IDs:
            print 'processing receiv ', r
            categATTM = []
            for j in xrange(len(propaFileName)): # loop propagation manner
                if j>0:
                    [ATTMspectrum, ATTMtotal] = calc_ATTM(Mij[j][r][2:10], CCF)  
                else:
                    [ATTMspectrum, ATTMtotal] = calc_ATTM(Mij[j][r][1:9], CCF)  
                categATTM.append([ATTMspectrum, ATTMtotal]) # leave the spectrum as empty to save some calculation
            affectedRcvATTM.append(categATTM)        
        print 'arrangin and writing out'
        receiverCategATTMdict = {ID:affectedRcvATTM[n] for n, ID in enumerate(IDs)}
        f = open('for_map\\receiver_categATTMdict_IDs_Leq_'+self.saveTag+'.pkl', 'wb')
        pickle.dump([receiverCategATTMdict, IDs], f)
        f.close()
        return receiverCategATTMdict



class PolylineShape():
    def __init__(self, polylineShapeFile):
        self.dbf = dbflib.open(polylineShapeFile)
        self.shp = shapelib.ShapeFile(polylineShapeFile)
        self.PLNUM = self.shp.info()[0]
        print 'shape info: -> ', self.shp.info()
    
    def write_point_shape_out(self, pointShapeFileName, pointsList):
        w2shp = shapelib.create(pointShapeFileName, shapelib.SHPT_POINT)
        w2dbf = dbflib.create(pointShapeFileName)
        w2dbf.add_field('ID', dbflib.FTInteger, 10, 0) # create 3 field for the ID and x, y coordinate
        w2dbf.add_field('x', dbflib.FTDouble, 16, 2)
        w2dbf.add_field('y', dbflib.FTDouble, 16, 2)
        i = 0
        for pts in pointsList:
            for pt in pts:
                shpObj = shapelib.SHPObject(shapelib.SHPT_POINT, i, [[pt]])
                w2shp.write_object(i, shpObj)
                w2dbf.write_record(i, {'ID':i})
                w2dbf.write_record(i, {'x':pt[0]})
                w2dbf.write_record(i, {'y':pt[1]})
                i += 1
        w2shp.close()
        w2dbf.close()
    
    def polyShape_to_points_by_eqDist(self, distance):
        allPoints = []
        for m in xrange(self.PLNUM):
            print 'processing polyline %d' %m
            # read shape object and feed function polyline_to_points_by_eqDist
            pobj = self.shp.read_object(m)
            polyline = pobj.vertices()[0]
            points = self.polyline_to_points_by_eqDist(polyline, distance)
            allPoints.append(points)
        self.shp.close()
        return allPoints # a list of points for different polyline segments as: [[(x1, y1), (x2, y2)...], [(x1, y1), (x2, y2)...]]
        
    def polyline_to_points_by_eqDist(self, polyline, distance):
        ''' Split the polyline to points by distance.
            polyline = [(x1, y1), (x2, y2), (x3, y3)...] is a list of points
            distance = number. The distance between two adjacent points. The
            distance is along the polyline NOT the direct distance from two points
        '''
        i = 0
        dni = 0.0 #distance from the inserted pt to the vertice of polyline
        points = [] # coordinate of the splitted points as (x,y)
        points.append(polyline[0][0:2])
        while i<len(polyline)-1:
            distij = dist2D(polyline[i], polyline[i+1])
            dni = dni+distij
            
            # find the inserting points
            ix = polyline[i][0]
            iy = polyline[i][1]
            if dni>=distance:
                ix = ix+(distance-(dni-distij))*(polyline[i+1][0]-polyline[i][0])/distij
                iy = iy+(distance-(dni-distij))*(polyline[i+1][1]-polyline[i][1])/distij
                points.append((ix, iy))
                ixiy2next = dist2D((ix, iy), polyline[i+1]) # updating the distance from the inserting points to the next vertice of the polyline
                while ixiy2next>distance:
                    ix = ix+distance*(polyline[i+1][0]-polyline[i][0])/distij
                    iy = iy+distance*(polyline[i+1][1]-polyline[i][1])/distij
                    points.append((ix, iy))
                    ixiy2next = dist2D((ix,iy), polyline[i+1])
                dni = ixiy2next # updating again
            i = i+1
        return points # point list as: [(x1, y1), (x2, y2), (x3, y3)...]
           

def split_polyline_to_points(lineShapeFileName, distance, newPointShapeFileName):
    ''' main function to split a polyline to points
        lineShapeFileName = shape file name of the poly line
        distance = distance between two points
        newPointShapeFileName = file name of the generated points
    '''
    pobj = PolylineShape(lineShapeFileName)
    pointlist = pobj.polyShape_to_points_by_eqDist(distance)
    pobj.write_point_shape_out(newPointShapeFileName, pointlist)


class PolyShape_PointShape(PolylineShape):
    def __init__(self, lineShapeFileName, newPointShapeFileName):
        PolylineShape.__init__(self,lineShapeFileName)
        self.IDlabel = 'ID'
        self.frStr = ['63', '125', '250', '500', '1K', '2K', '4K', '8K']
        self.nSHPFile = newPointShapeFileName
        self.den = ['DAY', 'EVE', 'NGT']
        self.fieldInfo = [(1, self.IDlabel, 16, 0)] # 1 -> integer
        fieldInfoDAY = [(2, 'LDAY_'+v, 16, 1) for v in self.frStr]  # 2 -> double
        fieldInfoEVE = [(2, 'LEVE_'+v, 16, 1) for v in self.frStr]  # 2 -> double
        fieldInfoNGT = [(2, 'LNGT_'+v, 16, 1) for v in self.frStr]  # 2 -> double
        self.newFieldNames = [fieldInfoDAY, fieldInfoEVE, fieldInfoNGT]
        self.powerFields = [['EDAY_'+v for v in self.frStr], \
                            ['EEVE_'+v for v in self.frStr], \
                            ['ENI_'+v for v in self.frStr]]
        
        
    def assign_power(self, ptDistance):
        # show the field info
        print 'filed in the polyline shape file: '
        for n in range(self.dbf.field_count()):
            print self.dbf.field_info(n)     
            
        for d in range(3): # day eveing and night
            newFieldName = self.newFieldNames[d]
            powerField = self.powerFields[d]
            den = self.den[d]
            ptDBF = dbflib.create(self.nSHPFile+'_'+den) 
            for fi in self.fieldInfo+newFieldName:
                ptDBF.add_field(fi[1], fi[0], fi[2], fi[3])
            ptSHP = shapelib.create(self.nSHPFile+'_'+den, shapelib.SHPT_POINT)
            b = 0
            for r in range(self.dbf.record_count()):
                pobj = self.shp.read_object(r)
                polyline = pobj.vertices()[0]
                rec = np.asarray([self.dbf.read_record(r)[fdn] for fdn in powerField])
                
                points = self.polyline_to_points_by_eqDist(polyline, ptDistance)
                try:
                    L7octave = 10.*np.log10(10**(0.1*rec)/len(points))
                except:
                    L7octave = np.asarry([0.0]*len(rec))            
                for pt in points:
                    print ('writing point %d (%0.1f, %0.1f)' %(b, pt[0], pt[1]))
                    shpObj = shapelib.SHPObject(shapelib.SHPT_POINT, b, [[pt]])
                    ptSHP.write_object(b, shpObj)
                    
                    ptDBF.write_attribute(b, 0, b)
                    for rc in range(len(newFieldName)):
                        ptDBF.write_attribute(b, rc+1, L7octave[rc] )
                    b += 1
            ptDBF.close()
            ptSHP.close()
    
    def form_record(self, ID, fieldName, L7octave):
            fieldnames = [self.IDlabel] + fieldName
            fieldvalues = [ID] + list(L7octave)
            record = {fieldnames[n]:fieldvalues[n] for n in xrange(len(fieldnames))}
            return record
        
def assign_lvl_from_polyLine_to_points(lineShapeFileName, ptDistance, newPointShapeFileName):
    ''' main function to split a polylie to points and assign the SPL
    '''
    ppobj = PolyShape_PointShape(lineShapeFileName, newPointShapeFileName)
    ppobj.assign_power(ptDistance)



  