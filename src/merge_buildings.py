import shapelib, dbflib
import shutil
import numpy as np

def dist2D(pti, pt2):
    return np.sqrt((pti[0]-pt2[0])**2.+(pti[1]-pt2[1])**2.)
    
class DonutsBuildings():
    def __init__(self, buildingSHPFile):
        self.buildingSHPFile = buildingSHPFile
        self.shp = shapelib.ShapeFile(buildingSHPFile)
        self.shpNum, self.shpType = self.shp.info()[0], self.shp.info()[1]
        
    def length_of_polygon(self, vertices):
        vts = vertices[0]
        if vts[0]==vts[-1]:
            vts.pop()
        dists = [dist2D(vts[i], vts[i+1]) for i in range(len(vts)-1)]
        return sum(dists)
        
    def write_shp(self, newShape):            
        tempshp = shapelib.create('temp', self.shpType)
        for b in range(self.shpNum):
            print ('\t merging -> %d' %b)
            shpobj = self.shp.read_object(b) #(5, i, [[(), (),()], []])
            if len(shpobj.vertices())>1:
                shpobj = shapelib.SHPObject(self.shpType, b, [shpobj.vertices()[0]])
            tempshp.write_object(b, shpobj)
        tempshp.close()
        shutil.copy(self.buildingSHPFile+'.dbf', 'temp'+'.dbf')
        tempshp.close()        
        
        self.clear_small_polygons('temp', newShape, 20.)
    
    def clear_small_polygons(self, shapeFile, newShapeFile, lenLimit):
        shp = shapelib.ShapeFile(shapeFile)
        dbf = dbflib.open(shapeFile)
        newSHP = shapelib.create(newShapeFile, self.shpType)
        newDBF = dbflib.create(newShapeFile)
        for f in range(dbf.field_count()):
            fi = dbf.field_info(f)
            newDBF.add_field(fi[1], fi[0], fi[2], fi[3])
        bb = 0
        for b in range(shp.info()[0]):
            sobj = shp.read_object(b)
            rec = dbf.read_record(b)
            if self.length_of_polygon(sobj.vertices()) > lenLimit:
                shpobj = shapelib.SHPObject(self.shpType, bb, [sobj.vertices()[0]])
                newSHP.write_object(bb, shpobj)
                newDBF.write_record(bb, rec)
                bb += 1
        shp.close()
        dbf.close()
        newSHP.close()
        newDBF.close()
            
if __name__=='__main__':
    bld = DonutsBuildings('test-building')
bld.write_shp('test-merge')