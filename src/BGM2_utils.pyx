# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 15:57:22 2013
tools to calculate BGM2
@author: wgwei
"""
import shapelib, dbflib
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, sin, cos, atan, abs, log10, asin
from line_x_poly import line_x_poly, intersection

def point_inside_polygon(double x, double y, poly):
    cdef int i, n        
    cdef double plx, ply, p2x, p2y, xinters
    
    inside =False
    n = len(poly)
    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside    

cdef int point_inside_rectangle(double ptx, double pty, double xi, double yi,\
             double x2, double y2, double x3, double y3, double x4, double y4):
    if ptx<min(xi, x2, x3, x4):
        return 0
    elif ptx>max(xi, x2, x3, x4):
        return 0
    elif pty<min(yi, y2, y3, y4):
        return 0
    elif pty>max(yi, y2, y3, y4):
        return 0
    else:
        return 1
    
def get_rectangular_buffer(double pix, double piy, double p2x, double p2y, double distance=30.):
    cdef double xmin, xmax, ymin, ymax
    xmin = min(pix, p2x)-distance
    xmax = max(pix, p2x)+distance
    ymin = min(piy, p2y)-distance
    ymax = max(p2y, p2y)+distance   
    
    return xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax
    
def get_polygon_buffer(double pix, double piy, double p2x, double p2y, double distance=30.):
    ''' offset a distance to create a rectangular buffer
        return a list with all the points in count-clock direction
    '''
    cdef double xi, yi, x2, y2, theta, phi, deltax, deltay, xmin, xmax, ymin, ymax
    phi  = atan(abs((pix-p2x)/(piy-p2y)))
    theta = 3.14159/2-phi
    deltax = distance/sin(theta)
    deltay = distance/sin(phi)
    xmin = min(pix, p2x)-distance
    xmax = max(pix, p2x)+distance
    ymin = min(piy, p2y)-distance
    ymax = max(p2y, p2y)+distance    
        
    if (p2y-piy)/(p2x-pix)<0:        
        return [(xmin, ymax), (xmin, ymax-deltay), (xmax-deltax, ymin),(xmax, ymin), (xmax, ymin+deltay), (xmin+deltax, ymax), (xmin, ymax)]
    elif (p2y-piy)/(p2x-pix)>0: 
        return [(xmin, ymin), (xmin+deltax, ymin), (xmax, ymax-deltay),(xmax, ymax), (xmax-deltax, ymax), (xmin, ymin+deltay), (xmin, ymin)]
    else:
        return [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)]

def link_rcv_src_building(receiverObjs, sourceObjs, buildingObjs):
    ''' link every building to the receivers and sources IDs which is checked
        whether the centriod point of the building inside the rectangular buffer
        or not
    '''
    cdef int N, M, n, m, q, Q
    cdef double rx, ry, sx, sy
    N = len(receiverObjs)
    M = len(sourceObjs)
    Q = len(buildingObjs)
    linkRSB = [[[] for _ in xrange(M)] for _ in xrange(N)]
    for n in xrange(N):
        print 'linking receiver -> %d' %n
        rx = receiverObjs[n].vertix[0]
        ry = receiverObjs[n].vertix[1]
        for m in xrange(M):
            sx = sourceObjs[m].vertix[0]
            sy = sourceObjs[m].vertix[1]
            recBuffer = get_polygon_buffer(rx, ry, sx, sy, distance=50.)
            for q in xrange(Q):
                inBuffer = point_inside_polygon(buildingObjs[q].centroidx,buildingObjs[q].centroidy, recBuffer)
                if inBuffer:
                    linkRSB[n][n].append(buildingObjs[q].num)
    return linkRSB

cdef double dist_pt_to_line(double ptx, double pty, double xi, double yi, double x2, double y2):
    return abs((x2-xi)*(yi-pty)-(xi-ptx)*(y2-yi))/sqrt((x2-xi)**2.+(y2-yi)**2.)
    
def get_blocked_building_num(double sourcex, double sourcey, double receiverx, double receivery, buildingObjs, double distLim=50.0):
    buildingNum = []
    cdef double xmin, xmax, ymin, ymax, x, y
    xmin = min(sourcex, receiverx)-distLim
    xmax = max(sourcex, receiverx)+distLim
    ymin = min(sourcey, receivery)-distLim
    ymax = max(sourcey, receivery)+distLim   
    for bobj in buildingObjs:
        x, y = bobj.centroidx, bobj.centroidy
        distanceToLine = dist_pt_to_line(x, y, sourcex, sourcey, receiverx, receivery)
        if distanceToLine<distLim:
            inside = point_inside_rectangle(x,y,xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax)
            if inside:
                buildingNum.append(bobj.num)
    return buildingNum
    

def get_pti_pt2_ext_pt_by_dist(double ptix, double ptiy,double pt2x, double pt2y, double distance):
    ''' calculate the point on the extension line from pti to pt2
        at distance=distance
    '''
    cdef double disti2, theta
    disti2 = sqrt((ptix-pt2x)**2+(ptiy-pt2y)**2)
    theta = atan(abs((pt2y-ptiy)/(pt2x-ptix)))
    return pt2x+distance*cos(theta)*np.sign(pt2x-ptix), pt2y+distance*sin(theta)*np.sign(pt2y-ptiy)
    
def calc_intersectPt(int ptNum, polygonPtList):
    cdef np.ndarray polygonx = np.zeros(ptNum, dtype=np.double)
    cdef np.ndarray polygony = np.zeros(ptNum, dtype=np.double)
    cdef int m
    for m in xrange(ptNum):
        polygonx[m], polygony[m] = polygonPtList[m][0], polygonPtList[m][1]
    return polygonx, polygony

def dist2Dpyx(double pix, double piy, double p2x, double p2y):
    return sqrt((pix-p2x)**2+(piy-p2y)**2)
cdef double dist2D(double pix, double piy, double p2x, double p2y):
    return sqrt((pix-p2x)**2+(piy-p2y)**2)
 

def distance_to_building(double sourcex, double sourcey,  double receiverx, double receivery, buildingObjs, vailableIDs):
    ''' calculate the distance to the closest building '''
    cdef int n, m, M
    cdef double ds, dr, dsr, dsi, dri
    dsr = sqrt((sourcex-receiverx)**2+(sourcey-receivery)**2) # distance source to receiver
    ds = dsr # distance from source to building. Initialize it to a big number
    dr = dsr # distance from receiver to building
    
    # initialize it to a minus number. if it's not updated, we'll know there
    # is no obstacles between source and receiver
    recs = -1 #record number
    recr = -1
    for n in vailableIDs:
        bobj = buildingObjs[n]
        intersecPt = line_x_poly(sourcex,sourcey, receiverx, receivery, bobj.polygonx, bobj.polygony)
        if intersecPt: # not empty
            for isec in intersecPt:
                dsi = sqrt((sourcex-isec[0])**2+(sourcey-isec[1])**2)
                if dsi<ds:
                    recs = n
                ds = min(dsi, ds)                                     
                dri = sqrt((receiverx-isec[0])**2+(receivery-isec[1])**2)
                if dri<dr:
                    recr = n
                dr = min(dri, dr)
    return ds, dr, recs, recr
    
class Building():
    def __init__(self, relativeHeight, polygonPtList, int ID, int num):
        ''' relativeHeight =  a float/double number
            polygonPtList = [(x1, y1), (x2, y2), ....,(x1, y1)]
            ID = an integer
        '''
        cdef double x, y
        self.relativeHeight = relativeHeight
        self.polygonPtList = polygonPtList
        self.ID = ID
        ptNum = len(polygonPtList)
        self.polygonx, self.polygony = calc_intersectPt(ptNum, polygonPtList)
        x,y =np.mean(self.polygonx), np.mean(self.polygony)
        self.centroidx, self.centroidy = x, y
        self.num = num

def packBuildingToPKL(buildingShapeFile, IDrelativeHeight='RELHEIGHT', ID='ID'):
    cdef int N, n
    shp = shapelib.ShapeFile(buildingShapeFile)
    dbf = dbflib.open(buildingShapeFile)
    N = shp.info()[0]
    buildingObjs = []
    for n in xrange(N):
        sobj = shp.read_object(n)
        adbf = dbf.read_record(n)
        buildingObjs.append(Building(adbf[IDrelativeHeight], sobj.vertices()[0], adbf[ID], n))
    return buildingObjs
              
              
              
##### ######################### ######################### ################33
################   optimize the model part   ##########################333333   

cdef double LdiffOverLff_simplified(double waveLen, double rs, double rr, double thetas, double thetar, double w):
    cdef double pSquare, Ys, Yr, Ygr, Ysm
    L = rs + rr + w
    Ys = sqrt(3)*(cos(2.0/3.0*thetas)-0.5) * sqrt(2*rs*(w+rr)/(waveLen*L))
    Yr = sqrt(3)*(cos(2.0/3.0*thetar)-0.5) * sqrt(2*rr*(w+rs)/(waveLen*L))
    if Ys > Yr:
        Ygr = Ys   # greater Y
        Ysm = sqrt(3.0)*(cos(2.0/3.0*thetar)-0.5) * sqrt(2*rr*w/(waveLen*(w+rr)))   # smaller Y
    else:
        Ygr = Yr   # greater Y
        Ysm = sqrt(3.0)*(cos(2.0/3.0*thetas)-0.5) * sqrt(2*rs*w/(waveLen*(w+rs)))   # smaller Y
    pSquare = ((0.37/(Ygr+0.37))**2)*((0.37/(Ysm+0.37))**2)        
    return 10.0*log10(pSquare)
        
cdef np.ndarray flatRoof_level_NoGround(np.ndarray waveLens, double rs, double rr, double thetas, double thetar, double w):
    cdef np.ndarray LdNoGround = np.zeros(len(waveLens), dtype = np.double)
    for n, waveLen in enumerate(waveLens):
        LdNoGround[n]  = LdiffOverLff_simplified(waveLen, rs, rr, thetas, thetar, w)
    return LdNoGround
        
cdef np.ndarray flatRoof_level_ground(double sposx, double sposy, double rposx, double rposy, double N1pointx, double N1pointy, double N2pointx, double N2pointy, double phis, double phir, np.ndarray waveLens):
    cdef double rs, rr, ts, tr, sposImagx, sposImagy, rsImag, \
                tsImag, rposImagx,  rposImagy, rrImag, trImag, \
                D1, D2, D3, D4, w, v
    cdef np.ndarray LdNoGround1, LdNoGround2, LdNoGround3, LdNoGround4, LdGround
    cdef int n
    
    LdNoGround1 = np.zeros(len(waveLens), dtype = np.double)
    LdNoGround2 = LdNoGround1
    LdNoGround3 = LdNoGround1
    LdNoGround4 = LdNoGround1
    LdGround = LdNoGround1
    w = abs(N2pointx - N2pointx)
    rr = sqrt((rposx-N2pointx)**2.+(rposy-N2pointy)**2.)
    rs = sqrt((sposx-N1pointx)**2.+(sposy-N1pointy)**2.)
    ts = asin(abs(N1pointx-sposx)/rs)
    tr = asin(abs(rposx-N2pointx)/rr)
    
    sposImagx = sposx
    sposImagy = -sposy    
    rsImag = sqrt((sposImagx-N1pointx)**2.+(sposImagy-N1pointy)**2.)
    tsImag = asin(abs(sposImagx-N1pointx)/rsImag)
    rposImagx = rposx
    rposImagy = -rposy
    rrImag = sqrt((rposImagx-N2pointx)**2.+(rposImagy-N2pointy)**2.)
    trImag = asin(abs(rposx-N2pointx)/rrImag)
    
    LdNoGround1 = flatRoof_level_NoGround(waveLens, rs, rr, ts, tr, w)
    D1 = rs + rr + w
    D2 = rsImag + rr + w
    LdNoGround2 = flatRoof_level_NoGround(waveLens, rsImag, rr, tsImag, tr, w) + 20.*log10(D1/D2)
    D3 = rs + rrImag + w
    LdNoGround3 = flatRoof_level_NoGround(waveLens, rs, rrImag, ts, trImag, w) + 20.*log10(D1/D3)
    D4 = rsImag + rrImag + w
    LdNoGround4 = flatRoof_level_NoGround(waveLens, rsImag, rrImag, tsImag, trImag, w) + 20.*log10(D1/D4)
    
    for n in xrange(len(waveLens)):
        v = 10.*log10(10.0**(0.1*LdNoGround1[n]) + 10.0**(0.1*LdNoGround2[n]) + 10.0**(0.1*LdNoGround3[n]) + 10.0**(0.1*LdNoGround4[n]))
        LdGround[n] = v
    return LdGround
    
def fitModel_refl(double h1, double h2, double Hs, double Hr, double Wi, \
            double sHi, double rHi, double distSR,\
             double N1pointx, double N1pointy, double N2pointx, double N2pointy,\
             double Ws, double Wr, double rs, double rr, double sposx, double sposy, \
             double rposx, double rposy, double phis, double phir, np.ndarray waveLength,):
    """ h1 = Hi-sHi
        h2 = Hi-rHi
        sHi = source height
        rHi = receiver height
        srDist is the distance between source and receiver
        N1point is building vertix in source canyon (x, y)
        N2point is building vertex in receiver canyon (x, y)
    """       
    cdef int N = len(waveLength)
    cdef double C, AroofCan, Ainter
    cdef np.ndarray C1s, C1r, C3s, C3r, Lhs, Lhr, AcanE, Abar
    Lhs = -40.*np.ones(N, dtype=np.double)
    Lhr = Lhs
    C1s = np.zeros(N, dtype=np.double)
    C1r = C1s
    C3s = C1s
    C3r = C1s
    AcanE = C1s
    if sHi==0:
        sHi = 0.01
    if rHi==0:
        rHi = 0.01
    if N1pointy!=N2pointy:
        N1pointy = (N1pointy+N2pointy)/2.0
        N2pointy = N1pointy
    
    Abar = flatRoof_level_ground(0.0, sHi, rposx, rposy, N1pointx, N1pointy, N2pointx, N2pointy, phis, phir, waveLength)
    
    # reflection coefficient alpha=0.97, beta=0.97
    C = 1.5*Ws+Wi+1.5*Wr
    C1s = 1./(3.31*cos(phir)*np.sqrt(rr/waveLength)+1.)**2.
    C1r = 1./(3.31*cos(phis)*np.sqrt(rs/waveLength)+1.)**2.
    C3s = 3.31*h1*np.sqrt(Wi/waveLength)+0.5*Ws+rr+Wi
    C3r = 3.31*h2*np.sqrt(Wi/waveLength)+0.5*Wr+rs+Wi
        
    if (Hs-sHi)/h1>=1.:
        Lhs = np.zeros(N, dtype=np.double)
    if (Hs-sHi)/h1>0.333333 : # 1./3.
        Lhs = -6.17*(1-(Hs-sHi)/h1)*(1-1.37*np.log10(np.sqrt(waveLength*Ws)/Wi))
        if Lhs.any()>0:
            Lhs = np.zeros(N, dtype=np.double)

    if (Hr-rHi)/h2>=1.:
        Lhr = np.zeros(N, dtype=np.double)
    if (Hr-rHi)/h2>0.333333 :
        Lhr = -6.17*(1-(Hr-rHi)/h2)*(1-1.37*np.log10(np.sqrt(waveLength*Wr)/Wi))
        if Lhr.any()>0:
            Lhr = np.zeros(N, dtype=np.double)

    AcanE = 10.5*np.log10(12.46*C1s* 0.97**2. * distSR**2. * 10.**(0.1*Lhs)/((C3s+Ws)**2.)\
        + 22.24*C1r* 0.97**2. * distSR**2. * 10.**(0.1*Lhr)/((C3r+Wr)**2.)\
        + 0.05*0.97**2.*0.97**2. * distSR**2. * 10.**(0.1*Lhs)*10.**(0.1*Lhr) / ((3.31*h1/np.sqrt(waveLength)+C)*(3.31*h2/np.sqrt(waveLength)+C)))
#========================================================================================================           

    if Hs!=0 or Hr!=0:
        AroofNocan = 0.0
        if Hs!=0 and Hr!=0:
            AroofCan = 5.
        else:
            AroofCan = 2.5
    else:
        AroofNocan = 0.27*Abar+2.9
        AroofCan = 0.
    Ainter = distSR/100.0
    if Ainter>5.0:
        Ainter = 5.0
    return [Abar, AcanE, Ainter, AroofNocan, AroofCan] 
   

def cal_ds_dr_h(double v1x, double v1y, double v2x, double v2y, double v3x, double v3y, double phis, double phir, double sht, double rht):
    ''' calculate ds and dr in Jens's scattering model
        v1 is point 1 (x, y) usually is source position
        v2 is the intersection (x, y)
        v3 is the receiver position (x,y)
        phis and phir are angles defined in pierce defraction equation
    '''    
    cdef double distV1V2, distV2V3, distV1V3, angle1, angle2, angle123, h, dr, ds
    distV1V2 = sqrt((v2y-v1y)**2+(v2x-v1x)**2)
    distV2V3 = sqrt((v3y-v2y)**2+(v3x-v2x)**2)
    distV1V3 = sqrt((v3y-v1y)**2+(v3x-v1x)**2)
    try:
        angle1 = atan(abs((v2x-v1x)/(v2y-v1y)))
    except ZeroDivisionError:
        angle1 = np.pi/2-0.1 # assign a angle close to pi/2
    try:
        angle2 = atan(abs((v3x-v2x)/(v3y-v2y)))
    except ZeroDivisionError:      
        angle2 = np.pi/2-0.1 # assign a angle close to pi/2
    angle123 = angle1+angle2
    h = distV1V2*distV2V3*sin(angle123)/distV1V3    
    dr = sqrt(abs(distV2V3**2-h**2))
    ds = sqrt(abs(distV1V2**2-h**2))
    if h<=0 or h!=h or ds<=0 or ds!=ds or dr<=0 or dr!=dr:
        assert(False)
    return [ds, dr, h]   
     
   
def fitModel_scatter(double Cvsq, double Ctsq, double sht, double rht, np.ndarray fr,\
        double srDist, double p_sourcex, double p_sourcey, double N1pointx, double N1pointy,\
        double p_receiverx, double p_receivery, double N2pointx, double N2pointy,\
         double phis, double phir, double Ws, double Wr, double height):
    ''' test Jens' scattering model    
    '''
    cdef int n, N
    cdef double H0, d0, f0, epslon, tb1, tb2, tb3, vb1, vb2, vb3, ds, dr, h, ct, cv, gema2, AscatCanyon
    cdef np.ndarray AscatBarCt, AscatBarCv, AscatBar, Ascat, mixAttenuation
    N = len(fr)
    AscatBarCt = np.zeros(N, dtype = np.double)
    AscatBarCv = AscatBarCt
    AscatBar = AscatBarCt
    Ascat = AscatBarCt
    H0 = 10.0
    d0 = 10.0
    f0 = 1000.0
    epslon = 0.000004
    tb1 = -49.6+10*log10(Ctsq)
    tb2 = 11.5
    tb3 = -13.1
    vb1 = -52.8+10*log10(Cvsq)
    vb2 = 11.3
    vb3 = -17.1    
    V2 = intersection(p_sourcex, p_sourcey, N1pointx, N1pointy, p_receiverx, p_receivery, N2pointx, N2pointy)
    [ds, dr, h] = cal_ds_dr_h(p_sourcex, p_sourcey, V2[0], V2[1], p_receiverx, p_receivery, phis, phir, sht, rht)    
    for n in xrange(N):
        ct = tb1 + tb2*log10(h/d0) + tb3*log10((h**2)/(dr*ds)+(h**2)*epslon) + 3.33*log10(fr[n]/f0)      
        cv = vb1 + vb2*log10(h/d0) + vb3*log10((h**2)/(dr*ds)+(h**2)*epslon) + 3.33*log10(fr[n]/f0)
        AscatBar[n] = 10*log10(10**(0.1*ct) + 10**(0.1*cv))   
    if height==0:
        print 'height: ', height
        assert(False)        
    if Ws==False and Wr==False:
        AscatCanyon = 0.0
    elif Ws==False and Wr!=False:
        gema2 = 2.0*height/Wr
        AscatCanyon = 7.0 + gema2*log10(height/H0)
    elif Ws!=False and Wr==False:
        gema2 = 2.0*height/Ws
        AscatCanyon = 7.0 + gema2*log10(height/H0)
    elif Ws!=False and Wr!=False:
        # should be careful: when both Ws and Wr are small, but height is big. then the
        # canyon term can become very big!!! That could be the main reason why this model
        # could overestimate the level when implimment it in the GIS system
        gema2 = 2.0*height*(1.0/Ws+1.0/Wr)
        AscatCanyon = 14.0 + gema2*log10(height/H0)        
    Ascat = AscatBar+AscatCanyon     
    mixAttenuation =  Ascat    
    return mixAttenuation     



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
