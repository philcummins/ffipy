#!/usr/bin/env python
#import dcs_python
#import kiklib
import math as mt
import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 

def kiksyn(nt,dt,ib,ic,str,dip,rak,z,az,p,g,t12,
           velmod_filename=None,pad=False,verbose=False):
    """
        Wrapper for Kikuchi body wave synthetic routine
    Input:
    nt,dt            length and sampling of output time series
    ib,ic            ib=1 for P; ic=1,2,3 for Z,NS,EW, resp.
                     ib=2 for SV; ic=1,2 for Z,radial, resp.
                     ib=3 for SH; (ic ignored)
    p, g             slowness and geometrical spreading (from raypgeom)
    z, az            event depth and source-to-recever azimuth
    str,dip,rak      fault geometry
    t12              half-duration of symmetrical trangle stf
    velmod_filename  velocity model file name
    Output
    syn              real array of length nt
    """
# From Hong Kie's ffi program:
#      read(5,*) ms,t1,t2
#        read(1,*) az(js),az2(js),del(js),p(js),g(js),ix0(js)
#        read(1,*) im(js),ib(js),ic(js)
#     -- read velocity model --------------------------------------------
#      read(2,'(a40)') dsn
#      read(2,*) tqp,tqs,nl,(vp(l),vs(l),den(l),dep(l),l=1,nl)
#      read(2,*) nl1,(vp1(l),vs1(l),den1(l),dep1(l),l=1,nl1)
#      read(2,*) nl2,(vp2(l),vs2(l),den2(l),dep2(l),l=1,nl2)
    
    lines = open(velmod_filename,'r').readlines()
    buf = []
    for line in lines[1:]: buf.extend(line.split())
    tqp,tqs = [float(x) for x in buf[0:2]]
    # Pack velocity model
    nl = int(buf[2])
    nl1 = int(buf[nl*4+3])
    nl2 = int(buf[(nl+nl1+1)*4])
    vmod=np.zeros((4,(nl+nl1+nl2)))
    vmod[0:4,       0:nl      ] = np.array([float(x) for x in buf[3:nl*4+3]]).reshape((nl,4)).transpose()
    vmod[0:4,    nl:nl+nl1    ] = np.array([float(x) for x in buf[4*(nl+1):4*(nl+nl1+1)]]).reshape((nl1,4)).transpose()
    vmod[0:4,nl+nl1:nl+nl1+nl2] = np.array([float(x) for x in buf[4*(nl+nl1+1)+1:4*(nl+nl1+nl2+1)+1]]).reshape((nl2,4)).transpose()
    # Pad nt to power of 2 if needed
    if mt.log(nt,2) % 1 != 0.:
        nt1 = int(mt.pow(2,int(mt.log(nt,2))+1))
        if verbose: print 'nt=%d raised to next power of two (%d)' % (nt,nt1)
    else:
        nt1 = nt
    if pad:
        nt1 *=2
    # Now call synthetic generating routine
    syn = kiklib.kiksyn(nt1,dt,ib,ic,az,p,g,str,dip,rak,z,t12,
                        tqp,tqs,vmod,[nl,nl1,nl2],velmod_filename!=None)
    return syn[0:nt]

def raypgeom(edep,ir0,ip):
    """
    Wrapper for Kikuchi raypath geometry routine
    """
    gfact,pp0,vp0 = kiklib.raypgeom(edep,ir0,ip)

    return gfact,pp0,vp0

def fault_disp(elon,elat,edep,strk,dip,lnt,wdt,disl1,disl2,
               rlon,rlat,alp=0.5,utm=False,dmax=-1.):
    """
    Wrapper around fault_slip from dcs_subs
    """
    edsp,ndsp,zdsp = dcs_python.fault_disp(alp,elon,elat,edep,strk,dip,\
                                           lnt,wdt,disl1,disl2,rlon,rlat,utm,dmax)
    return edsp,ndsp,zdsp

#!/usr/bin/env python


def stdaz(slat,slon,rlat,rlon,unit):

    pib2 = 0.5*mt.pi
    ecc  = .003367
    re   = 6378.388
    ecc=0.
    re=6371.
    ec1  = (1.-ecc)*(1.-ecc)
    #************ Calculate epicentral distances, azimuths, etc. ************
    aa  = ec1*mt.sin(slat*DEGR)     # Calculate event geocentric coordinates
    bb  =  mt.cos(slat*DEGR)
    glats = mt.atan2(aa, bb)
    d1  =  mt.sin(glats)
    sps = d1*d1
    d1  =  mt.cos(glats)
    cps = d1*d1
    rs  = re*mt.sqrt(ec1/((1.-ecc*cps)*(1.-ecc*cps) + ecc*ecc*sps*cps))
    trs = pib2 - glats
    prs = slon * DEGR
    ass  =  mt.sin(trs) * mt.cos(prs)    # Calculate direction cosines for source
    bs  =  mt.sin(trs) * mt.sin(prs)
    cs  =  mt.cos(trs)
    ds  =  mt.sin(prs)
    es  = -mt.cos(prs)
    gs  =  mt.cos(trs) * mt.cos(prs)
    hs  =  mt.cos(trs) * mt.sin(prs)
    os  = -mt.sin(trs)
    # Calculate station geocentric coordinates
    glatr = mt.atan2(ec1*mt.sin(rlat*DEGR) , mt.cos(rlat*DEGR))
    d1  = mt.sin(glatr)
    spr = d1 * d1
    d1  = mt.cos(glatr)
    cpr = d1 * d1
    rr  = re*mt.sqrt(ec1/((1.-ecc*cpr)*(1.-ecc*cpr) + ecc*ecc*spr*cpr))
    trr = pib2 - glatr
    prr = rlon * DEGR
    ar  =  mt.sin(trr)*mt.cos(prr)    # Calculate direction cosines for receiver
    br  =  mt.sin(trr)*mt.sin(prr)
    cr  =  mt.cos(trr)
    dr  =  mt.sin(prr)
    er  = -mt.cos(prr)
    gr  =  mt.cos(trr) * mt.cos(prr)
    hr  =  mt.cos(trr) * mt.sin(prr)
    pr  = -mt.sin(trr)
    cosdr  = ass * ar + bs * br + cs * cr      # Calculate distance
    deltar = mt.acos(cosdr)                                    # in radians 
    deltak = deltar * .5 * (rr + rs)                      # in kilometers
    delta  = deltar / DEGR                                 # and degrees 
    szs = ds * ar + es * br                      # Calculate azimuth 
    czs = gs * ar + hs * br + os * cr
    szr = dr * ass + er * bs
    czr = gr * ass + hr * bs + pr * cs
    # azima - azimuth to source from station 
    # bazim - backazimuth to station from source
    # cazim - azimuth of wavefront at array
    if szr == 0.:
        azima = 0. 
        bazim = 180. 
    else:
        bazim = mt.atan2(-szs, -czs) / DEGR
        azima = mt.atan2(-szr, -czr) / DEGR
    if bazim <   0.:
        bazim = bazim+360.
    if azima <   0.:
        azima = azima+360.
    cazim = azima + 180.
    if cazim > 360.:
        cazim = cazim-360.
    # Return station distance and azimuth for any value of iop 
    dist = delta
    if unit == 'km':
        dist = deltak
    elif unit == 'deg':
        dist = delta
    elif unit == 'rad':
        dist = deltar
    return(dist, bazim, azima)


def project_b_ac(alon,alat,blon,blat,clon,clat):
    global DEGR
    # Projection of point b onto the line a-c, where a-b-c are taken
    # (counter?)clockwise
    dip = 0.

    #  (ac,azac,bazim) = stdaz(alat,alon,clat,clon)
    csac = mt.cos((90.-alat)*DEGR)*mt.cos((90.-clat)*DEGR) + \
           mt.sin((90.-alat)*DEGR)*mt.sin((90.-clat)*DEGR)* \
           mt.cos((alon-clon)*DEGR)
    ac   = mt.acos(csac)/DEGR
    azac = -mt.asin(mt.sin((alon-clon)*DEGR)*mt.sin((90.-clat)*DEGR)/ \
                 mt.sin(ac*DEGR))/DEGR
    bd = dist_b_ac(alon,alat,blon,blat,clon,clat)	
    if bd >= 0.:
        angl = azac+90.
    else:
        angl = azac-90.
    if angl < -360:
        angl = angl+360
    if angl >  360:
        angl = angl-360
    (dlon,dlat) = prj_point(blon,blat,abs(bd),angl)
    dist = dist_b_ac(alon,alat,dlon,dlat,clon,clat)	
#    dz = bd*mt.tan(dip*DEGR)*DEGK
    return(dlon,dlat)

def prj_point(lon,lat,dist,azim):
    global DEGR
    # Projection of point (lon,lat) dist degrees along an azimuth
    # Returns (lon,lat) of projected point.

    if azim == 0.:
        azm = 0
        delt = 90.-lat-dist
    else:
        csa = mt.cos(dist*DEGR)*mt.cos((90.-lat)*DEGR) + \
              mt.sin(dist*DEGR)*mt.sin((90.-lat)*DEGR)* \
              mt.cos(azim*DEGR)
        delt = mt.acos(csa)/DEGR
        azm = -mt.asin(mt.sin(azim*DEGR)*mt.sin(dist*DEGR)/ \
                    mt.sin(delt*DEGR))/DEGR
    lat = 90.-delt  
    lon = lon-azm
    if lon < -180:
        lon = lon+360
    if lon >  180:
        lon = lon-360
    return(lon,lat)


def dist_b_ac(alon,alat,blon,blat,clon,clat):
    global DEGR
    # Distance (in deg) of point b from line a-c, where a-b-c 
    # are taken clockwise
    csab = mt.cos((90.-alat)*DEGR)*mt.cos((90.-blat)*DEGR) +\
           mt.sin((90.-alat)*DEGR)*mt.sin((90.-blat)*DEGR)* \
           mt.cos((alon-blon)*DEGR)
    ab   = mt.acos(csab)/DEGR
    azab = -mt.asin(mt.sin((alon-blon)*DEGR)*mt.sin((90.-blat)*DEGR)/ \
                 mt.sin(ab*DEGR))/DEGR
    csac = mt.cos((90.-alat)*DEGR)*mt.cos((90.-clat)*DEGR) + \
           mt.sin((90.-alat)*DEGR)*mt.sin((90.-clat)*DEGR)* \
           mt.cos((alon-clon)*DEGR)
    ac   = mt.acos(csac)/DEGR
    azac = -mt.asin(mt.sin((alon-clon)*DEGR)*mt.sin((90.-clat)*DEGR)/ \
                 mt.sin(ac*DEGR))/DEGR
    A = azac - azab            # angle between a-c and a-b
    dist = mt.asin(mt.sin(A*DEGR)*mt.sin(ab*DEGR))/DEGR
    return(dist)

def dist_lonlat(lon1,lat1,lon2,lat2):
    global DEGR
# Use Haversine formula?
    a = mt.sin(0.5*(lat1-lat2)*DEGR)
    b = mt.sin(0.5*(lon1-lon2)*DEGR)
    dst0 = 2.*mt.asin(mt.sqrt(a*a+b*b*mt.cos(lat1*DEGR)*mt.cos(lat2*DEGR)))/DEGR
    dist = mt.acos(mt.cos((90.-lat1)*DEGR)*mt.cos((90.-lat2)*DEGR) +
                mt.sin((90.-lat1)*DEGR)*mt.sin((90.-lat2)*DEGR)*
                mt.cos((lon1-lon2)*DEGR))/DEGR
#    print '??? ',dst0,dist
    return (dist)



D2K  = 6371.*mt.pi/180.
DEGR = mt.pi/180.
Verbose = False

def grid(etg0,hypo_lat,hypo_lon,hypo_dep,cntr_lat,cntr_lon,strike,dip,rake,\
         length,width,nsfx,nsfy):
    '''Caluclat a lat,lon,depth grid for subdfualts on a rectangular fault plane.
    etg0                          an (arbitrary) event tag
    hypo_lat,hypo_lon,hypo_dep    Hypocenterl information
    cntr_lat,cntr_lon             Centroid lat,lon
    strike,dip,rake               Fault orientation (usually from cmt)
    length,width                  Total fault length and width
    nsfx,nsfy                     No. of subfaults along strike (x) and down dip (y)'''

    (dist, bazim, azima) = stdaz(hypo_lat,hypo_lon,cntr_lat,cntr_lon,"deg")
    if Verbose:
        print "subfaults:\n\thypo_lat,hypo_lon = %g,%g, cntr_lat,cntr_lon=%g,%g" % \
              (hypo_lat,hypo_lon,cntr_lat,cntr_lon)
        print "\tCentroid %7.2f km at strike %6.2f deg from hypocenter" % \
              (dist*D2K,bazim)
    # bazim is the direction of the centroid from the hypocenter
    if mt.sin((strike-bazim)*DEGR) < 0:
        (tmp_lon,tmp_lat) = prj_point(hypo_lon,hypo_lat,length/D2K,strike)
    else:
        (tmp_lon,tmp_lat) = prj_point(hypo_lon,hypo_lat,length/D2K,strike+180.)

    (tmp_lon,tmp_lat) = project_b_ac(hypo_lon,hypo_lat,cntr_lon,cntr_lat,\
				     tmp_lon,tmp_lat)
    # (tmp_lon,tmp_lat) is the projection of the centroid along strike from hypo
    (distx) = dist_lonlat(hypo_lon,hypo_lat,tmp_lon,tmp_lat)
    if mt.cos((strike-bazim)*DEGR) < 0: distx = -distx
    if Verbose:
        print "\tAlong-strike epicentre-centroid distance is %7.2f km" % (distx*D2K)

    (disty) = dist_lonlat(cntr_lon,cntr_lat,tmp_lon,tmp_lat)
    if mt.sin((strike-bazim)*DEGR) < 0: disty = -disty
    if Verbose:
        print "\tUp-dip epicentre-centroid distance is %7.2f km" % (disty*D2K)

    e_strk_x  = 0.5*length-distx*D2K
    e_dip_y   = 0.5*width-disty*D2K/mt.cos(dip*DEGR)
    top_depth = hypo_dep - (width-e_dip_y)*mt.sin(dip*DEGR)
    if top_depth < 0:
        # The shallowest would keep the centre of subfault below 4.0
        print 'Adjusting top of fault from %g to 4.0 km depth' % top_depth
        e_dip_y = e_dip_y + (2.0-top_depth)/mt.sin(dip*DEGR)
    if Verbose:
        print "\tEpicenter in fault coords is (%5.2f,%5.2f)" % (e_strk_x,e_dip_y)

    ffi_io = True
    gmt_io = True
    # Open output files
    if ffi_io: ffi_file = open('subfaults.grid','w')
    if gmt_io: gmt_file = open(etg0+'_subfaults.xyz','w')

    cosdp = mt.cos(dip*DEGR)
    # Project origin of fault coordinate system
    (olon,olat) = prj_point(hypo_lon,hypo_lat,-e_strk_x/D2K,strike)
    (olon,olat) = prj_point(olon,olat,e_dip_y*cosdp/D2K,90.+strike)
    odep = hypo_dep + e_dip_y*mt.sin(dip*DEGR)
    
    dsfx  = length/nsfx
    dsfy  = width/nsfy
    nptx = nsfx
    npty = nsfy
    iex= int (e_strk_x/dsfx) + 1
    iey= nsfy - int (e_dip_y/dsfy)
    source_subfault = (iex-1)*nsfy + iey
    if Verbose:
        print "\tolon=%g, olat=%g, odep=%g" % (olon,olat,odep)
        print "\tnsfx=%d, dsfx=%g, nsfy=%d, dsfy=%g\n" % (nsfx,dsfx,nsfy,dsfy)
        print "\t%d,%d Source is in subfault %d\n" % (iex,iey,source_subfault)
    if iex <= 0 or iex > nsfx or iey > nsfy or iey <= 0:
        print 'Epicenter in subfault %d,%d - grid() failed!\n' % (iex,iey)
        return np.zeros((0))
    if ffi_io:
        ffi_file.write("%d  0                   # hypocenter subfault, tapering\n" %\
                       source_subfault)
    fault_grid = np.zeros((nsfx*nsfy,7))
    isf = 0
    for ix in range(1,nsfx+1):
        sfx=(ix-1)*dsfx
        for iy in range(1,nsfy+1):
            sfy=width-iy*dsfy
            if gmt_io: gmt_file.write(">>\n")
            titl = "subfault ix=%d,iy=%d, npts = %dx%d" % (ix,iy,nptx,npty)
            # GMT output
            (rlon,rlat) = prj_point(olon,olat,sfy*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,sfx/D2K,strike)
            dep = odep - sfy*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            lt0 = lat
            ln0 = lon
            dp0 = dep
            (rlon,rlat) = prj_point(olon,olat,(sfy+dsfy)*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,sfx/D2K,strike)
            dep = odep - (sfy+dsfy)*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            (rlon,rlat) = prj_point(olon,olat,(sfy+dsfy)*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,(sfx+dsfx)/D2K,strike)
            dep = odep - (sfy+dsfy)*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            (rlon,rlat) = prj_point(olon,olat,sfy*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,(sfx+dsfx)/D2K,strike)
            dep = odep - sfy*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (ln0,lt0,dp0))
            #
            (rlon,rlat) = prj_point(olon,olat,(sfy+0.5*dsfy)*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,(sfx+0.5*dsfx)/D2K,strike)
            dep = odep - (sfy+0.5*dsfy)*mt.sin(dip*DEGR)
            fault_grid[isf,0:8] = [lon,lat,dep,-(360.-strike),dip,rake,dsfx*dsfy*1.e6]
            isf += 1
            if ffi_io: ffi_file.write("%9.4f %9.4f %7.2f %6.1f %6.1f %6.1f %8.2f 0\n" % \
                                      (lon,lat,dep,-(360.-strike),dip,rake,dsfx*dsfy))
        
    if ffi_io: ffi_file.close()
    if gmt_io: gmt_file.close()
    return fault_grid,sbfs


def fault_grid(etg0,hypo_lat,hypo_lon,hypo_dep,xdrct,h_top,strike,dip,rake,\
    length,width,nsfx,nsfy,Verbose=False,ffi_io=True,gmt_io=True):
    '''
    Caluclate a lat,lon,depth grid for subdfaults on a rectangular fault plane.
    etg0                          an (arbitrary) event tag
    hypo_lat,hypo_lon,hypo_dep    Hypocenterl information
    xdrct                         Rupture directivity: -1 for ani-strike,
                                                        0 for bilateral
                                                        1 for along strike
    h_top                          Depth to top of fault
    strike,dip,rake               Fault orientation (usually from cmt)
    length,width                  Total fault length and width
    nsfx,nsfy                     No. of subfaults along strike (x) and down dip (y)
    gmt_io                        A file of grid coordinates in gmt formnat will be written
    ffi_io                        A file of grid parameters in ffi formnat will be written
    Returns
          (grid,sbfs), where
                    grid: A numpy array containing center lon,lat,dep,-(360.-strike),dip,rake,
                          dsfx*dsfy (in km^2) for each subfault.
                    sbfs: A list of 4-element lists for each subfault, with each of the 4 elements
                          consisting of float triplets containing (lon,lat,depth) for each subfault
                          corner
      *** len(sbfs) == 0 indicates error return (E.g., hypo_dep lies above or below the fault) ***
'''
    grid = []
    sbfs = []

    # Sanity check on event depth
    h_bottom = h_top+width*mt.sin(dip*DEGR)
    if hypo_dep <= h_top:
        print 'fault_grid: Hypocentre depth (%g) above top of fault (%g) - Failed!' % (hypo_dep,h_top)
        return np.array(grid),np.array(sbfs)
    elif hypo_dep >= h_bottom:
        print 'fault_grid: Hypocentre depth (%g) below bottom of fault (%g) - Failed!' % (hypo_dep,h_bottom)
        return np.array(grid),np.array(sbfs)
    # Caluclate fault coordinates of epicentre
    dx = length/float(nsfx)
    xdist = 0.5*(length - xdrct*(length-dx))
    ydist = (h_bottom-hypo_dep)/mt.sin(dip*DEGR)
    # Calculate fault centroid
    cosdp = mt.cos(dip*DEGR)
    (tmp_lon,tmp_lat) = prj_point(hypo_lon,hypo_lat,(0.5*width-ydist)*cosdp/D2K,strike-90.)
    (cntr_lon,cntr_lat) = prj_point(tmp_lon,tmp_lat,(0.5*length-xdist)/D2K,strike)
    # Project origin of fault coordinate system
    (olon,olat) = prj_point(cntr_lon,cntr_lat,-0.5*length/D2K,strike)
    (olon,olat) = prj_point(olon,olat,0.5*width*cosdp/D2K,90.+strike)
    odep = h_bottom
    
    # Open output files
    if ffi_io: ffi_file = open(etg0+'_subfaults.grid','w')
    if gmt_io: gmt_file = open(etg0+'_subfaults.xyz','w')
    
    dsfx  = length/nsfx
    dsfy  = width/nsfy
    nptx = nsfx
    npty = nsfy
    iex= int (xdist/dsfx) + 1
    iey= nsfy - int (ydist/dsfy)
    source_subfault = (iex-1)*nsfy + iey
    if Verbose:
        print "fault_grid:\n\thypo_lat,hypo_lon,hypo_dep\t=%7.2f,%7.2f,%g\n\tcntr_lat,cntr_lon\t\t=%7.2f,%7.2f" % \
              (hypo_lat,hypo_lon,hypo_dep,cntr_lat,cntr_lon)
        print "\tolat,olon,odep\t\t\t=%7.2f,%7.2f,%6.2f" % (olat,olon,odep)
        print "\tnsfx=%d, dsfx=%g, nsfy=%d, dsfy=%g\n" % (nsfx,dsfx,nsfy,dsfy)
        print '\tepicenter x,y = %g,%g (%d,%d)' % (xdist,ydist,iex,iey)
        print "\tSource is in subfault %d\n" % (source_subfault)
    if iex <= 0 or iex > nsfx or iey > nsfy or iey <= 0:
        print 'Epicenter in subfault %d,%d - grid() failed!\n' % (iex,iey)
        return np.array(grid),np.array(sbfs)
    if ffi_io:
        ffi_file.write("%d  0                   # hypocenter subfault, tapering\n" %\
                       source_subfault)
    isf = 0
    for ix in range(1,nsfx+1):
        sfx=(ix-1)*dsfx
        for iy in range(1,nsfy+1):
            sbfs.append([])
            sfy=width-iy*dsfy
            if gmt_io: gmt_file.write(">>\n")
            titl = "subfault ix=%d,iy=%d, npts = %dx%d" % (ix,iy,nptx,npty)
            # GMT output
            (rlon,rlat) = prj_point(olon,olat,sfy*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,sfx/D2K,strike)
            dep = odep - sfy*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            sbfs[-1].append((lon,lat,dep))
            lt0 = lat
            ln0 = lon
            dp0 = dep
            (rlon,rlat) = prj_point(olon,olat,(sfy+dsfy)*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,sfx/D2K,strike)
            dep = odep - (sfy+dsfy)*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            sbfs[-1].append((lon,lat,dep))
            (rlon,rlat) = prj_point(olon,olat,(sfy+dsfy)*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,(sfx+dsfx)/D2K,strike)
            dep = odep - (sfy+dsfy)*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            sbfs[-1].append((lon,lat,dep))
            (rlon,rlat) = prj_point(olon,olat,sfy*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,(sfx+dsfx)/D2K,strike)
            dep = odep - sfy*mt.sin(dip*DEGR)
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (lon,lat,dep))
            sbfs[-1].append((lon,lat,dep))
            if gmt_io: gmt_file.write("%7.3f %7.3f %5.2f\n" % (ln0,lt0,dp0))
            #
            (rlon,rlat) = prj_point(olon,olat,(sfy+0.5*dsfy)*cosdp/D2K,strike-90.)
            (lon,lat)   = prj_point(rlon,rlat,(sfx+0.5*dsfx)/D2K,strike)
            dep = odep - (sfy+0.5*dsfy)*mt.sin(dip*DEGR)
            grid.append([lon,lat,dep,-(360.-strike),dip,rake,dsfx*dsfy*1.e6])
            isf += 1
            if ffi_io: ffi_file.write("%9.4f %9.4f %7.2f %6.1f %6.1f %6.1f %8.2f 0\n" % \
                                      (lon,lat,dep,-(360.-strike),dip,rake,dsfx*dsfy))
        
    if ffi_io: ffi_file.close()
    if gmt_io: gmt_file.close()
    return np.array(grid),np.array(sbfs)
    
    
def stacov(bx1,(elat,elon), stas, lats, lons, filename = None):
    #~ lines = open('fort.1','r')
    #~ stas = []
    #~ lats = []
    #~ lons = []
    #~ for line in lines:
        #~ if 'ref' in line:
            #~ stas.append(line.split()[0])
            #~ continue
        #~ if len(stas) > len(lats):
            #~ lat,lon = [float(x) for x in line.split()[0:2]]
            #~ lats.append(lat)
            #~ lons.append(lon)
    map = Basemap(projection='ortho',lat_0=elat,lon_0=elon,
              resolution='l',area_thresh=10000.)
    x,y = map(lons,lats)
    map.scatter(x,y,s=60,c='blue',marker='o',edgecolors='none',zorder=10)
    for i in range(0,len(stas)):
        sx=x[i]
        sy=y[i]
        #plt.text(sx,sy,stas[i],fontsize=12,color='black')
#        map.drawcoastlines()
        map.drawcountries()
        map.fillcontinents(color='coral')
        map.drawmapboundary()
        #map.drawmeridians(np.arange(0,360,30))
        #map.drawparallels(np.arange(-90,90,30))
    x1,y1 = map(elon,elat)
    plt.plot(x1,y1,'c*',markersize=40)   
    if filename:
        plt.savefig(filename, dpi=100  ,transparent=True, bbox_inches=0 )
