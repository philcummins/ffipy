#!/usr/bin/env python
import getopt
import operator
import sys
import math as mt
import numpy as np
import pylab as pl
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.colors import rgb2hex
from matplotlib.patches import Polygon
from matplotlib import mpl
from matplotlib import pyplot as plt
from obspy.core import read
import os
import csv
from obspy.imaging.beachball import Beach

#print os.environ['PYTHONPATH']
#sys.path.append( os.environ['PYTHONPATH'])
sys.path.append( os.environ['HOME'] + '/Seismology/src/ffi/python')
from ffi import fault_disp
from netCDF4 import Dataset

plrs = {'KIP':'d','MOO':'c','WAKE':'c','GUMO':'c','MIDW':'c','NWAO':'c','DAV':'c',
        'VNDA':'c','CASY':'c','MAJO':'c','TATO':'c','BTDF':'c','PET':'c','BJT':'c',
        'MBWA':'c','COR':'d','COLA':'c','ANMO':'d','PAYG':'d','PTCN':'d'}

def usage():
    print 'ffiplot.py'


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hd:o:vf:p", ["help", "dstmax=", "output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            return(0)
        else:
            assert False, "unhandled option"
#    if len(args) != 1:
#        usage()
#        return(-1)
    #####
    ffi_dir = '.'

    # Read epicentre from i_inv
    iinv_lines = open(ffi_dir+'/'+'i_inv','r').readlines()
    elon,elat,edep = [float(x) for x in iinv_lines[1].split()[0:3]]
#
    io = int(iinv_lines[15].split()[0])
    strike = float(iinv_lines[16].split()[3])
    if strike < 0.: strike += 360.
    dip    = float(iinv_lines[16].split()[4])
    rake   = float(iinv_lines[16].split()[4])
    print 'elon,elat = %g,%g; strike,dip = %g,%g' % (elon,elat,strike,dip)

    
    sbfs,top,nx,ddx,ny1,ny2,ddy = load_subfaults(ffi_dir+'/'+'Tonga.xyz')
    print ny1,ny2

    lnmn = sbfs[:,:,0].min()
    lnmx = sbfs[:,:,0].max()
    ltmn = sbfs[:,:,1].min()
    ltmx = sbfs[:,:,1].max()
    print lnmn,lnmx,ltmn,ltmx
    length = nx*ddx
    width = ny1*ddy
    bottom = top + width/mt.sin(dip*mt.pi/180.)
    print '%d subfaults read, ddy=%g, ddx=%g, ny=%d, nx=%d, length=%g, width=%g, top=%g, bottom=%g ' % \
        (len(sbfs),ddy,ddx,ny1,nx,length,width,top,bottom)
        #sys.exit(0)
    # Set up xy subfaults
    sfxys = []
    x = 0. 
    y = width-ddy
    for sbf in sbfs:
        sfxys.append([])
        sfxys[-1].append((x    ,y    ))
        sfxys[-1].append((x    ,y+ddy))
        sfxys[-1].append((x+ddx,y+ddy))
        sfxys[-1].append((x+ddx,y    ))
        sfxys[-1].append((x    ,y    ))
        y -= ddy
        if y < -0.1:
            y = width-ddy
            x += ddx
    # Load the rupture file and map to subfaults
    rpf_filename = 'rup-final'
    rpls = open(ffi_dir+'/'+rpf_filename,'r').readlines()
    slpmax = 0.
    for line in rpls:
        tmp  = line.split()
        slpmax = max(slpmax,0.01*float(tmp[6]))

    #slpmax *= 2
    slps = []
    slpt = []
    rlt = []
    rln = []
    stk = []
    dep = []
    rke=[]
    dpts = []
    dips = []
    cmap = pl.cm.hot_r # use reeversed 'hot' colormap
    norm = mpl.colors.Normalize(vmin=0., vmax=slpmax)
    for line in rpls:
        (lon,lat,dep,strike,dip,rake,slip,x,y) = [float(x) for x in line.split()] 
        rln.append(lon)
        rlt.append(lat)
        dpts.append(dep)
        stk.append(strike)
        dips.append(dip)
        rke.append(rake)
        slpt.append(slip)
        rake = rke[-1]
        slip = 0.01*slpt[-1]
        for i in range(0,len(sbfs)):
            if InsidePolygon(sbfs[i][:,0:2],Point(lon,lat)):
            # Set color
            # calling colormap with value between 0 and 1 returns
            # rgba value.  
                color = cmap(slip/slpmax)[:3]
                slps.append((i,color))
                break
        if i == len(sbfs)-1:
            print 'Failed to find subfault for lon,lat = %g,%g' % (lon,lat)
    print 'mapped %d rupture points, maximum slip is %g metres' % (len(slps),slpmax)


    #
    figprops = dict(figsize=(15., 7.5 / 1.618))
    #, dpi=128)                                          # Figure properties
    adjustprops = dict(left=0.05, bottom=0.05, right=0.9, top=0.93, wspace=0.4, hspace=0.4)      # Subplot properties
   
    fig = pl.figure(**figprops)                                                              # New figure
                                                                                             #fig.subplots_adjust(**adjustprops)                                                          # Tunes the subplot layout

    ax = fig.add_subplot(1, 3, 3)
    # Map Figure
    # Set map limits and plot the map
    print lnmn,lnmx,ltmn,ltmx
    xpnd = 1.
    lllon = 0.5*(lnmn+lnmx)-xpnd*(lnmx-lnmn)
    lllat = 0.5*(ltmn+ltmx)-xpnd*(ltmx-ltmn)
    urlon = 0.5*(lnmn+lnmx)+xpnd*(lnmx-lnmn)
    urlat = 0.5*(ltmn+ltmx)+xpnd*(ltmx-ltmn)
    m = Basemap(lllon,lllat,urlon,urlat,projection='merc',resolution='h')
    m.drawmapboundary(fill_color='aqua') 
    # fill continents, set lake color same as ocean color. 
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    if urlon - lllon < 5.:
        dlon = 1.0
    else:
        dlon =0.5
    m.drawparallels(np.arange(-90.,90.,0.5),labels=[1,0,0,0],labelstyle='+/-')
    m.drawmeridians(np.arange(-180.,180.,dlon),labels=[0,0,0,1],labelstyle='+/-')
    # Convert subfaults to map coordinates
    #ax = pl.gca()
    # First plot the polygon boundaries
    sfply = []
    for sbf in sbfs:
        sfply.append([])
        for i in range(0,len(sbf)):
            sfply[-1].append(m(sbf[i][0],sbf[i][1]))
        poly = Polygon(sfply[-1],facecolor='none',edgecolor='black',linewidth=0.1)
        ax.add_patch(poly)
    # Then plot the polygon faces
    for slp in slps:
        poly =  Polygon(sfply[slp[0]],facecolor=rgb2hex(slp[1]),edgecolor='blue')
        ax.add_patch(poly)

    (x,y) = m(elon,elat)
    plt.plot(x,y,'c*',markersize=20,zorder=41)

    (x,y) = m(-174.33, -20.02,)
    #plt.plot(x,y,'c*',markersize=20,zorder=41)

    (x,y,) = m(-174.00,-20.19)
    #plt.plot(x,y,'c*',markersize=20,zorder=41)

    #
    ax.annotate("IFP",
            xy=m(-174.041, -20.546), xycoords='data',fontsize=16,
            xytext=(0.45,0.1), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3")
            )
    ax.annotate("MFP",
            xy=m(-174.433,-19.530), xycoords='data',fontsize=16,
            xytext=(0.075,0.8), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3")
            )
    #
    if strike < 0.: strike += 360
    ax.text(-0.075, 0.95, 'c', transform=ax.transAxes,
            fontsize=20, fontweight='regular', va='top')
    # get axes bounds.
    pos = ax.get_position()
    l, b, w, h = pos.bounds

    # create axes instance for colorbar on right.
    cax = pl.axes([l+w+0.01, b+0.02, 0.02, 0.9*h])
    # draw colorbar on right.
    cbar = mpl.colorbar.ColorbarBase(ax=cax,cmap=cmap,norm=norm,orientation='vertical')
    cbar.set_label('metres slip')
    #

    #
    #bx1 = fig.add_subplot(1, 3, 2)
    #
    #plotzdisp(bx1,'okd_in')
    #stacov(bx1,elon,elat,)
    #ffi_fpslipplot(bx1,sfxys,length,width,bottom,top)
    #ex=50.
    #ey=27.2
    #fpslipxy(bx1,sbfs,slps,bottom,top,nx,ddx,ny1,ny2,ddy,ex=ex,ey=ey)

    #
    cx = fig.add_subplot(1, 3, 2)
    cx.text(-0.075, 0.95, 'b', transform=cx.transAxes,
            fontsize=24, fontweight='regular', va='top')
    azs = stacov(cx,elon,elat,edep)
    azs['COLA'] = 500.
    #plotsac(cx,['R', 'L'])
    #
    dx = fig.add_subplot(1, 3, 1)
    plotsac(dx,['P', 'S'],order=azs)
    dx.text(-0.075, 0.95, 'a', transform=dx.transAxes,
            fontsize=24, fontweight='regular', va='top')

    pl.savefig('ffi_'+os.path.abspath('.').split('/')[-1]+'.pdf')
    pl.show()

#
def plotzdisp(ax,okdin):
    lns = open(okdin,'r').readlines()
    elon=[]
    elat=[]
    edep=[]
    strk=[]
    dip=[]
    lng=[]
    wdt=[]
    disl1=[]
    disl2 = []
    for ln in lns[1:]:
        [la,lo,de,st,di,ln,wd,d1,d2,d3] = [float(x) for x in ln.split()]
        elat.append(la)
        elon.append(lo)
        edep.append(de)
        strk.append(st)
        dip.append(di)
        lng.append(ln)
        wdt.append(wd)
        disl1.append(d1)
        disl2.append(d2)
    elat=np.array(elat)
    elon=np.array(elon)
    edep=np.array(edep)
    strk=np.array(strk)
    dip=np.array(dip)
    lng=np.array(lng)
    wdt=np.array(wdt)
    disl1=np.array(disl1)
    disl2=np.array(disl2)
    rlon = np.arange(lllon,urlon+.05,.001)
    rlat = np.arange(lllat,urlat+.05,.001)
    x,y = np.meshgrid(rlon,rlat)
    u,v,w = fault_disp(elon, elat, edep, strk, dip, lng, wdt, disl1,\
                       disl2, x.flatten(),y.flatten())
    z=w.reshape((len(rlat),len(rlon)))
    #
    nc = Dataset('Flores_eq.nc','w',format='NETCDF3_CLASSIC')
    nc.title = 'Flores vertical displacement'
    nc.source = ''
    nc.createDimension('side',2)
    nc.createDimension('xysize',len(rlon)*len(rlat))
    y_range = nc.createVariable('y_range','d', ('side',))
    y_range.units = 'degrees'
    y_range[:] = [rlat[0],rlat[-1]]
    x_range = nc.createVariable('x_range','d', ('side',))
    x_range.units = 'degrees'
    x_range[:] = [rlon[0],rlon[-1]]
    z_range = nc.createVariable('z_range','d', ('side',))
    z_range.units = 'meters'
    z_range[:] = [w.min(),w.max()]

    spacing = nc.createVariable('spacing','d',('side',))
    spacing[:] = [rlon[1]-rlon[0],rlat[1]-rlat[0]]
    dimension = nc.createVariable('dimension','i',('side',))
    dimension[:] = [len(rlon),len(rlat)]
    grid_data = nc.createVariable('z','f', ('xysize',))
    grid_data.scale_factor = np.array([1.])
    grid_data.add_offset = np.array([0.])
    grid_data.node_offset = np.array([0])
    q = np.flipud(z)
    q = q.flatten()
    grid_data[:] = q.astype('float32')
    nc.close()
    #
    z=w.reshape((len(rlat),len(rlon)))
    m = Basemap(lllon,lllat,urlon,urlat,projection='merc',resolution='h')
    m.drawmapboundary(fill_color='aqua') 
# fill continents, set lake color same as ocean color. 
    xs,ys = m(x,y)
    plt.contourf(xs,ys,z)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('metres uplit/subsidence')
#        map.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color='coral')
    m.drawmapboundary()
    if urlon - lllon < 5.:
        dlon = 1.0
    else:
        dlon =0.5
    m.drawparallels(np.arange(-90.,90.,0.5),labels=[1,0,0,0],labelstyle='+/-')
    m.drawmeridians(np.arange(-180.,180.,dlon),labels=[0,0,0,1],labelstyle='+/-')    



# Determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This fuction
# returns True or False.  The algorithm is called
# "Ray Casting Method".

# determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs.

def ffi_fpslipplot(bx1,sfxys,length,width,top,bottom):
    bx1.set_xlim(0, length)
    bx1.set_ylim(0, width)
    bx1.set_ylabel('Along-dip Distance (km)')
    bx1.set_xlabel('Along-strike Distance (km)')
    # First plot the polygon boundaries
    for sbf in sfxys:
        poly = Polygon(sbf,facecolor='none',edgecolor='black')
        bx1.add_patch(poly)
    # Then plot the polygon faces
    for slp in slps:
        poly =  Polygon(sfxys[slp[0]],facecolor=rgb2hex(slp[1]),edgecolor='blue')
        bx1.add_patch(poly)
    bx2 = bx1.twinx()
    bx2.set_ylim(bottom,top)
    bx2.set_ylabel('Depth(km)')

def fpslipxy(bx1,sbfs,slps,top,bottom,nx,ddx,ny1,ny2,ddy,ex=None,ey=None):
    bx1.set_xlim(0, nx*ddx)
    bx1.set_ylim(0, (ny1+ny2)*ddy)
    bx1.set_ylabel('Up-dip distance (km)')
    bx1.set_xlabel('Along-strike Distance (km)')
    if ex != None:
        plt.plot(ex,ey,'c*',markersize=20,zorder=41)


    # First plot the polygon boundaries
    x = 0.
    y = (ny1-1)*ddy
    for i in range(0,nx*(ny1+ny2)):
        if i==nx*ny1:
            y = (ny1+ny2-1)*ddy
            x = 0.
        sfxy = [(x,y), (x,y+ddy), (x+ddx,y+ddy), (x+ddx,y), (x,y)]
        poly = Polygon(sfxy,facecolor='none',edgecolor='black')
        bx1.add_patch(poly)
        y -= ddy
        if i<nx*ny1:
            if y < -0.1:
                y = (ny1-1)*ddy
                x += ddx
        else:
            if y < ny1*ddy:
                y = (ny1+ny2-1)*ddy
                x += ddx
    # Then plot the polygon faces
    x = 0.
    y = (ny1-1)*ddy
    for i in range(0,nx*(ny1+ny2)):
        if i==nx*ny1:
            y = ny1*ddy
            x = (nx-1)*ddx
        slp = slps[i]
        sfxy = [(x,y), (x,y+ddy), (x+ddx,y+ddy), (x+ddx,y), (x,y)]
        poly = Polygon(sfxy,facecolor=rgb2hex(slp[1]),edgecolor='blue')
        bx1.add_patch(poly)
        if i<nx*ny1:
            y -= ddy
            if y < -0.1:
                y = (ny1-1)*ddy
                x += ddx
        else:
            y += ddy
            if y > (ny1+ny2-1)*ddy:
                y = ny1*ddy
                x -= ddx
    plt.plot([0.,nx*ddx],[ny1*ddy,ny1*ddy],linewidth=5,color='black')
    #    poly =  Polygon(sfxys[slp[0]],facecolor=rgb2hex(slp[1]),edgecolor='blue')
    #bx2 = bx1.twinx()
    #bx2.set_ylim(bottom,top)
    #bx2.set_ylabel('Depth(km)')

def stacov(bx1,elon,elat,edep):
    lines = open('fort.1','r').readlines()
    sites = {}
    for line in open('iris.site','r'):
        sta = line.split()[0]
        lat,lon = [float(x) for x in line.split()[3:5]]
        sites[sta] = (lat,lon)
    stas = []
    lats = []
    lons = []
    azs  = {}
    ps   = []
    clrs = []
    for i in range(0,len(lines)):
        line = lines[i]
        sta = line.split()[0]
        if sta in stas or not 'ref' in line:
            continue
        ib = int(lines[i+2].split()[1])
        if ib < 5:
            az,baz,gcarc,p = [float(x) for x in lines[i+1].split()[0:4]]
            azs[sta] = az
            ps.append(p)
        if not sites.has_key(sta):
            print 'Station %s not found in site table' % sta
            continue
        if not plrs.has_key(sta):
            clrs.append('gray')
        elif plrs[sta] == 'c':
            clrs.append('yellow')
        else:
            clrs.append('cyan')            
        stas.append(sta)
        lat,lon = sites[sta]
        #print sta,lat,lon,az,baz,gcarc
        lats.append(lat)
        lons.append(lon)
    map = Basemap(projection='ortho',lat_0=elat,lon_0=elon,
              resolution='l',area_thresh=10000.)
    x,y = map(lons,lats)
    map.scatter(x,y,s=25,c=clrs,marker='o',edgecolors='black',zorder=10)
    for i in range(0,len(stas)):
        sx=x[i]
        sy=y[i]
        if stas[i] == 'MOO':
            plt.text(sx,sy-800000,stas[i],fontsize=12,color='black')
        else:
            plt.text(sx,sy,stas[i],fontsize=12,color='black')
    map.fillcontinents(color='coral')
    map.drawmapboundary()
    map.drawmeridians(np.arange(0,360,30))
    map.drawparallels(np.arange(-90,90,30))
    x1,y1 = map(elon,elat)
    rad = 3000000
    fm = [194., 18., 90.]
    bx1.add_collection(Beach(fm, size=100, width=2.*rad, xy=(x1, y1), facecolor='red', linewidth=.6,zorder=45))
    fm = [359.,80,90.]
    bx1.add_collection(Beach(fm, size=100, width=2.*rad, xy=(x1, y1), facecolor='blue', linewidth=.6,alpha=0.5,zorder=47))
    r_src = 6371.-edep
    x2 = []
    y2 = []
    for i in range(0,len(stas)):
        xi = mt.asin(min(1.,7.95*ps[i]))
        xi = mt.sin(0.5*xi)
        x2.append(x1+rad*xi*mt.sin(azs[stas[i]]*mt.pi/180.))
        y2.append(y1+rad*xi*mt.cos(azs[stas[i]]*mt.pi/180.))
    map.scatter(x2,y2,s=25,c=clrs,marker='o',edgecolors='black',zorder=50)
    bx1.annotate("IFP",xy=(0.50,0.275), xycoords='axes fraction',fontsize=16,
            xytext=(0.55,0.15), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",zorder=1)
            )
    bx1.annotate("MFP",xy=(0.565,0.735), xycoords='axes fraction',fontsize=16,
            xytext=(0.70,0.75), textcoords='axes fraction',arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",zorder=2)
            )
    return azs

def plotsac(ax,phases,ncol=4,order=None):
    sgrms = {}
    i = 0
    for phase in phases:
        for file in os.listdir('.'):
            if file.endswith(phase+'.syn'): 
                sta = file.split('.')[0]
                if order != None and not order.has_key(sta):
                    print 'plotsac: sta %s not found in order{}' % sta
                    cobtinue
                staph = file.split('.')[0]+'.'+phase
                if order != None:
                    sgrms[staph] = order[sta] 
                else:
                    sgrms[staph] = i
                    i += 1

    # Order the seismograms
    staphas = sorted(sgrms.iteritems(), key=operator.itemgetter(1))

    tlen = 0.
    for staph,az in staphas:
        sac = read(staph+'.obs')[0]
        tlen = max(tlen,sac.stats.npts*sac.stats.delta)
        sac = read(staph+'.syn')[0]
        tlen = max(tlen,sac.stats.npts*sac.stats.delta)

    nwf = len(staphas)
    xmx = (ncol+1)*tlen
    ymx = nwf/ncol
    print 'ncol=%d, tlen = %g, nwf=%d, xmx=%g, ymx=%g' % (ncol,tlen,nwf,xmx,ymx)

    ax.set_xlim(0,xmx)
    ax.set_ylim(0.,ymx)
    #pl.title(' & '.join(phases)+' Waveforms')
    yoff = ymx-.5
    xoff = 0.25*tlen
    for staph,az in staphas:
            obs = read(staph+'.obs')[0]
            syn = read(staph+'.syn')[0]
            depmin = min(obs.data.min(),syn.data.min())
            depmax = max(obs.data.max(),syn.data.max())
            deprng = depmax-depmin
            if False:
                rngobs = obs.data.max()-obs.data.min()
                rngsyn = syn.data.max()-syn.data.min()
            else:
                rngobs = 2.*deprng
                rngsyn = 2.*deprng
                #depmid = 0.5*(depmin+depmax)
            time = xoff+obs.stats.delta*np.arange(0,obs.stats.npts)
            data = yoff+(obs.data-obs.data.mean())/deprng
            #        print sta,len(time),len(yoff+(obs.data-depmid)/deprng)
            #pl.plot(time,yoff+(obs.data-obs.data.mean())/rngobs,color='blue')
            pl.plot(time,yoff+obs.data/rngobs,color='blue')
            pl.hold(True)
            pl.plot(time,yoff+syn.data/rngobs,color='red')
            pl.text(xoff-0.1*tlen,yoff+0.25,staph,fontsize=9)
            yoff -= 1.
            if (yoff <= 0.):
                yoff = ymx-0.5
                xoff += 1.25*tlen
    ax.set_yticklabels([])

def point_in_poly(x,y,poly):

    n = len(poly)
    inside =False

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

## 
 
class Point:
    def __init__(self, h, v):
        self.h = h
        self.v = v
 
def InsidePolygon(polygon, p):
   angle = 0.0
   n = len(polygon)
 
   for i, (h, v) in enumerate(polygon):
      p1 = Point(h - p.h, v - p.v)
      h, v = polygon[(i + 1) % n]
      p2 = Point(h - p.h, v - p.v)
      angle += Angle2D(p1.h, p1.v, p2.h, p2.v);
 
   if abs(angle) < mt.pi:
      return False
 
   return True
 
def Angle2D(x1, y1, x2, y2):
   theta1 = mt.atan2(y1, x1)
   theta2 = mt.atan2(y2, x2)
   dtheta = theta2 - theta1
   while dtheta > mt.pi:
      dtheta -= 2.0 * mt.pi
   while dtheta < -mt.pi:
      dtheta += 2.0 * mt.pi
 
   return dtheta

def load_subfaults(filename):
    """ Load the subfault coordinates """
    sbfls = open(filename,'r').readlines()
    sbfs = []
    nx = 1
    ny1 = 1
    ny2 = 0
    dep_prev = 0.
    for line in sbfls:
        if line.startswith('>>'):
            sbfs.append([])
            continue
        (lon,lat,dep) = [float(x) for x in line.split()]
        if len(sbfs) == 1 and len(sbfs[0]) == 1:
            ddy = 111.195*mt.sqrt((lon-sbfs[0][0][0])**2 +
                                    (lat-sbfs[0][0][1])**2 +
                                       ((dep-sbfs[0][0][2])/111.195)**2)
            ddy = float(int(ddy+.1))
            dip1 = 180.*mt.asin(abs(dep-sbfs[-1][0][2])/ddy)/mt.pi
            dip1 = float(int(dip1+.5))
        elif len(sbfs) == 1 and len(sbfs[0]) == 2:
            ddx = 111.195*mt.sqrt((lon-sbfs[0][1][0])**2 +
                                    (lat-sbfs[0][1][1])**2 +
                                       ((dep-sbfs[0][1][2])/111.195)**2)
            ddx = float(int(ddx+.1))
        if len(sbfs[-1]) == 1:
            dip2 = 180.*mt.asin(abs(dep-sbfs[-1][0][2])/ddy)/mt.pi
            dip2 = float(int(dip2+.5))
        if len(sbfs[-1]) == 1:
            if dep < dep_prev:
                if dip2 == dip1:
                    nx += 1
                    ny1 = 1
                else:
                    ny2 = 1
            else:
                if dip2 == dip1:
                    ny1 += 1
                else:
                    ny2 += 1
            dep_prev = dep
        if len(sbfs) == 1 and len(sbfs[0]) == 1: top = dep
        sbfs[-1].append((lon,lat,dep))
    return np.array(sbfs),top,nx,ddx,ny1,ny2,ddy

###
if __name__ == "__main__":
    sys.exit(main())
