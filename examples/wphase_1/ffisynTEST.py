#!/home/roberto/bin/epy
import bpfilter as bp
import glob
import sys
import os
import numpy as np
from ffi import fault_grid
from obspy.core import read, Trace, Stream 
import DataUtils as  DU
from obspy.core.util.geodetics import locations2degrees, gps2DistAzimuth
import matplotlib.pyplot as plt
from obspy.taup.taup import getTravelTimes
from scipy.signal import triang, convolve
from mpl_toolkits.axes_grid1 import host_subplot, make_axes_locatable
from numpy.random import rand

def main(argv=sys.argv): 
    
    #Earth's parameters 
    #~ beta = 4.e3 #m/s 
    #~ rho = 3.e3 #kg/m^3 
    #~ mu = rho*beta*beta
    
    PLotSt = ["IU.TRQA.00.LHZ",
             "IU.LVC.00.LHZ",
             "II.NNA.00.LHZ",
              "IU.RAR.00.LHZ"]
             
             
    #PlotSubf = [143, 133, 123, 113, 103, 93,
     #           83, 73, 63, 53]
    PlotSubf = [6,3]

    
    
    #Set rup_vel = 0 to have a point source solution
    RupVel = 2.1 #Chilean eq from Lay et al
    t_h     = 10. # Half duration for each sf  
    noiselevel = 0.0# L1 norm level of noise
    mu =40e9
    #W-Phase filter 
    corners = 4.
    fmin = 0.001
    fmax = 0.005
    
    ### Data from Chilean 2010 EQ (Same as W phase inv.) 
    strike = 18.
    dip    = 18.
    rake   = 104. # 109.
    
    rakeA = rake + 45.
    rakeB = rake - 45.
    
    
    ### Fault's grid parameters
    nsx   = 21 #Number of sf along strike
    nsy   = 11 #Number of sf along dip
    flen  = 600. #Fault's longitude [km] along strike
    fwid  = 300. #Fault's longitude [km] along dip
    direc = 0    #Directivity 0 = bilateral
    Min_h = 10.  #Min depth of the fault
    
    
    ### Derivated parameters:
    nsf = nsx*nsy
    sflen = flen/float(nsx)         
    sfwid = fwid/float(nsy)
    swp = [1, 0, 2] # useful to swap (lat,lon, depth)  
    mindist = flen*fwid # minimun dist to the hypcen (initializing)
    
    ###Chessboard
    #weight = np.load("RealSol.npy") 
    weight = np.zeros(nsf)
    weight[::2] = 1 
    #weight[::2] = 1 
    #~ weight[10]=15
    #~ weight[5001]=10
    #~ weight[3201]=2
    
    
    
    ## Setting dirs and reading files.
    GFdir = "/home/roberto/data/GFS/"
    workdir = os.path.abspath(".")+"/"
    datadir = workdir + "DATA/"
    tracesfilename = workdir + "goodtraces.dat"
    tracesdir = workdir + "WPtraces/"
    
    try:
        reqfilename    = glob.glob(workdir + '*.syn.req')[0]
    except IndexError:   
        print "There is not *.syn.req file in the dir"
        sys.exit()
    
    basename = reqfilename.split("/")[-1][:-4]
    
    if not os.path.exists(tracesfilename): 
        print tracesfilename, "does not exist."
        exit()
    
    if not os.path.exists(datadir):
            os.makedirs(datadir)
    
    if not os.path.exists(tracesdir):
            os.makedirs(tracesdir)
 
    tracesfile = open(tracesfilename)    
    reqfile =  open(reqfilename)    
    
    trlist = readtraces(tracesfile)
    eqdata = readreq(reqfile)    

    tracesfile.close()
    reqfile.close()   
    
    ####Hypocentre from
    ### http://earthquake.usgs.gov/earthquakes/eqinthenews/2010/us2010tfan/    
    cmteplat = -35.91#-35.85#-36.03#-35.83
    cmteplon = -72.73#-72.72#-72.83# -72.67
    cmtepdepth= 35.
    eq_hyp = (cmteplat,cmteplon,cmtepdepth)
    
    
      ############
    

    # Defining the sf system
    grid, sblt = fault_grid('CL-2010',cmteplat,cmteplon,
                            cmtepdepth, direc,
                            Min_h, strike, dip, rake, flen,fwid ,nsx,nsy,
                            Verbose=False,ffi_io=True,gmt_io=True)
    
    print ('CL-2010',cmteplat,cmteplon,
                            cmtepdepth, direc,
                            Min_h, strike, dip, rake, flen,fwid ,nsx,nsy)
    print grid[0][1]
    #sys.exit()
    #This calculation is inside of the loop
    #~ NP = [strike, dip, rake]
    #~ M = np.array(NodalPlanetoMT(NP))  
    #~ Mp = np.sum(M**2)/np.sqrt(2)    
     
    #############################################################################
    ######Determining the sf closest to the hypocentre:    
    min_Dist_hyp_subf = flen *fwid
    for subf in range(nsf):
        sblat   = grid[subf][1]
        sblon   = grid[subf][0]
        sbdepth = grid[subf][2]              
        sf_hyp =  (sblat,sblon, sbdepth)        
        Dist_hyp_subf = hypo2dist(eq_hyp,sf_hyp)
        if Dist_hyp_subf < min_Dist_hyp_subf:
            min_Dist_hyp_subf = Dist_hyp_subf
            min_sb_hyp = sf_hyp
            hyp_subf = subf
    ####Determining trimming times:    
    test_tr = read(GFdir + "H003.5/PP/GF.0001.SY.LHZ.SAC")[0]
    t0 = test_tr.stats.starttime
    TrimmingTimes = {}   # Min. Distace from the fault to each station. 
    A =0
    for trid in trlist:     
        metafile = workdir + "DATA/" + "META." + trid + ".xml"
        META = DU.getMetadataFromXML(metafile)[trid]
        stlat = META['latitude']
        stlon = META['longitude'] 
        dist =   locations2degrees(min_sb_hyp[0],min_sb_hyp[1],\
                                   stlat,stlon) 
        parrivaltime = getTravelTimes(dist,min_sb_hyp[2])[0]['time']        
        ta = t0 + parrivaltime
        tb = ta + round(15.*dist) 
        TrimmingTimes[trid] = (ta, tb)
        
    
    ###########################

      
    
    DIST = []
    # Ordering the stations in terms of distance
    for trid in trlist: 
        metafile = workdir + "DATA/" + "META." + trid + ".xml"
        META = DU.getMetadataFromXML(metafile)[trid]
        lat = META['latitude']
        lon = META['longitude']
        trdist = locations2degrees(cmteplat,
                                   cmteplon,lat,lon) 
        DIST.append(trdist)   

    DistIndex = lstargsort(DIST)
    trlist = [trlist[i] for i in DistIndex]
  
    stdistribution = StDistandAzi(trlist, eq_hyp , workdir + "DATA/")
    StDistributionPlot(stdistribution)
    #exit()
    #Main loop
   

 

        
    for subf in range(nsf):
        print subf
        sflat   = grid[subf][1]
        sflon   = grid[subf][0]           
        sfdepth = grid[subf][2]
        #~ strike = grid[subf][3] #+ 360.
        #~ dip    = grid[subf][4]
        #~ rake   = grid[subf][5] #     
        NP = [strike, dip, rake]  
        NPA = [strike, dip, rakeA]
        NPB = [strike, dip, rakeB]        


        
        M = np.array(NodalPlanetoMT(NP))   
        MA = np.array(NodalPlanetoMT(NPA)) 
        MB = np.array(NodalPlanetoMT(NPB)) 
        #Time delay is calculated as the time in which 
        #the rupture reach the subfault
            
        sf_hyp = (sflat, sflon, sfdepth) 
        Dist_ep_subf = hypo2dist(eq_hyp,sf_hyp)
        
        if Dist_ep_subf < mindist:
            mindist = Dist_ep_subf
            minsubf = subf
        
                
        if RupVel == 0:
            t_d = eqdata['time_shift']
        else:
            t_d = round(Dist_ep_subf/RupVel) #-59.
       
        print sflat, sflon, sfdepth
        # Looking for the best depth dir:
        depth = []
        depthdir = []
        for file in os.listdir(GFdir):
            if file[-2:] == ".5":
                depthdir.append(file)
                depth.append(float(file[1:-2]))            
        BestDirIndex = np.argsort(abs(sfdepth\
                                  - np.array(depth)))[0]      
        hdir = GFdir + depthdir[BestDirIndex] + "/"     
        
        ###

        SYN = np.array([])
        SYNA = np.array([])
        SYNB = np.array([])
        for trid in trlist:     
            
            metafile = workdir + "DATA/" + "META." + trid + ".xml"
            META = DU.getMetadataFromXML(metafile)[trid]
            lat = META['latitude']
            lon = META['longitude']  
            
            #Subfault loop               
            #GFs Selection:
            ##Change to folloing loop
            
            dist = locations2degrees(sflat,sflon,lat,lon)                                
            azi =  -np.pi/180.*gps2DistAzimuth(lat,lon,
                       sflat,sflon)[2] 
            trPPsy,  trRRsy, trRTsy,  trTTsy = \
                                       GFSelectZ(hdir,dist)          
            
            
 
            
            trROT =  MTrotationZ(azi, trPPsy,  trRRsy, trRTsy,  trTTsy) 
            orig = trROT[0].stats.starttime  
            dt = trROT[0].stats.delta                       

            trianglen = 2.*int(t_h/dt)-1.
            FirstValid = int(trianglen/2.) + 1 # to delete
            window = triang(trianglen)
            window /= np.sum(window)
            #window = np.array([1.])
            
      
            
            
            parrivaltime = getTravelTimes(dist,sfdepth)[0]['time']
            
            t1 = TrimmingTimes[trid][0] - t_d
            t2 = TrimmingTimes[trid][1] - t_d
            
            
            
            for trR in trROT:
                trR.data *= 10.**-21 ## To get M in Nm                   
                trR.data -= trR.data[0]
                AUX1 = len(trR)
                trR.data = convolve(trR.data,window,mode='valid') 
                AUX2 = len(trR)
                mean = np.mean(np.hstack((trR.data[0]*np.ones(FirstValid),\
                               trR.data[:60./trR.stats.delta*1.-FirstValid+1])))
                #mean = np.mean(trR.data[:60])
                trR.data -= mean      
                trR.data = bp.bandpassfilter(trR.data,len(trR), trR.stats.delta,\
                                             corners , 1 , fmin, fmax)  
                t_l = dt*0.5*(AUX1 - AUX2)                             
                trR.trim(t1-t_l,t2-t_l, pad=True, fill_value=trR.data[0])  #We lost t_h due to the convolution        
            


                   
            #~ for trR in trROT:
                #~ trR.data *= 10.**-23 ## To get M in Nm                   
                #~ trR.data -= trR.data[0]
 
                #~ trR.data = convolve(trR.data,window,mode='same') 

                #~ #mean = np.mean(np.hstack((trR.data[0]*np.ones(FirstValid),\
                               #~ #trR.data[:60./trR.stats.delta*1.-FirstValid+1])))
                #~ mean = np.mean(trR.data[:60])
                #~ trR.data -= mean      
                #~ trR.data = bp.bandpassfilter(trR.data,len(trR), trR.stats.delta,\
                                             #~ corners , 1 , fmin, fmax)  
                            
                #~ trR.trim(t1,t2,pad=True, fill_value=trR.data[0])     
           
            trROT = np.array(trROT)  
            syn  =  np.dot(trROT.T,M) 
            synA =  np.dot(trROT.T,MA)
            synB =  np.dot(trROT.T,MB)
            
            SYN = np.append(SYN,syn)  
            SYNA = np.append(SYNA,synA)
            SYNB = np.append(SYNB,synB)
            
            
        print np.shape(A), np.shape(np.array([SYN]))    
        if subf == 0: 
            A = np.array([SYN])
            AA = np.array([SYNA])
            AB = np.array([SYNB])
        else:
            A = np.append(A,np.array([SYN]),0)    
            AA = np.append(AA,np.array([SYNA]),0)
            AB = np.append(AB,np.array([SYNB]),0)
            
            
            
    AC = np.vstack((AA,AB))
    print np.shape(AC)
    print np.shape(weight)
    B = np.dot(A.T,weight)
    stsyn = Stream()
    n = 0
    Ntraces= {}
    for trid in trlist: 
        spid = trid.split(".")        
        print trid
        NMIN = 1. + (TrimmingTimes[trid][1] - TrimmingTimes[trid][0]) / dt
        Ntraces[trid] = (n,NMIN + n)
        trsyn = Trace(B[n:NMIN+n])   
        n += NMIN        
        trsyn.stats.network = spid[0]
        trsyn.stats.station = spid[1]
        trsyn.stats.location = spid[2]
        trsyn.stats.channel = spid[3] 
        trsyn = AddNoise(trsyn,level = noiselevel)
        #trsyn.stats.starttime = 
        stsyn.append(trsyn)
        
       
    stsyn.write(workdir+"WPtraces/" + basename + ".decov.trim.mseed",
                 format="MSEED")           
                
    #####################################################    
    # Plotting:
    #####################################################
    #we are going to reflect the y axis later, so:
    print minsubf
    hypsbloc = [minsubf / nsy , -(minsubf % nsy) - 2]

    #Creating the strike and dip axis:
    StrikeAx= np.linspace(0,flen,nsx+1)
    DipAx= np.linspace(0,fwid,nsy+1)
    DepthAx = DipAx*np.sin(np.pi/180.*dip) + Min_h    
    hlstrike = StrikeAx[hypsbloc[0]] + sflen*0.5
        
    hldip = DipAx[hypsbloc[1]] + sfwid*0.5 
    hldepth = DepthAx[hypsbloc[1]] + sfwid*0.5*np.sin(np.pi/180.*dip)
       
    StrikeAx = StrikeAx - hlstrike
    DipAx =     DipAx   - hldip
 

    
    XX, YY = np.meshgrid(StrikeAx, DepthAx)
    XX, ZZ = np.meshgrid(StrikeAx, DipAx )

   
    sbarea = sflen*sfwid
    
    SLIPS = weight.reshape(nsx,nsy).T#[::-1,:]
    SLIPS /= mu*1.e6*sbarea
    
    ######Plot:#####################
    plt.figure()
    ax = host_subplot(111)
    im = ax.pcolor(XX, YY, SLIPS, cmap="jet")    
    ax.set_ylabel('Depth [km]')       
    ax.set_ylim(DepthAx[-1],DepthAx[0])  
    
    # Creating a twin plot 
    ax2 = ax.twinx()
    #im2 = ax2.pcolor(XX, ZZ, SLIPS[::-1,:], cmap="Greys") 
    im2 = ax2.pcolor(XX, ZZ, SLIPS[::-1,:], cmap="jet")    
    
    ax2.set_ylabel('Distance along the dip [km]')
    ax2.set_xlabel('Distance along the strike [km]')    
    ax2.set_ylim(DipAx[0],DipAx[-1])
    ax2.set_xlim(StrikeAx[0],StrikeAx[-1])       
                         
                         
    ax.axis["bottom"].major_ticklabels.set_visible(False) 
    ax2.axis["bottom"].major_ticklabels.set_visible(False)
    ax2.axis["top"].set_visible(True)
    ax2.axis["top"].label.set_visible(True)
    
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Slip [m]") 
    ax2.plot([0], [0], '*', ms=225./(nsy+4))
    ax2.set_xticks(ax2.get_xticks()[1:-1])
    #ax.set_yticks(ax.get_yticks()[1:])
    #ax2.set_yticks(ax2.get_yticks()[:-1])
    

    
    #########Plotting the selected traces:
    nsp = len(PLotSt) * len(PlotSubf)
    plt.figure(figsize=(13,11))
    plt.title("Synthetics for rake = " + str(round(rake)))
    mindis = []
    maxdis = []
    for i, trid in enumerate(PLotSt):   
        x = np.arange(0,Ntraces[trid][1]-Ntraces[trid][0],
                      dt)
        for j, subf in enumerate(PlotSubf):
            y = A[subf, Ntraces[trid][0]:Ntraces[trid][1]]
            if j == 0:
                yy = y
            else:
                yy = np.vstack((yy,y))        
        maxdis.append(np.max(yy))
        mindis.append(np.min(yy))
        
    

    for i, trid in enumerate(PLotSt):   
        x = np.arange(0,Ntraces[trid][1]-Ntraces[trid][0],
                      dt)

        for j, subf in enumerate(PlotSubf):
            y = A[subf, Ntraces[trid][0]:Ntraces[trid][1]]
            plt.subplot2grid((len(PlotSubf), len(PLotSt)),
                              (j, i))                                
            plt.plot(x,y, linewidth=2.5)
            if j == 0:
                plt.title(trid)
            fig = plt.gca()            
            fig.axes.get_yaxis().set_ticks([])
            fig.set_ylabel(str(subf),rotation=0)
            fig.set_xlim((x[0],x[-1]))
            fig.set_ylim((mindis[i],maxdis[i]))
            if subf != PlotSubf[-1]:
                fig.axes.get_xaxis().set_ticks([])

    
    plt.show()
    
    
    
    
def readtraces(tracesfile):
    trlist = []
    while 1:
        line = tracesfile.readline().rstrip('\r\n')
        if not line: break       
        trlist.append(line.split()[0])           
    return trlist
    
def readreq(reqfile):    
    parameters={}    
    while 1: 
        line = reqfile.readline().rstrip('\r\n')
        if not line: break 
        line = line.split()
        key = line[0]
        val = line[1:]
        parameters[key] = val 
    
    parameters['eplat'] = float(parameters['eplat'][0])
    parameters['eplon'] = float(parameters['eplon'][0])
    parameters['epdepth'] = float(parameters['epdepth'][0])
    parameters['origin_time'] = parameters['origin_time'][0]
    parameters['before_p'] = float(parameters['before_p'][0])
    parameters['after_p'] = float(parameters['after_p'][0])
    parameters['half_duration'] = float(parameters['half_duration'][0])
    parameters['time_shift'] = float(parameters['time_shift'][0])
    parameters['maxdist'] = float(parameters['maxdist'][0])
    
    mt = ('mrr', 'mtt', 'mpp', 'mtp', 'mrt', 'mrp' )
    
    if all (k in parameters for k in mt):
        parameters['mt'] = [ float(parameters['mrr'][0]),\
                             float(parameters['mtt'][0]),\
                             float(parameters['mpp'][0]),\
                             float(parameters['mrt'][0]),\
                             float(parameters['mrp'][0]),\
                             float(parameters['mtp'][0]),\
                           ]
   
    return parameters    

def NodalPlanetoMT(NP):
    #M = [rr, tt, pp, rt, rp, tp]
    #NP = [Strike, dip, rake]

    ToRad = np.pi/180.
    
    phi   = ToRad *NP[0]
    delt  = ToRad *NP[1]
    lam   = ToRad *NP[2]

    Mtt =  -(np.sin(delt)*np.cos(lam)*np.sin(2*(phi))\
           + np.sin(2*delt)*np.sin(lam)*np.sin(phi)*np.sin(phi))# = Mxx   
    Mpp =    np.sin(delt)*np.cos(lam)*np.sin(2*(phi))\
           - np.sin(2*delt)*np.sin(lam)*np.cos(phi)*np.cos(phi)#  = Myy
    Mrr =    np.sin(2*delt)*np.sin(lam)                        #  = Mzz  
    Mtp =    -(np.sin(delt)*np.cos(lam)*np.cos(2*(phi)) +\
               0.5*np.sin(2*delt)*np.sin(lam)*np.sin(2*phi))  #    =-Mxy
    Mrp =   np.cos(delt)*np.cos(lam)*np.sin(phi)\
            - np.cos(2*delt)*np.sin(lam)*np.cos(phi)            # =-Myz 
    Mrt =   -(np.cos(delt)*np.cos(lam)*np.cos(phi) +\
             np.cos(2*delt)*np.sin(lam)*np.sin(phi))            # = Mxz
    
    M =[Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
    return M
def lstargsort(seq):
   #http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
   #by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)      
    

def GFSelectZ(hdir,dist): 
    dist_str = str(int(dist*10.)/2*2+1)#Some GFs have only odd dists.
    #dist_str = str(int(dist*10.))#Some GFs have only odd dists.
    dist_form = dist_str.zfill(4)    
     
    ## Loading files
    trPP = read(hdir + "PP/GF." + dist_form + ".SY.LHZ.SAC"  )[0]    
    #trPP.data -= trPP.data[0]


    trRR = read(hdir + "RR/GF." + dist_form + ".SY.LHZ.SAC"  )[0] 
    #trRR.data -= trRR.data[0]


    trRT = read(hdir + "RT/GF." + dist_form + ".SY.LHZ.SAC"  )[0] 
    #trRT.data -= trRT.data[0]



    trTT = read(hdir + "TT/GF." + dist_form + ".SY.LHZ.SAC"  )[0] 
    #trTT.data -= trTT.data[0]   
    return trPP, trRR, trRT,  trTT  
    
    
def MTrotationZ(azi, trPPsy,  trRRsy, trRTsy,  trTTsy):
    #Creating traces for the rotated syn:
        
        trRRsy_rot =   trRRsy.copy()
        trPPsy_rot =   trPPsy.copy()
        trTTsy_rot =   trTTsy.copy()
        trTPsy_rot =   trTTsy.copy()
        trRTsy_rot =   trRTsy.copy()
        trRPsy_rot =   trPPsy.copy()
                               
        #Obtaining the rotated synthetics:

        
        sinp = np.sin(azi)
        cosp = np.cos(azi)
        
        trRRsy_rot.data =   trRRsy.data           
        trPPsy_rot.data =   sinp*sinp*trTTsy.data+cosp*cosp*trPPsy.data
        trTTsy_rot.data =   cosp*cosp*trTTsy.data+sinp*sinp*trPPsy.data
        trTPsy_rot.data =   2.*sinp*cosp*(trTTsy.data-trPPsy.data)
        trRTsy_rot.data =   cosp*trRTsy.data
        trRPsy_rot.data =   sinp*trRTsy.data  
        return [trRRsy_rot, trTTsy_rot,  trPPsy_rot, trRTsy_rot,trRPsy_rot, trTPsy_rot]    
        #return [trRRsy_rot, trPPsy_rot,  trTTsy_rot, trTPsy_rot,trRTsy_rot, trRPsy_rot] 

def hypo2dist(hyploc1,hyploc2,r=6371.):
    ''' Given  2 tuples with hypcenter coordinates
        (lat, lon, depth) returns the distance in KM
        between them. '''
    r1 = r-hyploc1[2]
    r2 = r-hyploc2[2]    
    theta = np.pi/180.*locations2degrees(hyploc1[0],hyploc1[1],\
                                         hyploc2[0],hyploc2[1])   
    ##Law of Cosines:    
    l = np.sqrt(r1**2.+r2**2.-2.*r2*r1*np.cos(theta))
    return l


def AddNoise(trace, level=0.1):
    '''Given a obspytrace return a noisy trace.
       level is the noise ratio (In the sense of the L1 norm)'''
    
        
      
    noise = np.mean(np.abs(trace.data)) * level\
                    * (rand(len(trace.data))*4.-2)
    noisy = trace.copy() 
    noisy.data = trace.data + noise
    
    #noise = 4*level*(rand(trace.stats.npts)-0.5)*trace.data    
    #noisy = trace.copy()
    #noisy.data = noise + trace.data
   
    
    return noisy
    
    
def StDistandAzi(traces, hyp, dir):
    '''
    Given a list with st ids, a tuple with 
    (lat, lon, depth) of the hypocentre and 
    thw directory with the xml metafiles; it returns 
    dictionary with the distance and the azimuth of
    each station       
       
    '''
    Stdistribution = {}
    for trid in traces:     
            
        metafile = dir + "META." + trid + ".xml"
        META = DU.getMetadataFromXML(metafile)[trid]
        lat = META['latitude']
        lon = META['longitude']            
        dist = locations2degrees(hyp[0],hyp[1],lat,lon)                                
        azi =  -np.pi/180.*gps2DistAzimuth(lat,lon,
                   hyp[0],hyp[1])[2]
        Stdistribution[trid] = {'azi' : azi, 'dist' : dist} 

    
    #~ fig = plt.figure()  
    #~ plt.rc('grid', color='#316931', linewidth=1, linestyle='-')
    #~ plt.show()
    
    return Stdistribution
        
def StDistributionPlot(Stdistribution):
    ''' 
        Given a dict cointaning azi and dist for each stations
        it makes a plot of the distribution
    ''' 
        
    fig = plt.figure(figsize=(10,10))
    plt.rc('grid', color='#316931', linewidth=1, linestyle='-')  
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='#d5de9c')
   
    

    for trid, val in Stdistribution.iteritems():
        idsplit = trid.split(".")
        stname = idsplit[1]       
        ax.plot([val['azi']],[val['dist']],'o', ms=15)
        ax.annotate(stname,xytext=(val['azi'], val['dist']),
                    xy = (val['azi'], val['dist']))
    ax.set_rgrids(ax.get_yticks(),  angle=135)
    ticks = [str(tick) + r'$^o$' for tick in ax.get_yticks()]
    print ticks
    ax.set_yticklabels(ticks)          
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(1)
    ax.set_thetagrids([0, 90, 180, 270], labels = 
                       ['N', 'E', 'S', 'W' ]  )
       




    
if __name__ == "__main__":
    main() 
