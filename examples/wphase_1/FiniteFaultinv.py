#!/usr/bin/env python
from runave import runave
import bpfilter as bp
import sys
import os
from obspy.core import read, Trace, Stream, UTCDateTime
from obspy.taup.taup import getTravelTimes
from obspy.core.util.geodetics import locations2degrees, gps2DistAzimuth
import numpy as np
import DataUtils as DU
from obspy.iris import Client
import matplotlib.pyplot as plt 
from scipy.signal import triang, convolve, detrend
from scipy.linalg import lstsq, pinv2
from obspy.imaging.beachball import Beachball 
from obspy.signal.array_analysis import cosine_taper
from ffi import fault_grid
from scipy.optimize import nnls
from mpl_toolkits.axes_grid1 import host_subplot, make_axes_locatable
from matplotlib.ticker import FormatStrFormatter





def main(argv=sys.argv): 

    global eplat, eplon, epdepth, orig
    
    
    GFdir = "/home/roberto/data/GFS/"
    beta = 4.e3 #m/s 
    rho = 3.e3 #kg/m^3 
    mu = rho*beta*beta
    mu =40e9
    
    Lbdm0min = 1e-26*np.array([125.])
    Lbdsmooth = 1e-26*np.array([100.])
    
    #~ Lbdm0min = 1e-26*np.linspace(60.,500,40)
    #~ Lbdsmooth = 1e-26*np.linspace(60.,500,40)#*0.5

    corners = 4.
    fmin = 0.001
    fmax = 0.005
    
    ### Data from Chilean 2010 EQ (Same as W phase inv.) 
    strike = 18.
    dip    = 18.
    rake   = 104. # 109.
    #rake = 45.
    
    rakeA = rake + 45.
    rakeB = rake - 45.
    ####################
    nsx = 21
    nsy = 11
    Min_h = 10.
    flen  = 600. #Fault's longitude [km] along strike
    fwid  = 300. #Fault's longitude [km] along dip
    sflen = flen/float(nsx)
    sfwid = fwid/float(nsy)    
    swp = [1, 0, 2]
    nsf = nsx*nsy
    ###################
    
    t_h     = 10.
    MISFIT = np.array([])
    #RUPVEL = np.arange(1.0, 5.0, 0.05)
    RupVel = 2.1 # Best fit
    #RupVel = 2.25 #From Lay et al.
    
    
    #for RupVel in RUPVEL:
    print "****************************"
    print RupVel
    print "****************************"
    NP = [strike, dip, rake]
    NPA = [strike, dip, rakeA]
    NPB = [strike, dip, rakeB]
    
    M  = np.array(NodalPlanetoMT(NP))  
    MA = np.array(NodalPlanetoMT(NPA)) 
    MB = np.array(NodalPlanetoMT(NPB)) 
    
    Mp = np.sum(M**2)/np.sqrt(2)
    
    
    #############
        #Loading req file and EQparameters
    parameters={}    
    with open(argv[1],'r') as file:
        for line in file:
            line = line.split()
            key = line[0]
            val = line[1:]
            parameters[key] = val 
    #~ cmteplat = float(parameters['eplat'][0])
    #~ cmteplon = float(parameters['eplon'][0])
    #~ cmtepdepth=float(parameters['epdepth'][0])    
    orig = UTCDateTime(parameters['origin_time'][0])
    
    ####Hypocentre from
    ### http://earthquake.usgs.gov/earthquakes/eqinthenews/2010/us2010tfan/    
    cmteplat = -35.91#-35.85#-36.03#-35.83
    cmteplon = -72.73#-72.72#-72.83# -72.67
    cmtepdepth= 35.
    eq_hyp = (cmteplat,cmteplon,cmtepdepth)
      ############

    
    
    grid, sblt = fault_grid('CL-2010',cmteplat,cmteplon,cmtepdepth,0, Min_h,\
                             strike, dip, rake, flen, fwid, nsx, nsy,                            Verbose=False, ffi_io=True, gmt_io=True)
                             
    print ('CL-2010',cmteplat,cmteplon,cmtepdepth,0, Min_h,\
                             strike, dip, rake, flen, fwid, nsx, nsy,\
                            )
    print grid[0][1]
    #sys.exit()
    #############
    #Loading files and setting dirs:
    
    inputfile =  os.path.abspath(argv[1]) 
    if not os.path.exists(inputfile): print inputfile, "does not exist."; exit() 
    workdir = "/".join(inputfile.split("/")[:-1]) 
    basename = inputfile.split("/")[-1][:-4]
    if workdir[-1] != "/": workdir += "/"
    
    try :
        os.mkdir(workdir+"WPinv")
    except OSError:
        pass#print "Directory WPtraces already exists. Skipping"
    
    trfile = open(workdir+"goodtraces.dat")
    trlist = []
    #Loading Good traces files:
    while 1:
        line = trfile.readline().rstrip('\r\n')
        if not line: break        
        trlist.append(line.split()[0])        
    trfile.close()
    #############
    
    # Reading traces:    
    st = read(workdir+"WPtraces/" + basename + ".decov.trim.mseed")  
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
    print hyp_subf,  min_sb_hyp,  min_Dist_hyp_subf    
    
    
    ####Determining trimming times:
    
    test_tr = read(GFdir + "H003.5/PP/GF.0001.SY.LHZ.SAC")[0]
    t0 = test_tr.stats.starttime
    TrimmingTimes = {}   # Min. Distace from the fault to each station. 
    A =0
    for trid in trlist:     
        tr = st.select(id=trid)[0]
        metafile = workdir + "DATA/" + "META." + tr.id + ".xml"
        META = DU.getMetadataFromXML(metafile)[tr.id]
        stlat = META['latitude']
        stlon = META['longitude'] 
        dist =   locations2degrees(min_sb_hyp[0],min_sb_hyp[1],\
                                   stlat,stlon) 
        parrivaltime = getTravelTimes(dist,min_sb_hyp[2])[0]['time']        
        ta = t0 + parrivaltime
        tb = ta + round(15.*dist) 
        TrimmingTimes[trid] = (ta, tb)
        
    
    ##############################################################################
     
    #####

    DIST = []
    # Ordering the stations in terms of distance
    for trid in trlist: 
        metafile = workdir + "DATA/" + "META." + trid + ".xml"
        META = DU.getMetadataFromXML(metafile)[trid]
        lat = META['latitude']
        lon = META['longitude']
        trdist = locations2degrees(cmteplat,cmteplon,lat,lon) 
        DIST.append(trdist)   

    
    DistIndex = lstargsort(DIST)
    
    if len(argv) == 3:
        trlist = [argv[2]]
        OneStation = True
    else:
        trlist = [trlist[i] for i in DistIndex]
        OneStation = False
        
   ##### 

    client = Client()
    ObservedDisp = np.array([])   
    gridlat = []
    gridlon = []
    griddepth = []
    sbarea = []
    mindist = flen*fwid # min distance hyp-subfault 
    

    ##########Loop for each subfault
    for subf in range(nsf):
        print "**********"
        print subf
        eplat   = grid[subf][1]
        eplon   = grid[subf][0]           
        epdepth = grid[subf][2]
        
        ## Storing the subfault's location centered in the hypcenter 
        gridlat.append(eplat-cmteplat)
        gridlon.append(eplon-cmteplon)
        griddepth.append(epdepth)
        
        
        strike = grid[subf][3] #+ 360.
        dip    = grid[subf][4]
        rake   = grid[subf][5] #     
        NP = [strike, dip, rake]
        
        M = np.array(NodalPlanetoMT(NP))   

        
                #Calculating the time dalay:
            
        sf_hyp = (eplat,eplon, epdepth)        
        Dist_ep_subf = hypo2dist(eq_hyp,sf_hyp)
        t_d = round(Dist_ep_subf/RupVel) #-59.
        print eplat,eplon, epdepth
    
        #t_d  = 0.
        
    
        # Determining depth dir:
        depth = []
        depthdir = []
        for file in os.listdir(GFdir):
            if file[-2:] == ".5":
                depthdir.append(file)
                depth.append(float(file[1:-2]))            
        BestDirIndex = np.argsort(abs(epdepth-np.array(depth)))[0]      
        hdir = GFdir + depthdir[BestDirIndex] + "/"   
        # hdir is the absolute path to the closest deepth. 
        
        
        
        SYN = np.array([])
        SYNA = np.array([])
        SYNB = np.array([])

        
        #Main loop :
        for trid in trlist:  
                       
            tr = st.select(id=trid)[0]    
            metafile = workdir + "DATA/" + "META." + tr.id + ".xml"
            META = DU.getMetadataFromXML(metafile)[tr.id]
            lat = META['latitude']
            lon = META['longitude']    
            trPPsy,  trRRsy, trRTsy,  trTTsy = \
                                   GFSelectZ(lat,lon,hdir) 
            
            tr.stats.delta = trPPsy.stats.delta
            azi =   -np.pi/180.*gps2DistAzimuth(lat,lon,\
                                               eplat,eplon)[2]
            trROT = MTrotationZ(azi, trPPsy,  trRRsy, trRTsy,  trTTsy)        
                        
            
                    #Triangle 
            dt = trROT[0].stats.delta          
            trianglen = 2.*t_h/dt-1.
            window = triang(trianglen)
            window /= np.sum(window)
            #window = np.array([1.])
            
            FirstValid = int(trianglen/2.) + 1
            dist =   locations2degrees(eplat,eplon,lat,lon) 
            parrivaltime = getTravelTimes(dist,epdepth)[0]['time']
            
            t1 = TrimmingTimes[trid][0] - t_d
            t2 = TrimmingTimes[trid][1] - t_d
            
   
            #~ t1 = trROT[0].stats.starttime + parrivaltime- t_d
            #~ t2 = t1+ round(MinDist[tr.id]*15. )
                           
           
            N = len(trROT[0])
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
                trR.data = bp.bandpassfilter(trR.data,len(trR), trR.stats.                                              delta,
                corners , 1 , fmin, fmax)  
                t_l = dt*0.5*(AUX1 - AUX2)                             
                trR.trim(t1-t_l,t2-t_l, pad=True, fill_value=trR.data[0])  #We lost t_h due to the convolution                  
            
            
         
            #~ for trR in trROT:
                #~ trR.data *= 10.**-23 ## To get M in Nm
                #~ trR.data -= trR.data[0]                
                #~ trR.data = convolve(trR.data,window,mode='same')
                #~ # mean = np.mean(np.hstack((trR.data[0]*np.ones(FirstValid),\
                             #~ # trR.data[:60./trR.stats.delta*1.-FirstValid+1])))
                #~ mean = np.mean(trR.data[:60])               
                #~ trR.data -= mean
                #~ trR.data = bp.bandpassfilter(trR.data,len(trR), trR.stats.delta,\
                                             #~ corners ,1 , fmin, fmax)
                #~ trR.trim(t1,t2,pad=True, fill_value=trR.data[0])  

            nmin = min(len(tr.data),len(trROT[0].data))             
            tr.data = tr.data[:nmin]
            for trR in trROT:
                trR.data = trR.data[:nmin]
              
                
             #############            
            trROT = np.array(trROT)  
            syn  =  np.dot(trROT.T,M) 
            synA =  np.dot(trROT.T,MA)
            synB =  np.dot(trROT.T,MB)
            
            SYN = np.append(SYN,syn)  
            SYNA = np.append(SYNA,synA)
            SYNB = np.append(SYNB,synB)
            
            if subf == 0 : ObservedDisp =  np.append(ObservedDisp,tr.data,0) 
            
  
        sbarea.append(grid[subf][6])
   
        print np.shape(A), np.shape(np.array([SYN]))
        if subf == 0: 
            A  = np.array([SYN])
            AA = np.array([SYNA])
            AB = np.array([SYNB])
        else:
            A = np.append(A,np.array([SYN]),0) 
            AA = np.append(AA,np.array([SYNA]),0)
            AB = np.append(AB,np.array([SYNB]),0)
        
    
    
    #Full matrix with the two rake's component
    AC = np.vstack((AA,AB))

#MISFIT = np.array([])
########## Stabilizing the solution:         

#### Moment minimization:
#~ constraintD  = np.zeros(nsf)
#~ ObservedDispcons = np.append(ObservedDisp,constraintD)
#~ for lbd in Lbd:
    #~ constraintF  = lbd*np.eye(nsf,nsf)         
    #~ Acons = np.append(A,constraintF,1)   
    #~ print np.shape(Acons.T), np.shape(ObservedDispcons)
    #~ R = nnls(Acons.T,ObservedDispcons)
    #~ M = R[0]
    #~ #M = np.zeros(nsf)
    #~ #M[::2] = 1
    #~ fit = np.dot(A.T,M)
    #~ misfit = 100.*np.sum(np.abs(fit-ObservedDisp))\
             #~ /np.sum(np.abs(ObservedDisp))
    
    #~ MISFIT = np.append(MISFIT,misfit)
#~ plt.figure()
#~ plt.plot(Lbd,MISFIT)
#~ ###########################################
#~ ### Smoothing:
#~ constraintF_base = SmoothMatrix(nsx,nsy)
#~ constraintD  = np.zeros(np.shape(constraintF_base)[0])
#~ ObservedDispcons = np.append(ObservedDisp,constraintD)
#~ for lbd in Lbd:
    #~ constraintF  = lbd*constraintF_base   
    #~ Acons = np.append(A,constraintF.T,1)   
    #~ #print np.shape(Acons.T), np.shape(ObservedDispcons)
    #~ R = nnls(Acons.T,ObservedDispcons)
    #~ M = R[0]
    #~ fit = np.dot(A.T,M)
    #~ misfit = 100.*np.sum(np.abs(fit-ObservedDisp))\
             #~ /np.sum(np.abs(ObservedDisp))
    #~ print lbd, misfit
    #~ MISFIT = np.append(MISFIT,misfit)
#~ ###########################################    
###########################################
#~ ##### Moment Minimization (including rake projections):
#~ constraintD  = np.zeros(2*nsf)
#~ ObservedDispcons = np.append(ObservedDisp,constraintD)
#~ for lbd in Lbd:
    #~ constraintF  = lbd*np.eye(2*nsf,2*nsf)         
    #~ ACcons = np.append(AC,constraintF,1)   
    #~ print np.shape(ACcons.T), np.shape(ObservedDispcons)
    #~ R = nnls(ACcons.T,ObservedDispcons)
    #~ M = R[0]
    #~ fit = np.dot(AC.T,M)
    #~ misfit = 100.*np.sum(np.abs(fit-ObservedDisp))\
             #~ /np.sum(np.abs(ObservedDisp))        
    #~ MISFIT = np.append(MISFIT,misfit)  
    #~ M = np.sqrt(M[:nsf]**2+M[nsf:]**2)

##############################################
### Smoothing (including rake projections):
#~ constraintF_base = SmoothMatrix(nsx,nsy)
#~ Nbase = np.shape(constraintF_base)[0]
#~ constraintD  = np.zeros(2*Nbase)
#~ constraintF_base_big = np.zeros((2*Nbase, 2*nsf))
#~ constraintF_base_big[:Nbase,:nsf]= constraintF_base
#~ constraintF_base_big[Nbase:,nsf:]= constraintF_base 
#~ ObservedDispcons = np.append(ObservedDisp,constraintD)
#~ for lbd in Lbd:
    #~ constraintF  = lbd*constraintF_base_big   
    #~ ACcons = np.append(AC,constraintF.T,1)   
    #~ #print np.shape(Acons.T), np.shape(ObservedDispcons)
    #~ R = nnls(ACcons.T,ObservedDispcons)
    #~ M = R[0]
    #~ fit = np.dot(AC.T,M)
    #~ misfit = 100.*np.sum(np.abs(fit-ObservedDisp))\
             #~ /np.sum(np.abs(ObservedDisp))
    #~ print lbd, misfit
    #~ MISFIT = np.append(MISFIT,misfit)
#~ M = np.sqrt(M[:nsf]**2+M[nsf:]**2)    
###########################################    
#~ ##### Moment Minimization and  Smoothing
  #~ #### (including rake projections):
    #~ mom0 = []
    #~ constraintF_base = SmoothMatrix(nsx,nsy)
    #~ Nbase = np.shape(constraintF_base)[0]
    #~ constraintDsmoo = np.zeros(2*Nbase)
    #~ constraintDmin  = np.zeros(2*nsf)
    #~ constraintF_base_big = np.zeros((2*Nbase, 2*nsf))
    #~ constraintF_base_big[:Nbase,:nsf]= constraintF_base
    #~ constraintF_base_big[Nbase:,nsf:]= constraintF_base 
    #~ ObservedDispcons = np.concatenate((ObservedDisp,
                                  #~ constraintDmin,
                             #~ constraintDsmoo  ))    
   
    #~ for lbdm0 in Lbdm0min:
        #~ constraintFmin  = lbdm0*np.eye(2*nsf,2*nsf)
        #~ for lbdsm in Lbdsmooth:              
            #~ constraintFsmoo  = lbdsm*constraintF_base_big 
            #~ ACcons = np.hstack((AC, constraintFmin, constraintFsmoo.T))   
            #~ print lbdm0, lbdsm
            #~ R = nnls(ACcons.T,ObservedDispcons)
            #~ M = R[0]
            #~ fit = np.dot(AC.T,M)
            #~ misfit = 100.*np.sum(np.abs(fit-ObservedDisp))\
                     #~ /np.sum(np.abs(ObservedDisp))        
            #~ MISFIT = np.append(MISFIT,misfit) 
            #~ MA = M[:nsf]
            #~ MB = M[nsf:]
            #~ M = np.sqrt(MA**2+MB**2)
            #~ mom0.append(np.sum(M))
    ##############################################

    # Rotation to the rake's conventional angle:  
    #MB, MA = Rot2D(MB,MA,-rakeB)
    print np.shape(M), np.shape(A.T)
    R = nnls(A.T,ObservedDisp)
    M = R[0]
    
    #~ M = np.zeros(nsf)
    #~ M[::2] = 1 
    fit = np.dot(A.T,M)
    MA = M
    MB = M
    
    np.save("RealSol", M)
      
    nm0 = np.size(Lbdm0min) 
    nsmth = np.size(Lbdsmooth)
    #~ plt.figure()
    #~ plt.pcolor(1./Lbdsmooth, 1./Lbdm0min,MISFIT.reshape(nm0,nsmth))
    #~ plt.xlabel(r'$1/ \lambda_{2}$', fontsize = 24)
    #~ plt.ylabel(r'$1/ \lambda_{1}$',fontsize = 24 )
    #~ plt.ylim((1./Lbdm0min).min(),(1./Lbdm0min).max() )
    #~ plt.ylim((1./Lbdsmooth).min(),(1./Lbdsmooth).max() )
    #~ cbar = plt.colorbar()
    #~ cbar.set_label("Misfit %")
    #~ print np.shape(Lbdm0min), np.shape(mom0)
    
    #~ plt.figure()
    #~ CS = plt.contour(1./Lbdsmooth, 1./Lbdm0min,MISFIT.reshape(nm0,nsmth) )
    #~ plt.xlabel(r'$1/ \lambda_{2}$', fontsize = 24)
    #~ plt.ylabel(r'$1/ \lambda_{1}$',fontsize = 24 )
    #~ plt.clabel(CS, inline=1, fontsize=10)
    #~ plt.title('Misfit')
    
    
    
    #~ plt.figure()
    #~ plt.plot(1./Lbdm0min,MISFIT)
    #~ plt.xlabel(r'$1/ \lambda_{2}$', fontsize = 24)
    #~ plt.ylabel("Misfit %")
    #~ plt.figure()
    #~ plt.plot(Lbdm0min,mom0)
    #~ plt.ylabel(r'$M_0\, [Nm]$', fontsize = 24)
    #~ plt.xlabel(r'$\lambda_{M0}$', fontsize = 24)

   
    misfit = 100.*np.sum(np.abs(fit-ObservedDisp))/np.sum(np.abs(ObservedDisp))
    print "Residual: ", 1000.*R[1]
    print misfit
    
    
    #SLIP = M*Mp/mu/(1.e6*np.array(sbarea))
    
    sbarea = sflen*sfwid
    SLIP = M/(mu*1.e6*sbarea)
    SLIP = SLIP.reshape(nsx,nsy).T[::-1]
    moment = M.reshape(nsx,nsy).T[::-1]
    
    plt.figure(figsize = (13,5))
    plt.plot(fit,'b' ,label="Fit")
    plt.plot(ObservedDisp,'r',label="Observed")
    plt.xlabel("Time [s]")
    plt.ylabel("Displacement [m]")
    plt.legend()
    
    
    np.set_printoptions(linewidth=1000,precision=3)
    print "***********"
    print sbarea
    print SLIP
    print np.mean(SLIP)
    print "Moment:"
    print np.sum(M)
 

    ### SLIPS Distribution (as the synthetics) :
    SLIPS = M.reshape(nsx,nsy).T
    SLIPS /=  mu*1.e6*sbarea
    
    
    #~ #########Ploting slip distribution:
    #~ #we are going to reflect the y axis later, so:
    hypsbloc = [hyp_subf / nsy , -(hyp_subf % nsy) - 2]

    #Creating the strike and dip axis:
    StrikeAx= np.linspace(0,flen,nsx+1)
    DipAx= np.linspace(0,fwid,nsy+1)
    DepthAx = DipAx*np.sin(np.pi/180.*dip) + Min_h
    print DepthAx
    hlstrike = StrikeAx[hypsbloc[0]] + sflen*0.5
    #we are going to reflect the axis later, so:
    hldip = DipAx[hypsbloc[1]] + sfwid*0.5 
    hldepth = DepthAx[hypsbloc[1]] + sfwid*0.5*np.sin(np.pi/180.*dip)    
    
    StrikeAx = StrikeAx - hlstrike
    DipAx =     DipAx   - hldip
    
    XX, YY = np.meshgrid(StrikeAx, DepthAx)
    XX, ZZ = np.meshgrid(StrikeAx, DipAx )  

   ######Plot: (Old colormap: "gist_rainbow_r")
    plt.figure(figsize = (13,6))
    ax = host_subplot(111)
    im = ax.pcolor(XX, YY, SLIPS, cmap="jet")    
    ax.set_ylabel('Depth [km]')       
    ax.set_ylim(DepthAx[-1],DepthAx[0])  
    
    # Creating a twin plot 
    ax2 = ax.twinx()
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

    
    #~ ### Rake plot:
    plt.figure(figsize = (13,6))
    fig = host_subplot(111)
    XXq, ZZq = np.meshgrid(StrikeAx[:-1]+sflen, DipAx[:-1]+sfwid )
    Q = plt.quiver(XXq,ZZq, MB.reshape(nsx,nsy).T[::-1,:]/(mu*1.e6*sbarea), 
                    MA.reshape(nsx,nsy).T[::-1,:]/(mu*1.e6*sbarea),
                    SLIPS[::-1,:],
                units='xy',scale = 0.5  ,  linewidths=(2,), 
                edgecolors=('k'), headaxislength=5  )
    fig.set_ylim([ZZq.min()-80,ZZq.max()+80])
    fig.set_xlim([XXq.min()-20, XXq.max()+20 ])
    fig.set_ylabel('Distance along dip [km]') 
    fig.set_xlabel('Distance along the strike [km]') 
    
    fig2 = fig.twinx()
    fig2.set_xlabel('Distance along the strike [km]') 
    
    fig.axis["bottom"].major_ticklabels.set_visible(False) 
    fig.axis["bottom"].label.set_visible(False)
    fig2.axis["top"].set_visible(True)
    fig2.axis["top"].label.set_visible(True)
    fig2.axis["right"].major_ticklabels.set_visible(False)

    divider = make_axes_locatable(fig)
    cax = divider.append_axes("bottom", size="5%", pad=0.1)
    cb = plt.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Slip [m]") 
    

    
    
    plt.show()

    
        #############
    #~ print np.shape(MISFIT),  np.shape(RUPVEL)
    #~ plt.figure()
    #~ plt.plot(RUPVEL,MISFIT)
    #~ plt.xlabel("Rupture Velocity [km/s]")
    #~ plt.ylabel("Misfit %")
    #~ plt.show()
     

    print np.shape(MB.reshape(nsx,nsy).T)
    print np.shape(ZZ)

        



def GFSelectZ(stlat,stlon,hdir): 
    dist = locations2degrees(eplat,eplon,stlat,stlon) 
    
    
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

    
    
def CosTaper(win):
    taper = 0.5*(1.-np.cos(np.pi*np.linspace(0.,1.,win)))
    #trace[:win] *= taper

    return taper    
    
    
    
     
def lstargsort(seq):
   #http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
   #by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)   
    
    
    
def MTrotationZ(azi, trPPsy,  trRRsy, trRTsy,  trTTsy):
    #The result follow the order of NodalPlanetoMT
    
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



def degrees2kilometer(deg,r=6371.):
    return deg*r*np.pi/180.   
    
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
    
def SmoothMatrix(nsx,nsy):
    ''' Given the vertical and horizontal number of subfaults it 
        return a matrix F with the smoothing conditions. This matrix 
        should be append the Synthetics matrix before the inversion'''
        
    nsf =  nsx*nsy
    F = np.zeros((1,nsf))
    i = 0
    for subA in range(nsf):
        for subB in range(subA, nsf):
            if neighbour(subA, subB,(nsx,nsy)):
            #if subA >= 18 and subB >= 15:
                row = np.zeros((1,nsf))
                row[0, subA] = 1.
                row[0, subB] = -1.
                F = np.append(F, row, 0 )
                i += 1
    F = F[1:,:]
    
    return F
        
def neighbour(subi,subj,shape):
    '''Given 2 subfaults numbers and the shape of the subfaults'
       arrange it returns whether the subfaults are neighbours or 
       not. shape should be a tuple'''     
    
    nb = False
    
    indexi = [subi / shape[1], subi % shape[1]]
    indexj = [subj / shape[1], subj % shape[1]]
    
    #Diagonal neighbours cond
    nbconddiag = abs(indexi[0]-indexj[0]) == 1 and\
                abs(indexi[1]-indexj[1]) == 1
             
    # Direct neighbours con
    nbconddir = (abs(indexi[0]-indexj[0]) == 0 and\
                 abs(indexi[1]-indexj[1]) == 1) or\
                (abs(indexi[1]-indexj[1]) == 0 and\
                 abs(indexi[0]-indexj[0]) == 1)
    if (nbconddiag or nbconddir):
        nb = True
       
      
    return nb 
def Rot2D(x,y,angle):
    '''
    Given 2 vector x,y (same length) it performs a couterclockwise
    2d-rotation in a angle "angle" (degrees) for earch pair of elements       
    '''
    ang = -angle*np.pi/180.
    xr =  x*np.cos(ang) - y*np.sin(ang)
    yr =  x*np.sin(ang) + y*np.cos(ang)
    return xr,yr

    
if __name__ == "__main__":
    main()