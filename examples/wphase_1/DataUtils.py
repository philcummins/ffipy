import bpfilter as bp
from obspy.core import read, Trace, Stream, UTCDateTime
import matplotlib.pyplot as plt 
import numpy as np
from obspy.xseed import Parser
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.signal.rotate import rotate_NE_RT
from obspy.taup.taup import getTravelTimes,locations2degrees
import os
from scipy.integrate import cumtrapz
from obspy.signal.filter import bandpass
from obspy.signal.invsim import paz2AmpValueOfFreqResp
from scipy.signal import boxcar,triang, convolve, detrend
from subprocess import PIPE, Popen


def get_azimuths(sp):
 ''' Return a Dictionary with 
        id, azimuth
        for each station in a Parser file '''
 blk = sp.blockettes
 AZI = {}
 nst = 0
 for i in range(len(blk[50])):
    nchan = blk[50][i].number_of_channels
    Netcod = blk[50][i].network_code
    stcod  = blk[50][i].station_call_letters
    for j in range(nst,nst+nchan): 
      Locid = blk[52][j].location_identifier
      Chanid= blk[52][j].channel_identifier
      azi = blk[52][j].azimuth
      stid =  '.'.join([Netcod, stcod, Locid, Chanid])
      AZI[stid]= azi  
      nst += 1 
    
 return AZI

 
 
def rot_12_NE(st, AZI): 
 '''Perform a  12 -> NE rotation:'''       
 
 st2 = st.select(channel="??1")
 for tr in st2:
  id1 = tr.id 
  id2 = '.'.join([tr.stats.network,tr.stats.station,tr.stats.location,tr.stats.channel[:-1]+"2"])
  tr1 = tr
  tr2 = st.select(id=id2)[0]
  timeA = max(tr1.stats.starttime,tr2.stats.starttime)
  timeB = min(tr1.stats.endtime,tr2.stats.endtime)  
  tr1.trim(timeA,timeB)
  tr2.trim(timeA,timeB)   
  azi = AZI[id1] * np.pi/180.
  C1 = tr1.data.copy() 
  C2 = tr2.data.copy()
  tr1.data = C1*np.cos(azi)-C2*np.sin(azi)
  tr2.data = C1*np.sin(azi)+C2*np.cos(azi)
 return st
 
 

def PAZfromBLK(sp):
  
 ''' 
 Return a  dictionary which contains the channel id and a 
 tuple with the poles and zeros of Parser file
 ''' 

 blk = sp.blockettes
 PAZ = {}
 nst = 0
 k=0

# Getting rid of 'empty' blockettes  
 for i in range(len(blk[53])):
  
    Empty= blk[53][i].number_of_complex_poles==0 and\
  	 blk[53][i].number_of_complex_zeros==0
  
    if  Empty:
      blk[53][i]= "empty blockette"
 blk[53][:] = (A for A in  blk[53] if A != "empty blockette")   


#Building the dictionary
 for i in range(len(blk[50])):
    nchan = blk[50][i].number_of_channels
    Netcod = blk[50][i].network_code
    stcod  = blk[50][i].station_call_letters
    for j in range(nst,nst+nchan): 
      Locid = blk[52][j].location_identifier
      Chanid= blk[52][j].channel_identifier
      azi = blk[52][j].azimuth
      stid =  '.'.join([Netcod, stcod, Locid, Chanid])
      RP = blk[53][j].real_pole
      IP = blk[53][j].imaginary_pole
      RZ = blk[53][j].real_zero
      IZ = blk[53][j].imaginary_zero

      RP=np.array(RP)
      IP=np.array(IP)
      RZ=np.array(RZ)
      IZ=np.array(IZ) 

      poles= RP + 1j*IP
      zeros= RZ + 1j*IZ
      PAZ[stid]= (poles, zeros)  
     
      nst+=1
      
 return PAZ 
 
def PAZfromRESP(respfile):

 polesid = "B053F15-18"
 zeroesid= "B053F10-13"
 gainid =  "B053F07"
 sensid = "B058F04"
 

 cmdp =  "cat " + respfile + " | grep " + polesid + " |awk '{print $3, $4 }' "
 cmdz =  "cat " + respfile + " | grep " + zeroesid + " |awk '{print $3, $4 }' " 
 cmdg =  "cat " + respfile + " | grep " + gainid + " |awk '{print $5 }' "
 cmds =  "cat " + respfile + " | grep " + sensid + " |awk '{print $3 }' "
 
 #Popen(cmd, shell=True,close_fds=True, stdout=PIPE).stdout

 p = Popen(cmdp, shell=True,close_fds=True, stdout=PIPE).stdout
 z = Popen(cmdz, shell=True,close_fds=True, stdout=PIPE).stdout
 g = Popen(cmdg, shell=True,close_fds=True, stdout=PIPE).stdout
 s = Popen(cmds, shell=True,close_fds=True, stdout=PIPE).stdout
 poles = np.loadtxt(p).view(complex).reshape(-1)
 zeros= np.loadtxt(z).view(complex).reshape(-1) 
 gain =  np.array(np.loadtxt(g))
 #Some respfiles have a 2nd A0 = 1.(We need the 1st one)
 if gain.size > 1. : gain = gain[0]  
 gain = np.float(gain)
 
 sens = np.loadtxt(s)[-1] #The last one is the overall sensitivity 
 
 PAZ = {}
 PAZ['poles'] = poles
 PAZ['zeros'] = zeros
 PAZ['gain']  = gain
 PAZ['sensitivity'] = sens
 
 p.close()
 z.close()
 g.close()
 s.close() 
 
 
 return PAZ

def RTdeconv(data, om0, h, G, dt,corners=4, baselinelen=60., taperlen= 10.\
             , fmin = 0.001, fmax = 0.005 ):
 '''  Using the coef. of one station computes 
     the displacement trace.  
     Input: 
     data: Raw data to be deconvoluted
     om0, h, G: Natural frequency, damping constant, sensitivity of the seismometer 
     dt : sampling rate (sec)
     corners : order of the Butterworth filter to be applied.      
     baselinelen : time (in Sec ) to consider to define the baseline.
     taperlen:  Percentage of the trace to be tapered at the beginning.
     fmin, fmax: Corner frequencies of the Butterworth filter to be applied.
     Output
     Time serie with the deconvoluted displacement trace     
               '''

 baseline = np.mean(data[:int(baselinelen/dt)])
 nwin = int(10.*len(data)/100.)
 
 if len(data) < 2*nwin: nwin = len(data)
 
 data -= baseline 
 taper = 0.5*(1.-np.cos(np.pi*np.linspace(0.,1.,nwin)))
 data[:nwin]*= taper  
 #data[-nwin:]*= taper[::-1]
  
 datap1 = data[1:] 
 datap2 = data[2:] 
 
 ## Acceleration + double integration
 c0 = 1./G/dt
 c1 = -2.*(1. + h*om0*dt)/G/dt
 c2 = (1. + 2.*h*om0*dt + dt*dt*om0*om0)/G/dt 
 
 accel = np.zeros(len(data))
 aux = c2*datap2 + c1*datap1[:-1] + c0*data[:-2] 
 accel[2:] = np.cumsum(aux) 
 
 accel = bp.bandpassfilter(accel,len(accel),dt ,corners , 1 , fmin,fmax)  

 vel = np.zeros(len(data)) 
 vel[1:] = 0.5*dt*np.cumsum(accel[:-1]+accel[1:])
 
 dis = np.zeros(len(data))
 dis[1:] = 0.5*dt*np.cumsum(vel[:-1]+vel[1:])
 return dis 
 



def testRESP(respfile):

 TransFuncsid = "B053F03"
 Respinid= "B053F05"
 

 cmdtrans =  "cat " + respfile + " | grep " + TransFuncsid + " |awk '{print $5 }' "
 cmdResp =  "cat " + respfile + " | grep " +  Respinid + " |awk '{print $6 }' "
 
 t = Popen(cmdtrans, shell=True,close_fds=True, stdout=PIPE).stdout
 r = Popen(cmdResp , shell=True,close_fds=True, stdout=PIPE).stdout
 
 trans = t.readline()[:-1] # We are checking just in the first line with the pattern
 Resp  = r.readline()[:-1]
 
 t.close()
 r.close()
 
 return (trans, Resp)
 
def getCOEFF(respfile):

  paz=PAZfromRESP(respfile)
  poles  = paz['poles']
  zeroes = paz['zeros']
  poles = poles[np.argsort(np.abs(poles))]
  C = np.abs(poles[2:].prod()/zeroes[2:].prod())
  om0 =  np.sqrt(np.abs(poles[0]*poles[1]))
  G = C*paz['sensitivity']/paz['gain'] 
  
  #print C, paz['sensitivity'], paz['gain'] 
  h = np.real((-poles[0]-poles[1])/2./om0)
 
 #That's wrong! In case TF=B, one has to divide omega by 2*pi 
  trans, Resp = testRESP(respfile)
  if trans == "B": om0 = om0*2.*np.pi # transforming hz to rad*hz
  if h > .9: G, om0, h = np.NaN, np.NaN, np.NaN
  
  return (om0, h, G)
  
def getCOEFFfit(sens,freq,resp):

  
 
  X = (2.*np.pi*freq)**2
  Y = X*X*(1./resp/resp-1.)
  fit = np.polyfit(X,Y,1)
  
  G = sens
  om0 = fit[1]**(.25)
  h = np.sqrt(fit[0]/4./np.sqrt(fit[1]) + 0.5)
  return (om0, h, G)
  
def getMetadataFromXML(file):
  ''' 
  Given a stationXML file this function returns a dictionary containing every station id and:
  Latitude, Longitude, Type of Transfer fuction, Elevation, Azimuth, Dip, poles, zeros, 
  gain (normalization factor) and sensitivity. 
  '''
  from xml.dom import minidom  
  xmlfile = minidom.parse(file)

  stlist = xmlfile.getElementsByTagName('Station')
  metadata= {}
  for st in stlist:
      netcode = st.attributes["net_code"].value
      stcode  = st.attributes["sta_code"].value
      chlist = st.getElementsByTagName('Channel')
      for ch in chlist:
          loccode =  ch.attributes["loc_code"].value.strip()
          chcode  =  ch.attributes["chan_code"].value
          stid    = '.'.join([netcode,stcode,loccode, chcode])
          azi = float(ch.getElementsByTagName('Azimuth')[0].firstChild.data)    
          lat = float(ch.getElementsByTagName('Lat')[0].firstChild.data)
          lon = float(ch.getElementsByTagName('Lon')[0].firstChild.data)
          tf = ch.getElementsByTagName('PzTransferFunctionType')[0].firstChild.data
          tf = TransferFunctionType(tf)
          Elev = float(ch.getElementsByTagName('Elevation')[0].firstChild.data)
          dip = float(ch.getElementsByTagName('Dip')[0].firstChild.data)
          gain= float(ch.getElementsByTagName('NormalizationFactor')[0].firstChild.data)
          sens = float(ch.getElementsByTagName('SensitivityValue')[0].firstChild.data)       
        
        
          polelist = ch.getElementsByTagName('Pole')
          poles = []        
          for pole in polelist:
              real = pole.getElementsByTagName('Real')[0].firstChild.data
              real = float(real)
              imag =  pole.getElementsByTagName('Imaginary')[0].firstChild.data
              imag = float(imag)
              poles.append(real + 1j*imag)
            
            
          zerolist = ch.getElementsByTagName('Zero')
          zeros = []       
          for zero in zerolist:
              real = zero.getElementsByTagName('Real')[0].firstChild.data
              real = float(real)
              imag =  zero.getElementsByTagName('Imaginary')[0].firstChild.data
              imag = float(imag)
              zeros.append(real + 1j*imag) 
        
          metadata[stid] = { 'latitude'       : lat,
                             'longitude'      : lon,
                             'elevation'      : Elev,
                             'transfer_function': tf,
                             'azimuth'        : azi,
                             'dip'            : dip,
                             'zeros'          : zeros,
                             'poles'          : poles,
                             'gain'           : gain,
                             'sensitivity'    : sens            
                           }
       
  return metadata
  
  
  
def LstArgSort(seq):
    ''' Sorts a list and return the arguments '''
   #http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
   #by ubuntu
    return sorted(range(len(seq)), key=seq.__getitem__)   
 
 
def TransferFunctionType(transfertext):
    if transfertext == "LAPLACE (RAD/SEC)":
        Type = "A"
    elif transfertext == "LAPLACE (HZ)":
        Type = "B"
    else:
        Type = "Unknown"   
        
    return Type
    


def replace_pattern(file, pattern, subst):
    ''' Using a temporary file it replaces a regular expresion in a 
        given file.
    '''
    
    from tempfile import mkstemp
    from shutil import move
    from os import remove, close
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file)
    #Move new file
    move(abs_path, file)
    
def Rot2D(x,y,angle):
    '''
    Given 2 vector x,y (same length) it performs a couterclockwise
    2d-rotation in a angle "angle" (degrees) for earch pair of elements       
    '''
    ang = -angle*np.pi/180.
    xr =  x*np.cos(ang) - y*np.sin(ang)
    yr =  x*np.sin(ang) + y*np.cos(ang)
    return xr,yr
