c      program testpath
c      real c(11),t(11)
c      nt=11
c      data t /20,30,40,50,60,80,100,125,150,175,200/
c      call readhrv()
c      ela=20.
c      elo=18.
c      sla=34.
c      slo=-118.
c      call phvpath(ela,elo,sla,slo,t,nt,1,c)
c      do i=1,nt
c	print*, t(i),c(i)
c      enddo
c      end

      subroutine phvpath(evla,evlo,stla,stlo,per,nper,ityp,c,u,q)
c
c     -- computes path averaged phasevelocities for the Harvard model --
c
      real per(*),c(*),u(*),q(*)
      real ttime(4096)
      real*8 ela,elo,sla,slo,baz,az,dist,dsegment
      real*8 grla,grlo,graz
      ela=evla
      elo=evlo
      sla=stla
      slo=stlo
      if(nper .gt.4096) then
	print*,'phvpath: nper is too large'
	stop
      endif
      do j=1,nper
	ttime(j)=0.0
      enddo
c     -- compute path averaged phasevelocity ---------------------------
c
c     -- dr is the target spacing for the segments ---------------------
      dr=1.0
      call thetasr(ela,elo,sla,slo,baz,az,dist)
      nsegments=dist/dr
      dsegment=dist/nsegments
c     print*,nsegments,dsegment
      aseg=0.0
      do i=0,nsegments
	call gettp(ela,elo,az,dsegment*i,grla,grlo,graz)
	call phv_hrvt(grla,grlo,per,nper,ityp,c,u,q)
	do j=1,nper
	  if(i.eq.0 .or.i.eq.nsegments) then
	    ttime(j)=ttime(j)+0.5*111.*dsegment/c(j)
	    aseg=aseg+0.5*111.*dsegment
	  else
	    ttime(j)=ttime(j)+111.*dsegment/c(j)
	    aseg=aseg+111.*dsegment
	  endif
        enddo
      enddo
c     print*,111.*dist,aseg/11.
c     print*,grla,grlo
c     print*,sla,slo

      do j=1,nper
	c(j)=111.*dist/ttime(j)
      enddo
      return
      end
c     ------------------------------------------------------------------

      subroutine thetasr(tts,pps,ttr,ppr,bazm,azm,dist)
c     ------------------------------------------------------------------
c     -- back azimuth, azimuth and distance between source and receiver:
c     -- source coordinate (ts,ps,), receiver coordinate (tr,pr); ------
c     -- back azimuth (azimuth of source from receiver); ---------------
c     -- azimuth (azimuth of receiver from source);  -------------------
c     -- They are all in degrees, measured clockwise from north. -------
      implicit real*8 (a-h,o-z)
      real*8,intent(in):: tts,pps,ttr,ppr
      real*8,intent(out):: bazm,azm,dist
      double precision x,xd,dph,dis,sang,cang,rdph,rps,rts,rtr
      if(tts.eq.ttr.and.pps.eq.ppr) then
	 dist=0.0
	 azm=0.0
	 bazm=0.0
	 return
      endif
        
      if(pps.lt.0.0) then
         ps=pps+360.
        else
         ps=pps
        endif
      if(ppr.lt.0.0) then
         pr=ppr+360.
      else
         pr=ppr
        endif
      ts=90.-tts
      tr=90.-ttr

      conv=3.141592653589793/180.0d0
      rps=ps*conv
      rts=ts*conv
      rpr=pr*conv
      rtr=tr*conv
      xd=pr-ps
      if(xd .gt. -360.0 .and. xd .le. -180.0) then
        icase=1
        dph=xd+360.0
      elseif(xd .gt. -180.0 .and. xd .le. 0.0) then
        icase=2
        dph= -xd
      elseif(xd .gt. 0.0 .and. xd .le. 180.0) then
        icase=3
        dph=xd
      elseif(xd .gt. 180.0 .and. xd .le. 360.0 ) then
        icase=4
        dph = 360.0 - xd
      endif

      rdph=dph*conv
      x=dcos(rts)*dcos(rtr) + dsin(rts)*dsin(rtr)*dcos(rdph)
      dis=dacos(x)
      dist=dis/conv

      if(dis.ne.0.0d0) then
         sang=dsin(rts)*dsin(rdph)/dsin(dis)
         cang=(dcos(rts)-x*dcos(rtr))/(dsin(dis)*dsin(rtr))
         bazm=datan2(sang,cang)/conv
         if(icase .eq. 2 .or. icase .eq. 4) bazm= -bazm

         sang=dsin(rtr)*dsin(rdph)/dsin(dis)
         cang=(dcos(rtr)-x*dcos(rts))/(dsin(rts)*dsin(dis))
         azm=datan2(sang,cang)/conv
         if(icase .eq. 1 .or. icase .eq. 3) azm = -azm
c       -- change of convention: clockwise from north
         azm = -azm
         bazm = -bazm
      else
	 azm=0.0
	 bazm=0.0
      endif
      return
      end
c
      subroutine gettp(tthe,pphi,azm,dis,athe,aphi,aazm)
c     -- for a given (the,phi) and azm and dis, calculates (athe,aphi) -
c     -- and azimuth azm (actaully baz) at the new point.  -------------
c     -- Note that if azm=0.0, this routine causes zero divide. --------
      								
      double precision  tthe,the,phi,pphi,azm,dis,athe,aphi,aazm
      double precision  rthe,rphi,rdis,ct,cd,st,sd,razm,cb,ctp,stp
      double precision  sb,sc,sa,cc,cang,ca,aang
      the=90.-tthe
      phi=pphi
      if(phi.lt.0.0) phi=phi+360.
      conv=3.141592653589793/180.0d0
      rthe=the*conv
      rphi=phi*conv
      rdis=dis*conv
      ct=dcos(rthe)
      cd=dcos(rdis)
      st=dsin(rthe)
      sd=dsin(rdis)
      razm=azm*conv
      cb=dcos(razm)
      
      ctp=ct*cd+st*sd*cb
      athe=dacos(ctp)/conv
      if(athe .lt. 0.0) athe=athe+180.
      stp=dsin(athe*conv)
      
      sb=dsin(razm)
      sc=st*sb/stp
      sa=sd*sb/stp
      
      if(the .eq.90.0) then
        cc = -cd*sc*cb/sb
        cang=datan2(sc,cc)/conv
        ca= -cb*cc+sb*sc*cd
        aang=datan2(sa,ca)/conv
        aphi = phi + aang
      else
        if(the .ne. 90.0) then
          ca=(ctp*st-sd*cb)/stp/ct
          aang=datan2(sa,ca)/conv
          cc=(ct*sd-st*cd*cb)/stp
          cang=datan2(sc,cc)/conv
          aphi=phi+aang
        endif
      endif
      if(aphi .lt. 0.0) aphi= aphi+360.0
      if(aphi .gt. 360.0) aphi=aphi-360.0
      aazm=360. - cang
      if(aazm .gt. 360.) aazm=aazm - 360.0
      if(aazm .lt. 0.0 ) aazm=aazm + 360.0
      athe=90.-athe
      return
      end
c
c      program testhrv
c      parameter(Ldim=40)
c
c      real*8 clm(Ldim+1,Ldim+1,18), slm(Ldim+1,Ldim+1,18),tlm(18),norm
c      real*8 lat,lon,dv,phv_hrv
c      real T(14), c(20)
c      data T /10,20,30,40,50,60,70,80,90,100,125,150,175,200/
c      common /hrvcoef/ clm,slm,tlm,norm
c
c      lat=30.
c      lon=95.
c      nt=14
c
c      call readhrv()
c      call phv_hrvt(lat,lon,T,nt,5,c)
c      do i=1,nt
c        write(*,*) T(i),c(i)
c      enddo
c      end
c
c
c
c
c     -- subroutines to read the harvard surface wave data sets --------

      subroutine readhrv()

c     -- reads the Harvard model parameters ----------------------------
c     ------------------------------------------------------------------
      implicit none

      integer*4 im, M, il, L, Ldim, maxL, minL, taper, ip, np, idx
      integer*4 nprem
      parameter(Ldim=40)
      common /hrvcoef/ clm,slm,tlm,norm

      real*8 norm, junk, pi, lat, lon, filt
      real*8 clm(Ldim+1,Ldim+1,18), slm(Ldim+1,Ldim+1,18),tlm(18)

      character*80 coeffile(18),coefpath,file, xyzfile, ajunk
      data tlm /35.,37.,40.,45.,50.,60.,75.,100.,150.,
     1          35.,37.,40.,45.,50.,60.,75.,100.,150./
      data coeffile /'L35_SH.txt','L37_SH.txt','L40_SH.txt',
     1               'L45_SH.txt','L50_SH.txt','L60_SH.txt',
     2               'L75_SH.txt','L100_SH.txt','L150_SH.txt',
     3               'R35_SH.txt','R37_SH.txt','R40_SH.txt',
     4               'R45_SH.txt','R50_SH.txt','R60_SH.txt',
     5               'R75_SH.txt','R100_SH.txt','R150_SH.txt'/

c      common /hrvcoef/ clm,slm,tlm,norm
      coefpath='/d/atws/1/hongkie/Share/Harvard-phv'
      pi = 3.1415927d0

c     ------------------------------------------------------------------

      np=18
      minL=0
      maxL=40
      taper=0
      norm=1./pi

 10   format(a80)

      idx=index(coefpath,' ')-1
      do 20 ip=1,np

	file=coefpath(1:idx)//'/'//coeffile(ip)
        open(unit=10,file=file)
	rewind 10
        do L = 0, maxL
          if ((L .gt. taper).and.(taper .gt. 0)) then
            filt = exp(log(0.01)*(L-taper)**2 / (maxL-taper)**2)
          else 
            filt = 1.0
          end if
          il = L+1
          do M = 0, L
            im = M+1
	    if (m.eq.0) then
               read(10,*) junk, junk, clm(il,im,ip)
	    else
               read(10,*) junk, junk, clm(il,im,ip), slm(il,im,ip)
	    endif
            if (L .lt. minL) then
               clm(il,im,ip) = 0.0d0
               slm(il,im,ip) = 0.0d0
            else
               clm(il,im,ip) = clm(il,im,ip) * filt
               slm(il,im,ip) = slm(il,im,ip) * filt
            end if
          end do
        end do
	close(10)
20    continue
      return
      end
c     ------------------------------------------------------------------

c     ------------------------------------------------------------------

      subroutine phv_hrvt(lat,lon,T,nt,wtype,c,u,q)

      implicit none

      integer*4 im,M,il,L,Ldim,maxL,minL,ip1,ip2, wtype,it0,ilm,i,j,nt
      integer*4 ip,nprem(12)

      parameter(Ldim=40)

      real*8 norm, junk, pi, lat, lon, filt
      real*8 dlon, dlat, dphi, phi, theta, dtheta
      real*8 Plm, A1,A2,plgndr,A(9)
      real*8 clm(Ldim+1,Ldim+1,18), slm(Ldim+1,Ldim+1,18)
      real*8 cmphi, smphi,tlm(18),costh,conv,dinterp0
      real*4 T(*),c(*),u(*),q(*),cpr,ainterp
      real*4 tprem(50,12),cprem(50,12),uprem(50,12),qprem(50,12)
      real*4 fprem(100)
      common /hrvcoef/clm,slm,tlm,norm
      data nprem /36,36,35,34,31,34,27,33,25,32,23,31/
      data (tprem(i,1),i=1,36) /
     1   926.93,819.20,737.40,619.86,538.23,477.47,430.07,375.28,333.40,
     2   300.18,265.12,237.46,215.04,192.35,173.98,156.09,141.54,127.65,
     3   114.78,104.27, 94.54, 85.67, 77.66, 70.49, 63.65, 57.68, 52.16,
     4    47.39, 43.06, 39.01, 35.44, 32.21, 29.25, 26.57, 24.14, 21.92/
      data (cprem(i,1),i=1,36) /
     1     6.64,  6.52,  6.39,  6.15,  5.95,  5.78,  5.64,  5.47,  5.34,
     2     5.23,  5.12,  5.03,  4.96,  4.90,  4.84,  4.79,  4.75,  4.72,
     3     4.68,  4.65,  4.63,  4.60,  4.58,  4.56,  4.54,  4.52,  4.50,
     4     4.48,  4.46,  4.43,  4.40,  4.37,  4.32,  4.27,  4.21,  4.15/
      data (uprem(i,1),i=1,36) /
     1     5.83,  5.53,  5.30,  4.99,  4.79,  4.66,  4.57,  4.48,  4.43,
     2     4.41,  4.39,  4.38,  4.38,  4.38,  4.38,  4.38,  4.38,  4.38,
     3     4.38,  4.38,  4.38,  4.37,  4.36,  4.36,  4.34,  4.32,  4.30,
     4     4.26,  4.22,  4.16,  4.09,  3.99,  3.88,  3.76,  3.63,  3.51/
      data (qprem(i,1),i=1,36) /
     1   205.42,195.63,187.07,173.23,162.92,155.23,149.48,143.38,139.31,
     2   136.54,134.14,132.67,131.79,131.22,131.04,131.18,131.59,132.36,
     3   133.52,134.97,136.91,139.44,142.71,146.88,152.64,160.02,170.19,
     4   183.45,201.83,228.57,265.17,314.95,378.90,449.41,514.09,561.10/
      data (tprem(i,3),i=1,35) /
     1   808.99,694.86,630.72,571.27,519.32,438.55,381.68,339.77,307.16,
     2   269.47,240.79,218.26,194.87,176.74,159.04,142.59,129.59,117.36,
     3   106.06, 95.77, 86.50, 78.21, 70.82, 64.26, 58.07, 52.67, 47.70,
     4    43.19, 39.13, 35.50, 32.25, 29.26, 26.53, 24.06, 21.84/
      data (cprem(i,3),i=1,35) /
     1    32.99, 16.46, 14.10, 12.74, 11.86, 10.74,  9.99,  9.43,  8.99,
     2     8.49,  8.11,  7.80,  7.47,  7.19,  6.90,  6.61,  6.37,  6.15,
     3     5.94,  5.77,  5.61,  5.47,  5.36,  5.26,  5.16,  5.08,  5.01,
     4     4.94,  4.88,  4.83,  4.78,  4.74,  4.71,  4.68,  4.66/
      data (uprem(i,3),i=1,35) /
     1     2.58,  5.36,  6.28,  6.85,  7.10,  6.98,  6.61,  6.33,  6.16,
     2     5.98,  5.80,  5.62,  5.37,  5.15,  4.93,  4.74,  4.63,  4.56,
     3     4.51,  4.48,  4.46,  4.45,  4.43,  4.42,  4.41,  4.40,  4.38,
     4     4.37,  4.35,  4.35,  4.35,  4.37,  4.40,  4.43,  4.45/
      data (qprem(i,3),i=1,35) /
     1   259.95,252.84,249.65,246.27,242.16,232.14,223.24,217.33,213.16,
     2   207.31,200.41,192.37,181.00,170.41,159.67,150.75,145.23,141.58,
     3   139.55,138.70,138.52,138.54,138.34,137.59,135.88,133.10,128.98,
     4   123.57,117.22,110.58,104.37, 99.02, 94.96, 92.27, 90.67/
      data (tprem(i,5),i=1,31) /
     1   457.07,402.40,363.14,323.88,289.21,260.90,228.72,204.83,186.14,
     2   166.54,151.08,135.72,121.36,110.15, 99.92, 90.82, 81.94, 74.20,
     3    67.36, 60.86, 55.20, 49.98, 45.23, 40.93, 37.06, 33.59, 30.49,
     4    27.64, 25.07, 22.76, 20.67/
      data (cprem(i,5),i=1,31) /
     1    58.39, 18.09, 14.70, 13.01, 12.04, 11.37, 10.61, 10.02,  9.56,
     2     9.07,  8.69,  8.31,  7.95,  7.65,  7.35,  7.05,  6.74,  6.46,
     3     6.22,  6.01,  5.82,  5.66,  5.51,  5.39,  5.28,  5.19,  5.12,
     4     5.06,  5.00,  4.95,  4.90/
      data (uprem(i,5),i=1,31) /
     1     1.29,  4.61,  6.09,  7.15,  7.54,  7.41,  6.97,  6.65,  6.43,
     2     6.22,  6.07,  5.91,  5.70,  5.46,  5.16,  4.88,  4.68,  4.58,
     3     4.54,  4.50,  4.47,  4.43,  4.41,  4.42,  4.44,  4.47,  4.49,
     4     4.50,  4.49,  4.46,  4.43/
      data (qprem(i,5),i=1,31) /
     1   203.93,213.10,223.40,234.49,240.96,241.68,237.47,231.54,225.21,
     2   217.49,211.51,205.53,197.40,186.62,172.75,159.46,148.64,141.93,
     3   137.78,134.88,132.78,131.03,129.63,128.68,128.20,127.96,127.50,
     4   126.26,123.84,120.38,116.67/
      data (tprem(i,7),i=1,27) /
     1   312.16,277.22,250.21,222.50,197.85,178.31,158.79,143.74,128.78,
     2   116.84,105.40, 94.88, 85.39, 76.96, 69.66, 62.98, 57.01, 51.62,
     3    46.70, 42.42, 38.54, 34.92, 31.73, 28.82, 26.17, 23.78, 21.60/
      data (cprem(i,7),i=1,27) /
     1    85.49, 19.25, 15.24, 13.33, 12.26, 11.51, 10.73, 10.13,  9.56,
     2     9.14,  8.73,  8.35,  8.01,  7.71,  7.41,  7.10,  6.78,  6.49,
     3     6.23,  6.03,  5.85,  5.69,  5.55,  5.42,  5.30,  5.20,  5.13/
      data (uprem(i,7),i=1,27) /
     1      .92,  4.39,  5.96,  7.20,  7.55,  7.16,  6.69,  6.52,  6.40,
     2     6.28,  6.10,  5.92,  5.79,  5.59,  5.27,  4.89,  4.63,  4.54,
     3     4.53,  4.53,  4.51,  4.45,  4.40,  4.37,  4.38,  4.42,  4.49/
      data (qprem(i,7),i=1,27) /
     1   215.72,227.62,237.10,244.53,242.68,232.77,222.67,219.46,218.73,
     2   216.96,212.51,207.71,204.65,198.65,184.69,165.59,150.38,141.44,
     3   136.78,134.48,132.97,131.38,129.47,127.41,125.92,125.82,127.47/
      data (tprem(i,9),i=1,25) /
     1   232.31,209.26,185.41,165.24,147.17,132.95,119.26,106.63, 96.78,
     2    87.56, 79.03, 71.32, 64.53, 58.61, 53.14, 48.04, 43.61, 39.47,
     3    35.77, 32.46, 29.43, 26.69, 24.24, 22.00, 19.96/
      data (cprem(i,9),i=1,25) /
     1   114.88, 20.14, 14.89, 13.09, 12.09, 11.36, 10.66, 10.01,  9.51,
     2     9.05,  8.66,  8.31,  8.00,  7.72,  7.42,  7.09,  6.77,  6.48,
     3     6.23,  6.03,  5.85,  5.69,  5.55,  5.42,  5.31/
      data (uprem(i,9),i=1,25) /
     1      .67,  3.98,  5.86,  7.20,  7.44,  7.07,  6.76,  6.47,  6.27,
     2     6.18,  6.13,  6.00,  5.80,  5.57,  5.21,  4.81,  4.61,  4.55,
     3     4.55,  4.55,  4.52,  4.47,  4.43,  4.41,  4.41/
      data (qprem(i,9),i=1,25) /
     1   220.86,220.88,222.84,229.27,232.24,231.39,229.30,222.95,217.81,
     2   218.00,219.86,216.59,208.80,199.50,183.95,163.61,149.39,141.08,
     3   136.61,134.10,132.52,131.61,131.45,132.00,132.74/
      data (tprem(i,11),i=1,23) /
     1   186.79,169.56,151.29,134.51,122.05,109.82, 98.84, 88.74, 80.52,
     2    72.94, 66.20, 59.73, 54.08, 49.10, 44.48, 40.26, 36.59, 33.15,
     3    30.09, 27.28, 24.74, 22.47, 20.41/
      data (cprem(i,11),i=1,23) /
     1   142.87, 20.53, 15.12, 13.23, 12.38, 11.57, 10.80, 10.14,  9.65,
     2     9.22,  8.83,  8.43,  8.09,  7.80,  7.53,  7.23,  6.90,  6.58,
     3     6.32,  6.10,  5.92,  5.76,  5.61/
      data (uprem(i,11),i=1,23) /
     1      .51,  3.78,  5.74,  7.35,  7.61,  6.97,  6.60,  6.57,  6.54,
     2     6.35,  6.06,  5.87,  5.80,  5.73,  5.49,  4.97,  4.59,  4.52,
     3     4.55,  4.57,  4.55,  4.51,  4.44/
      data (qprem(i,11),i=1,23) /
     1   199.53,204.46,218.20,236.65,240.47,230.57,224.36,229.06,233.58,
     2   229.39,218.49,211.49,208.64,205.10,196.20,174.06,154.05,143.09,
     3   137.81,135.43,134.77,134.48,133.68/
      data (tprem(i,2),i=1,36) /
     1   963.19,811.83,707.46,633.60,536.94,473.27,426.19,374.07,335.83,
     2   297.67,268.43,239.58,216.42,193.84,175.39,157.75,143.23,129.53,
     3   116.87,105.36, 95.01, 85.77, 77.55, 70.27, 63.83, 57.81, 52.55,
     4    47.73, 43.37, 39.29, 35.69, 32.42, 29.39, 26.72, 24.26, 22.05/
      data (cprem(i,2),i=1,36) /
     1     6.39,  6.57,  6.66,  6.65,  6.48,  6.27,  6.06,  5.78,  5.54,
     2     5.27,  5.06,  4.84,  4.68,  4.54,  4.43,  4.34,  4.27,  4.20,
     3     4.15,  4.11,  4.07,  4.04,  4.02,  4.00,  3.98,  3.97,  3.96,
     4     3.95,  3.94,  3.93,  3.91,  3.90,  3.89,  3.87,  3.84,  3.81/
      data (uprem(i,2),i=1,36) /
     1     7.88,  7.55,  6.95,  6.24,  5.25,  4.81,  4.53,  4.19,  3.93,
     2     3.72,  3.61,  3.57,  3.57,  3.60,  3.62,  3.65,  3.68,  3.70,
     3     3.72,  3.74,  3.76,  3.78,  3.80,  3.82,  3.83,  3.84,  3.84,
     4     3.84,  3.83,  3.81,  3.79,  3.76,  3.71,  3.65,  3.57,  3.47/
      data (qprem(i,2),i=1,36) /
     1   347.36,342.06,337.39,332.76,322.05,307.19,288.74,259.06,232.42,
     2   205.36,186.63,170.73,159.60,149.71,142.07,135.03,129.54,124.81,
     3   121.15,118.79,117.90,118.58,120.91,124.99,130.91,139.42,150.39,
     4   164.94,183.87,209.18,241.14,282.08,334.54,396.59,468.32,542.30/
      data (tprem(i,4),i=1,34) /
     1   852.62,729.78,657.01,555.77,465.46,391.38,336.05,299.53,263.58,
     2   236.16,214.18,190.87,172.52,154.56,138.11,125.41,113.72,103.09,
     3    93.46, 84.76, 76.92, 69.89, 63.15, 57.22, 52.01, 47.18, 42.76,
     4    38.79, 35.23, 31.95, 28.96, 26.26, 23.84, 21.67/
      data (cprem(i,4),i=1,34) /
     1    10.43,  9.97,  9.37,  8.47,  8.19,  8.18,  8.22,  8.10,  7.79,
     2     7.53,  7.33,  7.11,  6.93,  6.73,  6.51,  6.32,  6.12,  5.93,
     3     5.75,  5.59,  5.45,  5.33,  5.22,  5.12,  5.05,  4.98,  4.91,
     4     4.86,  4.80,  4.76,  4.71,  4.67,  4.63,  4.60/
      data (uprem(i,4),i=1,34) /
     1     8.71,  6.95,  5.44,  6.14,  7.72,  8.40,  8.22,  6.36,  5.92,
     2     5.82,  5.75,  5.64,  5.49,  5.27,  5.00,  4.77,  4.59,  4.47,
     3     4.40,  4.37,  4.36,  4.35,  4.36,  4.37,  4.37,  4.37,  4.36,
     4     4.34,  4.32,  4.30,  4.29,  4.29,  4.30,  4.31/
      data (qprem(i,4),i=1,34) /
     1   271.09,291.88,345.67,379.45,378.34,365.26,293.43,165.99,155.62,
     2   155.80,156.79,157.16,155.54,150.90,143.33,135.82,128.78,123.25,
     3   119.47,117.14,115.89,115.30,115.01,114.70,114.08,112.87,110.91,
     4   108.19,104.93,101.39, 98.03, 95.21, 93.11, 91.78/
      data (tprem(i,6),i=1,34) /
     1   805.03,725.06,594.94,536.21,448.68,388.78,344.84,308.56,273.42,
     2   244.36,220.81,192.96,174.06,156.64,139.57,126.20,113.55,102.05,
     3    91.90, 83.12, 75.52, 68.31, 61.99, 56.03, 50.85, 46.11, 41.85,
     4    38.05, 34.49, 31.32, 28.38, 25.79, 23.44, 21.29/
      data (cprem(i,6),i=1,34) /
     1    14.21, 12.27, 10.35,  9.95,  9.39,  8.95,  8.60,  8.37,  8.37,
     2     8.40,  8.43,  8.47,  8.36,  8.11,  7.86,  7.64,  7.42,  7.20,
     3     6.97,  6.74,  6.50,  6.27,  6.06,  5.88,  5.73,  5.58,  5.45,
     4     5.33,  5.22,  5.12,  5.05,  4.98,  4.93,  4.88/
      data (uprem(i,6),i=1,34) /
     1     6.10,  5.19,  7.18,  7.42,  7.07,  6.67,  6.48,  7.63,  8.64,
     2     8.73,  8.73,  8.69,  6.45,  6.32,  6.15,  5.98,  5.78,  5.55,
     3     5.26,  4.96,  4.73,  4.61,  4.57,  4.56,  4.52,  4.45,  4.37,
     4     4.32,  4.33,  4.37,  4.41,  4.43,  4.42,  4.38/
      data (qprem(i,6),i=1,34) /
     1   415.44,380.16,237.93,211.61,188.24,176.18,174.29,257.84,390.01,
     2   406.98,411.52,408.68,188.13,179.05,171.04,165.61,161.39,157.52,
     3   152.41,146.12,140.64,136.93,135.11,134.08,133.08,131.47,129.09,
     4   126.50,124.37,123.07,122.14,120.99,119.02,116.04/
      data (tprem(i,8),i=1,33) /
     1   903.98,705.62,545.46,447.52,392.21,354.65,310.44,273.44,242.39,
     2   216.92,196.14,177.92,159.47,144.48,128.41,115.57,104.23, 93.85,
     3    84.46, 76.10, 68.71, 62.31, 56.48, 51.06, 46.33, 42.07, 38.16,
     4    34.60, 31.42, 28.51, 25.89, 23.48, 21.33/
      data (cprem(i,8),i=1,33) /
     1    17.71, 16.21, 16.31, 16.26, 15.70, 13.28, 11.21, 10.10,  9.44,
     2     9.00,  8.68,  8.49,  8.51,  8.53,  8.54,  8.55,  8.26,  7.97,
     3     7.71,  7.46,  7.24,  7.02,  6.78,  6.51,  6.24,  6.00,  5.81,
     4     5.66,  5.53,  5.41,  5.31,  5.21,  5.12/
      data (uprem(i,8),i=1,33) /
     1     6.95, 16.88, 16.26, 15.79,  5.54,  5.28,  5.54,  6.06,  6.38,
     2     6.50,  6.51,  8.59,  8.68,  8.67,  8.66,  8.63,  6.13,  5.99,
     3     5.85,  5.71,  5.56,  5.29,  4.89,  4.53,  4.37,  4.38,  4.45,
     4     4.50,  4.49,  4.46,  4.41,  4.39,  4.38/
      data (qprem(i,8),i=1,33) /
     1   366.63, 90.91, 89.99, 89.92,275.53,263.60,249.88,236.20,222.38,
     2   210.03,200.05,401.73,421.66,423.24,424.67,423.73,214.05,199.72,
     3   189.49,183.77,179.41,171.59,158.83,145.69,137.14,132.54,130.66,
     4   130.33,130.49,130.57,130.21,129.19,127.78/
      data (tprem(i,10),i=1,32) /
     1   707.90,580.62,488.05,438.67,380.67,331.83,294.39,258.76,232.82,
     2   211.16,186.25,165.67,149.04,132.64,119.81,107.51, 96.32, 87.29,
     3    78.55, 70.85, 64.15, 57.86, 52.48, 47.55, 43.19, 39.18, 35.50,
     4    32.26, 29.28, 26.59, 24.11, 21.87/
      data (cprem(i,10),i=1,32) /
     1    37.70, 27.58, 23.43, 20.28, 16.18, 16.08, 16.00, 14.73, 12.74,
     2    11.49, 10.48,  9.86,  9.42,  9.01,  8.68,  8.56,  8.57,  8.57,
     3     8.29,  8.01,  7.75,  7.48,  7.23,  6.99,  6.74,  6.49,  6.25,
     4     6.04,  5.85,  5.69,  5.54,  5.41/
      data (uprem(i,10),i=1,32) /
     1    11.63, 12.85, 13.21,  4.95, 15.50, 15.34, 15.22,  5.97,  5.71,
     2     6.07,  6.55,  6.73,  6.71,  6.55,  6.36,  8.64,  8.63,  6.44,
     3     6.28,  6.03,  5.77,  5.53,  5.36,  5.14,  4.85,  4.64,  4.56,
     4     4.51,  4.47,  4.44,  4.38,  4.33/
      data (qprem(i,10),i=1,32) /
     1   355.05,434.15,480.14,290.18, 89.80, 89.68, 89.60,302.40,268.22,
     2   260.27,258.94,258.82,256.42,248.45,235.76,426.37,427.17,242.99,
     3   239.53,226.95,210.27,194.36,183.88,172.65,158.43,146.77,140.44,
     4   137.94,136.89,135.17,131.67,127.88/
      data (tprem(i,12),i=1,31) /
     1   583.50,478.17,420.25,369.91,332.15,283.64,240.56,212.95,187.61,
     2   166.78,150.58,134.10,121.09,108.61, 98.37, 88.53, 79.75, 72.29,
     3    65.32, 59.19, 53.77, 48.77, 44.27, 40.23, 36.51, 33.07, 30.00,
     4    27.21, 24.73, 22.46, 20.42/
      data (cprem(i,12),i=1,31) /
     1    45.74, 33.49, 21.17, 19.68, 18.54, 16.60, 15.85, 15.04, 13.77,
     2    12.31, 11.31, 10.47,  9.87,  9.33,  8.94,  8.61,  8.58,  8.58,
     3     8.34,  8.10,  7.88,  7.64,  7.38,  7.13,  6.87,  6.60,  6.34,
     4     6.12,  5.92,  5.74,  5.59/
      data (uprem(i,12),i=1,31) /
     1    28.52,  2.73, 13.15, 12.69, 11.84,  8.26, 15.08,  9.39,  7.36,
     2     6.40,  6.52,  6.49,  6.34,  6.34,  6.43,  6.45,  8.62,  8.61,
     3     6.35,  6.30,  6.04,  5.69,  5.44,  5.23,  4.89,  4.62,  4.58,
     4     4.51,  4.41,  4.43,  4.47/
      data (qprem(i,12),i=1,31) /
     1    90.75,317.86,489.18,502.49,506.32,418.20, 89.48,385.70,345.56,
     2   284.69,268.85,253.87,240.42,234.78,238.17,242.73,428.08,428.30,
     3   248.95,251.80,242.51,222.23,202.12,187.62,171.50,156.04,146.18,
     4   138.72,132.87,130.79,131.43/

      pi = dacos(-1.0d0)
      conv=pi/180.d0
      maxL=40

      do j=1,9
	A(j)=0.0d0
      enddo
      do j=1,nprem(wtype)
	fprem(j)=1./tprem(j,wtype)
      enddo


      phi = lon * conv
      theta=(90.-lat)*conv
      costh=dcos(theta)

c     -- set love or rayleigh ------------------------------------------
      if (wtype.eq.1) then
	it0=0
      elseif (wtype.eq.2) then
	it0=9
      endif

c     ------------------------------------------------------------------
      if(wtype .lt. 3) then
        do L = 0, maxL
          do M = 0, L
            il = L+1
            im = M+1
	    Plm=plgndr(l,m,costh)
            cmphi = dcos(dble(M)*phi)
            smphi = dsin(dble(M)*phi)
	    do ip=1,9
              A(ip)=A(ip)+
     1	          Plm*(clm(il,im,ip+it0)*cmphi+slm(il,im,ip+it0)*smphi)
            end do
          end do
        end do
      endif

c     -- determine the period brackets ---------------------------------
      do j=1,nt
	cpr=ainterp(fprem,cprem(1,wtype),1./T(j),nprem(wtype))
	u(j)=ainterp(fprem,uprem(1,wtype),1./T(j),nprem(wtype))
	q(j)=ainterp(fprem,qprem(1,wtype),1./T(j),nprem(wtype))
        if(wtype .lt. 3) then
          c(j)=dinterp0(tlm,A,1.d0*T(j),9)
          c(j)=cpr*(1.+c(j)/100.)
        else
          c(j)=cpr
	endif
      enddo

      return
      end
c     ------------------------------------------------------------------
      real function ainterp(x,y,x0,nx)
 
c     linear interpolation/extrapolation. x should be monotonous
c     otherwise it should work fine.
c 
c                             Hong Kie Thio, January, 1996
      real x(*),y(*)
      j=0
      isign=1
      if(x(1).gt.x(nx)) isign=-1
      do 10, i=1,nx
	if(isign*(x0-x(i)).gt..0)j=j+1
10    continue
      if(j.ge.nx) then
  	 ainterp=y(nx)+(y(nx)-y(nx-1))*(x0-x(nx))/(x(nx)-x(nx-1))
      else
	 if(j.eq.0) j=1
         ainterp=y(j)+(y(j+1)-y(j))*(x0-x(j))/(x(j+1)-x(j))
      endif
      return
      end

      real*8 function dinterp(x,y,x0,nx)
 
c     linear interpolation/extrapolation. x should be monotonous
c     otherwise it should work fine.
c 
c                             Hong Kie Thio, January, 1996
      real*8 x(*),y(*),x0
      j=0
      isign=1
      if(x(1).gt.x(nx)) isign=-1
      do 10, i=1,nx
	if(isign*(x0-x(i)).gt..0)j=j+1
10    continue
      if(j.ge.nx) then
  	 dinterp=y(nx)+(y(nx)-y(nx-1))*(x0-x(nx))/(x(nx)-x(nx-1))
      else
	 if(j.eq.0) j=1
         dinterp=y(j)+(y(j+1)-y(j))*(x0-x(j))/(x(j+1)-x(j))
      endif
      return
      end


      real function ainterp0(x,y,x0,nx)
 
c     linear interpolation only. x should be monotonous
c     otherwise it should work fine.
c 
c                             Hong Kie Thio, January, 1996
      real x(*),y(*)
      j=0
      isign=1
      if(x(1).gt.x(nx)) isign=-1
      do 10, i=1,nx
	if(isign*(x0-x(i)).gt..0)j=j+1
10    continue
      if(j.ge.nx) then
  	 ainterp0=y(nx)
      else if(j.eq.0) then
         ainterp0=y(1)
      else
         ainterp0=y(j)+(y(j+1)-y(j))*(x0-x(j))/(x(j+1)-x(j))
      endif
      return
      end

      real*8 function dinterp0(x,y,x0,nx)
 
c     linear interpolation only. x should be monotonous
c     otherwise it should work fine.
c 
c                             Hong Kie Thio, January, 1996
      real*8 x(*),y(*),x0
      j=0
      isign=1
      if(x(1).gt.x(nx)) isign=-1
      do 10, i=1,nx
	if(isign*(x0-x(i)).gt..0)j=j+1
10    continue
      if(j.ge.nx) then
  	 dinterp0=y(nx)
      else if(j.eq.0) then
         dinterp0=y(1)
      else
         dinterp0=y(j)+(y(j+1)-y(j))*(x0-x(j))/(x(j+1)-x(j))
      endif
      return
      end

      function ainterpl(xx,yy,xx0,nx)
 
c     logarithmic interpolation of function y(x(i))
c     where both x should be mononous and both x and y should
c     be positive. datapoint outside the interval x(1)-x(nx)
c     are extrapolated, using the nearest dy/dx.
c     
c                                 Hong Kie Thio, January 1996

      real xx(*),yy(*),x(10000),y(10000)
      if(xx0.le..0.or.nx.gt.10000) then
	 print*,'Negative x0 or too many points, no can do'
	 stop
      else
         do 5, i=1,nx
	    if(xx(i).le..0 .or. yy(i).le. .0) then
	       print*,'Negative x or y, cannot take log',xx(i),yy(i)
	       stop
            else
	       x(i)=log10(xx(i))
	       y(i)=log10(yy(i))
	    endif
5        continue
      endif
      x0=log10(xx0)
      ainterpl=ainterp(x,y,x0,nx)
      if(ainterpl.lt.-36.) then
	 ainterpl=-36.
	 print*,'Underflow, set to 1e-36'
      elseif(ainterpl.gt.36) then
	 ainterpl=36.
	 print*,'Overflow, set to 1e36'
      endif
      ainterpl=10.**ainterpl
      return
      end
      real*8 function plgndr(l,m,x)
      implicit real*8 (a-h,o-z)
      if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1.)pause 'bad arguments'
      pmm=1.0d0
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      do 15 i=l-m+1,l+m,1
        plgndr=plgndr/sqrt(i*1.)
15    continue
      anorm=dsqrt(((2*l+1)/(4*3.1415927d0)))
      plgndr=plgndr*anorm
      return
      end
