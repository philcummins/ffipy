c*********************************************************************
c      read(5,*) ms,t1,t2
c        read(1,*) az(js),az2(js),del(js),p(js),g(js),ix0(js)
c        read(1,*) im(js),ib(js),ic(js)
c     -- read velocity model --------------------------------------------
c      read(2,'(a40)') dsn
c      read(2,*) tqp,tqs,nl,(vp(l),vs(l),den(l),dep(l),l=1,nl)
c      read(2,*) nl1,(vp1(l),vs1(l),den1(l),dep1(l),l=1,nl1)
c      read(2,*) nl2,(vp2(l),vs2(l),den2(l),dep2(l),l=1,nl2)
c  f = ds(i)*g(i) , with ds9i) = 1
c si should be a complex array of length nt and elements (1.,0.)
c*********************************************************************
      subroutine kiksyn(sy,nt,dt,ib,ic,str,dip,rak,z,az,p,g,
     +     tqp,tqs,vmod,nl,nl12)
      PARAMETER (NL0=10)
      real,intent(out):: sy(nt)
      real,intent(in):: vmod(12,nl)
      real,intent(in):: dt,str,dip,rak,z,az,p,g,tqp,tqs
      integer,intent(in):: nl12(2),ib,ic
      integer nls(3)
      complex*8 sqp(nt),sqs(nt),si(nt)
      real so(nt),gr(nt)
      COMMON /STR0/NL ,VP (NL0),VS (NL0),DEN (NL0),DEP (NL0)
      COMMON /STR1/NL1,VP1(NL0),VS1(NL0),DEN1(NL0),DEP1(NL0)
      COMMON /STR2/NL2,VP2(NL0),VS2(NL0),DEN2(NL0),DEP2(NL0)
c     Set velocity model common block
      nl1 = nl12(1)
      nl2 = nl12(2)
      do 100 i=1,nl
         vp(i)   = vmod(1,i)
         vs(i)   = vmod(2,i)
         den(i)  = vmod(3,i)
         dep(i)  = vmod(4,i)
         vp1(i)  = vmod(5,i)
         vs1(i)  = vmod(6,i)
         den1(i) = vmod(7,i)
         dep1(i) = vmod(8,i)
         vp2(i)  = vmod(9,i)
         vs2(i)  = vmod(10,i)
         den2(i) = vmod(11,i)
         dep2(i) = vmod(12,i)
 100  enddo
      df = 1./(dt*nt)
      call qfm(sqp,nt,tqp,df)
      call qfm(sqs,nt,tqs,df)
c         -- compute source time function -------------------------------
c Assum triangle ms=1
      ms = 1
      call stime(so,nt,dt,ms,t1,t2)
c     -- compute body waves -------------------------------------
c     dummy instrument response
      do 200 i=1,nt
 200     si(i) = cmplx(1.,0.)
      call bodyw(gr,nt,dt,ib,ic,str,dip,rak,z,az,p,g,sqp,sqs,si)
c     -- convolve green's functions with the source time --------
      call conv(gr,so,sy,mtg,mtg,mtg)
      return
      end
c*********************************************************************
      real function delay(az,sai,da,dr,dz,v)
*===================================================*
c     -- computes delay times with respect to an origin -----------------
c     -- assuming plane wave --------------------------------------------
c    az   source to receiver azimuth
c    sai  product of horizontal slowness and p or s velocity at source
c    da   azimuth of subfault centroid from that of hypocenter subfault
c    dr   distance of   ''  
c    dz   depth of      ''
c    v    p or s velocity at source depth
*===================================================*

      conv=3.1415927/180.

      dx=dr*sin(da*conv)
      dy=dr*cos(da*conv)
      cai=sqrt(1.-sai*sai)
      zr=cai
      yr=cos(az*conv)*sai
      xr=sin(az*conv)*sai

      delay=(xr*dx+yr*dy+zr*dz)/v

      return
      end
c*********************************************************************
      SUBROUTINE BODYW(X,N,DT,IB,IC,F1,D1,A1,H,AZ,P,FC0,ZQP,ZQS,ZI)
*===================================================*
*  Calculate synthetic waveforms                    *
*    Type of body wave  IB= 1/2/3/4: P/SV/SH/PP     *
*    Component          IC= 1/2/3: UD/NS/EW (IB=1/4)*
*                       IC= 1/2  : UD/HR    (IB=2)  *
*                       IC= any  : SH       (IB=3)  *
*============================ Ver.900715 ===========*
* Modification:                                     *
* 1) An isotropic component of M.T. is added        *
*         for d1 > 360 (degree)        -900531      *
* 2) Source layer # is determined in this subroutine*
*      independently from a reference point -900715 *
*===================================================*
      PARAMETER (ND=4096,NL0=10,RADIUS=6371,PI=3.141593)
      IMPLICIT COMPLEX*8 (Z)
      DIMENSION X(N),Z(ND),ZI(ND),ZQP(ND),ZQS(ND)
     -,  ZR0(ND),ZRPU(ND),ZRPD(ND),ZRSU(ND),ZRSD(ND),ZRPP(ND),ZDM(ND)
      COMMON /STR0/NL ,VP (NL0),VS (NL0),DEN (NL0),DEP (NL0)
      COMMON /STR1/NL1,VP1(NL0),VS1(NL0),DEN1(NL0),DEP1(NL0)
      COMMON /STR2/NL2,VP2(NL0),VS2(NL0),DEN2(NL0),DEP2(NL0)
         DF=1/(DT*N)
         DW=DF*2*PI
         TL=DT*N
* < Source layer # >
       hl=0.
       do 15 l=1,nl-1
	 hl=hl+dep(l)
	 dh=h-hl
15     if(dh.lt.0.) goto 16
16      ll=l
	hl=hl-dep(ll)
* < Radiation pattern >
      CALL RADP(AZ-F1,D1,A1,P,PD,PU,SVD,SVU,SHD,SHU,VP(LL),VS(LL))
* < Structure effects: Near-source & near-reciever >
         DH=H-HL
       IF(IB.EQ.3) GOTO 1
       CALL REFL(ZRPU,ZRPD,ZRSU,ZRSD,DW,N,P,VP,VS,DEN,DEP,NL,LL,IB,TR1)
       IF(IC.EQ.1)
     - CALL CNVR(ZDM,ZR0,DW,N,P,VP1,VS1,DEN1,DEP1,NL1,IB,TR2)
       IF(IC.NE.1)
     - CALL CNVR(ZR0,ZDM,DW,N,P,VP1,VS1,DEN1,DEP1,NL1,IB,TR2)
*   PP-reflector
       CALL REFL(ZRPP,ZDM,ZDM,ZDM
     -              ,DW,N,P,VP2,VS2,DEN2,DEP2,NL2,NL2,IB,TR0)
              GOTO 2
1     CALL REFLSH(ZRSU,ZRSD,DW,N,P,VS,DEN,DEP,NL,LL,TR1)
      CALL CNVRSH(ZR0,DW,N,P,VS1,DEN1,DEP1,NL1,TR2)
2      CONTINUE
* < Delay for depth phase >
          YS=SQRT(1/VS(LL)**2-P**2)
          TS=YS*DH
       IF(IB.EQ.1.OR.IB.EQ.4) THEN
          YP=SQRT(1/VP(LL)**2-P**2)
          TP=YP*DH
       IF(IB.EQ.1)              DELY=TR1+TR2-TP
* an additional delay of 10 sec is put for PP wave:
       IF(IB.EQ.4)              DELY=TR1+TR2-TP-10.
       ELSE
                                DELY=TR1+TR2-TS
       ENDIF
      DO 3 I=1,N/2
           W=DW*(I-1)
*--------------------------------------------------------------
* P or PP wave
      IF(IB.EQ.1.OR.IB.EQ.4) THEN
        FC=FC0/(4*PI*DEN(NL)*VP(NL)**3)
       Z(I) =      ZRPD(I)* PD*EXP(CMPLX(0.,+W*TP))
*   Exclude pP & sP phases outside the time window
      IF(LL.EQ.NL.AND.TP*2.0.GE.TL) GOTO 11
           Z(I) = Z(I)+ZRPU(I)* PU*EXP(CMPLX(0.,-W*TP))
     -                -ZRSD(I)*SVD*EXP(CMPLX(0.,+W*TS))
     -                -ZRSU(I)*SVU*EXP(CMPLX(0.,-W*TS))
11      Z(I) = Z(I)*ZQP(I)
       IF(IB.EQ.1) GOTO 50
*  PP-reflector & additional Q & Hilbert-transform
         Z(I) = Z(I)*ZRPP(I)*ZQP(I)*CMPLX(0.,1.)
* SV-wave
      ELSEIF(IB.EQ.2) THEN
        FC=FC0/(4*PI*DEN(NL)*VS(NL)**3)
        Z(I) =      ZRSD(I)*SVD*EXP(CMPLX(0.,+W*TS))
*   Exclude pS & sS phases outside the time window
      IF(LL.EQ.NL.AND.TS*2.0.GE.TL) GOTO 21
           Z(I) = Z(I)-ZRPU(I)* PU*EXP(CMPLX(0.,-W*TP))
     -                -ZRPD(I)* PD*EXP(CMPLX(0.,+W*TP))
     -                +ZRSU(I)*SVU*EXP(CMPLX(0.,-W*TS))
21      Z(I) = Z(I)*ZQS(I)
* SH-wave
      ELSEIF(IB.EQ.3) THEN
        FC=FC0/(4*PI*DEN(NL)*VS(NL)**2*VS(LL))
        Z(I) =      ZRSD(I)*SHD*EXP(CMPLX(0.,+W*TS))
*   Exclude sS phase outside the time window
      IF(LL.EQ.NL.AND.TS*2.0.GE.TL) GOTO 31
           Z(I) = Z(I)+ZRSU(I)*SHU*EXP(CMPLX(0.,-W*TS))
31      Z(I) = Z(I)*ZQS(I)
      ENDIF
*--------------------------------------------------------------
50      Z(I)=FC*Z(I)*ZR0(I)*ZI(I)*EXP(CMPLX(0.,W*DELY))
      IF(I.EQ.1) GOTO 3
       Z(N+2-I)=CONJG(Z(I))
3     CONTINUE
         Z(N/2+1)=0
      CALL CFFT(Z,N,1)
      DO 5 I=1,N
5     X(I)=Z(I)*DF
      END
c     program test
c     real x(100)
c     read*, dt
c     read*,mode,t1,t2
c     n=10
c     call stime(x,n,dt,mode,t1,t2)
c     a=0.0
c     do i=2,n
c        a=a+0.5*(x(i)+x(i-1))*dt
c     enddo
c     do i=1,n
c        print*,(i-1)*dt,x(i)
c     enddo
c     print*,a
c     end

c*********************************************************************
c*********************************************************************
      SUBROUTINE STIME(X,N,DT,MODE,T1,T2)
*   Source time function
      PARAMETER(PI=3.1415926)
      REAL X(N)
      DO 1 I=1,N
1     X(I)=0
           IT1=T1/DT+1.5
           IT2=T2/DT+1.5
           IT3=(T1+T2)/DT+1.5

* mode=0 : Impulsive time function
      IF(MODE.EQ.0) THEN
           X(1)=1/DT
           RETURN

* mode=1 : Ramp function       if T2 < T1
*          Symmetric trapezoid if T2 >= T1
      ELSEIF(MODE.EQ.1) THEN
           IF(T2.GE.T1) GOTO 23
             DO 21 I=2,IT1
21          X(I)=(I-1)*DT/T1
             DO 22 I=IT1+1,N
22          X(I)=1.0
          RETURN
23    CONTINUE
         if((t1+t2) .lt. dt) then
            x(2)=1./dt
            return
         else
            a=0
            DO 2 I=2,IT1
              X(I)=(I-1)*DT/(T1*T2)
2             a=a+0.5*(x(i)+x(i-1))*dt
            DO 3 I=IT1+1,IT2
              X(I)=1/T2
3             a=a+0.5*(x(i)+x(i-1))*dt
            DO 4 I=IT2+1,IT3
              X(I)=(T1+T2-(I-1)*DT)/(T1*T2)
4             a=a+0.5*(x(i)+x(i-1))*dt
           a=a+0.5*(x(it3+1)+x(it3))*dt
           do I=1,IT3
              x(i)=x(i)/a
           enddo
           RETURN
         endif

* mode=2 : Ramp function with cosine taper if T2 < T1
*          Symmetric trapezoid             if T2 >= T1
      ELSEIF(MODE.EQ.2) THEN
           IF(T2.GE.T1) GOTO 53
             DO 51 I=2,IT1
51         X(I)=(1-COS(PI*(I-1)*DT/T1))*.5
             DO 52 I=IT1+1,N
52          X(I)=1.0
          RETURN
53    CONTINUE
         DO 5 I=2,IT1
5          X(I)=(1-COS(PI*(I-1)*DT/T1))/(2*T2)
         DO 6 I=IT1+1,IT2
6          X(I)=1/T2
         DO 7 I=IT2+1,IT3
7          X(I)=(1-COS(PI*(T1+T2-(I-1)*DT)/T1))/(2*T2)
           RETURN

* mode=3 : Triangle
      ELSEIF(MODE.EQ.3) THEN
         DO 8 I=2,IT1
8          X(I)=2*(I-1)*DT/(T1*T2)
         DO 9 I=IT1+1,IT2
9          X(I)=2*(T2-(I-1)*DT)/(T2-T1)/T2
          RETURN
       ENDIF
        END
c*********************************************************************
c*********************************************************************
      SUBROUTINE QFm(Z,N,TQ,DF)
*  < Q-filter >
C   modification is made for imaginary part
C    following Prof. Kanamori's suggestion  -01/06/29
      COMPLEX Z(N),Z1
C     FN=DF*N/2
      pi=3.141593
      S2=3.*(2.*pi/tq)
      DO 1 I=2,N/2
      F=DF*(I-1)
C     Z1=CMPLX(0.,2*F*TQ)*LOG(CMPLX(0.,F/FN))
      Z1=CMPLX(0.,2*F*TQ)*LOG(CMPLX(0.,F/s2))
      Z(I)=EXP(Z1)
1     Z(N+2-I)=CONJG(Z(I))
      Z(1)=1
      Z(N/2+1)=0
      END
c*********************************************************************
c*********************************************************************
      SUBROUTINE REFL(ZRPU,ZRPD,ZRSU,ZRSD,DW,N,
     -               P,VP,VS,DEN,DEP,NL,L,IB,TR)
*  < Haskell(1953,BSSA;1962;JGR)'s matrix >
*      IB = 1/2 for P/SV
*      P  = Ray parameter: sin(ih)/v
*  --------------------------------------
*   For IB=1/4(incident P)
*     ZRPU = delta(l)' /delta(nl)''
*     ZRPD = delta(l)''/delta(nl)''
*     ZRSU =2omega(l)' /delta(nl)''
*     ZRSD =2omega(l)''/delta(nl)''
*   For IB=2(incident SV)
*     ZRPU = delta(l)' /2omega(nl)''
*     ZRPD = delta(l)''/2omega(nl)''
*     ZRSU = omega(l)' / omega(nl)''
*     ZRSD = omega(l)''/ omega(nl)''
*  --------------------------------------
*   delta,omega are coefficients of <<potentials>>
*  --------------------------------------
*     TR = P-wave travel time from l to nl layer
*  --------------------------------------
*    Critical incidence angle of SV is considered -980827
*  --------------------------------------
      IMPLICIT COMPLEX*8 (Z)
       REAL VP(NL),VS(NL),DEN(NL),DEP(NL)
       DIMENSION ZRPU(N),ZRPD(N),ZRSU(N),ZRSD(N)
     -,  ZEL(4,4),ZE(4,4),ZA(4,4),ZAA(4,4),ZA1(4,4),ZJ(4,4),ZJL(4,4)
           CALL CLEAR(ZEL,4)
           CALL CLEAR(ZE,4)
           CALL CLEAR(ZA,4)
           CALL CLEAR(ZAA,4)
           CALL CLEAR(ZA1,4)
           CALL CLEAR(ZJ,4)
           CALL CLEAR(ZJL,4)
           P2=P**2
           TR=0
       DO 10 M=L,NL-1
10        TR=TR+DEP(M)*SQRT(1/VP(M)**2-P2)
* E-1 matrix for l layers
           GM=2*(P*VS(L))**2
           VP2=1/VP(L)**2
           Y1=SQRT(VP2-P2)
           RA=Y1/P
           VS2=1/VS(L)**2
           Y2=SQRT(VS2-P2)
           RB=Y2/P
          ZEL(1,1)=-2*(VS(L)/VP(L))**2
          ZEL(1,3)=1/(DEN(L)*VP(L)**2)
          ZEL(2,2)=(GM-1)/(VP(L)**2*RA*P2)
          ZEL(2,4)=1/(DEN(L)*VP(L)**2*RA)
          ZEL(3,1)=(GM-1)/(GM*RB)
          ZEL(3,3)=-P2/(DEN(L)*GM*RB)
          ZEL(4,2)=1
          ZEL(4,4)=P2/(DEN(L)*GM)
* E-1 matrix for nl layers
           GM=2*(P*VS(NL))**2
           VP2=1/VP(NL)**2
C980827          Y1=SQRT(VP2-P2)
           if(vp2.ge.p2) ZY1=SQRT(VP2-P2)
           if(vp2.lt.p2) ZY1=cmplx(0.,-SQRT(p2-VP2))
C980827          RA=Y1/P
           ZRA=ZY1/P
           VS2=1/VS(NL)**2
           Y2=SQRT(VS2-P2)
           RB=Y2/P
          ZE(1,1)=-2*(VS(NL)/VP(NL))**2
          ZE(1,3)=1/(DEN(NL)*VP(NL)**2)
C980827         ZE(2,2)=(GM-1)/(VP(NL)**2*RA*P2)
          ZE(2,2)=(GM-1)/(VP(NL)**2*ZRA*P2)
C980827         ZE(2,4)=1/(DEN(NL)*VP(NL)**2*RA)
          ZE(2,4)=1/(DEN(NL)*VP(NL)**2*ZRA)
          ZE(3,1)=(GM-1)/(GM*RB)
          ZE(3,3)=-P2/(DEN(NL)*GM*RB)
          ZE(4,2)=1
          ZE(4,4)=P2/(DEN(NL)*GM)
*
      DO 100 I=1,N/2
         W=(I-1)*DW
         DO 1 J1=1,4
         DO 1 J2=1,4
1        ZAA(J1,J2)=0
         DO 2 J=1,4
2        ZAA(J,J)=1
        IF(L.EQ.1) CALL PROD(ZEL,ZAA,ZJL,4)
       DO 110 M=1,NL-1
           VP2=1/VP(M)**2
           Y1=SQRT(VP2-P2)
           PM=Y1*W*DEP(M)
           CPM=COS(PM)
           ZSPM=CMPLX(0.,SIN(PM))
           RA=Y1/P
           RMC2=DEN(M)/P2
         IF(VS(M).NE.0.) THEN
           VS2=1/VS(M)**2
           Y2=SQRT(VS2-P2)
           GM=2*(P*VS(M))**2
           QM=Y2*W*DEP(M)
           CQM=COS(QM)
           ZSQM=CMPLX(0.,SIN(QM))
           RB=Y2/P
* A-matrix
          ZA(1,1)=GM*CPM-(GM-1)*CQM
          ZA(1,2)=(GM-1)/RA*ZSPM+GM*RB*ZSQM
          ZA(1,3)=-(1/RMC2)*(CPM-CQM)
          ZA(1,4)=(1/RMC2)*(ZSPM/RA+RB*ZSQM)
          ZA(2,1)=-(GM*RA*ZSPM+(GM-1)/RB*ZSQM)
          ZA(2,2)=-(GM-1)*CPM+GM*CQM
          ZA(2,3)=(1/RMC2)*(RA*ZSPM+ZSQM/RB)
          ZA(2,4)=ZA(1,3)
          ZA(3,1)=RMC2*GM*(GM-1)*(CPM-CQM)
          ZA(3,2)=RMC2*( (GM-1)**2/RA*ZSPM+GM**2*RB*ZSQM)
          ZA(3,3)=ZA(2,2)
          ZA(3,4)=ZA(1,2)
          ZA(4,1)=RMC2*(GM**2*RA*ZSPM+(GM-1)**2/RB*ZSQM)
          ZA(4,2)=ZA(3,1)
          ZA(4,3)=ZA(2,1)
          ZA(4,4)=ZA(1,1)
* Water layer
         ELSE
          ZA(1,1)=1
          ZA(1,2)=-ZSPM/RA
          ZA(1,3)=-(1/RMC2)*CPM
          ZA(1,4)=(1/RMC2)*ZSPM/RA
          ZA(2,1)=0
          ZA(2,2)=CPM
          ZA(2,3)=RA/RMC2*ZSPM
          ZA(2,4)=ZA(1,3)
          ZA(3,1)=0
          ZA(3,2)=RMC2/RA*ZSPM
          ZA(3,3)=ZA(2,2)
          ZA(3,4)=ZA(1,2)
          ZA(4,1)=0
          ZA(4,2)=0
          ZA(4,3)=0
          ZA(4,4)=0
       ENDIF
*
       DO 101 J1=1,4
       DO 101 J2=1,4
101       ZA1(J1,J2)=ZAA(J1,J2)
        CALL PROD(ZA,ZA1,ZAA,4)
* J-matrix for l layers
        IF(M.EQ.L-1) CALL PROD(ZEL,ZAA,ZJL,4)
110     CONTINUE
*
* J-matrix for nl layers
         CALL PROD(ZE,ZAA,ZJ,4)
*
         ZJ1=ZJ(4,2)-ZJ(3,2)
         ZJ2=ZJ(3,1)-ZJ(4,1)
         ZJ3=ZJ(2,2)-ZJ(1,2)
         ZJ4=ZJ(1,1)-ZJ(2,1)
            ZDET=ZJ4*ZJ1-ZJ3*ZJ2
* Propagator coefficients
       IF(IB.EQ.1.OR.IB.EQ.4) THEN
         ZRPU(I)=(ZJ1*(ZJL(1,1)+ZJL(2,1))+ZJ2*(ZJL(1,2)+ZJL(2,2)))/ZDET
         ZRPD(I)=(ZJ1*(ZJL(1,1)-ZJL(2,1))+ZJ2*(ZJL(1,2)-ZJL(2,2)))/ZDET
       ZRSU(I)=2*(ZJ1*(ZJL(3,1)+ZJL(4,1))+ZJ2*(ZJL(3,2)+ZJL(4,2)))/ZDET
       ZRSD(I)=2*(ZJ1*(ZJL(4,1)-ZJL(3,1))+ZJ2*(ZJL(4,2)-ZJL(3,2)))/ZDET
       ELSE
       ZRPU(I)=(ZJ3*(ZJL(1,1)+ZJL(2,1))+ZJ4*(ZJL(1,2)+ZJL(2,2)))/ZDET/2
       ZRPD(I)=(ZJ3*(ZJL(1,1)-ZJL(2,1))+ZJ4*(ZJL(1,2)-ZJL(2,2)))/ZDET/2
         ZRSU(I)=(ZJ3*(ZJL(3,1)+ZJL(4,1))+ZJ4*(ZJL(3,2)+ZJL(4,2)))/ZDET
         ZRSD(I)=(ZJ3*(ZJL(4,1)-ZJL(3,1))+ZJ4*(ZJL(4,2)-ZJL(3,2)))/ZDET
       ENDIF
161    IF(I.EQ.1) GOTO 100
         ZRPU(N+2-I)=CONJG(ZRPU(I))
         ZRPD(N+2-I)=CONJG(ZRPD(I))
         ZRSU(N+2-I)=CONJG(ZRSU(I))
         ZRSD(N+2-I)=CONJG(ZRSD(I))
100   CONTINUE
         ZRPU(N/2+1)=0
         ZRPD(N/2+1)=0
         ZRSU(N/2+1)=0
         ZRSD(N/2+1)=0
      END
c*********************************************************************
c*********************************************************************
      SUBROUTINE REFLSH(ZRSU,ZRSD,DW,N,P,VS,DEN,DEP,NL,L,TR)
*  < Haskell(1953,BSSA;1960,JGR) matrix for SH wave >
*      P  = Ray parameter: sin(ih)/v
*  --------------------------------------
*     ZRSU = v(l)' / v(nl)''
*     ZRSD = v(l)''/ v(nl)''
*  v is coefficient of <<displacement>>
*  --------------------------------------
*     TR = S-wave travel time from l to nl layer
*  --------------------------------------
      IMPLICIT COMPLEX*8 (Z)
       REAL VS(NL),DEN(NL),DEP(NL)
       DIMENSION ZRSU(N),ZRSD(N)
     -  ,ZA(2,2),ZAA(2,2),ZA1(2,2),ZAL(2,2)
           CALL CLEAR(ZA,2)
           CALL CLEAR(ZAA,2)
           CALL CLEAR(ZA1,2)
           CALL CLEAR(ZAL,2)
           P2=P**2
           TR=0
       DO 10 M=L,NL-1
          IF(VS(M).EQ.0.) GOTO 10
          TR=TR+DEP(M)*SQRT(1/VS(M)**2-P2)
10     CONTINUE
      DO 100 I=1,N/2
         W=(I-1)*DW
         DO 1 J1=1,2
         DO 1 J2=1,2
1        ZAA(J1,J2)=0
         DO 2 J=1,2
2        ZAA(J,J)=1
        IF(L.NE.1) GOTO 3
         ZAL(1,1)=1.
         ZAL(2,2)=1.
         ZAL(1,2)=0.
         ZAL(2,1)=0.
3      DO 110 M=1,NL-1
        IF(VS(M).NE.0.) THEN
           VS2=1/VS(M)**2
           Y2=SQRT(VS2-P2)
           RGB=DEN(M)*VS(M)**2*Y2/P
           QM=Y2*W*DEP(M)
           CQM=COS(QM)
           ZSQM=CMPLX(0.,SIN(QM))
* A-matrix
          ZA(1,1)=CQM
          ZA(1,2)=1/RGB*ZSQM
          ZA(2,1)=RGB*ZSQM
          ZA(2,2)=CQM
* Water layer
         ELSE
          ZA(1,1)=1
          ZA(1,2)=0
          ZA(2,1)=0
          ZA(2,2)=0
       ENDIF
*
       DO 101 J1=1,2
       DO 101 J2=1,2
101       ZA1(J1,J2)=ZAA(J1,J2)
        CALL PROD(ZA,ZA1,ZAA,2)
* J-matrix for l layers
        IF(M.NE.L-1) GOTO 110
       DO 102 J1=1,2
       DO 102 J2=1,2
102       ZAL(J1,J2)=ZAA(J1,J2)
110     CONTINUE
*
          VS2=1/VS(NL)**2
          Y2=SQRT(VS2-P2)
          RGB=DEN(NL)*VS(NL)**2*Y2/P
          ZDET=ZAA(1,1)+ZAA(2,1)/RGB
          VS2=1/VS(L)**2
          Y2=SQRT(VS2-P2)
          RGBL=DEN(L)*VS(L)**2*Y2/P
* Propagator coefficients
         ZRSU(I)=(ZAL(1,1)-ZAL(2,1)/RGBL)/ZDET
         ZRSD(I)=(ZAL(1,1)+ZAL(2,1)/RGBL)/ZDET
161    IF(I.EQ.1) GOTO 100
         ZRSU(N+2-I)=CONJG(ZRSU(I))
         ZRSD(N+2-I)=CONJG(ZRSD(I))
100   CONTINUE
         ZRSU(N/2+1)=0
         ZRSD(N/2+1)=0
      END
c*********************************************************************
c*********************************************************************
      SUBROUTINE CNVR(ZR1,ZR2,DW,N,P,VP,VS,DEN,DEP,NL,IB,TR)
*  < Displacement-to-displacement conversion near the free surface >
*      P  = Ray parameter: sin(ih)/v
*  --------------------------------------
*   For IB=1/4(incident P)
*      ZR1 = u(0)/uP(nl) (P to Horizontal component)
*      ZR2 =-w(0)/uP(nl) (P to Vertical component(up))
*   For IB=2  (incident SV)
*      ZR1 = u(0)/uS(nl) (SV to Horizontal component)
*      ZR2 =-w(0)/uS(nl) (SV to Vertical component(up))
*  --------------------------------------
*      TR = P-wave travel time from nl-th layer to top
*  --------------------------------------
*      Critical incidence angle of SV is considered  -980827
*  --------------------------------------
      IMPLICIT COMPLEX*8 (Z)
       REAL VP(NL),VS(NL),DEN(NL),DEP(NL)
      DIMENSION ZR1(N),ZR2(N),ZE(4,4),ZA(4,4),ZAA(4,4),ZA1(4,4),ZJ(4,4)
           CALL CLEAR(ZE,4)
           CALL CLEAR(ZA,4)
           CALL CLEAR(ZAA,4)
           CALL CLEAR(ZA1,4)
           CALL CLEAR(ZJ,4)
           P2=P**2
           TR=0
       DO 10 M=1,NL-1
10        TR=TR+DEP(M)*SQRT(1/VP(M)**2-P2)
* E-1 matrix for nl layers
           GM=2*(P*VS(NL))**2
           VP2=1/VP(NL)**2
C980827          Y1=SQRT(VP2-P2)
            if(vp2.ge.p2) ZY1=SQRT(VP2-P2)
            if(vp2.lt.p2) ZY1=cmplx(0.,-SQRT(p2-VP2))
C980827          RA=Y1/P
           ZRA=ZY1/P
           VS2=1/VS(NL)**2
           Y2=SQRT(VS2-P2)
           RB=Y2/P
          ZE(1,1)=-2*(VS(NL)/VP(NL))**2
          ZE(1,3)=1/(DEN(NL)*VP(NL)**2)
C980827         ZE(2,2)=(GM-1)/(VP(NL)**2*RA*P2)
          ZE(2,2)=(GM-1)/(VP(NL)**2*ZRA*P2)
C980827         ZE(2,4)=1/(DEN(NL)*VP(NL)**2*RA)
          ZE(2,4)=1/(DEN(NL)*VP(NL)**2*ZRA)
          ZE(3,1)=(GM-1)/(GM*RB)
          ZE(3,3)=-P2/(DEN(NL)*GM*RB)
          ZE(4,2)=1
          ZE(4,4)=P2/(DEN(NL)*GM)
*
      DO 100 I=1,N/2
         W=(I-1)*DW
         DO 1 J1=1,4
         DO 1 J2=1,4
1        ZAA(J1,J2)=0
         DO 2 J=1,4
2        ZAA(J,J)=1
       DO 110 M=1,NL-1
           VP2=1/VP(M)**2
           Y1=SQRT(VP2-P2)
           PM=Y1*W*DEP(M)
           CPM=COS(PM)
           ZSPM=CMPLX(0.,SIN(PM))
           RA=Y1/P
           RMC2=DEN(M)/P**2
         IF(VS(M).NE.0.) THEN
          VS2=1/VS(M)**2
           Y2=SQRT(VS2-P2)
           GM=2*(P*VS(M))**2
           QM=Y2*W*DEP(M)
           CQM=COS(QM)
           ZSQM=CMPLX(0.,SIN(QM))
           RB=Y2/P
* A-matrix
          ZA(1,1)=GM*CPM-(GM-1)*CQM
          ZA(1,2)=(GM-1)/RA*ZSPM+GM*RB*ZSQM
          ZA(1,3)=-(1/RMC2)*(CPM-CQM)
          ZA(1,4)=(1/RMC2)*(ZSPM/RA+RB*ZSQM)
          ZA(2,1)=-(GM*RA*ZSPM+(GM-1)/RB*ZSQM)
          ZA(2,2)=-(GM-1)*CPM+GM*CQM
          ZA(2,3)=(1/RMC2)*(RA*ZSPM+ZSQM/RB)
          ZA(2,4)=ZA(1,3)
          ZA(3,1)=RMC2*GM*(GM-1)*(CPM-CQM)
          ZA(3,2)=RMC2*( (GM-1)**2/RA*ZSPM+GM**2*RB*ZSQM)
          ZA(3,3)=ZA(2,2)
          ZA(3,4)=ZA(1,2)
          ZA(4,1)=RMC2*(GM**2*RA*ZSPM+(GM-1)**2/RB*ZSQM)
          ZA(4,2)=ZA(3,1)
          ZA(4,3)=ZA(2,1)
          ZA(4,4)=ZA(1,1)
* Water layer
         ELSE
          ZA(1,1)=1
          ZA(1,2)=-ZSPM/RA
          ZA(1,3)=-(1/RMC2)*CPM
          ZA(1,4)=(1/RMC2)*ZSPM/RA
          ZA(2,1)=0
          ZA(2,2)=CPM
          ZA(2,3)=RA/RMC2*ZSPM
          ZA(2,4)=ZA(1,3)
          ZA(3,1)=0
          ZA(3,2)=RMC2/RA*ZSPM
          ZA(3,3)=ZA(2,2)
          ZA(3,4)=ZA(1,2)
          ZA(4,1)=0
          ZA(4,2)=0
          ZA(4,3)=0
          ZA(4,4)=0
       ENDIF
*
       DO 101 J1=1,4
       DO 101 J2=1,4
101       ZA1(J1,J2)=ZAA(J1,J2)
        CALL PROD(ZA,ZA1,ZAA,4)
110     CONTINUE
*
* J-matrix for nl layers
         CALL PROD(ZE,ZAA,ZJ,4)
*
         ZJ1=ZJ(4,2)-ZJ(3,2)
         ZJ2=ZJ(3,1)-ZJ(4,1)
         ZJ3=ZJ(2,2)-ZJ(1,2)
         ZJ4=ZJ(1,1)-ZJ(2,1)
            ZDET=ZJ4*ZJ1-ZJ3*ZJ2
* Conversion at the free surface
         IF(IB.EQ.1.OR.IB.EQ.4) THEN
*   for P-wave
           ZR1(I)=-ZJ1*2/(ZDET*VP(NL)*P)
           ZR2(I)=+ZJ2*2/(ZDET*VP(NL)*P)
         ELSE
*   for S-wave
           ZR1(I)=ZJ3/(ZDET*VS(NL)*P)
           ZR2(I)=ZJ4/(ZDET*VS(NL)*P)
         ENDIF
       IF(I.EQ.1) GOTO 100
         ZR1(N+2-I)=CONJG(ZR1(I))
         ZR2(N+2-I)=CONJG(ZR2(I))
100   CONTINUE
         ZR1(N/2+1)=0
         ZR2(N/2+1)=0
      END
c*********************************************************************
c*********************************************************************
      SUBROUTINE CNVRSH(ZRC,DW,N,P,VS,DEN,DEP,NL,TR)
*  < Displacement-to-displacement conversion near the free surface >
*      P  = Ray parameter: sin(ih)/v
*  --------------------------------------
*      ZRC = v(0)/vSH(nl) (SH to Horizontal component)
*  --------------------------------------
*      TR = S-wave travel time from nl-th layer to top
*  --------------------------------------
      IMPLICIT COMPLEX*8 (Z)
       REAL VS(NL),DEN(NL),DEP(NL)
       DIMENSION ZRC(N),ZA(2,2),ZAA(2,2),ZA1(2,2)
           CALL CLEAR(ZA,2)
           CALL CLEAR(ZAA,2)
           CALL CLEAR(ZA1,2)
           P2=P**2
           TR=0
       DO 10 M=1,NL-1
          IF(VS(M).EQ.0.) GOTO 10
          TR=TR+DEP(M)*SQRT(1/VS(M)**2-P2)
10     CONTINUE
      DO 100 I=1,N/2
         W=(I-1)*DW
         DO 1 J1=1,2
         DO 1 J2=1,2
1        ZAA(J1,J2)=0
         DO 2 J=1,2
2        ZAA(J,J)=1
3      DO 110 M=1,NL-1
           VS2=1/VS(M)**2
           Y2=SQRT(VS2-P2)
           RGB=DEN(M)*VS(M)**2*Y2/P
           QM=Y2*W*DEP(M)
           CQM=COS(QM)
           ZSQM=CMPLX(0.,SIN(QM))
* A-matrix
          ZA(1,1)=CQM
          ZA(1,2)=1/RGB*ZSQM
          ZA(2,1)=RGB*ZSQM
          ZA(2,2)=CQM
*
       DO 101 J1=1,2
       DO 101 J2=1,2
101       ZA1(J1,J2)=ZAA(J1,J2)
        CALL PROD(ZA,ZA1,ZAA,2)
110     CONTINUE
*
          VS2=1/VS(NL)**2
          Y2=SQRT(VS2-P2)
          RGB=DEN(NL)*VS(NL)**2*Y2/P
          ZDET=ZAA(1,1)+ZAA(2,1)/RGB
         ZRC(I)=2/ZDET
161    IF(I.EQ.1) GOTO 100
         ZRC(N+2-I)=CONJG(ZRC(I))
100   CONTINUE
         ZRC(N/2+1)=0
      END
c*********************************************************************
c*********************************************************************
      SUBROUTINE RADP(AZ,D0,A0,P,PD,PU,SVD,SVU,SHD,SHU,VP,VS)
*  < Radiation pattern >
*  -----------------------------------
*        PD,PU for P,p
*        SD,SU for S,s(SV&SH)
*     Polarity
*        SH : Clockwise
*        SV : SH x P
*  -----------------------------------
*  May 29,1990
*    adding radiation pattern for isotropic moment-tensor
*        for D0 > 360
*  -----------------------------------
      IF(D0.GT.360.) THEN
         PD=1
         PU=1
         SVD=0
         SVU=0
         SHD=0
         SHU=0
      RETURN
      ENDIF
      RAD=.0174533
      DR0=D0*RAD
      AR0=A0*RAD
      THET=AZ*RAD
      SIH=P*VP
      IF(SIH.GE.1.) CIH=0.
      IF(SIH.LT.1.) CIH=SQRT(1-SIH**2)
      SIH2=P*VS
      CIH2=SQRT(1-SIH2**2)
C     P2=P*P
C   *** C PARAMETERS
      C1=SIH**2
      C2=-2.*SIH*CIH
      C3=2-3*SIH**2
C    *** SV PARAMETERS
      SV1=SIH2*CIH2
      SV2=2*SIH2**2-1
      SV3=-3*SIH2*CIH2
C    *** SH PARAMETERS
      SH1=SIH2
      SH2=CIH2
C    ***** A PARAMETERS
      THE2 =THET*2.
      DR02=DR0*2.
      STH=SIN(THET)
      CTH=COS(THET)
      STH2=SIN(THE2 )
      CTH2=COS(THE2 )
      SD0=SIN(DR0)
      CD0=COS(DR0)
      SD02=SIN(DR02)
      CD02=COS(DR02)
      CA0=COS(AR0)
      SA0=SIN(AR0)
      A1=STH2*CA0*SD0+.5*CTH2*SA0*SD02
      A2=CTH*CA0*CD0-STH*SA0*CD02
      A3=.5*SA0*SD02
      A4=CTH2*CA0*SD0-.5*STH2*SA0*SD02
      A5=STH*CA0*CD0+CTH*SA0*CD02
*
      PD =A1*C1+A2*C2+A3*C3
      PU =A1*C1-A2*C2+A3*C3
      SVD =A1*SV1+A2*SV2+A3*SV3
      SVU =-A1*SV1+A2*SV2-A3*SV3
      SHD =A4*SH1+A5*SH2
      SHU =A4*SH1-A5*SH2
      END
c*********************************************************************
c*********************************************************************
      SUBROUTINE CONV(X,Y,Z,L,M,N)
**********************************
*       Z = X * Y                *
**********************************
      DIMENSION X(L),Y(M),Z(N)
      DO 1 I=1,N
      Z(I)=0.
      I0=MIN0(L,I)
      DO 1 J=1,I0
      I1=I-J+1
      IF(I1.LE.0) GOTO 1
      IF(I1.GT.M) GOTO 1
      Z(I)=Z(I)+X(J)*Y(I1)
    1 CONTINUE
      END
c*********************************************************************
c*********************************************************************
        SUBROUTINE PROD(ZA,ZB,ZC,N)
* Product of matrix
        COMPLEX*8 ZA(N,N),ZB(N,N),ZC(N,N)
       DO 1 J1=1,N
       DO 1 J2=1,N
        ZC(J1,J2)=0
         DO 1 J=1,N
1       ZC(J1,J2)=ZC(J1,J2)+ZA(J1,J)*ZB(J,J2)
       END
c*********************************************************************
c*********************************************************************
      SUBROUTINE CFFT(X,N,ID)
* <  FFT for complex variables >
      IMPLICIT COMPLEX*8 (X-Z)
      DIMENSION X(N)
       PI=SIGN(3.141593,ID*1.)
        N2=N/2
        L=N
1      NI=L/2
        Z=CEXP(CMPLX(0.,PI/NI))
        ZK=1
      DO 3 K=1,NI
       DO 2 M=1,N/L
          J=K+2*(M-1)*NI
          XJ=X(J)
          X(J)=XJ+X(J+NI)
2         X(J+NI)=(XJ-X(J+NI))*ZK
3      ZK=ZK*Z
        L=NI
      IF(L.GT.1) GOTO 1
      DO 10 I=1,N
         J=I-1
         M=N2
         L=1
5       L=L+M*(J-2*(J/2))
         M=M/2
         J=J/2
       IF(J.GE.1) GOTO 5
       IF(I.GT.L) GOTO 10
         XI=X(I)
         X(I)=X(L)
         X(L)=XI
10    CONTINUE
      END
c*********************************************************************
c*********************************************************************
       SUBROUTINE CLEAR(Z,N)
       COMPLEX Z(N,N)
       DO 1 I=1,N
       DO 1 J=1,N
1      Z(I,J)=0
       END
c*********************************************************************
c*********************************************************************
      SUBROUTINE raypgeom(GFACT,PP0,H,VP00,ir0,ip)
*  <<  DELTA VERSUS IH CURVE  >>
      IMPLICIT REAL*8 (A-H)
      IMPLICIT REAL*8 (O-Z)
      parameter(dh1=20.0, ndep1=145, dh=4.0, nddh=5)
      DIMENSION DELTA(500),V(1000),Q(200),R(1000)
      DIMENSION VV(200),TIH(500),vs(200)
      real*8,intent(out):: VP00,GFACT(115),PP0(115)
      integer,intent(in):: ir0,ip
      real*8,intent(in):: H


      data vv/5.80000, 6.50000, 8.04118, 8.04353, 8.04588, 8.04823, 
     +8.06389, 8.11945, 8.17500, 8.23056, 8.28612, 8.35475, 
     +8.42775, 8.50075, 8.57375, 8.64675, 8.71975, 8.79275, 
     +8.86575, 8.93875, 9.01175, 9.41040, 9.47760, 9.54480, 
     +9.61200, 9.67920, 9.74640, 9.81360, 9.88080, 9.94800, 
     +10.01520, 10.08240, 10.14960, 10.80329, 10.85645, 10.90961, 
     +10.96277, 11.01593, 11.06484, 11.10063, 11.13607, 11.17116, 
     +11.20591, 11.24032, 11.27439, 11.30813, 11.34154, 11.37464, 
     +11.40742, 11.43989, 11.47205, 11.50392, 11.53549, 11.56676, 
     +11.59776, 11.62847, 11.65891, 11.68907, 11.71897, 11.74861, 
     +11.77800, 11.80713, 11.83601, 11.86466, 11.89307, 11.92125, 
     +11.94920, 11.97693, 12.00444, 12.03174, 12.05884, 12.08573, 
     +12.11243, 12.13894, 12.16526, 12.19139, 12.21736, 12.24314, 
     +12.26876, 12.29422, 12.31953, 12.34468, 12.36968, 12.39454, 
     +12.41926, 12.44385, 12.46831, 12.49266, 12.51688, 12.54099, 
     +12.56499, 12.58889, 12.61269, 12.63640, 12.66002, 12.68356, 
     +12.70702, 12.73041, 12.75373, 12.77698, 12.80018, 12.82332, 
     +12.84641, 12.86947, 12.89248, 12.91546, 12.93841, 12.96134, 
     +12.98424, 13.00714, 13.03002, 13.05290, 13.07578, 13.09867, 
     +13.12157, 13.14448, 13.16742, 13.19038, 13.21337, 13.23640, 
     +13.25946, 13.28257, 13.30574, 13.32896, 13.35223, 13.37558, 
     +13.39899, 13.42248, 13.44604, 13.46970, 13.49344, 13.51728, 
     +13.54121, 13.56526, 13.58941, 13.61367, 13.63806, 13.65756, 
     +13.66217, 13.66679, 13.67141, 13.67603, 13.68064, 13.68526, 
     +13.68988,  55*0.0/
      data vs/3.36000, 3.75000, 4.47353, 4.48059, 4.48765, 4.49471, 
     + 4.50100, 4.50500, 4.50900, 4.51300, 4.51700, 4.54810, 
     + 4.58290, 4.61770, 4.65250, 4.68730, 4.72210, 4.75690, 
     + 4.79170, 4.82650, 4.86130, 5.10180, 5.14420, 5.18660, 
     + 5.22900, 5.27140, 5.31380, 5.35620, 5.39860, 5.44100, 
     + 5.48340, 5.52580, 5.56820, 5.96297, 6.01487, 6.06676, 
     + 6.11866, 6.17056, 6.21337, 6.22876, 6.24397, 6.25900, 
     + 6.27385, 6.28853, 6.30304, 6.31738, 6.33155, 6.34556, 
     + 6.35941, 6.37310, 6.38664, 6.40002, 6.41326, 6.42634, 
     + 6.43928, 6.45208, 6.46474, 6.47726, 6.48965, 6.50191, 
     + 6.51404, 6.52604, 6.53791, 6.54967, 6.56130, 6.57282, 
     + 6.58423, 6.59553, 6.60672, 6.61780, 6.62878, 6.63966, 
     + 6.65044, 6.66112, 6.67172, 6.68222, 6.69264, 6.70297, 
     + 6.71322, 6.72339, 6.73348, 6.74350, 6.75345, 6.76333, 
     + 6.77314, 6.78289, 6.79258, 6.80221, 6.81179, 6.82131, 
     + 6.83078, 6.84020, 6.84958, 6.85891, 6.86820, 6.87746, 
     + 6.88668, 6.89587, 6.90503, 6.91416, 6.92327, 6.93236, 
     + 6.94142, 6.95047, 6.95951, 6.96854, 6.97755, 6.98656, 
     + 6.99557, 7.00458, 7.01358, 7.02260, 7.03162, 7.04064, 
     + 7.04969, 7.05874, 7.06782, 7.07691, 7.08603, 7.09517, 
     + 7.10434, 7.11354, 7.12277, 7.13204, 7.14135, 7.15070, 
     + 7.16009, 7.16953, 7.17902, 7.18856, 7.19815, 7.20780, 
     + 7.21751, 7.22728, 7.23712, 7.24702, 7.25700, 7.26574, 
     + 7.27071, 7.27568, 7.28064, 7.28561, 7.29058, 7.29554, 
     + 7.30051, 55*0.0/
      data q/600.0, 600.0, 400.0, 80.0, 90., 100., 120., 150.,
     +200., 200., 200., 200., 200., 200., 200., 200.,
     +200., 200., 200., 200., 200., 200., 200., 200.,
     +200., 200., 200., 200., 200., 200., 250., 250.,
     +250., 250., 250., 250., 250., 250., 250., 250.,
     +300., 300., 300., 300., 300., 3000., 3000., 3000.,
     +3000., 3000., 6000., 6000., 6000., 6000., 6000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +10000., 10000., 10000., 10000., 10000., 10000., 10000., 10000.,
     +2000., 2000., 2000., 2000., 2000., 200., 200., 200.,
     +200., 200., 200., 200., 200., 200., 200., 200.,
     +100., 100., 100., 100., 100., 100., 100., 100.,
     +60., 60., 54*0.0/

      NDEP=NDDH*(NDEP1-1)
      DO 11 I=2,NDEP1
         DO 11 J=1,NDDH
            IJ=NDDH*(I-2)+J
	    if(ip.eq.1) then
              V(IJ)=VV(I-1)+(VV(I)-VV(I-1))/DH1*DH*FLOAT(J)
	    else
              V(IJ)=Vs(I-1)+(Vs(I)-Vs(I-1))/DH1*DH*FLOAT(J)
	    endif
   11 CONTINUE
C
      RAD=3.1415926/180.
      L=H/DH+1.
      VP00=V(L)
      ir=1

*  Moho-to-Moho ray-path
      if(ir0.eq.1) then
         L=11
	 ir=L
      endif

      R0=6370.
      RH=R0-H
      DO 9 I=1,NDEP
    9    R(I)=R0-DH*FLOAT(I)
C
      DO 100 J=1,180
         TH=RAD*FLOAT(J)*0.5
         P=RH*SIN(TH)/V(L)
         PA=P*V(NDEP)/R(NDEP)
         DEL=0.
         IF(PA.GT.1.) THEN
            R1=RH
            R2=R(L)
            I=L
   10       PA=P*V(I)/R1
            IF(PA.GE.1.) GOTO 20
            DEL1=DASIN(PA)
            PA=P*V(I)/R2
            IF(PA.GE.1.) THEN
               DEL2=3.141593/2.
               DEL=DEL+DEL2-DEL1
               GOTO 20
            ELSE
               DEL2=DASIN(PA)
               DEL=DEL+DEL2-DEL1
               R1=R(I)
               I=I+1
               R2=R(I)
               GOTO 10
            END IF
   20       DEL=DEL*2.
            R2=RH
            I=L-1
            IF(I.EQ.ir-1) GOTO 30
            R1=R(I)
   25       DEL1=DASIN(P*V(I)/R1)
            DEL2=DASIN(P*V(I)/R2)
            DEL=DEL+DEL2-DEL1
            R2=R1
            I=I-1
            IF(I.EQ.ir-1) GOTO 30
            R1=R(I)
            GOTO 25
   30       DEL=DEL+DASIN(P*V(ir)/R0)-DASIN(P*V(ir)/R2)
         else
   90       DEL=3.141593
	 endif
  100 DELTA(J)=DEL
C
      DO 200 J= 5,125
         DEL=RAD*FLOAT(J)
         I=1
  110    CONTINUE
         IF(DEL.GT.DELTA(I)) THEN
            M=I-1
            GOTO 120
         ELSE
            I=I+1
            GOTO 110
         END IF
  120    TIH(J)=M+(DELTA(M)-DEL)/(DELTA(M)-DELTA(M+1))
         TIH(J)=TIH(J)*0.5
  200 CONTINUE
      DO 300 J=20,115
C        DTIH=(2.*TIH(J-4)+TIH(J-2)-TIH(J+2)-2.*TIH(J+4))/20.
         DTIH=0.
         DO 301 KK=-15,15
  301    DTIH=DTIH+KK*TIH(J-KK)
         DTIH=DTIH/2480
         SH=SIN(TIH(J)*RAD)
         PP0(J)=SH/V(L)
         PP=V(ir)*PP0(J)
         CS0=SQRT(1.-PP*PP)
         SDEL=SIN(FLOAT(J)*RAD)
         GG=V(L)*SH*DTIH/V(ir)/SDEL/CS0
C        G=SQRT(GG)/R0*1.E5
  300    GFACT(J)=SQRT(GG)/R0*1.E5
      END
