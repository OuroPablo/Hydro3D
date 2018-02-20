!=======================================================================
!			Pablo Ouro Barba
!			Cardiff 2013-2014
!=======================================================================
!######################################################################
       Double precision function spline_33(r)	!Gil JCP 2013
!######################################################################
      implicit none 
       Double precision, intent(in) :: r

      IF (r .le. -3.0D0) THEN
        spline_33 = 0.D0
      ELSE IF ((r.gt. -3.0D0).and.(r.le. -2.0D0)) THEN
        spline_33 = 
     &   -29.d0/7560.d0*r**7 -5.d0/72.d0*r**6 -21.d0/40.d0*r**5
     &   -17.d0/8.d0*r**4    -39.d0/8.d0*r**3 -243.d0/40.d0*r**2 
     &   -27.d0/8.d0*r       -81.d0/280.d0
      ELSE IF ((r.gt. -2.0D0).and.(r.le. -1.0D0)) THEN
        spline_33 = 
     &    17.d0/1512.d0*r**7 +1.d0/8.d0*r**6   +13.d0/24.d0*r**5
     &    +79.d0/72.d0*r**4  +65.d0/72.d0*r**3 +7.d0/120.d0*r**2 
     &   + 13.d0/72.d0*r     +1447.d0/2520.d0
      ELSE IF ((r.gt. -1.0D0).and.(r.le. 0.0D0)) THEN
        spline_33 = 
     &   -11.d0/756.d0*r**7  -1.d0/18.d0*r**6 
     &   +7.d0/36.d0*r**4    -29.d0/60.d0*r**2 
     &   +691.d0/1260.d0
      ELSE IF ((r.gt. 0.0D0).and.(r.le. 1.0D0)) THEN
        spline_33 = 
     &   11.d0/756.d0*r**7  -1.d0/18.d0*r**6 
     &   +7.d0/36.d0*r**4   -29.d0/60.d0*r**2 
     &   +691.d0/1260.d0
      ELSE IF ((r.gt. 1.0D0).and.(r.le. 2.0D0)) THEN
        spline_33 = 
     &    -17.d0/1512.d0*r**7 +1.d0/8.d0*r**6   -13.d0/24.d0*r**5
     &    +79.d0/72.d0*r**4   -65.d0/72.d0*r**3 +7.d0/120.d0*r**2 
     &    -13.d0/72.d0*r      +1447.d0/2520.d0
      ELSE IF ((r.gt. 2.0D0).and.(r.le. 3.0D0)) THEN
        spline_33 = 
     &    29.d0/7560.d0*r**7 -5.d0/72.d0*r**6 +21.d0/40.d0*r**5
     &   -17.d0/8.d0*r**4    +39.d0/8.d0*r**3 -243.d0/40.d0*r**2 
     &   +27.d0/8.d0*r       -81.d0/280.d0
	ELSE IF (r .gt. 3.0D0) THEN
	  spline_33 = 0.D0
      END IF
      RETURN
      End Function
!######################################################################
       Double precision function phi_r1smth(r)
!######################################################################
      implicit none 
       Double precision, intent(in) :: r
       Double precision :: PI,abr
       PI = 4.D0*DATAN(1.D0)
       abr=SQRT(r*r)
      IF (abr.ge.1.5) THEN
        phi_r1smth = 0.0
      ELSE IF ((abr.lt.1.5).and.(abr.ge.0.5)) THEN
        phi_r1smth = 9./8.-3.*abr/2+abr**2/2
      ELSE IF ((abr.lt.0.5).and.(abr.ge.0.0)) THEN
        phi_r1smth = 3./4.-abr**2
      END IF
      RETURN
      End Function
!######################################################################
      Double precision function phi_r2smth(r)
!######################################################################
      implicit none 
       Double precision, intent(in) :: r
       Double precision :: PI
       PI = 4.D0*DATAN(1.D0)
      IF (r.le.-2.5) THEN
        phi_r2smth = 0.0
      ELSE IF ((r.ge.-2.5).and.(r.le.-1.5)) THEN
        phi_r2smth= -1./8./PI*(-5.*PI-2.*PI*r+4.*sin(PI/4.*(-2.*r-1.)))
      ELSE IF ((r.ge.-1.5).and.(r.le.0.0)) THEN
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(-2.*r+1.))
     &                           -2.*sin(PI/4.*(-2.*r-1.)))
      ELSE IF ((r.ge.0.0).and.(r.le.1.5)) THEN
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(2.*r+1.))
     &                           -2.*sin(PI/4.*(2.*r-1.)))
      ELSE IF ((r.ge.1.5).and.(r.le.2.5)) THEN
        phi_r2smth= -1./8./PI*(-5.*PI+2.*PI*r+4.*sin(PI/4.*(2.*r-1.)))
      ELSE IF (r.ge.2.5) THEN
        phi_r2smth = 0.0
      END IF

      RETURN
      End Function
!######################################################################
       Double precision function phi_r3(r)
!######################################################################
      implicit none 
       Double precision, intent(in) :: r
      IF (r.le.-1.5) THEN
        phi_r3 = 0.0
      ELSE IF ((r.ge.-1.5).and.(r.le.-0.5)) THEN
        phi_r3 = 1.0/6.0*(5.0+3.0*r-sqrt(-3.0*(1.0+r)**2+1.0))
      ELSE IF ((r.ge.-0.5).and.(r.le.0.0)) THEN
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      ELSE IF ((r.ge.0.0).and.(r.le.0.5)) THEN
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      ELSE IF ((r.ge.0.5).and.(r.le.1.5)) THEN
        phi_r3 = 1.0/6.0*(5.0-3.0*r-sqrt(-3.0*(1.0-r)**2+1.0))
      ELSE IF (r.ge.1.5) THEN
        phi_r3 = 0.0
      END IF

      RETURN
      End Function
!######################################################################
       Double precision function phi_r3smth2(r)
!######################################################################
      implicit none 
       Double precision, intent(in) :: r
       Double precision :: PI
       PI = 4.D0*DATAN(1.D0)
      IF (r.le.-2.0) THEN
        phi_r3smth2 = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r3smth2 = 55.0/48.0 - sqrt(3.0)*pi/108.0 + 13.0*r/12.0
     &+ r**2/4.0 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
     &+ sqrt(3.0)/36.0*ASIN(sqrt(3.0)/2.0*(-2.0*r-3.0))
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r3smth2 = 17.0/48.0 + sqrt(3.0)*pi/108.0 - r/4.0
     &- r**2/4.0 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
     &- sqrt(3.0)/12.0*ASIN(sqrt(3.0)/2.0*(-2.0*r-1.0))
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r3smth2 = 17.0/48.0 + sqrt(3.0)*pi/108.0 + r/4.0
     &- r**2/4.0 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
     &- sqrt(3.0)/12.0*ASIN(sqrt(3.0)/2.0*(2.0*r-1.0))
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r3smth2 = 55.0/48.0 - sqrt(3.0)*pi/108.0 - 13.0*r/12.0
     &+ r**2/4.0 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
     &+ sqrt(3.0)/36.0*ASIN(sqrt(3.0)/2.0*(2.0*r-3.0))
      ELSE IF (r.ge.2.0) THEN
        phi_r3smth2 = 0.0
      END IF

      RETURN
      End Function
!######################################################################
       Double precision function phi_r3smth(r)
!######################################################################
      implicit none 
      Double precision, intent(in) :: r
      Double precision :: scal1,scal2,scal3,scal4,scal5,scal6

      scal1 = 1.095450017660160615261  ! 55.0/48.0 - sqrt(3.0)*pi/108.0
      scal2 = 1.083333333333333333333  ! 13.0/12.0
      scal3 = 0.4045499823398393847387 ! 17.0/48.0 + sqrt(3.0)*pi/108.0
      scal4 = 0.0481125224324688137091 ! sqrt(3.0)/36.0
      scal5 = 0.1443375672974064411273 ! sqrt(3.0)/12.0
      scal6 = 0.8660254037844386467637 ! sqrt(3.0)/2.0

      IF (r.le.-2.d0) THEN
        phi_r3smth = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r3smth = scal1 + scal2*r + scal4*ASIN(scal6*(-2.0*r-3.0))
     &+ 0.25*r**2 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r3smth = scal3 - r/4.0 - scal5*ASIN(scal6*(-2.0*r-1.0))
     &- 0.25*r**2 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r3smth = scal3 + r/4.0 - scal5*ASIN(scal6*(2.0*r-1.0))
     &- 0.25*r**2 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r3smth = scal1 - scal2*r + scal4*ASIN(scal6*(2.0*r-3.0))
     &+ 0.25*r**2 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
      ELSE IF (r.ge.2.0) THEN
        phi_r3smth = 0.0
      END IF

      RETURN
      End Function
!######################################################################
       Double precision function phi_r4(r)
!######################################################################
      implicit none 
       Double precision, intent(in) :: r
      IF (r.le.-2.0) THEN
        phi_r4 = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r4 = 1.0/8.0*(5.0+2.0*r-sqrt(-7.0-12.0*r-4.0*r**2))
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r4 = 1.0/8.0*(3.0+2.0*r+sqrt(1.0-4.0*r-4.0*r**2))
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r4 = 1.0/8.0*(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*r**2))
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r4 = 1.0/8.0*(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*r**2))
      ELSE IF (r.ge.2.0) THEN
        phi_r4 = 0.0
      END IF

      RETURN
      End Function
!######################################################################
       Double precision function phi_r4smth(r)
!######################################################################
      implicit none
       Double precision, intent(in) :: r
       Double precision :: PI,ar
       PI = 4.D0*DATAN(1.D0)
	ar=abs(r)
      IF ((ar.ge.0.0).and.(ar.le.0.5)) THEN
	phi_r4smth = 3.0/8.0+ PI/32.0 - r*r/4.0
      ELSE IF ((ar.ge.0.5).and.(ar.le.1.5)) THEN
	phi_r4smth = 1.0/4.0 + (1.0-ar)/8.0*SQRT(-2.0+8.0*ar -4.0*r*r)
     & -1.0/8.0*ASIN(sqrt(2.0)*(ar-1.0))
      ELSE IF ((ar.ge.1.5).and.(ar.le.2.5)) THEN
	phi_r4smth = 17.0/16.0 - PI/64.0 - 3.0*ar/4.0+ r*r/8.0
     & + (ar-2.0)/16.0 * SQRT(-14.0+ 16.0*ar - 4.0*r*r)
     & + 1.0/16.0*ASIN(sqrt(2.0)*(ar-2.0))
      ELSE IF (ar.ge.2.5) THEN
        phi_r4smth = 0.0
      END IF
      RETURN
      End Function
!######################################################################
       Double precision function dh(dx,dy,dz,xij,yij,zij,Xl,Yl,Zl,order)
!######################################################################
      implicit none 
       Double precision,    intent(in) :: dx,dy,dz,xij,yij,zij,Xl,Yl,Zl
       INTEGER, intent(in) :: order
       Double precision :: phi_r2smth,phi_r3smth,phi_r4
       Double precision :: phi_r3,phi_r1smth,phi_r4smth,spline_33

      SELECT CASE (order)

         CASE (1)
      dh =  phi_r1smth((xij-Xl)/dx) 
     & * phi_r1smth((yij-Yl)/dy) * phi_r1smth((zij-Zl)/dz)
         CASE (2)
      dh =   phi_r2smth((xij-Xl)/dx) 
     & * phi_r2smth((yij-Yl)/dy) * phi_r2smth((zij-Zl)/dz)
         CASE (3)
       dh =   phi_r3smth((xij-Xl)/dx) 
     &   * phi_r3smth((yij-Yl)/dy) * phi_r3smth((zij-Zl)/dz)
         CASE (4)
       dh =   phi_r4smth((xij-Xl)/dx) 
     &   * phi_r4smth((yij-Yl)/dy) * phi_r4smth((zij-Zl)/dz)
         CASE (5)
       dh =   phi_r3((xij-Xl)/dx) 
     &   * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
         CASE (6)
      dh =   phi_r4((xij-Xl)/dx)
     &   * phi_r4((yij-Yl)/dy) * phi_r4((zij-Zl)/dz)
         CASE (7)
      dh =   spline_33((xij-Xl)/dx)
     &   * spline_33((yij-Yl)/dy) * spline_33((zij-Zl)/dz)
      Case default

         print*, '===ERROR==='
         print*, ' order of delta function is not selected '
          STOP

      End Select

      RETURN
      End Function
!######################################################################
!######################################################################
!		MOVING LEAST-SQUARE methodology
! Implemented by Pablo (2015)
! Algorithm from Ramirez and Nogueira (University of A Coruna)
!######################################################################
!######################################################################
	subroutine ShapeFunction_MLS(UorV,L,ib)
!######################################################################
         use vars
	 use multidata
         use mpi
	 use imb
         implicit none
         INTEGER :: order_mls,nbase,MAXCELL,MAXNEIG
	 LOGICAL :: vanella_mls
         INTEGER, intent(in) :: UorV,L,ib
         INTEGER :: I,J,K,ipface,inodehalo,nn,neignum
         INTEGER:: MVefx(126),MVefy(126),MVefz(126)
         Double precision :: vmaxdist,xN,yN,zN,h0,xv2,yv2,zv2
    	 Double precision::  Xf(126),Yf(126),Zf(126),hx,hy,hz
         Double precision :: offsety,offsetx,offsetz,ratioHKx1,ratioHKx2
         Double precision :: ratioHKy1,ratioHKy2,ratioHKz1,ratioHKz2
         Double precision :: kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,dist_ff
	 Double precision :: parametros(9),W(126),difx(126)
	 Double precision :: dify(126),difz(126)
         Double precision :: dX2(126),dY2(126),dZ2(126)
         Double precision :: FF(126),dFFX(126),dFFY(126),dFFZ(126)
         Double precision :: ddFFXX(126),ddFFXY(126),ddFFYY(126)
         Double precision :: ddFFXZ(126),ddFFZZ(126),ddFFYZ(126)
	 Double precision :: nmls,k_mls
         INTEGER :: MLS_neignum(maxnodeIBS)

!"""  Data for MLS

	MAXCELL=10000 ; MAXNEIG=65 
	MLS_neignum=MAXNEIG
	order_mls=1
	k_mls=5.0d+00

        if (order_mls .eq. 1) nbase=4
        if (order_mls .eq. 2) nbase=10
        if (order_mls .eq. 3) nbase=20
!""""

	Do i=1,MAXNEIG
 	  MVefx(i)=0 ; MVefy(i)=0 ; MVefz(i)=0 
	Enddo
	nn=0  ; nxl=1.9999999

	IF (UorV.eq.1) then
       DO I = 1, dom(ib)%ttc_i 
       IF ( dom(ib)%x(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &      dom(ib)%x(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx))  GOTO 100
        DO J = 1, dom(ib)%ttc_j
       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy))  GOTO 101
          DO K = 1, dom(ib)%ttc_k
       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 102
		nn=nn+1
	IF (nn.gt.MAXNEIG) write(6,*)'More Neighbours than MAXNEIG!',nn
	    MVefx(nn)=I   ;  MVefy(nn)=J	;  MVefz(nn)=K
           Xf(nn)=dom(ib)%x(i)  ; Yf(nn)=dom(ib)%yc(j)
	    Zf(nn)=dom(ib)%zc(k)
102 	CONTINUE
         END DO
101 	CONTINUE
        END DO
100 	CONTINUE
       END DO
	ENDIF
	IF (UorV.eq.2) then
       DO I = 1, dom(ib)%ttc_i
       IF ( dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &      dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx) )  GOTO 200 
        DO J = 1, dom(ib)%ttc_j
       IF (dom(ib)%y(j) .gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &     dom(ib)%y(j) .lt.(nodey_loc(L)-nxl*dom(ib)%dy) )  GOTO 201
          DO K = 1, dom(ib)%ttc_k
	IF( dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz) )  GOTO 202
		nn=nn+1
	IF (nn.gt.MAXNEIG)
     &       write(6,*)'More Neighbours than MAXNEIG!',nn,nxl
		MVefx(nn)=I   ;  MVefy(nn)=J	;  MVefz(nn)=K
              Xf(nn)=dom(ib)%xc(i)  
	       Yf(nn)=dom(ib)%y(j)
	       Zf(nn)=dom(ib)%zc(k)
202 	CONTINUE
         END DO
201 	CONTINUE
        END DO
200 	CONTINUE
       END DO
	ENDIF
	IF (UorV.eq.3) then
       DO I = 1, dom(ib)%ttc_i 
       IF ( dom(ib)%xc(i).gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
     &      dom(ib)%xc(i).lt.(nodex_loc(L)-nxl*dom(ib)%dx)) GOTO 300
        DO J = 1, dom(ib)%ttc_j
       IF ( dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
     &      dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy)) GOTO 301
          DO K = 1, dom(ib)%ttc_k
       IF ( dom(ib)%z(k) .gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
     &      dom(ib)%z(k) .lt.(nodez_loc(L)-nxl*dom(ib)%dz) ) GOTO 302
		nn=nn+1
	IF (nn.gt.MAXNEIG)
     &       write(6,*)'More Neighbours than MAXNEIG!',nn,nxl
		MVefx(nn)=I   ;  MVefy(nn)=J	;  MVefz(nn)=K
              Xf(nn)=dom(ib)%xc(i)  
	       Yf(nn)=dom(ib)%yc(j)
	       Zf(nn)=dom(ib)%z(k)
302 	CONTINUE
         END DO
301 	CONTINUE
        END DO
300 	CONTINUE
       END DO
	ENDIF

	neignum=nn
	MLS_neignum(L)=neignum	
!--------------------------------------------------------
        vmaxdist= 0.d+00   
    	 do j=1, neignum
	  difx(j)= 0.d+0 ;    dify(j)= 0.d+0;    difz(j)= 0.d+00
	 enddo

    	  do j=1, neignum
	 dist_ff=DSQRT((nodex_loc(L)-Xf(j))**2+(nodey_loc(L)-Yf(j))**2
     &		+(nodez_loc(L)-Zf(j))**2)
	    if(dist_ff.gt.vmaxdist) vmaxdist= dist_ff
	     	difx(j)=Xf(j)-nodex_loc(L)
		dify(j)=Yf(j)-nodey_loc(L)
		difz(j)=Zf(j)-nodez_loc(L)
	  enddo

	h0=1.0d+00
	hx= h0*maxval(DABS(difx(1:neignum)))
	hy= h0*maxval(DABS(dify(1:neignum)))
	hz= h0*maxval(DABS(difz(1:neignum)))

	parametros(1)= k_mls
	parametros(2)= k_mls
	parametros(3)= k_mls
	parametros(4)= k_mls
	parametros(5)= k_mls
	parametros(6)= k_mls
	parametros(7)= 0.0d+00
	parametros(8)= 0.0d+00
	parametros(9)= 0.0d+00

	ratioHKx1= parametros(1)
	ratioHKx2= parametros(2)
	ratioHKy1= parametros(3)
	ratioHKy2= parametros(4)
	ratioHKz1= parametros(5)
	ratioHKz2= parametros(6)
	offsetx  = parametros(7)
	offsety  = parametros(8)
	offsetz  = parametros(9)

	h0= 1.0d+00
	kwx1= h0*ratioHKx1
	kwx2= h0*ratioHKx2
	kwy1= h0*ratioHKy1
	kwy2= h0*ratioHKy2
	kwz1= h0*ratioHKz1
	kwz2= h0*ratioHKz2

	call  kerneltododer(W,dX2,dY2,dZ2,difx,dify,difz,hx,hy,hz,
     &  neignum,kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,MAXNEIG)

        call forma3Dder(FF,W,dX2,dY2,dZ2,neignum,Xf,Yf,Zf
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),hx,hy,hz,
     & dFFX,dFFY,dFFZ,ddFFXX,ddFFXY,ddFFXZ,ddFFYY,ddFFYZ,ddFFZZ,
     &  MAXNEIG,nbase)

	  do j=1, MLS_neignum(L)
	   IF (UorV.eq.1) then
	   	dh1_loc(L,j)=FF(j)  	;I_nr_U(L,j)=MVefx(j)
	    	J_nr_U(L,j)=MVefy(j)	;K_nr_U(L,j)=MVefz(j)
	   ELSEIF(UorV.eq.2) then
	    	dh2_loc(L,j)=FF(j)	;I_nr_V(L,j)=MVefx(j)
		J_nr_V(L,j)=MVefy(j)	;K_nr_V(L,j)=MVefz(j)
	   ELSEIF(UorV.eq.3) then
	   	dh3_loc(L,j)=FF(j)	;I_nr_W(L,j)=MVefx(j)
		J_nr_W(L,j)=MVefy(j)	;K_nr_W(L,j)=MVefz(j)
	   ENDIF
          enddo

	  IF(UorV.eq.1) KmaxU(L)=MLS_neignum(L)
	  IF(UorV.eq.2) KmaxV(L)=MLS_neignum(L)
	  IF(UorV.eq.3) KmaxW(L)=MLS_neignum(L)

	return
	end
!######################################################################
      subroutine kerneltododer(W,dX1,dY1,dZ1,difx,dify,difz,hx,hy,hz,
     &  neignum,kwx1,kwx2,kwy1,kwy2,kwz1,kwz2,MAXNEIG)
!######################################################################
        use mpi
	implicit none
	Double precision, intent(in) :: hx,hy,hz,kwx1,kwx2,kwy1
	Double precision, intent(in) :: kwy2,kwz1,kwz2
	Double precision, intent(in) :: difx(MAXNEIG)
	Double precision, intent(in) ::dify(MAXNEIG)
	Double precision, intent(in) :: difz(MAXNEIG)
	INTEGER, intent(in):: neignum,MAXNEIG
        INTEGER :: I,J
        Double precision:: cero, coef,cx1,cx2,cy1,cy2,cz1,cz2,PI
        Double precision:: dmx,dmy,dmz,Wx,Wy,Wz,dWx,dWy,dWz
        Double precision:: coef1x,coef2x,coefd1x,coefd2x,coefd3x
        Double precision:: coef1y,coef2y,coefd1y,coefd2y,coefd3y
        Double precision:: coef1z,coef2z,coefd1z,coefd2z,coefd3z
	Double precision:: sx(MAXNEIG),sy(MAXNEIG),sz(MAXNEIG)
	Double precision:: W(MAXNEIG),dX1(MAXNEIG),dY1(MAXNEIG),dZ1(MAXNEIG)

        cero=1.d-08  ;   coef=1.d+00  ;    PI = 4.D0*DATAN(1.D0)

        do j=1, neignum
	  sx(j)=0.0d+00   ;  sy(j)=0.0d+00;  sz(j)=0.0d+00
	enddo

! Primero cogemos las distancias que nos interesan
      do j=1, neignum
	sx(j)=dabs(difx(j)) ; sy(j)=dabs(dify(j)); sz(j)=dabs(difz(j))
	dmx=2.0d+00*hx 	; dmy=2.0d+00*hy 	;  dmz=2.0d+00*hz
	cx1= dmx/kwx1  	;  cx2= dmx/kwx2	
	cy1= dmy/kwy1	;  cy2= dmy/kwy2	
	cz1= dmz/kwz1	;  cz2= dmz/kwz2	

	if (difx(j).ge.0.0d+00)then
	  coef1x= exp(-(sx(j)/cx1)**2) - exp(-(dmx/cx1)**2)
	  coef2x= 1.0d+00 - exp(-(dmx/cx1)**2)
	  coefd1x=exp(-(sx(j)/cx1)**2)
	  coefd2x=coef2x
	  coefd3x=-2.0d+00*sx(j)/(cx1*cx1)
	else
	  coef1x= exp(-(sx(j)/cx2)**2) - exp(-(dmx/cx2)**2)
	  coef2x= 1.0d+00 - exp(-(dmx/cx2)**2)
	  coefd1x=exp(-(sx(j)/cx2)**2)
	  coefd2x=coef2x
	  coefd3x=-2.0d+00*sx(j)/(cx2*cx2)
	endif
	if (dify(j).ge.0.0d+00)then
	  coef1y= exp(-(sy(j)/cy1)**2) - exp(-(dmy/cy1)**2)
	  coef2y= 1.0d+00 - exp(-(dmy/cy1)**2)	
	  coefd1y=exp(-(sy(j)/cy1)**2)
	  coefd2y=coef2y
	  coefd3y=-2.0d+00*sy(j)/(cy1*cy1)
	else
	  coef1y= exp(-(sy(j)/cy2)**2) - exp(-(dmy/cy2)**2)
	  coef2y= 1.0d+00 - exp(-(dmy/cy2)**2)
	  coefd1y=exp(-(sy(j)/cy2)**2)
	  coefd2y=coef2y
	  coefd3y=-2.0d+00*sy(j)/(cy2*cy2)
	endif
	if (difz(j).ge.0.0d+00)then
	  coef1z= exp(-(sz(j)/cz1)**2) - exp(-(dmz/cz1)**2)
	  coef2z= 1.0d+00 - exp(-(dmz/cz1)**2)	
	  coefd1z=exp(-(sz(j)/cz1)**2)
	  coefd2z=coef2z
	  coefd3z=-2.0d+00*sz(j)/(cz1*cz1)
	else
	  coef1z= exp(-(sz(j)/cz2)**2) - exp(-(dmz/cz2)**2)
	  coef2z= 1.0d+00 - exp(-(dmz/cz2)**2)
	  coefd1z=exp(-(sz(j)/cz2)**2)
	  coefd2z=coef2z
	  coefd3z=-2.0d+00*sz(j)/(cz2*cz2)
	endif

	Wx= coef1x/coef2x ; Wy= coef1y/coef2y ;	Wz= coef1z/coef2z

	dWx= coefd3x*coefd1x/coefd2x
	dWy= coefd3y*coefd1y/coefd2y
	dWz= coefd3z*coefd1z/coefd2z

	W(j)= Wx*Wy*Wz 

	if(sx(j).gt.cero)then
	  dX1(j)= Wy*Wz*dWx*(difx(j)/sx(j))
	 else
	  dX1(j)= Wz*Wy*dWx*(difx(j))  
	endif
	if(sy(j).gt.cero)then
	  dY1(j)= Wx*Wz*dWy*(dify(j)/sy(j))
	 else
	  dY1(j)= Wx*Wz*dWy*(dify(j))
	endif
	if(sz(j).gt.cero)then
	  dZ1(j)= Wx*Wy*dWz*(difz(j)/sz(j))
	 else
	  dZ1(j)= Wx*Wy*dWz*(difz(j)) 
	endif
      enddo !j-loop

	return
	end
!######################################################################
      subroutine forma3Dder(FF,W2,dX2,dY2,dZ2,neignum,Xmls,Ymls,Zmls
     &  ,xgau,ygau,zgau,hx,hy,hz,
     &   dFFX,dFFY,dFFZ,ddFFXX,ddFFXY,ddFFXZ,ddFFYY,ddFFYZ,ddFFZZ,
     &   MAXNEIG,nbase)
!######################################################################
        implicit none
        INTEGER, intent(in) :: neignum,MAXNEIG,nbase
        Double precision, intent(in) :: hx,hy,hz,xgau,ygau,zgau
        Double precision, intent(in) :: W2(MAXNEIG)
	Double precision, intent(in) :: Xmls(MAXNEIG),Ymls(MAXNEIG)
	Double precision, intent(in) :: Zmls(MAXNEIG)
	Double precision, intent(in) :: dX2(MAXNEIG),dY2(MAXNEIG)
	Double precision, intent(in) :: dZ2(MAXNEIG)
        INTEGER :: I,J,K,ifila,icolumna,ibase,nelim,icol
        Double precision :: xg,yg,zg,wmax,wmin,valor
        Double precision :: FF(MAXNEIG),dFFX(MAXNEIG),dFFY(MAXNEIG)
        Double precision :: dFFZ(MAXNEIG)
        Double precision :: ddFFXX(MAXNEIG),ddFFXY(MAXNEIG)
        Double precision :: ddFFYY(MAXNEIG)
        Double precision :: ddFFXZ(MAXNEIG),ddFFYZ(MAXNEIG)
        Double precision :: ddFFZZ(MAXNEIG)
	Double precision :: A(nbase,nbase),vIA(nbase,nbase)
        Double precision :: eye(nbase,nbase)
        Double precision :: C(nbase,MAXNEIG),PtC(MAXNEIG,MAXNEIG)
        Double precision :: auxiliar(MAXNEIG,nbase)
	Double precision :: auxiliar1(nbase),auxiliar2(nbase)
        Double precision :: sol(nbase,nbase),wsv(10),vsv(10,10),coef

c Inicializamos la matriz C
        do i=1,nbase
	 do j=1, MAXNEIG
	  C(i,j)=0.d+00
	 enddo
	enddo
c Construimos la matriz A

      do i=1,nbase
	do j=1,nbase
	 A(i,j)= 0.d+00
	 vIA(i,j)=0.d+00
	 if(i.eq.j)then
	  eye(i,j)= 1.d+00
	 else
	  eye(i,j)= 0.d+00
	 endif
	enddo
      enddo

c Construimos un vector auxiliar para hacer los calculos
      do k=1,neignum
       xg=(Xmls(k)-xgau)/hx ; yg=(Ymls(k)-ygau)/hy ;zg=(Zmls(k)-zgau)/hz

      if(nbase.eq.4) then
      	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
      elseif(nbase.eq.10) then
	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
	auxiliar(k, 5)= xg*xg
	auxiliar(k, 6)= xg*yg
	auxiliar(k, 7)= xg*zg
	auxiliar(k, 8)= yg*yg
	auxiliar(k, 9)= yg*zg
	auxiliar(k,10)= zg*zg
      elseif(nbase.eq.20) then
	auxiliar(k, 1)= 1.d+00
	auxiliar(k, 2)= xg
	auxiliar(k, 3)= yg
	auxiliar(k, 4)= zg
	auxiliar(k, 5)= xg*xg
	auxiliar(k, 6)= xg*yg
	auxiliar(k, 7)= xg*zg
	auxiliar(k, 8)= yg*yg
	auxiliar(k, 9)= yg*zg
	auxiliar(k,10)= zg*zg
	auxiliar(k,11)= xg*xg*xg
	auxiliar(k,12)= xg*xg*yg
	auxiliar(k,13)= xg*xg*zg
	auxiliar(k,14)= xg*yg*zg
	auxiliar(k,15)= yg*yg*yg
	auxiliar(k,16)= yg*yg*xg
	auxiliar(k,17)= yg*yg*zg
	auxiliar(k,18)= zg*zg*zg
	auxiliar(k,19)= zg*zg*xg
	auxiliar(k,20)= zg*zg*yg
       endif

       do ifila= 1,nbase  
	do icolumna= 1,nbase
     	 A(ifila,icolumna)= A(ifila,icolumna)+
     &         W2(k)*auxiliar(k,ifila)*auxiliar(k,icolumna)
	enddo
       enddo
      enddo



c Calculamos la inversa de A
c Para calcularla resolvemos 6 sistemas de ecuaciones
c Primero triangularizamos la matriz A y hacemos las operaciones por filas
c sobre la matriz identidad
      do i=1,nbase-1
c Comprobamos que el pivote no sea nulo
        if(A(i,i).lt.1.d-06)then
	  write(6,*)' Pivote nulo',hx,hy,hz,i
	  pause
	endif
       do j=i+1,nbase
        coef=-A(j,i)/A(i,i)
          do k=i,nbase
            A(j,k)=A(j,k)+coef*A(i,k)
          enddo
c Aplicamos la operación a las columnas de la identidad
	    do icol=1,nbase
            eye(j,icol)= eye(j,icol)+coef*eye(i,icol)
	    enddo
        enddo
      enddo
      do icol=1,nbase
	do i=1,nbase
	 auxiliar1(i)= eye(i,icol)
	enddo
         call Back_Substitution(nbase,A,auxiliar1,auxiliar2)
	do i=1,nbase
	 vIA(i,icol)= auxiliar2(i)
	enddo
       enddo


c We calculate C, which is A^-1*B
      do i=1,nbase
	do k=1, neignum
 	 C(i,k)= 0.d+00
	  do j=1,nbase  
	  C(i,k)= C(i,k) + vIA(i,j)*auxiliar(k,j)
	  enddo
	 C(i,k)= C(i,k)*W2(k)
	enddo
      enddo
c We calculate I-Pt*C
      do i=1,neignum
	do j=1,neignum

	valor= 0.d+00
	do ibase=1,nbase
		valor= valor + C(ibase,j)*auxiliar(i,ibase)
	enddo

		PtC(i,j)=-valor
		if(i.eq.j)PtC(i,j)=PtC(i,j)+1.d+00
	  enddo
	enddo

      if(nbase.eq.4) then
        do k=1,neignum
	   FF(k)= C(1,k)
	enddo
        do j=1,neignum
	 dFFX(j)=C(2,j)/hx
	 dFFY(j)=C(3,j)/hy
	 dFFZ(j)=C(4,j)/hz
	 do k=1,neignum
           dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
	 enddo
	enddo
      elseif(nbase.eq.10) then
        do k=1,neignum
 	  FF(k)= C(1,k)
	enddo
        do j=1,neignum
	  dFFX(j)=C(2,j)/hx
	  dFFY(j)=C(3,j)/hy
	 dFFZ(j)=C(4,j)/hz
            do k=1,neignum
             dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
             dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
             dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
            enddo
	enddo
        do j=1,neignum
	  ddFFXY(j)=        C(5,j)/(hx*hy)
	  ddFFXZ(j)=        C(6,j)/(hx*hz)
	  ddFFYZ(j)=        C(7,j)/(hy*hz)
	  ddFFXX(j)= 2.d+00*C(8,j)/(hx*hx)  
	  ddFFYY(j)= 2.d+00*C(9,j)/(hy*hy)
	  ddFFZZ(j)=2.d+00*C(10,j)/(hz*hz)
        enddo
      elseif(nbase.eq.20) then
        do k=1,neignum
	  FF  (k)= C(1,k)
	enddo
        do j=1,neignum
	  dFFX(j)=C(2,j)/hx
	  dFFY(j)=C(3,j)/hy
	  dFFZ(j)=C(4,j)/hZ
	 do k=1,neignum
           dFFX(j)=dFFX(j)+dX2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFY(j)=dFFY(j)+dY2(k)*C(1,k)*PtC(k,j)/W2(k)
           dFFZ(j)=dFFZ(j)+dZ2(k)*C(1,k)*PtC(k,j)/W2(k)
	 enddo
	enddo
         do j=1,neignum
	  ddFFXX(j)= 2.d+00*C( 5,j)/(hx*hx)
	  ddFFXY(j)=        C( 6,j)/(hx*hy)
	  ddFFXZ(j)=        C( 7,j)/(hx*hz)
	  ddFFYY(j)= 2.d+00*C( 8,j)/(hy*hy)
	  ddFFYZ(j)=        C( 9,j)/(hy*hz)
	  ddFFZZ(j)= 2.d+00*C(10,j)/(hz*hz)
         enddo
	endif

	return
	end
!######################################################################
	subroutine Back_Substitution(ndim,A,b,x)
!######################################################################
	INTEGER :: i
        INTEGER, intent(in) :: ndim
        Double precision :: A(ndim,ndim),b(ndim),x(ndim)

      x(ndim)=b(ndim)/A(ndim,ndim)
      do i=ndim-1,1,-1
        suma=0
          do j=i+1,ndim
            suma=suma+A(i,j)*x(j)
          enddo
        x(i)=(b(i)-suma)/A(i,i)
      enddo
	return
	end

