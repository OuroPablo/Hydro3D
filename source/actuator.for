!######################################################################
      SUBROUTINE  act_line_geom(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
      DOUBLE PRECISION:: PI,an,d1,Dtot,alph,dummy
      INTEGER      :: L,nin,I,K,nlay,strlen,n_act,BladeN,fileact
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::r_in,c_in,Pit_in

       PI = 4.D0*DATAN(1.D0)

         write(char_block,'(i3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='geom_ActL_'//TRIM(ADJUSTL(char_block))//'.dat'
         open (unit=2, file=gridfile)
!----      Load the file and proceed to interpolate     ----------
 	open(unit=1, name=filepoints(numIB))
	   read(1,*)n_act
	   read(1,*)
	   allocate(r_in(n_act),c_in(n_act),Pit_in(n_act))
	 DO L=1,n_act
	   read(1,*)r_in(L),c_in(L),Pit_in(L)
       ENDDO
      close (1)

	 BladeN=3				!!!*****
	 nin=INT((r_in(n_act)-r_in(1))/dzm)+1

	 r_act(1+maxnodeIBS)  =r_in(1)+0.19999d0*dzm
	 c_act(1+maxnodeIBS)  =c_in(1)
       Pit_act(1+maxnodeIBS)=Pit_in(1) 
	 nodex(numIB,1)=0.D0
	 nodey(numIB,1)=0.D0
	 nodez(numIB,1)=r_act(1)

	 DO L = 2,nin
		r_act(L)  = r_act(L-1)+dzm
	  Do K =1,n_act-1
	   if(r_act(L+maxnodeIBS).le.r_act(1+maxnodeIBS)) then
	    	c_act(L+maxnodeIBS)  = c_in(K)
		Pit_act(L+maxnodeIBS)=Pit_in(K)
	    	nodex(numIB,L)=0.D0
	   	nodey(numIB,L)=0.D0
	   	nodez(numIB,L)=r_act(L+maxnodeIBS)
	   endif
	   if(r_act(L+maxnodeIBS).gt.r_in(K) 
     & .and.r_act(L+maxnodeIBS).le.r_in(K+1)) then
	  d1=DABS(r_in(K)-r_act(L+maxnodeIBS)) ;  dtot=r_in(K+1)-r_in(K) 
	  c_act(L+maxnodeIBS)  = c_in(K)  +(c_in(K+1)-c_in(K))*d1/dtot
	  Pit_act(L+maxnodeIBS)= Pit_in(K)+(Pit_in(K+1)-Pit_in(K))*d1/dtot
	  nodex(numIB,L)=0.D0
	  nodey(numIB,L)=0.D0
	  nodez(numIB,L)=r_act(L+maxnodeIBS)
		WRITE(6,*)L,r_act(L+maxnodeIBS)
     &   ,c_act(L+maxnodeIBS),Pit_act(L+maxnodeIBS)
	   endif
 	  enddo
	 ENDDO

!Multiple blades generation
	  write(6,*)'Total # points of one blade:',nin,BladeN*nin	
	  alph=2.D0*PI/BladeN  ; k=nin+maxnodeIBS
	DO L=1,BladeN-1
	 DO I=1,nin
	  K=K+1
	  nodex(numIB,K)= nodex(numIB,nin*(L-1)+I)
	  nodey(numIB,K)= nodey(numIB,nin*(L-1)+I)*DCOS(alph)
     & 		     +nodez(numIB,nin*(L-1)+I)*DSIN(alph)
	  nodez(numIB,K)=-nodey(numIB,nin*(L-1)+I)*DSIN(alph)
     & 		     +nodez(numIB,nin*(L-1)+I)*DCOS(alph) 
	  r_act(K)	 = r_act(I+maxnodeIBS)
	  c_act(K)   = c_act(I+maxnodeIBS)
	  Pit_act(K) = Pit_act(I+maxnodeIBS)
	 ENDDO	
	ENDDO
	 maxnode=0  ; nodes(numIB)=nin*BladeN
	 Cxor(numIB)=Cx(numIB);Cyor(numIB)=Cy(numIB);Czor(numIB)=Cz(numIB)	
	 maxnode = max(maxnode,nodes(numIB))
   	   DO i=1,nodes(numIB)		
	    nodexlocal(numIB,i)=nodex(numIB,i)+1.d-7  !Set this as the local coordinates
	    nodeylocal(numIB,i)=nodey(numIB,i)+1.d-7 
	    nodezlocal(numIB,i)=nodez(numIB,i)+1.d-7 
	    nodex(numIB,i)=nodex(numIB,i)+Cxor(numIB)
	    nodey(numIB,i)=nodey(numIB,i)+Cyor(numIB)
	    nodez(numIB,i)=nodez(numIB,i)+Czor(numIB)
	   enddo   

	 write(2,*)'VARIABLES=x,y,z,Radius,Chord,Pitch'
	  do L=1,nodes(numIB)
	   write(2,89)nodex(numIB,L),nodey(numIB,L),nodez(numIB,L)
     &  ,r_act(L+maxnodeIBS),c_act(L+maxnodeIBS),Pit_act(L+maxnodeIBS)
	  enddo
	 close(2) 

	   deallocate(r_in,c_in,Pit_in)

!Open file to output structural forces:
	   write(char_block,'(i3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='FEM_ActL_'//TRIM(ADJUSTL(char_block))//'.dat'
	   fileact=1200+numIB
         open (unit=fileact, file=gridfile)
	   write(fileact,*)'Title = Structural forces Actuator Turbine'	
	   write(fileact,*)'Variables=Deg,T1,T2,T3,Q1,Q2,Q3'
     &  ,',BR1,BR2,BR3,BT1,BT2,BT3,MY1,MY2,MY3'
                           
   87 FORMAT (a,i2,a,i6)
   89 FORMAT (10f13.6)
      RETURN
      end
!######################################################################
      SUBROUTINE ActuatorLine_Initial
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
	  call MPI_BCAST(c_act,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Pit_act,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(r_act,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
	END SUBROUTINE
!######################################################################
      SUBROUTINE ActuatorLine(M,L,ib)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: M,L,ib
	integer :: fileact
      DOUBLE PRECISION :: Utheta,Vrel,Vrot,phi,alpeff,Lift,Drag,F1,GVAL

!actuator line - HATT: forced velocities based on rotational speed
!Relative velocity in the plane of rotation, i.e. y-z for HATT
		Utheta= radsin(M)*R0_loc(L) - (U_Beta2_loc(L)*DCOS(rads(M)
     &  	+alpha0_loc(L)) + U_Beta3_loc(L)*DSIN(rads(M)+alpha0_loc(L))) 

!Module of the effective velocity
		Vrel=DSQRT(U_Beta1_loc(L)**2+Utheta**2)
!Effective angle of attack	
		alpeff=ATAN(U_Beta1_loc(L)/Utheta)
!Angle of attack, Pit_act= local pitch angle. Both in DEGREES
		phi = alpeff*180.D0/3.1416D0 - Pit_act(L)
!Calculate the lift and drag coefficients depending on the angle and Re
		call Stallard(phi)
!Compute equivalent hydrodynamic forces:
		Lift=0.5d0*1.d0*Cl_act*c_act(L)*dom(ib)%dz*Vrel**2
		Drag=0.5d0*1.d0*Cd_act*c_act(L)*dom(ib)%dz*Vrel**2		
!Corrections from Shen et al.
		gval=dexp(-0.125d0*
     & 	     (3.d0*radsin(M)*r_act(nodes(M)/3)/ubulk-21.d0))+0.1d0
		F1=2.d0/3.1416D0*
     & 	   DACOS(dexp(-gval*(3.d0*(r_act(nodes(M)/3)-r_act(L))
     &       /(2.d0*r_act(L)*sin(alpeff)))))
		Lift=Lift*F1 ; Drag=Drag*F1
!Compute the projection of the forces over X-Y-Z coordinates
		FX1_loc(L)=-(Drag*dsin(alpeff)+Lift*dcos(alpeff))
		FX2_loc(L)= (Drag*dcos(alpeff)-Lift*dsin(alpeff))
     &  *Dcos(rads(M)+alpha0_loc(L))
		FX3_loc(L)=-(Drag*dcos(alpeff)-Lift*dsin(alpeff))
     &  *Dsin(rads(M)+alpha0_loc(L))
!Normalise by the area so the ibm subroutine can be used (Wu and Porte-Agel '11): - CHECK THIS STEP - 
		FX1_loc(L)=FX1_loc(L)/(dom(ib)%dz*c_act(L))
		FX2_loc(L)=FX2_loc(L)/(dom(ib)%dz*c_act(L))
		FX3_loc(L)=FX3_loc(L)/(dom(ib)%dz*c_act(L))

!	IF(L.eq.10)THEN
!	write(6,'(a,i2,3f12.5)')' L :',L,phi,U_Beta1_loc(L),Utheta
!	write(6,'(a,i2,4f12.5)')' L :',L,phi
!     & ,1-(U_Beta1_loc(L)/UBULK)
!     &  ,(Drag*dcos(alpeff)-Lift*dsin(alpeff))*1000.d0
!	write(6,*)Utheta,radsin(M)*R0_loc(L),F1
!	write(6,*)Pit_act(L),alpeff*180.D0/3.1416D0,phi
!	write(6,*)Cl_act,Cd_act,c_act(L)
!	write(6,*)Lift,Drag,dzm,Vrel
!	write(6,*)U_Beta1_loc(L),U_Beta2_loc(L),U_Beta3_loc(L)
!	write(6,*)phi,Cl_act,c_act(L),dom(ib)%dz,Vrel**2
!	write(6,*)FX2_loc(L)*dcos(rads(M))-FX3_loc(L)*dsin(rads(M))
!	write(6,*)Drag*dsin(alpeff)+Lift*dcos(alpeff)
!	write(6,*)(-Drag*dcos(alpeff)+Lift*dsin(alpeff))
!	write(6,*)' --- '
!	ENDIF


   88 FORMAT (1f12.6,5e15.7)
      RETURN
      END SUBROUTINE 
!######################################################################
      SUBROUTINE ActuatorLine_FEM(numIB)
!######################################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      integer, intent (in) :: numIB
	integer ::M,L,fileact
	DOUBLE PRECISION :: PI,alpharads,RtoHub
	DOUBLE PRECISION :: T1_ACT,T2_ACT,T3_ACT,Q1_ACT,Q2_ACT,Q3_ACT
	DOUBLE PRECISION :: BR1_ACT,BR2_ACT,BR3_ACT,BT1_ACT,BT2_ACT,BT3_ACT
	DOUBLE PRECISION :: Myaw1,Myaw2,Myaw3,rho
	DOUBLE PRECISION :: F_X_MEAN(100),PHI_MEAN(100),F_Q_MEAN(100),FQ
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile
	INTEGER:: strlen

 	IF(itime.eq. itime_start) then
		F_X_MEAN=0.D0 ;PHI_MEAN=0.D0 ;F_Q_MEAN=0.D0 
	ENDIF

	PI = 4.D0*DATAN(1.D0) ; M=numIB
	alpharads=rads(M)*180.d0/PI ;	rho=1000.d0		
	T1_ACT=0.D0 	; T2_ACT=0.D0 	; T3_ACT=0.D0 
	Q1_ACT=0.D0 	; Q2_ACT=0.D0 	; Q3_ACT=0.D0 
	BR1_ACT=0.D0	; BR2_ACT=0.D0 	; BR3_ACT=0.D0 
	BT1_ACT=0.D0 	; BT2_ACT=0.D0 	; BT3_ACT=0.D0 
	Myaw1=0.D0 		; Myaw2=0.D0 	; Myaw3=0.D0 

! Time-averaged coefficients

!Open file to output mean coefficients:
	   write(char_block,'(i3)') numIB
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='MeanCoefs_ActL_'//TRIM(ADJUSTL(char_block))//'.dat'
	   fileact=1400+numIB
         open (unit=fileact, file=gridfile)
	   write(fileact,*)'Variables="r/R","Fx","Fq"'
	 do L=1,nodes(M)
	   if(L.le.nodes(M)/3) THEN

	 F_X_MEAN(L)=(F_X_MEAN(L)*(ITIME-1)+FX1(M,L)*dzm**2)/ITIME
	 FQ=(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    -FX3(M,L)*dsin(rads(M)+alpha0(M,L)))
     &     *R0(M,L)*dzm**2*radsin(M)
	 F_Q_MEAN(L)=(F_Q_MEAN(L)*(ITIME-1)+FQ)/ITIME

	write(fileact,'(3E18.6)')
     &    r_act(L)/r_act(nodes(M)/3),F_X_MEAN(L)*1000.d0,F_Q_MEAN(L)*1000.d0

	endif;enddo
	close(fileact)

! Structural forces:
	 do L=1,nodes(M)
		RtoHub=R0(M,L)-r_act(1)
	   if(L.le.nodes(M)/3) THEN
		T1_ACT=T1_ACT+FX1(M,L)*dzm*c_act(L)*rho
		Q1_ACT=Q1_ACT+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)*rho
		BR1_ACT=BR1_ACT+FX1(M,L)*dzm*c_act(L)*RtoHub*rho
		BT1_ACT=BT1_ACT+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dcos(rads(M)+alpha0(M,L))
		Myaw1=Myaw1+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dsin(rads(M)+alpha0(M,L))
	   ELSE IF(L.GT.nodes(M)/3 .and. L.le.nodes(M)*2/3 ) then
		T2_ACT=T2_ACT+FX1(M,L)*dzm*c_act(L)*rho
		Q2_ACT=Q2_ACT+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)*rho	
		BR2_ACT=BR2_ACT+FX1(M,L)*dzm*c_act(L)*RtoHub*rho
		BT2_ACT=BT2_ACT+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dcos(rads(M)+alpha0(M,L))
		Myaw2=Myaw2+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dsin(rads(M)+alpha0(M,L))
	   ELSE IF(L.GT.nodes(M)*2/3 .and. L.le.nodes(M)) then
		T3_ACT=T3_ACT+FX1(M,L)*dzm*c_act(L)*rho
		Q3_ACT=Q3_ACT+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &    	-FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)*rho
		BR3_ACT=BR3_ACT+FX1(M,L)*dzm*c_act(L)*RtoHub*rho
		BT3_ACT=BT3_ACT+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dcos(rads(M)+alpha0(M,L))
		Myaw3=Myaw3+FX1(M,L)*dzm*c_act(L)
     &		*R0(M,L)*rho*dsin(rads(M)+alpha0(M,L))
	   ENDIF	
	 end do
	fileact=1200+numIB
	write(fileact,88)alpharads,T1_ACT,T2_ACT,T3_ACT,Q1_ACT,Q2_ACT,Q3_ACT
     & ,BR1_ACT,BR2_ACT,BR3_ACT,BT1_ACT,BT2_ACT,BT3_ACT
     & ,Myaw1,Myaw2,Myaw3

   88 FORMAT (1f12.6,18e15.6)
	END SUBROUTINE
!######################################################################
      SUBROUTINE Stallard(phi)
!######################################################################
      use imb
      Double precision, intent(in) :: phi
      Double precision :: Cl_1,Cl_2,d1,dtot,Cd_1,Cd_2,ang1,ang2

	 ang1=0.d0 ;ang2=0.d0;Cl_1=0.d0;Cd_1=0.d0;Cl_2=0.d0;Cd_2=0.d0

      IF ((phi.gt.-30.d0).and.(phi.le.-20.d0)) THEN
	 ang1 = -30.d0 ; ang2 = -20.d0
	 Cl_1=-0.70 ; Cd_1= 0.48 ; Cl_2 = -0.50 ; Cd_2 = 0.25
      ELSE IF ((phi.gt.-20.d0).and.(phi.le.-10.d0)) THEN
	 ang1 = -20.d0 ; ang2 = -10.d0
	 Cl_1=-0.50 ; Cd_1= 0.25 ; Cl_2 = -0.30 ; Cd_2 = 0.13

      ELSE IF ((phi.gt.-10.d0).and.(phi.le.-7.d0)) THEN
	 ang1 = -10.d0 ; ang2 = -7.d0
	 Cl_1=-0.30 ; Cd_1= 0.13 ; Cl_2 = -0.40 ; Cd_2 = 0.12

      ELSE IF ((phi.gt.-7.d0).and.(phi.le.-6.d0)) THEN
	 ang1 = -7.d0 ; ang2 = -6.d0
	 Cl_1=-0.40 ; Cd_1= 0.12 ; Cl_2 = -0.51 ; Cd_2 = 0.11

      ELSE IF ((phi.gt.-6.d0).and.(phi.le.0.d0)) THEN
	 ang1 = -6.d0 ; ang2 = 0.d0
	 Cl_1=-0.51 ; Cd_1= 0.11 ; Cl_2 = 0.24 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.0.d0).and.(phi.le.1.d0)) THEN
	 ang1 = 0.d0 ; ang2 = 1.d0
	 Cl_1= 0.24 ; Cd_1= 0.06 ; Cl_2 = 0.37 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.1.d0).and.(phi.le.2.d0)) THEN
	 ang1 = 1.d0 ; ang2 = 2.d0
	 Cl_1= 0.37 ; Cd_1= 0.06 ; Cl_2 = 0.49 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.2.d0).and.(phi.le.4.d0)) THEN
	 ang1 = 2.d0 ; ang2 = 4.d0
	 Cl_1= 0.49 ; Cd_1= 0.06 ; Cl_2 = 0.74 ; Cd_2 = 0.06

      ELSE IF ((phi.gt.4.d0).and.(phi.le.8.d0)) THEN
	 ang1 = 4.d0 ; ang2 = 8.d0
	 Cl_1= 0.74 ; Cd_1= 0.06 ; Cl_2 = 1.24 ; Cd_2 = 0.10

      ELSE IF ((phi.gt.8.d0).and.(phi.le.9.d0)) THEN
	 ang1 = 8.d0 ; ang2 = 9.d0
	 Cl_1= 1.24 ; Cd_1= 0.10 ; Cl_2 = 1.19 ; Cd_2 = 0.12

      ELSE IF ((phi.gt.9.d0).and.(phi.le.10.d0)) THEN
	 ang1 = 9.d0 ; ang2 = 10.d0
	 Cl_1= 1.19 ; Cd_1= 0.12 ; Cl_2 = 1.11 ; Cd_2 = 0.14

      ELSE IF ((phi.gt.10.d0).and.(phi.le.12.d0)) THEN
	 ang1 = 10.d0 ; ang2 = 12.d0
	 Cl_1= 1.11 ; Cd_1= 0.14 ; Cl_2 = 0.97 ; Cd_2 = 0.23

      ELSE IF ((phi.gt.12.d0).and.(phi.le.14.d0)) THEN
	 ang1 = 12.d0 ; ang2 = 14.d0
	 Cl_1= 0.97 ; Cd_1= 0.23 ; Cl_2 = 0.91; Cd_2 = 0.29

      ELSE IF ((phi.gt.14.d0).and.(phi.le.16.d0)) THEN
	 ang1 = 14.d0 ; ang2 = 16.d0
	 Cl_1= 0.91 ; Cd_1= 0.29 ; Cl_2 = 0.89 ; Cd_2 = 0.35

      ELSE IF ((phi.gt.16.d0).and.(phi.le.18.d0)) THEN
	 ang1 = 16.d0 ; ang2 = 18.d0
	 Cl_1= 0.89 ; Cd_1= 0.35 ; Cl_2 = 0.91 ; Cd_2 = 0.41

      ELSE IF ((phi.gt.18.d0).and.(phi.le.20.d0)) THEN
	 ang1 = 18.d0 ; ang2 = 20.d0
	 Cl_1= 0.91 ; Cd_1= 0.41 ; Cl_2 = 0.94 ; Cd_2 = 0.46

      ELSE IF ((phi.gt.20.d0).and.(phi.le.25.d0)) THEN
	 ang1 = 20.d0 ; ang2 = 25.d0
	 Cl_1= 0.94 ; Cd_1= 0.46 ; Cl_2 = 1.06 ; Cd_2 = 0.59

      ELSE IF ((phi.gt.25.d0).and.(phi.le.30.d0)) THEN
	 ang1 = 25.d0 ; ang2 = 30.d0
	 Cl_1= 1.06 ; Cd_1= 0.59 ; Cl_2 = 1.16 ; Cd_2 = 0.74

      ELSE IF ((phi.gt.30.d0).and.(phi.le.35.d0)) THEN
	 ang1 = 30.d0 ; ang2 = 35.d0
	 Cl_1= 1.16 ; Cd_1= 0.74 ; Cl_2 = 1.19 ; Cd_2 = 0.89

      ELSE IF ((phi.gt.35.d0).and.(phi.le.40.d0)) THEN
	 ang1 = 35.d0 ; ang2 = 40.d0
	 Cl_1= 1.19 ; Cd_1= 0.89 ; Cl_2 = 1.16 ; Cd_2 = 1.01

      ELSE IF ((phi.gt.40.d0).and.(phi.le.45.d0)) THEN
	 ang1 = 40.d0 ; ang2 = 45.d0
	 Cl_1= 1.16; Cd_1= 1.01 ; Cl_2 = 1.11 ; Cd_2 = 1.07

      ELSE IF ((phi.gt.45.d0).and.(phi.le.50.d0)) THEN
	 ang1 = 45.d0 ; ang2 = 50.d0
	 Cl_1= 1.11; Cd_1= 1.07 ; Cl_2 = 1.10 ; Cd_2 = 1.07

      ELSE IF ((phi.gt.50.d0)) THEN
	 ang1 = 50.d0 ; ang2 = 100.d0
	 Cl_1= 1.11; Cd_1= 1.07 ; Cl_2 = 1.10 ; Cd_2 = 1.07
      END IF

	 d1=DABS( ang1 - phi) ;  dtot = ang2 - ang1
	 Cl_act = Cl_1 + (Cl_2 - Cl_1)*d1/dtot
	 Cd_act = Cd_1 + (Cd_2 - Cd_1)*d1/dtot

      RETURN
      End
