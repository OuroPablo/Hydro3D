!#############################################################################
	subroutine SEM
!#############################################################################
	USE IFPORT
	USE module_SEM
	USE vars
	USE multidata
	USE mpi
	IMPLICIT NONE
	DOUBLE PRECISION :: VOL,Ly,Lz,ENNE,XMIN,XMAX,YMIN,TI
	DOUBLE PRECISION :: YMAX,ZMIN,ZMAX,PI,U0,HU,RIZ,RDIVz,UAVE(3)
	DOUBLE PRECISION :: SIGMA_VALUE,NE_SEM,MAXVSEM,MINVSEM,rrandom
	INTEGER :: DIVY,DIVZ,IY,IZ,II,N,I,J,M,IT,IGLOBAL
	INTEGER,allocatable,dimension(:)::lsy,lsz,ley,lez
!THE BOX DIMENSIONS ARE DIFINED AS [XLENGHT] * [Ly] * [Lz]
![DIVX], [DIVY] & [DIVZ] ARE THE NUMBERS OF SPACIAL POINTS.
	PI = 2*ACOS(0.0D0) ; Ly  = yen-yst ; Lz= zen-zst 
	U0 = ubulk 	; TI=TI_SEM

	ALLOCATE(iddom(jdom*kdom),ljdom(jdom*kdom),lkdom(jdom*kdom))
	ALLOCATE(elemyst(jdom*kdom),elemzst(kdom*jdom))
	ALLOCATE(elemyen(jdom*kdom),elemzen(kdom*jdom))
	ALLOCATE(ley(jdom),lez(kdom),lsy(jdom),lsz(kdom))

	 ljdom=0 	; lkdom=0 	;  DIVY = 0 	; DIVZ = 0
	 elemyst=0 	; elemyen=0	;  elemzst=0 	; elemzen=0
	 I=0 ;ley=0 ; lez=0	;  lsy=0 		; lsz=0

!ID OF THE BLOCK
	 DO N=1,kdom ; DO J=1,jdom
	 I=I+1 	
	  iddom(I)=idom*(J-1)+(jdom*idom*(N-1))
	 ENDDO ; ENDDO
!DIVISION ON THE Y AND Z DIRECTIONS 
        lsy(1)=1
	  ley(1)=lsy(1)+NINT((ycor(iddom(1),2)-ycor(iddom(1),1))/g_dy)+1
	 Do J=2,jdom
	  lsy(J)=ley(J-1)+1
	  ley(J)=lsy(J)+NINT((ycor(iddom(J),2)-ycor(iddom(J),1))/g_dy)+1
	  DIVY = DIVY + NINT((ycor(iddom(J),2)-ycor(iddom(J),1))/g_dy) +1
	 Enddo
        lsz(1)=1
	  lez(1)=lsz(1)+NINT((zcor(iddom(1),2)-zcor(iddom(1),1))/g_dz)+1
	 DO N=1,kdom 
	  lsz(N)=lez(N-1)+1 
	  lez(N)=lsz(N)+NINT((zcor(iddom(N),2)-zcor(iddom(N),1))/g_dz)+1
	  DIVZ = DIVZ + NINT((zcor(iddom(N),2)-zcor(iddom(N),1))/g_dz) +1
	 ENDDO
!THE DIVISIONS FOR EACH OF THE DOMAINS IS DETERMINED:
	 I=0
	 DO N=1,kdom  ;  DO J=1,jdom
	 I=I+1 	
	  elemyst(I) = lsy(J)  ;  elemzst(I) = lsz(N)
	  elemyen(I) = ley(J)  ;  elemzen(I) = lez(N)
	  ljdom(I) = elemyen(I) -  elemyst(I) + 1	
	  lkdom(I) = elemzen(I) -  elemzst(I) + 1 
	 enddo ; enddo

	write(6,*)'Divisions :',DIVY,DIVZ, ' Number of BLOCKS  :',jdom,kdom	
	write(6,'(10i4)')lsy(:)
	write(6,'(10i4)')ley(:)-lsy(:)
	write(6,'(10i4)')lsz(:)
	write(6,'(10i4)')lez(:)-lsz(:)
!THE AMPLITUDE OF THE VORTICES--> turbulence integral length scale 
	allocate (SIGMA(3,DIVY,DIVZ))

!    ******     -> change numbers accordingly      ******  !
	SIGMA(1,:,:) = 0.252D0
	SIGMA(2,:,:) = 0.1485d0
	SIGMA(3,:,:) = 0.1125d0 

!THE NUMBER OF EDDIES IS THE INLET SURFACE DIVIDED BY THE SURFACE OF EACH TURBULENT SPOT
![N](INTEGER) AND [ENNE](DOUBLE PRECISION) REPRESENT THE NUMBER OF EDDIES
	NE_SEM = Ly/Sigma(2,1,1)*Lz/Sigma(3,1,1)

	N = INT(NE_SEM)	;  ENNE = REAL(N)

 	if(myrank.eq.0) write(6,*)'The number of SEM eddies is :',N

!  [REYNOLDS(6)] IS A VECTOR WITH THE SIX ELEMENTS OF REYNOLDS STRESSES.
!  |REYNOLDS(1)  REYNOLDS(2)  REYNOLDS(4)|
!  |REYNOLDS(2)  REYNOLDS(3)  REYNOLDS(5)|
!  |REYNOLDS(4)  REYNOLDS(5)  REYNOLDS(6)|
	REYNOLDS=(/(TI*U0)**2, 0.0D0,(TI*U0)**2, 0.0D0, 0.0D0,(TI*U0)**2/)  ![M/S]

!ALLOCATION OF THE EDDIES VECTOR.
![VSEM(DIVX,DIVY,DIVZ,3)] IS THE INSTANTANEOUS VELOCITY VECTOR IN THE POINT WITH X,Y,Z COMPONENTS
![X_EDDY(3,N)] IS THE N-TH EDDY LOCATION X,Y,Z
![EPSILO(3,N)] IS THE N-TH EDDY INTENSITY IN X,Y,Z
![MOLT(3,N)]   IS THE MATRIX PRODUCT [R(3,3)] * [EPSILO(3,N)]
![Ksem(N)]     IS A VECTOR USED TO CHECK WHITCH EDDY IS OUTSIDE THE BOX AFTER THE CONVECTION

	ALLOCATE(Vsem(DIVY,DIVZ,3),Usem(DIVY,DIVZ,3),Ksem(N))
	ALLOCATE(X_EDDY(3,N),EPSILO(3,N),MOLT(3,N))
	PRINT *, "SEM total time interval = ", DT * ITMAX_SEM, " (s)"

	XMIN=1.D8 ; XMAX=1.D-8 ; YMIN=1.D8 ; YMAX=1.D-8
	ZMIN=1.D8 ; ZMAX=1.D-8

!DEFINITION OF EDDY LENGTH SCALE AND INITIAL VELOCITY FIELD
	DO IY=1,DIVY
	  DO IZ=1,DIVZ
	   Usem(IY,IZ,:)=(/ U0 ,0.D0,0.D0/)
	   UAVE(:) = UAVE(:) + Usem(IY,IZ,:)
!CALCULATION OF THE EDDY BOX PARAMETERS
	XMIN=MIN(0.0D0-SIGMA(1,IY,IZ),XMIN);XMAX=MAX(0.0D0+SIGMA(1,IY,IZ),XMAX)
	YMIN=MIN(0.0D0-SIGMA(2,IY,IZ),YMIN);YMAX=MAX(Ly   +SIGMA(2,IY,IZ),YMAX)
	ZMIN=MIN(0.0D0-SIGMA(3,IY,IZ),ZMIN);ZMAX=MAX(Lz   +SIGMA(3,IY,IZ),ZMAX)
	  END DO
	END DO

	UAVE = UAVE / (DIVY * DIVZ)
	VOL  = (XMAX - XMIN) * (YMAX - YMIN) * (ZMAX - ZMIN)

	call RANDOM_SEED()
!GENERATION OF THE EDDY LOCATION INSIDE THE BOX AND INITIALIZATION OF THE [Ksem] VECTOR
	DO II=1,N

		call RANDOM_NUMBER(rrandom)
	  X_EDDY(1,II) = (XMAX - XMIN) * rrandom + XMIN
	  X_EDDY(2,II) = (YMAX - YMIN) * rrandom + YMIN
	  X_EDDY(3,II) = (ZMAX - ZMIN) * rrandom + ZMIN
	  Ksem(II) = 0
!INITIALIZATION OF THE INTENSITIES. FOR EVERY DIRECTION THE AVERAGE INTENSITY VALUE IS CALCULATED AND
!IT IS FORCED TO BE LOWER THAN THE [VLIM] VALUE
		call RANDOM_NUMBER(rrandom)
	  EPSILO(1,II) = (rrandom*2.0D0 - 1.0D0)
	  EPSILO(2,II) = (rrandom*2.0D0 - 1.0D0)
	  EPSILO(3,II) = (rrandom*2.0D0 - 1.0D0)
	END DO
!INITIALIZATION OF THE [R(3,3)] MATRIX WITH THE CHOLENSKY DECOMPOSITION
!OF THE REYNOLDS STRESS TENSOR
	R = 0.d0
	R(1,1) = DSQRT(REYNOLDS(1))
	R(2,1) = REYNOLDS(2) / R(1,1)
	R(2,2) = DSQRT(REYNOLDS(3) - R(2,1)*R(2,1))
	R(3,1) = REYNOLDS(4) / R(1,1)
	R(3,2) = (REYNOLDS(5) - R(2,1)*R(3,1)) / R(2,2)
	R(3,3) = DSQRT(REYNOLDS(6) - R(3,1)*R(3,1) - R(3,2)*R(3,2))
!BEGINNING OF TIME ITERATIONS
	DO IT=1,ITMAX_SEM			!PARALLELISE THIS LOOP
!PRINTINGS OF GLOBAL VELOCITY AND OF CONVECTION VELOCITY
	 Do I=1,jdom*kdom
		  IGLOBAL=500+I
	   WRITE (FILEGLOBAL,'(A13,i4.4,a1,I5.5,A4)')
     &		'inflow/Inlet_',iddom(I),'_',IT,'.dat'
	   OPEN(UNIT=IGLOBAL, FILE=FILEGLOBAL,STATUS="UNKNOWN", ACTION="WRITE")
	   WRITE (IGLOBAL,*)'Variables=up,vp,wp'
	   WRITE (IGLOBAL,*)
     &  'zone  i=',ljdom(I),', j=',lkdom(I),', k=1 f=point'
	 Enddo

	  MOLT = MATMUL(R,EPSILO) ! !aij*epsij   MATRIX MULTIPLICATION
	  MAXVSEM=0.D0 ; MINVSEM=1000000.D0
!------BEGINNING OF SPATIAL ITERATION
	  DO IY = 1,DIVY
	    DO IZ = 1,DIVZ
!X_POINT = GRID POINT COORDINATES
	   X_POINT=(/0.0D0,IY*Ly/DIVY+YMIN,IZ*Lz/DIVZ+ZMIN/)
	   Vsem(IY,IZ,:)=(/ 0.d0, 0.d0,  0.d0 /)
!------------BEGINNING OF EDDIES ITERATIONS
		DO II=1,N
		TEMP(:) = DABS(X_POINT(:) - X_EDDY(:,II))
	  IF (TEMP(1).LT.SIGMA(1,IY,IZ) .AND. 
     &      TEMP(2).LT.SIGMA(2,IY,IZ) .AND.
     &      TEMP(3).LT.SIGMA(3,IY,IZ)) THEN
	        Vsem(IY,IZ,:)=Vsem(IY,IZ,:)+MOLT(:,II)
     &        *(DSQRT(1.5D0)**3.0D0)*DSQRT(VOL)
     &        *(1.0D0-TEMP(1)/SIGMA(1,IY,IZ))
     &        *(1.0D0-TEMP(2)/SIGMA(2,IY,IZ))
     &        *(1.0D0-TEMP(3)/SIGMA(3,IY,IZ))
     &        /DSQRT(SIGMA(1,IY,IZ)*SIGMA(2,IY,IZ)*SIGMA(3,IY,IZ)) 
		 END IF
		END DO
! Instananeous velocity=mean velocity + SEM fluctuation velocity
	  Vsem(IY,IZ,:) = Vsem(IY,IZ,:) / DSQRT(ENNE) 
	  if(Vsem(IY,IZ,1).gt.MAXVSEM) MAXVSEM=Vsem(IY,IZ,1)
	  if(Vsem(IY,IZ,1).le.MINVSEM) MINVSEM=Vsem(IY,IZ,1)
	    END DO
	  END DO

	 IF(MOD(IT,50).eq.0) write(6,'(a,i6,a8,3f13.6)')
     &   'ITERATION ',IT,' Time(s)',(IT-1)*dt,MAXVSEM,MINVSEM
!------END OF SPATIAL ITERATIONS
	      Do I=1,jdom*kdom
	        IGLOBAL=500+I
		   Do M=elemzst(I),elemzen(I)
		    Do J=elemyst(I),elemyen(I)
                  WRITE(500+I,'(3E15.6)')Vsem(J,M,:)
	  	  enddo ;enddo 
	    	 CLOSE (UNIT=IGLOBAL)
	  	Enddo
!--------BEGINNING OF EDDIES CONVECTION ITERATIONS
	  DO II=1,N
!RE-CALCULATION OF THE EDDIES POSITION. IF ANY EDDY GOES BEYOND THE BOX LIMITS
!IT IS RESTARTED AT THE SURFACE FACING THE EXIT.
	    X_EDDY(1,II) = X_EDDY(1,II) + UAVE(1) * dt
	    X_EDDY(2,II) = X_EDDY(2,II) + UAVE(2) * dt
	    X_EDDY(3,II) = X_EDDY(3,II) + UAVE(3) * dt
!AFTER THE EDDIES CONVECTIONS ARE NECESSARY SOME TESTS TO CHECK IF ANY EDDY IS NOW OUTSIDE
!THE SEM BOX DEFINED EARLIER
!WHEN A EDDY IS RE-GENERATE THE Ksem(I) FACTOR ASSUME THE 1 VALUE. THIS VALUE IS USED LATER TO
!GENERATE A NEW INTENSITY FOR THE NEW EDDY

	    IF (X_EDDY(1,II) > XMAX) THEN
		call RANDOM_NUMBER(rrandom)
		X_EDDY(1,II) = XMIN
		X_EDDY(2,II) = (YMAX - YMIN) * rrandom + YMIN
		X_EDDY(3,II) = (ZMAX - ZMIN) * rrandom + ZMIN
		Ksem(II) = 1
	    ELSE IF (X_EDDY(1,II) < XMIN) THEN
		call RANDOM_NUMBER(rrandom)
		X_EDDY(1,II) = XMAX
		X_EDDY(2,II) = (YMAX - YMIN) * rrandom + ZMIN
		X_EDDY(3,II) = (ZMAX - ZMIN) * rrandom + ZMIN
		Ksem(II) = 1
	    ELSE IF (X_EDDY(2,II) > YMAX) THEN
		call RANDOM_NUMBER(rrandom)
		X_EDDY(1,II) = (XMAX - XMIN) * rrandom + XMIN
		X_EDDY(2,II) = YMIN
		X_EDDY(3,II) = (ZMAX - ZMIN) * rrandom + ZMIN
		Ksem(II) = 1
	    ELSE IF (X_EDDY(2,II) < YMIN) THEN
		call RANDOM_NUMBER(rrandom)
		X_EDDY(1,II) = (XMAX - XMIN) * rrandom + XMIN
		X_EDDY(2,II) = YMAX
		X_EDDY(3,II) = (ZMAX - ZMIN) * rrandom + ZMIN
		Ksem(II) = 1
	    ELSE IF  (X_EDDY(3,II) > ZMAX) THEN
		call RANDOM_NUMBER(rrandom)
		X_EDDY(1,II) = (XMAX - XMIN) * rrandom + XMIN
		X_EDDY(2,II) = (YMAX - YMIN) * rrandom + ZMIN
		X_EDDY(3,II) = ZMIN 
		Ksem(II) = 1
	    ELSE IF (X_EDDY(3,II) < ZMIN) THEN
		call RANDOM_NUMBER(rrandom)
		X_EDDY(1,II) = (XMAX - XMIN) * rrandom + XMIN
		X_EDDY(2,II) = (YMAX - YMIN) * rrandom + ZMIN
		X_EDDY(3,II) = ZMAX
		Ksem(II) = 1
	    END IF
!INTENSITY GENERATION FOR THE RE-CREATED EDDIES. WE ARE USING THE Ksem FACTOR AS EXPLAINED FATOR.
	    IF (Ksem(II)== 1) THEN
		call RANDOM_NUMBER(rrandom)
		EPSILO(3,II) = (rrandom*2.0D0 - 1.0D0)
		EPSILO(2,II) = (rrandom*2.0D0 - 1.0D0)
		EPSILO(1,II) = (rrandom*2.0D0 - 1.0D0)
	    END IF
	    Ksem(II) = 0
	  END DO
!---------END OF EDDIES CONVECTIONS ITERATIONS
	END DO
!----END OF TIME ITERATIONS

	END SUBROUTINE 
