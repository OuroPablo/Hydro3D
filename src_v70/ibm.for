!#############################################################
      SUBROUTINE IMB_INITIAL
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      DOUBLE PRECISION :: PI,revoltime,anst,iteratime
      INTEGER      :: L,I,M,strlen,maxn
      CHARACTER*8  :: char_block
      CHARACTER*31 :: gridfile
      CHARACTER*80 :: dummyline
        PI = 4.D0*DATAN(1.D0)

	master=0 ! 0 is going to be always the master processor

       open (unit=1, file='in_geom.cin')
       read (1,*) dummyline
       read (1,*) ibm_out_forces
       read (1,*) yangcase
       read (1,*) mdfsteps
       read (1,*) bodynum

!Allocate variables that all MPI need to know:
      allocate(nodes(bodynum),rotating(bodynum),rads(bodynum))
	allocate(imb_shape(bodynum),imbnumber(bodynum))
	allocate(radsin(bodynum),filepoints(bodynum))
	allocate(turax(bodynum),reddelta(bodynum)) 
	allocate(ibturbine(bodynum)) 
	nodes=0 	; rotating=.false. 	; rads=0.d0
	imb_shape=1	; imbnumber=1		; ibturbine=.false.
	radsin=0.d0	; turax =1			; reddelta=1.d0      

	IF (myrank.ne.master) RETURN
!Allocate variables only needed by the master:
	allocate(Cx(bodynum),Cxor(bodynum),Cy(bodynum),Cyor(bodynum))
	allocate(Cz(bodynum),Czor(bodynum),pitch(bodynum))
	allocate(R(bodynum),l2norm(bodynum))	
	allocate(xaero(bodynum),yaero(bodynum),zaero(bodynum))
      allocate(cmax(bodynum),axis(bodynum),nscatter(bodynum))
	allocate(linfin(bodynum),zini(bodynum),zend(bodynum))        
	
	Cx=0.d0 	; Cxor=0.d0 ; Cy=0.d0 
	Cyor=0.d0 	; Cz=0.d0 	; Czor=0.d0 
	pitch=0.d0 	; R=0.d0 	; l2norm=0.d0 	
	xaero=0.d0 	; yaero=0.d0; zaero=0.d0
	cmax=1	; axis=1	; linfin=0.d0 ; zini=0.d0 ; zend=0.d0
	nscatter=0	

        DO M=1,bodynum
		read (1,*) dummyline 
		read (1,*) imb_shape(M)
		read (1,*) linfin(M),zini(M),zend(M)		
		read (1,*) Cx(M),Cy(M),Cz(M)
		read (1,*) R(M)
		read (1,*) cmax(M)
		read (1,*) axis(M)
		read (1,*) filepoints(M)
		read (1,*) rotating(M)
		read (1,*) reddelta(M)			
		read (1,*) dummyline  		!--- Turbine parameters: 
		read (1,*) ibturbine(M)				
		read (1,*) turax(M)				
		read (1,*) xaero(M),yaero(M),zaero(M)		
		read (1,*) pitch(M)		
		read (1,*) imbnumber(M)
		read (1,*) radsin(M)		

	   if(.not.ibturbine(M)) then
	    xaero(M)=0.d0 ; yaero(M)=0.d0 ; zaero(M)=0.d0 	
	    pitch(M)=0.d0 ; turax(M)=1    ; imbnumber(M)=1
	   endif
	   if(turax(M).eq.2 ) then !HATT
	    xaero(M)=0.d0 ; yaero(M)=0.d0 ; zaero(M)=0.d0 		
	    pitch(M)=0.d0 ; imbnumber(M)=1 	
	   endif
	   if (.not.rotating(M)) radsin(M)=0.d0

!Actuator line:
	   if(ibturbine(M) .and. turax(M).eq.3 .and. M.eq.1) then
	    allocate(r_act(5000),c_act(5000),Pit_act(5000))
	    r_act=0.d0 ; c_act=0.d0 ;Pit_act=0.d0 
	   endif
!Write possible combinations of movements/body types that doesn't work
      End do !M
      close (1)
	WRITE(6,*)' '
	WRITE(6,*)'===========  Immersed Boundary Details  =========='

!CALCULATE TO WHICH BLOCK IS THE CENTRE OF THE BODIES

	 maxn=1999000	!maximum number of Lagrangian allowed

       allocate (nodex(bodynum,maxn),nodexlocal(bodynum,maxn))
	 allocate (nodey(bodynum,maxn),nodeylocal(bodynum,maxn))
	 allocate (nodez(bodynum,maxn),nodezlocal(bodynum,maxn))

      nodex = 0.d0	; nodey = 0.d0 	  ; nodez = 0.d0  
      nodexlocal = 0.d0	; nodeylocal = 0.d0 ; nodezlocal = 0.d0 
	maxnodeIBS=0
	
      dxm=g_dx/rdivmax ; dym=g_dy/rdivmax ; dzm=g_dz/rdivmax	!Minimum grid sizes
	 
	write(6,*)'Largest rdivmax  :',rdivmax	
	write(6,'(a,3e12.4)')'Smallest gridsize: ',dxm,dym,dzm	 
		
	Do M=1,bodynum
        IF (imb_shape(M).eq.1) call imb_square(M)
        IF (imb_shape(M).eq.2) call imb_cylinder(M)
        IF (imb_shape(M).eq.3) call imb_cube(M)
	  IF (imb_shape(M).eq.4) call imb_sphere(M)
        IF (imb_shape(M).eq.11)call imb_cone(M)
        IF (imb_shape(M).eq.12)call imb_dune(M)
        IF (imb_shape(M).eq.13)call imb_hemisphere(M)
 	  IF (imb_shape(M).eq.5 .and. turax(M).le.2) call imb_file(M)		!07_2017
 	  IF (imb_shape(M).eq.5 .and. turax(M).eq.3) call act_line_geom(M)
   	maxnodeIBS=maxnodeIBS+nodes(M)
	IF (maxnodeIBS.gt.maxn)
     & write(6,*)'Too many ib points, change maxn in imb.for'
	IF (maxnodeIBS.gt.maxn) STOP
	Enddo

!The velocity and force vectors/matrix are allocated:
      allocate (U_Beta1(bodynum,maxnodeIBS),U_Beta2(bodynum,maxnodeIBS))
      allocate (U_Beta3(bodynum,maxnodeIBS),imb_block(maxnodeIBS))
      allocate (alpha0(bodynum,maxnodeIBS),R0(bodynum,maxnodeIBS))
!Initiate all these variables:
	U_Beta1=0.d0  ; U_Beta2=0.d0  ; U_Beta3=0.d0
	R0=0.d0       ; alpha0=0.d0	  ; imb_block=0

	if (myrank.eq.master) then
       allocate (FX1(bodynum,maxnodeIBS)) ; FX1=0.D0 
       allocate (FX2(bodynum,maxnodeIBS)) ; FX2=0.D0 
       allocate (FX3(bodynum,maxnodeIBS)) ; FX3=0.D0 
       allocate (FX1M(bodynum,maxnodeIBS)) ; FX1M=0.D0 
       allocate (FX2M(bodynum,maxnodeIBS)) ; FX2M=0.D0 
       allocate (FX3M(bodynum,maxnodeIBS)) ; FX3M=0.D0 
	endif
!Local angle and radius:
	 call imb_alpha0	

	if (ibm_out_forces.eq.1) then
!CREATE OUTPUT FILES FOR THE IMMERSED BOUNDARIES:
	L=0
	DO M=1,bodynum
	 DO i=1,imbnumber(M)
             L=L+1 ; forcefilej=399+L
	 IF (ibturbine(M)) then !Rotating VATT
         write(char_block,'(i2)') L
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(2-strlen))//TRIM(ADJUSTL(char_block))
         gridfile='F_Blade_'//TRIM(ADJUSTL(char_block))//'.dat'
         open (unit=forcefilej, file=gridfile, status="unknown",
     &	action="write")
	if (turax(M).eq.1) then
         write (forcefilej,*)'Variables=Deg,Fx,Fy,T1,T2,N1,N2,M'
	else if (turax(M).eq.2) then
         write (forcefilej,*)'Variables=Deg,Fx,Fy,Fz,Ft'
	else if (turax(M).eq.3) then
         write (forcefilej,*)'Variables=Deg,Cx,Cp,T,P'
	   goto 363
	endif
	     IF(I.EQ.1) then
		lambda=radsin(M)*R(M)/ubulk
		sigma=imbnumber(M)*1.d0/(R(M)*2*PI)
		revoltime=2.d0*PI/radsin(M)
		anst=radsin(M)*dt*180.0/PI
		iteratime=(360.d0/anst)
	write(6,*)'     '	
	write(6,'(a,i2,a)')   ' ****** Turbine   ',M,'rotating details,'
	write(6,'(a,f12.4)')  '        TSR     : ',lambda
	write(6,'(a,f12.4,a)')'        Solidity: ',sigma*100,' %'
	write(6,'(a,f12.4,a)')'Time per revolut: ',revoltime,' sec'
	write(6,'(a,f12.4,a)')'Iter per revolut: ',iteratime,' it'
	write(6,'(a,f12.4,a)')'     Angle step : ',anst,'deg/iteration'
	    ENDIF
	  ENDIF
363	CONTINUE	!09-2017

	 IF (.not.ibturbine(M).and. .not.rotating(L)) then !Rotating VATT
         write(char_block,'(i2)') L
         strlen=LEN(TRIM(ADJUSTL(char_block)))
         char_block=REPEAT('0',(2-strlen))//TRIM(ADJUSTL(char_block))
         if(imb_shape(M).eq.1) then
         gridfile='F_Squ_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif        
         if(imb_shape(M).eq.2) then
         gridfile='F_Cyl_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif
         if(imb_shape(M).eq.3) then
         gridfile='F_Cub_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif
         if(imb_shape(M).eq.4) then
         gridfile='F_Sph_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif     
         if(imb_shape(M).eq.5 .and. .not.rotating(L)) then
         gridfile='F_Bod_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif 	      
         if(imb_shape(M).ge.11) then
         gridfile='F_Con_'//TRIM(ADJUSTL(char_block))//'.dat'
     	 endif 
         open (unit=forcefilej, file=gridfile)
           write (forcefilej,*)'Variables=CTIME,Fx,Fy,Fz'
	 ENDIF

        if(imb_shape(M)==5 .and..not.ibturbine(M) .and.rotating(L)) then
           write(char_block,'(i2)') L
           strlen=LEN(TRIM(ADJUSTL(char_block)))
           char_block=REPEAT('0',(2-strlen))//TRIM(ADJUSTL(char_block))
           gridfile='F_Body_'//TRIM(ADJUSTL(char_block))//'.dat'
           open (unit=forcefilej, file=gridfile)
           write (forcefilej,*)'Variables=Deg,Fx,Fy,Fz,Ft'
     	   endif       
    
	 Enddo !i

	ENDDO !M

	WRITE(6,*)' '
	WRITE(6,*)'Total # of IB POINTS.........',maxnodeIBS
	WRITE(6,*)' '

	endif !If the IBM force output

      RETURN
      end
!#############################################################
      SUBROUTINE imb_alpha0	
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER  :: M,L,iii,K,strlen
      DOUBLE PRECISION :: PI
      CHARACTER*8  :: char_block2
      CHARACTER*31 :: gridfile

       PI = 4.D0*DATAN(1.D0)
      
	Do M=1,bodynum

	IF(imb_shape(M).eq.5) then

!!!! TURBINES  !!!!!!!!!!
        write(char_block2,'(I3)') M
         strlen=LEN(TRIM(ADJUSTL(char_block2)))
         char_block2=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(char_block2))
         gridfile='g_angles_'//TRIM(ADJUSTL(char_block2))//'.dat'
         open (unit=2, file=gridfile)
	write(2,*)'variables=x,y,z,al0,R0'

	 if(.not.ibturbine(M)) then  !Not a turbine

           do L=1,nodes(M)    			!!!!!assumed z-axis...
            alpha0(M,L)=datan(nodexlocal(M,L)/nodeylocal(M,L))
            R0(M,L)=dsqrt(nodexlocal(M,L)**2+nodeylocal(M,L)**2)
         if(nodexlocal(M,L).gt.0.d0 .and. nodeylocal(M,L).gt.0.d0)then !1st quarter. Alpha>0
          alpha0(M,L)=alpha0(M,L)
         endif
         if(nodexlocal(M,L).le.0.d0 .and. nodeylocal(M,L).ge.0.d0)then !2nd quarter. Alpha<0
          alpha0(M,L)=2*PI+alpha0(M,L)
         endif 
         if(nodexlocal(M,L).lt.0.d0 .and. nodeylocal(M,L).lt.0.d0)then !3rd quarter. Alpha>0 
          alpha0(M,L)=PI+alpha0(M,L)
         endif 
         if(nodexlocal(M,L).gt.0.d0 .and. nodeylocal(M,L).lt.0.d0)then !4th quarter. Alpha>0
          alpha0(M,L)=PI+alpha0(M,L)
         endif
            enddo

	 else


	IF (turax(M).eq.1) then		! Vertical Axis Turbine
	   K=nodes(M)/imbnumber(M)
          do L=1,K 
          alpha0(M,L)=-datan(nodexlocal(M,L)/(nodeylocal(M,L)+R(M))) !!!!!!!! 30 Aug Pablo
          R0(M,L)=dsqrt((nodexlocal(M,L))**2+(nodeylocal(M,L)+R(M))**2)
          enddo
          Do iii=1,imbnumber(M)-1
           do L=1,K
             alpha0(M,L+K*iii)=(2.D0*PI/imbnumber(M))*iii+alpha0(M,L)
	 if(alpha0(M,L+K*iii).gt.2.d0*PI) 
     &	 alpha0(M,L+K*iii)=alpha0(M,L+K*iii)-2.d0*PI
             R0(M,L+K*iii)= R0(M,L)
           enddo
          Enddo
        ENDIF

	IF (turax(M).eq.2) then	! Horizontal Axis Turbine
          do L=1,nodes(M) 
          alpha0(M,L)=datan(nodeylocal(M,L)/(nodezlocal(M,L)))
          if(nodeylocal(M,L).ge.0.d0 .and. nodezlocal(M,L).ge.0.d0) then !1st quarter. Alpha>0
         	alpha0(M,L)=alpha0(M,L)
          endif           
          if(nodeylocal(M,L).lt.0.d0 .and. nodezlocal(M,L).gt.0.d0) then !2nd quarter. Alpha<0
          	alpha0(M,L)=2.D0*PI+alpha0(M,L)
          endif 
          if(nodeylocal(M,L).lt.0.d0 .and. nodezlocal(M,L).lt.0.d0) then !3rd quarter. Alpha>0
         	 alpha0(M,L)=PI+alpha0(M,L)
          endif
          if(nodeylocal(M,L).gt.0.d0 .and. nodezlocal(M,L).lt.0.d0) then !4th quarter. Alpha<0
          	alpha0(M,L)=PI+alpha0(M,L)
          endif
          R0(M,L)=dsqrt((nodeylocal(M,L))**2+(nodezlocal(M,L))**2)                               
          enddo
        ENDIF


	IF (turax(M).eq.3) then		! ACTUATOR LINE
	   K=nodes(M)/imbnumber(M)
          do L=1,K 
          alpha0(M,L)=0.d0
          R0(M,L)=r_act(L)
          enddo
          Do iii=1,imbnumber(M)-1
           do L=1,K
             alpha0(M,L+K*iii)=(2.D0*PI/imbnumber(M))*iii+alpha0(M,L)
	    if(alpha0(M,L+K*iii).gt.2.d0*PI) 
     &	 alpha0(M,L+K*iii)=alpha0(M,L+K*iii)-2.d0*PI
             R0(M,L+K*iii)= R0(M,L)
           enddo
          Enddo
        ENDIF
                 
	 endif

         do L=1,nodes(M)
	    write(2,'(5f12.6)')nodex(M,L),nodey(M,L),nodez(M,L)
     &		,alpha0(M,L)*180/3.1416,R0(M,L)
	   enddo
       close(2)

	endif

	Enddo ! M

   88 FORMAT (i5)
   89 FORMAT (5f12.5)
        RETURN
        END SUBROUTINE
!######################################################################
      SUBROUTINE PartLoc_Initial
!######################################################################
      use vars
      use imb
      use mpi
      use multidata
      implicit none
	INTEGER :: K,L,N,numIBslv

       call MPI_BCAST(maxnodeIBS,1,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)		!Total # IB points
       call MPI_BCAST(nodes,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)		!# IB points of each body
       call MPI_BCAST(imb_shape,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)		!IB shape of each body
       call MPI_BCAST(imbnumber,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)		!# IB bodies of each body
       call MPI_BCAST(turax,bodynum,MPI_INTEGER,
     &  master,MPI_COMM_WORLD,ierr)		!# axis of rotation for turbines     
       call MPI_BCAST(reddelta,bodynum,MPI_DOUBLE_PRECISION,
     &  master,MPI_COMM_WORLD,ierr) 	!Reduction factor
       call MPI_BCAST(radsin,bodynum,MPI_DOUBLE_PRECISION,
     &  master,MPI_COMM_WORLD,ierr)		!Rotational velocity of each IB body
       call MPI_BCAST(rotating,bodynum,MPI_LOGICAL,
     &  master,MPI_COMM_WORLD,ierr)		!If the body rotates
       call MPI_BCAST(ibturbine,bodynum,MPI_LOGICAL,
     &  master,MPI_COMM_WORLD,ierr)		!If the body is a turbine
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	nIBslv=maxnodeIBS

!Now each processor will store the # of markers assigned (nIBslv)		
      allocate (kmaxU(nIBslv),kmaxV(nIBslv),kmaxW(nIBslv))
      allocate (nodex_loc(nIBslv),U_Beta1_loc(nIBslv),FX1_loc(nIBslv))
      allocate (nodey_loc(nIBslv),U_Beta2_loc(nIBslv),FX2_loc(nIBslv))
      allocate (nodez_loc(nIBslv),U_Beta3_loc(nIBslv),FX3_loc(nIBslv))
      allocate (alpha0_loc(maxnodeIBS),R0_loc(maxnodeIBS))
	allocate (imb_block_loc(maxnodeIBS),lag_bod_loc(maxnodeIBS))
	allocate (imbinblock_loc(num_domains),rott_loc(maxnodeIBS))
      allocate (FX1NF(nIBslv),FX2NF(nIBslv),FX3NF(nIBslv))

	  kmaxU=0 		 ; kmaxV=0 			; kmaxW=0	         		  
	  nodex_loc=0.d0 	 ; nodey_loc=0.d0 	; nodez_loc=0.d0 
	  U_Beta1_loc=0.d0 ; U_Beta2_loc =0.d0 	; U_Beta3_loc =0.d0 
	  FX1_loc =0.d0    ; FX2_loc =0.d0    	; FX3_loc =0.d0 
	  imbinblock_loc=0 ; imb_block_loc=0  	; lag_bod_loc = 0
 	  alpha0_loc =0.d0 ; R0_loc =0.d0 		; rott_loc=2
	  FX1NF=0.D0	 ; FX2NF=0.D0		; FX3NF=0.D0
								
!nxl is the length of the kernel used for the delta functions. 
!nl is the number of neighbours considered
	if (yangcase.eq.2 .or. yangcase.eq.4.or. yangcase.eq.7) then
	 nxl=2.499999d0 	; nl=125
	endif
	if (yangcase.eq.3 .or. yangcase.eq.6) then
	 nxl=1.999999d0  	; nl=64 
	endif
	if (yangcase.eq.1 .or. yangcase.eq.5) then
	 nxl=1.499999d0  	; nl=27
	endif

       allocate (dh1_loc(nIBslv,nl),dh2_loc(nIBslv,nl))
       allocate (dh3_loc(nIBslv,nl))
	 allocate (I_nr_U(nIBslv,nl),J_nr_U(nIBslv,nl))
	 allocate (I_nr_V(nIBslv,nl),J_nr_V(nIBslv,nl))
	 allocate (I_nr_W(nIBslv,nl),J_nr_W(nIBslv,nl))
	 allocate (K_nr_U(nIBslv,nl),K_nr_V(nIBslv,nl))
	 allocate (K_nr_W(nIBslv,nl))

	  dh1_loc=0.d0	; dh2_loc=0.d0	; dh3_loc=0.d0
	  I_nr_U=0		; J_nr_u=0 		; K_nr_U=0
	  I_nr_V=0 		; J_nr_V=0 		; K_nr_V=0
	  I_nr_W=0 		; J_nr_W=0 		; K_nr_W=0
                                         
	RETURN
	END
!######################################################################
      SUBROUTINE IB_previous
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      INTEGER :: K
	Do K=1,bodynum

!Distribute properties of the actuator line turbines
	if(itime==itime_start .and. imb_shape(K)==5 .and. turax(K)==3) then
	     CALL ActuatorLine_Initial(K)			!Turbine via Actuator Line model
	endif

	 IF (imb_shape(K).eq.5 .and. rotating(K)) then
	   if(ibturbine(K)) then
 		      rads(K)=radsin(K)*CTIME			!Prescribed rotational speed

	   else
       	rads(K)=-radsin(K)*3.1416D0/180.D0*DSIN(0.1983*CTIME)	!Pitching airfoil case
	   endif
		 call imb_moved(K) 		 		!In shapes.for
	   if(ibturbine(K) .and. turax(K).eq.1) call imb_moved_shades(K)  !In shapes.for
	 ENDIF
	Enddo	

		Call PartLoc  					!Asign markers to domains and processors 
!           if (itime.eq.itime_start)   	call imb_vel_to_zero		

	IF(itime.eq.itime_start) then				!Generate delta func for steady bodies
	 if(myrank.eq.master)write(6,*)'Delta functions initiating'
	   Call Deltah 
	 if(myrank.eq.master)write(6,*)'Delta functions generated'
	ENDIF

      END SUBROUTINE
!######################################################################
      SUBROUTINE IBM
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      INTEGER :: NF,M,L,N
!call exchange subroutines to fill the ghost cells with Ustar values
	call exchange(11) ; call exchange(22) ; call exchange(33)

	 FX1NF=0.d0;FX2NF=0.d0;FX3NF=0.d0

	DO NF =1,mdfsteps+1	!MDF loops. +1 as the default loop for IB is 0 MDF loops
	  call imb_openmp
	  rott_loc=2 !This is set to avoid to calculate dh again when NF>1
	ENDDO

!		call imb_averaging !Calculate mean IB force values
	if (mod(itime,2) .eq. 0 .and. ibm_out_forces.eq.1) then
		call caldrag	!Output of IB forces
		if(ibturbine(1) .and. turax(1).eq.2) then
		 call imb_FEM		!Structural loadings of all blades
		 call imb_FEM_oneblade !Structural loadings for a single blade divided into sections
	      endif
	endif

      if (mod(itime,n_out).eq.0) call imb_pressure		!returns coordinates and forces of markers every nout 
!update the velocities at the ghost cell values	
	call exchange(11) ; call exchange(22) ; call exchange(33)
      END SUBROUTINE
!######################################################################
      SUBROUTINE PartLoc 
!######################################################################
      use vars
      use imb
      use mpi
      use multidata
      implicit none
      INTEGER 	:: M,L,ii,nxdom,nydom,nzdom,tnm,N,nx,ny,nz
      DOUBLE PRECISION :: lxdom(idom+1),lydom(jdom+1),lzdom(kdom+1)
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::nodex_mas,nodey_mas
	DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::nodez_mas
	INTEGER,ALLOCATABLE,DIMENSION(:)::lag_bod_mas

       IF(myrank.eq.master)THEN
!Variables used by the master than involve all the markers	
	allocate(lag_bod_mas(maxnodeIBS),nodex_mas(maxnodeIBS))
	allocate(nodey_mas(maxnodeIBS),nodez_mas(maxnodeIBS))
	
	lxdom=0 		; lydom=0 		; lzdom=0
	lxdom(1)=xst 	; lydom(1)=yst 	; lzdom(1)=zst
	do N=2,idom+1
	  lxdom(N)=(xcor(N-2,2)-xcor(N-2,1))+lxdom(N-1)
	enddo
	do N=2,jdom+1
	  lydom(N)=(ycor((N-2)*idom,2)-ycor((N-2)*idom,1))+lydom(N-1)
	enddo
	do N=2,kdom+1
	 lzdom(N)=
     &	 (zcor((N-2)*idom*jdom,2)-zcor((N-2)*idom*jdom,1))+lzdom(N-1)
	enddo

	 imbinblock_loc=0 ; imbinblk=0  	!# Points in each block
	 imb_block=0  				!Block id to which every particle belongs
	 ii=0  ;	tnm=0		
	do M=1,bodynum  			!Perform this operation to all IB bodies
	    DO L=1,nodes(M)		!Analyze all IB poins of the body.
		ii=ii+1
	  Do nx=1,idom 
	   if( (nodex(M,L)-1.d-11).gt.lxdom(nx) .and. 
     &		(nodex(M,L)-1.d-11).le.lxdom(nx+1) )THEN
		nxdom=nx-1
	 	GOTO 490
	   endif 	
	  Enddo
490	CONTINUE
	  Do ny=1,jdom 
	   if( (nodey(M,L)-1.d-11).gt.lydom(ny) .and. 
     &		(nodey(M,L)-1.d-11).le.lydom(ny+1) )THEN
		nydom=ny-1
	 	GOTO 491
	   endif 	
	  Enddo
491	CONTINUE
	  Do nz=1,kdom 
	   if( (nodez(M,L)-1.d-11).gt.lzdom(nz) .and. 
     &		(nodez(M,L)-1.d-11).le.lzdom(nz+1) )THEN
		nzdom=nz-1
	 	GOTO 492
	   endif 	
	  Enddo
492	CONTINUE
		imb_block(ii)=idom*jdom*nzdom+idom*nydom+nxdom
		imbinblk(imb_block(ii)+1)=imbinblk(imb_block(ii)+1)+1
	    ENDDO
	enddo 

	do L=1,num_domains		!Check in all the domains
	 tnm=tnm+imbinblk(L)
	 IF (itime.eq.itime_start .AND. imbinblk(L).ne.0)  
     &	   write(6,*)'Dom,#markrs',L-1,imbinblk(L),tnm  
	 imbinblock_loc(L)=imbinblk(L) !New variable for all the other MPI
	enddo
!Warning if some point is not assigned to some domain
	if(tnm.lt.maxnodeIBS) 
     &  write(6,*)'Some Lagrangian are not assigned to a domain!!!CHECK'	


  	ii=0
	 DO M=1,bodynum
	    DO L=1,nodes(M)
	      ii=ii+1
		imb_block_loc(ii)=imb_block(ii)
		lag_bod_loc(ii)=M 
	      nodex_loc(ii)=nodex(M,L) ; nodey_loc(ii)=nodey(M,L)
	      nodez_loc(ii)=nodez(M,L)

	     IF (itime.eq.itime_start) then !THIS IS DONE ONCE
		R0_loc(ii)=R0(M,L) ; alpha0_loc(ii)=alpha0(M,L)
	     ENDIF 
		rott_loc(ii)=1 	!Moving Lagrangian
		IF(.not.rotating(M))rott_loc(ii)=2  	!Static Lagrangian
	  ENDDO
	 Enddo
	ENDIF !master
	
!Now the properties of the markers owned by each processor is scattered
!We can't use BCAST because then we share all the vectors (+memory)
	IF(itime.eq.itime_start) then
	  call MPI_BCAST(lag_bod_loc,maxnodeIBS,MPI_INTEGER,
     &			master,MPI_COMM_WORLD,ierr)	!# of the body to which the Lag is.
        call MPI_BCAST(alpha0_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(R0_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
	ENDIF

        call MPI_BCAST(rott_loc,maxnodeIBS,MPI_INTEGER,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(imbinblock_loc,num_domains,MPI_INTEGER,
     & 		master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(imb_block_loc,maxnodeIBS,MPI_INTEGER,
     & 		master,MPI_COMM_WORLD,ierr)
	
        call MPI_BCAST(nodex_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     & 		master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nodey_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nodez_loc,maxnodeIBS,MPI_DOUBLE_PRECISION,
     &			master,MPI_COMM_WORLD,ierr)

	if(myrank.eq.master) then
	 deallocate(nodex_mas,nodey_mas,nodez_mas)
	endif

      RETURN
      END
!######################################################################
      SUBROUTINE Deltah
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      DOUBLE PRECISION :: dh,dhtotal
      INTEGER :: I,J,L,ib,K,nnnmls
	 nnnmls=0
	IF(nnnmls.eq.0) then

	Do ib=1,nbp  !Loop through all the blocks of one processor
       if (imbinblock_loc(dom_id(ib)+1).ne.0) then !IF THERE ARE NO POINTS IN THE BLOCK
      Do L = 1,maxnodeIBS !investigate all the IB points
      	nl=0 ;dhtotal=0.d0
	IF(imb_block_loc(L).EQ.dom_id(ib)) THEN !If the IB point is not in the present block
	IF(rott_loc(L).eq.2) THEN	!If the Lagrangian is dynamic:exit
!NEIGHBOURS FOR THE U-GRID
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%x(i) .ge.(nodex_loc(L)-nxl*dom(ib)%dx) .AND.
     &     dom(ib)%x(i) .lt.(nodex_loc(L)+nxl*dom(ib)%dx)) THEN
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).ge.(nodey_loc(L)-nxl*dom(ib)%dy) .AND.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)+nxl*dom(ib)%dy)) THEN
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).ge.(nodez_loc(L)-nxl*dom(ib)%dz) .AND.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)+nxl*dom(ib)%dz)) THEN 
!nl indicates the number of the neighbour and dh1 the delta functions value.
	nl=nl+1
 	  dh1_loc(L,nl)=dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase) 
!The index of the neighbours number nl to the Lagrangian L are:    
	  I_nr_U(L,nl)=I ;  J_nr_U(L,nl)=J ;  K_nr_U(L,nl)=K
	  dhtotal=dhtotal+dh1_loc(L,nl)
	  if(dhtotal.ge.0.99999) goto 876
	    ENDIF
            END DO
	   ENDIF
           END DO
	  ENDIF
          END DO
        dh1_loc(L,nl)=dh1_loc(L,nl)/dhtotal
876	continue       
	kmaxU(L)=nl !# of neighbours of the Lagrangian L
!NEIGHBOURS FOR THE V-GRID
      	nl=0 ;dhtotal=0.d0
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).ge.(nodex_loc(L)-nxl*dom(ib)%dx) .AND.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)+nxl*dom(ib)%dx)) THEN
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%y(j) .ge.(nodey_loc(L)-nxl*dom(ib)%dy) .AND.
     &     dom(ib)%y(j) .lt.(nodey_loc(L)+nxl*dom(ib)%dy)) THEN
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).ge.(nodez_loc(L)-nxl*dom(ib)%dz) .AND.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)+nxl*dom(ib)%dz) ) THEN 
	nl=nl+1
 	  dh2_loc(L,nl)=dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  
	 I_nr_V(L,nl)=I ;  J_nr_V(L,nl)=J ;  K_nr_V(L,nl)=K
	 dhtotal=dhtotal+dh2_loc(L,nl)
	  if(dhtotal.ge.0.99999) goto 877
	    ENDIF
            END DO
	   ENDIF
           END DO
	  ENDIF
          END DO
        dh2_loc(L,nl)=dh2_loc(L,nl)/dhtotal
877	continue              
	kmaxV(L)=nl 
!NEIGHBOURS FOR THE W-GRID
      	nl=0 ;dhtotal=0.d0
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).ge.(nodex_loc(L)-nxl*dom(ib)%dx) .AND.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)+nxl*dom(ib)%dx) )  THEN 
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).ge.(nodey_loc(L)-nxl*dom(ib)%dy) .AND.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)+nxl*dom(ib)%dy) )  THEN 
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%z(k) .ge.(nodez_loc(L)-nxl*dom(ib)%dz) .AND.
     &     dom(ib)%z(k) .lt.(nodez_loc(L)+nxl*dom(ib)%dz) )  THEN 
	nl=nl+1
 	  dh3_loc(L,nl)=dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%YC(J),dom(ib)%Z(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  
	 I_nr_W(L,nl)=I ;  J_nr_W(L,nl)=J ;  K_nr_W(L,nl)=K
	 dhtotal=dhtotal+dh3_loc(L,nl)
	  if(dhtotal.ge.0.99999) goto 878
	    ENDIF
            END DO
	   ENDIF
           END DO
	  ENDIF
          END DO
        dh3_loc(L,nl)=dh3_loc(L,nl)/dhtotal
878	 continue                     
	kmaxW(L)=nl 	
	ENDIF
	ENDIF
        Enddo
       ENDIF  
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ELSE !MLS IS USED        
	Do ib=1,nbp  
     	if(myrank.eq.master)	 write(6,*)'Starting MLS in block#',ib
          if (imbinblock_loc(dom_id(ib)+1).eq.0) GOTO 601 !IF THERE ARE NO POINTS IN THE BLOCK
          Do L = 1,maxnodeIBS
	   IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 701 
	   IF(rott_loc(L).ne.2) GOTO 701	!Rotating
	   call ShapeFunction_MLS(1,L,ib)
	   call ShapeFunction_MLS(2,L,ib)
	   call ShapeFunction_MLS(3,L,ib)
701	CONTINUE   
        Enddo
     	if(myrank.eq.master)	 write(6,*)'Ended MLS in block #',ib
       ENDDO	!ib blocks
601	 CONTINUE
	ENDIF	 !Delta functions or MLS

	 call MPI_BARRIER(MPI_COMM_WORLD,IERR)	
!	 if(myrank.eq.master)write(6,*)'Delta functions generated'
      RETURN
      END
!######################################################################
      SUBROUTINE imb_openmp
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      INTEGER :: I,J,K,L,M,ib,nt,nnmls,KK,iii
      DOUBLE PRECISION :: dh,dhtotal
      DOUBLE PRECISION :: PI,UIB_loc,VIB_loc,WIB_loc,A1,A2,R1,R2,fbeta

	U_Beta1_loc=0.d0 ; U_Beta2_loc=0.d0 ; U_Beta3_loc=0.d0	!Pablo
	Do ib=1,nbp

      if (imbinblock_loc(dom_id(ib)+1).ne.0) THEN

	if (maxnodeIBS.le.100) nt = 1
	if (maxnodeIBS.gt.100) nt = OMP_threads

	call OMP_SET_NUM_THREADS(nt)

	nnmls=2
	IF (nnmls.eq.1) then
!$OMP parallel DEFAULT(SHARED)PRIVATE(I,J,K,L,nl)
!$OMP DO SCHEDULE(DYNAMIC,50)
      Do L = 1,maxnodeIBS
	  IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 901 
        IF(rott_loc(L).eq.1 ) then
	  call ShapeFunction_MLS(1,L,ib) !u-velocity
	  call ShapeFunction_MLS(2,L,ib) !v-velocity
	  call ShapeFunction_MLS(3,L,ib) !w-velocity
        ENDIF
       Do nl=1,KmaxU(L)
	  I=I_nr_U(L,nl) ;  J=J_nr_U(L,nl) ;  K=K_nr_U(L,nl)
        U_Beta1_loc(L)=U_Beta1_loc(L)+dom(ib)%USTAR(I,J,K)*dh1_loc(L,nl)      
       Enddo 
       Do nl=1,KmaxV(L)
        I=I_nr_V(L,nl) ;  J=J_nr_V(L,nl) ;  K=K_nr_V(L,nl)
        U_Beta2_loc(L)=U_Beta2_loc(L)+dom(ib)%VSTAR(I,J,K)*dh2_loc(L,nl)
       Enddo
       Do nl=1,KmaxW(L)
	 I=I_nr_W(L,nl) ;  J=J_nr_W(L,nl) ;  K=K_nr_W(L,nl)
        U_Beta3_loc(L)=U_Beta3_loc(L)+dom(ib)%WSTAR(I,J,K)*dh3_loc(L,nl)
       Enddo
901	 CONTINUE   
      Enddo
!$OMP end DO
!$OMP END PARALLEL	

	else
!Using delta functions. Interpolation of U-velocities
!$OMP parallel DEFAULT(SHARED)PRIVATE(I,J,K,L,nl)
!$OMP DO SCHEDULE(DYNAMIC,50)
      Do L = 1,maxnodeIBS
	 nl=0 
	IF(imb_block_loc(L).eq.dom_id(ib)) THEN
	 IF( rott_loc(L).eq.1 )then
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%x(i) .ge.(nodex_loc(L)-nxl*dom(ib)%dx) .and.
     &     dom(ib)%x(i) .lt.(nodex_loc(L)+nxl*dom(ib)%dx)) THEN
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).ge.(nodey_loc(L)-nxl*dom(ib)%dy) .and.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)+nxl*dom(ib)%dy)) THEN
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).ge.(nodez_loc(L)-nxl*dom(ib)%dz) .and.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)+nxl*dom(ib)%dz)) THEN 
	 nl=nl+1
        dh1_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%X(I),dom(ib)%YC(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  

        U_Beta1_loc(L)=U_Beta1_loc(L)+dom(ib)%USTAR(I,J,K)*dh1_loc(L,nl)

        KmaxU(L)=nl ; I_nr_U(L,nl)=I ; J_nr_U(L,nl)=J ; K_nr_U(L,nl)=K    

		if (yangcase.eq.2 .and. nl.ge.125) GOTO 700
		if (yangcase.eq.4 .and. nl.ge.125) GOTO 700
		if (yangcase.eq.7 .and. nl.ge.125) GOTO 700
		if (yangcase.eq.3 .and. nl.ge.64)  GOTO 700
		if (yangcase.eq.6 .and. nl.ge.64)  GOTO 700
		if (yangcase.eq.1 .and. nl.ge.27)  GOTO 700
		if (yangcase.eq.5 .and. nl.ge.27)  GOTO 700
 
	    ENDIF
            END DO
	   ENDIF
           END DO
	  ENDIF
          END DO
       if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'
        
	ELSE
	  Do nl=1,KmaxU(L)
	  I=I_nr_U(L,nl) ;  J=J_nr_U(L,nl) ;  K=K_nr_U(L,nl)
        U_Beta1_loc(L)=U_Beta1_loc(L)+dom(ib)%USTAR(I,J,K)*dh1_loc(L,nl)      !FMA 
  	  Enddo  	 
	 ENDIF
700 	continue
	endif
      Enddo
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC,50)
      Do L = 1,maxnodeIBS
	IF(imb_block_loc(L).eq.dom_id(ib)) THEN
	 nl=0 
	IF( rott_loc(L).eq.1 )then
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).ge.(nodex_loc(L)-nxl*dom(ib)%dx) .and.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)+nxl*dom(ib)%dx)) then
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%y(j) .ge.(nodey_loc(L)-nxl*dom(ib)%dy) .and.
     &     dom(ib)%y(j) .lt.(nodey_loc(L)+nxl*dom(ib)%dy)) then
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).ge.(nodez_loc(L)-nxl*dom(ib)%dz) .and.
     &     dom(ib)%zc(k).lt.(nodez_loc(L)+nxl*dom(ib)%dz)) then 
	 nl=nl+1
        dh2_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%Y(J),dom(ib)%ZC(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase) 		!June 2015
        U_Beta2_loc(L)=U_Beta2_loc(L)+dom(ib)%VSTAR(I,J,K)*dh2_loc(L,nl)         
        KmaxV(L)=nl ; I_nr_V(L,nl)=I ; J_nr_V(L,nl)=J ; K_nr_V(L,nl)=K

	if (yangcase.eq.2 .and. nl.ge.125) GOTO 701
	if (yangcase.eq.4 .and. nl.ge.125) GOTO 701
	if (yangcase.eq.7 .and. nl.ge.125) GOTO 701
	if (yangcase.eq.3 .and. nl.ge.64)  GOTO 701
	if (yangcase.eq.6 .and. nl.ge.64)  GOTO 701
	if (yangcase.eq.1 .and. nl.ge.27)  GOTO 701
	if (yangcase.eq.5 .and. nl.ge.27)  GOTO 701
 
	    ENDIF
            END DO
	   ENDIF
           END DO
	  ENDIF
          END DO
       if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'         
	ELSE
	 Do nl=1,KmaxV(L)
	  I=I_nr_V(L,nl) ;  J=J_nr_V(L,nl) ;  K=K_nr_V(L,nl)
        U_Beta2_loc(L)=U_Beta2_loc(L)+dom(ib)%VSTAR(I,J,K)*dh2_loc(L,nl)
  	 Enddo
	ENDIF
701	continue
	endif
      Enddo
!$OMP END DO
!$OMP DO SCHEDULE(DYNAMIC,50)
      Do L = 1,maxnodeIBS
	IF(imb_block_loc(L).eq.dom_id(ib)) THEN
	 nl=0 
	IF( rott_loc(L).eq.1 )then
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).ge.(nodex_loc(L)-nxl*dom(ib)%dx) .and.
     &     dom(ib)%xc(i).lt.(nodex_loc(L)+nxl*dom(ib)%dx)) then 
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).ge.(nodey_loc(L)-nxl*dom(ib)%dy) .and.
     &     dom(ib)%yc(j).lt.(nodey_loc(L)+nxl*dom(ib)%dy)) then
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%z(k) .ge.(nodez_loc(L)-nxl*dom(ib)%dz) .and.
     &     dom(ib)%z(k) .lt.(nodez_loc(L)+nxl*dom(ib)%dz)) then
	 nl=nl+1
        dh3_loc(L,nl)= dh(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz,
     &  dom(ib)%XC(I),dom(ib)%YC(J),dom(ib)%Z(K)
     & ,nodex_loc(L),nodey_loc(L),nodez_loc(L),yangcase)  		!June 2015
        U_Beta3_loc(L)=U_Beta3_loc(L)+dom(ib)%WSTAR(I,J,K)*dh3_loc(L,nl)
        KmaxW(L)=nl ;  I_nr_W(L,nl)=I ; J_nr_W(L,nl)=J ; K_nr_W(L,nl)=K

	if (yangcase.eq.2 .and. nl.ge.125) GOTO 702		!March 2016
	if (yangcase.eq.4 .and. nl.ge.125) GOTO 702
	if (yangcase.eq.7 .and. nl.ge.125) GOTO 702
	if (yangcase.eq.3 .and. nl.ge.64)  GOTO 702
	if (yangcase.eq.6 .and. nl.ge.64)  GOTO 702
	if (yangcase.eq.1 .and. nl.ge.27)  GOTO 702
	if (yangcase.eq.5 .and. nl.ge.27)  GOTO 702
 
	    ENDIF
            END DO
	   ENDIF
           END DO
	  ENDIF
          END DO
        if (nl.eq.0) write(6,*)L,'nl is equal to 0!!'
       ELSE
	  Do nl=1,KmaxW(L)
	  I=I_nr_W(L,nl) ;  J=J_nr_W(L,nl) ;  K=K_nr_W(L,nl)
        U_Beta3_loc(L)=U_Beta3_loc(L)+dom(ib)%WSTAR(I,J,K)*dh3_loc(L,nl)
  	  Enddo
	 ENDIF
702	continue
	endif
      Enddo
!$OMP end DO
!$OMP END PARALLEL	
	  endif
	 ENDIF
	Enddo !ib-loop
!#################   SUBROUTINE calfl   #################################
	 FX1_loc = 0.d0     ; FX2_loc = 0.d0	; FX3_loc = 0.d0
 
	DO ib=1,nbp
	IF (imbinblock_loc(dom_id(ib)+1).NE.0) THEN !No points within the block

!$OMP parallel DEFAULT(SHARED)PRIVATE(L)
!$OMP DO SCHEDULE(DYNAMIC,50)
        Do L = 1,maxnodeIBS
	  IF(imb_block_loc(L).EQ.dom_id(ib)) THEN              !REMOVE GOTO
	   UIB_loc = 0.d0; VIB_loc = 0.d0; WIB_loc = 0.d0
		M=lag_bod_loc(L) 

 	  IF (imb_shape(M).eq.5 .and. turax(M).eq.3) THEN	!Pablo 09/2017
		call ActuatorLine(M,L,ib)
	  ELSE
	   IF(ibturbine(M)) then 			!If it is a turbine
	    IF (turax(M).eq.1) then			! Vertical Axis Turbine
	      UIB_loc=-radsin(M)*R0_loc(L)*dcos(rads(M)+alpha0_loc(L))
	      VIB_loc=-radsin(M)*R0_loc(L)*dsin(rads(M)+alpha0_loc(L))
	      WIB_loc= 0.d0	
	    ENDIF	     
	    IF (turax(M).eq.2) then 		! Horizontal Axis Turbine
	      UIB_loc= 0.d0			
	      VIB_loc= radsin(M)*R0_loc(L)*dcos(rads(M)+alpha0_loc(L))
	      WIB_loc=-radsin(M)*R0_loc(L)*dsin(rads(M)+alpha0_loc(L))
	    ENDIF	     
	   ENDIF
	  IF (.not.ibturbine(M).AND. imb_shape(M)==5.and.rotating(M)) then	! from file but moving
!	      UIB_loc=-radsin(M)*R0_loc(L)*dcos(rads(M)+alpha0_loc(L)) !Pitching airfoil simu
!	      VIB_loc=-radsin(M)*R0_loc(L)*dsin(rads(M)+alpha0_loc(L))
!	      WIB_loc= 0.d0	
!		A1=-radsin(M)*3.1416D0/180.D0*DSIN(0.1983*(CTIME-dt))
!		A2=-radsin(M)*3.1416D0/180.D0*DSIN(0.1983*CTIME)
		R1=R0_loc(L)*DSIN(A1+alpha0_loc(L))
		R2=R0_loc(L)*DSIN(A2+alpha0_loc(L))
		UIB_loc=(R2-R1)/dt

		R1=R0_loc(L)*DCOS(A1+alpha0_loc(L))
		R2=R0_loc(L)*DCOS(A2+alpha0_loc(L))
		VIB_loc=(R2-R1)/dt
		WIB_loc=0.D0
	    ENDIF

        FX1_loc(L)=UIB_loc-U_Beta1_loc(L)
        FX2_loc(L)=VIB_loc-U_Beta2_loc(L)
        FX3_loc(L)=WIB_loc-U_Beta3_loc(L)  
	
	endif!actuator line

	   ENDIF
	  ENDDO  
!$OMP end DO
!$OMP END PARALLEL
                         
	ENDIF
       Enddo !ib-loop	

!Accumulate the forces as MDF can be performed
	   Do L=1,maxnodeIBS
	    FX1NF(L)=FX1_loc(L)+FX1NF(L)
	    FX2NF(L)=FX2_loc(L)+FX2NF(L)
	    FX3NF(L)=FX3_loc(L)+FX3NF(L)
	   ENDDO
!!################   SUBROUTINE distfbeta   ###################
	Do ib=1,nbp
       if(imbinblock_loc(dom_id(ib)+1).NE.0) THEN
!$OMP parallel DEFAULT(SHARED)PRIVATE(L)
!$OMP DO SCHEDULE(DYNAMIC,50)
      Do L = 1,maxnodeIBS
	IF(imb_block_loc(L).EQ.dom_id(ib)) THEN      !REMOVE GOTO
	 Do nl=1,KmaxU(L)
          I=I_nr_U(L,nl) ;  J=J_nr_U(L,nl) ;  K=K_nr_U(L,nl)
          fbeta = FX1_loc(L)*dh1_loc(L,nl)*reddelta(lag_bod_loc(L))    
          dom(ib)%USTAR(I,J,K) = dom(ib)%USTAR(I,J,K) + fbeta 
  	 Enddo	
	 Do nl=1,KmaxV(L)
	   I=I_nr_V(L,nl) ;  J=J_nr_V(L,nl);  K=K_nr_V(L,nl)
          fbeta = FX2_loc(L)*dh2_loc(L,nl)*reddelta(lag_bod_loc(L))   
          dom(ib)%VSTAR(I,J,K) = dom(ib)%VSTAR(I,J,K) + fbeta
  	 Enddo
	 Do nl=1,KmaxW(L)
	   I=I_nr_W(L,nl) ;  J=J_nr_W(L,nl);  K=K_nr_W(L,nl)
          fbeta = FX3_loc(L)*dh3_loc(L,nl)*reddelta(lag_bod_loc(L))      
          dom(ib)%WSTAR(I,J,K) = dom(ib)%WSTAR(I,J,K) + fbeta
  	 Enddo
	ENDIF
      End do 
!$OMP end DO
!$OMP END PARALLEL
	ENDIF
	Enddo !ib-loop 
	
      RETURN
      END
!######################################################################
      SUBROUTINE caldrag
!######################################################################
	use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER :: N,M,L,j,iii,inipts,finpts,totalpoints
      DOUBLE PRECISION::PI,fx_loc,fy_loc,fz_loc,alpharads
      DOUBLE PRECISION::ft_loc,fn_loc,ft2_loc,fn2_loc,M_loc
      DOUBLE PRECISION::sumvel,l1normU,l1normV
	DOUBLE PRECISION,allocatable,dimension(:)::FX1_mas,FX2_mas
	DOUBLE PRECISION,allocatable,dimension(:)::FX3_mas
      PI = 4.D0*DATAN(1.D0)
      allocate (FX1_mas(maxnodeIBS),FX2_mas(maxnodeIBS))
      allocate (FX3_mas(maxnodeIBS))
	FX1_mas =0.d0 ; FX2_mas =0.d0 ; FX3_mas =0.d0 
	FX1_loc =0.d0 ; FX2_loc =0.d0 ; FX3_loc =0.d0 
        
! To gather the forces ordered by the marker global index L
! the forces calculated by the masters are sent to the MASTER PROC
	   Do L=1,maxnodeIBS
	    FX1_loc(L)=FX1NF(L) ;FX2_loc(L)=FX2NF(L) ; FX3_loc(L)=FX3NF(L)
	   ENDDO
        call MPI_ALLREDUCE (FX1_loc,FX1_mas,maxnodeIBS,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE (FX2_loc,FX2_mas,maxnodeIBS,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE (FX3_loc,FX3_mas,maxnodeIBS,
     &      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Now the forces are transmitted, the resultants are calcualted by the MASTER
	IF (Myrank.EQ.master) THEN
	 FX1=0.d0;FX2=0.d0;FX3=0.d0 
	N=0 
	Do M=1,bodynum
	 Do L=1,nodes(M)
	  N=N+1
	  FX1(M,L)=FX1_mas(N) ; FX2(M,L)=FX2_mas(N) ; FX3(M,L)=FX3_mas(N)
	 enddo
	enddo

	J=0
	Do M = 1,bodynum
	IF (imb_shape(M).eq.5 .and. rotating(M)) then
	  Do iii=1,imbnumber(M)
		J=J+1 ; forcefilej=399+J    
 	  fx_loc = 0.d0   ; fy_loc = 0.d0 ; fz_loc = 0.d0 
	  ft_loc = 0.d0	; fn_loc = 0.d0 ; ft2_loc = 0.d0 ; fn2_loc = 0.d0
	  totalpoints=nodes(M)/imbnumber(M)
	  inipts=(iii-1)*totalpoints+1 
	  finpts=iii*totalpoints
	  alpharads= (rads(M)+(iii-1)*2.*PI/imbnumber(M))*180.d0/PI
!!!!!
	if (turax(M).eq.2) then !HATT
	 do L=1,nodes(M)
	  fx_loc=fx_loc+FX1(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt	
	  fy_loc=fy_loc+FX2(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt 	
	  fz_loc=fz_loc+FX3(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt
	  ft_loc=ft_loc+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))-
     &   	FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)/dt
     &  	*dxm*dym*dzm*reddelta(lag_bod_loc(L)) 	
	 end do
       write(forcefilej,88)alpharads,fx_loc,fy_loc,fz_loc,ft_loc
!!!!!
	ELSE IF (turax(M).eq.1) then 	!VATT
	M_loc=0.d0
	 do L=inipts,finpts
	  fx_loc=fx_loc+FX1(M,L)*dxm*dym*dzm/dt!*reddelta(lag_bod_loc(L))
	  fy_loc=fy_loc+FX2(M,L)*dxm*dym*dzm/dt  	
!	  fz_loc = fz_loc+FX3(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt
	  ft_loc=ft_loc+(FX1(M,L)*dcos(rads(M)+alpha0(M,L))+
     &     FX2(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)/dt*dxm*dym*dzm
	  ft2_loc=ft2_loc+(FX1(M,L)*dcos(rads(M)+alpha0(M,L))-
     &     FX2(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)/dt*dxm*dym*dzm
	  fn_loc=fn_loc+(FX1(M,L)*dsin(rads(M)+alpha0(M,L))+
     &     FX2(M,L)*dcos(rads(M)+alpha0(M,L)))/dt*dxm*dym*dzm
	  fn2_loc=fn2_loc+(FX1(M,L)*dsin(rads(M)+alpha0(M,L))-
     &     FX2(M,L)*dcos(rads(M)+alpha0(M,L)))/dt*dxm*dym*dzm
	   M_loc=M_loc+
     &     (FX1(M,L)*(R0(M,L)*dcos(alpha0(M,L))-
     &                R(M)*dcos(rads(M)+(iii-1)*2.*PI/imbnumber(M)))+
     &      FX2(M,L)*(R0(M,L)*dsin(alpha0(M,L))-
     &                R(M)*dsin(rads(M)+(iii-1)*2.*PI/imbnumber(M))))
     &	/dt*dxm*dym*dzm
	 end do

         write(forcefilej,88)alpharads,fx_loc,fy_loc
     &	,ft_loc,ft2_loc,fn_loc,fn2_loc,M_loc
	ENDIF

	  Enddo !iii-loop

!!!!! Actuator line - Pablo 09/2017
	IF (turax(M).eq.3) then 	!AL
		fx_loc = 0.d0   ; ft_loc = 0.d0 ;forcefilej=400   
	  	alpharads=rads(M)*180.d0/PI 
	 do L=1,nodes(M)
	  fx_loc=fx_loc+ FX1(M,L)*dzm*c_act(L)	
	  ft_loc=ft_loc+(FX2(M,L)*dcos(rads(M)+alpha0(M,L))
     &         -FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)*dzm*c_act(L)
	 end do
	call ActuatorLine_FEM(M)
!For the AL plot the coefficients of thrust and power:
       write(forcefilej,88)alpharads
     &  ,2.d0*fx_loc/(PI*0.135**2*ubulk**2)
     &  ,2.d0*ft_loc*radsin(M)/(PI*0.135**2*ubulk**3)
     &  ,fx_loc,ft_loc*radsin(M)
       write(6,88)alpharads,2.d0*fx_loc/(PI*0.135**2*ubulk**2)
     &  ,2.d0*ft_loc*radsin(M)/(PI*0.135**2*ubulk**3)
	ENDIF
!!!!!
	ELSE
  	  fx_loc = 0.d0   ; fy_loc = 0.d0 ; fz_loc = 0.d0 
	  J=J+1 ;	 forcefilej=399+J    
	  do L=1,nodes(M)	
	   fx_loc=fx_loc+FX1(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt	
	   fy_loc=fy_loc+FX2(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt 	
	   fz_loc=fz_loc+FX3(M,L)*dxm*dym*dzm*reddelta(lag_bod_loc(L))/dt
	  end do
         write(forcefilej,88) CTIME,fx_loc,fy_loc,fz_loc
	ENDIF

       End do !M loop

!l2-norm is calculated in reference to the final velocitiy field
!	sumvel=0.d0 ; l1normU=0.d0 ; l1normV=0.d0 ; l2norm=0.d0
!	DO M=1,bodynum
!	 Do L=1,nodes(M)
!    	   sumvel=sumvel+(FX1(M,L)**2+FX2(M,L)**2+FX3(M,L)**2)
!	   l1normU=l1normU+(ABS(FX1(M,L)))!+ABS(FX2(M,L)+ABS(FX3(M,L))
!	   l1normV=l1normV+(ABS(FX2(M,L)))!+ABS(FX2(M,L)+ABS(FX3(M,L))
!	 End Do
!	 l2norm(M)=DSQRT(sumvel/nodes(M)) 
!	 l1normU=l1normU/nodes(M)	;  l1normV=l1normV/nodes(M)
!	End do
!	 if(bodynum.eq.1) write(757,'(4F20.7)')
!     &       CTIME,l2norm(1),l1normU,l1normV
!	 if(bodynum.ge.2) write(757,'(3F20.7)')CTIME,l2norm(1),l2norm(2)

	ENDIF !master

	deallocate (FX1_mas,FX2_mas,FX3_mas)

   88 FORMAT (14f15.7)
   98 FORMAT (14E15.7)
      RETURN
      END SUBROUTINE caldrag
!#############################################################
      SUBROUTINE imb_pressure
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER      :: ib,M,L,I,J,K,strlen,strlen2,Geom_Time1    
      CHARACTER*8  :: char_block,ibnum
      CHARACTER*31 :: gridfile1

	IF (myrank.eq.master) then

	Do M=1,bodynum
	
	IF (rotating(M)) then
         write(ibnum,'(I2)') M
           strlen2=LEN(TRIM(ADJUSTL(ibnum)))
           ibnum=REPEAT('0',(2-strlen2))//TRIM(ADJUSTL(ibnum))
           write(char_block,'(I6)') itime
           strlen=LEN(TRIM(ADJUSTL(char_block)))
           char_block=REPEAT('0',(6-strlen))//TRIM(ADJUSTL(char_block))
           gridfile1='PnFint_'//TRIM(ADJUSTL(char_block))//
     &	   '_'//TRIM(ADJUSTL(ibnum))//'.plt'
         open (unit=Geom_Time1, file=gridfile1)
	   write(Geom_Time1,*)'TITLE = "Forces at IB body"'
	   write(Geom_Time1,*)'VARIABLES=x,y,z,Fx,Fy,Fz'
	   write(Geom_Time1,*)'zone i= ',nodes(M),'DATAPACKING = POINT'
	  do L=1,nodes(M)
	   write (Geom_Time1,89) nodex(M,L),nodey(M,L),nodez(M,L),
     &  FX1(M,L)*dxm*dym*dzm*reddelta(M)/dt*1000.d0,
     &  FX2(M,L)*dxm*dym*dzm*reddelta(M)/dt*1000.d0,
     &  FX3(M,L)*dxm*dym*dzm*reddelta(M)/dt*1000.d0
	  enddo
	 close(Geom_Time1)
	
	ENDIF !rotating

	Enddo !M

	ENDIF

   89 FORMAT (17e14.5)
        RETURN
        END SUBROUTINE	
!#############################################################
      SUBROUTINE imb_FEM
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      INTEGER          :: M,L,K,numblades  
      DOUBLE PRECISION :: dh,alpharads,PI,x1,x2,x3,Rhub
      DOUBLE PRECISION :: Mxtot(3),Mqtot(3),Mzytot(3),Mzxtot(3)
	DOUBLE PRECISION,allocatable,dimension(:)::Mx,Mq,Tx,Tq
	DOUBLE PRECISION,allocatable,dimension(:)::R0_x,R0_y,R0_z

	IF (myrank.eq.master) then
      PI = 4.D0*DATAN(1.D0)
	Do M=1,1!bodynum
	 numblades=3
	  if(itime.eq.itime_start+1) then
         open (unit=1414, file='FEM_1.dat')
	   write(1414,*)'TITLE = "Forces at IB. Fx(kN), Mom (Nm)"'
	   write(1414,*)'VARIABLES=Deg,X1,X2,X3',
     &  ',BR1,BR2,BR3,BT1,BT2,BT3,TX1,TX2,TX3,TQ1,TQ2,TQ3'
	  endif

!Distances of each marker
	 allocate(R0_y(nodes(M)),R0_z(nodes(M)),R0_x(nodes(M)))
	   Rhub=0.0475d0 ; R0_x=0.d0 ; R0_y=0.d0 ; R0_z=0.d0 ;
	   do L=1,nodes(M)/numblades
!		R0_z(L)= ABS(R0(M,L)*cos(alpha0(M,L)))-Rhub	!Always>0 
		R0_z(L)= R0(M,L)-Rhub	!Pablo 07_2017
		R0_y(L)= R0(M,L)*sin(alpha0(M,L))
		R0_x(L)=nodexlocal(M,L)-Rhub-0.035*0.25-0.0048525

	 if(R0_z(L).le.0.d0) write(6,*)'FEM Ro_z',L,'Very small',R0_z(L)
	   enddo
	   do L=1,nodes(M)/numblades
		R0_x(L+nodes(M)/numblades)=R0_x(L) 
		R0_x(L+2*nodes(M)/numblades)=R0_x(L)
		R0_y(L+nodes(M)/numblades)=R0_y(L) 
		R0_y(L+2*nodes(M)/numblades)=R0_y(L)
		R0_z(L+nodes(M)/numblades)=R0_z(L) 
		R0_z(L+2*nodes(M)/numblades)=R0_z(L)	
	   enddo

	 allocate(Mx(nodes(M)),Tx(nodes(M)),Tq(nodes(M)))
	 allocate(Mq(nodes(M)))

	 Mx=0.d0 ; Mq=0.d0 ; Tx=0.d0 ; Tq=0.d0 
	 Mxtot=0.d0 ; Mqtot=0.d0 ; Mzytot=0.d0 ; Mzxtot=0.d0
	 x1=0.d0 ; x2=0.d0 ; x3=0.d0 
!Thrust produced at each blade:
	  do L=1,nodes(M)/numblades
	   x1=x1+FX1(M,L)*dxm*dym*dzm/dt*1000.d0
	   x2=x2+FX1(M,L+1*nodes(M)/numblades)*dxm*dym*dzm/dt*1000.d0
	   x3=x3+FX1(M,L+2*nodes(M)/numblades)*dxm*dym*dzm/dt*1000.d0
	  enddo

!Bending moments:
!Flapwise bending moment:
	  do L=1,nodes(M)
	   Mx(L)=FX1(M,L)*R0_z(L)*dxm*dym*dzm*reddelta(M)/dt*1000.d0
	  enddo
	  do k=1,numblades
	    do L=1,nodes(M)/numblades
	 	Mxtot(k)=Mxtot(k)+Mx(L+nodes(M)/numblades*(k-1))
	    enddo
	  enddo

!Edgewise or radial bending moment:
	  do L=1,nodes(M)
	   Mq(L)=(FX2(M,L)*dcos(rads(M)+alpha0(M,L))-
     &   	FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0_z(L)
     &      *dxm*dym*dzm*reddelta(M)/dt*1000.d0
	  enddo
	  do k=1,numblades
	    do L=1,nodes(M)/numblades
		Mqtot(k)=Mqtot(k)+Mq(L+nodes(M)/numblades*(k-1))
	    enddo
	  enddo

!Torsional moments:
	  do L=1,nodes(M)
	   Tx(L)=(FX1(M,L)*dxm*dym*dzm*reddelta(M)/dt)*R0_y(L)*1000.d0
	   Tq(L)=(FX2(M,L)*dcos(rads(M)+alpha0(M,L))-
     &   	    FX3(M,L)*dsin(rads(M)+alpha0(M,L)))
     &          *dxm*dym*dzm*reddelta(M)/dt*R0_x(L)*1000.d0
	  enddo
	 do k=1,numblades
	    do L=1,nodes(M)/numblades
	 	Mzxtot(k)=Mzxtot(k)+Tx(L+nodes(M)/numblades*(k-1))
		Mzytot(k)=Mzytot(k)+Tq(L+nodes(M)/numblades*(k-1))
	    enddo
	 enddo

	  alpharads=rads(M)*180.d0/PI

	  write (1414,89)alpharads,x1,x2,x3,
     &  Mxtot(1),Mxtot(2),Mxtot(3),Mqtot(1),Mqtot(2),Mqtot(3),
     &  Mzxtot(1),Mzxtot(2),Mzxtot(3),Mzytot(1),Mzytot(2),Mzytot(3)
	
	 deallocate(Mx,Mq,Tx,Tq,R0_x,R0_y,R0_z)
	  
	Enddo !M
	ENDIF

   89 FORMAT (16e15.5)
        RETURN
        END SUBROUTINE	
!#############################################################
      SUBROUTINE imb_FEM_oneblade
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
	INTEGER,parameter	:: numsections=5,numblades=3
      INTEGER          :: M,L,K 
      DOUBLE PRECISION :: alpharads,Rhub,PI
      DOUBLE PRECISION :: fx,fx_tot,ft,ft_tot,Mx,Mq,Mx_tot,Mq_tot
	DOUBLE PRECISION,dimension(numsections)::Mxtot,Mqtot
	DOUBLE PRECISION,dimension(numsections)::Xtot,Qtot

	IF (myrank.eq.master) then
      PI = 4.D0*DATAN(1.D0)

	 Do M=1,1!bodynum

	 if(itime.eq.itime_start+1) then
	  if (M.eq.1) then
         open (unit=1414, file='FEM_1.dat')
	   write(1414,*)'TITLE = "Forces at IB#1. Fx(kN), Mom (Nm)"'
	   write(1414,*)'VARIABLES=Deg,T1,T2,T3,T4,T5,Q1,Q2,Q3,Q4,Q5,',
     &  'MX1,MX2,MX3,MX4,MX5,MQ1,MQ2,MQ3,MQ4,MQ5'
	  else
         open (unit=1415, file='FEM_2.dat')
	   write(1415,*)'TITLE = "Forces at IB#2. Fx(kN), Mom (Nm)"'
	   write(1415,*)'VARIABLES=Deg,T1,T2,T3,T4,T5,Q1,Q2,Q3,Q4,Q5,',
     &  'MX1,MX2,MX3,MX4,MX5,MQ1,MQ2,MQ3,MQ4,MQ5'
	  endif
	 endif
		Rhub=0.0500d0 

	 Mx=0.d0 ; Mq=0.d0 ; Mx_tot=0.d0 ; Mq_tot=0.d0
       fx=0.d0 ; ft=0.d0 ; fx_tot=0.d0 ; ft_tot=0.d0
	 Mxtot=0.d0 ; Mqtot=0.d0 ; Xtot=0.d0 ; Qtot=0.d0

!Thrust produced at each blade:
	 do L=1,nodes(M)/numblades
         fx=FX1(M,L)*dxm*dym*dzm*reddelta(1)/dt*1000.d0
         fx_tot=fx_tot+fx
           if (L.ge.1 .and. L.lt.nodes(M)/numblades/numsections)
     &     Xtot(1)=Xtot(1)+fx

           if     (L.ge.nodes(M)/numblades*1/numsections 
     &  .and. L.lt.nodes(M)/numblades*2/numsections) Xtot(2)=Xtot(2)+fx

           if     (L.ge.nodes(M)/numblades*2/numsections 
     &  .and. L.lt.nodes(M)/numblades*3/numsections) Xtot(3)=Xtot(3)+fx

           if     (L.ge.nodes(M)/numblades*3/numsections 
     &  .and. L.lt.nodes(M)/numblades*4/numsections) Xtot(4)=Xtot(4)+fx

           if     (L.ge.nodes(M)/numblades*4/numsections 
     &  .and. L.le.nodes(M)/numblades*5/numsections) Xtot(5)=Xtot(5)+fx

	 enddo

!Torque produced at each blade:
         do L=1,nodes(M)/numblades
          ft=(FX2(M,L)*dcos(rads(M)+alpha0(M,L))-
     &       FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*R0(M,L)/dt
     &       *dxm*dym*dzm*reddelta(lag_bod_loc(L))*1000.d0   
          ft_tot=ft_tot+ft
           if (L.ge.1 .and. L.lt.nodes(M)/numblades/numsections)
     &     Qtot(1)=Qtot(1)+ft

           if  (L.ge.nodes(M)/numblades*1/numsections 
     &  .and. L.lt.nodes(M)/numblades*2/numsections)  Qtot(2)=Qtot(2)+ft

           if  (L.ge.nodes(M)/numblades*2/numsections 
     &  .and. L.lt.nodes(M)/numblades*3/numsections)  Qtot(3)=Qtot(3)+ft

           if  (L.ge.nodes(M)/numblades*3/numsections 
     &  .and. L.lt.nodes(M)/numblades*4/numsections)  Qtot(4)=Qtot(4)+ft

           if  (L.ge.nodes(M)/numblades*4/numsections 
     &  .and. L.le.nodes(M)/numblades*5/numsections)  Qtot(5)=Qtot(5)+ft
         enddo


!Flapwise bending moment:
         do L=1,nodes(M)/numblades
	    Mx=FX1(M,L)*(R0(M,L)-Rhub)*dxm*dym*dzm*reddelta(M)/dt*1000.d0
          Mx_tot=Mx_tot+Mx

           if (L.ge.1 .and. L.lt.nodes(M)/numblades/numsections)
     &    Mxtot(1)=Mxtot(1)+Mx

           if  (L.ge.nodes(M)/numblades*1/numsections 
     & .and. L.lt.nodes(M)/numblades*2/numsections) Mxtot(2)=Mxtot(2)+Mx

           if  (L.ge.nodes(M)/numblades*2/numsections 
     & .and. L.lt.nodes(M)/numblades*3/numsections) Mxtot(3)=Mxtot(3)+Mx

           if  (L.ge.nodes(M)/numblades*3/numsections 
     & .and. L.lt.nodes(M)/numblades*4/numsections) Mxtot(4)=Mxtot(4)+Mx

           if  (L.ge.nodes(M)/numblades*4/numsections 
     & .and. L.le.nodes(M)/numblades*5/numsections) Mxtot(5)=Mxtot(5)+Mx
         enddo

!Edgewise or radial bending moment:
         do L=1,nodes(M)/numblades
	    Mq=(FX2(M,L)*dcos(rads(M)+alpha0(M,L))-
     &   	  FX3(M,L)*dsin(rads(M)+alpha0(M,L)))*(R0(M,L)-Rhub)
     &        *dxm*dym*dzm*reddelta(M)/dt*1000.d0
          Mq_tot=Mq_tot+Mq

           if (L.ge.1 .and. L.lt.nodes(M)/numblades/numsections)
     &    Mqtot(1)=Mqtot(1)+Mq

           if  (L.ge.nodes(M)/numblades*1/numsections 
     &  .and. L.lt.nodes(M)/numblades*2/numsections) Mqtot(2)=Mqtot(2)+Mq

           if  (L.ge.nodes(M)/numblades*2/numsections 
     & .and. L.lt.nodes(M)/numblades*3/numsections) Mqtot(3)=Mqtot(3)+Mq

           if  (L.ge.nodes(M)/numblades*3/numsections 
     & .and. L.lt.nodes(M)/numblades*4/numsections) Mqtot(4)=Mqtot(4)+Mq

           if  (L.ge.nodes(M)/numblades*4/numsections 
     & .and. L.le.nodes(M)/numblades*5/numsections) Mqtot(5)=Mqtot(5)+Mq
         enddo

	  alpharads=DABS(rads(M))*180.d0/PI

	 if (M.eq.1) then
	  write (1414,89)alpharads,Xtot(:),Qtot(:),Mxtot(:),Mqtot(:)
	 else
	  write (1415,89)alpharads,Xtot(:),Qtot(:),Mxtot(:),Mqtot(:)
	 endif

	Enddo !M
	ENDIF

   89 FORMAT (1f15.5,25e15.6)
        RETURN
        END SUBROUTINE	
!######################################################################
      SUBROUTINE imb_vel_to_zero
!######################################################################
      use vars
      use mpi
      use multidata
      use imb
      implicit none
      DOUBLE PRECISION :: dh,umin,xmax,xmin,ymax,ymin,zmax,zmin
      INTEGER :: I,J,K,L,ib

	xmax=0.d0;  xmin=10000.d0  ;ymax=0.d0;  ymin=10000.d0   
	zmax=0.d0;  zmin=10000.d0   

	umin=1.d-0
	
	Do ib=1,nbp  !Loop through all the blocks of one processor

       if (imbinblock_loc(dom_id(ib)+1).eq.0) GOTO 600 !IF THERE ARE NO POINTS IN THE BLOCK

      Do L = 1,maxnodeIBS !investigate all the IB points
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 700 !If the IB point is not in the present block
	xmax=max(xmax,nodex_loc(L)) ;	xmin=min(xmin,nodex_loc(L))
	ymax=max(ymax,nodey_loc(L)) ;	ymin=min(ymin,nodey_loc(L))
	zmax=max(zmax,nodez_loc(L)) ;	zmin=min(zmin,nodez_loc(L))
700 	CONTINUE
       Enddo !L

      Do L = 1,maxnodeIBS !investigate all the IB points
	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 705 !If the IB point is not in the present block
!NEIGHBOURS FOR THE U-GRID
          DO I = 1, dom(ib)%ttc_i 
       IF (dom(ib)%xc(i).ge.xmin .and. dom(ib)%xc(i) .le.xmax) THEN
           DO J = 1, dom(ib)%ttc_j 
       IF (dom(ib)%yc(j).ge.ymin .and. dom(ib)%yc(j) .le.ymax) THEN
            DO K = 1, dom(ib)%ttc_k 
       IF (dom(ib)%zc(k).ge.zmin .and. dom(ib)%zc(k) .le.zmax) THEN
!nl indicates the number of the neighbour and dh1 the delta functions value.
!The index of the neighbours number nl to the Lagrangian L are:    
	   dom(ib)%U(I,J,K)=umin ; dom(ib)%V(I,J,K)=umin
	   dom(ib)%W(I,J,K)=umin
	ENDIF
            END DO
	ENDIF
           END DO
	ENDIF
          END DO
705 	CONTINUE
       Enddo !L

!      Do L = 1,maxnodeIBS !investigate all the IB points
!	IF(imb_block_loc(L).ne.dom_id(ib)) GOTO 700 !If the IB point is not in the present block
!NEIGHBOURS FOR THE U-GRID
!          DO I = 1, dom(ib)%ttc_i 
!       IF (dom(ib)%xc(i) .gt.(nodex_loc(L)+nxl*dom(ib)%dx) .or.
!     &     dom(ib)%xc(i) .lt.(nodex_loc(L)-nxl*dom(ib)%dx)) GOTO 210
!           DO J = 1, dom(ib)%ttc_j 
!       IF (dom(ib)%yc(j).gt.(nodey_loc(L)+nxl*dom(ib)%dy) .or.
!     &     dom(ib)%yc(j).lt.(nodey_loc(L)-nxl*dom(ib)%dy)) GOTO 211
!            DO K = 1, dom(ib)%ttc_k 
!       IF (dom(ib)%zc(k).gt.(nodez_loc(L)+nxl*dom(ib)%dz) .or.
!     &     dom(ib)%zc(k).lt.(nodez_loc(L)-nxl*dom(ib)%dz)) GOTO 212 
!nl indicates the number of the neighbour and dh1 the delta functions value.
!The index of the neighbours number nl to the Lagrangian L are:    
!	   dom(ib)%U(I,J,K)=umin
!	   dom(ib)%V(I,J,K)=umin
!	   dom(ib)%U(I,J,K)=umin
!212         CONTINUE
!            END DO
!211         CONTINUE
!           END DO
!210         CONTINUE
!          END DO
!700 	CONTINUE
!       Enddo !L
600	CONTINUE   
	Enddo !ib
      RETURN
      END
!#############################################################
      SUBROUTINE imb_averaging
!#############################################################
      use vars
      use multidata
      use imb
      use mpi
      implicit none
      double precision    :: facp1,facm1,facp2,facm2
      INTEGER      :: ib,M,L,I,J,K,strlen,strlen2,Geom_Time1    
      CHARACTER*8  :: char_block,ibnum
      CHARACTER*31 :: gridfile1

	IF (myrank.eq.master) then

	
!.....For first order moments
        if (ctime.ge.t_start_averaging1) then
           facp1 = 1./( ((ctime-t_start_averaging1)/dt) + 1)
        else
           facp1 = 1.
        endif
           facm1 = 1. - facp1
!.....For second order moments
        if (ctime.ge.t_start_averaging2) then
           facp2 = 1./( ((ctime-t_start_averaging2)/dt) + 1)
        else
           facp2 = 1.
        endif
           facm2 = 1. - facp2

	Do M=1,bodynum ; do L=1,nodes(M)
	   FX1M(M,L)=facm1*FX1M(M,L)+facp1*FX1(M,L)
	   FX2M(M,L)=facm1*FX2M(M,L)+facp1*FX2(M,L)
	   FX3M(M,L)=facm1*FX3M(M,L)+facp1*FX3(M,L)
	ENDDO; enddo

! Write out the mean values
      if (mod(itime,n_out).eq.0) then
	Do M=1,bodynum

         write(ibnum,'(I2)') M
           strlen2=LEN(TRIM(ADJUSTL(ibnum)))
           ibnum=REPEAT('0',(2-strlen2))//TRIM(ADJUSTL(ibnum))
           write(char_block,'(I6)') itime
           strlen=LEN(TRIM(ADJUSTL(char_block)))
           char_block=REPEAT('0',(6-strlen))//TRIM(ADJUSTL(char_block))
           gridfile1='PnFmean_'//TRIM(ADJUSTL(char_block))//
     &	   '_'//TRIM(ADJUSTL(ibnum))//'.plt'
         open (unit=Geom_Time1, file=gridfile1)
	   write(Geom_Time1,*)'TITLE = "Mean Forces at IBs"'
	   write(Geom_Time1,*)'VARIABLES=x,y,z,FxM,FyM,FzM'
	   write(Geom_Time1,*)'zone i= ',nodes(M),'DATAPACKING = POINT'
	  do L=1,nodes(M)
	   write (Geom_Time1,89) nodex(M,L),nodey(M,L),nodez(M,L),
     & 	FX1M(M,L)*dxm*dym*dzm*1000.D0/dt,
     &  	FX2M(M,L)*dxm*dym*dzm*1000.D0/dt,
     &  	FX3M(M,L)*dxm*dym*dzm*1000.D0/dt
	  enddo
	 close(Geom_Time1)
	Enddo !M

	endif

	ENDIF !master

   89 FORMAT (17e14.5)
        RETURN
        END SUBROUTINE	
