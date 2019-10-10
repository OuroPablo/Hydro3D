!=======================================================================
!    			LAGRANGIAN PARTICLE TRACKING
!			Bruño Fraga Bugallo
!			Cardiff 2013-2015
!=======================================================================
!##########################################################################
        subroutine alloc_pt
!##########################################################################
        use vars
        use mpi
        use multidata
 	  use vars_pt
        implicit none

        integer l,out_cnt   
	  logical,allocatable,dimension(:):: out_pt			

!	wtime_refresh = MPI_WTIME ( ) 

	allocate(out_pt(np))

	out_pt=.false.

	out_cnt = 0
	
!	do ib=1,nbp
!      	is = dom(ib)%isp
!      	ie = dom(ib)%iep
!      	js = dom(ib)%jsp
!      	je = dom(ib)%jep
!       	ks = dom(ib)%ksp
!      	ke = dom(ib)%kep

	do l=1,np

!	Comprobar si permanece en el dominio
      	if ((xp_pt(l).le.xst).or.(xp_pt(l)
     &.ge.xen)) then
       		out_pt(l) = .TRUE.
	 		out_cnt = out_cnt + 1
	 		goto 50
      	elseif ((yp_pt(l).le.yst).or.(yp_pt(l)
     &.ge.yen)) then
       		out_pt(l) = .TRUE.
	 		out_cnt = out_cnt + 1
	 		goto 50
      	elseif ((zp_pt(l).lt.zst).or.(zp_pt(l)
     &.ge.zen)) then
       		out_pt(l) = .TRUE.
	 		out_cnt = out_cnt + 1
	 		goto 50
      	end if

50	continue

		if (.not.out_pt(l)) then
      		xpold(l-out_cnt)=xp_pt(l)
           		ypold(l-out_cnt)=yp_pt(l)
           		zpold(l-out_cnt)=zp_pt(l)
           		uopold(l-out_cnt)=uop_pt(l)
           		vopold(l-out_cnt)=vop_pt(l)
           		wopold(l-out_cnt)=wop_pt(l)
			dp_old(l-out_cnt)=dp_pt(l)
		endif

      	enddo
!	enddo

	np = np - out_cnt

	if (out_cnt.gt.0) then
		write(202,*)'Time:',ntime,'Removing',out_cnt,'particles'
     &,'. Total:',np
	endif

		deallocate (out_pt)
      	deallocate (xp_pt,yp_pt,zp_pt)
      	deallocate (uop_pt,vop_pt,wop_pt)
		deallocate (dp_pt)
		deallocate (Fu,Fv,Fw)

         	allocate (xp_pt(np),yp_pt(np),zp_pt(np))
         	allocate (uop_pt(np),vop_pt(np),wop_pt(np))
		allocate (dp_pt(np))
		allocate (Fu(np),Fv(np),Fw(np))

	do l=1,np
			xp_pt(l) = xpold(l)
           		yp_pt(l) = ypold(l)
           		zp_pt(l) = zpold(l)
           		uop_pt(l)= uopold(l)
           		vop_pt(l)= vopold(l)
           		wop_pt(l)= wopold(l)
			dp_pt(l) = dp_old(l)
	enddo

	deallocate (xpold,ypold,zpold,uopold,vopold,wopold,dp_old)

      allocate(xpold(np),ypold(np),zpold(np))
    	allocate(uopold(np),vopold(np),wopold(np)) 
	allocate(dp_old(np))

	do l=1,np
			xpold(l)=xp_pt(l)
           		ypold(l)=yp_pt(l)
           		zpold(l)=zp_pt(l)
           		uopold(l)=uop_pt(l)
           		vopold(l)=vop_pt(l)
           		wopold(l)=wop_pt(l)
			dp_old(l)=dp_pt(l)
	enddo	

!	  wtime_refresh = MPI_WTIME ( ) - wtime_refresh
!	  write(203,*)'reallocate parts: ',wtime_refresh	

      return
      end subroutine


!##########################################################################
        subroutine release_pt
!##########################################################################
        use vars
        use mpi
        use multidata
 	  use vars_pt
        implicit none

        integer l
	  integer X2a,X2b,X2c,X2d
        real X,X3,X4,X5,X6	
	  double precision random_number_normal

!	  wtime_release = MPI_WTIME ( )

	CALL RANDOM_SEED

!	open(35,file='LPT.cin')

!       	read(35,*)
!       	read(35,*)
!       	read(35,*) ptnr
!       	read(35,*) Dp,random_dp,sigma
!       	read(35,*) Dw
!		read(35,*)
!		read(35,*) random
!		div = 10./(0.5*Dw)	
       	np = np + ptnr

	  	write(202,*) ntime,'Releasing',ptnr,'new particles.', 
     &			' Total:',np

      	deallocate (xp_pt,yp_pt,zp_pt)
      	deallocate (uop_pt,vop_pt,wop_pt)
		deallocate (dp_pt)
		deallocate (Fu,Fv,Fw)

         	allocate (xp_pt(np),yp_pt(np),zp_pt(np))
         	allocate (uop_pt(np),vop_pt(np),wop_pt(np))
		allocate (dp_pt(np))
		allocate (Fu(np),Fv(np),Fw(np))

         do l=1,(np-ptnr)

         	xp_pt(l)=xpold(l)
         	yp_pt(l)=ypold(l)
         	zp_pt(l)=zpold(l)

         	uop_pt(l)=uopold(l)
         	vop_pt(l)=vopold(l)
         	wop_pt(l)=wopold(l)    

		dp_pt(l)=dp_old(l)

         end do

         deallocate (xpold,ypold,zpold)
         allocate(xpold(np),ypold(np),zpold(np))
         deallocate (uopold,vopold,wopold)
         allocate (uopold(np),vopold(np),wopold(np))
         deallocate (dp_old)
         allocate (dp_old(np))

	   if (random) then

!         	read(35,*)xp,yp,zp,uop,vop,wop

!	write(401,*)'A Release',
!     &	np-ptnr+1,xp_pt(np-ptnr+1),yp_pt(np-ptnr+1),
!     &	zp_pt(np-ptnr+1),out_pt(np-ptnr+1),' '

         	do l=(np-ptnr+1),np
			CALL RANDOM_NUMBER(X)
			X2a=int(X*10.)
			X3=X/10.	!1-10cm
			CALL RANDOM_NUMBER(X)
			X2b=int(X*10.)
			X4=X/100.	!1-10mm
			CALL RANDOM_NUMBER(X)
			X2c=int(X*10.)
			X5=X/1000.	!0.1-1mm
			CALL RANDOM_NUMBER(X)
			X2d=int(X*10.)
			X6=X/10000.	!0.01-0.1mm
			xp_pt(l)=xp+((-1.0)**X2a*X3					!til 4 cm (2 cm radius both sides)
     &		+(-1.0)**X2b*X4+(-1.0)**X2c*X5+(-1.0)**X2d*X6)/div		

			CALL RANDOM_NUMBER(X)
			X2a=int(X*10.)
			X3=X/10.	!1-10cm
			CALL RANDOM_NUMBER(X)
			X2b=int(X*10.)
			X4=X/100.	!1-10mm
			CALL RANDOM_NUMBER(X)
			X2c=int(X*10.)
			X5=X/1000.	!0.1-1mm
			CALL RANDOM_NUMBER(X)
			X2d=int(X*10.)
			X6=X/10000.	!0.01-0.1mm
			yp_pt(l)=yp+((-1.0)**X2a*X3					!til 4 cm (2 cm radius both sides)
     &		+(-1.0)**X2b*X4+(-1.0)**X2c*X5+(-1.0)**X2d*X6)/div	

!			CALL RANDOM_NUMBER(X)
!			X4=X/100.	!1-10mm
!			CALL RANDOM_NUMBER(X)
!			X2c=int(X*10.)
!			X5=X/1000.	!0.1-1mm
!			CALL RANDOM_NUMBER(X)
!			X2d=int(X*10.)
!			X6=X/10000.	!0.01-0.1mm
			zp_pt(l)=zp! + (X4							!til 5 mm (only up side)
!     &		+(-1.0)**X2c*X5+(-1.0)**X2d*X6)/2.	

			if (random_dp) then						!random dimeter of particle

			!CALL RANDOM_NUMBER(X)
			!X2c=int(X*10.)
			!X5=X/1000.	!0.1-1mm
			!CALL RANDOM_NUMBER(X)
			!X2d=int(X*10.)
			!X6=X/10000.	!0.01-0.1mm
			!dp_pt(l)= Dp + ((-1.0)**X2c*X5+(-1.0)**X2d*X6)
			dp_pt(l)= random_number_normal(Dp,sigma)

			else

				dp_pt(l) = Dp

			endif

			uop_pt(l)=uop
			vop_pt(l)=vop
			wop_pt(l)=wop

         	end do
	
		else										!if release isnt random, then Dp is constant
       		do l=np-ptnr+1,np
      			read(35,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &					uop_pt(l),vop_pt(l),wop_pt(l)
				xpold(l)=xp_pt(l)
				ypold(l)=yp_pt(l)
				zpold(l)=zp_pt(l)
				uopold(l)=uop_pt(l)
				vopold(l)=uop_pt(l)
				wopold(l)=uop_pt(l)
				dp_pt(l) = Dp
      		 end do
		endif	

!         close(35)

      return
      end subroutine
C#############################################################
      SUBROUTINE INIT_PARTICLE
C#############################################################
        use multidata
        use mpi
        use vars   
   	  use vars_pt

	implicit none	

	integer l,ib
	integer X2a,X2b,X2c,X2d
      real X,X3,X4,X5,X6	
	double precision :: Dw,random_number_normal
      integer :: i,j,k,tti,ttj,ttk

	CALL RANDOM_SEED

       open(10,file='in_LPT.cin')
	 read(10,*) PSIcell,order
       read(10,*) tsteps_pt,tsnr
       read(10,*) ptnr
       read(10,*) Dp,random_dp,sigma
       read(10,*) Dw
	 read(10,*) DF
	 read(10,*) random
	 div = 10./(0.5*Dw)

	 np=ptnr

       if (LRESTART) then
       	open(20,file='final_particle.dat')
       	read(20,*) np
       	read(20,*) cnt_pt
       end if

       allocate (xp_pt(np),yp_pt(np),zp_pt(np))
       allocate (uop_pt(np),vop_pt(np),wop_pt(np))
       allocate (xpold(np),ypold(np),zpold(np))
       allocate (uopold(np),vopold(np),wopold(np))
	 allocate (dp_pt(np),dp_old(np))
	 allocate (Fu(np),Fv(np),Fw(np))

	 allocate(ptsinproc(nprocs))

	 IF (myrank.eq.0) then

	 write(6,*)'random Dp=',random_dp

       if (LRESTART) then

       do l=1,np
       	read(20,*) xp_pt(l),yp_pt(l),zp_pt(l)
       	read(20,*) uop_pt(l),vop_pt(l),wop_pt(l)
		read(20,*) dp_pt(l)
		xpold(l)=xp_pt(l)
		ypold(l)=yp_pt(l)
		zpold(l)=zp_pt(l)
		uopold(l)=uop_pt(l)
		vopold(l)=vop_pt(l)
		wopold(l)=wop_pt(l)
		dp_old(l)=dp_pt(l)
		Fu(l)=0 ; Fv(l)=0 ; Fw(l)=0
       end do
       do ib=1,nbp
      	tti=dom(ib)%ttc_i
            ttj=dom(ib)%ttc_j
            ttk=dom(ib)%ttc_k
            do k=1,ttk
               do j=1,ttj
                  do i=1,tti
     				dom(ib)%uoo(i,j,k) = dom(ib)%u(i,j,k)
     				dom(ib)%voo(i,j,k) = dom(ib)%v(i,j,k)
     				dom(ib)%woo(i,j,k) = dom(ib)%w(i,j,k)
			enddo
               end do
            end do
	enddo

       elseif (.not.LRESTART) then
	  	write(202,*) 'First release:',np,'new particles. Total:',np

		if (random) then

         	read(10,*)xp,yp,zp,uop,vop,wop

!       	read(10,*) xp_pt(1),yp_pt(1),zp_pt(1),
!     &		uop_pt(1),vop_pt(1),wop_pt(1)
!		xpold(1) = xp_pt(1)
!		ypold(1) = yp_pt(1)
!		zpold(1) = zp_pt(1)
!		uopold(1)= uop_pt(1)
!		vopold(1)= vop_pt(1)
!		wopold(1)= wop_pt(1)
!		dp_pt(1) = Dp
!		dp_old(1)= Dp
!		Fu(1)=0 ; Fv(1)=0 ; Fw(1)=0

		if (np.gt.1) then

	       do l=1,np
			CALL RANDOM_NUMBER(X)
			X2a=int(X*10.)
			X3=X/10.	!1-10cm
			CALL RANDOM_NUMBER(X)
			X2b=int(X*10.)
			X4=X/100.	!1-10mm
			CALL RANDOM_NUMBER(X)
			X2c=int(X*10.)
			X5=X/1000.	!0.1-1mm
			CALL RANDOM_NUMBER(X)
			X2d=int(X*10.)
			X6=X/10000.	!0.01-0.1mm
			xp_pt(l)=xp+((-1.0)**X2a*X3				!til 4 cm (2 cm radius both sides)
     &		+(-1.0)**X2b*X4+(-1.0)**X2c*X5+(-1.0)**X2d*X6)/div		

			CALL RANDOM_NUMBER(X)
			X2a=int(X*10.)
			X3=X/10.	!1-10cm
			CALL RANDOM_NUMBER(X)
			X2b=int(X*10.)
			X4=X/100.	!1-10mm
			CALL RANDOM_NUMBER(X)
			X2c=int(X*10.)
			X5=X/1000.	!0.1-1mm
			CALL RANDOM_NUMBER(X)
			X2d=int(X*10.)
			X6=X/10000.	!0.01-0.1mm
			yp_pt(l)=yp+((-1.0)**X2a*X3				!til 4 cm (2 cm radius both sides)
     &		+(-1.0)**X2b*X4+(-1.0)**X2c*X5+(-1.0)**X2d*X6)/div	

!			CALL RANDOM_NUMBER(X)
!			X4=X/100.	!1-10mm
!			CALL RANDOM_NUMBER(X)
!			X2c=int(X*10.)
!			X5=X/1000.	!0.1-1mm
!			CALL RANDOM_NUMBER(X)
!			X2d=int(X*10.)
!			X6=X/10000.	!0.01-0.1mm
			zp_pt(l)=zp! + (X4						!til 5 mm (only up side)
!     &		+(-1.0)**X2c*X5+(-1.0)**X2d*X6)/2.	

			if (random_dp) then						!random dimeter of particle

			!CALL RANDOM_NUMBER(X)
			!X2c=int(X*10.)
			!X5=X/1000.	!0.1-1mm
			!CALL RANDOM_NUMBER(X)
			!X2d=int(X*10.)
			!X6=X/10000.	!0.01-0.1mm
			!dp_pt(l)= Dp + ((-1.0)**X2c*X5+(-1.0)**X2d*X6)
			dp_pt(l)= random_number_normal(Dp,sigma)

			else

				dp_pt(l) = Dp

			endif
			uop_pt(l)=uop
			vop_pt(l)=vop
			wop_pt(l)=wop

			xpold(l)=xp_pt(l)
			ypold(l)=yp_pt(l)
			zpold(l)=zp_pt(l)
			uopold(l)=uop_pt(l)
			vopold(l)=vop_pt(l)
			wopold(l)=wop_pt(l)
			dp_old(l)=dp_pt(l)
			Fu(l)=0 ; Fv(l)=0 ; Fw(l)=0
!			write(202,*)l,X,xp_pt(l),yp_pt(l),zp_pt(l)
	       end do

		endif
		else
       		do l=1,np
      			read(10,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &					uop_pt(l),vop_pt(l),wop_pt(l)
				xpold(l)=xp_pt(l)
				ypold(l)=yp_pt(l)
				zpold(l)=zp_pt(l)
				uopold(l)=uop_pt(l)
				vopold(l)=uop_pt(l)
				wopold(l)=uop_pt(l)
				dp_pt(l) = Dp						!if release isnt random, then Dp is constant
				dp_old(l)=dp_pt(l)
				Fu(l)=0 ; Fv(l)=0 ; Fw(l)=0
      		 end do
		endif		
 		call TECPARTICLE(0)
      end if										!restart

	endif											!myrank

      RETURN
      END SUBROUTINE

C **********************************************************************
      SUBROUTINE TECPLOT(num_output)
C **********************************************************************

        use multidata
        use mpi
        use vars   
   	  use vars_pt

	implicit none	

      integer strlen,i,j,k,ib,ni,nj,nk,ii,idfile,num_output
	integer is,ie,js,je,ks,ke
      character(LEN=18) filename
      character(LEN=3) b_str,c_str
	double precision u_cn,v_cn,w_cn,p_cn,S_cn,k_cn,eps_cn,vis_cn
   		
!        if (LRESTART) KK1=KK+KK2
!        if (LRESTART.eq..false.) KK1=KK

	do ib=1,nbp

	  idfile=600+dom_id(ib)

        write(b_str,'(I3)') num_output
        strlen=LEN(TRIM(ADJUSTL(b_str)))
        b_str=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(b_str)) ! e.g. "001"
        write(c_str,'(I3)') dom_id(ib)
        strlen=LEN(TRIM(ADJUSTL(c_str)))
        c_str=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(c_str)) ! e.g. "001"

      filename='tecout_'//b_str//'_'//c_str//'.dat'

      OPEN (UNIT=idfile, FILE=filename)

      WRITE (idfile,*) 'TITLE = ', '"Eulerian field"'
      WRITE (idfile,"(A)")'VARIABLES = "X","Y","Z","U","V","W","P","S",
     &"k","eps","vis"'


        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        ni=ie-(is-1)+1
        nj=je-(js-1)+1
        nk=ke-(ks-1)+1

      !WRITE(idfile,*)'ZONE T="','id:',dom_id(ib),'it:',ntime,'"'
      WRITE(idfile,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime
	WRITE(idfile,*)'I=',ni,', J=',nj,', K=',nk,'F=POINT'

		do k=ks-1,ke
		do j=js-1,je
		do i=is-1,ie

                 u_cn  =0.25*(dom(ib)%u(i,j,k)+
     &dom(ib)%u(i,j+1,k)+dom(ib)%u(i,j,k+1)+
     &dom(ib)%u(i,j+1,k+1))

                 v_cn  =0.25*(dom(ib)%v(i,j,k)+
     &dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j,k+1)+
     &dom(ib)%v(i+1,j,k+1))

                 w_cn  =0.25*(dom(ib)%w(i,j,k)+
     &dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j+1,k)+
     &dom(ib)%w(i+1,j+1,k)) 

                 p_cn  =0.125*(dom(ib)%p(i,j,k)+
     &dom(ib)%p(i+1,j,k)    +dom(ib)%p(i,j+1,k)+
     &dom(ib)%p(i+1,j+1,k)  +dom(ib)%p(i,j,k+1)+
     &dom(ib)%p(i+1,j,k+1)  +dom(ib)%p(i,j+1,k+1)+
     &dom(ib)%p(i+1,j+1,k+1))
                 S_cn  =0.125*(dom(ib)%S(i,j,k)+
     &dom(ib)%S(i+1,j,k)    +dom(ib)%S(i,j+1,k)+
     &dom(ib)%S(i+1,j+1,k)  +dom(ib)%S(i,j,k+1)+
     &dom(ib)%S(i+1,j,k+1)  +dom(ib)%S(i,j+1,k+1)+
     &dom(ib)%S(i+1,j+1,k+1))
                 k_cn  =0.125*(dom(ib)%ksgs(i,j,k)+
     &dom(ib)%ksgs(i+1,j,k)    +dom(ib)%ksgs(i,j+1,k)+
     &dom(ib)%ksgs(i+1,j+1,k)  +dom(ib)%ksgs(i,j,k+1)+
     &dom(ib)%ksgs(i+1,j,k+1)  +dom(ib)%ksgs(i,j+1,k+1)+
     &dom(ib)%ksgs(i+1,j+1,k+1))
                 eps_cn  =0.125*(dom(ib)%eps(i,j,k)+
     &dom(ib)%eps(i+1,j,k)    +dom(ib)%eps(i,j+1,k)+
     &dom(ib)%eps(i+1,j+1,k)  +dom(ib)%eps(i,j,k+1)+
     &dom(ib)%eps(i+1,j,k+1)  +dom(ib)%eps(i,j+1,k+1)+
     &dom(ib)%eps(i+1,j+1,k+1))
                 vis_cn  =0.125*(dom(ib)%vis(i,j,k)+
     &dom(ib)%vis(i+1,j,k)    +dom(ib)%vis(i,j+1,k)+
     &dom(ib)%vis(i+1,j+1,k)  +dom(ib)%vis(i,j,k+1)+
     &dom(ib)%vis(i+1,j,k+1)  +dom(ib)%vis(i,j+1,k+1)+
     &dom(ib)%vis(i+1,j+1,k+1))

      write (idfile,'(11e14.6)') dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k)
     & ,u_cn,v_cn,w_cn,p_cn,S_cn,k_cn,eps_cn,vis_cn

		enddo
		enddo
		enddo
!		write (90,*) dom(ib)%isp,dom(ib)%iep,
!     & dom(ib)%jsp,dom(ib)%jep,dom(ib)%ksp,dom(ib)%kep
	end do



      close (idfile)

!   88 FORMAT (10F15.8)

      END SUBROUTINE


C **********************************************************************
      SUBROUTINE TECPARTICLE(num_output)
C **********************************************************************
C
  
        use multidata
        use mpi
        use vars   
   	  use vars_pt

	implicit none	

      integer l,strlen,ib,num_output
      character(LEN=80) filename,filename2
      character(LEN=3) b_str

        write(b_str,'(I3)') num_output
        strlen=LEN(TRIM(ADJUSTL(b_str)))
        b_str=REPEAT('0',(3-strlen))//TRIM(ADJUSTL(b_str)) ! e.g. "001"

        filename='tecout_'//b_str//'_pt.dat'

        OPEN (UNIT=95, FILE=TRIM(ADJUSTL(filename)))

      WRITE (95,*) 'TITLE = ', '"Lagrangian field"'
      WRITE (95,"(A)")'VARIABLES = "X","Y","Z","U<sub>Lag<\sub>","V<sub>
     &Lag<\sub>","W<sub>Lag<\sub>","D<sub>Lag<\sub>","F<sub>u","F<sub>v"
     &,"F<sub>w"'
      WRITE (95,*) 'ZONE T= "', 'id:',ntime,'np:',np,'"'	

        do l=1,np
      	WRITE (95,*) xp_pt(l),yp_pt(l),zp_pt(l)
     &			,uop_pt(l),vop_pt(l),wop_pt(l)
     &			,dp_pt(l),Fu(l),Fv(l),Fw(l)
        end do
	  close (95)

!   88 FORMAT (10F15.8)

      END SUBROUTINE
!##########################################################################      
      subroutine particle_tracking

!	Calculates particles' velocities and the resulting source terms

!##########################################################################
        use multidata
        use mpi
        use vars   
   	  use vars_pt
        use omp_lib, only : omp_get_num_threads,
     &                    omp_get_thread_num,
     &                    omp_set_num_threads
      implicit none

      integer i,j,k,l
      integer ib,is,ie,js,je,ks,ke  
	integer nt,m
	integer iballs_u,iballe_u,jballs_u,jballe_u,kballs_u,kballe_u
	integer iballs_v,iballe_v,jballs_v,jballe_v,kballs_v,kballe_v
	integer iballs_w,iballe_w,jballs_w,jballe_w,kballs_w,kballe_w
      real REp,rho_p,rx,ry,rz
      double precision a,b,c,wx,wy,wz,Cd,ao,bo,co
	double precision dwdy,dvdz,dudz,dwdx,dvdx,dudy
	double precision dh,delta
	double precision Vcell,Vp!,Vball
!      double precision, allocatable, dimension(:):: Fpu,Fpv,Fpw
      double precision, allocatable, dimension(:):: up_pt,vp_pt,wp_pt
      double precision, allocatable, dimension(:):: ui_pt,vi_pt,wi_pt
      double precision, allocatable, dimension(:):: uoi_pt,voi_pt,woi_pt
	integer,allocatable,dimension(:)::  ip,jp,kp,ipu,jpv,kpw
	integer,dimension(nprocs) :: strider

      allocate (ui_pt(np_loc),vi_pt(np_loc),wi_pt(np_loc))
      allocate (uoi_pt(np_loc),voi_pt(np_loc),woi_pt(np_loc))
	allocate (ip(np_loc),jp(np_loc),kp(np_loc))
	allocate (ipu(np_loc),jpv(np_loc),kpw(np_loc))
	allocate (up_pt(np_loc),vp_pt(np_loc),wp_pt(np_loc))
      allocate (Fpu(np_loc),Fpv(np_loc),Fpw(np_loc))

	rho_p = 1.25d0
	dens = 1000.d0

	if (np_loc.le.100) nt = 1
	if (np_loc.gt.100) nt = OMP_threads

	call OMP_SET_NUM_THREADS(nt)

	SELECT CASE (order)
         	CASE (1)
			m = 1	!1.5d0
         	CASE (2)
			m = 2	!2.5d0
         	CASE (3)
			m = 2	!2.d0
         	CASE (4)
			m = 2	!2.5d0
         	CASE (5)
			m = 1	!1.5d0
         	CASE (6)
			m = 2	!2.d0
	end select

!loop in domains
	do ib=1,nbp

		Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz	

!computational domain limits (one-axis index)	
      	is = dom(ib)%isp
      	ie = dom(ib)%iep
      	js = dom(ib)%jsp
      	je = dom(ib)%jep
       	ks = dom(ib)%ksp
      	ke = dom(ib)%kep

!loop in particles
!$OMP 	PARALLEL DEFAULT (SHARED), PRIVATE(i,j,k,l,
!$OMP&	iballs_u,iballe_u,jballs_u,jballe_u,kballs_u,kballe_u,
!$OMP&	iballs_v,iballe_v,jballs_v,jballe_v,kballs_v,kballe_v,
!$OMP&	iballs_w,iballe_w,jballs_w,jballe_w,kballs_w,kballe_w,
!$OMP&	REp,rx,ry,rz,Vp,delta,
!$OMP&	a,b,c,ao,bo,co,Cd,wx,wy,wz,
!$OMP&	dwdy,dvdz,dudz,dvdx,dudy,dwdx)

!$OMP DO SCHEDULE (DYNAMIC,1)
      do l=1,np_loc

	IF (id(l).eq.dom_id(ib)) then							!particle belongs to THIS block
	
!      if ((xp_loc(l).lt.dom(ib)%x(is-1))  .or.
!     & 	(xp_loc(l).gt.dom(ib)%x(ie))  .or.
!     &	(yp_loc(l).lt.dom(ib)%y(js-1)).or.
!     &      (yp_loc(l).gt.dom(ib)%y(je))  .or.
!     &	(zp_loc(l).lt.dom(ib)%z(ks-1)).or.
!     &	(zp_loc(l).gt.dom(ib)%z(ke)))
		
!		goto 400
!	endif

	Vp = 3.1416*dp_loc(l)**3.d0/6.d0

!loop in cells
!      do k=ks,ke ; do j=js,je ; do i=is,ie
!checking in what cell is the part centre
!      if ((xp_loc(l).ge.dom(ib)%x(i-1)).and.
!     & (xp_loc(l).le.dom(ib)%x(i))) then

!      if ((yp_loc(l).ge.dom(ib)%y(j-1)).and.
!     & (yp_loc(l).le.dom(ib)%y(j))) then

!      if ((zp_loc(l).ge.dom(ib)%z(k-1)).and.
!     & (zp_loc(l).le.dom(ib)%z(k))) then 

	ip(l)=INT((xp_loc(l)-dom(ib)%x(is-1)-1.d-12)/dom(ib)%dx)+1+pl
	jp(l)=INT((yp_loc(l)-dom(ib)%y(js-1)-1.d-12)/dom(ib)%dy)+1+pl
	kp(l)=INT((zp_loc(l)-dom(ib)%z(ks-1)-1.d-12)/dom(ib)%dz)+1+pl

!	write(myrank+700,*)'======================================='
!     	write(myrank+700,*) 'proc and dom:',myrank,'->',dom_id(ib)
!     	write(myrank+700,*) 'tstep:',itime
!     	write(myrank+700,*)'ip',ip(l),jp(l),kp(l)


!	goto 100

!      endif ; endif ; endif
!	enddo ; enddo ; enddo

!100	continue

!locate the u,v and w nodes

      if (xp_loc(l).gt.dom(ib)%xc(ip(l))) then

	ipu(l) = ip(l)
!	ipux(l) = ip(l) - 1

      elseif (xp_loc(l).le.dom(ib)%xc(ip(l))) then

	ipu(l) = ip(l) - 1
!	ipux(l) = ip(l)

      endif

      if (yp_loc(l).gt.dom(ib)%yc(jp(l))) then

	jpv(l) = jp(l)
!	jpvy(l) = jp(l) - 1

      elseif (yp_loc(l).le.dom(ib)%yc(jp(l))) then

	jpv(l) = jp(l) - 1
!	jpvy(l) = jp(l)

      end if

      if (zp_loc(l).gt.dom(ib)%zc(kp(l))) then

	kpw(l) = kp(l)
!	kpwz(l) = kp(l) - 1

      elseif (zp_loc(l).le.dom(ib)%zc(kp(l))) then

	kpw(l) = kp(l) - 1
!	kpwz(l) = kp(l)

      end if		

!     	write(myrank+700,*)'ipu',ipu(l),jpv(l),kpw(l)

	rx = max(dp_loc(l),dom(ib)%dx)
	ry = max(dp_loc(l),dom(ib)%dy)
	rz = max(dp_loc(l),dom(ib)%dz)

	!Ball

	IF (order.eq.3.or.order.eq.6) then
	if (ipu(l).eq.ip(l)) then
		iballs_u = ipu(l) - 1 * NINT(rx/dom(ib)%dx)
		iballe_u = ipu(l) + m * NINT(rx/dom(ib)%dx)
		iballs_v = ip(l) - 1 * NINT(rx/dom(ib)%dx)
		iballe_v = ip(l) + m * NINT(rx/dom(ib)%dx)
		iballs_w = ip(l) - 1 * NINT(rx/dom(ib)%dx)
		iballe_w = ip(l) + m * NINT(rx/dom(ib)%dx)
	else
		iballs_u = ipu(l) - m * NINT(rx/dom(ib)%dx)
		iballe_u = ipu(l) + 1 * NINT(rx/dom(ib)%dx)
		iballs_v = ip(l) - m * NINT(rx/dom(ib)%dx)
		iballe_v = ip(l) + 1 * NINT(rx/dom(ib)%dx)
		iballs_w = ip(l) - m * NINT(rx/dom(ib)%dx)
		iballe_w = ip(l) + 1 * NINT(rx/dom(ib)%dx)
	endif
	if (jpv(l).eq.jp(l)) then
		jballs_u = jp(l) - 1 * NINT(ry/dom(ib)%dy)
		jballe_u = jp(l) + m * NINT(ry/dom(ib)%dy)
		jballs_v = jpv(l) - 1 * NINT(ry/dom(ib)%dy)
		jballe_v = jpv(l) + m * NINT(ry/dom(ib)%dy)
		jballs_w = jp(l) - 1 * NINT(ry/dom(ib)%dy)
		jballe_w = jp(l) + m * NINT(ry/dom(ib)%dy)
	else
		jballs_u = jp(l) - m * NINT(ry/dom(ib)%dy)
		jballe_u = jp(l) + 1 * NINT(ry/dom(ib)%dy)
		jballs_v = jpv(l) - m * NINT(ry/dom(ib)%dy)
		jballe_v = jpv(l) + 1 * NINT(ry/dom(ib)%dy)
		jballs_w = jp(l) - m * NINT(ry/dom(ib)%dy)
		jballe_w = jp(l) + 1 * NINT(ry/dom(ib)%dy)
	endif
	if (kpw(l).eq.kp(l)) then
		kballs_u = kp(l) - 1 * NINT(rz/dom(ib)%dz)
		kballe_u = kp(l) + m * NINT(rz/dom(ib)%dz)
		kballs_v = kp(l) - 1 * NINT(rz/dom(ib)%dz)
		kballe_v = kp(l) + m * NINT(rz/dom(ib)%dz)
		kballs_w = kpw(l) - 1 * NINT(rz/dom(ib)%dz)
		kballe_w = kpw(l) + m * NINT(rz/dom(ib)%dz)
	else
		kballs_u = kp(l) - m * NINT(rz/dom(ib)%dz)
		kballe_u = kp(l) + 1 * NINT(rz/dom(ib)%dz)
		kballs_v = kp(l) - m * NINT(rz/dom(ib)%dz)
		kballe_v = kp(l) + 1 * NINT(rz/dom(ib)%dz)
		kballs_w = kpw(l) - m * NINT(rz/dom(ib)%dz)
		kballe_w = kpw(l) + 1 * NINT(rz/dom(ib)%dz)
	endif
	ELSE	
		iballs_u = ipu(l) - m * NINT(rx/dom(ib)%dx)
		iballe_u = ipu(l) + m * NINT(rx/dom(ib)%dx)
		jballs_u = jp(l) - m * NINT(ry/dom(ib)%dy)
		jballe_u = jp(l) + m * NINT(ry/dom(ib)%dy)
		kballs_u = kp(l) - m * NINT(rz/dom(ib)%dz)
		kballe_u = kp(l) + m * NINT(rz/dom(ib)%dz)

		iballs_v = ip(l) - m * NINT(rx/dom(ib)%dx)
		iballe_v = ip(l) + m * NINT(rx/dom(ib)%dx)
		jballs_v = jpv(l) - m * NINT(ry/dom(ib)%dy)
		jballe_v = jpv(l) + m * NINT(ry/dom(ib)%dy)
		kballs_v = kp(l) - m * NINT(rz/dom(ib)%dz)
		kballe_v = kp(l) + m * NINT(rz/dom(ib)%dz)

		iballs_w = ip(l) - m * NINT(rx/dom(ib)%dx)
		iballe_w = ip(l) + m * NINT(rx/dom(ib)%dx)
		jballs_w = jp(l) - m * NINT(ry/dom(ib)%dy)
		jballe_w = jp(l) + m * NINT(ry/dom(ib)%dy)
		kballs_w = kpw(l) - m * NINT(rz/dom(ib)%dz)
		kballe_w = kpw(l) + m * NINT(rz/dom(ib)%dz)
	ENDIF

	iballs_u = max(iballs_u,1)
	iballe_u = min(iballe_u,dom(ib)%ttc_i)
	jballs_u = max(jballs_u,1)
	jballe_u = min(jballe_u,dom(ib)%ttc_j)
	kballs_u = max(kballs_u,1)
	kballe_u = min(kballe_u,dom(ib)%ttc_k)
	iballs_v = max(iballs_v,1)
	iballe_v = min(iballe_v,dom(ib)%ttc_i)
	jballs_v = max(jballs_v,1)
	jballe_v = min(jballe_v,dom(ib)%ttc_j)
	kballs_v = max(kballs_v,1)
	kballe_v = min(kballe_v,dom(ib)%ttc_k)
	iballs_w = max(iballs_w,1)
	iballe_w = min(iballe_w,dom(ib)%ttc_i)
	jballs_w = max(jballs_w,1)
	jballe_w = min(jballe_w,dom(ib)%ttc_j)
	kballs_w = max(kballs_w,1)
	kballe_w = min(kballe_w,dom(ib)%ttc_k)

!     	write(myrank+700,*)'ip_u',iballs_u,iballe_u
!     	write(myrank+700,*)'ip_w',iballs_w,iballe_w
!     	write(myrank+700,*)'kp_u',kballs_u,kballe_u
!     	write(myrank+700,*)'kp_w',kballs_w,kballe_w


!	if (m.ne.2.d0) then
!	Vball = (iballe-iballs+1)*(jballe-jballs+1)*
!     &	  (kballe-kballs+1)*Vcell
!	else	
!	Vball = (iballe-iballs)*(jballe-jballs)*
!     &	  (kballe-kballs)*Vcell
!	endif

		uoi_pt(l) = 0.d0
		voi_pt(l) = 0.d0
		woi_pt(l) = 0.d0

	do i=iballs_u,iballe_u 
	do j=jballs_u,jballe_u 
	do k=kballs_u,kballe_u

		uoi_pt(l) = uoi_pt(l) + dom(ib)%uoo(i,j,k)*
     &	dh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j)
     &	,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

	enddo ; enddo ; enddo

	do i=iballs_v,iballe_v 
	do j=jballs_v,jballe_v 
	do k=kballs_v,kballe_v

		voi_pt(l) = voi_pt(l) + dom(ib)%voo(i,j,k)*
     &	dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j)
     &	,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

	enddo ; enddo ; enddo

	do i=iballs_w,iballe_w 
	do j=jballs_w,jballe_w 
	do k=kballs_w,kballe_w

		woi_pt(l) = woi_pt(l) + dom(ib)%woo(i,j,k)*
     &	dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j)
     &	,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

	enddo ; enddo ; enddo

!     	write(myrank+700,*)'uoi',uoi_pt(l),voi_pt(l),woi_pt(l)

		up_pt(l) = uop_loc(l)
		vp_pt(l) = vop_loc(l)
		wp_pt(l) = wop_loc(l)
	
		ui_pt(l) = 0.d0
		vi_pt(l) = 0.d0
		wi_pt(l) = 0.d0

	do i=iballs_u,iballe_u 
	do j=jballs_u,jballe_u 
	do k=kballs_u,kballe_u

     		delta = dh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j)
     &	,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

		ui_pt(l) = ui_pt(l) + dom(ib)%ustar(i,j,k) * delta

	enddo ; enddo ; enddo

	do i=iballs_v,iballe_v 
	do j=jballs_v,jballe_v 
	do k=kballs_v,kballe_v

		delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j)
     &	,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

		vi_pt(l) = vi_pt(l) + dom(ib)%vstar(i,j,k) * delta

	enddo ; enddo ; enddo

	do i=iballs_w,iballe_w 
	do j=jballs_w,jballe_w 
	do k=kballs_w,kballe_w

		delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j)
     &	,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

		wi_pt(l) = wi_pt(l) + dom(ib)%wstar(i,j,k) * delta

	enddo ; enddo ; enddo

!     	write(myrank+700,*)'ui',ui_pt(l),vi_pt(l),wi_pt(l)

	!Slip vel components 
      a = up_pt(l)-ui_pt(l)
      b = vp_pt(l)-vi_pt(l)
	c = wp_pt(l)-wi_pt(l)

!     	write(myrank+700,*)'wslip',a,b,c

       REp = dp_loc(l)* (sqrt((uop_loc(l)-ui_pt(l))**2.d0+(vop_loc(l)
     &  	-vi_pt(l))**2.d0+(wop_loc(l)-wi_pt(l))**2.d0))/(1.d0/Re)

        if (REp.le.800) Cd = 24.d0*(1.d0+0.15d0*(REp**0.687d0))/REp
        if (REp.gt.800) Cd = 0.44d0

	!Vorticity calculation	
	dwdy=(dom(ib)%wstar(ip(l),jp(l)+1,kpw(l))
     &	-dom(ib)%wstar(ip(l),jp(l)-1,kpw(l)))
     &	/(dom(ib)%y(jp(l)+1)-dom(ib)%y(jp(l)-1))	
	dvdz=(dom(ib)%vstar(ip(l),jpv(l),kp(l)+1)
     &	-dom(ib)%vstar(ip(l),jpv(l),kp(l)-1))
     &	/(dom(ib)%z(kp(l)+1)-dom(ib)%z(kp(l)-1))		
	wx = dwdy-dvdz
	dudz=(dom(ib)%ustar(ipu(l),jp(l),kp(l)+1)
     &	-dom(ib)%ustar(ipu(l),jp(l),kp(l)-1))
     &	/(dom(ib)%z(kp(l)+1)-dom(ib)%z(kp(l)-1))	
	dwdx=(dom(ib)%wstar(ip(l)+1,jp(l),kpw(l))
     &	-dom(ib)%wstar(ip(l)-1,jp(l),kpw(l)))
     &	/(dom(ib)%x(ip(l)+1)-dom(ib)%x(ip(l)-1))		
	wy = dudz-dwdx
	dvdx=(dom(ib)%vstar(ip(l)+1,jpv(l),kp(l))
     &	-dom(ib)%vstar(ip(l)-1,jpv(l),kp(l)))
     &	/(dom(ib)%x(ip(l)+1)-dom(ib)%x(ip(l)-1))	
	dudy=(dom(ib)%ustar(ipu(l),jp(l)+1,kp(l))
     &	-dom(ib)%ustar(ipu(l),jp(l)-1,kp(l)))
     &	/(dom(ib)%y(jp(l)+1)-dom(ib)%y(jp(l)-1))			
	wz = dvdx-dudy

      up_pt(l) = uop_loc(l) + dt * (3.0d0*((ui_pt(l)-uoi_pt(l))/dt)	!Buoyancy, stress
     &  -(3.0d0/(2.0d0*dp_loc(l)))*Cd*sqrt(a**2.d0+b**2.d0+c**2.d0)*a	!Added Mass and drag
     &  -2.0d0*0.53d0*(b*wz-c*wy))								!Lift


      vp_pt(l) = vop_loc(l) + dt* (3.0d0*((vi_pt(l)-voi_pt(l))/dt)
     &  -(3.0d0/(2.0d0*dp_loc(l)))*Cd*sqrt(a**2.d0+b**2.d0+c**2.d0)*b
     &  -2.0d0*0.53d0*(c*wx-a*wz))


      wp_pt(l) = wop_loc(l) + dt* (2.0d0*9.81d0+3.0d0*
     &	((wi_pt(l)-woi_pt(l))/dt)	
     &-(3.0d0/(2.0d0*dp_loc(l)))*Cd*sqrt(a**2.d0+b**2.d0+c**2.d0)*c		
     &-2.0d0*0.53d0*(a*wy-b*wx))


!     	write(myrank+700,*)'up',up_pt(l),vp_pt(l),wp_pt(l)

	!New slip vel components 
      a = up_pt(l)-ui_pt(l)
      b = vp_pt(l)-vi_pt(l)
	c = wp_pt(l)-wi_pt(l)

	!Fluid stresses
!	Fsu(l) = dens*3.14d0*dp_loc(l)**3.d0*((ui_pt(l)-uoi_pt(l))/dt)/6.0d0
!	Fsv(l) = dens*3.14d0*dp_loc(l)**3.d0*((vi_pt(l)-voi_pt(l))/dt)/6.0d0
!	Fsw(l) = dens*3.14d0*dp_loc(l)**3.d0*((wi_pt(l)-woi_pt(l))/dt)/6.0d0

	!Added mass
!	Fau(l) = -dens*3.14d0*(dp_loc(l)**3.d0)*(a-ao)/(12.0d0*dt)
!	Fav(l) = -dens*3.14d0*(dp_loc(l)**3.d0)*(b-bo)/(12.0d0*dt)
!	Faw(l) = -dens*3.14d0*(dp_loc(l)**3.d0)*(c-co)/(12.0d0*dt)

	!Drag
!	Fdu(l) = -dens*3.14d0*(dp_loc(l)**2.d0)*Cd*(sqrt(a**2.d0+b**2.d0
!     &	+c**2.d0))*a/8.0d0
!	Fdv(l) = -dens*3.14d0*(dp_loc(l)**2.d0)*Cd*(sqrt(a**2.d0+b**2.d0
!     &	+c**2.d0))*b/8.0d0
!	Fdw(l) = -dens*3.14d0*(dp_loc(l)**2.d0)*Cd*(sqrt(a**2.d0+b**2.d0
!     &	+c**2.d0))*c/8.0d0

	
	!Lift
!      Flu(l) = -0.53d0*rho_p*3.14d0*(dp_loc(l)**3.d0)*(b*wz-c*wy)
!     &	/6.0d0
!      Flv(l) = -0.53d0*rho_p*3.14d0*(dp_loc(l)**3.d0)*(c*wx-a*wz)
!     &	/6.0d0
!      Flw(l) = -0.53d0*rho_p*3.14d0*(dp_loc(l)**3.d0)*(a*wy-b*wx)
!     &	/6.0d0

	!Interphase Force (bubble->liquid)

!	Fpu(l) = -(Fau(l) + Fdu(l) + Flu(l))
!	Fpv(l) = -(Fav(l) + Fdv(l) + Flv(l))
!	Fpw(l) = -(Faw(l) + Fdw(l) + Flw(l))

	if (.not.DF) then

      Fpu(l) = -(3.0d0*((ui_pt(l)-uoi_pt(l))/dt)	
     &  -(3.0d0/(2.0d0*dp_loc(l)))*Cd*sqrt(a**2.d0+b**2.d0+c**2.d0)*a  
     &  -2.0d0*0.53d0*(b*wz-c*wy))				

      Fpv(l) = -(3.0d0*((vi_pt(l)-voi_pt(l))/dt)
     &  -(3.0d0/(2.0d0*dp_loc(l)))*Cd*sqrt(a**2.d0+b**2.d0+c**2.d0)*b  
     &  -2.0d0*0.53d0*(c*wx-a*wz))

      Fpw(l) = -(3.0d0*((wi_pt(l)-woi_pt(l))/dt)	
     &-(3.0d0/(2.0d0*dp_loc(l)))*Cd*sqrt(a**2.d0+b**2.d0+c**2.d0)*c  			
     &-2.0d0*0.53d0*(a*wy-b*wx))

	else

	Fpu(l) = (a/dt) * (rho_p/dens)!- Fsu(l)/(Vp*dens)
	Fpv(l) = (b/dt) * (rho_p/dens) !- Fsv(l)/(Vp*dens)
	Fpw(l) = (c/dt) * (rho_p/dens) !- (Fsw(l) - Fgw(l))/(Vp*dens)
!	write(201,*) 'Fpw',Fpw(l)

	endif

!     	write(myrank+700,*)'Fp',Fpu(l),Fpv(l),Fpw(l)

!$OMP CRITICAL
	if (PSIcell) then

           	dom(ib)%ustar(ipu(l),jp(l),kp(l)) = 
     &	dom(ib)%ustar(ipu(l),jp(l),kp(l)) + dt * alfapr * Fpu(l) *
     &      Vp/Vcell

           	dom(ib)%vstar(ip(l),jpv(l),kp(l)) = 
     &	dom(ib)%vstar(ip(l),jpv(l),kp(l)) + dt * alfapr * Fpv(l) *
     &      Vp/Vcell

           	dom(ib)%wstar(ip(l),jp(l),kpw(l)) = 
     &	dom(ib)%wstar(ip(l),jp(l),kpw(l)) + dt * alfapr * Fpw(l) *
     &      Vp/Vcell
	else
		do i=iballs_u,iballe_u 
		do j=jballs_u,jballe_u 
		do k=kballs_u,kballe_u

		delta = dh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j)				
     &	,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

           	dom(ib)%ustar(i,j,k) = 
     &	dom(ib)%ustar(i,j,k) + dt * alfapr * Fpu(l) * delta *
     &      Vp/Vcell

		enddo ; enddo ; enddo


		do i=iballs_v,iballe_v 
		do j=jballs_v,jballe_v 
		do k=kballs_v,kballe_v

     		delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j)
     &	,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

           	dom(ib)%vstar(i,j,k) = 
     &	dom(ib)%vstar(i,j,k) + dt * alfapr * Fpv(l) * delta *
     &      Vp/Vcell

		enddo ; enddo ; enddo


		do i=iballs_w,iballe_w 
		do j=jballs_w,jballe_w 
		do k=kballs_w,kballe_w

	      delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j)
     &	,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

           	dom(ib)%wstar(i,j,k) = 
     &	dom(ib)%wstar(i,j,k) + dt * alfapr * Fpw(l) * delta *
     &      Vp/Vcell

		enddo ; enddo ; enddo

	endif

!$OMP END CRITICAL

!	Actualizar velocidad paso previo
	      uop_loc(l) = up_pt(l)
      	vop_loc(l) = vp_pt(l)
      	wop_loc(l) = wp_pt(l)

!	Actualizar posicion de particula
      	xp_loc(l)=xp_loc(l)+up_pt(l)*dt
      	yp_loc(l)=yp_loc(l)+vp_pt(l)*dt
      	zp_loc(l)=zp_loc(l)+wp_pt(l)*dt

!     		write(myrank+700,*)'xp_loc',xp_loc(l),yp_loc(l),zp_loc(l)

	ENDIF	  !if the particle belongs to the block

      end do  !end of loop in particles
!$OMP ENDDO
!$OMP END PARALLEL


      end do	!end loop in domains

!		deallocate (Fsu,Fau,Fdu,Flu)
!		deallocate (Fsv,Fav,Fdv,Flv)
!		deallocate (Fgw,Fsw,Faw,Fdw,Flw)
      	deallocate (ui_pt,vi_pt,wi_pt,uoi_pt,voi_pt,woi_pt)
		deallocate (ip,jp,kp,ipu,jpv,kpw)
      	deallocate (up_pt,vp_pt,wp_pt)
		deallocate (id)
		deallocate (dp_loc)
!		deallocate (ipux,jpvy,kpwz)
!		deallocate (dh_acumu,dh_acumv,dh_acumw)

      return
      end subroutine

!##########################################################################      
      subroutine final_LPT

!	sends backp(l) to master processor 

!##########################################################################

        use multidata
        use mpi
   	  use vars_pt

      implicit none

	integer,dimension(nprocs) :: strider
	integer s

	strider(1) = 0
	do s=2,nprocs
		strider(s) = ptsinproc(s-1) + strider(s-1)
	enddo

!	do s=1,nprocs
!     		write(myrank+800,*)'proc',s,'At',ntime
!     		write(myrank+800,*)'pts',ptsinproc(s)
!     		write(myrank+800,*)'strider',strider(s)
!	enddo

       call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_GATHERV(xp_loc,np_loc,MPI_DOUBLE_PRECISION,xp_pt
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(yp_loc,np_loc,MPI_DOUBLE_PRECISION,yp_pt
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
        call MPI_GATHERV(zp_loc,np_loc,MPI_DOUBLE_PRECISION,zp_pt
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

        call MPI_GATHERV(uop_loc,np_loc,MPI_DOUBLE_PRECISION,uop_pt
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
        call MPI_GATHERV(vop_loc,np_loc,MPI_DOUBLE_PRECISION,vop_pt
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(wop_loc,np_loc,MPI_DOUBLE_PRECISION,wop_pt
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

        call MPI_GATHERV(Fpu,np_loc,MPI_DOUBLE_PRECISION,Fu
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
        call MPI_GATHERV(Fpv,np_loc,MPI_DOUBLE_PRECISION,Fv
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHERV(Fpw,np_loc,MPI_DOUBLE_PRECISION,Fw
     &,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

!	if (myrank.eq.0) then
!		do l=1,np
!     		write(myrank+800,*)'xp',xp_pt(l),yp_pt(l),zp_pt(l)
!		enddo
!	endif

	if (np_loc.gt.0) then
      	deallocate (xp_loc,yp_loc,zp_loc)
      	deallocate (uop_loc,vop_loc,wop_loc)
      	deallocate (Fpu,Fpv,Fpw)
	endif

      return
      end subroutine
!######################################################################
      SUBROUTINE MPI_pt
!	Distributes the particles among the processors
!######################################################################
      use vars
      use mpi
      use multidata
	use vars_pt

      implicit none

      integer :: l,o,ii,nxdom,nydom,nzdom,ntot
	real :: lx,ly,lz
	integer,allocatable,dimension(:)::	lpt_proc,lpt_block
	integer,allocatable,dimension(:)::	id_MPI,id_MPI_loc
      double precision, allocatable, dimension(:):: X_MPI,Y_MPI,Z_MPI
      double precision, allocatable, dimension(:):: X_MPI_loc,Y_MPI_loc
      double precision, allocatable, dimension(:):: Z_MPI_loc
      double precision, allocatable, dimension(:):: U_MPI,V_MPI,W_MPI
      double precision, allocatable, dimension(:):: U_MPI_loc,V_MPI_loc
      double precision, allocatable, dimension(:):: W_MPI_loc
      double precision, allocatable, dimension(:):: dp_MPI,dp_MPI_loc

	call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      call MPI_BCAST(np,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	ntot=np*nprocs

	allocate(X_MPI(ntot),Y_MPI(ntot),Z_MPI(ntot))
	allocate(U_MPI(ntot),V_MPI(ntot),W_MPI(ntot))
	allocate(dp_MPI(ntot),id_MPI(ntot))
	allocate(X_MPI_loc(np),Y_MPI_loc(np),Z_MPI_loc(np))
	allocate(U_MPI_loc(np),V_MPI_loc(np),W_MPI_loc(np))
	allocate(dp_MPI_loc(np),id_MPI_loc(np))
	allocate(lpt_block(np),lpt_proc(np))

      IF (Myrank.eq.0) THEN

	X_MPI = -1.d12

	lx=xen-xst
	ly=yen-yst
	lz=zen-zst

	ptsinproc=0  ; lpt_block=0  ; lpt_proc=0

		do l=1,np

			nxdom=INT((xp_pt(l)-xst-1d-12)/(lx/idom))				!Domain where point belongs in i

			nydom=INT((yp_pt(l)-yst-1d-12)/(ly/jdom))				!Domain where point belongs in j

			nzdom=INT((zp_pt(l)-zst-1d-12)/(lz/kdom))				!Domain where point belongs in k

			lpt_block(l)=idom*jdom*(nzdom)+idom*(nydom)
     &					+nxdom					!Domain where point belongs (0 to num_dom)

			lpt_proc(l)=dom_ad(lpt_block(l))+1 				!Processor where point belongs (1 to num_procs+1)

!			ptsinblk(lpt_block(l)+1)=ptsinblk(lpt_block(l)+1)+1	!array which tells how many points are in each domain

			ptsinproc(lpt_proc(l))=ptsinproc(lpt_proc(l))+1		!array which tells how many points are in each processor


!		write(6,*)lpt_block(l)
!     		write(6,*)xp_pt(l)
!     		write(6,*)(lx/idom)*(nxdom)
!		write(6,*)'======================================='

		enddo

		do o=1,nprocs
!		write(6,*)'1',o,ptsinproc(o)
	  	   ii=1
		   DO l=1,np
			if(lpt_proc(l).eq.(o)) then
				X_MPI(ii+np*(o-1))=xp_pt(l)				!superarrays
				Y_MPI(ii+np*(o-1))=yp_pt(l)
				Z_MPI(ii+np*(o-1))=zp_pt(l)

				U_MPI(ii+np*(o-1))=uop_pt(l)
				V_MPI(ii+np*(o-1))=vop_pt(l)
				W_MPI(ii+np*(o-1))=wop_pt(l)

				dp_MPI(ii+np*(o-1))=dp_pt(l)

				id_MPI(ii+np*(o-1))=lpt_block(l)

				ii=ii+1
			endif
		   ENDDO
		enddo

	ENDIF

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_BCAST(ptsinproc,nprocs,MPI_INTEGER,0,
     &			MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(ptsinproc,1,MPI_INTEGER,np_loc,1,
     &		MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!	if (np_loc.eq.0) then
!		np_loc=1
!	endif


!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(X_MPI,np,MPI_DOUBLE_PRECISION,X_MPI_loc,		!location
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(Y_MPI,np,MPI_DOUBLE_PRECISION,Y_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(Z_MPI,np,MPI_DOUBLE_PRECISION,Z_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(U_MPI,np,MPI_DOUBLE_PRECISION,U_MPI_loc,		!velocities (of prev tstep)
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(V_MPI,np,MPI_DOUBLE_PRECISION,V_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(W_MPI,np,MPI_DOUBLE_PRECISION,W_MPI_loc,
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(dp_MPI,np,MPI_DOUBLE_PRECISION,dp_MPI_loc,	!diameter
     &		np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(id_MPI,np,MPI_INTEGER,id_MPI_loc,			!block to which they belong
     &		np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


	if (np_loc.gt.0) then

      allocate (xp_loc(np_loc),yp_loc(np_loc),zp_loc(np_loc))
      allocate (uop_loc(np_loc),vop_loc(np_loc),wop_loc(np_loc))
	allocate (dp_loc(np_loc))
!	allocate (Fsu(np_loc),Fau(np_loc),Fdu(np_loc),Flu(np_loc))
!	allocate (Fsv(np_loc),Fav(np_loc),Fdv(np_loc),Flv(np_loc))
!	allocate (Faw(np_loc),Fdw(np_loc),Flw(np_loc))
!	allocate (Fgw(np_loc),Fsw(np_loc))
!	allocate (ipux(np_loc),jpvy(np_loc),kpwz(np_loc))
	allocate (id(np_loc))
!	allocate (dh_acumu(np_loc),dh_acumv(np_loc),dh_acumw(np_loc))


		ii=0
		do l=1,np
			if (X_MPI_loc(l).ne.-1.d12) then	
				ii=ii+1		
				xp_loc(l)=X_MPI_loc(l)
				yp_loc(l)=Y_MPI_loc(l)
				zp_loc(l)=Z_MPI_loc(l)

				uop_loc(l)=U_MPI_loc(l)
				vop_loc(l)=V_MPI_loc(l)
				wop_loc(l)=W_MPI_loc(l)

				dp_loc(l)=dp_MPI_loc(l)

				id(l)=id_MPI_loc(l)
			endif
		enddo
		if (ii.ne.np_loc) then
			write(6,*)'MPI ERROR in proc:',myrank
			write(6,*)'np_loc=',np_loc,'=/=',ii
			stop
		endif

	endif

		deallocate(lpt_proc,lpt_block)
		deallocate(X_MPI,Y_MPI,Z_MPI,U_MPI,V_MPI,W_MPI)
		deallocate(X_MPI_loc,Y_MPI_loc,Z_MPI_loc)
		deallocate(U_MPI_loc,V_MPI_loc,W_MPI_loc)
		deallocate(dp_MPI,dp_MPI_loc,id_MPI,id_MPI_loc)

      RETURN
      END
