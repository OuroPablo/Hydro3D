!!=======================================================================
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

		if (out_pt(l).eq..false.) then
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
