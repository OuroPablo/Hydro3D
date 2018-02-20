!=======================================================================
!    			LAGRANGIAN PARTICLE TRACKING
!			Bruño Fraga Bugallo
!			Cardiff 2013-2015
!=======================================================================
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

       open(10,file='LPT.cin')
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

       elseif (LRESTART.eq..false.) then
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



