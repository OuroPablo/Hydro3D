!##########################################################################
      MODULE var_rough
!##########################################################################

        real hd,sigma,d50
        integer :: maxk,iblock,rough_blockno,rough_proc,iall,jall
        integer, allocatable,dimension (:) :: rough_block
        real,allocatable,dimension (:,:) :: zrough_all
        real,allocatable, dimension(:,:):: zbp
        
        type multi_rough
       
        real, pointer, dimension(:,:,:) :: rough
        real, pointer, dimension(:,:)   :: z_rough
 
        end type multi_rough

        type (multi_rough), pointer, dimension(:) :: rough_dom

	END MODULE
!##########################################################################
      SUBROUTINE INIT_ROUGH
!##########################################################################
      use vars
      use mpi
      use var_rough
      use multidata
      implicit none
      integer :: ib,proc_no,L,N,ni,nj,nk
      integer ::i,j,k,nx,ny,istart,jstart
      CHARACTER*31 :: gridfile
      character (LEN=80)   :: tecfile
      character (LEN=4)    :: char_block

      gridfile='rough_info.cin'
      open (unit=1, file=gridfile)  
      read (1,*)     
      read (1,*) d50
      read (1,*) hd
      read (1,*) sigma

! generate the whole bed
      iall=(xen-xst)/g_dx+2*pl
      jall=(yen-yst)/g_dy+2*pl
      nx=(iall-2*pl)/idom !for each subdomain
      ny=(jall-2*pl)/jdom
      rough_blockno=idom*jdom !only for the bottom layer
       
      allocate(zrough_all(iall,jall))
      allocate(zbp(iall,jall))

      if(myrank.eq.0) then
 
        CALL Roughness_Function

!write the total roughness files      
      open (UNIT=460,file='totalrough.plt')
      WRITE (460,*) 'TITLE = ',' Roughness Function'
      WRITE (460,*) 'VARIABLES = "X", "Y", "Zb", "Zbed" '

      write (460,*)'zone ', ' i=',iall-2*pl,', ',
     &  ' j=',jall-2*pl,' f=point'		     
      do j=pl+1,jall-pl
        do i=pl+1,iall-pl    
       write (460,'(4e14.6)') g_dx*(0.5+real(i-pl-1)),
     & g_dy*(0.5+real(j-pl-1)),
     & zrough_all(i,j),zrough_all(i,j)
        end do
      end do    
      CLOSE(460)
       
!write the seperate roughness files      
      DO L=1,rough_blockno !in bottom layer domains
      i=mod((L-1),idom)
      j=int((L-1)/idom)
      
      istart=i*nx+1
      jstart=j*ny+1
      
      write(char_block,'(I4)') (L-1)
      tecfile='rough'//TRIM(ADJUSTL(char_block))//'.dat'
      open (UNIT=455,file=TRIM(tecfile))		     
      do j=jstart,jstart+2*pl+ny-1
        do i=istart,istart+2*pl+nx-1     
          write (455,*) zrough_all(i,j)
        end do
      end do    
      CLOSE(455)

      END DO !L  
      end if !myrank 0

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       allocate (rough_dom(nbp))

!ALL processor read the roughness files       
      DO ib = 1, nbp
       ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j
       

       allocate (rough_dom(ib)%z_rough(ni,nj))
	rough_dom(ib)%z_rough=0.d0   
    
       N=mod(dom_id(ib),rough_blockno)

      write(char_block,'(I4)') N
      tecfile='rough'//TRIM(ADJUSTL(char_block))//'.dat'
      open (UNIT=455,file=TRIM(tecfile))
        do j=1,nj      
        do i=1,ni
        read(455,*)rough_dom(ib)%z_rough(i,j)
!        WRITE(6,*) rough_dom(ib)%z_rough(i,j)
        end do
        end do    
      CLOSE(455)
      END DO !ib

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!determine parameter rough 1 or 0 &maxk
       DO ib = 1, nbp
       ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
       allocate(rough_dom(ib)%rough(ni,nj,nk))
       rough_dom(ib)%rough=0.0
       maxk=0
      
       DO L=1,rough_blockno
         IF (dom_id(ib).eq.(L-1))  THEN
       do j=1,nj      
       do i=1,ni
       do k=1,nk

         if (dom(ib)%zc(k).lt.rough_dom(ib)%z_rough(i,j)) then
          
	   rough_dom(ib)%rough(i,j,k)=1.0
           maxk=MAX(maxk,k)
           
         end if
       end do
       end do
       end do 

          END IF
       end do !	L

       end do !ib


      RETURN
      END SUBROUTINE
!##########################################################################	
	Subroutine Roughness_Function
!##########################################################################
      use vars
      use mpi
      use var_rough
      use multidata
      implicit none
    

      integer :: ni, nj,nk
 
      real :: b
      real    :: xistep,yjstep,zdelta,zbav, rms
      real    :: random_number_normal2,maxelev,minelev,maxelev2,minelev2
      integer :: i,j,k,l,istep,jstep,ii1,ii2,ii,icount
      character (LEN=80)   :: tecfile, cdummy
      character (LEN=4)    :: char_block
      real random, xicount

     
      d50=d50/hd

	xistep= (d50/g_dx)
	yjstep= (d50/g_dy)

        istep=NINT(xistep)
	jstep=NINT(yjstep)

	if (istep.eq.0) istep=1
        if (jstep.eq.0) jstep=1

      istep=1
      jstep=1

      call RANDOM_SEED
    
      zrough_all=0.
      maxelev =-100.
      minelev = 100.
 
	do i = pl+1,iall-pl, istep
        do j = pl+1,jall-pl, jstep
	zrough_all(i,j)=random_number_normal2(0.0,sigma)*d50
                           
	maxelev=MAX(maxelev,zrough_all(i,j))
	minelev=MIN(minelev,zrough_all(i,j))
        
        do l=1,istep         !keep those cube same bed
        do k=1,jstep    
        if(((i+l-1).le.iall).and.((j+k-1).le.jall)) then                  
        zrough_all(i+l-1,j+k-1)=zrough_all(i,j)
        end if
        end do 
        end do           

        end do
	end do

!  shift bed upwards so that the troughs are at zero
      maxelev2=-100.
      minelev2=100.
      zbav=0.0   

      do i = pl+1,iall-pl
        do j = pl+1,jall-pl
          zrough_all(i,j)= zrough_all(i,j)+ABS(minelev)+0.5*g_dz !move 0.5deltaz upwards        
          maxelev2=MAX(maxelev2,zrough_all(i,j))
          minelev2=MIN(minelev2,zrough_all(i,j))
          zbav=zbav+zrough_all(i,j)
        end do
      end do

      zbav=zbav/(iall-2*pl)/(jall-2*pl)
 
! calculate corresponding values
      rms=0.0
      do i = pl+1,iall-pl
        do j = pl+1,jall-pl       
	  zbp(i,j)=(zrough_all(i,j)-zbav)
          rms=rms+zbp(i,j)**2 
  	end do
      end do
      rms=sqrt(rms/(iall-2*pl)/(jall-2*pl))

       do i = 1+pl,ni-pl
        do j = 1+pl,nj-pl
          zbp(i,j)=zbp(i,j)/(sigma*d50)
         !WRITE (688,*) i,j,'zbp:',zbp(i,j)
        end do
       end do

! determine ghost cell
      do i=1,pl
      do j = pl+1,jall-pl
      zrough_all(i,j)=zrough_all(i+iall-2*pl,j)
      end do
      end do

      do i=iall-pl+1,iall
      do j = pl+1,jall-pl
      zrough_all(i,j)=zrough_all(i-iall+2*pl,j)
      end do
      end do

      do j=1,pl
      do i = pl+1,iall-pl
      zrough_all(i,j)=zrough_all(i,j+jall-2*pl)
      end do
      end do

      do j=jall-pl+1,jall
      do i = pl+1,iall-pl
      zrough_all(i,j)=zrough_all(i,j-jall+2*pl)
      end do
      end do

! open write files for recording roughness information
      open (UNIT=688,file='info_rough.dat')
      WRITE (688,*) '==================================='
      WRITE (688,*) 'roughness diameter:  ',d50      
      WRITE (688,*) 'maximum bed elevation:  ', maxelev2
      WRITE (688,*) 'minimum bed elevation:  ', minelev2
      WRITE (688,*) 'mean bed elevation:  ',zbav 
      write (688,*) 'rms of bed elevation:  ', rms 

      close (688)

      RETURN
      END Subroutine

!##########################################################################	
      Subroutine Rough_velocity
!##########################################################################
      use vars
      use mpi
      use var_rough     
      use multidata

      IMPLICIT NONE
      INTEGER :: I,J,K,L,ib,nip,njp,N
      CHARACTER*31 :: gridfile

      DO ib = 1, nbp

       nip=dom(ib)%ttc_i; njp=dom(ib)%ttc_j

       DO L=1,rough_blockno

         IF (dom_id(ib).eq.(L-1))  THEN
       
        do k=dom(ib)%ksp,maxk           !2,maxk
          do j=dom(ib)%jsp,dom(ib)%jep   !2,njp-1
            do i=dom(ib)%isp,dom(ib)%iep    !2,nip-1
              if (rough_dom(ib)%rough(i,j,k).eq.1.0) then
                dom(ib)%ustar(i-1,j,k)=0.0
                dom(ib)%vstar(i,j-1,k)=0.0
                dom(ib)%wstar(i,j,k-1)=0.0
                dom(ib)%ustar(i,j,k)=0.
                dom(ib)%vstar(i,j,k)=0.
                dom(ib)%wstar(i,j,k)=0.
              end if
            end do
          end do
        end do

         END IF

       END DO

      END DO


      RETURN
      END SUBROUTINE
!##########################################################################
      Subroutine Rough_restart
!##########################################################################

       use vars
       use mpi
       use var_rough
       use multidata

       IMPLICIT NONE
       integer i,j,k,ib,n
       integer :: ni, nj,nk,L
       real dummy
       character (LEN=80)   :: tecfile,charac
       character (LEN=4)    :: char_block
      
      rough_blockno=idom*jdom !only for the bottom layer
      allocate (rough_dom(nbp))

      DO ib = 1, nbp
       ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
       allocate (rough_dom(ib)%z_rough(ni,nj))
       rough_dom(ib)%z_rough=0.d0

       N=mod(dom_id(ib),rough_blockno)
       
      write(char_block,'(I4)') N
      tecfile='rough'//TRIM(ADJUSTL(char_block))//'.dat'
      open (UNIT=455,file=TRIM(tecfile))
        do j=1,nj      
        do i=1,ni
        read(455,*)rough_dom(ib)%z_rough(i,j) 
        !WRITE(6,*) rough_dom(ib)%z_rough(i,j)     
        end do
        end do    
      CLOSE(455)
      END DO !ib

!determine parameter rough 1 or 0 &maxk
       DO ib = 1, nbp
       ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
       allocate(rough_dom(ib)%rough(ni,nj,nk))
       rough_dom(ib)%rough=0.0
       maxk=0
      
       DO L=1,rough_blockno
         IF (dom_id(ib).eq.(L-1))  THEN
       do j=1,nj      
       do i=1,ni
       do k=1,nk
         if (dom(ib)%zc(k).lt.rough_dom(ib)%z_rough(i,j)) then
          
	   rough_dom(ib)%rough(i,j,k)=1.0
           maxk=MAX(maxk,k)
           
         end if
       end do
       end do
       end do 
       
          END IF
       end do !	L

       end do !ib

       RETURN
       END SUBROUTINE         
!##########################################################################
C******************ADD NOISE****************************************************
C*******************************************************************************

!-----------------------------------------------------------------------
      FUNCTION random_number_normal2(mean,sigma) RESULT( fn_val )
!         Generate random numbers
!         with a normal distribution with given mean and standard deviaton.
!
!         Generate a random normal deviate using the polar method.
!         Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!                    normal variables', Siam Rev., vol.6, 260-264, 1964.
!------------------------------------------------------------------------
      IMPLICIT NONE
      REAL  :: fn_val
      REAL  :: mean,sigma
!
!.... Local variables
!
      REAL            :: u, sum
      REAL, SAVE      :: v, sln
      LOGICAL, SAVE   :: second = .FALSE.
      REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

      IF (second) THEN
!
!...... If second, use the second random number generated on last call
!
        second = .false.
        fn_val = v*sln

      ELSE
!
!...... First call; generate a pair of random normals
!
        second = .true.
        DO
          CALL RANDOM_NUMBER( u )
          CALL RANDOM_NUMBER( v )
          u = SCALE( u, 1 ) - one
          v = SCALE( v, 1 ) - one
!
!.........vsmall added to prevent LOG(zero) / zero
!
          sum = u*u + v*v + vsmall
          IF(sum < one) EXIT
        END DO
        sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
        fn_val = u*sln
      END IF
!
!.....set mean and standard deviation
!
      fn_val = fn_val * sigma + mean
!
      RETURN
      END FUNCTION random_number_normal2

C*******************************************************************************
C*******************************************************************************
