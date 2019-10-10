!##########################################################################
      MODULE var_rough
!##########################################################################

        real hd,sigma
        integer :: maxk,iblock,rough_blockno,rough_proc
        integer, allocatable,dimension (:) :: rough_block

        type multi_rough
        real d50
        real, pointer, dimension(:,:,:) :: rough,irough
        real, pointer, dimension(:,:)   :: zbp,z_rough
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
      real :: d50_dummy
      CHARACTER*31 :: gridfile

      gridfile='in_rough_info.cin'
      open (unit=1, file=gridfile)      
      read (1,*) d50_dummy,hd,sigma,rough_blockno,maxk

      allocate (rough_block(1:rough_blockno),rough_dom(nbp))

      DO N=1,rough_blockno    
       read (1,*) rough_block(N)
      END DO          

      DO ib = 1, nbp
        ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k

        allocate(rough_dom(ib)%rough(ni,nj,nk))
        allocate(rough_dom(ib)%irough(ni,nj,nk))
        allocate(rough_dom(ib)%zbp(ni,nj)) 
        allocate(rough_dom(ib)%z_rough(ni,nj))
        rough_dom(ib)%d50 = d50_dummy

        DO L=1,rough_blockno
          IF (dom_id(ib).eq.rough_block(L)) CALL Roughness_Function(ib)
        END DO      

      END DO
    
      RETURN
      END SUBROUTINE
!##########################################################################	
	Subroutine Roughness_Function(ib)
!##########################################################################
      use vars
      use mpi
      use var_rough
      use multidata
      implicit none
      integer, intent(in) :: ib
      integer :: ni, nj,nk
      real aa(301), A(301),b
      real    :: xistep,yjstep,zdelta,zbav, rms
      real    :: random_number_normal2,maxelev,minelev,maxelev2,minelev2
      integer :: i,j,k,istep,jstep,ii1,ii2,ii,icount!,ib
      character (LEN=80)   :: tecfile, cdummy
      character (LEN=4)    :: char_block
      real 	  :: random, xicount

      ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k

     	rough_dom(ib)%d50=rough_dom(ib)%d50/1000.0/hd

	xistep= (rough_dom(ib)%d50/dom(ib)%dx)
	yjstep= (rough_dom(ib)%d50/dom(ib)%dy)

      istep=NINT(xistep) ; jstep=NINT(yjstep)

	if (istep.eq.0) istep=1
      if (jstep.eq.0) jstep=1

      istep=1
      jstep=1

      call RANDOM_SEED
      
      rough_dom(ib)%z_rough=10.
      maxelev =-100.
      minelev = 100.

	do i = 1,ni, istep
        do j = 1, nj, jstep
	rough_dom(ib)%z_rough(i,j)= random_number_normal2(0.0,sigma)
     &                                *rough_dom(ib)%d50
	maxelev=MAX(maxelev,rough_dom(ib)%z_rough(i,j))
	minelev=MIN(minelev,rough_dom(ib)%z_rough(i,j))
        enddo
	end do

!  shift bed upwards so that the troughs are at zero
      maxelev2=-100.
      minelev2=100.

!  modify boundaries
      if (dom(ib)%iprev.lt.0) then
        do i = 1,pl
          do j = 1, nj
            do k = pl,1,-1
              if (rough_dom(ib)%z_rough(i,j).gt.5.0) Then
        rough_dom(ib)%z_rough(i,j)=rough_dom(ib)%z_rough(k+1,j)
              end if
            end do
          enddo
        end do
      end if 

!  periodic boundaries
      do i = ni,ni
        do j = 1, nj
!           z_rough(i,j)=z_rough(1,j)
        enddo
      end do

      if (dom(ib)%jprev.lt.0) then
        do i = 1,ni
          do j = pl, 1, -1
            if (rough_dom(ib)%z_rough(i,j).gt.5.0) 
     & rough_dom(ib)%z_rough(i,j)=rough_dom(ib)%z_rough(i,j+1)
          enddo
        end do    
      end if

      if (dom(ib)%jnext.lt.0) then
        do i = 1,ni
          do j =  nj-pl+1,nj
            if (rough_dom(ib)%z_rough(i,j).gt.5.0) 
     & rough_dom(ib)%z_rough(i,j)=rough_dom(ib)%z_rough(i,j-1)
          enddo
	  end do
      end if

      zbav=0.0

      write(char_block,'(I4)') dom_id(ib)
      tecfile='maxk'//TRIM(ADJUSTL(char_block))//'.dat'
      open (UNIT=688,file=TRIM(tecfile))
   
      do i = 1,ni
        do j = 1, nj
          rough_dom(ib)%z_rough(i,j)= rough_dom(ib)%z_rough(i,j)+
     & ABS(minelev)       
          b=rough_dom(ib)%z_rough(i,j)/rough_dom(ib)%d50
          WRITE (688,*) i,j,'b:',b
          maxelev2=MAX(maxelev2,rough_dom(ib)%z_rough(i,j))
          minelev2=MIN(minelev2,rough_dom(ib)%z_rough(i,j))
          zbav=zbav+rough_dom(ib)%z_rough(i,j)
        end do
      end do

      zbav=zbav/(ni*nj) 

      WRITE (688,*) 'mycpu#:',myrank,'myblock#:',dom_id(ib)
      WRITE (688,*) '==================================='
      WRITE (688,*) 'maximum bed elevation:  ', maxelev2
      WRITE (688,*) 'minimum bed elevation:  ', minelev2
      WRITE (688,*) 'mean bed elevation:  ',zbav
      WRITE (688,*) 'roughness diameter:  ',rough_dom(ib)%d50       

C determine roughness geometry function
      rms=0.0
      do i = 1,ni
        do j = 1, nj         
	    rough_dom(ib)%zbp(i,j)=(rough_dom(ib)%z_rough(i,j)-zbav)
          rms=rms+rough_dom(ib)%zbp(i,j)**2 
  	  enddo
      end do
      rms=rms/(ni-nj)

      do i = 1,ni
        do j = 1, nj
          rough_dom(ib)%zbp(i,j)=rough_dom(ib)%zbp(i,j)/
     & (sigma*rough_dom(ib)%d50)
!	write (78,69) 'mycpu#:',myrank,i,j,rough_dom(ib)%zbp(i,j)
        enddo
      end do

      aa(1)=-3.0
      do ii=2,301
	  aa(ii)=aa(ii-1)+0.02      
	  xicount=0
	  do i =1, ni
	    do j=1,nj
            if (rough_dom(ib)%zbp(i,j).lt. aa(ii)) xicount=xicount+1
	    end do
	  end do
	  A(ii)=xicount/(ni*nj)
      end do

   88 FORMAT(a,i4,2F15.6)
   69 FORMAT(a,i4,i4,i4,1F15.6)

! determine porosity 
      rough_dom(ib)%rough=0.0
	rough_dom(ib)%irough=0
	maxk=0

      do i=1,ni
	  do j=1,nj
          do k=1,nk

!	      zdelta=dom(ib)%zc(k)-dom(ib)%zc(1)

!	      if (zdelta.lt.rough_dom(ib)%z_rough(i,j)) then
	      if (dom(ib)%zc(k).lt.rough_dom(ib)%z_rough(i,j)) then
	        rough_dom(ib)%rough(i,j,k)=1.0
	        rough_dom(ib)%irough(i,j,k)=1
	        maxk=MAX(maxk,k)
	      end if

          end do
  	  end do
      end do

       
      write (688,*) 'maximum value of k is:  ', maxk
      close (688)

      write(char_block,'(I4)') dom_id(ib)
      tecfile='rough'//TRIM(ADJUSTL(char_block))//'.plt'

      open (UNIT=455,file=TRIM(tecfile))		     
      WRITE (455,*) 'TITLE = ',' Roughness Function'
      WRITE (455,*) 'VARIABLES = "X", "Y", "Z", "ROUGH", "ZR" '

      write (455,*)'zone ', ' i=',ni,', ',
     &  ' j=',nj,', k= ',nk,' f=point'

      do k=1,nk
        do j=1,nj
          do i=1,ni
            write (455,'(5e14.6)') dom(ib)%xc(i),
     &dom(ib)%yc(j),dom(ib)%zc(k),rough_dom(ib)%rough(i,j,k)
     &,rough_dom(ib)%z_rough(i,j)
          end do
        end do
      end do

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

         IF (dom_id(ib).eq.rough_block(L))  THEN
       
        do k=dom(ib)%ksp,maxk+pl-1    !2,maxk
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
       integer i,j,k,ib
       integer :: ni, nj,nk,L
       real dummy
       character (LEN=80)   :: tecfile,charac
       character (LEN=4)    :: char_block

       DO ib = 1, nbp

      ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k

       DO L=1,rough_blockno

       IF (dom_id(ib).eq.rough_block(L))  THEN

       write(char_block,'(I4)') dom_id(ib)
       tecfile='rough'//TRIM(ADJUSTL(char_block))//'.plt'

       open (UNIT=455,file=TRIM(tecfile))
       do i=1,3
         READ(455,*) charac
       end do

       do k=1,nk
         do j=1,nj
           do i=1,ni
             READ (455,'(5e14.6)') dummy,dummy,dummy,
     &rough_dom(ib)%rough(i,j,k),rough_dom(ib)%z_rough(i,j)
           end do
         end do
       end do

       END IF
       END DO
      
       END DO	

       RETURN
       END SUBROUTINE         

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
