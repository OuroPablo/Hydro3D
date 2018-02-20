!######################################################################
      subroutine bound_LSM(op)
!######################################################################
      use vars
      use LSM
      use multidata

      implicit none
      integer :: I,J,K,np,nq,nr,op,ib,ipl
      double precision, pointer, dimension(:,:,:) :: FI

        do ib=1,nbp

       select case (op)
        Case (14)  
          fi => dom(ib)%PHI
        Case (15)  
          fi => dom(ib)%PHI_REINIT
        Case (16)  
          fi => dom(ib)%phi_new
        Case (17)  
          fi => dom(ib)%phi_init
        Case (18)  
          fi => dom(ib)%dens
        Case (19)  
          fi => dom(ib)%mu
       end select 

          np = dom(ib)%ttc_i; nq = dom(ib)%ttc_j; nr = dom(ib)%ttc_k  
!
!..... Boundary Conditions for PHI, zero gradient..
!
!=== EAST ===
      IF (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) THEN   
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl
        DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl
	  DO ipl = 1,pl
	    FI (dom(ib)%iep+ipl,j,k) = FI(dom(ib)%iep+ipl-1,j,k)

!Fix the level set funtion at the outlet of the domain:
	  if(op.le.17) then ! level set functions
         if (dom(ib)%bc_east.eq.2 .or.dom(ib)%bc_east.eq.21) then
                 if (dom(ib)%zc(k).lt.length)   then
	    FI(dom(ib)%isp-ipl,j,k) = 1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).gt.length)   then
	    FI(dom(ib)%isp-ipl,j,k) = -1.0*abs(dom(ib)%zc(k)-length)
                 else if (dom(ib)%zc(k).eq.length)  then
	    FI(dom(ib)%isp-ipl,j,k) = 0.0
                 end if  
  	  endif; endif

	  ENDDO
        END DO  
       END DO
      END IF
!
!=== WEST ===
      IF (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) THEN
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl
        DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl 
	  DO ipl = 1,pl
	    FI(dom(ib)%isp-ipl,j,k) = FI(dom(ib)%isp-ipl+1,j,k)
	  ENDDO
        END DO
       END DO
      END IF
!
!=== BOTTOM ===
      IF (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) THEN 
       DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl
	  DO ipl = 1,pl
	    FI (i,j,dom(ib)%ksp-ipl) = FI(i,j,dom(ib)%ksp-ipl+1)
	  ENDDO
        END DO
       END DO
      END IF
!
!=== TOP ===   
      IF (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) THEN 
       DO j = dom(ib)%jsp-pl, dom(ib)%jep+pl
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl  
	  DO ipl = 1,pl
	    FI (i,j,dom(ib)%kep+ipl) = FI(i,j,dom(ib)%kep+ipl-1)
	  ENDDO
        END DO
       END DO   
      END IF
!
!=== SOUTH ===
      IF (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) THEN
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl   
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl    
	  DO ipl = 1,pl
	    FI (i,dom(ib)%jsp-ipl,k) = FI(i,dom(ib)%jsp-ipl+1,k)
	  ENDDO
        END DO
       END DO
      END IF
!
!=== NORTH ===  
      IF (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) THEN 
       DO k = dom(ib)%ksp-pl, dom(ib)%kep+pl   
        DO i = dom(ib)%isp-pl, dom(ib)%iep+pl     
	  DO ipl = 0,pl
	    FI (i,dom(ib)%jep+ipl,k) = FI(i,dom(ib)%jep+ipl-1,k)
	  ENDDO
        END DO
       END DO   
      END IF

      end do

	call exchange(op) 

      RETURN
      end subroutine bound_LSM
!#####################################################################
