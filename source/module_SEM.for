!#############################################################################
	module module_SEM
!#############################################################################
	SAVE
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: Vsem ,Usem
	DOUBLE PRECISION,DIMENSION(:,:),  ALLOCATABLE :: X_EDDY,EPSILO,MOLT
	DOUBLE PRECISION,DIMENSION(:,:,:),  ALLOCATABLE :: SIGMA
	DOUBLE PRECISION:: R(3,3),X_POINT(3),REYNOLDS(6),TEMP(3)
	DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE ::  Ksem
	CHARACTER*44 :: FILEGLOBAL
	INTEGER,ALLOCATABLE:: elemyst(:),elemyen(:),elemzst(:),elemzen(:)
	INTEGER,ALLOCATABLE:: iddom(:),ljdom(:),lkdom(:)
	end
