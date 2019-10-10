!#############################################################################
	module module_SEM
!#############################################################################
	SAVE
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE :: Vsem ,Usem
	DOUBLE PRECISION,DIMENSION(:,:),  ALLOCATABLE :: X_EDDY,EPSILO,MOLT
	DOUBLE PRECISION,DIMENSION(:,:,:),  ALLOCATABLE :: SIGMA
	DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE ::  Ksem
	CHARACTER*44 :: FILEGLOBAL
	INTEGER,ALLOCATABLE:: elemyst(:),elemyen(:),elemzst(:),elemzen(:)
	INTEGER,ALLOCATABLE:: iddom(:),ljdom(:),lkdom(:)
	end
