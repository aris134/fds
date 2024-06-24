PROGRAM LOOP3D
    ! This program is the MPI-OpenMP GPU version of the original OpenMP-based loop3d program
    ! Code modifications were made to enable the mini-app to run across multiple compute nodes
    USE OMP_LIB
    IMPLICIT NONE
    include 'mpif.h'
    ! MPI declarations
    INTEGER :: IERR, MY_RANK, NUM_TASKS, REQUEST

    ! Miscellaneous declarations
    INTEGER :: ISTEP
    INTEGER :: RK, NBR_RK
    INTEGER, PARAMETER :: HALO_DEPTH = 1
    INTEGER, PARAMETER :: NUM_TIME_STEPS = 1
    INTEGER, PARAMETER :: BOX_SIZE = 6
    INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
    INTEGER, PARAMETER :: IBAR = 1024, JBAR = 1024, KBAR =1024
    INTEGER, PARAMETER :: NEDGE = 12
    REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB
    REAL(EB) :: FLUX_AVGS(3)
    INTEGER :: MY_BBOX(6)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: HALO_SEND_SIZES, HALO_SEND_DISPLS, SEND_COORD_DISPLS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: HALO_RECV_SIZES, HALO_RECV_DISPLS, RECV_COORD_DISPLS
    INTEGER :: TOTAL_SEND_PTS, TOTAL_RECV_PTS
    INTEGER :: NUM_NBRS ! number of ranks that neighbor this one
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NBR_RANKS ! list of ranks that neighbor this one
    
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: GX, GY, GZ, RHO_0
    REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:) :: X, Y, Z, RDXN, RDYN, RDZN, RDX, RDY, RDZ
    REAL(EB), POINTER, DIMENSION(:) :: XP, YP, ZP, RDXNP, RDYNP, RDZNP, RDXP, RDYP, RDZP
    REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: U, V, W, D, RHO, WORK1, WORK2, WORK3, WORK4, WORK5, WORK6
    REAL(EB), POINTER, DIMENSION(:,:,:) :: DP, RHOP, UU, VV, WW, OMY, OMX, OMZ, TXZ, TXY, TYZ
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: OMX_SEND, OMY_SEND, OMZ_SEND, TXZ_SEND, TXY_SEND, TYZ_SEND, U_SEND, V_SEND, W_SEND
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: OMX_RECV, OMY_RECV, OMZ_RECV, TXZ_RECV, TXY_RECV, TYZ_RECV, U_RECV, V_RECV, W_RECV
    INTEGER, ALLOCATABLE, DIMENSION(:) :: SEND_COORDS, RECV_COORDS
    REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU, FVX, FVY, FVZ
    INTEGER, ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: CELL_INDEX
    INTEGER, POINTER, DIMENSION(:,:,:) :: CC ! pointer to CELL_INDEX
    INTEGER, ALLOCATABLE, DIMENSION (:) :: CELL_INDEX_SEND, CELL_INDEX_RECV
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_TASK_BBOXES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_TASK_EXT_BBOXES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NBR_EXT_BBOXES
    REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
                DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
                DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
                VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
                RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,T_SIM_START,T_SIM_END,&
                LOOP_TIME,MIN_LOOP_TIME,MAX_LOOP_TIME,AVG_LOOP_TIME,&
                T_SETUP_START,T_SETUP_END,SETUP_TIME,MAX_SETUP_TIME,TOTAL_TIME
    INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,ICMAX,IE,MAX_EDGE,GLOB_MAX_IC,GLOB_MAX_EDGE,NT
    INTEGER :: GLOB_MIN_I, GLOB_MIN_J, GLOB_MIN_K, GMINIP1, GMINJP1, GMINKP1, GLOB_MAX_I, GLOB_MAX_J, GLOB_MAX_K
    INTEGER :: GLOBAL_MIN(3)
    INTEGER :: GLOBAL_MAX(3)
    INTEGER :: LOC_MIN_I, LOC_MIN_J, LOC_MIN_K, LOC_MAX_I, LOC_MAX_J, LOC_MAX_K, LOC_MAX_I_P1, LOC_MAX_J_P1, LOC_MAX_K_P1, &
               LOC_MIN_I_P1, LOC_MIN_J_P1, LOC_MIN_K_P1
    INTEGER :: LOC_MIN_I_M1, LOC_MIN_J_M1, LOC_MIN_K_M1
    INTEGER :: LOC_MIN_HALO_I, LOC_MIN_HALO_J, LOC_MIN_HALO_K, &
               LOC_MAX_HALO_I, LOC_MAX_HALO_J, LOC_MAX_HALO_K
    CHARACTER(LEN=50) :: FILENAME

TYPE CELL_TYPE
   INTEGER :: EDGE_INDEX(NEDGE)=0
END TYPE CELL_TYPE
TYPE(CELL_TYPE), ALLOCATABLE, DIMENSION(:) :: CELL

TYPE EDGE_TYPE
   REAL(EB) :: OMEGA(-2:2)=-1.E6_EB
   REAL(EB) :: TAU(-2:2)=-1.E6_EB
END TYPE EDGE_TYPE
TYPE(EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE
    
    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUM_TASKS, IERR)

    T_SETUP_START = MPI_WTIME()

    ! Write out Starting:
    !$OMP PARALLEL
    !$OMP MASTER
    !$ NT = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP BARRIER
    !$OMP END PARALLEL
    
    IF (MY_RANK==0) THEN
        WRITE(FILENAME,'(A,I4,A,I2,A,I4,A)') 'loop3d_',IBAR,'GRID_',NT,'THR_',INT(NUM_TASKS),'RANKS_GPU_MPI_V1.txt'
        WRITE(*,*) 'Starting Loop3D, out file: ',TRIM(FILENAME)
        WRITE(*,*) 'Number of detected devices: ',OMP_GET_NUM_DEVICES()
        OPEN(UNIT=10,FILE=TRIM(FILENAME),STATUS='UNKNOWN')
        WRITE(10,*) 'Number of devices=',OMP_GET_NUM_DEVICES()
        WRITE(10,*) 'Starting Loop3D'
        WRITE(10,*) 'IBAR=',IBAR,' JBAR=',JBAR,' KBAR=',KBAR,' OMP_NUM_THREADS=',NT
    END IF

    GLOB_MIN_I = 0
    GLOB_MIN_J = 0
    GLOB_MIN_K = 0
    GLOBAL_MIN(1) = GLOB_MIN_I
    GLOBAL_MIN(2) = GLOB_MIN_J
    GLOBAL_MIN(3) = GLOB_MIN_K
    GMINIP1 = GLOB_MIN_I + 1
    GMINJP1 = GLOB_MIN_J + 1
    GMINKP1 = GLOB_MIN_K + 1
    GLOB_MAX_I = IBAR
    GLOB_MAX_J = JBAR
    GLOB_MAX_K = KBAR
    GLOBAL_MAX(1) = GLOB_MAX_I
    GLOBAL_MAX(2) = GLOB_MAX_J
    GLOBAL_MAX(3) = GLOB_MAX_K

    ! Decompose the domain
    ALLOCATE(ALL_TASK_BBOXES(NUM_TASKS * BOX_SIZE)) ! AoS format: ..., imin_p, jmin_p, kmin_p, imax_p, jmax_p, kmax_p, imin_(p+1), ...
    ALLOCATE(ALL_TASK_EXT_BBOXES(NUM_TASKS * BOX_SIZE)) ! all task bboxes extended by halo depth
    CALL DECOMPOSE_DOMAIN(ALL_TASK_BBOXES, GLOB_MIN_I, GLOB_MIN_J, GLOB_MIN_K, GLOB_MAX_I, GLOB_MAX_J, GLOB_MAX_K, BOX_SIZE, NUM_TASKS, MY_RANK)
    ! extract the local task bounds from ALL_TASK_BBOXES. These are INCLUSIVE endpoints of the rank's local domain bounding box
    LOC_MIN_I = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE)
    LOC_MIN_J = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 1)
    LOC_MIN_K = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 2)
    LOC_MAX_I = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 3)
    LOC_MAX_J = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 4)
    LOC_MAX_K = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 5)
    MY_BBOX(1) = LOC_MIN_I
    MY_BBOX(2) = LOC_MIN_J
    MY_BBOX(3) = LOC_MIN_K
    MY_BBOX(4) = LOC_MAX_I
    MY_BBOX(5) = LOC_MAX_J
    MY_BBOX(6) = LOC_MAX_K
    ! Debugging: print out rank's bounding box coordinates
    ! WRITE(*,*) 'Rank=',MY_RANK,' LOC_MIN_I=',LOC_MIN_I,' LOC_MIN_J=',LOC_MIN_J,' LOC_MIN_K=',LOC_MIN_K,' LOC_MAX_I=',LOC_MAX_I,' LOC_MAX_J=',LOC_MAX_J,' LOC_MAX_K=',LOC_MAX_K

    ! Determine the neighboring ranks, and neighboring rank bounding boxes...
    CALL DETERMINE_NUM_NEIGHBORS(MY_RANK, GLOBAL_MIN, GLOBAL_MAX, MY_BBOX, NUM_TASKS, NUM_NBRS, ALL_TASK_BBOXES, ALL_TASK_EXT_BBOXES, HALO_DEPTH, BOX_SIZE)
    ALLOCATE(NBR_EXT_BBOXES(NUM_NBRS * BOX_SIZE))
    ALLOCATE(NBR_RANKS(NUM_NBRS))
    ALLOCATE(HALO_SEND_SIZES(NUM_NBRS))
    ALLOCATE(HALO_SEND_DISPLS(NUM_NBRS))
    ALLOCATE(HALO_RECV_SIZES(NUM_NBRS))
    ALLOCATE(HALO_RECV_DISPLS(NUM_NBRS))
    CALL GET_NBR_BBOXES(MY_RANK, MY_BBOX, ALL_TASK_EXT_BBOXES, NBR_EXT_BBOXES, NBR_RANKS, NUM_TASKS, BOX_SIZE)

    ! Debugging: print nbr extended bboxes
    !DO RK=1,NUM_NBRS
    !    NBR_RK = NBR_RANKS(RK)
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_RANK=',NBR_RK
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_EXT_BBOX(',(RK-1) * BOX_SIZE,')=',NBR_EXT_BBOXES((RK-1) * BOX_SIZE)
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_EXT_BBOX(',(RK-1) * BOX_SIZE+1,')=',NBR_EXT_BBOXES((RK-1) * BOX_SIZE+1)
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_EXT_BBOX(',(RK-1) * BOX_SIZE+2,')=',NBR_EXT_BBOXES((RK-1) * BOX_SIZE+2)
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_EXT_BBOX(',(RK-1) * BOX_SIZE+3,')=',NBR_EXT_BBOXES((RK-1) * BOX_SIZE+3)
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_EXT_BBOX(',(RK-1) * BOX_SIZE+4,')=',NBR_EXT_BBOXES((RK-1) * BOX_SIZE+4)
    !    WRITE(*,*) 'MY_RANK=',MY_RANK,' NBR_EXT_BBOX(',(RK-1) * BOX_SIZE+5,')=',NBR_EXT_BBOXES((RK-1) * BOX_SIZE+5)
    !END DO

    ! Set the halo maximum values
    IF (LOC_MAX_I + HALO_DEPTH <= GLOB_MAX_I) THEN
        LOC_MAX_HALO_I = LOC_MAX_I + HALO_DEPTH
    ELSE
        LOC_MAX_HALO_I = GLOB_MAX_I
    END IF

    IF (LOC_MAX_J + HALO_DEPTH <= GLOB_MAX_J) THEN
        LOC_MAX_HALO_J = LOC_MAX_J + HALO_DEPTH
    ELSE
        LOC_MAX_HALO_J = GLOB_MAX_J
    END IF

    IF (LOC_MAX_K + HALO_DEPTH <= GLOB_MAX_K) THEN
        LOC_MAX_HALO_K = LOC_MAX_K + HALO_DEPTH
    ELSE
        LOC_MAX_HALO_K = GLOB_MAX_K
    END IF

    ! Set the halo minimum values
    IF (LOC_MIN_I - HALO_DEPTH >= GLOB_MIN_I) THEN
        LOC_MIN_HALO_I = LOC_MIN_I - HALO_DEPTH
    ELSE
        LOC_MIN_HALO_I = GLOB_MIN_I
    END IF

    IF (LOC_MIN_J - HALO_DEPTH >= GLOB_MIN_J) THEN
        LOC_MIN_HALO_J = LOC_MIN_J - HALO_DEPTH
    ELSE
        LOC_MIN_HALO_J = GLOB_MIN_J
    END IF

    IF (LOC_MIN_K - HALO_DEPTH >= GLOB_MIN_K) THEN
        LOC_MIN_HALO_K = LOC_MIN_K - HALO_DEPTH
    ELSE
        LOC_MIN_HALO_K = GLOB_MIN_K
    END IF

    ! Set +1 offsets
    LOC_MAX_I_P1 = LOC_MAX_I + 1
    LOC_MAX_J_P1 = LOC_MAX_J + 1
    LOC_MAX_K_P1 = LOC_MAX_K + 1

    IF (LOC_MIN_I == GLOB_MIN_I) THEN
        LOC_MIN_I_P1 = LOC_MIN_I + 1
    ELSE
        LOC_MIN_I_P1 = LOC_MIN_I
    END IF

    IF (LOC_MIN_J == GLOB_MIN_J) THEN
        LOC_MIN_J_P1 = LOC_MIN_J + 1
    ELSE
        LOC_MIN_J_P1 = LOC_MIN_J
    END IF

    IF (LOC_MIN_K == GLOB_MIN_K) THEN
        LOC_MIN_K_P1 = LOC_MIN_K + 1
    ELSE
        LOC_MIN_K_P1 = LOC_MIN_K
    END IF

    ! Pre-setup MPI comm
    CALL DETERMINE_HALO_SEND_SIZES(HALO_SEND_SIZES, TOTAL_SEND_PTS, NBR_EXT_BBOXES, NUM_NBRS, BOX_SIZE)
    ! Exchange message sizes
    CALL EXCHANGE_SIZES(HALO_SEND_SIZES, HALO_RECV_SIZES, NUM_NBRS)
    CALL COMPUTE_TOTAL_RECV_PTS(TOTAL_RECV_PTS, HALO_RECV_SIZES, NUM_NBRS)
    ! Fill the halo send and receive displacements
    CALL DETERMINE_HALO_DISPLACEMENTS(HALO_SEND_SIZES, HALO_SEND_DISPLS, NUM_NBRS)
    CALL DETERMINE_HALO_DISPLACEMENTS(HALO_RECV_SIZES, HALO_RECV_DISPLS, NUM_NBRS)

    ALLOCATE(SEND_COORDS(3*TOTAL_SEND_PTS)) ! AoS : I0,J0,K0,I1,J1,K1,...
    ALLOCATE(RECV_COORDS(3*TOTAL_RECV_PTS))

    !WRITE(*,*) 'RANK=',MY_RANK,',LOC_MIN_HALO_I=',LOC_MIN_HALO_I,' LOC_MIN_HALO_J=',LOC_MIN_HALO_J,' LOC_MIN_HALO_K=',LOC_MIN_HALO_K
    !WRITE(*,*) 'RANK=',MY_RANK,',LOC_MAX_HALO_I=',LOC_MAX_HALO_I,' LOC_MAX_HALO_J=',LOC_MAX_HALO_J,' LOC_MAX_HALO_K=',LOC_MAX_HALO_K
    !WRITE(*,*) 'RANK=',MY_RANK,',LOC_MIN_I_P1=',LOC_MIN_I_P1,' LOC_MIN_J_P1=',LOC_MIN_J_P1,' LOC_MIN_K_P1=',LOC_MIN_K_P1

    ! Allocate vars in CPU:
    ALLOCATE(X(LOC_MIN_HALO_I:LOC_MAX_HALO_I), Y(LOC_MIN_HALO_J:LOC_MAX_HALO_J), Z(LOC_MIN_HALO_K:LOC_MAX_HALO_K))
    ALLOCATE(RDXN(LOC_MIN_HALO_I:LOC_MAX_HALO_I), RDYN(LOC_MIN_HALO_J:LOC_MAX_HALO_J), RDZN(LOC_MIN_HALO_K:LOC_MAX_HALO_K))
    ALLOCATE(RDX(LOC_MIN_HALO_I:LOC_MAX_HALO_I), RDY(LOC_MIN_HALO_J:LOC_MAX_HALO_J), RDZ(LOC_MIN_HALO_K:LOC_MAX_HALO_K))
    ALLOCATE(U(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(V(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(W(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(MU(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(D(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(RHO(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(WORK1(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(WORK2(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(WORK3(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(WORK4(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(WORK5(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(WORK6(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(DP(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(RHOP(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(RHO_0(LOC_MIN_HALO_K:LOC_MAX_HALO_K),GX(LOC_MIN_HALO_I:LOC_MAX_HALO_I),GY(LOC_MIN_HALO_I:LOC_MAX_HALO_I))
    ALLOCATE(GZ(LOC_MIN_HALO_I:LOC_MAX_HALO_I))
    ALLOCATE(FVX(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(FVY(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(FVZ(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    ALLOCATE(CELL_INDEX(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
    
    ! Initialize (leave halo uninitialized as it will be addressed during MPI comm):
    DO I=LOC_MIN_I,LOC_MAX_I
       X(I)    = REAL(I,EB)
       RDXN(I) = 1.0_EB
       RDX(I)  = 1.0_EB
    ENDDO
    DO J=LOC_MIN_J,LOC_MAX_J
       Y(J)    = REAL(J,EB)
       RDYN(J) = 1.0_EB
       RDY(J)  = 1.0_EB
    ENDDO
    DO K=LOC_MIN_K,LOC_MAX_K
       Z(K)    = REAL(K,EB)
       RDZN(K) = 1.0_EB
       RDZ(K)  = 1.0_EB
    ENDDO

    XP => X
    YP => Y
    ZP => Z
    RDXNP => RDXN
    RDXP => RDX
    RDYNP => RDYN
    RDYP => RDY
    RDZNP => RDZN
    RDZP => RDZ

    ! 1. Exchange coordinates
    ALLOCATE(SEND_COORD_DISPLS(NUM_NBRS))
    ALLOCATE(RECV_COORD_DISPLS(NUM_NBRS))
    CALL COMPUTE_COORDINATE_DISPLS(NUM_NBRS,HALO_SEND_SIZES,HALO_RECV_SIZES,SEND_COORD_DISPLS,RECV_COORD_DISPLS)
    CALL EXCHANGE_COORDINATES(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_RECV_SIZES,&
                              SEND_COORDS,RECV_COORDS,XP,YP,ZP,RDXNP,RDXP,RDYNP,RDYP,RDZNP,RDZP,MY_RANK,NBR_RANKS,&
                              SEND_COORD_DISPLS,RECV_COORD_DISPLS)

    ! Cell Index, CELL and EDGE:
    ! Again, leave halo initialization for MPI comm
    CELL_INDEX = 0
    IC = 0
    DO K=LOC_MIN_K_P1,LOC_MAX_K
       DO J=LOC_MIN_J_P1,LOC_MAX_J
          DO I=LOC_MIN_I_P1,LOC_MAX_I
            IF( .NOT. (ANY( K==(/GMINKP1,GLOB_MAX_K/) ) .OR. ANY( J==(/GMINJP1,GLOB_MAX_J/) ) .OR. ANY( I==(/GMINIP1,GLOB_MAX_I/) )) ) CYCLE
             ! Obtain the same IC for a given I-J-K as in the serial case
             IC = IC_VAL(I,J,K,GLOB_MAX_I,GLOB_MAX_J,GLOB_MAX_K)
             CELL_INDEX(I,J,K) = IC
          ENDDO
       ENDDO
    ENDDO
    ! For now let's have an MPI_Allreduce here so that each rank has the same size CELL
    ! This is because we're essentially keeping a copy of the same data structure on each rank (for now)
    CALL MPI_ALLREDUCE(IC,GLOB_MAX_IC,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
    ALLOCATE(CELL(0:GLOB_MAX_IC))

    ! 2. Exchange cell index
    CC => CELL_INDEX
    ALLOCATE(CELL_INDEX_SEND(TOTAL_SEND_PTS))
    ALLOCATE(CELL_INDEX_RECV(TOTAL_RECV_PTS))

    CALL EXCHANGE_CELL_INDEX(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,&
                            HALO_RECV_DISPLS,RECV_COORDS,RECV_COORD_DISPLS,CELL_INDEX_SEND,CELL_INDEX_RECV,CC,MY_RANK)

    IC = 0
    MAX_EDGE=-1
    DO K=LOC_MIN_K_P1,LOC_MAX_K
       DO J=LOC_MIN_J_P1,LOC_MAX_J
          DO I=LOC_MIN_I_P1,LOC_MAX_I
            IF( .NOT. (ANY( K==(/GMINKP1,GLOB_MAX_K/) ) .OR. ANY( J==(/GMINJP1,GLOB_MAX_J/) ) .OR. ANY( I==(/GMINIP1,GLOB_MAX_I/) )) ) CYCLE
             ! Obtain the same IC for a given I-J-K as in the serial case
             IC = IC_VAL(I,J,K,GLOB_MAX_I,GLOB_MAX_J,GLOB_MAX_K)
             CELL(IC)%EDGE_INDEX(1)  = CELL_INDEX(I  ,J  ,K  ) + 1
             CELL(IC)%EDGE_INDEX(2)  = CELL_INDEX(I+1,J  ,K  ) + 1
             CELL(IC)%EDGE_INDEX(3)  = CELL_INDEX(I  ,J  ,K  ) + 2
             CELL(IC)%EDGE_INDEX(4)  = CELL_INDEX(I  ,J+1,K  ) + 2
             CELL(IC)%EDGE_INDEX(5)  = CELL_INDEX(I  ,J  ,K  ) + 3
             CELL(IC)%EDGE_INDEX(6)  = CELL_INDEX(I  ,J  ,K-1) + 3
             CELL(IC)%EDGE_INDEX(7)  = CELL_INDEX(I  ,J  ,K  ) + 4
             CELL(IC)%EDGE_INDEX(8)  = CELL_INDEX(I  ,J+1,K  ) + 4
             CELL(IC)%EDGE_INDEX(9)  = CELL_INDEX(I  ,J  ,K  ) + 5
             CELL(IC)%EDGE_INDEX(10) = CELL_INDEX(I  ,J  ,K-1) + 5 
             CELL(IC)%EDGE_INDEX(11) = CELL_INDEX(I  ,J  ,K  ) + 6
             CELL(IC)%EDGE_INDEX(12) = CELL_INDEX(I+1,J  ,K  ) + 6
             DO IE=1,NEDGE
                MAX_EDGE = MAX(MAX_EDGE,CELL(IC)%EDGE_INDEX(IE))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! For now let's have an MPI_Allreduce here so that each rank has the same size EDGE
    ! This is because we're essentially keeping a copy of the same data structure on each rank
    CALL MPI_ALLREDUCE(MAX_EDGE,GLOB_MAX_EDGE,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)
    ALLOCATE(EDGE(0:GLOB_MAX_EDGE))

    DO IE=1,MAX_EDGE,2
        EDGE(IE)%OMEGA = 1.5E-4_EB
        EDGE(IE)%TAU   = 2.5E-4_EB
    ENDDO
    
    RHO   = 1.19_EB
    RHO_0 = 1.19_EB
    D     = 0.0015_EB
    MU    = 0.0019_EB
    WORK1 = 0.0_EB; WORK2 = 0.0_EB; WORK3 = 0.0_EB; WORK4 = 0.0_EB; WORK5 = 0.0_EB; WORK6 = 0.0_EB
    GX(:) = 0.0_EB; GY(:) = 0.0_EB; GZ(:) = 1.0_EB

    ! U, V, W: determine velocity components over owned points
    DO K=LOC_MIN_K,LOC_MAX_K
       DO J=LOC_MIN_J,LOC_MAX_J
          DO I=LOC_MIN_I,LOC_MAX_I
             ! Some Trig functions:
             U(I,J,K) =  SIN(X(I))*COS(Y(J))*COS(Z(K))
             V(I,J,K) = -COS(X(I))*SIN(Y(J))*COS(Z(K))
             W(I,J,K) =  COS(X(I))*COS(Y(J))*SIN(Z(K))
          ENDDO
       ENDDO
    ENDDO

    UU => U
    VV => V
    WW => W
    DP => D
    RHOP => RHO
    TXY => WORK1
    TXZ => WORK2
    TYZ => WORK3
    OMX => WORK4
    OMY => WORK5
    OMZ => WORK6

    ! 3. Exchange velocities
    ALLOCATE(U_SEND(TOTAL_SEND_PTS))
    ALLOCATE(V_SEND(TOTAL_SEND_PTS))
    ALLOCATE(W_SEND(TOTAL_SEND_PTS))
    ALLOCATE(U_RECV(TOTAL_RECV_PTS))
    ALLOCATE(V_RECV(TOTAL_RECV_PTS))
    ALLOCATE(W_RECV(TOTAL_RECV_PTS))

    CALL EXCHANGE_VELOCITIES(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,&
                             HALO_RECV_DISPLS,U_SEND,V_SEND,W_SEND,U_RECV,V_RECV,W_RECV,&
                             RECV_COORDS,RECV_COORD_DISPLS,UU,VV,WW,&
                             MY_RANK,LOC_MIN_HALO_I,LOC_MIN_HALO_J,LOC_MIN_HALO_K,LOC_MAX_HALO_I,LOC_MAX_HALO_J,LOC_MAX_HALO_K)

    ! Compute Tau OMG:

    ! Compute vorticity & stress terms
    DO K=LOC_MIN_K,LOC_MAX_K
        DO J=LOC_MIN_J,LOC_MAX_J
           DO I=LOC_MIN_I,LOC_MAX_I
              DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
              DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
              DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
              DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
              DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
              DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
              OMX(I,J,K) = DWDY - DVDZ
              OMY(I,J,K) = DUDZ - DWDX
              OMZ(I,J,K) = DVDX - DUDY
              MUX = 0.25_EB*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
              MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
              MUZ = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
              TXY(I,J,K) = MUZ*(DVDX + DUDY)
              TXZ(I,J,K) = MUY*(DUDZ + DWDX)
              TYZ(I,J,K) = MUX*(DVDZ + DWDY)
           ENDDO
        ENDDO
     ENDDO

    ! 4. Exchange w and t
    ALLOCATE(OMX_SEND(TOTAL_SEND_PTS)) ! TOTAL_SEND_PTS ACROSS ALL NBRS. SEND OFFSETS IN HALO_SEND_SIZES
    ALLOCATE(OMY_SEND(TOTAL_SEND_PTS))
    ALLOCATE(OMZ_SEND(TOTAL_SEND_PTS))
    ALLOCATE(TXY_SEND(TOTAL_SEND_PTS))
    ALLOCATE(TXZ_SEND(TOTAL_SEND_PTS))
    ALLOCATE(TYZ_SEND(TOTAL_SEND_PTS))

    ALLOCATE(OMX_RECV(TOTAL_RECV_PTS))
    ALLOCATE(OMY_RECV(TOTAL_RECV_PTS))
    ALLOCATE(OMZ_RECV(TOTAL_RECV_PTS))
    ALLOCATE(TXY_RECV(TOTAL_RECV_PTS))
    ALLOCATE(TXZ_RECV(TOTAL_RECV_PTS))
    ALLOCATE(TYZ_RECV(TOTAL_RECV_PTS))

    CALL EXCHANGE_OMEGA_TAU(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,&
                            HALO_RECV_DISPLS,RECV_COORDS,RECV_COORD_DISPLS,OMX_SEND,OMX_RECV,OMY_SEND,OMY_RECV,OMZ_SEND,&
                            OMZ_RECV,TXZ_SEND,TXZ_RECV,TXY_SEND,TXY_RECV,TYZ_SEND,TYZ_RECV,&
                            OMX,OMY,OMZ,TXZ,TXY,TYZ)

    

    !$OMP TARGET ENTER DATA MAP(TO:UU(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1)) &
    !$OMP MAP(TO:VV(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                   &
    !$OMP MAP(TO:WW(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                   &            
    !$OMP MAP(TO:OMX(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(TO:OMY(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(TO:OMZ(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(TO:TXY(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(TO:TXZ(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(TO:TYZ(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(TO:CC(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                   &
    !$OMP MAP(TO:CELL(0:GLOB_MAX_IC)) MAP(TO:EDGE(0:GLOB_MAX_EDGE))                                                         &
    !$OMP MAP(TO:RHOP(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                 &
    !$OMP MAP(TO:RDX(LOC_MIN_HALO_I:LOC_MAX_HALO_I)) MAP(TO:RDY(LOC_MIN_HALO_J:LOC_MAX_HALO_J))                             &
    !$OMP MAP(TO:RDZ(LOC_MIN_HALO_K:LOC_MAX_HALO_K)) MAP(TO:RDXN(LOC_MIN_HALO_I:LOC_MAX_HALO_I))                            &
    !$OMP MAP(TO:RDYN(LOC_MIN_HALO_J:LOC_MAX_HALO_J)) MAP(TO:RDZN(LOC_MIN_HALO_K:LOC_MAX_HALO_K))                           &
    !$OMP MAP(TO:RHO_0(LOC_MIN_HALO_K:LOC_MAX_HALO_K)) MAP(TO:GX(LOC_MIN_HALO_I:LOC_MAX_HALO_I))                            &
    !$OMP MAP(TO:GY(LOC_MIN_HALO_I:LOC_MAX_HALO_I)) MAP(TO:GZ(LOC_MIN_HALO_I:LOC_MAX_HALO_I))                               &
    !$OMP MAP(TO:MU(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                   &
    !$OMP MAP(TO:DP(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))
                            
    T_SETUP_END = MPI_WTIME()
    SETUP_TIME = T_SETUP_END - T_SETUP_START
    CALL MPI_REDUCE(SETUP_TIME, MAX_SETUP_TIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, IERR)

    IF (MY_RANK==0) WRITE(*,*) 'BEGINNING SIM LOOP...'
    T_SIM_START = MPI_WTIME()
    SIM_LOOP: DO ISTEP = 1, NUM_TIME_STEPS
       ! WRITE(*,*) 'TIME STEP=',ISTEP
       CALL LOOP3D_OMP_GPU()
       ! CALL LOOP3D_SERIAL_CPU()
    END DO SIM_LOOP
    T_SIM_END = MPI_WTIME()
    CALL MPI_Barrier(MPI_COMM_WORLD, IERR)
    IF (MY_RANK==0) WRITE(*,*) 'FINISHED SIM LOOP...'
    !$OMP TARGET EXIT DATA MAP(FROM:FVX(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1)) &
    !$OMP MAP(FROM:FVY(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))                  &
    !$OMP MAP(FROM:FVZ(LOC_MIN_HALO_I:LOC_MAX_I_P1,LOC_MIN_HALO_J:LOC_MAX_J_P1,LOC_MIN_HALO_K:LOC_MAX_K_P1))

    ! Global reduction to compute the mean flux terms
    FLUX_AVGS(1) = SUM(FVX(LOC_MIN_I_P1:LOC_MAX_I,LOC_MIN_J_P1:LOC_MAX_J,LOC_MIN_K_P1:LOC_MAX_K))
    FLUX_AVGS(2) = SUM(FVY(LOC_MIN_I_P1:LOC_MAX_I,LOC_MIN_J_P1:LOC_MAX_J,LOC_MIN_K_P1:LOC_MAX_K))
    FLUX_AVGS(3) = SUM(FVZ(LOC_MIN_I_P1:LOC_MAX_I,LOC_MIN_J_P1:LOC_MAX_J,LOC_MIN_K_P1:LOC_MAX_K))

    CALL MPI_REDUCE(FLUX_AVGS, FLUX_AVGS, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERR)

    LOOP_TIME = T_SIM_END - T_SIM_START
    CALL MPI_REDUCE(LOOP_TIME, MIN_LOOP_TIME, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_REDUCE(LOOP_TIME, MAX_LOOP_TIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_REDUCE(LOOP_TIME, AVG_LOOP_TIME, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERR)

    IF (MY_RANK == 0) THEN
        FLUX_AVGS(1) = FLUX_AVGS(1) / (GLOB_MAX_I*GLOB_MAX_J*GLOB_MAX_K)
        FLUX_AVGS(2) = FLUX_AVGS(2) / (GLOB_MAX_I*GLOB_MAX_J*GLOB_MAX_K)
        FLUX_AVGS(3) = FLUX_AVGS(3) / (GLOB_MAX_I*GLOB_MAX_J*GLOB_MAX_K)
        AVG_LOOP_TIME = AVG_LOOP_TIME / NUM_TASKS
        WRITE(10,*) 'Setup Time=',MAX_SETUP_TIME
        WRITE(10,*) 'Loop Time(avg,min,max)=',AVG_LOOP_TIME,MIN_LOOP_TIME,MAX_LOOP_TIME
        WRITE(10,*) 'Total Time=',MAX_LOOP_TIME+MAX_SETUP_TIME
        WRITE(10,*) 'mean FVX =',FLUX_AVGS(1)
        WRITE(10,*) 'mean FVY =',FLUX_AVGS(2)
        WRITE(10,*) 'mean FVZ =',FLUX_AVGS(3)
        WRITE(10,*) 'Ending Loop3D'
        CLOSE(10)
        WRITE(*,*) 'Loop3D done.'
    END IF

    CALL MPI_FINALIZE(IERR)

    CONTAINS
    
    SUBROUTINE LOOP3D_OMP_GPU()
       
       ! Compute x-direction flux term FVX

!$OMP PARALLEL PRIVATE(WP,WM,VP,VM,UP,UM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXZP,TXZM,TXYP,TXYM,TYZP,TYZM, &
!$OMP& IC,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,RRHO,DUDX,DVDY,DWDZ,VTRM)

!$OMP DO SCHEDULE(STATIC) PRIVATE(WOMY, VOMZ, TXXP, TXXM, DTXXDX, DTXYDY, DTXZDZ)
        DO K=LOC_MIN_K_P1,LOC_MAX_K
            DO J=LOC_MIN_J_P1,LOC_MAX_J
               DO I=LOC_MIN_I,LOC_MAX_I
                  WP    = WW(I,J,K)   + WW(I+1,J,K)
                  WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
                  VP    = VV(I,J,K)   + VV(I+1,J,K)
                  VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
                  OMYP  = OMY(I,J,K)
                  OMYM  = OMY(I,J,K-1)
                  OMZP  = OMZ(I,J,K)
                  OMZM  = OMZ(I,J-1,K)
                  TXZP  = TXZ(I,J,K)
                  TXZM  = TXZ(I,J,K-1)
                  TXYP  = TXY(I,J,K)
                  TXYM  = TXY(I,J-1,K)
                  IC    = CC(I,J,K)
                  IEYP  = CELL(IC)%EDGE_INDEX(8)
                  IEYM  = CELL(IC)%EDGE_INDEX(6)
                  IEZP  = CELL(IC)%EDGE_INDEX(12)
                  IEZM  = CELL(IC)%EDGE_INDEX(10)
                  IF (EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
                     OMYP = EDGE(IEYP)%OMEGA(-1)
                     TXZP = EDGE(IEYP)%TAU(-1)
                  ENDIF
                  IF (EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
                     OMYM = EDGE(IEYM)%OMEGA( 1)
                     TXZM = EDGE(IEYM)%TAU( 1)
                  ENDIF
                  IF (EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
                     OMZP = EDGE(IEZP)%OMEGA(-2)
                     TXYP = EDGE(IEZP)%TAU(-2)
                  ENDIF
                  IF (EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
                     OMZM = EDGE(IEZM)%OMEGA( 2)
                     TXYM = EDGE(IEZM)%TAU( 2)
                  ENDIF
                  WOMY  = WP*OMYP + WM*OMYM
                  VOMZ  = VP*OMZP + VM*OMZM
                  RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
                  DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
                  DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
                  TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
                  DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
                  DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
                  TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
                  DTXXDX= RDXN(I)*(TXXP-TXXM)
                  DTXYDY= RDY(J) *(TXYP-TXYM)
                  DTXZDZ= RDZ(K) *(TXZP-TXZM)
                  VTRM  = DTXXDX + DTXYDY + DTXZDZ
                  FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO NOWAIT
         
         ! Compute y-direction flux term FVY
         
         !$OMP DO SCHEDULE(STATIC) PRIVATE(WOMX, UOMZ, TYYP, TYYM, DTXYDX, DTYYDY, DTYZDZ)
         DO K=LOC_MIN_K_P1,LOC_MAX_K
            DO J=LOC_MIN_J,LOC_MAX_J
               DO I=LOC_MIN_I_P1,LOC_MAX_I
                  UP    = UU(I,J,K)   + UU(I,J+1,K)
                  UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
                  WP    = WW(I,J,K)   + WW(I,J+1,K)
                  WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
                  OMXP  = OMX(I,J,K)
                  OMXM  = OMX(I,J,K-1)
                  OMZP  = OMZ(I,J,K)
                  OMZM  = OMZ(I-1,J,K)
                  TYZP  = TYZ(I,J,K)
                  TYZM  = TYZ(I,J,K-1)
                  TXYP  = TXY(I,J,K)
                  TXYM  = TXY(I-1,J,K)
                  IC    = CC(I,J,K)
                  IEXP  = CELL(IC)%EDGE_INDEX(4)
                  IEXM  = CELL(IC)%EDGE_INDEX(2)
                  IEZP  = CELL(IC)%EDGE_INDEX(12)
                  IEZM  = CELL(IC)%EDGE_INDEX(11)
                  IF (EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
                     OMXP = EDGE(IEXP)%OMEGA(-2)
                     TYZP = EDGE(IEXP)%TAU(-2)
                  ENDIF
                  IF (EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
                     OMXM = EDGE(IEXM)%OMEGA( 2)
                     TYZM = EDGE(IEXM)%TAU( 2)
                  ENDIF
                  IF (EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
                     OMZP = EDGE(IEZP)%OMEGA(-1)
                     TXYP = EDGE(IEZP)%TAU(-1)
                  ENDIF
                  IF (EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
                     OMZM = EDGE(IEZM)%OMEGA( 1)
                     TXYM = EDGE(IEZM)%TAU( 1)
                  ENDIF
                  WOMX  = WP*OMXP + WM*OMXM
                  UOMZ  = UP*OMZP + UM*OMZM
                  RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
                  DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
                  DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
                  TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
                  DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
                  DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
                  TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
                  DTXYDX= RDX(I) *(TXYP-TXYM)
                  DTYYDY= RDYN(J)*(TYYP-TYYM)
                  DTYZDZ= RDZ(K) *(TYZP-TYZM)
                  VTRM  = DTXYDX + DTYYDY + DTYZDZ
                  FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*RHO_0(K) - VTRM)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO NOWAIT
         
         ! Compute z-direction flux term FVZ
         
         !$OMP DO SCHEDULE(STATIC) PRIVATE(UOMY, VOMX, TZZP, TZZM, DTXZDX, DTYZDY, DTZZDZ)
         DO K=LOC_MIN_K,LOC_MAX_K
            DO J=LOC_MIN_J_P1,LOC_MAX_J
               DO I=LOC_MIN_I_P1,LOC_MAX_I
                  UP    = UU(I,J,K)   + UU(I,J,K+1)
                  UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
                  VP    = VV(I,J,K)   + VV(I,J,K+1)
                  VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
                  OMYP  = OMY(I,J,K)
                  OMYM  = OMY(I-1,J,K)
                  OMXP  = OMX(I,J,K)
                  OMXM  = OMX(I,J-1,K)
                  TXZP  = TXZ(I,J,K)
                  TXZM  = TXZ(I-1,J,K)
                  TYZP  = TYZ(I,J,K)
                  TYZM  = TYZ(I,J-1,K)
                  IC    = CC(I,J,K)
                  IEXP  = CELL(IC)%EDGE_INDEX(4)
                  IEXM  = CELL(IC)%EDGE_INDEX(3)
                  IEYP  = CELL(IC)%EDGE_INDEX(8)
                  IEYM  = CELL(IC)%EDGE_INDEX(7)
                  IF (EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
                     OMXP = EDGE(IEXP)%OMEGA(-1)
                     TYZP = EDGE(IEXP)%TAU(-1)
                  ENDIF
                  IF (EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
                     OMXM = EDGE(IEXM)%OMEGA( 1)
                     TYZM = EDGE(IEXM)%TAU( 1)
                  ENDIF
                  IF (EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
                     OMYP = EDGE(IEYP)%OMEGA(-2)
                     TXZP = EDGE(IEYP)%TAU(-2)
                  ENDIF
                  IF (EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
                     OMYM = EDGE(IEYM)%OMEGA( 2)
                     TXZM = EDGE(IEYM)%TAU( 2)
                  ENDIF
                  UOMY  = UP*OMYP + UM*OMYM
                  VOMX  = VP*OMXP + VM*OMXM
                  RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
                  DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
                  DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
                  TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
                  DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
                  DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
                  TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
                  DTXZDX= RDX(I) *(TXZP-TXZM)
                  DTYZDY= RDY(J) *(TYZP-TYZM)
                  DTZZDZ= RDZN(K)*(TZZP-TZZM)
                  VTRM  = DTXZDX + DTYZDY + DTZZDZ
                  FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(RHO_0(K)+RHO_0(K+1)) - VTRM)
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO NOWAIT
         
         !$OMP END PARALLEL
    END SUBROUTINE LOOP3D_OMP_GPU

    SUBROUTINE EXCHANGE_OMEGA_TAU(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,&
                            HALO_RECV_DISPLS,RECV_COORDS,RECV_COORD_DISPLS,OMX_SEND,OMX_RECV,OMY_SEND,OMY_RECV,OMZ_SEND,&
                            OMZ_RECV,TXZ_SEND,TXZ_RECV,TXY_SEND,TXY_RECV,TYZ_SEND,TYZ_RECV,&
                            OMX,OMY,OMZ,TXZ,TXY,TYZ)
        
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: NBR_EXT_BBOXES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_DISPLS (:)
        INTEGER, INTENT(IN) :: HALO_RECV_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_RECV_DISPLS (:)
        INTEGER, INTENT(IN) :: RECV_COORDS (:)
        INTEGER, INTENT(IN) :: RECV_COORD_DISPLS (:)
        REAL(EB), INTENT(IN OUT) :: OMX_SEND (:)
        REAL(EB), INTENT(IN OUT) :: OMY_SEND (:)
        REAL(EB), INTENT(IN OUT) :: OMZ_SEND (:)
        REAL(EB), INTENT(IN OUT) :: TXZ_SEND (:)
        REAL(EB), INTENT(IN OUT) :: TXY_SEND (:)
        REAL(EB), INTENT(IN OUT) :: TYZ_SEND (:)
        REAL(EB), INTENT(IN OUT) :: OMX_RECV (:)
        REAL(EB), INTENT(IN OUT) :: OMY_RECV (:)
        REAL(EB), INTENT(IN OUT) :: OMZ_RECV (:)
        REAL(EB), INTENT(IN OUT) :: TXZ_RECV (:)
        REAL(EB), INTENT(IN OUT) :: TXY_RECV (:)
        REAL(EB), INTENT(IN OUT) :: TYZ_RECV (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: OMX (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: OMY (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: OMZ (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: TXZ (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: TXY (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: TYZ (:,:,:)
        INTEGER :: RK, POINT, COORD_OFFSET, RANK_OFFSET, RECV_PTS, I, J, K

        DO RK=1, NUM_NBRS
            RANK_OFFSET = HALO_SEND_DISPLS(RK)
            POINT = 0
            DO K=LOC_MIN_K,LOC_MAX_K
                DO J=LOC_MIN_J,LOC_MAX_J
                    DO I=LOC_MIN_I,LOC_MAX_I
                        IF (IS_PT_INSIDE(I, J, K, NBR_EXT_BBOXES, RK-1, BOX_SIZE)) THEN
                            OMX_SEND(RANK_OFFSET + POINT) = OMX(I,J,K)
                            OMY_SEND(RANK_OFFSET + POINT) = OMY(I,J,K)
                            OMZ_SEND(RANK_OFFSET + POINT) = OMZ(I,J,K)
                            TXZ_SEND(RANK_OFFSET + POINT) = TXZ(I,J,K)
                            TXY_SEND(RANK_OFFSET + POINT) = TXY(I,J,K)
                            TYZ_SEND(RANK_OFFSET + POINT) = TYZ(I,J,K)
                            POINT = POINT + 1
                        END IF
                    END DO
                END DO
            END DO
        END DO

        ! Exchange each of the velocity components along the halo regions
        CALL EXCHANGE_DOUBLE_DATA(OMX_SEND,OMX_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(OMY_SEND,OMY_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(OMZ_SEND,OMZ_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(TXY_SEND,TXY_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(TXZ_SEND,TXZ_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(TYZ_SEND,TYZ_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)

        ! Set up the communicated velocities
        DO RK=1, NUM_NBRS
            RANK_OFFSET = HALO_RECV_DISPLS(RK)
            COORD_OFFSET = RECV_COORD_DISPLS(RK)
            RECV_PTS = HALO_RECV_SIZES(RK)
            DO POINT=0,RECV_PTS-1
                I = RECV_COORDS(COORD_OFFSET + POINT*3)
                J = RECV_COORDS(COORD_OFFSET + POINT*3 + 1)
                K = RECV_COORDS(COORD_OFFSET + POINT*3 + 2)
                OMX(I,J,K) = OMX_RECV(RANK_OFFSET + POINT)
                OMY(I,J,K) = OMY_RECV(RANK_OFFSET + POINT)
                OMZ(I,J,K) = OMZ_RECV(RANK_OFFSET + POINT)
                TXZ(I,J,K) = TXZ_RECV(RANK_OFFSET + POINT)
                TXY(I,J,K) = TXY_RECV(RANK_OFFSET + POINT)
                TYZ(I,J,K) = TYZ_RECV(RANK_OFFSET + POINT)
            END DO
        END DO

    END SUBROUTINE EXCHANGE_OMEGA_TAU

    SUBROUTINE EXCHANGE_CELL_INDEX(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,&
                                  HALO_RECV_DISPLS,RECV_COORDS,RECV_COORD_DISPLS,CELL_INDEX_SEND,CELL_INDEX_RECV,CC,MY_RANK)
    
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: NBR_EXT_BBOXES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_DISPLS (:)
        INTEGER, INTENT(IN) :: HALO_RECV_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_RECV_DISPLS (:)
        INTEGER, INTENT(IN) :: RECV_COORDS (:)
        INTEGER, INTENT(IN) :: RECV_COORD_DISPLS (:)
        INTEGER, INTENT(IN) :: MY_RANK
        INTEGER, POINTER, INTENT(IN OUT) :: CC (:,:,:) ! pointer to cell index
        INTEGER, INTENT(IN OUT) :: CELL_INDEX_SEND (:)
        INTEGER, INTENT(IN OUT) :: CELL_INDEX_RECV (:)
        INTEGER :: RK, POINT, COORD_OFFSET, RANK_OFFSET, RECV_PTS, I, J, K, DUMMY

        DO RK=1, NUM_NBRS
            RANK_OFFSET = HALO_SEND_DISPLS(RK)
            POINT = 0
            DO K=LOC_MIN_K,LOC_MAX_K
                DO J=LOC_MIN_J,LOC_MAX_J
                    DO I=LOC_MIN_I,LOC_MAX_I
                        IF (IS_PT_INSIDE(I, J, K, NBR_EXT_BBOXES, RK-1, BOX_SIZE)) THEN
                            CELL_INDEX_SEND(RANK_OFFSET + POINT) = CC(I,J,K)
                            POINT = POINT + 1
                        END IF
                    END DO
                END DO
            END DO
        END DO
        
        CALL EXCHANGE_INTEGER_DATA(CELL_INDEX_SEND,CELL_INDEX_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)


        ! Setup received cell index 
        DO RK=1, NUM_NBRS
            RANK_OFFSET = HALO_RECV_DISPLS(RK)
            COORD_OFFSET = RECV_COORD_DISPLS(RK)
            RECV_PTS = HALO_RECV_SIZES(RK)
            DO POINT=0,RECV_PTS-1
                I = RECV_COORDS(COORD_OFFSET + POINT*3)
                J = RECV_COORDS(COORD_OFFSET + POINT*3 + 1)
                K = RECV_COORDS(COORD_OFFSET + POINT*3 + 2)
                CC(I,J,K) = CELL_INDEX_RECV(RANK_OFFSET + POINT)
            END DO
        END DO

    END SUBROUTINE EXCHANGE_CELL_INDEX

    SUBROUTINE EXCHANGE_VELOCITIES(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,&
                                   HALO_RECV_DISPLS,U_SEND,V_SEND,W_SEND,U_RECV,V_RECV,W_RECV,&
                                   RECV_COORDS,RECV_COORD_DISPLS,UU,VV,WW,&
                                   MY_RANK,LOC_MIN_HALO_I,LOC_MIN_HALO_J,LOC_MIN_HALO_K,LOC_MAX_HALO_I,LOC_MAX_HALO_J,LOC_MAX_HALO_K)

        INTEGER, INTENT(IN) :: MY_RANK,LOC_MIN_HALO_I,LOC_MIN_HALO_J,LOC_MIN_HALO_K,LOC_MAX_HALO_I,LOC_MAX_HALO_J,LOC_MAX_HALO_K
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: NBR_EXT_BBOXES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_DISPLS (:)
        INTEGER, INTENT(IN) :: HALO_RECV_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_RECV_DISPLS (:)
        INTEGER, INTENT(IN) :: RECV_COORDS (:)
        INTEGER, INTENT(IN) :: RECV_COORD_DISPLS (:)
        REAL(EB), INTENT(IN OUT) :: U_SEND (:)
        REAL(EB), INTENT(IN OUT) :: V_SEND (:)
        REAL(EB), INTENT(IN OUT) :: W_SEND (:)
        REAL(EB), INTENT(IN OUT) :: U_RECV (:)
        REAL(EB), INTENT(IN OUT) :: V_RECV (:)
        REAL(EB), INTENT(IN OUT) :: W_RECV (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: UU (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: VV (:,:,:)
        REAL(EB), POINTER, INTENT(IN OUT) :: WW (:,:,:)
        INTEGER :: RK, POINT, COORD_OFFSET, RANK_OFFSET, RECV_PTS, I, J, K

        DO RK=1, NUM_NBRS
            RANK_OFFSET = HALO_SEND_DISPLS(RK)
            POINT = 0
            DO K=LOC_MIN_K,LOC_MAX_K
                DO J=LOC_MIN_J,LOC_MAX_J
                    DO I=LOC_MIN_I,LOC_MAX_I
                        IF (IS_PT_INSIDE(I, J, K, NBR_EXT_BBOXES, RK-1, BOX_SIZE)) THEN
                            U_SEND(RANK_OFFSET + POINT) = UU(I,J,K)
                            V_SEND(RANK_OFFSET + POINT) = VV(I,J,K)
                            W_SEND(RANK_OFFSET + POINT) = WW(I,J,K)
                            POINT = POINT + 1
                        END IF
                    END DO
                END DO
            END DO
        END DO
        ! Exchange each of the velocity components along the halo regions
        CALL EXCHANGE_DOUBLE_DATA(U_SEND,U_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(V_SEND,V_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        CALL EXCHANGE_DOUBLE_DATA(W_SEND,W_RECV,NUM_NBRS,HALO_SEND_SIZES,HALO_SEND_DISPLS,HALO_RECV_SIZES,HALO_RECV_DISPLS)
        ! Set up the communicated velocities
        DO RK=1, NUM_NBRS
            RANK_OFFSET = HALO_RECV_DISPLS(RK)
            COORD_OFFSET = RECV_COORD_DISPLS(RK)
            RECV_PTS = HALO_RECV_SIZES(RK)
            DO POINT=0,RECV_PTS-1
                I = RECV_COORDS(COORD_OFFSET + POINT*3)
                J = RECV_COORDS(COORD_OFFSET + POINT*3 + 1)
                K = RECV_COORDS(COORD_OFFSET + POINT*3 + 2)
                ! Error checker. This statement should never print.
                IF ((I < LOC_MIN_HALO_I .OR. I > LOC_MAX_HALO_I) .OR. (J < LOC_MIN_HALO_J .OR. J > LOC_MAX_HALO_J) .OR. (K < LOC_MIN_HALO_K .OR. K > LOC_MAX_HALO_K)) THEN 
                    WRITE(*,*) 'RANK ',MY_RANK,' RECEIVED AN INVALID COORDINATE (I,J,K): ',I,',',J,',',K,' FROM NEIGHBOR RANK=',NBR_RANKS(RK)
                END IF
                UU(I,J,K) = U_RECV(RANK_OFFSET + POINT)
                VV(I,J,K) = V_RECV(RANK_OFFSET + POINT)
                WW(I,J,K) = W_RECV(RANK_OFFSET + POINT)
            END DO
        END DO

    END SUBROUTINE EXCHANGE_VELOCITIES

    SUBROUTINE COMPUTE_COORDINATE_DISPLS(NUM_NBRS,HALO_SEND_SIZES,HALO_RECV_SIZES,SEND_COORD_DISPLS,RECV_COORD_DISPLS)
        
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: HALO_SEND_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_RECV_SIZES (:)
        INTEGER, INTENT(OUT) :: SEND_COORD_DISPLS (:)
        INTEGER, INTENT(OUT) :: RECV_COORD_DISPLS (:)
        INTEGER :: RK, TOTAL_SUM, XYZ_DIM

        TOTAL_SUM = 1
        XYZ_DIM = 3
        DO RK=1,NUM_NBRS
            SEND_COORD_DISPLS(RK) = TOTAL_SUM
            TOTAL_SUM = TOTAL_SUM + XYZ_DIM * HALO_SEND_SIZES(RK)
        END DO

        TOTAL_SUM = 1
        DO RK=1,NUM_NBRS
            RECV_COORD_DISPLS(RK) = TOTAL_SUM
            TOTAL_SUM = TOTAL_SUM + XYZ_DIM * HALO_RECV_SIZES(RK)
        END DO

    END SUBROUTINE COMPUTE_COORDINATE_DISPLS
    
    SUBROUTINE EXCHANGE_COORDINATES(NUM_NBRS,NBR_EXT_BBOXES,HALO_SEND_SIZES,HALO_RECV_SIZES,SEND_COORDS,RECV_COORDS,&
                                    XP,YP,ZP,RDXNP,RDXP,RDYNP,RDYP,RDZNP,RDZP,MY_RANK,NBR_RANKS,&
                                    SEND_COORD_DISPLS,RECV_COORD_DISPLS)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: NBR_RANKS (:)
        INTEGER, INTENT(IN) :: NBR_EXT_BBOXES (:)
        INTEGER, INTENT(IN) :: HALO_SEND_SIZES (:)
        INTEGER, INTENT(IN) :: HALO_RECV_SIZES (:)
        INTEGER, INTENT(IN OUT) :: SEND_COORDS (:)
        INTEGER, INTENT(IN OUT) :: RECV_COORDS (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: XP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: YP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: ZP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: RDXNP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: RDXP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: RDYNP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: RDYP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: RDZNP (:)
        REAL(EB), POINTER, INTENT(IN OUT) :: RDZP (:)
        INTEGER, INTENT(IN) :: MY_RANK
        INTEGER, INTENT(IN) :: SEND_COORD_DISPLS (:)
        INTEGER, INTENT(IN) :: RECV_COORD_DISPLS (:)
        INTEGER :: RK, POINT, COORD_OFFSET, RECV_PTS, I, J, K, XYZ_DIM

        XYZ_DIM = 3

        DO RK=1, NUM_NBRS
            COORD_OFFSET = SEND_COORD_DISPLS(RK)
            POINT = 0
            DO K=LOC_MIN_K,LOC_MAX_K
                DO J=LOC_MIN_J,LOC_MAX_J
                    DO I=LOC_MIN_I,LOC_MAX_I
                        IF (IS_PT_INSIDE(I, J, K, NBR_EXT_BBOXES, RK-1, BOX_SIZE)) THEN
                            SEND_COORDS(COORD_OFFSET + POINT*3 + 0) = I
                            SEND_COORDS(COORD_OFFSET + POINT*3 + 1) = J
                            SEND_COORDS(COORD_OFFSET + POINT*3 + 2) = K
                            POINT = POINT + 1
                        END IF
                    END DO
                END DO
            END DO
        END DO

        ! Exchange the point coordinates
        CALL EXCHANGE_INTEGER_DATA(SEND_COORDS,RECV_COORDS,NUM_NBRS,XYZ_DIM*HALO_SEND_SIZES,SEND_COORD_DISPLS,XYZ_DIM*HALO_RECV_SIZES,RECV_COORD_DISPLS)

        ! Setup received coordinate data
        DO RK=1, NUM_NBRS
            COORD_OFFSET = RECV_COORD_DISPLS(RK)
            RECV_PTS = HALO_RECV_SIZES(RK)
            DO POINT=0,RECV_PTS-1
                I = RECV_COORDS(COORD_OFFSET + POINT*3)
                J = RECV_COORDS(COORD_OFFSET + POINT*3 + 1)
                K = RECV_COORDS(COORD_OFFSET + POINT*3 + 2)
                XP(I) = REAL(I,EB)
                RDXNP(I) = 1.0_EB
                RDXP(I) = 1.0_EB
                YP(J) = REAL(J,EB)
                RDYNP(J) = 1.0_EB
                RDYP(J) = 1.0_EB
                ZP(K) = REAL(K,EB)
                RDZNP(K) = 1.0_EB
                RDZP(K) = 1.0_EB
            END DO
        END DO

    END SUBROUTINE EXCHANGE_COORDINATES

    SUBROUTINE EXCHANGE_DOUBLE_DATA(SEND_DATA, RECV_DATA, NUM_NBRS,SEND_SIZES,SEND_DISPLS,&
                                    RECV_SIZES,RECV_DISPLS)
        REAL(EB), INTENT(IN) :: SEND_DATA (:)
        REAL(EB), INTENT(IN OUT) :: RECV_DATA (:)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: SEND_SIZES (:)
        INTEGER, INTENT(IN) :: SEND_DISPLS (:)
        INTEGER, INTENT(IN) :: RECV_SIZES (:)
        INTEGER, INTENT(IN) :: RECV_DISPLS (:)
        INTEGER :: SEND_REQUEST(NUM_NBRS)
        INTEGER :: RECV_REQUEST(NUM_NBRS)
        INTEGER :: RK, NBR_RANK
        INTEGER :: TAG

        TAG = 323

        DO RK=1,NUM_NBRS
            NBR_RANK = NBR_RANKS(RK)
            CALL MPI_IRECV(RECV_DATA(RECV_DISPLS(RK)),RECV_SIZES(RK),MPI_DOUBLE_PRECISION,&
                                     NBR_RANK,TAG,MPI_COMM_WORLD,RECV_REQUEST(RK),IERR)
        END DO

        DO RK=1,NUM_NBRS
            NBR_RANK = NBR_RANKS(RK)
            CALL MPI_ISEND(SEND_DATA(SEND_DISPLS(RK)),SEND_SIZES(RK),MPI_DOUBLE_PRECISION,&
                                     NBR_RANK,TAG,MPI_COMM_WORLD,SEND_REQUEST(RK),IERR)
        END DO

        CALL MPI_WAITALL(NUM_NBRS,SEND_REQUEST,MPI_STATUSES_IGNORE,IERR)
        CALL MPI_WAITALL(NUM_NBRS,RECV_REQUEST,MPI_STATUSES_IGNORE,IERR)

    END SUBROUTINE EXCHANGE_DOUBLE_DATA

    SUBROUTINE EXCHANGE_INTEGER_DATA(SEND_DATA, RECV_DATA, NUM_NBRS, SEND_SIZES, SEND_DISPLS,&
                                     RECV_SIZES, RECV_DISPLS)
        INTEGER, INTENT(IN) :: SEND_DATA (:)
        INTEGER, INTENT(IN OUT) :: RECV_DATA (:)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: SEND_SIZES (:) ! need to multiply by 3 with send coords
        INTEGER, INTENT(IN) :: SEND_DISPLS (:) ! need to multiply by 3 with send coords
        INTEGER, INTENT(IN) :: RECV_SIZES (:)
        INTEGER, INTENT(IN) :: RECV_DISPLS (:)
        INTEGER :: SEND_REQUEST(NUM_NBRS)
        INTEGER :: RECV_REQUEST(NUM_NBRS)
        INTEGER :: RK, NBR_RANK
        INTEGER :: TAG

        TAG = 222

        DO RK=1,NUM_NBRS
            NBR_RANK = NBR_RANKS(RK)
            CALL MPI_IRECV(RECV_DATA(RECV_DISPLS(RK)),RECV_SIZES(RK),MPI_INTEGER,&
                                     NBR_RANK,TAG,MPI_COMM_WORLD,RECV_REQUEST(RK),IERR)
        END DO

        DO RK=1,NUM_NBRS
            NBR_RANK = NBR_RANKS(RK)
            CALL MPI_ISEND(SEND_DATA(SEND_DISPLS(RK)),SEND_SIZES(RK),MPI_INTEGER,&
                                     NBR_RANK,TAG,MPI_COMM_WORLD,SEND_REQUEST(RK),IERR)
        END DO

        CALL MPI_WAITALL(NUM_NBRS,SEND_REQUEST,MPI_STATUSES_IGNORE,IERR)
        CALL MPI_WAITALL(NUM_NBRS,RECV_REQUEST,MPI_STATUSES_IGNORE,IERR)


    END SUBROUTINE EXCHANGE_INTEGER_DATA

    SUBROUTINE EXCHANGE_SIZES(SEND_SIZES, RECV_SIZES, NUM_NBRS)
        INTEGER, INTENT(IN) :: SEND_SIZES (:)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN OUT) :: RECV_SIZES (:)
        INTEGER :: SEND_REQUEST(NUM_NBRS)
        INTEGER :: RECV_REQUEST(NUM_NBRS)
        INTEGER :: TAG
        INTEGER :: RK
        INTEGER :: NBR_RANK
        INTEGER :: IERR

        TAG = 123

        ! Post receives
        DO RK=1,NUM_NBRS
            NBR_RANK = NBR_RANKS(RK)
            CALL MPI_IRECV(RECV_SIZES(RK), 1, MPI_INTEGER, NBR_RANK, TAG, MPI_COMM_WORLD, RECV_REQUEST(RK), IERR)
        END DO

        ! Post sends
        DO RK=1,NUM_NBRS
            NBR_RANK = NBR_RANKS(RK)
            CALL MPI_ISEND(SEND_SIZES(RK), 1, MPI_INTEGER, NBR_RANK, TAG, MPI_COMM_WORLD, SEND_REQUEST(RK), IERR)
        END DO

        CALL MPI_WAITALL(NUM_NBRS,SEND_REQUEST,MPI_STATUSES_IGNORE,IERR)
        CALL MPI_WAITALL(NUM_NBRS,RECV_REQUEST,MPI_STATUSES_IGNORE,IERR)

    END SUBROUTINE EXCHANGE_SIZES

    SUBROUTINE COMPUTE_TOTAL_RECV_PTS(TOTAL_RECV_PTS, HALO_RECV_SIZES, NUM_NBRS)
        INTEGER, INTENT(IN OUT) :: TOTAL_RECV_PTS
        INTEGER, INTENT(IN) :: HALO_RECV_SIZES (:)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER :: RK
        
        TOTAL_RECV_PTS = 0

        DO RK=1,NUM_NBRS
            TOTAL_RECV_PTS = TOTAL_RECV_PTS + HALO_RECV_SIZES(RK)
        END DO

    END SUBROUTINE COMPUTE_TOTAL_RECV_PTS
    
    SUBROUTINE DETERMINE_HALO_SEND_SIZES(HALO_SEND_SIZES, TOTAL_SEND_PTS, NBR_EXT_BBOXES, NUM_NBRS, BOX_SIZE)
        INTEGER, INTENT(IN OUT) :: HALO_SEND_SIZES (:)
        INTEGER, INTENT(IN OUT) :: TOTAL_SEND_PTS
        INTEGER, INTENT(IN) :: NBR_EXT_BBOXES (:)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER :: RK, I, J, K
        
        HALO_SEND_SIZES = 0
        TOTAL_SEND_PTS = 0

        DO K=LOC_MIN_K,LOC_MAX_K
            DO J=LOC_MIN_J,LOC_MAX_J
                DO I=LOC_MIN_I,LOC_MAX_I
                    DO RK=1,NUM_NBRS
                        IF (IS_PT_INSIDE(I, J, K, NBR_EXT_BBOXES, RK-1, BOX_SIZE)) THEN
                            HALO_SEND_SIZES(RK) = HALO_SEND_SIZES(RK) + 1
                            TOTAL_SEND_PTS = TOTAL_SEND_PTS + 1
                        END IF
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE DETERMINE_HALO_SEND_SIZES

    SUBROUTINE DETERMINE_HALO_DISPLACEMENTS(HALO_SIZES, HALO_DISPLS, NUM_NBRS)
        INTEGER, INTENT(IN) :: HALO_SIZES (:)
        INTEGER, INTENT(IN OUT) :: HALO_DISPLS (:)
        INTEGER, INTENT(IN) :: NUM_NBRS
        INTEGER :: TOTAL_PTS  ! fortran being indexed beginning at 1
        INTEGER :: RK

        TOTAL_PTS = 1

        DO RK=1,NUM_NBRS
            HALO_DISPLS(RK) = TOTAL_PTS
            TOTAL_PTS = TOTAL_PTS + HALO_SIZES(RK)
        END DO

    END SUBROUTINE DETERMINE_HALO_DISPLACEMENTS

    SUBROUTINE GET_NBR_BBOXES(MY_RANK, MY_BBOX, ALL_TASK_EXT_BBOXES, NBR_EXT_BBOXES, NBR_RANKS, NUM_TASKS, BOX_SIZE)
        INTEGER, INTENT(IN) :: MY_RANK
        INTEGER, INTENT(IN) :: MY_BBOX (:)
        INTEGER, INTENT(IN) :: ALL_TASK_EXT_BBOXES (:)
        INTEGER, INTENT(INOUT) :: NBR_EXT_BBOXES (:)
        INTEGER, INTENT(INOUT) :: NBR_RANKS (:)
        INTEGER, INTENT(IN) :: NUM_TASKS
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER :: NBR_EXT_BOX (BOX_SIZE)
        INTEGER :: RK, NBR_COUNT

        NBR_COUNT = 1
        DO RK=0,NUM_TASKS-1
            IF (RK == MY_RANK) CYCLE
            NBR_EXT_BOX(1) = ALL_TASK_EXT_BBOXES(RK * BOX_SIZE)
            NBR_EXT_BOX(2) = ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 1)
            NBR_EXT_BOX(3) = ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 2)
            NBR_EXT_BOX(4) = ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 3)
            NBR_EXT_BOX(5) = ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 4)
            NBR_EXT_BOX(6) = ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 5)
            IF (NBR_EXT_BOX(1) > MY_BBOX(4)) CYCLE
            IF (NBR_EXT_BOX(4) < MY_BBOX(1)) CYCLE
            IF (NBR_EXT_BOX(2) > MY_BBOX(5)) CYCLE
            IF (NBR_EXT_BOX(5) < MY_BBOX(2)) CYCLE
            IF (NBR_EXT_BOX(3) > MY_BBOX(6)) CYCLE
            IF (NBR_EXT_BOX(6) < MY_BBOX(3)) CYCLE
            NBR_EXT_BBOXES((NBR_COUNT-1) * BOX_SIZE) = NBR_EXT_BOX(1)
            NBR_EXT_BBOXES((NBR_COUNT-1) * BOX_SIZE + 1) = NBR_EXT_BOX(2)
            NBR_EXT_BBOXES((NBR_COUNT-1) * BOX_SIZE + 2) = NBR_EXT_BOX(3)
            NBR_EXT_BBOXES((NBR_COUNT-1) * BOX_SIZE + 3) = NBR_EXT_BOX(4)
            NBR_EXT_BBOXES((NBR_COUNT-1) * BOX_SIZE + 4) = NBR_EXT_BOX(5)
            NBR_EXT_BBOXES((NBR_COUNT-1)* BOX_SIZE + 5) = NBR_EXT_BOX(6)
            NBR_RANKS(NBR_COUNT) = RK
            NBR_COUNT = NBR_COUNT + 1
        END DO
    END SUBROUTINE GET_NBR_BBOXES

    SUBROUTINE DETERMINE_NUM_NEIGHBORS(MY_RANK, GLOBAL_MIN, GLOBAL_MAX, MY_BBOX, NUM_TASKS, NUM_NBRS, ALL_TASK_BBOXES, ALL_TASK_EXT_BBOXES, HALO_DEPTH, BOX_SIZE)
        INTEGER, INTENT(IN) :: MY_RANK
        INTEGER, INTENT(IN) :: GLOBAL_MIN (:)
        INTEGER, INTENT(IN) :: GLOBAL_MAX (:)
        INTEGER, INTENT(IN) :: MY_BBOX (:)
        INTEGER, INTENT(IN) :: NUM_TASKS
        INTEGER, INTENT(INOUT) :: NUM_NBRS
        INTEGER, INTENT(IN) :: ALL_TASK_BBOXES (:)
        INTEGER, INTENT(INOUT) :: ALL_TASK_EXT_BBOXES (:)
        INTEGER, INTENT(IN) :: HALO_DEPTH
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER :: NBR_EXT_BOX (BOX_SIZE)
        INTEGER :: RK

        NUM_NBRS = 0

        DO RK=0,NUM_TASKS-1
            IF (RK == MY_RANK) CYCLE
            NBR_EXT_BOX(1) = ALL_TASK_BBOXES(RK * BOX_SIZE)
            NBR_EXT_BOX(2) = ALL_TASK_BBOXES(RK * BOX_SIZE + 1)
            NBR_EXT_BOX(3) = ALL_TASK_BBOXES(RK * BOX_SIZE + 2)
            NBR_EXT_BOX(4) = ALL_TASK_BBOXES(RK * BOX_SIZE + 3)
            NBR_EXT_BOX(5) = ALL_TASK_BBOXES(RK * BOX_SIZE + 4)
            NBR_EXT_BOX(6) = ALL_TASK_BBOXES(RK * BOX_SIZE + 5)
            IF (NBR_EXT_BOX(1) - HALO_DEPTH >= GLOBAL_MIN(1)) THEN
                NBR_EXT_BOX(1) = NBR_EXT_BOX(1) - HALO_DEPTH
            END IF
            IF (NBR_EXT_BOX(2) - HALO_DEPTH >= GLOBAL_MIN(2)) THEN
                NBR_EXT_BOX(2) = NBR_EXT_BOX(2) - HALO_DEPTH
            END IF
            IF (NBR_EXT_BOX(3) - HALO_DEPTH >= GLOBAL_MIN(3)) THEN
                NBR_EXT_BOX(3) = NBR_EXT_BOX(3) - HALO_DEPTH
            END IF
            IF (NBR_EXT_BOX(4) + HALO_DEPTH <= GLOBAL_MAX(1)) THEN
                NBR_EXT_BOX(4) = NBR_EXT_BOX(4) + HALO_DEPTH
            END IF
            IF (NBR_EXT_BOX(5) + HALO_DEPTH <= GLOBAL_MAX(2)) THEN
                NBR_EXT_BOX(5) = NBR_EXT_BOX(5) + HALO_DEPTH
            END IF
            IF (NBR_EXT_BOX(6) + HALO_DEPTH <= GLOBAL_MAX(3)) THEN
                NBR_EXT_BOX(6) = NBR_EXT_BOX(6) + HALO_DEPTH
            END IF
            ALL_TASK_EXT_BBOXES(RK * BOX_SIZE) = NBR_EXT_BOX(1)
            ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 1) = NBR_EXT_BOX(2)
            ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 2) = NBR_EXT_BOX(3)
            ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 3) = NBR_EXT_BOX(4)
            ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 4) = NBR_EXT_BOX(5)
            ALL_TASK_EXT_BBOXES(RK * BOX_SIZE + 5) = NBR_EXT_BOX(6)
            IF (NBR_EXT_BOX(1) > MY_BBOX(4)) CYCLE
            IF (NBR_EXT_BOX(4) < MY_BBOX(1)) CYCLE
            IF (NBR_EXT_BOX(2) > MY_BBOX(5)) CYCLE
            IF (NBR_EXT_BOX(5) < MY_BBOX(2)) CYCLE
            IF (NBR_EXT_BOX(3) > MY_BBOX(6)) CYCLE
            IF (NBR_EXT_BOX(6) < MY_BBOX(3)) CYCLE
            NUM_NBRS = NUM_NBRS + 1
        END DO

    END SUBROUTINE DETERMINE_NUM_NEIGHBORS

    SUBROUTINE DECOMPOSE_DOMAIN(LOC_BOXES, MIN_I, MIN_J, MIN_K, MAX_I, MAX_J, MAX_K, BOX_SIZE, NUM_TASKS, MY_RANK)
        INTEGER, INTENT(INOUT) :: LOC_BOXES (:) ! AoS
        INTEGER, INTENT(IN) :: MIN_I, MIN_J, MIN_K, MAX_I, MAX_J, MAX_K
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER, INTENT(IN) :: NUM_TASKS, MY_RANK
        INTEGER :: MAX_BOX_IND, NEW_BOX_IND, NEW_BOX_1_IND, NEW_BOX_2_IND, DIM, DIM_VALUE, MAX_DIM, MAX_DIM_VALUE
        INTEGER :: RUNNING_BOX_TOTAL
        INTEGER :: ONE_THIRD, MID_PT
        INTEGER :: TOTAL_NUM_PROCS
        INTEGER :: PROC
        INTEGER :: MAX_ELE
        REAL(EB) :: LOC_VOLUMES(NUM_TASKS)

        TOTAL_NUM_PROCS = 1
        RUNNING_BOX_TOTAL = TOTAL_NUM_PROCS
        LOC_VOLUMES = 0
        
        ! initialize the first box to the global box coordinates...
        LOC_BOXES(0) = MIN_I
        LOC_BOXES(1) = MIN_J
        LOC_BOXES(2) = MIN_K
        LOC_BOXES(3) = MAX_I
        LOC_BOXES(4) = MAX_J
        LOC_BOXES(5) = MAX_K

        ! if NUM_TASKS is a power of 3, then first break that down by 3
        IF (MOD(NUM_TASKS, 3) == 0) THEN
            MAX_BOX_IND = 0
            MAX_DIM = -1
            MAX_DIM_VALUE = -1
            DO DIM=0,2
                DIM_VALUE = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) - LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM) + 1
                IF (DIM_VALUE > MAX_DIM_VALUE) THEN
                    MAX_DIM_VALUE = DIM_VALUE
                    MAX_DIM = DIM
                END IF
            END DO
            NEW_BOX_1_IND = RUNNING_BOX_TOTAL
            NEW_BOX_2_IND = RUNNING_BOX_TOTAL + 1
            DO DIM=0,2
                IF (DIM == MAX_DIM) THEN
                    ONE_THIRD = (LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) - LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM)) / 3
                    ! original box gets one third. [ --- maxbox --- ] --> [ maxbox/3 | newbox1 | newbox2]
                    LOC_BOXES(NEW_BOX_2_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM)
                    LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM) + ONE_THIRD
                    LOC_BOXES(NEW_BOX_1_IND * BOX_SIZE + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) + 1
                    LOC_BOXES(NEW_BOX_1_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(NEW_BOX_1_IND * BOX_SIZE + DIM) + ONE_THIRD
                    LOC_BOXES(NEW_BOX_2_IND * BOX_SIZE + DIM) = LOC_BOXES(NEW_BOX_1_IND * BOX_SIZE + 3 + DIM) + 1
                ELSE
                    LOC_BOXES(NEW_BOX_1_IND * BOX_SIZE + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM)
                    LOC_BOXES(NEW_BOX_1_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM)
                    LOC_BOXES(NEW_BOX_2_IND * BOX_SIZE + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM)
                    LOC_BOXES(NEW_BOX_2_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM)
                END IF
            END DO
            ! compute box volumes
            LOC_VOLUMES(MAX_BOX_IND+1) = VOLUME(LOC_BOXES, MAX_BOX_IND, BOX_SIZE)
            LOC_VOLUMES(NEW_BOX_1_IND+1) = VOLUME(LOC_BOXES, NEW_BOX_1_IND, BOX_SIZE)
            LOC_VOLUMES(NEW_BOX_2_IND+1) = VOLUME(LOC_BOXES, NEW_BOX_2_IND, BOX_SIZE)
            TOTAL_NUM_PROCS = TOTAL_NUM_PROCS + 2
            RUNNING_BOX_TOTAL = TOTAL_NUM_PROCS
        END IF

        DO PROC=TOTAL_NUM_PROCS,NUM_TASKS-1
            ! find the box with the biggest volume to split
            ! subtract 1 because we're using zero-indexing
            MAX_BOX_IND = MAXVAL(MAXLOC(LOC_VOLUMES)) - 1
            MAX_DIM = -1
            MAX_DIM_VALUE = -1
            DO DIM=0,2 ! find the longest dim of
                DIM_VALUE = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) - LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM) + 1
                IF (DIM_VALUE > MAX_DIM_VALUE) THEN
                    MAX_DIM_VALUE = DIM_VALUE
                    MAX_DIM = DIM
                END IF
            END DO
            ! split the largest dim and then update
            NEW_BOX_IND = RUNNING_BOX_TOTAL
            DO DIM=0,2
                IF (DIM == MAX_DIM) THEN
                    MID_PT = (LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) + LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM)) / 2
                    LOC_BOXES(NEW_BOX_IND * BOX_SIZE + DIM) = MID_PT + 1
                    LOC_BOXES(NEW_BOX_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM)
                    LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM) = MID_PT
                ELSE
                    LOC_BOXES(NEW_BOX_IND * BOX_SIZE + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + DIM)
                    LOC_BOXES(NEW_BOX_IND * BOX_SIZE + 3 + DIM) = LOC_BOXES(MAX_BOX_IND * BOX_SIZE + 3 + DIM)
                END IF
            END DO
            LOC_VOLUMES(MAX_BOX_IND+1) = VOLUME(LOC_BOXES, MAX_BOX_IND, BOX_SIZE)
            LOC_VOLUMES(NEW_BOX_IND+1) = VOLUME(LOC_BOXES, NEW_BOX_IND, BOX_SIZE)
            RUNNING_BOX_TOTAL = RUNNING_BOX_TOTAL + 1
        END DO

    END SUBROUTINE DECOMPOSE_DOMAIN
    
    ! RK_M1 : rank - 1.
    LOGICAL FUNCTION IS_PT_INSIDE(I, J, K, NBR_EXT_BBOXES, RK_M1, BOX_SIZE)
        INTEGER, INTENT(IN) :: I, J, K, RK_M1, BOX_SIZE
        INTEGER, INTENT(IN) :: NBR_EXT_BBOXES (:)
        INTEGER :: I_MIN, J_MIN, K_MIN, I_MAX, J_MAX, K_MAX
        
        IS_PT_INSIDE = .TRUE.
        
        I_MIN = NBR_EXT_BBOXES(RK_M1 * BOX_SIZE)
        J_MIN = NBR_EXT_BBOXES(RK_M1 * BOX_SIZE + 1)
        K_MIN = NBR_EXT_BBOXES(RK_M1 * BOX_SIZE + 2)
        I_MAX = NBR_EXT_BBOXES(RK_M1 * BOX_SIZE + 3)
        J_MAX = NBR_EXT_BBOXES(RK_M1 * BOX_SIZE + 4)
        K_MAX = NBR_EXT_BBOXES(RK_M1 * BOX_SIZE + 5)

        IF (I < I_MIN .OR. I > I_MAX) IS_PT_INSIDE = .FALSE.
        IF (J < J_MIN .OR. J > J_MAX) IS_PT_INSIDE = .FALSE.
        IF (K < K_MIN .OR. K > K_MAX) IS_PT_INSIDE = .FALSE.

    END FUNCTION IS_PT_INSIDE

    REAL(EB) FUNCTION VOLUME(LOC_BOXES, BOX_IND, BOX_SIZE)
        INTEGER, INTENT(IN) :: LOC_BOXES (:)
        INTEGER, INTENT(IN) :: BOX_IND
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER :: X_MIN, Y_MIN, Z_MIN, X_MAX, Y_MAX, Z_MAX
        INTEGER :: X_LEN, Y_LEN, Z_LEN
        INTEGER :: X_OFFSET, Y_OFFSET, Z_OFFSET
        
        X_OFFSET = 0
        Y_OFFSET = 0
        Z_OFFSET = 0
         
        X_MIN = LOC_BOXES(BOX_IND * BOX_SIZE)
        Y_MIN = LOC_BOXES(BOX_IND * BOX_SIZE + 1)
        Z_MIN = LOC_BOXES(BOX_IND * BOX_SIZE + 2)
        
        X_MAX = LOC_BOXES(BOX_IND * BOX_SIZE + 3)
        Y_MAX = LOC_BOXES(BOX_IND * BOX_SIZE + 4)
        Z_MAX = LOC_BOXES(BOX_IND * BOX_SIZE + 5)

        IF (X_MIN > 0) X_OFFSET = 1
        IF (Y_MIN > 0) Y_OFFSET = 1
        IF (Z_MIN > 0) Z_OFFSET = 1

        X_LEN = X_MAX - X_MIN + X_OFFSET
        Y_LEN = Y_MAX - Y_MIN + Y_OFFSET
        Z_LEN = Z_MAX - Z_MIN + Z_OFFSET

        VOLUME = X_LEN * Y_LEN * Z_LEN
    
    END FUNCTION VOLUME

    INTEGER FUNCTION IC_VAL(I,J,K,IBAR,JBAR,KBAR)
        INTEGER, INTENT(IN) :: I, J, K, IBAR, JBAR, KBAR
        INTEGER :: VAL_1, VAL_2, AREA, PERIM, SIDE3, I_OFFSET
        I_OFFSET=0
        AREA = IBAR * JBAR
        PERIM = 2*IBAR + 2*(JBAR-2)
        SIDE3 = IBAR + 2*(JBAR-2)
        
        VAL_1 = I + (J-1)*IBAR
         
        IF (K == 1) THEN
            IC_VAL = VAL_1
        ELSE IF (K > 1 .AND. K < KBAR) THEN
            IF (J==1) THEN
                VAL_2 = I
            ELSE IF (J < JBAR) THEN
                IF (I > 1) THEN
                   I_OFFSET = 2
                ELSE
                   I_OFFSET = 1
                END IF 
                VAL_2 = IBAR + (J-1)*2 - 2 + I_OFFSET
            ELSE
                VAL_2 = SIDE3 + I
            END IF
                IC_VAL = AREA + VAL_2 + (K-2)*PERIM
        ELSE
            IC_VAL = AREA + (KBAR-2)*PERIM + VAL_1
        END IF
            
    END FUNCTION IC_VAL

    SUBROUTINE LOOP3D_SERIAL_CPU()
    
        ! Compute x-direction flux term FVX
            
        !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
        DO K=LOC_MIN_K_P1,LOC_MAX_K
           DO J=LOC_MIN_J_P1,LOC_MAX_J
              DO I=LOC_MIN_I,LOC_MAX_I
                 WP    = WW(I,J,K)   + WW(I+1,J,K)
                 WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
                 VP    = VV(I,J,K)   + VV(I+1,J,K)
                 VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
                 OMYP  = OMY(I,J,K)
                 OMYM  = OMY(I,J,K-1)
                 OMZP  = OMZ(I,J,K)
                 OMZM  = OMZ(I,J-1,K)
                 TXZP  = TXZ(I,J,K)
                 TXZM  = TXZ(I,J,K-1)
                 TXYP  = TXY(I,J,K)
                 TXYM  = TXY(I,J-1,K)
                 IC    = CELL_INDEX(I,J,K)
                 IEYP  = CELL(IC)%EDGE_INDEX(8)
                 IEYM  = CELL(IC)%EDGE_INDEX(6)
                 IEZP  = CELL(IC)%EDGE_INDEX(12)
                 IEZM  = CELL(IC)%EDGE_INDEX(10)
                 IF (EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
                    OMYP = EDGE(IEYP)%OMEGA(-1)
                    TXZP = EDGE(IEYP)%TAU(-1)
                 ENDIF
                 IF (EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
                    OMYM = EDGE(IEYM)%OMEGA( 1)
                    TXZM = EDGE(IEYM)%TAU( 1)
                 ENDIF
                 IF (EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
                    OMZP = EDGE(IEZP)%OMEGA(-2)
                    TXYP = EDGE(IEZP)%TAU(-2)
                 ENDIF
                 IF (EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
                    OMZM = EDGE(IEZM)%OMEGA( 2)
                    TXYM = EDGE(IEZM)%TAU( 2)
                 ENDIF
                 WOMY  = WP*OMYP + WM*OMYM
                 VOMZ  = VP*OMZP + VM*OMZM
                 RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
                 DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
                 DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
                 TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
                 DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
                 DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
                 TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
                 DTXXDX= RDXN(I)*(TXXP-TXXM)
                 DTXYDY= RDY(J) *(TXYP-TXYM)
                 DTXZDZ= RDZ(K) *(TXZP-TXZM)
                 VTRM  = DTXXDX + DTXYDY + DTXZDZ
                 FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
              ENDDO
           ENDDO
        ENDDO
        
        ! Compute y-direction flux term FVY

        !$OMP TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
        DO K=LOC_MIN_K_P1,LOC_MAX_K
           DO J=LOC_MIN_J,LOC_MAX_J
              DO I=LOC_MIN_I_P1,LOC_MAX_I
                 UP    = UU(I,J,K)   + UU(I,J+1,K)
                 UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
                 WP    = WW(I,J,K)   + WW(I,J+1,K)
                 WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
                 OMXP  = OMX(I,J,K)
                 OMXM  = OMX(I,J,K-1)
                 OMZP  = OMZ(I,J,K)
                 OMZM  = OMZ(I-1,J,K)
                 TYZP  = TYZ(I,J,K)
                 TYZM  = TYZ(I,J,K-1)
                 TXYP  = TXY(I,J,K)
                 TXYM  = TXY(I-1,J,K)
                 IC    = CELL_INDEX(I,J,K)
                 IEXP  = CELL(IC)%EDGE_INDEX(4)
                 IEXM  = CELL(IC)%EDGE_INDEX(2)
                 IEZP  = CELL(IC)%EDGE_INDEX(12)
                 IEZM  = CELL(IC)%EDGE_INDEX(11)
                 IF (EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
                    OMXP = EDGE(IEXP)%OMEGA(-2)
                    TYZP = EDGE(IEXP)%TAU(-2)
                 ENDIF
                 IF (EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
                    OMXM = EDGE(IEXM)%OMEGA( 2)
                    TYZM = EDGE(IEXM)%TAU( 2)
                 ENDIF
                 IF (EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
                    OMZP = EDGE(IEZP)%OMEGA(-1)
                    TXYP = EDGE(IEZP)%TAU(-1)
                 ENDIF
                 IF (EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
                    OMZM = EDGE(IEZM)%OMEGA( 1)
                    TXYM = EDGE(IEZM)%TAU( 1)
                 ENDIF
                 WOMX  = WP*OMXP + WM*OMXM
                 UOMZ  = UP*OMZP + UM*OMZM
                 RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
                 DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
                 DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
                 TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
                 DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
                 DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
                 TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
                 DTXYDX= RDX(I) *(TXYP-TXYM)
                 DTYYDY= RDYN(J)*(TYYP-TYYM)
                 DTYZDZ= RDZ(K) *(TYZP-TYZM)
                 VTRM  = DTXYDX + DTYYDY + DTYZDZ
                 FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*RHO_0(K) - VTRM)
              ENDDO
           ENDDO
        ENDDO
        
        ! Compute z-direction flux term FVZ
        
        !$OMP TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
        DO K=LOC_MIN_K,LOC_MAX_K
           DO J=LOC_MIN_J_P1,LOC_MAX_J
              DO I=LOC_MIN_I_P1,LOC_MAX_I
                 UP    = UU(I,J,K)   + UU(I,J,K+1)
                 UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
                 VP    = VV(I,J,K)   + VV(I,J,K+1)
                 VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
                 OMYP  = OMY(I,J,K)
                 OMYM  = OMY(I-1,J,K)
                 OMXP  = OMX(I,J,K)
                 OMXM  = OMX(I,J-1,K)
                 TXZP  = TXZ(I,J,K)
                 TXZM  = TXZ(I-1,J,K)
                 TYZP  = TYZ(I,J,K)
                 TYZM  = TYZ(I,J-1,K)
                 IC    = CELL_INDEX(I,J,K)
                 IEXP  = CELL(IC)%EDGE_INDEX(4)
                 IEXM  = CELL(IC)%EDGE_INDEX(3)
                 IEYP  = CELL(IC)%EDGE_INDEX(8)
                 IEYM  = CELL(IC)%EDGE_INDEX(7)
                 IF (EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
                    OMXP = EDGE(IEXP)%OMEGA(-1)
                    TYZP = EDGE(IEXP)%TAU(-1)
                 ENDIF
                 IF (EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
                    OMXM = EDGE(IEXM)%OMEGA( 1)
                    TYZM = EDGE(IEXM)%TAU( 1)
                 ENDIF
                 IF (EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
                    OMYP = EDGE(IEYP)%OMEGA(-2)
                    TXZP = EDGE(IEYP)%TAU(-2)
                 ENDIF
                 IF (EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
                    OMYM = EDGE(IEYM)%OMEGA( 2)
                    TXZM = EDGE(IEYM)%TAU( 2)
                 ENDIF
                 UOMY  = UP*OMYP + UM*OMYM
                 VOMX  = VP*OMXP + VM*OMXM
                 RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
                 DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
                 DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
                 TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
                 DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
                 DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
                 TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
                 DTXZDX= RDX(I) *(TXZP-TXZM)
                 DTYZDY= RDY(J) *(TYZP-TYZM)
                 DTZZDZ= RDZN(K)*(TZZP-TZZM)
                 VTRM  = DTXZDX + DTYZDY + DTZZDZ
                 FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(RHO_0(K)+RHO_0(K+1)) - VTRM)
              ENDDO
           ENDDO
        ENDDO    
        
        END SUBROUTINE LOOP3D_SERIAL_CPU


    END PROGRAM LOOP3D
    
    
    
    
    


