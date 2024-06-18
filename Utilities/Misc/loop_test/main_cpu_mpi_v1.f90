PROGRAM LOOP3D

    USE OMP_LIB
    IMPLICIT NONE
    include 'mpif.h'
    
    ! MPI declarations
    INTEGER :: IERR, MY_RANK, NUM_TASKS, REQUEST

    ! Miscellaneous declarations
    INTEGER :: ISTEP
    INTEGER, PARAMETER :: NUM_TIME_STEPS = 100
    INTEGER, PARAMETER :: BOX_SIZE = 6
    INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
    INTEGER, PARAMETER :: IBAR = 256, JBAR = 256, KBAR = 256
    INTEGER, PARAMETER :: NEDGE = 12
    INTEGER, PARAMETER :: IBP1 = IBAR+1, JBP1 = JBAR+1, KBP1 = KBAR+1
    REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB
    REAL(EB) :: FLUX_AVGS(3)
    
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: GX, GY, GZ, RHO_0, X, Y, Z, RDXN, RDYN, RDZN, RDX, RDY, RDZ
    REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: U, V, W, D, RHO, WORK1, WORK2, WORK3, WORK4, WORK5, WORK6
    REAL(EB), POINTER, DIMENSION(:,:,:) :: DP, RHOP, UU, VV, WW, OMY, OMX, OMZ, TXZ, TXY, TYZ
    REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU, FVX, FVY, FVZ
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: CELL_INDEX
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ALL_TASK_BBOXES
    REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
                DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
                DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
                VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
                RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,T_NOW,T_END,MY_TIME,MIN_TIME,MAX_TIME,AVG_TIME
    INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,ICMAX,IE,MAX_EDGE,MAX_EDGE_P1,NT
    INTEGER :: GLOB_MIN_I, GLOB_MIN_J, GLOB_MIN_K, GLOB_MAX_I, GLOB_MAX_J, GLOB_MAX_K
    INTEGER :: LOC_MIN_I, LOC_MIN_J, LOC_MIN_K, LOC_MAX_I, LOC_MAX_J, LOC_MAX_K, LOC_MAX_I_P1, LOC_MAX_J_P1, LOC_MAX_K_P1, &
               LOC_MIN_I_P1, LOC_MIN_J_P1, LOC_MIN_K_P1
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

    ! Write out Starting:
    !$OMP PARALLEL
    !$OMP MASTER
    !$ NT = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP BARRIER
    !$OMP END PARALLEL
    
    WRITE(*,*) 'Hello from Rank=',MY_RANK

    IF (MY_RANK==0) THEN
        WRITE(FILENAME,'(A,I3,A,I3,A,I3,A)') 'loop3d_',IBAR,'GRID_',NT,'THR_',NUM_TIME_STEPS,'STEPS_CPU_MPI_V1.txt'
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
    GLOB_MAX_I = IBAR
    GLOB_MAX_J = JBAR
    GLOB_MAX_K = KBAR

    ! Decompose the domain
    ALLOCATE(ALL_TASK_BBOXES(NUM_TASKS * BOX_SIZE)) ! AoS format: ..., imin_p, jmin_p, kmin_p, imax_p, jmax_p, kmax_p, imin_(p+1), ...
    CALL DECOMPOSE_DOMAIN(ALL_TASK_BBOXES, GLOB_MIN_I, GLOB_MIN_J, GLOB_MIN_K, GLOB_MAX_I, GLOB_MAX_J, GLOB_MAX_K, BOX_SIZE, NUM_TASKS)
    ! extract the local task bounds from ALL_TASK_BBOXES. These are INCLUSIVE endpoints of the rank's local domain bounding box
    LOC_MIN_I = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE)
    LOC_MIN_J = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 1)
    LOC_MIN_K = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 2)
    LOC_MAX_I = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 3)
    LOC_MAX_J = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 4)
    LOC_MAX_K = ALL_TASK_BBOXES(MY_RANK * BOX_SIZE + 5)
    WRITE(*,*) 'Rank=',MY_RANK,' LOC_MIN_I=',LOC_MIN_I,' LOC_MIN_J=',LOC_MIN_J,' LOC_MIN_K=',LOC_MIN_K,' LOC_MAX_I=',LOC_MAX_I,' LOC_MAX_J=',LOC_MAX_J,' LOC_MAX_K=',LOC_MAX_K

    LOC_MIN_I_P1 = LOC_MIN_I + 1
    LOC_MIN_J_P1 = LOC_MIN_J + 1
    LOC_MIN_K_P1 = LOC_MIN_K + 1
    LOC_MAX_I_P1 = LOC_MAX_I + 1
    LOC_MAX_J_P1 = LOC_MAX_J + 1
    LOC_MAX_K_P1 = LOC_MAX_K + 1
    !TODO: setup halos
    ! Allocate vars in CPU:
    
    ALLOCATE(X(LOC_MIN_I:LOC_MAX_I), Y(LOC_MIN_J:LOC_MAX_J), Z(LOC_MIN_K:LOC_MAX_K))
    ALLOCATE(RDXN(LOC_MIN_I:LOC_MAX_I), RDYN(LOC_MIN_J:LOC_MAX_J), RDZN(LOC_MIN_K:LOC_MAX_K))
    ALLOCATE(RDX(LOC_MIN_I:LOC_MAX_I), RDY(LOC_MIN_J:LOC_MAX_J), RDZ(LOC_MIN_K:LOC_MAX_K))
    ALLOCATE(U(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(V(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(W(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(MU(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(D(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(RHO(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(WORK1(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(WORK2(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(WORK3(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(WORK4(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(WORK5(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(WORK6(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(DP(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(RHOP(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(RHO_0(LOC_MIN_K:LOC_MAX_K),GX(LOC_MIN_I:LOC_MAX_I),GY(LOC_MIN_J:LOC_MAX_J))
    ALLOCATE(GZ(LOC_MIN_K:LOC_MAX_K))
    ALLOCATE(FVX(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(FVY(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(FVZ(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    ALLOCATE(CELL_INDEX(LOC_MIN_I:LOC_MAX_I_P1,LOC_MIN_J:LOC_MAX_J_P1,LOC_MIN_K:LOC_MAX_K_P1))
    
    ! Initialize:
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
    
    ! Cell Index, CELL and EDGE:
    CELL_INDEX = 0
    IC = 0
    DO K=LOC_MIN_K_P1,LOC_MAX_K
       DO J=LOC_MIN_J_P1,LOC_MAX_J
          DO I=LOC_MIN_I_P1,LOC_MAX_I
             IF( .NOT. (ANY( K==(/LOC_MIN_K_P1,LOC_MAX_K/) ) .OR. ANY( J==(/LOC_MIN_J_P1,LOC_MAX_J/) ) .OR. ANY( I==(/LOC_MIN_I_P1,LOC_MAX_I/) )) ) CYCLE
             IC = IC + 1
             CELL_INDEX(I,J,K) = IC
          ENDDO
       ENDDO
    ENDDO
    ALLOCATE(CELL(0:IC))

    IC = 0
    MAX_EDGE=-1
    DO K=LOC_MIN_K_P1,LOC_MAX_K
       DO J=LOC_MIN_J_P1,LOC_MAX_J
          DO I=LOC_MIN_I_P1,LOC_MAX_I
             IF( .NOT. (ANY( K==(/LOC_MIN_K_P1,LOC_MAX_K/) ) .OR. ANY( J==(/LOC_MIN_J_P1,LOC_MAX_J/) ) .OR. ANY( I==(/LOC_MIN_I_P1,LOC_MAX_I/) )) ) CYCLE
             IC = IC + 1
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
    ALLOCATE(EDGE(0:MAX_EDGE))
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
    
    ! U, V, W:
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
    
    ! Compute Tau OMG:
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
    
    IF (MY_RANK==0) WRITE(*,*) 'BEGINNING SIM LOOP...'
    T_NOW = MPI_WTIME()
    SIM_LOOP: DO ISTEP = 1, NUM_TIME_STEPS
       ! WRITE(*,*) 'TIME STEP=',ISTEP
       CALL LOOP3D_OMP_CPU()
    END DO SIM_LOOP
    T_END = MPI_WTIME()
    CALL MPI_Barrier(MPI_COMM_WORLD, IERR)
    IF (MY_RANK==0) WRITE(*,*) 'FINISHED SIM LOOP...'
    
    ! Global reduction to compute the mean flux terms
    FLUX_AVGS(1) = SUM(FVX(LOC_MIN_I_P1:LOC_MAX_I,LOC_MIN_J_P1:LOC_MAX_J,LOC_MIN_K_P1:LOC_MAX_K))
    FLUX_AVGS(2) = SUM(FVY(LOC_MIN_I_P1:LOC_MAX_I,LOC_MIN_J_P1:LOC_MAX_J,LOC_MIN_K_P1:LOC_MAX_K))
    FLUX_AVGS(3) = SUM(FVZ(LOC_MIN_I_P1:LOC_MAX_I,LOC_MIN_J_P1:LOC_MAX_J,LOC_MIN_K_P1:LOC_MAX_K))

    CALL MPI_REDUCE(FLUX_AVGS, FLUX_AVGS, 3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERR)
    
   !TODO calculate min, max, avg time
    MY_TIME = T_END - T_NOW
    CALL MPI_REDUCE(MY_TIME, MIN_TIME, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_REDUCE(MY_TIME, MAX_TIME, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, IERR)
    CALL MPI_REDUCE(MY_TIME, AVG_TIME, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, IERR)


    IF (MY_RANK == 0) THEN
        FLUX_AVGS(1) = FLUX_AVGS(1) / (GLOB_MAX_I*GLOB_MAX_J*GLOB_MAX_K)
        FLUX_AVGS(2) = FLUX_AVGS(2) / (GLOB_MAX_I*GLOB_MAX_J*GLOB_MAX_K)
        FLUX_AVGS(3) = FLUX_AVGS(3) / (GLOB_MAX_I*GLOB_MAX_J*GLOB_MAX_K)
        AVG_TIME = AVG_TIME / NUM_TASKS
        WRITE(10,*) 'Time(avg,min,max)=',AVG_TIME,MIN_TIME,MAX_TIME
        WRITE(10,*) 'mean FVX =',FLUX_AVGS(1)
        WRITE(10,*) 'mean FVY =',FLUX_AVGS(2)
        WRITE(10,*) 'mean FVZ =',FLUX_AVGS(3)
        WRITE(10,*) 'Ending Loop3D'
        CLOSE(10)
        WRITE(*,*) 'Loop3D done.'
    END IF

    CALL MPI_FINALIZE(IERR)

    CONTAINS
    
    SUBROUTINE LOOP3D_OMP_CPU()
       
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
         !$OMP END DO NOWAIT
         
         !$OMP END PARALLEL
    END SUBROUTINE LOOP3D_OMP_CPU

    SUBROUTINE DECOMPOSE_DOMAIN(LOC_BOXES, MIN_I, MIN_J, MIN_K, MAX_I, MAX_J, MAX_K, BOX_SIZE, NUM_TASKS)
        INTEGER, INTENT(IN OUT) :: LOC_BOXES (:) ! AoS
        INTEGER, INTENT(IN) :: MIN_I, MIN_J, MIN_K, MAX_I, MAX_J, MAX_K
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER, INTENT(IN) :: NUM_TASKS
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
            LOC_VOLUMES(MAX_BOX_IND) = VOLUME(LOC_BOXES, MAX_BOX_IND, BOX_SIZE)
            LOC_VOLUMES(NEW_BOX_1_IND) = VOLUME(LOC_BOXES, NEW_BOX_1_IND, BOX_SIZE)
            LOC_VOLUMES(NEW_BOX_2_IND) = VOLUME(LOC_BOXES, NEW_BOX_2_IND, BOX_SIZE)
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
            LOC_VOLUMES(MAX_BOX_IND) = VOLUME(LOC_BOXES, MAX_BOX_IND, BOX_SIZE)
            LOC_VOLUMES(NEW_BOX_IND) = VOLUME(LOC_BOXES, NEW_BOX_IND, BOX_SIZE)
            RUNNING_BOX_TOTAL = RUNNING_BOX_TOTAL + 1
        END DO

    END SUBROUTINE DECOMPOSE_DOMAIN
    
    REAL(EB) FUNCTION VOLUME(LOC_BOXES, BOX_IND, BOX_SIZE)
        INTEGER, INTENT(IN) :: LOC_BOXES (:)
        INTEGER, INTENT(IN) :: BOX_IND
        INTEGER, INTENT(IN) :: BOX_SIZE
        INTEGER :: X_LEN, Y_LEN, Z_LEN

        X_LEN = LOC_BOXES(BOX_IND * BOX_SIZE + 3) - LOC_BOXES(BOX_IND * BOX_SIZE) + 1
        Y_LEN = LOC_BOXES(BOX_IND * BOX_SIZE + 4) - LOC_BOXES(BOX_IND * BOX_SIZE + 1) + 1
        Z_LEN = LOC_BOXES(BOX_IND * BOX_SIZE + 5) - LOC_BOXES(BOX_IND * BOX_SIZE + 2) + 1
        
        VOLUME = X_LEN * Y_LEN * Z_LEN
    
    END FUNCTION VOLUME

    END PROGRAM LOOP3D
    
    
    
    
    


