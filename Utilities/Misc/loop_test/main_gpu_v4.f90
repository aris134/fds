PROGRAM LOOP3D

    USE OMP_LIB
    IMPLICIT NONE
    
    ! Miscellaneous declarations
    INTEGER :: ISTEP
    INTEGER, PARAMETER :: NUM_TIME_STEPS = 100
    INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
    INTEGER, PARAMETER :: IBAR = 256, JBAR = 256, KBAR = 256
    INTEGER, PARAMETER :: NEDGE = 12
    INTEGER, PARAMETER :: IBP1 = IBAR+1, JBP1 = JBAR+1, KBP1 = KBAR+1
    REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB
    
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: GX, GY, GZ, RHO_0, X, Y, Z, RDXN, RDYN, RDZN, RDX, RDY, RDZ
    REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: U, V, W, D, RHO, WORK1, WORK2, WORK3, WORK4, WORK5, WORK6
    REAL(EB), POINTER, DIMENSION(:,:,:) :: DP, RHOP, UU, VV, WW, OMY, OMX, OMZ, TXZ, TXY, TYZ
    REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU, FVX, FVY, FVZ
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: CELL_INDEX
    REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
                DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
                DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
                VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
                RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,T_NOW,T_END
    INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,ICMAX,IE,MAX_EDGE,NT
    CHARACTER(LEN=50) :: FILENAME
    
    INTEGER, ALLOCATABLE, DIMENSION(:) :: CELL
    
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: OMEGA_ARR
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: TAU_ARR
    
    ! Write out Starting:
    !$OMP PARALLEL
    !$OMP MASTER
    !$ NT = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP BARRIER
    !$OMP END PARALLEL
    
    WRITE(FILENAME,'(A,I3,A,I3,A,I3,A)') 'loop3d_',IBAR,'GRID_',NT,'THR_',NUM_TIME_STEPS,'STEPS_GPU_V4.txt'
    WRITE(*,*) 'Starting Loop3D, out file: ',TRIM(FILENAME)
    WRITE(*,*) 'Number of detected devices: ',OMP_GET_NUM_DEVICES()
    OPEN(UNIT=10,FILE=TRIM(FILENAME),STATUS='UNKNOWN')
    WRITE(10,*) 'Number of devices=',OMP_GET_NUM_DEVICES()
    WRITE(10,*) 'Starting Loop3D'
    WRITE(10,*) 'IBAR=',IBAR,' JBAR=',JBAR,' KBAR=',KBAR,' OMP_NUM_THREADS=',NT
    
    ! Allocate vars in CPU:
    ALLOCATE(X(0:IBAR), Y(0:JBAR), Z(0:KBAR), RDXN(0:IBAR), RDYN(0:JBAR), RDZN(0:KBAR), RDX(0:IBAR), RDY(0:JBAR), RDZ(0:KBAR))
    ALLOCATE(U(0:IBP1,0:JBP1,0:KBP1), V(0:IBP1,0:JBP1,0:KBP1), W(0:IBP1,0:JBP1,0:KBP1), MU(0:IBP1,0:JBP1,0:KBP1) )
    ALLOCATE(D(0:IBP1,0:JBP1,0:KBP1), RHO(0:IBP1,0:JBP1,0:KBP1), WORK1(0:IBP1,0:JBP1,0:KBP1))
    ALLOCATE(WORK2(0:IBP1,0:JBP1,0:KBP1), WORK3(0:IBP1,0:JBP1,0:KBP1), WORK4(0:IBP1,0:JBP1,0:KBP1))
    ALLOCATE(WORK5(0:IBP1,0:JBP1,0:KBP1), WORK6(0:IBP1,0:JBP1,0:KBP1), DP(0:IBP1,0:JBP1,0:KBP1), RHOP(0:IBP1,0:JBP1,0:KBP1))
    ALLOCATE(RHO_0(0:KBAR), GX(0:IBAR), GY(0:IBAR), GZ(0:IBAR))
    ALLOCATE(FVX(0:IBP1,0:JBP1,0:KBP1), FVY(0:IBP1,0:JBP1,0:KBP1), FVZ(0:IBP1,0:JBP1,0:KBP1))
    ALLOCATE(CELL_INDEX(0:IBP1,0:JBP1,0:KBP1))
    
    ! Initialize:
    DO I=0,IBAR
       X(I)    = REAL(I,EB)
       RDXN(I) = 1.0_EB
       RDX(I)  = 1.0_EB
    ENDDO
    DO J=0,JBAR
       Y(J)    = REAL(J,EB)
       RDYN(J) = 1.0_EB
       RDY(J)  = 1.0_EB
    ENDDO
    DO K=0,KBAR
       Z(K)    = REAL(K,EB)
       RDZN(K) = 1.0_EB
       RDZ(K)  = 1.0_EB
    ENDDO
    
    ! Cell Index, CELL and EDGE:
    CELL_INDEX = 0
    IC = 0
    DO K=1,KBAR
       DO J=1,JBAR
          DO I=1,IBAR
             IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
             IC = IC + 1
             CELL_INDEX(I,J,K) = IC
          ENDDO
       ENDDO
    ENDDO
    ICMAX = IC
    ALLOCATE(CELL(0:NEDGE*ICMAX))
    
    IC = 0
    MAX_EDGE=-1
    DO K=1,KBAR
       DO J=1,JBAR
          DO I=1,IBAR
             IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
             IC = IC + 1
             CELL(IC * NEDGE + 1) = CELL_INDEX(I  ,J  ,K  ) + 1
             CELL(IC * NEDGE + 2) = CELL_INDEX(I+1,J  ,K  ) + 1
             CELL(IC * NEDGE + 3) = CELL_INDEX(I  ,J  ,K  ) + 2
             CELL(IC * NEDGE + 4) = CELL_INDEX(I  ,J+1,K  ) + 2
             CELL(IC * NEDGE + 5) = CELL_INDEX(I  ,J  ,K  ) + 3
             CELL(IC * NEDGE + 6) = CELL_INDEX(I  ,J  ,K-1) + 3
             CELL(IC * NEDGE + 7) = CELL_INDEX(I  ,J  ,K  ) + 4
             CELL(IC * NEDGE + 8) = CELL_INDEX(I  ,J+1,K  ) + 4
             CELL(IC * NEDGE + 9) = CELL_INDEX(I  ,J  ,K  ) + 5
             CELL(IC * NEDGE + 10) = CELL_INDEX(I  ,J  ,K-1) + 5
             CELL(IC * NEDGE + 11) = CELL_INDEX(I  ,J  ,K  ) + 6
             CELL(IC * NEDGE + 12) = CELL_INDEX(I+1,J  ,K  ) + 6
             DO IE=1,NEDGE
                MAX_EDGE = MAX(MAX_EDGE,CELL(IC * NEDGE + IE))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ALLOCATE(OMEGA_ARR(0:4*MAX_EDGE))
    ALLOCATE(TAU_ARR(0:4*MAX_EDGE))

    OMEGA_ARR(0:4*MAX_EDGE) = -1.E6_EB
    TAU_ARR(0:4*MAX_EDGE) = -1.E6_EB

    DO IE=1,MAX_EDGE,2
       DO I=0,3
           OMEGA_ARR(IE * 4 + I) = 1.5E-4_EB
           TAU_ARR(IE * 4 + I) = 2.5E-4_EB
       END DO
    ENDDO
    
    RHO   = 1.19_EB
    RHO_0 = 1.19_EB
    D     = 0.0015_EB
    MU    = 0.0019_EB
    WORK1 = 0.0_EB; WORK2 = 0.0_EB; WORK3 = 0.0_EB; WORK4 = 0.0_EB; WORK5 = 0.0_EB; WORK6 = 0.0_EB
    GX(:) = 0.0_EB; GY(:) = 0.0_EB; GZ(:) = 1.0_EB
    
    ! U, V, W:
    DO K=0,KBAR
       DO J=0,JBAR
          DO I=0,IBAR
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
    
    DO K=0,KBAR
       DO J=0,JBAR
          DO I=0,IBAR
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
    
    

    !$OMP TARGET ENTER DATA MAP(TO:UU(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:VV(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:WW(0:IBP1,0:JBP1,0:KBP1))    &
    !$OMP MAP(TO:OMX(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:OMY(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:OMZ(0:IBP1,0:JBP1,0:KBP1))                   &
    !$OMP MAP(TO:TXY(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:TXZ(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:TYZ(0:IBP1,0:JBP1,0:KBP1))                   &
    !$OMP MAP(TO:CELL_INDEX(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:CELL(0:NEDGE*ICMAX)) MAP(TO:OMEGA_ARR(0:4*MAX_EDGE))                    &
    !$OMP MAP(TO:TAU_ARR(0:4*MAX_EDGE)) MAP(TO:RHOP(0:IBP1,0:JBP1,0:KBP1)) MAP(TO:RDX(0:IBAR)) MAP(TO:RDY(0:JBAR))                &
    !$OMP MAP(TO:RDZ(0:KBAR)) MAP(TO:RDXN(0:IBAR)) MAP(TO:RDYN(0:JBAR)) MAP(TO:RDZN(0:KBAR)) MAP(TO:RHO_0(0:KBAR))                &
    !$OMP MAP(TO:GX(0:IBAR)) MAP(TO:GY(0:IBAR)) MAP(TO:GZ(0:IBAR)) MAP(TO:MU(0:IBP1,0:JBP1,0:KBP1))                               &
    !$OMP MAP(TO:DP(0:IBP1,0:JBP1,0:KBP1))                                                                                        
    T_NOW = OMP_GET_WTIME()
    SIM_LOOP: DO ISTEP = 1, NUM_TIME_STEPS
       CALL LOOP3D_OMP_GPU()
    END DO SIM_LOOP
    T_END = OMP_GET_WTIME()
    !$OMP TARGET EXIT DATA MAP(FROM:FVX(0:IBP1,0:JBP1,0:KBP1),FVY(0:IBP1,0:JBP1,0:KBP1),FVZ(0:IBP1,0:JBP1,0:KBP1))
    
    
    WRITE(10,*) 'Time=',T_END-T_NOW
    WRITE(10,*) 'mean FVX =',SUM(FVX(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
    WRITE(10,*) 'mean FVY =',SUM(FVY(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
    WRITE(10,*) 'mean FVZ =',SUM(FVZ(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
    WRITE(10,*) 'Ending Loop3D'
    CLOSE(10)
    WRITE(*,*) 'Loop3D done.'
    
    CONTAINS
    
    SUBROUTINE LOOP3D_OMP_GPU()
       
       ! Compute x-direction flux term FVX
       
       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
       DO K=1,KBAR
          DO J=1,JBAR
             DO I=0,IBAR
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
                IEYP = CELL(IC * NEDGE + 8)
                IEYM = CELL(IC * NEDGE + 6)
                IEZP = CELL(IC * NEDGE + 12)
                IEZM = CELL(IC * NEDGE + 10)
                IF (OMEGA_ARR(IEYP * 4)>-1.E5_EB) THEN
                   OMYP = OMEGA_ARR(IEYP * 4)
                   TXZP = TAU_ARR(IEYP * 4)
                ENDIF
                IF (OMEGA_ARR(IEYM * 4 + 1)>-1.E5_EB) THEN
                   OMYM = OMEGA_ARR(IEYM * 4 + 1)
                   TXZM = TAU_ARR(IEYM * 4 + 1)
                ENDIF
                IF (OMEGA_ARR(IEZP * 4 + 2)>-1.E5_EB) THEN
                   OMZP = OMEGA_ARR(IEZP * 4 + 2)
                   TXYP = TAU_ARR(IEZP * 4 + 2)
                ENDIF
                IF (OMEGA_ARR(IEZM * 4 + 3)>-1.E5_EB) THEN
                   OMZM = OMEGA_ARR(IEZM * 4 + 3)
                   TXYM = TAU_ARR(IEZM * 4 + 3)
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
       DO K=1,KBAR
          DO J=0,JBAR
             DO I=1,IBAR
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
                IEXP = CELL(IC * NEDGE + 4)
                IEXM = CELL(IC * NEDGE + 2)
                IEZP = CELL(IC * NEDGE + 12)
                IEZM = CELL(IC * NEDGE + 11)
                IF (OMEGA_ARR(IEXP * 4 + 2)>-1.E5_EB) THEN
                   OMXP = OMEGA_ARR(IEXP * 4 + 2)
                   TYZP = TAU_ARR(IEXP * 4 + 2)
                ENDIF
                IF (OMEGA_ARR(IEXM * 4 + 3)>-1.E5_EB) THEN
                   OMXM = OMEGA_ARR(IEXM * 4 + 3)
                   TYZM = TAU_ARR(IEXM * 4 + 3)
                ENDIF
                IF (OMEGA_ARR(IEZP * 4)>-1.E5_EB) THEN
                   OMZP = OMEGA_ARR(IEZP * 4)
                   TXYP = TAU_ARR(IEZP * 4)
                ENDIF
                IF (OMEGA_ARR(IEZM * 4 + 1)>-1.E5_EB) THEN
                   OMZM = OMEGA_ARR(IEZM * 4 + 1)
                   TXYM = TAU_ARR(IEZM * 4 + 1)
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
       DO K=0,KBAR
          DO J=1,JBAR
             DO I=1,IBAR
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
                IEXP = CELL(IC * NEDGE + 4)
                IEXM = CELL(IC * NEDGE + 3)
                IEYP = CELL(IC * NEDGE + 8)
                IEYM = CELL(IC * NEDGE + 7)
                IF (OMEGA_ARR(IEXP * 4)>-1.E5_EB) THEN
                   OMXP = OMEGA_ARR(IEXP * 4)
                   TYZP = TAU_ARR(IEXP * 4)
                ENDIF
                IF (OMEGA_ARR(IEXM * 4 + 1)>-1.E5_EB) THEN
                   OMXM = OMEGA_ARR(IEXM * 4 + 1)
                   TYZM = TAU_ARR(IEXM * 4 + 1)
                ENDIF
                IF (OMEGA_ARR(IEYP * 4 + 2)>-1.E5_EB) THEN
                   OMYP = OMEGA_ARR(IEYP * 4 + 2)
                   TXZP = TAU_ARR(IEYP * 4 + 2)
                ENDIF
                IF (OMEGA_ARR(IEYM * 4 + 3)>-1.E5_EB) THEN
                   OMYM = OMEGA_ARR(IEYM * 4 + 3)
                   TXZM = TAU_ARR(IEYM * 4 + 3)
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
       
       !$OMP TASKWAIT
       END SUBROUTINE LOOP3D_OMP_GPU
    
    END PROGRAM LOOP3D
    
    
    
    
    

