    module ModBingham
    contains 
    
    Subroutine ESM_BINGHAM(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, &
      DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME,NPROPS,PROPS,NUMBEROFPHASES,NTENS)


      implicit double precision (a-h, o-z) 
      integer :: NTENS, NSTATEV, NADDVAR, NPROPS, NPT, NOEL, IDSET, NUMBEROFPHASES
      double precision :: EUNLOADING, PLASTICMULTIPLIER
      CHARACTER*80 CMNAME   
      DIMENSION STRESS(NTENS), DSTRAN(NTENS),STATEV(NSTATEV),ADDITIONALVAR(NADDVAR),PROPS(NPROPS)
      

!---Local variables required in standard UMAT
        integer :: IStep, TimeStep
        double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
        double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
        double precision, dimension(:), allocatable :: stran
        double precision, dimension(:), allocatable :: time
        double precision, dimension(:), allocatable :: predef
        double precision, dimension(:), allocatable :: dpred    
        double precision, dimension(:), allocatable :: coords
        double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
        double precision, dimension(:,:), allocatable :: drot
        double precision, dimension(:,:), allocatable :: dfgrd0
        double precision, dimension(:,:), allocatable :: dfgrd1
        double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
        double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
        double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
        double precision :: pnewdt, dtime, temp, dtemp, celent
        double precision :: Value ! auxiliary variable holding any real valued number
        double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation  
              
  
        integer :: ndi, nshr, layer, kspt, kstep, kinc     

!---Local variables defned by the user
! e.g. integer :: var_local	  
!---User can define here additional variables     
        
        allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
              coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )
    
!Initialization
        Eunloading = 0.0
        PlasticMultiplier = 0.0
          
!Rename additional variables
        Porosity = AdditionalVar(1)
        WaterPressure = AdditionalVar(2)
        WaterPressure0 = AdditionalVar(3)
        GasPressure = AdditionalVar(4)
        GasPressure0 = AdditionalVar(5)
        DegreeSaturation = AdditionalVar(6)
        time(1) = AdditionalVar(7)   !TotalRealTime
        time(2) = AdditionalVar(8)   !OverallTotalTime
        dtime = AdditionalVar(9)     !TimeIncrement
        IStep = AdditionalVar(10)    
        TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1   
!Call the UMAT
          call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, & 
           dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
           dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      
!---Definition of Eunloading -> required to define the max time step
      Eunloading = props(5)
!---Always define this value to run the simulation
    
! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
    
          


        return

    end subroutine ESM_BINGHAM
    
     SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)


      implicit double precision (a-h, o-z) 

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV), &
      DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS), &
      STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1), &
      PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


        ! ARGUMENTS:
        !          I/O  TYPE
        !  PROPS    I   R()  : LIST WITH MODEL PARAMETERS
        !  DSTRAN   I   R()  : STRAIN INCREMENT
        !  DDSDDE   O   R(,) : MATERIAL STIFFNESS MATRIX
        !  STRESS  I/O  R()  : STRESSES
        !  STATEV  I/O  R()  : STATE VARIABLES
 
        !---  LOCAL VARIABLES
        ! Elasto-viscoplastic return mapping (Perzyna, n = 1). Yield criterion: sqrt(J2) <= YIELDSTRESS.
        ! Stationary shear flow gives tau = YIELDSTRESS + VISCOSITY * shear rate, for any DTIME.
        INTEGER :: I
        double precision, dimension(NTENS) :: DEV ! deviatoric stress
        double precision :: EEL, ENU, YIELDSTRESS, VISCOSITY, BULKMODULUSLIQUID, LIQUIDPRESSURECAVITATION, ONE, TWO, G
        double precision :: D1, D2
        double precision :: DSTRANVOL, LIQUIDPRESSURE0, LIQUIDPRESSURE
        double precision :: QTRIAL, QRETURN, A


        EEL = PROPS(1)   
        ENU = PROPS(2)
        YIELDSTRESS=PROPS(3)
        VISCOSITY = PROPS(4)
        BULKMODULUSLIQUID=PROPS(5)
        LIQUIDPRESSURECAVITATION=PROPS(6)
        ONE = 1.0D0
        TWO = 2.0D0
        G = EEL/(TWO*(ONE + ENU))

        ! MATERIAL STIFFNESS MATRIX (ELASTIC)
        D1 = BULKMODULUSLIQUID + 4.0D0 * G / 3.0D0
        D2 = BULKMODULUSLIQUID - TWO * G / 3.0D0
        DDSDDE = 0.0D0
        DDSDDE(1:3,1:3) = D2
        DO I = 1, 3
          DDSDDE(I,I) = D1
        END DO
        DO I = 4, NTENS
          DDSDDE(I,I) = G
        END DO

        ! VOLUMETRIC PART (TENSION POSITIVE), CUT OFF AT THE CAVITATION PRESSURE.
        ! Only the updated pressure is checked: a cavitated point can be compressed again.
        ! The cut-off does not change the deviatoric stress.
        DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
        LIQUIDPRESSURE0 = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.0D0
        LIQUIDPRESSURE = MIN(LIQUIDPRESSURE0 + BULKMODULUSLIQUID * DSTRANVOL, LIQUIDPRESSURECAVITATION)

        ! DEVIATORIC ELASTIC PREDICTOR (DSTRAN(4:NTENS) ARE ENGINEERING SHEAR STRAINS)
        DEV(1:3) = STRESS(1:3) - LIQUIDPRESSURE0 + TWO * G * (DSTRAN(1:3) - DSTRANVOL / 3.0D0)
        DEV(4:NTENS) = STRESS(4:NTENS) + G * DSTRAN(4:NTENS)
        QTRIAL = SQRT( 0.5D0 * SUM(DEV(1:3)**2) + SUM(DEV(4:NTENS)**2) ) ! sqrt(J2)

        ! VISCOPLASTIC CORRECTOR: BACKWARD EULER OF (VISCOPLASTIC SHEAR RATE) = (Q - YIELDSTRESS) / VISCOSITY
        IF (QTRIAL .GT. YIELDSTRESS) THEN
          IF (VISCOSITY .GT. 0.0D0) THEN
            A = G * DTIME / VISCOSITY
            QRETURN = (QTRIAL + A * YIELDSTRESS) / (ONE + A)
          ELSE
            QRETURN = YIELDSTRESS ! inviscid limit
          END IF
          DEV = DEV * (QRETURN / QTRIAL)
        END IF

        STRESS(1:3) = LIQUIDPRESSURE + DEV(1:3)
        STRESS(4:NTENS) = DEV(4:NTENS)

     END SUBROUTINE UMAT     
     
    end module ModBingham 