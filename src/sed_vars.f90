! This module declares the main variables used by the sediment transport
! model program.

! By: Luis Morales

! Release: Jun-2017, GIWS

module sed_vars

    implicit none

    ! Variables from MESH (TEMPORALY INITIALIZED HERE !!!)
    ! integer, parameter :: ngru = 5

    !integer, parameter :: NA = 24, ntsteps = 24, yCount = 4, xCount = 13
    !real, parameter :: DELX = 1000., DELY = 1000. !DELT = 30*60
    !real, dimension(NA) :: imperCellAre
    !integer, dimension(NA) :: ipos, jpos
!   real, dimension(yCount, xCount) :: rank = reshape( (/ 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	19.0,	20.0,	21.0,	11.0,	18.0,	0.0, &
!                                                0.0,	1.0,	3.0,	0.0,	12.0,	15.0,	16.0,	17.0,	6.0,	22.0,	23.0,	24.0,	0.0, &
!                                                0.0,	4.0,	7.0,	10.0,	13.0,	14.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	25.0, &
!                                                2.0,	5.0,	8.0,	9.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0 /), &
!                                                shape(rank), order=(/2,1/))

!    real, dimension(yCount, xCount) :: next = reshape( (/ 0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	20.0,	21.0,	23.0,	24.0,	24.0,	0.0, &
!                                                0.0,	7.0,	7.0,    0.0,	14.0,	16.0,	17.0,	20.0,	22.0,	23.0,	24.0,	25.0,	0.0, &
!                                                0.0,	7.0,	10.0,	13.0,	14.0,	16.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0, &
!                                                5.0,	8.0,	10.0,	10.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0 /), &
!                                                shape(next), order=(/2,1/))
!    real, dimension(yCount, xCount) :: dummy
!    character(len=8), dimension(yCount, xCount) :: dummyC

!    integer, dimension (yCount,xCount) :: cell = reshape( (/ 0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0, &
!                                                0,	1,	1,	1,	1,	1,	0,	1,	0,	1,	1,	1,	0, &
!                                                0,	1,	1,	1,	0,	1,	1,	0,	1,	0,	1,	1,	0, &
!                                                0,	0,	0,	0,	0,	0,	1,	0,	1,	0,	0,	0,	0 /), &
!                                                shape(cell), order=(/2,1/))

!    real, dimension (yCount,xCount) :: precip = reshape( (/ 0.1,	0.1,	0.1,	0.1,	0.2,	0.1,	0.3,	0.1,	0.0,	0.0,	0.1,	0.2,	0.1, &
!                                                0.3,	0.2,	0.6,	0.4,	0.3,	1.0,	0.2,	0.1,	0.2,	0.1,	0.1,	0.1,	0.4, &
!                                                0.4,	1.0,	1.0,	0.6,	0.4,	0.9,	1.0,	0.5,	1.0,	0.0,	1.0,	0.1,	0.2, &
!                                                0.1,	0.1,	0.3,	0.4,	0.6,	0.5,	1.0,	0.4,	1.0,	0.6,	0.3,	0.4,	0.5 /), &
!                                                shape(cell), order=(/2,1/))


    !real, dimension(NA) :: waterDepth, precip, evapotran, waterSslope, flowWidth, flowVelocity, discharge
    !real, dimension(NA) :: waterDepthB, waterSslopeB, flowWidthB, flowVelocityB, dischargeB

    ! Paramenters contained in 'MESH_drainage_database.r2c' file. See fileL:
    ! 'H:\MESH\r1064\Driver\MESH_Driver\shd_variables.f90'

    ! Path to folder with input files
    character(len=80), parameter :: filenamePR = "PREC_H.r2c", &
                                    filenameEV = "EVAP_H.r2c", &
                                    filenameRBM = "rbm_input.r2c", &
                                    filenameWR = "WR_RUNOVF_H.r2c"
    character(len=80) :: casename, casedir, OUTFIELDfolder, dateTime,  filename, &
                        dateIn, date1PR, date1EV, date1WR,  date2PR, date2EV, date2WR, &
                        date1RBM1, date2RBM1, date1RBM2, date2RBM2, date1RBM3, date2RBM3, &
                        date1RBM4, date2RBM4, date1RBM5, date2RBM5, date1RBM6, date2RBM6
    character(len=150) :: MESHdir, filepath, filepathOUTFIELD, filepathRBM
    integer :: DELTac
    type GridParams
         !* NA: Total number of grids in the basin. [-]
        integer :: NA = 0

        !* NAA: Total number of contributing grids in the basin. [-]
        integer :: NAA = 0

        real :: &
            xOrigin = 0.0, &
            yOrigin = 0.0

        real :: AL = 0.0

        integer :: &
            xCount = 0, &
            yCount = 0

        real :: &
            xDelta = 0.0, &
            yDelta = 0.0

        !> RANK: Rank of the grid by order of elevation and streamflow direction. [-]
        real, dimension(:,:), allocatable :: rank

        !> NEXT: Rank of the grid that is immediately downstream of this grid. [-]
        real, dimension(:,:), allocatable :: next

        !* ELEV: Elevation of the grid. [m]
        real, dimension(:,:), allocatable :: ELEV

        !* AREA: Area of the grid. [m^2]
        real, dimension(:,:), allocatable :: AREA

        !* SLOPE_CHNL: Channel slope
        real, dimension(:,:), allocatable :: SLOPE_CHNL

        !* CHNL_LEN: Channel length (may exist in other files as 'chnllength', 'rl', or 'ch_length')
        real, dimension(:,:), allocatable :: CHNL_LEN
    end type GridParams

    type outputPointInfo
        character(len=150) :: outputDir
        integer :: nGrPoint
    end type outputPointInfo
    type(outputPointInfo) :: opin

    type, extends(outputPointInfo) :: outputFieInfo
        integer :: nField
    end type outputFieInfo
    type (outputFieInfo) :: ofin

    type outGrPoint
        integer :: grNumb, idsedclass
        character(len=15) :: varname, sedname
    end type outGrPoint
    type(outGrPoint), dimension(:), allocatable :: ogrp, ofie



    type GridSize
        real :: DELX, DELY
    end type GridSize

    type IterDate
        integer :: year, month, day, jday, hour, mins, secs
    end type IterDate
    type(IterDate) :: startDate, stopDate, currDate

    type(GridParams) :: shd
    type(GridSize), dimension(:), allocatable :: grs
    integer :: NA, yCount, xCount
    !integer, parameter :: ntsteps = 24
    !real :: DELX, DELY
    real, dimension(:), allocatable:: imperCellAre
    integer, dimension(:), allocatable :: ipos, jpos
    integer, dimension(:,:), allocatable :: rank, next
    real, dimension(:,:), allocatable  :: dummy, variMat, mat1PR, mat1EV, mat1WR, &
                                        mat2PR, mat2EV, mat2WR, mat1RBM1, mat2RBM1, &
                                        mat1RBM2, mat2RBM2, mat1RBM3, mat2RBM3, &
                                        mat1RBM4, mat2RBM4, mat1RBM5, mat2RBM5, &
                                        mat1RBM6, mat2RBM6
    integer, dimension(:,:), allocatable  :: dummyI
    character(len=8), dimension(:,:), allocatable  :: dummyC

    type meteorologicalVariables
        real :: precip, evapotran
    end type meteorologicalVariables
    !type(meteorologicalVariables), dimension(NA) :: mv
    type(meteorologicalVariables), dimension(:), allocatable :: mv
    character(len=10) :: metVarName


    type :: overlandFlowHydraulics
        real :: slope, discharge, velocity, width, depth
    end type overlandFlowHydraulics
    !type(overlandFlowHydraulics), dimension(NA) :: ofh, ofhB
    type(overlandFlowHydraulics), dimension(:), allocatable :: ofh, ofhB

    type, extends(overlandFlowHydraulics) :: reachHydraulics
        real :: length
    end type reachHydraulics
    !type(reachHydraulics), dimension(NA) :: rh
    type(reachHydraulics), dimension(:), allocatable :: rh
    real, dimension(:,:), allocatable :: DISC, INFL, QDIF, DEPT, WIDT, VELO


    ! stability parameters in the finite difference scheme
    real :: theta, theta_r, phi

    ! Define constant
    real, parameter :: gravi = 9.81
    real, parameter :: pi = 3.14159265359
    real, parameter :: rhow = 1000.0
    real, parameter :: vis = 1.2e-6         ! (m^2/s)
    !real, parameter :: theta = 0.65, theta_r = 0.5 , phi = 0.5 ! stability parameters in the finite difference scheme
    integer, parameter :: nsedpar = 14 ! Number of possible sediment particle diameters.
    integer, parameter :: nchbedly = 2 ! Number of channel bed layers
    real, parameter :: parDLim = 0.062 ! Particle diamenter limit (mm) between fine and non-fine sediments.
    integer, parameter :: DELTout = 3600 ! Time step at which met and hydro vars are read from MESH.

    ! Dummy variables
    integer :: i, j, k, ios
    integer :: iosPR = 0, iosEV = 0, iosWR = 0, iosRBM

    ! Input file units
    !integer, parameter :: unitpar = 1
    integer, parameter :: unitParam = 10
    integer, parameter :: unitSoVeCha = 20
    integer, parameter :: unitGrid = 30
    integer, parameter :: unitDrainDB = 40
    integer, parameter :: unitPR = 60
    integer, parameter :: unitEV = 80
    integer, parameter :: unitWR = 100
    integer, parameter :: unitRBM = 120
    !integer, parameter :: unitOUTFIELDprecip = 50
    !integer, parameter :: unitOUTFIELDevap = 60
    !integer, parameter :: unit_rbm_input = 70
    !integer, parameter :: unitOUTFIELDwr_runovf = 80

    ! Declaring input data variables
    type :: soilParti
        character(len=15) :: name
        real :: minD, maxD, meanD
    end type soilParti
    type(soilParti), dimension(nsedpar) :: sp

    integer :: overlFlowCapaMethod, instreamFlowCapaMethod
    real :: FPCRIT ! FPCRIT=Maximun sustainable concentration of fine particles (m3/m3)
    integer :: DELT


    type :: gridCellAttrib
        real :: iniSoilDepth, sedFluxBondC, cellCanopyCov, cellGroundCov
        integer :: soilType
        character(len=8) :: cellVege
    end type gridCellAttrib
    type(gridCellAttrib), dimension(:), allocatable  :: gca

    type :: channelBebAttrib
        real :: iniThickBed
        integer :: soilTypeBed, soilTypeBank
    end type channelBebAttrib
    type(channelBebAttrib), dimension(:), allocatable  :: cba

    integer :: NSOIL
    type soilAttrib
        real, dimension(nsedpar) :: frac
        real :: density, porosity, soilDetach, overlandDetach, chanBankDetach
    end type soilAttrib
    type(soilAttrib), dimension(:), allocatable :: sa

    integer :: NVEGE
    type :: vegeAttrib
        character(len=9) :: name
        real :: fallHeight, dropDiam, percDrip
    end type vegeAttrib
    type(vegeAttrib), dimension(:), allocatable :: va

    type :: chbedlayer
        real :: thick, D16, D50, D84, D99
        real,dimension(nsedpar) :: frac
    end type chbedlayer

    real :: threChanCon, maxBedThick, ratioStress
    !integer :: NBEDSOIL
    type :: chanBedSoilAttrib
        type(chbedlayer),dimension(nchbedly) :: ly
        real, dimension(nsedpar) :: frac
        real :: density, porosity, thick
    end type chanBedSoilAttrib
    !type(chanBedSoilAttrib), dimension(:), allocatable :: cbsa


    ! Set-up variables
    type, extends(soilAttrib)  :: soilCellAttr
        real :: diameter
    end type  soilCellAttr
    type(soilCellAttr), dimension(:), allocatable :: sca, bsca
    type(vegeAttrib), dimension(:), allocatable :: vca
    type(chanBedSoilAttrib), dimension(:), allocatable :: cbsca, cbscaB

!    type, extends(chanBedSoilAttrib) :: chanBedSoilCellAttrib
!        real, dimension(14) :: frac
!        real :: density, porosity
!    end type chanBedSoilAttrib
!    type(chanBedSoilCellAttrib), dimension(:), allocatable :: cbsca




    !
!    type :: Kr_soil
!        real, parameter :: Kr_clay, Kr_siltyClay, Kr_siltyClayLoam, &
!                           Kr_silt, Kr_siltLoam, Kr_loam, Kr_sandyLoam, &
!                           Kr_sand
!    end type Kr_soil
!    type(Kr_soil) :: KrSoils


    ! soil type
    !type :: soil
    !    character(len=9) :: name
    !    real :: diameter, density, surfPorosity, soilDetach, overlandDetach
    !end type soil

    ! Leaf drip from vegetation
    !type :: vegeta
    !    character(len=9) :: name
    !    real :: fallHeight, dropDiam, percDrip
    !end type vegeta

!    ! Vegetation data in each cell
!    type :: cellVegeta
!        character(len=8) :: soilname, vegename
!        real :: dropDiam, percDrip, fallHeight, soilDetach, groundCov, canopyCov
!    end type cellVegeta
!
!    ! Vegetation data in each cell
!    type :: cellVegeta
!        character(len=8) :: soilname, vegename
!        real :: dropDiam, percDrip, fallHeight, soilDetach, groundCov, canopyCov
!    end type cellVegeta


    ! Values at cell edges
    type :: cellNeighbor
        integer :: east, north, west, south
    end type cellNeighbor

    ! Fluxes at cell edges
    type :: fluxNeighbor
        real :: east, north, west, south
    end type fluxNeighbor

    type, extends(fluxNeighbor) :: inOutCellFlux
        integer :: nfluxes
        integer :: s_east, s_north, s_west, s_south
    end type inOutCellFlux

    type :: cellSedAtt
        real :: pot_Dz, ava_Dz, Dz, SD
        real, dimension(nsedpar) :: C
    end type cellSedAtt

    type :: inOutNode
        integer :: reachOut, NreachIn
        integer, dimension(:), allocatable :: reachIn
    end type inOutNode



    type :: inStreamRoutingReach
        real, dimension(nsedpar) :: Vs, G, Vs_up, G_up, bke, C
    end type inStreamRoutingReach

    type :: inStreamRoutingNode
        real, dimension(nsedpar) :: gs, Dz
        !real :: SD
    end type inStreamRoutingNode




!    !Soil detachment by rainfall impact
!    type :: TrainSoilEr
!        real :: soilVul     ! Soil vulverability to rainfall detach.
!        real :: damEffecSw  ! Damping efectiveness of surface water.
!    end type TrainSoilEr
!


    ! Define general input parameters
    !integer :: nsoilType, nvegeType
    !type(soil), allocatable :: soilChar(:)
    !type(vegeta), allocatable :: vegeChar(:)
    !character(len=8), allocatable :: cellSoil(:)
    !character(len=8), allocatable :: cellVege(:)
    real, allocatable :: D_R(:), D_F(:), G(:)
    !type(soil), allocatable :: cellSoilChar(:)
    !type(vegeta), allocatable :: cellVegeChar(:)
    type(cellNeighbor), allocatable :: cn(:)
    type(fluxNeighbor), allocatable ::  G_out(:), waterDepth_edge(:), waterSslope_edge(:), flowVelocity_edge(:), discharge_edge(:),&
                                        waterDepthB_edge(:), waterSslopeB_edge(:), flowVelocityB_edge(:), dischargeB_edge(:)
    type(inOutCellFlux), dimension(:), allocatable :: FB, F
    type(cellSedAtt), dimension(:), allocatable :: csaB, csa

    ! In-stream routing
    type(inOutNode), dimension(:), allocatable :: ion
    type(inStreamRoutingReach), dimension(:), allocatable :: isrr, isrrB
    type(inStreamRoutingNode), dimension(:), allocatable :: isrn

    !real, allocatable :: C(:)


!
!    ! Soil detachment by rainfall impact
!    type(TrainSoilEr), allocatable :: rainSoilEr(:)
!    real, allocatable, rainDropEro(:,:,:)




    !type(cellEdges), dimension (yCount,xCount) :: waterDepth_edge, waterSslope_edge, flowVelocity_edge, discharge_edge





end module sed_vars
