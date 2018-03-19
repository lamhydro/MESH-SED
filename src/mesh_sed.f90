program mesh_sed

    ! Call modules
    !use sed_tools
    use sed_vars
    use sed_config
    use sed_input
    use sed_output
    use sed_rainDropDetachment
    use sed_overlandFlowDetachment
    use sed_overlandFlowTransCapa
    use sed_hillslopeRouting
    use sed_inStreamRouting
    use sed_chanBankErosion

    ! Declare internal variables
    !integer :: titer

    print *,'!------------------------------------------------------!'
    print *,'! STARTING : MESH-SED                                  !'
    print *,'!------------------------------------------------------!'


    !> Read the study case name
    call mesh_sed_case

    ! Read 'MESH_sed_options.r2c'
    !call read_sed_options

    !> Read 'MESH_sed_param.ini'
    call read_sed_param

    !> Read 'MESH_drainage_database.r2c'
    call read_MESH_drainage_database

    !> Allocate memory space of important variables
    call sed_allocate_var

    !> Get (i,j) position of RANK cells and neighbor cells.
    call sed_config_init

    !> Read 'MESH_sed_gridCell.ini'
    call read_sed_gridCell

    !> Read 'MESH_sed_soilVegChan.ini'
    call read_sed_soilVegChan





    !print*, cn(24)%east, cn(24)%north, cn(24)%west, cn(24)%south
    !print *,'!------------------------------------------------------!'
    !print *,'! Reading ', 'MESH_parameters_sediment.ini'
    !print *,'!------------------------------------------------------!'
    !call read_MESH_parameters_sediment

!    do i=1,nsoilType
!        print *, soilChar(i)%name,  soilChar(i)%diameter, soilChar(i)%density, soilChar(i)%surfPorosity, soilChar(i)%soilDetach
!    end do

!    do i = 1, nvegeType
!        print *, vegeChar(i)%name,  vegeChar(i)%fallHeight, vegeChar(i)%dropDiam, vegeChar(i)%percDrip
!    end do

!    do i = 1, NA
!        !write(*,'(13(a9))') cellVege(k)
!        print*, cellVege(i)
!    end do

    !do i=1,yCount
    !    print *,  cell(i,:)
    !end do


    ! Initialize some variables
    !> Set the x and y dimension of each grid cell.
    call setGridCellSize

    !> Set land, channel and bank soil and vegetation attributes for each active cell
    call setSoilAndVegeParamInCell()
    !do i = 1, NA
        !do j = 1,xCount
        !print *,  cellSoilChar(i)%name
        !end do
    !end do

    !> Set the inflow and outflow reaches for each active cell
    call setInAndOutStreamToNode()
    !do i=1,NA
        !if (allocated(ion(i)%reachIn)) print*, ion(i)%reachIn
        !print*, ion(i)%reachIn(1)
        !print*, ion(i)%reachOut
    !end do

    !> Set initial conditions in the cell
    call sed_initial_conditions



    ! Opening Meteo. variables files
    !filepathOUTFIELD = trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/"
    !filenamePR = "PREC_H.r2c"
    !call openFile(unitOUTFIELDprecip, filepath, filename) ! Should be 'PREC_H' for hourly PREC
    !filenameEV = "EVAP_H.r2c"
    !call openFile(unitOUTFIELDevap, filepath, filename) ! Should be 'EVAP_H' for hourly EVAP

    ! Opening rbm_input file for in-stream flow variables
    !filepathRBM = trim(ADJUSTL(MESHdir))
    !filenameRBM = "rbm_input.r2c"
    !call openFile(unit_rbm_input, filepath, filename)

    ! Opening WR_RUNOVF_H.r2c for overland flow variables
    !filepath = trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/"
    !filenameWR = "WR_RUNOVF_H.r2c"
    !call openFile(unitOUTFIELDwr_runovf, filepath, filename)


    !> Open the files to write time series at different grid points
    call openFilesForGrPointTS()

    ! MAIN LOOP: ITERATION THROUGH TIME
10  format('"',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',I2.2,':',I2.2,'.000"')
    do

        call currDate_update
        !print *, currDate%year, ' ', currDate%month,' ',  currDate%day, ' ',currDate%jday, &
        !        ' ', currDate%hour, ' ',currDate%mins, ' ', currDate%secs

        !write (dateIn, '( '"',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',I2.2,':',I2.2,'.000"' )') currDate%year, &
        write (dateIn, 10) currDate%year, &
        currDate%month,  currDate%day,  currDate%hour, currDate%mins, currDate%secs
        print *, 'DATE & TIME: ', dateIn

        !dateIn = '"'//"2002/01/30      1:15:20.000"//'"'
        !call parseDateTime(dateIn, currDate)
        !stop


        !> READING MET DATA
        call readMetData

        !> READING OVERLAND FLOW DATA
        call readOverLandFlowData

        !> READING IN-STREAM FLOW DATA
        call readInStreamFlowData

        !-----------------------------------------------------
        ! HILLSLOPE TRANSPORT
        !-----------------------------------------------------
        !> Raindrop detachment
        call rainDropDetachCell()
        !do i=1,NA
            !write(*,'(13(F15.2))')  D_R(i,:)
        !    print*, mv(i)%precip, ' ', D_R(i)
        !end do

        !> Overland flow detachment
        call overlandFlowDetachCell()
        !do i=1,NA
            !write(*,'(13(F15.2))')  D_F(i,:)
        !    pr

        ! Calculate variables at the cell edges
         call varsAtCellEdge()

        ! Hillslope routing
        !- 1) Estimate transport capacity at the outlet of the cell
        !print*,
        !print*, '----1----'
        call overlandFlowTransCapa_outCell()
        !do i=1,NA
            !write(*,'(13(F15.2))')  G(i,:)
        !    print*, i, G_out(i)%east, G_out(i)%north, G_out(i)%west, G_out(i)%south
        !end do

        ! Initilize variables
        call sed_before_vars()

        !- 2) Estimate potential sediment concentration
        !print*,
        !print*, '----2----'
        call potentialSedConc()
        !do i=1,NA
            !write(*,'(13(F15.2))')  G(i,:)
            !print*, i, F(i)%nfluxes, F(i)%east, F(i)%north, F(i)%west, F(i)%south !C(i)
        !    print*, csa(i)%C
        !end do

        !- 3) Estimate potential change in surface elevation
        !print*,
        !print*, '----3----'
        call potentialChanSufEle()
        !do i=1,NA
        !    print*, csa(i)%pot_Dz
        !end do

        !- 4) Estimate available change in surface elevation
        !print*,
        !print*, '----4----'
        call availableChanSufEle()
        !do i=1,NA
        !    print*, csa(i)%ava_Dz
        !end do

        !- 5) Compare pot_Dz vs. ava_Dz
        !print*,
        !print*, '----5----'
        call cellConAndChanSufEle()
        !do i=1,NA
        !    print*, csa(i)%Dz, csa(i)%SD
        !end do

        !-----------------------------------------------------
        ! IN-STREAM TRANSPORT
        !-----------------------------------------------------
        ! Instream routing
        ! - Longitudinal sediment velocity
        call longSedVelocity_grid

        ! - Sediment inflow to channel
        ! -- Channel bank erosion
        call grid_bankErosion
        ! -- Sediment inflows
        call flowSedInput

        ! - In-stream sediment routing
        call inStreamRouting
        !do i=1,NA
        !    print*, isrr(i)%G
        !end do
        ! - Update frac diameters in active and parent layers
        call setFracDiame

        !-----------------------------------------------------
        ! CHECKIG THE MASS BALANCE
        !-----------------------------------------------------


        !-----------------------------------------------------
        ! WRITE OUT MODEL PRONOSTIC VARIABLES
        !-----------------------------------------------------
        call writeGrPointTS




        if (stopExec()) exit

    end do

    call closeFilesForGrPointTS()

    print *,'!------------------------------------------------------!'
    print *,'! ENDING : MESH-SED                                    !'
    print *,'!------------------------------------------------------!'


!    ! MAIN LOOP: ITERATION THROUGH TIME
!    titer = 1
!    do
!
!
!        !> READING MET DATA
!        call readMetData
!
!        !> READING OVERLAND FLOW DATA
!        call readOverLandFlowData
!
!
!        !> READING IN-STREAM FLOW DATA
!        call readInStreamFlowData
!
!
!
!
!
!
!
!!        ! The data set-up in this loop must come from the hydrological model
!!        do i = 1, NA
!!        !do i =1,yCount
!!        !    do j = 1,xCount
!!
!!            ! Met. data
!!            !precip(i) = 10. !rand()*500.              ! Precipitation in mm/h
!!            !evapotran(i) = precip(i)*0.6/3.6e6 ! Evapotranspiration in m/s
!!
!!            ! Met. data
!!            !mv(i)%precip = 10. !rand()*500.              ! Precipitation in mm/h
!!            !mv(i)%evapotran = mv(i)%precip*0.6/3.6e6 ! Evapotranspiration in m/s
!!
!!            ! Overland flow data
!!            !waterDepth(i) = 10. !rand()*20             ! Water depth in mm
!!            !waterSslope(i) = 0.01               ! Water surface slope ~ bottom slope. It is read once!!
!!            !flowWidth(i) = sqrt(DELX**2 + DELY**2)                 ! Flow width = cell width (m). It is read once!!
!!            !flowVelocity(i) = 0.01                ! Water flow velocity (m/s).
!!            !discharge(i) = waterDepth(i)*1.e-3*flowWidth(i)*flowVelocity(i) ! Discharge in the cell in m3/s
!!
!!            ! Overland flow data
!!            ofh(i)%depth = 10. !rand()*20             ! Water depth in mm
!!            ofh(i)%slope = 0.01               ! Water surface slope ~ bottom slope. It is read once!!
!!            ofh(i)%width = sqrt(DELX**2 + DELY**2)                 ! Flow width = cell width (m). It is read once!!
!!            ofh(i)%velocity = 0.01                ! Water flow velocity (m/s).
!!            ofh(i)%discharge = ofh(i)%depth*1.e-3*ofh(i)%width*ofh(i)%velocity ! Discharge in the cell in m3/s
!!
!!            ! Instream flow data
!!            rh(i)%slope = 0.01 ! In-stream water surface slope. It is aproximate to the channel bottom slope.
!!            rh(i)%velocity = 0.02! Channel flow velocity (m/s)
!!            rh(i)%width = 10. ! Channel width (m)
!!            rh(i)%depth = 1.! Channel water depth (m)
!!            rh(i)%discharge =rh(i)%velocity*(rh(i)%width*rh(i)%depth) ! In-stream discharge (m3/s). Assume a rectagular section
!!            rh(i)%length = DELX ! Channel reach length
!!
!!        !    end do
!!        !end do
!!        end do
!!        do i = 1, NA
!!            !write(*,'(13(f10.2))')  waterDepth(i,:)
!!            print*, waterDepth(i)
!!        end do
!
!
!
!
!        !-----------------------------------------------------
!        ! HILLSLOPE TRANSPORT
!        !-----------------------------------------------------
!
!        ! Raindrop detachment
!        call rainDropDetachCell()
!        !do i=1,NA
!            !write(*,'(13(F15.2))')  D_R(i,:)
!        !    print*, mv(i)%precip, ' ', D_R(i)
!        !end do
!
!
!        ! Overland flow detachment
!        call overlandFlowDetachCell()
!        !do i=1,NA
!            !write(*,'(13(F15.2))')  D_F(i,:)
!        !    print*, ofh(i)%depth, ' ',ofh(i)%slope, ' ', D_F(i)
!        !end do
!
!    !    call overlandFlowTransCapaCell()
!    !    do i=1,NA
!    !        !write(*,'(13(F15.2))')  G(i,:)
!    !        print*, G(i)
!    !    end do
!
!        ! Calculate variables at the cell edges
!         call varsAtCellEdge()
!
!        ! Hillslope routing
!        !- 1) Estimate transport capacity at the outlet of the cell
!        print*,
!        print*, '----1----'
!        call overlandFlowTransCapa_outCell()
!        !do i=1,NA
!            !write(*,'(13(F15.2))')  G(i,:)
!        !    print*, i, G_out(i)%east, G_out(i)%north, G_out(i)%west, G_out(i)%south
!        !end do
!        DELTac = 0
!        do while (DELTac < DELTout)
!            ! Initilize variables
!            call sed_before_vars()
!
!
!
!            !- 2) Estimate potential sediment concentration
!            print*,
!            print*, '----2----'
!            call potentialSedConc()
!            do i=1,NA
!                !write(*,'(13(F15.2))')  G(i,:)
!                !print*, i, F(i)%nfluxes, F(i)%east, F(i)%north, F(i)%west, F(i)%south !C(i)
!                print*, csa(i)%C
!            end do
!
!            !- 3) Estimate potential change in surface elevation
!            print*,
!            print*, '----3----'
!            call potentialChanSufEle()
!            do i=1,NA
!                print*, csa(i)%pot_Dz
!            end do
!
!            !- 4) Estimate available change in surface elevation
!            print*,
!            print*, '----4----'
!            call availableChanSufEle()
!            do i=1,NA
!                print*, csa(i)%ava_Dz
!            end do
!
!            !- 5) Compare pot_Dz vs. ava_Dz
!            print*,
!            print*, '----5----'
!            call cellConAndChanSufEle()
!            do i=1,NA
!                print*, csa(i)%Dz, csa(i)%SD
!            end do
!
!            DELTac = DELTac + DELT
!
!        end do
!
!
!        !-----------------------------------------------------
!        ! IN-STREAM TRANSPORT
!        !-----------------------------------------------------
!        ! Instream routing
!        ! - Longitudinal sediment velocity
!        call longSedVelocity_grid
!
!        ! - Sediment inflow to channel
!        ! -- Channel bank erosion
!        call grid_bankErosion
!        ! -- Sediment inflows
!        call flowSedInput
!
!        ! - In-stream sediment routing
!        call inStreamRouting
!        do i=1,NA
!            print*, isrr(i)%G
!        end do
!        ! - Update frac diameters in active and parent layers
!        call setFracDiame
!
!        print*,'-------'
!        print*,'Time iter: ', titer
!        print*,'-------'
!
!        if (titer>2) stop
!        titer = titer + 1
!        !stop
!
!    end do
!
!    ! Closing files containing outfields and rbm_input
!    call closeFile(unitOUTFIELDprecip)
!    call closeFile(unitOUTFIELDevap)
!    call closeFile(unitOUTFIELDwr_runovf)
!    call closeFile(unit_rbm_input)


        !print*, '>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        !do i=1,NA

        !    print*, csaB%C(i)
        !end do

!    call hydroAtCellEdges()
!    do i=1,NA
!        !write(*,'(13(F10.2))')   waterDepth_edge(i,:)%south
!        print*,
!    end do

    !do i = 1, ngru
    !    print*, rainSoilEr(i)%soilVul, rainSoilEr(i)%damEffecSw
    !end do



end program mesh_sed
