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
    use sed_massBalance

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

    !> Read 'MESH_sed_reservoir.ini'
    call read_sed_reservoir





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

    !> Update time step
    call currDate_update
    write (dateIn, 10) currDate%year, &
        currDate%month,  currDate%day,  currDate%hour, currDate%mins, currDate%secs
    print *, 'DATE & TIME: ', dateIn

    !> Open and read header of MET DATA input files
    !call getMat1_MESH_OUTFIELD(filepathOUTFIELD, filenamePR, unitPR, mat1PR, date1PR, mat2PR, date2PR, dateIn, iosPR)
    !print *, date1PR, date2PR
    !call readHeader_MESH_OUTFIELD(filepathOUTFIELD, filenameEV, unitEV, mat1EV, date1EV, dateIn, iosEV)
    !call getMat1_MESH_OUTFIELD(filepathOUTFIELD, filenameEV, unitEV, mat1EV, date1EV, mat2EV, date2EV, dateIn, iosEV)
    !print *, date1EV, date2EV
    !stop

    !> Open and read header of OVERLAND FLOW DATA input file
    !call getMat1_MESH_OUTFIELD(filepathOUTFIELD, filenameWR, unitWR, mat1WR, date1WR, mat2WR, date2WR, dateIn, iosWR)

    !> Open and read header of IN-STREAM FLOW DATA
!    call getMat1_rbm_input(filepathRBM, filenameRBM, unitRBM, mat1RBM1, date1RBM1, mat1RBM2, date1RBM2, &
!                                     mat1RBM3, date1RBM3, mat1RBM4, date1RBM4, mat1RBM5, date1RBM5, mat1RBM6, date1RBM6, &
!                                     mat2RBM1, date2RBM1, mat2RBM2, date2RBM2, mat2RBM3, date2RBM3, mat2RBM4, date2RBM4, &
!                                     mat2RBM5, date2RBM5, mat2RBM6, date2RBM6, dateIn, iosRBM)

    call getMat1_sed_input(sedi%filepathSED, sedi%filenameSED, sedi%unitSED, mat1SED1, date1SED1, mat1SED2, date1SED2, &
                                     mat1SED3, date1SED3, mat1SED4, date1SED4, mat1SED5, date1SED5, mat1SED6, date1SED6, &
                                     mat1SED7, date1SED7, mat1SED8, date1SED8, mat1SED9, date1SED9, mat1SED10, date1SED10, &
                                     mat1SED11, date1SED11, &
                                     mat2SED1, date2SED1, mat2SED2, date2SED2, mat2SED3, date2SED3, mat2SED4, date2SED4, &
                                     mat2SED5, date2SED5, mat2SED6, date2SED6, mat2SED7, date2SED7, mat2SED8, date2SED8, &
                                     mat2SED9, date2SED9, mat2SED10, date2SED10, mat2SED11, date2SED11, &
                                     dateIn, sedi%iosSED)

!    print *,   trim(date1RBM1), ' ', trim(date2RBM1)
!    print *,   trim(date1RBM2), ' ', trim(date2RBM2)
!    print *,   trim(date1RBM3), ' ', trim(date2RBM3)
!    print *,   trim(date1RBM4), ' ', trim(date2RBM4)
!    print *,   trim(date1RBM5), ' ', trim(date2RBM5)
!    print *,   trim(date1RBM6), ' ', trim(date2RBM6)


    !> Open the files to write time series at different grid points
    call openFilesForGrPointTS_2()

    !> Open the files to write time series of fields (griddded data)
    !call openFilesForFieldTS()

    ! Save initial variables  for next time-step calculation
    !call sed_before_vars()

    ! MAIN LOOP: ITERATION THROUGH TIME
    !10  format('"',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',I2.2,':',I2.2,'.000"')
10  format(I4,'/',I2.2,'/',I2.2,1X,I2.2,':',I2.2,':',I2.2,'.000')
    do

        !call currDate_update
        !print *, currDate%year, ' ', currDate%month,' ',  currDate%day, ' ',currDate%jday, &
        !        ' ', currDate%hour, ' ',currDate%mins, ' ', currDate%secs
        !write (dateIn, '( '"',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',I2.2,':',I2.2,'.000"' )') currDate%year, &
        !write (dateIn, 10) currDate%year, &
        !currDate%month,  currDate%day,  currDate%hour, currDate%mins, currDate%secs
        !print *, 'DATE & TIME: ', dateIn

        !dateIn = '"'//"2002/01/30      1:15:20.000"//'"'
        !call parseDateTime(dateIn, currDate)
        !stop

        !-----------------------------------------------------
        ! READING INPUT DATA
        !-----------------------------------------------------
        !> READ sed_input.r2c (MESH OUTPUT FILE)
        call readMESHoutputData()

        !-----------------------------------------------------
        ! HILLSLOPE LOADS
        !-----------------------------------------------------
        call hillSlopeLoad_grid()

        !-----------------------------------------------------
        ! RIVER BANK EROSION LOADS
        !-----------------------------------------------------
        call bankEroLoad_grid()

        !-----------------------------------------------------
        ! RIVER BED EROSION(RESUSPENSION) LOADS
        !-----------------------------------------------------
        call streamBottomResLoad_grid()

        !-----------------------------------------------------
        ! MASS BALANCE AT EACH INDIVIDUAL CELL TO ESTIMATE NEW
        ! CONCENTRATIONS
        !-----------------------------------------------------
        call massBalance2()
        !print *, 'L_bank: ', sum(cmb(848)%L_bank(:)), 'L_hill: ', sum(cmb(848)%L_hill(:)), 'L_in: ',&
        !            sum(cmb(848)%L_in(:)), 'L_res: ', sum(cmb(848)%L_res(:)), rh(848)%depth, &
        !            sum(cmb(848)%C(:)*cbsca(848)%frac(:)), &
        !            sum(cmb(848)%C_pot(:)*cbsca(848)%frac(:)), rh(848)%discharge

		if (isnan(sum(cmb(848)%L_bank(:)))) stop '"L_bank" is a NAN'
		if (isnan(sum(cmb(848)%L_hill(:)))) stop '"L_hill" is a NAN'
		if (isnan(sum(cmb(848)%L_in(:)))) stop '"L_in" is a NAN'
		if (isnan(sum(cmb(848)%L_res(:)))) stop '"L_res" is a NAN'
		if (isnan(sum(cmb(848)%C_pot(:)*cbsca(848)%frac(:)))) stop '"C_pot" is a NAN'
        !print *, 'L_bank: ', sum(cmb(1000)%L_bank(:)), 'L_hill: ', sum(cmb(1000)%L_hill(:)), 'L_in: ',&
        !            sum(cmb(1000)%L_in(:)), 'L_res: ', sum(cmb(1000)%L_res(:)), rh(1000)%depth, &
        !            sum(cmb(1000)%C(:)*cbsca(1000)%frac(:)), &
        !            sum(cmb(1000)%C_pot(:)*cbsca(1000)%frac(:)), rh(1000)%discharge

        !print *, sum(cmb(999)%C(:)*cbsca(999)%frac(:)), rh(999)%velocity, sum(cmb(999)%L_out(:)), &
        !        sum(cmb(1000)%C(:)*cbsca(1000)%frac(:)), rh(1000)%velocity, sum(cmb(1000)%L_out(:))
        !print *, sum(cbsca(1325)%frac(:)), sum(bsca(1325)%frac(:))
        !stop


        !-----------------------------------------------------
        ! WRITE OUT MODEL PRONOSTIC VARIABLES
        !-----------------------------------------------------
        !call writeGrPointTS_2
        call writeHourAveEstimate
        !call writeFieldTS

        ! Save variables at the current time step for next time-step calculation
        !call sed_before_vars()


        if (stopExec()) exit

        !-----------------------------------------------------
        ! UPDATE AND PRINT DATE AND TIME
        !-----------------------------------------------------
        call currDate_update
        write (dateIn, 10) currDate%year, &
        currDate%month,  currDate%day,  currDate%hour, currDate%mins, currDate%secs
        print *, 'DATE & TIME: ', dateIn

    end do

    !> Clossing MET DATA files
    !call close_MESH_OUTFIELD(unitPR)
    !call close_MESH_OUTFIELD(unitEV)

    !> Clossing OVERLAND FLOW DATA file
    !call close_MESH_OUTFIELD(unitWR)

    !> Clossing IN-STREAM FLOW DATA file
    !call close_MESH_OUTFIELD(unitRBM)

    !> Clossing sed_input.r2c file
    call close_MESH_OUTFIELD(sedi%unitSED)


    !> Closing the files to write time series at different grid points
    call closeFilesForGrPointTS()

    !> Closing the files to write time series of fields (griddded data)
    !call closeFilesForFieldTS()





end program mesh_sed

