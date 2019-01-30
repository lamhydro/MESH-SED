!------------------------------------------------------------------------------
! Sediment transport in cold region catchments: the MESH-SED model
!------------------------------------------------------------------------------
!
!> @brief
!> MODULE: sed_config
!>
!> @detail This module contains subroutine to allocate variables, set-up soil
!> and vegetation grid-cell characteristics, save state variables for the early
!> time step, estimate characteristics sediment particle diameters, set-up
!> model initial conditions, set-up the grid cell connectivity, control time step
!> and model execution termination.
!>
!> @author Luis Morales (LAM), GIWS & GWF.
!> - July, 2017
!> @date January, 2019-LAM
!> - Documenting the code
!> @todo
!---------------------------------------------------------------------------

module sed_config
    use sed_vars

    implicit none

    contains

        ! Allocating variables
        !---------------------------------------------------------------------------
        !> @brief
        !> Allocate variables
        !>
        !> @detail Allocate all static and dynamic variables that change in time.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine sed_allocate_var

            implicit none

            ! Allocate input variables
            !allocate(gca%iniSoilDepth(NA), gca%sedFluxBondC(NA), gca%cellCanopyCov(NA),&
            !        gca%cellGroundCov(NA), gca%soilType(NA), gca%cellVege(NA))

            !allocate(cba%iniThickBed(NA), cba%soilTypeBed(NA), cba%soilTypeBank(NA))

            allocate(imperCellAre(NA), ipos(NA), jpos(NA), dummy(yCount, xCount), mat1PR(yCount, xCount), &
                    mat2PR(yCount, xCount), mat1EV(yCount, xCount), mat2EV(yCount, xCount), &
                    mat1WR(yCount, xCount), mat2WR(yCount, xCount), variMat(yCount, xCount), &
                    dummyC(yCount,xCount), dummyI(yCount,xCount), grs(NA), DISC(yCount,xCount), INFL(yCount,xCount), &
                    QDIF(yCount,xCount), DEPT(yCount,xCount), WIDT(yCount,xCount), VELO(yCount,xCount) )

            allocate(sedi%DISC(yCount,xCount), sedi%DEPT(yCount,xCount), sedi%WIDT(yCount,xCount), &
                    sedi%CHLE(yCount,xCount),  sedi%CHSL(yCount,xCount), sedi%VELO(yCount,xCount), &
                    sedi%PREP(yCount,xCount), sedi%EVAP(yCount,xCount), sedi%OFDE(yCount,xCount), &
                    sedi%LASL(yCount,xCount), sedi%CEWI(yCount,xCount))


!            allocate(mat1RBM1(yCount, xCount),mat2RBM1(yCount, xCount), mat1RBM2(yCount, xCount), &
!                     mat2RBM2(yCount, xCount), mat1RBM3(yCount, xCount), mat2RBM3(yCount, xCount), &
!                     mat1RBM4(yCount, xCount), mat2RBM4(yCount, xCount), mat1RBM5(yCount, xCount), &
!                     mat2RBM5(yCount, xCount), mat1RBM6(yCount, xCount), mat2RBM6(yCount, xCount) )

            allocate(mat1SED1(yCount, xCount),mat2SED1(yCount, xCount), mat1SED2(yCount, xCount), &
                     mat2SED2(yCount, xCount), mat1SED3(yCount, xCount), mat2SED3(yCount, xCount), &
                     mat1SED4(yCount, xCount), mat2SED4(yCount, xCount), mat1SED5(yCount, xCount), &
                     mat2SED5(yCount, xCount), mat1SED6(yCount, xCount), mat2SED6(yCount, xCount), &
                     mat1SED7(yCount, xCount), mat2SED7(yCount, xCount), mat1SED8(yCount, xCount), &
                     mat2SED8(yCount, xCount), mat1SED9(yCount, xCount), mat2SED9(yCount, xCount), &
                     mat1SED10(yCount, xCount), mat2SED10(yCount, xCount), mat1SED11(yCount, xCount), &
                     mat2SED11(yCount, xCount) )


            allocate(mv(NA), ofh(NA), ofhB(NA), rh(NA))

            allocate(gca(NA), cba(NA))

            allocate(sca(NA), vca(NA), cbsca(NA), cbscaB(NA), bsca(NA))


            ! Allocate variables
            !allocate(cellSoilChar(NA), cellVegeChar(NA))
            allocate(D_R(NA), D_F(NA))
            allocate(cn(NA), G_out(NA), F(NA), waterDepth_edge(NA),&
                    waterSslope_edge(NA), flowVelocity_edge(NA), discharge_edge(NA))
            allocate(FB(NA), waterDepthB_edge(NA),&
                    waterSslopeB_edge(NA), flowVelocityB_edge(NA), dischargeB_edge(NA))

            allocate(csa(NA), csaB(NA))
            !allocate(csaB%C(NA), csaB%pot_Dz(NA), csaB%ava_Dz(NA))

            allocate(ion(NA), isrr(NA), isrrB(NA))
            allocate(isrn(NA))
            !allocate(bke)

        end subroutine sed_allocate_var

        !---------------------------------------------------------------------------
        !> @brief
        !> Set up soil and vegetation attributes to each cell.
        !>
        !> @detail Assign to each active grid cell soil and vegetation attributes.
        !> Characteristics are assigned for overland soil, channel bank soil and
        !> channel bed soil.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine setSoilAndVegeParamInCell()
        ! Save in each grid-cell all the parameters needed to estimate
        ! rain-drop detachment (D_R)

            use sed_tools
            implicit none

            !allocate(cellSoilChar(yCount,xCount), cellVegeChar(yCount,xCount))
            !allocate(cellSoilChar(NA), cellVegeChar(NA))

            do i = 1, NA
            !do i = 1,yCount
            !    do j = 1,xCount
                    !if (cell(i,j) == 1) then

                        ! For soil characteristics
!                        do k = 1, nsoilType
!                            if (trim(cellSoil(i)) == trim(soilChar(k)%name)) then
!                                cellSoilChar(i)%name = soilChar(k)%name
!                                cellSoilChar(i)%diameter = soilChar(k)%diameter
!                                cellSoilChar(i)%density = soilChar(k)%density
!                                cellSoilChar(i)%surfPorosity = soilChar(k)%surfPorosity
!                                cellSoilChar(i)%soilDetach = soilChar(k)%soilDetach
!                                cellSoilChar(i)%overlandDetach = soilChar(k)%overlandDetach
!                            end if
!                        end do

                        ! For soil characteristics
                        sca(i)%frac = sa(gca(i)%soilType)%frac
                        sca(i)%diameter = fracDiame(sca(i)%frac, sp%meanD, 0.50) !sum(sca(i)%frac*sp%meanD)
                        !print*, 'herere', sca(i)%diameter
                        sca(i)%density  = sa(gca(i)%soilType)%density
                        !print *, 'den ', i, gca(i)%soilType, sca(i)%density
                        sca(i)%porosity = sa(gca(i)%soilType)%porosity
                        sca(i)%soilDetach       = sa(gca(i)%soilType)%soilDetach
                        !sca(i)%overlandDetach   = sa(gca(i)%soilType)%overlandDetach

                        ! For river bank soil characteristics
                        bsca(i)%frac = sa(cba(i)%soilTypeBank)%frac
                        bsca(i)%diameter = fracDiame(bsca(i)%frac, sp%meanD, 0.50) !sum(sca(i)%frac*sp%meanD)
                        bsca(i)%density   = sa(cba(i)%soilTypeBank)%density
                        bsca(i)%chanBankDetach   = sa(cba(i)%soilTypeBank)%chanBankDetach


                        ! For vegetation characteristics
!                        do k = 1, nvegeType
!                            if (trim(cellVege(i)) == trim(vegeChar(k)%name)) then
!                                cellVegeChar(i)%name = vegeChar(k)%name
!                                cellVegeChar(i)%dropDiam = vegeChar(k)%dropDiam
!                                cellVegeChar(i)%percDrip = vegeChar(k)%percDrip
!                                cellVegeChar(i)%fallHeight = vegeChar(k)%fallHeight
!                            end if
!                        end do
!                        print*, 'here', cellVegeChar(i)%name

                        ! For vegetation characteristics
                        do k = 1, NVEGE
                            if (trim(gca(i)%cellVege) == trim(va(k)%name)) then
                                vca(i)%name = va(k)%name
                                vca(i)%dropDiam = va(k)%dropDiam
                                vca(i)%percDrip = va(k)%percDrip
                                vca(i)%fallHeight = va(k)%fallHeight
                                exit
                            end if
                        end do
                        !print*, 'here2', vca(i)%name

                        ! For channel-bed soil characteristics
                        cbsca(i)%density = sa(cba(i)%soilTypeBed)%density
                        cbsca(i)%porosity = sa(cba(i)%soilTypeBed)%porosity


                    !end if
                !end do
            !end do
            end do
        end subroutine setSoilAndVegeParamInCell

        !---------------------------------------------------------------------------
        !> @brief
        !> Save the value of state variables for a prior time step.
        !>
        !> @detail In finite difference methods, the estate variables are calculated
        !> \f$(X^{t+1})\f$ sometimes based on the previous state \f$(X^{t-1})\f$ and
        !> the current state \f$(X^{t-1})\f$.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine sed_before_vars()
            implicit none

            ! Hydro variables at the centre
            !waterDepthB = waterDepth
            !waterSslopeB = waterSslope
            !flowVelocityB = flowVelocity
            !dischargeB = discharge

            ! Hydro variables at the centre
            ofhB = ofh
            !ofhB%depth = ofh%depth
            !ofhB%slope = ofh%slope
            !ofhB%velocity = ofh%velocity
            !ofhB%discharge = ofh%discharge


            ! Hydro variables at the edge
            waterDepthB_edge = waterDepth_edge
            waterSslopeB_edge = waterDepth_edge
            flowVelocityB_edge = flowVelocity_edge
            dischargeB_edge = discharge_edge

            ! Sed variables at the centre
            csaB = csa
            !csaB%C = csa%C
            !csaB%pot_Dz = csa%pot_Dz
            !csaB%ava_Dz = csa%ava_Dz
            !csaB%SD = csa%SD


            ! Sed variables at the edge
            FB = F

            ! Instream sediment routing variables
            isrrB = isrr
            cbscaB = cbsca

        end subroutine sed_before_vars

        !---------------------------------------------------------------------------
        !> @brief
        !> Estimate some characteristic diameters.
        !>
        !> @detail For each soil type at each active grid-cell, this subroutine
        !> estimate \f$ D_{16}\f$, \f$ D_{50}\f$, \f$ D_{84}\f$ and \f$ D_{99}\f$
        !> from the particle size density distribution.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine setFracDiame()
            use sed_tools
            implicit none


            do i = 1, NA
                do j = 1, nchbedly
                    ! These diam. need to be calculated from the distribution
                    cbsca(i)%ly(j)%D16 = fracDiame(cbsca(i)%ly(j)%frac, sp%meanD, 0.16)
                    cbsca(i)%ly(j)%D50 = fracDiame(cbsca(i)%ly(j)%frac, sp%meanD, 0.50)
                    cbsca(i)%ly(j)%D84 = fracDiame(cbsca(i)%ly(j)%frac, sp%meanD, 0.84)
                    cbsca(i)%ly(j)%D99 = fracDiame(cbsca(i)%ly(j)%frac, sp%meanD, 0.99)
                end do
            end do

        end subroutine setFracDiame



        !---------------------------------------------------------------------------
        !> @brief
        !> Set up initial conditions.
        !>
        !> @detail Set up initial values such as the path to the sediment transport
        !> project, the startDate, the initial values of transport capacity, the
        !> the initial sediment concentration, the initial depth of loose soil,
        !> the initial thickness of active and parent bed layers and their sediment
        !> class composition.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine sed_initial_conditions()

            implicit none

            !> Set directory paths
            !filepathOUTFIELD = trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/"
            !> filepathRBM = trim(ADJUSTL(MESHdir))
            sedi%filepathSED = trim(ADJUSTL(MESHdir))

            !> Set current date as start date
            currDate = startDate

            ! Set inital transport capacity
            F(:)%east = 0.; F(:)%north = 0.; F(:)%west = 0.; F(:)%south = 0.
            do i = 1, NA

                ! Set initial concentration
                csa(i)%C(:) = 0.

                ! Set initial depth of loose soil (m)
                csa(i)%SD = gca(i)%iniSoilDepth

                ! Set initial thickness and particle fractions of active and parent channel bed layers
                ! This values must be D99 or 10 mm for sandy river beds.
                do j = 1, nchbedly
                    cbsca(i)%ly(j)%frac = sa(cba(i)%soilTypeBed)%frac
                    cbsca(i)%ly(j)%thick = cba(i)%iniThickBed/nchbedly
                end do
                cbsca(i)%thick = sum(cbsca(i)%ly(:)%thick)

!                cbsca(i)%frac_aly = cbsa(cba(i)%soilTypeBed)%frac   ! Proportion of sed. in the active bed layer
!                cbsca(i)%frac_ply = cbsa(cba(i)%soilTypeBed)%frac   ! Proportion of sed. in the parent bed layer
!                cbsca(i)%thick_aly = cba(i)%iniThickBed             ! Active bed layer thickness
!                cbsca(i)%thick_ply = cba(i)%iniThickBed             ! Parent bed layer thickness
!                cbsca(i)%thick = cbsca(i)%thick_aly + cbsca(i)%thick_ply ! Total layer thickness
!                cbsca(i)%D99_ply = 0.01 ! Sed. diam. for which 99% of the orig. parent lay. sed. part. are finer.


                ! Initial thickness of channel bed material (m)
                !do k = 1, nsedpar
                !    isrn(i)%SD(k) = cbsca(i)%frac(k)*cba(i)%iniThickBed
                !end do

            end do

            ! Set fraction diameters
            call setFracDiame()

        end subroutine sed_initial_conditions


        !---------------------------------------------------------------------------
        !> @brief
        !> Set up some variables related with the model configuration.
        !>
        !> @detail Set up the position (i,j) of active cells according to rank and
        !> the cell neighbours to each active cell
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine sed_config_init

            implicit none

            ! Allocating memory for variables
            !allocate(soilVul(ngru), damEffecSw(ngru))
            !allocate(rainSoilEr(ngru))

            ! Set position (i,j) of active cells according to rank
            do k = 1,NA
                do i = 1,yCount
                    do j = 1, xCount
                        if (rank(i,j) == k) then
                            ipos(k) = i
                            jpos(k) = j
                            !nextV(k) = next(i,j)
                        end if
                    end do
                end do
            end do


            ! Get the neighbor cells to each cell
            cn(:)%east = 0; cn(:)%south = 0; cn(:)%west = 0; cn(:)%north = 0 ! Inizialization
            do k = 1, NA
                ! East
                if (jpos(k)+1 <= xCount) then
                    if (rank(ipos(k),jpos(k)+1) /= 0.) then
                        cn(k)%east = rank(ipos(k),jpos(k)+1)
                    end if
                end if
                ! North
                if (ipos(k)+1 <= yCount) then
                    if (rank(ipos(k)+1,jpos(k)) /= 0.) then
                        cn(k)%south = rank(ipos(k)+1,jpos(k))
                    end if
                end if
                ! West
                if (jpos(k)-1 >= 1) then
                    if (rank(ipos(k),jpos(k)-1) /= 0.) then
                        cn(k)%west = rank(ipos(k),jpos(k)-1)
                    end if
                end if
                ! South
                if (ipos(k)-1 >= 1) then
                    if (rank(ipos(k)-1,jpos(k)) /= 0.) then
                        cn(k)%north = rank(ipos(k)-1,jpos(k))
                    end if
                end if

            end do

        end subroutine sed_config_init

        !---------------------------------------------------------------------------
        !> @brief
        !> Get the in and out reaches for each active cell.
        !>
        !> @detail For each active cell, get the upstream reaches that flow into the
        !> cell and the downstream receiver cell.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine setInAndOutStreamToNode()
        !> For eachs node (cell) find the elements (reaches) that flow into and the element(reach)
        !> that flow out

            implicit none

            !> Get the number of inreaches for each active cell
            do k = 1, NA
                !ion(k)%reachOut = next(ipos(k),jpos(k))

                i = 0
                do j = 1,NA
                    if (next(ipos(j),jpos(j)) == k) i = i + 1
                end do
                !print*, i, allocated(ion(k)%reachIn)
                if (i > 0) allocate(ion(k)%reachIn(i))
                ion(k)%NreachIn = i

            end do

            !> Set the inreaches and the outreach for each active cell
            do k = 1, NA
                ion(k)%reachOut = next(ipos(k),jpos(k))

                i = 1
                do j = 1,NA
                    if (next(ipos(j),jpos(j)) == k) then
                        ion(k)%reachIn(i) = rank(ipos(j),jpos(j))
                        i = i + 1
                    end if
                end do

            end do

        end subroutine setInAndOutStreamToNode


        !---------------------------------------------------------------------------
        !> @brief
        !> Set the x and y dimension of each grid cell.
        !>
        !> @detail Set the x and y dimension of each grid cell assuming square
        !> cells.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo this need to be checked as grid cells are irregular and not squares
        !---------------------------------------------------------------------------
        subroutine setGridCellSize
            implicit none
            do i = 1, NA
                grs(i)%DELX = sqrt(shd%AREA(ipos(i),jpos(i)))
                grs(i)%DELY = sqrt(shd%AREA(ipos(i),jpos(i)))
            end do
        end subroutine setGridCellSize

        !---------------------------------------------------------------------------
        !> @brief
        !> Set up the time step
        !>
        !> @detail Set up time step at each iteration so, year, month, day, jday,
        !> hours, minutes and seconds are updated here.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !> @todo
        !---------------------------------------------------------------------------
        subroutine currDate_update()
            use sed_tools
            implicit none

            integer :: nti, i1

            nti = DELT/60

            if (nti < 1) nti = 1

            !> Increment the current time-step.
            do i1 = 1,nti
                !currDate%secs = currDate%secs + DELT! increment the current time step
                currDate%secs = currDate%secs + DELT/nti! increment the current time step
                if (currDate%secs == 60) then
                    currDate%secs = 0
                    currDate%mins = currDate%mins + 1
                    if (currDate%mins == 60) then
                        currDate%mins = 0
                        currDate%hour = currDate%hour + 1
                        if (currDate%hour == 24) then
                            currDate%hour = 0
                            currDate%jday = currDate%jday + 1
                            if (currDate%jday >= 366) then
                                if (mod(currDate%year, 400) == 0) then !LEAP YEAR
                                    if (currDate%jday == 367) then
                                        currDate%jday = 1
                                        currDate%year = currDate%year + 1
                                    end if
                                else if (mod(currDate%year, 100) == 0) then !NOT A LEAP YEAR
                                    currDate%jday = 1
                                    currDate%year = currDate%year + 1
                                else if (mod(currDate%year, 4) == 0) then !LEAP YEAR
                                    if (currDate%jday == 367) then
                                        currDate%jday = 1
                                        currDate%year = currDate%year + 1
                                    end if
                                else !NOT A LEAP YEAR
                                    currDate%jday = 1
                                    currDate%year = currDate%year + 1
                                end if
                            end if
                        end if
                    end if
                end if
             end do

!            !> Year:
!            if (old_year /= currDate%year) then
!                ic%count_year = ic%count_year + 1
!            end if
!
!            !> Julian day:
!            if (old_jday /= currDate%jday) then
!                ic%count_jday = ic%count_jday + 1
!                ic%ts_daily = 0
!            end if

            !> Determine the current month and day.
            call julian2monthday(currDate%jday, currDate%year, currDate%month, currDate%day)

!            !> Month:
!            if (old_month /= currDate%month) then
!                ic%count_month = ic%count_month + 1
!            end if
!
!            !> Hourly:
!            if (old_hour /= currDate%hour) then
!                ic%count_hour = ic%count_hour + 1
!                ic%ts_hourly = 0
!            end if
!
!            !> Minutes:
!            if (old_mins /= currDate%mins) then
!                ic%count_mins = ic%count_mins + 1
!                if (currDate%mins == 0 .or. currDate%mins == 30) then
!                    ic%ts_halfhourly = 0
!                end if
!            end if

    !debug: Print the now.
    !print *, "now: Y JD M D HR"
    !print *, currDate_year, currDate_jday, currDate_month, currDate_day, currDate_hour

    !debug: Print count.
    !print *, "count: Y JD M D HR"
    !print *, ic%count_year, ic%count_jday, ic%count_month, ic%count_day, ic%count_hour

            !> Update time-step counters.
!            ic%ts_daily = ic%ts_daily + 1
!            ic%ts_hourly = ic%ts_hourly + 1
!            ic%ts_halfhourly = ic%ts_halfhourly + 1
!            ic%ts_count = ic%ts_count + 1

        end subroutine currDate_update

        !---------------------------------------------------------------------------
        !> @brief
        !> Set up the time step
        !>
        !> @detail Set up time step at each iteration so, year, month, day, jday,
        !> hours, minutes and seconds are updated here.
        !>
        !> @author Luis Morales (LAM), GIWS & GWF.
        !> - July, 2017
        !> @date January, 2019-LAM
        !> - Documenting the code
        !>
        !> @param[in] currDate% Current date (year julian-day hour minutes seconds)
        !> @param[in] stopDate% Date (year julian-day hour minutes seconds) at which
        !> the program stop
        !> @return stopExec Logical value (.true. or .false.) to stop program
        !> execution.
        !> @todo
        !---------------------------------------------------------------------------
        logical function stopExec()
            implicit none

!            yearCheck = currDate%year >= stopDate%year
!            jdayCheck = currDate%jday >= stopDate%jday
!            hourCheck = currDate%hour >= stopDate%hour
!            minsCheck = currDate%mins >= stopDate%mins
!            secsCheck = currDate%secs >= stopDate%secs

           if((currDate%year >= stopDate%year) .and. &
            (currDate%jday >= stopDate%jday) .and. &
            (currDate%hour >= stopDate%hour) .and. &
            (currDate%mins >= stopDate%mins) .and. &
            (currDate%secs >= stopDate%secs)) then

                stopExec = .true.
            else
                stopExec = .false.
            end if
        end function stopExec



end module sed_config
