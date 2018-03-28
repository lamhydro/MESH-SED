!> Define new variables for the sediment transport model
module sed_input
    use sed_vars

    implicit none

    contains

	!subroutine initSedConditions()

	!end subroutine initSedConditions


        !> Read the name of the input directory of the study case.
        !> and set the study case directory.
        subroutine mesh_sed_case

            implicit none

            call getarg(1,casename)

            if(len_trim(casename) == 0) then
                print *, "Please provide study case name on command line"
                print *, "Stopping..."
                stop
            end if

            ! Name of the study case directory
            casedir = './' // trim(ADJUSTL(casename)) // '/'

        end subroutine mesh_sed_case


!        subroutine read_sed_options
!            implicit none
!
!            ! Open the file
!            open(unit = unitParam,file=trim(ADJUSTL(casedir)) // "MESH_sed_options.ini",status="old",action="read",iostat=ios)
!            if (ios /= 0) then
!                close(unit = unitParam)
!                print *, "The file ", "MESH_sed_options.ini", " could not be opened."
!                print *, "Stopping..."
!                stop
!
!            else
!                print *, 'READING: MESH_sed_options.ini '
!            end if
!
!            read(unit = unitParam,fmt = *)
!            read(unit = unitParam,fmt = *) MESHdir
!            read(unit = unitParam,fmt = *)
!            read(unit = unitParam,fmt = *) OUTFIELDfolder
!
!            close(unit = unitParam)
!
!        end subroutine read_sed_options


        !> Read "MESH_sed_param.ini" file. Estimate the mean diameter of
        !> each sediment class
        subroutine read_sed_param
            use sed_tools
            implicit none


            !> Set sediment particle diameters
            ! Representative sediment and soil particle diameters (mm) (Wentworth (1922)). Taken from UCL geog.
            sp%name = (/ 'clay          ', 'veryFineSilt  ', 'fineSilt      ', 'mediumSilt    ', &
                         'coarseSilt    ', 'veryFineSand  ', 'fineSand      ', 'mediumSand    ', &
                         'coarseSand    ', 'veryCoarseSand', 'granule       ', 'pebble        ', &
                         'cobble        ', 'boulder       ' /)
            sp%minD = (/ 0.00006, 0.0039, 0.0078, 0.0156, 0.0313, 0.0625, 0.125, 0.25, &
                          0.5, 1., 2., 4., 64., 256. /)
            sp%maxD = (/ 0.0039, 0.0078, 0.0156, 0.0313, 0.0625, 0.125, 0.25, &
                          0.5, 1., 2., 4., 64., 256., 4096. /)
            sp%meanD = (sp%minD + sp%maxD)*0.5

            ! Open the file
            open(unit = unitParam,file=trim(ADJUSTL(casedir)) // "MESH_sed_param.ini", &
                status="old",action="read",iostat=ios)
            if (ios /= 0) then
                close(unit = unitParam)
                print *, "The file ", "MESH_sed_param.ini", " could not be opened."
                print *, "Stopping..."
                stop

            else
                print *, 'READING: MESH_sed_param.ini '
            end if

!            read(unit = unitParam,fmt = *)
!            do i = 1, nsedpar
!                read(unit = unitParam, fmt = *) sp(i)%name, sp(i)%minD, sp(i)%maxD
!                sp(i)%meanD = (sp(i)%minD + sp(i)%maxD)*0.5
!                !print*, sp(i)%name, sp(i)%minD, sp(i)%maxD
!            end do
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) MESHdir
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) OUTFIELDfolder
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) theta, theta_r, phi
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) overlFlowCapaMethod
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) instreamFlowCapaMethod
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) DELT
            !> Simulation starting and stopping dates.
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = '(5i4)') startDate%year, startDate%jday, startDate%hour, startDate%mins, startDate%secs
            read(unit = unitParam,fmt = '(5i4)') stopDate%year, stopDate%jday, stopDate%hour, stopDate%mins, stopDate%secs
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) FPCRIT
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) opin%outputDir
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) opin%nGrPoint
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *)
            allocate(ogrp(opin%nGrPoint))
            !print *, opin%nGrPoint
            do i = 1, opin%nGrPoint
                !read(unit = unitParam, fmt = '(i8,a15,a15)') ogrp(i)%grNumb, ogrp(i)%varname, ogrp(i)%sedname
                !print *, i
                read(unit = unitParam, fmt = '(i10, a15, a15)') ogrp(i)%grNumb, ogrp(i)%varname, ogrp(i)%sedname
                ogrp(i)%idsedclass = getIdSedClass(ogrp(i)%sedname)
                !print *, ogrp(i)%grNumb, '  ',ogrp(i)%varname, '  ', ogrp(i)%sedname, ' ', ogrp(i)%idsedclass
            end do
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) ofin%outputDir
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *) ofin%nField
            read(unit = unitParam,fmt = *)
            read(unit = unitParam,fmt = *)
            allocate(ofie(ofin%nField))
            do i = 1, ofin%nField
                read(unit = unitParam, fmt = '(a15, a15)') ofie(i)%varname, ofie(i)%sedname
                ofie(i)%idsedclass = getIdSedClass(ofie(i)%sedname)
            end do



            !stop


            !> Determine the current month and day.
            call julian2monthday(startDate%jday, startDate%year, startDate%month, startDate%day)
            call julian2monthday(stopDate%jday, stopDate%year, stopDate%month, stopDate%day)

            !print *, startDate%year, startDate%jday, startDate%hour, startDate%mins, startDate%secs
            !print *, stopDate%year, stopDate%jday, stopDate%hour, stopDate%mins, stopDate%secs

            close(unit = unitParam)
            !stop
            !print*,  instreamFlowCapaMethod
            !print*, sp%meanD

        end subroutine read_sed_param

        !> Get the position of a a sed. particle class within the vector of classes.
        integer function getIdSedClass(sedname)
            implicit none
            character(len=15), intent(in) :: sedname
            integer :: i

            do i = 1, nsedpar
                if (trim(sp(i)%name) == trim(sedname)) then
                    getIdSedClass = i
                    exit
                end if
            end do
        end function getIdSedClass


        !> Read "MESH_sed_gridCell.ini" file.
        subroutine read_sed_gridCell
            implicit none


            !real, dimension(yCount, xCount) :: tempoR
            !integer, dimension(yCount, xCount) :: tempoI
            !character(len=8), dimension(yCount, xCount) :: tempoC

            ! Open the file
            open(unit = unitGrid,file=trim(ADJUSTL(casedir)) // "MESH_sed_gridCell.ini",status="old",action="read",iostat=ios)
            if (ios /= 0) then
                close(unit = unitGrid)
                print *, "The file ", "MESH_sed_gridCell.ini", " could not be opened."
                print *, "Stopping..."
                stop

            else
                print *, 'READING: MESH_sed_gridCell.ini '
            end if

            ! Reading grid cell attributes
            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummy(i,:) !tempoR(i,:)
            end do
            do i = 1, NA
                gca(i)%iniSoilDepth = dummy(ipos(i),jpos(i)) !tempoR(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummy(i,:) !tempoR(i,:)
            end do
            do i = 1, NA
                gca(i)%cellGroundCov = dummy(ipos(i),jpos(i)) !tempoR(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummyI(i,:) !tempoI(i,:)
            end do
            do i = 1, NA
                gca(i)%soilType = dummyI(ipos(i),jpos(i)) !tempoI(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummy(i,:) !tempoR(i,:)
            end do
            do i = 1, NA
                gca(i)%sedFluxBondC = dummy(ipos(i),jpos(i)) !tempoR(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummyC(i,:) !tempoC(i,:)
            end do
            do i = 1, NA
                gca(i)%cellVege = dummyC(ipos(i),jpos(i)) !tempoC(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummy(i,:) !tempoR(i,:)
            end do
            do i = 1, NA
                gca(i)%cellCanopyCov = dummy(ipos(i),jpos(i)) !tempoR(ipos(i),jpos(i))
            end do

            ! Reading channel bed attributes
            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummyI(i,:) !tempoI(i,:)
            end do
            do i = 1, NA
                cba(i)%soilTypeBed = dummyI(ipos(i),jpos(i)) !tempoI(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummyI(i,:) !tempoI(i,:)
            end do
            do i = 1, NA
                cba(i)%soilTypeBank = dummyI(ipos(i),jpos(i)) !tempoI(ipos(i),jpos(i))
            end do

            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummy(i,:) !tempoR(i,:)
            end do
            do i = 1, NA
                cba(i)%iniThickBed = dummy(ipos(i),jpos(i)) !tempoR(ipos(i),jpos(i))
            end do

            close(unit = unitGrid)

        end subroutine read_sed_gridCell


        !> Read "MESH_sed_soilVegChan.ini" file.
        subroutine read_sed_soilVegChan
            implicit none

            ! Open the file
            open(unit = unitSoVeCha,file=trim(ADJUSTL(casedir)) // "MESH_sed_soilVegChan.ini",status="old",action="read",iostat=ios)
            if (ios /= 0) then
                close(unit = unitSoVeCha)
                print *, "The file ", "MESH_sed_soilVegChan.ini", " could not be opened."
                print *, "Stopping..."
                stop

            else
                print *, 'READING: MESH_sed_soilVegChan.ini '
            end if

            ! Read soil attributes
            read(unit = unitSoVeCha,fmt = *)
            read(unit = unitSoVeCha, fmt = *) NSOIL
            allocate(sa(NSOIL))
            do j = 1,NSOIL
                read(unit = unitSoVeCha, fmt = *) (sa(j)%frac(i), i=1,nsedpar)
                read(unit = unitSoVeCha, fmt = *) sa(j)%density, sa(j)%porosity, &
                                                  sa(j)%soilDetach, sa(j)%overlandDetach, &
                                                  sa(j)%chanBankDetach
                !print*, sa(j)%frac
                !print*, sa(j)%density, sa(j)%porosity, &
                !                                  sa(j)%soilDetach, sa(j)%overlandDetach, &
                !                                  sa(j)%chanBankDetach
            end do

            ! Read vegetation atributes
            read(unit = unitSoVeCha, fmt = *)
            read(unit = unitSoVeCha, fmt = *) NVEGE
            allocate(va(NVEGE))
            read(unit = unitSoVeCha, fmt = *) (va(i)%name, i=1,NVEGE)
            read(unit = unitSoVeCha, fmt = *) (va(i)%fallHeight, i=1,NVEGE)
            read(unit = unitSoVeCha, fmt = *) (va(i)%dropDiam, i=1,NVEGE)
            read(unit = unitSoVeCha, fmt = *) (va(i)%percDrip, i=1,NVEGE)

            ! Read channel-bed soil attributes
            read(unit = unitSoVeCha,fmt = *)
            read(unit = unitSoVeCha,fmt = *) threChanCon, maxBedThick, ratioStress
!            read(unit = unitSoVeCha, fmt = *) NBEDSOIL
!            allocate(cbsa(NBEDSOIL))
!            do j = 1,NBEDSOIL
!                read(unit = unitSoVeCha, fmt = *) (cbsa(j)%frac(i), i=1,nsedpar)
!                read(unit = unitSoVeCha, fmt = *) cbsa(j)%density, cbsa(j)%porosity
!
!                !print*, cbsa(j)%frac
!                !print*, cbsa(j)%density, cbsa(j)%porosity
!            end do

            close(unit = unitSoVeCha)

        end subroutine read_sed_soilVegChan

        !> Read "MESH_drainage_database.r2c" file. This is a input file for
        !> MESH.
        subroutine read_MESH_drainage_database
            implicit none

            character(len=80) :: header
            integer :: i, j
            integer, dimension(30) :: attOrder = 0
            character(len=20), dimension(30) :: attName = ' '


            ! Open the file
            open(unit = unitDrainDB, &
                file=trim(ADJUSTL(MESHdir)) // "MESH_drainage_database.r2c", &
                status="old", &
                action="read", &
                iostat=ios)

            if (ios /= 0) then
                close(unit = unitDrainDB)
                print *, "The file ", "MESH_drainage_database.r2c", " could not be opened."
                print *, "Stopping..."
                stop

            else
                print *, 'READING: MESH_drainage_database.r2c '
            end if
            i = 1
            do while (header /= ':EndHeader')
                read(unitDrainDB,'(A80)') header

                if (header(1:19) == ':NominalGridSize_AL') then
                    read(header(20:),*) shd%AL
                    !print *, shd%AL
                end if

                if (header(1:16) == ':TotalNumOfGrids') then
                    read(header(17:),*) shd%NAA
                    !print *, shd%NAA
                end if

                if (header(1:16) == ':NumGridsInBasin') then
                    read(header(17:),*) shd%NA
                    !print *, shd%NA
                end if

                if (header(1:8) == ':xOrigin') then
                    read(header(9:),*) shd%xOrigin
                    !print *, shd%xOrigin
                end if
                if (header(1:8) == ':yOrigin') then
                    read(header(9:),*) shd%yOrigin
                    !print *, shd%yOrigin
                end if

                if (header(1:7) == ':xCount') then
                    read(header(8:),*) shd%xCount
                    !print *, shd%xCount
                end if
                if (header(1:7) == ':yCount') then
                    read(header(8:),*) shd%yCount
                    !print *, shd%yCount
                end if

                if (header(1:7) == ':xDelta') then
                    read(header(8:),*) shd%xDelta
                    !print *, shd%xDelta
                end if
                if (header(1:7) == ':yDelta') then
                    read(header(8:),*) shd%yDelta
                    !print *, shd%yDelta
                end if
                if (header(1:14) == ':AttributeName') then
                    read(header(16:17),*) attOrder(i)
                    !print *, attOrder(i)
                    read(header(18:len_trim(header)),*) attName(i)
                    !print *, attName(i)
                    i = i + 1
                end if

            end do


            ! Assigning values to similar variables
            NA = shd%NA
            yCount = shd%yCount
            xCount = shd%xCount
            !DELX = shd%AL
            !DELY = shd%AL

            allocate(shd%rank(yCount,xCount), shd%next(yCount,xCount), &
                shd%ELEV(yCount,xCount), shd%AREA(yCount,xCount), &
                shd%SLOPE_CHNL(yCount,xCount), shd%CHNL_LEN(yCount,xCount))

            do j = 1,maxval(attOrder)
                if (trim(attName(j)) == 'Rank') then
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *) shd%rank(i,:)
                    end do
                else if (trim(attName(j)) == 'Next') then
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *) shd%next(i,:)
                    end do
                else if (trim(attName(j)) == 'ChnlSlope') then
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *) shd%SLOPE_CHNL(i,:)
                    end do
                else if (trim(attName(j)) == 'Elev') then
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *) shd%ELEV(i,:)
                    end do
                else if (trim(attName(j)) == 'ChnlLength') then
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *) shd%CHNL_LEN(i,:)
                    end do
                else if (trim(attName(j)) == 'GridArea') then
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *) shd%AREA(i,:)
                    end do
                else
                    do i = 1, yCount
                        read(unit = unitDrainDB, fmt = *)
                    end do
                end if
            end do
            allocate(rank(yCount,xCount), next(yCount,xCount))
            rank(:,:) = nint(shd%rank(:,:))
            !read(unitDrainDB, *) ((shd%next(i, j), i = 1, yCount), j = 1, xCount)
            !do i = 1, yCount
            !    read(unit = unitDrainDB, fmt = *) shd%next(i,:)
            !end do
            next(:,:) = nint(shd%next(:,:))

            close(unit = unitDrainDB)

        end subroutine read_MESH_drainage_database





        !> Read a MESH gridded output file (e.g. 'PREC_H.r2c')
        subroutine read_MESH_OUTFIELD(unitOUTFIELD, metF, dateT)
            use sed_tools
            implicit none

            integer, intent(in) :: unitOUTFIELD
            real, dimension(:,:), intent(out) :: metF
            character(len=80), intent(out) :: dateT

            integer :: i, ios = 0
            character(len=80) :: keyword, line, string1, string2, delim
            !integer :: FRAME_NO1, FRAME_NO2

            delim = '\t'
            do while(ios.eq.0)
              read(unit=unitOUTFIELD, FMT='((A))', iostat=ios) line
              !print *,'idx: ', index(line, '"'), '  ',len_trim(line)
              !print *, trim(line(index(line, '"'):len_trim(line)))
              dateT = trim(line(index(line, '"'):len_trim(line)))
              call split_string(line, string1, string2, delim)
              !print *, trim(string1), ' ' ,trim(string2)
              !read(unit=unitOUTFIELD, *) keyword, FRAME_NO1, FRAME_NO2, TIME
              keyword = trim(string1)
              if(keyword .eq. ':Frame') then
                  do i = 1, yCount
                      read(unit = unitOUTFIELD, fmt = *) metF(i,:)
                  end do
                  read(unit = unitOUTFIELD, fmt = *)

                  exit
              end if

            end do

        end subroutine read_MESH_OUTFIELD

        !> Read the 'rbm_input.r2c' file produced in MESH.
        subroutine read_rbm_input(unit_rbm_input, DISC, INFL, QDIF, DEPT, WIDT, VELO, dateT)
            use sed_tools
            implicit none

            integer, intent(in) :: unit_rbm_input
            real, dimension(:,:), intent(out) :: DISC, INFL, QDIF, DEPT, WIDT, VELO
            character(len=80), intent(out) :: dateT
            integer :: i, ios = 0
            character(len=80) :: keyword, line, string1, string2, delim
            !integer :: FRAME_NO1, FRAME_NO2

            delim = '\t'
            do while(ios.eq.0)
              read(unit=unit_rbm_input, FMT='((A))', iostat=ios) line
              !print *,'idx: ', index(line, '"'), '  ',len_trim(line)
              !print *, trim(line(index(line, '"'):len_trim(line)))
              dateT = trim(line(index(line, '"'):len_trim(line)))
              call split_string(line, string1, string2, delim)
              !print *, trim(string1), ' ' ,trim(string2)
              !read(unit=unitOUTFIELD, *) keyword, FRAME_NO1, FRAME_NO2, TIME
              keyword = trim(string1)
              if(keyword .eq. ':Frame') then
                  do i = 1, yCount
                      read(unit = unit_rbm_input, fmt = *) DISC(i,:)
                  end do
                  read(unit = unit_rbm_input, fmt = *)
                  read(unit = unit_rbm_input, fmt = *)
                  do i = 1, yCount
                      read(unit = unit_rbm_input, fmt = *) INFL(i,:)
                  end do
                  read(unit = unit_rbm_input, fmt = *)
                  read(unit = unit_rbm_input, fmt = *)
                  do i = 1, yCount
                      read(unit = unit_rbm_input, fmt = *) QDIF(i,:)
                  end do
                  read(unit = unit_rbm_input, fmt = *)
                  read(unit = unit_rbm_input, fmt = *)
                  do i = 1, yCount
                      read(unit = unit_rbm_input, fmt = *) DEPT(i,:)
                  end do
                  read(unit = unit_rbm_input, fmt = *)
                  read(unit = unit_rbm_input, fmt = *)
                  do i = 1, yCount
                      read(unit = unit_rbm_input, fmt = *) WIDT(i,:)
                  end do
                  read(unit = unit_rbm_input, fmt = *)
                  read(unit = unit_rbm_input, fmt = *)
                  do i = 1, yCount
                      read(unit = unit_rbm_input, fmt = *) VELO(i,:)
                  end do
                  read(unit = unit_rbm_input, fmt = *)

                  exit
              end if

            end do
        end subroutine read_rbm_input



        !> Read meteoroligcal variables from opened file
        subroutine readMetData
            implicit none

            !>- Read precipitation data
            !call read_MESH_OUTFIELD(unitOUTFIELDprecip, dummy, dateTime)
            call interp_MESH_OUTFIELD(filepathOUTFIELD, filenamePR, dummy, dateIn)
            do i = 1, NA
                mv(i)%precip = dummy(ipos(i),jpos(i))
            end do
            !>- Read precipitation data
            !call read_MESH_OUTFIELD(unitOUTFIELDevap, dummy, dateTime)
            call interp_MESH_OUTFIELD(filepathOUTFIELD, filenameEV, dummy, dateIn)
            do i = 1, NA
                mv(i)%evapotran = dummy(ipos(i),jpos(i))
            end do
        end subroutine readMetData

        !> Read overland flow data
        subroutine readOverLandFlowData
            implicit none

            !call read_MESH_OUTFIELD(unitOUTFIELDwr_runovf, dummy, dateTime)
            call interp_MESH_OUTFIELD(filepathOUTFIELD, filenameWR, dummy, dateIn)
            do i = 1, NA
                ofh(i)%discharge = (dummy(ipos(i),jpos(i))/1000)*shd%AREA(ipos(i),jpos(i))/3600 ! Discharge in the cell in m3/s. 1h=3600s
                ofh(i)%depth = (0.34*ofh(i)%discharge**0.341)*1000.             ! Water depth in mm According to Allen et al, 1994
                ofh(i)%slope = shd%SLOPE_CHNL(ipos(i),jpos(i))                ! Water surface slope.
                ofh(i)%width = sqrt( shd%AREA(ipos(i),jpos(i)) )                 ! Flow width = cell width (m).
                if (ofh(i)%depth /= 0.) then
                    ofh(i)%velocity = ofh(i)%discharge/((ofh(i)%depth/1000)*ofh(i)%width)               ! Water flow velocity (m/s).
                else
                    ofh(i)%velocity = 0.
                end if

            end do
        end subroutine readOverLandFlowData


        subroutine readInStreamFlowData
            implicit none

            !call read_rbm_input(unit_rbm_input, DISC, INFL, QDIF, DEPT, WIDT, VELO, dateTime)
            call interp_rbm_input(filepathRBM, filenameRBM, DISC, dateIn, 1)
            call interp_rbm_input(filepathRBM, filenameRBM, INFL, dateIn, 2)
            call interp_rbm_input(filepathRBM, filenameRBM, QDIF, dateIn, 3)
            call interp_rbm_input(filepathRBM, filenameRBM, DEPT, dateIn, 4)
            call interp_rbm_input(filepathRBM, filenameRBM, WIDT, dateIn, 5)
            call interp_rbm_input(filepathRBM, filenameRBM, VELO, dateIn, 6)
            do i = 1, NA
                rh(i)%slope = shd%SLOPE_CHNL(ipos(i),jpos(i)) ! In-stream water surface slope. It is aproximate to the channel bottom slope.
                rh(i)%velocity = VELO(ipos(i),jpos(i)) ! Channel flow velocity (m/s)
                rh(i)%width =  WIDT(ipos(i),jpos(i)) ! Channel width (m)
                rh(i)%depth = DEPT(ipos(i),jpos(i)) ! Channel water depth (m)
                rh(i)%discharge = DISC(ipos(i),jpos(i)) ! In-stream discharge (m3/s). Assume a rectagular section
                rh(i)%length = shd%CHNL_LEN(ipos(i),jpos(i)) ! Channel reach length
            end do
        end subroutine readInStreamFlowData


        !> Inerpolate matrix in MESH OUTFIELDS files
        subroutine interp_MESH_OUTFIELD(filepath, filename, matF, inDate)
            !use sed_tools
            implicit none

            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath
            real, dimension(:,:), intent(out) :: matF
            character(len=80), intent(in) :: inDate

            integer :: i, ios = 0
            character(len=80) :: line, date1, date2
            type(IterDate) :: dateC, date11, date22
            integer :: fileunit = 100
            real, dimension(yCount,xCount)  :: mat1, mat2
            real :: secsC, secs11, secs22, dsecs11 !, dsecs22


            !print *,  trim(ADJUSTL(filepath)) // trim(ADJUSTL(filename))
            open(unit = fileunit, &
                !file=trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/" // &
                 !trim(ADJUSTL(metVar)) //  ".r2c", &
                file =  trim(ADJUSTL(filepath)) // trim(ADJUSTL(filename)), &
                status="old", &
                action="read", &
                iostat=ios)

            if (ios /= 0) then
                close(unit = fileunit)
                print *, "The file ", trim(ADJUSTL(filename)), " could not be opened."
                print *, "Stopping..."
                stop

            else
                !print *, 'OPENING: ', trim(ADJUSTL(filename))
            end if

            !call parseDateComponent(inDate, dateC)
            call parseDateTime(inDate, dateC)
            secsC = dateC%hour*3600 + dateC%mins*60 + dateC%secs
            !print *, dateC%year, ' ', dateC%month, ' ', dateC%day, ' ', dateC%hour, &
            !        dateC%mins, ' ', dateC%secs

            !> Read the file header
            do
                read(unit=fileunit, FMT='((A))', iostat=ios) line
                if (trim(line) .eq. ':endHeader') exit
            end do

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date1 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                  read(unit = fileunit, fmt = *, iostat=ios) mat1(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            do while(ios.eq.0)
              !if (begi==0) read(unit=fileunit, FMT='((A))', iostat=ios) line
              !if (trim(line) .eq. ':endHeader') begi = 1
              !if (begi==1) then
                !> Record date and matrix of date before

!                if (date1==inDate) then
!                    matF = mat1
!                    print *, 'date1 --', date1
!                    exit
!                end if

                !> Record date and matrix of date before
                read(unit=fileunit, FMT='((A))', iostat=ios) line
                date2 = trim(line(index(line, '"'):len_trim(line)))
                do i = 1, yCount
                      read(unit = fileunit, fmt = *, iostat=ios) mat2(i,:)
                end do
                read(unit = fileunit, fmt = *, iostat=ios)
!                if (date2==inDate) then
!                    matF = mat2
!                    print *, 'date2 --', date2
!                    exit
!                end if

                !call parseDateComponent(date1, date11)
                call parseDateTime(date1, date11)
                !call parseDateComponent(date2, date22)
                call parseDateTime(date2, date22)

                if (date11%year == dateC%year .or. date22%year == dateC%year) then
                    if (date11%month == dateC%month .or. date22%month == dateC%month) then
                        if (date11%day == dateC%day .or. date22%day == dateC%day) then

                            !if (date11%hour <= dateC%hour .and. date22%hour >= dateC%hour) then
                            secs11 = date11%hour*3600 + date11%mins*60 + date11%secs
                            secs22 = date22%hour*3600 + date22%mins*60 + date22%secs
                            if (secs11 <= secsC .and. secs22 >= secsC) then
                                !if (date11%mins <= dateC%mins .and. date22%mins >= dateC%mins) then
                                !    if (date11%secs <= dateC%secs .and. date22%secs >= dateC%secs) then
                                        dsecs11 = secsC - secs11
                                        !dsecs22 = secs22 - secsC

                                        !if (dsecs11 < dsecs22) then
                                            !matF(:,:) = mat1(:,:)
                                        matF = (mat2-mat1)*(dsecs11/3600) + mat1
                                        !print *,'date1 ', date1
                                        !print *,'date2 ', date2
                                        !print *, 'heeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
                                        !else
                                        !    !matF(:,:) = mat2(:,:)
                                        !end if
                                        exit
                                !    end if
                                !end if
                            end if
                            if (dateC%hour == 0)then
                                dsecs11 = secsC
                                matF = (mat2-mat1)*(dsecs11/3600) + mat1
                                !print *,'date1 ', date1
                                !print *,'date2 ', date2
                                exit
                            end if
                        end if
                    end if
                end if
                date1=date2
                mat1=mat2

            end do

            close(unit = fileunit)

        end subroutine interp_MESH_OUTFIELD


        !> Inerpolate matrix in rbm_input_file
        subroutine interp_rbm_input(filepath, filename, matF, inDate, varOrder)
            !use sed_tools
            implicit none

            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath
            real, dimension(:,:), intent(out) :: matF
            character(len=80), intent(in) :: inDate
            integer, intent(in) :: varOrder

            integer :: i, ios = 0
            character(len=80) :: line, date1, date2
            type(IterDate) :: dateC, date11, date22
            integer :: fileunit = 100, nvars = 6
            real, dimension(yCount,xCount)  :: mat1, mat2
            real :: secsC, secs11, secs22, dsecs11 !, dsecs22


            !print *,  trim(ADJUSTL(filepath)) // trim(ADJUSTL(filename))
            open(unit = fileunit, &
                !file=trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/" // &
                 !trim(ADJUSTL(metVar)) //  ".r2c", &
                file =  trim(ADJUSTL(filepath)) // trim(ADJUSTL(filename)), &
                status="old", &
                action="read", &
                iostat=ios)

            if (ios /= 0) then
                close(unit = fileunit)
                print *, "The file ", trim(ADJUSTL(filename)), " could not be opened."
                print *, "Stopping..."
                stop

            else
                !print *, 'OPENING: ', trim(ADJUSTL(filename))
            end if

            !call parseDateComponent(inDate, dateC)
            call parseDateTime(inDate, dateC)
            secsC = dateC%hour*3600 + dateC%mins*60 + dateC%secs

            ! Read header
            do
                read(unit=fileunit, FMT='((A))', iostat=ios) line
                if (trim(line) .eq. ':endHeader') exit
            end do

            !> Read lines to omit
            do i = 1, (yCount+2)*(varOrder-1)
                  read(unit=fileunit, FMT='((A))', iostat=ios) line
            end do

            !> Read initial block
            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date1 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                  read(unit = fileunit, fmt = *, iostat=ios) mat1(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)


            do while(ios.eq.0)

                !> Read lines to omit
                do i = 1, (yCount+2)*(nvars -1)
                  read(unit=fileunit, FMT='((A))', iostat=ios) line
                end do

                !> Record date and matrix of date before
                read(unit=fileunit, FMT='((A))', iostat=ios) line
                date2 = trim(line(index(line, '"'):len_trim(line)))
                do i = 1, yCount
                      read(unit = fileunit, fmt = *, iostat=ios) mat2(i,:)
                end do
                read(unit = fileunit, fmt = *, iostat=ios)

                !call parseDateComponent(date1, date11)
                call parseDateTime(date1, date11)
                !call parseDateComponent(date2, date22)
                call parseDateTime(date2, date22)

                if (date11%year == dateC%year .or. date22%year == dateC%year) then
                    if (date11%month == dateC%month .or. date22%month == dateC%month) then
                        if (date11%day == dateC%day .or. date22%day == dateC%day) then
                            secs11 = date11%hour*3600 + date11%mins*60 + date11%secs
                            secs22 = date22%hour*3600 + date22%mins*60 + date22%secs
                            if (secs11 <= secsC .and. secs22 >= secsC) then

                                dsecs11 = secsC - secs11
                                matF = (mat2-mat1)*(dsecs11/3600) + mat1
                                exit

                            end if
                            if (dateC%hour == 0)then
                                dsecs11 = secsC
                                matF = (mat2-mat1)*(dsecs11/3600) + mat1
                                exit
                            end if
                        end if
                    end if
                end if
                date1=date2
                mat1=mat2

            end do

            close(unit = fileunit)

        end subroutine interp_rbm_input


        !> Parse a stringdate (e.g. "2000/01/01 01:10:05.000") into datet
        subroutine parseDateComponent(stringdate, datet)
            use sed_tools
            implicit none

            character(len=80), intent(in) :: stringdate
            type(IterDate), intent(out) :: datet

            character(len=80) :: delim, string1, string2

            delim = '\t'
            call split_string(stringdate, string1, string2, delim)
            read(string1(2:5),'(i5)') datet%year != nint(trim(string1(2:5)))
            read(string1(7:8),'(i5)') datet%month != nint(trim(string1(7:8)))
            read(string1(10:11),'(i5)') datet%day != nint(trim(string1(10:11)))
            read(string2(1:2),'(i5)') datet%hour != nint(trim(string2(1:2)))
            read(string2(4:5),'(i5)') datet%mins != nint(trim(string2(4:5)))
            read(string2(7:8),'(i5)') datet%secs != nint(trim(string2(7:8)))
            !datet%year = trim( string1(1:index(dateC, '/')-1) )
            !date = trim(dateC(index(dateC, '"')+1:len_trim(dateC)))
            !time = trim(timeC(1:len_trim(timeC)-1))
            !hour = trim(dateC(index(dateC, '"')

        end subroutine

        subroutine parseDateTime(stringDateTime, datet)
            use sed_tools
            implicit none
            character(len=80), intent(in) :: stringDateTime
            character(len=80) :: delim, datee, time, year, string, string1,  month, day, hour, mins, secs
            type(IterDate), intent(out) :: datet

            delim = ' '
            call split_string(stringDateTime, datee, time, delim)

            ! For date
            delim = '/'
            call split_string(datee, year, string, delim)
            year = adjustl(trim(year(index(year, '"')+1:len_trim(year))))
            delim = '/'
            call split_string(string, month, day, delim)
            !month = adjustl(month)
            !day = adjustl(day)

            ! For time
            time=adjustl(time) !"2002/1/3  1:00:00.000"
            !print *, time
            delim = ':'
            call split_string(time, hour, string, delim)
            !hour = adjustl(hour)
            !print *, 'H ',hour
            !print *, string
            delim = ':'
            call split_string(string, mins, string1, delim)
            !mins = adjustl(mins)
            !print *, 'MM ',mins
            delim = '.'
            call split_string(string1, secs, string, delim)
            !secs = adjustl(secs)

            !print *, year, ' ', month, ' ', day, ' ', hour, ' ', mins, ' ', secs


            read(year,'(i5)') datet%year != nint(trim(year))
            !print *, 'hereer1'
            read(month,'(i5)') datet%month != nint(trim(month))
            !print *, 'hereer2'
            read(day,'(i5)') datet%day != nint(trim(day))
            !print *, 'hereer3'
            read(hour,'(i5)') datet%hour != nint(trim(hour))
            read(mins,'(i5)') datet%mins != nint(trim(mins))
            read(secs,'(i5)') datet%secs != nint(trim(secs))


            !print *, 'S ',secs

            !print *, datet%year, ' ', datet%month, ' ', datet%day, ' ', datet%hour, ' ', datet%mins, ' ', datet%secs
            !print *, year, ' ', month, ' ', day, ' ', hour, ' ', mins, ' ', secs

        end subroutine parseDateTime



!        subroutine read_MESH_parameters_sediment
!
!
!            !use sed_vars
!
!            implicit none
!
!            ! Allocating memory for variables
!            !allocate(soilVul(ngru), damEffecSw(ngru))
!            !allocate(rainSoilEr(ngru))
!
!
!            ! Open the file
!            open(unit=unitpar,file="MESH_parameters_sediment.ini",status="old",action="read",iostat=ios)
!            if (ios /= 0) then
!                print *, "The file ", "MESH_parameters_sediment.dat", " does not exist."
!                print *, "Stopping..."
!                stop
!            end if
!
!            ! Read soil data
!            read(unit = unitpar, fmt = *) overlFlowCapaMethod
!            read(unit = unitpar,fmt = *)
!            read(unit = unitpar,fmt = *)
!            read(unit = unitpar, fmt = *) nsoilType
!            allocate(soilChar(nsoilType))
!            read(unit = unitpar, fmt = *) (soilChar(i)%name, i=1,nsoilType)
!            read(unit = unitpar, fmt = *) (soilChar(i)%diameter, i=1,nsoilType)
!            read(unit = unitpar, fmt = *) (soilChar(i)%density, i=1,nsoilType)
!            read(unit = unitpar, fmt = *) (soilChar(i)%surfPorosity, i=1,nsoilType)
!            read(unit = unitpar, fmt = *) (soilChar(i)%soilDetach, i=1,nsoilType)
!            read(unit = unitpar, fmt = *) (soilChar(i)%overlandDetach, i=1,nsoilType)
!
!            read(unit = unitpar, fmt = *)
!            !allocate(cellSoil(yCount,xCount))
!            allocate(cellSoil(NA))
!            !read(unit = unitpar, fmt = *) dummyC
!            do i = 1, yCount
!                read(unit = unitpar, fmt = *) dummyC(i,:)
!            end do
!            do k = 1, NA
!                !print*, ipos(k),jpos(k)
!                cellSoil(k) = dummyC(ipos(k),jpos(k))
!            end do
!
!    !        k = 1
!    !        do i=1,yCount
!    !            do j = 1, xCount
!    !                if(rank(i,j)>0) then
!    !                    cellSoil(k)
!    !                    k = k + 1
!    !                end if
!    !                read(unit = unitpar, fmt = *) cellSoil(i,:)
!    !            end do
!    !        end do
!
!            read(unit = unitpar, fmt = *)
!            allocate(cellGroundCov(NA))
!            do i = 1, yCount
!                read(unit = unitpar, fmt = *) dummy(i,:)
!            end do
!            do k = 1, NA
!                cellGroundCov(k) = dummy(ipos(k),jpos(k))
!            end do
!    !        allocate(cellGroundCov(yCount,xCount))
!    !        do i=1,yCount
!    !            read(unit = unitpar, fmt = *) cellGroundCov(i,:)
!    !        end do
!
!            ! Read vegetation
!            read(unit = unitpar, fmt = *)
!            read(unit = unitpar, fmt = *) nvegeType
!            allocate(vegeChar(nvegeType))
!            read(unit = unitpar, fmt = *) (vegeChar(i)%name, i=1,nvegeType)
!            read(unit = unitpar, fmt = *) (vegeChar(i)%fallHeight, i=1,nvegeType)
!            read(unit = unitpar, fmt = *) (vegeChar(i)%dropDiam, i=1,nvegeType)
!            read(unit = unitpar, fmt = *) (vegeChar(i)%percDrip, i=1,nvegeType)
!
!            read(unit = unitpar, fmt = *)
!            allocate(cellVege(NA))
!            do i = 1, yCount
!                read(unit = unitpar, fmt = *) dummyC(i,:)
!            end do
!            do k = 1, NA
!                cellVege(k) = dummyC(ipos(k),jpos(k))
!            end do
!    !
!    !        do i=1,yCount
!    !            read(unit = unitpar, fmt = *) cellVege(i,:)
!    !        end do
!
!            read(unit = unitpar, fmt = *)
!            allocate(cellCanopyCov(NA))
!            do i = 1, yCount
!                read(unit = unitpar, fmt = *) dummy(i,:)
!            end do
!            do k = 1, NA
!                cellCanopyCov(k) = dummy(ipos(k),jpos(k))
!            end do
!    !        allocate(cellCanopyCov(yCount,xCount))
!    !        do i=1,yCount
!    !            read(unit = unitpar, fmt = *) cellCanopyCov(i,:)
!    !        end do
!
!
!            !read(unit = unitpar, fmt = *) (soilVul(i), i=1,ngru)
!            !read(unit = unitpar, fmt = *) (soilVul(i), i=1,ngru)
!            !read(unit = unitpar, fmt = *) (damEffecSw(i), i=1,ngru)
!
!            ! Close the file
!            close(unit = unitpar)
!
!
!
!        end subroutine read_MESH_parameters_sediment



end module sed_input
