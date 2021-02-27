!------------------------------------------------------------------------------
! Sediment transport in cold region catchments: the MESH-SED model
!------------------------------------------------------------------------------
!
!> @brief
!> MODULE: sed_input
!>
!> @detail This module contains subroutines and functions to read the input
!> required by MESH-SED
!>
!> @author Luis Morales (LAM), GIWS & GWF.
!> - July, 2017
!> @date January, 2019-LAM
!> - Documenting the code
!> @todo There are some unused subroutines/functions that can be commented
!> in order ignore them in the compilation.
!---------------------------------------------------------------------------
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
!            read(unit = unitParam,fmt = *)
!            read(unit = unitParam,fmt = *) OUTFIELDfolder
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

            getIdSedClass = 0
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
                !print *, i, ipos(i),jpos(i), gca(i)%soilType
            end do
            !stop

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
            
            read(unit = unitGrid,fmt = *)
            do i = 1, yCount
                read(unit = unitGrid, fmt = *) dummyI(i,:) !tempoR(i,:)
            end do
            do i = 1, NA
                riverType(i) = dummyI(ipos(i),jpos(i)) !tempoR(ipos(i),jpos(i))
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
            read(unit = unitSoVeCha,fmt = *) ratioStress

            ! Read river type
            read(unit = unitSoVeCha,fmt = *)
            read(unit = unitSoVeCha, fmt = *) NRIVT
            allocate(rtp(NRIVT))
            do j = 1,NRIVT
                read(unit = unitSoVeCha, fmt = *) rtp(j)%a, rtp(j)%b
            end do

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

        !> Read "MESH_sed_reservoir.ini" file that characterize reservoir in the basin
        subroutine read_sed_reservoir
            implicit none

            ! Open the file
            open(unit = unitReser,file=trim(ADJUSTL(casedir)) // "MESH_sed_reservoir.ini", &
                status="old",action="read",iostat=ios)
            if (ios /= 0) then
                close(unit = unitParam)
                print *, "The file ", "MESH_sed_reservoir.ini", " could not be opened."
                print *, "Stopping..."
                stop

            else
                print *, 'READING: MESH_sed_reservoir.ini '
            end if
            read(unit = unitReser,fmt = *)
            read(unit = unitReser,fmt = *) NRESE
            read(unit = unitReser,fmt = *)
            read(unit = unitReser,fmt = *)
            allocate(rese(NRESE))
            do i = 1, NRESE
                read(unit = unitReser, fmt = '(A15, I10, 2(F15.2))') rese(i)%name, rese(i)%rank, rese(i)%area, rese(i)%vol
                !print *, rese(i)%name, rese(i)%rank, rese(i)%area, rese(i)%vol
            end do

            close(unit = unitReser)
        end subroutine read_sed_reservoir


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

            !call interp_MESH_OUTFIELD(filepathOUTFIELD, filenamePR, dummy, dateIn)
            call interp2_MESH_OUTFIELD(unitPR, mat1PR, date1PR, mat2PR, date2PR, dummy, dateIn, iosPR)
            do i = 1, NA
                mv(i)%precip = dummy(ipos(i),jpos(i)) ! In mm/h
            end do

            !>- Read evapotranspiration data
            !call read_MESH_OUTFIELD(unitOUTFIELDevap, dummy, dateTime)

            !call interp_MESH_OUTFIELD(filepathOUTFIELD, filenameEV, dummy, dateIn)
            call interp2_MESH_OUTFIELD(unitEV, mat1EV, date1EV, mat2EV, date2EV, dummy, dateIn, iosEV)
            do i = 1, NA
                mv(i)%evapotran = dummy(ipos(i),jpos(i)) ! In mm/h
            end do
        end subroutine readMetData

        !> Read overland flow data
        subroutine readOverLandFlowData
            implicit none

            !call read_MESH_OUTFIELD(unitOUTFIELDwr_runovf, dummy, dateTime)

            !call interp_MESH_OUTFIELD(filepathOUTFIELD, filenameWR, dummy, dateIn)
            call interp2_MESH_OUTFIELD(unitWR, mat1WR, date1WR, mat2WR, date2WR, dummy, dateIn, iosWR)

            do i = 1, NA
                ofh(i)%discharge = (dummy(ipos(i),jpos(i))/1000)*shd%AREA(ipos(i),jpos(i))/3600 ! From mm per hour to m3/s. Discharge in the cell in m3/s. 1h=3600s
                ofh(i)%depth = (0.34*ofh(i)%discharge**0.341)*1000.             ! Water depth in mm According to Allen et al, 1994
                ofh(i)%slope = shd%SLOPE_CHNL(ipos(i),jpos(i))                ! Water surface slope.
                ofh(i)%width = sqrt( shd%AREA(ipos(i),jpos(i)) )                 ! Flow width = cell width (m).
                if (ofh(i)%depth /= 0.) then
                    ofh(i)%velocity = ofh(i)%discharge/((ofh(i)%depth/1000)*ofh(i)%width)               ! Water flow velocity (m/s).
                else
                    ofh(i)%velocity = 0.
                end if

            end do
            !print *, '->AREA', minval(ofh(:)%width)**2, ' ', maxval(ofh(:)%width)**2
!            print *, '->DIS (m3/s): ', minval(ofh(:)%discharge), ' ', maxval(ofh(:)%discharge)
!            print *, '->DE (mm): ', minval(ofh(:)%depth), ' ', maxval(ofh(:)%depth)
!            print *, '->SL (mm): ', minval(ofh(:)%slope), ' ', maxval(ofh(:)%slope)
!            print *, '->WI (mm): ', minval(ofh(:)%width), ' ', maxval(ofh(:)%width)
!            print *, '->VE (mm): ', minval(ofh(:)%velocity), ' ', maxval(ofh(:)%velocity)

        end subroutine readOverLandFlowData


!        subroutine readInStreamFlowData
!            implicit none
!
!            !call read_rbm_input(unit_rbm_input, DISC, INFL, QDIF, DEPT, WIDT, VELO, dateTime)
!
!            !call interp_rbm_input(filepathRBM, filenameRBM, DISC, dateIn, 1)
!            !call interp_rbm_input(filepathRBM, filenameRBM, INFL, dateIn, 2)
!            !call interp_rbm_input(filepathRBM, filenameRBM, QDIF, dateIn, 3)
!            !call interp_rbm_input(filepathRBM, filenameRBM, DEPT, dateIn, 4)
!            !call interp_rbm_input(filepathRBM, filenameRBM, WIDT, dateIn, 5)
!            !call interp_rbm_input(filepathRBM, filenameRBM, VELO, dateIn, 6)
!
!            call interp2_rbm_input(unitRBM, mat1RBM1, date1RBM1, mat1RBM2, date1RBM2, &
!                                     mat1RBM3, date1RBM3, mat1RBM4, date1RBM4, mat1RBM5, date1RBM5, mat1RBM6, date1RBM6, &
!                                     mat2RBM1, date2RBM1, mat2RBM2, date2RBM2, mat2RBM3, date2RBM3, mat2RBM4, date2RBM4, &
!                                     mat2RBM5, date2RBM5, mat2RBM6, date2RBM6, DISC, INFL, QDIF, DEPT, WIDT, &
!                                     VELO, dateIn, iosRBM)
!
!!            print *, '->DISC : ', sum(mat1RBM1), ' ',sum(DISC) ,' ', sum(mat2RBM1)
!!            print *, '->INFL : ', sum(mat1RBM2), ' ',sum(INFL) ,' ', sum(mat2RBM2)
!!            print *, '->QDIF : ', sum(mat1RBM3), ' ',sum(QDIF) ,' ', sum(mat2RBM3)
!!            print *, '->DEPT : ', sum(mat1RBM4), ' ',sum(DEPT) ,' ', sum(mat2RBM4)
!!            print *, '->WIDT : ', sum(mat1RBM5), ' ',sum(WIDT) ,' ', sum(mat2RBM5)
!!            print *, '->VELO : ', sum(mat1RBM6), ' ',sum(VELO) ,' ', sum(mat2RBM6)
!
!            do i = 1, NA
!                rh(i)%slope = shd%SLOPE_CHNL(ipos(i),jpos(i)) ! In-stream water surface slope. It is aproximate to the channel bottom slope.
!                rh(i)%velocity = VELO(ipos(i),jpos(i)) ! Channel flow velocity (m/s)
!                rh(i)%width =  WIDT(ipos(i),jpos(i)) ! Channel width (m)
!                rh(i)%depth = DEPT(ipos(i),jpos(i)) ! Channel water depth (m)
!                rh(i)%discharge = DISC(ipos(i),jpos(i)) ! In-stream discharge (m3/s). Assume a rectagular section
!                rh(i)%length = shd%CHNL_LEN(ipos(i),jpos(i)) ! Channel reach length
!                if (i == 1325 ) then
!                    print *, '-->', i, rh(i)%velocity, rh(i)%depth, 0.34*rh(i)%discharge**0.341, rh(i)%width, rh(i)%discharge
!                end if
!
!            end do
!!            print *, '->DISC : ', minval(rh(:)%discharge), ' ',maxval(rh(:)%discharge)
!!            print *, '->INFL : ', minval(INFL), ' ',maxval(INFL)
!!            print *, '->DEPT : ', minval(rh(:)%depth), ' ',maxval(rh(:)%depth)
!!            print *, '->WIDT : ', minval(rh(:)%width), ' ',maxval(rh(:)%width)
!!            print *, '->VELO : ', minval(rh(:)%velocity), ' ',maxval(rh(:)%velocity)
!!            print *, '->LENG : ', minval(rh(:)%length), ' ',maxval(rh(:)%length)
!
!        end subroutine readInStreamFlowData


        subroutine readMESHoutputData
            implicit none

            real :: mannDisc, mannVel

            call interp2_sed_input(mat1SED1, date1SED1, mat1SED2, date1SED2, &
                                     mat1SED3, date1SED3, mat1SED4, date1SED4, mat1SED5, date1SED5, mat1SED6, date1SED6, &
                                     mat1SED7, date1SED7, mat1SED8, date1SED8, mat1SED9, date1SED9, mat1SED10, date1SED10, &
                                     mat1SED11, date1SED11, &
                                     mat2SED1, date2SED1, mat2SED2, date2SED2, mat2SED3, date2SED3, mat2SED4, date2SED4, &
                                     mat2SED5, date2SED5, mat2SED6, date2SED6, mat2SED7, date2SED7, mat2SED8, date2SED8, &
                                     mat2SED9, date2SED9, &
                                     mat2SED10, date2SED10, mat2SED11, date2SED11, dateIn)

            do i = 1, NA
                !> Met variables

                !*  7.  Precipitation (mm h-1). Note: Accumulated over the time-step.
                if (sedi%PREP(ipos(i),jpos(i)) < 0.) then ! Check for negative values
                  mv(i)%precip = 0.
                else
                  mv(i)%precip = sedi%PREP(ipos(i),jpos(i))
                end if

                !*  8.  Evapotranspiration (mm h-1). Note: Of evapotranspiration accumulated over the time-step. UNITS????????????????
                if (sedi%EVAP(ipos(i),jpos(i)) < 0.) then ! Check for negative values
                  mv(i)%evapotran = 0.
                else
                  mv(i)%evapotran = sedi%EVAP(ipos(i),jpos(i))
                end if


                !> Overland flow hydraulic variables
                ofh(i)%depth = sedi%OFDE(ipos(i),jpos(i)) !*  9.  Overland water depth (mm). Note: Accumulated over the time-step.
                ofh(i)%slope = sedi%LASL(ipos(i),jpos(i)) !*  10. Surface slope (m m-1). SLOPE_INT isn't used in CLASS, so average slope from tiles in cell?
                ofh(i)%width = sedi%CEWI(ipos(i),jpos(i)) !*  11. Cell width (m).
                ofh(i)%discharge = (ofh(i)%depth*0.001/0.34)**(1/0.341) ! Overland flow discharge in m3/s. Eq. According to Allen et al, 1994
                mannDisc = (1/0.060)*(ofh(i)%depth*0.001*ofh(i)%width)*(ofh(i)%slope**0.5)* &
                            (ofh(i)%depth*0.001*ofh(i)%width/(2*ofh(i)%depth*0.001 + ofh(i)%width))**(2./3)
                if (ofh(i)%depth /= 0.) then
                    ofh(i)%velocity = ofh(i)%discharge/((ofh(i)%depth/1000.)*ofh(i)%width)               ! Water flow velocity (m/s).
                    mannVel = mannDisc/((ofh(i)%depth/1000.)*ofh(i)%width)
                else
                    ofh(i)%velocity = 0.
                    mannVel  = 0.
                end if
                !if (i == 1325 ) then
                    !print *, '-->', i, ofh(i)%depth, ofh(i)%discharge, mannDisc, ofh(i)%velocity, mannVel
                    !print *, i, mv(i)%precip, mv(i)%evapotran
                !end if


                !> In-stream variables
                !> @todo rh(i)%slope was to high so it was elevated to power of 2

                rh(i)%discharge = sedi%DISC(ipos(i),jpos(i)) !*  1.  Average flow (discharge) (m3 s-1). Note: Averaged over the time-step.
                rh(i)%length = sedi%CHLE(ipos(i),jpos(i)) !*  4.  Channel length (m).
                if (rh(i)%wbody == 0) then  !> For channels
                    !> rh(i)%width = rhwidthc*sedi%WIDT(ipos(i),jpos(i)) !*  3.  Channel width (m).
                    rh(i)%width = sedi%WIDT(ipos(i),jpos(i)) !*  3.  Channel width (m).
                    rh(i)%slope = 1.*( sedi%CHSL(ipos(i),jpos(i)) )**2 !*  5.  Channel slope (m m-1). slope = sqrt(SLOPE_CHNL)
                    rh(i)%velocity = rtpc(i)%a*sedi%VELO(ipos(i),jpos(i)) !*  6.  Stream velocity (m s-1). Take stream speed to be average flow (m3 s-1) divided by channel x-sec area (m2) (from rte_sub.f).
                    !> rh(i)%velocity = rtpc(i)%a*(rh(i)%discharge)**rtpc(i)%b !* Based on Leopold and
                    
                    !rh(i)%depth = sedi%DEPT(ipos(i),jpos(i)) !*  2.  Channel depth (m).
                    !rh(i)%depth = 0.34*rh(i)%discharge**0.341 !*  Channel depth (m). Estimated based on Allen et al, 1994 as MESH water depth is cte

                    !>  Channel depth (m). Estimated assuming rectangular section
                    if (rh(i)%velocity <= 1e-6  .OR. rh(i)%width <= 1e-6) then
                        rh(i)%depth = 0.
                    else
                        rh(i)%depth = rtpc(i)%b*rh(i)%discharge/(rh(i)%velocity * rh(i)%width)
                        !rh(i)%depth = rh(i)%discharge/(rh(i)%velocity * rh(i)%width)
                    end if

                else                 !> for reservoir
                    rh(i)%width = rh(i)%areaRe/rh(i)%length
                    rh(i)%slope = 0.
                    rh(i)%velocity = rh(i)%discharge/rh(i)%areaRe

                end if

                !if (i == 1325 ) then
                    !print *, '-->', i, rh(i)%velocity, rh(i)%depth, 0.34*rh(i)%discharge**0.341, rh(i)%width, rh(i)%discharge
                !end if

            end do
            !print *, mv(324)%precip, ofh(324)%discharge, ofh(324)%depth, rh(324)%discharge, rh(324)%depth


        end subroutine readMESHoutputData

        subroutine readHeader_MESH_OUTPUTFILES(filepath, filename, fileunit, ios)
            implicit none
            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath
            integer, intent(in) :: fileunit
            integer, intent(out) :: ios

            character(len=80) :: line
            !real :: secs11
            !integer :: ios = 0

            !> Open the file
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

            !> Read the header
            do
                read(unit=fileunit, FMT='((A))', iostat=ios) line
                if (trim(line) .eq. ':endHeader') exit
            end do
        end subroutine readHeader_MESH_OUTPUTFILES

        subroutine getMat1_MESH_OUTFIELD(filepath, filename, fileunit, mat1, date1, mat2, date2, inDate, ios)
            implicit none

            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath
            integer, intent(in) :: fileunit
            character(len=80), intent(in) :: inDate
            character(len=80), intent(out) :: date1, date2
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: mat1, mat2

            type(IterDate) :: dateC, date11
            character(len=80) :: line
            integer :: i

            call readHeader_MESH_OUTPUTFILES(filepath, filename, fileunit, ios)

            call parseDateTime(inDate, dateC)
            ! NOTE: MESH gridded output files go from 1H to 24H. Our program begin 0H and end at 23H.
            ! 1 is added below to adjust our hour accoridingly.
            dateC%hour = dateC%hour + 1
            !secsC = dateC%hour*3600 + dateC%mins*60 + dateC%secs

            do while(ios.eq.0)

                !> Read lines to omit
                !do i = 1, (yCount+2)*(varOrder-1)
                !      read(unit=fileunit, FMT='((A))', iostat=ios) line
                !end do

                read(unit=fileunit, FMT='((A))', iostat=ios) line
                date1 = trim(line(index(line, '"'):len_trim(line)))
                call parseDateTime(date1, date11)
                do i = 1, yCount
                    read(unit = fileunit, fmt = *, iostat=ios) mat1(i,:)
                end do
                read(unit = fileunit, fmt = *, iostat=ios)


                if (date11%year == dateC%year) then
                    if (date11%month == dateC%month) then
                        if (date11%day == dateC%day) then
                            if (date11%hour == dateC%hour) then
                                exit
                            end if
                        end if
                    end if
                end if
            end do

            call getMat2_MESH_OUTFIELD(fileunit, mat2, date2, ios)
        end subroutine getMat1_MESH_OUTFIELD




!        subroutine getMat1_MESH_OUTFIELD(filepath, filename, fileunit, mat1, date1, mat2, date2, inDate, ios)
!            implicit none
!
!            character(len=80), intent(in) :: filename
!            character(len=150), intent(in) :: filepath
!            integer, intent(in) :: fileunit
!            character(len=80), intent(in) :: inDate
!            character(len=80), intent(out) :: date1, date2
!            integer, intent(inout) :: ios
!            real, dimension(:,:), intent(out) :: mat1, mat2
!
!            type(IterDate) :: dateC, date11
!            character(len=80) :: line
!
!            call readHeader_MESH_OUTPUTFILES(filepath, filename, fileunit, ios)
!
!            call parseDateTime(inDate, dateC)
!            ! NOTE: MESH gridded output files go from 1H to 24H. Our program begin 0H and end at 23H.
!            ! 1 is added below to adjust our hour accoridingly.
!            dateC%hour = dateC%hour + 1
!            !secsC = dateC%hour*3600 + dateC%mins*60 + dateC%secs
!
!            do while(ios.eq.0)
!
!                read(unit=fileunit, FMT='((A))', iostat=ios) line
!                date1 = trim(line(index(line, '"'):len_trim(line)))
!                call parseDateTime(date1, date11)
!                do i = 1, yCount
!                    read(unit = fileunit, fmt = *, iostat=ios) mat1(i,:)
!                end do
!                read(unit = fileunit, fmt = *, iostat=ios)
!
!
!                if (date11%year == dateC%year) then
!                    if (date11%month == dateC%month) then
!                        if (date11%day == dateC%day) then
!                            if (date11%hour == dateC%hour) then
!                                exit
!                            end if
!                        end if
!                    end if
!                end if
!            end do
!
!            call getMat2_MESH_OUTFIELD(fileunit, mat2, date2, ios)
!        end subroutine getMat1_MESH_OUTFIELD


        subroutine getMat2_MESH_OUTFIELD(fileunit, mat2, date2, ios)
            implicit none
                integer, intent(in) :: fileunit
                character(len=80), intent(out) :: date2
                integer, intent(inout) :: ios
                real, dimension(:,:), intent(out) :: mat2

                character(len=80) :: line

                !> Read lines to omit
                !do i = 1, (yCount+2)*(varOrder-1)
                !      read(unit=fileunit, FMT='((A))', iostat=ios) line
                !end do

                read(unit=fileunit, FMT='((A))', iostat=ios) line
                date2 = trim(line(index(line, '"'):len_trim(line)))
                !call parseDateTime(date1, date11)
                do i = 1, yCount
                    read(unit = fileunit, fmt = *, iostat=ios) mat2(i,:)
                end do
                read(unit = fileunit, fmt = *, iostat=ios)
        end subroutine getMat2_MESH_OUTFIELD


        subroutine interp2_MESH_OUTFIELD(fileunit, mat1, date1, mat2, date2, matF, inDate, ios)
            implicit none

            integer, intent(in) :: fileunit
            character(len=80), intent(inout) :: date1, date2
            real, dimension(:,:), intent(inout) :: mat1, mat2
            character(len=80), intent(in) :: inDate
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: matF


            !real, dimension(yCount,xCount) :: matF
            type(IterDate) :: dateC, date22
            real :: secsC
            !real :: dsecs11

            call parseDateTime(inDate, dateC)
            secsC = dateC%mins*60 + dateC%secs

            !call parseDateTime(date1, date11)
            !secs11 = date11%mins*60 + date11%secs

            call parseDateTime(date2, date22)
            !secs22 = date22%mins*60 + date22%secs

            if (dateC%year == date22%year) then
                if (dateC%month == date22%month) then
                    if (dateC%day == date22%day) then
                        if (dateC%hour+1 == date22%hour) then
                            if (secsC >= 0.0) then
                                date1 = date2
                                mat1 = mat2
                                call getMat2_MESH_OUTFIELD(fileunit, mat2, date2, ios)
                            end if
                        end if
                    end if
                end if
            end if

            !dsecs11 = secsC - secs11
            !matF = (mat2-mat1)*(dsecs11/3600) + mat1
            matF = (mat2-mat1)*(secsC/3600) + mat1

        end subroutine interp2_MESH_OUTFIELD



        subroutine getMat1_rbm_input(filepath, filename, fileunit, mat11, date11, mat12, date12, &
                                     mat13, date13, mat14, date14, mat15, date15, mat16, date16, &
                                     mat21, date21, mat22, date22, mat23, date23, mat24, date24, &
                                     mat25, date25, mat26, date26, inDate, ios)
            implicit none

            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath
            integer, intent(in) :: fileunit
            character(len=80), intent(in) :: inDate
            character(len=80), intent(out) :: date11, date12, date13, date14, date15, date16, &
                                              date21, date22, date23, date24, date25, date26
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: mat11, mat12, mat13, mat14, mat15, mat16, &
                                                 mat21, mat22, mat23, mat24, mat25, mat26

            type(IterDate) :: dateC, date111
            character(len=80) :: line
            integer :: i

            call readHeader_MESH_OUTPUTFILES(filepath, filename, fileunit, ios)

            call parseDateTime(inDate, dateC)
            ! NOTE: MESH gridded output files go from 1H to 24H. Our program begin 0H and end at 23H.
            ! 1 is added below to adjust our hour accoridingly.
            dateC%hour = dateC%hour + 1
            !secsC = dateC%hour*3600 + dateC%mins*60 + dateC%secs

            do while(ios.eq.0)

                read(unit=fileunit, FMT='((A))', iostat=ios) line
                date11 = trim(line(index(line, '"'):len_trim(line)))
                call parseDateTime(date11, date111)
                do i = 1, yCount
                    read(unit = fileunit, fmt = *, iostat=ios) mat11(i,:)
                end do
                read(unit = fileunit, fmt = *, iostat=ios)


                if (date111%year == dateC%year) then
                    if (date111%month == dateC%month) then
                        if (date111%day == dateC%day) then
                            if (date111%hour == dateC%hour) then

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date12 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat12(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date13 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat13(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date14 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat14(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date15 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat15(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date16 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat16(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                exit
                            end if
                        end if
                    end if
                end if
            end do

            call getMat2_rbm_input(fileunit, mat21, date21, mat22, date22, mat23, date23, &
                                     mat24, date24, mat25, date25, mat26, date26, ios)
        end subroutine getMat1_rbm_input


        subroutine getMat1_sed_input(filepath, filename, fileunit, mat11, date11, mat12, date12, &
                                     mat13, date13, mat14, date14, mat15, date15, mat16, date16, &
                                     mat17, date17, mat18, date18, mat19, date19, mat110, date110, &
                                     mat111, date111, &
                                     mat21, date21, mat22, date22, mat23, date23, mat24, date24, &
                                     mat25, date25, mat26, date26, mat27, date27, mat28, date28, &
                                     mat29, date29, mat210, date210, mat211, date211, &
                                     inDate, ios)
            implicit none

            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath
            integer, intent(in) :: fileunit
            character(len=80), intent(in) :: inDate
            character(len=80), intent(out) :: date11, date12, date13, date14, date15, date16, &
                                              date17, date18, date19, date110, date111, &
                                              date21, date22, date23, date24, date25, date26, &
                                              date27, date28, date29, date210, date211
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: mat11, mat12, mat13, mat14, mat15, mat16, &
                                                 mat17, mat18, mat19, mat110, mat111, &
                                                 mat21, mat22, mat23, mat24, mat25, mat26, &
                                                 mat27, mat28, mat29, mat210, mat211

            type(IterDate) :: dateC, date1111
            character(len=80) :: line
            integer :: i

            call readHeader_MESH_OUTPUTFILES(filepath, filename, fileunit, ios)

            call parseDateTime(inDate, dateC)
            ! NOTE: MESH gridded output files go from 1H to 24H. Our program begin 0H and end at 23H.
            ! 1 is added below to adjust our hour accoridingly.
            dateC%hour = dateC%hour + 1
            !secsC = dateC%hour*3600 + dateC%mins*60 + dateC%secs

            do while(ios.eq.0)

                read(unit=fileunit, FMT='((A))', iostat=ios) line
                date11 = trim(line(index(line, '"'):len_trim(line)))
                call parseDateTime(date11, date1111)
                do i = 1, yCount
                    read(unit = fileunit, fmt = *, iostat=ios) mat11(i,:)
                end do
                read(unit = fileunit, fmt = *, iostat=ios)

                print *, date11
                if (date1111%year == dateC%year) then
                    if (date1111%month == dateC%month) then
                        if (date1111%day == dateC%day) then
                            if (date1111%hour == dateC%hour) then

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date12 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat12(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date13 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat13(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date14 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat14(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date15 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat15(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date16 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat16(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date17 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat17(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date18 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat18(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date19 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat19(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date110 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat110(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                read(unit=fileunit, FMT='((A))', iostat=ios) line
                                date111 = trim(line(index(line, '"'):len_trim(line)))
                                do i = 1, yCount
                                    read(unit = fileunit, fmt = *, iostat=ios) mat111(i,:)
                                end do
                                read(unit = fileunit, fmt = *, iostat=ios)

                                exit
                            end if
                        end if
                    end if
                end if
            end do

            call getMat2_sed_input(fileunit, mat21, date21, mat22, date22, mat23, date23, &
                                     mat24, date24, mat25, date25, mat26, date26, &
                                     mat27, date27, mat28, date28, mat29, date29, mat210, date210, &
                                     mat211, date211, ios)

        end subroutine getMat1_sed_input


        subroutine getMat2_rbm_input(fileunit, mat21, date21, mat22, date22, mat23, date23, &
                                     mat24, date24, mat25, date25, mat26, date26, ios)
            implicit none
            integer, intent(in) :: fileunit
            character(len=80), intent(out) :: date21, date22, date23, date24, date25, date26
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: mat21, mat22, mat23, mat24, mat25, mat26

            character(len=80) :: line

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date21 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat21(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date22 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat22(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date23 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat23(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date24 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat24(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date25 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat25(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date26 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat26(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

        end subroutine getMat2_rbm_input

        subroutine getMat2_sed_input(fileunit, mat21, date21, mat22, date22, mat23, date23, &
                                     mat24, date24, mat25, date25, mat26, date26, &
                                     mat27, date27, mat28, date28, mat29, date29, mat210, date210, &
                                     mat211, date211, ios)
            implicit none
            integer, intent(in) :: fileunit
            character(len=80), intent(out) :: date21, date22, date23, date24, date25, date26, &
                                              date27, date28, date29, date210, date211
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: mat21, mat22, mat23, mat24, mat25, mat26, &
                                                 mat27, mat28, mat29, mat210, mat211

            character(len=80) :: line

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date21 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat21(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date22 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat22(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date23 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat23(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date24 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat24(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date25 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat25(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date26 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat26(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date27 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat27(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date28 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat28(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date29 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat29(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date210 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat210(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

            read(unit=fileunit, FMT='((A))', iostat=ios) line
            date211 = trim(line(index(line, '"'):len_trim(line)))
            do i = 1, yCount
                read(unit = fileunit, fmt = *, iostat=ios) mat211(i,:)
            end do
            read(unit = fileunit, fmt = *, iostat=ios)

        end subroutine getMat2_sed_input


        subroutine interp2_rbm_input(fileunit, mat11, date11, mat12, date12, &
                                     mat13, date13, mat14, date14, mat15, date15, mat16, date16, &
                                     mat21, date21, mat22, date22, mat23, date23, mat24, date24, &
                                     mat25, date25, mat26, date26, matF1, matF2, matF3, matF4, matF5, &
                                     matF6, inDate, ios)
            implicit none

            integer, intent(in) :: fileunit
            character(len=80), intent(inout) :: date11, date12, date13, date14, date15, date16, &
                                                date21, date22, date23, date24, date25, date26
            real, dimension(:,:), intent(inout) :: mat11, mat12, mat13, mat14, mat15, mat16, &
                                                   mat21, mat22, mat23, mat24, mat25, mat26
            character(len=80), intent(in) :: inDate
            integer, intent(inout) :: ios
            real, dimension(:,:), intent(out) :: matF1, matF2, matF3, matF4, matF5, matF6


            !real, dimension(yCount,xCount) :: matF
            type(IterDate) :: dateC, date222
            real :: secsC
            !real :: dsecs11

            call parseDateTime(inDate, dateC)
            secsC = dateC%mins*60 + dateC%secs

            call parseDateTime(date21, date222)

            if (dateC%year == date222%year) then
                if (dateC%month == date222%month) then
                    if (dateC%day == date222%day) then
                        if (dateC%hour+1 == date222%hour) then
                            if (secsC >= 0.0) then
                                date11 = date21
                                mat11 = mat21

                                date12 = date22
                                mat12 = mat22

                                date13 = date23
                                mat13 = mat23

                                date14 = date24
                                mat14 = mat24

                                date15 = date25
                                mat15 = mat25

                                date16 = date26
                                mat16 = mat26

                                call getMat2_rbm_input(fileunit, mat21, date21, mat22, date22, mat23, date23, &
                                     mat24, date24, mat25, date25, mat26, date26, ios)
                            end if
                        end if
                    end if
                end if
            end if

            !dsecs11 = secsC - secs11
            !matF = (mat2-mat1)*(dsecs11/3600) + mat1
            matF1 = (mat21-mat11)*(secsC/3600) + mat11
            matF2 = (mat22-mat12)*(secsC/3600) + mat12
            matF3 = (mat23-mat13)*(secsC/3600) + mat13
            matF4 = (mat24-mat14)*(secsC/3600) + mat14
            matF5 = (mat25-mat15)*(secsC/3600) + mat15
            matF6 = (mat26-mat16)*(secsC/3600) + mat16

        end subroutine interp2_rbm_input


        subroutine interp2_sed_input(mat11, date11, mat12, date12, &
                                     mat13, date13, mat14, date14, mat15, date15, mat16, date16, &
                                     mat17, date17, mat18, date18, mat19, date19, mat110, date110, mat111, date111, &
                                     mat21, date21, mat22, date22, mat23, date23, mat24, date24, mat25, date25, &
                                     mat26, date26, mat27, date27, mat28, date28, mat29, date29, mat210, date210, &
                                     mat211, date211, inDate)
            implicit none

            !integer, intent(in) :: fileunit
            character(len=80), intent(inout) :: date11, date12, date13, date14, date15, date16, &
                                                date17, date18, date19, date110, date111, &
                                                date21, date22, date23, date24, date25, date26, &
                                                date27, date28, date29, date210, date211

            real, dimension(:,:), intent(inout) :: mat11, mat12, mat13, mat14, mat15, mat16, &
                                                   mat17, mat18, mat19, mat110, mat111, &
                                                   mat21, mat22, mat23, mat24, mat25, mat26, &
                                                   mat27, mat28, mat29, mat210, mat211

            character(len=80), intent(in) :: inDate
            !integer, intent(inout) :: ios
            !real, dimension(:,:), intent(out) :: matF1, matF2, matF3, matF4, matF5, matF6, &
            !                                     matF7, matF8, matF9, matF10, matF11


            !real, dimension(yCount,xCount) :: matF
            type(IterDate) :: dateC, date2222
            real :: secsC
            !real :: dsecs11

            call parseDateTime(inDate, dateC)
            secsC = dateC%mins*60 + dateC%secs

            call parseDateTime(date21, date2222)

            if (dateC%year == date2222%year) then
                if (dateC%month == date2222%month) then
                    if (dateC%day == date2222%day) then
                        if (dateC%hour+1 == date2222%hour) then
                            if (secsC >= 0.0) then
                                date11 = date21
                                mat11 = mat21

                                date12 = date22
                                mat12 = mat22

                                date13 = date23
                                mat13 = mat23

                                date14 = date24
                                mat14 = mat24

                                date15 = date25
                                mat15 = mat25

                                date16 = date26
                                mat16 = mat26

                                date17 = date27
                                mat17 = mat27

                                date18 = date28
                                mat18 = mat28

                                date19 = date29
                                mat19 = mat29

                                date110 = date210
                                mat110 = mat210

                                date111 = date211
                                mat111 = mat211

                                call getMat2_sed_input(sedi%unitSED, mat21, date21, mat22, date22, mat23, date23, &
                                     mat24, date24, mat25, date25, mat26, date26, &
                                     mat27, date27, mat28, date28, mat29, date29, mat210, date210, &
                                     mat211, date211, sedi%iosSED)
                            end if
                        end if
                    end if
                end if
            end if

            !dsecs11 = secsC - secs11
            !matF = (mat2-mat1)*(dsecs11/3600) + mat1
            sedi%DISC = (mat21-mat11)*(secsC/3600) + mat11
            sedi%DEPT = (mat22-mat12)*(secsC/3600) + mat12
            sedi%WIDT = (mat23-mat13)*(secsC/3600) + mat13
            sedi%CHLE = (mat24-mat14)*(secsC/3600) + mat14
            sedi%CHSL = (mat25-mat15)*(secsC/3600) + mat15
            sedi%VELO = (mat26-mat16)*(secsC/3600) + mat16
            sedi%PREP = (mat27-mat17)*(secsC/3600) + mat17
            sedi%EVAP = (mat28-mat18)*(secsC/3600) + mat18
            sedi%OFDE = (mat29-mat19)*(secsC/3600) + mat19
            sedi%LASL = (mat210-mat110)*(secsC/3600) + mat110
            sedi%CEWI = (mat211-mat111)*(secsC/3600) + mat111

        end subroutine interp2_sed_input



!        subroutine intep2_MESH_OUTFIELD(fileunit, mat1, date1, matF, inDate, ios)
!            implicit none
!            integer, intent(in) :: fileunit
!            real, dimension(:,:), intent(out) :: matF
!            real, dimension(:,:), intent(inout) :: mat1
!            character(len=80), intent(in) :: inDate
!            character(len=80), intent(inout) :: date1
!            integer :: ios
!
!            real, dimension(yCount,xCount) :: mat2
!            character(len=80) :: line, date2
!            type(IterDate) :: dateC, date11, date22
!            real :: secsC, secs11
!
!            call parseDateTime(inDate, dateC)
!            ! NOTE: MESH gridded output files go from 1H to 24H. Our program begin 0H and end at 23H.
!            ! 1 is added below to adjust our hour accoridingly.
!            dateC%hour = dateC%hour + 1
!            secsC = dateC%mins*60 + dateC%secs
!
!            call parseDateTime(date1, date11)
!            secs11 = date11%mins*60 + date11%secs
!
!            read(unit=fileunit, FMT='((A))', iostat=ios) line
!            date2 = trim(line(index(line, '"'):len_trim(line)))
!            call parseDateTime(date2, date22)
!            do i = 1, yCount
!                read(unit = fileunit, fmt = *, iostat=ios) mat2(i,:)
!            end do
!            read(unit = fileunit, fmt = *, iostat=ios)
!
!            dsecs11 = secsC - secs11
!            matF = (mat2-mat1)*(dsecs11/3600) + mat1
!            date1 = date2
!            mat1 = mat2
!
!        end subroutine intep2_MESH_OUTFIELD


!        subroutine intep2_MESH_OUTFIELD(fileunit, mat1, date1, matF, inDate, ios)
!            implicit none
!            integer, intent(in) :: fileunit
!            real, dimension(:,:), intent(out) :: matF
!            real, dimension(:,:), intent(inout) :: mat1
!            character(len=80), intent(in) :: inDate
!            character(len=80), intent(inout) :: date1
!            integer :: ios
!
!            real, dimension(yCount,xCount) :: mat2
!            character(len=80) :: line, date2
!            type(IterDate) :: dateC, date11, date22
!            real :: secsC, secs11, dsecs11
!
!
!            call parseDateTime(inDate, dateC)
!            ! NOTE: MESH gridded output files go from 1H to 24H. Our program begin 0H and end at 23H.
!            ! 1 is added below to adjust our hour accoridingly.
!            dateC%hour = dateC%hour + 1
!            secsC = dateC%mins*60 + dateC%secs
!
!            call parseDateTime(date1, date11)
!            secs11 = date11%mins*60 + date11%secs
!
!            do while(ios.eq.0)
!
!                read(unit=fileunit, FMT='((A))', iostat=ios) line
!                date2 = trim(line(index(line, '"'):len_trim(line)))
!                call parseDateTime(date2, date22)
!                do i = 1, yCount
!                    read(unit = fileunit, fmt = *, iostat=ios) mat2(i,:)
!                end do
!                read(unit = fileunit, fmt = *, iostat=ios)
!
!
!!                if (date22%year >= dateC%year) then
!!                    if (date22%month <= dateC%month) then
!!                        if (date22%day <= dateC%day) then
!
!                if (date11%year == dateC%year .or. date22%year == dateC%year) then
!                    if (date11%month == dateC%month .or. date22%month == dateC%month) then
!                        if (date11%day == dateC%day .or. date22%day == dateC%day) then
!                            !secs22 = date22%hour*3600 + date22%mins*60 + date22%secs
!                            if (dateC%hour + 1 == date22%hour + 1) then
!
!                                dsecs11 = secsC - secs11
!                                matF = (mat2-mat1)*(dsecs11/3600) + mat1
!
!                                exit
!                            end if
!!                            if (dateC%hour == 0)then
!!                                dsecs11 = secsC
!!                                matF = (mat2-mat1)*(dsecs11/3600) + mat1
!!                                !print *,'date1 ', date1
!!                                !print *,'date2 ', date2
!!                                exit
!!                            end if
!                        end if
!                    end if
!                end if
!                date1 = date2
!                mat1 = mat2
!            end do
!        end subroutine intep2_MESH_OUTFIELD


        subroutine close_MESH_OUTFIELD(fileunit)
            implicit none
            integer, intent(in) :: fileunit
            close(unit = fileunit)
        end subroutine close_MESH_OUTFIELD


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
