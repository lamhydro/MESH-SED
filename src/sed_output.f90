module sed_output
    use sed_vars
    use sed_tools

    implicit none

    contains


        subroutine openFilesForGrPointTS()
            implicit none

            integer :: fileunit
            character(len=80):: filename, grNumbStr

            do i =1,opin%nGrPoint
                fileunit = i+1000
                write (grNumbStr,"(I10)") ogrp(i)%grNumb
                filename = 'ts_'//trim(adjustl(grNumbStr)) // '_' // trim(ogrp(i)%varname) &
                            // '_' // trim(ogrp(i)%sedname) // '.out'
                !print *, filename
                call openFile(fileunit, opin%outputDir, filename)
            end do
        end subroutine openFilesForGrPointTS


        subroutine openFilesForFieldTS()
            implicit none

            integer :: fileunit
            character(len=80):: filename

            !> Initialize output matrix
            do k=1,yCount
                do j=1,xCount
                    variMat(k,j) = -99999.9
                end do
            end do

            do i =1,ofin%nField
                fileunit = i+500
                filename = 'field_'// trim(ofie(i)%varname) &
                            // '_' // trim(ofie(i)%sedname) // '.out'
                !print *, filename
                call openFile(fileunit, ofin%outputDir, filename)
            end do
        end subroutine openFilesForFieldTS


        subroutine writeGrPointTS()
            implicit none
            integer :: fileunit
            real :: variVal

            do i =1,opin%nGrPoint
                fileunit = i+1000
                if (trim(ogrp(i)%varname) == 'load') then
                    variVal = isrr(ogrp(i)%grNumb)%G(ogrp(i)%idsedclass)*cbsca(ogrp(i)%grNumb)%density ! multiply by density o get kg/s
                else ! for concentration
                    variVal = isrr(ogrp(i)%grNumb)%C(ogrp(i)%idsedclass)*cbsca(ogrp(i)%grNumb)%density * 1000. ! multiply by density and 1000 to get mg/l
                end if
                !print *, dateIn, variVal
                write(fileunit, FMT=*) trim(dateIn), variVal
            end do
        end subroutine writeGrPointTS

        subroutine writeFieldTS()
            implicit none
            integer :: fileunit,k

            do i = 1, ofin%nField
                fileunit = i+500
                write(fileunit, FMT=*) trim(dateIn)
                do k =1, NA
                    if (trim(ofie(i)%varname) == 'load') then
                        variMat(ipos(k), jpos(k)) = isrr(k)%G(ofie(i)%idsedclass)*cbsca(k)%density ! multiply by density o get kg/s
                    else ! for concentration
                        variMat(ipos(k), jpos(k)) = isrr(k)%C(ofie(i)%idsedclass)*cbsca(k)%density * 1000. ! multiply by density and 1000 to get mg/l
                    end if
                end do

                do k = 1,yCount
                    write(fileunit, '(*(F14.4))') (variMat(k,j) ,j=1,xCount)
                    !write(fileunit, *) (variMat(k,j), j=1,xCount)
                    !do j = 1, xCount
                    !    write(fileunit, FMT=*) trim(dateIn), variVal
                    !end do
                end do
            end do
        end subroutine writeFieldTS



        subroutine closeFilesForGrPointTS()
            implicit none
            integer :: fileunit

            do i =1,opin%nGrPoint
                fileunit = i+1000
                call closeFile(fileunit)
            end do
        end subroutine closeFilesForGrPointTS

        subroutine closeFilesForFieldTS()
            implicit none
            integer :: fileunit

            do i =1,ofin%nField
                fileunit = i+500
                call closeFile(fileunit)
            end do
        end subroutine closeFilesForFieldTS


end module sed_output
