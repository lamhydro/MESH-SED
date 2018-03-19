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



        subroutine writeGrPointTS()
            implicit none
            integer :: fileunit
            real :: variVal

            do i =1,opin%nGrPoint
                fileunit = i+1000
                if (trim(ogrp(i)%varname) == 'load') then
                    variVal = isrr(ogrp(i)%grNumb)%G(ogrp(i)%idsedclass)*cbsca(i)%density ! multiply by density o get kg/s
                else ! for concentration
                    variVal = isrr(ogrp(i)%grNumb)%C(ogrp(i)%idsedclass)*cbsca(i)%density * 1000. ! multiply by density and 1000 to get mg/l
                end if
                !print *, dateIn, variVal
                write(fileunit, FMT=*) trim(dateIn), variVal
            end do
        end subroutine writeGrPointTS



        subroutine closeFilesForGrPointTS()
            implicit none
            integer :: fileunit

            do i =1,opin%nGrPoint
                fileunit = i+1000
                call closeFile(fileunit)
            end do

        end subroutine closeFilesForGrPointTS
end module sed_output
