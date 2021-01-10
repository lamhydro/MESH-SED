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

        subroutine openFilesForGrPointTS_2()
            implicit none

            integer :: fileunit
            character(len=80):: filename, grNumbStr

            do i =1,opin%nGrPoint
                fileunit = i+1000
                write (grNumbStr,"(I10)") ogrp(i)%grNumb
                filename = 'ts_'//trim(adjustl(grNumbStr)) //  '.out'
                !print *, filename
                call openFile(fileunit, opin%outputDir, filename)

                write(fileunit, FMT='(1X, A, 10(",",A))') 'DateTime','L_in(kg/s)','L_bank(kg/s)','L_hill(kg/s)', &
                                                            'L_res(kg/s)','L_dep(kg/s)','L_out(kg/s)','C(kg/m3)','C_pot(kg/m3)', &
                                                            'Q(m3/s)','P(mm/h)'




            end do
        end subroutine openFilesForGrPointTS_2


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


!        subroutine writeGrPointTS_2()
!            implicit none
!            integer :: fileunit, nn
!            !real :: variVal
!
!            do i =1,opin%nGrPoint
!                fileunit = i+1000
!                !if (trim(ogrp(i)%varname) == 'load') then
!                !variVal = isrr(ogrp(i)%grNumb)%G(ogrp(i)%idsedclass)*cbsca(ogrp(i)%grNumb)%density ! multiply by density o get kg/s
!                !else ! for concentration
!                !    variVal = isrr(ogrp(i)%grNumb)%C(ogrp(i)%idsedclass)*cbsca(ogrp(i)%grNumb)%density * 1000. ! multiply by density and 1000 to get mg/l
!                !end if
!                !print *, dateIn, variVal
!                nn = ogrp(i)%grNumb
!                !write(fileunit, FMT=*) trim(dateIn), sum(cmb(nn)%L_bank(:)),  sum(cmb(nn)%L_hill(:)), &
!                !    sum(cmb(nn)%L_in(:)), rh(nn)%depth, sum(cmb(nn)%C(:)*cbsca(nn)%frac(:)), &
!                !    sum(cmb(nn)%C_pot(:)*cbsca(nn)%frac(:)), rh(nn)%discharge
!
!                write(fileunit, FMT='(1X, A, 7(",",1X, F16.7))') trim(dateIn), sum(cmb(nn)%L_bank(:)),  sum(cmb(nn)%L_hill(:)), &
!                    sum(cmb(nn)%L_in(:)), rh(nn)%depth, sum(cmb(nn)%C(:)*cbsca(nn)%frac(:)), &
!                    sum(cmb(nn)%C_pot(:)*cbsca(nn)%frac(:)), rh(nn)%discharge
!
!            end do
!        end subroutine writeGrPointTS_2

        subroutine writeGrPointTS_3()
            implicit none
            integer :: fileunit, nn
            !real :: variVal

            !print *, niter
            do i =1,opin%nGrPoint
                fileunit = i+1000
                !if (trim(ogrp(i)%varname) == 'load') then
                !variVal = isrr(ogrp(i)%grNumb)%G(ogrp(i)%idsedclass)*cbsca(ogrp(i)%grNumb)%density ! multiply by density o get kg/s
                !else ! for concentration
                !    variVal = isrr(ogrp(i)%grNumb)%C(ogrp(i)%idsedclass)*cbsca(ogrp(i)%grNumb)%density * 1000. ! multiply by density and 1000 to get mg/l
                !end if
                !print *, dateIn, variVal
                nn = ogrp(i)%grNumb
                !write(fileunit, FMT=*) trim(dateIn), sum(cmb(nn)%L_bank(:)),  sum(cmb(nn)%L_hill(:)), &
                !    sum(cmb(nn)%L_in(:)), rh(nn)%depth, sum(cmb(nn)%C(:)*cbsca(nn)%frac(:)), &
                !    sum(cmb(nn)%C_pot(:)*cbsca(nn)%frac(:)), rh(nn)%discharge

!                write(fileunit, FMT='(1X, A, 7(",",1X, F16.7))') trim(dateIn), &
!                                                                 cmb_h(nn)%L_bank/niter, &
!                                                                 cmb_h(nn)%L_hill/niter, &
!                                                                 cmb_h(nn)%L_in/niter, &
!                                                                 cmb_h(nn)%C/niter, &
!                                                                 cmb_h(nn)%C_pot/niter, &
!                                                                 cmb_h(nn)%rhQ/niter, &
!                                                                 cmb_h(nn)%rhH/niter

                write(fileunit, FMT='(1X, A, 10(",",F16.7))') trim(dateIn), &
                                                                 cmb_h(nn)%L_in/niter, &
                                                                 cmb_h(nn)%L_bank/niter, &
                                                                 cmb_h(nn)%L_hill/niter, &
                                                                 cmb_h(nn)%L_res/niter, &
                                                                 cmb_h(nn)%L_dep/niter, &
                                                                 cmb_h(nn)%L_out/niter, &
                                                                 cmb_h(nn)%C/niter, &
                                                                 cmb_h(nn)%C_pot/niter, &
                                                                 cmb_h(nn)%rhQ/niter, &
                                                                 mv(nn)%precip



            end do
        end subroutine writeGrPointTS_3

        subroutine writeHourAveEstimate()
            implicit none

            do i = 1, NA
                cmb_h(i)%L_in    = cmb_h(i)%L_in + sum(cmb(i)%L_in(:))
                cmb_h(i)%L_bank  = cmb_h(i)%L_bank + sum(cmb(i)%L_bank(:))
                cmb_h(i)%L_hill  = cmb_h(i)%L_hill + sum(cmb(i)%L_hill(:))
                cmb_h(i)%L_res   = cmb_h(i)%L_res   + sum(cmb(i)%L_res(:))
                cmb_h(i)%L_dep   = cmb_h(i)%L_dep   + sum(cmb(i)%L_dep(:))
                cmb_h(i)%L_out   = cmb_h(i)%L_out   + sum(cmb(i)%L_out(:))

                cmb_h(i)%C     = cmb_h(i)%C + sum(cmb(i)%C(:)*cbsca(i)%frac(:))
                cmb_h(i)%C_pot = cmb_h(i)%C_pot + sum(cmb(i)%C_pot(:)*cbsca(i)%frac(:))

                cmb_h(i)%rhQ = cmb_h(i)%rhQ + rh(i)%discharge
                cmb_h(i)%rhH = cmb_h(i)%rhH + rh(i)%depth
            end do
            niter = niter + 1

            if (currDate%mins == 0) then
                call writeGrPointTS_3
                cmb_h(:)%L_in = 0.; cmb_h(:)%L_bank = 0.; cmb_h(:)%L_hill = 0.; &
                cmb_h(:)%L_res = 0.; cmb_h(:)%L_dep = 0.; cmb_h(:)%L_out = 0.;&
                cmb_h(:)%C = 0.; cmb_h(:)%C_pot = 0.; cmb_h(:)%rhQ = 0.; cmb_h(:)%rhH = 0.

                niter = 0
            end if
        end subroutine writeHourAveEstimate

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
