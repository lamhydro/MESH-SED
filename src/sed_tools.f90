module sed_tools
    implicit none
    contains

        real function fracDiame(partiDistr, partiDiame, frac)
            implicit none

            real, dimension(:), intent(in) :: partiDistr, partiDiame
            real, intent(in) :: frac

            real, dimension(:), allocatable :: accPartiDistr
            integer :: i, nparti


            nparti = size(partiDistr)
            allocate(accPartiDistr(nparti))

            ! Accumulated particle distribution
            accPartiDistr(1) = partiDistr(1)
            do i = 2, nparti
                accPartiDistr(i) = accPartiDistr(i-1) + partiDistr(i)
            end do

            ! Interpolation of frac in the accum. parti. distri.
            fracDiame = linearInterpo(accPartiDistr,partiDiame,frac)

        end function fracDiame

        real function linearInterpo(x,y,xi)
            implicit none
            real, dimension(:), intent(in) :: x,y
            real, intent(in) :: xi

            integer :: n,i
            real :: dx, dy

            n = size(x)

            linearInterpo = 0.

            ! Check the limits for "xi"
            if (xi < x(1)) then
                !print*, 'Warning!, extrapolation, low limit'
                linearInterpo = y(1)
            else if (xi > x(n)) then
                !print*, 'Warning!, extrapolation, high limit'
                linearInterpo = y(n)
                !print*, xi, " is > 1, nothing to do."
                !stop
            else
                ! Lineal interpolation
                do i=1,n-1
                  if (x(i) <= xi .and. x(i+1) >= xi)  then
                   dx = x(i+1)-x(i)
                   dy = y(i+1)-y(i)
                   if (dx == 0) then
                    linearInterpo = y(i)
                   else
                    linearInterpo = y(i) + dy*(xi - x(i))/dx
                   end if
                   exit
                  end if
                end do
            end if

        end function linearInterpo

        ! split a string into 2 either side of a delimiter token
        SUBROUTINE split_string(instring, string1, string2, delim)
            CHARACTER(80) :: instring,delim
            CHARACTER(80),INTENT(OUT):: string1,string2
            INTEGER :: index

            instring = TRIM(instring)

            index = SCAN(instring,delim)
            string1 = instring(1:index-1)
            string2 = instring(index+1:)

        END SUBROUTINE split_string

        !> *****************************************************************
        !> Description: Convert Julian day to month and day in gregorian
        !> calendar given the Julian day and the year
        !> *****************************************************************
        subroutine Julian2MonthDay(jday,year,month,day)

            !> Input
            integer, intent(in) :: jday, year

            !> Output
            integer, intent(out) :: month, day

            !> Internals
            integer, parameter, dimension(12) :: &
                daysnl = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365], &
                daysyl = [31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]

            integer i, int_i, int_f

            do i = 2, 12

                if (leap_year(year) == 365) then
                    int_i = daysnl(i - 1)
                    int_f = daysnl(i)
                elseif (leap_year(year) == 366) then
                    int_i = daysyl(i - 1)
                    int_f = daysyl(i)
                end if

                if (jday <= 31) then
                    month = 1
                    day = jday
                    exit
                else
                    if ((jday > int_i) .and. (jday <= int_f)) then
                        month = i
                        day = jday - int_i
                        exit
                    elseif (jday == int_i) then
                        month = i - 1
                        if (leap_year(year) == 365) then
                            day = jday - daysnl(i - 1)
                            exit
                        elseif (leap_year(year) == 366) then
                            day = jday - daysyl(i - 1)
                            exit
                        end if
                    end if
                end if
            end do

        end subroutine !Julian2MonthDay

        !> *****************************************************************
        !> Description: Get the number of days in leap and normal years
        !> (365 or 366).
        !> *****************************************************************
        integer function leap_year(y) result(ndays)

            logical is_leap
            integer, intent(in) :: y

            is_leap = (mod(y, 4) == 0 .and. .not. mod(y, 100) == 0) .or. (mod(y, 400) == 0)

            if (is_leap) then
                ndays = 366
            else
                ndays = 365
            end if

        end function !leap_year


        !> Open a text file
        subroutine openFile(fileunit, filepath, filename)
            implicit none

            integer :: ios

            integer, intent(in) :: fileunit
            character(len=80), intent(in) :: filename
            character(len=150), intent(in) :: filepath

            !print *, trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/" // &
            !     trim(ADJUSTL(metVar)) // ".r2c"

            ! Open the file
            !print *,  trim(ADJUSTL(filepath)) // trim(ADJUSTL(filename))
            open(unit = fileunit, &
                !file=trim(ADJUSTL(MESHdir)) // trim(ADJUSTL(OUTFIELDfolder)) // "/" // &
                 !trim(ADJUSTL(metVar)) //  ".r2c", &
                file =  trim(ADJUSTL(filepath)) // trim(ADJUSTL(filename)), &
                status="replace", &
                action="readwrite", &
                iostat=ios)

            if (ios /= 0) then
                close(unit = fileunit)
                print *, "The file ", trim(ADJUSTL(filename)), " could not be opened."
                print *, "Stopping..."
                stop

            else
                !print *, 'OPENING: ', trim(ADJUSTL(filename))
            end if

        end subroutine openFile


        !> Close a text file
        subroutine closeFile(fileunit)
            implicit none

            integer, intent(in) :: fileunit

            ! Close the file
            close(unit = fileunit)

        end subroutine closeFile



end module sed_tools
