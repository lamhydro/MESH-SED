module sed_inStreamTransCapa
    use sed_overlandFlowDetachment
    implicit none

    contains

        function inStreamTransCapa_AckersWhite(rho, rhos, D, gravi, h, S, v, Vel, Q)
        ! Transport capacity of total loads (bed and suspended loads) by Ackers and White 1973
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - Vel: Water velocity (m/s)
        ! - h: Stream water depth (m)
        ! - S: Water surface slope
        ! - D: Sediment diameter (m). Authors advise D_35 size
        ! - rhos: Sediment density (kg/m^3)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - Q: Discharge (m^3 /s)
        ! Output:
        ! - inStreamTransCapa_AckersWhite: Volumetric sediment transport rate (m^3 s^-1)

            implicit none

            real, intent(in) :: rho, rhos, D, gravi, h, S, v, Vel, Q
            real :: inStreamTransCapa_AckersWhite, sg, Dgr, n, A, M, c, Vs, Fgr, Ggr

            ! Check Froude Number first before applying this procedure
    !        Fr = Vel/sqrt(gravi*h)
    !        if (Fr >= 0.8) then
    !            print *, "Froude number is >= 0.8. Ackers and White applies when Fr < 0.8"
    !            print *, "Change the methodology. Stopping..."
    !            stop
    !        end if

            ! Specific gravity
            sg = rhos/rho

            ! Dimensionless sediment diameter
            Dgr = D*( gravi*(sg-1)/(v**2) )**(1/3)

            ! Transition exponent (n), initial motion parameter (A), and
            ! coeff. (c) and exponent (M) in the sed. transport function
            if (Dgr > 60) then
                n = 0.
                A = 0.17
                M = 1.5
                c = 0.025
            else if(Dgr <= 60 .and. Dgr >= 1 ) then
                n = 1 - 0.56*log10(Dgr)
                A = 0.14 + 0.23/sqrt(Dgr)
                M = 9.66/Dgr + 1.34
                c = 10**(2.86*log10(Dgr) - (log10(Dgr))**2 - 3.53)
            else !else if(Dgr < 1) then
                n = 1
                A = 0.37
                M = 11
                c = 2.95e-4
                !print *, "Dgr is < 1. Ackers and White applies when Dgr >= 1."
                !print *, "Change the methodology. Stopping..."
                !stop

            end if

            ! Particle mobility
            ! - Shear velocity (m/s)
            Vs = shearVelocity(rho, gravi, h, S)

            Fgr = ( (Vs**n)/(sqrt(gravi*D*(sg-1))) )*( Vel /(sqrt(32.)*log10(10*h/D)) )**(1-n)

            ! Dimensionless sediment transport rate
            if (A < Fgr) then
                Ggr = c*((Fgr/A - 1)**M)
            else
                Ggr = 0.
            end if

            ! Volumetric sediment transport rate G (m^3 /s )
            inStreamTransCapa_AckersWhite = Q*Ggr*D*((Vel/Vs)**n)/h

        end function inStreamTransCapa_AckersWhite


        subroutine inStreamTransCapa_Day(rho, rhos, Dm, frac, D16, D50, D84, gravi, h, S, v, Vel, Q, Gm)
        ! Transport capacity of total loads (bed and suspended loads) by Day 1980 (a modifification of Ackers and White)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - rhos: Sediment density (kg/m^3)
        ! - Dm: Vector of sediment diameters (m)
        ! - frac: fractions ocupied by each sediment diameter in the sediment grain distribution.
        ! - D16: D_16 sediment diameter (m)
        ! - D50: D_50 sediment diameter (m)
        ! - D84: D_84 sediment diameter (m)
        ! - gravi: Gravity acce. (m/s^2)
        ! - Vel: Water velocity (m/s)
        ! - h: Water depth (m)
        ! - S: Water surface slope
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - Q: Discharge (m^3 /s)
        ! Output:
        ! - Gm: Vector of volumetric sediment transport rates for different sediment sizes (m^3 s^-1)

            implicit none

            real, intent(in) :: rho, rhos, D16, D50, D84, gravi, h, S, v, Vel, Q
            real, dimension(:), intent(in) :: Dm, frac
            real, dimension(:), intent(out) :: Gm

            real :: sg, DA, DgrA, A, Am, Dgrm,  n, M, c, Vs, Fgrm, Ggrm
            integer :: nm, i

            nm = size(Dm)

            ! - Shear velocity (m/s)
            Vs = shearVelocity(rho, gravi, h, S)

            ! Specific gravity
            sg = rhos/rho

            ! Critical diameter
            DA = 1.62*D50*((D84/D16)**(-0.28))

            ! Dimensionless grainsize for DA
            DgrA = DA * (gravi*(sg-1)/(v**2))**(1/3)

            ! Initial motion parameter for DA
            if (DgrA <= 60. .and. DgrA >= 1.) then
                A = 0.23/sqrt(DgrA) + 0.14
            else
                A = 0.17
            end if


            do i = 1, nm
                if (Dm(i) > 0.062*1.e-3) then   ! For coarse sediment
                    ! Initial motion parameter for Dm
                    Am = A*(0.4*((Dm(i)/DA)**(-0.5)) + 0.6)

                    ! Dimensionless sediment diameter for Dm
                    Dgrm = Dm(i)*( gravi*(sg-1)/(v**2) )**(1/3)

                    ! Transition exponent (n), and
                    ! coeff. (c) and exponent (M) in the sed. transport function
                    if (Dgrm > 60) then
                        n = 0.
                        M = 1.5
                        c = 0.025
                    else if(Dgrm <= 60 .and. Dgrm >= 1 ) then
                        n = 1 - 0.56*log10(Dgrm)
                        M = 9.66/Dgrm + 1.34
                        c = 10**(2.86*log10(Dgrm) - (log10(Dgrm))**2 - 3.53)
                    else !else if(Dgrm < 1) then
                        n = 1
                        M = 11
                        c = 2.95e-4
                        !print *, "Dgr is < 1. Ackers and White applies when Dgr >= 1."
                        !print *, "Change the methodology. Stopping..."
                        !stop

                    end if

                    ! Particle mobility
                    Fgrm = ( (Vs**n)/(sqrt(gravi*Dm(i)*(sg-1))) )*( Vel /(sqrt(32.)*log10(10*h/Dm(i))) )**(1-n)

                    ! Dimensionless sediment transport rate
                    if (Am < Fgrm) then
                        Ggrm = c*((Fgrm/Am - 1)**M)
                    else
                        Ggrm = 0.
                    end if

                    ! Volumetric sediment transport rate G (m^3 /s )
                    Gm(i) = frac(i)*Q*Ggrm*Dm(i)*((Vel/Vs)**n)/h
                end if

            end do

        end subroutine inStreamTransCapa_Day

        real function inStreamTransCapa_Day_1(rho, rhos, Dm, D16, D50, D84, gravi, h, S, v, Vel, Q)
        ! Transport capacity of total loads (bed and suspended loads) by Day 1980 (a modifification of Ackers and White)
        ! NOTE: This is the subroutine converted into function for one diameter once a time.
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - rhos: Sediment density (kg/m^3)
        ! - Dm: Sediment diameter (m)
        ! - D16: D_16 sediment diameter (m)
        ! - D50: D_50 sediment diameter (m)
        ! - D84: D_84 sediment diameter (m)
        ! - gravi: Gravity acce. (m/s^2)
        ! - Vel: Water velocity (m/s)
        ! - h: Water depth (m)
        ! - S: Water surface slope
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - Q: Discharge (m^3 /s)
        ! Output:
        ! - inStreamTransCapa_Day_1: Volumetric sediment transport rates for different sediment sizes (m^3 s^-1)

            implicit none

            real, intent(in) :: rho, rhos, Dm, D16, D50, D84, gravi, h, S, v, Vel, Q

            real :: sg, DA, DgrA, A, Am, Dgrm,  n, M, c, Vs, Fgrm, Ggrm

            ! - Shear velocity (m/s)
            Vs = shearVelocity(rho, gravi, h, S)

            ! Specific gravity
            sg = rhos/rho

            ! Critical diameter
            DA = 1.62*D50*((D84/D16)**(-0.28))

            ! Dimensionless grainsize for DA
            DgrA = DA * (gravi*(sg-1)/(v**2))**(1/3)

            ! Initial motion parameter for DA
            if (DgrA <= 60. .and. DgrA >= 1.) then
                A = 0.23/sqrt(DgrA) + 0.14
            else
                A = 0.17
            end if

            ! Initial motion parameter for Dm
            Am = A*(0.4*((Dm/DA)**(-0.5)) + 0.6)

            ! Dimensionless sediment diameter for Dm
            Dgrm = Dm*( gravi*(sg-1)/(v**2) )**(1/3)

            ! Transition exponent (n), and
            ! coeff. (c) and exponent (M) in the sed. transport function
            if (Dgrm > 60) then
                n = 0.
                M = 1.5
                c = 0.025
            else if(Dgrm <= 60 .and. Dgrm >= 1 ) then
                n = 1 - 0.56*log10(Dgrm)
                M = 9.66/Dgrm + 1.34
                c = 10**(2.86*log10(Dgrm) - (log10(Dgrm))**2 - 3.53)
            else !else if(Dgrm < 1) then
                n = 1
                M = 11
                c = 2.95e-4
                !print *, "Dgr is < 1. Ackers and White applies when Dgr >= 1."
                !print *, "Change the methodology. Stopping..."
                !stop

            end if

            ! Particle mobility
            Fgrm = ( (Vs**n)/(sqrt(gravi*Dm*(sg-1))) )*( Vel /(sqrt(32.)*log10(10*h/Dm)) )**(1-n)

            ! Dimensionless sediment transport rate
            if (Am < Fgrm) then
                Ggrm = c*((Fgrm/Am - 1)**M)
            else
                Ggrm = 0.
            end if

            !print *,'--', Q, Ggrm, Dm, Vel, Vs, n, h

            ! Volumetric sediment transport rate G (m^3 /s )
            inStreamTransCapa_Day_1 = Q*Ggrm*Dm*((Vel/Vs)**n)/h

        end function inStreamTransCapa_Day_1

end module sed_inStreamTransCapa
