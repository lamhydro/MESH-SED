module sed_chanBankErosion
    implicit none
    contains

        real function K(B,H)
        ! Proportionatily contant for rectangular channels.
        ! Input:
        ! - B: Effective channel width (m)
        ! - H: Flow depth (m)
        ! Output:
        ! - K: Proportionatily constant for rectangular channels.
            implicit none

            real, intent(in) :: B, H
            real :: ratio, a4, b4

            ratio = B/H

            if (ratio < 1) then
                a4 = 0.05; b4 = 0.41
            else if (ratio>=1 .and. ratio<2) then
                a4 = 0.24; b4 = 0.22
            else if (ratio>=2 .and. ratio<4) then
                a4 = 0.61; b4 = 0.035
            else
                a4 = 0.75; b4 = 0.0
            end if
            K = a4 + b4*ratio


        end function K

        real function bankFlowShearStress(rho, gravi, h, S, B)
        ! Flow shear stresses at channel bank(N/m^2)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - B: Effective channel width (m)
        ! Output:
        ! - bankFlowShearStress: Flow shear stresses at channel bank (N/m^2)
            use sed_overlandFlowDetachment
            implicit none
            real, intent(in) :: rho, gravi, h, S, B

            bankFlowShearStress = K(B,h)*flowShearStress(rho, gravi, h, S)

        end function bankFlowShearStress

        real function bankRateErosion(Kb, rho, gravi, h, S, D, v, rhos, B)
        ! Rate of erosion by channel flow at one of the two channel banks.
        ! Input:
        ! - Kb: Bank erodability coefficient (kg m^-2 s^-1)
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! - B: Effective channel width (m)
        ! - H: Flow depth (m)
        ! Output:
        ! - bankRateErosion: Rate of detachment of material per unit area of bank (kg m^-2 s^-1)
            use sed_overlandFlowDetachment
            implicit none
            real, intent(in) :: Kb, rho, gravi, h, S, D, v, rhos, B

            real :: tau_b, tau_bc

            tau_b = bankFlowShearStress(rho, gravi, h, S, B)
            tau_bc = critShearStress(rho, gravi, h, S, D, v, rhos)

            if (tau_b > tau_bc) then
                bankRateErosion = Kb * (tau_b/tau_bc - 1)
            else
                bankRateErosion = 0.
            end if

        end function bankRateErosion

        subroutine grid_bankErosion()
            use sed_vars
            implicit none

            real :: Eb

            do i = 1, NA
                Eb = 2*bankRateErosion(bsca(i)%chanBankDetach, rhow, gravi, &
                                    rh(i)%depth, rh(i)%slope, bsca(i)%diameter, &
                                    vis, bsca(i)%density, rh(i)%width)
                Eb = Eb * rh(i)%depth/bsca(i)%density    ! from kg m-2 s-1 to m3 m-1 s-1
                do j = 1, nsedpar
                    isrr(i)%bke(j) =  bsca(i)%frac(j) * Eb
                end do
            end do


        end subroutine grid_bankErosion





end module sed_chanBankErosion
