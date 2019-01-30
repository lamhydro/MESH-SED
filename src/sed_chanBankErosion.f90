!------------------------------------------------------------------------------
! Sediment transport in cold region catchments: the MESH-SED model
!------------------------------------------------------------------------------
!
! MODULE: sed_chanBankErosion
!
!> @author
!> Luis Morales, GIWS & GWF.
!
! DESCRIPTION:
!> Module to estimate erosion in channel banks
!
! REVISION HISTORY:
! 01 Jul 2018 - Initial Version
! TODO_25_Jan_2019 - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module sed_chanBankErosion
    implicit none
    contains

       !---------------------------------------------------------------------------
       !> @author
       !> Luis Morales, GIWS & GWF.
       !
       ! DESCRIPTION:
       !> Estimate the proportional constant \f$ K = a4 + b4\frac{B}{H} \f$ for rectangular channels.
       !
       ! REVISION HISTORY:
       ! TODO_25_Jan_2019 - TODO_describe_appropriate_changes - TODO_name
       !
       !> @param[in] B  Effective channel width (m)
       !> @param[in] h  Water depth (m)
       !> @return K Proportional constant for rectangular channels.
       !---------------------------------------------------------------------------
        real function K(B,h)
            implicit none

            real, intent(in) :: B, h
            real :: ratio, a4, b4

            ratio = B/h

            if (ratio < 1) then
                a4 = 0.05; b4 = 0.41
            else if (ratio>=1 .and. ratio<2) then
                a4 = 0.24; b4 = 0.22
            else if (ratio>=2 .and. ratio<4) then
                a4 = 0.61; b4 = 0.035
            else
                a4 = 0.75; b4 = 0.0
            end if
            K = a4 + b4*ratio   !> a4 y b4 are constants

        end function K

       !---------------------------------------------------------------------------
       !> @author
       !> Luis Morales, GIWS & GWF.
       !
       ! DESCRIPTION:
       !> Estimate the flow shear stresses at channel bank
       !> \f$ \tau_{b} = K\tau \f$
       !> where \f$ K \f$ is a proportional constant for rectangular channels
       !> and \f$ \tau \f$  the mean flow shear stress on the bed.
       !
       ! REVISION HISTORY:
       ! TODO_25_Jan_2019 - TODO_describe_appropriate_changes - TODO_name
       !
       !> @param[in] rho    Water density \f$ (kg/m^{3}) \f$
       !> @param[in] gravi  Gravity acce. \f$(m/s^{2} ) \f$
       !> @param[in] h  Water depth (m)
       !> @param[in] S  Water surface slope ~ bottom slope
       !> @param[in] B  Effective channel width (m)
       !> @return bankFlowShearStress   Flow shear stresses at channel bank in
       !> \f$ N/m^{2} \f$
       !---------------------------------------------------------------------------
        real function bankFlowShearStress(rho, gravi, h, S, B)
            use sed_overlandFlowDetachment
            implicit none
            real, intent(in) :: rho, gravi, h, S, B

            bankFlowShearStress = K(B,h)*flowShearStress(rho, gravi, h, S)

        end function bankFlowShearStress

       !---------------------------------------------------------------------------
       !> @author
       !> Luis Morales, GIWS & GWF.
       !
       ! DESCRIPTION:
       !> Estimate the rate of erosion by channel flow at one of the two channel banks as:
       !> \f$ E_b = \left\{ \begin{array}{ll} k_b = (\frac{\tau_{b}}{\tau_{bc}} - 1)  & \mbox{if } \tau_{b} > \tau_{bc} \\ 0 & \mbox{otherwise}  \end{array} \right. \f$
       !> in \f$ (kg m^{-2} s^{-1}) \f$.
       !
       ! REVISION HISTORY:
       ! TODO_25_Jan_2019 - TODO_describe_appropriate_changes - TODO_name
       !
       !> @param[in] Kb Bank erodability coefficient \f$ (kg m^{-2} s^{-1}) \f$
       !> @param[in] rho    Water density \f$ (kg/m^{3}) \f$
       !> @param[in] gravi  Gravity acce. \f$(m/s^{2} ) \f$
       !> @param[in] h  Water depth (m)
       !> @param[in] S  Water surface slope ~ bottom slope
       !> @param[in] D  Sediment diameter (m)
       !> @param[in] v  Kinematic viscosity of water \f$ (m^{2}/s) \f$
       !> @param[in] rhos   Sediment density \f$ (kg/m^{3}) \f$
       !> @param[in] B  Effective channel width (m)
       !> @return bankRateErosion   Rate of detachment of material per unit area of bank \f$ (kg m^{-2} s^{-1}) \f$.
       !---------------------------------------------------------------------------
        real function bankRateErosion(Kb, rho, gravi, h, S, D, v, rhos, B)
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

       !---------------------------------------------------------------------------
       !> @brief
       !> Estimate the bank erosion at each cell
       !>
       !> @detail Assuming a rectangular channel, the bank erosion is estimated at
       !> each channel for each sediment particle fraction and for both banks.
       !> The bank erosion is converted from
       !> \f$ kg m^{-2} s^{-1} \f$ to \f$ m^{3} m^{-1} s^{-1} \f$
       !>
       !> @author Luis Morales (LAM), GIWS & GWF.
       !> - July, 2017
       !> @date January, 2019-LAM
       !> - Documenting the code
       !> @todo are bank erosion estimates for each sediment particle size?
       !---------------------------------------------------------------------------
        subroutine grid_bankErosion()
            use sed_vars
            implicit none

            real :: Eb
            !> Loop for cell channels
            do i = 1, NA
                Eb = 2*bankRateErosion(bsca(i)%chanBankDetach, rhow, gravi, &
                                    rh(i)%depth, rh(i)%slope, bsca(i)%diameter, &
                                    vis, bsca(i)%density, rh(i)%width)
                Eb = Eb * rh(i)%depth/bsca(i)%density
                !> Loop for each sediment particle size
                do j = 1, nsedpar
                    isrr(i)%bke(j) =  bsca(i)%frac(j) * Eb
                end do
            end do


        end subroutine grid_bankErosion



end module sed_chanBankErosion
