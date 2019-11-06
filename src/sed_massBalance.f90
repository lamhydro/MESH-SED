!------------------------------------------------------------------------------
! Sediment transport in cold region catchments: the MESH-SED model
!------------------------------------------------------------------------------
!
! MODULE: sed_massBalance
!
!> @author
!> Luis Morales, GIWS & GWF.
!
! DESCRIPTION:
!> Module to estimate the mass balance through the river network
!
! REVISION HISTORY:
! 06 Sep 2019 - Initial Version
!------------------------------------------------------------------------------


module sed_massBalance
    use sed_chanBankErosion
    use sed_fineSedDepoAndResus

    implicit none
    contains

        real function depositedLoad(alpha, B, D, g, v, rhos, rho, h, S, C, l, parDLim, frac)
            implicit none
            real, intent(in) :: alpha, B, D, g, v, rhos, rho, h, S, C, l, parDLim, frac

            if (D > parDLim) then
                !depositedLoad = 0.
                depositedLoad = frac*sedSettlVelo(D*1.e-3, g, v, rhos, rho) * B * l * C !<  [kg s^-1]
                !print *, '->> coarse'
            else
                depositedLoad = frac*rateOfDeposition(alpha, B, D*1.e-3, g, v, rhos, rho, h, S, C/rhos)*l !< From [kg m^-1 s^-1] to [kg s^-1]
                !print *, '->> fine'
            end if

        end function depositedLoad

        real function outFlowLoad(frac, C, Q)
            implicit none
            real, intent(in) :: frac, C, Q
            outFlowLoad = frac * C * Q
        end function outFlowLoad

        real function sedMassBalanceDeri(C, Q, L_in, Vol, L_bank, L_hill, L_res, alpha, B, D, g, v, &
                                        rhos, rho, h, S, l, parDLim, frac)

            implicit none

            real, intent(in) :: C, Q, L_in, Vol, L_bank,  L_hill, L_res
            real, intent(in) :: alpha, B, D, g, v, rhos, rho, h, S, l, parDLim, frac
            real :: massIn, massOut, L_out, L_dep

            !> Deposited load
            L_dep = depositedLoad(alpha, B, D, g, v, rhos, rho, h, S, C, l, parDLim, frac)
            !if (D > parDLim) then
            !    L_dep = 0.
                !print *, '->> coarse'
            !else
            !    L_dep = frac*rateOfDeposition(alpha, B, D*1.e-3, g, v, rhos, rho, h, S, C/rhos)*l !< From [kg m^-1 s^-1] to [kg s^-1]
                !print *, '->> fine'
            !end if

            !> Outflow Load
            !L_out = frac * C * Q
            L_out = outFlowLoad(frac, C, Q)

            massIn  = L_in + L_bank + L_hill + L_res
            massOut = L_out + L_dep

            sedMassBalanceDeri = (massIn - massOut)/Vol
        end function sedMassBalanceDeri



        real function C_new_rk4(C, Q, L_in, Vol, L_bank, L_hill, L_res, alpha, B, D, g, v, rhos, rho, h, S, l, parDLim, dt, frac)
            implicit none

            real, intent(in) :: C, Q, L_in, Vol, L_bank, L_hill, L_res, alpha, B, D, g, v, rhos, rho, h, S, l, parDLim, frac
            integer, intent(in) :: dt
            real :: k1, k2, k3, k4

            !print *, 'C_new_rk4-------------------'

            k1 = sedMassBalanceDeri(C, Q, L_in, Vol, L_bank, L_hill,L_res,alpha, &
                                    B, D, g, v, rhos, rho, h, S, l, parDLim, frac)
            k2 = sedMassBalanceDeri(C+0.5*dt*k1, Q, L_in, Vol, L_bank, L_hill,L_res,alpha, &
                                    B, D, g, v, rhos, rho, h, S, l, parDLim, frac)
            k3 = sedMassBalanceDeri(C+0.5*dt*k2, Q, L_in, Vol, L_bank, L_hill,L_res,alpha, &
                                    B, D, g, v, rhos, rho, h, S, l, parDLim, frac)
            k4 = sedMassBalanceDeri(C+dt*k3, Q, L_in, Vol, L_bank, L_hill,L_res,alpha, &
                                    B, D, g, v, rhos, rho, h, S, l, parDLim, frac)


            C_new_rk4 = C + (1./6.)*(k1+2*k2+2*k3+k4)*dt

        end function C_new_rk4


        subroutine hillSlopeLoad_grid()
            use sed_vars
            use sed_overlandFlowTransCapa
            use sed_rainDropDetachment
            use sed_overlandFlowDetachment
            implicit none

            real :: totalHillSero, overFTcapa, aux

            !> Raindrop detachment
            call rainDropDetachCell()

            !> Overland flow detachment
            call overlandFlowDetachCell()


            !> Loads from hillslopes
            do i = 1, NA

                totalHillSero = (D_R(i)+D_F(i))*(grs(i)%DELX*grs(i)%DELY) !> From kg m^-2 s^-1 to kg s^-1
!                overFTcapa = min(overlandFlowTransCapa(overlFlowCapaMethod, rhow, sca(i)%density, &
!                                                    sca(i)%diameter/1000., &
!                                                    gravi, ofh(i)%velocity, &
!                                                    ofh(i)%depth/1000., &
!                                                    ofh(i)%slope, ofh(i)%width,vis), &
!                                                    FPCRIT*ofh(i)%discharge)
                overFTcapa = min(overlandFlowTransCapa(overlFlowCapaMethod, rhow, sca(i)%density, &
                                                    sca(i)%diameter/1000., &
                                                    gravi, ofh(i)%velocity, &
                                                    ofh(i)%depth/1000., &
                                                    ofh(i)%slope, ofh(i)%width,vis)*sca(i)%density, &
                                                    totalHillSero)


                aux = min(FPCRIT*ofh(i)%discharge*sca(i)%density,overFTcapa)
                !> Loop for each sediment particle size
                do j = 1, nsedpar
                    cmb(i)%L_hill(j) = sca(i)%frac(j)*aux
                end do

            end do
        end subroutine hillSlopeLoad_grid

        !> River bank erosion
        subroutine bankEroLoad_grid()
            use sed_vars
            use sed_chanBankErosion
            implicit none

            real :: Eb

            do i = 1, NA

                !> River bank erosion rate from \f$ (kg m^{-2} s^{-1}) \f$ to  \f$ (kg s^{-1}) \f$
                Eb = 2*bankRateErosion(bsca(i)%chanBankDetach*1.e-6, rhow, gravi, &
                                    rh(i)%depth, rh(i)%slope, bsca(i)%diameter*1e-3, &
                                    vis, bsca(i)%density, rh(i)%width)*(rh(i)%depth*rh(i)%length)

                !> Loop for each sediment particle size
                do j = 1, nsedpar
                    cmb(i)%L_bank(j) =  bsca(i)%frac(j) * Eb
                end do

            end do

        end subroutine bankEroLoad_grid

        subroutine streamBottomResLoad_grid()
            use sed_vars
            use sed_overlandFlowDetachment

            implicit none
            real :: res

            do i = 1, NA

                res = overlandFlowDetachment(rhow, gravi, rh(i)%depth,&
                            rh(i)%slope, bsca(i)%diameter/1000., vis, cbsca(i)%density, &
                            bsca(i)%chanBankDetach*1.e-6)

                !> Loop for each sediment particle size
                do j = 1, nsedpar
                    cmb(i)%L_res(j) =  cbsca(i)%frac(k) * res
                end do
            end do

        end subroutine streamBottomResLoad_grid


        !> Instream inflow loads from upstream
!        subroutine upstreamLoad_grid()
!            use sed_vars
!            implicit none
!
!            real :: L_in
!            integer :: rIn
!
!            do k = 1, nsedpar
!                do i = 1, NA
!                    L_in = 0.
!                    if (ion(i)%NreachIn > 0) then
!
!                        do j = 1, ion(i)%NreachIn
!                            rIn = ion(i)%reachIn(j)
!                            L_in = L_in + cmb(rIn)%C(k)*rh(rIn)%discharge
!                        end do
!
!                    end if
!                    cmb(i)%L_in(k) = L_in
!                end do
!            end do
!
!        end subroutine upstreamLoad_grid



        subroutine inStreamTransCapa_grid
            use sed_vars
            use sed_overlandFlowTransCapa
            use sed_inStreamRouting

            implicit none

            real :: Gc, wetArea, f_i

            !> Estimate the longitudinal sediment velocity (m/s)
            call longSedVelocity_grid

            do i = 1, NA
                wetArea =  rh(i)%depth*rh(i)%width
                do k = 1, nsedpar
                    if (cbsca(i)%frac(k) > 1.e-6) then
                        if (sp(k)%meanD > parDLim) then   ! For coarse sediment
                            !f_i = frac_fi(isrrB(i)%C(k),SFDELi, cbscaB(i)%ly(1)%frac(k)*cbscaB(i)%ly(1)%thick, dc)
                            f_i = 1.
                            if (instreamFlowCapaMethod == 1) then
                                Gc = f_i*overlandFlowTransCapa_engHan(rhow, cbsca(i)%density, &
                                                                sp(k)%meanD*1e-3, gravi, rh(i)%velocity, &
                                                                rh(i)%depth, rh(i)%slope, rh(i)%width)
                            else if (instreamFlowCapaMethod == 2) then
                                Gc = f_i*inStreamTransCapa_AckersWhite(rhow, cbsca(i)%density, &
                                                                sp(k)%meanD*1e-3, gravi, rh(i)%depth, rh(i)%slope, vis, &
                                                                rh(i)%velocity, rh(i)%discharge)
                            else if (instreamFlowCapaMethod == 3) then
                                Gc = f_i*inStreamTransCapa_Day_1(rhow, cbsca(i)%density, sp(k)%meanD*1e-3, &
                                                cbsca(i)%ly(1)%D16*1e-3, cbsca(i)%ly(1)%D50*1e-3, &
                                                cbsca(i)%ly(1)%D84*1e-3, gravi, rh(i)%depth, rh(i)%slope, &
                                                vis, rh(i)%velocity, rh(i)%discharge)
                            else
                                print *, "Wrong in-stream flow transport capacity method."
                                print *, "Stopping..."
                                stop
                            end if

                            !cmb(i)%C_pot(k) = cbsca(i)%density * sedCon(Gc, isrr(i)%Vs(k), wetArea) !> Potential concentration in kg/m3 based on transport capacity
                            cmb(i)%C_pot(k) = cbsca(i)%density * sedCon(Gc*cbsca(i)%frac(k), isrr(i)%Vs(k), wetArea) !> Potential concentration in kg/m3 based on transport capacity

                            !print *, 'Method: ', instreamFlowCapaMethod, cmb(i)%C_pot(k), Gc
                            !print *, 'Vs ', isrr(i)%Vs(k), '  V: ', rh(i)%velocity, &
                            !            sedCon(Gc*cbsca(i)%frac(k), isrr(i)%Vs(k), wetArea), cmb(i)%C_pot(k)
                            !stop
                        else ! For fine sediment

                            cmb(i)%C_pot(k) = cbsca(i)%frac(k)*cbsca(i)%density * FPCRIT

                        end if
                    else
                        cmb(i)%C_pot(k) = 0.
                    end if
                end do

            end do



        end subroutine inStreamTransCapa_grid





!        subroutine massBalance()
!                use sed_vars
!                implicit none
!
!                real :: delta, V, aux
!
!
!                !< Instream transport capacity
!                call inStreamTransCapa_grid()
!
!                do k = 1, nsedpar
!
!                    if (sp(k)%meanD > parDLim) then   ! For coarse sediment
!                        do i = 1, NA
!
!                            V = rh(i)%depth * rh(i)%width * rh(i)%length
!                            aux = C_new_rk4(cmb(i)%C(k), ofh(i)%discharge, cmb(i)%L_in(k), V, &
!                                        cmb(i)%L_bank(k),cmb(i)%L_hill(k), ratioStress, rh(i)%width, &
!                                        sp(k)%meanD, gravi, vis, cbsca(i)%density, rhow,rh(i)%depth, &
!                                        rh(i)%slope, rh(i)%length, parDLim, DELT)
!
!                                        !C_new_rk4(C, Q, L_in, Vol, L_bank, L_hill, alpha, B, D, g, v, rhos, rho, h, S, l, parDLim, dt)
!                            cmb(i)%C(k) = min(aux, cmb(i)%C_pot(k))
!
!                            !
!
!                            delta = cmb(i)%C(k)-cmb(i)%C_pot(k)
!
!                            !if delta < 0. then !< Erosion from bed
!                            !    isrn(i)%Dz(k) = isrn(i)%Dz(k) + (G_old - isrr(i)%G(k))/(rh(i)%width*isrr(i)%Vs(k)*(1-cbsca(i)%porosity))
!
!                            !else !< Deposition on bed
!                            !end if
!
!                            print *, 'herere'
!                            stop
!
!                        end do
!
!                    else !< For fine sediments
!                        do i = 1, NA
!
!                            V = rh(i)%depth * rh(i)%width * rh(i)%length
!                            aux = C_new_rk4(cmb(i)%C(k), ofh(i)%discharge, cmb(i)%L_in(k), V, &
!                                        cmb(i)%L_bank(k),cmb(i)%L_hill(k), ratioStress, rh(i)%width, &
!                                        sp(k)%meanD, gravi, vis, cbsca(i)%density, rhow, rh(i)%depth, &
!                                        rh(i)%slope, rh(i)%length, parDLim, DELT)
!                            cmb(i)%C(k) = min(aux, cmb(i)%C_pot(k))
!
!
!
!                            delta = cmb(i)%C(k)-cmb(i)%C_pot(k)
!
!                        end do
!
!
!                    end if
!
!
!!                    !> Instream inflow loads from upstream
!!                    L_in = 0.
!!                    if (ion(i)%NreachIn > 0) then
!!                        do j = 1, ion(i)%NreachIn
!!                            rIn = ion(i)%reachIn(j)
!!                            L_in = L_in + cmb(i)%L_in
!!                        end do
!!                    end if
!!
!!                    !do k = 1, nsedpar
!!
!!
!!
!!                    !isrr(i)%C(k) = isrrB(i)%C(k)+sedMassBalanceDeri*DELT
!!                    L_in = !> Load from upstream
!!                    L_sero = !> Load from instream erosion: bank erosion and bottom erosion
!!                    L_dep = !> Sediment deposition
!!                    L_hero = csa(i)%C(k)* !> Load from hillslope erosion
!!                    isrr(i)%C(k) = c_new_rk4(c, Qout, L_in, V, L_sero, L_dep, L_hero, dt)
!
!
!                !end do
!
!            end do
!
!        end subroutine massBalance

        subroutine massBalance2()
            use sed_vars
            implicit none

            real :: delta, V, C_new, L_in
            integer :: rIn


            !< Instream transport capacity
            call inStreamTransCapa_grid()

            do i = 1, NA
                V = rh(i)%depth * rh(i)%width * rh(i)%length
                do k = 1, nsedpar
                    !print *, '==========  ',  cbsca(i)%frac(k)
                    if (cbsca(i)%frac(k) > 1.e-6) then
                        !print *, '==========  ',  cbsca(i)%frac(k)

                        !> Estimate L_in first
                        L_in = 0
                        if (ion(i)%NreachIn > 0) then
                            do j = 1, ion(i)%NreachIn
                                rIn = ion(i)%reachIn(j)
                                L_in = L_in + cmb(rIn)%C(k) * cbsca(i)%frac(k) * rh(rIn)%discharge
                            end do
                        end if
                        cmb(i)%L_in(k) = L_in

                        !> Deposited load
                        cmb(i)%L_dep(k) = depositedLoad(ratioStress, rh(i)%width, sp(k)%meanD, gravi, vis,&
                                                        cbsca(i)%density, rhow, rh(i)%depth, rh(i)%slope, &
                                                        cmb(i)%C(k), rh(i)%length, parDLim, cbsca(i)%frac(k))

                        !> Outflow load
                        cmb(i)%L_out(k) = outFlowLoad(cbsca(i)%frac(k), cmb(i)%C(k), rh(i)%discharge)


                        C_new = C_new_rk4(cmb(i)%C(k), rh(i)%discharge, cmb(i)%L_in(k), V, &
                                    cmb(i)%L_bank(k),cmb(i)%L_hill(k),cmb(i)%L_res(k), ratioStress, rh(i)%width, &
                                    sp(k)%meanD, gravi, vis, cbsca(i)%density, rhow,rh(i)%depth, &
                                    rh(i)%slope, rh(i)%length, parDLim, DELT, cbsca(i)%frac(k))

                        if (C_new < 0.) then
                            C_new = cbsca(i)%frac(k)*cbsca(i)%density * FPCRIT
                        end if

!                        C_new = C_new_rk4(cmb(i)%C(k), rh(i)%discharge, L_in, V, &
!                                    cmb(i)%L_bank(k),cmb(i)%L_hill(k), ratioStress, rh(i)%width, &
!                                    sp(k)%meanD, gravi, vis, cbsca(i)%density, rhow,rh(i)%depth, &
!                                    rh(i)%slope, rh(i)%length, parDLim, DELT, cbsca(i)%frac(k))

                                    !C_new_rk4(C, Q, L_in, Vol, L_bank, L_hill, alpha, B, D, g, v, rhos, rho, h, S, l, parDLim, dt)

                        !> Deposited load
                        !cmb(i)%L_dep(k) = cbsca(i)%frac(k)*rateOfDeposition(alpha, B, D*1.e-3, g, v, rhos, rho, h, S, C/rhos)*l

                        cmb(i)%C(k) = min(C_new, cmb(i)%C_pot(k))

                        delta = cmb(i)%C(k)-cmb(i)%C_pot(k)
                    else
                        cmb(i)%C(k) = 0.

                    end if
                end do

            end do


        end subroutine massBalance2


end module sed_massBalance
