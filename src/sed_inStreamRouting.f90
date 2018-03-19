module sed_inStreamRouting
    use sed_overlandFlowDetachment
    use sed_inStreamTransCapa

    implicit none

    contains

        function bedShearVelocity(gravi, h, S)
        ! Bed shear velocity
        ! Input:
        ! - gravi: 9.8 (m s^-2)
        ! - h: Stream water depth (m)
        ! - S: Water surface slope
        ! Output
        ! - bedShearVelocity: Bed shear velocity (m/s)

            implicit none
            real, intent(in) :: gravi, h, S
            real :: bedShearVelocity

            bedShearVelocity = sqrt(gravi*h*S)

        end function bedShearVelocity


        function criBedShearVelocity(rho, gravi, h, S, D, v, rhos)
        ! Critical bed shear velocity
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Stream water depth (m)
        ! - S: Water surface slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! Output:
        ! - criBedShearVelocity: Critica bed shear velocity (m/s)

            implicit none

            real, intent(in) :: rho, gravi, h, S, D, v, rhos
            real :: criBedShearVelocity, tau_c

            tau_c = critShearStress(rho, gravi, h, S, D, v, rhos)

            criBedShearVelocity = (tau_c/rho)**0.5

        end function criBedShearVelocity


        function longSedVelocity(Vel, rho, gravi, h, S, D, v, rhos)
        ! Longitudinal sediment velocity (Phillips and Sutherland, 1985)
        ! Input:
        ! - Vel: Water velocity (m/s)
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Stream water depth (m)
        ! - S: Water surface slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! Output:
        ! - longSedVelocity: Longitudinal sediment velocity (m/s)

            implicit none

            real, intent(in) :: Vel, rho, gravi, h, S, D, v, rhos
            real :: longSedVelocity, Vb, Vcb

            if (D > 0.062*1e-3) then ! For particles bigger than Silt and Clay

                Vcb = criBedShearVelocity(rho, gravi, h, S, D, v, rhos)
                Vb = bedShearVelocity(gravi, h, S)
                if (Vb == 0.) then
                    longSedVelocity = 0.
                else
                    longSedVelocity = 8.5*Vb*sqrt(1-Vcb/Vb)
                end if

                if (longSedVelocity > Vel) longSedVelocity  = Vel
            else
                longSedVelocity = Vel
            end if

        end function longSedVelocity

        subroutine longSedVelocity_grid()
            use sed_vars
            implicit none

            do i = 1, NA
                do j = 1, nsedpar
                    !if (cbsca(i)%ly(1)%frac(j) > 0.) then
                        isrr(i)%Vs(j) = longSedVelocity(rh(i)%velocity, rhow, &
                                             gravi, rh(i)%depth, rh(i)%slope, sp(j)%meanD*1e-3, &
                                             vis, cbsca(i)%density)
                    !end if
                end do
            end do

        end subroutine longSedVelocity_grid


        subroutine flowSedInput()
            use sed_vars
            implicit none

            real :: gsIn, gsOut, gsIn_of, gsOut_of, gsIn_bf, gsOut_bf
            integer :: rIn, rOut

            do i = 1, NA

                do k = 1, nsedpar

                    ! Upstream of i
                    gsIn = 0
                    if (ion(i)%NreachIn > 0) then
                        do j = 1, ion(i)%NreachIn
                            rIn = ion(i)%reachIn(j)
                            ! From overland sediment flow
                            gsIn_of = 0.25*( csa(rIn)%C(k)*ofh(rIn)%discharge + &
                                        csaB(rIn)%C(k)*ofhB(rIn)%discharge )/rh(rIn)%length
                            ! From bank sediment flow
                            gsIn_bf = 0.25*(isrrB(rIn)%bke(k) + isrr(rIn)%bke(k))
                            gsIn = gsIn + gsIn_of + gsIn_bf

                        end do
                    end if

                    ! Downstream of i
                    rOut = ion(i)%reachOut
                    ! From overland sediment flow
                    gsOut_of = 0.25*( csa(rOut)%C(k)*ofh(rOut)%discharge + csaB(rOut)%C(k)*ofhB(rOut)%discharge )/rh(rOut)%length
                    ! From bank sediment flow
                    gsOut_bf = 0.25*(isrrB(rOut)%bke(k) + isrr(rOut)%bke(k))
                    gsOut = gsOut_of + gsOut_bf

                    isrn(i)%gs(k) = gsIn + gsOut
                    !print*, gs

                end do

            end do

        end subroutine flowSedInput


        ! Sediment concentration (m3/m3)
        real function sedCon(G, Vs, A)
            implicit none
            real, intent(in) :: G, Vs, A
            if (Vs == 0.) then
                sedCon = 0.
            else
                sedCon = G/(Vs*A)
            end if
        end function sedCon

        ! Fraction of non-fine sediment particles in suspension or transported
        real function frac_fi(FDELi,SFDELi, dci, dc)
            implicit none
            real, intent(in) :: FDELi,SFDELi, dci, dc
            real, dimension(3) :: vec
            vec(1) = 0.05
            if (SFDELi == 0.) then
                vec(2) = 0.
            else
                vec(2) = FDELi/SFDELi
            end if
            vec(3) = dci/dc
            frac_fi = maxval(vec)
        end function frac_fi



        subroutine inStreamRouting
            use sed_vars
            use sed_overlandFlowTransCapa
            use sed_fineSedDepoAndResus

            implicit none

            real :: ideltat, ideltal, length_up, factor, frac1, frac2, frac3, &
                    frac4, nume1, nume2, deno, Gc
            !real, dimension(nsedpar) :: G_Day
            integer :: rIn
            real :: SFDELi, f_i, dc, G_old !, Dz_dep

            ideltat = 1/DELT
            do i = 1, NA

                ! In-stream transp. capacity by Day, 1980
    !            if (instreamFlowCapaMethod == 3) then
    !                call inStreamTransCapa_Day(rhow, cbsca(i)%density, sp%meanD*1e-3, cbsca(i)%ly(1)%frac, &
    !                                        cbsca(i)%ly(1)%D16*1e-3, cbsca(i)%ly(1)%D50*1e-3, &
    !                                        cbsca(i)%ly(1)%D84*1e-3, gravi, rh(i)%depth, rh(i)%slope, &
    !                                        vis, rh(i)%velocity, rh(i)%discharge, G_Day)
    !            end if

                ! Total depth of non-fine sediment in the top bed layer and
                ! Acc. of vol. sed. conce. of non-fine sediment
                dc = 0.
                SFDELi = 0.
                do k = 1, nsedpar
                    if (sp(k)%meanD > parDLim) then   ! For coarse sediment
                        dc = dc + cbscaB(i)%ly(1)%frac(k)*cbscaB(i)%ly(1)%thick
                        SFDELi = SFDELi + isrrB(i)%C(k)
                    end if
                end do

                ! Calculating upstream variables
                isrr(i)%Vs_up = 0.
                isrr(i)%G_up = 0.
                length_up = 0
                if (ion(i)%NreachIn > 0) then
                    do j = 1, ion(i)%NreachIn
                        rIn = ion(i)%reachIn(j)
                        isrr(i)%Vs_up = isrr(i)%Vs_up + isrr(rIn)%Vs
                        isrr(i)%G_up = isrr(i)%G_up + isrr(rIn)%G
                        length_up = length_up + rh(rIn)%length
                    end do
                    isrr(i)%Vs_up = isrr(i)%Vs_up/ion(i)%NreachIn
                    length_up = length_up/ion(i)%NreachIn
                end if

                ideltal = 1/(0.5*(length_up+rh(i)%length))
                do k = 1, nsedpar
                    !if  (cbsca(i)%ly(1)%frac(k) > 0) then

                    if (isrrB(i)%Vs(k) ==0.) then
                        frac1 = 0
                    else
                        frac1 = isrrB(i)%G(k)/isrrB(i)%Vs(k)
                    end if
                    if (isrr(i)%Vs_up(k) ==0.) then
                        frac2 = 0
                    else
                        frac2 = isrr(i)%G_up(k)/isrr(i)%Vs_up(k)
                    end if
                    if (isrrB(i)%Vs_up(k) ==0.) then
                        frac3 = 0
                    else
                        frac3 = isrrB(i)%G_up(k)/isrrB(i)%Vs_up(k)
                    end if
                    if (sp(k)%meanD > parDLim) then   ! For coarse sediment
                        f_i = frac_fi(isrrB(i)%C(k),SFDELi, cbscaB(i)%ly(1)%frac(k)*cbscaB(i)%ly(1)%thick, dc)
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

                        isrr(i)%G(k) = Gc ! Transport capacity
                        if (isrr(i)%Vs(k) ==0.) then
                            frac4 = 0
                        else
                            frac4 = isrr(i)%G(k)/isrr(i)%Vs(k)
                        end if
                        nume1 = isrn(i)%gs(k)-ideltat*(phi*(frac4 - frac1) + (1-phi)*(frac2 - frac3))
                        nume2 = - ideltal*(theta_r*(isrr(i)%G(k)-isrr(i)%G_up(k)) + &
                                (1-theta_r)*(isrrB(i)%G(k)-isrrB(i)%G_up(k)))
                        deno = (1-cbsca(i)%porosity)*rh(i)%width/DELT
                        isrn(i)%Dz(k) = (nume1+nume2)/deno

                        if (isrn(i)%Dz(k) < 0.) then
                            if (isrn(i)%Dz(k) < -cbsca(i)%ly(1)%thick * cbscaB(i)%ly(1)%frac(k)) then
                                isrn(i)%Dz(k) = -cbsca(i)%ly(1)%thick * cbscaB(i)%ly(1)%frac(k)
                                factor = -(1-cbsca(i)%porosity)*rh(i)%width*isrn(i)%Dz(k)/DELT
                                nume1 = isrn(i)%gs(k)-ideltat*(-phi*frac1 + (1-phi)*(frac2 - frac3))
                                nume2 = factor - ideltal*(-theta_r*isrr(i)%G_up(k) + &
                                        (1-theta_r)*(isrrB(i)%G(k)-isrrB(i)%G_up(k)))
                                if (isrr(i)%Vs(k) == 0.) then
                                    deno = theta_r*ideltal
                                else
                                    deno = theta_r*ideltal + ideltat*phi/isrr(i)%Vs(k)
                                end if

                                isrr(i)%G(k) = (nume1 + nume2)/deno
                            end if
                        end if
                        isrr(i)%C(k) = sedCon(isrr(i)%G(k), isrr(i)%Vs(k), rh(i)%width*rh(i)%depth)
                        ! Induce deposition of hightly concetrate flow
                        if ( isrr(i)%C(k) > FPCRIT) then
                            isrr(i)%C(k) = FPCRIT
                            G_old = isrr(i)%G(k)
                            isrr(i)%G(k) = isrr(i)%C(k) * isrr(i)%Vs(k)*rh(i)%width*rh(i)%depth
                            isrn(i)%Dz(k) = isrn(i)%Dz(k) + (G_old - isrr(i)%G(k))*DELT/(rh(i)%width*rh(i)%length)
                        end if

                    else                            ! For fine sediment (silt and clay)

                        if (rh(i)%discharge == 0.) then
                            isrr(i)%G(k) = 0.
                            isrr(i)%C(k) = sedCon(isrr(i)%G(k), isrr(i)%Vs(k), rh(i)%width*rh(i)%depth)

                            if (isrr(i)%Vs(k) == 0.) then
                                frac4 = 0
                            else
                                frac4 = isrr(i)%G(k)/isrr(i)%Vs(k)
                            end if
                            nume1 = isrn(i)%gs(k)-ideltat*(phi*(frac4 - frac1) + (1-phi)*(frac2 - frac3))
                            nume2 = - ideltal*(theta_r*(isrr(i)%G(k)-isrr(i)%G_up(k)) + &
                                    (1-theta_r)*(isrrB(i)%G(k)-isrrB(i)%G_up(k)))
                            deno = (1-cbsca(i)%porosity)*rh(i)%width/DELT
                            isrn(i)%Dz(k) = (nume1+nume2)/deno
                        else
                            !Dz_dep = depthOfDeposition(ratioStress, rh(i)%width, sp(k)%meanD,
                            !                            gravi, vis, cbsca(i)%density, rhow,
                            !                            rh(i)%depth, rh(i)%slope, isrrB(i)%C(k), DELT)

                            isrn(i)%Dz(k) = -cbsca(i)%ly(1)%thick * cbsca(i)%ly(1)%frac(k)
                            factor = -(1-cbsca(i)%porosity)*rh(i)%width*isrn(i)%Dz(k)/DELT
                            nume1 = isrn(i)%gs(k)-ideltat*(-phi*frac1 + (1-phi)*(frac2 - frac3))
                            nume2 = factor - ideltal*(-theta_r*isrr(i)%G_up(k) + &
                                    (1-theta_r)*(isrrB(i)%G(k)-isrrB(i)%G_up(k)))
                            if (isrr(i)%Vs(k) == 0.) then
                                deno = theta_r*ideltal
                            else
                                deno = theta_r*ideltal + ideltat*phi/isrr(i)%Vs(k)
                            end if

                            isrr(i)%G(k) = (nume1 + nume2)/deno
                            isrr(i)%C(k) = sedCon(isrr(i)%G(k), isrr(i)%Vs(k), rh(i)%width*rh(i)%depth)



                            ! Induce deposition of hightly concetrate flow
                            if ( isrr(i)%C(k) > FPCRIT) then
                                isrr(i)%C(k) = FPCRIT
                                G_old = isrr(i)%G(k)
                                isrr(i)%G(k) = isrr(i)%C(k) * rh(i)%discharge
                                isrn(i)%Dz(k) = isrn(i)%Dz(k) + (G_old - isrr(i)%G(k))*DELT/(rh(i)%width*rh(i)%length)
                            end if

                        end if


                    end if

                    ! For numerical estability porpouses
                    if (isrr(i)%G(k) < 0) then
                        isrr(i)%G(k) = 0.
                        isrr(i)%C(k) = 0.
                    end if

                end do


                ! Updating channel bed layer (active and parent layers)
                cbsca(i)%thick = cbsca(i)%thick + sum(isrn(i)%Dz) ! Total layer thickness
                if (cbsca(i)%ly(2)%D99 <= cbsca(i)%thick) then
                    cbsca(i)%ly(1)%thick = cbsca(i)%ly(2)%D99
                else
                    cbsca(i)%ly(1)%thick = cbsca(i)%thick
                end if
                cbsca(i)%ly(2)%thick  = cbsca(i)%thick - cbsca(i)%ly(1)%thick

                if (cbsca(i)%ly(1)%thick == 0. .and. cbsca(i)%ly(2)%thick == 0.) then
                    cbsca(i)%ly(1)%frac = 0.
                    cbsca(i)%ly(2)%frac = 0.
                else
                    if (sum(isrn(i)%Dz) >= 0) then ! Net deposition
                        deno = cbscaB(i)%ly(1)%thick + sum(isrn(i)%Dz)
                        do k = 1, nsedpar
                            cbsca(i)%ly(1)%frac(k) = (cbscaB(i)%ly(1)%frac(k)*cbscaB(i)%ly(1)%thick + &
                                                    isrn(i)%Dz(k))/deno
                            cbsca(i)%ly(2)%frac(k) = (cbscaB(i)%ly(2)%frac(k)*cbscaB(i)%ly(2)%thick + &
                                                   cbsca(i)%ly(1)%frac(k)*(cbsca(i)%ly(2)%thick - cbscaB(i)%ly(2)%thick))/&
                                                   cbsca(i)%ly(2)%thick
                        end do

                    else                            ! Net erosion
                        do k = 1, nsedpar
                            cbsca(i)%ly(1)%frac(k) = (cbscaB(i)%ly(1)%frac(k)*cbscaB(i)%ly(1)%thick + &
                                                    isrn(i)%Dz(k) + cbscaB(i)%ly(2)%frac(k)*(cbsca(i)%ly(1)%thick - &
                                                    cbscaB(i)%ly(1)%thick - sum(isrn(i)%Dz)))/cbsca(i)%ly(1)%thick
                            cbsca(i)%ly(2)%frac(k) = (cbscaB(i)%ly(2)%frac(k)*cbscaB(i)%ly(2)%thick + &
                                                   cbscaB(i)%ly(2)%frac(k)*(cbsca(i)%ly(2)%thick - cbscaB(i)%ly(2)%thick))/&
                                                   cbsca(i)%ly(2)%thick
                        end do

                    end if
                 end if


            end do

        end subroutine inStreamRouting


end module sed_inStreamRouting
