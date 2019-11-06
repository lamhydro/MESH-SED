module sed_overlandFlowTransCapa

    implicit none


    contains



        function dimenShearStress(rho, gravi, h, S, D, rhos)
        ! Dimensionless shear stresses (Shield stress)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - rhos: Sediment density (kg/m^3)
        ! Output:
        ! - dimenShearStress: Dimensionless shear stresses
            use sed_overlandFlowDetachment

            implicit none

            real, intent(in) :: rho, gravi, h, S, D, rhos
            real :: dimenShearStress, tau

            ! Flow shear stresses
            tau = flowShearStress(rho, gravi, h, S)

            ! Dimensionless shear stress
            dimenShearStress = tau/(gravi*D*(rhos-rho))

        end function dimenShearStress


        function overlandFlowTransCapa_yalin(rho, rhos, D, W, gravi, h, S, v)
        ! Overland flow transport capacity by Yalin 1963
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - rhos: Sediment density (kg/m^3)
        ! - W: Width of the flow (m) (~with of the cell)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! Output:
        ! - overlandFlowTransCapa_yalin: Volumetric overland flow transport capacity (m^3 s^-1)

            use sed_overlandFlowDetachment

            implicit none

            real, intent(in) :: rho, rhos, D, W, gravi, h, S, v
            real :: overlandFlowTransCapa_yalin, F, Fc, delta, sg, a, Vs

            ! Dimensionless shear stress
            F = dimenShearStress(rho, gravi, h, S, D, rhos)

            ! Critical dimensionless shear stress
            Fc = critDimenShearStress(rho, gravi, h, S, D, v)
            !print*, F, Fc


            if (Fc > F) then
                !delta = 0.
                overlandFlowTransCapa_yalin = 0.
            else
                ! delta
                delta = F/Fc - 1

                ! Specific gravity
                sg = rhos/rho

                ! a
                a = 2.45*(sg**(-0.4))*sqrt(Fc)

                ! Shear velocity (m/s)
                Vs = shearVelocity(rho, gravi, h, S)

                ! Overland flow transport capacity by Yalin, 1963
                overlandFlowTransCapa_yalin = W * Vs * D * 0.635 * delta * ( 1 - (1/(a*delta)) * log(1+a*delta) )

            end if


        end function overlandFlowTransCapa_yalin

        function overlandFlowTransCapa_engHan(rho, rhos, D, gravi, V, h, S, W)
        ! Overland flow transport capacity by Engenlund-Hansen 1967
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - V = Water velocity (m/s)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - rhos: Sediment density (kg/m^3)
        ! - W: Width of the flow (m) (~with of the cell)
        ! Output:
        ! - overlandFlowTransCapa_engHan: Volumetric overland flow transport capacity (m^3 s^-1)
            implicit none

            real, intent(in) :: rho, rhos, D, gravi, V, h, S, W
            real :: overlandFlowTransCapa_engHan, sg


            ! Specific gravity
            sg = rhos/rho
            !print *, rhos, ' ', rho, ' ',  sg

            overlandFlowTransCapa_engHan = (0.05 * W * (V**2) * (h**1.5) * (S**1.5))/( ((sg-1)**2)*D*sqrt(gravi) )


        end function overlandFlowTransCapa_engHan


        function overlandFlowTransCapa(method, rho, rhos, D, gravi, V, h, S, W, vis)
            implicit none

            integer, intent(in) :: method
            real, intent(in) :: rho, rhos, D, gravi, V, h, S, W, vis
            real :: overlandFlowTransCapa


            ! By Yalin 1963
            if (method == 1) then
                overlandFlowTransCapa = overlandFlowTransCapa_yalin(rho, rhos, D, W, gravi, h, S, vis)
            ! By Engelund-Hansen 1967
            else if(method == 2) then
                overlandFlowTransCapa = overlandFlowTransCapa_engHan(rho, rhos, D, gravi, V, h, S, W)
            else
                print *
                print *, "Stopping... Choose a method for overland flow trans. capacity."
                print *
                stop

            end if


        end function overlandFlowTransCapa

        subroutine varsAtCellEdge()
            use sed_vars
            implicit none
            real :: d1, d2

            do i = 1, NA
                if (cn(i)%east /= 0) then
                    d1 = grs(i)%DELX*0.5
                    d2 = grs(cn(i)%east)%DELX*0.5
                    waterDepth_edge(i)%east = (ofh(i)%depth*d1 + ofh(cn(i)%east)%depth*d2)/(d1+d2)
                    waterSslope_edge(i)%east = (ofh(i)%slope*d1 + ofh(cn(i)%east)%slope*d2)/(d1+d2)
                    flowVelocity_edge(i)%east = (ofh(i)%velocity*d1 + ofh(cn(i)%east)%velocity*d2)/(d1+d2)
                    discharge_edge(i)%east = (ofh(i)%discharge*d1 + ofh(cn(i)%east)%discharge*d2)/(d1+d2)
                end if
                if (cn(i)%north /= 0) then
                    d1 = grs(i)%DELY*0.5
                    d2 = grs(cn(i)%north)%DELY*0.5
                    waterDepth_edge(i)%north = (ofh(i)%depth*d1 + ofh(cn(i)%north)%depth*d2)/(d1+d2)
                    waterSslope_edge(i)%north = (ofh(i)%slope*d1 + ofh(cn(i)%north)%slope*d2)/(d1+d2)
                    flowVelocity_edge(i)%north = (ofh(i)%velocity*d1 + ofh(cn(i)%north)%velocity*d2)/(d1+d2)
                    discharge_edge(i)%north = (ofh(i)%discharge*d1 + ofh(cn(i)%north)%discharge*d2)/(d1+d2)
                end if
                if (cn(i)%west /= 0) then
                    d1 = grs(i)%DELX*0.5
                    d2 = grs(cn(i)%west)%DELX*0.5
                    waterDepth_edge(i)%west = (ofh(i)%depth*d1 + ofh(cn(i)%west)%depth*d2)/(d1+d2)
                    waterSslope_edge(i)%west = (ofh(i)%slope*d1 + ofh(cn(i)%west)%slope*d2)/(d1+d2)
                    flowVelocity_edge(i)%west = (ofh(i)%velocity*d1 + ofh(cn(i)%west)%velocity*d2)/(d1+d2)
                    discharge_edge(i)%west = (ofh(i)%discharge*d1 + ofh(cn(i)%west)%discharge*d2)/(d1+d2)
                end if
                if (cn(i)%south /= 0) then
                    d1 = grs(i)%DELY*0.5
                    d2 = grs(cn(i)%south)%DELY*0.5
                    waterDepth_edge(i)%south = (ofh(i)%depth*d1 + ofh(cn(i)%south)%depth*d2)/(d1+d2)
                    waterSslope_edge(i)%south = (ofh(i)%slope*d1 + ofh(cn(i)%south)%slope*d2)/(d1+d2)
                    flowVelocity_edge(i)%south = (ofh(i)%velocity*d1 + ofh(cn(i)%south)%velocity*d2)/(d1+d2)
                    discharge_edge(i)%south = (ofh(i)%discharge*d1 + ofh(cn(i)%south)%discharge*d2)/(d1+d2)
                end if

            end do

        end subroutine varsAtCellEdge


        subroutine overlandFlowTransCapa_outCell()
            use sed_vars

            implicit none

    !        real :: waterDepth_east, waterSslope_east, flowVelocity_east
    !        real :: waterDepth_north, waterSslope_north, flowVelocity_north
    !        real :: waterDepth_west, waterSslope_west, flowVelocity_west
    !        real :: waterDepth_south, waterSslope_south, flowVelocity_south

            G_out(:)%east = 0.; G_out(:)%north = 0.; G_out(:)%west = 0.; G_out(:)%south = 0.
            do i = 1, NA
                if (cn(i)%east>i) then
    !                waterDepth_east = (waterDepth(i) + waterDepth(cn(i)%east))/2
    !                waterSslope_east = (waterSslope(i) + waterSslope(cn(i)%east))/2
    !                flowVelocity_east = (flowVelocity(i) + flowVelocity(cn(i)%east))/2
                    G_out(i)%east = min(overlandFlowTransCapa(overlFlowCapaMethod, rhow, sca(i)%density, &
                                                    sca(i)%diameter/1000., &
                                                    gravi, flowVelocity_edge(i)%east, &
                                                    waterDepth_edge(i)%east/1000., &
                                                    waterSslope_edge(i)%east, ofh(i)%width,vis), &
                                                    FPCRIT*discharge_edge(i)%east)
                end if
                if (cn(i)%north>i) then
    !                waterDepth_north = (waterDepth(i) + waterDepth(cn(i)%north))/2
    !                waterSslope_north = (waterSslope(i) + waterSslope(cn(i)%north))/2
    !                flowVelocity_north = (flowVelocity(i) + flowVelocity(cn(i)%north))/2
                    G_out(i)%north = min(overlandFlowTransCapa(overlFlowCapaMethod, rhow, sca(i)%density, &
                                                    sca(i)%diameter/1000., &
                                                    gravi, flowVelocity_edge(i)%north, &
                                                    waterDepth_edge(i)%north/1000., &
                                                    waterSslope_edge(i)%north, ofh(i)%width,vis),&
                                                    FPCRIT*discharge_edge(i)%north)
                end if
                if (cn(i)%west>i) then
    !                waterDepth_west = (waterDepth(i) + waterDepth(cn(i)%west))/2
    !                waterSslope_west = (waterSslope(i) + waterSslope(cn(i)%west))/2
    !                flowVelocity_west = (flowVelocity(i) + flowVelocity(cn(i)%west))/2
                    G_out(i)%west = min(overlandFlowTransCapa(overlFlowCapaMethod, rhow, sca(i)%density, &
                                                    sca(i)%diameter/1000., &
                                                    gravi, flowVelocity_edge(i)%west, &
                                                    waterDepth_edge(i)%west/1000., &
                                                    waterSslope_edge(i)%west, ofh(i)%width,vis), &
                                                    FPCRIT*discharge_edge(i)%west)
                end if
                if (cn(i)%south>i) then
    !                waterDepth_south = (waterDepth(i) + waterDepth(cn(i)%south))/2
    !                waterSslope_south = (waterSslope(i) + waterSslope(cn(i)%south))/2
    !                flowVelocity_south = (flowVelocity(i) + flowVelocity(cn(i)%south))/2
                    G_out(i)%south = min(overlandFlowTransCapa(overlFlowCapaMethod, rhow, sca(i)%density, &
                                                    sca(i)%diameter/1000., &
                                                    gravi, flowVelocity_edge(i)%south, &
                                                    waterDepth_edge(i)%south/1000., &
                                                    waterSslope_edge(i)%south, ofh(i)%width,vis), &
                                                    FPCRIT*discharge_edge(i)%south)
                end if

                !if i is not nextV
                    !
                    ! We estimate fluxes throug the boundaries
                !next(ipos,jpos)
            end do


        end subroutine overlandFlowTransCapa_outCell




    !    subroutine overlandFlowTransCapaCell()
    !    ! Estimate overland flow transport capacity for each active cell
    !
    !        use sed_vars
    !
    !        implicit none
    !
    !        !real :: dummy
    !
    !        !allocate(G(yCount,xCount))
    !        allocate(G(NA))
    !
    !        ! By Yalin 1963
    !        if (overlFlowCapaMethod == 1) then
    !            do i = 1, NA
    !            !do i = 1,yCount
    !            !    do j = 1,xCount
    !            !        if (cell(i,j) == 1) then
    !                        !dummy = dimenShearStress(rhow, gravi, waterDepth(i,j)/1000., waterSslope(i,j), &
    !                        !        cellSoilChar(i,j)%diameter/1000., cellSoilChar(i,j)%density)
    !                        !print*, dummy
    !                        G(i) = overlandFlowTransCapa_yalin(rhow, cellSoilChar(i)%density, &
    !                                cellSoilChar(i)%diameter/1000., flowWidth(i), &
    !                                gravi, waterDepth(i)/1000., waterSslope(i), vis)
    !             !       else
    !             !           G(i,j) = -9999.9
    !             !       end if
    !
    !             !   end do
    !            !end do
    !            end do
    !        ! By Engelund-Hansen 1967
    !        else if(overlFlowCapaMethod == 2) then
    !            do i = 1, NA
    !            !do i = 1,yCount
    !            !    do j = 1,xCount
    !            !        if (cell(i,j) == 1) then
    !                        G(i) = overlandFlowTransCapa_engHan(rhow, cellSoilChar(i)%density, &
    !                        cellSoilChar(i)%diameter/1000., gravi, flowVelocity(i), waterDepth(i)/1000., &
    !                        waterSslope(i), flowWidth(i))
    !            !        else
    !            !            G(i,j) = -9999.9
    !            !        end if
    !
    !            !    end do
    !            !end do
    !            end do
    !        else
    !            print *
    !            print *, "Stopping... Choose a method for overland flow trans. capacity."
    !            print *
    !            stop
    !
    !        end if
    !
    !	end subroutine overlandFlowTransCapaCell

end module sed_overlandFlowTransCapa
