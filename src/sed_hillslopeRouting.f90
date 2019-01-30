module sed_hillslopeRouting
    use sed_vars

    implicit none

    contains

        subroutine inOutCellFluxes()

            implicit none

            !real :: factor ! = 1/(rho*As), As = h*W. Convert G [kg/(m^2 s)] to [m^3/s]

            do i = 1, NA

                k = 0
                if (cn(i)%east /= 0.) then
                    !factor = 1/(cellSoilChar(i)%density*flowWidth(i)*waterDepth_edge(i)%east*1.e-3)
                    if (cn(i)%east > i) then
                        F(i)%east = G_out(i)%east
                        F(i)%s_east = 0
                    else
                        F(i)%east = G_out(cn(i)%east)%west
                        F(i)%s_east = 1
                    end if
                    k = k + 1
                end if
                if (cn(i)%north /= 0.) then
                    !factor = 1/(cellSoilChar(i)%density*flowWidth(i)*waterDepth_edge(i)%north*1.e-3)
                    if (cn(i)%north > i) then
                        F(i)%north = G_out(i)%north
                        F(i)%s_north = 0
                    else
                        F(i)%north = G_out(cn(i)%north)%south
                        F(i)%s_north = 1
                    end if
                    k = k + 1
                end if
                if (cn(i)%west /= 0.) then
                    !factor = 1/(cellSoilChar(i)%density*flowWidth(i)*waterDepth_edge(i)%west*1.e-3)
                    if (cn(i)%west > i) then
                        F(i)%west = G_out(i)%west
                        F(i)%s_west = 0
                    else
                        F(i)%west = G_out(cn(i)%west)%east
                        F(i)%s_west = 1
                    end if
                    k = k + 1
                end if
                if (cn(i)%south /= 0.) then
                    !factor = 1/(cellSoilChar(i)%density*flowWidth(i)*waterDepth_edge(i)%south*1.e-3)
                    if (cn(i)%south > i) then
                        F(i)%south = G_out(i)%south
                        F(i)%s_south = 0
                    else
                        F(i)%south = G_out(cn(i)%south)%north
                        F(i)%s_south = 1
                    end if
                    k = k + 1
                end if

                F(i)%nfluxes = k

                !C(i) = acumFrac/k

            end do

        end subroutine inOutCellFluxes


        subroutine potentialSedConc()

            implicit none

            !real :: discharge_east, discharge_north, discharge_west, discharge_south, acumFrac
            real :: FQeast, FQnorth, FQwest, FQsouth, aux

            call inOutCellFluxes()

            do i = 1, NA
                if (discharge_edge(i)%east == 0.) then
                    FQeast = 0.
                else
                    FQeast = F(i)%east/discharge_edge(i)%east
                end if
                if (discharge_edge(i)%north == 0.) then
                    FQnorth = 0.
                else
                    FQnorth = F(i)%north/discharge_edge(i)%north
                end if
                if (discharge_edge(i)%west == 0.) then
                    FQwest = 0.
                else
                    FQwest = F(i)%west/discharge_edge(i)%west
                end if
                if (discharge_edge(i)%south == 0.) then
                    FQsouth = 0.
                else
                    FQsouth = F(i)%south/discharge_edge(i)%south
                end if

                !print *, i, FQeast,' ' ,FQnorth, ' ',FQwest, ' ',FQsouth, ' ', F(i)

                aux = (FQeast + FQnorth + FQwest + FQsouth)/F(i)%nfluxes
                do j = 1, nsedpar
                    !C(i) = C_out(cn(i)%east)%east + C_out(cn(i)%north)%north + C_out(cn(i)%west)%west + C_out(cn(i)%south)%south
                    !if (C_out(i)%east == 0.) then
                    !    C_in(i)%east = C_out(cn(i)%east)%west
                    !end if



                    csa(i)%C(j) = sca(i)%frac(j)*aux

        !            C(i) = (F(i)%east/discharge_edge(i)%east + F(i)%north/discharge_edge(i)%north +&
        !                    F(i)%west/discharge_edge(i)%west + F(i)%south/discharge_edge(i)%south)/F(i)%nfluxes

        !            k = 0
        !            acumFrac = 0.
        !            if (cn(i)%east /= 0.) then
        !                discharge_east = (discharge(i) + discharge(cn(i)%east))/2
        !                if (cn(i)%east > i) then
        !                    acumFrac = acumFrac + G_out(i)%east/discharge_east
        !                else
        !                    acumFrac = acumFrac + G_out(cn(i)%east)%west/discharge_east
        !                end if
        !                k = k + 1
        !            end if
        !            if (cn(i)%north /= 0.) then
        !                discharge_north = (discharge(i) + discharge(cn(i)%north))/2
        !                if (cn(i)%north > i) then
        !                    acumFrac = acumFrac + G_out(i)%north/discharge_north
        !                else
        !                    acumFrac = acumFrac + G_out(cn(i)%north)%south/discharge_north
        !                end if
        !                k = k + 1
        !            end if
        !            if (cn(i)%west /= 0.) then
        !                discharge_west = (discharge(i) + discharge(cn(i)%west))/2
        !                if (cn(i)%west > i) then
        !                    acumFrac = acumFrac + G_out(i)%west/discharge_west
        !                else
        !                    acumFrac = acumFrac + G_out(cn(i)%west)%east/discharge_west
        !                end if
        !                k = k + 1
        !            end if
        !            if (cn(i)%south /= 0.) then
        !                discharge_south = (discharge(i) + discharge(cn(i)%south))/2
        !                if (cn(i)%south > i) then
        !                    acumFrac = acumFrac + G_out(i)%south/discharge_south
        !                else
        !                    acumFrac = acumFrac + G_out(cn(i)%south)%north/discharge_south
        !                end if
        !                k = k + 1
        !            end if

                    !C(i) = acumFrac/k
                end do
            end do

        end subroutine potentialSedConc


        subroutine potentialChanSufEle()

            implicit none

            real :: B, acc, term1, term2, term3, term4

            !B = DELT/(DELX*DELY)

            do i = 1, NA
                B = DELT/(grs(i)%DELX*grs(i)%DELY)
                acc = 0.
                term1 = -B*(theta*(F(i)%south-F(i)%north) + (1-theta)*(FB(i)%south - FB(i)%north))
                term2 = -B*(theta*(F(i)%east-F(i)%west) + (1-theta)*(FB(i)%east-FB(i)%west))
                term4 = (1-sca(i)%porosity)
                do j = 1, nsedpar
    !                csa(i)%pot_Dz =( -B*(theta*(F(i)%south-F(i)%north) + (1-theta)*(FB(i)%south - FB(i)%north)) - &
    !                            B*(theta*(F(i)%east-F(i)%west) + (1-theta)*(FB(i)%east-FB(i)%west)) - &
    !                            (ofh(i)%depth*1.e-3*csa(i)%C(j)- ofhB(i)%depth*1.e-3*csaB(i)%C(j)) )/(1-sca(i)%porosity)
                    term3 = -(ofh(i)%depth*1.e-3*csa(i)%C(j)- ofhB(i)%depth*1.e-3*csaB(i)%C(j))
                    acc = acc + (term1 + term2 + term3)/term4
                end do
                csa(i)%pot_Dz = acc
            end do

        end subroutine potentialChanSufEle


        subroutine availableChanSufEle()

            implicit none

            do i = 1, NA
                csa(i)%ava_Dz = -csa(i)%SD - DELT*(D_R(i)+D_F(i))/(sca(i)%density*(1-sca(i)%porosity))
            end do

        end subroutine availableChanSufEle


        subroutine cellConAndChanSufEle()

            implicit none

            real :: B, numeTerm1, numeTerm2, numeTerm3, deno

            !B = DELT/(DELX*DELY)

            do i = 1,NA
                B = DELT/(grs(i)%DELX*grs(i)%DELY)
                !do j = 1, nsedpar

                !if () then ! Deposition when it is positive
                if (csa(i)%pot_Dz >= 0. .or. csa(i)%pot_Dz > csa(i)%ava_Dz) then ! Erosion when is negative
                    csa(i)%Dz = csa(i)%pot_Dz
                    !csa%SD(i) =  csa%SD(i) + csa%Dz(i)
                else
                    csa(i)%Dz = csa(i)%ava_Dz

                    deno = ofh(i)%depth*1.e-3 + theta*B*((1-F(i)%s_south)*discharge_edge(i)%south -&
                                      (1-F(i)%s_north)*discharge_edge(i)%north +&
                                      (1-F(i)%s_east)*discharge_edge(i)%east   -&
                                      (1-F(i)%s_west)*discharge_edge(i)%west)

                    numeTerm2 = - B*(theta*(F(i)%s_south*F(i)%south-F(i)%s_north*F(i)%north)  &
                                +   (1-theta)*(FB(i)%s_south*FB(i)%south-FB(i)%s_north*FB(i)%north))
                    numeTerm3 = - B*(theta*(F(i)%s_east*F(i)%east-F(i)%s_west*F(i)%west)  &
                                +   (1-theta)*(FB(i)%s_east*FB(i)%east-FB(i)%s_west*FB(i)%west))
                    do j = 1, nsedpar
                        numeTerm1 = (ofhB(i)%depth*1.e-3*csaB(i)%C(j)-csa(i)%Dz*(1-sca(i)%porosity))

                        if (deno /= 0.) then
                            csa(i)%C(j) = sca(i)%frac(j)*( numeTerm1 + numeTerm2 + numeTerm3 )/deno
                        else
                            csa(i)%C(j) = 0.
                        end if
                        !> Se C = 0 when C<0. Guarantee numerical stability
                        if (csa(i)%C(j)< 0.) then
                            csa(i)%C(j) = 0.
                        end if
                    end do
                end if
                G_out(i)%east  = sum(csa(i)%C) * discharge_edge(i)%east
                G_out(i)%north = sum(csa(i)%C) * discharge_edge(i)%north
                G_out(i)%west  = sum(csa(i)%C) * discharge_edge(i)%west
                G_out(i)%south = sum(csa(i)%C) * discharge_edge(i)%south

                !end do


                ! Updating the loose soil layer depth
                if (csa(i)%SD+csa(i)%Dz < 0) then
                    csa(i)%SD = 0.
                else
                    csa(i)%SD =  csa(i)%SD + csa(i)%Dz
                end if

            end do


        end subroutine cellConAndChanSufEle



    !    subroutine hydroAtCellEdges()
    !        use sed_overlandFlowTransCapa
    !        use sed_vars
    !        implicit none

    !        do i = 1, NA
    !            if (cn(k)%east>i) then
    !                waterDepth_edge(i)%west = (waterDepth(i) + waterDepth(cn(k)%east))/2
    !                waterSslope_edge(i)%west = (waterSslope(i) + waterSslope(cn(k)%east))/2
    !                flowVelocity_edge(i)%west = (flowVelocity(i) + flowVelocity(cn(k)%east))/2
    !                G_east = overlandFlowTransCapa(rhow, cellSoilChar(i)%density, &
    !                                                cellSoilChar(i)%diameter/1000., &
    !                                                gravi, flowVelocity_edge(i)%west, &
    !                                                waterDepth_edge(i)%west/1000., &
    !                                                waterSslope_edge(i)%west, flowWidth(i),vis)
    !
    !                overlandFlowTransCapa_yalin(rhow, cellSoilChar(i)%density, &
    !                                cellSoilChar(i)%diameter/1000., flowWidth(i), &
    !                                gravi, waterDepth(i)/1000., waterSslope(i), vis)
    !
    !            end if
    !            if i is not nextV
    !                !
    !                ! We estimate fluxes throug the boundaries
    !            next(ipos,jpos)
    !        end do

    !        do i = 1,yCount
    !            do j = 1,xCount
    !
    !                ! West edge
    !                if (j+1 > xCount) then
    !                    waterDepth_edge(i,j)%west = waterDepth(i,j)
    !                    waterSslope_edge(i,j)%west = waterSslope(i,j)
    !                    flowVelocity_edge(i,j)%west = flowVelocity(i,j)
    !                else
    !                    waterDepth_edge(i,j)%west = (waterDepth(i,j)+waterDepth(i,j+1))/2
    !                    waterSslope_edge(i,j)%west = (waterSslope(i,j)+waterSslope(i,j+1))/2
    !                    flowVelocity_edge(i,j)%west = (flowVelocity(i,j)+flowVelocity(i,j+1))/2
    !                end if
    !
    !                ! East edge
    !                if (j-1 < 1) then
    !                    waterDepth_edge(i,j)%east = waterDepth(i,j)
    !                    waterSslope_edge(i,j)%east = waterSslope(i,j)
    !                    flowVelocity_edge(i,j)%east = flowVelocity(i,j)
    !                else
    !                    waterDepth_edge(i,j)%east = (waterDepth(i,j)+waterDepth(i,j-1))/2
    !                    waterSslope_edge(i,j)%east = (waterSslope(i,j)+waterSslope(i,j-1))/2
    !                    flowVelocity_edge(i,j)%east = (flowVelocity(i,j)+flowVelocity(i,j-1))/2
    !                end if
    !
    !                ! South edge
    !                if (i+1 > yCount) then
    !                    waterDepth_edge(i,j)%south = waterDepth(i,j)
    !                    waterSslope_edge(i,j)%south = waterSslope(i,j)
    !                    flowVelocity_edge(i,j)%south = flowVelocity(i,j)
    !                else
    !                    waterDepth_edge(i,j)%south = (waterDepth(i,j)+waterDepth(i+1,j))/2
    !                    waterSslope_edge(i,j)%south = (waterSslope(i,j)+waterSslope(i+1,j))/2
    !                    flowVelocity_edge(i,j)%south = (flowVelocity(i,j)+flowVelocity(i+1,j))/2
    !                end if
    !
    !                ! North edge
    !                if (i-1 < 1) then
    !                    waterDepth_edge(i,j)%north = waterDepth(i,j)
    !                    waterSslope_edge(i,j)%north = waterSslope(i,j)
    !                    flowVelocity_edge(i,j)%north = flowVelocity(i,j)
    !                else
    !                    waterDepth_edge(i,j)%north = (waterDepth(i,j)+waterDepth(i-1,j))/2
    !                    waterSslope_edge(i,j)%north = (waterSslope(i,j)+waterSslope(i-1,j))/2
    !                    flowVelocity_edge(i,j)%north = (flowVelocity(i,j)+flowVelocity(i-1,j))/2
    !                end if
    !
    !                !end if

    !            end do
    !        end do


    !    end subroutine hydroAtCellEdges



end module sed_hillslopeRouting

