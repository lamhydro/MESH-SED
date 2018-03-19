module sed_rainDropDetachment

    implicit none
    contains

        function sqrMomRain(I)
        ! Squared momentum of rainfall (M_R). Rainfall erosivity
        ! Input:
        ! - I: Rainfall intensity in mm/h
        !Output
        ! - sqrMomRain: Squared momentum for rain ((kg m s^-1)^2 m^-2 s^-1)
            implicit none

            real :: sqrMomRain, alpha, beta
            real, intent(in) :: I

            ! Look up 'alpha' and 'betha' parameter in table A.1
            if (I >= 0 .and. I < 10) then
                alpha = 3214.9
                beta = 1.6896
            else if (I >= 10 .and. I < 50) then
                alpha = 583.4
                beta = 1.5545
            else if (I >= 50 .and. I < 100) then
                alpha = 133.1
                beta = 1.4242
            else if (I >= 100 .and. I < 250) then
                alpha = 29.9
                beta = 1.2821
            ! WARNING for I >= 250 there are not data in the literature.
            ! we assume they are the same than for the above range.
            else
                alpha = 29.9
                beta = 1.2821
            end if

            ! Estimate the squared momentum
            sqrMomRain = alpha*(I/3.6e6)**beta

        end function sqrMomRain



        function sqrMomLeafDrip(d,drip,draina,rho, pi, X, gravi)
        ! Squared momentum of leaf drip (M_D). Leaf drip erosivity
        ! Input:
        ! - d: Leaf drip diameter (mm)
        ! - drip: Proportion of drainage which fall as leaf drip (0-1)
        ! - draina: Canopy drainage (Evapotranspiration) (m/s)
        ! - rho: Water density (kg/m^3)
        ! - pi: Constant = 3.1416
        ! - X: Fall height of leaf drip (m)
        ! - gravi: Gravity acce. (m/s^2)
        !Output
        ! - sqrMomLeafDrip: Squared momentum for leaf drip ((kg m s^-1)^2 m^-2 s^-1)
            implicit none

            real, intent(in) :: d,drip,draina,rho, pi, X, gravi
            real :: sqrMomLeafDrip, leafDripVol, v

            ! Leaf drip fall velocity (v) (m/s)
            v = leafDripFallVel(rho,d,X,gravi,pi)

            leafDripVol = pi*((d*1e-3)**3)/6
            sqrMomLeafDrip = (drip*draina*(rho*v*leafDripVol)**2)/leafDripVol

        end function sqrMomLeafDrip



        function leafDripFallVel(rho,d,X,gravi,pi)
        ! Leaf drip fall velocity (v)
        ! Input:
        ! - d: Leaf drip diameter (mm)
        ! - rho: Water density (kg/m^3)
        ! - pi: Constant = 3.1416
        ! - X: Fall height of leaf drip (m)
        ! - gravi: Gravity acce. (m/s^2)
        !Output
        ! - leafDripFallVel: Leaf drip fall velocity (m/s)
            implicit none

            real, intent(in) :: rho, d, X, gravi, pi
            real :: leafDripFallVel, leafDripMass, beta

            ! Friction constant (beta) (kg/m)
            beta = frictionConst(d,rho,X,pi)

            leafDripMass = rho*pi*((d*1e-3)**3)/6
            leafDripFallVel = sqrt((leafDripMass*gravi/beta)*(1-exp(-2*X*beta/leafDripMass)))

        end function leafDripFallVel



        function frictionConst(d,rho,X,pi)
        ! Friction constant (beta)
        ! Input:
        ! - d: Leaf drip diameter (mm)
        ! - rho: Water density (kg/m^3)
        ! - pi: Constant = 3.1416
        ! - X: Fall height of leaf drip (m)
        !Output
        ! - frictionConst: Friction constant (kg/m)
            implicit none

            real, intent(in) :: d,rho,X,pi
            real :: frictionConst, leafDripMass


            leafDripMass = rho*pi*((d*1e-3)**3)/6
            frictionConst = 0.
            if (d <= 3.3) then
                frictionConst = leafDripMass/(2200*d*1e-3)
            end if
            if (d > 3.3 .and. X < 7.5) then
                frictionConst = leafDripMass/(1640*d*1e-3 + 1.93)
            end if
            if (d > 3.3 .and. X >= 7.5) then
                frictionConst = leafDripMass/(660*d*1e-3 + 5.14)
            end if

        end function frictionConst



        function waterDepthCorrFact(I,h,d)
        ! Water depth correction factor (Fw)
        ! Input:
        ! - I: Rainfall intensity (mm/h)
        ! - h: Water depth (mm)
        ! - d: Leaf drop diameter (mm)
        ! Output:
        ! - waterDepthCorrFact: Water depth correction factor

            implicit none

            real, intent(in) :: I,h,d
            real :: waterDepthCorrFact, dm

            if (I > 0.) then
                ! Median rain-drop diameter (m)
                dm = 1.24e-3*I**0.182
                dm = dm*1000. ! from (m) to (mm)
            else
                ! dm is equal to leaf-drop diameter (just for certain time after the rainfall event)
                dm = d
            end if

            if (h > dm) then
                waterDepthCorrFact = exp(1-h/dm)
            else
                waterDepthCorrFact = 1.
            end if


        end function waterDepthCorrFact


        function rainDropDetach(I, d, drip, draina, rho, pi, X, gravi, h, K_R, C_G, C_C)
        ! Rain-drop detachment (D_R) (kg m^-2 s^-1)
        ! Input:
        ! - I: Rainfall intensity (mm/h)
        ! - d: Leaf drip diameter (mm)
        ! - drip: Proportion of drainage which fall as leaf drip (0-1)
        ! - draina: Canopy drainage (Evapotranspiration) (m/s)
        ! - rho: Water density (kg/m^3)
        ! - pi: Constant = 3.1416
        ! - X: Fall height of leaf drip (m)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (mm)
        ! - K_R: Raindrop soil detachment coefficient (J^-1 = kg^-1 m^-2 s^2)
        ! - C_G: Proportion of soil covered by ground cover (0-1)
        ! - C_C: Proportion of soil covered by canopy cover (0-1)
        !Output
        ! - rainDropDetach: Soil detached by raindrop impac (kg  m^-2 s^-1)

            implicit none

            real, intent(in) :: I, d, drip, draina, rho, pi, X, gravi, h, K_R, C_G, C_C
            real :: rainDropDetach, M_R, M_D, Fw


            M_R = sqrMomRain(I)
            M_D = sqrMomLeafDrip(d, drip, draina, rho, pi, X, gravi)
            Fw = waterDepthCorrFact(I, h, d)

            rainDropDetach = K_R*Fw*(1 - C_G)*( (1 - C_C)*M_R + M_D )

        end function rainDropDetach





        subroutine rainDropDetachCell()
        ! Estimate rain-drop detachment for each active cell

            use sed_vars

            implicit none

            !real :: dummy, dummy1, dummy2, dummy3, dummy4

            !allocate(D_R(yCount,xCount))
            !allocate(D_R(NA))

            do i = 1, NA
            !do i = 1,yCount
            !    do j = 1,xCount
            !        if (cell(i,j) == 1) then
                        !dummy = sqrMomRain(precip(i,j))
                        !dummy1 = leafDripFallVel(rhow, cellRainDetaParam(i,j)%dropDiam, cellRainDetaParam(i,j)%fallHeight, gravi, pi)
                        !dummy2 = sqrMomLeafDrip(cellRainDetaParam(i,j)%dropDiam, cellRainDetaParam(i,j)%percDrip, &
                        !            evapotran(i,j), rhow, pi, cellRainDetaParam(i,j)%fallHeight, gravi)
                        !dummy3 = waterDepthCorrFact(precip(i,j),waterDepth(i,j),cellRainDetaParam(i,j)%dropDiam)
                        !print *, dummy, dummy1, dummy2, dummy3
                        !write(*,'(2(F15.2))') sqrMomRain(precip(i,j)), dummy1

                        !D_R(i) = rainDropDetach(precip(i), cellVegeChar(i)%dropDiam, &
                        ! cellVegeChar(i)%percDrip, evapotran(i), rhow, pi, &
                        ! cellVegeChar(i)%fallHeight, gravi, waterDepth(i), &
                        ! cellSoilChar(i)%soilDetach, cellGroundCov(i), &
                        ! cellCanopyCov(i))

                        D_R(i) = rainDropDetach(mv(i)%precip, vca(i)%dropDiam, &
                         vca(i)%percDrip, mv(i)%evapotran, rhow, pi, &
                         vca(i)%fallHeight, gravi, ofh(i)%depth, &
                         sca(i)%soilDetach, gca(i)%cellGroundCov, &
                         gca(i)%cellCanopyCov)

             !       else
             !           D_R(i,j) = -9999.9
             !       end if

             !   end do
            !end do
            end do

        end subroutine rainDropDetachCell




    !	subroutine rainfallDropEro
    !        implicit none
    !
    !        allocate(rainDropEro(NA,ngru))
    !
    !        do i = 1,NA
    !            do j = 1, ngru
    !                do k = 1, nsedClass
    !                    rainDropEro(i,j,k) =
    !                    print *, sedClass(i)%name,  sedClass(i)%diameter, sedClass(i)%density, sedClass(i)%setVelocity
    !                end do
    !            end do
    !        end do
    !
    !
    !	end subroutine rainfallDropEro

end module sed_rainDropDetachment
