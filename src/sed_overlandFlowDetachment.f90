module sed_overlandFlowDetachment

    implicit none

    contains

        real function flowShearStress(rho, gravi, h, S)
        ! Flow shear stresses (N/m^2)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! Output:
        ! - flowShearStress: Flow shear stresses (N/m^2)
            implicit none

            real, intent(in) :: rho, gravi, h, S

            flowShearStress = rho*gravi*h*S

        end function flowShearStress

        real function shearVelocity(rho, gravi, h, S)
        ! Shear velocity (m/s)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! Output:
        ! - shearVelocity: Shear velocity (m/s)
            implicit none

            real, intent(in) :: rho, gravi, h, S
            real :: tau

            ! Flow shear stresses
            tau = flowShearStress(rho, gravi, h, S)

            ! Shear velocity
            shearVelocity = sqrt(tau/rho)

        end function shearVelocity

        real function critDimenShearStress(rho, gravi, h, S, D, v)
        ! Critical dimensionless shear stresses
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! Output:
        ! - critDimenShearStress: Critical dimensionless shear stresses
            implicit none

            real, intent(in) :: rho, gravi, h, S, D, v
            real ::  Vs, Rey, a, b

            ! Shear velocity (m/s)
            Vs = shearVelocity(rho, gravi, h, S)

            ! Particle Reynolds number
            Rey = max(0.03,Vs*D/v)

            ! Critical dimensionles shear stress
            a = 0.10; b = -0.30
            if (Rey >= 0.03 .and. Rey <= 1.0) then
                a = 0.10; b = -0.30
            else if (Rey > 1.0 .and. Rey <= 6.0) then
                a = 0.10; b = -0.62
            else if (Rey > 6.0 .and. Rey <= 30.0) then
                a = 0.033; b = 0.0
            else if (Rey > 30.0 .and. Rey <= 135.0) then
                a = 0.013; b = 0.28
            else if (Rey > 135.0 .and. Rey <= 400.0) then
                a = 0.030; b = 0.10
            else if (Rey > 400.0) then
                a = 0.056; b = 0.0
            end if
            !print *, 'Rey:', Rey, a , b
            critDimenShearStress = a*(Rey**b)

        end function critDimenShearStress

        real function critShearStress(rho, gravi, h, S, D, v, rhos)
        ! Critical shear stresses for erosion (N/m^2)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! Output:
        ! - critShearStress: Critical shear stresses for erosion  (N/m^2)
            implicit none

            real, intent(in) :: rho, gravi, h, S, D, v, rhos
            real :: Fc

            Fc = critDimenShearStress(rho, gravi, h, S, D, v)
            !print *, 'FC :', Fc, rhos, rho, gravi, D

            ! Critical shear stress
            critShearStress = Fc*(rhos - rho)*gravi*D

        end function critShearStress

        real function overlandFlowDetachment(rho, gravi, h, S, D, v, rhos, K_F)
        ! Overland flow detachment (kg m^-2 s^-1)
        ! Input:
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! - K_F: Overland flow detachment coefficient (kg m^-2 s^-1)
        ! Output:
        ! - overlandFlowDetachment: Overland flow detachment (kg m^-2 s^-1)
            implicit none

            real, intent(in) :: rho, gravi, h, S, D, v, rhos, K_F
            real :: tau, tau_c

            ! Flow shear stresses
            tau = flowShearStress(rho, gravi, h, S)

            ! Critical shear stress
            tau_c = critShearStress(rho, gravi, h, S, D, v, rhos)

            ! Overland flow detachment
            if (tau > tau_c) then
                overlandFlowDetachment = K_F*(tau/tau_c - 1)
            else
                overlandFlowDetachment = 0.
            end if

        end function overlandFlowDetachment


        subroutine overlandFlowDetachCell()
        ! Estimate overland flow detachment for each active cell

            use sed_vars

            implicit none


            !allocate(D_F(yCount,xCount))
            !allocate(D_F(NA))

            do i = 1, NA
            !do i = 1,yCount
            !    do j = 1,xCount
                    !if (cell(i,j) == 1) then
                        !D_F(i) = overlandFlowDetachment(rhow, gravi, waterDepth(i)/1000.,&
                        !            waterSslope(i), cellSoilChar(i)%diameter/1000., vis, cellSoilChar(i)%density, &
                        !            cellSoilChar(i)%overlandDetach)
                        !print *, 'Here: ', rhow, gravi, ofh(i)%depth/1000.,&
                        !            ofh(i)%slope, sca(i)%diameter/1000., vis, sca(i)%density, &
                        !            sca(i)%overlandDetach

                        !D_F(i) = overlandFlowDetachment(rhow, gravi, ofh(i)%depth/1000.,&
                        !            ofh(i)%slope, sca(i)%diameter/1000., vis, sca(i)%density, &
                        !            sca(i)%overlandDetach*1.e-6)

                        D_F(i) = overlandFlowDetachment(rhow, gravi, ofh(i)%depth,&
                                    ofh(i)%slope, sca(i)%diameter, vis, sca(i)%density, &
                                    sca(i)%overlandDetach)
                    !else
                    !    D_F(i,j) = -9999.9
                    !end if

            !    end do
            !end do
            end do

        end subroutine overlandFlowDetachCell


end module sed_overlandFlowDetachment

