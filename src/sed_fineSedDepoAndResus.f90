module sed_fineSedDepoAndResus
    implicit none
    contains

        real function sedSettlVelo(D, g, v, rhos, rho)
        ! Settling velocity of fine sediment particle
        ! Input:
        ! - D: Particle diameter (m)
        ! - g: Gravity accel. (m s^-2)
        ! - v: Kinematic viscosity (m^2 s^-1)
        ! - rhos: Sediment particle density (kg m^-3)
        ! - rho: Water density (kg m^-3)
        ! Output:
        ! - sedSettlVelo: Settling velocity (m s^-1)
            implicit none
            real, intent(in) :: D, g, v, rhos, rho

            sedSettlVelo = ((D**2)*g/(18*v))*((rhos/rho) - 1)

        end function sedSettlVelo


        real function rateOfDeposition(alpha, B, D, g, v, rhos, rho, h, S, C)
        ! Rate of deposition of fine sediment
        ! Input:
        ! - alpha: Relationship of critical shear stresses for desposition and ini. of motion
        ! - Channel width (m)
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! - C: Concentration of fine particles in suspension (m3/m3)
        ! Output:
        ! - rateOfDeposition: Critical shear stresses for erosion  (Kg s^-1 m^-1)
            use sed_overlandFlowDetachment
            implicit none
            real, intent(in) :: alpha, B, D, g, v, rhos, rho, h, S, C
            real :: tau, tau_dc, Ws

            tau_dc = alpha* critShearStress(rho, g, h, S, D, v, rhos)

            if (tau_dc == 0.) then
                rateOfDeposition = 0.
            else
                tau = flowShearStress(rho, g, h, S)
                if (tau < tau_dc) then
                    Ws = sedSettlVelo(D, g, v, rhos, rho)
                    rateOfDeposition = B*Ws*C*rhos*(1-(tau/tau_dc))
                else
                    rateOfDeposition = 0.
                end if
            end if

        end function rateOfDeposition

        real function depthOfDeposition(alpha, B, D, g, v, rhos, rho, h, S, C, Dt)
        ! Depth of sediment deposited
        ! Input:
        ! - alpha: Relationship of critical shear stresses for desposition and ini. of motion
        ! - Channel width (m)
        ! - rho: Water density (kg/m^3)
        ! - gravi: Gravity acce. (m/s^2)
        ! - h: Water depth (m)
        ! - S: Water surface slope ~ bottom slope
        ! - D: Sediment diameter (m)
        ! - v: Kinematic viscosity of water (m^2/s)
        ! - rhos: Sediment density (kg/m^3)
        ! - C: Concentration of fine particles in suspension (m3/m3)
        ! Output:
        ! - depthOfDeposition: Depth of sediment deposited (m)
            implicit none
            real, intent(in) :: alpha, B, D, g, v, rhos, rho, h, S, C, Dt
            real :: Ed

            Ed = rateOfDeposition(alpha, B, D, g, v, rhos, rho, h, S, C)
            depthOfDeposition  = Ed * Dt/(rhos * B)

        end function  depthOfDeposition


end module sed_fineSedDepoAndResus
