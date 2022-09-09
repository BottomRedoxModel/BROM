!-----------------------------------------------------------------------
! brom_injection is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Anfisa Berezina
!-----------------------------------------------------------------------
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:
   module fabm_niva_brom_injection
!
! !DESCRIPTION:
!
! brom_injection parameterizes chemical/ biogeochemical transfomration of an injected matter.
! brom_injection consists of 1 state variables ( in nitrogen units):
! - waste - is a compound that can be can be decayed and /or oxidized to inorganic nutrients 
!   with a strong oxidizing agent under acidic condition (oxygen).
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Anfisa Berezina
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_injection
!     Variable identifiers
! variables defined in this modules
      type (type_state_variable_id)        :: id_waste
      type (type_state_variable_id)        :: id_oxy, id_nut, id_pom, id_dom
! variables defined in this modules
      type (type_state_variable_id)        :: id_po4, id_nh4, id_dic
! global dependencies
      type (type_dependency_id)            :: id_par, id_temp, id_salt, id_depth, id_thickness
      type (type_horizontal_dependency_id) :: id_long
      type (type_diagnostic_variable_id)   :: id_waste_diss, id_waste_miner, id_waste_decomp
!     Model parameters
       real(rk) :: r_waste_diss,r_waste_miner,r_waste_decomp
       real(rk) :: Wwaste, Bu, beta_da, Tda
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the Brom_Injection model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_injection), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
    real(rk),parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
  ! waste
   call self%get_parameter(self%r_waste_decomp, 'r_waste_decomp', '1/d', &
           'Specific rate of waste decomposition',default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_waste_diss, 'r_waste_diss', '1/d', &
           'Specific rate of waste decay',default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_waste_miner, 'r_waste_miner', '1/d', &
           'Specific rate of waste mineralization',default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Wwaste,         'Wwaste',        'm/s', &
           'vertical velocity of Waste (<0 for sinking)', default=-1.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Bu,             'Bu',            'nd',  &
           'Burial coeficient for lower boundary',      default=0.25_rk)
   call self%get_parameter(self%beta_da,        'beta_da',         'nd', 'Coefficient for dependence of mineralization on t ', default=20._rk)
   call self%get_parameter(self%Tda,            'Tda',             'nd', 'Coefficient for dependence of mineralization on t ', default=13._rk)
   ! Register state variables
   call self%register_state_variable(self%id_waste,'Waste','mmol/m**3', &
           'Waste   compound ', 0.0_rk, minimum=0.0_rk, vertical_movement=self%Wwaste)
   ! Register link to external variables
   call self%register_state_dependency(self%id_oxy,'Oxy','mmol/m**3','OXY')
   call self%register_state_dependency(self%id_pom,'POM','mmol/m**3','POM')
   call self%register_state_dependency(self%id_nut,'NUT','mmol/m**3','NUT')
   call self%register_state_dependency(self%id_dom,'DOM','mmol/m**3','DOM')
   !call self%register_state_dependency(self%id_po4,'PO4','mmol/m**3','PO4')
   !call self%register_state_dependency(self%id_nh4,'NH4','mmol/m**3','NH4')
   !call self%register_state_dependency(self%id_dic,'DIC','mmol/m**3','DIC')
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_waste_decomp,'waste_decomp', &
           'mmol/m**3/d',  'waste_decomp,  Decomposition of waste',&
            output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_waste_diss,'waste_diss', &
           'mmol/m**3/d',  'waste_diss,  Decay of waste',&
            output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_waste_miner,'waste_miner', &
           'mmol/m**3/d',  'waste_miner,  Mineralization of waste with oxygen',&
            output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
!   call self%register_dependency(self%id_depth,standard_variables%depth)
!   call self%register_dependency(self%id_thickness,standard_variables%cell_thickness)
!   call self%register_dependency(self%id_long, standard_variables%longitude)
   
  !!!       call self%register_dependency(self%id_bedstress,standard_variables%bottom_stress)
!     call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
!    call self%register_dependency(self%id_pres,standard_variables%pressure)
!   call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
!   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_injection),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) ::  waste,  oxy,  pom,  nut, po4, nh4, dom, dic
   real(rk) :: dwaste, dpom, ddom
 ! external parameters   
   real(rk) :: t, depth,long, thickness
 ! Rates of biogeochemical processes
   real(rk) :: waste_decomp           ! decay mineralization (1/d)
   real(rk) :: waste_diss            ! decay mineralization (1/d)
   real(rk) :: waste_miner            ! oxic mineralization of waste (1/d)
!   real(rk) :: waste_diss_denitr        ! suboxic mineralization of waste (denitrification) (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_waste,waste)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_pom,pom)
   _GET_(self%id_nut,nut)
   _GET_(self%id_dom,dom)
   !_GET_(self%id_po4,po4)
   !_GET_(self%id_nh4,nh4)
   !_GET_(self%id_dic,dic)

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,t)                ! temperature
!!   _GET_(self%id_depth,depth)           ! depth, m
!!   _GET_HORIZONTAL_(self%id_long,long)  ! longitude for "dx", m
!!   _GET_(self%id_thickness,thickness)   ! cell thickness, "Hz", m

!--------------------------------------------------------------
! Oxic mineralization of waste depends on T
!!!!!!!!!!!!!!   waste_decomp = self%r_waste_decomp*waste
   waste_decomp = self%r_waste_decomp*waste*(1.0_rk+self%beta_da*yy(self%tda,t))
   waste_miner = self%r_waste_miner*waste*(1.0_rk+self%beta_da*yy(self%tda,t))
   waste_diss = self%r_waste_diss*waste*(1.0_rk+self%beta_da*yy(self%tda,t))
!   waste_diss = self%r_waste_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*waste
! Suboxic mineralization depends on T,O2,NO3/NO2
   !waste_diss_denitr = self%r_pom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
   !                        * (0.5-0.5*tanh(self%O2LimC-oxy)) &
   !                        * (1-tanh(1.-nut))*pom
! Mineralization of OM, ammonification and growth of NUT

!if (waste.gt.0.01) write (*,*) waste, waste_decomp
! Now we can summarize processes and write state variables sink/sources:
!--------------------------------------------------------------
! OXY
!--------------------------------------------------------------
! Changes of OXY due to OM production and decay!
   dwaste = -waste_decomp-waste_miner-waste_diss
   dpom   =  waste_decomp
   ddom   =  waste_diss
!--------------------------------------------------------------

!derivatives for FABM
   _SET_ODE_(self%id_waste,dwaste)
   _SET_ODE_(self%id_pom,  dpom)
   _SET_ODE_(self%id_dom,  ddom)
   _SET_ODE_(self%id_nut,  waste_miner) ! case OXYDEP
!   _SET_ODE_(self%id_nh4,  waste_miner)! case FABM
!   _SET_ODE_(self%id_po4,  waste_miner/16.0_rk)! case FABM
!   _SET_ODE_(self%id_dic,  waste_miner*6.625_rk)! case FABM
   _SET_ODE_(self%id_oxy,  -6.625*waste_miner)


   ! Export diagnostic variables
_SET_DIAGNOSTIC_(self%id_waste_decomp,waste_decomp)
_SET_DIAGNOSTIC_(self%id_waste_miner,waste_miner)
_SET_DIAGNOSTIC_(self%id_waste_diss,waste_diss)
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!

! !INPUT PARAMETERS:
   class (type_niva_brom_injection),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: waste

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_waste,waste)

   ! BURYING into the sediments, mmol/m2/s (sinking rates "Wxxx" are in m/s and positive upward)
   _SET_BOTTOM_EXCHANGE_(self%id_waste,self%Bu*self%Wwaste*waste)

   _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------
  !
  ! Original author(s): Hans Burchard, Karsten Bolding
  ! DESCRIPTION:
  ! This is a squared Michaelis-Menten type of limiter
  !
  real(rk) function yy(a,x)
    real(rk),intent(in):: a,x

    yy=x**2._rk/(a**2._rk+x**2._rk)
  end function yy
end module
