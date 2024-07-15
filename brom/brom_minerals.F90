!-----------------------------------------------------------------------
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_minerals
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_minerals
!     Variable identifiers
      type (type_state_variable_id)        :: id_Mg, id_K, id_Br, id_Na, id_Cl, id_Ca_u, id_CO3, id_SO4
      type (type_state_variable_id)        :: id_Na2CaSO42, id_K2Ca2MgSO44, id_CaMgCO32, id_CaSO4
      type (type_state_variable_id)        :: id_KCl, id_MgSO4
      
!      type (type_bottom_state_variable_id) :: id_NaCl_bot
      type (type_dependency_id)            :: id_temp
            
!      type (type_diagnostic_variable_id)            :: id_psu
      type (type_diagnostic_variable_id)            :: id_Na2CaSO42_sat
      type (type_diagnostic_variable_id)            :: id_K2Ca2MgSO44_sat
      type (type_diagnostic_variable_id)            :: id_CaMgCO32_sat
      type (type_diagnostic_variable_id)            :: id_CaSO4_sat
      type (type_diagnostic_variable_id)            :: id_MgSO4_sat
      type (type_diagnostic_variable_id)            :: id_KCl_sat
              real(rk) :: K_CaSO4
              real(rk) :: K_CaSO4_diss
              real(rk) :: K_CaSO4_form
              real(rk) :: K_MgSO4
              real(rk) :: K_MgSO4_diss
              real(rk) :: K_MgSO4_form
              real(rk) :: K_KCl
              real(rk) :: K_KCl_diss
              real(rk) :: K_KCl_form
              real(rk) :: K_Na2CaSO42
              real(rk) :: K_Na2CaSO42_diss
              real(rk) :: K_Na2CaSO42_form
      real(rk) :: w_NaCl      
      real(rk) :: porosity

   contains
      procedure :: initialize
      procedure :: do
 !     procedure :: do_bottom
   end type

!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the BROM equilibrium constant model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_minerals), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !!!!!!!!!!!!!!!!!!!!!!!!!   37.3 mol**2/L**2 =37.3*10e6*10e6 mmol**2/m**2
   call self%get_parameter(self%K_CaSO4,     'K_CaSO4',     'mmol**2/m**6', 'CaSO4 (gypsum) conditional equilibrium constant',default=2.6e5_rk) 
   call self%get_parameter(self%K_CaSO4_diss,'K_CaSO4_diss','mmol/m**3/d',  'Specific rate of dissolution of CaSO4',default=1._rk)
   call self%get_parameter(self%K_CaSO4_form,'K_CaSO4_form','mmol/m**3/d',  'Specific rate of formation of CaSO4',default=1._rk)
   call self%get_parameter(self%K_MgSO4,     'K_MgSO4',     'mmol**2/m**6', 'MgSO4 (?) conditional equilibrium constant',default=2.6e5_rk) 
   call self%get_parameter(self%K_MgSO4_diss,'K_MgSO4_diss','mmol/m**3/d',  'Specific rate of dissolution of MgSO4',default=1._rk)
   call self%get_parameter(self%K_MgSO4_form,'K_MgSO4_form','mmol/m**3/d',  'Specific rate of formation of MgSO4',default=1._rk)
   call self%get_parameter(self%K_KCl,     'K_KCl',     'mmol**2/m**6', 'KCl (sylvite) conditional equilibrium constant',default=2.6e5_rk) 
   call self%get_parameter(self%K_KCl_diss,'K_KCl_diss','mmol/m**3/d',  'Specific rate of dissolution of KCl',default=1._rk)
   call self%get_parameter(self%K_KCl_form,'K_KCl_form','mmol/m**3/d',  'Specific rate of formation of KCl',default=1._rk)
   call self%get_parameter(self%K_Na2CaSO42,     'K_Na2CaSO42',     'mmol**2/m**6', 'Na2CaSO42 (glauberite) conditional equilibrium constant',default=2.6e5_rk) 
   call self%get_parameter(self%K_Na2CaSO42_diss,'K_Na2CaSO42_diss','mmol/m**3/d',  'Specific rate of dissolution of Na2CaSO42')
   call self%get_parameter(self%K_Na2CaSO42_form,'K_Na2CaSO42_form','mmol/m**3/d',  'Specific rate of formation of Na2CaSO42')
   call self%get_parameter(self%w_NaCl,     'w_NaCl',     'm/d',          'Sedimentation rate for particulate NaCl',default=3._rk)

   call self%register_state_variable(self%id_Ca_u,   'Ca_u',   'mmol/m**3','Ca_u', minimum=0.0_rk)
   call self%register_state_variable(self%id_Mg,   'Mg',   'mmol/m**3','Mg', minimum=0.0_rk)
   call self%register_state_variable(self%id_K,    'K',    'mmol/m**3','K',  minimum=0.0_rk)
   call self%register_state_variable(self%id_Br,   'Br',   'mmol/m**3','Br', minimum=0.0_rk)
   call self%register_state_variable(self%id_Na2CaSO42,  'Na2CaSO42',  'mmol/m**3','Na2CaSO42 glauberite', minimum=0.0_rk,vertical_movement=-self%w_NaCl/86400._rk)
!   call self%register_state_variable(self%id_K2Ca2MgSO44,'K2Ca2MgSO44','mmol/m**3','K2Ca2MgSO44 polyhalite',minimum=0.0_rk,vertical_movement=-self%w_NaCl/86400._rk)
!   call self%register_state_variable(self%id_CaMgCO32,   'CaMgCO32',   'mmol/m**3','CaMgCO32 dolomite',    minimum=0.0_rk,vertical_movement=-self%w_NaCl/86400._rk)
   call self%register_state_variable(self%id_CaSO4,      'CaSO4',      'mmol/m**3','CaSO4 gypsum',         minimum=0.0_rk,vertical_movement=-self%w_NaCl/86400._rk)
   call self%register_state_variable(self%id_MgSO4,      'MgSO4',      'mmol/m**3','MgSO4 ',         minimum=0.0_rk,vertical_movement=-self%w_NaCl/86400._rk)
   call self%register_state_variable(self%id_KCl,        'KCl',        'mmol/m**3','KCl sylvite',          minimum=0.0_rk,vertical_movement=-self%w_NaCl/86400._rk)

   call self%register_state_dependency(self%id_Cl,  'Cl',  'mmol/m**3','Cl')
   call self%register_state_dependency(self%id_Na,  'Na',  'mmol/m**3','Na')
   call self%register_state_dependency(self%id_SO4,  'SO4',  'mmol/m**3','SO4')

!  call self%register_diagnostic_variable(self%id_psu,     'salt' ,   '-','salinity in PSU', standard_variable=standard_variables%practical_salinity)
   call self%register_diagnostic_variable(self%id_CaSO4_sat,'CaSO4_sat','-','CaSO4_saturation')
   call self%register_diagnostic_variable(self%id_MgSO4_sat,'MgSO4_sat','-','MgSO4_saturation')
   call self%register_diagnostic_variable(self%id_KCl_sat,'KCl_sat','-','KCl_saturation')
   call self%register_diagnostic_variable(self%id_Na2CaSO42_sat,'Na2CaSO42_sat','-','Na2CaSO42_saturation')

   call self%register_dependency(self%id_temp,standard_variables%temperature)
!Specify that rates are per day (default: per second)
    self%dt = 86400._rk
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
   class (type_niva_brom_minerals),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) :: temp, Na, Mg, Ca_u, K, Cl, SO4
   real(rk) :: CaSO4, Om_CaSO4,CaSO4_prec,CaSO4_diss
   real(rk) :: MgSO4, Om_MgSO4,MgSO4_prec,MgSO4_diss
   real(rk) :: KCl,   Om_KCl,  KCl_prec,  KCl_diss
   real(rk) :: Na2CaSO42,   Om_Na2CaSO42,  Na2CaSO42_prec,  Na2CaSO42_diss
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature - not used yet

   _GET_(self%id_Ca_u,Ca_u)
   _GET_(self%id_Mg,Mg)
   _GET_(self%id_K,K)
   _GET_(self%id_Na,Na)
   _GET_(self%id_Cl,Cl)
   _GET_(self%id_SO4,SO4)
   _GET_(self%id_CaSO4,CaSO4)
   _GET_(self%id_MgSO4,MgSO4)
   _GET_(self%id_KCl,KCl)
   _GET_(self%id_Na2CaSO42,Na2CaSO42)

! CaSO4 saturation state (NB all units mmol/m**3)
    Om_CaSO4=Ca_u*SO4/(self%K_CaSO4)*f_t_Ksp(temp)
   _SET_DIAGNOSTIC_(self%id_CaSO4_sat,Om_CaSO4)
   CaSO4_prec=max(0.0_rk,self%K_CaSO4_form*max(0._rk,(Om_CaSO4-1._rk)))
   CaSO4_diss=max(0.0_rk,self%K_CaSO4_diss*max(0._rk,(1._rk-Om_CaSO4))*CaSO4)
   
! MgSO4 saturation state (NB all units mmol/m**3)
    Om_MgSO4=Mg*SO4/(self%K_MgSO4)*f_t_Ksp(temp)
   _SET_DIAGNOSTIC_(self%id_MgSO4_sat,Om_MgSO4)
   MgSO4_prec=max(0.0_rk,self%K_MgSO4_form*max(0._rk,(Om_MgSO4-1._rk)))
   MgSO4_diss=max(0.0_rk,self%K_MgSO4_diss*max(0._rk,(1._rk-Om_MgSO4))*MgSO4)

! KCl saturation state (NB all units mmol/m**3)
    Om_KCl=K*Cl/(self%K_KCl)
   _SET_DIAGNOSTIC_(self%id_KCl_sat,Om_KCl)*f_t_Ksp(temp)
   KCl_prec=max(0.0_rk,self%K_KCl_form*max(0._rk,(Om_KCl-1._rk)))
   KCl_diss=max(0.0_rk,self%K_KCl_diss*max(0._rk,(1._rk-Om_KCl))*KCl)

! Na2CaSO42 glauberite. saturation state (NB all units mmol/m**3)
!      CaSO4 +2Na + SO4 -> Na2Ca(SO4)2
    Om_Na2CaSO42=Na*Na*CaSO4*SO4/(self%K_Na2CaSO42)*f_t_Ksp(temp)
   _SET_DIAGNOSTIC_(self%id_Na2CaSO42_sat,Om_Na2CaSO42)
   Na2CaSO42_prec=max(0.0_rk,self%K_Na2CaSO42_form*max(0._rk,(Om_Na2CaSO42-1._rk)))
   Na2CaSO42_diss=max(0.0_rk,self%K_Na2CaSO42_diss*max(0._rk,(1._rk-Om_Na2CaSO42))*Na2CaSO42)

! resulting changes
   _SET_ODE_(self%id_Na,   -2*Na2CaSO42_prec+2*Na2CaSO42_diss)
   _SET_ODE_(self%id_Ca_u,   -CaSO4_prec+CaSO4_diss)
   _SET_ODE_(self%id_Mg,   -MgSO4_prec+MgSO4_diss)
   _SET_ODE_(self%id_SO4,  -CaSO4_prec+CaSO4_diss-Na2CaSO42_prec+Na2CaSO42_diss-MgSO4_prec+MgSO4_diss)
   _SET_ODE_(self%id_K,    -KCl_prec+KCl_diss)
   _SET_ODE_(self%id_Cl,   -KCl_prec+KCl_diss)

   _SET_ODE_(self%id_CaSO4,+CaSO4_prec-CaSO4_diss-Na2CaSO42_prec+Na2CaSO42_diss)
   _SET_ODE_(self%id_MgSO4,+MgSO4_prec-MgSO4_diss)
   _SET_ODE_(self%id_KCl,  +KCl_prec-KCl_diss)   
   _SET_ODE_(self%id_Na2CaSO42, +Na2CaSO42_prec-Na2CaSO42_diss)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
!   Dependence of solubility product on temperature  !
  elemental real(rk) function f_t_Ksp(temperature) 
    real(rk),intent(in):: temperature

    f_t_Ksp = exp(2.86623-571.361/(temperature+272.15) &
           +77.128/(temperature+272.15)/(temperature+272.15))/2.57

  end function f_t_Ksp
!-----------------------------------------------------------------------
end module
