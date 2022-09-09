!-----------------------------------------------------------------------
! BROM2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM2 distribution.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_silicon
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_silicon
    !all descriptions are in initialize subroutine
    !variables allocated here  
    type(type_state_variable_id)       :: id_Sipart
    !state variables dependencies
    type(type_state_variable_id)       :: id_Si
    !diagnostic variables
    type (type_diagnostic_variable_id) :: id_Si_dissolution, id_Si_clay_miner, id_Si_precip
    !Model parameters
    !sinking
    real(rk):: Wsed
    !specific rates of biogeochemical processes
    !----Si---------!
    real(rk):: K_sipart_diss, K_sipart_diss_limit
    real(rk):: K_sipart_to_minerals, Si_diss_max
  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_silicon), intent(inout), target :: self
    integer,                        intent(in)            :: configunit

    !-----Model parameters------
    !Sinking
    call self%get_parameter(&
         self%Wsed,'Wsed','[1/day]',&
         'Rate of sinking of detritus (POP, POML)',&
         default=5.00_rk)

    !Specific rates of biogeochemical processes
    !----Si-------!
    call self%get_parameter(&
         self%K_sipart_diss, 'K_sipart_diss', '[1/day]',&
         'Si dissollution rate constant',&
         default=0.10_rk)
    
    call self%get_parameter(&
         self%K_sipart_diss_limit, 'K_sipart_diss_limit', '[uM]',&
         ' Si particulate maximum concentration, that is biogenic and can be dissolved',&
         default=30.0_rk)
    call self%get_parameter(&
         self%K_sipart_to_minerals, 'K_sipart_to_minerals',&
         '[1/day]','Si_part transformation into minerals not modeled here',&
         default=0.10_rk)
    call self%get_parameter(&
         self%Si_diss_max, 'Si_diss_max',&
         'mmol/m**3','Max conc. of diss Si for precipitation',&
         default=1500.0_rk)    

    !Register state variables
    call self%register_state_variable(&
         self%id_Sipart, 'Sipart', 'mmol/m**3','Si Particulate',&
         minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
    !Register state dependencies
    call self%register_state_dependency(&
         self%id_Si, 'Si', 'mmol/m**3','Si')
    !Register diagnostics
    call self%register_diagnostic_variable(&
         self%id_Si_dissolution,'Si_dissolution','mmol/m**3/d',&
         'Biogeonic particulate Si dissolution',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Si_clay_miner,'Si_clay_miner','mmol/m**3/d',&
         'Si_part transformation into minerals not modeled here',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Si_precip,'Si_precip','mmol/m**3/d',&
         'Precipitation of dissolved Si at high concentrations',&
         output=output_time_step_integrated)
    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_silicon),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: Sipart, Si
    !increments
    real(rk):: d_Si,d_Sipart
    !processes
    real(rk):: Si_dissolution, Si_clay_miner, Si_precip

    _LOOP_BEGIN_
      !Retrieve current variable values
      !state
      _GET_(self%id_Si,Si)
      _GET_(self%id_Sipart,Sipart)
! biogeonic silicate dissolution possible for small concentrations of Sipart (i,e, biogenic)
      Si_dissolution = self%K_sipart_diss*Sipart*(1._rk-0.5_rk*(1._rk+tanh(Sipart-self%K_sipart_diss_limit)))
! Formation of minerals with Al etc. (DeMaster, 2003, Treatise of Geochemistry, vol.7)
    if (Sipart.gt.1000.0_rk) then
      Si_clay_miner = self%K_sipart_to_minerals*(Sipart-1000.0_rk)
    else
      Si_clay_miner = 0.0_rk
    endif
! Precipitation of dissolved Si at high concentrations (Strakhov, 1978)
    if (Si.gt.self%Si_diss_max) then
      Si_precip= 0.1*(Si-self%Si_diss_max)
    else
      Si_precip= 0.0_rk
    endif
      !Set increments
      d_Si = Si_dissolution-Si_precip
      d_Sipart = -Si_dissolution-Si_clay_miner+Si_precip
      
    _SET_ODE_(self%id_Sipart,d_Sipart)
    _SET_ODE_(self%id_Si,d_Si)

    _SET_DIAGNOSTIC_(self%id_Si_dissolution,Si_dissolution)
    _SET_DIAGNOSTIC_(self%id_Si_precip,Si_precip)
    _SET_DIAGNOSTIC_(self%id_Si_clay_miner,Si_clay_miner)

    _LOOP_END_
  end subroutine do
end module fabm_niva_brom_silicon
