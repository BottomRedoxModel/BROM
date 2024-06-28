!-----------------------------------------------------------------------
! BROM2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM2 distribution.
!-----------------------------------------------------------------------
!This is the module for the Nutrients state variables initialization

#include "fabm_driver.h"

module fabm_niva_brom_main_nutrients
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_main_nutrients
    !all descriptions are in the initialize subroutine
    !state variables
    type(type_state_variable_id):: id_NH4,id_NO2,id_NO3
    type(type_state_variable_id):: id_Si
    type(type_state_variable_id):: id_PO4
  contains
    procedure :: initialize
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_main_nutrients),intent(inout),target:: self
    integer,                              intent(in)          :: configunit

    !Register state variables
    call self%register_state_variable(&
         self%id_NH4,'NH4','mmol/m**3 N',&
         'ammonium',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_NO2,'NO2','mmol/m**3 N',&
         'nitrite',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_NO3,'NO3','mmol/m**3 N',&
         'nitrate',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_PO4,'PO4','mmol/m**3 P',&
         'phosphate',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_Si,'Si','mmol/m**3 Si',&
         'silicon',minimum=0.0_rk)
  end subroutine initialize
end module fabm_niva_brom_main_nutrients
