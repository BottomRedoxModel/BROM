!-----------------------------------------------------------------------
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Jorn Bruggeman
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_volumes
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_volumes
!     Variable identifiers
      type (type_state_variable_id)        :: id_V_air,id_V_wat,id_V_sed
      type (type_state_variable_id)        :: id_hz_m,id_z_m, id_area

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
   class (type_niva_brom_volumes), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!

   call self%register_state_variable(self%id_V_air,   'V_air', 'm**3','V air',   minimum=0.0_rk)
   call self%register_state_variable(self%id_V_wat,   'V_wat', 'm**3','V water',   minimum=0.0_rk)
   call self%register_state_variable(self%id_V_sed,   'V_sed', 'm**3','V sediment',   minimum=0.0_rk)
   call self%register_state_variable(self%id_hz_m,    'hz_m',  'm','hz, layer thickeness',   minimum=0.0_rk)
   call self%register_state_variable(self%id_z_m,     'z_m',   'm**3','z_m, mid-layer depth',   minimum=0.0_rk)
   call self%register_state_variable(self%id_area,    'area',  'm**2','area total of the layer',   minimum=0.0_rk)

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
   class (type_niva_brom_volumes),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) :: V_air, V_wat, V_sed
   real(rk) ::  hz_m, z_m, area
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _GET_(self%id_V_air,V_air)
   _GET_(self%id_V_wat,V_wat)
   _GET_(self%id_V_sed,V_sed)
   _GET_(self%id_hz_m,hz_m)
   _GET_(self%id_z_m,z_m)
   _GET_(self%id_area,area)
   
   _SET_ODE_(self%id_V_air,0.0_rk)
   _SET_ODE_(self%id_V_wat,0.0_rk)
   _SET_ODE_(self%id_V_sed,0.0_rk)
   _SET_ODE_(self%id_hz_m, 0.0_rk)
   _SET_ODE_(self%id_z_m,  0.0_rk)
   _SET_ODE_(self%id_area, 0.0_rk)
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC


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
!
! !INPUT PARAMETERS:
   class (type_niva_brom_volumes),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

end module
