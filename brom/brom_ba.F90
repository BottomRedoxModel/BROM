#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_ba
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman, Svetlana Pakhomova, Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_ba
!     Variable identifiers
      type (type_state_variable_id)        :: id_Ba, id_BaSO4, id_SO4      
      
      !---- Ba---------!
      real(rk) :: K_BaSO4=5.    ! BaSO4 equilibrium constant (Solubility Product Constant) (uM)=5  ( 5 umol/l, Wiki,09)  
      real(rk) :: K_BaSO4_form=1.4e-6 ! Specific rate of precipitation of BaSO4 from Ba with SO4 (1/day)=1.4e-6 (5x10-4 uM/yr,  Arndt,09)
      real(rk) :: K_BaSO4_diss=8.e-11 ! Specific rate of dissollution of BaSO4 to Ba and SO4  (1/day)=8.e-11 (3x10-8 1/yr, Arndt,09)
!     Model parameters
      real(rk) :: Wsed= 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)      
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
! !IROUTINE: Initialise the BROM equilibrium constant model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! This module described transformation of Ba, barium 
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_ba), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev, Svetlana Pakhomova
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Parameters, i.e. rate constants  
   call self%get_parameter(self%Wsed, 'Wsed', '[m/day]',  'Rate of sinking of detritus (POP, PON)',       default=5.00_rk)     
   call self%register_state_variable(self%id_Ba, 'Ba', 'mol/m**3','barium', minimum=0.0_rk)
   call self%register_state_variable(self%id_BaSO4, 'BaSO4', 'mol/m**3','barium sulphate', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   
     call self%register_state_dependency(self%id_SO4, 'SO4', 'mmol/m**3','SO4')

     self%dt = 86400.
   
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
   class (type_niva_brom_ba),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev 
!
! !LOCAL VARIABLES:
   real(rk) :: Ba, BaSO4 
   real(rk) :: SO4
   
   real(rk) ::  Om_BaSO4, baso4_prec, baso4_diss,ba_flux  
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _GET_(self%id_Ba,Ba) 
   _GET_(self%id_BaSO4,BaSO4)    
   _GET_(self%id_SO4, SO4)

!BaSO¤ formation/dissollution (Arndt,09)
          Om_BaSO4=SO4*Ba/(self%K_BaSO4) 
!% BaSO4 formation Ba2+ + SO42- -> BaSO4  (Arndt,09)
    baso4_prec = 0.0_rk !self%K_BaSO4_form*max(0._rk,(Om_BaSO4-1._rk))
!% BaSO4 dissollution  BaSO4  -> Ba2+ + SO42-   (Arndt,09)
    baso4_diss = 0.0_rk !self%K_BaSO4_diss*BaSO4*max(0._rk,(1._rk-Om_BaSO4))       

    _SET_ODE_(self%id_Ba,baso4_diss-baso4_prec+ba_flux)
    _SET_ODE_(self%id_BaSO4,-baso4_diss+baso4_prec)   
    _SET_ODE_(self%id_SO4,baso4_diss-baso4_prec)    

! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

end module