!-----------------------------------------------------------------------
! brom_acetate is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Elizaveta Protsenko, Alisa Ilinskaya
!-----------------------------------------------------------------------
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:
   module fabm_niva_brom_acetate
!
! !DESCRIPTION:
!
! brom_acetate parameterize  the Chemical Oxygen Demand, COD (https://en.wikipedia.org/wiki/Chemical_oxygen_demand) 
! OM mineralization, nitrification, and oxidation of reduced specied of S, Mn, Fe, present in suboxic conditions.
! brom_acetate consists of 1 state variables ( in oxygen-units):
! - acetate - is an organic compound CnHaObNc that can be can be fully oxidized to inorganic nutrients 
!   with a strong oxidizing agent under acidic condition (oxygen).
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): 
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_acetate
!     Variable identifiers
      type (type_state_variable_id)        :: id_oxy, id_dom, id_ac_tot
      type (type_dependency_id)            :: id_temp, id_salt
      type (type_diagnostic_variable_id)   :: id_ac_free_miner_rate, id_mg, id_ac_free, id_ac_Mg
!     Model parameters
      !---Organic matter mineralization---- !
       real(rk) :: r_ac_free_miner, Tda, beta_da, Wacetate, Bu, mg_s_ratio, Mg_ac_const
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
! !IROUTINE: Initialise the OXYDEP-COD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_acetate), intent(inout), target :: self
   integer,                      intent(in)            :: configunit
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
  ! acetate
   call self%get_parameter(self%mg_s_ratio, 'mg_s_ratio', '-', 'mg to s ratio',               default=0.10_rk)
   call self%get_parameter(self%r_ac_free_miner, 'r_ac_free_miner', '1/d', 'Specific rate of acetate mineralization',           default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%beta_da,        'beta_da',         'nd', 'Coefficient for dependence of mineralization on t ', default=20._rk)
   call self%get_parameter(self%Tda,            'Tda',             'nd', 'Coefficient for dependence of mineralization on t ', default=13._rk)
   call self%get_parameter(self%Wacetate,      'Wacetate',          'm/s', 'vertical velocity of acetate (<0 for sinking)',         default=-1.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Bu,             'Bu',             'nd',  'Burial coeficient for lower boundary',               default=0.25_rk)
   call self%get_parameter(self%Mg_ac_const,    'Mg_ac_const',      '-',  'Stability constant for Mg+acetate complex',       default=1.25_rk)
   ! Register state variables
!   call self%register_state_variable(self%id_ac_free,'ac_free','mmol/m**3','Acetate free', 0.0_rk, minimum=0.0_rk, vertical_movement=self%Wacetate)
   call self%register_state_variable(self%id_ac_tot,'ac_tot','mmol/m**3','Acetate total', 0.0_rk, minimum=0.0_rk, vertical_movement=self%Wacetate)
!   call self%register_state_variable(self%id_ac_Mg,'ac_Mg','mmol/m**3','Acetate+Mg', 0.0_rk, minimum=0.0_rk, vertical_movement=self%Wacetate)
   ! Register link to external variables
   call self%register_state_dependency(self%id_oxy,'Oxy','mmol/m**3','OXY')
!   call self%register_state_dependency(self%id_oxy,'DOM','mmol/m**3','DOM')
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_ac_free_miner_rate,'ac_free_miner_rate','mmol/m**3/d',  'ac_free_miner_rate,  Mineralization of acetate with oxygen',           &
                    output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_mg,'Mg','mmol/m**3',  'Mg, concentration', output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_ac_free,'ac_free','mmol/m**3',  'Acetate free', output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_ac_Mg,'ac_Mg','mmol/m**3',  'Acetate+Mg', output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
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
   class (type_niva_brom_acetate),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: ac_free, ac_tot, ac_Mg, oxy, dom, t, s

   real(rk) :: doxy, dac_free, dac_tot, dac_Mg
 ! Rates of biogeochemical processes
   real(rk) :: ac_free_miner_rate   ! oxic mineralization of acetate (1/d)
   real(rk) :: mg_ac_const   ! oxic mineralization of acetate (1/d)
   real(rk) :: id_ac_free           ! concentration of mineralized acetate
   real(rk) :: id_ac_mg             ! concentration of Mg+acetate complex
   real(rk) :: id_ac_tot            ! total concentration of acetate
   real(rk) :: Mg            ! Mg (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_ac_tot,ac_tot)

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,t)              ! temperature
   _GET_(self%id_salt,s)              ! salinity
!    Mg concentration
   Mg=self%mg_s_ratio*s

!--------------------------------------------------------------
!   Mg+acetate complex concentration 
!  we assume that ac_Mg= Mg*Ac_free*K and Ac_free = Ac_tot-Ac_Mg, so:
   ac_Mg= (self%mg_ac_const*Mg*ac_tot)/(1+self%mg_ac_const*Mg)
   ac_free= ac_tot - ac_Mg 

!--------------------------------------------------------------
! Oxic mineralization of acetate depends on T
   ac_free_miner_rate   = self%r_ac_free_miner*(1.+self%beta_da*yy(self%tda,t))*ac_free

!--------------------------------------------------------------
! Now we can summarize processes and write state variables sink/sources:
   doxy     = -2.0_rk*ac_free_miner_rate
   dac_tot  = -ac_free_miner_rate

!-------------------------------------------------------------
!derivatives for FABM
   _SET_ODE_(self%id_oxy, doxy)
   _SET_ODE_(self%id_ac_tot,dac_tot)

! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_ac_free_miner_rate,ac_free_miner_rate)   
   _SET_DIAGNOSTIC_(self%id_Mg,Mg)
   _SET_DIAGNOSTIC_(self%id_ac_free,ac_free)
   _SET_DIAGNOSTIC_(self%id_ac_Mg,ac_Mg)
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
   class (type_niva_brom_acetate),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: ac_tot

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_ac_tot,ac_tot)

   ! BURYING into the sediments, mmol/m2/s (sinking rates "Wxxx" are in m/s and positive upward)
   _SET_BOTTOM_EXCHANGE_(self%id_ac_tot,self%Bu*self%Wacetate*ac_tot)

   _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Saturation function squared
!
! !INTERFACE:
   real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !IN2PUT PARAMETERS:
   real(rk), intent(in)                :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)

   end function yy
!EOC

   end module fabm_niva_brom_acetate

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
