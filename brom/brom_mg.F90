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

module fabm_niva_brom_mg
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_mg
    !all descriptions are in initialize subroutine
    !variables allocated here
    type(type_state_variable_id)       :: id_MgOH2
    type(type_state_variable_id)       :: id_Mg, id_Alk
    !diagnostic dependencies
    type(type_dependency_id):: id_Hplus, id_salt
    !diagnostic variables
    type (type_diagnostic_variable_id) :: id_MgOH2_diss,  id_MgOH2_form, id_Om_MgOH2, id_MgOH2_mg_l

    !Model parameters
    !Sinking
    real(rk):: Wmg, Wmg_min, Wmg_max !, rate_diss_mult
    !specific rates of biogeochemical processes
    !----Mg---------!
    real(rk):: K_MgOH2_diss, K_MgOH2_form, K_MgOH2, MgOH2_diss_limit

  contains
    procedure :: initialize
    procedure :: do
    procedure :: get_vertical_movement
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_mg), intent(inout), target :: self
    integer,                        intent(in)            :: configunit

    !-----Model parameters------
    !Sinking
    call self%get_parameter(&
         self%Wmg,'Wmg','[1/day]',&
         'rate of sinking of particulate matter',&
         default=16.0_rk)

    call self%get_parameter(&
         self%Wmg_min,'Wmg_min','[1/day]',&
         'Min rate of sinking of particulate matter',&
         default=1.00_rk)

    call self%get_parameter(&
         self%Wmg_max,'Wmg_max','[1/day]',&
         'Max rate of sinking of particulate matter',&
         default=90.00_rk)

    !Specific rates of biogeochemical processes
    !----Mg-------!
    call self%get_parameter(&
         self%K_MgOH2,'K_MgOH2','[uM]',&
         'MgOH2 equilibrium constant (Solubility Product Constant)',&
         default=2510.0_rk)
    call self%get_parameter(&
         self%K_MgOH2_diss, 'K_MgOH2_diss', '[1/day]',&
         'MgOH2 dissollution rate constant',&
         default=0.10_rk)
    call self%get_parameter(&
         self%K_MgOH2_form, 'K_MgOH2_form', '[1/day]',&
         'MgOH2 formation rate constant',&
         default=0.10_rk)
    call self%get_parameter(&
         self%MgOH2_diss_limit, 'MgOH2_diss_limit', '[uM]',&
         ' Mg particulate maximum concentration',&
         default=30.0_rk)
!    call self%get_parameter(&
!         self%rate_diss_mult,'rate_diss_mult','[1/day]',&
!         'Multiplier for dissolution rate',&
!         default=1.00_rk)


    !Register state variables
    call self%register_state_variable(&
         self%id_MgOH2, 'MgOH2', 'mmol/m**3','Mg(OH)2',&
         minimum=0.0_rk,vertical_movement=-self%Wmg/86400._rk)
    call self%register_state_variable(&
         self%id_Mg, 'Mg', 'mmol/m**3','Mg dissolved additional',&
         minimum=0.0_rk,vertical_movement=0._rk)

    !Register state dependencies
    call self%register_dependency(self%id_Hplus,&
         'Hplus', 'mmol/m**3','H+ Hydrogen')
    call self%register_dependency(&
         self%id_salt,standard_variables%practical_salinity)

    !Register diagnostics
    call self%register_diagnostic_variable(&
         self%id_MgOH2_diss,'Diss_MgOH2','mmol/m**3/d',&
         'MgOH2 dissolution',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MgOH2_form,'Form_MgOH2','mmol/m**3/d',&
         'MgOH2 formation ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Om_MgOH2,'Om_MgOH2','n/d',&
         'Om solubility of MgOH2 ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MgOH2_mg_l,'MgOH2_mg_l','mg/l',&
         'Mg(OH)2 in mg/l ',&
         output=output_time_step_integrated)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_mg),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: MgOH2, Mg
    !diagnostic variables dependencies
    real(rk):: Hplus
    !increments
    real(rk):: d_Mg, d_MgOH2, d_Alk
    !processes
    real(rk):: MgOH2_diss, MgOH2_form,  Om_MgOH2
    real(rk):: OHminus
    !standard variables
    real(rk):: salt !,temp,pres

    _LOOP_BEGIN_
      !Retrieve current variable values
      !state
      _GET_(self%id_Mg,Mg)   !this is an "additional" Mg originated from added Mg(OH)2
      _GET_(self%id_MgOH2,MgOH2)
      _GET_(self%id_Hplus,Hplus)
!      _GET_(self%id_temp,temp)        ! temperature
      _GET_(self%id_salt,salt)        ! salinity


!!Mg(OH)2 (formation)/dissollution  Mg(OH)2 (s) <=> Mg2+ + 2 OH-
      OHminus = 10._rk**(-14._rk)/Hplus
!     Mg =  !from salt: 1 psu = 37.53 mg Mg/l
!     Mg = (53000.0_r/35.0rk)*salt !(Hawaii presentation, 2018) !uM
!     Mg = (0.054/35)*S ! in 1 l of SW=1.3 g =/24=0.054  mol
!  Solubility constant for MgOH2 precipitation/dissolution
      Om_MgOH2 = (((0.054_rk/35.0_rk)*salt+Mg/1000000._rk)*OHminus*OHminus)/(self%K_MgOH2) ! Mg converted into M from mmol/m3
!  MgOH2 dissolution  Mg(OH)2(s) -> Mg2+ + 2OH-   umol/d
      MgOH2_diss = max(24.0_rk,self%K_MgOH2_diss*83.346_rk*(MgOH2**(-0.524_rk))) &
                       *MgOH2*max(0._rk,(1._rk-Om_MgOH2)) !&
 !                      *self%rate_diss_mult
          if (MgOH2.lt.self%MgOH2_diss_limit) MgOH2_diss = 0.0_rk

!  MgOH2 formation    Mg2+ + 2OH- -> Mg(OH)2(s)
      MgOH2_form = self%K_MgOH2_form*max(0._rk,(Om_MgOH2-1._rk))
 !   endif
      !Set increments
      d_Mg    =  MgOH2_diss-MgOH2_form
      d_MgOH2 = -MgOH2_diss+MgOH2_form
      d_Alk = (&      !Alkalinity changes due to redox reactions:
             +2._rk*MgOH2_diss &  ! Mg(OH)2(s) -> Mg2+ + 2OH-
             -2._rk*MgOH2_form &  ! Mg2+ + 2OH- -> Mg(OH)2(s)
             )
    _SET_ODE_(self%id_MgOH2,d_MgOH2)
    _SET_ODE_(self%id_Mg,d_Mg)
    _SET_ODE_(self%id_Alk,d_Alk)

    _SET_DIAGNOSTIC_(self%id_MgOH2_diss,MgOH2_diss)
    _SET_DIAGNOSTIC_(self%id_MgOH2_form,MgOH2_form)
    _SET_DIAGNOSTIC_(self%id_Om_MgOH2,Om_MgOH2)
    _SET_DIAGNOSTIC_(self%id_MgOH2_mg_l,MgOH2*0.0171_rk)



    _LOOP_END_
  end subroutine do
    !
  !
  ! Increased manganese sinking via MnIV and MnIII oxides formation
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
     class (type_niva_brom_mg), intent(in) :: self
     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

     real(rk) :: Wmg
     real(rk) :: MgOH2
     real(rk) :: ranval

     _LOOP_BEGIN_
      _GET_(self%id_MgOH2,MgOH2)
        call random_number(ranval)   !as in duplicator.f90
  !    value = minimum + ranval*(maximum-minimum)
  !    call random(ranval)
         Wmg = -(self%Wmg_min + (self%Wmg_max-self%Wmg_min)*ranval)

      _ADD_VERTICAL_VELOCITY_(self%id_MgOH2, Wmg)

     _LOOP_END_

  end subroutine get_vertical_movement
end module

