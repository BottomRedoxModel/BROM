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

module fabm_niva_brom_nitrogen
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_nitrogen
    !all descriptions are in the initialize subroutine
    !state variables dependencies
    type(type_state_variable_id):: id_NH4,id_NO2,id_NO3
    type(type_state_variable_id):: id_O2
    type(type_state_variable_id):: id_POML,id_POMR,id_DOML,id_DOMR
    type(type_state_variable_id):: id_DIC,id_Alk
    type(type_state_variable_id):: id_PO4

!    type(type_diagnostic_variable_id):: id_DcPM_NOX,id_DcDM_NOX
    type(type_diagnostic_variable_id):: id_anammox!,id_DcPOMR_NOX, id_DcDOMR_NOX
    type(type_diagnostic_variable_id):: id_Nitrif1, id_Nitrif2
    type(type_diagnostic_variable_id):: id_DcPOML_NO3,id_DcDOML_NO3
    type(type_diagnostic_variable_id):: id_DcPOML_NO2,id_DcDOML_NO2
    type(type_diagnostic_variable_id):: id_DcPOMR_NO3,id_DcDOMR_NO3
    type(type_diagnostic_variable_id):: id_DcPOMR_NO2,id_DcDOMR_NO2
    type(type_diagnostic_variable_id):: id_DcTOM_NOX
    type(type_diagnostic_variable_id):: id_Denitr1,id_Denitr2
    !Model parameters
    !specific rates of biogeochemical processes
    !----N--------!
    real(rk):: K_nitrif1,K_nitrif2,O2s_nf
    real(rk):: K_annamox,O2s_dn
    real(rk):: K_omno_no3,K_omno_no2
    real(rk):: K_POML_NO3,K_POMR_NO3,K_POML_NO2,K_POMR_NO2
    real(rk):: K_DOML_NO3,K_DOMR_NO3,K_DOML_NO2,K_DOMR_NO2
    !---- Stoichiometric coefficients ----!
    real(rk):: r_c_n,r_n_p,r_n_no3,r_n_no2,s_OM_refr
  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_nitrogen), intent(inout), target :: self
    integer,                         intent(in)            :: configunit

    !----Model parameters----!
    !----N---------!
    call self%get_parameter(&
         self%K_nitrif1, 'K_nitrif1', '[1/day]',&
         'Spec.rate of 1st st. of nitrification',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_nitrif2, 'K_nitrif2', '[1/day]',&
         'Spec.rate of 2d st. of nitrification',&
         default=0.1_rk)
    call self%get_parameter(&
         self%K_POML_NO3, 'K_POML_NO3', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of POML',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_POML_NO2, 'K_POML_NO2', '[1/day]',&
         'Spec.rate of 2 stage of denitrif of POML',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_POMR_NO3, 'K_POMR_NO3', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of POMR',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_POMR_NO2, 'K_POMR_NO2', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of POMR',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_DOML_NO3, 'K_DOML_NO3', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of DOML',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_DOML_NO2, 'K_DOML_NO2', '[1/day]',&
         'Spec.rate of 2 stage of denitrif of DOML',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_DOMR_NO3, 'K_DOMR_NO3', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of DOMR',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_DOMR_NO2, 'K_DOMR_NO2', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of DOMR',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_omno_no3, 'K_omno_no3', '[uM N]',&
         'half sat. of no3 for OM denitr.',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_omno_no2, 'K_omno_no2', '[uM N]',&
         'half sat. of no2 for OM denitr.',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_annamox, 'K_annamox', '[1/day]',&
         'Spec.rate of Anammox',&
         default=0.8_rk)
    call self%get_parameter(&
         self%O2s_nf, 'O2s_nf', '[uM O]',&
         'half saturation for nitrification',&
         default=4.488_rk)
    call self%get_parameter(&
         self%O2s_dn, 'O2s_dn', '[uM O]',&
         'half saturation for denitrification',&
         default=10.0_rk)
    call self%get_parameter(&
         self%s_OM_refr, 's_OM_refr', '[uM N]',&
         'threshold of decay of refractory OM',&
         default=5.0_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(self%r_c_n,  'r_c_n',  '[-]','C[uM]/N[uM]',  default=8.0_rk)
    call self%get_parameter(self%r_n_p,  'r_n_p',  '[-]','N[uM]/P[uM]',  default=16.0_rk)
    call self%get_parameter(self%r_n_no3,'r_n_no3','[-]','N[uM]/NO3[uM]',default=0.075_rk)
    call self%get_parameter(self%r_n_no2,'r_n_no2','[-]','N[uM]/NO2[uM]',default=0.113_rk)

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3 N',&
         'ammonium')
    call self%register_state_dependency(&
         self%id_NO2,'NO2','mmol/m**3 N',&
         'nitrite')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3 N',&
         'nitrate')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3 C',&
         'total dissolved inorganic carbon',required=.false.)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3 O2',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_POML,'POML','mmol/m**3 N',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_POMR,'POMR','mmol/m**3 N',&
         'POM refractory')
    call self%register_state_dependency(&
         self%id_DOMR,'DOMR','mmol/m**3 N',&
         'DOM refractory')
    call self%register_state_dependency(&
         self%id_DOML,'DOML','mmol/m**3 N',&
         'dissolved organic nitrogen')
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3 P',&
         'phosphate',required=.false.)

    !Register diagnostic variables
    !call self%register_diagnostic_variable(self%id_DcPM_NOX,'DcPM_NOX','mmol/m**3 N',&
    !     'POML denitrification (1+2 stage)',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_DcPOMR_NOX,'DcPOMR_NOX','mmol/m**3 N',&
    !     'POMR denitrification (1+2 stage)',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_DcDOMR_NOX,'DcDOMR_NOX','mmol/m**3 N',&
    !     'DOMR denitrification (1+2 stage)',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_DcDM_NOX,'DcDM_NOX','mmol/m**3 N',&
    !     'DOML denitrification (1+2 stage)',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Nitrif1,'Nitrif1','mmol/m**3 N',&
         'Nitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Nitrif2,'Nitrif2','mmol/m**3 N',&
         'Nitrification 2 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_anammox,'Anammox','mmol/m**3 N',&
         'Anammox',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPOML_NO3,'DcPOML_NO3','mmol/m**3 N',&
         'POML denitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDOML_NO3,'DcDOML_NO3','mmol/m**3 N',&
         'DOML denitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPOML_NO2,'DcPOML_NO2','mmol/m**3 N',&
         'POML denitrification 2 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDOML_NO2,'DcDOML_NO2','mmol/m**3 N',&
         'DOML denitrification 2 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPOMR_NO3,'DcPOMR_NO3','mmol/m**3 N',&
         'POML denitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDOMR_NO3,'DcDOMR_NO3','mmol/m**3 N',&
         'DOML denitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcPOMR_NO2,'DcPOMR_NO2','mmol/m**3 N',&
         'POML denitrification 2 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDOMR_NO2,'DcDOMR_NO2','mmol/m**3 N',&
         'DOML denitrification 2 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcTOM_NOX,'DcTOM_NOX','mmol/m**3 N',&
         'Total OM denitrification all stages',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1,'Denitr1','mmol/m**3 N',&
         '(POML+DOML) denitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2,'Denitr2','mmol/m**3 N',&
         '(POML+DOML) denitrification 2 stage',output=output_time_step_integrated)

    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !check alkalinity
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_nitrogen),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: O2,POML,POMR,DOML,DOMR
    real(rk):: NO2,NO3,NH4
    !increments
    real(rk):: d_O2,d_DOML,d_POML,d_POMR,d_DOMR
    real(rk):: d_NO2,d_NO3,d_NH4
    real(rk):: d_Alk,d_DIC
    real(rk):: d_PO4
    !processes
!    real(rk):: DcPM_NOX,DcDM_NOX,DcPOMR_NOX,DcDOMR_NOX
    real(rk):: Nitrif1,Nitrif2,Anammox
    real(rk):: DcPOML_NO3,DcDOML_NO3,DcPOML_NO2,DcDOML_NO2
    real(rk):: DcPOMR_NO3,DcDOMR_NO3,DcPOMR_NO2,DcDOMR_NO2
    real(rk):: Denitr1,Denitr2,DcTOM_NOX

    _LOOP_BEGIN_
      !Retrieve current variable values
      !state
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_DOMR,DOMR)
      _GET_(self%id_NO2,NO2)
      _GET_(self%id_NO3,NO3)
      !solids
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_NH4,NH4)

      !N
      !Nitrification 1st stage: NH4+ + 1.5 O2 ->
      !                         NO2- + 2H+ + H2O (Canfield,2005)
      Nitrif1 = self%K_nitrif1*NH4*o2*0.5_rk*(1._rk+tanh(o2-self%O2s_nf))
      !Nitrification 2d stage: NO2- + 0.5 O2 -> NO3- (Canfield,2005)
      Nitrif2 = self%K_nitrif2*NO2*o2*0.5_rk*(1._rk+tanh(o2-self%O2s_nf))
      !in suboxic conditions
      !Anammox NO2- + NH4+ -> N2 + 2H2O (Canfield,2005)
      Anammox = self%K_annamox*NO2*NH4 &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn)))
      !OM denitrification (Richards, 1965)
      !(CH2O)106(NH3)16H3PO4 + 84.8HNO3 =
      ! 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4
      !POM and DOM denitrification (1st stage) (Anderson,1982)
      !1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2
      DcPOML_NO3 = self%K_POML_NO3*POML &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO3/(NO3+self%K_omno_no3)
      DcDOML_NO3 = self%K_DOML_NO3*DOML &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO3/(NO3+self%K_omno_no3)
      DcPOMR_NO3 = self%K_POMR_NO3*POMR &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
                  *NO3/(NO3+self%K_omno_no3)
      DcDOMR_NO3 = self%K_DOMR_NO3*DOMR &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
                  *NO3/(NO3+self%K_omno_no3)
      !POM and DOM denitrification (2d stage)
      !3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,1982)
      DcPOML_NO2 = self%K_POML_NO2*POML &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO2/(NO2+self%K_omno_no2)
      DcDOML_NO2 = self%K_DOML_NO2*DOML &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO2/(NO2+self%K_omno_no2)
      DcPOMR_NO2 = self%K_POMR_NO2*POMR &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
                  *NO2/(NO2+self%K_omno_no2)
      DcDOMR_NO2 = self%K_DOMR_NO2*DOMR &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
                  *NO2/(NO2+self%K_omno_no2)
      !Denitrification as consumpation of NOX and production of N2
      Denitr1 = self%r_n_no3*(DcPOMR_NO3+DcDOMR_NO3) 
      Denitr2 = self%r_n_no2*(DcPOMR_NO2+DcDOMR_NO2)
      !Summariazed OM mineralization in N (released NH4) units
      DcTOM_NOX = DcPOMR_NO3+DcPOMR_NO2+DcDOMR_NO3+DcDOMR_NO2

      !Set increments
   d_O2 = -1.5_rk*Nitrif1-0.5_rk*Nitrif2
      _SET_ODE_(self%id_O2,d_O2)
   d_DOML = -DcDOML_NO3-DcDOML_NO2
      _SET_ODE_(self%id_DOML,d_DOML)
   d_DOMR = DcDOML_NO3+DcDOML_NO2-DcDOMR_NO3-DcDOMR_NO2
      _SET_ODE_(self%id_DOMR,d_DOMR)
   d_POML = -DcPOML_NO3-DcPOML_NO2
      _SET_ODE_(self%id_POML,d_POML)
   d_POMR = DcPOML_NO3+DcPOML_NO2-DcPOMR_NO3-DcPOMR_NO2
      _SET_ODE_(self%id_POMR,d_POMR)
   d_NO2 = Nitrif1-Nitrif2-Anammox+self%r_n_no3*(DcDOMR_NO3+DcPOMR_NO3) &
          -self%r_n_no2*(DcPOMR_NO2+DcDOMR_NO2)
      _SET_ODE_(self%id_NO2,d_NO2)
   d_NO3 = Nitrif2-DcDOMR_NO3-self%r_n_no3*(DcDOMR_NO3+DcPOMR_NO3)
      _SET_ODE_(self%id_NO3,d_NO3)
   d_NH4 = DcDOML_NO3+DcDOML_NO2+DcPOML_NO3+DcPOML_NO2-Nitrif1-Anammox
      _SET_ODE_(self%id_NH4,d_NH4)
   d_DIC = (DcDOMR_NO3+DcDOMR_NO2+DcPOMR_NO3+DcPOMR_NO2)*self%r_c_n
      _SET_ODE_(self%id_DIC,d_DIC)
   d_PO4 = (DcDOML_NO3+DcDOML_NO2+DcPOML_NO3+DcPOML_NO2)/self%r_n_p
      _SET_ODE_(self%id_PO4,d_PO4)
   d_Alk = (&      !Alkalinity changes due to redox reactions:
             !NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O
             -2._rk*Nitrif1 &  !(Wolf-Gladrow, Zeebe, 2007)
             !3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 or
             !5 CH2O + 4 H+ + 4 NO3- -> 2 N2 + 5 CO2 + 7H2O
!             +1._rk*(DcPOML_NO2+DcDOML_NO2) &
!             + DcTOM_NOX &
             )
      _SET_ODE_(self%id_Alk,d_Alk)

      _SET_DIAGNOSTIC_(self%id_anammox,Anammox)
      _SET_DIAGNOSTIC_(self%id_DcPOML_NO3,DcPOML_NO3)
      _SET_DIAGNOSTIC_(self%id_DcDOML_NO3,DcDOML_NO3)
      _SET_DIAGNOSTIC_(self%id_DcPOML_NO2,DcPOML_NO2)
      _SET_DIAGNOSTIC_(self%id_DcDOML_NO2,DcDOML_NO2)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_NO3,DcPOMR_NO3)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_NO3,DcDOMR_NO3)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_NO2,DcPOMR_NO2)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_NO2,DcDOMR_NO2)
      _SET_DIAGNOSTIC_(self%id_DcTOM_NOX,DcTOM_NOX)
      _SET_DIAGNOSTIC_(self%id_Denitr1,Denitr1)
      _SET_DIAGNOSTIC_(self%id_Denitr2,Denitr2)
      _SET_DIAGNOSTIC_(self%id_Nitrif1,Nitrif1)
      _SET_DIAGNOSTIC_(self%id_Nitrif2,Nitrif2)

    _LOOP_END_
  end subroutine do
end module fabm_niva_brom_nitrogen
