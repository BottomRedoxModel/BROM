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

module fabm_niva_brom_methane
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_methane
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_CH4
    !state dependencies
    type(type_state_variable_id):: id_DIC,id_NH4,id_PO4
    type(type_state_variable_id):: id_DOML,id_POML,id_POMR,id_DOMR
    type(type_state_variable_id):: id_O2,id_NO3,id_SO4
    type(type_state_variable_id):: id_oxy,id_nut,id_pom,id_dom
    !global dependencies
    type (type_dependency_id)        :: id_salt
    !diagnostic variables by bacteria needed
    type(type_diagnostic_variable_id):: id_DcPOML_ch4,id_DcDOML_ch4
    type(type_diagnostic_variable_id):: id_DcPOMR_ch4,id_DcDOMR_ch4
    type(type_diagnostic_variable_id):: id_DcTOM_CH4, id_ch4_o2, id_ch4_so4
    
    logical   :: bgcmod_BROM ! added ice area limitation of air-sea flux, with new switch parameter use_aice.
    
    !Model parameters
    !specific rates of biogeochemical processes
    !Methane
    real(rk):: s_omso_o2,s_omso_no3,s_omch_so4,s_OM_refr
    real(rk):: K_DOML_ch4,K_POML_ch4,K_POMR_ch4,K_DOMR_ch4
    real(rk):: K_ch4_o2,K_ch4_so4
    !---- Stoichiometric coefficients ----!
    real(rk):: r_c_n, r_n_p

  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class(type_niva_brom_methane),intent(inout),target :: self
    integer,                      intent(in)           :: configunit

    !-----Model parameters------
    !Specific rates of biogeochemical processes
    !Methane
    call self%get_parameter(&
         self%s_omso_o2, 's_omso_o2', '[uM O]',&
         'threshold of o2 for OM sulfate reduction', default=25.0_rk)
    call self%get_parameter(&
         self%s_omso_no3, 's_omso_no3', '[uM N]',&
         'threshold of noX for OM sulfate reduction', default=5.0_rk)
    call self%get_parameter(&
         self%s_omch_so4, 's_omch_so4', '[uM S]',&
         'threshold of SO4 for CH4 production from OM', default=15000.0_rk)
    call self%get_parameter(&
         self%s_OM_refr, 's_OM_refr', '[uM N]',&
         'threshold of decay of refractory OM', default=5.0_rk)
    call self%get_parameter(&
         self%K_DOML_ch4, 'K_DOML_ch4', '[1/day]',&
         'Specific rate of CH4 production from DOML', default=0.00014_rk)
    call self%get_parameter(&
         self%K_POML_ch4, 'K_POML_ch4', '[1/day]',&
         'Specific rate of CH4 production from POML', default=0.00014_rk)
    call self%get_parameter(&
         self%K_POMR_ch4, 'K_POMR_ch4', '[1/day]',&
         'Specific rate of CH4 production from POMR', default=0.00014_rk)
    call self%get_parameter(&
         self%K_DOMR_ch4, 'K_DOMR_ch4', '[1/day]',&
         'Specific rate of CH4 production from DOMR', default=0.00014_rk)
    call self%get_parameter(&
         self%K_ch4_o2, 'K_ch4_o2', '[1/day]',&
         'Specific rate of oxidation of CH4 with O2', default=0.14_rk)
    call self%get_parameter(&
         self%K_ch4_so4, 'K_ch4_so4', '[1/day]',&
         'Specific rate of oxidation of CH4 with SO4', default=0.0000274_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(&
         self%r_c_n,   'r_c_n',  '[-]','C[uM]/N[uM]',default=6.625_rk)
    call self%get_parameter(&
         self%r_n_p,   'r_n_p',  '[-]','N[uM]/P[uM]',default=16.0_rk)
    !register state variables
    call self%register_state_variable(&
         self%id_CH4,'CH4','mmol/m**3','CH4',minimum=0.0_rk)
    !register state dependencies
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3','DIC')
    call self%register_state_dependency(&
         self%id_O2,'O2','mmol/m**3','dissolved oxygen')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3','nitrate')
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3','phosphate',required=.false.)
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3','ammonium',required=.false.)
    call self%register_state_dependency(&
         self%id_POML,'POML','mmol/m**3','particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_POMR,'POMR','mmol/m**3',&
         'particulate organic nitrogen',required=.false.)
    call self%register_state_dependency(&
         self%id_DOMR,'DOMR','mmol/m**3','DOMR',required=.false.)
    call self%register_state_dependency(&
         self%id_DOML,'DOML','mmol/m**3','dissolved organic nitrogen')
    call self%register_state_dependency(&
         self%id_SO4,'SO4','mmol/m**3','sulphate',required=.false.)
    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPOML_ch4,'DcPOML_ch4','mmol/m**3',&
         'CH4 production from POML',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_ch4,'DcDOML_ch4','mmol/m**3',&
         'CH4 production from DOML',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_ch4,'DcPOMR_ch4','mmol/m**3',&
         'CH4 production from POMR ',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_ch4,'DcDOMR_ch4','mmol/m**3',&
         'CH4 production from DOMR ',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcTOM_CH4,'DcTOM_CH4','mmol/m**3',&
         'Total OM mineralization with CH4 genesis',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ch4_o2,'ch4_o2','mmol/m**3/d',&
         'CH4 oxidation with oxygen ',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ch4_so4,'ch4_so4','mmol/m**3/d',&
         'CH4 anaeribic oxidation with SO4 ',output=output_time_step_integrated)
    ! Register environmental dependencies
    call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
    !Specify that rates are per day (default: per second)
    self%dt = 86400._rk
  end subroutine initialize
!
!
!
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_methane),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: DIC,SO4,NO3,NH4,PO4
    real(rk):: O2,CH4
    real(rk):: POML,POMR,DOML,DOMR
    !processes
    real(rk):: DcDOML_ch4,DcPOML_ch4,DcPOMR_CH4,DcDOMR_CH4
    real(rk):: DcTOM_CH4, ch4_o2,ch4_so4
    !LOCAL VARIABLES:
    real(rk):: salt
    !increments
    real(rk):: d_SO4,d_O2,d_CH4,d_NH4,d_DIC,d_PO4
    real(rk):: d_DOML,d_POML,d_POMR,d_DOMR

    _LOOP_BEGIN_
    
      !salinity      
      _GET_(self%id_salt,salt)         
      
      !state variables    
      _GET_(self%id_DOML,DOML) ! DOM in Oxydep
      _GET_(self%id_NO3,NO3)  ! NUT in Oxydep
      
  if (_AVAILABLE_(self%id_DOMR)) then      
      _GET_(self%id_DOMR,DOMR) 
  else 
      DOMR=0.0_rk
  endif
  
  if (_AVAILABLE_(self%id_SO4)) then
      _GET_(self%id_SO4,SO4) 
  else
      SO4 = salt*799.7_rk   ! convert from g Salt to g SO4 and then to uM : (2.649_rk/34.483_rk)/96.061_rk*10^6 
  endif 
  
  if (_AVAILABLE_(self%id_NH4)) then
      _GET_(self%id_NH4,NH4)
  else
      NH4 =0.0_rk
  endif
  
  if (_AVAILABLE_(self%id_PO4)) then
      _GET_(self%id_PO4,PO4)
  else
      PO4=0.0_rk
  endif
  
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_CH4,CH4)
      !solids
      _GET_(self%id_POML,POML)
      
  if (_AVAILABLE_(self%id_POMR)) then      
      _GET_(self%id_POMR,POMR)
  else
      POMR=0.0_rk
  endif    

  
      !CH4 production from DOML, POML, DOMR and POMR 
      !(CH2O)106(NH3)16H3PO4 -> 53 CO2 + 53 CH4 + 16 NH3 + H3PO4
      DcDOML_ch4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                *(1._rk-0.5_rk*(1._rk+tanh(SO4-self%s_omch_so4))) &
                * self%K_DOML_ch4*DOML
      DcPOML_ch4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                *(1._rk-0.5_rk*(1._rk+tanh(SO4-self%s_omch_so4))) &
                * self%K_POML_ch4*POML
      DcDOMR_CH4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                *(1._rk-0.5_rk*(1._rk+tanh(SO4-self%s_omch_so4))) &
                *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
                * self%K_DOMR_ch4*DOMR
      DcPOMR_CH4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                *(1._rk-0.5_rk*(1._rk+tanh(SO4-self%s_omch_so4))) &
                *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
                * self%K_POMR_ch4*POMR
      !total OM mineeralization 
      DcTOM_CH4=DcDOMR_CH4+DcPOMR_CH4
      !
      !CH4 oxidation with O2
      !CH4 + 2O2 = CO2 + 2H2O
      ch4_o2 = self%K_ch4_o2*CH4*O2
      !
      !CH4 anoxic oxidation with SO4
      !CH4 + SO42- + 2 H+  = CO2 + H2S + 2H2O
      ch4_so4 = self%K_ch4_so4*CH4*SO4

      !Set increments
   d_DOML = -DcDOML_ch4
      _SET_ODE_(self%id_DOML,d_DOML)
   d_POML = -DcPOML_ch4
      _SET_ODE_(self%id_POML,d_POML)
   d_DOMR = DcDOML_ch4-DcDOMR_ch4
      if (_AVAILABLE_(self%id_DOMR)) then 
          _SET_ODE_(self%id_DOMR,d_DOMR)
      endif
   d_POMR = DcPOML_ch4-DcPOMR_ch4
      if (_AVAILABLE_(self%id_POMR)) then 
         _SET_ODE_(self%id_POMR,d_POMR)
      endif
   d_O2 = -2._rk*ch4_o2
      _SET_ODE_(self%id_O2,d_O2)
   d_SO4 = -ch4_so4
     if (_AVAILABLE_(self%id_SO4)) then    
      _SET_ODE_(self%id_SO4,d_SO4)
     endif
   d_DIC = (DcDOMR_ch4+DcPOMR_ch4)*self%r_c_n+ch4_o2+ch4_so4
      _SET_ODE_(self%id_DIC,d_DIC)
   d_CH4 = 0.5_rk*(DcDOMR_ch4+DcPOMR_ch4)-ch4_o2-ch4_so4
      _SET_ODE_(self%id_CH4,d_CH4)      
   d_NH4 = DcDOML_ch4+DcPOML_ch4
     if (_AVAILABLE_(self%id_NH4)) then 
         _SET_ODE_(self%id_NH4,d_NH4)
     endif
   d_PO4 = (DcDOML_ch4+DcPOML_ch4)/self%r_n_p
     if (_AVAILABLE_(self%id_PO4)) then 
          _SET_ODE_(self%id_PO4,d_PO4)
     endif
      
      _SET_DIAGNOSTIC_(self%id_DcPOML_ch4,DcPOML_ch4)
      _SET_DIAGNOSTIC_(self%id_DcDOML_ch4,DcDOML_ch4)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_ch4,DcPOMR_ch4)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_ch4,DcDOMR_ch4)
      _SET_DIAGNOSTIC_(self%id_DcTOM_CH4,DcTOM_CH4)
      _SET_DIAGNOSTIC_(self%id_ch4_o2, ch4_o2)
      _SET_DIAGNOSTIC_(self%id_ch4_so4,ch4_so4)
      
    _LOOP_END_
  end subroutine do
end module
