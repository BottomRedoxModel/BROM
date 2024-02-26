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

module fabm_niva_brom_sulfur
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_sulfur
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_H2S
    type(type_state_variable_id):: id_S0,id_S2O3,id_SO4
    !state variables dependencies
    type(type_state_variable_id):: id_O2,id_NO3,id_NH4,id_PO4
    type(type_state_variable_id):: id_POML,id_POMR,id_DOML,id_DOMR
    type(type_state_variable_id):: id_DIC,id_Alk

    type(type_diagnostic_variable_id):: id_DcPOMR_S2O3,id_DcDOMR_S2O3,id_DcPOML_s2o3,id_DcDOML_s2o3
    type(type_diagnostic_variable_id):: id_DcPOMR_SO4, id_DcDOMR_SO4, id_DcPOML_so4, id_DcDOML_so4
    type(type_diagnostic_variable_id):: id_s2o3_no3,id_s0_no3,id_DcTOM_SOX
    type(type_diagnostic_variable_id):: id_s0_ox,id_s2o3_ox
    type(type_diagnostic_variable_id):: id_s0_disp,id_hs_ox
    type(type_diagnostic_variable_id):: id_hs_no3,id_s2o3_rd,id_so4_rd
    
    ! diagnostic dependencies
    type(type_dependency_id):: id_Wadd
    
    !Model parameters
    !specific rates of biogeochemical processes
    real(rk):: K_s0_disp,K_hs_ox,K_s0_ox,K_s0_no3
    real(rk):: K_POMR_so4,K_POMR_s2o3,K_POML_so4,K_POML_s2o3
    real(rk):: K_DOMR_so4,K_DOMR_s2o3,K_DOML_so4,K_DOML_s2o3
    real(rk):: K_s2o3_ox,K_s2o3_no3,K_hs_no3
    !---- Switches-------!
    real(rk):: s_omso_o2,s_omso_no3,s_OM_refr
    !---- Stoichiometric coefficients ----!
    real(rk):: r_c_n,r_n_p,r_n_s
    !sinking
    real(rk):: WS0 , WS0_tot
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
    class (type_niva_brom_sulfur), intent(inout), target :: self
    integer,                       intent(in)            :: configunit

    !-----Model parameters------
    !---- S---------!
    call self%get_parameter(&
         self%K_s0_disp,'K_s0_disp','[1/day]',&
         'Specific rate of S0 dispropotionation',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_hs_ox,'K_hs_ox','[1/day]',&
         'Specific rate of oxidation of H2S to S0 with O2',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_s0_ox,'K_s0_ox','[1/day]',&
         'Specific rate of oxidation of S0 with O2',&
         default=0.02_rk)
    call self%get_parameter(&
         self%K_s0_no3,'K_s0_no3','[1/day]',&
         'Specific rate of oxidation of S0 with NO3',&
         default=0.9_rk)
    call self%get_parameter(&
         self%K_s2o3_ox,'K_s2o3_ox','[1/day]',&
         'Specific rate of oxidation of S2O3 with O2',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_s2o3_no3,'K_s2o3_no3','[1/day]',&
         'Specific rate of oxidation of S2O3 with NO3',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_hs_no3, 'K_hs_no3', '[1/day]',&
         'Spec.rate of thiodenitrification.',&
         default=0.8_rk)
    call self%get_parameter(&
         self%K_POML_so4,'K_POML_so4','[1/day]',&
         'Specific rate of OM sulfate reduction with sulfate',&
         default=0.000005_rk)
    call self%get_parameter(&
         self%K_POMR_so4,'K_POMR_so4','[1/day]',&
         'Specific rate of POMR sulfate reduction with sulfate',&
         default=0.000005_rk)
    call self%get_parameter(&
         self%K_POML_s2o3,'K_POML_s2o3','[1/day]',&
         'Specific rate of OM sulfate reduction with thiosulfate',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_POMR_s2o3,'K_POMR_s2o3','[1/day]',&
         'Specific rate of POMR sulfate reduction with sulfate',&
         default=0.000005_rk)
    call self%get_parameter(&
         self%K_DOML_so4,'K_DOML_so4','[1/day]',&
         'Specific rate of DOML sulfate reduction with sulfate',&
         default=0.000005_rk)
    call self%get_parameter(&
         self%K_DOMR_so4,'K_DOMR_so4','[1/day]',&
         'Specific rate of DOMR sulfate reduction with sulfate',&
         default=0.000005_rk)
    call self%get_parameter(&
         self%K_DOML_s2o3,'K_DOML_s2o3','[1/day]',&
         'Specific rate of DOML sulfate reduction with thiosulfate',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_DOMR_s2o3,'K_DOMR_s2o3','[1/day]',&
         'Specific rate of DOMR sulfate reduction with thiosulfate',&
         default=0.000005_rk)
    !----Switches----!
    call self%get_parameter(&
         self%s_omso_o2, 's_omso_o2', '[uM O]',&
         'threshold of o2 for OM sulfate reduction',&
         default=25.0_rk)
    call self%get_parameter(&
         self%s_omso_no3, 's_omso_no3', '[uM N]',&
         'threshold of noX for OM sulfate reduction',&
         default=5.0_rk)
    call self%get_parameter(&
         self%s_OM_refr, 's_OM_refr', '[uM N]',&
         'threshold of decay of refractory OM',&
         default=5.0_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(&
         self%r_c_n,   'r_c_n',  '[-]',&
         'C[uM]/N[uM]',&
         default=6.625_rk)
    call self%get_parameter(&
         self%r_n_p,   'r_n_p',  '[-]',&
         'N[uM]/P[uM]',&
         default=16.0_rk)
    call self%get_parameter(&
         self%r_n_s,   'r_n_s',  '[-]',&
         'N[uM]/S[uM]',&
         default=0.302_rk)
    !----Sinking----!
    call self%get_parameter(&
         self%WS0,'WS0','1/day',&
         'S particles sinking rate',&
          default=0.0_rk)
    call self%get_parameter(&
         self%WS0_tot,'WS0_tot','[1/day]',&
         'Total accelerated sinking with absorbed Mn hydroxides',&
          default=0.0_rk)

    !Register state variables
    call self%register_state_variable(&
         self%id_H2S, 'H2S', 'mmol/m**3','H2S',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_S0 , 'S0',  'mmol/m**3','S0',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_S2O3,'S2O3','mmol/m**3','S2O3',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_SO4, 'SO4', 'mmol/m**3','SO4',&
         minimum=0.0_rk)

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3',&
         'total dissolved inorganic carbon',required=.false.)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3',&
         'phosphate',required=.false.)
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3',&
         'ammonium')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3',&
         'nitrate')
    call self%register_state_dependency(&
         self%id_POML,'POML','mmol/m**3',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_POMR,'POMR','mmol/m**3',&
         'POM refractory')
    call self%register_state_dependency(&
         self%id_DOMR,'DOMR','mmol/m**3',&
         'POM refractory')
    call self%register_state_dependency(&
         self%id_DOML,'DOML','mmol/m**3',&
         'dissolved organic nitrogen')

    !Register diagnostic variables
    !call self%register_diagnostic_variable(&
    !     self%id_DcPM_SO4,'DcPM_SO4','mmol/m**3',&
    !     'POM sulfatereduction (1+2 stage)',&
    !     output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_DcDM_SO4,'DcDM_SO4','mmol/m**3',&
    !     'DOM sulfatereduction (1+2 stage)',&
    !     output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_SO4,'DcPOMR_SO4','mmol/m**3',&
         'POMR sulfatereduction with SO4 ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_S2O3,'DcPOMR_S2O3','mmol/m**3',&
         'POMR sulfatereduction wirh S2O3',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_SO4,'DcDOMR_SO4','mmol/m**3',&
         'DOMR sulfatereduction with SO4 ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_S2O3,'DcDOMR_S2O3','mmol/m**3',&
         'DOMR sulfatereduction wirh S2O3',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOML_s2o3,'DcPOML_s2o3','mmol/m**3',&
         'POM sulfatereduction 2d stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_s2o3,'DcDOML_s2o3','mmol/m**3',&
         'DOM sulfatereduction 2d stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOML_SO4,'DcPOML_SO4','mmol/m**3',&
         'POML sulfatereduction 1st stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_SO4,'DcDOML_SO4','mmol/m**3',&
         'DOML sulfatereduction 1st stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcTOM_SOX,'DcTOM_SOX','mmol/m**3',&
         'Total OM sulfatereduction all stages',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_s2o3_no3,'s2o3_no3','mmol/m**3',&
         ' S2O3 with NO3 oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_s0_no3,'s0_no3','mmol/m**3',&
         'S0 with NO3 oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_s2o3_ox,'s2o3_ox','mmol/m**3',&
         'S2O3  with O2oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_s0_ox,'s0_ox','mmol/m**3',&
         'S0 with O2 oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_so4_rd,'so4_rd','mmol/m**3',&
         '(POM+DOM) sulfatereduction 1st stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_s2o3_rd,'s2o3_rd','mmol/m**3',&
         '(POM+DOM) sulfatereduction 2d stage ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_s0_disp,'s0_disp','mmol/m**3',&
         'S0 disproportionation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_hs_ox,'hs_ox','mmol/m**3',&
         'H2S with O2 oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_hs_no3,'hs_no3','mmol/m**3',&
         'H2S with NO3 oxidation',&
         output=output_time_step_integrated)
    !Register diagnostic dependencies
    call self%register_dependency(self%id_Wadd,'Wadd','[1/day]',&
         'Additional sinking velocity via Mn4 adsorptoin')

    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_sulfur),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: H2S,S0,S2O3,SO4
    !state dependencies
    real(rk):: O2
    real(rk):: POML,POMR,DOML,DOMR
    real(rk):: NO3,NH4,PO4
    !increments
    real(rk):: d_S2O3,d_SO4,d_S0,d_H2S
    real(rk):: d_O2,d_DOML,d_POML,d_POMR,d_DOMR
    real(rk):: d_NO3,d_NH4,d_DIC,d_PO4
    real(rk):: d_Alk
    !processess
    real(rk):: s0_disp,hs_ox,s0_ox,s0_no3,s2o3_ox,s2o3_no3,hs_no3
    real(rk):: DcPOML_so4,DcDOML_so4,DcPOML_s2o3,DcDOML_s2o3,s2o3_rd,so4_rd
    real(rk):: DcPOMR_SO4,DcPOMR_S2O3,DcDOMR_SO4,DcDOMR_S2O3
    !Summariazed OM mineralization in N units
    real(rk):: DcTOM_SOX

    _LOOP_BEGIN_
      !Retrieve current state variable values
      !state
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_DOMR,DOMR)
      _GET_(self%id_S2O3,S2O3)
      _GET_(self%id_SO4,SO4)
      _GET_(self%id_NO3,NO3)
      _GET_(self%id_NH4,NH4)
      _GET_(self%id_PO4,PO4)
      !solids
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)
      _GET_(self%id_S0,S0)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_H2S,H2S)

      !S
      !S0 disportionation: 4S0 + 3H2O -> 2H2S + S2O32= + 2H+
      s0_disp = self%K_s0_disp*S0
      !HS oxidation with O2: 2H2S + O2 -> 2S0 + 2H2O
      hs_ox = self%K_hs_ox*o2*H2S
      !S0 oxidation with O2: 2S0 + O2 + H2O -> S2O32= + 2H+
      s0_ox = self%K_s0_ox*o2*S0
      !S0 oxidation with NO3: 4S0 + 3NO3- + 7H2O -> 4SO4= + 3NH4+ + 2H+
      s0_no3 = self%K_s0_no3*NO3*S0
      !S2O3 oxidation with O2: S2O32= + 2O2 + 2OH- -> 2SO42= + H2O
      s2o3_ox = self%K_s2o3_ox*o2*S2O3
      !S2O3 oxidation with NO3: S2O3= + NO3- + 2H2O --> 2SO4= + NH4+
      s2o3_no3 = self%K_s2o3_no3*NO3*S2O3
      !Thiodenitrification: 3H2S + 4NO3- + 6OH- -> 3SO4= + 2N2 + 6H2O
      !(Volkov, 1984)
      hs_no3 = self%K_hs_no3*H2S*NO3
      !in anoxic conditions:
      !OM sulfatereduction (Boudreau, 1996)
      !(CH2O)106(NH3)16H3PO4 + 53SO42- = 106HCO3- + 16NH3 + H3PO4 + 53H2S
      !POM sulfatereduction (1st stage):
      DcPOML_so4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                 *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                 *self%K_POML_so4*SO4*POML
      !DOM sulfatereduction (1st stage):
      DcDOML_so4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                 *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                 *self%K_DOML_so4*SO4*DOML
      if (o2.gt.10._rk) then
        DcPOML_so4=0._rk
        DcDOML_so4=0._rk
      endif
      !POM sulfatereduction (2d stage):
      DcPOML_s2o3 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                  *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                  *self%K_POML_s2o3*S2O3*POML
      !DOML sulfatereduction (2d stage):
      DcDOML_s2o3 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                  *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                  *self%K_DOML_s2o3*S2O3*DOML
      if (o2.gt.10._rk) then
        DcPOML_s2o3=0._rk
        DcDOML_s2o3=0._rk
      endif
      so4_rd   = (DcPOML_so4+DcDOML_so4)/self%r_n_s  !in S units
      s2o3_rd  = (DcPOML_s2o3+DcDOML_s2o3)/self%r_n_s
      ! POMR mineralization with SO4 and S2O3
      DcPOMR_SO4  = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                  *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                  *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
                  *self%K_POMR_so4 *SO4 *POMR
      DcPOMR_S2O3 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                  *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                  *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
                  *self%K_POMR_s2o3*S2O3*POMR
      ! DOMR mineralization with SO4 and S2O3
      DcDOMR_SO4  = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                  *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                  *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
                  *self%K_DOMR_so4 *SO4 *DOMR
      DcDOMR_S2O3 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                  *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                  *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
                  *self%K_DOMR_s2o3*S2O3*DOMR
      !Summariazed total SO4 and S2O3 reduction for total OM in N units
      DcTOM_SOX=DcPOMR_SO4+DcPOMR_S2O3+DcDOMR_SO4+DcDOMR_S2O3

      !Set increments
   d_SO4 = hs_no3+0.5_rk*s2o3_ox+s0_no3+2._rk*s2o3_no3 &
           +(-DcPOMR_SO4-DcDOMR_SO4)/self%r_n_s 
      _SET_ODE_(self%id_SO4,d_SO4)
   d_S2O3 = 0.5_rk*s0_ox-s2o3_ox+0.25_rk*s0_disp-s2o3_no3 &
           +0.5_rk*((DcPOMR_SO4+DcDOMR_SO4-DcDOMR_S2O3-DcPOMR_S2O3))/self%r_n_s
      _SET_ODE_(self%id_S2O3,d_S2O3)
   d_S0 = hs_ox-s0_ox-s0_disp-s0_no3
      _SET_ODE_(self%id_S0,d_S0)
   d_H2S = -hs_ox+0.5_rk*s0_disp-hs_no3 &
           +(DcDOMR_S2O3+DcPOMR_S2O3)/self%r_n_s
      _SET_ODE_(self%id_H2S,d_H2S)
   d_O2 = -0.5_rk*hs_ox-0.5_rk*s0_ox-0.5_rk*s2o3_ox
      _SET_ODE_(self%id_O2,d_O2)
   d_DOML = -DcDOML_so4-DcDOML_s2o3
      _SET_ODE_(self%id_DOML,d_DOML)
   d_POML = -DcPOML_so4-DcPOML_s2o3
      _SET_ODE_(self%id_POML,d_POML)
   d_POMR=  DcPOML_so4+DcPOML_s2o3-DcPOMR_SO4-DcPOMR_S2O3
      _SET_ODE_(self%id_POMR,d_POMR)
   d_DOMR=  DcPOML_so4+DcPOML_s2o3-DcDOMR_SO4-DcDOMR_s2o3
      _SET_ODE_(self%id_DOMR,d_DOMR)
   d_NO3 = -1.6_rk*hs_no3-0.75_rk*s0_no3-s2o3_no3
      _SET_ODE_(self%id_NO3,d_NO3)
   d_NH4 = 0.75_rk*s0_no3+s2o3_no3+DcPOML_so4+DcDOML_so4+DcPOML_s2o3+DcDOML_s2o3
      _SET_ODE_(self%id_NH4,d_NH4)
   d_DIC =(DcDOMR_so4+DcDOMR_s2o3+DcPOMR_SO4+DcPOMR_S2O3)*self%r_c_n
      _SET_ODE_(self%id_DIC,d_DIC)
   d_PO4 = (DcPOML_so4+DcDOML_so4+DcPOML_s2o3+DcDOML_s2o3)/self%r_n_p
      _SET_ODE_(self%id_PO4,d_PO4)
   d_Alk = (&
             !(CH2O)106(NH3)16H3PO4 + 53SO42- =
             !106HCO3- + 16NH3 + H3PO4 + 53H2S (Boudreau, 1996)
!             +2._rk*(DcPM_SO4 +DcDM_SO4) &
             -0.5_rk*s0_disp & !4S0 + 3H2O -> 2H2S + S2O3-- + 2H+
             -1._rk*(-s0_ox) & !2S0 + O2 + H2O -> S2O3-- + 2H+
             -0.5_rk*s0_no3 &  !4S0 + 3NO3- + 7H2O -> 4SO4-- + 3NH4+ + 2H+
             -1._rk*s2o3_ox &  !S2O3-- + 2O2 + 2OH- -> 2SO4-- + H2O
             -0.4_rk*hs_no3 &  !5H2S + 8NO3- + 2OH+ -> 5SO4-- + 4N2 + 6H2O
             )
      _SET_ODE_(self%id_Alk,d_Alk)

      _SET_DIAGNOSTIC_(self%id_s2o3_no3,s2o3_no3)
      _SET_DIAGNOSTIC_(self%id_s0_no3,s0_no3)
      _SET_DIAGNOSTIC_(self%id_s2o3_ox,s2o3_ox)
      _SET_DIAGNOSTIC_(self%id_s0_ox,s0_ox)
      _SET_DIAGNOSTIC_(self%id_DcPOML_s2o3,DcPOML_s2o3)
      _SET_DIAGNOSTIC_(self%id_DcDOML_s2o3,DcDOML_s2o3)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_s2o3,DcPOMR_s2o3)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_s2o3,DcDOMR_s2o3)
      _SET_DIAGNOSTIC_(self%id_DcPOML_so4,DcPOML_so4)
      _SET_DIAGNOSTIC_(self%id_DcDOML_so4,DcDOML_so4)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_so4,DcPOMR_so4)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_so4,DcDOMR_so4)
      _SET_DIAGNOSTIC_(self%id_DcTOM_SOX,DcTOM_SOX)
      _SET_DIAGNOSTIC_(self%id_so4_rd,so4_rd)
      _SET_DIAGNOSTIC_(self%id_s2o3_rd,s2o3_rd)
      _SET_DIAGNOSTIC_(self%id_s0_disp,s0_disp)
      _SET_DIAGNOSTIC_(self%id_hs_ox,hs_ox)
      _SET_DIAGNOSTIC_(self%id_hs_no3,hs_no3)
    _LOOP_END_
  end subroutine do
  
    ! Set increased manganese sinking via MnIV and MnIII oxides formation
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
     class (type_niva_brom_sulfur), intent(in) :: self
     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
     
     real(rk) :: Wadd, WS0_tot
          
     _LOOP_BEGIN_
  
      _GET_(self%id_Wadd,Wadd)
     
      WS0_tot = self%WS0 + Wadd 
  
      _ADD_VERTICAL_VELOCITY_(self%id_S0, WS0_tot)
  
     _LOOP_END_
  
  end subroutine get_vertical_movement

end module fabm_niva_brom_sulfur
