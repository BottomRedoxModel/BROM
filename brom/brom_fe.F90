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

module fabm_niva_brom_fe
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_fe
    !all descriptions are in initialize subroutine
    type(type_state_variable_id):: id_Fe2,id_Fe3,id_FeS
    type(type_state_variable_id):: id_FeS2,id_FeCO3,id_Fe3PO42
    type(type_state_variable_id):: id_PO4_Fe3
    !state variables dependencies
    type(type_state_variable_id):: id_Mn2,id_Mn4,id_Mn3
    type(type_state_variable_id):: id_H2S
    type(type_state_variable_id):: id_S0,id_SO4
    type(type_state_variable_id):: id_O2
    type(type_state_variable_id):: id_POML,id_POMR,id_DOML,id_PO4,id_Si,id_NH4
    type(type_state_variable_id):: id_DIC,id_Alk,id_DOMR

    type(type_diagnostic_variable_id):: id_DcDOML_Fe,id_DcPOML_Fe,id_DcPOMR_Fe,id_DcDOMR_Fe
    type(type_diagnostic_variable_id):: id_DcTOM_Fe, id_fe_ox1
    type(type_diagnostic_variable_id):: id_fe_ox2,id_fe_ox3,id_feco3_diss
    type(type_diagnostic_variable_id):: id_feco3_form, id_feco3_ox
    type(type_diagnostic_variable_id):: id_fe_rd
    type(type_diagnostic_variable_id):: id_fe_p_compl,id_Kad_PO4
    !type(type_diagnostic_variable_id):: id_fe_p_diss
    !type(type_diagnostic_variable_id):: id_fe_si_compl
    type(type_diagnostic_variable_id):: id_fe3po42_diss,id_fe3po42_form,id_fe3po42_hs
    type(type_diagnostic_variable_id):: id_fes_form,id_fes_diss,id_fes_ox,id_feS2_form
    !diagnostic dependencies
    type(type_dependency_id):: id_Hplus
    type(type_dependency_id):: id_CO3
    type(type_dependency_id):: id_Wadd

    !Model parameters
    !sinking
    real(rk):: WFe , WFe_tot
    !specific rates of biogeochemical processes
    !---- Fe--------!
    real(rk):: K_fe_ox1,K_fe_ox2,K_fe_rd,K_fes,K_fes_form
    real(rk):: K_fes_diss,K_fes_ox,K_DOML_fe,K_POML_fe,K_POMR_fe,K_DOMR_fe
    real(rk):: K_fes2_form,K_fes2_ox,s_feox_fe2,s_ferd_fe3,K_feco3
    real(rk):: K_feco3_diss,K_feco3_form,K_feco3_ox, Sad_Fe3, K_PO4_Fe3
    real(rk):: K_fe3po42, K_fe3po42_diss, K_fe3po42_form, K_fe3po42_ox
    !---- S---------!
    real(rk):: K_ferd_hs
    !---- N--------!
    real(rk):: K_omno_no3, s_OM_refr
    !---- O--------!
    real(rk):: O2s_dn
    !---- Stoichiometric coefficients ----!
    real(rk):: r_n_p,r_c_n,r_fe_n
    !---- Partitioning coefficients ----!
    real(rk):: r_fe3_p,r_fe3_si
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
    class (type_niva_brom_fe), intent(inout), target :: self
    integer,                     intent(in)            :: configunit

    !-----Model parameters------
    !Sinking
    call self%get_parameter(&
         self%WFe,'WFe','1/day',&
         'Fe particles sinking rate',&
          default=7.0_rk)
    call self%get_parameter(&
         self%WFe_tot,'WFe_tot','[1/day]',&
         'Total accelerated sinking with absorbed Mn hydroxides',&
          default=7.0_rk)

    !Specific rates of biogeochemical processes
    !---- Fe---------!
    call self%get_parameter(&
         self%K_fe_ox1,'K_fe_ox1','[1/day]',&
         'Specific rate of oxidation of Fe2 to Fe3 with O2',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_fe_ox2,'K_fe_ox2','[1/day]',&
         'Specific rate of oxidation of Fe2 to Fe3 with MnO2',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_fe_rd,'K_fe_rd','[1/day]',&
         'Specific rate of reduction of Fe3 to Fe2 with H2S',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_fes,'K_fes','[uM]',&
         'FeS equilibrium constant (Solubility Product Constant)',&
         default=2510.0_rk)
    call self%get_parameter(&
         self%K_fes_form,'K_fes_form','[1/day]',&
         'Specific rate of precipitation of FeS from Fe2 with H2S',&
         default=5.e-5_rk)
    call self%get_parameter(&
         self%K_fes_diss,'K_fes_diss','[1/day]',&
         'Specific rate of dissollution of FeS to Fe2 and H2S',&
         default=1.e-6_rk)
    call self%get_parameter(&
         self%K_fes_ox,'K_fes_ox','[1/day]',&
         'Specific rate of oxidation of FeS with O2',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_DOML_fe,'K_DOML_fe','[1/day]',&
         'Specific rate of oxidation of DOML with Fe3',&
         default=0.00005_rk)
    call self%get_parameter(&
         self%K_POML_fe,'K_POML_fe','[1/day]',&
         'Specific rate of oxidation of POML with Fe3',&
         default=0.00001_rk)
    call self%get_parameter(&
         self%K_POMR_fe,'K_POMR_fe','[1/day]',&
         'Specific rate of oxidation of POMR with Fe3',&
         default=0.00001_rk)
    call self%get_parameter(&
         self%K_DOMR_fe,'K_DOMR_fe','[1/day]',&
         'Specific rate of oxidation of DOMR with Fe3',&
         default=0.00001_rk)
    call self%get_parameter(&
         self%K_fes2_form,'K_fes2_form','[1/day]',&
         'Specific rate of FeS2 formation by FeS oxidation by H2S',&
         default=0.000001_rk)
    call self%get_parameter(&
         self%K_fes2_ox,'K_fes2_ox','[1/uM/d]',&
         'Specific rate of pyrite oxidation by O2',&
         default=0.00044_rk)
    call self%get_parameter(&
         self%s_feox_fe2,'s_feox_fe2','[uM Fe]',&
         'threshold of Fe2 reduciton',&
         default=0.001_rk)
    call self%get_parameter(&
         self%s_ferd_fe3,'s_ferd_fe3','[uM Fe]',&
         'threshold of Fe3 reduciton',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_feco3,'K_feco3','[M]',&
         'Conditional equilibrium constant 1.8e-11',&
         default=1000.0_rk)
    call self%get_parameter(&
         self%K_feco3_diss,'K_feco3_diss','[1/day]',&
         'Specific rate of dissolution of FeCO3',&
         default=2.7e-7_rk)
    call self%get_parameter(&
         self%K_feco3_form,'K_feco3_form','[1/day]',&
         'Specific rate of formation of FeCO3',&
         default=2.7e-7_rk)
    call self%get_parameter(&
         self%K_feco3_ox,'K_feco3_ox','[1/day]',&
         'Specific rate of oxidation of FeCO3 with O2',&
         default=0.0027_rk)
    call self%register_diagnostic_variable(&
         self%id_Kad_PO4,'Kad_PO4','-',&
         'Kad_PO4',&
         output=output_time_step_integrated) ! conditional partitioning coef
   call self%get_parameter(self%K_fe3po42,    'K_fe3po42',      '[M]',      'Conditional equilibrium constant %  1.8e-11 ',     default=10000.0_rk)
   call self%get_parameter(self%K_fe3po42_diss, 'K_fe3po42_diss', '[1/day]', 'Specific rate of dissolution of Fe3PO42',   default=2.7e-7_rk)
   call self%get_parameter(self%K_fe3po42_form, 'K_fe3po42_form', '[1/day]', 'Specific rate of formation of Fe3PO42',    default=2.7e-7_rk)
   call self%get_parameter(self%K_fe3po42_ox, 'K_fe3po42_ox',   '[1/day]',  'Specific rate of oxidation of Fe3PO42 with O2',      default=0.0027_rk)
   call self%get_parameter(self%K_ferd_hs,'K_ferd_hs','[uM S]', 'half sat. of Fe reduction',       default=1.0_rk)
   call self%get_parameter(self%K_omno_no3, 'K_omno_no3', '[uM N]','half sat. of no3 for OM denitr.',   default=0.001_rk)
   call self%get_parameter(self%O2s_dn, 'O2s_dn', '[uM O]','half saturation for denitrification',     default=10.0_rk)
   call self%get_parameter(self%K_PO4_Fe3, 'K_PO4_Fe3', '[-]',  'partitioning coeff. for  PO4 on Fe3',   default=100000.0_rk)
   call self%get_parameter(self%Sad_Fe3,   'Sad_Fe3', '[-]',    'adsorbtion sites on Fe3',                default=0.01_rk)
       call self%get_parameter(&
         self%s_OM_refr, 's_OM_refr', '[uM N]',&
         'threshold of decay of refractory OM',&
         default=50.0_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(&
         self%r_n_p,   'r_n_p',  '[-]',&
         'N[uM]/P[uM]',&
         default=16.0_rk)
    call self%get_parameter(&
         self%r_c_n,   'r_c_n',  '[-]',&
         'C[uM]/N[uM]',&
         default=6.625_rk)
    call self%get_parameter(&
         self%r_fe_n,   'r_fe_n',  '[-]',&
         'Fe[uM]/N[uM]',&
         default=26.5_rk)
    call self%get_parameter(&
         self%r_fe3_p,  'r_fe3_p',  '[-]',&
         'Fe[uM]/P[uM] partitioning coeff. for Fe oxide',&
         default=2.7_rk)
    call self%get_parameter(&
         self%r_fe3_si,   'r_fe3_si',  '[-]',&
         'Fe[uM]/Si[uM] partitioning coeff. for Fe oxide',&
         default=2.7_rk)

    !Register state variables
    call self%register_state_variable(&
         self%id_Fe2, 'Fe2', 'mmol/m**3','Fe(II) dissolved',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_Fe3, 'Fe3', 'mmol/m**3','Fe(III) oxides',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_FeS, 'FeS', 'mmol/m**3','FeS',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_FeCO3, 'FeCO3', 'mmol/m**3','FeCO3',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_FeS2, 'FeS2', 'mmol/m**3','FeS2',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_Fe3PO42, 'Fe3PO42', 'mmol/m**3','Fe3PO42',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_PO4_Fe3, 'PO4_Fe3', 'mmol/m**3','PO4_Fe3 adsorbed',&
         minimum=0.0_rk)

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_Mn2, 'Mn2', 'mmol/m**3','Mn(II)')
    call self%register_state_dependency(&
         self%id_Mn4, 'Mn4', 'mmol/m**3','Mn(IV)')
    call self%register_state_dependency(&
         self%id_Mn3, 'Mn3', 'mmol/m**3','Mn(III)')
    call self%register_state_dependency(&
         self%id_H2S, 'H2S', 'mmol/m**3','H2S')
    call self%register_state_dependency(&
         self%id_S0 , 'S0',  'mmol/m**3','S0')
    call self%register_state_dependency(&
         self%id_SO4, 'SO4', 'mmol/m**3','SO4')
    call self%register_state_dependency(&
         self%id_Si, 'Si', 'mmol/m**3','Si')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3',&
         'total DIC')
    call self%register_state_dependency(&
         self%id_Alk,standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_PO4,'PO4','mmol/m**3',&
         'PO4')
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3',&
         'NH4')
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_POML,'POML','mmol/m**3',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_POMR,'POMR','mmol/m**3',&
         'POMR')
    call self%register_state_dependency(&
         self%id_DOMR,'DOMR','mmol/m**3',&
         'DOMR')
    call self%register_state_dependency(&
         self%id_DOML,'DOML','mmol/m**3',&
         'dissolved organic nitrogen')

    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPOML_Fe,'DcPOML_Fe','mmol/m**3',&
         'POML with Fe(III) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_Fe,'DcPOMR_Fe','mmol/m**3',&
         'POMR with Fe(III) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_Fe,'DcDOMR_Fe','mmol/m**3',&
         'DOMR with Fe(III) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_Fe,'DcDOML_Fe','mmol/m**3',&
         'DOML with Fe(III) mineralization  ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcTOM_Fe,'DcTOM_Fe','mmol/m**3',&
         'Total OM with Fe(III) mineralization  ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self% id_fe_ox1,'fe_ox1','mmol/m**3',&
         'Fe(II) with O2 oxidation ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self% id_fe_ox2,'fe_ox2','mmol/m**3',&
         'Fe(II)  with Mn(IV) oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self% id_fe_ox3,'fe_ox3','mmol/m**3',&
         'Fe(II)  with Mn(III) oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fe_rd,'fe_rd','mmol/m**3',&
         'Fe (III) with H2S reduction',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_feco3_diss,'feco3_diss','mmol/m**3',&
         'FeCO3 dissolusion',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_feco3_form,'feco3_form','mmol/m**3',&
         'FeCO3 formation ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_feco3_ox,'feco3_ox','mmol/m**3',&
         'FeCO3 oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fe_p_compl,'fe_p_compl','mmol/m**3',&
         'complexation of P with Fe(III)',&
         output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_fe_p_diss,'fe_p_diss','mmol/m**3',&
    !     'dissolution of complexation of P with Fe(III)',&
    !     output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_fe_si_compl,'fe_si_compl','mmol/m**3',&
    !     'complexation of Si with Fe(III)',&
    !     output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fes_form,'fes_form','mmol/m**3 day',&
         'FeS formation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fes_diss,'fes_diss','mmol/m**3 day',&
         'FeS dissolution',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fes_ox,'fes_ox','mmol/m**3 day',&
         'FeS oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_feS2_form,'feS2_form','mmol/m**3 day',&
         'FeS2 formation',&
         output=output_time_step_integrated)
!oxydation of FeCO3 is missed
    call self%register_diagnostic_variable(&
         self%id_fe3po42_diss,'fe3po42_diss','mmol/m**3',&
        'Fe3PO42 dissolusion',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fe3po42_form,'fe3po42_form','mmol/m**3',&
        'Fe3PO42 formation ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_fe3po42_hs,'fe3po42_hs','mmol/m**3',&
        'Fe3PO42 reaction with H2S ',           &
                output=output_time_step_integrated)
    !Register diagnostic dependencies
    call self%register_dependency(self%id_Hplus,&
         'Hplus', 'mmol/m**3','H+ Hydrogen')
    call self%register_dependency(self%id_CO3,&
      standard_variables%&
      mole_concentration_of_carbonate_expressed_as_carbon)
    call self%register_dependency(self%id_Wadd,'Wadd','[1/day]',&
         'Additional sinking velocity via Mn4 adsorptoin')

    !Specify that rates are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_fe),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: Fe2,Fe3,FeS,FeS2,FeCO3,Fe3PO42,PO4_Fe3
    !state dependencies
    real(rk):: O2,POML,POMR,DOML,DOMR
    real(rk):: Mn4,Mn3,H2S,NH4,PO4
    !diagnostic variables dependencies
    real(rk):: Hplus,CO3
    !increments
    real(rk):: d_Mn2,d_Mn4,d_Mn3
    real(rk):: d_Fe2,d_Fe3,d_FeS,d_FeS2,d_FeCO3,d_Fe3PO42,d_PO4_Fe3 
    real(rk):: d_SO4,d_S0,d_H2S,d_O2,d_DOML,d_POML,d_POMR,d_DOMR
    real(rk):: d_DIC,d_Si,d_PO4,d_NH4
    real(rk):: d_Alk
    !processes
    real(rk):: fe_ox1,fe_rd,fe_ox2,fe_ox3
    real(rk):: Om_FeS,fes_form,fes_diss,fes_ox
    real(rk):: fes2_form,fes2_ox,Om_FeCO3,feco3_form
    real(rk):: feco3_diss,feco3_ox
    real(rk):: fe_p_compl, Kad_PO4 !,fe_p_diss,fe_si_compl
    real(rk):: DcDOML_Fe,DcPOML_Fe,DcPOMR_Fe,DcDOMR_Fe,DcTOM_Fe
    real(rk):: Om_Fe3PO42,fe3po42_form
    real(rk):: fe3po42_diss,fe3po42_hs

    _LOOP_BEGIN_
      !Retrieve variable values
      !state
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_Fe2,Fe2)
      _GET_(self%id_PO4,PO4)
      _GET_(self%id_NH4,NH4)
      !solids
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)
      _GET_(self%id_DOMR,DOMR)
      _GET_(self%id_Mn4,Mn4)
      _GET_(self%id_Mn3,Mn3)
      _GET_(self%id_Fe3,Fe3)
      _GET_(self%id_FeS,FeS)
      _GET_(self%id_FeS2,FeS2)
      _GET_(self%id_FeCO3,FeCO3)
      _GET_(self%id_Fe3PO42,Fe3PO42)
      _GET_(self%id_PO4_Fe3,PO4_Fe3)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_H2S,H2S)
      !diagnostic
      _GET_(self%id_CO3,CO3)
      _GET_(self%id_Hplus,Hplus)

      !Fe
      !Fe2 oxidation1: 4Fe2+ + O2 + 10H2O -> 4Fe(OH)3 +8H+ (vanCappelen,96)
      fe_ox1 = 0.5_rk*(1._rk+tanh(Fe2-self%s_feox_fe2))*&
               self%K_fe_ox1*o2*Fe2
      !
      !Fe2 oxidation2: Fe2+ + MnO2 + 4H+ -> Fe3+ + Mn2+ + 2H2O (vanCappelen,96)
      !                2Fe2+ + MnO2 + 4H2O -> 2Fe(OH)3 + Mn2+ + 2H+ (Pakhomova, p.c.)
      fe_ox2 = 0.5_rk*(1._rk+tanh(Fe2-self%s_feox_fe2))*&
               0.5_rk*(1._rk+tanh(Mn4-self%s_feox_fe2))*&
               self%K_fe_ox2*Mn4*Fe2
      !Fe2 oxidation2: Fe2+ + Mn3+ 3H2O->  Fe(OH)3 + Mn2+ + 3H+ (Pakhomova, p.c.)
      fe_ox3 = 0.5_rk*(1._rk+tanh(Fe2-self%s_feox_fe2))*&
               0.5_rk*(1._rk+tanh(Mn3-self%s_feox_fe2))*&
               self%K_fe_ox2*Mn3*Fe2
      !
      !Fe3 reduction: 2Fe(OH)3 + HS- + 5H+ -> 2Fe2+ + S0 + 6H2O
      fe_rd = 0.5_rk*(1._rk+tanh(Fe3-self%s_ferd_fe3))*&
              self%K_fe_rd*Fe3*h2s/(h2s+self%K_ferd_hs)
      !
      !FeS formation/dissollution (Bektursunova,11)
      Om_FeS = H2S*Fe2/(self%K_fes*Hplus*1000000._rk)
      !
      !FeS formation Fe2+ + HS- -> FeS + H+ (Bektursunova,11)
      fes_form = self%K_fes_form*max(0._rk,(Om_FeS-1._rk))
      !
      !FeS dissollution FeS + H+ -> Fe2+ + HS (Bektursunova,11)
      fes_diss = self%K_fes_diss*FeS*max(0._rk,(1._rk-Om_FeS))
      !
      !FeS oxidation: FeS + 2.25O2 +H2O -> 0.5Fe2O3 + 2H+ +SO42-
      !(Soetaert,07) or FeS + 2O2 -> Fe2+ + SO42-(Bektursunova,11)
      fes_ox = self%K_fes_ox*FeS*O2
      !
      !Pyrite formation by FeS oxidation by H2S
      !FeS + H2S -> FeS2 + H2 (Rickard,97)
      fes2_form = self%K_fes2_form*H2S*FeS
      !
      !Pyrite oxidation by O2
      !FeS2 + 3.5 O2 + H2O = Fe2+ + 2SO42- + 2H+ (Wijsman,02)
      fes2_ox = self%K_fes2_ox*FeS2*o2
      !
      !FeCO3 precipitation/dissolution
      Om_FeCO3 = Fe2*CO3/(self%K_FeCO3)
      !
      !Fe2+ + CO3-- <-> FeCO3 (vanCappelen,96)
      feco3_form = self%K_feco3_form*max(0._rk,(Om_FeCO3-1._rk))
      feco3_diss = self%K_FeCO3_diss*FeCO3*max(0._rk,(1._rk-Om_FeCO3))
      !
      !FeCO3(s) + O2 + 2H2O = Fe2O3(s) + HCO3- + H+ (Morgan,05)
      !2FeCO3 + O2 + 4H2O = 2Fe(OH)3 + 2HCO3- + 6H+
      feco3_ox = self%K_feco3_ox*FeCO3*O2
      !
      !  Fe3(PO4)2 precipitation/dissolution
      !   Om_Fe3PO42=Fe2*PO4/(self%K_Fe3PO42*PO4)
        Om_Fe3PO42=Fe2*Fe2*Fe2*PO4*PO4/(self%K_fe3po42)
      ! ! 3Fe2+ + 2PO4--- <-> Fe3(PO4)2 (?):
        fe3po42_form=self%K_fe3po42_form*max(0._rk,(Om_Fe3PO42-1._rk))
        fe3po42_diss=self%K_Fe3PO42_diss*Fe3PO42*max(0._rk,(1._rk-Om_Fe3PO42))
      !%  4Fe3(PO4)2(s)  + 3O2  +   12H2O  =   6Fe2O3(s)  +   8PO4---  +   24H+ (?)=
        !fe3po42_ox=self%K_fe3po42_ox*Fe3PO42*O2
!%  Fe3(PO4)2(s)  + 3H2S = 3FeS(s)  +   2PO4---  +   6H+ (?)=
        fe3po42_hs=0._rk !Fe3PO42*H2S

      !(CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 ->
      ! 848HCO3-+ 424Fe2+ +318H2O +16NH3 +H3PO4 (Boudreau,1996) Fe units
      DcDOML_Fe = self%K_DOML_fe*DOML &
               *Fe3/(Fe3+self%K_omno_no3) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn)))
      DcPOML_Fe = self%K_POML_fe*POML &
               *Fe3/(Fe3+self%K_omno_no3) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn)))
      DcPOMR_Fe = self%K_POMR_fe*POMR&
               *Fe3/(Fe3+self%K_omno_no3) &
               *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn)))
      DcDOMR_Fe = self%K_DOMR_fe*DOMR&
               *Fe3/(Fe3+self%K_omno_no3) &
               *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn)))
      !!!complexation of P with Fe(III)
      !fe_p_compl = ((fe_ox1+fe_ox2+fe_ox3+fes_ox+feco3_ox)*PO4/(PO4+0.1) &
      !        -fe_rd-(DcDOML_Fe+DcPOML_Fe)*self%r_fe_n)/self%r_fe3_p 
!"      fe_p_compl = ((0.0006_rk*(fe_ox1+fe_ox2+fe_ox3+fes_ox+feco3_ox &
!"              -fe_rd-(DcDOMR_Fe+DcPOMR_Fe)*self%r_fe_n))/(1.e-9_rk/Hplus+0.06_rk*PO4))*PO4

 !!!!!! Sorption of PO4 on Fe oxides
     Kad_PO4=self%K_PO4_Fe3*self%Sad_Fe3*Fe3/(1.e-9_rk/Hplus+self%K_PO4_Fe3*PO4)
   fe_p_compl= Kad_PO4*(PO4+PO4_Fe3)/(1.0_rk+Kad_PO4)-PO4_Fe3
 !!!!  _SET_ODE_(self%id_Hg2_Fe3,hg2_fe3_compl)
      
!      fe_p_diss = 0.0_rk !((0.0006_rk*fe_rd)/(1.e-9_rk/Hplus+0.06_rk*PO4))*PO4
      !!!complexation of Si with Fe(III)
      !!fe_si_compl = (fe_rd-fe_ox1-fe_ox2+4._rk*DcDOML_Fe+4._rk*DcPOML_Fe)/&
      !!               self%r_fe3_si
!      fe_si_compl =  0.0_rk

      !Summariazed OM mineralization
      DcTOM_Fe = DcDOMR_Fe+DcPOMR_Fe

      !Set increments
   d_Mn2 = 0.5_rk*fe_ox2+fe_ox3
      _SET_ODE_(self%id_Mn2,d_Mn2)
   d_Mn3 = -fe_ox3
      _SET_ODE_(self%id_Mn3,d_Mn3)
   d_Mn4 = -0.5_rk*fe_ox2
      _SET_ODE_(self%id_Mn4,d_Mn4)
   d_Fe2 = -fe_ox1-fe_ox2-fe_ox3+fe_rd-fes_form+fes_diss &
               -feco3_form+feco3_diss-fe3po42_form+fe3po42_diss &
               +(DcPOMR_Fe+DcDOMR_Fe)*self%r_fe_n+feS2_ox
      _SET_ODE_(self%id_Fe2,d_Fe2)
   d_Fe3 = fe_ox1+fe_ox2+fe_ox3-fe_rd+fes_ox+feco3_ox&
              -(DcPOMR_Fe+DcDOMR_Fe)*self%r_fe_n
      _SET_ODE_(self%id_Fe3,d_Fe3)
   d_FeS = fes_form-fes_diss-fes_ox-feS2_form+fe3po42_hs
      _SET_ODE_(self%id_FeS,d_FeS)
   d_FeS2 = feS2_form-feS2_ox
      _SET_ODE_(self%id_FeS2,d_FeS2)
   d_FeCO3 = feco3_form-feco3_diss-feco3_ox
      _SET_ODE_(self%id_FeCO3,d_FeCO3)
  !d_Fe3PO42
      _SET_ODE_(self%id_Fe3PO42, fe3po42_form-fe3po42_diss-fe3po42_hs)
   d_SO4 = fes_ox+2._rk*feS2_ox
      _SET_ODE_(self%id_SO4,d_SO4)
   d_S0 = 0.5_rk*fe_rd
      _SET_ODE_(self%id_S0,d_S0)
   d_H2S = -0.5_rk*fe_rd-fes_form+fes_diss-feS2_form+fe3po42_hs
      _SET_ODE_(self%id_H2S,d_H2S)
   d_O2 = -0.25_rk*fe_ox1-2.25_rk*fes_ox-3.5_rk*feS2_ox+feco3_ox
      _SET_ODE_(self%id_O2,d_O2)
   d_DOML = -DcDOML_Fe
      _SET_ODE_(self%id_DOML,d_DOML)
   d_POML = -DcPOML_Fe
      _SET_ODE_(self%id_POML,d_POML)
   d_POMR = DcPOML_Fe-DcPOMR_Fe
      _SET_ODE_(self%id_POMR,d_POMR)
   d_DOMR = DcDOML_Fe-DcDOMR_Fe
      _SET_ODE_(self%id_DOMR,d_DOMR)
   d_DIC = (DcPOMR_Fe+DcDOMR_Fe)*self%r_c_n-feco3_form+feco3_diss+feco3_ox
      _SET_ODE_(self%id_DIC,d_DIC)
   d_Si = 0.0_rk !fe_si_compl
      _SET_ODE_(self%id_Si,d_Si)
   d_PO4 = (DcDOML_Fe+DcPOML_Fe)/self%r_n_p-fe_p_compl &
              -0.66_rk*fe3po42_form+0.66_rk*fe3po42_diss+0.66_rk*fe3po42_hs
      _SET_ODE_(self%id_PO4,d_PO4)
  !PO4 complexed with Fe3
      _SET_ODE_(self%id_PO4_Fe3,fe_p_compl)
   d_NH4 = DcDOML_Fe+DcPOML_Fe
      _SET_ODE_(self%id_NH4,d_NH4)
   d_Alk = (&                  !Alkalinity changes due to redox reactions:
             -2._rk*fe_ox1 &   !4Fe2+ + O2 +10H2O-> 4Fe(OH)3 +8H+
             -1._rk*fe_ox2 &   !2Fe2+ + MnO2 +4H2O -> 2Fe(OH)3 + Mn2+ +2H+
             -3._rk*fe_ox3 &   !Fe2+ + Mn3+ 3H2O->  Fe(OH)3 + Mn2+ + 3H+ (Pakhomova, p.c.)
             +2._rk*fe_rd &    !2Fe(OH)3 + HS- + 5H+ -> 2Fe2+ + S0 + 6H2O
             !(here and below d(AlK_H2S) is excluded, as give before)
             -1._rk*fes_form & !Fe2+ + H2S <-> FeS + H+
             +1._rk*fes_diss &
             -2._rk*fes_ox &  !FeS + 2.25O2 +H2O -> 0.5Fe2O3 + 2H+ +SO42-
             -2._rk*fes2_ox & !FeS2 + 3.5O2 + H2O -> Fe2+ + 2SO42- + 2H+
             -2._rk*feco3_form & !Fe2+ + CO3-- <-> FeCO3
             +2._rk*feco3_diss &
             !(CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 ->
             ! 848HCO3-+ 424Fe2+ +318H2O +16NH3 +H3PO4
             ! + 53._rk*(DcDOML_Fe+DcPOML_Fe) & !DcDOML_Fe is in N-units,i.e.848/16
             )
      _SET_ODE_(self%id_Alk,d_Alk)

      _SET_DIAGNOSTIC_(self%id_DcDOMR_Fe,DcDOMR_Fe)
      _SET_DIAGNOSTIC_(self%id_DcDOML_Fe,DcDOML_Fe)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_Fe,DcPOMR_Fe)
      _SET_DIAGNOSTIC_(self%id_DcPOML_Fe,DcPOML_Fe)
      _SET_DIAGNOSTIC_(self%id_DcTOM_Fe,DcTOM_Fe)
      _SET_DIAGNOSTIC_(self%id_fe_ox1,fe_ox1)
      _SET_DIAGNOSTIC_(self%id_fe_rd,fe_rd)
      _SET_DIAGNOSTIC_(self%id_fe_ox2,fe_ox2)
      _SET_DIAGNOSTIC_(self%id_fe_ox3,fe_ox3)
      _SET_DIAGNOSTIC_(self%id_feco3_diss,feco3_diss)
      _SET_DIAGNOSTIC_(self%id_feco3_ox,feco3_ox)
      _SET_DIAGNOSTIC_(self%id_feco3_form,feco3_form)
      _SET_DIAGNOSTIC_(self%id_fe3po42_diss,fe3po42_diss)
      _SET_DIAGNOSTIC_(self%id_fe3po42_form,fe3po42_form)
      _SET_DIAGNOSTIC_(self%id_fe3po42_hs,fe3po42_hs)
      _SET_DIAGNOSTIC_(self%id_fe_p_compl,fe_p_compl)
      _SET_DIAGNOSTIC_(self%id_Kad_PO4,Kad_PO4)
!      _SET_DIAGNOSTIC_(self%id_fe_p_diss,fe_p_diss)
!      _SET_DIAGNOSTIC_(self%id_fe_si_compl,fe_si_compl)
      _SET_DIAGNOSTIC_(self%id_fes_form,fes_form)
      _SET_DIAGNOSTIC_(self%id_fes_diss,fes_diss)
      _SET_DIAGNOSTIC_(self%id_fes_ox,fes_ox)
      _SET_DIAGNOSTIC_(self%id_feS2_form,feS2_form)
    _LOOP_END_
  end subroutine do
  
  ! Set increased manganese sinking via MnIV and MnIII oxides formation
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
     class (type_niva_brom_fe), intent(in) :: self
     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
     
     real(rk) :: Wadd, WFe_tot
          
     _LOOP_BEGIN_
  
      _GET_(self%id_Wadd,Wadd)
     
      WFe_tot = self%WFe + Wadd 
  
      _ADD_VERTICAL_VELOCITY_(self%id_Fe3, WFe_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_FeS, WFe_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_FeCO3, WFe_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_FeS2, WFe_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_Fe3PO42, WFe_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_PO4_Fe3, WFe_tot)
  
     _LOOP_END_
  
  end subroutine get_vertical_movement
end module fabm_niva_brom_fe
