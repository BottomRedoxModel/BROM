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

module fabm_niva_brom_manganese
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_manganese
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_Mn2,id_Mn3,id_Mn4,id_MnS,id_MnCO3,id_PO4_Mn3
    !state variables dependnecies
    type(type_state_variable_id):: id_H2S,id_PO4,id_NH4
    type(type_state_variable_id):: id_POML,id_POMR,id_DOML,id_DOMR
    type(type_state_variable_id):: id_DIC,id_Alk,id_S0,id_O2

    type(type_diagnostic_variable_id):: id_DcDOML_Mn4,id_DcPOML_Mn4,id_DcPOMR_Mn4,id_DcDOMR_Mn4
    type(type_diagnostic_variable_id):: id_DcDOML_Mn3,id_DcPOML_Mn3,id_DcPOMR_Mn3,id_DcDOMR_Mn3
    type(type_diagnostic_variable_id):: id_mn_ox2,id_mn_ox1,id_mn_rd1,id_DcTOM_MnX
    type(type_diagnostic_variable_id):: id_mn_rd2,id_mns_diss,id_mns_form,id_mns_ox
    type(type_diagnostic_variable_id):: id_mnco3_diss,id_mnco3_form
    type(type_diagnostic_variable_id):: id_Wadd

    !diagnostic dependencies
    type(type_dependency_id):: id_Hplus
    type(type_dependency_id):: id_CO3

    !Model parameters
    !sinking
    real(rk):: WMn,  Mn4WMn, WMn_tot
    !specific rates of biogeochemical processes
    !----Mn---------!
    real(rk):: K_mn_ox1,K_mn_ox2,K_mn_rd1,K_mn_rd2,K_mns
    real(rk):: K_mns_diss, K_mns_form, K_mns_ox
    real(rk):: K_mnco3,K_mnco3_diss,K_mnco3_form,K_mnco3_ox
    real(rk):: K_DOML_mn4,K_POML_mn4,K_POMR_mn4,K_DOMR_mn4
    real(rk):: K_DOML_mn3,K_POML_mn3,K_POMR_mn3,K_DOMR_mn3
    real(rk):: s_mnox_mn2,s_mnox_mn3,s_mnrd_mn4,s_mnrd_mn3,s_OM_refr
    !---- S---------!
    real(rk):: K_mnrd_hs
    !---- O--------!
    real(rk):: O2s_dn
    real(rk):: K_mnox_o2
    !---- Stoichiometric coefficients ----!
    real(rk):: r_n_p,r_c_n,r_mn_n
    !---- Partitioning coefficients ----!
    real(rk):: r_mn3_p
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
    class (type_niva_brom_manganese), intent(inout), target :: self
    integer,                          intent(in)            :: configunit
! !LOCAL VARIABLES:
     real(rk) :: EPS

    !-----Model parameters------
    !Sinking
    call self%get_parameter(&
         self%WMn,'WMn','[1/day]','Mn particles sinking rate',default=16.0_rk)
    call self%get_parameter(&
         self%Mn4WMn,'Mn4WMn','[1/day]','Mn4 threshold for sinking rate increase',default=16.0_rk)   
    call self%get_parameter(&
         self%WMn_tot,'WMn_tot','[1/day]','Total accelerated sinking of Mn hydroxides',default=16.0_rk)

    !Specific rates of biogeochemical processes
    !---- Mn---------!
    call self%get_parameter(&
         self%K_mn_ox1,'K_mn_ox1','[1/day]',&
         'Specific rate of oxidation of Mn2 to Mn3 with O2',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_mn_ox2,'K_mn_ox2','[1/day]',&
         'Specific rate of oxidation of Mn3 to Mn4 with O2',&
         default=0.1_rk)
    call self%get_parameter(&
         self%K_mn_rd1,'K_mn_rd1','[1/day]',&
         'Specific rate of reduction of Mn4 to Mn3 with H2S',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_mn_rd2,'K_mn_rd2','[1/day]',&
         'Specific rate of reduction of Mn3 to Mn2 with H2S',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_mns,'K_mns','[M]',&
         'Conditional equilibrium constant for MnS from Mn2 with H2S',&
         default=0.02_rk)
    call self%get_parameter(&
         self%K_mns_diss,'K_mns_diss','[1/day]',&
         'Specific rate of dissolution of MnS to Mn2 and H2S',&
         default=0.0001_rk)
    call self%get_parameter(&
         self%K_mns_form,'K_mns_form','[1/day]',&
         'Specific rate of formation of MnS from Mn2 with H2S',&
         default=0.00001_rk)
    call self%get_parameter(&
         self%K_mns_ox,'K_mns_ox','[1/day]',&
         'Specific rate of  MnS oxidation',&
         default=0.00001_rk)
    call self%get_parameter(&
         self%K_mnco3,'K_mnco3','[M]',&
         'Conditional equilibrium constant 1.8e-11 ',&
         default=1000.0_rk)
    call self%get_parameter(&
         self%K_mnco3_diss,'K_mnco3_diss','[1/day]',&
         'Specific rate of dissolution of MnCO3',&
         default=2.7e-7_rk)
    call self%get_parameter(&
         self%K_mnco3_form,'K_mnco3_form','[1/day]',&
         'Specific rate of formation of MnCO3',&
         default=2.7e-7_rk)
    call self%get_parameter(&
         self%K_mnco3_ox,'K_mnco3_ox','[1/day]',&
         'Specific rate of oxidation of MnCO3 with O2',&
         default=0.0027_rk)
    call self%get_parameter(&
         self%K_DOML_mn4,'K_DOML_mn4','[1/day]',&
         'Specific rate of oxidation of DOML with Mn4',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_DOMR_mn4,'K_DOMR_mn4','[1/day]',&
         'Specific rate of oxidation of DOMR with Mn4',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_POML_mn4,'K_POML_mn4','[1/day]',&
         'Specific rate of oxidation of POML with Mn4',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_POMR_mn4,'K_POMR_mn4','[1/day]',&
         'Specific rate of oxidation of POMR with Mn4',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_DOML_mn3,'K_DOML_mn3','[1/day]',&
         'Specific rate of oxidation of DOML with Mn3',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_DOMR_mn3,'K_DOMR_mn3','[1/day]',&
         'Specific rate of oxidation of DOMR with Mn3',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_POML_mn3,'K_POML_mn3','[1/day]',&
         'Specific rate of oxidation of POML with Mn3',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_POMR_mn3,'K_POMR_mn3','[1/day]',&
         'Specific rate of oxidation of POMR with Mn3',&
         default=0.001_rk)
    call self%get_parameter(&
         self%s_mnox_mn2,'s_mnox_mn2','[uM Mn]',&
         'threshold of Mn2 oxidation',&
         default=0.01_rk)
    call self%get_parameter(&
         self%s_mnox_mn3,'s_mnox_mn3','[uM Mn]',&
         'threshold of Mn3 oxidation',&
         default=0.01_rk)
    call self%get_parameter(&
         self%s_mnrd_mn4,'s_mnrd_mn4','[uM Mn]',&
         'threshold of Mn4 reduciton',&
         default=0.01_rk)
    call self%get_parameter(&
         self%s_mnrd_mn3,'s_mnrd_mn3','[uM Mn]',&
         'threshold of Mn3 reduciton',&
         default=0.01_rk)
    call self%get_parameter(&
         self%s_OM_refr, 's_OM_refr', '[uM N]',&
         'threshold of decay of refractory OM',&
         default=5.0_rk)
    !---- S--------!
    call self%get_parameter(&
         self%K_mnrd_hs,'K_mnrd_hs','[uM S]',&
         'half sat. of Mn reduction',&
         default=1.0_rk)
    !----O2--------!
    call self%get_parameter(&
         self%O2s_dn, 'O2s_dn', '[uM O2]',&
         'half saturation for denitrification',&
         default=10.0_rk)
    call self%get_parameter(&
         self%K_mnox_o2, 'K_mnox_o2', '[uM O2]',&
         'half sat. of Mn oxidation',&
         default=2.0_rk)
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
         self%r_mn_n,'r_mn_n','[-]',&
         'Mn[uM]/N[uM]',&
         default=13.25_rk)
    !---- Partitioning coefficients ----!
    call self%get_parameter(&
         self%r_mn3_p,'r_mn3_p','[-]',&
         'Mn[uM]/P[uM] partitioning coeff. for Mn(III)',&
         default=0.67_rk)
    ! for light
    call self%get_parameter(EPS, 'EPS', 'm^2/mg C', 'specific shortwave attenuation', default=2.208E-3_rk)     
    !Register state variables
    call self%register_state_variable(&
         self%id_Mn2, 'Mn2', 'mmol/m**3','Mn(II)',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_Mn3, 'Mn3', 'mmol/m**3','Mn(III)',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_Mn4, 'Mn4', 'mmol/m**3','Mn(IV)',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_MnS, 'MnS', 'mmol/m**3','MnS',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_MnCO3, 'MnCO3', 'mmol/m**3','MnCO3',&
         minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_PO4_Mn3, 'PO4_Mn3', 'mmol/m**3','PO4_Mn3, PO4 complexed with Mn(III)',&
         minimum=0.0_rk)

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_S0,'S0','mmol/m**3','S0')
    call self%register_state_dependency(&
         self%id_H2S,'H2S','mmol/m**3','H2S')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3',&
         'total dissolved inorganic carbon')
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3',&
         'PO4')
    call self%register_state_dependency(&
         self%id_nh4,'NH4','mmol/m**3',&
         'NH4')
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_POML,'POML','mmol/m**3',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_POMR,'POMR','mmol/m**3',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_DOML,'DOML','mmol/m**3',&
         'dissolved organic nitrogen')
    call self%register_state_dependency(&
         self%id_DOMR,'DOMR','mmol/m**3',&
         'particulate organic nitrogen')
    
    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPOML_Mn4,'DcPOML_Mn4','mmol/m**3/d',&
         'POML with Mn(IV) mineralization ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_Mn4,'DcDOML_Mn4','mmol/m**3/d',&
         'DOML with Mn(IV) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_Mn4,'DcPOMR_Mn4','mmol/m**3/d',&
         'POMR with Mn(IV) mineralization ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_Mn4,'DcDOMR_Mn4','mmol/m**3/d',&
         'DOMR with Mn(IV) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOML_Mn3,'DcPOML_Mn3','mmol/m**3/d',&
         'POML with Mn(III) mineralization ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_Mn3,'DcDOML_Mn3','mmol/m**3/d',&
         'DOML with Mn(III) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_Mn3,'DcPOMR_Mn3','mmol/m**3/d',&
         'POMR with Mn(III) mineralization ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_Mn3,'DcDOMR_Mn3','mmol/m**3/d',&
         'DOMR with Mn(III) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcTOM_MnX,'DcTOM_MnX','mmol/m**3/d',&
         'Total OM with Mn(III) and Mn(IV) mineralization',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mn_ox1,'mn_ox1','mmol/m**3',&
         'Mn(II) with O2 oxidation ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mn_ox2,'mn_ox2','mmol/m**3',&
         'Mn(III) with O2 oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mns_diss,'mns_diss','mmol/m**3',&
         'MnS dissolution',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mns_form,'mns_form','mmol/m**3',&
         'MnS formation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mns_ox,'mns_ox','mmol/m**3/d',&
         'MnS oxidation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mnco3_diss,'mnco3_diss','mmol/m**3',&
         'MnCO3 dissolusion',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mnco3_form,'mnco3_form','mmol/m**3',&
         'MnCO3 formation ',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mn_rd1,'mn_rd1','mmol/m**3',&
         'Mn(IV) with H2S reduction',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_mn_rd2,'mn_rd2','mmol/m**3',&
         'Mn(III) with H2S reduction',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Wadd,'Wadd','[1/day]',&
         'Additional sinking velocity via Mn4 adsorptoin',&
         source=source_get_vertical_movement)

    !Register diagnostic dependencies
    call self%register_dependency(self%id_Hplus,&
         'Hplus', 'mmol/m**3','H+ Hydrogen')
    call self%register_dependency(self%id_CO3,&
      standard_variables%&
      mole_concentration_of_carbonate_expressed_as_carbon)

   ! Register contribution to light extinction
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
        self%id_Mn4,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &                                       
        self%id_MnCO3,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,&
        self%id_MnS,scale_factor=EPS,include_background=.true.)    
    
    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_manganese),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: Mn2,Mn3,Mn4,MnS,MnCO3,PO4_Mn3
    !dependencies
    !state dependencies
    real(rk):: O2,POML,POMR,DOML,H2S,PO4,NH4,DOMR
    !diagnostic variables dependencies
    real(rk):: Hplus,CO3
    !increments
    real(rk):: d_Mn2,d_Mn3,d_Mn4,d_MnS,d_MnCO3,d_PO4_Mn3
    real(rk):: d_S0,d_H2S,d_O2,d_DOML,d_POML,d_POMR,d_DOMR
    real(rk):: d_DIC,d_PO4,d_Alk,d_NH4
    !processes
    real(rk):: mn_ox1,mn_ox2,mn_rd1,mn_rd2,Om_MnS,mns_form,mns_diss, mns_ox
    real(rk):: Om_MnCO3,mnco3_form,mnco3_diss,mnco3_ox
    real(rk):: DcDOML_Mn3, DcPOML_Mn3, DcDOMR_Mn3, DcPOMR_Mn3 
    real(rk):: DcPOML_Mn4, DcDOML_Mn4, DcPOMR_Mn4, DcDOMR_Mn4
    real(rk):: DcTOM_MnX
    !sinking
    real(rk):: WMn, Mn4WMn

    _LOOP_BEGIN_
      !Retrieve current (local) variable values.
      _GET_(self%id_Hplus,Hplus)
      _GET_(self%id_CO3,CO3)
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_PO4,PO4)
      _GET_(self%id_NH4,NH4)
      _GET_(self%id_Mn2,Mn2)
      _GET_(self%id_Mn3,Mn3)
      _GET_(self%id_PO4_Mn3,PO4_Mn3)
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)
      _GET_(self%id_DOMR,DOMR)
      _GET_(self%id_Mn4,Mn4)
      _GET_(self%id_MnS,MnS)
      _GET_(self%id_MnCO3,MnCO3)
      _GET_(self%id_O2,O2)
      _GET_(self%id_H2S,H2S)

      !Mn2 oxidation: 4Mn(2+)+O2+4H(+)->4Mn(3+)+2H2O (Canfield,2005)
      mn_ox1 = max(0._rk, 0.5_rk*(1._rk+tanh(Mn2-self%s_mnox_mn2))*&
               self%K_mn_ox1*Mn2*o2*o2/(o2+self%K_mnox_o2) )
      !Mn3 oxidation: 2Mn3+ + 0.5O2 + 3H20 -> 2MnO2 + 6H+ (Tebo,1997)
      mn_ox2 = max(0._rk, 0.5_rk*(1._rk+tanh(Mn3-self%s_mnox_mn3))*&
               self%K_mn_ox2*Mn3*o2/(o2+self%K_mnox_o2) )
      !Mn4 reduction: 2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0
      mn_rd1 = max(0._rk, 0.5_rk*(1._rk+tanh(Mn4-self%s_mnrd_mn4))*&
               self%K_mn_rd1*Mn4*h2s/(h2s+self%K_mnrd_hs) )
      !Mn3 reduction: 2Mn3+ + HS- -> 2Mn2+ + S0 + H+
      mn_rd2 = max(0._rk, 0.5_rk*(1._rk+tanh(Mn3-self%s_mnrd_mn3))*&
               self%K_mn_rd2*Mn3*h2s/(h2s+self%K_mnrd_hs) )
      !MnS formation/dissollution (dSED) Mn2+ + HS- = MnS(s) + H+
      Om_MnS = H2S*Mn2/(self%K_MnS*Hplus*1000000._rk)
      !Mn2+ + HS- -> MnS(s) + H+
      mns_form = self%K_mns_form*max(0._rk,(Om_MnS-1._rk))
      !MnS(s) + H+ -> Mn2+ + HS-
      mns_diss = self%K_MnS_diss*MnS*max(0._rk,(1._rk-Om_MnS))
      !MnS + 2O2 -> Mn(2+) + SO4(2-)
      mns_ox = self%K_mns_ox*MnS*O2
      !MnCO3 precipitation/dissolution
      Om_MnCO3 = Mn2*CO3/(self%K_MnCO3)
      !Mn2+ + CO3-- <-> MnCO3 (vanCappelen,96)
      mnco3_form = self%K_mnco3_form*max(0._rk,(Om_MnCO3-1._rk))
      mnco3_diss = self%K_mnco3_diss*MnCO3*max(0._rk,(1._rk-Om_MnCO3))
      !2MnCO3(s) + O2 + 2H2O = 2MnO2(s) + 2HCO3- + 2H+ (Morgan,05)
      mnco3_ox = self%K_mnco3_ox*MnCO3*O2
      !(CH2O)106(NH3)16(H3PO4) + 212MnO2 + 318CO2 + 106H2O ->
      !424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996) !in N units
      DcDOML_Mn3 = max(0._rk,self%K_DOML_mn3*DOML &    ! Mn(III)
               *Mn3/(Mn3+0.5_rk) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      DcPOML_Mn3 = max(0._rk,self%K_POML_mn3*POML &
               *Mn3/(Mn3+0.5_rk) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      DcDOMR_Mn3 = max(0._rk,self%K_DOMR_mn3*DOMR &
               *Mn3/(Mn3+0.5_rk) &
               *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      DcPOMR_Mn3 = max(0._rk,self%K_POMR_mn3*POMR &
               *Mn3/(Mn3+0.5_rk) &
               *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      
      DcDOML_Mn4 = max(0._rk, self%K_DOML_mn4*DOML &
               *Mn4/(Mn4+0.5_rk) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      DcPOML_Mn4 = max(0._rk, self%K_POML_mn4*POML &
               *Mn4/(Mn4+0.5_rk) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      DcPOMR_Mn4 = max(0._rk,self%K_POMR_mn4*POMR &
               *Mn4/(Mn4+0.5_rk) &
               *(0.5_rk*(1._rk+tanh((POMR-self%s_OM_refr)*0.1_rk))) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      DcDOMR_Mn4 = max(0._rk,self%K_DOMR_mn4*DOMR &
               *Mn4/(Mn4+0.5_rk) &
               *(0.5_rk*(1._rk+tanh((DOMR-self%s_OM_refr)*0.1_rk))) &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))))
      !Summariazed OM mineralization
      DcTOM_MnX = DcDOMR_Mn3+DcPOMR_Mn3+DcDOMR_Mn4+DcPOMR_Mn4
      !Set increments
      d_Mn2 = -mn_ox1+mn_rd2-mns_form+mns_diss-mnco3_form+mns_ox+&
                mnco3_diss+(DcDOMR_Mn3+DcPOMR_Mn3)*self%r_mn_n
      _SET_ODE_(self%id_Mn2,d_Mn2)
      d_Mn3 = mn_ox1-mn_ox2+mn_rd1-mn_rd2-&
              (DcDOMR_Mn4+DcPOMR_Mn4-DcDOMR_Mn3-DcPOMR_Mn3)*self%r_mn_n
      _SET_ODE_(self%id_Mn3,d_Mn3)
            !complexation of PO4 with Mn(III)
      d_PO4_Mn3 = - PO4_Mn3 &  ! old (last time step) PO4_Mn3
            + Mn3/self%r_mn3_p ! new PO4_Mn3 calculated for old Mn3
      if (d_PO4_Mn3.ge.PO4) d_PO4_Mn3 = PO4
      _SET_ODE_(self%id_PO4_Mn3,d_PO4_Mn3)
      d_Mn4 = mn_ox2-mn_rd1+mnco3_ox &
              +(-DcDOMR_Mn4-DcPOMR_Mn4)*self%r_mn_n
      _SET_ODE_(self%id_Mn4,d_Mn4)
      d_MnS = mns_form-mns_diss-mns_ox
      _SET_ODE_(self%id_MnS,d_MnS)
      d_MnCO3 = mnco3_form-mnco3_diss-mnco3_ox
      _SET_ODE_(self%id_MnCO3,d_MnCO3)
      d_S0 = 0.5_rk*mn_rd1+0.5_rk*mn_rd2
      _SET_ODE_(self%id_S0,d_S0)
      d_H2S = -0.5_rk*mn_rd1-0.5_rk*mn_rd2-mns_form+mns_diss
      _SET_ODE_(self%id_H2S,d_H2S)
      d_O2 = -0.25_rk*mn_ox1-0.25_rk*mn_ox2-0.5_rk*mnco3_ox-2._rk*mns_ox
      _SET_ODE_(self%id_O2,d_O2)
      d_DOML = -DcDOML_Mn4-DcDOML_Mn3
      _SET_ODE_(self%id_DOML,d_DOML)
      d_DOMR = DcDOML_Mn4+DcDOML_Mn3-DcDOMR_Mn4-DcDOMR_Mn3
      _SET_ODE_(self%id_DOMR,d_DOMR)
      d_POML = -DcPOML_Mn4-DcPOML_Mn3
      _SET_ODE_(self%id_POML,d_POML)
      d_POMR = DcPOML_Mn4+DcPOML_Mn3-DcPOMR_Mn4-DcPOMR_Mn3
      _SET_ODE_(self%id_POMR,d_POMR)
      d_DIC = (DcDOMR_Mn3+DcDOMR_Mn4+DcPOMR_Mn3+DcPOMR_Mn4)*self%r_c_n &
               -mnco3_form+mnco3_diss+mnco3_ox
      _SET_ODE_(self%id_DIC,d_DIC)
      d_PO4 = ((DcDOML_Mn4+DcPOML_Mn4+DcDOML_Mn3+DcPOML_Mn3)/self%r_n_p)-d_PO4_Mn3
      _SET_ODE_(self%id_PO4,d_PO4)
      d_NH4 = DcDOML_Mn4+DcPOML_Mn4+DcDOML_Mn3+DcPOML_Mn3
      _SET_ODE_(self%id_NH4,d_NH4)
      d_Alk = (&      !Alkalinity changes due to redox reactions:
             +1._rk*mn_ox1 &   !4Mn2+ + O2 + 4H+ -> 4Mn3+ + 2H2O
             -3._rk*mn_ox2 &   !2Mn3+ + 3H2O  + 0.5 O2 -> 2MnO2 + 6H+
             +3._rk*mn_rd1 &   !2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0
             -1._rk*mn_rd2 &   !2Mn3+ + HS- -> 2Mn2+ + S0 + H+
             !these 4 above, do we need it?
             -2._rk*mns_form & !Mn2+ + H2S <-> MnS + 2H+
             +2._rk*mns_diss &
             -2._rk*mnco3_form &!Mn2+ + CO3-- <-> MnCO3
             +2._rk*mnco3_diss &
             !(CH2O)106(NH3)16(H3PO4) + 212MnO2 + 318CO2 + 106H2O ->
             ! 424HCO3- + 212Mn2+ + 16NH3 + H3PO4
             ! +26.5_rk*(Dc_OM_Mn_total) & 
             ! DcDM_Mn is in N-units,i.e.424/16
             )
      _SET_ODE_(self%id_Alk,d_Alk)

      _SET_DIAGNOSTIC_(self%id_DcPOML_Mn3,DcPOML_Mn3)
      _SET_DIAGNOSTIC_(self%id_DcDOML_Mn3,DcDOML_Mn3)
      _SET_DIAGNOSTIC_(self%id_DcPOML_Mn4,DcPOML_Mn4)
      _SET_DIAGNOSTIC_(self%id_DcDOML_Mn4,DcDOML_Mn4)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_Mn3,DcPOMR_Mn3)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_Mn3,DcDOMR_Mn3)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_Mn4,DcPOMR_Mn4)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_Mn4,DcDOMR_Mn4)
      _SET_DIAGNOSTIC_(self%id_DcTOM_MnX,DcTOM_MnX)
      _SET_DIAGNOSTIC_(self%id_mn_ox1,mn_ox1)
      _SET_DIAGNOSTIC_(self%id_mn_ox2,mn_ox2)
      _SET_DIAGNOSTIC_(self%id_mn_rd1,mn_rd1)
      _SET_DIAGNOSTIC_(self%id_mn_rd2,mn_rd2)
      _SET_DIAGNOSTIC_(self%id_mns_form,mns_form)
      _SET_DIAGNOSTIC_(self%id_mns_diss,mns_diss)
      _SET_DIAGNOSTIC_(self%id_mns_ox,mns_ox)
      _SET_DIAGNOSTIC_(self%id_mnco3_diss,mnco3_diss)
      _SET_DIAGNOSTIC_(self%id_mnco3_form,mnco3_form)
    _LOOP_END_
  end subroutine do
  !
  !
  ! Set increased manganese sinking via MnIV and MnIII oxides formation
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
     class (type_niva_brom_manganese), intent(in) :: self
     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
     
     real(rk) :: Mn4, Wadd, WMn_tot
     
     _LOOP_BEGIN_
     
     _GET_(self%id_Mn4,Mn4)
     
     ! Calculate increased manganese sinking via MnIV and MnIII oxides formation

      Wadd = self%WMn*Mn4/(Mn4+self%Mn4WMn)
      
      WMn_tot = self%WMn + Wadd 

      _SET_DIAGNOSTIC_(self%id_Wadd,Wadd)
      
      _ADD_VERTICAL_VELOCITY_(self%id_Mn4, WMn_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_MnS, WMn_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_MnCO3, WMn_tot)
     _LOOP_END_

  end subroutine get_vertical_movement
end module fabm_niva_brom_manganese
