!-----------------------------------------------------------------------
! fabm_niva_brom_bio is
! free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_bio
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_bio
    !variables allocated here
    type(type_state_variable_id):: id_Phy,id_Het
    type(type_state_variable_id):: id_O2,id_POML,id_POMR,id_DOML,id_DOMR
    !state dependencies
    type(type_state_variable_id):: id_NO2,id_NO3,id_NH4,id_PO4
    type(type_state_variable_id):: id_Baae,id_Baan,id_Bhae,id_Bhan
    type(type_state_variable_id):: id_DIC,id_H2S,id_Si,id_Sipart,id_Alk

    type(type_diagnostic_variable_id):: id_DcPOML_O2,id_DcDOML_O2,id_DcPOMR_O2,id_DcDOMR_O2
    type(type_diagnostic_variable_id):: id_MortHet,id_Grazing,id_RespHet, id_DcTOM_O2
    type(type_diagnostic_variable_id):: id_GrazBhae,id_GrazBhan
    type(type_diagnostic_variable_id):: id_GrazBaae,id_GrazBaan
    type(type_diagnostic_variable_id):: id_GrazPhy,id_GrazPOP,id_GrazBact
    type(type_diagnostic_variable_id):: id_MortPhy,id_ExcrPhy,id_LimNH4
    type(type_diagnostic_variable_id):: id_LimN,id_GrowthPhy
    type(type_diagnostic_variable_id):: id_LimT,id_LimP,id_LimNO3,id_LimSi
    type(type_diagnostic_variable_id):: id_LimLight, id_N_fixation
    type(type_diagnostic_variable_id):: id_O2_rel_sat,id_O2_sat, id_POMTot,id_DOMTot, id_AOU

    type(type_dependency_id):: id_temp,id_salt,id_par,id_pres
    type(type_dependency_id):: id_Hplus
    type(type_dependency_id):: id_Wadd
    type (type_horizontal_dependency_id) :: id_windspeed, id_aice
    logical   :: use_aice ! added ice area limitation of air-sea flux, with new switch parameter use_aice.

    !Model parameters
    !specific rates of biogeochemical processes
    real(rk):: K_DOML_ox,K_POML_DOML,K_POML_ox,K_POMR_ox
    real(rk):: K_DOMR_ox,K_POMR_DOMR,Tda,beta_da,K_omox_o2
    !----Phy  ----------!
    real(rk):: K_phy_gro,k_Erlov,Iopt
    real(rk):: K_phy_mrt,K_phy_exc,LatLight
    integer :: phy_t_dependence ! select dependence on T: (1) ERGOM; (2) for Arctic; (3) ERSEM
    !----Het -----------!
    real(rk):: K_het_phy_gro,K_het_phy_lim,K_het_pom_gro,K_het_pom_lim,K_het_bac_gro
    real(rk):: K_het_res,K_het_mrt,Uz,Hz,limGrazBac, s_hmort_H2S
    !---- O2--------!
    !Upper boundary, for oxygen flux calculations
    real(rk):: pvel = 5._rk ! wind speed [m/s]
    real(rk):: a0 = 31.25_rk !oxygen saturation [uM]
    real(rk):: a1 = 14.603_rk !oxygen saturation [-]
    real(rk):: a2 = 0.4025_rk !oxygen saturation [1/degC]
    !---- N, P, Si--!
    real(rk):: K_nox_lim,K_nh4_lim,K_psi,K_nfix,K_po4_lim,K_si_lim
    !---- Sinking---!
    real(rk):: Wsed,Wphy,Whet ,Wsed_tot,Wphy_tot,Whet_tot
    !---- Stoichiometric coefficients ----!
    real(rk):: r_n_p, r_o_n, r_c_n, r_si_n      
    
  contains
    procedure :: initialize
    procedure :: do
    procedure :: do_surface
    procedure :: f_t
    procedure :: graz
    procedure :: get_vertical_movement

    end type
    contains

    subroutine initialize(self,configunit)
        class (type_niva_brom_bio), intent(inout), target :: self
        integer,                    intent(in)            :: configunit
 
 !LOCAL VARIABLES:      
   real(rk) :: EPS        

    call self%get_parameter(&
         self%LatLight,'LatLight','degree','Latitude',default=50.0_rk)
    call self%get_parameter(&
         self%K_DOML_ox,'K_DOML_ox','[1/day]',&
         'Specific rate of oxidation of DOML with O2',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_DOMR_ox,'K_DOMR_ox','[1/day]',&
         'Specific rate of oxidation of DOMR with O2',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_POML_ox,'K_POML_ox','[1/day]',&
         'Specific rate of oxidation of POML with O2',&
         default=0.002_rk)
    call self%get_parameter(&
         self%K_POMR_ox,'K_POMR_ox','[1/day]',&
         'Specific rate of oxidation of POMR with O2',&
         default=0.002_rk)
    call self%get_parameter(&
         self%K_POML_DOML, 'K_POML_DOML', '[1/day]',&
         'Specific rate of Autolysis of POML to DOML',&
         default=0.1_rk)
    call self%get_parameter(&
         self%K_POMR_DOMR, 'K_POMR_DOMR', '[1/day]',&
         'Specific rate of Autolysis of POMR to DOMR',&
         default=0.1_rk)
    call self%get_parameter(&
         self%beta_da,'beta_da','[1/day]',&
         'Temperature control coefficient for OM decay',&
         default=20.0_rk)
    call self%get_parameter(&
         self%K_omox_o2,'K_omox_o2','[uM]',&
         'half sat. of o2 for OM mineralization',&
         default=1.0_rk)
    call self%get_parameter(&
         self%Tda,'Tda','[1/day]',&
         'Temperature control coefficient for OM decay',&
         default=13.0_rk)

   call self%get_parameter(self%use_aice, 'use_aice', '', 'use ice area to limit air-sea flux',  default=.false.)
    
    !----Phy----------!
    call self%get_parameter(&
         self%K_phy_gro,'K_phy_gro','1/d','Maximum specific growth rate',&
         default=2.0_rk)
    call self%get_parameter(&
         self%k_Erlov,'k_Erlov', '1/m','Extinction coefficient',&
         default=0.05_rk)
    call self%get_parameter(&
         self%Iopt,'Iopt','Watts/m**2/h','Optimal irradiance',&
         default=25.0_rk)
    call self%get_parameter(&
         self%K_phy_mrt,'K_phy_mrt','1/d','Specific rate of mortality',&
         default=0.10_rk)
    call self%get_parameter(&
         self%K_phy_exc,'K_phy_exc','1/d','Specific rate of excretion',&
         default=0.01_rk)
    call self%get_parameter(&
         self%phy_t_dependence,'phy_t_dependence','-','T dependence fro Phy growth',&
         default=1)
    
    !----Het----------!
    call self%get_parameter(&
         self%K_het_phy_gro,'K_het_phy_gro','1/d',&
         'Max.spec. rate of grazing of Het on Phy',&
         default=1.0_rk)
    call self%get_parameter(&
         self%K_het_phy_lim,'K_het_phy_lim','nd',&
         'Half-sat.const.for grazing of Het on Phy for Phy/Het ratio',&
         default=1.1_rk)
    call self%get_parameter(&
         self%K_het_pom_gro,'K_het_pom_gro','mmol/m**3',&
         'Max.spec.rate of grazing of Het on POM',&
         default=0.70_rk)
    call self%get_parameter(&
         self%K_het_bac_gro,'K_het_bac_gro','mmol/m**3',&
         'Max.spec.rate of grazing of Het on POM',&
         default=0.70_rk)
    call self%get_parameter(&
         self%K_het_pom_lim,'K_het_pom_lim','nd',&
         'Half-sat.const.for grazing of Het on POM for POM/Het ratio',&
         default=0.2_rk)
    call self%get_parameter(&
         self%K_het_res,'K_het_res','1/d',&
         'Specific respiration rate',&
         default=0.02_rk)
    call self%get_parameter(&
         self%K_het_mrt,'K_het_mrt','1/d',&
         'Maximum specific rate of mortality of Het',&
         default=0.05_rk)
     call self%get_parameter(&
        self%s_hmort_H2S,'s_hmort_H2S','1/d',&
        'threshold of H2S for Het additional mortality',&
        default=0.05_rk)
    call self%get_parameter(&
         self%Uz,'Uz','nd',&
         'Food absorbency for Het',&
         default=0.5_rk)
    call self%get_parameter(&
         self%Hz,'Hz','nd',&
         'Ratio betw. diss. and part. excretes of Het',&
         default=0.5_rk)
    call self%get_parameter(&
         self%limGrazBac,'limGrazBac','mmol/m**3',&
         'Limiting parameter for bacteria grazing by Het',&
         default=2._rk)
    
    !----N---------------
    call self%get_parameter(&
         self%K_psi,'K_psi','[nd]',&
         'Strength of NH4 inhibition of NO3 uptake constant',&
         default=1.46_rk)
    call self%get_parameter(&
         self%K_nox_lim,'K_nox_lim','[mmol/m**3]',&
         'Half-sat.const.for uptake of NO3+NO2',&
         default=0.15_rk)
    call self%get_parameter(&
         self%K_nh4_lim,'K_nh4_lim','[mmol/m**3]',&
         'Half-sat.const.for uptake of NH4',&
         default=0.02_rk)
    call self%get_parameter(&
         self%K_nfix,'K_nfix','[1/d]',&
         'Max. specific rate of mitrogen fixation',&
         default=10._rk)
    
    !----P---------------
    call self%get_parameter(&
         self%K_po4_lim,'K_po4_lim','[mmol/m**3]',&
         'Half-sat. constant for uptake of PO4 by Phy',&
         default=0.02_rk)
    
    !----Si-------------
    call self%get_parameter(&
         self%K_si_lim,'K_si_lim','[mmol/m**3]',&
         'Half-sat. constant for uptake of Si by Phy',&
         default=0.02_rk)
    
    !----Sinking--------
    call self%get_parameter(&
         self%Wsed,'Wsed','[1/day]',&
         'Rate of sinking of detritus (POP, POM)',&
         default=5.00_rk)
    call self%get_parameter(&
         self%Wphy,'Wphy','[m/day]',&
         'Rate of sinking of Phy',&
         default=0.10_rk)
    call self%get_parameter(&
         self%Whet,'Whet','[m/day]',&
         'Rate of sinking of Het',&
         default=1.00_rk)
    call self%get_parameter(&
         self%Wphy_tot,'Wphy_tot','[1/day]',&
         'Total accelerated sinking with absorbed Mn hydroxides',&
          default=0.10_rk)
    call self%get_parameter(&
         self%Whet_tot,'Whet_tot','[1/day]',&
         'Total accelerated sinking with absorbed Mn hydroxides',&
          default=1.0_rk)
    call self%get_parameter(&
         self%Wsed_tot,'Wsed_tot','[1/day]',&
         'Total accelerated sinking with absorbed Mn hydroxides',&
          default=5.0_rk)
    
    !----Stoichiometric coefficients----!
    call self%get_parameter(self%r_n_p,'r_n_p','[-]','N[uM]/P[uM]', default=16.0_rk)
    call self%get_parameter(self%r_o_n,'r_o_n','[-]','O[uM]/N[uM]', default=6.625_rk)
    call self%get_parameter(self%r_c_n,'r_c_n','[-]','C[uM]/N[uM]', default=6.625_rk)
    call self%get_parameter(self%r_si_n,'r_si_n','[-]','Si[uM]/N[uM]', default=1.0_rk)
    
!   call self%get_parameter(self%transmodel, 'transmodel', 'na', 'Type of transport model', default=0.0_rk)
   call self%get_parameter(EPS, 'EPS', 'm^2/mg C', 'specific shortwave attenuation', default=2.208E-3_rk)     
    
    
    !Register state variables
    call self%register_state_variable(&
         self%id_Phy,'Phy','mmol/m**3','Phy',&
         minimum=0.0001_rk,initial_value=0.0001_rk)
    call self%register_state_variable(&
         self%id_Het,'Het','mmol/m**3','Het',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_POML,'POML','mmol/m**3','POML labile',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_POMR,'POMR','mmol/m**3','POMR semi-labile',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_DOML,'DOML','mmol/m**3','DOML labile',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_DOMR,'DOMR','mmol/m**3','DOMR semi-labile',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_O2,'O2','mmol/m**3','O2',minimum=0.0_rk)
   ! Register contribution to light extinction
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
        self%id_Phy,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &                                       
        self%id_Het,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,&
        self%id_POML,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,&
        self%id_POMR,scale_factor=EPS,include_background=.true.)
    
    !Register state dependencies
    call self%register_state_dependency(&
         self%id_PO4,'PO4','mmol/m**3','PO4')
    call self%register_state_dependency(&
         self%id_Si,'Si','mmol/m**3','Si')
    call self%register_state_dependency(&
         self%id_Sipart,'Sipart','mmol/m**3','Si particulate')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3','NO3')
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3','NH4')
    call self%register_state_dependency(&
         self%id_NO2,'NO2','mmol/m**3','NO2')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3','DIC')
    call self%register_state_dependency(&
         self%id_H2S,'H2S','mmol/m**3','H2S')
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_Baae,'Baae','mmol/m**3','aerobic autotrophic bacteria')
    call self%register_state_dependency(&
         self%id_Bhae,'Bhae','mmol/m**3','aerobic heterotrophic bacteria')
    call self%register_state_dependency(&
         self%id_Baan,'Baan','mmol/m**3','anaerobic aurotrophic bacteria')
    call self%register_state_dependency(&
         self%id_Bhan,'Bhan','mmol/m**3','anaerobic heterotrophic bacteria')
    
    !diagnostic dependency
    call self%register_dependency(&
         self%id_Hplus,'Hplus','mmol/m**3','H+ hydrogen')
    call self%register_dependency(self%id_Wadd,'Wadd','[1/day]',&
         'Additional sinking velocity via Mn4 adsorptoin')
    
    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPOML_O2,'DcPOML_O2','mmol/m**3',&
         'POML with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_O2,'DcPOMR_O2','mmol/m**3',&
         'POMR with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_O2,'DcDOMR_O2','mmol/m**3',&
         'DOMR with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_O2,'DcDOML_O2','mmol/m**3',&
         'DOM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortHet,'MortHet','mmol/m**3','Mortality of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Grazing,'Grazing','mmol/m**3','Grazing of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_RespHet,'RespHet','mmol/m**3','Respiration rate of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBhae,'GrazBhae','mmol/m**3','GrazBhae',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBhan,'GrazBhan','mmol/m**3','GrazBhan',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBaae,'GrazBaae','mmol/m**3','GrazBaae',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBaan,'GrazBaan','mmol/m**3','GrazBaan',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazPhy,'GrazPhy','mmol/m**3','GrazPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazPOP,'GrazPOP','mmol/m**3','GrazPOP',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazBact,'GrazBact','mmol/m**3','GrazBact',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortPhy,'MortPhy','mmol/m**3','MortPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ExcrPhy,'ExcrPhy','mmol/m**3','ExcrPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimNH4,'LimNH4','mmol/m**3','LimNH4',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimN,'LimN','mmol/m**3','LimN',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrowthPhy,'GrowthPhy','mmol/m**3','GrowthPhy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimT,'LimT','mmol/m**3','LimT',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimP,'LimP','mmol/m**3','LimP',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimNO3,'LimNO3','mmol/m**3','LimNO3',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimSi,'LimSi','mmol/m**3','LimSi',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimLight,'LimLight','mmol/m**3','LimLight',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_N_fixation,'N_fixation','mmol/m**3/d','N_fixation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_O2_rel_sat, &
         'O2_rel_sat','1','relative oxygen saturation', &
         standard_variable=standard_variables%fractional_saturation_of_oxygen)
    call self%register_diagnostic_variable(self%id_O2_sat, &
         'O2_sat','mmol O_2/m^3','oxygen saturation concentration')
    call self%register_diagnostic_variable(self%id_AOU, &
         'AOU','mmol O_2/m^3','Apparent Oxygen Utilization')
    call self%register_diagnostic_variable(&
         self%id_DcTOM_O2,'DcTOM_O2','mmol/m**3',&
         'Total OM_ oxidation with O2',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DOMTot,'DOMTot','mmol/m**3',&
         'DOMTot: refractory+labile',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_POMTot,'POMTot','mmol/m**3',&
         'POMTot: refractory+labile',output=output_time_step_integrated)
    
    !Register environmental dependencies
    call self%register_dependency(self%id_pres,standard_variables%pressure)
    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
    call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   if (self%use_aice) call self%register_horizontal_dependency(self%id_aice,type_horizontal_standard_variable(name='aice'))
    
    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
    

  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_bio),intent(in) :: self

        _DECLARE_ARGUMENTS_DO_
        real(rk):: temp,salt,pres,Iz
        real(rk):: NH4,NO2,NO3,PO4,Phy,Het,H2S,O2,Baae,Baan,Bhae,Bhan
        real(rk):: POML,POMR,DOML,DOMR,Si,Sipart,Alk,Hplus
        real(rk):: LimLight,LimT,LimP,LimNO3,LimNH4,LimN,LimSi
        real(rk):: GrowthPhy,MortPhy,ExcrPhy,dAlk,N_fixation
        real(rk):: GrazPhy,GrazPOP,GrazBaae,GrazBaan,GrazBhae
        real(rk):: GrazBhan,GrazBact,Grazing,RespHet,MortHet
        real(rk):: Autolysis_L,Autolysis_R,DcDOML_O2,DcPOML_O2,DcTOM_O2
        real(rk):: DcPOMR_O2,DcDOMR_O2,O2_sat,DOMTot,POMTot,AOU
        integer :: phy_t_dependence ! select dependence on T: (1) ERGOM; (2) for Arctic; (3) ERSEM
        !increments
        real(rk):: d_NO2,d_NO3,d_PO4,d_Si,d_DIC,d_O2,d_NH4
        real(rk):: d_Sipart,d_Phy,d_Het,d_Baae,d_Baan,d_Bhae,d_Bhan
        real(rk):: d_DOML,d_DOMR,d_POML,d_POMR,kf

        ! Enter spatial loops (if any)
        _LOOP_BEGIN_
            ! Retrieve current environmental conditions.
            _GET_(self%id_par,Iz) ! local photosynthetically active radiation
            _GET_(self%id_temp,temp) ! temperature
            _GET_(self%id_salt,salt) ! temperature
            _GET_(self%id_pres,pres) ! pressure in dbar
            ! Retrieve current (local) state variable values.
            !diagnostic
            _GET_(self%id_Hplus,Hplus)
            !state variables
            _GET_(self%id_NO2,NO2)
            _GET_(self%id_NO3,NO3)
            _GET_(self%id_PO4,PO4)
            _GET_(self%id_Si,Si)
            _GET_(self%id_Alk,Alk)
            _GET_(self%id_NH4,NH4)
            _GET_(self%id_H2S,H2S)
            _GET_(self%id_O2,O2)
            _GET_(self%id_Phy,Phy)
            _GET_(self%id_Het,Het)
            _GET_(self%id_Baae,Baae)
            _GET_(self%id_Baan,Baan)
            _GET_(self%id_Bhae,Bhae)
            _GET_(self%id_Bhan,Bhan)
            _GET_(self%id_POML,POML)
            _GET_(self%id_POMR,POMR)
            _GET_(self%id_DOML,DOML)
            _GET_(self%id_DOMR,DOMR)
            _GET_(self%id_Sipart,Sipart)

            !Phy
            !Influence of the Irradiance on photosynthesis
            LimLight = Iz/self%Iopt*exp(1._rk-Iz/self%Iopt)
            !Influence of Temperature on photosynthesis
            LimT = self%f_t(temp)
            !dependence of photosynthesis on P          
            LimP = yy(self%K_po4_lim*self%r_n_p,v_to_phy(PO4,Phy))
            !dependence of photosynthesis on Si
            LimSi = yy(self%K_si_lim/self%r_si_n,v_to_phy(Si,Phy))
            !dependence of photosynthesis on NO3+NO2
            LimNO3 = yy(self%K_nox_lim,v_to_phy(NO3+NO2,Phy))*&
            exp(-self%K_psi*v_to_phy(NH4,Phy)) 
            !dependence of photosynthesis on NH4
            LimNH4 = yy(self%K_nh4_lim,v_to_phy(NH4,Phy))*&
            (1._rk-exp(-self%K_psi*v_to_phy(NH4,Phy)))           
            !dependence of photosynthesis on N
            LimN = n_zero(min(1._rk,LimNO3+LimNH4))
            
            !Grouth of Phy (gross primary production in uM N)
            GrowthPhy = self%K_phy_gro*LimLight*LimT*min(LimP,LimN,LimSi)*n_zero(Phy)
            !Rate of mortality of phy
            MortPhy = Phy*(self%K_phy_mrt + thr_l(60._rk,O2,1._rk)*0.45_rk+&
                                            thr_l(20._rk,O2,1._rk)*0.45_rk) 
            !Excretion of phy
            ExcrPhy = self%K_phy_exc*Phy

            !Het
            !Grazing of Het on phy
            GrazPhy = self%K_het_phy_gro*Het*&
                yy(self%K_het_phy_lim,Phy/n_zero(Het))
            !Grazing of Het on detritus
            GrazPOP = self%K_het_pom_gro*Het*&
                yy(self%K_het_pom_lim,POML/n_zero(Het))
       
            !Grazing of Het on  bacteria            
            GrazBaae = 1.0_rk*self%graz(Baae,Het)            
            GrazBaan = 0.5_rk*self%graz(Baan,Het)
            GrazBhae = 1.0_rk*self%graz(Bhae,Het) 
            GrazBhan = 1.3_rk*self%graz(Bhan,Het) 
                        
            GrazBact =GrazBaae+GrazBaan+GrazBhae+GrazBhan
        
            !Total grazing of Het
            Grazing = GrazPhy+GrazPOP+GrazBact
            !Respiration of Het
            RespHet = self%K_het_res*Het*thr_h(20._rk,O2,1._rk)
            MortHet = Het*(self%K_het_mrt+thr_l(20._rk,O2,1._rk)*0.3_rk+&
            (0.5_rk+0.4_rk*tanh(H2S-self%s_hmort_H2S))*0.50_rk)

            !Nitrogen fixation described as appearence of NH4 available for
            !phytoplankton: N2 -> NH4 :
            N_fixation = self%K_nfix*LimP*&
            1._rk/(1._rk+((NO3+NO2+NH4)/n_zero(PO4)*16._rk)**4._rk)*GrowthPhy

            !POML and DOML (Savchuk, Wulff,1996)
            Autolysis_L = self%K_POML_DOML*POML
            Autolysis_R = self%K_POMR_DOMR*POMR
            
            !OM decay in N units for release of DIC and consumption of O2            
            !(CH2O)106(NH3)16H3PO4+106O2->106CO2+106H2O+16NH3+H3PO4
            
            kf = (O2/(O2+self%K_omox_o2))*&
                 (1._rk+self%beta_da*yy(self%tda,temp)) !koefficient   
            DcDOMR_O2 = self%K_DOMR_ox*DOMR*kf          
           
            DcDOML_O2 = self%K_DOML_ox*DOML*kf 
            DcDOMR_O2 = self%K_DOMR_ox*DOMR*kf 
            DcPOML_O2 = self%K_POML_ox*POML*kf            
            DcPOMR_O2 = self%K_POMR_ox*POMR*kf           
            DcTOM_O2  = DcPOMR_O2+DcDOMR_O2
      
      !components of temporal derivarives calculated in this module:
           d_POML = (-Autolysis_L-DcPOML_O2+MortPhy+MortHet+Grazing*&
                     (1._rk-self%Uz)*(1._rk-self%Hz)-GrazPOP)
            d_DOML = (Autolysis_L-DcDOML_O2+ExcrPhy+Grazing*(1._rk-self%Uz)*self%Hz)
            d_POMR = (DcPOML_O2-DcPOMR_O2-Autolysis_R)
            d_DOMR = (DcDOML_O2-DcDOMR_O2+Autolysis_R)
            d_NO2 = (-GrowthPhy*(LimNO3/LimN)*(n_zero(NO2)/n_zero(NO2+NO3)))
            d_NO3 = (-GrowthPhy*(LimNO3/LimN)*(n_zero(NO3)/n_zero(NO2+NO3)))
            d_PO4 = ((DcPOML_O2+DcDOML_O2-GrowthPhy+RespHet)/self%r_n_p)
            d_Si = ((-GrowthPhy+ExcrPhy)*self%r_si_n)            
            d_DIC = ((DcDOMR_O2+DcPOMR_O2-GrowthPhy+RespHet)*self%r_c_n)           
            d_O2 = ((-DcDOMR_O2-DcPOMR_O2+GrowthPhy-RespHet)*self%r_o_n)            
            d_NH4 = (DcPOML_O2+DcDOML_O2+RespHet+N_fixation-GrowthPhy*(LimNH4/LimN))            
            d_Sipart = ((MortPhy+GrazPhy)*self%r_si_n)
            d_Phy = (GrowthPhy-MortPhy-ExcrPhy-GrazPhy)
            d_Het = (self%Uz*Grazing-MortHet-RespHet)            
            d_Baae = -GrazBaae
            d_Baan = -GrazBaan
            d_Bhae = -GrazBhae
            d_Bhan = -GrazBhan 
            dAlk = 0.0_rk -d_PO4 -d_NO3 -d_NO2 +d_NH4           
                ! -1 mole per 1 mole of NO3- or NO2- or PO4-          
                ! +1 mole per 1 mole of NH4+ (Wollf-Gladrow, Zeebe,.. 2007)
            
            
            _SET_ODE_(self%id_POML,d_POML)            
            _SET_ODE_(self%id_DOML,d_DOML)            
            _SET_ODE_(self%id_POMR,d_POMR)            
            _SET_ODE_(self%id_DOMR,d_DOMR)            
            _SET_ODE_(self%id_NO2,d_NO2)            
            _SET_ODE_(self%id_NO3,d_NO3)            
            _SET_ODE_(self%id_PO4,d_PO4)            
            _SET_ODE_(self%id_Si,d_Si)            
            _SET_ODE_(self%id_DIC,d_DIC)            
            _SET_ODE_(self%id_O2,d_O2)            
            _SET_ODE_(self%id_NH4,d_NH4)            
            _SET_ODE_(self%id_Sipart,d_Sipart)            
            _SET_ODE_(self%id_Phy,d_Phy)            
            _SET_ODE_(self%id_Het,d_Het)            
            _SET_ODE_(self%id_Baae,d_Baae)
            _SET_ODE_(self%id_Baan,d_Baan)
            _SET_ODE_(self%id_Bhae,d_Bhae)            
            _SET_ODE_(self%id_Bhan,d_Bhan)            
            _SET_ODE_(self%id_Alk,dAlk)
      
      O2_sat = oxygen_saturation_concentration(temp,salt)
      POMTot=POML+POMR
      DOMTot=DOML+DOMR

      
      _SET_DIAGNOSTIC_(self%id_O2_sat,O2_sat)
      _SET_DIAGNOSTIC_(self%id_O2_rel_sat,max(0.0_rk,100.0_rk*O2/O2_sat))
      _SET_DIAGNOSTIC_(self%id_AOU,(O2_sat-O2))
      _SET_DIAGNOSTIC_(self%id_DcPOML_O2,DcPOML_O2)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_O2,DcPOMR_O2)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_O2,DcDOMR_O2)
      _SET_DIAGNOSTIC_(self%id_DcDOML_O2,DcDOML_O2)
      _SET_DIAGNOSTIC_(self%id_MortHet,MortHet)
      _SET_DIAGNOSTIC_(self%id_Grazing,Grazing)
      _SET_DIAGNOSTIC_(self%id_RespHet,RespHet)
      _SET_DIAGNOSTIC_(self%id_GrazBhae,GrazBhae)
      _SET_DIAGNOSTIC_(self%id_GrazBhan,GrazBhan)
      _SET_DIAGNOSTIC_(self%id_GrazBaae,GrazBaae)
      _SET_DIAGNOSTIC_(self%id_GrazBaan,GrazBaan)
      _SET_DIAGNOSTIC_(self%id_GrazPhy,GrazPhy)
      _SET_DIAGNOSTIC_(self%id_GrazPOP,GrazPOP)
      _SET_DIAGNOSTIC_(self%id_GrazBact,GrazBact)
      _SET_DIAGNOSTIC_(self%id_MortPhy,MortPhy)
      _SET_DIAGNOSTIC_(self%id_ExcrPhy,ExcrPhy)
      _SET_DIAGNOSTIC_(self%id_LimNH4,LimNH4)
      _SET_DIAGNOSTIC_(self%id_LimN,LimN)
      _SET_DIAGNOSTIC_(self%id_GrowthPhy,GrowthPhy)
      _SET_DIAGNOSTIC_(self%id_LimT,LimT)
      _SET_DIAGNOSTIC_(self%id_LimP,LimP)
      _SET_DIAGNOSTIC_(self%id_LimNO3,LimNO3)
      _SET_DIAGNOSTIC_(self%id_LimSi,LimSi)
      _SET_DIAGNOSTIC_(self%id_LimLight,LimLight)
      _SET_DIAGNOSTIC_(self%id_DcTOM_O2,DcTOM_O2)
      _SET_DIAGNOSTIC_(self%id_DOMTot,DOMTot)
      _SET_DIAGNOSTIC_(self%id_POMTot,POMTot)
      _SET_DIAGNOSTIC_(self%id_N_fixation,N_fixation)
    _LOOP_END_
  end subroutine do
 
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    class (type_niva_brom_bio),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_SURFACE_
    real(rk)                   :: O2, temp, salt, windspeed
    real(rk)                   :: Ox, Oa, TempT, Obe, Q_O2

    _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_O2,O2)
      _GET_(self%id_temp,temp) ! temperature
      _GET_(self%id_salt,salt) ! salinity
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)

      Ox = 1953.4_rk-128._rk*temp+3.9918_rk*temp*temp-&
           0.050091_rk*temp*temp*temp !(Wanninkoff, 1992)
      if (Ox>0._rk) then
        Oa = 0.028_rk*(windspeed**3._rk)*sqrt(400._rk/Ox)
      else
        Oa = 0._rk
      end if

      !Calculation of O2 saturation Obe according to UNESCO, 1986
      TempT = (temp+273.15_rk)/100._rk
      Obe = exp(-173.4292_rk+249.6339_rk/TempT+143.3483_rk*&
            log(TempT)-21.8492_rk*TempT+salt*(-0.033096_rk+&
            0.014259_rk*TempT-0.0017_rk*TempT*TempT)) !O2_sat
      Obe = Obe*1000._rk/22.4_rk !convert from ml/l into uM
      Q_O2 = windspeed*(Obe-O2) !After (Burchard et al., 2005)

      _SET_SURFACE_EXCHANGE_(self%id_O2,Q_O2)
    _HORIZONTAL_LOOP_END_
  end subroutine
  !
  !Phy temperature limiter
  !
  real(rk) function f_t(self,temperature)
    class (type_niva_brom_bio),intent(in) :: self
    real(rk)                  ,intent(in):: temperature

    real(rk):: bm
    real(rk):: cm
    real(rk):: t_0           !reference temperature
    real(rk):: temp_aug_rate !temperature augmentation rate
    real(rk):: q10       !Coefficient for uptake rate dependence on t
    real(rk):: t_upt_min !Low  t limit for uptake rate dependence on t
    real(rk):: t_upt_max !High t limit for uptake rate dependence on t



 
    if (self%phy_t_dependence == 1) then
      ! ERGOM
      bm = 0.12_rk
      cm = 1.4_rk
      f_t = exp(bm*temperature-cm)
    else if (self%phy_t_dependence == 2) then
      !for Arctic (Moore et al.,2002; Jin et al.,2008)
      t_0           = 0._rk
      temp_aug_rate = 0.0663_rk
      f_t = exp(temp_aug_rate*(temperature-t_0))
    else if (self%phy_t_dependence == 3) then
      ! ERSEM
      q10       = 2.0_rk
      t_upt_min = 10.0_rk
      t_upt_max = 32.0_rk
      f_t = q10**((temperature-t_upt_min)/10._rk)-&
            q10**((temperature-t_upt_max)/3._rk)
    end if
!   Some others:
 !  LimT     = 0.5(1+tanh((t-tmin)/smin)) (1-0.5(1+th((t-tmax)/smax))) !Smin= 15  Smax= 15  Tmin=  10 Tmax= 35   (Deb et al., .09)
 !  LimT     = exp(self%bm*temp-self%cm))        !Dependence on Temperature (used in (Ya,So, 2011) for Arctic)  
 !  LimT     = 1./(1.+exp(10.-temp))             !Dependence on Temperature (ERGOM for cya)
 !  LimT     = 1.-temp*temp/(temp*temp +12.*12.) !Dependence on Temperature (ERGOM for dia)
 !  LimT     = 2.**((temp- 10.)/10.) -2**((temp-32.)/3.) !(ERSEM)
 !  LimT     =q10*(T-20)/10 !Q10=1.88 (Gregoire, 2000)       
  end function f_t
  !
  !adapted from ERSEM
  !
  function oxygen_saturation_concentration(ETW,X1X) result(O2_sat)
    real(rk),                      intent(in) :: ETW,X1X
    real(rk)                                  :: O2_sat

    real(rk),parameter :: A1 = -173.4292_rk
    real(rk),parameter :: A2 = 249.6339_rk
    real(rk),parameter :: A3 = 143.3483_rk
    real(rk),parameter :: A4 = -21.8492_rk
    real(rk),parameter :: B1 = -0.033096_rk
    real(rk),parameter :: B2 = 0.014259_rk
    real(rk),parameter :: B3 = -0.0017_rk
    real(rk),parameter :: R = 8.3145_rk
    real(rk),parameter :: P = 101325_rk
    real(rk),parameter :: T = 273.15_rk

    ! volume of an ideal gas at standard temp (25C) and pressure (1 atm)
    real(rk),parameter :: VIDEAL = (R * 298.15_rk / P) *1000._rk

    real(rk)           :: ABT

    ! calc absolute temperature
    ABT = ETW + T

    ! calc theoretical oxygen saturation for temp + salinity
    ! From WEISS 1970 DEEP SEA RES 17, 721-735.
    ! units of ln(ml(STP)/l)
    O2_sat = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
            + A4 * (ABT/100._rk) &
            + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

    ! convert units to ml(STP)/l then to mMol/m3
    O2_sat = exp( O2_sat )
    O2_sat = O2_sat * 1000._rk / VIDEAL
  end function
    
    real(rk) function n_zero(var)
        real(rk),intent(in):: var
        n_zero = max(var,1.e-10_rk)
    end function n_zero
 
    real(rk) function graz(self,var,Het)
        class (type_niva_brom_bio),intent(in) :: self
        real(rk),intent(in):: var,Het
        !Het = n_zero(Het)
        graz = self%K_het_bac_gro*Het*yy(self%limGrazBac,var/n_zero(Het))   
    end function graz
    
    real(rk) function yy(a,x)
        ! Squared Michaelis-Menten type of limiter
        ! Original author(s): Hans Burchard, Karsten Bolding
        real(rk),intent(in):: a,x
        yy=x**2._rk/(a**2._rk+x**2._rk)
    end function yy
    
    
    real(rk) function v_to_phy(var,Phy)
    real(rk),intent(in):: var, Phy
        v_to_phy = n_zero(var)/n_zero(Phy)
    end function v_to_phy    

    real(rk) function thr_h(threshold_value,var_conc,koef)
        ! Threshold value for the reaction 
        ! koef 1 gives regular tgh function 
        ! 0.1 - smooth function 
        real(rk), intent(in) :: threshold_value,var_conc,koef
        thr_h = 0.5+0.5*tanh((var_conc-threshold_value)*koef)
    end function 
          
    real(rk) function thr_l(threshold_value,var_conc,koef)
        ! Threshold value for the reaction 
        real(rk), intent(in) :: threshold_value,var_conc,koef
        thr_l = 0.5-0.5*tanh((var_conc-threshold_value)*koef)
    end function 

  ! Set increased manganese sinking via MnIV and MnIII oxides formation
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
     class (type_niva_brom_bio), intent(in) :: self
     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
     
     real(rk) :: Wadd, Wphy_tot, Whet_tot, Wsed_tot
          
     _LOOP_BEGIN_
  
      _GET_(self%id_Wadd,Wadd)
     
      Wphy_tot = self%Wphy + 0.25_rk * Wadd
      Whet_tot = self%Whet + 0.5_rk * Wadd
      Wsed_tot = self%Wsed + Wadd
      
      _ADD_VERTICAL_VELOCITY_(self%id_Phy, Wphy_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_Het, Whet_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_POML, Wsed_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_POMR, Wsed_tot)
  
     _LOOP_END_
  
  end subroutine get_vertical_movement
end module fabm_niva_brom_bio
