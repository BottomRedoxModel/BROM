!-----------------------------------------------------------------------
! OXYDEP is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Jorn Bruggemann, Anfisa Berezina
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_niva_oxydep --- OXYDEP biogeochemical model based upon
! Yakushev et al, 2013  with  modifications by Jorn Bruggeman and Anfisa Berezina
! and adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module fabm_niva_oxydep
!
! !DESCRIPTION:
!
! OXYDEP targests on the silmplest possible way of parameterization of the oxygen  (DO) fate in changeable redox conditions.
! It has a simplified ecosystem, and simulates production of DO due to photosynthesis and consumation of DO for biota respiraion,
! OM mineralization, nitrification, and oxidation of reduced specied of S, Mn, Fe, present in suboxic conditions.
! OXYDEP consists of 6 state variables ( in N-units):
! - Phy - all the phototrophic organisms (phytoplankton and bacteria).
!   Phy grows due to photosynthesis, loses inorganic matter
!   due to respiraion, and loses organic matter in dissolved (DOM) and particulate (POM)
!   forms due to metabolism and mortality. Phy growth is limited by irradiance, temperature and NUT availability.
! - Het - heterotrophs, can consume Phy and POM,  produce DOM and POM and respirate NUT.
! - NUT - represents oxydized forms of nutrients (i.e. NO3 and NO2 for N),
!   that doesn't need additional  oxygen for nitrification.
! - DOM - is dissolved organic matter. DOM  includes all kinds of labile dissolved organic matter
!   and reduced forms of inorganic nutrients (i.e. NH4 and Urea for N).
! - POM - is particular organic matter (less labile than DOM). Temperature affects DOM and POM mineralization.
! - Oxy - is dissolved oxygen.
! For the details of  OxyDEP  implemented here (actually, a previous version withour Het), see (Yakushev et al, 2013)
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Jorn Bruggeman, Anfisa Berezina
!
	real(rk), parameter :: pi=3.141592653589793_rk
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_oxydep
!     Variable identifiers
      type (type_state_variable_id)        :: id_oxy,id_phy,id_het,id_nut,id_pom,id_dom
      type (type_state_variable_id)        :: id_dic, id_alk
    !global dependencies
      type (type_dependency_id)            :: id_par,id_temp, id_salt
      type (type_horizontal_dependency_id) :: id_windspeed, id_aice
      type (type_horizontal_dependency_id) :: id_lat
      type (type_global_dependency_id)     :: id_yday            
      logical   :: use_aice ! added ice area limitation of air-sea flux, with new switch parameter use_aice.
      
      type (type_diagnostic_variable_id)   :: id_MortHet,id_RespHet,id_GrazPhy,id_GrazPOM,id_GrowthPhy,id_MortPhy,id_ExcrPhy,id_RespPhy
      type (type_diagnostic_variable_id)   :: id_DOM_decay_ox,id_DOM_decay_denitr,id_POM_decay_ox,id_POM_decay_denitr
      type (type_diagnostic_variable_id)   :: id_LimT,id_LimP,id_LimN,id_LimLight, id_N_fixation, id_Autolys
!     Model parameters
      !----Phy -----------!
       real(rk) :: Max_uptake, Knut, r_phy_nut, r_phy_pom, r_phy_dom, r_phy_pom_anox
       real(rk) :: Iopt, alphaI, betaI, gammaD, O2_add_mor_phy
       integer  :: phy_t_dependence     ! dependence on T: (1) ERGOM; (2) for Arctic; (3) ERSEM 
       integer  :: phy_light_dependence ! dependence on light: (1) Steel; (2) Platt  
       integer  :: phy_light_daylength  ! daylength matters (1) or not (0)
      !----Het -----------!
       real(rk) :: r_phy_het, Kphy, r_pop_het, Kpop, r_het_nut, r_het_pom, Uz, Hz, r_het_pom_anox
      !----DOM, POM   --- !
       real(rk) :: r_pom_dom, r_pom_nut_oxy, r_pom_nut_nut, r_dom_nut_oxy, r_dom_nut_nut, Tda, beta_da
      ! real(rk) :: r_om_nut_sul  = 0.005    ! Specific rate of OM anoxic decay   [1/d]
      !----Oxy -----------!
       real(rk) :: O2_suboxic
      !  Lower boundary
       real(rk) :: Bu,Trel,b_ox,b_dom_ox,b_dom_anox,b_nut
      ! Upper boundary for oxygen flux calculations
       real(rk) :: l_bubble, s_bubble
      !  Stochiometric coefficients
       real(rk) :: NtoB, OtoN, NtoN, CtoN
      !---sinking-------------------------------!
       real(rk) :: Wphy, Whet, Wpom
       !--Type of transport model
       real(rk) :: transmodel
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
      procedure :: f_t
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OxyDep model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the OXYDEP namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_niva_oxydep), intent(inout), target :: self
   integer,                  intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev, Jorn Bruggemann, Anfisa Berezina
!
! !LOCAL VARIABLES:
   real(rk),parameter ::  d_per_s = 1.0_rk/86400.0_rk
   real(rk) :: EPS
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
  ! Phy
   call self%get_parameter(self%Max_uptake,'Max_uptake','1/d', 'Maximum nutrient uptake rate',                               default=5.0_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%Knut,      'Knut',      'nd',  'Half-sat.const. for uptake of NUT by Phy for NUT/Phy ratio', default=0.1_rk)
   call self%get_parameter(self%Iopt,      'Iopt',      'Watt/m2/h',              'Optimal irradiance',                   default=25.0_rk)
   call self%get_parameter(self%alphaI,    'alphaI',    '1/d/(Watt/m2)',          'Initial slope of PI-curve',            default=0.767_rk,scale_factor=d_per_s)
   call self%get_parameter(self%betaI,     'betaI',     '1/d/(Watt/m2)',          'Photoinhibition parameter',            default=0.013_rk,scale_factor=d_per_s)
   call self%get_parameter(self%phy_light_dependence, 'phy_light_dependence', '-',   'light dependence for Phy growth',      default=1)
   call self%get_parameter(self%phy_light_daylength,  'phy_light_daylength',  '-',   'influence of daylength at Phy growth', default=1)
   call self%get_parameter(self%gammaD,    'gammaD',    '-',                         'Adaptation to daylength parameter',    default=0.5_rk)
   call self%get_parameter(self%phy_t_dependence,'phy_t_dependence','-','T dependence for Phy growth',                       default=1)
   call self%get_parameter(self%r_phy_nut,     'r_phy_nut',     '1/d', 'Specific Phy  respiration rate',                     default=0.04_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_pom,     'r_phy_pom',     '1/d', 'Specific Phy rate of mortality',                     default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_dom,     'r_phy_dom',     '1/d', 'Specific Phy rate of excretion',                     default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_pom_anox, 'r_phy_pom_anox', '1/d', 'Specific additional Phy mortality in suboxic/anoxic conditions', default=0.4_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_het_pom_anox, 'r_het_pom_anox', '1/d', 'Specific additional Het mortality in suboxic/anoxic conditions', default=0.4_rk,scale_factor=d_per_s)
   call self%get_parameter(self%O2_add_mor_phy, 'O2_add_mor_phy', 'mmol/m3', 'Threshold O2 value for additional Phy mortality', default=20._rk)
  ! Het
   call self%get_parameter(self%r_phy_het, 'r_phy_het', '1/d', 'Max.spec. rate of grazing of Het on Phy',                    default=2.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Kphy,      'Kphy',      'nd',  'Half-sat.const.for grazing of Het on Phy for Phy/Het ratio ',default=0.1_rk)
   call self%get_parameter(self%r_pop_het, 'r_pop_het', '1/d', 'Max.spec. rate of grazing of Het on POM ',                   default=0.7_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Kpop,      'Kpop',      'nd',  'Half-sat.const.for grazing of Het on POM for POM/Het ratio ',default=2._rk)
   call self%get_parameter(self%r_het_nut, 'r_het_nut', '1/d', 'Specific Het respiration rate ',                             default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_het_pom, 'r_het_pom', '1/d', 'Specific Het mortality rate ',                               default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Uz,        'Uz',        'nd',  'Food absorbency for Het',                                    default=0.5_rk)
   call self%get_parameter(self%Hz,        'Hz',        'nd',  'Ratio betw. diss. and part. excretes of Het ',               default=0.5_rk)
  ! POM
   call self%get_parameter(self%r_pom_dom,     'r_pom_dom',     '1/d', 'Specific rate of POM decomposition ',                default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_pom_nut_oxy, 'r_pom_nut_oxy', '1/d', 'Specific rate of POM oxic decay  ',                  default=0.03_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_pom_nut_nut, 'r_pom_nut_nut', '1/d', 'Specific rate of POM denitrification  ',             default=0.006_rk,scale_factor=d_per_s)
  ! DOM
   call self%get_parameter(self%r_dom_nut_oxy, 'r_dom_nut_oxy', '1/d', 'Specific rate of DOM oxic decay  ',                  default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_dom_nut_nut, 'r_dom_nut_nut', '1/d', 'Specific rate of DOM denitrification  ',             default=0.002_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Tda,           'Tda',           'nd',  'Coefficient for dependence of mineralization on t ', default=13._rk)
   call self%get_parameter(self%beta_da,       'beta_da',       'nd',  'Coefficient for dependence of mineralization on t ', default=20._rk)
  ! O2
  !  Lower boundary
   call self%get_parameter(self%Bu,        'Bu',         'nd',      'Burial coeficient for lower boundary',                default=0.25_rk)
   call self%get_parameter(self%Trel,      'Trel',       's/m',     'Relaxation time for exchange with teh sediments',     default=1e5_rk)
   call self%get_parameter(self%O2_suboxic,    'O2_suboxic',     'mmol/m3', 'Limiting O2 value for oxic/suboxic switch',   default=40._rk)
   call self%get_parameter(self%b_ox,      'b_ox',       'mmol/m3', 'Oxy in the sediment',                                 default=0._rk)
   call self%get_parameter(self%b_dom_ox,  'b_dom_ox',   'mmol/m3', 'OM in the sediment (oxic conditions)',                default=2._rk)
   call self%get_parameter(self%b_dom_anox,'b_dom_anox', 'mmol/m3', 'OM in the sediment (anoxic conditions) ',             default=6._rk)
   call self%get_parameter(self%b_nut,     'b_nut',      'mmol/m3', 'NUT in the sediment',                                 default=0._rk)
  ! Upper boundary for oxygen flux calculations
   call self%get_parameter(self%l_bubble,   'l_bubble',   '[ - ]',   'large bubbles effect',                              default=1._rk)
   call self%get_parameter(self%s_bubble,   's_bubble',   '[ - ]',   'small bubbles effect',                              default=1._rk)
   call self%get_parameter(self%use_aice,   'use_aice',   '[ - ]',   'use ice area to limit air-sea flux',                default=.false.)
  !---sinking-------------------------------!
   call self%get_parameter(self%Wphy,     'Wphy',     'm/s', 'vertical velocity of Phy (<0 for sinking)',           default=-0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Whet,     'Whet',     'm/s', 'vertical velocity of het (<0 for sinking)',           default=-0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Wpom,     'Wpom',     'm/s', 'vertical velocity of POM (<0 for sinking)',           default=-1.0_rk,scale_factor=d_per_s)
  !  Stochiometric coefficients
   call self%get_parameter(self%NtoB,     'NtoB',         'uM(N)/mgWW/m3', 'N[uM]/BIOMASS [mg/m3]',                        default=0.016_rk)
   call self%get_parameter(self%OtoN,     'OtoN',         'uM(O)/uM(N)',   'Redfield (138/16) to NO3',                     default=8.625_rk)
   call self%get_parameter(self%CtoN,     'CtoN',         'uM(C)/uM(N)',   'Redfield (106/16) to NO3',                     default=6.625_rk)
   call self%get_parameter(self%NtoN,     'NtoN',         'uM(N)/uM(N)',   'Richards denitrification (84.8/16.)',          default=5.3_rk)
   
   call self%get_parameter(self%transmodel, 'transmodel', 'na', 'Type of transport model', default=0.0_rk)
   call self%get_parameter(EPS, 'EPS', 'm^2/mg C', 'specific shortwave attenuation', default=2.208E-3_rk)     

   ! Register state variables
   call self%register_state_variable(self%id_oxy,'Oxy','mmol/m**3','Oxy: Oxygen',  150.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_phy,'Phy','mmol/m**3','Phy: Autotrophs/Phytoplankton',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wphy)
   call self%register_state_variable(self%id_nut,'NUT','mmol/m**3','NUT: Nutrient',  1.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pom,'POM','mmol/m**3','POM: Particulate organic matter',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wpom)
   call self%register_state_variable(self%id_dom,'DOM','mmol/m**3','DOM: Dissolved organic matter',  0.1_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_het,'Het','mmol/m**3','Het: Heterotrophs/Zooplankton',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Whet)
   ! Register the contribution of all state variables to total nitrogen
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phy)
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_nut)
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_pom)
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_dom)
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_het)

   ! Register contribution to light extinction
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
        self%id_phy,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &                                       
        self%id_het,scale_factor=EPS,include_background=.true.)
   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,&
        self%id_pom,scale_factor=EPS,include_background=.true.)

   ! Register link to external DIC and/or Alk pool, if DIC variable name is provided in namelist. 
! if (_AVAILABLE_(self%id_dic)) then
   call self%register_state_dependency(self%id_dic,'DIC','mmol/m**3','total dissolved inorganic carbon',required=.false.)
! endif   
! if (_AVAILABLE_(self%id_alk)) then
   call self%register_state_dependency(self%id_alk,'Alk','mmol/m**3','total alkalinity',required=.false.)
! endif
   ! Register diagnostic variables
call self%register_diagnostic_variable(self%id_MortHet,'MortHet','mmol/m**3/d',  'MortHet,  Mortality of Het',           &
                    output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_RespHet,'RespHet','mmol/m**3/d',  'RespHet, Respiration rate of Het',     &
                output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazPhy,'GrazPhy','mmol/m**3/d',  'GrazPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazPOM,'GrazPOM','mmol/m**3/d',  'GrazPOM',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_MortPhy,'MortPhy','mmol/m**3/d',  'MortPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_ExcrPhy,'ExcrPhy','mmol/m**3/d',  'ExcrPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_RespPhy,'RespPhy','mmol/m**3/d',  'RespPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimN,'LimN','nd',  'LimN',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrowthPhy,'GrowthPhy','mmol/m**3/d', 'GrowthPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimT,'LimT','nd',  'LimT',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimLight,'LimLight','nd',  'LimLight',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_DOM_decay_ox,'DOM_decay_ox','mmol/m**3/d',  'DOM_decay_ox',             &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_DOM_decay_denitr,'DOM_decay_denitr','mmol/m**3/d',  'DOM_decay_denitr', &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_POM_decay_ox,'POM_decay_ox','mmol/m**3/d',  'POM_decay_ox',             &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_POM_decay_denitr,'POM_decay_denitr','mmol/m**3/d',  'POM_decay_denitr', &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_Autolys,'Autolys','mmol/m**3/d',  'Autolys', &
            output=output_time_step_integrated)
!call self%register_diagnostic_variable(self%id_N_fixation,'N_fixation','mmol/m**3/d',  'N_fixation',           &
!            output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   if (self%use_aice) call self%register_horizontal_dependency(self%id_aice,type_horizontal_standard_variable(name='aice'))
 !  call self%register_horizontal_dependency(self%id_day_length,type_horizontal_standard_variable(name='day_length'))
       ! Register environmental dependencies (lat,day length,day of year)
   call self%register_horizontal_dependency(self%id_lat,standard_variables%latitude)
   call self%register_global_dependency(self%id_yday,standard_variables%number_of_days_since_start_of_the_year)  
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
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) :: oxy, nut, pom, dom, phy, het, t
   real(rk) :: doxy, dnut, ddom, dpom, dphy, dhet
   real(rk) :: iopt, daylength, sundec, yday, lat
   real(rk) :: DIC, Alk, dAlk, dDIC
 ! Rates of biogeochemical processes
 ! Phy
   real(rk) :: GrowthPhy               ! Nutrient uptake rate (1/d)
   real(rk) :: Iz                      ! Irradiance at certain depth
   real(rk) :: LimLight                ! Photosynthesis dependencs on irradiance
   real(rk) :: LimT                    ! Photosynthesis dependencs on temperature
   real(rk) :: LimN                    ! Photosynthesis dependencs on nutrient
   real(rk) :: MortPhy                 ! Mortality of Phy (1/d)
   real(rk) :: ExcrPhy                 ! Excretion of Phy (1/d)
   real(rk) :: RespPhy                 ! Respiration of Phy (1/d)
 ! Het
   real(rk) :: GrazPhy, GrazPOM, RespHet, MortHet
 ! POM
   real(rk) :: Autolys                 ! Autolysis of POM to DOM (1/d)
   real(rk) :: POM_decay_ox            ! oxic mineralization of POM and ammonification (1/d)
   real(rk) :: POM_decay_denitr        ! suboxic mineralization of POM (denitrification) (1/d)
 ! DOM
   real(rk) :: DOM_decay_ox            ! oxic mineralization of POM and ammonification (1/d)
   real(rk) :: DOM_decay_denitr        ! suboxic mineralization of POM (denitrification) (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_phy,phy)
   _GET_(self%id_het,het)
   _GET_(self%id_pom,pom)
   _GET_(self%id_nut,nut)
   _GET_(self%id_dom,dom)
  if (_AVAILABLE_(self%id_dic)) then   
    _GET_(self%id_dic,DIC)
  else  
    DIC=0.0_rk
  endif
  if (_AVAILABLE_(self%id_alk)) then
   _GET_(self%id_alk,Alk)
  else  
   Alk=0.0_rk
  endif    
   ! Retrieve current environmental conditions.
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
   _GET_(self%id_temp,t)              ! temperature 
   _GET_HORIZONTAL_(self%id_lat,lat)  ! latitude in degreese  
   _GET_GLOBAL_(self%id_yday,yday)    ! Retrieve time since beginning of the year

!--------------------------------------------------------------
! Phy
!--------------------------------------------------------------
   ! Growth of Phy and uptake of NUT   
   if (self%phy_light_dependence == 1) then   !Dependence on Irradiance 
     LimLight = Iz/self%Iopt*exp(1-Iz/self%Iopt)                 ! Steel     
   else if (self%phy_light_dependence == 2) then
     LimLight = (1._rk-exp(-self%alphaI*Iz/(self%Max_uptake))) & ! units are 1/s !
                *exp(-self%betaI*Iz/(self%Max_uptake)) ! (Platt et al., 1980) 
   endif  
   if (self%phy_light_daylength == 1) then  ! correction for adaptation to day length (cf. Gilstad and Sakshaug 1990)
   ! Earth declination (converted to radians)
     sundec= (23.45_rk*sin((360*(284+yday)/365)*pi/180._rk))*pi/180._rk
   ! Calculate day length as the day' share (i.e. /24)
     daylength = acos(min(1._rk,max(-1._rk,-tan(lat*pi/180._rk)*tan(sundec))))*180_rk/pi*2._rk/15._rk/24._rk
     LimLight = LimLight*(self%gammaD + 0.5_rk)/(self%gammaD + daylength)
   endif
  
      LimT     = self%f_t(t)                         !Dependence on Temperature
      LimN     = yy(self%Knut,nut/(max(0.0001,phy))) !Dependence on nutrients

     GrowthPhy = self%Max_uptake*LimLight*LimT*LimN*phy

   ! Respiraion of Phy and increase of NUT
     RespPhy=self%r_phy_nut*phy

   ! Methabolism of Phy and increase of DOM
     ExcrPhy=self%r_phy_dom*phy

   ! Mortality of Phy and increase of POM
     MortPhy = phy*(self%r_phy_pom &
   ! Additional  mortaliny in low oxygen conditions
              + F_subox(oxy,self%O2_add_mor_phy)*self%r_phy_pom_anox)

!--------------------------------------------------------------
! Het
!--------------------------------------------------------------
   ! Grazing of Het on Phy
    GrazPhy = het*self%r_phy_het*yy(self%Kphy,(max(0.0_rk,phy-0.01))/max(het,0.0001))

   ! Grazing of Het on POM
    GrazPOM = self%r_pop_het*het*yy(self%Kpop,(max(0.0_rk,pom-0.01))/max(het,0.0001))

   ! Respiraion of Het and increase of NUT
    RespHet = self%r_het_nut*het 

   ! Mortality of Het and increase of POM
    MortHet = het*(self%r_het_pom+F_subox(oxy,self%O2_suboxic)*self%r_het_pom_anox) 

!--------------------------------------------------------------
! POM
!--------------------------------------------------------------
   ! Decomposition of POM and increase of DOM
    Autolys = self%r_pom_dom*pom

   ! Oxic mineralization of POM and ammonification depend on T
    POM_decay_ox   = self%r_pom_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*pom
   ! Suboxic mineralization of OM (denitrification and anammox)
   ! denitrification: (CH2O)106(NH3)16H3PO4 + 84.8HNO3 -> 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4 (Richards, 1965)
   ! anammox: NO2- + NH4+ -> N2 + 2H2O (Canfield,2005)
    POM_decay_denitr = self%r_pom_nut_nut*pom &           ! depends on NUT (NO3+NO2) and DOM (NH4+Urea+"real"DON)
                      *(1.+self%beta_da*yy(self%tda,t)) & ! depends on T
                      *F_subox(oxy,self%O2_suboxic) &     ! starts when O2<O2_suboxic
                      *F_ox(nut,0.01_rk)                  ! stops at NUT<0.01 

!--------------------------------------------------------------
! DOM
!--------------------------------------------------------------
! Oxic mineralization of DOM and ammonification depend on T
    DOM_decay_ox   = self%r_dom_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*dom
! Suboxic mineralization of OM (denitrification and anammox)
    DOM_decay_denitr = self%r_dom_nut_nut*dom &           ! depends on NUT (NO3+NO2) and DOM (NH4+Urea+"real"DON)
                      *(1.+self%beta_da*yy(self%tda,t)) & ! depends on T
                      *F_subox(oxy,self%O2_suboxic) &     ! starts when O2<O2_suboxic
                      *F_ox(nut,0.01_rk)                  ! stops at NUT<0.01 

   ! Now we can summarize processes and write state variables sink/sources:

!--------------------------------------------------------------
! Oxy
!--------------------------------------------------------------
   ! Changes of Oxy due to OM production and decay!
    doxy = -self%OtoN*(POM_decay_ox + DOM_decay_ox &
   !        + (POM_decay_denitr + DOM_decay_denitr) & !denitrification doesn't change oxygen
           - GrowthPhy + RespPhy + RespHet) &
   ! additional consumption of Oxy due to oxidation of reduced froms of S,Mn,Fe etc.
   ! in suboxic conditions equals consumption for NH4 oxidation (Yakushev et al, 2008 in Nauka Kubani)
           - F_subox(oxy,self%O2_suboxic)*self%OtoN*DOM_decay_ox  ! starts when O2<O2_suboxic
!--------------------------------------------------------------
! NUT
!--------------------------------------------------------------
   dnut = -GrowthPhy+RespPhy+RespHet+POM_decay_ox+DOM_decay_ox-self%NtoN*(POM_decay_denitr+DOM_decay_denitr)
            ! Decrease of NUT (as NO3+NO2) due to denitrification
   ddom = ExcrPhy+Autolys-DOM_decay_ox+(GrazPhy+GrazPOM)*(1.-self%Uz)*self%Hz+POM_decay_denitr 
            ! Denitrification of "real" DOM into NH4 (DOM_decay_denitr) will not change state variable DOM 
   dpom = MortPhy-Autolys-POM_decay_ox-POM_decay_denitr+(GrazPhy+GrazPOM)*(1.-self%Uz)*(1.-self%Hz)-GrazPOM+MortHet
   dphy = GrowthPhy-RespPhy-ExcrPhy-MortPhy-GrazPhy
   dhet = self%Uz*(GrazPhy+GrazPOM)-MortHet-RespHet
!--------------------------------------------------------------

!derivatives for FABM
   _SET_ODE_(self%id_oxy,doxy)
   _SET_ODE_(self%id_nut,dnut)
   _SET_ODE_(self%id_dom,ddom)
   _SET_ODE_(self%id_pom,dpom)
   _SET_ODE_(self%id_phy,dphy)
   _SET_ODE_(self%id_het,dhet)
   
   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
 if (_AVAILABLE_(self%id_dic)) then
     dDIC=self%CtoN*(-GrowthPhy+RespPhy+RespHet+POM_decay_ox+DOM_decay_ox)
   _SET_ODE_(self%id_dic,dDIC)       
 endif
   ! If an externally maintained Alk pool is present, change the Alk 
 if (_AVAILABLE_(self%id_alk)) then    
     dAlk=( +GrowthPhy+self%NtoN*(POM_decay_denitr+DOM_decay_denitr)&   ! Alk increase if NO3 consumes
      -RespPhy-RespHet-POM_decay_ox & ! Alk decrease if NO3 produces
      -2.0_rk*DOM_decay_ox & !nitrification leads to 2 X Alk decrease
      +POM_decay_denitr    & !denitrification leads to 1 X Alk decrease (Wolf-Gladrow,2007)
      +ExcrPhy+Autolys+(GrazPhy+GrazPOM)*(1.-self%Uz)*self%Hz &
          )                 ! Alk increases with an increase of NH4+DOM (i.e. we assume here DOM=orgAlk)
   _SET_ODE_(self%id_alk,dAlk)
endif
   ! Export diagnostic variables
_SET_DIAGNOSTIC_(self%id_MortHet,MortHet)
_SET_DIAGNOSTIC_(self%id_RespHet,RespHet)
_SET_DIAGNOSTIC_(self%id_GrazPhy,GrazPhy)
_SET_DIAGNOSTIC_(self%id_GrazPOM,GrazPOM)
_SET_DIAGNOSTIC_(self%id_MortPhy,MortPhy)
_SET_DIAGNOSTIC_(self%id_RespPhy,RespPhy)
_SET_DIAGNOSTIC_(self%id_ExcrPhy,ExcrPhy)
_SET_DIAGNOSTIC_(self%id_LimN,LimN)
_SET_DIAGNOSTIC_(self%id_LimT,LimT)
_SET_DIAGNOSTIC_(self%id_LimLight,LimLight)
_SET_DIAGNOSTIC_(self%id_GrowthPhy,GrowthPhy)
_SET_DIAGNOSTIC_(self%id_DOM_decay_ox,DOM_decay_ox)
_SET_DIAGNOSTIC_(self%id_DOM_decay_denitr,DOM_decay_denitr)
_SET_DIAGNOSTIC_(self%id_POM_decay_ox,POM_decay_ox)
_SET_DIAGNOSTIC_(self%id_POM_decay_denitr,POM_decay_denitr)
_SET_DIAGNOSTIC_(self%id_Autolys,Autolys)

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
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
!  Oxygen air-water flux. Adopted from ERSEM

! !INPUT PARAMETERS:
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

! !LOCAL VARIABLES:
   real(rk)                   :: O2, temp, salt, windspeed, aice
   real(rk)                   :: ko2o, TempT, Osat, Sc
   real(rk):: Q_O2, Qs, Qb, Qi ! gas fluxes: total, surface, large bubbles, small bubbles
   real(rk):: u_fric, Cd, kp, kc, delta_P ! bubbles related
   
   _HORIZONTAL_LOOP_BEGIN_
    _GET_(self%id_oxy,O2)
    _GET_(self%id_temp,temp)              ! temperature
    _GET_(self%id_salt,salt)              ! salinity
    _GET_HORIZONTAL_(self%id_windspeed,windspeed)
    if (self%use_aice) then
      _GET_HORIZONTAL_(self%id_aice,aice)
    end if
    
!The total flux (Qnet) as sum of flux through the water surface (Qs),  large bubbles surface (Qb) 
!      and from collapsing small bubbles (Qi). (McNeil, D Assaro 2007) Qnet = Qs + Qb + Qi
    
!Qs followintg NERSEM
!Sc = 1953.4_rk-128._rk*temp+3.9918_rk*temp*temp-0.050091_rk*temp*temp*temp !(Wanninkoff, 1992) old Sc
! formulation for the Schmidt number in NERSEM for O2 following Wanninkhof 2014
    Sc = 1920.4_rk - 135.6_rk*temp + 5.2122_rk*temp**2._rk - 0.10939_rk*temp**3 + 0.00093777_rk*temp**4._rk   

! Calculation of O2 saturation Osat according to UNESCO, 1986
   TempT = (temp+273.15_rk)/100._rk
   Osat = exp(-173.4292_rk+249.6339_rk/TempT+143.3483_rk*log(TempT)-21.8492_rk*TempT &
             +salt*(-0.033096_rk+0.014259_rk*TempT-0.0017_rk*TempT*TempT)) !Osat
   Osat = Osat*1000._rk/22.4_rk  ! convert from ml/l into uM

   ko2o = 0.251_rk*windspeed**2.0_rk*(sc/660._rk)**(-0.5_rk) ! Wanninkhof 2014

   Qs = ko2o*(Osat-O2)*0.24/86400._rk ! 0.24 is to convert from [cm/h] to [m/day]
!   Qs = windspeed*(Osat-O2)/86400._rk !After (Burchard et al., 2005)
   
 ! ! Qb and Qi following (Emerson, Bushinski, 2016)
          if(windspeed.le.11._rk) Cd=0.0012  ! drag coefficient 
          if(windspeed.gt.11._rk.and.windspeed.lt.20._rk) Cd=(0.49_rk+0.065_rk*windspeed)*0.001_rk
          if(windspeed.ge.20._rk) Cd=0.0018
          u_fric= 0.034*windspeed*sqrt(Cd) ! friction velocity
!          kp= 5.5_rk*(u_fric**2.76_rk)*((Sc/660._rk)**(-0.67_rk)) ! mass transfer coeff. for large bubbles
          kp= 5.5_rk*(u_fric**2.76_rk)*((Sc/660._rk)**(-0.67_rk)) ! mass transfer coeff. for large bubbles
          kc= 5.56_rk*(u_fric**3.86_rk) ! mass transfer coeff. for small bubbles
          delta_P= 1.52_rk*(u_fric**1.06_rk) !fractional increase in pressure in larg bubbles
      Qb= kp*((1._rk+delta_P)*Osat-O2)  !mmol/(m2xs)
      Qi= kc*Osat      !mmol/(m2xs)
   
    if (self%l_bubble.lt.0.5_rk.or.windspeed.lt.10._rk) then 
        Qb=0.0_rk   
    endif     
    if (self%s_bubble.lt.0.5_rk.or.windspeed.lt.10._rk) then 
        Qi=0.0_rk  
    endif        
   
   Q_O2 = Qs+Qb+Qi !mmol/(m2xs)
   if (self%use_aice) then
       Q_O2 = (1.0_rk-aice)*Q_O2 !Limit flux to area fraction not covered by ice
   end if

  _SET_SURFACE_EXCHANGE_(self%id_oxy,Q_O2)

_HORIZONTAL_LOOP_END_

   end subroutine

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
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: oxy,pom,dom,phy,het,nut
   real(rk) :: transmodel

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_nut,nut)
   _GET_(self%id_phy,phy)
   _GET_(self%id_het,het)
   _GET_(self%id_pom,pom)
   _GET_(self%id_dom,dom)
   _GET_(self%id_oxy,oxy)

   ! BURYING into the sediments, mmol/m2/s (sinking rates "Wxxx" are in m/s and positive upward)
   _SET_BOTTOM_EXCHANGE_(self%id_pom,self%Bu*self%Wpom*pom)
   _SET_BOTTOM_EXCHANGE_(self%id_phy,self%Bu*self%Wphy*phy)
   _SET_BOTTOM_EXCHANGE_(self%id_het,self%Bu*self%Whet*het)

   ! we use here the relaxation condition with relaxation time Trel
    if (self%transmodel.ge.0) then
      
       ! DOWNWARD flux of OXY dependent on redox conditions:
       !---  in suboxic conditions oxy is additionally consumed due to its depletion in the sediments
       _SET_BOTTOM_EXCHANGE_(self%id_oxy,(F_ox(oxy,self%O2_suboxic)*(self%b_ox-oxy)+F_subox(oxy,self%O2_suboxic)*(0._rk-oxy))/self%Trel)
       
       ! CHANGEABLE flux of NUT (NO3+NO2) dependent on redox conditions:
       !--- in oxic conditions nut is released from the oxygenated sediments
       !--- in suboxic conditions nut is consumed by the sediments due to denitrification in the sediments
       _SET_BOTTOM_EXCHANGE_(self%id_nut,(F_ox(oxy,self%O2_suboxic)*(self%b_nut-nut)+F_subox(oxy,self%O2_suboxic)*(0._rk-nut))/self%Trel)

       ! UPWARD flux of DOM (DON+NH4) dependent on redox conditions:
       !---  in suboxic conditions dom is additionally released from the sediments (doubles)
       _SET_BOTTOM_EXCHANGE_(self%id_dom,(F_ox(oxy,self%O2_suboxic)*(self%b_dom_ox-dom)+F_subox(oxy,self%O2_suboxic)*(2._rk*self%b_dom_ox-dom))/self%Trel)

    endif
   _HORIZONTAL_LOOP_END_

   end subroutine

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: A dependence for a process possible at high concntrations 
!
! !INTERFACE:
   real(rk) function F_ox(conc,threshold)
!
! !DESCRIPTION:
! A dependence for processes possible at higher than "threhold" concntrations 
!   
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: conc,threshold
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev, Anfisa Berezina
!
!EOP
!
!BOC
   F_ox=(0.5_rk+0.5_rk*tanh(conc-threshold))
   
   end function F_ox
!EOC
!-----------------------------------------------------------------------
   
!-----------------------------------------------------------------------
!BOP
! !IROUTINE: A dependence for a process possible at low concntrations 
!
! !INTERFACE:
   real(rk) function F_subox(conc,threshold)
!
! !DESCRIPTION:
! A dependence for processes possible at lower than "threhold" concntrations 
!   
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: conc,threshold
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev, Anfisa Berezina
!
!EOP
!
!BOC
   F_subox=(0.5_rk-0.5_rk*tanh(conc-threshold))
   
   end function F_subox
!EOC

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
!
!BOC
   yy=x**2/(a**2+x**2)

   end function yy
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Phy photosynthesis temperature limitation
  
  real(rk) function f_t(self,temperature)
    class (type_niva_oxydep),  intent(in):: self
    real(rk)                ,  intent(in):: temperature

    if (self%phy_t_dependence == 1) then
      ! ERGOM
      !   bm = 0.12_rk
      !   cm = 1.4_rk
      ! f_t = exp(bm*temperature-cm)
      f_t = exp(0.12_rk*temperature-1.4_rk)
    else if (self%phy_t_dependence == 2) then
      !for Arctic (Moore et al.,2002; Jin et al.,2008)
      !   t_0           = 0._rk      !reference temperature
      !   temp_aug_rate = 0.0663_rk  !temperature augmentation rate
      ! f_t = exp(temp_aug_rate*(temperature-t_0))
      f_t = exp(0.0663_rk*(temperature-0._rk))
    else if (self%phy_t_dependence == 3) then
      !ERSEM
      !   q10       = 2.0_rk    !Coefficient for uptake rate dependence on t
      !   t_upt_min = 10.0_rk   !Low  t limit for uptake rate dependence on t
      !   t_upt_max = 32.0_rk   !High t limit for uptake rate dependence on t
      ! f_t = q10**((temperature-t_upt_min)/10._rk)-&
      !       q10**((temperature-t_upt_max)/3._rk)
      f_t = max(0.0_rk, 2.0_rk**((temperature-10.0_rk)/10._rk) &
                      - 2.0_rk**((temperature-32.0_rk)/3._rk))
    end if
!   Some others:
 !      LimT     = self%q10**((t-self%t_upt_min)/10.) - self%q10**((t-self%t_upt_max)/3.) !Dependence on Temperature (ERSEM)
 !      LimT     = 0.5(1+tanh((t-tmin)/smin)) (1-0.5(1+th((t-tmax)/smax))) !Smin= 15  Smax= 15  Tmin=  10 Tmax= 35 !Dependence on Temperature   (Deb et al., .09)
 !      LimT     = exp(self%bm*temp-self%cm))        !Dependence on Temperature (used in (Ya,So, 2011) for Arctic)  
 !      LimT     = 1./(1.+exp(10.-temp))             !Dependence on Temperature (ERGOM for cya)
 !      LimT     = 1.-temp*temp/(temp*temp +12.*12.) !Dependence on Temperature (ERGOM for dia)
 !      LimT     = 2.**((temp- 10.)/10.) -2**((temp-32.)/3.) !(ERSEM)
 !      LimT     =q10*(T-20)/10 !Q10=1.88 (Gr., 2000)       
  end function f_t
   
end module fabm_niva_oxydep
