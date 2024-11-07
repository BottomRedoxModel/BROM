#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_partit
!
! !DESCRIPTION:
!   This module equilibrates total dissolved Substance 
!   between free dissolved form (Ci_free) and 
!   Subastance partitioned with living orgamisms (Ci_phy and Ci_het), 
!   particulate (Ci_POM) and dissolved (Ci_DOM) organic
!   matter, using partitioning approach. The rate of equilibration 
!   is assumed to be fast process compared with the typical 
!   model timestep
!
! !USES:PON

   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_partit
!     Variable identifiers
      type (type_state_variable_id)        :: id_Ci_free, id_Ci_phy, id_Ci_het, id_Ci_POM, id_Ci_DOM, id_Ci_miner! , id_Ci_tot 
      type (type_state_variable_id)        :: id_Phy, id_Het, id_Baae, id_Bhae, id_Baan, id_Bhan, id_NH4, id_Sipart, id_O2
      type (type_state_variable_id)        :: id_Mn4, id_FeS, id_FeS2, id_POML, id_DOML, id_POMR, id_DOMR
      type (type_dependency_id)            :: id_temp, id_par, id_depth
      type (type_dependency_id)            :: id_Hplus,  id_Wadd      
      type (type_diagnostic_variable_id)   :: id_Ci_tot_diss, id_Ci_tot_part        

      real(rk) :: Wsed, Wphy, Whet, Wm
      real(rk) :: Iopt, O2_suboxic
      real(rk) :: Kow_bio, Kow_pom, Kow_dom
      real(rk) :: K_biodegrad, K_biodegrad_anae, K_hydrolysis, K_photolysis
   contains
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement      
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
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_partit), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev
!
!EOP
!-----------------------------------------------------------------------
!BOC
! Parameters, i.e. rate constants    
   call self%get_parameter(self%K_biodegrad,      'K_biodegrad', '[1/day]','rate of biodegradation, 1/d', default=00000.0_rk)
   call self%get_parameter(self%K_biodegrad_anae, 'K_biodegrad_anae', '[1/day]','rate of biodegradation in anaerobic cond., 1/d', default=00000.0_rk)
   call self%get_parameter(self%K_hydrolysis, 'K_hydrolysis','[1/day]','rate of hydrolysis, 1/d',     default=00000.0_rk)
   call self%get_parameter(self%K_photolysis, 'K_photolysis','[1/day]','rate of photolysis, 1/d',     default=00000.0_rk)
   call self%get_parameter(self%Kow_bio, 'Kow_bio','[-]','partitioning coeff. BIO/water ', default=100000.0_rk)
   call self%get_parameter(self%Kow_pom, 'Kow_pom','[-]','partitioning coeff. POM/water ', default=100000.0_rk)
   call self%get_parameter(self%Kow_dom, 'Kow_dom','[-]','partitioning coeff. DOM/water ', default=100000.0_rk)
   call self%get_parameter(self%Wphy, 'Wphy', '[m/day]',  'Rate of sinking of Phy',                 default=0.10_rk)
   call self%get_parameter(self%Whet, 'Whet', '[m/day]',  'Rate of sinking of Het',                 default=5.00_rk) 
   call self%get_parameter(self%Wsed, 'Wsed', '[m/day]',  'Rate of sinking of detritus (POP, POM)', default=5.00_rk) 
   call self%get_parameter(self%Wm,   'Wm','   [m/day]',  'Rate of accelerated sinking of metals',  default=7.0_rk)
   call self%get_parameter(self%O2_suboxic, 'O2_suboxic', 'mmol/m3', 'Threshold O2 value for oxic/suboxic switch', default=40._rk)
   call self%get_parameter(self%Iopt,       'Iopt',       'Watts/m**2/h', 'Optimal irradiance',   default=25.0_rk)

   
    !Register state variables 
   call self%register_state_variable(self%id_Ci_free,  'Ci_free', 'mol/m**3', 'Ci_free', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ci_phy, 'Ci_phy', 'mol/m**3','Ci_phy', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ci_het, 'Ci_het', 'mol/m**3','Ci_het', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ci_POM, 'Ci_POM', 'mol/m**3','Ci_POM', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ci_DOM, 'Ci_DOM', 'mol/m**3','Ci_DOM', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ci_miner, 'Ci_miner', 'mol/m**3','Ci_miner', minimum=0.0_rk)

    !Register state dependencies 
   call self%register_state_dependency(self%id_Phy, 'Phy', 'mmol/m**3','Phy')
   call self%register_state_dependency(self%id_Het, 'Het', 'mmol/m**3','Het')
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','aerobic autotrophic bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','aerobic heterotrophic bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','anaerobic aurotrophic bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','anaerobic heterotrophic bacteria')
   call self%register_state_dependency(self%id_POML,'POML','mmol/m**3','labile POM')
   call self%register_state_dependency(self%id_DOML,'DOML','mmol/m**3','labile DOM')
   call self%register_state_dependency(self%id_POMR,'POMR','mmol/m**3','refractory POM')
   call self%register_state_dependency(self%id_DOMR,'DOMR','mmol/m**3','refractory DOM')
   call self%register_state_dependency(self%id_O2,  'O2',    'mmol/m**3', 'O2')

   call self%register_dependency(self%id_Hplus,'Hplus','mmol/m**3','H+ hydrogen')
   call self%register_dependency(self%id_Wadd,'Wadd','[1/day]',   'Additional sinking velocity via Mn4 adsorptoin')

    !Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_Ci_tot_diss,'Ci_tot_diss','mol/m**3','Ci_tot_diss',output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_Ci_tot_part,'Ci_tot_part','mol/m**3','Ci_tot_part',output=output_time_step_integrated)   

   ! Specify that are rates computed in this module are per day (default: per second)
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
   class (type_niva_brom_partit),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: Ci_free, Ci_phy, Ci_het, Ci_POM, Ci_DOM, Ci_miner ! , Ci_tot 
   real(rk) :: Ci_tot_diss,Ci_tot_part
   real(rk) :: Phy, Het, Baae, Bhae, Baan, Bhan, POMR, DOMR, POML, DOML, O2
   real(rk) :: dCi_free, dCi_phy, dCi_het, dCi_POM, dCi_DOM, dCi_miner 
   real(rk) :: temp, depth, Iz
   real(rk) :: Iopt, O2_suboxic
   real(rk) :: dCi_tot_diss, dCi_tot_part
   real(rk) :: pol_phy  ! "new" pollutant in phytoplankton and bacteria,"ng?"/l
   real(rk) :: pol_het     ! "new" pollutant in heterotrophs,"ng?"/l
   real(rk) :: pol_dom     ! "new" pollutant in DOM, "ng?"/l
   real(rk) :: pol_pom     ! "new" pollutant in POM, "ng?"/l
   real(rk) :: pol_free     ! "new" pollutant in dissolved INORGANIC,"ng?"/l, i.e. not partitioned  
   
   real(rk) :: Ci_total   ! total pollutant, "ng?"/l
   
   real(rk) :: sha_phy ! % share of polutant in phytoplankton and bacteria
   real(rk) :: sha_het    ! % share of polutant in heterotrophs
   real(rk) :: sha_pom    ! % share of polutant in POM
   real(rk) :: sha_dom    ! % share of polutant in DOM
   real(rk) :: sha_miner  ! % share of polutant in clay (Sipart)  
   real(rk) :: sha_free   ! % share of 'free' polutant 
   real(rk) :: sha_poten_partit ! % share of potential partitioning based on volumes and part.coeffs 
   real(rk) :: uMn2lip=0.0084   !coeff.to transfer POM (umol/l N)->(g DryWeight/l) 
   real(rk) :: decay_total      !total rate of decay of Ci (d-1)
!====================================================
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
  ! Environment
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_depth,depth)            ! depth
   _GET_(self%id_par,Iz)                 ! local photosynthetically active radiation

   _GET_(self%id_Ci_free,Ci_free)   
   _GET_(self%id_Ci_phy,Ci_phy)   
   _GET_(self%id_Ci_het,Ci_het)   
   _GET_(self%id_Ci_POM,Ci_POM)   
   _GET_(self%id_Ci_DOM,Ci_DOM)   
   _GET_(self%id_Ci_miner,Ci_miner)   
 !  _GET_(self%id_Ci_tot,Ci_tot)
   
   _GET_(self%id_Phy,Phy)    
   _GET_(self%id_Het,Het)    
   _GET_(self%id_Baae,Baae)
   _GET_(self%id_Bhae,Bhae)    
   _GET_(self%id_Baan,Baan)    
   _GET_(self%id_Bhan,Bhan)    
!   _GET_(self%id_NH4,NH4)    
   _GET_(self%id_POML,POML)    
   _GET_(self%id_DOML,DOML)
   _GET_(self%id_POMR,POMR)    
   _GET_(self%id_DOMR,DOMR) 
   _GET_(self%id_O2,O2)
   
! Let's first assume that all the polutant is dissolved INORGANIC...
        Ci_total = Ci_free+ Ci_phy+ Ci_het+ Ci_POM+ Ci_DOM + Ci_miner  !total amount of pollutant 

! We assume that density of organic matter is the same as that of 
!  water, i.e. 1 g=1 ml and operate with weight units to caluclate 
!  the shares of pollutant partitioning medias:
       if((Phy+Baae+Bhae+Baan+Bhan)<=0._rk) then 
        sha_phy=0._rk 
       else
        sha_phy=uMn2lip/1000._rk*(Phy+Baae+Bhae+Baan+Bhan)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif 

       if(Het<=0._rk) then 
        sha_het=0._rk 
       else
        sha_het=uMn2lip/1000._rk*Het  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif 
       
       if((POML+POMR)<=0._rk) then 
        sha_pom=0._rk 
       else
        sha_pom=uMn2lip/1000._rk*(POML+POMR)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif     
       
       if((DOML+DOMR)<=0._rk) then 
        sha_dom=0.0_rk 
       else
        sha_dom=uMn2lip/1000._rk*(DOML+DOMR)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif    

      sha_free = 1._rk-sha_phy-sha_het-sha_pom-sha_dom ! i.e Volume(weight in [kg]) of 1l of water minus volumes of org. and part. forms 
!! sum of potential shared due to volumes AND part.coeffs.  Kow_water=1. needed for correct units is excluded

!! Now we calculate "new" concentrations, "pol_XXX" from total Ci_total using calculated shares.
!!       The free Ci  left as dissolved free
   pol_free = Ci_total * sha_free / (sha_free + self%Kow_bio*sha_phy &
           + self%Kow_bio*sha_het + self%Kow_pom*sha_pom + self%Kow_dom*sha_dom)
!! subst.conc. partitioned to phy and bact
   pol_phy=max(0.0_rk,self%Kow_bio*pol_free*sha_phy/sha_free)
!! subst.conc. partitioned to het
   pol_het = max(0.0_rk,self%Kow_bio*pol_free*sha_het/sha_free)
!! subst.conc. partitioning to POM
   pol_pom = max(0.0_rk,self%Kow_pom*pol_free*sha_pom/sha_free)           
!! subst.conc. partitioning to DOM
   pol_dom = max(0.0_rk,self%Kow_dom*pol_free*sha_dom/sha_free)

! difference betweeen new and old state variable, needed for FABM 
    dCi_free    = -Ci_free    +pol_free
    dCi_phy = -Ci_phy +pol_phy
    dCi_het    = -Ci_het    +pol_het
    dCi_POM    = -Ci_POM    +pol_pom
    dCi_DOM    = -Ci_DOM    +pol_dom
    dCi_miner  =  0.0_rk

! Add to the calculated difference decrease of Ci_xxxxx due to biogegradaion, hydrolysis and photolysis:     
!    decay_total = ( &
!         self%K_biodegrad-(self%K_biodegrad-self%K_biodegrad_anae)*thr_l(15._rk,O2,1._rk) & ! biodegradation f(O2)
!       + self%K_hydrolysis &                          ! hydropysis
!       + self%K_photolysis*Iz/25.*exp(1._rk-Iz/25.))  !photolysis f(light)

    dCi_free = dCi_free - Ci_free *(self%K_biodegrad-(self%K_biodegrad-self%K_biodegrad_anae)*thr_l(self%O2_suboxic,O2,1._rk) &
                                       +self%K_photolysis*Iz/self%Iopt*exp(1._rk-Iz/self%Iopt) + self%K_hydrolysis) 
    dCi_phy   = dCi_phy - Ci_phy  * self%K_biodegrad-(self%K_biodegrad-self%K_biodegrad_anae)*thr_l(self%O2_suboxic,O2,1._rk)
    dCi_het   = dCi_het - Ci_het  * self%K_biodegrad-(self%K_biodegrad-self%K_biodegrad_anae)*thr_l(self%O2_suboxic,O2,1._rk)
    dCi_POM   = dCi_POM - Ci_POM  * self%K_biodegrad-(self%K_biodegrad-self%K_biodegrad_anae)*thr_l(self%O2_suboxic,O2,1._rk)
    dCi_DOM   = dCi_DOM - Ci_DOM  *(self%K_biodegrad-(self%K_biodegrad-self%K_biodegrad_anae)*thr_l(self%O2_suboxic,O2,1._rk) &
                                       +self%K_photolysis*Iz/self%Iopt*exp(1._rk-Iz/self%Iopt))                 
    dCi_miner =  0.0_rk

    dCi_tot_diss  =  dCi_free +dCi_DOM
    dCi_tot_part  =  dCi_phy +dCi_het +dCi_POM +dCi_miner
    
   _SET_ODE_(self%id_Ci_free,  dCi_free)
   _SET_ODE_(self%id_Ci_phy,   dCi_phy)
   _SET_ODE_(self%id_Ci_het,   dCi_het)
   _SET_ODE_(self%id_Ci_POM,   dCi_POM)
   _SET_ODE_(self%id_Ci_DOM,   dCi_DOM)
   _SET_ODE_(self%id_Ci_miner, dCi_miner)

   _SET_DIAGNOSTIC_(self%id_Ci_tot_diss, pol_free+pol_dom)   
   _SET_DIAGNOSTIC_(self%id_Ci_tot_part, pol_phy+pol_het+pol_pom)      
   
! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC
!-------------------------------------------------------------------------
real(rk) function thr_h(threshold_value,var_conc,koef)
   ! Threshold value for the reaction 
   ! koef 1 gives regular tgh function 
   ! 0.1 - smooth function 
   real(rk), intent(in) :: threshold_value,var_conc,koef
   thr_h = 0.5+0.5*tanh((var_conc-threshold_value)*koef)
end function    
!-------------------------------------------------------------------------
real(rk) function thr_l(threshold_value,var_conc,koef)
! Threshold value for the reaction 
real(rk), intent(in) :: threshold_value,var_conc,koef
thr_l = 0.5-0.5*tanh((var_conc-threshold_value)*koef)
end function 
!-------------------------------------------------------------------------
  ! Set increased manganese sinking via MnIV and MnIII oxides formation
   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_niva_brom_partit), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      
      real(rk) :: Wadd, Wphy_tot, Whet_tot, Wsed_tot, Wm_tot

      _LOOP_BEGIN_
   
      _GET_(self%id_Wadd,Wadd)
      
       Wphy_tot = self%Wphy + 0.25_rk * Wadd
       Whet_tot = self%Whet + 0.5_rk * Wadd
       Wsed_tot = self%Wsed + Wadd
       Wm_tot = self%Wm + Wadd
 
      _ADD_VERTICAL_VELOCITY_(self%id_Ci_phy, Wphy_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_Ci_het, Whet_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_Ci_POM, Wsed_tot)
      _ADD_VERTICAL_VELOCITY_(self%id_Ci_miner, Wsed_tot)
  

      _LOOP_END_
   
   end subroutine get_vertical_movement
end module