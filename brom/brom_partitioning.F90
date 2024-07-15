#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_partitioning
!
! !DESCRIPTION:
!   This module equilibrates total dissolved Substance 
!   between free dissolved form (Subst_dis) and 
!   Subastance partitioned with living orgamisms (Subst_biota), 
!   particulate (Subst_POM) and dissolved (Subst_DOM) organic
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
   type,extends(type_base_model),public :: type_niva_brom_partitioning
!     Variable identifiers
      type (type_state_variable_id)        :: id_Subst_dis, id_Subst_biota, id_Subst_POM, id_Subst_DOM, id_Subst_miner! , id_Subst_tot 
      type (type_state_variable_id)        :: id_Phy, id_Het, id_Baae, id_Bhae, id_Baan, id_Bhan, id_NH4, id_Sipart
      type (type_state_variable_id)        :: id_Mn4, id_FeS, id_FeS2, id_POML, id_DOML, id_POMR, id_DOMR
      type (type_dependency_id)            :: id_temp, id_par
      type (type_diagnostic_variable_id)   :: id_Subst_tot_diss, id_Subst_tot_part        

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
! !IROUTINE: Initialise the BROM equilibrium constant model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_partitioning), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev
!
!EOP
!-----------------------------------------------------------------------
!BOC

!    call self%register_state_variable(self%id_Subst_tot, 'Subst_tot', 'mol/m**3','Subst_tot', minimum=0.0_rk)

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
   call self%register_state_dependency(self%id_Subst_dis,  'Subst_dis', 'mol/m**3', 'Subst_dis') !,'Subst dissolved', minimum=0.0_rk)
   call self%register_state_dependency(self%id_Subst_biota, 'Subst_biota', 'mol/m**3','Subst_biota') !, minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_dependency(self%id_Subst_POM, 'Subst_POM', 'mol/m**3','Subst_POM') !, minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_dependency(self%id_Subst_DOM, 'Subst_DOM', 'mol/m**3','Subst_DOM') !, minimum=0.0_rk)
   call self%register_state_dependency(self%id_Subst_miner, 'Subst_miner', 'mol/m**3','Subst_miner') !, minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

    !Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_Subst_tot_diss,'Subst_tot_diss','mol/m**3','Subst_tot_diss')
   call self%register_diagnostic_variable(self%id_Subst_tot_part,'Subst_tot_part','mol/m**3','Subst_tot_part')   

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
   class (type_niva_brom_partitioning),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: Subst_dis, Subst_biota, Subst_POM, Subst_DOM, Subst_miner ! , Subst_tot 
   real(rk) :: Subst_tot_diss,Subst_tot_part
   real(rk) :: Phy, Het, Baae, Bhae, Baan, Bhan, POMR, DOMR, POML, DOML
   real(rk) :: dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, dSubst_miner !, dSubst_tot
   real(rk) :: dSubst_tot_diss, dSubst_tot_part
   real(rk) :: pol_bio  ! pollutant in BIOTA,"ng?"/l
   real(rk) :: pol_dom  ! pollutant in DOM, "ng?"/l
   real(rk) :: pol_pom  ! pollutant in POM, "ng?"/l
   
   real(rk) :: Subst_total   ! total pollutant, "ng?"/l
   real(rk) :: Subst_free     ! pollutant in dissolved INORGANIC,"ng?"/l, i.e. not partitioned  
   
   real(rk) :: sha_bio ! % share of polutant in living organisms
   real(rk) :: sha_pom ! % share of polutant in POM
   real(rk) :: sha_dom ! % share of polutant in DOM
   real(rk) :: sha_miner ! % share of polutant in clay (Sipart)  
   real(rk) :: sha_free ! % share of 'free' polutant 
   real(rk) :: uMn2lip=0.0084 !coeff.to transfer POM (umol/l N)->(g DryWeight/l) 
   !real(rk) :: rho_FeS= 5.90E7 !    # Density of FeS [mmolFe/m3] (default = 5.90E7 mmolFe/m3)
   !real(rk) :: rho_FeS2=4.17E7 !    # Density of FeS2 [mmolFe/m3] (default = 4.17E7 mmolFe/m3)
   !real(rk) :: rho_Mn4= 5.78E7 !    # Density of Mn4 [mmolMn/m3] (default = 5.78E7 mmolMn/m3)   
!====================================================

!  Ni
  real(rk) :: Kow_bio = 39811.  ! part.coeff. BIO/water (Allisson, 2005, Table 1, in L/kg)
  real(rk) :: Kow_pom = 39811.  ! part.coeff. POM/water (Allisson, 2005, Table 1, in L/kg)
  real(rk) :: Kow_dom = 125893. ! part.coeff. DOM/water (Allisson, 2005, Table 1, in L/kg)
! /Ni   
!====================================================
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment

   _GET_(self%id_Subst_dis,Subst_dis)   
   _GET_(self%id_Subst_biota,Subst_biota)   
   _GET_(self%id_Subst_POM,Subst_POM)   
   _GET_(self%id_Subst_DOM,Subst_DOM)   
   _GET_(self%id_Subst_miner,Subst_miner)   
 !  _GET_(self%id_Subst_tot,Subst_tot)
   
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
   
! Let's first assume that all the polutant is dissolved INORGANIC...
        Subst_total = Subst_dis+ Subst_biota+ Subst_POM+ Subst_DOM + Subst_miner  !total amount of pollutant 

! We assume that density of organic matter is the same as that of 
!  water, i.e. 1 g=1 ml and operate with weight units to caluclate 
!  the shares of pollutant partitioning medias:
       if((Phy+Het+Baae+Bhae+Baan+Bhan)<=0.) then 
        sha_bio=0. 
       else
        sha_bio=uMn2lip/1000.*(Phy+Het+Baae+Bhae+Baan+Bhan)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif 
       
       if((POML+POMR)<=0.) then 
        sha_pom=0. 
       else
        sha_pom=uMn2lip/1000.*(POML+POMR)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif     
       
       if((DOML+DOMR)<=0.) then 
        sha_dom=0. 
       else
        sha_dom=uMn2lip/1000.*(DOML+DOMR)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif    

      sha_free = 1.-sha_bio-sha_pom-sha_dom ! i.e Volume(weight in [kg]) of 1l of water minus volumes of org. and part. forms 
!! The free subst.conc. left as dissolved free
      Subst_free = Subst_total* sha_free /(sha_free +Kow_bio*sha_bio              &
     &            +Kow_pom*sha_pom +Kow_dom*sha_dom) ! Kow_water=1. needed for correcn units is excluded
!! subst.conc. partitioned to biota
     pol_bio=max(0.,Kow_bio*Subst_free*sha_bio/sha_free)
!! subst.conc. partitioning to POM
     pol_pom=max(0.,Kow_pom*Subst_free*sha_pom/sha_free)
!! subst.conc. partitioning to DOM
     pol_dom=max(0.,Kow_dom*Subst_free*sha_dom/sha_free)

! difference betweeen new and old state variable, needed for FABM
    dSubst_dis   = -Subst_dis   +Subst_free
    dSubst_biota = -Subst_biota +pol_bio
    dSubst_POM   = -Subst_POM   +pol_pom
    dSubst_DOM   = -Subst_DOM   +pol_dom
    dSubst_miner = 0.0_rk

    dSubst_tot_diss  = dSubst_dis+dSubst_DOM
    dSubst_tot_part  = dSubst_POM+dSubst_miner
    
   _SET_ODE_(self%id_Subst_dis,  dSubst_dis)
   _SET_ODE_(self%id_Subst_biota, dSubst_biota)
   _SET_ODE_(self%id_Subst_POM, dSubst_POM)
   _SET_ODE_(self%id_Subst_DOM, dSubst_DOM)
   _SET_ODE_(self%id_Subst_miner, dSubst_miner)

   _SET_DIAGNOSTIC_(self%id_Subst_tot_diss, dSubst_tot_diss)   
   _SET_DIAGNOSTIC_(self%id_Subst_tot_part, dSubst_tot_part)      
   
! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

end module