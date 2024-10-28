#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_ni
!
! !DESCRIPTION:
! This module described transformation of Ni connceted with 
! formation/dissolution of minerals and adsorption on solids.
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY: version 25.03.2017
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_ni

!     Variable identifiers

! variables defined in other modules
   type (type_state_variable_id)        :: id_H2S, id_O2, id_Mn4, id_Fe3, id_FeS, id_FeS2
   type (type_state_variable_id)        :: id_Baae,id_Bhae,id_Baan,id_Bhan
   type (type_state_variable_id)        :: id_Phy, id_Het, id_POMR,id_DOMR, id_SO4
! diagnostic dependences defined in other modules
   type (type_dependency_id):: id_Hplus
! variables defined in this modules
   type (type_state_variable_id)        :: id_Ni, id_NiS, id_Ni_biota, id_Ni_POM, id_Ni_DOM
   type (type_state_variable_id)        :: id_Ni_Mn4, id_Ni_Fe3,id_Ni_FeS, id_Ni_FeS2
! diagnostic dependences defined in this modules
   type (type_diagnostic_variable_id)   :: id_Ni_tot, id_Ni_tot_diss
   type (type_diagnostic_variable_id)   :: id_NiS_diss, id_NiS_form, id_NiS_ox
   type (type_diagnostic_variable_id)   :: id_ni_mn4_compl, id_ni_fe3_compl,id_ni_fes_compl
   type (type_diagnostic_variable_id)   :: id_ni_fes2_compl, id_ni_dom_compl, id_ni_pom_compl, id_ni_bio_compl
   type (type_diagnostic_variable_id)   :: id_Kad_Mn4, id_Kad_Fe3, id_Kad_FeS, id_Kad_FeS2, id_Kad_DOM, id_Kad_POM, id_Kad_bio
! Hg transformation coefficients
   real(rk) ::  K_NiS, K_NiS_form, K_NiS_diss, K_NiS_ox
   real(rk) ::  KNi_Mn4, Sad_Mn4, KNi_Fe3, Sad_Fe3, KNi_FeS, Sad_FeS
   real(rk) ::  KNi_FeS2, Sad_FeS2, KNi_DOM, Sad_DOM, KNi_POM, Sad_POM, KNi_bio, Sad_bio
   real(rk) ::  r_fes_ni, r_fes2_ni, r_mn4_ni, r_fe3_ni 
   real(rk) ::  Kow_bio_Ni, Kow_pom_Ni, Kow_dom_Ni, K_relax
   real(rk) ::  Wsed, Wphy, Wm

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
! This module described transformation of Ni connceted with 
! formation/dissolution of minerals and adsorption on solids.
!
! !INPUT PARAMETERS:
   class (type_niva_brom_ni), intent(inout), target :: self
   integer,                     intent(in)          :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev 
!
!EOP
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------   
! Parameters, i.e. rate constants  
   call self%get_parameter(self%K_NiS,      'K_NiS',    '[uM]',     'K_NiS equilibrium constant (Solubility Product Constant)', default=2510.0_rk)
   call self%get_parameter(self%K_NiS_form, 'K_NiS_form', '[1/day]','Specific rate of precipitation of NiS from Ni with H2S',   default=5.e-5_rk)
   call self%get_parameter(self%K_NiS_diss, 'K_NiS_diss', '[1/day]','Specific rate of dissollution of NiS to Ni and H2S',       default=1.e-6_rk)
   call self%get_parameter(self%K_NiS_ox,   'K_NiS_ox',   '[1/day]','Specific rate of oxidation of NiS',               default=0.0_rk) 
   call self%get_parameter(self%KNi_Mn4,  'KNi_Mn4', ' [-]',  'partitioning coeff. for Ni on Mn4',   default=100000.0_rk)
   call self%get_parameter(self%KNi_Fe3,  'KNi_Fe3', ' [-]',  'partitioning coeff. for Ni on Fe3',   default=100000.0_rk)
   call self%get_parameter(self%KNi_FeS,  'KNi_FeS',  '[-]',  'partitioning coeff. for Ni on FeS',   default=100000.0_rk)
   call self%get_parameter(self%KNi_FeS2, 'KNi_FeS2', '[-]',  'partitioning coeff. for Ni on FeS2',   default=100000.0_rk)
   call self%get_parameter(self%KNi_DOM, 'KNi_DOM', '[-]',  'partitioning coeff. for Ni on DOM',   default=100000.0_rk)
   call self%get_parameter(self%KNi_POM, 'KNi_POM', '[-]',  'partitioning coeff. for Ni on POM',   default=100000.0_rk)
   call self%get_parameter(self%KNi_bio, 'KNi_bio', '[-]',  'partitioning coeff. for Ni on bio',   default=100000.0_rk)
   call self%get_parameter(self%Sad_Mn4,   'Sad_Mn4',  '[-]',   'adsorbtion sites on Mn4',                default=0.01_rk)
   call self%get_parameter(self%Sad_Fe3,   'Sad_Fe3',  '[-]',   'adsorbtion sites on Fe3',                default=0.01_rk)
   call self%get_parameter(self%Sad_FeS,   'Sad_FeS',  '[-]',   'adsorbtion sites on FeS',                default=0.01_rk)
   call self%get_parameter(self%Sad_FeS2,  'Sad_FeS2', '[-]',   'adsorbtion sites on FeS2',               default=0.01_rk)
   call self%get_parameter(self%Sad_DOM,  'Sad_DOM', '[-]',   'adsorbtion sites on DOM',               default=0.01_rk)
   call self%get_parameter(self%Sad_POM,  'Sad_POM', '[-]',   'adsorbtion sites on POM',               default=0.01_rk)
   call self%get_parameter(self%Sad_bio,  'Sad_bio', '[-]',   'adsorbtion sites on bio',               default=0.01_rk)
   call self%get_parameter(self%Kow_bio_Ni, 'Kow_bio_Ni','[-]','partitioning koeff. for biota for Ni', default=100000.0_rk)
   call self%get_parameter(self%Kow_POM_Ni, 'Kow_POM_Ni','[-]','partitioning koeff. for POM for Ni',   default=100000.0_rk)
   call self%get_parameter(self%Kow_DOM_Ni, 'Kow_DOM_Ni','[-]','partitioning koeff. for DOM for Ni',   default=100000.0_rk)
   call self%get_parameter(self%K_relax,'K_relax','[-]','relaxation koeff.  for partitioning and adsorption', default=0.01_rk)

   call self%get_parameter(self%Wsed, 'Wsed', '[m/day]',  'Rate of sinking of detritus (POP, POM)',       default=5.00_rk)     
   call self%get_parameter(self%Wphy, 'Wphy', '[m/day]',  'Rate of sinking of Phy',                       default=0.10_rk)
   call self%get_parameter(self%Wm,   'Wm','   [m/day]',  'Rate of accelerated sinking of metals',        default=7.0_rk)
 
   !---- Ni---------!      
   call self%register_state_variable(self%id_Ni,      'Ni',    'mmol/m**3','Ni',         minimum=0.0_rk)
   call self%register_state_variable(self%id_Ni_biota,'Ni_biota','mmol/m**3','Ni_biota', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_Ni_POM,  'Ni_POM','mmol/m**3','Ni_POM', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Ni_DOM,  'Ni_DOM','mmol/m**3','Ni_DOM', minimum=0.0_rk)
   call self%register_state_variable(self%id_NiS,     'NiS',  'mmol/m**3','NiS',       minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_Ni_Mn4,  'Ni_Mn4','mmol/m**3','Ni_Mn4', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Ni_Fe3,  'Ni_Fe3','mmol/m**3','Ni_Fe3', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Ni_FeS,  'Ni_FeS','mmol/m**3','Ni_FeS', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Ni_FeS2, 'Ni_FeS2','mmol/m**3','Ni_FeS2', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)

   call self%register_state_dependency(self%id_Mn4, 'Mn4',   'mmol/m**3', 'Mn4')
   call self%register_state_dependency(self%id_Fe3, 'Fe3',   'mmol/m**3', 'Fe3')
   call self%register_state_dependency(self%id_FeS, 'FeS',   'mmol/m**3', 'FeS')
   call self%register_state_dependency(self%id_FeS2,'FeS2',  'mmol/m**3', 'FeS2')
   call self%register_state_dependency(self%id_H2S, 'H2S',   'mmol/m**3', 'H2S')
   call self%register_state_dependency(self%id_O2,  'O2',    'mmol/m**3', 'O2')
   call self%register_state_dependency(self%id_SO4, 'SO4',   'mmol/m**3',  'SO4')
   call self%register_state_dependency(self%id_Phy, 'Phy', 'mmol/m**3','Phy')
   call self%register_state_dependency(self%id_Het, 'Het', 'mmol/m**3','Het')
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_POMR,'POMR','mmol/m**3','POM refr')
   call self%register_state_dependency(self%id_DOMR,'DOMR','mmol/m**3','DOM refr')

    !call self%register_diagnostic_variable(self%id_ni_mn4_compl,'ni_mn4_compl','mmol/m**3/d',&
    !     'Ni adsorption on Mn4',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_ni_fes_compl,'ni_fes_compl','mmol/m**3/d',&
    !     'Ni adsorpion on FeS',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_ni_fes2_compl,'ni_fes2_compl','mmol/m**3/d',&
    !     'Ni adsorpion on FeS2',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(self%id_Ni_free,'Ni_free','mmol/m**3',&
    !     'Ni2+ free dissolved',output=output_time_step_integrated) ! dissolved inorganic
    call self%register_diagnostic_variable(self%id_Ni_tot_diss,'Ni_tot_diss','mmol/m**3',&
         'Ni total dissolved',output=output_time_step_integrated) ! dissolved total 
    call self%register_diagnostic_variable(self%id_Ni_tot,'Ni_total','mmol/m**3',&
         'Ni total',output=output_time_step_integrated) ! dissolved and particulate 
    call self%register_diagnostic_variable(self%id_NiS_form,'NiS_form','mmol/m**3',&
         'NiS formation rate',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_NiS_diss,'NiS_diss','mmol/m**3/d',&
         'NiS dissolution rate',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_NiS_ox,'NiS_ox','mmol/m**3',&
         'NiS oxidation rate',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Kad_Mn4,'Kad_Mn4','-',&
         'Kad_Mn4', output=output_time_step_integrated) !  sorption coeff. on Mn4
    call self%register_diagnostic_variable(self%id_Kad_Fe3,'Kad_Fe3','-',&
         'Kad_Fe3', output=output_time_step_integrated) !  sorption coeff. on Fe3
    call self%register_diagnostic_variable(self%id_Kad_FeS,'Kad_FeS','-',&
         'Kad_FeS', output=output_time_step_integrated) !  sorption coeff. on FeS
    call self%register_diagnostic_variable(self%id_Kad_FeS2,'Kad_FeS2','-',&
         'Kad_FeS2', output=output_time_step_integrated) !  sorption coeff. on FeS2
    call self%register_diagnostic_variable(self%id_Kad_DOM,'Kad_DOM','-',&
         'Kad_DOM', output=output_time_step_integrated) !  sorption coeff. on DOM
    call self%register_diagnostic_variable(self%id_Kad_POM,'Kad_POM','-',&
         'Kad_POM', output=output_time_step_integrated) !  sorption coeff. on POM
    call self%register_diagnostic_variable(self%id_Kad_bio,'Kad_bio','-',&
         'Kad_bio', output=output_time_step_integrated) !  sorption coeff. on biota
    call self%register_diagnostic_variable(self%id_ni_mn4_compl,'ni_mn4_compl','-',&
         'ni_mn4_compl', output=output_time_step_integrated) !  sorption on Mn4
    call self%register_diagnostic_variable(self%id_ni_fe3_compl,'ni_fe3_compl','-',&
         'ni_fe3_compl', output=output_time_step_integrated) !  sorption on Fe3
    call self%register_diagnostic_variable(self%id_ni_fes_compl,'ni_fes_compl','-',&
         'ni_fes_compl', output=output_time_step_integrated) !  sorption on FeS
    call self%register_diagnostic_variable(self%id_ni_fes2_compl,'ni_fes2_compl','-',&
         'ni_fes2_compl', output=output_time_step_integrated) !  sorption on FeS2
    call self%register_diagnostic_variable(self%id_ni_dom_compl,'ni_dom_compl','-',&
         'ni_dom_compl', output=output_time_step_integrated) !  sorption on DOM
    call self%register_diagnostic_variable(self%id_ni_pom_compl,'ni_pom_compl','-',&
         'ni_pom_compl', output=output_time_step_integrated) !  sorption on POM
    call self%register_diagnostic_variable(self%id_ni_bio_compl,'ni_bio_compl','-',&
         'ni_bio_compl', output=output_time_step_integrated) !  sorption on biota
    
 !Register diagnostic dependencies
    call self%register_dependency(self%id_Hplus,'Hplus', 'mmol/m**3','H+ Hydrogen')
    
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
   class (type_niva_brom_ni),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY: 
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev 
!
! !LOCAL VARIABLES:
   real(rk) ::  Ni, NiS, Ni_biota, Ni_POM, Ni_DOM
   real(rk) ::  Ni_Mn4, Ni_Fe3, Ni_FeS, Ni_FeS2
   real(rk) ::  ni_mn4_compl, ni_fe3_compl, ni_fes_compl, ni_fes2_compl, ni_dom_compl, ni_pom_compl, ni_bio_compl
   real(rk) ::  Ni_tot, Ni_tot_diss
   real(rk) ::  temp, O2, Mn4, Fe3, FeS, FeS2, depth, H2S, SO4 
   real(rk) ::  Phy, Het, Bhae, Baae, Bhan, Baan, POMR, DOMR 
   real(rk) ::  Om_NiS, NiS_form, NiS_diss, NiS_ox
   real(rk) ::  dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM
   real(rk) ::  dNi, dNiS, Kad_Mn4, Kad_Fe3, Kad_FeS, Kad_FeS2, Kad_DOM, Kad_POM, Kad_bio
  !diagnostic variables dependencies
   real(rk):: Hplus
   !EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

    ! Our own state variables
    _GET_(self%id_Ni_biota,Ni_biota)
    _GET_(self%id_Ni_POM,Ni_POM)
    _GET_(self%id_Ni_DOM,Ni_DOM)
    _GET_(self%id_Ni,Ni)
    _GET_(self%id_NiS,NiS)
    _GET_(self%id_Ni_Mn4,Ni_Mn4)
    _GET_(self%id_Ni_Fe3,Ni_Fe3)
    _GET_(self%id_Ni_FeS,Ni_FeS)
    _GET_(self%id_Ni_FeS2,Ni_FeS2)
   ! other modules state variables
    _GET_(self%id_H2S,H2S)
    _GET_(self%id_FeS,FeS)
    _GET_(self%id_FeS2,FeS2)
    _GET_(self%id_Fe3,Fe3)
    _GET_(self%id_Mn4,Mn4)
    _GET_(self%id_O2,O2)
    _GET_(self%id_Phy,Phy) 
    _GET_(self%id_Het,Het)
    _GET_(self%id_POMR,POMR)
    _GET_(self%id_DOMR,DOMR)
    _GET_(self%id_Baae,Baae)
    _GET_(self%id_Bhae,Bhae)
    _GET_(self%id_Baan,Baan)
    _GET_(self%id_Bhan,Bhan)
    _GET_(self%id_SO4,SO4)  
    !other modules diagnostics
    _GET_(self%id_Hplus,Hplus)
!-----------------------------------------------------------------
! NiS  formation/dissollution (REF1) is calculated for Ni_free excuding aborbed Ni 
      Om_NiS=H2S*Ni/(self%K_NiS)
!% FeS formation Fe2+ + HS- -> FeS + H+ (Bektursunova, 11)
    NiS_form=self%K_NiS_form*max(0._rk,(Om_NiS-1._rk)) !*1.e-20
!% FeS dissollution FeS + H+ -> Fe2+ + HS (Bektursunova, 11)
    NiS_diss=self%K_NiS_diss*NiS*max(0._rk,(1._rk-Om_NiS)) !*1.e-20
!% HgS oxydation  HgS + 2O2 -> Hg2+ + SO42-  ()
    NiS_ox=self%K_NiS_ox*NiS*O2/(O2+2._rk) 
!    NiS_ox=self%K_NiS_ox*NiS*thr_h(1.0_rk,O2,1._rk) 
!    if(O2.lt.0.5_rk) NiS_ox=0.0_rk

!!!  ! partitioning betweeen dissolved Hg(II) and OM
!!!    call partit (Ni, Ni_biota, Ni_POM, Ni_DOM, &
!!!                 Phy, Het, Baae, Bhae, Baan, Bhan, POMR, DOMR, &
!!!                 dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, &
!!!                 self%Kow_bio_Ni,self%Kow_POM_Ni,self%Kow_DOM_Ni)
!!!    
!!!    dSubst_dis=self%K_relax*dSubst_dis
!!!    dSubst_biota=self%K_relax*dSubst_biota
!!!    dSubst_POM=self%K_relax*dSubst_POM
!!!    dSubst_DOM=self%K_relax*dSubst_DOM
!!!    
!!!   _SET_ODE_(self%id_Ni_biota,dSubst_biota)
!!!   _SET_ODE_(self%id_Ni_POM,dSubst_POM)
!!!   _SET_ODE_(self%id_Ni_DOM,dSubst_DOM)
!!!!   _SET_ODE_(self%id_Ni_free,0.0)

!! Sorption of Ni on Mn oxides
       Kad_Mn4=self%KNi_Mn4*self%Sad_Mn4*Mn4/(Hplus*1000000._rk+self%KNi_Mn4*Ni)
     ni_mn4_compl = self%K_relax*Kad_Mn4*(Ni+Ni_Mn4)/(1.0_rk+Kad_Mn4)-Ni_Mn4
   _SET_ODE_(self%id_Ni_Mn4,ni_mn4_compl)

!! Sorption of Ni on Fe3
       Kad_Fe3=self%KNi_Fe3*self%Sad_Fe3*Fe3/(Hplus*1000000._rk+self%KNi_Fe3*Ni)
     ni_fe3_compl = self%K_relax*Kad_Fe3*(Ni+Ni_Fe3)/(1.0_rk+Kad_Fe3)-Ni_Fe3
   _SET_ODE_(self%id_Ni_Fe3,ni_fe3_compl)

!! Sorption of Ni on FeS
       Kad_FeS=self%KNi_FeS*self%Sad_FeS*FeS/(Hplus*1000000._rk+self%KNi_FeS*Ni)
     ni_fes_compl = self%K_relax*Kad_FeS*(Ni+Ni_FeS)/(1.0_rk+Kad_FeS)-Ni_FeS
   _SET_ODE_(self%id_Ni_FeS,ni_fes_compl)

!! Sorption of Ni on FeS2
       Kad_FeS2=self%KNi_FeS2*self%Sad_FeS2*FeS2/(Hplus*1000000._rk+self%KNi_FeS2*Ni)
     ni_fes2_compl = self%K_relax*Kad_FeS2*(Ni+Ni_FeS2)/(1.0_rk+Kad_FeS2)-Ni_FeS2
   _SET_ODE_(self%id_Ni_FeS2,ni_fes2_compl)
   
!! Sorption of Ni on DOMR
       Kad_DOM=self%KNi_DOM*self%Sad_DOM*DOMR/(Hplus*1000000._rk+self%KNi_DOM*Ni)
     ni_dom_compl = self%K_relax*Kad_DOM*(Ni+Ni_DOM)/(1.0_rk+Kad_DOM)-Ni_DOM
   _SET_ODE_(self%id_Ni_DOM,ni_dom_compl)

   !! Sorption of Ni on POMR
       Kad_POM=self%KNi_POM*self%Sad_POM*POMR/(Hplus*1000000._rk+self%KNi_POM*Ni)
     ni_pom_compl = self%K_relax*Kad_POM*(Ni+Ni_POM)/(1.0_rk+Kad_POM)-Ni_POM
   _SET_ODE_(self%id_Ni_POM,ni_pom_compl)
   
   !! Sorption of Ni on biota
       Kad_bio=self%KNi_bio*self%Sad_bio*(Phy+Het+Baae+Bhae+Baan+Bhan)/(Hplus*1000000._rk+self%KNi_bio*Ni)
     ni_bio_compl = self%K_relax*Kad_bio*(Ni+Ni_biota)/(1.0_rk+Kad_bio)-Ni_biota
   _SET_ODE_(self%id_Ni_biota,ni_bio_compl)   
   


   _SET_ODE_(self%id_Ni, NiS_diss-NiS_form+NiS_ox-ni_dom_compl-ni_pom_compl-ni_bio_compl-ni_fe3_compl-ni_fes_compl-ni_fes2_compl-ni_mn4_compl)  !+dSubst_dis  ! "Ni" includes Ni free and Ni adsorped
   _SET_ODE_(self%id_NiS,-NiS_diss+NiS_form-NiS_ox)
   _SET_ODE_(self%id_H2S,NiS_diss-NiS_form)
   _SET_ODE_(self%id_O2,-NiS_ox)
   _SET_ODE_(self%id_SO4, NiS_ox)

!      _SET_DIAGNOSTIC_(self%id_Ni_Mn4,Ni_Mn4)
!      _SET_DIAGNOSTIC_(self%id_Ni_FeS,Ni_FeS)
!      _SET_DIAGNOSTIC_(self%id_Ni_FeS2,Ni_FeS2)
!      _SET_DIAGNOSTIC_(self%id_Ni_free,Ni_free)
      _SET_DIAGNOSTIC_(self%id_Ni_tot,Ni+Ni_Mn4+Ni_Fe3+Ni_FeS+Ni_FeS2+Ni_DOM+Ni_POM+Ni_biota) !
      _SET_DIAGNOSTIC_(self%id_Ni_tot_diss,Ni+Ni_DOM)
      _SET_DIAGNOSTIC_(self%id_NiS_form,NiS_form)
      _SET_DIAGNOSTIC_(self%id_NiS_diss,NiS_diss)
      _SET_DIAGNOSTIC_(self%id_NiS_ox,NiS_ox)

! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
   subroutine partit (Subst_dis,Subst_biota, Subst_POM, Subst_DOM, &
                      Phy, Het, Baae, Bhae, Baan, Bhan, POMR, DOMR, &
                      dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, &
                      Kow_bio, Kow_pom, Kow_dom)
   ! !LOCAL VARIABLES:
   real(rk) :: Subst_dis, Subst_biota, Subst_POM, Subst_DOM,  Subst_tot 
   real(rk) :: Phy, Het, Baae, Bhae, Baan, Bhan, POMR, DOMR
   real(rk) :: dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, dSubst_tot
   real(rk) :: dSubst_tot_diss, dSubst_tot_part
   real(rk) :: pol_bio  ! pollutant in BIOTA,"ng?"/l
   real(rk) :: pol_dom  ! pollutant in DOM, "ng?"/l
   real(rk) :: pol_pom  ! pollutant in POM, "ng?"/l
   
   real(rk) :: Subst_total   ! total pollutant, "ng?"/l
   real(rk) :: Subst_free     ! pollutant in dissolved INORGANIC,"ng?"/l, i.e. not partitioned  
   
   real(rk) :: sha_bio ! % share of polutant in living organisms
   real(rk) :: sha_pom ! % share of polutant in POM
   real(rk) :: sha_dom ! % share of polutant in DOM
   real(rk) :: sha_free ! % share of 'free' polutant 
   real(rk) :: uMn2lip=0.0084 !coeff.to transfer POM (umol/l N)->(g DryWeight/l) 
   real(rk) :: rho_FeS= 5.90E7 !    # Density of FeS [mmolFe/m3] (default = 5.90E7 mmolFe/m3)
   real(rk) :: rho_FeS2=4.17E7 !    # Density of FeS2 [mmolFe/m3] (default = 4.17E7 mmolFe/m3)
   real(rk) :: rho_Mn4= 5.78E7 !    # Density of Mn4 [mmolMn/m3] (default = 5.78E7 mmolMn/m3)   
   real(rk) :: rho_Fe3= 5.90E7 !    # Density of Fe3 [mmolFe/m3] (default = 5.90E7 mmolFe/m3)

!====================================================

!  Hg
  real(rk) :: Kow_bio != 100000.  ! part.coeff. BIO/water (Allisson, 2005, Table 1, in L/kg)
  real(rk) :: Kow_pom != 100000.  ! part.coeff. POM/water (Allisson, 2005, Table 1, in L/kg)
  real(rk) :: Kow_dom != 100000. ! part.coeff. DOM/water  (Allisson, 2005, Table 1, in L/kg)
! /Hg   
!====================================================
!EOP 
!---
   
! Let's first assume that all the polutant is dissolved INORGANIC...
        Subst_total = Subst_dis+ Subst_biota+ Subst_POM+ Subst_DOM ! total amount of pollutant 

! We assume that density of organic matter is the same as that of 
!  water, i.e. 1 g=1 ml and operate with weight units to caluclate 
!  the shares of pollutant partitioning medias:
       if((Phy+Het+Baae+Bhae+Baan+Bhan)<=0.) then 
        sha_bio=0. 
       else
        sha_bio=uMn2lip/1000.*(Phy+Het+Baae+Bhae+Baan+Bhan)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif 
       
       if(POMR<=0.) then 
        sha_pom=0. 
       else
        sha_pom=uMn2lip/1000.*POMR  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif     
       
       if(DOMR<=0.) then 
        sha_dom=0. 
       else
        sha_dom=uMn2lip/1000.*DOMR  ! Volume(weight in kg, g->kg=/1000) of BIO
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

    end subroutine partit

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

end module