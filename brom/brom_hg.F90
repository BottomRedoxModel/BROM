#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_hg
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_hg
!     Variable identifiers
! variables defined in other modules
   type (type_state_variable_id)        :: id_H2S, id_O2, id_Mn4, id_Fe3, id_Baae,id_Bhae,id_Baan,id_Bhan
   type (type_state_variable_id)        :: id_Phy, id_Het, id_POML,id_POMR, id_DON, id_SO4
! global dependencies
   type (type_dependency_id)            :: id_temp,id_par,id_depth
! diagnostic dependences defined in other modules
   type (type_dependency_id) :: id_Hplus
! variables defined in this modules
   type (type_state_variable_id)        :: id_Hg0, id_Hg2,id_MeHg,id_HgS
   type (type_state_variable_id)        :: id_MeHg_biota, id_MeHg_POM, id_MeHg_DOM, id_MeHg_Fe3, id_MeHg_Mn4, id_MeHg_free
   type (type_state_variable_id)        :: id_Hg2_biota, id_Hg2_POM, id_Hg2_DOM, id_Hg2_Fe3, id_Hg2_Mn4, id_Hg2_free
! diagnostic dependences defined in this modules
   type (type_diagnostic_variable_id)   :: id_Hg2_tot, id_Hg2_tot_diss
   type (type_diagnostic_variable_id)   :: id_MeHg_tot, id_MeHg_tot_diss, id_MeHg_procent_diss
   type (type_diagnostic_variable_id)   :: id_Hg_tot, id_Hg_tot_diss, id_Kad_Hg2, id_Kad_MeHg, id_Om_HgS
   type (type_diagnostic_variable_id)   :: id_hg2_fe3_compl, id_hg2_mn4_compl, id_mehg_fe3_compl, id_mehg_mn4_compl
! Hg transformation coefficients
   real(rk) :: K_hg2_mehg , K_mehg_hg2, K_hg2_hg0, K_hg0_hg2
   real(rk) :: K_mehg_irr_degr ! photo-degradation
   real(rk) :: K_hg0_irr_ox ! photo-oxidation of hg0    
   real(rk) :: K_hg2_irr_red ! photo-reduction of Hg2      
   real(rk) :: K_HgS, K_hgs_form, K_hgs_ox, K_hgs_diss, K_mehg_h2s, O2s_nf      
   real(rk) :: KHg2_Fe3, Sad_Fe3, KMeHg_Fe3, r_fe3_mehg,  r_mn4_mehg, r_fe3_hg2, r_mn4_hg2, r_fe_n 
   real(rk) :: Kow_bio_Hg2, Kow_pom_Hg2, Kow_dom_Hg2, Kow_bio_MeHg, Kow_pom_MeHg, Kow_dom_MeHg,K_relax
   real(rk) :: Wsed, Wphy, Wm            

   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
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
! Module describes transformation of mercury (Hg)
!
! !INPUT PARAMETERS:
   class (type_niva_brom_hg), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev
!
!EOP
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------   
! Parameters, i.e. rate constants
   call self%get_parameter(self%K_HgS,          'K_HgS',          ' - ', 'HgS equilibrium constant',       default=0.0_rk)     
   call self%get_parameter(self%K_hg2_mehg,     'K_hg2_mehg',     ' - ', 'Coef. of mercury methilation',   default=0.0_rk)  
   call self%get_parameter(self%K_mehg_hg2,     'K_mehg_hg2',     ' - ', 'Coef. of demethylation of MeHg', default=0.0_rk) 
   call self%get_parameter(self%K_hg2_hg0,      'K_hg2_hg0',      ' - ', 'Coef. biotic reduction of Hg2',  default=0.0_rk)    
   call self%get_parameter(self%K_hg0_hg2,      'K_hg0_hg2',      ' - ', 'Coef. of dark oxidation of Hg0', default=0.0_rk)    
   call self%get_parameter(self%K_mehg_irr_degr,'K_mehg_irr_degr',' - ', 'photo-degradation of MeHg',      default=0.0_rk) 
   call self%get_parameter(self%K_hg0_irr_ox,   'K_hg0_irr_ox',   ' - ', 'photo-oxidation of hg0',         default=0.0_rk)    
   call self%get_parameter(self%K_hgs_form,     'K_hgs_form',     ' - ', 'Formation of HgS',               default=0.0_rk) 
   call self%get_parameter(self%K_hgs_ox,       'K_hgs_ox',       ' - ', 'Oxidation of HgS',               default=0.0_rk) 
   call self%get_parameter(self%K_hgs_diss,     'K_hgs_diss',     ' - ', 'Dissolution of HgS',             default=0.0_rk)  
   call self%get_parameter(self%K_hg2_irr_red,  'K_hg2_irr_red',  ' - ', 'photo-reduction of Hg2 ',        default=0.0_rk)    
   call self%get_parameter(self%K_mehg_h2s,     'K_mehg_h2s',     ' - ', 'Coef. of reduction of MeHg',     default=0.0_rk)

   call self%get_parameter(self%r_fe3_mehg,'r_fe3_mehg', '[-]', 'Fe3[uM]/MeHg[uM] partitioning coeff. for Fe3',  default=833.0_rk)
   call self%get_parameter(self%r_mn4_mehg,'r_mn4_mehg',  '[-]','MnO2[uM]/MeHg[uM] partitioning coeff. for MnO2',default=833.0_rk)
   call self%get_parameter(self%r_fe3_hg2, 'r_fe3_hg2', '[-]',  'Fe3[uM]/Hg2[uM] partitioning coeff. for Fe3',   default=833.0_rk)
   call self%get_parameter(self%r_mn4_hg2, 'r_mn4_hg2',  '[-]', 'MnO2[uM]/Hg2[uM] partitioning coeff. for MnO2', default=833.0_rk)
   call self%get_parameter(self%r_fe_n,    'r_fe_n',  '[-]',    'fe/n for OM',                            default=1.0_rk) 
   call self%get_parameter(self%KHg2_Fe3,  'KHg2_Fe3', '[-]',   'partitioning coeff. for  Hg2 on Fe3',    default=100000.0_rk)
   call self%get_parameter(self%KMeHg_Fe3, 'KMeHg_Fe3', '[-]',  'partitioning coeff. for  MeHg on Fe3',   default=100000.0_rk)
   call self%get_parameter(self%Sad_Fe3,   'Sad_Fe3', '[-]',    'adsorbtion sites on Fe3',                default=0.01_rk)
   call self%get_parameter(self%Kow_bio_Hg2, 'Kow_bio_Hg2','[-]','partitioning koeff. for biota for Hg2', default=100000.0_rk)
   call self%get_parameter(self%Kow_POM_Hg2, 'Kow_POM_Hg2','[-]','partitioning koeff. for POM for Hg2',   default=100000.0_rk)
   call self%get_parameter(self%Kow_DOM_Hg2, 'Kow_DOM_Hg2','[-]','partitioning koeff. for DOM for Hg2',   default=100000.0_rk)
   call self%get_parameter(self%Kow_bio_MeHg,'Kow_bio_MeHg','[-]','partitioning koeff. for biota for MeHg', default=100000.0_rk)
   call self%get_parameter(self%Kow_POM_MeHg,'Kow_POM_MeHg','[-]','partitioning koeff. for POM for MeHg', default=100000.0_rk)
   call self%get_parameter(self%Kow_DOM_MeHg,'Kow_DOM_MeHg','[-]','partitioning koeff. for DOM for MeHg', default=100000.0_rk)
   call self%get_parameter(self%K_relax,'K_relax','[-]','relaxation koeff.  for partitioning and adsorption', default=0.01_rk)

   call self%get_parameter(self%Wsed, 'Wsed', '[m/day]',  'Rate of sinking of detritus (POM)',       default=5.00_rk)     
   call self%get_parameter(self%Wphy, 'Wphy', '[m/day]',  'Rate of sinking of Phy',                       default=0.10_rk)
   call self%get_parameter(self%Wm,   'Wm','   [m/day]',  'Rate of accelerated sinking of metals',        default=7.0_rk)
   call self%get_parameter(self%O2s_nf, 'O2s_nf', '[uM O]','half saturation for nitrification',default=4.488_rk)

   call self%register_state_variable(self%id_Hg0, 'Hg0', 'mmol/m**3','Hg(0)', minimum=0.0_rk)
   call self%register_state_variable(self%id_Hg2, 'Hg2', 'mmol/m**3','Hg(II)',minimum=0.0_rk)
   call self%register_state_variable(self%id_MeHg,'MeHg','mmol/m**3','MeHg',  minimum=0.0_rk)
   call self%register_state_variable(self%id_HgS, 'HgS', 'mmol/m**3','HgS',   minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_MeHg_biota,'MeHg_biota','mmol/m**3','MeHg_biota', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_MeHg_POM,  'MeHg_POM','mmol/m**3','MeHg_POM', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_MeHg_DOM,  'MeHg_DOM','mmol/m**3','MeHg_DOM', minimum=0.0_rk)
   call self%register_state_variable(self%id_MeHg_Mn4,  'MeHg_Mn4','mmol/m**3','MeHg_Mn4', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_MeHg_Fe3,  'MeHg_Fe3','mmol/m**3','MeHg_Fe3', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_MeHg_free, 'MeHg_free','mmol/m**3','MeHg_free', minimum=0.0_rk)
   call self%register_state_variable(self%id_Hg2_biota,'Hg2_biota','mmol/m**3','Hg2_biota', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_Hg2_POM,  'Hg2_POM','mmol/m**3','Hg2_POM', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Hg2_DOM,  'Hg2_DOM','mmol/m**3','Hg2_DOM', minimum=0.0_rk)
   call self%register_state_variable(self%id_Hg2_Mn4,  'Hg2_Mn4','mmol/m**3','Hg2_Mn4', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Hg2_Fe3,  'Hg2_Fe3','mmol/m**3','Hg2_Fe3', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Hg2_free, 'Hg2_free','mmol/m**3','Hg2_free', minimum=0.0_rk)

   call self%register_state_dependency(self%id_O2,  'O2',    'mmol/m**3', 'O2')
   call self%register_state_dependency(self%id_Mn4, 'Mn4',   'mmol/m**3', 'Mn4')
   call self%register_state_dependency(self%id_Fe3, 'Fe3',   'mmol/m**3', 'Fe3')
   call self%register_state_dependency(self%id_H2S, 'H2S',   'mmol/m**3',  'H2S')
   call self%register_state_dependency(self%id_SO4, 'SO4',   'mmol/m**3',  'SO4')
   call self%register_state_dependency(self%id_Phy, 'Phy', 'mmol/m**3','Phy')
   call self%register_state_dependency(self%id_Het, 'Het', 'mmol/m**3','Het')
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_POML,'POML','mmol/m**3','particulate organic nitrogen')
   call self%register_state_dependency(self%id_POMR,'POMR','mmol/m**3','POM refractory')
   call self%register_state_dependency(self%id_DON,'DON','mmol/m**3','dissolved organic nitrogen')
   

    !call self%register_diagnostic_variable(&
    !     self%id_Hg2_free,'Hg2_free','mmol/m**3',&
    !     'Hg2+ free + part on biota',&
    !     output=output_time_step_integrated) ! dissolved inorganic
    call self%register_diagnostic_variable(&
         self%id_Hg2_tot_diss,'Hg2_tot_diss','mmol/m**3',&
         'Hg(II) total dissolved',&
         output=output_time_step_integrated) ! dissolved organic and inorganic
    call self%register_diagnostic_variable(&
         self%id_Hg2_tot,'Hg2_total','mmol/m**3',&
         'Hg(II) total',&
         output=output_time_step_integrated) ! dissolved and particulate
    call self%register_diagnostic_variable(&
         self%id_Om_HgS,'Om_HgS','-',&
         'Om_HgS saturation',&
         output=output_time_step_integrated) ! saturation of HgS

!    call self%register_diagnostic_variable(&
!         self%id_MeHg_free,'MeHg_free','mmol/m**3',&
!         'MeHg+ free + part on biota',&
!         output=output_time_step_integrated) ! dissolved inorganic
    call self%register_diagnostic_variable(&
         self%id_MeHg_tot_diss,'MeHg_tot_diss','mmol/m**3',&
         'MeHg total dissolved',&
         output=output_time_step_integrated) ! dissolved organic and inorganic
    call self%register_diagnostic_variable(&
         self%id_MeHg_tot,'MeHg_total','mmol/m**3',&
         'MeHg total',&
         output=output_time_step_integrated) ! dissolved and particulate 
    call self%register_diagnostic_variable(&
         self%id_Hg_tot_diss,'Hg_tot_diss','mmol/m**3',&
         'Hg total dissolved',&
         output=output_time_step_integrated) ! dissolved  
    call self%register_diagnostic_variable(&
         self%id_Hg_tot,'Hg_total','mmol/m**3',&
         'Hg total',&
         output=output_time_step_integrated) ! dissolved and particulate 
    call self%register_diagnostic_variable(&
         self%id_MeHg_procent_diss,'MeHg_procent_diss','%',&
         'MeHg_procent_diss',&
         output=output_time_step_integrated) ! share of MeHg 
    call self%register_diagnostic_variable(&
         self%id_Kad_Hg2,'Kad_Hg2','-',&
         'Kad_Hg2',&
         output=output_time_step_integrated) ! conditional partitioning coef
    call self%register_diagnostic_variable(&
         self%id_Kad_MeHg,'Kad_MeHg','-',&
         'Kad_MeHg',&
         output=output_time_step_integrated) ! conditional partitioning coef
    call self%register_diagnostic_variable(&
         self%id_hg2_fe3_compl,'hg2_fe3_compl','mmol/m**3',&
         'hg2_fe3_compl',&
         output=output_time_step_integrated) ! adsorbtion of Hg2 on Fe3 
    call self%register_diagnostic_variable(&
         self%id_hg2_mn4_compl,'hg2_mn4_compl','mmol/m**3',&
         'hg2_mn4_compl',&
         output=output_time_step_integrated) ! adsorbtion of Hg2 on Mn4    
    call self%register_diagnostic_variable(&
         self%id_mehg_fe3_compl,'mehg_fe3_compl','mmol/m**3',&
         'mehg_fe3_compl',&
         output=output_time_step_integrated) ! adsorbtion of MeHg on Fe3  
    call self%register_diagnostic_variable(&
         self%id_mehg_mn4_compl,'mehg_mn4_compl','mmol/m**3',&
         'mehg_mn4_compl',&
         output=output_time_step_integrated) ! adsorbtion of MeHg on Mn4 
 !Register diagnostic dependencies
   call self%register_dependency(self%id_Hplus,&
         'Hplus', 'mmol/m**3','H+ Hydrogen')
   call self%register_dependency(self%id_temp, standard_variables%temperature) 
   call self%register_dependency(self%id_depth,standard_variables%pressure)  
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
!   call self%register_dependency(self%id_Izt,'Izt','W/m2','downwelling_photosynthetic_radiative_flux')
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
! !DESCRIPTION: This module descibes biogeochemical transformation of mercury (Hg) in the seawater
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_hg),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman, Svetlana Pakhomova, Evgeniy Yakushev
!
! !LOCAL VARIABLES:
   real(rk) ::  temp, O2, Mn4, Fe3, depth
   real(rk) ::  Hg0, Hg2, MeHg, HgS, Iz, H2S, Om_HgS, Hg2_flux, SO4 
   real(rk) ::  Phy, Het, Bhae, Baae, Bhan, Baan, DON, POML, POMR
   real(rk) ::  Hg2_biota, Hg2_POM, Hg2_DOM
   real(rk) ::  Hg2_Mn4, Hg2_Fe3, Hg2_free, Hg2_tot, Hg2_tot_diss
   real(rk) ::  MeHg_biota, MeHg_POM, MeHg_DOM
   real(rk) ::  MeHg_Mn4, MeHg_Fe3, MeHg_free, MeHg_tot, MeHg_tot_diss
   real(rk) ::  Hg_tot, Hg_tot_diss, MeHg_procent_diss
   real(rk) ::  hgs_form, hgs_diss, hgs_ox, hg2_hg0, hg0_hg2, mehg_h2s
   real(rk) ::  hg0_irr_ox, hg2_irr_red, mehg_irr_degr, hg2_mehg, mehg_hg2
   real(rk) ::  hg2_fe3_compl, hg2_mn4_compl, mehg_fe3_compl, mehg_mn4_compl
   real(rk) ::  dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM
   real(rk) ::  dHg2, dHg0, dMeHg, Kad_Hg2, Kad_MeHg
  !diagnostic variables dependencies
   real(rk):: Hplus
   
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_depth,depth)            ! depth
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation

    ! Our own state variables
    _GET_(self%id_Hg0,Hg0)
    _GET_(self%id_Hg2,Hg2)
    _GET_(self%id_Hg2_biota,Hg2_biota) 
    _GET_(self%id_Hg2_POM,Hg2_POM) 
    _GET_(self%id_Hg2_DOM,Hg2_DOM) 
    _GET_(self%id_Hg2_free,Hg2_free)
    _GET_(self%id_Hg2_Fe3,Hg2_Fe3)
    _GET_(self%id_Hg2_Mn4,Hg2_Mn4)
    _GET_(self%id_MeHg,MeHg) 
    _GET_(self%id_MeHg_biota,MeHg_biota) 
    _GET_(self%id_MeHg_POM,MeHg_POM) 
    _GET_(self%id_MeHg_DOM,MeHg_DOM) 
    _GET_(self%id_MeHg_free,MeHg_free)
    _GET_(self%id_MeHg_Fe3,MeHg_Fe3) 
    _GET_(self%id_MeHg_Mn4,MeHg_Mn4)
    _GET_(self%id_HgS,HgS)
   ! other modules state variables
    _GET_(self%id_Phy,Phy) 
    _GET_(self%id_Het,Het)
    _GET_(self%id_POML,POML)
    _GET_(self%id_POMR,POMR)
    _GET_(self%id_DON,DON)
    _GET_(self%id_Baae,Baae)   
    _GET_(self%id_Bhae,Bhae)
    _GET_(self%id_Baan,Baan)   
    _GET_(self%id_Bhan,Bhan)
    _GET_(self%id_H2S,H2S)
    _GET_(self%id_SO4,SO4)
    _GET_(self%id_Mn4,Mn4)
    _GET_(self%id_Fe3,Fe3)
    _GET_(self%id_O2,O2)  
    !diagnostic
    _GET_(self%id_Hplus,Hplus)
!-----------------------------------------------------------------
    ! Hg species (Knigthes 2008)
!% Hg0 bioreduction  Hg0 -> Hg2+  ()
    hg0_hg2=self%K_hg0_hg2*Hg0         
!% Hg2 biooxydation  Hg2+ + 0.5O2 + 2H+-> Hg0 + H2O   ()
    hg2_hg0=self%K_hg2_hg0*Hg2*0.5*(1.+tanh(o2-self%O2s_nf))  
!% Hg2 methylation Hg2+  -> MeHg   ()
    hg2_mehg=self%K_hg2_mehg*Hg2*0.5*(1.+tanh(Bhan-50.0_rk))
!% MeHg demethylation MeHg  -> Hg2+   ()
    mehg_hg2=(self%K_mehg_hg2/10.0_rk + self%K_mehg_hg2*0.5*(1.+tanh(Bhan-5.0_rk)))*MeHg 
!% MeHg reduction MeHg + H2S -> Hg2+ + ???   ()
    mehg_h2s=self%K_mehg_h2s*MeHg*(0.5_rk+0.5_rk*tanh(H2S+50.0_rk))
!% HgS saturarion state
    Om_HgS=H2S*Hg2/(self%K_HgS) 
!% HgS formation Hg2+ + H2S -> HgS + 2H+ ()
    hgs_form=max(0._rk,self%K_hgs_form*max(0._rk,(Om_HgS-1._rk)))
    if (Hg2<0.000001.or.H2S<0.01) hgs_form=0._rk
!% HgS dissolution  HgS + 2H+ -> Hg2+ + H2S   ()
    hgs_diss=self%K_hgs_diss*HgS*max(0._rk,(1._rk-Om_HgS))
!    if (HgS<0.000001) hgs_diss=0._rk 
!% HgS oxydation  HgS + 2O2 -> Hg2+ + SO42-  ()
    hgs_ox=self%K_hgs_ox*HgS*(0.5_rk+0.5_rk*tanh(O2+1.0_rk))     
    !hgs_form = 0.0_rk
    !hgs_diss = 0.0_rk
    !hgs_ox = 0.0_rk
!% Hg2 photo reduction  Hg2+ -> Hg0   ()
    hg2_irr_red=self%K_hg2_irr_red*Hg2*Iz/25.*exp(1._rk-Iz/25.)
!% Hg0 photo oxydation  Hg0 -> Hg2+   ()
    hg0_irr_ox=self%K_hg0_irr_ox*Hg0*Iz/25.*exp(1._rk-Iz/25.)
!% MeHg photo degradation MeHg  -> Hg0   
    mehg_irr_degr=self%K_mehg_irr_degr*MeHg*Iz/25.*exp(1._rk-Iz/25.)
  !_______
  ! Hg(II)
  !
  ! partitioning betweeen dissolved Hg(II) and OM
    call partit (Hg2, Hg2_biota, Hg2_POM, Hg2_DOM, &
                 Phy, Het, Baae, Bhae, Baan, Bhan, POML, DON, &
                 dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, &
                 self%Kow_bio_Hg2,self%Kow_POM_Hg2,self%Kow_DOM_Hg2)
    
    dSubst_dis=self%K_relax*dSubst_dis
    dSubst_biota=self%K_relax*dSubst_biota
    dSubst_POM=self%K_relax*dSubst_POM
    dSubst_DOM=self%K_relax*dSubst_DOM
    
   _SET_ODE_(self%id_Hg2_biota,dSubst_biota)
   _SET_ODE_(self%id_Hg2_POM,dSubst_POM)
   _SET_ODE_(self%id_Hg2_DOM,dSubst_DOM)
   _SET_ODE_(self%id_Hg2_free,0.0)

 !! Sorption of Hg(II) on Mn oxides
     hg2_mn4_compl = 0.0_rk 
   _SET_ODE_(self%id_Hg2_Mn4,hg2_mn4_compl)

 !! Sorption of Hg(II) on Fe oxides
       Kad_Hg2=self%KHg2_Fe3*self%Sad_Fe3*Fe3/(Hplus*1000000._rk+self%KHg2_Fe3*Hg2)
     hg2_fe3_compl= self%K_relax*Kad_Hg2*(Hg2+Hg2_Fe3)/(1.0_rk+Kad_Hg2)-Hg2_Fe3
   _SET_ODE_(self%id_Hg2_Fe3,hg2_fe3_compl)

    dHg2= hg0_hg2-hg2_hg0-hg2_mehg+mehg_hg2 &
         +hg0_irr_ox-hg2_irr_red+mehg_h2s &
         -hgs_form+hgs_diss+hgs_ox &
         +dSubst_dis-hg2_fe3_compl-hg2_mn4_compl
   _SET_ODE_(self%id_Hg2, dHg2)
   _SET_ODE_(self%id_HgS,hgs_form-hgs_diss-hgs_ox)
  ! Hg2 related diagnostics
    Hg2_tot_diss= Hg2+Hg2_DOM
   _SET_DIAGNOSTIC_(self%id_Hg2_tot_diss,Hg2_tot_diss)
    Hg2_tot     = Hg2+Hg2_Fe3+Hg2_biota+Hg2_POM+Hg2_DOM !+Hg2_Mn4
   _SET_DIAGNOSTIC_(self%id_Hg2_tot,Hg2_tot)
   
  !________
  ! MeHg
   !
  ! partitioning betweeen dissolved MeHg and OM
    call partit (MeHg, MeHg_biota, MeHg_POM, MeHg_DOM, &
                 Phy, Het, Baae, Bhae, Baan, Bhan, POML, DON, &
                 dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, &
                 self%Kow_bio_MeHg,self%Kow_POM_MeHg,self%Kow_DOM_MeHg)
        
    dSubst_dis=self%K_relax*dSubst_dis
    dSubst_biota=self%K_relax*dSubst_biota
    dSubst_POM=self%K_relax*dSubst_POM
    dSubst_DOM=self%K_relax*dSubst_DOM
    
   _SET_ODE_(self%id_MeHg_biota,dSubst_biota)
   _SET_ODE_(self%id_MeHg_POM,dSubst_POM)
   _SET_ODE_(self%id_MeHg_DOM,dSubst_DOM)
   _SET_ODE_(self%id_MeHg_free,0.0)

 !! Sorption of Hg(II) on Mn oxides
     mehg_mn4_compl = 0.0_rk 
   _SET_ODE_(self%id_MeHg_Mn4,mehg_mn4_compl)
    
 !! Sorption of MeHg on Fe oxides  
     Kad_MeHg=self%KMeHg_Fe3*self%Sad_Fe3*Fe3/(Hplus*1000000._rk+self%KMeHg_Fe3*MeHg)
    mehg_fe3_compl= self%K_relax*Kad_MeHg*(MeHg+MeHg_Fe3)/(1.0_rk+Kad_MeHg)-MeHg_Fe3
   _SET_ODE_(self%id_MeHg_Fe3,mehg_fe3_compl)

    dMeHg=hg2_mehg-mehg_hg2-mehg_irr_degr-mehg_h2s &
         +dSubst_dis-mehg_fe3_compl-mehg_mn4_compl
   _SET_ODE_(self%id_MeHg,dMeHg)
  ! MeHg related diagnostics
    MeHg_tot_diss=MeHg+MeHg_DOM
    _SET_DIAGNOSTIC_(self%id_MeHg_tot_diss,MeHg_tot_diss)
    MeHg_tot=MeHg+MeHg_Fe3+MeHg_biota+MeHg_POM+MeHg_DOM !+MeHg_Mn4
    _SET_DIAGNOSTIC_(self%id_MeHg_tot,MeHg_tot)

  !_______
  ! Hg0
    dHg0=-hg0_hg2+hg2_hg0 &
          - hg0_irr_ox+hg2_irr_red+mehg_irr_degr
   _SET_ODE_(self%id_Hg0,dHg0)
   
    Hg_tot_diss=MeHg_tot_diss+Hg2_tot_diss+Hg0
    
    Hg_tot=1.0_rk*(HgS +MeHg_tot+Hg2_tot+Hg0)
    
    MeHg_procent_diss=100.0_rk*MeHg_tot_diss/Hg_tot_diss
 !! other modules variables
   _SET_ODE_(self%id_Baan,0.0_rk)   
   _SET_ODE_(self%id_Bhan,0.0_rk)
   _SET_ODE_(self%id_Fe3,0.0_rk)
   _SET_ODE_(self%id_Mn4,0.0_rk)
   _SET_ODE_(self%id_H2S,-hgs_form+hgs_diss)
   _SET_ODE_(self%id_SO4,+hgs_ox)
   _SET_ODE_(self%id_O2,-hg2_hg0-hgs_ox)
 !! this module diagnostics
      _SET_DIAGNOSTIC_(self%id_Hg_tot_diss,Hg_tot_diss)
      _SET_DIAGNOSTIC_(self%id_Hg_tot,Hg_tot)
      _SET_DIAGNOSTIC_(self%id_Om_HgS,Om_HgS)
      _SET_DIAGNOSTIC_(self%id_Kad_Hg2,Kad_Hg2)
      _SET_DIAGNOSTIC_(self%id_Kad_MeHg,Kad_MeHg)
      _SET_DIAGNOSTIC_(self%id_MeHg_procent_diss,MeHg_procent_diss)
      _SET_DIAGNOSTIC_(self%id_hg2_fe3_compl,hg2_fe3_compl)
      _SET_DIAGNOSTIC_(self%id_hg2_mn4_compl,hg2_mn4_compl)
      _SET_DIAGNOSTIC_(self%id_mehg_fe3_compl,mehg_fe3_compl)
      _SET_DIAGNOSTIC_(self%id_mehg_mn4_compl,mehg_mn4_compl)
 ! not in use
 !     _SET_DIAGNOSTIC_(self%id_Hg2_free,Hg2_free)
 !     _SET_DIAGNOSTIC_(self%id_MeHg_free,MeHg_free)


! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC
!-----------------------------------------------------------------------
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
! Sea water Hg(0) exchange.   
   
! !INPUT PARAMETERS:
   class (type_niva_brom_hg),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: pCO2w, xk, Ox, Q_Hg0, Hg0
   real(rk) :: temp, Kc0, salt
   real(rk) :: Sc, TK, fwind !PML
   real(rk) :: Hg0_air
   real(rk) :: windspeed

   _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_temp,temp)              ! temperature
      _GET_(self%id_Hg0,Hg0)              ! temperature
      
      TK=(temp + 273.15)
      windspeed=5.
      Hg0_air=5.0e-7 !5.0e-8 ! 0.01/200./1000. !convert from ng/l into mmol/m3

! PML
! calculate the Scmidt number and unit conversions
      Sc = 2073.1_rk-125.62_rk*temp+3.6276_rk*temp**2._rk-0.043219_rk*&
           temp**3.0_rk
      fwind = (0.222_rk*windspeed**2_rk+0.333_rk*windspeed)*&
              (Sc/660._rk)**(-0.5_rk)
      fwind=fwind*24._rk/100._rk !convert to m/day
! flux depends on the difference in partial pressures, wind and henry
! here it is rescaled to mmol/m2/d
!          flux = fwind * HENRY * ( PCO2A - PCO2W ) * dcf      

      Q_Hg0= fwind * (Hg0_air- max(0e0,Hg0))

      _SET_SURFACE_EXCHANGE_(self%id_Hg0,Q_Hg0)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
   
   
   

!-----------------------------------------------------------------------
   subroutine partit (Subst_dis,Subst_biota, Subst_POM, Subst_DOM, &
                      Phy, Het, Baae, Bhae, Baan, Bhan, POML, DON, &
                      dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, &
                      Kow_bio, Kow_pom, Kow_dom)
   ! !LOCAL VARIABLES:
   real(rk) :: Subst_dis, Subst_biota, Subst_POM, Subst_DOM,  Subst_tot 
   real(rk) :: Phy, Het, Baae, Bhae, Baan, Bhan, POML, DON
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
   !real(rk) :: rho_FeS= 5.90E7 !    # Density of FeS [mmolFe/m3] (default = 5.90E7 mmolFe/m3)
   !real(rk) :: rho_FeS2=4.17E7 !    # Density of FeS2 [mmolFe/m3] (default = 4.17E7 mmolFe/m3)
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
       
       if(POML<=0.) then 
        sha_pom=0. 
       else
        sha_pom=uMn2lip/1000.*POML  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif     
       
       if(DON<=0.) then 
        sha_dom=0. 
       else
        sha_dom=uMn2lip/1000.*DON  ! Volume(weight in kg, g->kg=/1000) of BIO
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
end module