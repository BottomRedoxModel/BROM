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

module fabm_niva_brom_calcium
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_calcium
    !all descriptions are in initialize subroutine
    type(type_state_variable_id):: id_CaCO3
    !state variables dependencies
    type(type_state_variable_id):: id_DIC,id_Alk

    type(type_diagnostic_variable_id):: id_Om_Ca,id_Om_Ar,id_Ca
    type(type_diagnostic_variable_id):: id_CaCO3_form,id_CaCO3_diss

    !standard variables dependencies
    type(type_dependency_id):: id_temp,id_salt,id_pres
    !diagnostic variables dependencies
    type(type_dependency_id):: id_CO3
    type(type_dependency_id):: id_Wadd

    !Model parameters
    real(rk):: K_caco3_diss,K_caco3_form
    real(rk):: WCa, WCa_tot
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
    class(type_niva_brom_calcium),intent(inout),target :: self
    integer,                      intent(in)           :: configunit

    !-----Model parameters------
        !Sinking
    call self%get_parameter(&
         self%WCa,'WCa','[1/day]',&
         'Rate of sinking of detritus (POP, POM)',&
         default=5.00_rk)
    call self%get_parameter(&
         self%WCa_tot,'WCa_tot','[1/day]',&
         'Total accelerated sinking with absorbed Mn hydroxides',&
          default=5.0_rk)
    !Calcium
    call self%get_parameter(&
         self%K_caco3_diss, 'K_caco3_diss', '[1/day]',&
         'CaCO3 dissollution rate constant',&
         default=3.0_rk)
    call self%get_parameter(&
         self%K_caco3_form, 'K_caco3_form', '[1/day]',&
         'CaCO3 precipitation rate constant',&
         default=0.0001_rk)

    !registering variables
    !state variables
    call self%register_state_variable(&
         self%id_CaCO3, 'CaCO3', 'mmol/m**3','CaCO3',&
         minimum=0.0_rk)

    !registering dependencies
    !state
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3',&
         'total dissolved inorganic carbon',required=.false.)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)

    !diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_Ca,'Ca','mmol/m**3','Ca++')
    call self%register_diagnostic_variable(&
         self%id_Om_Ca,'Om_Ca','-','CaCO3-Calcite saturation')
    call self%register_diagnostic_variable(&
         self%id_Om_Ar,'Om_Ar','-','CaCO3-Aragonite saturation')

    call self%register_diagnostic_variable(&
         self%id_CaCO3_form,'CaCO3_form','-','CaCO3 formation')
    call self%register_diagnostic_variable(&
         self%id_CaCO3_diss,'CaCO3_diss','-','CaCO3 dissolution')
    
    !Register environmental dependencies
    call self%register_dependency(self%id_temp,&
         standard_variables%temperature)
    call self%register_dependency(self%id_salt,&
         standard_variables%practical_salinity)
    call self%register_dependency(self%id_pres,&
         standard_variables%pressure)
    !diagnostic
    call self%register_dependency(self%id_CO3,'CO3','mmol/m**3','CO3--')
    call self%register_dependency(self%id_Wadd,'Wadd','[1/day]',&
         'Additional sinking velocity via Mn4 adsorptoin')

    !self%dt = 86400 states that all rates and cross-boundary fluxes
    !(arguments to _SET_ODE_, _SET_SURFACE_EXCHANGE_, etc.) will need
    !to be divided [by FABM] by 86400 to arrive at units in per second
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !---------------------------------------------------
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_calcium),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !standard variables
    real(rk):: temp,salt,pres
    !state variables
    real(rk):: CaCO3
    !diagnostic variables
    real(rk):: Om_Ca,Om_Ar,Ca
    real(rk):: co3
    !increments
    real(rk):: d_DIC,d_CaCO3,d_Alk
    !processes
    real(rk):: CaCO3_form,CaCO3_diss
    !auxiliary
    real(rk):: K_Cal,K_Ara

    _LOOP_BEGIN_
      ! Environment
      _GET_(self%id_pres,pres)        ! pressure in dbar
      _GET_(self%id_temp,temp)        ! temperature
      _GET_(self%id_salt,salt)        ! salinity
      !diagnostic variables
      _GET_(self%id_CO3,CO3)
      !state variable
      _GET_(self%id_CaCO3,CaCO3)

      !
      call CaCO3solub(temp,salt,0.1_rk*pres,&
                      Ca,K_Cal,K_Ara)
      Om_Ca=(co3/1000000._rk)*Ca/K_Cal !Saturation (Omega) for calcite
      Om_Ar=(co3/1000000._rk)*Ca/K_Ara !Saturation (Omega) for aragonite
      !Ca
      !CaCO3 precipitation/dissolution (Luff et al., 2001)
      caco3_form = self%K_caco3_form*max(0._rk & !Ca2+ + CO32- -> CaCO3
                  ,(Om_Ar-1._rk))
      caco3_diss = CaCO3*self%K_caco3_diss & !CaCO3 -> Ca2+ + CO32-
                  *(max(0._rk,(1._rk-Om_Ar)))**4.5_rk
      !DIC
      d_DIC = -caco3_form+caco3_diss
      _SET_ODE_(self%id_DIC,d_DIC)
      !Calcium
      d_CaCO3 = caco3_form-caco3_diss
      _SET_ODE_(self%id_CaCO3,d_CaCO3)
      !Alkalinity changes due to redox reactions:
      d_Alk = (&
             -2._rk*caco3_form & !Ca2+ + CO32- -> CaCO3
             +2._rk*caco3_diss & !CaCO3 -> Ca2+ + CO32-
             )
      _SET_ODE_(self%id_Alk,d_Alk)
      
      _SET_DIAGNOSTIC_(self%id_Ca,Ca)
      _SET_DIAGNOSTIC_(self%id_Om_Ca,Om_Ca)
      _SET_DIAGNOSTIC_(self%id_Om_Ar,Om_Ar)
      _SET_DIAGNOSTIC_(self%id_caco3_diss,caco3_diss)
      _SET_DIAGNOSTIC_(self%id_caco3_form,caco3_form)
    _LOOP_END_
  contains
    !
    !CaCO3 solubility EYA 2010-08-17
    !
    subroutine CaCO3solub(temp,salt,Pbar,Ca,KCal,KAra)
      !*********************************************************************
      !SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
      !Inputs: WhichKs%, Sal, TempCi, Pdbari, TCi, pHi, Kc1, Kc2
      !Outputs: OmegaCa, OmegaAr
      !This calculates omega, the solubility ratio, for calcite and
      !aragonite.
      !This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
      !      where Ksp is the solubility product (either KCa or KAr).
      !*********************************************************************
      !These are from:
      !Mucci, Alphonso, The solubility of calcite and aragonite in seawater
      !      at various salinities, temperatures, and one atmosphere total
      !      pressure, American Journal of Science 283:781-799, 1983.
      !Ingle, S. E., Solubility of calcite in the ocean,
      !      Marine Chemistry 3:301-319, 1975,
      !Millero, Frank, The thermodynamics of the carbonate system in
      !seawater,
      !      Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
      !Ingle et al, The solubility of calcite in seawater at atmospheric
      !pressure and 35%o salinity, Marine Chemistry 1:295-307, 1973.
      !Berner, R. A., The solubility of calcite and aragonite in seawater in
      !      atmospheric pressure and 34.5%o salinity, American Journal of
      !      Science 276:713-730, 1976.
      !Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
      !Culberson,C.H. and Pytkowicz,R.M., Effect of pressure on carbonic
      !acid,
      !      boric acid, and the pHi of seawater, Limnology and Oceanography
      !      13:403-417, 1968.
      !*********************************************************************
      real(rk),intent(in) :: temp, salt, Pbar
      real(rk),intent(out):: Ca, KCal, KAra
      real(rk):: tempK, logKCal, logKAra, RT, &
                 deltaVKCal,KappaKCal,lnKCalfac,deltaVKAra, &
                 KappaKAra, lnKArafac, RGasConstant

      Ca      = 0.02128_rk/40.087_rk*(salt/1.80655_rk)!in mol/kg-SW
      tempK   = temp+273.15_rk

      !CalciteSolubility:
      !Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
      logKCal = -171.9065_rk-0.077993_rk*tempK+2839.319_rk/tempK &
                +71.595_rk*LOG10(tempK) &
                +(-0.77712_rk+0.0028426_rk*TempK &
                +178.34_rk/tempK)*salt**(0.5_rk) &
                -0.07711_rk*salt+0.0041249_rk*salt**(1.5_rk)

      !check for T=25C, S=35psu: logKcal = 6.3693
      KCal = 10._rk**(logKCal)!in (mol/kg-SW)^2

      !AragoniteSolubility:
      !Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
      logKAra = -171.945_rk-0.077993_rk*tempK+2903.293_rk/tempK &
                +71.595_rk*LOG10(tempK) &
                +(-0.068393_rk+0.0017276_rk*tempK &
                +88.135_rk/tempK)*salt**(0.5_rk) &
                -0.10018_rk*salt+0.0059415_rk*salt**(1.5_rk)
      !check for T=25C, S=35psu: logKcal = 6.1883
      KAra = 10.**(logKAra)!in (mol/kg-SW)^2

      !PressureCorrectionForCalcite:
      !Ingle, Marine Chemistry 3:301-319, 1975
      !same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
      !has typos (-.5304, -.3692, and 10^3 for Kappa factor)
      RGasConstant = 83.1451_rk
      deltaVKCal = -48.76_rk+0.5304_rk*temp
      KappaKCal  = (-11.76_rk+0.3692_rk*temp)/1000._rk
      lnKCalfac  = (-deltaVKCal+0.5_rk*KappaKCal*Pbar) &
                   *Pbar/(RGasConstant*tempK)
      KCal       = KCal*exp(lnKCalfac)

      !PressureCorrectionForAragonite:
      !Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
      !same as Millero, GCA 1995 except for typos (-.5304, -.3692,
      !and 10^3 for Kappa factor)
      deltaVKAra = deltaVKCal+2.8_rk
      KappaKAra  = KappaKCal
      lnKArafac  = (-deltaVKAra+0.5_rk*KappaKAra*Pbar) &
                   *Pbar/(RGasConstant*tempK)
      KAra       = KAra*exp(lnKArafac)
    end subroutine CaCO3solub
  end subroutine do

  ! Set increased manganese sinking via MnIV and MnIII oxides formation
  subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
     class (type_niva_brom_calcium), intent(in) :: self
     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
     
     real(rk) :: Wadd, WCa_tot
          
     _LOOP_BEGIN_
  
      _GET_(self%id_Wadd,Wadd)
     
      WCa_tot = self%WCa + Wadd
  
      _ADD_VERTICAL_VELOCITY_(self%id_CaCO3, WCa_tot)
!      _ADD_VERTICAL_VELOCITY_(self%id_CaCO3, self%WCa)
  
     _LOOP_END_

  end subroutine get_vertical_movement
end module
