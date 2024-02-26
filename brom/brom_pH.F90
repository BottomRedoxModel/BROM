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

module fabm_niva_brom_ph
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_ph
    !state variables dependencies
    type(type_state_variable_id):: id_DIC,id_Alk
    type(type_state_variable_id):: id_PO4,id_Si,id_NH4,id_H2S,id_SO4
    !all descriptions are in the initialize subroutine
    type(type_diagnostic_variable_id) :: id_pH,id_Hplus
    !standard variables dependencies
    type(type_dependency_id):: id_temp,id_salt
    !diagnosric variables dependencies
    type(type_dependency_id):: id_Kc1,id_Kc2,id_Kb
    type(type_dependency_id):: id_Kp1,id_Kp2,id_Kp3,id_KSi
    type(type_dependency_id):: id_Knh4,id_Kh2s,id_kso4,id_kflu
    type(type_dependency_id):: id_Kw,id_tot_free
    logical :: Alk_from_Sal ! switch parameter for calculation of total Alk from  (true/false)
    logical :: SO4_from_Sal ! switch parameter for calculation of SO4 from salinity (true/false)
    logical :: Sal_from_Alk ! switch parameter for calculation of salinity from total Alk (true/false)
    
  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class(type_niva_brom_ph),intent(inout),target :: self
    integer,                 intent(in)           :: configunit

    call self%get_parameter(self%Alk_from_Sal,'Alk_from_Sal', &
        '[ - ]','calculate Alk from Salinity', default=.false.)
    call self%get_parameter(self%SO4_from_Sal,'SO4_from_Sal', &
        '[ - ]','calculate SO4 from Salinity', default=.false.)
    call self%get_parameter(self%Sal_from_Alk,'Sal_from_Alk', &
        '[ - ]','calculate Salinity from Alk', default=.false.)
  
    !Register state dependencies
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3',&
         'total dissolved inorganic carbon',required=.false.)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3','phosphate',required=.false.)
    call self%register_state_dependency(&
         self%id_Si, 'Si', 'mmol/m**3','silicate',required=.false.)
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3','ammonium',required=.false.)
    call self%register_state_dependency(&
         self%id_H2S,'H2S','mmol/m**3','hydrogen sulfide',required=.false.)
    call self%register_state_dependency(&
         self%id_SO4,'SO4','mmol/m**3','sulphate',required=.false.)

    !Register diagnostic variables
    call self%register_diagnostic_variable(self%id_ph,&
         'pH','-','pH',standard_variable=&
         standard_variables%ph_reported_on_total_scale)
    call self%register_diagnostic_variable(&
         self%id_Hplus, 'Hplus', 'M','H+ Hydrogen')

    !Register standard dependencies
    call self%register_dependency(&
         self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(&
         self%id_temp,standard_variables%temperature)

    !Register diagnostic dependencies
    call self%register_dependency(&
         self%id_Kc1,'Kc1','-','[H+][HCO3-]/[H2CO3]')
    call self%register_dependency(&
         self%id_Kc2,'Kc2','-','[H+][CO3--]/[HCO3-]')
    call self%register_dependency(&
         self%id_Kb, 'Kb','-','[H+][B(OH)4-]/[B(OH)3]')
    call self%register_dependency(&
         self%id_Kp1,'Kp1','-','[H+][H2PO4-]/[H3PO4]')
    call self%register_dependency(&
         self%id_Kp2,'Kp2','-','[H+][HPO4--]/[H2PO4-]')
    call self%register_dependency(&
         self%id_Kp3,'Kp3','-','[H+][PO4---]/[HPO4--]')
    call self%register_dependency(&
         self%id_KSi,'KSi','-','[H+][H3SiO4-]/[Si(OH)4]')
    call self%register_dependency(&
         self%id_Knh4,'Knh4','-','[H+][NH3]/[NH4]')
    call self%register_dependency(self%id_Kh2s,'Kh2s','-','[H+][HS-]/[H2S]')
    call self%register_dependency(&
         self%id_kso4,'kso4','-','[H+][HSO4-]/[H2SO4]')
    call self%register_dependency(self%id_kflu,'kflu','-','[H+][F-]/[HF]')
    call self%register_dependency(&
         self%id_Kw, 'Kw','-','[H+][OH-]/H2O')
    call self%register_dependency(&
         self%id_tot_free,'tot_free','-','ratio H_Tot/H_free')    
    
  end subroutine initialize
  !
  !does we need to recalculate concentrations from mmol/m3 to mol/kg?
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_ph),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    real(rk):: temp,salt
    real(rk):: DIC,Alk,PO4,Si,NH4,H2S,SO4
    real(rk):: Kc1,Kc2,Kb,Kp1,Kp2,Kp3,KSi
    real(rk):: Knh4,Kh2s,kso4,kflu,Kw,tot_free
    real(rk):: H_,pH
    real(rk):: Bt,flutot

    _LOOP_BEGIN_
      ! Environment
      _GET_(self%id_temp,temp) !temperature
      _GET_(self%id_salt,salt) !salinity
      ! state variables
      _GET_(self%id_Alk,Alk)
      _GET_(self%id_DIC,DIC)
      _GET_(self%id_PO4,PO4)
      _GET_(self%id_Si,Si)
      _GET_(self%id_SO4,SO4)
      !gases
      _GET_(self%id_NH4,NH4)
      _GET_(self%id_H2S,H2S)
      ! Equilibrium constants
      _GET_(self%id_Kc1,  Kc1)
      _GET_(self%id_Kc2,  Kc2)
      _GET_(self%id_Kb,   Kb)
      _GET_(self%id_Kp1,  Kp1)
      _GET_(self%id_Kp2,  Kp2)
      _GET_(self%id_Kp3,  Kp3)
      _GET_(self%id_KSi,  KSi)
      _GET_(self%id_Knh4, Knh4)
      _GET_(self%id_Kh2s, Kh2s)
      _GET_(self%id_kso4, kso4)
      _GET_(self%id_kflu, kflu)
      _GET_(self%id_Kw,   Kw)
      _GET_(self%id_tot_free, tot_free)

      !returns total alkalinity as a function of salinity if prescribed
       if (self%Alk_from_Sal) then      
      !References:
      ! (Lee et al., 2006) for N.Atlantic mmol/kg
       !  Alk = 2305._rk + 53.97_rk*(salt-35._rk) + 2.74_rk*(salt-35._rk)*(salt-35._rk) &
       !      - 1.16_rk*(temp-20._rk)-0.040_rk*(temp -20_rk)*(temp -20_rk)
       !  Alk = salt*59.88_rk+278.79_rk  ! ???
         Alk = salt*0.068_rk*1000._rk   !Murray, 2014
       endif  
       
       ! calculation of salinity from total Alk (to calculate AB, AF and SO4 ot influenced by salinity) 
       if (self%Sal_from_Alk) then
         salt = Alk/0.068_rk/1000._rk   !Murray, 2014
       endif  
           
      !returns total borate concentration in mol/kg-SW
      !References:
      ! Uppstrom (1974),cited by Dickson et al.(2007,chapter 5,p 10)
      ! Millero (1982) cited in Millero (1995)
      !pH scale  : N/A
      Bt = 0.000416_rk*(salt/35._rk)
      Bt = Bt*(1027._rk/1000._rk)*1.e6_rk !mmol/m3
      
      !returns total fluoride concentration in mol/kg-SW
      !References: Culkin (1965) (???)
      !pH scale  : N/A
      flutot = 0.000068_rk*(salt/35._rk)
      flutot = flutot*(1027._rk/1000._rk)*1.e6_rk !mmol/m3
      
      !returns total sulfate concentration in mol/kg-SW
       if(self%SO4_from_Sal) then
      !References:
      ! Morris A.W. and Riley, J.P. (1966) quoted in Dickson et al. (2007)
      ! pH scale  : N/A
         SO4 = (0.1400_rk/96.062_rk)*(salt/1.80655_rk)
         SO4 = SO4*(1027._rk/1000._rk)*1.e6_rk !mmol/m3
       endif  

!-----------------------------------------------------------------
      H_ = ph_solver(Alk, DIC, Bt, &
                     PO4, Si, NH4, H2S, SO4, flutot, &
                     kc1, kc2, kb, kp1, kp2, kp3, ksi, Knh4, Kh2s, &
                     kso4, kflu, kw, tot_free)
      pH = -LOG10(H_)

      _SET_DIAGNOSTIC_(self%id_pH,pH)
      _SET_DIAGNOSTIC_(self%id_Hplus,H_)
    _LOOP_END_
  end subroutine do
  !
  !
  !
  real(rk) function ph_solver(alktot, dictot, bortot, &
          po4tot, siltot, nh4tot, hstot, so4tot, flutot, &
          kc1, kc2, kb, kp1, kp2, kp3, ksi, kn, khs, &
          kso4, kflu, kw, ph_scale) result(ans)
    !adopted from:
    !Munhoven, G.: Mathematics of the total alkalinity-pH equation
    ! pathway to robust and universal solution algorithms:
    ! the SolveSAPHE package v1.0.1, Geosci. Model Dev., 6, 1367
    ! 1388, doi:10.5194/gmd-6-1367-2013, 2013.
    real(rk), intent(in):: alktot !total alkalinity
    real(rk), intent(in):: dictot !total dissolved
                                  !inorganic carbon
    real(rk), intent(in):: bortot !total boron
    real(rk), intent(in):: po4tot !total po4
    real(rk), intent(in):: siltot !total silicate
    real(rk), intent(in):: nh4tot !total nh3
    real(rk), intent(in):: hstot  !total hs
    real(rk), intent(in):: so4tot !total so4
    real(rk), intent(in):: flutot !total fluor
    !dissociation constants mol/kg-SW in total pH scale
    real(rk), intent(in):: kc1 !1st of carbonic acid
    real(rk), intent(in):: kc2 !2nd of carbonic acid
    real(rk), intent(in):: kb  !boric acid
    real(rk), intent(in):: kp1, kp2, kp3 !phosphoric acid
    real(rk), intent(in):: ksi !silicic acid
    real(rk), intent(in):: kn  !Ammonia
    real(rk), intent(in):: khs !Hydrogen sulfide
    real(rk), intent(in):: kso4
    real(rk), intent(in):: kflu
    real(rk), intent(in):: kw  !water
    !ph_scale factor - ratio H_Tot/H_free
    real(rk), intent(in):: ph_scale

    !discriminant of the main equation R([H+])=0
    !r - value of main eq., dr - its derivative
    real(rk):: delta, r, dr
    ![H+] starting value and etc.
    real(rk):: initial_h, h, h_prev
    !bracketing bounds for alk
    real(rk):: alkinf, alksup
    !bracketing bounds for [H+]
    real(rk):: h_min, h_max
    !required for convergence test
    real(rk):: h_fac
    !auxiliary
    real(rk):: iter, absmin
    logical:: exitnow

    !calculate initial [H+] value
    call initial_h_do(alktot, dictot, bortot, &
                      kc1, kc2, kb, initial_h)

    !calculate bracketing bounds for alk
    alkinf = -po4tot-so4tot-flutot
    alksup = dictot+dictot+bortot+ &
             po4tot+po4tot+siltot+ &
             nh4tot+hstot

    !calculate discriminant for lower bound
    delta = (alktot-alkinf)**2+4._rk*kw/ph_scale
    !calculate lower bound
    if (alktot >= alkinf) then
      h_min = 2._rk*kw/(alktot-alkinf+sqrt(delta))
    else
      h_min = ph_scale*(-(alktot-alkinf)+sqrt(delta))/2._rk
    end if

    !calculate discriminant for upper bound
    delta = (alktot-alksup)**2+4._rk*kw/ph_scale
    !calculate upper bound
    if (alktot <= alksup) then
      h_max = ph_scale*(-(alktot-alksup)+sqrt(delta))/2._rk
    else
      h_max = 2._rk*kw/(alktot-alksup+sqrt(delta))
    end if

    !main algorithm
    h = max(min(h_max, initial_h), h_min)
    iter = 0
    absmin = huge(1._rk)
    do
      h_prev = h
      !return value (r) of main equation for given [H+]
      !and its derivative (dr)
      !call r_calc(h, alktot, dictot, bortot, po4tot, &
      !    siltot, nh4tot, hstot, so4tot, flutot, &
      !    kc1, kc2, kb, kp1, kp2, kp3, ksi, kn, khs, &
      !    kw,kso4, kflu, ph_scale, r, dr)
      call r_calc_old(h, alktot, dictot, bortot, po4tot, &
                      siltot, nh4tot, hstot, &
                      kc1, kc2, kb, kp1, kp2, kp3, ksi, kn, khs, &
                      kw, r, dr)
      !adapt bracketing interval
      if (r > 0._rk) then
        h_min = h_prev
      else if (r < 0._rk) then
        h_max = h_prev
      else
        !h is root
        exit
      end if

      iter = iter + 1
      if (abs(r) >= 0.5_rk*absmin) then
      !bisection method pH-Alk space
        h = sqrt(h_max*h_min)
        !required for convergence test
        h_fac = (h-h_prev)/h_prev
      else
        h_fac = -r/(dr*h_prev)
        if(abs(h_fac) > 1.0_rk) then
        !Newton-Raphson at pH-Alk space
          h = h_prev*exp(h_fac)
        else
        !Newton-Raphson at H-Alk space
          h = h_prev+h_fac*h_prev
        end if
        !boundaries check
        !if not succeed perform bisection step
        if (h<h_min) then
          h = sqrt(h_prev*h_min)
          h_fac = (h-h_prev)/h_prev
        end if
        if (h>h_max) then
          h = sqrt(h_prev*h_max)
          h_fac = (h-h_prev)/h_prev
        end if
      end if

      absmin = min(abs(r), absmin)
      exitnow = (abs(h_fac) < 1.e-8_rk)
      if (exitnow) exit
    end do

    ans = h
  end function ph_solver
  !
  !
  !
  subroutine initial_h_do(alktot, dictot, bortot, &
          kc1, kc2, kb, initial_h)
    !adopted from:
    !Munhoven, G.: Mathematics of the total alkalinity-pH equation
    ! pathway to robust and universal solution algorithms:
    ! the SolveSAPHE package v1.0.1, Geosci. Model Dev., 6, 1367
    ! 1388, doi:10.5194/gmd-6-1367-2013, 2013.
    !calculates initial value for [H+]
    real(rk), intent(in):: alktot
    real(rk), intent(in):: dictot
    real(rk), intent(in):: bortot
    real(rk), intent(in):: kc1
    real(rk), intent(in):: kc2
    real(rk), intent(in):: kb
    real(rk), intent(out):: initial_h

    real(rk):: dic_alk, bor_alk !auxiliary variables
    real(rk):: a0, a1, a2 !coefficients of the cubic polynomial
    real(rk):: d, d_sqrt !discriminant and square of it
    real(rk):: h_min !local minimum

    if (alktot <= 0._rk) then
      initial_h = 1.e-3_rk
    else if (alktot >= (2._rk*dictot+bortot)) then
      initial_h = 1.e-10_rk
    else
      dic_alk = dictot/alktot
      bor_alk = bortot/alktot

      !coefficients of the cubic polynomial equation
      ![H+]*3+a2[H+]*2+a1[H+]+a0=0
      !derived from alkalinity-pH equation
      !for given alkalinity, dictot, bortot
      a2 = kb*(1._rk-bor_alk)+kc1*(1._rk-dic_alk)
      a1 = kc1*kb*(1._rk-bor_alk-dic_alk) &
          +kc2*(1._rk-(dic_alk+dic_alk))
      a0 = kc2*kb*(1._rk-bor_alk-(dic_alk+dic_alk))

      !discriminant of the quadratic equation
      !(derivative of the cubic equation)
      !for the minimum close to the root
      d = a2*a2-3._rk*a1

      if (d > 0._rk) then
        d_sqrt = sqrt(d)
        if (a2 < 0) then
          h_min = (-a2+d_sqrt)/3._rk
        else
          h_min = -a1/(a2+d_sqrt)
        end if

        initial_h = h_min+sqrt(-(a0+h_min*(a1+h_min*(a2+h_min))) &
                    /d_sqrt)
      else
        initial_h = 1.e-7_rk
      end if
    end if
  end subroutine initial_h_do
  !
  !
  !
  subroutine r_calc(h, alktot, dictot, bortot, po4tot, &
          siltot, nh4tot, hstot, so4tot, flutot, &
          kc1, kc2, kb, kp1, kp2, kp3, ksi, kn, khs, &
          kso4, kflu, kw, ph_scale, r, dr)
    !adopted from:
    !Munhoven, G.: Mathematics of the total alkalinity-pH equation
    ! pathway to robust and universal solution algorithms:
    ! the SolveSAPHE package v1.0.1, Geosci. Model Dev., 6, 1367
    ! 1388, doi:10.5194/gmd-6-1367-2013, 2013.
    !return value of main equation for given [H+]
    !and value of its derivative
    real(rk), intent(in):: h ![H+]
    real(rk), intent(in):: alktot !total alkalinity
    real(rk), intent(in):: dictot !total dissolved
                                  !inorganic carbon
    real(rk), intent(in):: bortot !total boron
    real(rk), intent(in):: po4tot !total po4
    real(rk), intent(in):: siltot !total silicate
    real(rk), intent(in):: nh4tot !total nh3
    real(rk), intent(in):: hstot  !total hs
    real(rk), intent(in):: so4tot !total so4
    real(rk), intent(in):: flutot !total fluor
    !dissociation constants
    real(rk), intent(in):: kc1 !1st of carbonic acid
    real(rk), intent(in):: kc2 !2nd of carbonic acid
    real(rk), intent(in):: kb  !boric acid
    real(rk), intent(in):: kp1, kp2, kp3 !phosphoric acid
    real(rk), intent(in):: ksi !silicic acid
    real(rk), intent(in):: kn  !Ammonia
    real(rk), intent(in):: khs !Hydrogen sulfide
    real(rk), intent(in):: kso4
    real(rk), intent(in):: kflu
    real(rk), intent(in):: kw  !water
    !ph_scale_factor
    real(rk), intent(in):: ph_scale
    !output variables
    real(rk), intent(out):: r, dr

    real(rk):: dic1, dic2, dic, dddic, ddic
    real(rk):: bor1, bor2, bor, ddbor, dbor
    real(rk):: po4_1, po4_2, po4, ddpo4, dpo4
    real(rk):: sil1, sil2, sil, ddsil, dsil
    real(rk):: nh4_1, nh4_2, nh4, ddnh4, dnh4
    real(rk):: h2s_1, h2s_2, h2s, ddh2s, dh2s
    real(rk):: so4_1, so4_2, so4, ddso4,dso4
    real(rk):: flu_1, flu_2, flu, ddflu, dflu
    real(rk):: wat

    !H2CO3 - HCO3 - CO3
    dic1 = 2._rk*kc2 + h*       kc1
    dic2 =       kc2 + h*(      kc1 + h)
    dic  = dictot * (dic1/dic2)
    !B(OH)3 - B(OH)4
    bor1 =       kb
    bor2 =       kb + h
    bor  = bortot * (bor1/bor2)
    !H3PO4 - H2PO4 - HPO4 - PO4
    po4_1 = 3._rk*kp3 + h*(2._rk*kp2 + h* kp1)
    po4_2 =       kp3 + h*(      kp2 + h*(kp1 + h))
    po4   = po4tot * (po4_1/po4_2 - 1._rk) ! Zero level of H3PO4 = 1
    !H4SiO4 - H3SiO4
    sil1 =       ksi
    sil2 =       ksi + h
    sil  = siltot * (sil1/sil2)
    !NH4 - NH3
    nh4_1 =       kn
    nh4_2 =       kn + h
    nh4   = nh4tot * (nh4_1/nh4_2)
    !H2S - HS
    h2s_1 =       khs
    h2s_2 =       khs + h
    h2s   = hstot * (h2s_1/h2s_2)
    !HSO4 - SO4
    so4_1 =       kso4
    so4_2 =       kso4 + h
    so4   = so4tot * (so4_1/so4_2 - 1._rk)
    !HF - F
    flu_1 =       kflu
    flu_2 =       kflu + h
    flu   = flutot * (flu_1/flu_2 - 1._rk)
    !H2O - OH
    wat   = kw/h - h/ph_scale

    r = dic + bor + po4 + sil &
      + nh4 + h2s + so4 + flu &
      + wat - alktot

    !H2CO3 - HCO3 - CO3
    dddic = kc1*kc2 + h*(4._rk*kc2+h*kc1)
    ddic  = -dictot*(dddic/dic2**2)
    !B(OH)3 - B(OH)4
    ddbor = kb
    dbor  = -bortot*(ddbor/bor2**2)
    !H3PO4 - H2PO4 - HPO4 - PO4
    ddpo4 = kp2*kp3 + h*(4._rk*kp1*kp3 &
          + h*(9._rk*kp3 + kp1*kp2     &
          + h*(4._rk*kp2               &
          + h*       kp1)))
    dpo4   = -po4tot * (ddpo4/po4_2**2)
    !H4SiO4 - H3SiO4
    ddsil = ksi
    dsil  = -siltot * (ddsil/sil2**2)
    !NH4 - NH3
    ddnh4 = kn
    dnh4  = -nh4tot * (ddnh4/nh4_2**2)
    !H2S - HS
    ddh2s = khs
    dh2s  = -hstot * (ddh2s/h2s_2**2)
    !HSO4 - SO4
    ddso4 = kso4
    dso4  = -so4tot * (ddso4/so4_2**2)
    !HF - F
    ddflu = kflu
    dflu  = -flutot * (ddflu/flu_2**2)

    dr = ddic + dbor + dpo4 + dsil &
       + dnh4 + dh2s + dso4 + dflu &
       - kw/h**2 - 1._rk/ph_scale
  end subroutine r_calc
  !
  !
  !
  subroutine r_calc_old(H_, Alk, Ct, Bt, Pt, &
          Sit, NHt, H2St, &
          kc1, kc2, kb, kp1, kp2, kp3, ksi, Knh4, Kh2s, &
          kw, r, dr)
    !return value of main equation for given [H+]
    !and value of its derivative
    real(rk), intent(in):: H_ ![H+]
    real(rk), intent(in):: Alk !total alkalinity
    real(rk), intent(in):: Ct !total dissolved
                                  !inorganic carbon
    real(rk), intent(in):: Bt !total boron
    real(rk), intent(in):: Pt !total po4
    real(rk), intent(in):: Sit !total silicate
    real(rk), intent(in):: NHt !total nh3
    real(rk), intent(in):: H2St  !total hs
    !dissociation constants
    real(rk), intent(in):: kc1 !1st of carbonic acid
    real(rk), intent(in):: kc2 !2nd of carbonic acid
    real(rk), intent(in):: kb  !boric acid
    real(rk), intent(in):: kp1, kp2, kp3 !phosphoric acid
    real(rk), intent(in):: ksi !silicic acid
    real(rk), intent(in):: Knh4  !Ammonia
    real(rk), intent(in):: Kh2s !Hydrogen sulfide
    real(rk), intent(in):: kw  !water
    !output variables
    real(rk), intent(out):: r, dr

    real(rk):: T1, T2, T12, K12p, K123p, HKR123p

    T1  = H_/Kc1
    T2  = H_/Kc2
    T12 = (1e0+T2+T1*T2)

    K12p  = Kp1*Kp2
    K123p = Kp1*Kp2*Kp3
    HKR123p = 1e0/(((H_+Kp1)*H_+K12p)*H_+K123p)

    !"Alk"=[HCO3-]+2[CO3--]
    !i.e.=([H]/Kc2+2.)*Ct/(1+[H]/Kc2+[H]*[H]/(Kc1*Kc2)) carbonate alkalinity
    r = Ct*(2.+H_/Kc2)/(1.+H_/Kc2+H_/Kc1*H_/Kc2) &
    ![B(OH)4-] i.e.= Btot*Kb/(Kb+[H+]) boric alkalinity
    + Bt*Kb/(Kb+H_) &
    ![OH-] i.e.= Kw/[H+]
    + Kw/H_ &
    ![H+]
    - H_ &
    !Alk_tot
    - Alk &
    ![HPO4--]+2.*[PO4---]-[H3PO4-] i.e. phosphoric alkalinity
    + Pt*((Kp1*Kp2-H_*H_)*H_+2e0*Kp1*Kp2*Kp3) &
    / (((H_+Kp1)*H_+Kp1*Kp2)*H_+Kp1*Kp2*Kp3) &
    ![H3SiO4-] i.e.=Sit*KSi/(KSi+[H+])! silicate alkalinity
    + Sit*KSi/(KSi+H_) &
    ![HS-] i.e.=[H2St]*Kh2s1/(Kh2s1+[H+]) hydrogen sulfide alkalinity
    + H2St*Kh2s/(Kh2s+H_) &
    ![NH3] i.e.=NHt*Knh4/(Knh4+[H+]) ammonia alkalinity
    + NHt*Knh4/(Knh4+H_)

    !d(AH)/d[H+]
    dr = Ct*(T2-4e0*T1-T12)/(Kc2*T12*T12) &
    - Bt*Kb/(Kb+H_)**2. &
    - Kw/H_**2.-1. &
    - Pt*((((Kp1*H_+4e0*K12p)*H_ &
    + Kp1*K12p+9e0*K123p)*H_ &
    + 4e0*Kp1*K123p)*H_+K12p*K123p)*HKR123p*HKR123p &
    - Sit*KSi/(KSi+H_)**2. &
    - H2St*Kh2s/(Kh2s+H_)**2. &
    - NHt*Knh4/(Knh4+H_)**2.
  end subroutine r_calc_old
end module
