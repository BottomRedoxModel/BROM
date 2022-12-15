!-----------------------------------------------------------------------
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_bubble
    use fabm_types
    implicit none
    private
    type,extends(type_base_model),public :: type_niva_brom_bubble
   !all descriptions are in initialize subroutine
   !variables allocated here  
    type(type_state_variable_id)         :: id_Bubble, id_DIC, id_CH4 
   !state variables dependencies
   !global dependencies
    type (type_dependency_id)            :: id_temp, id_salt, id_pres, id_depth
   !diagnostic variables
    type(type_diagnostic_variable_id)    :: id_P_B, id_flux1bubble, id_Bubble_dissolution, id_r_bub
   ! diagnostic dependencies
    type(type_dependency_id):: id_pCO2
   !Model parameters
    real(rk):: W_float !floating rate
    real(rk):: R, P_A, pi, g, TK, N_bub, sigma, Henry, Df, what_gas !K_Bubble_diss, 
    
    contains
        procedure :: initialize
        procedure :: do
    end type
    
    contains
    !
  subroutine initialize(self,configunit)
    class (type_niva_brom_bubble), intent(inout), target :: self
    integer,                       intent(in)            :: configunit

   !-----Model parameters------
   !Sinking
    call self%get_parameter(self%W_float,&
      'W_float','[m/day]','Rate of floating of bubble',  default=5.00_rk)
   !Specific rates of biogeochemical processes
!    call self%get_parameter(self%K_Bubble_diss,&
!      'K_Bubble_diss', '[1/d]','K_Bubble_diss', default=100.0_rk)   
    call self%get_parameter(self%R,&
      'R', '[(n m)/(TK Mol)]','R', default=8.314_rk)   
    call self%get_parameter(self%P_A,&
      'P_A', '[pa]', 'P_A', default=101325.0_rk)   
    call self%get_parameter(self%pi,&
      'pi', '[ - ]','pi', default=3.14159_rk)   
    call self%get_parameter(self%g,&
      'g', '[m/s2]', 'gravity acceleration', default=9.8_rk)   
    call self%get_parameter(self%TK,&
      'TK', '[deg Kelvin]', 'Celsius to Kelvin', default=273.15_rk)
    call self%get_parameter(self%sigma,&
      'sigma', '[n/m]','sigma', default=72.86E-3_rk)
    call self%get_parameter(self%N_bub,&
      'N_bub', '[ - ]', 'Number of bubbles', default=100000._rk)   
    call self%get_parameter(self%what_gas,&
      'what_gas', '[ - ]', 'what gas in bubbles', default=100000._rk)
    call self%get_parameter(self%Henry,&
      'Henry', '[ - ]', 'Henry constant', default=100000._rk)   
    call self%get_parameter(self%Df,&
      'Df', '[ m2/s ]', 'CO2 diffusivity', default=0.000001_rk)   
   !Register state variables
    call self%register_state_variable(self%id_Bubble,&
      'Bubble', 'mmol/m**3','Bubble',minimum=0.0_rk,vertical_movement=-self%W_float/86400._rk)
   !Register state dependencies
    call self%register_state_dependency(self%id_DIC, &
      'DIC', 'mmol/m**3','DIC')
    call self%register_state_dependency(self%id_CH4, &
      'CH4', 'mmol/m**3','CH4')
   !Register diagnostic variables 
    call self%register_diagnostic_variable(self%id_r_bub,&
       'r_bub', 'm','radius of the bubble', output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_P_B,&
       'P_B','m','pressure inside the bubbles', output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_flux1bubble,&
       'flux1bubble','uM/d','1 bubble flux', output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Bubble_dissolution,&
       'Bubble_dissolution','uM/d','Bubble_dissolution', output=output_time_step_integrated)
   !Register global diagnostic dependencies
    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_dependency(self%id_pres,standard_variables%pressure) !in (db)
    call self%register_dependency(self%id_depth,standard_variables%depth)   !in (m)
    call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   !Register diagnostic dependencies
    call self%register_dependency(self%id_pCO2,'pCO2','ppm','partial pressure of CO2')
        
   !Specify that rates are per day (default: per second)
     self%dt = 86400._rk
  end subroutine initialize
    !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_bubble),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_
   !state variables
    real(rk):: Bubble, DIC, CH4
   !diagnostic variables
    real(rk):: r_bub, P_B, flux1bubble,Bubble_dissolution
   !increments and misc
     real(rk):: d_Bubble, Vol_tot, sigma_tau, sigma, Sc, k_B
     real(rk):: temp, salt, pressure, depth, pCO2, Henry_corr
   _LOOP_BEGIN_
   !Retrieve current variable values
   !state variables
    _GET_(self%id_Bubble,Bubble)
    _GET_(self%id_DIC,DIC)
    _GET_(self%id_CH4,CH4)
    _GET_(self%id_temp,temp) ! temperature
    _GET_(self%id_salt,salt) ! temperature
    _GET_(self%id_pres,pressure) ! pressure in dbar
    _GET_(self%id_pCO2,pCO2)

    depth=pressure-10.0_rk

    ! volume of all bubbles, 10.3 converts m into pa; 0.001 converts mmol into mol
    Vol_tot=Bubble*0.001_rk*self%R*(self%TK+temp)/(self%P_A+depth*self%P_A/10.3_rk)

    ! radius of the bubble:
    r_bub =(Vol_tot/self%N_bub*3.0_rk/4.0_rk/self%pi)**(1.0_rk/3.0_rk)

    ! pressure sigma-tau (kg/m3):
    call svan(salt, temp, depth, sigma_tau)

!real(rk):: Kh_theta,coef_temp,koef_salt,pressure,kh_corrected
!Henry = 3.3e-4 mol/m3/Pa
!kh_theta = 3.3e-4  * 101.325 !mol/Atm 
!coef_temp = 2400 ![K]
!koef_salt = 10**(-0.5*0.127) !0.127 L/Mole
!pressure = 1 + depth/10 !in Atm
!kh_corrected = kh_theta * exp(coef_temp*(1/abs_temp - 1/298.15))* 1.e6 * pressure * koef_salt !Micromol/l

   ! Henry constant correction for temperature (K), pressure (atm) and "salinity=35psu" (Mol):
     Henry_corr=self%Henry*exp(2400._rk*(1._rk/(self%TK+temp) - 1._rk/298.15_rk))* 1.e6 * (1._rk + depth/10._rk) * 10**(-0.5*0.127)

   ! gas-water transfer rate kB (m/s):
    if (r_bub.ge.0.0025) then
      k_B =0.065*(self%Df**0.5_rk)
    else
      k_B =1.13_rk*(0.2_rk/(0.45_rk+40*r_bub))*(self%Df**0.5_rk) 
    endif

    if(r_bub.gt.0.000001) then
    ! pressure inside one bubble (Pa):
        P_B = self%P_A + (1000._rk+sigma_tau)*self%g*depth + 2*self%sigma/r_bub
    ! one bubble dissolution:
      flux1bubble = 4._rk*self%pi*(r_bub)**2._rk *k_B &
        *(P_B*self%Henry - 101325._rk/1000000.0_rk*pCO2*self%Henry) & !convert pCO2 from uatm into Pa
        *86400.0_rk ! convert to per day
    else
      P_B = self%P_A + (1000._rk+sigma_tau)*self%g*depth
      flux1bubble = 0.0_rk
    endif

    ! all bubbles dissolution 
    Bubble_dissolution = self%N_bub*flux1bubble
    d_Bubble = Bubble_dissolution

    _SET_ODE_(self%id_Bubble,-d_Bubble)
    if (self%what_gas.lt.1.0_rk) then
    _SET_ODE_(self%id_DIC, d_Bubble)
    _SET_ODE_(self%id_CH4, 0.0_rk)
    else
    _SET_ODE_(self%id_DIC, 0.0_rk)
    _SET_ODE_(self%id_CH4, d_Bubble)
    endif
    _SET_DIAGNOSTIC_(self%id_r_bub,r_bub)
    _SET_DIAGNOSTIC_(self%id_P_B,P_B)
    _SET_DIAGNOSTIC_(self%id_flux1bubble,flux1bubble)
    _SET_DIAGNOSTIC_(self%id_Bubble_dissolution,Bubble_dissolution)
    
    _LOOP_END_
  end subroutine do

!=======================================================================================================================
    subroutine svan(s, t, po, sigma)
!/*
!c------specific volume anomaly based on 1980 equation of state for
!c      seawater and 1978 practical salinity scale
!c      pressure          PO     decibars
!c      temperature        T     degree celsius (IPTS-68)
!c      salinitty          S     (PSS-78)
!c      spec.vol.anom.  SVAN     1.0E-8 m**3/Kg
!c      density anom.   SIGMA    Kg/m**3
!*/

    real(rk)  p,sig,sr,r1,r2,r3,s,t,po,sigma
    real(rk)  a,b,c,d,e,a1,b1,aw,bw,k,ko,kw,k35,v350p,sva,dk,gam,pk,dr35p,dvan

    real(rk) r3500, r4 ,dr350
    data  r3500 /1028.1063/, r4/4.8314E-4/,dr350/28.106331/

    p=po/10.
    sr=sqrt(abs(s))

    r1= ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*t &
           -9.095290E-3)*t+6.793952E-2)*t-28.263737
    r2= (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t &
          +8.24493E-1
    r3= (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3

    sig=(r4*s + r3*sr + r2)*s +r1

    v350p=1.0/r3500
    sva=-sig*v350p/(r3500+sig)
    sigma= sig + dr350

    if (p.eq.0.0) return

    e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
       bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
    b = bw + e*s

    d= 1.91075E-4
    c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
    aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
    a = (d*sr + c)*s + aw

    b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
    a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
    kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t &
           +148.4206)*t-1930.06
    ko = (b1*sr + a1)*s + kw

    dk = (b*p+a)*p+ko
    k35 = (5.03217E-5*p+3.359406)*p+21582.27
    gam=p/k35
    pk=1.0-gam
    sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))

    v350p= v350p*pk
    dr35p=gam/v350p
    dvan= sva/(v350p*(v350p+sva))
    sigma = dr350 + dr35p -dvan
    return
    end subroutine svan
!=======================================================================================================================
end module fabm_niva_brom_bubble
