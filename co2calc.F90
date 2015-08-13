MODULE co2calc

  !-----------------------------------------------------------------------------
  !   based upon OCMIP2 co2calc
  !
  !   CVS:$Id: co2calc.F90 941 2006-05-12 21:36:48Z klindsay $
  !   CVS:$Name$
  !-----------------------------------------------------------------------------

  use BGC_parms

#ifdef CCSMCOUPLED
   !*** ccsm
  USE shr_vmath_mod
#endif

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !-----------------------------------------------------------------------------

  PRIVATE
  PUBLIC :: co2calc_1point, comp_CO3terms, comp_co3_sat_vals

  !-----------------------------------------------------------------------------
  !   module parameters
  !-----------------------------------------------------------------------------

  real (BGC_r8), parameter ::  &
     c0 = 0.0_BGC_r8,   &
     c1 = 1.0_BGC_r8,   &
     c2 = 2.0_BGC_r8,   &
     c3 = 3.0_BGC_r8,   &
     c10 = 10.0_BGC_r8,   &
     c1000 = 1000.0_BGC_r8,   &
     p5 = 0.5_BGC_r8,   &
     p001 = 0.001_BGC_r8

!  these need to be passed in
  real (BGC_r8), parameter ::  &
!    rho_sw    = 4.1_BGC_r8/3.996_BGC_r8,    & ! density of salt water (g/cm^3)
     rho_sw    = 1.026_BGC_r8,    & ! density of salt water (g/cm^3) from SHR_CONST
     T0_Kelvin = 273.15_BGC_r8             ! zero point for Celsius

  !-----------------------------------------------------------------------------
  !   The current setting of xacc, a tolerance critera, will result in co2star
  !   being accurate to 3 significant figures (xx.y). Making xacc bigger will
  !   result in faster convergence also, but this is not recommended (xacc of
  !   10**-9 drops precision to 2 significant figures).
  !-----------------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: xacc = 1e-10_BGC_r8
  INTEGER(KIND=BGC_i4), PARAMETER :: max_bracket_grow_it = 3
  INTEGER(KIND=BGC_i4), PARAMETER :: maxit = 100

  REAL(KIND=BGC_r8), PARAMETER :: salt_min = 0.1_BGC_r8
  REAL(KIND=BGC_r8), PARAMETER :: dic_min  = salt_min / 35.0_BGC_r8 * 1944.0_BGC_r8
  REAL(KIND=BGC_r8), PARAMETER :: alk_min  = salt_min / 35.0_BGC_r8 * 2225.0_BGC_r8

  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------

  REAL(KIND=BGC_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
       kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
       bt, st, ft, dic, ta, pt, sit

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE co2calc_1point(depth, locmip_k1_k2_bug_fix, lcomp_co3_coeffs, &
       temp, salt, dic_in, ta_in, pt_in, sit_in, phlo, phhi, ph, xco2_in, atmpres, &
       co2star, dco2star, pCO2surf, dpco2)

    !---------------------------------------------------------------------------
    !   SUBROUTINE co2calc_row
    !
    !   PURPOSE : Calculate delta co2*, etc. from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    LOGICAL(KIND=BGC_log), INTENT(IN) :: locmip_k1_k2_bug_fix
    LOGICAL(KIND=BGC_log), INTENT(IN) :: lcomp_co3_coeffs
    REAL(KIND=BGC_r8), INTENT(IN) :: &
         depth,    & ! depth (meters)
         temp,     & ! temperature (degrees C)
         salt,     & ! salinity (PSU)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         xco2_in,  & ! atmospheric mole fraction CO2 in dry air (ppmv)
         atmpres     ! atmospheric pressure (atmosphere)

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(INOUT) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(OUT) :: &
         ph,       & ! computed ph values, for initial guess on next time step
         co2star,  & ! CO2*water (nmol/cm^3)
         dco2star, & ! delta CO2 (nmol/cm^3)
         pco2surf, & ! oceanic pCO2 (ppmv)
         dpco2       ! Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4) :: i
    INTEGER(KIND=BGC_i4) :: k

    REAL(KIND=BGC_r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass,  & ! (mmol/m^3) -> (mol/kg)
         co2starair,   & ! co2star saturation
         htotal2

    REAL(KIND=BGC_r8) :: &
         press_bar,    & ! pressure at z=depth (bars)
         xco2,         & ! atmospheric CO2 (atm)
         htotal,       & ! free concentration of H ion
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_BGC_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    k = 1

    !---------------------------------------------------------------------------
    !   compute thermodynamic CO3 coefficients
    !---------------------------------------------------------------------------

!  below is from POP ref_pressure
    press_bar = 0.059808_BGC_r8*(exp(-0.025_BGC_r8*depth) - c1)     &
            + 0.100766_BGC_r8*depth + 2.28405e-7_BGC_r8*depth**2

    IF (lcomp_co3_coeffs) THEN
       CALL comp_co3_coeffs( k, press_bar, temp, salt, k0, k1, k2, ff, &
                            k1_k2_pH_tot=locmip_k1_k2_bug_fix)
    END IF

    !---------------------------------------------------------------------------
    !   compute htotal
    !---------------------------------------------------------------------------

    CALL comp_htotal(k, temp, dic_in, ta_in, pt_in, sit_in, &
                     k1, k2, phlo, phhi, htotal)

    !---------------------------------------------------------------------------
    !   convert xco2 from uatm to atm
    !---------------------------------------------------------------------------

    xco2 = xco2_in * 1e-6_BGC_r8

    !---------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    !
    !   Compute co2starair
    !---------------------------------------------------------------------------

    htotal2 = htotal ** 2
    co2star = dic(1) * htotal2 / &
         (htotal2 + k1 * htotal + k1 * k2)
    co2starair = xco2 * ff * atmpres
    dco2star = co2starair - co2star
    ph = -LOG10(htotal)

    !---------------------------------------------------------------------
    !   Add two output arguments for storing pCO2surf
    !   Should we be using K0 or ff for the solubility here?
    !---------------------------------------------------------------------

    pCO2surf = co2star / ff
    dpCO2    = pCO2surf - xco2 * atmpres

    !---------------------------------------------------------------------
    !   Convert units of output arguments
    !   Note: pCO2surf and dpCO2 are calculated in atm above.
    !---------------------------------------------------------------------

    co2star  = co2star * mass_to_vol
    dco2star = dco2star * mass_to_vol

    pCO2surf = pCO2surf * 1e6_BGC_r8
    dpCO2    = dpCO2 * 1e6_BGC_r8

  END SUBROUTINE co2calc_1point

  !*****************************************************************************

  SUBROUTINE comp_CO3terms(k, depth, lcomp_co3_coeffs, temp, salt, &
       dic_in, ta_in, pt_in, sit_in, phlo, phhi, ph, H2CO3, HCO3, CO3)

    !---------------------------------------------------------------------------
    !   SUBROUTINE comp_CO3terms
    !
    !   PURPOSE : Calculate H2CO3, HCO3, CO3 from
    !             total alkalinity, total CO2, temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4), INTENT(IN) :: k
    LOGICAL(KIND=BGC_log), INTENT(IN) :: lcomp_co3_coeffs
    REAL(KIND=BGC_r8), INTENT(IN) :: &
         depth,    & ! depth (meters)
         temp,     & ! temperature (degrees C)
         salt,     & ! salinity (PSU)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in      ! inorganic silicate (nmol/cm^3)

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(INOUT) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(OUT) :: &
         pH,         & ! computed ph values, for initial guess on next time step
         H2CO3,      & ! Carbonic Acid Concentration
         HCO3,       & ! Bicarbonate Ion Concentration
         CO3           ! Carbonate Ion Concentration

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4) :: i

    REAL(KIND=BGC_r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass,  & ! (mmol/m^3) -> (mol/kg)
         htotal2, denom

    REAL(KIND=BGC_r8) :: &
         htotal,       & ! free concentration of H ion
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_BGC_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    !------------------------------------------------------------------------
    !   compute thermodynamic CO3 coefficients
    !------------------------------------------------------------------------

    IF (lcomp_co3_coeffs) THEN
       CALL comp_co3_coeffs(k, depth, temp, salt, k0, k1, k2, ff, k1_k2_pH_tot=.true.)
    END IF

    !------------------------------------------------------------------------
    !   compute htotal
    !------------------------------------------------------------------------

    CALL comp_htotal(k, temp, dic_in, &
                     ta_in, pt_in, sit_in, k1, k2, &
                     phlo, phhi, htotal)

    !------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49-51)
    !------------------------------------------------------------------------

    htotal2  = htotal ** 2
    denom    = c1 / (htotal2 + k1 * htotal + k1 * k2)
    H2CO3 = dic(1) * htotal2 * denom
    HCO3  = dic(1) * k1 * htotal * denom
    CO3   = dic(1) * k1 * k2 * denom
    ph    = -LOG10(htotal)

    !------------------------------------------------------------------
    !   Convert units of output arguments
    !------------------------------------------------------------------

    H2CO3 = H2CO3 * mass_to_vol
    HCO3  = HCO3 * mass_to_vol
    CO3   = CO3 * mass_to_vol

  END SUBROUTINE comp_CO3terms

  !*****************************************************************************

  SUBROUTINE comp_co3_coeffs(k, depth, temp, salt, sk0, sk1, sk2, sff, k1_k2_pH_tot)

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4), INTENT(IN) :: k
    REAL(KIND=BGC_r8), INTENT(IN) :: &
         depth,    & ! depth (meters)
         temp,     & ! temperature (degrees C)
         salt        ! salinity (PSU)
    LOGICAL(KIND=BGC_log), INTENT(IN) :: k1_k2_pH_tot

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

!maltrud these are scalar versions--need to copy from array(1) due to shr_vmath
    REAL(KIND=BGC_r8), INTENT(OUT) :: &
         sk0,sk1,sk2,     & ! equilibrium constants for CO2 species
         sff              ! fugacity of CO2

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), dimension(1) :: &  ! need to be arrays for shr_vmath
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2


    INTEGER(KIND=BGC_i4) :: i

    REAL(KIND=BGC_r8) :: &
         press_bar       ! pressure at level k [bars]

    REAL(KIND=BGC_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
         salt_lim,     & ! bounded salt
         tk,           & ! temperature (K)
         is,           & ! ionic strength
         scl,          & ! chlorinity
         tk100, tk1002, invtk, dlogtk, is2, sqrtis, &
         s2, sqrts, s15, invRtk, arg, &
         deltaV,Kappa,lnKfac,Kfac, & ! pressure correction terms
         log_1_m_1p005em3_s, &
         log_1_p_tot_sulfate_div_ks

    !---------------------------------------------------------------------------

!   press_bar = ref_pressure(k)
!  below is from POP ref_pressure
    press_bar = 0.059808_BGC_r8*(exp(-0.025_BGC_r8*depth) - c1)     &
            + 0.100766_BGC_r8*depth + 2.28405e-7_BGC_r8*depth**2

    !---------------------------------------------------------------------------
    !   Calculate all constants needed to convert between various
    !   measured carbon species. References for each equation are
    !   noted in the code.  Once calculated, the constants are stored
    !   and passed in the common block "const". The original version
    !   of this code was based on the code by Dickson in Version 2 of
    !   "Handbook of Methods for the Analysis of the Various Parameters
    !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
    !   p25-26).
    !   Derive simple terms used more than once
    !---------------------------------------------------------------------------

    salt_lim = max(salt,salt_min)
    tk       = T0_Kelvin + temp
    tk100    = tk * 1e-2_BGC_r8
    tk1002   = tk100 * tk100
    invtk    = c1 / tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(tk, dlogtk, 1)
#else
    dlogtk   = LOG(tk)
#endif
    invRtk   = (c1 / 83.1451_BGC_r8) * invtk

    is       = 19.924_BGC_r8 * salt_lim / (c1000 - 1.005_BGC_r8 * salt_lim)
    is2      = is * is
#ifdef CCSMCOUPLED
    CALL shr_vmath_sqrt(is, sqrtis, 1)
    CALL shr_vmath_sqrt(salt_lim, sqrts, 1)
#else
    sqrtis   = SQRT(is)
    sqrts    = SQRT(salt_lim)
#endif
    s2       = salt_lim * salt_lim
    scl      = salt_lim / 1.80655_BGC_r8

    arg = c1 - 0.001005_BGC_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(arg, log_1_m_1p005em3_s, 1)
#else
    log_1_m_1p005em3_s = LOG(arg)
#endif

    !---------------------------------------------------------------------------
    !   f = k0(1-pH2O)*correction term for non-ideality
    !   Weiss & Price (1980, Mar. Chem., 8, 347-359;
    !                 Eq 13 with table 6 values)
    !---------------------------------------------------------------------------

    arg = -162.8301_BGC_r8 + 218.2968_BGC_r8 / tk100 + &
          90.9241_BGC_r8 * (dlogtk + LOG(1e-2_BGC_r8)) - 1.47696_BGC_r8 * tk1002 + &
          salt_lim * (.025695_BGC_r8 - .025225_BGC_r8 * tk100 + 0.0049867_BGC_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ff, 1)
#else
    ff = EXP(arg)
#endif
    sff = ff(1)

    !---------------------------------------------------------------------------
    !   K0 from Weiss 1974
    !---------------------------------------------------------------------------

    arg = 93.4517_BGC_r8 / tk100 - 60.2409_BGC_r8 + 23.3585_BGC_r8 * (dlogtk + LOG(1e-2_BGC_r8)) + &
          salt_lim * (.023517_BGC_r8 - 0.023656_BGC_r8 * tk100 + 0.0047036_BGC_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k0, 1)
#else
    k0 = EXP(arg)
#endif
    sk0 = k0(1)

    !---------------------------------------------------------------------------
    !   k1 = [H][HCO3]/[H2CO3]
    !   k2 = [H][CO3]/[HCO3]
    !   if k1_k2_pH_tot == .true., then use
    !      Lueker, Dickson, Keeling (2000) using Mehrbach et al. data on total scale
    !   otherwise, use
    !      Millero p.664 (1995) using Mehrbach et al. data on seawater scale
    !      this is only present to be consistent w/ OCMIP2 code
    !      it should not be used for new runs
    !      the only reason to use it is to be compatible with prior
    !      long spun up runs that had used it
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 3633.86_BGC_r8 * invtk - 61.2172_BGC_r8 + &
             9.67770_BGC_r8 * dlogtk - 0.011555_BGC_r8 * salt_lim + &
             0.0001152_BGC_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 3670.7_BGC_r8 * invtk - 62.008_BGC_r8 + &
             9.7944_BGC_r8 * dlogtk - 0.0118_BGC_r8 * salt_lim + &
             0.000116_BGC_r8 * s2
    END IF
    arg = -LOG(c10) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1, 1)
#else
    k1 = EXP(arg)
#endif
    sk1 = k1(1)

    IF (k > 1) THEN
       deltaV = -25.5_BGC_r8 + 0.1271_BGC_r8 * temp
       Kappa  = (-3.08_BGC_r8 + 0.0877_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k1 = k1 * Kfac
    END IF

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 471.78_BGC_r8 * invtk + 25.9290_BGC_r8 - &
             3.16967_BGC_r8 * dlogtk - 0.01781_BGC_r8 * salt_lim + 0.0001122_BGC_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 1394.7_BGC_r8 * invtk + 4.777_BGC_r8 - &
             0.0184_BGC_r8 * salt_lim + 0.000118_BGC_r8 * s2
    END IF
    arg = -LOG(c10) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2, 1)
#else
    k2 = EXP(arg)
#endif
    sk2 = k2(1)

    IF (k > 1) THEN
       deltaV = -15.82_BGC_r8 - 0.0219_BGC_r8 * temp
       Kappa  = (1.13_BGC_r8 - 0.1475_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k2 = k2 * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   kb = [H][BO2]/[HBO2]
    !   Millero p.669 (1995) using data from Dickson (1990)
    !   CO2SYS states that this in on total pH scale
    !   pressure correction from Millero 1979, p. 1657
    !      omitting salinity contribution
    !---------------------------------------------------------------------------

    arg = (-8966.90_BGC_r8 - 2890.53_BGC_r8 * sqrts - &
           77.942_BGC_r8 * salt_lim + 1.728_BGC_r8 * salt_lim * sqrts - &
           0.0996_BGC_r8 * s2) * invtk + &
          (148.0248_BGC_r8 + 137.1942_BGC_r8 * sqrts + 1.62142_BGC_r8 * salt_lim) + &
          (-24.4344_BGC_r8 - 25.085_BGC_r8 * sqrts - 0.2474_BGC_r8 * salt_lim) * dlogtk + &
          0.053105_BGC_r8 * sqrts * tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kb, 1)
#else
    kb = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -29.48_BGC_r8 + (0.1622_BGC_r8 - 0.002608_BGC_r8 * temp) * temp
       Kappa  = -2.84_BGC_r8 * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kb = kb * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   k1p = [H][H2PO4]/[H3PO4]
    !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4576.752_BGC_r8 * invtk + 115.525_BGC_r8 - &
          18.453_BGC_r8 * dlogtk + &
          (-106.736_BGC_r8 * invtk + 0.69171_BGC_r8) * sqrts + &
          (-0.65643_BGC_r8 * invtk - 0.01844_BGC_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1p, 1)
#else
    k1p = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -14.51_BGC_r8 + (0.1211_BGC_r8 - 0.000321_BGC_r8 * temp) * temp
       Kappa  = (-2.67_BGC_r8 + 0.0427_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k1p = k1p * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   k2p = [H][HPO4]/[H2PO4]
    !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -8814.715_BGC_r8 * invtk + 172.0883_BGC_r8 - &
          27.927_BGC_r8 * dlogtk + &
          (-160.340_BGC_r8 * invtk + 1.3566_BGC_r8) * sqrts + &
          (0.37335_BGC_r8 * invtk - 0.05778_BGC_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2p, 1)
#else
    k2p = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -23.12_BGC_r8 + (0.1758_BGC_r8 - 0.002647_BGC_r8 * temp) * temp
       Kappa  = (-5.15_BGC_r8 + 0.09_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k2p = k2p * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   k3p = [H][PO4]/[HPO4]
    !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -3070.75_BGC_r8 * invtk - 18.141_BGC_r8 + &
          (17.27039_BGC_r8 * invtk + 2.81197_BGC_r8) * sqrts + &
          (-44.99486_BGC_r8 * invtk - 0.09984_BGC_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k3p, 1)
#else
    k3p = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -26.57_BGC_r8 + (0.202_BGC_r8 - 0.003042_BGC_r8 * temp) * temp
       Kappa  = (-4.08_BGC_r8 + 0.0714_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       k3p = k3p * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   ksi = [H][SiO(OH)3]/[Si(OH)4]
    !   Millero p.671 (1995) using data from Yao and Millero (1995)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !      apply boric acid values
    !---------------------------------------------------------------------------

    arg = -8904.2_BGC_r8 * invtk + 117.385_BGC_r8 - &
          19.334_BGC_r8 * dlogtk + &
          (-458.79_BGC_r8 * invtk + 3.5913_BGC_r8) * sqrtis + &
          (188.74_BGC_r8 * invtk - 1.5998_BGC_r8) * is + &
          (-12.1652_BGC_r8 * invtk + 0.07871_BGC_r8) * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ksi, 1)
#else
    ksi = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -29.48_BGC_r8 + (0.1622_BGC_r8 - 0.002608_BGC_r8 * temp) * temp
       Kappa  = -2.84_BGC_r8 * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       ksi = ksi * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   kw = [H][OH]
    !   Millero p.670 (1995) using composite data
    !   following DOE Handbook, 0.015 substracted from constant to
    !   approximately convert from SWS pH scale to total pH scale
    !   pressure correction from Millero 1983
    !      note that deltaV coeffs in Millero 1995 are those actually
    !      freshwater deltaV coeffs from Millero 1983
    !---------------------------------------------------------------------------

    arg = -13847.26_BGC_r8 * invtk + 148.9652_BGC_r8 - 23.6521_BGC_r8 * dlogtk + &
          (118.67_BGC_r8 * invtk - 5.977_BGC_r8 + 1.0495_BGC_r8 * dlogtk) * sqrts - &
          0.01615_BGC_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kw, 1)
#else
    kw = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -20.02_BGC_r8 + (0.1119_BGC_r8 - 0.001409_BGC_r8 * temp) * temp
       Kappa  = (-5.13_BGC_r8 + 0.0794_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kw = kw * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   ks = [H][SO4]/[HSO4], free pH scale
    !   Dickson (1990, J. chem. Thermodynamics 22, 113)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4276.1_BGC_r8 * invtk + 141.328_BGC_r8 - 23.093_BGC_r8 * dlogtk + &
          (-13856.0_BGC_r8 * invtk + 324.57_BGC_r8 - 47.986_BGC_r8 * dlogtk) * sqrtis + &
          (35474.0_BGC_r8 * invtk - 771.54_BGC_r8 + 114.723_BGC_r8 * dlogtk) * is - &
          2698.0_BGC_r8 * invtk * is * sqrtis + &
          1776.0_BGC_r8 * invtk * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ks, 1)
#else
    ks = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -18.03_BGC_r8 + (0.0466_BGC_r8 + 0.000316_BGC_r8 * temp) * temp
       Kappa  = (-4.53_BGC_r8 + 0.09_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       ks = ks * Kfac
    END IF

    !---------------------------------------------------------------------
    !   kf = [H][F]/[HF]
    !   Dickson and Riley (1979) -- change pH scale to total
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------

    arg = c1 + (0.1400_BGC_r8 / 96.062_BGC_r8) * (scl) / ks
#ifdef CCSMCOUPLED
       CALL shr_vmath_log(arg, log_1_p_tot_sulfate_div_ks, 1)
#else
    log_1_p_tot_sulfate_div_ks = LOG(arg)
#endif
    arg = 1590.2_BGC_r8 * invtk - 12.641_BGC_r8 + 1.525_BGC_r8 * sqrtis + &
          log_1_m_1p005em3_s + log_1_p_tot_sulfate_div_ks
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kf, 1)
#else
    kf = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -9.78_BGC_r8 - (0.009_BGC_r8 + 0.000942_BGC_r8 * temp) * temp
       Kappa  = (-3.91_BGC_r8 + 0.054_BGC_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
       Kfac = EXP(lnKfac)
#endif
       kf = kf * Kfac
    END IF

    !---------------------------------------------------------------------
    !   Calculate concentrations for borate, sulfate, and fluoride
    !   bt : Uppstrom (1974)
    !   st : Morris & Riley (1966)
    !   ft : Riley (1965)
    !---------------------------------------------------------------------

    bt = 0.000232_BGC_r8 / 10.811_BGC_r8 * scl
    st = 0.14_BGC_r8 / 96.062_BGC_r8 * scl
    ft = 0.000067_BGC_r8 / 18.9984_BGC_r8 * scl

  END SUBROUTINE comp_co3_coeffs

  !*****************************************************************************

  SUBROUTINE comp_htotal(k, temp, dic_in, ta_in, pt_in, sit_in, &
                         k1, k2, phlo, phhi, htotal)

    !---------------------------------------------------------------------------
    !   SUBROUTINE comp_htotal
    !
    !   PURPOSE : Calculate htotal from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4), INTENT(IN) :: k
    REAL(KIND=BGC_r8), INTENT(IN) :: &
         temp,     & ! temperature (degrees C)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         k1,k2       ! equilibrium constants for CO2 species

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(INOUT) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(OUT) :: &
         htotal      ! free concentration of H ion

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4) :: i

    REAL(KIND=BGC_r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass     ! (mmol/m^3) -> (mol/kg)

    REAL(KIND=BGC_r8) :: &
         x1, x2          ! bounds on htotal for solver

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_BGC_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass
    !---------------------------------------------------------------------------

    dic  = max(dic_in,dic_min) * vol_to_mass
    ta   = max(ta_in,alk_min)  * vol_to_mass
    pt   = max(pt_in,c0)       * vol_to_mass
    sit  = max(sit_in,c0)      * vol_to_mass

    x1 = c10 ** (-phhi)
    x2 = c10 ** (-phlo)

    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------

    CALL drtsafe_row( k, k1, k2, x1, x2, xacc, htotal)

  END SUBROUTINE comp_htotal

  !*****************************************************************************

  SUBROUTINE drtsafe_row(k, k1, k2, x1, x2, xacc, soln)

    !---------------------------------------------------------------------------
    !   Vectorized version of drtsafe, which was a modified version of
    !      Numerical Recipes algorithm.
    !   Keith Lindsay, Oct 1999
    !
    !   Algorithm comment :
    !      Iteration from Newtons method is used unless it leaves
    !      bracketing interval or the dx is > 0.5 the previous dx.
    !      In that case, bisection method is used.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4), INTENT(IN) :: k
    REAL(KIND=BGC_r8), INTENT(IN) :: k1, k2
    REAL(KIND=BGC_r8), INTENT(IN) :: xacc

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(INOUT) :: x1, x2

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(OUT) :: soln

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    LOGICAL(KIND=BGC_log) :: leave_bracket, dx_decrease, mask
    INTEGER(KIND=BGC_i4) ::  i, it
    REAL(KIND=BGC_r8) :: temp
    REAL(KIND=BGC_r8) :: xlo, xhi, flo, fhi, f, df, dxold, dx

    !---------------------------------------------------------------------------
    !   bracket root at each location and set up first iteration
    !---------------------------------------------------------------------------

    it = 0

    DO
       CALL talk_row(k1, k2, x1, flo, df)
       CALL talk_row(k1, k2, x2, fhi, df)

       mask = (flo > c0 .AND. fhi > c0) .OR. &
              (flo < c0 .AND. fhi < c0)

       IF (.NOT. mask) EXIT

       it = it + 1

       IF (it > max_bracket_grow_it) THEN
!         CALL shr_sys_abort('bounding bracket for pH solution not found')
       END IF

       dx = sqrt(x2 / x1)
       x2 = x2 * dx
       x1 = x1 / dx
    END DO

    IF (flo .LT. c0) THEN
       xlo = x1
       xhi = x2
    ELSE
       xlo = x2
       xhi = x1
       temp = flo
       flo = fhi
       fhi = temp
    END IF
    soln = p5 * (xlo + xhi)
    dxold = ABS(xlo - xhi)
    dx = dxold

    CALL talk_row(k1, k2, soln, f, df)

    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

    mask = .true.
    DO it = 1,maxit
       leave_bracket = ((soln - xhi) * df - f) * &
            ((soln - xlo) * df - f) .GE. 0
       dx_decrease = ABS(c2 * f) .LE. ABS(dxold * df)
       IF (leave_bracket .OR. .NOT. dx_decrease) THEN
          dxold = dx
          dx = p5 * (xhi - xlo)
          soln = xlo + dx
          IF (xlo .EQ. soln) mask = .FALSE.
       ELSE
          dxold = dx
          dx = -f / df
          temp = soln
          soln = soln + dx
          IF (temp .EQ. soln) mask = .FALSE.
       END IF
       IF (ABS(dx) .LT. xacc) mask = .FALSE.

       IF (.NOT. mask) RETURN

       CALL talk_row(k1, k2, soln, f, df)

       IF (f .LT. c0) THEN
          xlo = soln
          flo = f
       ELSE
          xhi = soln
          fhi = f
       END IF

    END DO ! iteration loop

#ifdef CCSMCOUPLED
!   CALL shr_sys_abort('lack of convergence in drtsafe_row')
#endif

  END SUBROUTINE drtsafe_row

  !*****************************************************************************

  SUBROUTINE talk_row(k1, k2, x, fn, df)

    !---------------------------------------------------------------------------
    !   This routine computes TA as a function of DIC, htotal and constants.
    !   It also calculates the derivative of this function with respect to
    !   htotal. It is used in the iterative solution for htotal. In the call
    !   "x" is the input value for htotal, "fn" is the calculated value for
    !   TA and "df" is the value for dTA/dhtotal.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(IN) :: k1, k2, x

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(OUT) :: fn, df

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4) :: i

    REAL(KIND=BGC_r8) :: &
         x1, x1_r, x2, x2_r, x3, k12, k12p, k123p, &
         a, a_r, a2_r, da, b, b_r, b2_r, db, c, c_r, &
         kb_p_x1_r, ksi_p_x1_r, c1_p_c_ks_x1_r_r, c1_p_kf_x1_r_r

    !---------------------------------------------------------------------------

    x1 = x
    x1_r = c1 / x1
    x2 = x1 * x1
    x2_r = x1_r * x1_r
    x3 = x2 * x1
    k12 = k1 * k2
    k12p = k1p(1) * k2p(1)
    k123p = k12p * k3p(1)
    a = x3 + k1p(1) * x2 + k12p * x1 + k123p
    a_r = c1 / a
    a2_r = a_r * a_r
    da = c3 * x2 + c2 * k1p(1) * x1 + k12p
    b = x2 + k1 * x1 + k12
    b_r = c1 / b
    b2_r = b_r * b_r
    db = c2 * x1 + k1
    c = c1 + st(1) / ks(1)
    c_r = c1 / c
    kb_p_x1_r = c1 / (kb(1) + x1)
    ksi_p_x1_r = c1 / (ksi(1) + x1)
    c1_p_c_ks_x1_r_r = c1 / (c1 + c * ks(1) * x1_r)
    c1_p_kf_x1_r_r = c1 / (c1 + kf(1) * x1_r)

    !---------------------------------------------------------------------
    !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
    !---------------------------------------------------------------------

    fn = k1 * dic(1) * x1 * b_r &
         + c2 * dic(1) * k12 * b_r &
         + bt(1) * kb(1) * kb_p_x1_r &
         + kw(1) * x1_r &
         + pt(1) * k12p * x1 * a_r &
         + c2 * pt(1) * k123p * a_r &
         + sit(1) * ksi(1) * ksi_p_x1_r &
         - x1 * c_r &
         - st(1) * c1_p_c_ks_x1_r_r &
         - ft(1) * c1_p_kf_x1_r_r &
         - pt(1) * x3 * a_r &
         - ta(1)

    !---------------------------------------------------------------------
    !   df = d(fn)/dx
    !---------------------------------------------------------------------

    df = k1 * dic(1) * (b - x1 * db) * b2_r &
         - c2 * dic(1) * k12 * db * b2_r &
         - bt(1) * kb(1) * kb_p_x1_r * kb_p_x1_r &
         - kw(1) * x2_r &
         + (pt(1) * k12p * (a - x1 * da)) * a2_r &
         - c2 * pt(1) * k123p * da * a2_r &
         - sit(1) * ksi(1) * ksi_p_x1_r * ksi_p_x1_r &
         - c1 * c_r &
         - st(1) * c1_p_c_ks_x1_r_r * c1_p_c_ks_x1_r_r * (c * ks(1) * x2_r) &
         - ft(1) * c1_p_kf_x1_r_r * c1_p_kf_x1_r_r * kf(1) * x2_r &
         - pt(1) * x2 * (c3 * a - x1 * da) * a2_r

  END SUBROUTINE talk_row

  !*****************************************************************************

  SUBROUTINE comp_co3_sat_vals(k, depth, temp, salt, co3_sat_calc, co3_sat_arag)

    !---------------------------------------------------------------------------
    !   SUBROUTINE comp_co3_sat_vals
    !
    !   PURPOSE : Calculate co3 concentration at calcite and aragonite saturation
    !             from temp, salinity (s), press
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4), INTENT(IN) :: k
    REAL(KIND=BGC_r8), INTENT(IN) :: &
         depth,    & ! depth (meters)
         temp,     & ! temperature (degrees C)
         salt        ! salinity (PSU)

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=BGC_r8), INTENT(OUT) :: &
         co3_sat_calc,&! co3 concentration at calcite saturation
         co3_sat_arag  ! co3 concentration at aragonite saturation

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=BGC_i4) :: i, j

    REAL(KIND=BGC_r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         press_bar       ! pressure at level k [bars]

    REAL(KIND=BGC_r8), dimension(1) :: &  !  need to be arrays to use shr_vmath
         salt_lim,     & ! bounded salt
         tk,           & ! temperature (K)
         log10tk,invtk,sqrts,s15,invRtk,arg,&
         K_calc,       & ! thermodynamic constant for calcite
         K_arag,       & ! thermodynamic constant for aragonite
         deltaV,Kappa, & ! pressure correction terms
         lnKfac,Kfac,  & ! pressure correction terms
         inv_Ca          ! inverse of Calcium concentration (mol/kg)

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_BGC_r8 * rho_sw

    !---------------------------------------------------------------------------

!   press_bar = ref_pressure(k)
!  below is from POP ref_pressure
    press_bar = 0.059808_BGC_r8*(exp(-0.025_BGC_r8*depth) - c1)     &
            + 0.100766_BGC_r8*depth + 2.28405e-7_BGC_r8*depth**2

       salt_lim = max(salt,salt_min)
       tk       = T0_Kelvin + temp
#ifdef CCSMCOUPLED
       CALL shr_vmath_log(tk, log10tk, 1)
#else
       log10tk  = LOG(tk)
#endif
!maltrud dividing by non-shared log(10) here.  i guess since it is a scalar?
       log10tk  = log10tk/LOG(c10)
       invtk    = c1 / tk
       invRtk   = (c1 / 83.1451_BGC_r8) * invtk

#ifdef CCSMCOUPLED
       CALL shr_vmath_sqrt(salt_lim, sqrts, 1)
#else
       sqrts    = SQRT(salt_lim)
#endif
       s15      = sqrts * salt_lim

       !------------------------------------------------------------------------
       !   Compute K_calc, K_arag, apply pressure factor
       !   Mucci, Amer. J. of Science 283:781-799, 1983 & Millero 1979
       !------------------------------------------------------------------------

       arg = -171.9065_BGC_r8 - 0.077993_BGC_r8 * tk + 2839.319_BGC_r8 * invtk + 71.595_BGC_r8 * log10tk + &
             (-0.77712_BGC_r8 + 0.0028426_BGC_r8 * tk + 178.34_BGC_r8 * invtk) * sqrts - &
             0.07711_BGC_r8 * salt_lim + 0.0041249_BGC_r8 * s15
       arg = LOG(c10) * arg
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(arg, K_calc, 1)
#else
       K_calc = EXP(arg)
#endif

       IF (k > 1) THEN
          deltaV = -48.76_BGC_r8 + 0.5304_BGC_r8 * temp
          Kappa  = (-11.76_BGC_r8 + 0.3692_BGC_r8 * temp) * p001
          lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
          Kfac = EXP(lnKfac)
#endif
          K_calc = K_calc * Kfac
       END IF

       arg = -171.945_BGC_r8 - 0.077993_BGC_r8 * tk + 2903.293_BGC_r8 * invtk + 71.595_BGC_r8 * log10tk + &
            (-0.068393_BGC_r8 + 0.0017276_BGC_r8 * tk + 88.135_BGC_r8 * invtk) * sqrts - &
            0.10018_BGC_r8 * salt_lim + 0.0059415_BGC_r8 * s15
       arg = LOG(c10) * arg
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(arg, K_arag, 1)
#else
       K_arag = EXP(arg)
#endif

       IF (k > 1) THEN
          deltaV = deltaV + 2.8_BGC_r8
          lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, 1)
#else
          Kfac = EXP(lnKfac)
#endif
          K_arag = K_arag * Kfac
       END IF

    !------------------------------------------------------------------
    !   Compute CO3 concentration at calcite & aragonite saturation
    !------------------------------------------------------------------

    inv_Ca(1) = (35.0_BGC_r8 / 0.01028_BGC_r8) / salt_lim(1) 
    co3_sat_calc = K_calc(1) * inv_Ca(1)
    co3_sat_arag = K_arag(1) * inv_Ca(1)

    !------------------------------------------------------------------
    !   Convert units of output arguments
    !------------------------------------------------------------------

    co3_sat_calc = co3_sat_calc * mass_to_vol
    co3_sat_arag = co3_sat_arag * mass_to_vol

  END SUBROUTINE comp_co3_sat_vals

  !*****************************************************************************

END MODULE co2calc
