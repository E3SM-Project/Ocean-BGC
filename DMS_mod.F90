!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module DMS_mod

!BOP
! !MODULE: DMS_mod
!
! !DESCRIPTION:
!
 !------------------------------------------------------------------------------
 !   Multispecies trace gas processing routine based on Chu, Elliott, Maltrud,
 !   Erickson and company papers > 2000, on references therein,
 !   and also on assorted supplementary materials as indicated in comments.
 !   Designed for insertion into POP and CCSM along with DML ecosys_mod
 !   Trace gases defined roughly as those of concentration order nanomolar or less
 !   This excludes CO2 and O2 which are handled inside NCAR ecodynamics
 !------------------------------------------------------------------------------
 !------------------------------------------------------------------------------
 !   variables/subroutines/functions used from other modules
 !   The following are called upon extensively in tracegas, and so appear at
 !   the module level. The use statements for variables that are only needed
 !   locally are dealt with at the module subprogram level.
 !------------------------------------------------------------------------------

! !REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! !USES:


! !INPUT PARAMETERS:
  !-----------------------------------------------------------------------------
  !   include tracegas and ecosystem parameters
  !   all variables from these modules have a parm_ prefix
  !-----------------------------------------------------------------------------

   use DMS_parms

   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      DMS_tracer_cnt,     &
      DMS_init,           &
      DMS_SurfaceFluxes,  &
      DMS_SourceSink

!-----------------------------------------------------------------------
!  module variables
!-----------------------------------------------------------------------

   integer (DMS_i4), parameter :: &
      DMS_tracer_cnt = 14

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: DMS_init
! !INTERFACE:

 subroutine DMS_init(DMS_indices)

! !DESCRIPTION:
!  Initialize tracegas tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(DMS_indices_type), intent(inout) :: DMS_indices

! !INPUT/OUTPUT PARAMETERS:

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   DMS_indices%short_name(DMS_indices%dms_ind)='DMS'
   DMS_indices%long_name(DMS_indices%dms_ind)='DiMethyl Sulfide'

   DMS_indices%short_name(DMS_indices%dmsp_ind)='DMSP'
   DMS_indices%long_name(DMS_indices%dmsp_ind)='Dimethylsulfoniopropionate'

   DMS_indices%short_name(DMS_indices%no3_ind)='NO3'
   DMS_indices%long_name(DMS_indices%no3_ind)='Dissolved Inorganic Nitrate'

   DMS_indices%short_name(DMS_indices%doc_ind)='DOC'
   DMS_indices%long_name(DMS_indices%doc_ind)='Dissolved Organic Carbon'

   DMS_indices%short_name(DMS_indices%zooC_ind)='zooC'
   DMS_indices%long_name(DMS_indices%zooC_ind)='Zooplankton Carbon'

   DMS_indices%short_name(DMS_indices%spChl_ind) = 'spChl'
   DMS_indices%long_name(DMS_indices%spChl_ind)  = ' Small Phytoplankton Chlorophyll'

   DMS_indices%short_name(DMS_indices%diatChl_ind) = 'diatChl'
   DMS_indices%long_name(DMS_indices%diatChl_ind)  = ' Diatom Chlorophyll'

   DMS_indices%short_name(DMS_indices%diazChl_ind) = 'diazChl'
   DMS_indices%long_name(DMS_indices%diazChl_ind)  = ' Diazotroph Chlorophyll'

   DMS_indices%short_name(DMS_indices%phaeoChl_ind) = 'phaeoChl'
   DMS_indices%long_name(DMS_indices%phaeoChl_ind)  = 'Phaeocystis Chlorophyll'

   DMS_indices%short_name(DMS_indices%spC_ind) = 'spC'
   DMS_indices%long_name(DMS_indices%spC_ind)  = ' Small Phytoplankton Carbon'

   DMS_indices%short_name(DMS_indices%diatC_ind) = 'diatC'
   DMS_indices%long_name(DMS_indices%diatC_ind)  = ' Diatom Carbon'

   DMS_indices%short_name(DMS_indices%diazC_ind) = 'diazC'
   DMS_indices%long_name(DMS_indices%diazC_ind)  = ' Diazotroph Carbon'

   DMS_indices%short_name(DMS_indices%phaeoC_ind) = 'phaeoC'
   DMS_indices%long_name(DMS_indices%phaeoC_ind)  = 'Phaeocystis Carbon'

   DMS_indices%short_name(DMS_indices%spCaCO3_ind) = 'spCaCO3'
   DMS_indices%long_name(DMS_indices%spCaCO3_ind)  = ' Small Phytoplankton Calcium Carbonate'

   DMS_indices%units(:)                       = 'mmol/m^3'

!-----------------------------------------------------------------------
!EOC

 end subroutine DMS_init

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: DMS_SourceSink
! !INTERFACE:

 subroutine DMS_SourceSink(DMS_indices, DMS_input, DMS_forcing, DMS_output, DMS_diagnostic_fields, &
                           numLevelsMax, numColumnsMax, numColumns)
                           

! !DESCRIPTION:
!  Compute time derivatives for tracegas state variables
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(DMS_indices_type),     intent(in ) :: DMS_indices
  type(DMS_input_type),       intent(in ) :: DMS_input
  type(DMS_forcing_type),     intent(in ) :: DMS_forcing

  integer (DMS_i4) :: numLevelsMax, numColumnsMax, numColumns

! !OUTPUT PARAMETERS:

  type(DMS_output_type),      intent(inout) :: DMS_output
  type(DMS_diagnostics_type), intent(inout) :: DMS_diagnostic_fields

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    real(DMS_r8) :: &
         totalChl,       & ! total chlorophyll from all phytoplankton
         PAR_out,        & ! photosynthetically available radiation (W/m^2)
         PAR_in,         & ! photosynthetically available radiation (W/m^2)
         KPARdz,         & ! PAR adsorption coefficient (non-dim)
         PAR_avg,        & ! average PAR over mixed layer depth (W/m^2)
         UV_out,         & ! generic UV radiation (W/m^2)
         UV_in,          & ! generic UV radiation (W/m^2)
         KUVdz,          & ! UV adsorption coefficient (non-dim)
         UV_avg            ! average UV (W/m^2)

    !--------------------------------------------------------------------------------
    !   DMS and DMSP are the dissolved forms carried as tracers.
    !   The community would refer to them as DMS/Pd.
    !   The analytically particulate form DMSPp is for the moment
    !   defined as computationally local in order to minimize expense.
    !   Averaging and output will be added as necessary.
    !--------------------------------------------------------------------------------

    real(DMS_r8), allocatable, dimension(:,:) :: & !mmol molecule/m^3
         DMS_loc,        & ! local copy of model DMS
         DMSP_loc,       & ! local copy of model DMSP
         NO3_loc,        & ! local copy of model NO3
         DOC_loc,        & ! local copy of model DOC
         zooC_loc,       & ! local copy of model zooC
         spC_loc,        & ! local copy of model spC
         spCaCO3_loc,    & ! local copy of model spCaCO3
         diatC_loc,      & ! local copy of model diatC
         diazC_loc,      & ! local copy of model diazC
         phaeoC_loc,     & ! local copy of model phaeoC
         spChl_loc,      & ! local copy of model spChl
         diatChl_loc,    & ! local copy of model diatChl
         diazChl_loc,    & ! local copy of model diazChl
         phaeoChl_loc      ! local copy of model phaeoChl

    real(DMS_r8) :: &
         Fcocco            ! Fraction of small carbon associated with calcite

    !--------------------------------------------------------------------------------
    !   Conversion to nitrogen currency as the base for initial sulfur simulations.
    !   Links will most often be made to models such as Vallina et al. 2008,
    !   which tend to track organisms as N.
    !   Thus the CCSM standard carbon quantities are transposed.
    !   The trace gas system really begins here.
    !   Note that the DML phyto-classes are often subdivided,
    !   e.g. into specialists such as phaeocystis or nonS producers cyanobacteria.
    !   Ordering of the organisms in the code is intended to reflect
    !   prioritization of this decomposition process.
    !   E.g. coccolithophorids are partitioned prior to the cyano
    !   but after phaeocystis due to its traditional local dominance.
    !   Diatom, diazotroph and zooplanktonic bins are already sulfur appropriate
    !   and so are not directly altered here.
    !   Heterotrophic, recycling bacteria will be decoupled from phaeocystis
    !   which is colonial and generates antibiotics.
    !--------------------------------------------------------------------------------

    real(DMS_r8) :: & !all mmol N/m^3
         diatN_loc,       & ! analog to diatC_loc
         phaeoN_loc,      & ! nitrogen associated with phaeocystis
         coccoN_loc,      & ! nitrogen associated with coccolithophores
         cyanoN_loc,      & ! nitrogen associated with cyanobacteria
         eukarN_loc,      & ! analog to spC_loc but remove specialists
         diazN_loc,       & ! analog to diazC_loc
         phytoN_loc,      & ! total phytoplanktonic nitrogen excluding phaeo
         zooN_loc           ! analog to zooC_loc

    !--------------------------------------------------------------------------------
    !   Further conversion to sulfur quantities.
    !   These may be thought of as particulate DMSP by class
    !   and so may be zero as for example in the case of prokaryotes.
    !   A simple renaming is included at the end of the list
    !   to align with DMS community convention.
    !   Phaeocystis is pulled out of the total again,
    !   but this time because it is not subject to grazing.
    !   An even, steady state leakage from colonies is assumed.
    !--------------------------------------------------------------------------------

    real(DMS_r8) :: & !all mmol S/m^3
         diatS_loc,       & ! DMSP in diatoms
         phaeoS_loc,      & ! DMSP in phaeocystis
         coccoS_loc,      & ! DMSP in coccolithophores
         cyanoS_loc,      & ! DMSP in cyanobacteria
         eukarS_loc,      & ! DMSP in smalls remaining
         diazS_loc,       & ! DMSP in diazotrophs
         phytoS_loc,      & ! total phytoplanktonic sulfur excluding phaeo
         zooS_loc           ! total zooplanktonic sulfur

    !--------------------------------------------------------------------------------
    !   The mechanism traces to reviews such as Kiene et al. 2000 and Simo 2004
    !   supplemented by detail taken from low dimensionality models,
    !   examples being Lefevre et al. 2002, Vallina et al. 2008, Toole et al. 2008.
    !   Best descriptions to date appear in our reports on piston velocity testing.
    !   Channels are similar but not identical to those incorporated into
    !   the first (SciDAC) coupled ocean-atmosphere chemistry runs in CCSM.
    !   Key points include:
    !   Phytoplanktonic intracellular contents consistent with Stefels 2000.
    !   Zooplanktonic dependence built into the release constant k
    !   which relaxes coupling to chlorophyll and works toward a summer phase lag.
    !   Sunda et al. (2002) cell internal oxidant stress under high irradiance
    !   treated ad hoc by enhancing picoeukaryotic cell content
    !   in proportion to an arbitrary chlorophyll decrement.
    !   This situation arose because ultraviolet calculations were running late
    !   relative to the recent CODiM intercomparison.
    !   Real UV penetration will soon be computed per Chu et al. CO simulations.
    !   Direct summer exudation remains a speculative explanation for summer peaks.
    !   It may be thought of as closely related to the decrement approach
    !   but is not actually simulated here in a direct sense.
    !   DMS yield from bacterial processing of DMSP is given a T dependence.
    !   The concept is that microbial sulfur demand will be higher
    !   in cold, nutrient rich waters. See Kiene et al. 2000 for concepts involved.
    !   True sulfur utilization must clearly be incorporated very soon.
    !   This can begin with application of a cycling time to diagnosed densities,
    !   but may eventually entail addition of a dynamic bacterial module.
    !   Microbial consumption is rendered 2nd order per the above reviews,
    !   but in particular based on the equatorial data of Kiene and Bates 1990.
    !   Densities depend on free (non-colonial) phytoplankton distributions.
    !   Injection scale is a dial for dealing simultaneously with uncertainties in
    !   disruption rate and intracellular DMSP content.
    !   Optimization issues are thus focused upon average release.
    !   This collapses their dimensionality considerably.
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   These kinetically important quantities have geographic dependence.
    !--------------------------------------------------------------------------------

    real(DMS_r8) :: &
         k_S_p,        & ! 1st order constant for DMSP release from phyto (1/sec)
         yield,        & ! fraction of DMSP converion to DMS
         B_diagnosed,  & ! local bacterial density (mmol N/m^3)
         j_dms           ! overall dms photolysis rate (1/sec)

    real(DMS_r8) :: &
         T_ind             ! upward temperature index

    real(DMS_r8) :: &
         Cocco_frac,       & ! local fraction of smalls as coccos
         Cyano_frac,       & ! local fraction of smalls as cyanos
         Eukar_frac          ! local fraction of remaining smalls

    real(DMS_r8) :: &
         Sp_dec            ! small phytoplanktonic decrement

    real(DMS_r8) :: &
         Stress_fac        ! local up regulation

    !--------------------------------------------------------------------------------
    !   Surface temperature fixed moving down the column
    !   for several purposes within the fuzzy and binary logic 
    !--------------------------------------------------------------------------------

    real(DMS_r8) :: &
         SST_loc           ! local surface temperature (C)

    !--------------------------------------------------------------------------------
    !   Average cell contents have geographic dependence for zooplantkon
    !   because they are weighted averages over food items.
    !--------------------------------------------------------------------------------

    real(DMS_r8) :: &
         Rs2n_zoo        ! S/N assuming weighted average of consumable content

    !--------------------------------------------------------------------------------
    !   Begin declaration of source sink terms.
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   By and large the structure of the mechanism is reflected.
    !   For example in the present configuration,
    !   most sulfur release from phytoplankton takes the form DMSP.
    !   In anticipation of increased mechanistic complexity
    !   certain likely future permutations are included, e.g. exudation.
    !   The concepts here are clarity and motivation.
    !   The ports are simply readied for later connections.
    !   A few omissions will be obvious.
    !   There are to date no indications that DMSP photolyzes in the column.
    !   Note that sea-air transfer is time split into other routines.
    !--------------------------------------------------------------------------------

    real(DMS_r8) :: &
         dms_s_exu,    & ! DMS source from exudation (mmol S/m^3/sec)
         dms_s_dmsp,   & ! DMS source by conversion of DMSP (mmol S/m^3/sec)
         dms_s           ! DMS source total (mmol S/m^3/sec)

    real(DMS_r8) :: &
         dms_r_B,      & ! DMS removal by bacteria (mmol S/m^3/sec)
         dms_r_phot,   & ! DMS removal by photolysis (mmol S/m^3/sec)
         dms_r_bkgnd,  & ! DMS removal low level in thermocline
         dms_r           ! DMS removal total (mmol S/m^3/sec)

    real(DMS_r8) :: &
         dmsp_s_phaeo, & ! DMSP source from phaeo (mmol S/m^3/sec)
         dmsp_s_nonphaeo,   & ! DMSP source from other phytoplankton (mmol S/m^3/sec)
         dmsp_s_zoo,   & ! DMSP source from zooplankton (mmol S/m^3/sec)
         dmsp_s          ! DMSP source total (mmol S/m^3/sec)

   real(DMS_r8) :: &
         dmsp_r_B,       & ! DMSP removal by bacteria (mmol S/m^3/sec)
         dmsp_r_bkgnd,   & ! DMSP removal low level in thermocline
         dmsp_r            ! DMSP removal total (mmol S/m^3/sec)

   real(DMS_r8) :: &
      work

   integer(DMS_i4) ::  &
      column, kmax, k

   integer (DMS_i4) :: &
      dms_ind,         & ! DMS index
      dmsp_ind,        & ! DMSP index
      no3_ind,         & ! DMSP index
      doc_ind,         & ! DMSP index
      zooC_ind,        & ! DMSP index
      spC_ind,         & ! DMSP index
      spCaCO3_ind,     & ! DMSP index
      diatC_ind,       & ! DMSP index
      diazC_ind,       & ! DMSP index
      phaeoC_ind,      & ! DMSP index
      spChl_ind,       & ! DMSP index
      diatChl_ind,     & ! DMSP index
      diazChl_ind,     & ! DMSP index
      phaeoChl_ind       ! DMSP index

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  initialize  all tendencies to zero
!-----------------------------------------------------------------------

   DMS_output%DMS_tendencies = 0.0_DMS_r8

!-----------------------------------------------------------------------
!  allocate local copies of tracers
!-----------------------------------------------------------------------

   allocate(DMS_loc(numLevelsMax,numColumns))
   allocate(DMSP_loc(numLevelsMax,numColumns))
   allocate(NO3_loc(numLevelsMax,numColumns))
   allocate(DOC_loc(numLevelsMax,numColumns))
   allocate(spC_loc(numLevelsMax,numColumns))
   allocate(spCaCO3_loc(numLevelsMax,numColumns))
   allocate(diatC_loc(numLevelsMax,numColumns))
   allocate(diazC_loc(numLevelsMax,numColumns))
   allocate(phaeoC_loc(numLevelsMax,numColumns))
   allocate(spChl_loc(numLevelsMax,numColumns))
   allocate(diatChl_loc(numLevelsMax,numColumns))
   allocate(diazChl_loc(numLevelsMax,numColumns))
   allocate(phaeoChl_loc(numLevelsMax,numColumns))
   allocate(zooC_loc(numLevelsMax,numColumns))

!-----------------------------------------------------------------------
!  assign indices.  this is not necessary but results in fewer
!    differences between original and new code.
!-----------------------------------------------------------------------

   no3_ind        = DMS_indices%no3_ind
   doc_ind        = DMS_indices%doc_ind
   zooC_ind       = DMS_indices%zooC_ind
   spC_ind        = DMS_indices%spC_ind
   diatC_ind      = DMS_indices%diatC_ind
   diazC_ind      = DMS_indices%diazC_ind
   phaeoC_ind     = DMS_indices%phaeoC_ind
   spChl_ind      = DMS_indices%spChl_ind
   diatChl_ind    = DMS_indices%diatChl_ind
   diazChl_ind    = DMS_indices%diazChl_ind
   phaeoChl_ind   = DMS_indices%phaeoChl_ind
   spCaCO3_ind    = DMS_indices%spCaCO3_ind

   dms_ind        = DMS_indices%dms_ind
   dmsp_ind       = DMS_indices%dmsp_ind

!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   setup_loop: do column = 1, numColumns

   kmax = DMS_input%number_of_active_levels(column)
   if (kmax < 1) cycle setup_loop

   do k = 1, kmax

    !---------------------------------------------------------------------------
    !   create local copies of requisite ecotracers
    !   treat negative values as zero and apply mask to locals
    !---------------------------------------------------------------------------

   NO3_loc(k,column)       = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,no3_ind))
   DOC_loc(k,column)       = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,doc_ind))
   zooC_loc(k,column)      = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,zooC_ind))
   spC_loc(k,column)       = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,spC_ind))
   diatC_loc(k,column)     = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,diatC_ind))
   diazC_loc(k,column)     = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,diazC_ind))
   phaeoC_loc(k,column)    = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,phaeoC_ind))
   spChl_loc(k,column)     = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,spChl_ind))
   diatChl_loc(k,column)   = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,diatChl_ind))
   diazChl_loc(k,column)   = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,diazChl_ind))
   phaeoChl_loc(k,column)  = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,phaeoChl_ind))
   spCaCO3_loc(k,column)   = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,spCaCO3_ind))

   DMS_loc(k,column)       = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,dms_ind))
   DMSP_loc(k,column)      = max(0.0_DMS_r8, DMS_input%DMS_tracers(k,column,dmsp_ind))

   end do  !  end of setup k loop

   enddo  setup_loop  !  end of setup column loop

!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   column_loop: do column = 1, numColumns

   kmax = DMS_input%number_of_active_levels(column)
   if (kmax < 1) cycle column_loop

!-----------------------------------------------------------------------
!  various k==1 initializations
!-----------------------------------------------------------------------

   SST_loc = DMS_forcing%SST(column)

   PAR_out = max (0.0_DMS_r8, DMS_forcing%ShortWaveFlux_surface(column))
   PAR_out = PAR_out*f_qsw_par_DMS

! UV is 1% of PAR
   UV_out = PAR_out*0.01_DMS_r8

!-----------------------------------------------------------------------
!  loop over levels
!-----------------------------------------------------------------------

   do k = 1, kmax

    !---------------------------------------------------------------------------
    !  Baseline phytoplanktonic sulfur release rate constant is here adjusted
    !  by local zooplanktonic densities normalized to an average.
    !  This has the effect of partially decoupling sulfur from chlorophyll.
    !  Carbon is used in this case since conversions would simply cancel.
    !  Mortality modulates the decoupling but may be zeroed.
    !  Note that phaeocystis release is treated as an independent leak
    !  adjusted in fact to control grazing.
    !  This is normally given the base rate but could be as slow as senescence.
    !---------------------------------------------------------------------------

   k_S_p = k_S_p_base * (mort + (zooC_loc(k,column)/0.3_DMS_r8))

   UV_in  = UV_out
    
   KUVdz = (0.01e-2_DMS_r8 * DOC_loc(k,column) + 0.04e-4_DMS_r8) * DMS_input%cell_thickness(k,column)

   UV_out = UV_in * exp(-KUVdz)
   UV_avg = UV_in * (1.0_DMS_r8 - exp(-KUVdz)) / KUVdz
    
   PAR_in = PAR_out
   totalChl = spChl_loc(k,Column) + diatChl_loc(k,Column) + diazChl_loc(k,Column)   &
            + phaeoChl_loc(k,Column)
   work = max(totalChl, 0.02_DMS_r8)
   if (work < 0.13224_DMS_r8) then
      KPARdz = 0.000919_DMS_r8*(work**0.3536_DMS_r8)
   else
      KPARdz = 0.001131_DMS_r8*(work**0.4562_DMS_r8)
   end if

   KPARdz = KPARdz * DMS_input%cell_thickness(k,column)

   PAR_out = PAR_in * exp(-KPARdz)
   PAR_avg = PAR_in * (1.0_DMS_r8 - exp(-KPARdz)) / KPARdz

    !---------------------------------------------------------------------------
    !   Momentarily we follow Gabric et al. 1993 who speculate implicitly that
    !   DMS photolytic wavelengths attenuate as PAR.
    !   This is confirmed in part by the review Mopper and Kieber 2002,
    !   But modernization will be required per Toole et al. 2004 and related.
    !   Photolysis per unit intensity is based by Gabric on
    !   Brimblecombe and Shooter 1986 results for wavelengths unknown.
    !---------------------------------------------------------------------------

    j_dms = j_dms_perI * PAR_avg

    !---------------------------------------------------------------------------
    !   Fcocco is the fraction of sp organic matter in coccolithophores
    !   as imported from the driver DML ecology.
    !   In the S model this takes precedence over all but phaeocystis.
    !---------------------------------------------------------------------------

    Fcocco = spCaCO3_loc(k,column) / (spC_loc(k,column) + epsC)
    if (Fcocco > 0.4_DMS_r8) Fcocco = 0.4_DMS_r8

    Cocco_frac = Fcocco

    !--------------------------------------------------------------------------
    !   Fuzzy logical segregation into relevant subclasses
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !   Further distribution of the DML small phytoplanktonic biomass
    !   is now undertaken per requirements of the sulfur cycle.
    !   Phaeocystis is considered dominant and ultimately trumps other autotrophs.
    !--------------------------------------------------------------------------

    T_ind  = (SST_loc - T_lo)/(T_hi - T_lo)

    if (T_ind <= 0.0_DMS_r8) T_ind = 0.0_DMS_r8
    if (T_ind >= 1.0_DMS_r8) T_ind = 1.0_DMS_r8

    Cyano_frac  = (T_ind * (Max_cyano_frac - Min_cyano_frac)) + Min_cyano_frac
    Cyano_frac  = (1.0_DMS_r8 - Cocco_frac) * Cyano_frac

    Eukar_frac = 1.0_DMS_r8 - Cocco_frac - Cyano_frac

    !--------------------------------------------------------------------------
    !  Convert to nitrogen currency plus distribute into new classes
    !--------------------------------------------------------------------------

    diatN_loc    = R * diatC_loc(k,column)
    phaeoN_loc   = R * phaeoC_loc(k,column)
    coccoN_loc   = Cocco_frac * R * spC_loc(k,column)
    cyanoN_loc   = Cyano_frac * R * spC_loc(k,column)
    eukarN_loc   = Eukar_frac * R * spC_loc(k,column)
    diazN_loc    = R*diazC_loc(k,column)
    zooN_loc     = R*zooC_loc(k,column)

    !--------------------------------------------------------------------------
    !  Collect noncolonial nitrogen
    !  (swang) since phaeocystis is grazed in BEC, include it in grazing-dmsp production too
    !--------------------------------------------------------------------------

    phytoN_loc = diatN_loc + coccoN_loc + cyanoN_loc + eukarN_loc + diazN_loc  &
    			+ phaeoN_loc 

    !--------------------------------------------------------------------------
    !  Gyre oxidant stress now exerts upregulation per chlorophyll decrement.
    !  Functional form maintains round figure parameters and powers,
    !  but was essentially determined by offline, postprocessing trial and error.
    !  Real physiology must be incorporated in short order.
    !--------------------------------------------------------------------------

    Sp_dec = (Sp_ref - spChl_loc(k,column))/Sp_ref

    if (Sp_dec <= 0.0_DMS_r8) Sp_dec = 0.0_DMS_r8
    if (Sp_dec >= 1.0_DMS_r8) Sp_dec = 1.0_DMS_r8

    Stress_fac = 1.0_DMS_r8 + Stress_mult* Sp_dec * Sp_dec

    if (Stress_fac >= 10.0_DMS_r8) Stress_fac = 10.0_DMS_r8

    !--------------------------------------------------------------------------
    !  Yields for bacterial conversion of DMSP also determined by fuzzy logic.
    !  Since Phaeocystis is colonial and produces antimicrobials,
    !  its habitat constitutes an exception at very low S demand.
	!  (swang) assume an optimal temperautre for cryoprotention 
    !--------------------------------------------------------------------------

    yield = (T_ind * (Max_yld - Min_yld)) + Min_yld

    if (SST_loc < T_cryo_hi .and. SST_loc > T_cryo_lo) yield = 0.5_DMS_r8 
    if (SST_loc < -1.0_DMS_r8) yield = 0.25_DMS_r8 

    !--------------------------------------------------------------------------
    !  Phytoplanktonic sulfur content determined.
    !  This is where oxidant stress is in fact applied.
    !--------------------------------------------------------------------------

    diatS_loc  = Rs2n_diat  * diatN_loc
    phaeoS_loc = Rs2n_phaeo * phaeoN_loc
    coccoS_loc = Rs2n_cocco * coccoN_loc
    cyanoS_loc = Rs2n_cyano * cyanoN_loc
    eukarS_loc = Rs2n_eukar * eukarN_loc * Stress_fac
    diazS_loc  = Rs2n_diaz   * diazN_loc

    !--------------------------------------------------------------------------
    !  Collect noncolonials
    ! (swang) assume only a fraction (40%) of phaeo S contribute to dsmp release when grazed
    !--------------------------------------------------------------------------

    phytoS_loc = diatS_loc + coccoS_loc + cyanoS_loc + eukarS_loc + diazS_loc &
                + G_phaeo_S * phaeoS_loc  

    !--------------------------------------------------------------------------
    !  Weight the zooplanktonic reduced sulfur content per the various
    !  phytoplankton types available unprotected in the column as food.
    !  Observe that plant nitrogen may total zero under some circumstances
    !  and in this case it is assumed that all concentrations are equally small.
    !  Note that phaeocystis has been segregated entirely.
    !  In colonial form it exhibits multiple grazing inhibition strategies.
    !--------------------------------------------------------------------------

    if (phytoN_loc > 0.0_DMS_r8) then
    Rs2n_zoo = (Rs2n_diat   * diatN_loc  + &
    			 G_phaeo_S * Rs2n_phaeo   * phaeoN_loc  + &
                 Rs2n_cocco * coccoN_loc + &
                 Rs2n_cyano * cyanoN_loc + &
                 Rs2n_eukar * eukarN_loc * Stress_fac + &
                 Rs2n_diaz   * diazN_loc)/phytoN_loc
    else
! (swang) include 1group phaeo 
    Rs2n_zoo = (Rs2n_diat + Rs2n_cocco + Rs2n_cyano + Rs2n_eukar + Rs2n_diaz &
    			+ Rs2n_phaeo)/6.0
    end if

    zooS_loc = Rs2n_zoo * zooN_loc

    !---------------------------------------------------------------------------
    !   Given a nitrogen distribution diagnose B which is bacterial N.
    !   Phaeocystis excluded because it is colonial and exudes antibiotics.
    !   Given a turnover time B could serve as the basis for S demand calcs.
    !   However our tendency is to move directly to bacterial populations.
    !   Second order DMS removal is based on reports of biotic uptake presented
    !   by Kiene and Bates 1991 for the eastern equatorial Pacific.
    !---------------------------------------------------------------------------

    B_diagnosed = B_preexp*(phytoN_loc**B_exp)

    !-------------------------------------------------------------------------
    !   Construction of kinetic terms for the sulfur cycle
    !-------------------------------------------------------------------------

    dms_s_dmsp   = yield * k_conv * DMSP_loc(k,column)
    dms_s        = dms_s_dmsp

    dms_r_B      = k_S_B * B_diagnosed * DMS_loc(k,column)
    dms_r_phot   = j_dms * DMS_loc(k,column)
    dms_r_bkgnd  = k_bkgnd * DMS_loc(k,column)
    dms_r        = dms_r_B + dms_r_phot + dms_r_bkgnd

    dmsp_s_phaeo = inject_scale * k_S_p_base * phaeoS_loc 
    dmsp_s_nonphaeo   = inject_scale * k_S_p * phytoS_loc
    dmsp_s_zoo   = inject_scale * k_S_z * zooS_loc
    dmsp_s       = dmsp_s_phaeo + dmsp_s_nonphaeo + dmsp_s_zoo

    dmsp_r_B     = k_conv * DMSP_loc(k,column)
    dmsp_r_bkgnd = k_bkgnd * DMSP_loc(k,column)
    dmsp_r       = dmsp_r_B + dmsp_r_bkgnd

    DMS_output%DMS_tendencies(k,column,dms_ind)  = dms_s - dms_r
    DMS_output%DMS_tendencies(k,column,dmsp_ind) = dmsp_s - dmsp_r
! all other tendencies were initialized to 0

! DMS source terms
       DMS_diagnostic_fields%diag_DMS_S_DMSP(k,column) = dms_s_dmsp
       DMS_diagnostic_fields%diag_DMS_S_TOTAL(k,column) = dms_s

! DMS removal terms
       DMS_diagnostic_fields%diag_DMS_R_B(k,column) = dms_r_B
       DMS_diagnostic_fields%diag_DMS_R_PHOT(k,column) = dms_r_phot
       DMS_diagnostic_fields%diag_DMS_R_BKGND(k,column) = dms_r_bkgnd
       DMS_diagnostic_fields%diag_DMS_R_TOTAL(k,column) = dms_r

! DMSP source terms
       DMS_diagnostic_fields%diag_DMSP_S_PHAEO(k,column) = dmsp_s_phaeo
       DMS_diagnostic_fields%diag_DMSP_S_NONPHAEO(k,column) = dmsp_s_nonphaeo
       DMS_diagnostic_fields%diag_DMSP_S_ZOO(k,column) = dmsp_s_zoo
       DMS_diagnostic_fields%diag_DMSP_S_TOTAL(k,column) = dmsp_s

! DMSP removal terms
       DMS_diagnostic_fields%diag_DMSP_R_B(k,column) = dmsp_r_B
       DMS_diagnostic_fields%diag_DMSP_R_BKGND(k,column) = dmsp_r_bkgnd
       DMS_diagnostic_fields%diag_DMSP_R_TOTAL(k,column) = dmsp_r

! fractional compositions
       DMS_diagnostic_fields%diag_Cyano_frac(k,column) = Cyano_frac
       DMS_diagnostic_fields%diag_Cocco_frac(k,column) = Cocco_frac
       DMS_diagnostic_fields%diag_Eukar_frac(k,column) = Eukar_frac

! sulfur content
       DMS_diagnostic_fields%diag_diatS(k,column) = diatS_loc
       DMS_diagnostic_fields%diag_diatN(k,column) = diatN_loc
       DMS_diagnostic_fields%diag_phytoN(k,column) = phytoN_loc
       DMS_diagnostic_fields%diag_coccoS(k,column) = coccoS_loc
       DMS_diagnostic_fields%diag_cyanoS(k,column) = cyanoS_loc
       DMS_diagnostic_fields%diag_eukarS(k,column) = eukarS_loc
       DMS_diagnostic_fields%diag_diazS(k,column) = diazS_loc
       DMS_diagnostic_fields%diag_phaeoS(k,column) = phaeoS_loc
       DMS_diagnostic_fields%diag_zooS(k,column) = zooS_loc

! other
       DMS_diagnostic_fields%diag_zooCC(k,column) = zooC_loc(k,column)
       DMS_diagnostic_fields%diag_RSNzoo(k,column) = Rs2n_zoo

   enddo ! k loop

   enddo column_loop ! i loop

!-----------------------------------------------------------------------
!EOC

 end subroutine DMS_SourceSink

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: DMS_SurfaceFluxes
! !INTERFACE:

 subroutine DMS_SurfaceFluxes(DMS_indices, DMS_input, DMS_forcing,   &
                              DMS_flux_diagnostic_fields, &
                              numColumnsMax, numColumns)

! !DESCRIPTION:
!  Compute surface fluxes for tracegas tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(DMS_indices_type), intent(in )   :: DMS_indices
  type(DMS_input_type),   intent(in )   :: DMS_input
  type(DMS_forcing_type), intent(inout) :: DMS_forcing

  integer (DMS_i4) :: numColumnsMax, numColumns

! !OUTPUT PARAMETERS:

  type(DMS_flux_diagnostics_type), intent(inout) :: DMS_flux_diagnostic_fields

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (DMS_r8) :: &
      seaSurfaceTemp,   &
      seaSurfaceSalt,   &
      seaSurfaceDMS,   &
      xkw,              &
      xkw_ice,          &
      SCHMIDT_DMS,       &
      DMSSAT_1atm,       &
      pistonVel_DMS,     &
      DMSSAT,           &
      WIND_SPEED,       &
      FW92,         & ! apportionment to Wanninkhof 1992
      FLM86,        & ! apportionment to Liss and Merlivat 1986
      XKW_W92,      & ! the Wanninkhof limit
      XKW_LM86      ! the Liss and Merlivat limit

   real (DMS_r8) :: scalar_temp

   integer (DMS_i4) :: &
      column            ! index for looping over columns

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------

    real(DMS_r8), parameter :: &
         a  = 0.31_DMS_r8,     &    ! W92
         e1 = 0.17_DMS_r8,     &    ! LM86 from here 
         e2 = 2.85_DMS_r8,     &
         e3 = 0.612_DMS_r8,    &    
         e4 = 5.9_DMS_r8,      &    
         e5 = 26.79_DMS_r8,    &
         e6 = 0.612_DMS_r8           

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  compute DMS flux
!-----------------------------------------------------------------------

   if (DMS_forcing%lcalc_DMS_gas_flux) then

!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   column_loop: do column = 1, numColumns

      seaSurfaceDMS  = max(0.0_DMS_r8, DMS_input%DMS_tracers(1,column,DMS_indices%dms_ind))
      seaSurfaceTemp = DMS_forcing%SST(column)
      seaSurfaceSalt = DMS_forcing%SSS(column)

      if (DMS_forcing%iceFraction(column) < 0.0_DMS_r8) DMS_forcing%iceFraction(column) = 0.0_DMS_r8
      if (DMS_forcing%iceFraction(column) > 1.0_DMS_r8) DMS_forcing%iceFraction(column) = 1.0_DMS_r8

         SCHMIDT_DMS = SCHMIDT_DMS_singleValue(seaSurfaceTemp)
!        SCHMIDT_DMS = min(SCHMIDT_DMS,1.e20)
!        SCHMIDT_DMS = max(SCHMIDT_DMS,1.)

! convert to m/s
      WIND_SPEED = sqrt(abs(DMS_forcing%windSpeedSquared10m(column)))*0.01_DMS_r8

      XKW_W92 = &
         a*((660.0_DMS_r8/SCHMIDT_DMS)**0.500_DMS_r8)*WIND_SPEED*WIND_SPEED
      XKW_LM86 = &
         e2*((600.0_DMS_r8/SCHMIDT_DMS)**0.500_DMS_r8)*(WIND_SPEED - 3.6_DMS_r8) &
       + e3*((600.0_DMS_r8/SCHMIDT_DMS)**0.667_DMS_r8)

      if (WIND_SPEED < 3.6_DMS_r8) xkw = XKW_W92
      if ((WIND_SPEED >= 3.6_DMS_r8) .and. (WIND_SPEED < 5.6_DMS_r8)) then
         FLM86 = 0.5_DMS_r8*(WIND_SPEED - 3.6_DMS_r8)
         FW92  = 1.0_DMS_r8 - FLM86
         xkw = FW92*XKW_W92 + FLM86*XKW_LM86
      end if
      if (WIND_SPEED >= 5.6_DMS_r8) xkw = XKW_LM86

      xkw = xkw/3600.0_DMS_r8 ! conversion to cm/s
      xkw_ice = (1.0_DMS_r8 - DMS_forcing%iceFraction(column)) * xkw

         DMSSAT_1atm = DMSSAT_singleValue(seaSurfaceTemp, seaSurfaceSalt)

         pistonVel_DMS = xkw_ice * SQRT(660.0_DMS_r8 / SCHMIDT_DMS)
         DMSSAT = DMS_forcing%surfacePressure(column) * DMSSAT_1atm
         DMS_forcing%netFlux(column, DMS_indices%dms_ind) = pistonVel_DMS * (DMSSAT - seaSurfaceDMS)
         DMS_forcing%netFlux(column, DMS_indices%dmsp_ind) = 0.0_DMS_r8

         DMS_flux_diagnostic_fields%diag_DMS_IFRAC(column) = DMS_forcing%iceFraction(column)
         DMS_flux_diagnostic_fields%diag_DMS_XKW(column) = xkw_ice
         DMS_flux_diagnostic_fields%diag_DMS_ATM_PRESS(column) = DMS_forcing%surfacePressure(column)
         DMS_flux_diagnostic_fields%diag_DMS_PV(column) = pistonVel_DMS
         DMS_flux_diagnostic_fields%diag_DMS_SCHMIDT(column) = SCHMIDT_DMS
         DMS_flux_diagnostic_fields%diag_DMS_SAT(column) = DMSSAT
         DMS_flux_diagnostic_fields%diag_DMS_SURF(column) = seaSurfaceDMS
         DMS_flux_diagnostic_fields%diag_DMS_WS(column) = WIND_SPEED

       enddo column_loop

    endif  ! lflux_gas_dms

!-----------------------------------------------------------------------
!EOC

 end subroutine DMS_SurfaceFluxes

!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_DMS
! !INTERFACE:

 function SCHMIDT_DMS_singleValue(SST)

! !DESCRIPTION:
!---------------------------------------------------------------------------
!   Compute Schmidt number in seawater as function of SST
!   where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!   ref : Kettle and Andreae 2000
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (DMS_r8), intent(in) :: SST

! !OUTPUT PARAMETERS:

   real (DMS_r8) :: SCHMIDT_DMS_singleValue

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   coefficients in expansion
    !---------------------------------------------------------------------------

    real(DMS_r8), parameter :: &
         a = 2674.0_DMS_r8, &
         b = 147.12_DMS_r8, &
         c = 3.726_DMS_r8, &
         d = 0.038_DMS_r8

!-----------------------------------------------------------------------

   SCHMIDT_DMS_singleValue = a + SST * (-b + SST * (c + SST * (-d)))

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_DMS_singleValue

!*****************************************************************************
!BOP
! !IROUTINE: DMSSAT_singleValue
! !INTERFACE:

 function DMSSAT_singleValue(SST, SSS)

! !DESCRIPTION:
!
!   Sat functions normally compute for a given molecule
!   a sea surface saturation concentration estimate at 1 atm total pressure
!   in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
!   in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
!   For DMS the assumption is made per Kettle and Andreae 2000
!   that the atmospheric concentration is negligible.
!   Henrys Law is thus not accounted here but passes for temperature and
!   salinity influences are preserved pending.
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   real (DMS_r8), intent(in) :: &
      SST, & ! sea surface temperature (C)
      SSS    ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

    real (DMS_r8) :: DMSSAT_singleValue

!EOP
!BOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!---------------------------------------------------------------------------
!   Units should lead to mmol/m^3 for saturation level
!---------------------------------------------------------------------------

   DMSSAT_singleValue = 0.0_DMS_r8

!-----------------------------------------------------------------------
!EOC

 end function DMSSAT_singleValue

!-----------------------------------------------------------------------

 end module DMS_mod
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
