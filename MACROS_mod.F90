!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module MACROS_mod

!BOP
! !MODULE: MACROS_mod
!
! !DESCRIPTION:
!
!------------------------------------------------------------------------------
!   Quick first cut at global macromolecules/surfactants, due to AMCE deadlines
!   Will retain much of the DMS template material on this pass,
!   potentially for later use.
!   For scientific and `technical details see Elliott et al. 2014
!                                             Burrows et al. 2014
!                                             Ogunro et al. 2015
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
  !   include macromolecular and ecosystem parameters
  !   all variables from these modules have a parm_ prefix
  !-----------------------------------------------------------------------------

   use MACROS_parms

   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      MACROS_tracer_cnt,     &
      MACROS_init,           &
      MACROS_SourceSink

!-----------------------------------------------------------------------
!  module variables
!-----------------------------------------------------------------------

   integer (MACROS_i4), parameter :: &
      MACROS_tracer_cnt = 8

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: MACROS_init
! !INTERFACE:

 subroutine MACROS_init(MACROS_indices)

! !DESCRIPTION:
!  Initialize macromolecules module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(MACROS_indices_type), intent(inout) :: MACROS_indices

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

   MACROS_indices%short_name(MACROS_indices%prot_ind)='PROT'
   MACROS_indices%long_name(MACROS_indices%prot_ind)='Proteins'

   MACROS_indices%short_name(MACROS_indices%poly_ind)='POLY'
   MACROS_indices%long_name(MACROS_indices%poly_ind)='Polysaccharides'

   MACROS_indices%short_name(MACROS_indices%lip_ind)='LIP'
   MACROS_indices%long_name(MACROS_indices%lip_ind)='Lipids'

   MACROS_indices%short_name(MACROS_indices%zooC_ind)='zooC'
   MACROS_indices%long_name(MACROS_indices%zooC_ind)='Zooplankton Carbon'

   MACROS_indices%short_name(MACROS_indices%spC_ind) = 'spC'
   MACROS_indices%long_name(MACROS_indices%spC_ind)  = ' Small Phytoplankton Carbon'

   MACROS_indices%short_name(MACROS_indices%diatC_ind) = 'diatC'
   MACROS_indices%long_name(MACROS_indices%diatC_ind)  = ' Diatom Carbon'

   MACROS_indices%short_name(MACROS_indices%diazC_ind) = 'diazC'
   MACROS_indices%long_name(MACROS_indices%diazC_ind)  = ' Diazotroph Carbon'

   MACROS_indices%short_name(MACROS_indices%phaeoC_ind) = 'phaeoC'
   MACROS_indices%long_name(MACROS_indices%phaeoC_ind)  = 'Phaeocystis Carbon'

   MACROS_indices%units(:)                       = 'mmol/m^3'

!-----------------------------------------------------------------------
!EOC

 end subroutine MACROS_init

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: MACROS_SourceSink
! !INTERFACE:

 subroutine MACROS_SourceSink(MACROS_indices, MACROS_input, MACROS_output, MACROS_diagnostic_fields, &
                           numLevelsMax, numColumnsMax, numColumns)
                           

! !DESCRIPTION:
!  Compute time derivatives for tracegas state variables
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  type(MACROS_indices_type),     intent(in ) :: MACROS_indices
  type(MACROS_input_type),       intent(in ) :: MACROS_input

  integer (MACROS_i4) :: numLevelsMax, numColumnsMax, numColumns

! !OUTPUT PARAMETERS:

  type(MACROS_output_type),      intent(inout) :: MACROS_output
  type(MACROS_diagnostics_type), intent(inout) :: MACROS_diagnostic_fields

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   The strategy from here is to substitute:
    !   Proteins for DMS, Polysaccacharides for DMPS, then adding lipids as 3rd C type.
    !   Maintain the ecosystem structure in this first cut.  
    !   For example, phaeocystis may exude polysaccharides.
    !   Dissolved organosulfur becomes dissolved and speciated biomacromolecular carbon. 
    !   Note that units are now carbon atoms.
    !--------------------------------------------------------------------------------

    real(MACROS_r8), allocatable, dimension(:,:) :: & !mmol carbon/m^3
         Prot_loc,       & ! local copy of model proteins
         Poly_loc,       & ! local copy of model polysaccharide
         Lip_loc,        & ! local copy of model lipids
         zooC_loc,       & ! local copy of model zooC
         spC_loc,        & ! local copy of model spC
         diatC_loc,      & ! local copy of model diatC
         diazC_loc,      & ! local copy of model diazC
         phaeoC_loc        ! local copy of model phaeoC

    !--------------------------------------------------------------------------------
    !   Macro production requires a sum of phytoplanktonic carbon.
    !   SW states phaeo is grazed so it is included.
    !   Early compile also shows column intermediates required.
    !   These are analogs of e.g. diatN_ from the DMS template
    !--------------------------------------------------------------------------------

    real(MACROS_r8) :: & !all mmol carbon/m^3
         spCk_loc,        & ! column test -ugly, there must be a better way 
         diatCk_loc,      & ! column test 
         diazCk_loc,      & ! column test 
         phaeoCk_loc,     & ! column test
         phytoC_loc         ! total phytoplanktonic carbon

    !--------------------------------------------------------------------------------
    !   Quantities for conversion to sulfur appeared here but are now irrelevant  
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   The biomacromolecular distribution mechanism is unique to the present model
    !   but that is in fact why/how it can be so simple.
    !   The primary reference will hopefully be Ogunro et al. 2015 in Biogeochemistry.
    !   But this one is still under revision.
    !   Very similar information is available in ocean sections of Burrows et al. 2014
    !   and physical chemical theory is laid out in Elliott et al. 2014.      
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   Only a few quantities have geographic dependence in early tests.
    !--------------------------------------------------------------------------------

    real(MACROS_r8) :: &
         k_C_p           ! 1st order constant for dissolved C release from phyto (1/sec)

    !--------------------------------------------------------------------------------
    !   Begin declaration of source sink terms.
    !--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------
    !   We will preview a macromolecular mechanism with the sort of kinetics
    !   imagined in Ogunro et al. 2015 but superseded there by a link to semis.
    !   Three separate time evolutions will be permitted then later normalized. 
    !   See O15 Table 1 for an explanation of the rate constants.
    !   Return and align with DOC as the representative for BEC semilabiles only
    !   if and when ACME deadlines can be met.     
    !--------------------------------------------------------------------------------

    real(MACROS_r8) :: &
         prot_s_disr,    & ! Protein source from cell disruption (mmol C/m^3/sec)
         poly_s_disr,    & ! Polysaccharide source from cell disruption (mmol C/m^3/sec)
         lip_s_disr        ! Lipid source from cell disruption (mmol C/m^3/sec)

    real(MACROS_r8) :: &
         prot_r_bac,    & ! Protein removal by heterotrophic bacteria (mmol C/m^3/sec)
         poly_r_bac,    & ! Polysacch removal by heterotrophic bacteria (mmol C/m^3/sec)
         lip_r_bac        ! Lipid removal by heterotrophic bacteria (mmol C/m^3/sec)

    real(MACROS_r8) :: &
         prot_s,        & ! Protein removal total (mmol C/m^3/sec)
         poly_s,        & ! Polysacch removal total (mmol C/m^3/sec)
         lip_s            ! Lipid removal total (mmol C/m^3/sec)

    real(MACROS_r8) :: &
         prot_r,        & ! Protein removal total (mmol C/m^3/sec)
         poly_r,        & ! Polysacch removal total (mmol C/m^3/sec)
         lip_r            ! Lipid removal total (mmol C/m^3/sec)

   real(MACROS_r8) :: &
      work

   integer(MACROS_i4) ::  &
      column, kmax, k

   integer (MACROS_i4) :: &
      prot_ind,        & ! MACROS index
      poly_ind,        & ! MACROS index
      lip_ind,         & ! MACROS index
      zooC_ind,        & ! MACROS index
      spC_ind,         & ! MACROS index
      diatC_ind,       & ! MACROS index
      diazC_ind,       & ! MACROS index
      phaeoC_ind         ! MACROS index

!-----------------------------------------------------------------------
!  initialize  all tendencies to zero
!-----------------------------------------------------------------------

   MACROS_output%MACROS_tendencies = 0.0_MACROS_r8

!-----------------------------------------------------------------------
!  allocate local copies of tracers
!-----------------------------------------------------------------------

   allocate(prot_loc(numLevelsMax,numColumns))
   allocate(poly_loc(numLevelsMax,numColumns))
   allocate(lip_loc(numLevelsMax,numColumns))
   allocate(spC_loc(numLevelsMax,numColumns))
   allocate(diatC_loc(numLevelsMax,numColumns))
   allocate(diazC_loc(numLevelsMax,numColumns))
   allocate(phaeoC_loc(numLevelsMax,numColumns))
   allocate(zooC_loc(numLevelsMax,numColumns))

!-----------------------------------------------------------------------
!  assign indices.  this is not necessary but results in fewer
!  differences between original and new code.
!-----------------------------------------------------------------------

   zooC_ind       = MACROS_indices%zooC_ind
   spC_ind        = MACROS_indices%spC_ind
   diatC_ind      = MACROS_indices%diatC_ind
   diazC_ind      = MACROS_indices%diazC_ind
   phaeoC_ind     = MACROS_indices%phaeoC_ind

   prot_ind       = MACROS_indices%prot_ind
   poly_ind       = MACROS_indices%poly_ind
   lip_ind        = MACROS_indices%lip_ind
 
!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   setup_loop: do column = 1, numColumns

   kmax = MACROS_input%number_of_active_levels(column)
   if (kmax < 1) cycle setup_loop

   do k = 1, kmax

    !---------------------------------------------------------------------------
    !   create local copies of requisite ecotracers
    !   treat negative values as zero and apply mask to locals
    !---------------------------------------------------------------------------

   zooC_loc(k,column)      = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,zooC_ind))
   spC_loc(k,column)       = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,spC_ind))
   diatC_loc(k,column)     = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,diatC_ind))
   diazC_loc(k,column)     = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,diazC_ind))
   phaeoC_loc(k,column)    = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,phaeoC_ind))

   prot_loc(k,column)       = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,prot_ind))
   poly_loc(k,column)       = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,poly_ind))
   lip_loc(k,column)        = max(0.0_MACROS_r8, MACROS_input%MACROS_tracers(k,column,lip_ind))

   end do  !  end of setup k loop

   enddo  setup_loop  !  end of setup column loop

!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   column_loop: do column = 1, numColumns

   kmax = MACROS_input%number_of_active_levels(column)
   if (kmax < 1) cycle column_loop

!-----------------------------------------------------------------------
!  loop over levels
!-----------------------------------------------------------------------

   do k = 1, kmax

    !---------------------------------------------------------------------------
    !  Baseline phytoplanktonic carbon release rate constate is here adjusted
    !  by local zooplanktonic densities normalized to an average.
    !  This has the effect of partially decoupling macromolecules from chlorophyll.
    !  Mortality modulates the decoupling but may be zeroed.
    !---------------------------------------------------------------------------

   k_C_p = k_C_p_base * (mort + (zooC_loc(k,column)/zooC_avg))

    !--------------------------------------------------------------------------
    !  Convert to column intermediate
    !--------------------------------------------------------------------------

    spCk_loc      = spC_loc(k,column)
    diatCk_loc    = diatC_loc(k,column)
    phaeoCk_loc   = phaeoC_loc(k,column)
    diazCk_loc    = diazC_loc(k,column)

    !--------------------------------------------------------------------------
    !  Compute total phytoplanktonic carbon.
    !  No need to distinguish ecostructural details at this stage.
    !  Recall that ACME deadlines loom.
    !--------------------------------------------------------------------------

    phytoC_loc  = diatCk_loc + phaeoCk_loc + spCk_loc + diazCk_loc
    
    !-------------------------------------------------------------------------
    !   Construction of kinetic terms for the macromolecules
    !-------------------------------------------------------------------------

    prot_s_disr  = inject_scale*f_prot*k_C_p*phytoC_loc
    poly_s_disr  = inject_scale*f_poly*k_C_p*phytoC_loc
    lip_s_disr   = inject_scale*f_lip *k_C_p*phytoC_loc

    prot_r_bac   = k_prot_bac*prot_loc(k,column)
    poly_r_bac   = k_poly_bac*poly_loc(k,column)
    lip_r_bac    = k_lip_bac *lip_loc(k,column)

    prot_s = prot_s_disr
    poly_s = poly_s_disr
    lip_s  = lip_s_disr

    prot_r = prot_r_bac
    poly_r = poly_r_bac
    lip_r  = lip_r_bac

    MACROS_output%MACROS_tendencies(k,column,prot_ind) = prot_s - prot_r
    MACROS_output%MACROS_tendencies(k,column,poly_ind) = poly_s - poly_r
    MACROS_output%MACROS_tendencies(k,column,lip_ind)  = lip_s  - lip_r
 
! all other tendencies were initialized to 0

! MACROS source terms -keep diagnostics ASAP to start

       MACROS_diagnostic_fields%diag_PROT_S_TOTAL(k,column) = prot_s
       MACROS_diagnostic_fields%diag_POLY_S_TOTAL(k,column) = poly_s
       MACROS_diagnostic_fields%diag_LIP_S_TOTAL(k,column)  = lip_s

       MACROS_diagnostic_fields%diag_PROT_R_TOTAL(k,column) = prot_r
       MACROS_diagnostic_fields%diag_POLY_R_TOTAL(k,column) = poly_r
       MACROS_diagnostic_fields%diag_LIP_R_TOTAL(k,column)  = lip_r

   enddo ! k loop

   enddo column_loop ! i loop

!-----------------------------------------------------------------------
!EOC

 end subroutine MACROS_SourceSink

!***********************************************************************
!  All gas transfer stripped since the macromolecules aren't... gases
!***********************************************************************

 end module MACROS_mod
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
