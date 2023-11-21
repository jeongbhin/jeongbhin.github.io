!****************************************************************************
!****************************************************************************
module constants

  ! Physical or arithmetic constants

  real(kind=8), parameter :: pii    = 3.14159265d0
  real(kind=8), parameter :: twopi  = 2.0d0*pii
  real(kind=8), parameter :: pii43  = 4.0d0*pii/3.0d0
  real(kind=8), parameter :: third = 1.0d0/3.0d0
  
  real(kind=8), parameter :: xmp = 1.672622d-24 ! Proton mass in grams
  real(kind=8), parameter :: amu = 1.660531d-24 ! Atomic mass unit in grams
                                       ! (amu is 1/12 of mass of carbon atom)
  real(kind=8), parameter :: xme = 9.109383d-28 ! Electron mass in grams
  real(kind=8), parameter :: qcgs= 4.803205d-10 ! Proton charge in cgs units:
                                                ! _E_lectro_S_tatic _U_nities
  real(kind=8), parameter :: xkb  = 1.3806488d-16    ! Boltzmann's constant (cgs)
  real(kind=8), parameter :: ccgs = 2.99792458d10   ! speed of light (cm/s)
  real(kind=8), parameter :: xh   = 6.626070d-27     ! Planck constant (erg-sec)
  real(kind=8), parameter :: xh_bar = xh/(2.0d0*pii)
  
  real(kind=8), parameter :: rm_prot = xmp*(ccgs**2) ! proton rest mass en. [erg]
  real(kind=8), parameter :: rm_elec = xme*(ccgs**2) ! electron rest mass en. [erg]
  
  real(kind=8), parameter :: ergtev = 6.242d11  ! conversion ergs to eV
  real(kind=8), parameter :: etkev  = 6.242d08  ! conversion ergs to keV
  real(kind=8), parameter :: etmev  = 6.242d05  ! conversion ergs to MeV
  real(kind=8), parameter :: evterg = 1.602d-12 ! conversion eV to ergs
  real(kind=8), parameter :: xkevte = 1.602d-09 ! conversion keV to ergs
  real(kind=8), parameter :: xmevte = 1.602d-06 ! conversion MeV to ergs
  real(kind=8), parameter :: xjansky  = 1.0d-23 ! Jansky in erg cm^(-2)
  
  real(kind=8), parameter :: degtrd = twopi/360.0  ! Degrees to radians
  real(kind=8), parameter :: radtdg = 1.0d0/degtrd ! Radians to degrees
  real(kind=8), parameter :: sec_p_yr = 3.156d07   ! secs in a year
  real(kind=8), parameter :: pc_to_cm = 3.084d18   ! cm in a parsec

  ! Set the minimum value of theta such that cos(theta) can be
  !  distinguished from 1.0 in whichever precision
  real(kind=8), parameter :: min_tht8 = 1.d1 * sqrt(2.d0*epsilon(1.d0))
  
  real(kind=8), parameter :: B_CMB0 = 3.27d-06 ! Equivalent B field[G] to CMB
                                               !   energy density at a
                                               !   redshift of 0
  real(kind=8), parameter :: T_CMB0 = 2.725d0  ! Temperature[K] of CMB at a
                                               !   redshift of 0
end module constants


!****************************************************************************
!****************************************************************************
module parameters
  
  ! Parameters governing array sizes.  All of these should be interpreted as
  !  maximum values for their respective quantity; the code will quite
  !  happily run if not all the array is used.
  
  integer, parameter :: na_p = 100000  ! Max # of particles at each pcut
  integer, parameter :: na_i = 5       ! Max # of different ion species
  integer, parameter :: na_g = 110     ! Max # of elements in grid arrays
  integer, parameter :: na_c = 100     ! Max # elements in pcut array
  integer, parameter :: na_its = 50    ! Max # of iterations
  
  integer, parameter :: psd_max=200 ! Max # of bins usable for phase space
                                    !  distribution and associated calcs.
                                    ! Applies to both momentm and angular
                                    !  dimensions
     
  integer, parameter :: num_therm_bins=150 ! # of bins in thermal dist
  
  integer, parameter :: na_cr = 10 * na_p  ! Max # of thermal particle 
                                           !   crossings to hold before
                                           !   writing to scratch file
  
  real(kind=8), parameter :: beta_rel_fl = 0.02d0   ! Cutoffs for non-rel vs
  real(kind=8), parameter :: en_rel_pt   = 0.005d0  !   rel eqs, for fluid &
                                                    !   for particles
  
end module parameters


!****************************************************************************
!****************************************************************************
module controls
  use parameters, only: na_g, na_i, na_c
  
  ! Variables here are runtime constants, not compile-time constants
  ! Most of these are set during the call to subroutine data_input
  
  integer :: n_ions   ! Number of ion species
  integer :: n_itrs   ! Number of iterations to perform
  integer :: n_pcuts  ! Number of momentum levels to split at
  integer :: n_tcuts  ! Number of times to use in particle tracking
  integer :: n_xspec  ! Number of locations to take spectra
  integer :: i_in_distr  ! Form of input distribution. 1 = thermal,
                         !   2 = delta function, 3 = other (not implemented)
  integer :: n_pts_inj, n_pts_pcut, n_pts_pcut_hi  ! Number of CRs to use
                                                   !   during pcuts
  integer :: n_old_skip, n_old_profs, n_old_per_prof  ! Used when reading in
                                                      !   old profiles
  integer :: psd_bins_per_dec_mom,  &!& Number of bins per decade to use in
             psd_bins_per_dec_tht    !    PSD in mom and ang dimensions
  integer :: psd_lin_cos_bins,      &!& Number of linear-spaced (cos) bins,
             psd_log_tht_decs        !    # of log-spaced (theta) decades
  
  
  logical :: l_sc_elec  ! Are electrons included as a separate species?
  logical :: l_inj_wt   ! How to assign initial weights to injected particle
                        !   distribution.  "T" = all particles have equal
                        !   weights.  "F" = all bins have equal weights.
  logical :: l_oblique  ! Whether shock is parallel to UpS B field
  logical :: l_use_prp  ! Use a probability of return plane DwS rather than
                        !   a free escape boundary
  logical :: dont_shock   ! Force compression ratio of shock to 1.0, i.e.
                          !   eliminates shock entirely.  For testing
                          !   purposes only.
  logical :: dont_scatter   ! Eliminate calls to scattering subroutine.
                            !   For testing purposes only.
  logical :: dont_DSA   ! Prevent particles from moving DwS -> UpS, i.e.
                        !   entering DSA process.  For testing purposes only.
  logical :: do_smoothing   ! Allow smoothing of shock profile between
                            !   iterations.
  logical :: do_prof_fac_damp   ! Increase old profile weight at later
                                !   iterations, damping effect of smoothing
  logical :: do_fast_push   ! Allow fast transport to specified UpS location
                            !   instead of propagating from x_start_grid
  logical :: do_rad_losses   ! Calculate radiative losses for electrons
  logical :: do_old_prof   ! Read in an old profile rather than starting
                           !   with fresh UM shock
  logical :: do_retro   ! Use retro time calculation for DwS CRs beyond PRP;
                        !   ignored if not tracking CR age
  logical :: use_custom_frg   ! Set dependence of del_tht_max on gyroradius
                              !   in subr "scattering"
  logical :: use_custom_epsB  ! Set a custom function for btot_grid in subr
                              !   setup_profile
  logical :: do_multi_dNdps   ! Write dN/dp for each iteration to a separate
                              !   file?  If "false", only last one is kept
  logical :: do_tcuts   ! Track particle ages during propagation?
  

!--! Floats related to profile in/out speeds
  real(kind=8) :: u_Z, beta_Z, gam_Z  ! Shock speed (cm/s, u/c, Lorentz
                                      !   factor) in shock frame
  real(kind=8) :: r_comp, rRH         ! Compression ratio (u_Z/u_2 in shock
                                      !   frame), and Rankine-Hugoniot value
  real(kind=8) :: u_2, beta_2, gam_2  ! DwS speeds in shock frame
  real(kind=8) :: theta_u2            ! If in oblique shock, angle[deg] btwn
                                      !   u_2 and shock normal
  real(kind=8) :: mach_sonic,     &!& ! Sonic and Alfven Mach numbers for
                  mach_alfven         !   shock
  
!--! Floats related to other UpS/DwS conditions
  real(kind=8) :: bmag_Z        ! Far UpS magnetic field strength[Gauss]
  real(kind=8) :: bmag_2        ! Far DwS magnetic field strength[Gauss]
  real(kind=8) :: theta_BZ      ! Angle[deg] btwn UpS B and shock normal
  real(kind=8) :: theta_B2      ! Angle[deg] btwn DwS B and shock normal
  real(kind=8) :: tZ_elec  ! Far UpS electron temperature[K]; ignored if
                           !   electrons are separate particle species
  real(kind=8) :: en_inj   ! Injection energy[keV] if UpS dist is delta
                           !   function.  Ignored otherwise.
  real(kind=8) :: emin_therm_fac  ! Factor setting minimum PSD energy if
                                  !   UpS dist is thermal
  
!--! Far UpS fluxes
  real(kind=8) :: flux_px_UpS, &!& ! Far UpS momentum flux, x-component
                  flux_pz_UpS, &!& ! Far UpS momentum flux, z-component
                  flux_en_UpS      ! Far UpS energy flux
  
!--! Floats related to the grid
  real(kind=8) :: rg0  ! Conversion factor between "code" and cgs distance.
                       !   Set in "BMAGZ".  Equal to gyroradius of proton
                       !   with v = u_Z
  real(kind=8) :: x_grid_start_rg,   &!& Start and stop positions[rg0] for
                  x_grid_stop_rg      !    the grid
  real(kind=8) :: x_grid_start,      &!& Start and stop positions[cm] for
                  x_grid_stop         !    the grid
  real(kind=8) :: x_fast_stop_rg   ! Location[rg0] at which fast push ends
  
  real(kind=8) :: en_pcut_hi  ! Energy separating n_pts_pcut & n_pts_pcut_hi;
                              !   actual pcut depends on CR species
!--! Floats related to scattering
  real(kind=8) :: eta_mfp   ! Ratio of MFP to gyroradius. = 1 is Bohm
                            !   diffusion
  real(kind=8) :: xn_per_coarse,   &!& Number of divisions of gyro period
                  xn_per_fine       !    when using coarse/fine scattering
  
!--! Floats related to dN/dp normalization and Cosmic Microwave Background
  real(kind=8) :: jet_dist_kpc, redshift  ! Distance between observer and
                                          !   source, & associated redshift;
                                          !   needed to handle CMB
  real(kind=8) :: jet_rad_pc       ! Jet radius and two measures for surface
  real(kind=8) :: jet_sph_frac,   &!&  area
                  jet_open_ang_deg
  
!--! Floats related to limits on CR energy
  real(kind=8) :: age_max  ! Maximum CR age[sec] in explosion frame
  real(kind=8) :: Emax_keV, Emax_keV_per_aa, pmax_cgs   ! Limits on CR energy
  real(kind=8) :: feb_UpS, feb_DwS   ! UpS and DwS locations[cm] of FEBs
  
!--! Floats related to electrons
  real(kind=8) :: p_elec_crit,   &!& Momentum[cgs] and Lorentz factor below
                  gam_elec_crit   !    which electrons have constant MFP
  real(kind=8) :: en_xfer_frac   ! Fraction of UpS (shock frame) ion energy
                                 !   transferred to electrons at 1st shock
                                 !   crossing
  
!--! Floats related to increasing B field
  real(kind=8) :: bturb_comp_frac  ! Degree to which Bturb compression is
                                   !   used in propagation & rad losses
  real(kind=8) :: bfield_amp       ! Additional amplification of B field
                                   !   beyond Bturb compression; has no
                                   !   effect if bturb_comp_frac = 0
  
!--! Floats related to smoothing
  real(kind=8) :: prof_wt_fac   ! Weighting factor for old profile when
                                !   smoothing. > 1 favors old prof
  real(kind=8) :: smooth_mom_en_fac   ! Degree to use momentum/energy eqs
                                      !   when smoothing.  0 for all mom.,
                                      !   1 for all energy
  real(kind=8) :: smooth_press_flux_psd_fac   ! Degree to use flux eqs (=0)
                                              !   PSD (=1) to find pressure
                                              !   when smoothing
  real(kind=8) :: x_art_start_rg, x_art_scale   ! If using artif. smoothing,
                                                !   start and scale factors
  
!--! Mass, charge, UpS temp and UpS density for particle species
  real(kind=8), dimension(na_i) :: aa_ion = 0.d0,  zz_ion   = 0.d0,       &!&
                                   tZ_ion = 0.d0,  denZ_ion = 0.d0
!--! Injection rate for each particle species (0=no DSA, 1=thermal leakage)
  real(kind=8), dimension(na_i) :: inj_fracs = 1.d0

!--! Cutoff momenta[aa*m_pc] btwn pcuts.  Converted to cgs units in-run
  real(kind=8), dimension(na_c) :: pcuts_in
  

  real(kind=8), dimension(na_g) :: x_spec ! Locations where particle spectra
                                          !   should be calculated.  Used
                                          !   independently of PSD
  
!--! Array to hold cutoff times for particle tracking
  real(kind=8), dimension(na_c) :: tcuts
  
end module controls


!****************************************************************************
!****************************************************************************
module psd_vars
  use parameters, only: na_g, psd_max
  
    ! Related to total number of bins along each dimension; bins_per_dec
    !   held in module controls as they're read in from the input file
  integer :: num_psd_mom_bins, num_psd_tht_bins
  
  real(kind=8) :: psd_mom_min, psd_cos_fine, del_cos, psd_tht_min
  
  real(kind=8), dimension(0:psd_max) :: psd_mom_bounds, psd_tht_bounds
  
  real(kind=8), dimension(0:psd_max, 0:psd_max, na_g) :: psd

end module psd_vars


!****************************************************************************
!****************************************************************************
module grid_vars
  use parameters, only: na_g
  
  ! The grid *boundaries* below refer to the left-hand side of the particular
  !  zone.  Other variables (uxsk, etc.) are the value within that zone.

  integer :: n_grid   ! Number of grid zones
  
  real(kind=8), dimension(0:na_g) :: x_grid_cm, &!& Grid boundaries, in cgs
                                     x_grid_rg   !    and rg0 units
  
  real(kind=8), dimension(0:na_g) :: uxsk_grid,   &!& x and z components of
                                     uzsk_grid,   &!&   bulk flow velocity in
                                     utot_grid     !    shock frame, and tot
  real(kind=8), dimension(0:na_g) :: gam_sf_grid, &!& Bulk flow Lorentz
                                     gam_ef_grid   !    factor, shock and 
                                                   !    explosion frames
  real(kind=8), dimension(0:na_g) :: beta_ef_grid  !  Bulk flow speed/c,
                                                   !    explosion frame
  
  real(kind=8), dimension(0:na_g) :: btot_grid   ! Total magnetic field
  real(kind=8), dimension(0:na_g) :: theta_grid  ! Angle[rad] of B to shock
                                                 !   normal (i.e. x axis)
  real(kind=8), dimension(0:na_g) :: epsB_grid   ! Energy density fraction in
                                                 !   magnetic field
  
end module grid_vars


!****************************************************************************
!****************************************************************************
module iteration_vars

  use parameters, only: na_g, na_i, na_its, psd_max, na_p, na_c
  
  integer :: i_itr
  
!--! Arrays used to track the flux of particles on (or off) the grid
  real(kind=8), dimension(na_g) :: pxx_flux, pxz_flux, en_flux
  real(kind=8), dimension(na_i) :: esc_flux
  real(kind=8), dimension(na_i, na_its) :: px_esc_feb, en_esc_feb
  real(kind=8), dimension(0:psd_max, na_i) :: esc_en_eff, esc_num_eff

!--! Arrays for tracking pressure & other thermo quantities
  real(kind=8), dimension(na_g) :: press_psd_par, press_psd_perp, en_den_psd
  real(kind=8), dimension(na_g, 2) :: gam_adiab_grid

!--! Arrays for holding thermal distribution information; they're set at
   !   the start of the run, but included here because of chance they could
   !   change due to fast push
  integer, dimension(na_i) :: n_pts_MB
  real(kind=8), dimension(na_p, na_i) :: ptot_inj, wt_inj
  
!--! Arrays for holding information about particle counts and spectra at
   !   various tcuts
  real(kind=8), dimension(na_c, na_i) :: wt_coupled
  real(kind=8), dimension(0:psd_max, na_c, na_i) :: spec_coupled
  
end module iteration_vars


!****************************************************************************
!****************************************************************************
module species_vars
  use parameters, only: na_p, na_cr, psd_max, na_g
  
  ! Variables here refer to the species currently being propagated within
  !  the shock structure
  !
  ! See top of loop_pcut for explanation of particle arrays
  !
  ! zz_in is particle charge, in units of proton charge; not currently used,
  !   but declared in case ionization gains or losses are ever implemented
  
  integer :: i_ion
  
  real(kind=8) :: o_o_mc
  
  integer, dimension(na_p) :: i_grid_in
  real(kind=8), dimension(na_p) :: xwt_in, ptpf_in, pbpf_in, x_PT_cm_in,  &!&
       zz_in
  
!--! Arrays will hold crossing data for thermal particles; px and pt are
   !   shock frame values
  integer :: n_cr_count
  integer, dimension(na_g) :: num_crossings
  integer, dimension(na_cr) :: therm_grid
  real(kind=8), dimension(na_cr) :: therm_pxsk, therm_ptsk, therm_xwt
  
  
!--! Spectra at x_spec locations
  real(kind=8), dimension(0:psd_max, na_g) :: spec_sf, spec_pf
  
  
!--! Escaping spectra UpS and DwS from shock; 2-D arrays store angular
   !   information
  real(kind=8), dimension(0:psd_max) :: esc_spec_feb_UpS, esc_spec_feb_DwS
  real(kind=8), dimension(0:psd_max, 0:psd_max) :: esc_psd_feb_UpS,       &!&
       esc_psd_feb_DwS

end module species_vars


!****************************************************************************
!****************************************************************************
module pcut_vars
  use parameters, only: na_p
  
  ! Particle arrays with the suffix "_sav" are filled with particles to be
  !  carried over to the next pcut.  During subroutine next_pcut the values
  !  are transferred with any necessary replications to arrays with the
  !  suffix "_new"
  
  ! i_grid:     Current grid zone number of particle
  ! tcut:       Next time at which tracking should occur
  ! l_DwS:      Whether the particle has been downstream
  ! l_inj:      Whether the particle has been injected into DSA, i.e. has
  !               made at least one UpS -> DwS -> UpS cycle
  ! xwt:        Particle weight
  ! ptpf, pbpf: Total momentum (plasma frame) of particle, and the component
  !               parallel to the magnetic field
  ! x_PT_cm:    Current position of particle
  ! xn_per_period: Number of steps to divide gyroperiod into; more means
  !               finer scattering
  ! zz:         Particle charge, in units of proton charge
  ! prp_x_cm:   Location of DwS probability of return plane
  ! acctime_sec: Total elapsed time (in explosion frame, not shock or plasma)
  !               since crossing shock for first time
  ! phi_rad:    Phase angle of gyro-period
  !DOLATER: the variable l_DwS is superfluous.  After pcut #1 all particles
  !  are guaranteed to have been downstream, and before pcut #1 all
  !  particles are guaranteed NOT to have been downstream.  Set at the top
  !  of loop_pcut instead of carrying around the extra memory.
  
  integer :: i_cut, n_pts_use
  
  logical, dimension(na_p) :: l_save     ! Whether or not to save particle
                                         !  for next pcut
  
  integer, dimension(na_p) :: i_grid_sav, tcut_sav
  logical, dimension(na_p) :: l_DwS_sav, l_inj_sav
  real(kind=8), dimension(na_p) :: xwt_sav, ptpf_sav, pbpf_sav,           &!&
       x_PT_cm_sav, xn_per_sav, zz_sav, prp_x_cm_sav, acctime_sec_sav,    &!&
       phi_rad_sav
  
  integer, dimension(na_p) :: i_grid_new, tcut_new
  logical, dimension(na_p) :: l_DwS_new, l_inj_new
  real(kind=8), dimension(na_p) :: xwt_new, ptpf_new, pbpf_new,           &!&
       x_PT_cm_new, xn_per_new, zz_new, prp_x_cm_new, acctime_sec_new,    &!&
       phi_rad_new
  
end module pcut_vars


!****************************************************************************
!****************************************************************************
module randnum
  implicit none
  
  integer :: iseed_in
  integer :: x=123456789, y=362436069, z=521288629, w=916191069
!$omp threadprivate(x, y, z, w)

contains
   
   subroutine kiss(rand_out)
      integer :: i
      real(kind=8) :: rand_out
! The  KISS (Keep It Simple Stupid) random number generator.
! http://www.fortran.com/kiss.f90 . Rather modified from original.

      do while(.true.)
         x = 69069 * x + 1327217885
         y = m (m (m (y, 13), - 17), 5)
         z = 18000 * iand (z, 65535) + ishft (z, - 16)
         w = 30903 * iand (w, 65535) + ishft (w, - 16)
         i = x + y + ishft (z, 16) + w
         
         rand_out = i*2.33d-10 + 0.5d0
         if((rand_out .gt. 0.d0) .and. (rand_out .lt. 1.d0)) return
      enddo
   
   end subroutine
   
   function m(k, n)
         integer :: m, k, n
         m = ieor (k, ishft (k, n) )
   end function m
   
   subroutine kisset (ix, iy, iz, iw)
      integer :: ix, iy, iz, iw
      x = ix
      y = iy
      z = iz
      w = iw
   end subroutine kisset

end module randnum


!****************************************************************************
!****************************************************************************
module debug
  use parameters
  
  real(kind=8) :: running_taken, running_given, running_given_2, p_avg, en_avg
  real(kind=8), dimension(na_g, na_i) :: en_den = 0.d0, therm_en_den = 0.d0
  
  real(kind=8) :: den_pf
  real(kind=8), dimension(na_g) :: zone_vol, en_pool
end module


!****************************************************************************
!****************************************************************************
! Monte Carlo program begins here
!****************************************************************************
!****************************************************************************

program MC

  use constants
  use parameters
  use controls
  use psd_vars
  use grid_vars
  use iteration_vars
  use species_vars
  use pcut_vars
  use randnum
  !DEBUGLINE
  use debug
  
  implicit none

  ! Loop indices and miscellaneous quantities
integer :: i, i_tmp, i_grid_feb, n_print_pt, i_prt, n_pts_max, i_shock,   &!&
     n_saved, n_pts_target, iseed_mod, i_fin, helix_count
real(kind=8) :: gam_adiab_2_RH, rand, wt_running, pcut_prev, jet_dist_Mpc
  ! For tracking time
real(kind=8) :: t_start, t_end, run_time
  ! Related to the phase space distribution
real(kind=8) :: psd_tht_fine, Emin_keV, aa_min, gam, aa_max, psd_mom_max
  ! Related to radiative losses
real(kind=8) :: elec_mass_ratio, elec_rm, sigma_T, rad_loss_fac, B_CMBz
  ! For holding escaping fluxes
real(kind=8) :: px_esc_flux_UpS_tot, en_esc_flux_UpS_tot, sum_press_DwS,  &!&
     sum_KEden_DwS, px_esc_UpS, en_esc_UpS
real(kind=8), dimension(na_its) :: px_esc_flux_UpS, en_esc_flux_UpS
  ! EXPECTED escaping fluxes
real(kind=8) :: q_esc_cal_en_avg, q_esc_cal_px_avg
real(kind=8), dimension(na_its) :: q_esc_cal_en, q_esc_cal_px
  ! For energy transfer from ions to electrons
integer :: i_start, i_stop, n_split
real(kind=8) :: z_max, z_curr, elec_wt_fac, frac_to_keep, gampf_i,        &!&
     gampf_f, ptpf_f, scale_fac, en_to_Xfer
real(kind=8), dimension(na_g) :: eps_target, en_Xfer_pool, en_recv_pool
  ! Related to the current particle species
real(kind=8) :: aa, zz, p_pcut_hi, pmax_cutoff
real(kind=8), dimension(na_c) :: pcuts_use
  ! Properties of the individual particle
integer :: i_grid, i_return, i_grid_old, tcut_curr
logical :: l_DwS, l_inj
real(kind=8) :: xwt, ptpf, pbpf, x_PT_cm, xn_per, prp_x_cm, acctime_sec,  &!&
     phi_rad, y_PT_cm, z_PT_cm, gam_ptpf, gyro_denom, gyro_rad_tot_cm,    &!&
     p_perp_b_pf, gyro_rad_cm, ptpf_old, gyro_period_sec, x_PT_old,       &!&
     y_PT_old, z_PT_old, phi_rad_old, gyro_rad_tmp, vel
  ! Properties of the bulk flow
real(kind=8) :: uxsk, uzsk, utot, gam_usf, gam_uef, beta_uef, bmag, btht, &!&
     b_sin_th, b_cos_th, uxsk_old, uzsk_old, utot_old, gam_usf_old,       &!&
     b_sin_old, b_cos_old, ptsk, pxsk, pysk, pzsk, pbsk, p_perp_b_sk,     &!&
     gam_ptsk
  ! Other variables used in loop_pt
integer :: i_reason, nc_unit
logical :: first_iter, keep_looping, do_therm_pt_formatted, do_prob_ret,  &!&
     lose_pt, loop_again
real(kind=8) :: B_CMB_loc, B_tot_sq, dp_synch, t_step, x_move_bpar,       &!&
     vpt_o_u2, L_diff
  ! Arrays that hold dN/dp's & associated quantities
real(kind=8), dimension(na_g) :: zone_pop
real(kind=8), dimension(0:psd_max,0:psd_max,na_g) :: d2N_dpdcos_ef
real(kind=8), dimension(0:psd_max,na_g,3) :: dNdp_therm,                  &!&
     dNdp_therm_pvals, dNdp_cr
real(kind=8), dimension(na_its) :: gam_adiab_DwS
  ! End of iteration writeout
integer :: n_avg
real(kind=8) :: px_esc_avg, en_esc_avg, gamma_DwS_esc, gamma_UpS
  
  
!--! Start the wall clock for this run
  call get_time(t_start)
  
  
!--! Get input, control variables, etc.
  call data_input(u_Z, gam_Z, beta_Z, n_ions, aa_ion, zz_ion, tZ_ion,     &!&
     denZ_ion, l_sc_elec, tZ_elec, i_in_distr, en_inj, l_inj_wt, Emax_keV,&!&
     Emax_keV_per_aa, pmax_cgs, eta_mfp, bmag_Z, rg0, theta_BZ, l_oblique,&!&
     x_grid_start_rg, x_grid_stop_rg, feb_UpS, feb_DwS, l_use_prp,        &!&
     n_xspec, x_spec, n_itrs, xn_per_coarse, xn_per_fine, n_pts_inj,      &!&
     n_pts_pcut, n_pts_pcut_hi, en_pcut_hi, n_pcuts, pcuts_in, dont_shock,&!&
     dont_scatter, dont_DSA, do_smoothing, prof_wt_fac, do_prof_fac_damp, &!&
     smooth_mom_en_fac, smooth_press_flux_psd_fac, r_comp, rRH,           &!&
     do_old_prof, n_old_skip, n_old_profs, n_old_per_prof, age_max,       &!&
     do_retro, do_fast_push, x_fast_stop_rg, x_art_start_rg, x_art_scale, &!&
     p_elec_crit, gam_elec_crit, do_rad_losses, jet_rad_pc, jet_sph_frac, &!&
     jet_open_ang_deg, jet_dist_kpc, redshift, en_xfer_frac,              &!&
     bturb_comp_frac, bfield_amp, gam_adiab_2_RH, psd_bins_per_dec_mom,   &!&
     psd_bins_per_dec_tht, psd_lin_cos_bins, psd_log_tht_decs, u_2,       &!&
     beta_2, gam_2, bmag_2, theta_B2, theta_u2, use_custom_frg,           &!&
     emin_therm_fac, do_multi_dNdps, do_tcuts, n_tcuts, tcuts, inj_fracs, &!&
     use_custom_epsB)
     
     
!--! Set quantities related to the phase space distribution, including the
   !   bins
  psd_cos_fine = 1.d0  -  2.d0 / real(psd_lin_cos_bins+1)
  psd_tht_fine = acos(psd_cos_fine)
  psd_tht_min  = psd_tht_fine / 10.d0**psd_log_tht_decs
  
  if( i_in_distr .eq. 1 ) then
    ! Set minimum PSD energy using thermal distribution for UpS plasma
    Emin_keV = huge(1.d0)
    do i = 1, n_ions
      if( (xkb * tZ_ion(i) * etkev) .lt. Emin_keV ) then
        Emin_keV = xkb * tZ_ion(i) * etkev
      endif
    enddo
    
    ! Allow for a few extra zones below the thermal peak
    Emin_keV = emin_therm_fac * Emin_keV
  
  else if( i_in_distr .eq. 2 ) then
    ! Set minimum PSD energy using delta-function dist for UpS plasma;
    !   allow for a few extra zones below the location of the distribution
    Emin_keV = 0.2 * en_inj
    
  endif
  
  ! Determine minimum momentum associated with the given energy, which will
  !   occur for the lightest particle species.  Use a cutoff of 0.1% of the
  !   rest-mass energy for rel/non-rel calculation.
  aa_min = minval( aa_ion(1:n_ions) )
  if( (Emin_keV*xkevte) .lt. (1.d-3*aa_min*rm_prot) ) then
    psd_mom_min = sqrt( 2.d0 * aa_min*xmp * Emin_keV * xkevte )
  else
    gam         = 1.d0  +  (Emin_keV*xkevte)/(aa_min*rm_prot)
    psd_mom_min = aa_min * xmp * ccgs * sqrt( gam**2 - 1.d0 )
  endif
  
  ! Now find the maximum momentum for the PSD (this will be adjusted due to
  !   SF->PF Lorentz transformation).  How to actually calculate it depends
  !   on the user-specified maximum energy "ENMAX"
  aa_max = maxval( aa_ion(1:n_ions) )
  if( Emax_keV .gt. 0.d0 ) then
    gam         = 1.d0  +  (Emax_keV*xkevte)/(aa_max*rm_prot)
    psd_mom_max = aa_max * xmp * ccgs * sqrt( gam**2 - 1.d0 )
  
  else if( Emax_keV_per_aa .gt. 0.d0 ) then
    gam         = 1.d0  +  (Emax_keV_per_aa*xkevte/rm_prot)
    psd_mom_max = aa_max * xmp * ccgs * sqrt( gam**2 - 1.d0 )
  
  else if( pmax_cgs .gt. 0.d0 ) then
    psd_mom_max = pmax_cgs
  
  else
    ! Something has gone very wrong.
    write(*,"(2A)") 'Max CR energy not set in data_input, so can not ',   &!&
        'set PSD bins.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  ! Adjust max momentum based on a SF->PF Lorentz transform
  psd_mom_max = 2.d0 * gam_Z * psd_mom_max
  
  call set_psd_bins(psd_mom_min, psd_mom_max, psd_bins_per_dec_mom,       &!&
     psd_bins_per_dec_tht, psd_lin_cos_bins, psd_cos_fine, psd_tht_min,   &!&
     num_psd_mom_bins, num_psd_tht_bins, del_cos, psd_mom_bounds,         &!&
     psd_tht_bounds)
  
  
!--! Set up the computational grid
  call setup_grid(x_grid_start, x_grid_stop, n_grid, x_grid_rg, x_grid_cm)
  
  
!--! Check x_spec data to make sure it falls within the grid, and add it to
   !   the output file
  if( n_xspec .gt. 0 ) then
    do i = 1, n_xspec
      if( (x_spec(i) .lt. x_grid_start) .or.                              &!&
          (x_spec(i) .gt. x_grid_stop) ) then
        write(*,"(A,I0,2A)") 'ERROR: x_spec position ',i,' falls ',       &!&
            'outside grid start/stop bounds.'
        write(*,"(A)") 'Stopping program now.'
        stop
      endif
      
      write(9,"(3X,2A,I2,2(3X,ES10.3E2))") 'x position[rg0,pc] for ',     &!&
          'spectrum calculation:  ',i,x_spec(i)/rg0, x_spec(i)/pc_to_cm
    enddo
    
    write(9,*)
  endif
  
  
!--! Get grid zone numbers for the location of the UpS FEB
   do i = 1, n_grid
     if( x_grid_cm(i) .gt. feb_UpS ) then
       i_grid_feb = i - 1
       exit
     endif
   enddo
  
  
!--! Because redshift will be needed to compute radiative losses (it affects
   !   both the energy and density of CMB photons), calculate it here.  If
   !   redshift was provided during data_input, cosmo_calc will return
   !   distance instead
   ! Cosmo_calc expects distance in megaparsecs, so convert from value read
   !   in during data_input
  jet_dist_Mpc = jet_dist_kpc * 1.d-3
  call cosmo_calc(jet_dist_Mpc, redshift) 
  
  
!--! Set a handful of constants related to radiative losses.  elec_rm will
   !   not be the same as rm_elec in module "constants" unless (a) there are
   !   electrons in the run, and (b) they are true electrons, with aa = 
   !   xme/xmp.  (elec_rm may in fact = rm_prot, but in that case it won't
   !   be used because radiative losses will never be calculated.)
   ! To convert rad_loss_fac to dp/dt, multiply by p^2*B^2, both in cgs.
   ! Prefactor of rad_loss_fac comes from average over pitch and is Eq (16)
   !   of Sturner+ (1997) [1997ApJ...490..619S].  Note the extra factor of
   !   c in the denominator, because code tracks dp/dt, not dE/dt as given
   !   in Sturner+ (1997).
  elec_mass_ratio = minval( aa_ion(1:n_ions) )
  elec_rm = elec_mass_ratio * rm_prot
  sigma_T = 8.d0*third*pii * (qcgs**2 / elec_rm)**2
  rad_loss_fac = 4.d0*third * ccgs * sigma_T                              &!&
                /  (ccgs**3 * (xmp*elec_mass_ratio)**2 * 8.d0*pii)
  B_CMBz  = B_CMB0 * (1.d0 + redshift)**2
  
  
!--! Zero out total escaping fluxes and calculate the far UpS fluxes
  px_esc_flux_UpS_tot = 0.d0
  en_esc_flux_UpS_tot = 0.d0
  px_esc_flux_UpS(:)  = 0.d0
  en_esc_flux_UpS(:)  = 0.d0
  call upstream_fluxes(flux_px_UpS, flux_pz_UpS, flux_en_UpS)
  
  
!--! Determine upstream Mach numbers (sonic & Alfven)
  call upstream_machs(mach_sonic, mach_alfven)
  
  
!--! Set up the initial shock profile, or read it in from a file
  if( .not. do_old_prof ) then
    call setup_profile(uxsk_grid, uzsk_grid, utot_grid, gam_sf_grid,      &!&
       beta_ef_grid, gam_ef_grid, btot_grid, theta_grid, epsB_grid, bmag_2)
  else
    call read_old_prof(n_old_skip, n_old_profs, n_old_per_prof, x_grid_rg,&!&
       x_grid_cm, uxsk_grid, uzsk_grid, utot_grid, gam_sf_grid,           &!&
       beta_ef_grid, gam_ef_grid, btot_grid, epsB_grid, theta_grid,       &!&
       n_grid, u_Z, gam_Z, rg0, r_comp, rRH, beta_Z, bmag_Z, u_2, beta_2, &!&
       gam_2, theta_u2, bmag_2, theta_BZ, theta_B2, flux_px_UpS,          &!&
       flux_pz_UpS, flux_en_UpS)
    ! Must set far UpS and DwS limits manually, since they won't be read in
    !   from the file
    x_grid_rg(0)        = -1.d30
    x_grid_rg(n_grid+1) =  1.d30
    x_grid_cm(0)        = -1.d30 * rg0
    x_grid_cm(n_grid+1) =  1.d30 * rg0
  endif
  
  
!--! Find the location of the shock
  do i = 1, n_grid
    if( (x_grid_rg(i) .eq. 0.d0) .or.                                     &!&
       ((x_grid_rg(i) .lt. 0.d0) .and. (x_grid_rg(i+1) .gt. 0.d0)) ) then
      i_shock = i
      exit
    endif
  enddo
  
  
!--! How frequently will the code print during initial particle propagation?
  if( n_pts_inj .lt. 500 ) then
    n_print_pt = 25
  else
    n_print_pt = 500
  endif
  
  
!--! The random number seed depends on the specific particle.  Compute a
   !   necessary quantity to determine that seed
  n_pts_max = max(n_pts_pcut, n_pts_pcut_hi)
  
  
!--! Because electrons might have different Monte Carlo weights than protons,
   !   protons, set that ratio here
   ! WARNING: assumes electrons are last ion species of input file
  if( l_sc_elec ) then
    elec_wt_fac = 1.d0 / denZ_ion(n_ions)
  else
    elec_wt_fac = 0.d0
  endif
  
  
!--! Print a bunch of data about the run to screen/file
  call print_input(n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,          &!&
     num_psd_mom_bins, num_psd_tht_bins, n_xspec, n_pcuts, n_grid, rRH,   &!&
     r_comp, u_Z, beta_Z, gam_Z, u_2, beta_2, gam_2, denZ_ion, bmag_Z,    &!&
     bmag_2, theta_BZ, theta_B2, theta_u2, aa_ion, tZ_ion, tZ_elec,       &!&
     mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine, feb_UpS,        &!&
     feb_DwS, rg0, age_max, en_pcut_hi, do_fast_push, bturb_comp_frac)


  !**************************************************************************
  !**************************************************************************
  !  Main computational loops
  !
  !  Loops are set in the following sequence of nesting
  !   1) loop_itr:   i_itr    Iteration number
  !   2) loop_ion:   i_ion    Particle species number
  !   3) loop_pcut:  i_cut    Particle splitting (pcut) number
  !   4) loop_pt:    i_prt    Individual particle number
  !
  !
  ! Start of loop over iterations
  !**************************************************************************
  loop_itr: do i_itr = 1, n_itrs
  
  !--! Zero out numerous quantities that will be modified over the course of
     !   this iteration.  Minimally positive number is used to prevent errors
     !   when taking logarithms later
    pxx_flux(:) = 1.d-99
    pxz_flux(:) = 1.d-99
    en_flux(:)  = 1.d-99
    
    esc_spec_feb_UpS(:) = 1.d-99
    esc_spec_feb_DwS(:) = 1.d-99
    
    press_psd_par(:)  = 1.d-99
    press_psd_perp(:) = 1.d-99
    en_den_psd(:)     = 1.d-99
    
    wt_coupled(:,:)   = 1.d-99
    
    
  !--! Additionally, set/reset scalar quantities that will change
    sum_press_DwS = 1.d-99
    sum_KEden_DwS = 1.d-99
    
    en_esc_UpS    = 1.d-99
    px_esc_UpS    = 1.d-99
    
    
  !--! To facilitate energy transfer from ions to electrons, calculate here
     !   the target energy density fraction for electrons at each grid zone,
     !   and zero out the pool of plasma-frame energy that will be taken from
     !   ions and donated to electrons
     ! Per Ardaneh+ (2015) [2015arXiv150705374A], eps_elec is proportional to
     !   sqrt(eps_b).  eps_b is itself roughly proportional to density**2 --
     !   B**2 is proportional to z**2 (z being the density compression
     !   factor) -- so eps_elec should vary roughly linearly with density.
    z_max = gam_Z * beta_Z / (gam_2 * beta_2)
    do i = 1, n_grid
      z_curr = gam_Z * u_Z / (gam_sf_grid(i) * uxsk_grid(i))
      
      if( uxsk_grid(i) .eq. u_Z ) then
        eps_target(i) = 0.d0
      else
        eps_target(i) = en_xfer_frac * (z_curr - 1.d0) / (z_max - 1.d0)
      endif
    enddo
    
    en_Xfer_pool(:) = 0.d0
    en_recv_pool(:) = 0.d0
    eps_target(n_grid+1:) = 0.d0
    !DEBUGLINE
    en_den(:,:) = 0.d0
    therm_en_den(:,:) = 0.d0
    
    
    !************************************************************************
    !  Start of loop over particle species
    !
    !  First species always protons.  If electrons present they MUST be last
    !    species.
    !  Each species has mass number "aa" and charge number "zz" in units of
    !    proton mass and charge, respectively.
    !************************************************************************
    loop_ion: do i_ion = 1, n_ions
      
      aa = aa_ion(i_ion)
      zz = zz_ion(i_ion)
      
      o_o_mc = 1.d0 / (aa*xmp * ccgs)
      
      
    !--! At the start of each ion, print a glyph to the screen
      write(*,"(2(A,I0))") '*** Iteration # ', i_itr, '   species # ', i_ion
      
      
    !--! Determine the pcut at which to switch from low-E particle counts to
       !   high-E particle counts.  Recall that en_pcut_hi has units of keV
       !   per aa, so when dividing by particle mass the factor of aa is
       !   already present in the denominator.
       ! Also set the maximum momentum cutoff based on the values given in
       !   keyword "ENMAX"
      if( (en_pcut_hi * xkevte / rm_prot) .lt. en_rel_pt ) then
        p_pcut_hi = sqrt( 2.d0 * en_pcut_hi * xkevte / rm_prot )
      else
        p_pcut_hi = aa * xmp * ccgs                                       &!&
                   * sqrt( (en_pcut_hi*xkevte/rm_prot + 1.d0)**2  -  1.d0 )
      endif
      
      if( Emax_keV .gt. 0.d0 ) then
        gam         = 1.d0  +  (Emax_keV*xkevte)/(aa*rm_prot)
        pmax_cutoff = aa*xmp * ccgs * sqrt( gam**2 - 1.d0 )
      else if( Emax_keV_per_aa .gt. 0.d0 ) then
        gam         = 1.d0  +  (Emax_keV_per_aa*xkevte/rm_prot)
        pmax_cutoff = aa*xmp * ccgs * sqrt( gam**2 - 1.d0 )
      else if( pmax_cgs .gt. 0.d0 ) then
        pmax_cutoff = pmax_cgs
      else
        ! Something has gone very wrong.
        write(*,"(2A)") "Max CR energy not set in data_input, so can't ", &!&
            'set pmax_cutoff.'
        write(*,"(A)") 'Stopping program now.'
        stop
      endif
      
      
    !--! Zero out the phase space distributions and set variables related to
       !   tracking thermal particles
      psd(:,:,:)           = 1.d-99
      esc_psd_feb_UpS(:,:) = 1.d-99
      esc_psd_feb_DwS(:,:) = 1.d-99
      n_cr_count       = 0
      num_crossings(:) = 0
      therm_grid(:) = 0
      therm_pxsk(:) = 0.d0
      therm_ptsk(:) = 0.d0
      therm_xwt(:)  = 0.d0
      
      ! In addition to initializing the phase space distribution, open the
      !  scratch (i.e. temporary) file to which we will write information
      !  about thermal particle grid crossings
      nc_unit = 2525
      do_therm_pt_formatted = .false.   !DOLATER: add keyword for this?
      if( do_therm_pt_formatted ) then
        open(unit=nc_unit,status="replace",form="formatted",&!&
             file="mc_crossings.dat")
      else
        open(unit=nc_unit,status="scratch",form="unformatted")
      endif
      
      
    !--! To maintain identical results between OpenMP and serial 
       !   versions, set RNG seed based on current iteration/ion/pcut/
       !   particle number
      iseed_mod =  (i_itr - 1)*n_ions  +  (i_ion - 1)
      call kisset(iseed_in -   iseed_mod, iseed_in - 2*iseed_mod,     &!&
                  iseed_in - 3*iseed_mod, iseed_in - 4*iseed_mod)
      
      
    !--! Initialize the particle populations that will be propagated through
       !   the shock structure
      call init_pop(do_fast_push, i_in_distr, i_ion, aa, n_pts_use,       &!&
         xwt_in, ptpf_in, pbpf_in, x_PT_cm_in, i_grid_in, pxx_flux,       &!&
         pxz_flux, en_flux)
      
      
    !--! Assign the various particle properties to the population
      xwt_new(1:n_pts_use)     = xwt_in(1:n_pts_use)
      ptpf_new(1:n_pts_use)    = ptpf_in(1:n_pts_use)
      pbpf_new(1:n_pts_use)    = pbpf_in(1:n_pts_use)
      x_PT_cm_new(1:n_pts_use) = x_PT_cm_in(1:n_pts_use)
      i_grid_new(1:n_pts_use)  = i_grid_in(1:n_pts_use)
      
      l_DwS_new(1:n_pts_use)       = .false.
      l_inj_new(1:n_pts_use)       = .false.
      xn_per_new(1:n_pts_use)      = xn_per_fine
      prp_x_cm_new(1:n_pts_use)    = x_grid_stop
      acctime_sec_new(1:n_pts_use) = 0.d0
      tcut_new(1:n_pts_use)        = 1
      
      do i_prt = 1, n_pts_use
        call kiss(rand)
        phi_rad_new(i_prt) = 2.d0*pii*rand
      enddo
      
      wt_running = xwt_in(1)  ! Weight of remaining particles, printed after
                              !   each pcut; note that this will not be
                              !   correct for all particles if they were
                              !   originally created so each thermal bin
                              !   would have equal weight
      
      
    !--! When using OpenMP, the array en_Xfer_pool can't be conditionally
       !   assigned shared or reduction status, so it can't be used for both
       !   the ion and electron loops.  To get around this, use one array to
       !   hold the donated energy, and another to hold the received energy.
      en_recv_pool(:) = en_Xfer_pool(:)
      
      
    !--! The array of pcuts read in by data_input has units momentum/mc.
       !   Convert to momentum for this species
      pcuts_use(1:n_pcuts) = pcuts_in(1:n_pcuts)  *  aa*xmp*ccgs
      
      !**********************************************************************
      !  Start of loop over pcuts
      !
      !  Particle splitting and resetting of n_pts_use is handled by the call
      !    to "new_pcut" at the end of each iteration of loop_pcut.
      !**********************************************************************
      loop_pcut: do i_cut = 1, n_pcuts
        
        
        ! Initialize all of the *_sav arrays to help prevent bleedover
        !   between pcuts or ion species
        l_save(:) = .false.  ! Whole array must be initialized in case number
                             !   of particles changes from pcut to pcut
        xwt_sav(:)         = 0.d0
        ptpf_sav(:)        = 0.d0
        pbpf_sav(:)        = 0.d0
        x_PT_cm_sav(:)     = 0.d0
        i_grid_sav(:)      = 0
        l_DwS_sav(:)       = .false.
        l_inj_sav(:)       = .false.
        xn_per_sav(:)      = 0.d0
        prp_x_cm_sav(:)    = 0.d0
        acctime_sec_sav(:) = 0.d0
        phi_rad_sav(:)     = 0.d0
        tcut_sav(:)        = 0
        
        ! A separate variable tracks the number of finished particles, so
        !   that race conditions can be avoided in OMP mode
        i_fin = 0
        
        ! For high-energy electrons in a strong magnetic field, need to know
        !   previous cutoff momentum for calculating new PRP downstream
        if( i_cut .gt. 1 ) pcut_prev = pcuts_use(i_cut-1)
          
        
        !********************************************************************
        !  Start of loop over particles
        !
        !  Quick explanation of particle properties.  All units cgs where
        !    a unit exists.
        !
        !  xwt: weighting value (used in momentum splitting)
        !  ptpf: total plasma frame momentum
        !  pbpf: momentum along B field in plasma frame.  NOT along x-axis
        !    unless upstream orientation is parallel
        !  x_PT_cm: current particle position
        !  i_grid: current grid zone number
        !  l_DwS: whether particle has been downstream
        !  l_inj: whether particle has been back upstream (i.e. is a CR)
        !  xn_per: time steps per gyro period, delta t = T_g/xn_per
        !  prp_x_cm: DwS position of PRP; adjusted to allow all particles,
        !    regardless of momentum, to isotropize before reaching
        !  acctime_sec: time since crossing shock for first time; not
        !    started until l_DwS = .true.
        !  phi_rad: phase angle of gyration rel. to z axis
        !  tcut_curr: next time at which particle tracking takes place
        !********************************************************************
!$omp parallel do default(none), schedule(dynamic,1), num_threads(4),     &!&
!
! OMP NOTE: Parameters and true constants don't need to be explicitly
!           listed below.  They are shared by default.
!
!$omp shared(n_pts_use, n_pts_max, i_itr, i_ion, o_o_mc, zz, aa,          &!&
!$omp        pmax_cutoff, rad_loss_fac, pcuts_use, i_cut, i_fin,          &!&
!$omp        n_print_pt, B_CMBz, tcuts),                                  &!&
!  Shared control variables
!$omp shared(n_ions, n_pcuts, iseed_in, l_sc_elec, en_xfer_frac, i_shock, &!&
!$omp        elec_wt_fac, dont_scatter, feb_UpS, age_max, do_rad_losses,  &!&
!$omp        dont_DSA, u_2, eta_mfp, nc_unit, do_therm_pt_formatted,      &!&
!$omp        i_grid_feb, feb_DwS, p_elec_crit, gam_elec_crit, pcut_prev,  &!&
!$omp        do_tcuts, inj_fracs, use_custom_epsB),                       &!&
!  Shared particle arrays
!$omp shared(xwt_new, ptpf_new, pbpf_new, x_PT_cm_new, i_grid_new,        &!&
!$omp        l_DwS_new, l_inj_new, xn_per_new, prp_x_cm_new,              &!&
!$omp        acctime_sec_new, phi_rad_new, tcut_new),                     &!&
!$omp shared(l_save, xwt_sav, ptpf_sav, pbpf_sav, x_PT_cm_sav, i_grid_sav,&!&
!$omp        l_DwS_sav, l_inj_sav, xn_per_sav, prp_x_cm_sav,              &!&
!$omp        acctime_sec_sav, phi_rad_sav, tcut_sav),                     &!&
!  Shared grid variables
!$omp shared(n_grid, btot_grid, uxsk_grid, uzsk_grid, utot_grid,          &!&
!$omp        gam_sf_grid, gam_ef_grid, beta_ef_grid, theta_grid,          &!&
!$omp        eps_target, en_Xfer_pool, en_recv_pool, x_grid_stop),        &!&
!  Shared flux/particle tracking variables
!$omp shared(n_cr_count, psd, esc_psd_feb_DwS, esc_psd_feb_UpS,           &!&
!$omp        pxx_flux, pxz_flux, en_flux, num_crossings, en_esc_UpS,      &!&
!$omp        px_esc_UpS, sum_press_DwS, sum_KEden_DwS, esc_flux,          &!&
!$omp        px_esc_feb, en_esc_feb, esc_en_eff, esc_num_eff, spec_sf,    &!&
!$omp        spec_pf),                                                    &!&
!  Private variables that are properties of the specific particle
!$omp private(iseed_mod, xwt, ptpf, pbpf, x_PT_cm, i_grid, i_grid_old,    &!&
!$omp         l_DwS, l_inj, xn_per, prp_x_cm, acctime_sec, phi_rad,       &!&
!$omp         y_PT_cm, z_PT_cm, gam_ptpf, gyro_denom, gyro_rad_tot_cm,    &!&
!$omp         gyro_period_sec, p_perp_b_pf, gyro_rad_cm, ptsk, pxsk,      &!&
!$omp         pysk, pzsk, pbsk, p_perp_b_sk, gam_ptsk, x_PT_old, i_start, &!&
!$omp         i_stop, n_split, ptpf_old, y_PT_old, z_PT_old, phi_rad_old, &!&
!$omp         tcut_curr),                                                 &!&
!  Private properties related to the particle's position in the bulk flow
!$omp private(uxsk, uzsk, utot, gam_usf, gam_uef, beta_uef, bmag, btht,   &!&
!$omp         b_sin_th, b_cos_th, uxsk_old, uzsk_old, utot_old,           &!&
!$omp         gam_usf_old, b_sin_old, b_cos_old),                         &!&
!  Additional private variables needed for the loop
!$omp private(helix_count, keep_looping, first_iter, i_return, i_reason,  &!&
!$omp         frac_to_keep, gampf_i, gampf_f, ptpf_f, scale_fac,          &!&
!$omp         en_to_Xfer, B_CMB_loc, B_tot_sq, dp_synch, t_step,          &!&
!$omp         x_move_bpar, rand, vpt_o_u2, L_diff, do_prob_ret,           &!&
!$omp         gyro_rad_tmp, lose_pt, vel, loop_again)
        loop_pt: do i_prt = 1, n_pts_use
          
        !--! To maintain identical results between OpenMP and serial 
           !   versions, set RNG seed based on current iteration/ion/pcut/
           !   particle number
          iseed_mod =  (i_itr - 1)*n_ions*n_pcuts*n_pts_max               &!&
                     + (i_ion - 1)       *n_pcuts*n_pts_max               &!&
                     + (i_cut - 1)               *n_pts_max               &!&
                     +  i_prt
          call kisset(iseed_in +   iseed_mod, iseed_in + 2*iseed_mod,     &!&
                      iseed_in + 3*iseed_mod, iseed_in + 4*iseed_mod)
          
          
        !--! Reset the counter for number of times through the main loop
          helix_count = 0
          
    
        !--! Get the properties of the particle we're about to treat
          xwt         = xwt_new(i_prt)
          ptpf        = ptpf_new(i_prt)
          pbpf        = pbpf_new(i_prt)
          x_PT_cm     = x_PT_cm_new(i_prt)
          i_grid      = i_grid_new(i_prt)
          i_grid_old  = i_grid  ! Needed for energy transfer
          l_DwS       = l_DwS_new(i_prt)
          l_inj       = l_inj_new(i_prt)
          xn_per      = xn_per_new(i_prt)
          prp_x_cm    = prp_x_cm_new(i_prt)
          acctime_sec = acctime_sec_new(i_prt)
          phi_rad     = phi_rad_new(i_prt)
          tcut_curr   = tcut_new(i_prt)
          
          y_PT_cm     = 0.d0  ! Not currently tracked, but could be added
          z_PT_cm     = 0.d0  !   in at later date
          
          gam_ptpf    = sqrt( 1.d0  +  (ptpf * o_o_mc)**2 )
          
          ! Constant that will be used repeatedly during loop
          gyro_denom = 1.d0 / (zz*qcgs * btot_grid(i_grid))
          if( use_custom_epsB .and. (x_PT_cm .gt. x_grid_stop) ) then
            gyro_denom = gyro_denom * ( x_PT_cm / x_grid_stop )**0.25d0
          endif
          
          ! Gyroradius assuming all motion is perpendicular to B field;
          !   pitch-angle-correct gyroradius is gyro_rad_cm
          gyro_rad_tot_cm = ptpf * ccgs * gyro_denom
          
          ! Gyroperiod in seconds
          gyro_period_sec = twopi * gam_ptpf * aa*xmp * ccgs * gyro_denom
          
          
        !--! Get the properties of the grid zone the particle's in
          uxsk     = uxsk_grid(i_grid)
          uzsk     = uzsk_grid(i_grid)
          utot     = utot_grid(i_grid)
          gam_usf  = gam_sf_grid(i_grid)
          gam_uef  = gam_ef_grid(i_grid)
          beta_uef = beta_ef_grid(i_grid)
          bmag     = btot_grid(i_grid)
          btht     = theta_grid(i_grid)
          
          b_sin_th = sin(btht)
          b_cos_th = cos(btht)
          
          
          !------------------------------------------------------------------
          ! Helix loop: The following loop is a large-scale restructuring of
          !   the original code's loops 5702 and 2001.
          ! The original code was slightly clearer in intent, but the goto
          !   statements made it impossible to parallelize for CPU computing;
          !   this code has sacrificed some of the clarity in exchange for
          !   potential speed gains down the road.
          !
          !    ORIGINAL                       THIS CODE
          !  5702 continue                 first_iter = .true.
          !
          !   Code Block 1                 do while( keep_looping )
          !                                  if(first_iter OR condition 1)
          !  2001 continue                     first_iter = .false.
          !                                    Code Block 1
          !   Code Block 2                   else
          !                                    Code Block 3
          !  if(condition 1) goto 5702       endif
          !
          !   Code Block 3                   Code Block 2
          !                                enddo
          !  goto 2001
          !
          ! Code Block 1: start of helix loop; minor setup
          ! Code Block 2: particle movement; check on DwS PRP return
          ! Code Block 3: grid zone changes; non-DwS escape (exception: if no
          !   scattering); radiative losses; momentum splitting; scattering
          !
          ! loop_helix in this code is the old loop 2001.  The first part of
          !   loop 5702 has been subsumed into the if statement, which
          !   determines whether Code Block 3 would be run or whether the
          !   code would cycle back to run Code Block 1 again.
          ! Following every combination of logical possibilities should
          !   result in the same path through the two loops.
          !------------------------------------------------------------------
          keep_looping = .true.
          first_iter   = .true.
          i_return     = -1     ! Control variable for "condition 1":
                                !  -1: default
                                !   0: particle escapes
                                !   1: particle returns from DwS via convection
                                !   2: particle didn't enter return calcs
          i_reason   = 0
          lose_pt    = .false.
          
          loop_helix: do while( keep_looping )
            
            ! Track number of times through the main loop.  This will only
            !   be needed for electrons at high energies when using
            !   radiative losses
            helix_count = helix_count + 1
            
            
            if( first_iter .or. (i_return .eq. 1) ) then
            
              first_iter = .false.
              
              
              ! Code Block 1: start of helix loop; minor setup
              !--------------------------------------------------------------
              ! Calculate momentum perpendicular to magnetic field...
              if( ptpf .lt. abs(pbpf) ) then
                p_perp_b_pf = 1.d-6 * ptpf
                pbpf        = sign( sqrt( ptpf**2 - p_perp_b_pf**2 ), pbpf )
                !CHECKTHIS: does this *ever* happen?!
                write(*,"(A)") 'Warning: ptpf < pbpf at top of loop_helix'
              else
                p_perp_b_pf = sqrt( ptpf**2  -  pbpf**2 )
              endif
              
              ! ...and turn that into a gyroradius
              gyro_rad_cm = p_perp_b_pf * ccgs * gyro_denom
              !--------------------------------------------------------------
              ! End of Code Block 1
              
              
            else
            
              
              ! Code Block 3: grid zone changes; non-DwS escape (exception:
              !   if no scattering); radiative losses; momentum splitting;
              !   scattering
              !--------------------------------------------------------------
            !--! Store values from the previous run through loop_helix; only
               !   values needed by subroutine Xform_p_psp are kept here
              uxsk_old    = uxsk
              uzsk_old    = uzsk
              utot_old    = utot
              gam_usf_old = gam_usf
              b_sin_old   = b_sin_th
              b_cos_old   = b_cos_th
              
              
            !--! Get new values for this run through loop_helix
              uxsk     = uxsk_grid(i_grid)
              uzsk     = uzsk_grid(i_grid)
              utot     = utot_grid(i_grid)
              gam_usf  = gam_sf_grid(i_grid)
              gam_uef  = gam_ef_grid(i_grid)
              beta_uef = beta_ef_grid(i_grid)
              bmag     = btot_grid(i_grid)
              btht     = theta_grid(i_grid)
              b_sin_th = sin(btht)
              b_cos_th = cos(btht)
              
              if( use_custom_epsB .and. (x_PT_cm .gt. x_grid_stop) ) then
                bmag = btot_grid(n_grid) * ( x_grid_stop / x_PT_cm )**0.25
              endif
              
              gyro_denom = 1.d0 / (zz*qcgs * bmag)
              
              
            !--! If particle crossed a velocity gradient, find its new shock
               !   frame properties
              if( uxsk .ne. uxsk_old ) then
                call Xform_p_PSP(aa, ptpf, pbpf, p_perp_b_pf, gam_ptpf,   &!&
                   phi_rad, uxsk_old, uzsk_old, utot_old, gam_usf_old,    &!&
                   b_cos_old, b_sin_old, uxsk, uzsk, utot, gam_usf,       &!&
                   b_cos_th, b_sin_th, ptsk, pxsk, pysk, pzsk, pbsk,      &!&
                   p_perp_b_sk, gam_ptsk)
              
                gyro_rad_cm = p_perp_b_pf * ccgs * gyro_denom
                gyro_rad_tot_cm = ptpf * ccgs * gyro_denom
              endif
              
              
          !--! Energy transfer from ions to electrons in the UpS region; only
             !   bother with this if
             !   (1) Electrons are a separate species and energy transfer is
             !     enabled
             !   (2) Particles are still on their first trip to the DwS
             !     region of the shock structure
             !   (3) Particles have entered a new grid zone, and energy
             !     should be either subtracted or added
            if( l_sc_elec .and. (en_xfer_frac .gt. 0.d0) ) then
            
              if( (.not. l_inj) .and. (x_PT_old .le. 0.d0) ) then
              
                if( i_grid_old .ne. i_grid ) then
                  
                  i_start = i_grid_old
                  i_stop  = min(i_grid, i_shock)
                  
                  ! Subtract energy from the ions and add it to the pool of
                  !   energy for this grid zone.
                  if( (aa .ge. 1.d0) .and.                                &!&
                      (maxval(eps_target(i_start+1:i_stop)) .gt. 0.d0) ) then
                    
                    ! Subtract energy based on the difference between current
                    !   and previous values of eps_elec
                    frac_to_keep = (1.d0 - eps_target(i_stop))            &!&
                                  / (1.d0 - eps_target(i_start))
                    gampf_i = sqrt( 1.d0  +  (ptpf*o_o_mc)**2 )
                    gampf_f = frac_to_keep*(gampf_i - 1.d0)  +  1.d0
                    
                    ! Split the energy equally among all grid zones crossed
                    !   during this scattering step; the donated energy is
                    !   weighted by xwt.
                    n_split = count( eps_target(i_start+1:i_stop) .gt. 0.d0 )
!$omp critical
                    do i = i_start+1, i_stop
                      if( eps_target(i) .gt. 0.d0 ) then
                        en_Xfer_pool(i) =  en_Xfer_pool(i)                &!&
                                     +  (gampf_i - gampf_f) * aa          &!&
                                       * rm_prot * xwt / float( n_split )
                      endif
                    enddo
!$omp end critical
                    
                    ! Calculate the new momentum based on the new energy, and
                    !   rescale components accordingly
                    ptpf_f    = aa*xmp * ccgs * sqrt( gampf_f**2 - 1.d0 )
                    scale_fac = ptpf_f / ptpf
                   
                    pbpf = pbpf * scale_fac
                    p_perp_b_pf = p_perp_b_pf * scale_fac
                    ptpf = ptpf_f
                    gam_ptpf = gampf_f
                    
                  else if( maxval(en_recv_pool(i_start+1:i_stop)) .gt. 0.d0 )then
                    
                    ! For electrons, add pooled energy.  Include energy from
                    !   all cells electron crossed in this scattering step.
                    ! Also modify the amount of energy to reflect the number
                    !   of electrons each MC particle represents
                    en_to_Xfer = sum( en_recv_pool(i_start+1:i_stop) )
                    en_to_Xfer = en_to_Xfer * elec_wt_fac
                    
                    gampf_i = sqrt( 1.d0  +  (ptpf*o_o_mc)**2 )
                    gampf_f = gampf_i  +  en_to_Xfer / (aa*rm_prot)
                    
                    ! Calculate the new momentum based on the new energy, and
                    !   rescale components accordingly
                    ptpf_f    = aa*xmp * ccgs * sqrt( gampf_f**2 - 1.d0 )
                    scale_fac = ptpf_f / ptpf
                    
                    pbpf = pbpf * scale_fac
                    p_perp_b_pf = p_perp_b_pf * scale_fac
                    ptpf = ptpf_f
                    gam_ptpf = gampf_f
                    
                  endif  ! check on particle species
                
                !--! Since the plasma-frame momenta have changed, recalculate
                   !   the shock-frame momenta                 
                  call Xform_p_PS(aa, pbpf, p_perp_b_pf, gam_ptpf,        &!&
                     phi_rad, uxsk, uzsk, utot, gam_usf, b_cos_th,        &!&
                     b_sin_th, ptsk, pxsk, pzsk, gam_ptsk)
                  
                endif  ! check on grid zone
                
              endif  ! check on whether particles are thermal and UpS
            
            endif  ! check on l_sc_elec and en_xfer_frac
            
            
            !--! Particle escape: DwS with scattering disabled
              if( dont_scatter .and. (x_PT_cm .gt. 10.d0*gyro_rad_cm) ) then
                i_return = 0
                i_reason = 1
                
                keep_looping = .false.
                cycle loop_helix
              endif
              
              
            !--! Particle escape: pmax (note that effects are same as for
               !   escape at UpS FEB)
              if( ptpf .gt. pmax_cutoff ) then
                ! Transform plasma frame momentum into shock frame to test
                !   there also
                call Xform_p_PS(aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad, &!&
                   uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th, ptsk,   &!&
                   pxsk, pzsk, gam_ptsk)
                
                if( ptsk .gt. pmax_cutoff ) then
                  i_reason = 2
                
                  keep_looping = .false.
                  cycle loop_helix
                endif  ! check on ptsk
              endif  ! check on ptpf
              
              
            !--! Particle escape: UpS FEB (note that effects are same as for
               !   escape by pmax)
              if( l_inj .and. (x_PT_cm .lt. feb_UpS) ) then
                i_reason = 2
              
                keep_looping = .false.
                cycle loop_helix
              endif
              
              
            !--! Particle escape: age_max
              if( (age_max .gt. 0.d0) .and. (acctime_sec .gt. age_max) ) then
                i_reason = 3
                
                keep_looping = .false.
                cycle loop_helix
              endif
              
              
            !--! Particle escape: transverse motion
               !DOLATER: implement new keyword related to transverse motion
!comm               if( yz_pos_max .gt. 0.d0 ) then
!comm                 trans_dist = sqrt( y_PT**2  +  z_PT**2 )
!comm                 
!comm                 if( trans_dist .gt. yz_pos_max ) then
!comm                   i_reason = 2
!comm                   
!comm                   keep_looping = .false.
!comm                   cycle loop_helix                
!comm                 endif
!comm               endif
              
              
            !--! Synchrotron/ICCMB losses for electrons
              if( do_rad_losses .and. (aa .lt. 1.d0) ) then
                
                ! Store previous ptpf value
                ptpf_old = ptpf
                
                ! Compute effective magnetic field for radiative losses. When
                !   squared to get energy density, one factor of gam_uef in
                !   B_CMB_loc represents increased number density, while the
                !   other is increased energy per photon.
                B_CMB_loc = B_CMBz * gam_uef
                
                B_tot_sq = bmag**2 + B_CMB_loc**2
                
                
                ! Note that here dp_synch is actually dp/p.  If this value is
                !   too large we will directly integrate from p_i to get p_f,
                !   since the discrete approach would result in too high a
                !   loss in a single time step
                dp_synch = rad_loss_fac * B_tot_sq * ptpf * t_step
                
                ! Correction to make sure electrons don't lose too much
                !   energy in a single time step
                if( dp_synch .gt. 1.d-2 ) then
                  ptpf = ptpf / (1.d0 + dp_synch)
                else
                  ptpf = ptpf * (1.d0 - dp_synch) ! Put second factor of ptpf
                                                  !   back into dp_synch
                endif
                
                ! Catch electrons that have somehow lost all their energy in
                !   a single time step
                if(ptpf .le. 0.d0) then
                  ptpf = 1.d-99
                  pbpf = 1.d-99
                  p_perp_b_pf = 1.d-99
                  gam_ptpf = 1.d0
                  
                  i_reason = 4
                  
                  keep_looping = .false.
                  cycle loop_helix
                endif
                
                ! Recalculate gam_ptpf since that's needed elsewhere
                gam_ptpf = sqrt( (ptpf*o_o_mc)**2  +  1.d0 )
                
                ! Modify components of ptpf due to losses
                pbpf        = pbpf * ptpf/ptpf_old
                p_perp_b_pf = p_perp_b_pf * ptpf/ptpf_old
                
                ! Also recalculate gyroradii
                gyro_rad_tot_cm = ptpf * ccgs * gyro_denom
                gyro_rad_cm     = p_perp_b_pf * ccgs * gyro_denom
              endif
              
              
            !--!
            !--! -- SCATTERING -- SCATTERING -- SCATTERING -- SCATTERING --
            !--!
              if( dont_scatter ) then
                ! Do nothing
              else
                call scattering(aa, gyro_denom, ptpf, gam_ptpf, xn_per,   &!&
                   pbpf, p_perp_b_pf, phi_rad, gyro_period_sec)
              endif
              
              
            !--! Update acceleration time in explosion frame, so convert
               !   t_step from plasma frame
               ! Only start the clock once a particle has crossed the shock
               !   for the first time
              if( l_DwS ) then
                acctime_sec = acctime_sec + t_step*gam_uef
                
                if( do_tcuts ) then
                  if( acctime_sec .ge. tcuts(tcut_curr) ) then
                    call tcut_track(tcut_curr, xwt, ptpf)
                    
                    tcut_curr = tcut_curr + 1
                  endif
                endif
              endif
              
              
            !--! Remove ions at splitting momentum
              if( (ptpf .gt. pcuts_use(i_cut)) .and. l_DwS ) then
                l_save(i_prt) = .true.
                
                xwt_sav(i_prt)         = xwt
                ptpf_sav(i_prt)        = ptpf
                pbpf_sav(i_prt)        = pbpf
                x_PT_cm_sav(i_prt)     = x_PT_cm
                i_grid_sav(i_prt)      = i_grid
                l_DwS_sav(i_prt)       = l_DwS
                l_inj_sav(i_prt)       = l_inj
                xn_per_sav(i_prt)      = xn_per
                if( x_PT_cm .lt. prp_x_cm ) then   ! Ensure particle is
                  prp_x_cm_sav(i_prt)  = prp_x_cm   !   within PRP
                else
                  prp_x_cm_sav(i_prt)  = x_PT_cm * 1.1d0
                endif
                acctime_sec_sav(i_prt) = acctime_sec
                phi_rad_sav(i_prt)     = phi_rad
                tcut_sav(i_prt)        = tcut_curr
                
                keep_looping = .false.
                cycle loop_helix
              endif
              
              
            !--! Shift between coarse/fine xn_per; right now the code uses
               !   only one value of xn_per everywhere in the ptot/x_PT phase
               !   space, so this subroutine does nothing.
              call xn_per_shift(xn_per, x_PT_cm, gyro_rad_tot_cm)
              !--------------------------------------------------------------
              ! End of Code Block 3
              
              
            endif  ! check on i_return
            
            
            ! Code Block 2: particle movement; fluxes; DwS escape/return
            !---------------------------------------------------------------
          !--! Store old x/y/z/phi positions
            x_PT_old    = x_PT_cm
            y_PT_old    = y_PT_cm
            z_PT_old    = z_PT_cm
            phi_rad_old = phi_rad
            
            
          !--! Find time step
            t_step = gyro_period_sec / xn_per
            
            
            ! Odd little loop that only matters if DSA has been disabled per
            !   the input file
            !----------------------------------------------------------------
            loop_again = .true.
            loop_no_DSA: do while( loop_again )
            
            !--! Update position/phase angle.
               ! Phase angle is easy: add appropriate amount and make sure
               !   result is in [0,2pi)
               ! For position, use inverse Lorentz transformations where
               !   primed frame is the plasma frame, moving to the right
               !   along the x axis, and the unprimed frame is the
               !   (stationary) shock frame.  Take frames to be coincident at
               !   t = t' = 0.  Then
               !
               !      x  =  gam_usf * (x' + uxsk*t'),
               !
               !   where t' is t_step and x' is distance moved in plasma
               !   frame (i.e. x_move_bpar below).
               ! Remember to take gyration about magnetic field into account
              phi_rad = phi_rad  +  twopi/xn_per
              if( phi_rad .ge. twopi ) phi_rad = phi_rad - twopi
              if( phi_rad .lt. 0.d0  ) phi_rad = phi_rad + twopi
              
              x_move_bpar = pbpf * t_step / (gam_ptpf * aa*xmp)
              
              x_PT_cm = x_PT_old                                          &!&
                       +  gam_usf * ( x_move_bpar * b_cos_th              &!&
                                     -  gyro_rad_cm * b_sin_th            &!&
                                       * (cos(phi_rad) - cos(phi_rad_old))&!&
                                     + uxsk * t_step )
              
              ! Don't care about transverse motion unless it's being tracked
              !DOLATER: add new keyword for tracking transverse motion
!comm               if( yz_pos_max .gt. 0.d0 ) then
!comm                 
!comm                 y_PT_cm = y_PT_old  +  gyro_rad_cm * ( sin(phi_rad)       &!&
!comm                                                       - sin(phi_rad_old) )
!comm                 
!comm                 z_PT_cm = z_PT_old  +  gam_usf * ( x_move_bpar * b_sin_th   &!&
!comm                       - gyro_rad_cm*b_cos_th*(cos(phi_rad)-cos(phi_rad_old))&!&
!comm                                                   + uzsk * t_step )
!comm               endif
              
              
            !--!DOLATER: set flags for energy transfer here?  Search original
               !   code for "out of ions"
              
              
            !--! Reflect particles under two conditions:
               !   (1) they've been DwS and are crossing back UpS
               !   (2) either DSA is disabled or they fail an injection check
               ! Cycle again if particles get reflected
               if( (x_PT_cm .le. 0.d0) .and. (x_PT_old .gt. 0.d0) .and.   &!&
                   (.not. l_inj) ) then
                 if( dont_DSA .or. (inj_fracs(i_ion) .lt. 1.d0) ) then
                   
                   call kiss(rand)
                   if( dont_DSA .or. (rand .gt. inj_fracs(i_ion)) ) then
                     ! Reflect particles with negative pbpf; randomize phase
                     !   of the rest
                     if( pbpf .lt. 0.d0 ) then
                       pbpf = -1.d0 * pbpf
                     else
                       call kiss(rand)
                       phi_rad = rand * twopi
                     endif
                   else
                     loop_again = .false.
                   endif
                 
                 else
                   
                   loop_again = .false.
                   
                 endif
               else
               
                 loop_again = .false.
                 
               endif
            enddo loop_no_DSA
            !----------------------------------------------------------------
            ! DSA injection prevented if specified
            
            
          !--! If particle crosses shock going UpS -> DwS
            if( (x_PT_old .lt. 0.d0) .and. (x_PT_cm .ge. 0.d0) ) then
              
              l_DwS = .true.
              
              ! Ensure the downstream region is sufficiently long to allow
              !   particle to isotropize.  This simple equation is decidedly
              !   non-trivial, and comes from two assumptions:
              !  1) the particle's diffusion coefficient D may be described
              !   by D = 1/3 * eta_mfp * r_g * v_pt (by default; a different
              !   f(r_g) may be specified as desired in place of eta*r_g)
              !  2) the relation between the diffusion coefficient D and the
              !   diffusion length L is L = D/<u>, where <u> is the average
              !   speed of diffusion.  Assuming isotropic particles in the
              !   DwS frame, <u> = u_2 since the average thermal *velocity*
              !   of the population is 0.
              !DOLATER: include f(r_g) in place of eta*r_g to allow for
              !  arbitrary diffusion
              vpt_o_u2 = ptpf/(aa*xmp*gam_ptpf * u_2)
              L_diff = third * eta_mfp * gyro_rad_tot_cm * vpt_o_u2
              
              if( prp_x_cm .lt. L_diff ) prp_x_cm = L_diff
              
            endif
            
            
          !--! Test for injection into the acceleration process
            if( l_DwS .and. (x_PT_cm .lt. 0.d0) ) l_inj = .true.
            
            
          !--! Calculate fluxes due to this motion; also locates new grid
             !   zone
            call all_flux(i_prt, aa, pbpf, p_perp_b_pf, ptpf, gam_ptpf,   &!&
               phi_rad, xwt, i_grid_old, i_grid, uxsk, uzsk, utot,        &!&
               gam_usf, b_cos_th, b_sin_th, x_PT_cm, x_PT_old, l_inj,     &!&
               nc_unit, do_therm_pt_formatted, i_grid_feb, pxx_flux,      &!&
               pxz_flux, en_flux, en_esc_UpS, px_esc_UpS, spec_sf,        &!&
               spec_pf, n_cr_count, num_crossings, psd)
            
            
          !--! Downstream escape/return test; assume we'll call subroutine
             !   prob_return unless told otherwise by particle location/info
            do_prob_ret = .true.
            
            if( (feb_DwS .gt. 0.d0) .and. (x_PT_cm .gt. feb_DwS) ) then
              
              ! Particle flagged as escaping downstream
              i_return    = 0
              do_prob_ret = .false.
              
            else if( x_PT_cm .gt. 1.1d0*prp_x_cm ) then
            
              ! The following simple equation is decidedly non-trivial, and
              !   comes from two assumptions:
              !  1) the particle's diffusion coefficient D may be described
              !   by D = 1/3 * eta_mfp * r_g * v_pt (by default; a different
              !   f(r_g) may be specified as desired in place of eta*r_g)
              !  2) the relation between the diffusion coefficient D and the
              !   diffusion length L is L = D/<u>, where <u> is the average
              !   speed of diffusion.  Assuming isotropic particles in the
              !   DwS frame, <u> = u_2 since the average thermal *velocity*
              !   of the population is 0.
              !DOLATER: include f(r_g) in place of eta*r_g to allow for
              !  arbitrary diffusion
              if( (aa .lt. 1.d0) .and. (ptpf .lt. p_elec_crit) ) then
                vpt_o_u2     = p_elec_crit / (aa*xmp*gam_elec_crit * u_2)
                gyro_rad_tmp = p_elec_crit * ccgs * gyro_denom
                L_diff       = third * eta_mfp * gyro_rad_tmp * vpt_o_u2
              else
                vpt_o_u2 = ptpf / (aa*xmp*gam_ptpf * u_2)
                L_diff   = third * eta_mfp * gyro_rad_tot_cm * vpt_o_u2
              endif
              
              if( x_PT_cm .gt. (6.91d0 * L_diff) ) then
                i_return    = 0
                do_prob_ret = .false.
              endif
              
            endif
            
            if( do_prob_ret ) then
              call prob_return(rad_loss_fac, B_CMBz, x_PT_old, aa, zz,    &!&
                 gyro_denom, i_return, x_PT_cm, prp_x_cm, ptpf, gam_ptpf, &!&
                 pbpf, p_perp_b_pf, acctime_sec, phi_rad, lose_pt,        &!&
                 helix_count, pcut_prev, xwt, tcut_curr)
            endif
            
            ! If particle escaped DwS, handle final calculations here
            if( i_return .eq. 0 ) then
              
              if( (gam_ptpf - 1.d0) .lt. en_rel_pt ) then
                vel = ptpf / (aa*xmp)
              else
                vel = ptpf / (aa*xmp * gam_ptpf)
              endif
              
!$omp atomic
              sum_press_DwS = sum_press_DwS  +  third * ptpf * vel * xwt
!$omp atomic
              sum_KEden_DwS = sum_KEden_DwS  +  (gam_ptpf - 1.d0)         &!&
                                               * aa*rm_prot * xwt
              
              i_reason = 1
              if( lose_pt ) i_reason = 4
              
              keep_looping = .false.
              cycle loop_helix
            endif
            !----------------------------------------------------------------
            ! End of Code Block 2
            
          enddo loop_helix
          !------------------------------------------------------------------
          ! End of loop moving/tracking particles on/off grid
          
          
          if( .not. l_save(i_prt) ) then
            call particle_finish(aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad,&!&
              uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th, i_reason,    &!&
              xwt, esc_psd_feb_DwS, esc_psd_feb_UpS, esc_flux, px_esc_feb,&!&
              en_esc_feb, esc_en_eff, esc_num_eff)
          endif
          
          
        !--! Particle counting
          if( (i_cut .eq. 1) .and.                                        &!&
              ( (i_prt .eq. 1) .or. (mod(i_prt, n_print_pt) .eq. 0) ) ) then
            write(*,"(10X,A,I6)") 'particle = ', i_prt
          endif
          
!$omp atomic
          i_fin = i_fin + 1
          if( mod(i_fin, 16) .eq. 0 ) then
            call print_progress_bar(i_fin, n_pts_use)
          endif
          
              
        enddo loop_pt
!$omp end parallel do
        !********************************************************************
        ! Conclusion of particle loop
        
        
        n_saved = count( l_save )
        
        
        call get_time(t_end)
        run_time = t_end - t_start
        
        write(*,"(2(A,I2),A,I3,2ES10.2E2,2(A,I6),2ES10.2E2)")             &!&
           ' itr=', i_itr, ' ion=', i_ion, ' icut=', i_cut,               &!&
           pcuts_in(i_cut), pcuts_use(i_cut)/(xmp*ccgs), '  n_sav=',      &!&
           n_saved, '/', n_pts_use, wt_running, run_time
        write(9,"(2(A,I2),A,I3,2ES10.2E2,2(A,I6),2ES10.2E2)")             &!&
           ' itr=', i_itr, ' ion=', i_ion, ' icut=', i_cut,               &!&
           pcuts_in(i_cut), pcuts_use(i_cut)/(xmp*ccgs), '  n_sav=',      &!&
           n_saved, '/', n_pts_use, wt_running, run_time
        
        
      !--! If no particles saved, don't bother with remaining pcuts
        if( n_saved .eq. 0 ) exit loop_pcut
        
        
      !--! Prepare population for next pcut
        if( pcuts_use(i_cut) .lt. p_pcut_hi ) then
          n_pts_target = n_pts_pcut
        else
          n_pts_target = n_pts_pcut_hi
        endif
        
        call new_pcut(n_pts_target, n_saved, l_save, i_grid_sav,          &!&
           l_DwS_sav, l_inj_sav, xwt_sav, ptpf_sav, pbpf_sav, x_PT_cm_sav,&!&
           xn_per_sav, prp_x_cm_sav, acctime_sec_sav, phi_rad_sav,        &!&
           tcut_sav, n_pts_use, i_grid_new, l_DwS_new, l_inj_new,         &!&
           xwt_new, ptpf_new, pbpf_new, x_PT_cm_new, xn_per_new,          &!&
           prp_x_cm_new, acctime_sec_new, phi_rad_new, tcut_new, wt_running)
        
        
      enddo loop_pcut
      !**********************************************************************
      ! Conclusion of pcuts loop
      
      write(*,*)
      write(9,*)
      
      
    !--! Obtain the dN/dp arrays for the current species: one 2-D array for
       !   each grid zone.
      call get_normalized_dNdp(dNdp_therm, dNdp_therm_pvals, dNdp_cr,     &!&
         nc_unit, do_therm_pt_formatted, zone_pop)
      
      
    !--! Get pressure (both components) and kinetic energy density everywhere
       !   in the shock profile
      call thermo_calcs(num_crossings, n_cr_count, therm_grid, therm_pxsk,&!&
         therm_ptsk, therm_xwt, nc_unit, do_therm_pt_formatted, psd,      &!&
         zone_pop, press_psd_par, press_psd_perp, en_den_psd)
      
      
    !--! Transform just the electron PSD into the explosion frame, since we
       !   will need it shortly to calculate inverse Compton emission
      call get_dNdp_2D(nc_unit, do_therm_pt_formatted, zone_pop,          &!&
         d2N_dpdcos_ef)
      
      
    !--! Output the spectra associated with this ion species
      !DOLATER: spectrum_plot
      
      
    !--! Print out escaping particle population for this species
      call print_dNdp_esc(esc_psd_feb_UpS, esc_psd_feb_DwS)
      
      
    enddo loop_ion
    !DEBUGLINE
    do i = 1, n_grid
      write(888,"(2i4,15es11.3e2)") i_itr, i, &!&
!         x_grid_rg(i), &!&
         therm_en_den(i,1:n_ions)/zone_vol(i),&!&
         therm_en_den(i,n_ions)/max(sum(therm_en_den(i,1:n_ions)),1.d-99),&!&
         sum(therm_en_den(i,1:n_ions)/zone_vol(i)), &!&
         en_den(i,1:n_ions)/zone_vol(i),&!&
         en_den(i,n_ions)/max(sum(en_den(i,1:n_ions)),1.d-99),&!&
         sum(en_den(i,1:n_ions)/zone_vol(i))
    enddo
    call print_plot_vals(888)
    !************************************************************************
    ! Conclusion of species loop
    
    
  !--! Compute the escaping flux for this iteration
    px_esc_flux_UpS(i_itr) = px_esc_UpS / flux_px_UpS
    en_esc_flux_UpS(i_itr) = en_esc_UpS / flux_en_UpS
    
    
  !--! Compute the adiabatic index everywhere on grid now that all particles
     !   are accounted for.  First store value from previous iteration, since
     !   both the pre- and post-iteration adiabatic indices will be used in
     !   the profile smoothing subroutine smooth_grid
    if( i_itr .eq. 1) then
      do i = 1, n_grid
        if( x_grid_cm(i) .le. 0.d0 ) then
          gam_adiab_grid(i,1) = 5.d0*third  ! #assumecold
        else
          gam_adiab_grid(i,1) = gam_adiab_2_RH
        endif
      enddo
    else
      gam_adiab_grid(:,1) = gam_adiab_grid(:,2)
    endif
    
    gam_adiab_grid(:,2) = 1.d0  +  (press_psd_par(:)+press_psd_perp(:))   &!&
                                  / en_den_psd(:)
    
    where( en_den_psd(:) .eq. 1.d-99 ) gam_adiab_grid(:,2) = 1.d-99
    
    
  !--! Also compute the adiabatic index of particles that were lost DwS
    gam_adiab_DwS(i_itr) = 1.d0 + sum_press_DwS/sum_KEden_DwS


  !--! Calculate expected escaping fluxes, now that adiabatic index is known
     !   far DwS.  Also average them so that the smoothing subroutine treats
     ! calculated and actual escaping fluxes identically
    call q_esc_calcs(gam_adiab_DwS(i_itr), q_esc_cal_px(i_itr),           &!&
       q_esc_cal_en(i_itr))
    n_avg      = min(i_itr, 4)
    q_esc_cal_px_avg = 0.d0
    q_esc_cal_en_avg = 0.d0
    do i = 1, n_avg
      i_tmp = i_itr + 1 - i
      q_esc_cal_px_avg = q_esc_cal_px_avg  +  q_esc_cal_px(i_tmp)
      q_esc_cal_en_avg = q_esc_cal_en_avg  +  q_esc_cal_en(i_tmp)
    enddo
    q_esc_cal_px_avg = q_esc_cal_px_avg / float(n_avg)
    q_esc_cal_en_avg = q_esc_cal_en_avg / float(n_avg)
    
    
  !--! In runs with OpenMP, race conditions may cause roundoff error at the
     !   least significant decimal place of the **_flux variables.  This
     !   causes slight variations in the smoothing algorithm across different
     !   runs.  Correct for that here by rounding off the last two decimal
     !   places (the second place is a fudge factor).
    do i = 1, n_grid
      
      ! pxx_flux
      call round(pxx_flux(i), 13)
      
      ! pxz_flux
      ! Commented out because it's not necessary for smoothing parallel
      !   shocks
!comm       call round(pxz_flux(i), 13)
      
      ! en_flux
      call round(en_flux(i), 13)
      
    enddo
    
    
  !--! Output grid data for this iteration, and smooth the grid for the next
     !   iteration
    if( l_oblique ) then
      write(*,"(2A)") 'ERROR: grid smoothing not coded for oblique shocks.'
      write(*,"(2A)") 'Stopping program now.'
      stop
    else
      call smooth_grid_par(i_itr, i_shock, n_grid, x_grid_rg, x_grid_cm,  &!&
         gam_adiab_grid, uzsk_grid, theta_grid, press_psd_par,            &!&
         press_psd_perp, flux_px_UpS, flux_en_UpS, gam_adiab_2_RH,        &!&
         q_esc_cal_px_avg, q_esc_cal_en_avg, pxx_flux, en_flux, uxsk_grid,&!&
         gam_sf_grid, btot_grid, utot_grid, gam_ef_grid, beta_ef_grid,    &!&
         epsB_grid)
    endif
    
    
  !--! Compute average escaping flux over last four iterations and write to
     !   file; average over all iterations if i_itr < 4.
    n_avg      = min(i_itr, 4)
    px_esc_avg = 0.d0
    en_esc_avg = 0.d0
    do i = 1, n_avg
      i_tmp = i_itr + 1 - i
      px_esc_avg = px_esc_avg  +  px_esc_flux_UpS(i_tmp)
      en_esc_avg = en_esc_avg  +  en_esc_flux_UpS(i_tmp)
    enddo
    px_esc_avg = px_esc_avg / float(n_avg)
    en_esc_avg = en_esc_avg / float(n_avg)
    
    write(9,"(2A)") ' Parallel shock q_esc from Double et al (2004) ',    &!&
        'equations:'
    write(9,"(A,ES10.3E2)") '     Esc. energy flux/UpS    = ',            &!&
        q_esc_cal_en_avg
    write(9,"(A,ES10.3E2)") '     Esc. momentum flux/UpS  = ',            &!&
        q_esc_cal_px_avg
    en_esc_flux_UpS(i_itr) = max(en_esc_flux_UpS(i_itr), 1.0d-99)
    en_esc_avg             = max(en_esc_avg,             1.0d-99)
    px_esc_flux_UpS(i_itr) = max(px_esc_flux_UpS(i_itr), 1.0d-99)
    px_esc_avg             = max(px_esc_avg,             1.0d-99)
    write(9,"(A,I0,2(A,ES9.2E2))") ' Esc. en flux FEB/UpS  for i_itr = ', &!&
        i_itr, ':   en esc = ', en_esc_flux_UpS(i_itr),                   &!&
        '   Avg. esc en  = ', en_esc_avg
    write(9,"(A,I0,2(A,ES9.2E2))") ' Esc. pxx flux FEB/UpS for i_itr = ', &!&
        i_itr, ':  pxx esc = ', px_esc_flux_UpS(i_itr),                   &!&
        '   Avg. esc pxx = ', px_esc_avg
    
    if( q_esc_cal_px_avg .eq. 0.d0 ) then
      write(9,"(A)") ' Avg q_px_MC/q_px_cal N/A, b/c q_px_cal = 0'
    else
      write(9,"(A,ES9.2E2)") ' Avg q_px_MC/q_px_cal. = ',                 &!&
          px_esc_avg/q_esc_cal_px_avg
    endif
    if( q_esc_cal_en_avg .eq. 0.d0 ) then
      write(9,"(A)") ' Avg q_en_MC/q_en_cal N/A, b/c q_en_cal = 0'
    else
      write(9,"(A,ES9.2E2)") ' Avg q_en_MC/q_en_cal. = ',                 &!&
          en_esc_avg/q_esc_cal_en_avg
    endif
    write(9,*)
    
    
  !--! Compute various adiabatic indices and write them out
    n_avg      = min(i_itr, 4)
    gamma_DwS_esc = 0.d0
    do i = 1, n_avg
      i_tmp = i_itr + 1 - i
      
      gamma_DwS_esc = gamma_DwS_esc  +  gam_adiab_DwS(i_tmp)
    enddo
    gamma_DwS_esc = gamma_DwS_esc / float(n_avg)
    gamma_UpS     = 5.d0 * third ! #assumecold
    
    write(9,"(A,I0)") ' Iteration #', i_itr
    write(9,"(2(A,ES11.3E2))") '   r_comp = ', r_comp, '      r_RH = ',rRH
    write(9,"(A,ES11.3E2)") '   Adiab index for far UpS particles = ',    &!&
        gamma_UpS
    write(9,"(A,ES11.3E2)") '   Adiab index for DwS PRP particles = ',    &!&
        gamma_DwS_esc
    write(9,"(A,ES11.3E2)") '   Adiab index from R-H relations    = ',    &!&
        gam_adiab_2_RH
    write(9,*)
    
    
  !--! If tcut tracking was enabled, print out the results here
    if( do_tcuts ) call tcut_print()
    
    
  enddo loop_itr
  !**************************************************************************
  ! Conclusion of iteration loop
  
  
!--! End the run
  call get_time(t_end)
  
  run_time = t_end - t_start
  
  write(*,*)
  write(*,"(2(A,ES12.2E2),A)") ' Finished.  Run time = ', run_time,       &!&
      ' sec, ', run_time/60.d0, ' min'
  
  write(9,*)
  write(9,"(2(A,ES12.2E2),A)") ' Finished.  Run time = ', run_time,       &!&
      ' sec, ', run_time/60.d0, ' min'
  write(9,*)
  
stop
end program


!****************************************************************************
!****************************************************************************
subroutine all_flux(i_prt, aa, pbpf, p_perp_b_pf, ptpf, gam_ptpf, phi_rad,&!&
     xwt, i_grid_old, i_grid, uxsk, uzsk, utot, gam_usf, b_cos_th,        &!&
     b_sin_th, x_PT_cm, x_PT_old, l_inj, nc_unit, do_therm_pt_formatted,  &!&
     i_grid_feb, pxx_flux, pxz_flux, en_flux, en_esc_UpS, px_esc_UpS,     &!&
     spec_sf, spec_pf, n_cr_count, num_crossings, psd)

! Tracks particle flux due to motion on the grid.  Also finds number of new
!   grid zone
!
! Input arguments:
!#DOLATER
! Output arguments:
!#DOLATER
! Input/output arguments:
!#DOLATER

use constants
use parameters, only: na_g, psd_max, en_rel_pt, na_cr
use controls, only: n_xspec, x_spec, feb_UpS, l_oblique, gam_Z, u_Z
use grid_vars, only: n_grid, x_grid_cm
  ! The therm_*** arrays can be pulled from the module because they are
  !    inherently thread-safe.  Only n_cr_count and num_crossings need
  !    protection from race conditions, and so need to be included explicitly
  !    in the arguments of all_flux.
use species_vars, only: therm_grid, therm_pxsk, therm_ptsk, therm_xwt

implicit none

  ! Input arguments
integer, intent(in) :: i_prt, nc_unit, i_grid_feb
logical, intent(in) :: l_inj, do_therm_pt_formatted
real(kind=8), intent(in) :: aa, pbpf, p_perp_b_pf, ptpf, gam_ptpf,        &!&
     phi_rad, xwt, uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th, x_PT_cm,&!&
     x_PT_old
  ! Output arguments
  ! Input/output arguments
integer, intent(inout) :: i_grid, i_grid_old, n_cr_count
integer, dimension(na_g), intent(inout) :: num_crossings
real(kind=8), intent(inout) :: px_esc_UpS, en_esc_UpS
real(kind=8), dimension(na_g), intent(inout) :: pxx_flux, pxz_flux, en_flux
real(kind=8), dimension(0:psd_max, na_g), intent(inout) :: spec_sf, spec_pf
real(kind=8), dimension(0:psd_max, 0:psd_max, na_g), intent(inout) :: psd

  ! Local variables
real(kind=8), parameter :: spike_away = 1000.d0  ! Max value for 1/cosine
integer :: i_pt, j_th, i, i_pt_pf, j_th_pf
real(kind=8) :: ptsk, pxsk, pzsk, gam_ptsk, pt_o_px_sk, o_o_vxsk,         &!&
     o_o_pxsk, pt_o_px_pf, en_flux_add, flux_wt_fac


!--! Very early check to see if particle crossed a grid zone boundary
  i_grid_old = i_grid
  
  if( x_PT_cm .gt. x_PT_old ) then
  
  !--! Loop over grid zones until we find the next boundary the particle
     !   hasn't crossed
    do i = i_grid+1, n_grid+1
    !--! Did particle stop before reaching this boundary?
      if( x_grid_cm(i) .gt. x_PT_cm ) then  ! Next grid zone boundary is DwS
        i_grid = i - 1                      !   of particle location, so
        exit                                !   i_grid should be zone just
      endif                                 !   UpS
    enddo
  
  else
  
  !--! Loop over grid zones until we find the next boundary the particle
     !   hasn't crossed
    do i = i_grid, 1, -1
    
    !--! Did particle stop before reaching this boundary?
      if( x_grid_cm(i) .le. x_PT_cm ) then  ! Next grid zone boundary is UpS
        i_grid = i                          !   of particle location, so
        exit                                !   i_grid should be this zone
      endif
    enddo
  
  endif
  
  ! Don't bother with *any* of the other computations if
  !  (1) Grid zone hasn't changed, and
  !  (2) There's no possibility of crossing an intra-grid detector
  ! #DOLATER: remove extra return statement, folding rest of computation into
  !   "else" block of if-then
  if( (i_grid .eq. i_grid_old) .and. (i_grid .gt. i_grid_feb) .and.       &!&
      (n_xspec .eq. 0) ) then
    return
  endif
  
  
!--! Convert plasma frame momentum to shock frame; determine a few values 
   !   that will be reused during the call to all_flux
  call Xform_p_PS(aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad, uxsk, uzsk,   &!&
     utot, gam_usf, b_cos_th, b_sin_th, ptsk, pxsk, pzsk, gam_ptsk)
  
  if( ptsk .gt. abs(pxsk*spike_away) ) then
    pt_o_px_sk = spike_away
    o_o_vxsk   = spike_away/uxsk  ! Minimum shock frame velocity is a small
                                  !   fraction of local bulk flow speed
  else
    o_o_pxsk   = 1.d0 / pxsk
    pt_o_px_sk = ptsk * o_o_pxsk
    o_o_vxsk   = gam_ptsk * aa*xmp * o_o_pxsk
  endif
  
  if( ptpf .gt. abs(pbpf*spike_away) ) then
    pt_o_px_pf = spike_away
  else
    pt_o_px_pf = abs(ptpf / pbpf)
  endif
  
  ! Kinetic energy only; rest mass energy NOT included
  if( (gam_ptsk - 1.d0) .gt. en_rel_pt ) then
    en_flux_add = (gam_ptsk - 1.d0) * aa*rm_prot * xwt
  else
    en_flux_add = ptsk**2 / (2.d0 * aa*xmp) * xwt
  endif
  
  
!--! Calculate spectrum at x_spec locations if needed
  if( n_xspec .gt. 0 ) then
    call get_psd_bins(pxsk, ptsk, .false., i_pt,    j_th)
    call get_psd_bins(pbpf, ptpf, .false., i_pt_pf, j_th_pf)
        
    do i = 1, n_xspec
      if( ( (x_PT_old .lt. x_spec(i)) .and. (x_PT_cm .ge. x_spec(i)) ).or.&!&
          ( (x_PT_cm .le. x_spec(i)) .and. (x_PT_old .gt. x_spec(i)) ) ) then
        
      !--! Spectrum in shock frame
!$omp atomic
        spec_sf(i_pt, i) = spec_sf(i_pt,i)  +  xwt * pt_o_px_sk
        
        
        !DOLATER: if the shock is oblique, need to replace pbpf with pxpf
        !  below, as well as changing the outputs from Xform_p_PS to include
        !  pxpf at all
        if( l_oblique ) then
          write(*,"(2A)") 'ERROR in all_flux: cannot calculate ',         &!&
              'flux_wt_fac when pbpf != pxpf'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
      !--! Spectrum in plasma frame; flux_wt_fac corresponds to vxpf/vxsk and
         !   measures relative likelihood of crossing in the plasma frame
         !   given a known crossing in the shock frame
        flux_wt_fac = abs(pbpf/pxsk) * (gam_ptsk/gam_ptpf)
        
!$omp atomic
        spec_pf(i_pt_pf, i) = spec_pf(i_pt_pf, i)                         &!&
                             +  xwt * pt_o_px_pf * flux_wt_fac
      endif
    enddo
  endif
  
  
!--! Check if particle has already been injected into acceleration process;
   !   if it has get the PSD bins it will be placed into
  if( l_inj ) call get_psd_bins(pxsk, ptsk, .true., i_pt, j_th)
  
  
  ! Main loop for all_flux: depending on motion of particle, travel UpS or
  !   DwS and update fluxes across all zone boundaries the particle crossed
  ! WARNING: the particle quantity "xwt" is the fraction of the far UpS
  !   density each particle represents.  However, the actual flux is
  !        gam_Z * u_Z * den_Z,
  !   which means that the flux contribution of each particle must be
  !   increased by a factor of gam_Z*u_Z.
  !--------------------------------------------------------------------------
!--! Downstream first
  if( x_PT_cm .gt. x_PT_old ) then
  
  !--! Loop over grid zones boundaries that the particle has crossed
    do i = i_grid_old+1, i_grid
      
    !--! Update fluxes
!$omp atomic
      pxx_flux(i) = pxx_flux(i)  +      pxsk *xwt * gam_Z*u_z
!$omp atomic
      pxz_flux(i) = pxz_flux(i)  +  abs(pzsk)*xwt * gam_Z*u_Z
      
!$omp atomic
      en_flux(i)  = en_flux(i)   +  en_flux_add   * gam_Z*u_Z
      
      
    !--! Update PSD, or add to list of tracked thermal particles
      if( l_inj ) then
        
!$omp atomic
        psd(i_pt, j_th, i) = psd(i_pt, j_th, i)  +  xwt * abs(o_o_vxsk)
!$omp end atomic
        
      else
        
        ! Particle is a thermal particle.  Increment crossing counter and
        !   attempt to add entry to the various arrays.  If arrays already
        !   full, write out crossing data to scratch file for tracking them.
        !   Then update array holding number of crossings.
!$omp critical
        if( n_cr_count .lt. na_cr ) then
          n_cr_count = n_cr_count + 1
          therm_grid(n_cr_count) = i
          therm_pxsk(n_cr_count) = pxsk
          therm_ptsk(n_cr_count) = ptsk
          therm_xwt(n_cr_count)  = xwt * abs(o_o_vxsk)
        else
          ! Need to write to scratch file, formatted or otherwise
          if( do_therm_pt_formatted ) then
            write(nc_unit,"(I3,I5,4ES20.12E2)") i, i_prt, pxsk, ptsk,     &!&
                xwt * abs(o_o_vxsk)
          else
            write(nc_unit) i, i_prt, pxsk, ptsk, xwt*abs(o_o_vxsk)
          endif
          
        endif  ! check on n_cr_count
!$omp end critical
        
!$omp atomic
        num_crossings(i) = num_crossings(i) + 1
        
      endif  ! check on l_inj
      
    enddo  ! loop over grid zones
    
    
  else
!--! Particle has moved upstream
  
  !--! Loop over grid zones boundaries that the particle has crossed
    do i = i_grid_old, i_grid+1, -1
      
    !--! Is particle UpS of free escape boundary after crossing DwS?  Then
       !   Don't count flux contributions to grid zones UpS from FEB
      if( l_inj .and. (i .le. i_grid_feb) ) cycle
      
    !--! Update fluxes; note the minus signs force pxx_flux to increase (b/c
       !   particle is moving UpS and pxsk < 0) and force en_flux to decrease
!$omp atomic
      pxx_flux(i) = pxx_flux(i)  -      pxsk *xwt * gam_Z*u_Z
!$omp atomic
      pxz_flux(i) = pxz_flux(i)  +  abs(pzsk)*xwt * gam_Z*u_Z
      
!$omp atomic
      en_flux(i)  = en_flux(i)   -  en_flux_add   * gam_Z*u_Z
      
      
    !--! Update PSD, or add to list of tracked thermal particles
      if( l_inj ) then
        
!$omp atomic
        psd(i_pt, j_th, i) = psd(i_pt, j_th, i)  +  xwt * abs(o_o_vxsk)
!$omp end atomic
        
      else
        
        ! Particle is a thermal particle.  Increment crossing counter and
        !   attempt to add entry to the various arrays.  If arrays already
        !   full, write out crossing data to scratch file for tracking them.
        !   Then update array holding number of crossings.
!$omp critical
        if( n_cr_count .lt. na_cr ) then
          n_cr_count = n_cr_count + 1
          therm_grid(n_cr_count) = i
          therm_pxsk(n_cr_count) = pxsk
          therm_ptsk(n_cr_count) = ptsk
          therm_xwt(n_cr_count)  = xwt * abs(o_o_vxsk)
        else
          ! Need to write to scratch file, formatted or otherwise
          if( do_therm_pt_formatted ) then
            write(nc_unit,"(I3,I5,4ES20.12E2)") i, i_prt, pxsk, ptsk,      &!&
                xwt * abs(o_o_vxsk)
          else
            write(nc_unit) i, i_prt, pxsk, ptsk, xwt*abs(o_o_vxsk)
          endif
          
        endif  ! check on n_cr_count
!$omp end critical
        
!$omp atomic
        num_crossings(i) = num_crossings(i) + 1
        
      endif  ! check on l_inj
      
    enddo  ! loop over grid zones
    
  endif  ! check on direction of motion
  !--------------------------------------------------------------------------
  ! Finished updating fluxes for crossed grid boundaries  
  
  
!--! One final task: update tracker for escaping flux at FEB if needed; don't
   !   forget that flux must be rescaled
!$omp critical
  if( l_inj .and. (x_PT_cm .lt. feb_UpS) .and. (x_PT_old .ge. feb_UpS) ) then
    en_esc_UpS = en_esc_UpS  +  en_flux_add * gam_Z*u_Z
    px_esc_UpS = px_esc_UpS  -  pxsk * xwt  * gam_Z*u_Z
  endif
!$omp end critical

return
end subroutine all_flux


!****************************************************************************
!****************************************************************************
subroutine calc_DwS(l_oblique, bmag_Z, r_comp, beta_Z, beta_2, gam_2,     &!&
     bmag_2, theta_B2, theta_u2)
! Commented out argument list is for eventual oblique extension
!comm (l_oblique, theta_BZ, bmag_Z, r_comp, beta_Z, gam_Z,   &!&
!comm n_ions, aa_ion, zz_ion, denZ_ion, tZ_ion, l_sc_elec, tZ_elec, beta_2,&!&
!comm gam_2, bmag_2, theta_B2, theta_u2)

! Uses the Rankine-Hugoniot jump conditions to calculate the downstream
!   conditions for a test particle shock.  Big difference between this
!   subroutine and calc_rRH is that we already know what the DwS speed is,
!   courtesy of r_comp in the input.
!
! Input arguments:
!  1) l_oblique: whether the shock is oblique
!  2) theta_BZ: angle[deg] between UpS magnetic field and shock normal
!  3) bmag_Z: far UpS magnetic field strength[Gauss]
!  4) r_comp: compression ratio of shock
!  5) beta_Z: UpS bulk fluid speed, over c
!  6) gam_Z: Lorentz factor associated with beta_Z
!  7) n_ions: number of different ion species
!  8) aa_ion: array of particle species' atomic mass numbers
!  9) zz_ion: array of particle species' charge numbers
!  10) denZ_ion: array of far UpS densities for particle species
!  11) tZ_ion: array of particle species' far UpS temperatures
!  12) l_sc_elec: flag for whether electrons are a separate species
!  13) tZ_elec: if electrons are not a separate species, this is their far
!    UpS temperature
! Output arguments
!  1) beta_2: (total) bulk fluid speed DwS
!  2) gam_2: Lorentz factor associated with beta_2
!  3) bmag_2: DwS magnetic field strength[Gauss]
!  4) theta_B2: angle[deg] between DwS magnetic field and shock normal
!  5) theta_u2: angle[deg] between DwS fluid velocity and shock normal

use constants
use parameters, only: na_i, beta_rel_fl

implicit none

  ! Input arguments
logical, intent(in) :: l_oblique
real(kind=8), intent(in) :: bmag_Z, r_comp, beta_Z
  ! Output arguments
real(kind=8), intent(out) :: beta_2, gam_2, bmag_2, theta_B2, theta_u2

  ! Local variables
logical :: l_nonrel


  !--------------------------------------------------------------------------
  !  Four possibilities for R-H relations: nonrel/rel and parallel/oblique.
  !    Determine which of the four to use.  Cutoff for nonrel/rel is set in
  !    module 'controls'
  !--------------------------------------------------------------------------
  if( beta_Z .lt. beta_rel_fl ) then
    l_nonrel = .true.
  else
    l_nonrel = .false.
  endif
  
  
  !--------------------------------------------------------------------------
  !  Possibility 1: Parallel at any shock speed
  !--------------------------------------------------------------------------
  if( .not. l_oblique ) then
    
    beta_2   = beta_Z / r_comp
    gam_2    = 1.d0 / sqrt( 1.d0 - beta_2**2 )
    theta_u2 = 0.d0
    
    bmag_2   = bmag_Z
    theta_B2 = 0.d0
  
      
  !--------------------------------------------------------------------------
  !  Possibility 2: Oblique at any shock speed.  Not currently supported by
  !    code, but included here in case code is extended in future.
  !--------------------------------------------------------------------------
  else
    write(*,"(2A)") 'ERROR in calc_DwS: not implemented for oblique ',    &!&
        'shocks yet.'
    write(*,"(2A)") "If this ever changes, don't forget to update the",   &!&
        ' shock profile in subroutine "setup_profile".'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif

return
end subroutine calc_DwS

!****************************************************************************
!****************************************************************************
subroutine calc_rRH(beta_Z, gam_Z, n_ions, aa_ion, zz_ion, denZ_ion,      &!&
     tZ_ion, l_sc_elec, tZ_elec, rRH, l_oblique, gam_adiab_2_RH)

! Uses the Rankine-Hugoniot jump conditions to calculate the compression
!   ratio for a shock assuming test-particle conditions.  In other words,
!   (1) sharp shock, (2) negligible/no DSA, and (3) no escaping flux.
! Additionally assumes that the inflowing plasma has non-rel thermal speeds
!   to make UpS adiabatic index exactly 5/3.
!
! Input arguments:
!  1) beta_Z, gam_Z: shock speed and Lorentz factor
!  2) n_ions: number of different ion species
!  3) aa_ion: array of particle species' atomic mass numbers
!  4) zz_ion: array of particle species' charge numbers
!  4) denZ_ion: array of far UpS densities for particle species
!  5) tZ_ion: array of particle species' far UpS temperatures
!  6) l_sc_elec: flag for whether electrons are a separate species
!  7) tZ_elec: if electrons are not a separate species, this is their far
!    UpS temperature
!  8) l_oblique: controls whether to use parallel or oblique formulations of
!    R-H relations
! Output arguments:
!  1) rRH: Rankine-Hugoniot compression ratio
!  2) gam_adiab_2_RH: ratio of specific heats (adiabatic index) for DwS
!    region, assuming r_comp = rRH

use constants
use parameters, only: na_i, beta_rel_fl

implicit none

  ! Input arguments
integer, intent(in) :: n_ions
logical, intent(in) :: l_sc_elec, l_oblique
real(kind=8), intent(in) :: beta_Z, gam_Z, tZ_elec
real(kind=8), dimension(na_i), intent(in) :: aa_ion, zz_ion, denZ_ion, tZ_ion
  ! Output arguments
real(kind=8), intent(out) :: rRH, gam_adiab_2_RH

  ! Local variables
integer, parameter :: maxitrs = 10000
integer :: i, j
logical :: l_nonrel
real(kind=8), parameter :: target_err = 1.d-6
real(kind=8) :: press_Z, rho_Z, den_elec, gam_sph, c_s, M_Z, sqrt_term,   &!&
     w_Z, UpS_mom_flux, UpS_num_flux, gam_2_min, w_fac_max, rm_avg,       &!&
     p2_max_A, p2_max_B, p2_max_C, p2_max_sq, p2_max, p2_guess,           &!&
     p2_guess_p, del_p2_guess, p2_o_mc, P_fac, w_fac, gam_2, gambeta_2,   &!&
     F_p2_guess, F_p2_guess_p, Fprime_p2_guess, p2_guess_next, err_curr,  &!&
     p2_found, e_fac, phi_fac, beta_2
real(kind=8), dimension(na_i) :: rm_ion, den_rel_ion


  !--------------------------------------------------------------------------
  !  Four possibilities for R-H relations: nonrel/rel and parallel/oblique.
  !    Determine which of the four to use.  Cutoff for nonrel/rel is set in
  !    module 'controls'
  !--------------------------------------------------------------------------
  if( beta_Z .lt. beta_rel_fl ) then
    l_nonrel = .true.
  else
    l_nonrel = .false.
  endif
  
  
  !--------------------------------------------------------------------------
  !  Possibility 1: Nonrel, parallel
  !    Solution comes from Ellison (1985) [1985JGR....90...29E].  Uses far
  !    UpS Mach number to calculate rRH.
  !--------------------------------------------------------------------------
  if( l_nonrel  .and.  (.not. l_oblique) ) then
    
    ! Calculate thermal pressure of far upstream gas
    press_Z  = 0.d0
    rho_Z    = 0.d0
    den_elec = 0.d0
    do j = 1, n_ions
      press_Z  =  press_Z  +  denZ_ion(j) * xkb*tZ_ion(j)
      rho_Z    =  rho_Z    +  denZ_ion(j) * xmp*aa_ion(j)
      
      if( aa_ion(j) .ge. 1.d0 ) den_elec  =  den_elec + denZ_ion(j)*zz_ion(j)
    enddo
    
    ! If electrons were not a separate species, add them in here
    if( .not. l_sc_elec ) then
      press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
      rho_Z    =  rho_Z    +  den_elec * xme
    endif
    
    ! Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
    !   to calculate the far UpS sound speed and Mach number
    ! #assumecold
    gam_sph = 5.d0 / 3.d0
    c_s     = sqrt( gam_sph * press_Z / rho_Z )
    M_Z     = beta_Z * ccgs / c_s
    
    ! Finally, use Equation (11) from Ellison (1985) to calculate rRH.  Note
    !   that q = 0 here b/c we assume no escaping flux
    sqrt_term = 9.d0  *  ( 1.d0 - 1.d0/M_Z**2 )**2
    sqrt_term = sqrt( sqrt_term )
    
    rRH = 8.d0 / ( 5.d0 + 3.d0/M_Z**2 - sqrt_term )
    
    ! In non-rel case, downstream adiabatic index is pegged to 5/3
    gam_adiab_2_RH = 5.d0/3.d0
    
    
  !--------------------------------------------------------------------------
  !  Possibility 2: Relativistic, parallel
  !   Solution comes from Ellison+ (1990) [1991ApJ...378..214E].  Uses
  !     relativistic Rankine-Hugoniot relations.  See that paper for details
  !     of equations and associated quantities.  Briefly,
  !        R-H1:             g0 * n0 * b0  =  g2 * n2 * b2
  !        R-H2:  g0^2 * w0 * b0^2  +  P0  =  g2^2 * w2 * b2^2  +  P2
  !        R-H3:           g0^2 * w0 * b0  =  g2^2 * w2 * b2
  !     where
  !        w    = E_rm  +  E_ke  +  P,  <--- enthalpy as total energy density
  !                                           + pressure
  !        E_rm = n * m * c^2           <--- rest mass energy density
  !        E_ke = n * m * c^2 * e(p)    <---  kinetic  energy density, with
  !           e(p) =  sqrt( 1 + (p/mc)^2 )  -  1
  !        P    = 1/3 * n * p * v       <--- pressure
  !
  !   Assumes that downstream particle distributions are delta-functions.
  !     Solves for p2 using Newton's method, then works backwards to rRH.
  !--------------------------------------------------------------------------
  else if( (.not. l_nonrel)  .and.  (.not. l_oblique) ) then
    
    ! Calculate thermal pressure of far upstream gas
    press_Z  = 0.d0
    rho_Z    = 0.d0
    den_elec = 0.d0
    do j = 1, n_ions
      press_Z  =  press_Z  +  denZ_ion(j) * xkb*tZ_ion(j)
      rho_Z    =  rho_Z    +  denZ_ion(j) * xmp*aa_ion(j)
      
      if( aa_ion(j) .ge. 1.d0 ) den_elec = den_elec                       &!&
                                          +  denZ_ion(j) * zz_ion(j)
      
      ! Calculate two quantities to be used during loop to find rRH: the
      !   rest mass-energy of each species, and the density relative to
      !   protons
      rm_ion(j)      = xmp * aa_ion(j) * ccgs**2
      den_rel_ion(j) = denZ_ion(j) / denZ_ion(1)
    enddo
    
    ! If electrons were not a separate species, add them in here
    if( .not. l_sc_elec ) then
      press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
      rho_Z    =  rho_Z    +  den_elec * xme
    endif
    
    
    ! Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
    !   to calculate the far UpS enthalpy
    ! #assumecold
    gam_sph = 5.d0 / 3.d0
    w_Z     = rho_Z * ccgs**2  +  gam_sph/(gam_sph - 1.d0) * press_Z
    
    ! Calculate the far UpS momentum flux
    UpS_mom_flux = gam_Z**2 * w_Z * beta_Z**2  +  press_Z
    UpS_num_flux = gam_Z * denZ_ion(1) * beta_Z ! Protons only here; not
                                                !  strictly correct but
                                                !  appropriate for later use
    
    
    ! Now use Newton's method to determine the downstream momentum that
    !   satisfies the R-H relations.
    ! Assumptions: (1) momentum distribution functions are delta-functions
    !   rather than thermal, and (2) any non-proton species have p /propto m
    !------------------------------------------------------------------------
    ! We know a priori that the R-H compression ratio will be between 3.0
    !   and 4.0.  Use a compression ratio of 4.5 as an upper bound, which
    !   sets a minimum value for gam_2 and in turn an upper limit on the
    !   downstream momentum.
    gam_2_min = 1.d0 / sqrt( 1.d0  -  (beta_Z/4.5d0)**2 )
    w_fac_max = gam_Z * w_Z / (denZ_ion(1) * gam_2_min)
    
    ! In the following quadratic equation, we'll use w_fac_max/rm_avg in the
    !   coefficients.  So calculate rm_avg here
    rm_avg = 0.d0
    do j = 1, n_ions
      rm_avg = rm_avg  +  rm_ion(j) * den_rel_ion(j)
    enddo
    ! Handle possibility that electrons aren't included self-consistently
    if( .not. l_sc_elec ) then
      rm_avg = rm_avg  +  rm_elec * den_elec/denZ_ion(1)
    endif
    
    p2_max_A  = 16.d0/9.d0
    p2_max_B  = 8.d0/3.d0 - (w_fac_max/rm_avg)**2
    p2_max_C  = 1.d0 - (w_fac_max/rm_avg)**2
    p2_max_sq = (-p2_max_B + sqrt(p2_max_B**2 - 4.d0*p2_max_A*p2_max_C))  &!&
               / (2.d0 * p2_max_A)
    p2_max    = sqrt(p2_max_sq) * xmp*ccgs

    ! Initial guess for downstream proton momentum, as well as nearby
    !   location to use for finite difference approximation
    p2_guess   = p2_max / 1.001d0
    p2_guess_p = p2_max
    del_p2_guess = p2_guess_p - p2_guess
    
    
    do i = 1, maxitrs
      
      !------------------------------------------
      ! Calculate downstream momentum fluxes associated with both p2_guess
      !   and p2_guess_p
      !  1) Calculate enthalpy, which is summed over all particle species
      !  2) Calculate pressure similarly
      !------------------------------------------
      p2_o_mc = p2_guess / (xmp * ccgs)  ! Protons only here because of
                                         !   how the math in P_fac & w_fac
                                         !   works out
      
      P_fac = 0.d0
      w_fac = 0.d0
      do j = 1, n_ions
        ! Pressure with proton density factored out
        P_fac = P_fac  +  den_rel_ion(j) * third * rm_ion(j)              &!&
                         * p2_o_mc**2 / sqrt( 1.d0  +  p2_o_mc**2 )
        
        ! Enthalpy with proton density factored out
        w_fac = w_fac  +  den_rel_ion(j) * rm_ion(j)                      &!&
                         * ( sqrt( 1.d0 + p2_o_mc**2 )                    &!&
                            + third * p2_o_mc**2 / sqrt(1.d0 + p2_o_mc**2 ) )
      enddo
      
      ! Handle possibility that electrons aren't included self-consistently
      if( .not. l_sc_elec ) then
        P_fac = P_fac  +  (den_elec / denZ_ion(1)) * rm_elec              &!&
                        * third * p2_o_mc**2 / sqrt( 1.d0 + p2_o_mc**2 )
        w_fac = w_fac  +  (den_elec / denZ_ion(1)) * rm_elec              &!&
                        * ( sqrt( 1.d0 + p2_o_mc**2 )                     &!&
                           + third * p2_o_mc**2 / sqrt( 1.d0 + p2_o_mc**2 ) )
      endif

      gam_2     = gam_Z * w_Z / (denZ_ion(1) * w_fac)  ! Using only proton
                                                       !   density is correct
      gambeta_2 = sqrt( gam_2**2 - 1.d0 )
      
      F_p2_guess = gambeta_2 * UpS_num_flux * w_fac                       &!&
                 +  UpS_num_flux * P_fac / gambeta_2  -  UpS_mom_flux
      
      !------------------------------------------
      
      p2_o_mc = p2_guess_p / (xmp * ccgs)  ! Protons only here because of
                                           !   how the math in P_fac & w_fac
                                           !   works out
      
      P_fac = 0.d0
      w_fac = 0.d0
      do j = 1, n_ions
        ! Pressure with proton density factored out
        P_fac = P_fac  +  den_rel_ion(j) * third * rm_ion(j)              &!&
                         * p2_o_mc**2 / sqrt( 1.d0  +  p2_o_mc**2 )
        
        ! Enthalpy with proton density factored out
        w_fac = w_fac  +  den_rel_ion(j) * rm_ion(j)                      &!&
                         * ( sqrt( 1.d0 + p2_o_mc**2 )                    &!&
                            + third * p2_o_mc**2 / sqrt(1.d0 + p2_o_mc**2 ) )
      enddo
      
      ! Handle possibility that electrons aren't included self-consistently
      if( .not. l_sc_elec ) then
        P_fac = P_fac  +  (den_elec / denZ_ion(1)) * rm_elec              &!&
                        * third * p2_o_mc**2 / sqrt( 1.d0 + p2_o_mc**2 )
        w_fac = w_fac  +  (den_elec / denZ_ion(1)) * rm_elec              &!&
                        * ( sqrt( 1.d0 + p2_o_mc**2 )                     &!&
                           + third * p2_o_mc**2 / sqrt( 1.d0 + p2_o_mc**2 ) )
      endif
      
      gam_2     = gam_Z * w_Z / (denZ_ion(1) * w_fac)  ! Using only proton
                                                       !   density is correct
      gambeta_2 = sqrt( gam_2**2 - 1.d0 )
      
      F_p2_guess_p = gambeta_2 * UpS_num_flux * w_fac                     &!&
                   +  UpS_num_flux * P_fac / gambeta_2  -  UpS_mom_flux
      !------------------------------------------
      
      
      ! Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
      Fprime_p2_guess = (F_p2_guess_p - F_p2_guess) / del_p2_guess
      
      
      ! Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
      p2_guess_next = p2_guess  -  F_p2_guess/Fprime_p2_guess
      
      
      ! Relative change in this step
      err_curr = (p2_guess_next - p2_guess) / p2_guess
      
      ! If the relative change is small enough, we've found our solution and
      !   can exit the loop; otherwise return for another cycle
      if( abs(err_curr) .lt. target_err ) then
        p2_found = p2_guess_next
        exit
      endif
      
      
      ! Make sure new value for p2_guess is less than p2_max.  Use a weighted
      !   average for p2_guess to converge faster.
      if( p2_guess_next*1.001d0 .ge. p2_max ) then
        p2_guess_p   = 0.2d0 * (p2_guess  +  4.d0*p2_max)
        p2_guess     = p2_guess_p / 1.001d0
        del_p2_guess = p2_guess_p - p2_guess
      else
        p2_guess     = p2_guess_next
        p2_guess_p   = 1.001d0 * p2_guess
        del_p2_guess = p2_guess_p - p2_guess
      endif
      
    enddo
    
    
    ! Did we hit the maximum number of iterations without finding the
    !   flux-conserving solution?
    if( i .ge. maxitrs ) then
      write(*,"(A)") 'ERROR in calc_rRH: Newton method did not find solution'
      write(*,"(A)") 'Stopping program now.'
      stop
    endif
    
    
    ! Calculate the compression ratio beta_Z/beta_2 associated with p2_found
    p2_o_mc = p2_found / (xmp * ccgs)  ! Protons only here because of how the
                                       !   math in w_fac works out
    
    ! Pressure, internal energy, and enthalpy with proton density factored
    !   out
    P_fac = 0.d0
    e_fac = 0.d0
    w_fac = 0.d0
    do j = 1, n_ions
      P_fac = P_fac  +  den_rel_ion(j) * third * rm_ion(j)                &!&
                       * p2_o_mc**2 / sqrt( 1.d0  +  p2_o_mc**2 )
        
      e_fac = e_fac  +  den_rel_ion(j) * rm_ion(j)                        &!&
                       * ( sqrt( 1.d0 + p2_o_mc**2 ) - 1.d0 )
        
      w_fac = w_fac  +  den_rel_ion(j) * rm_ion(j)                        &!&
                       * ( sqrt( 1.d0 + p2_o_mc**2 )                      &!&
                          + third * p2_o_mc**2 / sqrt(1.d0 + p2_o_mc**2 ) )
    enddo
    
    ! Handle possibility that electrons aren't included self-consistently
    if( .not. l_sc_elec ) then
      P_fac = P_fac  +  (den_elec / denZ_ion(1)) * rm_elec                &!&
                      * third * p2_o_mc**2 / sqrt( 1.d0 + p2_o_mc**2 )
      e_fac = e_fac  +  (den_elec / denZ_ion(1)) * rm_elec                &!&
                      * ( sqrt( 1.d0 + p2_o_mc**2 ) - 1.d0 )
      w_fac = w_fac  +  (den_elec / denZ_ion(1)) * rm_elec                &!&
                      * ( sqrt( 1.d0 + p2_o_mc**2 )                       &!&
                         + third * p2_o_mc**2 / sqrt( 1.d0 + p2_o_mc**2 ) )
    endif
    
    ! Calculate adiabatic index downstream
    phi_fac = e_fac/P_fac
    gam_adiab_2_RH = (phi_fac + 1.d0)/phi_fac
    
    ! Finally, get downstream speed and compression ratio
    gam_2  = gam_Z * w_Z / (denZ_ion(1) * w_fac)  ! Using only proton
                                                  !   density is correct
    beta_2 = sqrt( 1.d0 - 1.d0/gam_2**2 )
    
    
    rRH = beta_Z/beta_2
    !------------------------------------------------------------------------
    ! rRH found using Newton's method
    
  !--------------------------------------------------------------------------
  !  Possibility 3: Oblique at any shock speed.  Not currently supported by
  !    code, but included here in case code is extended in future.
  !--------------------------------------------------------------------------
  else
    write(*,"(2A)") 'ERROR in calc_rRH: not implemented for oblique ',    &!&
        'shocks yet.'
    write(*,"(2A)") "If this ever changes, don't forget to update the",   &!&
        ' shock profile in subroutine "setup_profile".'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif

return
end subroutine calc_rRH


!****************************************************************************
!****************************************************************************
subroutine data_input(u_Z, gam_Z, beta_Z, n_ions, aa_ion, zz_ion, tZ_ion, &!&
     denZ_ion, l_sc_elec, tZ_elec, i_in_distr, en_inj, l_inj_wt, Emax_keV,&!&
     Emax_keV_per_aa, pmax_cgs, eta_mfp, bmag_Z, rg0, theta_BZ, l_oblique,&!&
     x_grid_start_rg, x_grid_stop_rg, feb_UpS, feb_DwS, l_use_prp,        &!&
     n_xspec, x_spec, n_itrs, xn_per_coarse, xn_per_fine, n_pts_inj,      &!&
     n_pts_pcut, n_pts_pcut_hi, en_pcut_hi, n_pcuts, pcuts_in, dont_shock,&!&
     dont_scatter, dont_DSA, do_smoothing, prof_wt_fac, do_prof_fac_damp, &!&
     smooth_mom_en_fac, smooth_press_flux_psd_fac, r_comp, rRH,           &!&
     do_old_prof, n_old_skip, n_old_profs, n_old_per_prof, age_max,       &!&
     do_retro, do_fast_push, x_fast_stop_rg, x_art_start_rg, x_art_scale, &!&
     p_elec_crit, gam_elec_crit, do_rad_losses, jet_rad_pc, jet_sph_frac, &!&
     jet_open_ang_deg, jet_dist_kpc, redshift, en_xfer_frac,              &!&
     bturb_comp_frac, bfield_amp, gam_adiab_2_RH, psd_bins_per_dec_mom,   &!&
     psd_bins_per_dec_tht, psd_lin_cos_bins, psd_log_tht_decs, u_2,       &!&
     beta_2, gam_2, bmag_2, theta_B2, theta_u2, use_custom_frg,           &!&
     emin_therm_fac, do_multi_dNdps, do_tcuts, n_tcuts, tcuts, inj_fracs, &!&
     use_custom_epsB)

use constants
use parameters
use randnum

implicit none

! Read in data.  Unlike in previous versions of the code, data MUST be piped
!  in from an input file using the syntax
!     program.exe < mc_in.dat
! Entering data with keyboard is no longer an option!
!
! NOTE: to add new variables, take these steps
!  1) Select keyword that does not conflict with an existing keyword (the
!    keywords are listed in alpha order to help with this)
!  2) Create an "isset_*****" flag, initialized to "false"
!  3) Add appropriate code in both the read-in and the default sections of
!    this subroutine
!  4) Add the "isset_*****" flag to the long if check at the end of the
!    subroutine


  ! Output arguments
integer, intent(out) :: n_ions, i_in_distr, n_xspec, n_itrs, n_pts_inj,   &!&
     n_pts_pcut, n_pts_pcut_hi, n_pcuts, n_old_skip, n_old_profs,         &!&
     n_old_per_prof, psd_bins_per_dec_mom, psd_bins_per_dec_tht,          &!&
     psd_lin_cos_bins, psd_log_tht_decs, n_tcuts
logical, intent(out) :: l_sc_elec, l_inj_wt, l_oblique, l_use_prp,        &!&
     dont_shock, dont_scatter, dont_DSA, do_smoothing, do_prof_fac_damp,  &!&
     do_old_prof, do_retro, do_fast_push, do_rad_losses, use_custom_frg,  &!&
     do_multi_dNdps, do_tcuts, use_custom_epsB
real(kind=8), intent(out) :: u_Z, gam_Z, beta_Z, tZ_elec, en_inj,         &!&
     Emax_keV, Emax_keV_per_aa, pmax_cgs, eta_mfp, bmag_Z, rg0, theta_BZ, &!&
     x_grid_start_rg, x_grid_stop_rg, feb_UpS, feb_DwS,                   &!&
     xn_per_coarse, xn_per_fine, en_pcut_hi, prof_wt_fac,                 &!&
     smooth_mom_en_fac, smooth_press_flux_psd_fac, r_comp, age_max,       &!&
     x_fast_stop_rg, x_art_start_rg, x_art_scale, p_elec_crit,            &!&
     gam_elec_crit, jet_rad_pc, jet_sph_frac, jet_open_ang_deg,           &!&
     jet_dist_kpc, redshift, en_xfer_frac, bturb_comp_frac, bfield_amp,   &!&
     gam_adiab_2_RH, rRH, u_2, beta_2, gam_2, bmag_2, theta_B2, theta_u2, &!&
     emin_therm_fac
real(kind=8), dimension(na_i) :: aa_ion, zz_ion, tZ_ion, denZ_ion, inj_fracs
real(kind=8), dimension(na_g) :: x_spec
real(kind=8), dimension(na_c) :: pcuts_in, tcuts

  ! Local variables
character(len=5) :: keyword
integer :: iseed, i, i_no_shock, i_no_scatter, i_no_DSA, i_smooth,        &!&
     i_damp_prof_wt_fac, i_read_old_prof, i_retro, i_fast_push,           &!&
     i_rad_losses, i_custom_frg, i_multi_dNdps, i_ions, i_custom_epsB
real(kind=8) :: den_elec, pcut_curr, tcut_curr, Emax_eff, pmax_eff,       &!&
     en_elec_crit_keV, en_elec_crit_rm
real(kind=8), dimension(2) :: febup, febdw, jetfr
real(kind=8), dimension(3) :: skspd, enmax

  ! Variables to flag when quantities are set
logical ::  &!&
     isset_AGEMX = .false., isset_ARTSM = .false., isset_BAMPF = .false., &!&
     isset_BTRBF = .false., isset_BMAGZ = .false., isset_EMNFP = .false., &!&
     isset_ENINJ = .false., isset_ENMAX = .false., isset_ENXFR = .false., &!&
     isset_FEBDW = .false., isset_FEBUP = .false., isset_FPUSH = .false., &!&
     isset_FPSTP = .false., isset_GYFAC = .false., isset_INDST = .false., &!&
     isset_INJWT = .false., isset_ISEED = .false., isset_JETDS = .false., &!&
     isset_JETFR = .false., isset_JETRD = .false., isset_NIONS = .false., &!&
     isset_NITRS = .false., isset_NODSA = .false., isset_NORAD = .false., &!&
     isset_NOSCT = .false., isset_NOSHK = .false., isset_NPTHI = .false., &!&
     isset_NPTLO = .false., isset_NSPEC = .false., isset_OLDIN = .false., &!&
     isset_OLDDT = .false., isset_PCUTS = .false., isset_RCOMP = .false., &!&
     isset_RETRO = .false., isset_SKSPD = .false., isset_SMIWT = .false., &!&
     isset_SMMOE = .false., isset_SMPFP = .false., isset_SMVWT = .false., &!&
     isset_SMSHK = .false., isset_TELEC = .false., isset_THTBZ = .false., &!&
     isset_XGDDW = .false., isset_XGDUP = .false., isset_XNPER = .false., &!&
     isset_PSDBD = .false., isset_PSDTB = .false., isset_NWFRG = .false., &!&
     isset_EMNFC = .false., isset_RDSHF = .false., isset_DNDPS = .false., &!&
     isset_TCUTS = .false., isset_INJFR = .false., isset_NWEPB = .false.


open(unit=9,status='unknown',file='mc_out.dat')

! Read the variable at the start of the input file
read(5,"(A5)",advance='no') keyword


! Loop over input file, reading in quantities as long as the current keyword
!  isn't "ENDIN", which signifies the end of input data.
!----------------------------------------------------------------------------
do while( keyword .ne. "ENDIN" )
  
  select case( keyword )
    !--------------------------------------------
    case( "AGEMX" )
        read(5,*) age_max
        isset_AGEMX = .true.
        !DOLATER?: force do_retro = true if age_max > 0?
        
        if( age_max .le. 0.d0 ) age_max = -1.d0
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, age_max, 'Maximum ',   &!&
            'allowed CR age (sec, explosion frame).  Ignored if <= 0'
    
    !--------------------------------------------
    case( "ARTSM" )
        read(5,*) x_art_start_rg, x_art_scale
        isset_ARTSM = .true.
        
        write(9,"(A5,2X,2ES10.3E2,5X,2A)") keyword, x_art_start_rg,       &!&
            x_art_scale, 'Artif. smoothing: start position[rg0], scale ', &!&
            'factor.  Ignored if 1st input >= 0'
    
    !--------------------------------------------
    case( "BAMPF" )
        read(5,*) bfield_amp
        isset_BAMPF = .true.
        
        if( bfield_amp .lt. 1.d0 ) then
          write(*,"(A)") 'ERROR in "BAMPF": must be >= 1.d0'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        if( isset_BTRBF .and. (bturb_comp_frac .eq. 0.d0) .and.           &!&
                (bfield_amp .gt. 1.d0) ) then
          write(*,"(2A)") 'ERROR in "BAMPF": bfield_amp > 1 has no ',     &!&
              'effect if "BTRBF" = 0'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, bfield_amp,            &!&
            'Amplification factor to use for B field.  1.0 means no amp'
    
    !--------------------------------------------
    case( "BTRBF" )
        read(5,*) bturb_comp_frac
        isset_BTRBF = .true.
        
        if( (bturb_comp_frac .lt. 0.d0) .or.                              &!&
            (bturb_comp_frac .gt. 1.d0) ) then
          write(*,"(A)") 'ERROR IN "BTRBF": must be in [0,1]'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        if( isset_BAMPF .and. (bfield_amp .gt. 1.d0) .and.                &!&
           (bturb_comp_frac .eq. 0.d0) ) then
          write(*,"(2A)") 'ERROR in "BTRBF": bfield_amp > 1 has no ',     &!&
              'effect if "BTRBF" = 0'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES9.2E2,5X,2A)") keyword, bturb_comp_frac,        &!&
            'Fraction of compressed B turb to use in scattering & ',      &!&
            'losses.  Must be in [0,1].  0 means ignore'
    
    !--------------------------------------------
    case( "BMAGZ" )
        read(5,*) bmag_Z
        isset_BMAGZ = .true.
        
        if( .not. isset_SKSPD ) then
          
          write(*,"(2A)") 'ERROR: to calculate rg0, u0 must be set in ',  &!&
              '"SKSPD" prior to "BMAGZ".  Modify ordering.'
          write(*,"(A)") 'Stopping program now.'
          stop
          
        else
          
          ! rg0 below is the gyroradius of a proton whose speed is u_Z that
          !   is gyrating in a field of strength bmag_Z
          ! Note that this formula is relativistically correct
          rg0 = (gam_Z * xmp * ccgs**2 * beta_Z) / (qcgs * bmag_Z)
          
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, bmag_Z, 'Far UpS ',    &!&
            'magnetic field in Gauss'
    
    !--------------------------------------------
    case( "DNDPS" )
        read(5,*) i_multi_dNdps
        isset_DNDPS = .true.
        
        if( i_multi_dNdps .eq. 66 ) then
          do_multi_dNdps = .true.
        else
          i_multi_dNdps = 0
          do_multi_dNdps = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_multi_dNdps, 'Enter "66"', &!&
            ' to write a separate dNdp for each iteration'
    
    !--------------------------------------------
    case( "EMNFC" )
        read(5,*) emin_therm_fac
        isset_EMNFC = .true.
        
        write(9,"(A5,2X,ES10.3E2,5X,3A)") keyword, emin_therm_fac, 'Min ',&!&
            'PSD limit (factor below thermal peak) if "INDST" = 1.  ',    &!&
            'Ignored otherwise.'
    
    !--------------------------------------------
    case( "EMNFP" )
        read(5,*) en_elec_crit_keV
        isset_EMNFP = .true.
        
        ! If needed, convert input energy[keV] to momentum and Lorentz factor
        if( en_elec_crit_keV .gt. 0 ) then
          en_elec_crit_rm = en_elec_crit_keV * xkevte / rm_elec
          
          ! Different forms for nonrel and rel momenta
          if( en_elec_crit_rm .lt. 1.d-2 ) then
            p_elec_crit   = (xme*ccgs) * sqrt( en_elec_crit_rm * 2.d0 )
            gam_elec_crit = 1.d0
          else
            p_elec_crit   = (xme*ccgs)                                    &!&
                           * sqrt( (en_elec_crit_rm + 1.d0)**2 - 1.d0 )
            gam_elec_crit = en_elec_crit_rm + 1.d0
          endif
          
        else
          en_elec_crit_keV = -1.d0
          p_elec_crit      = -1.d0
          gam_elec_crit    = -1.d0
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, en_elec_crit_keV,      &!&
            'Kinetic energy[keV] below which electrons have constant ',   &!&
            'mfp.  Ignored if <= 0'
    
    !--------------------------------------------
    case( "ENINJ" )
        read(5,*) en_inj
        isset_ENINJ = .true.
        
        write(9,"(A5,2X,ES10.3E2,5X,3A)") keyword, en_inj, 'Injection ',  &!&
            'energy[keV] and min PSD limit if "INDST" = 2.  Ignored ',    &!&
            'otherwise.'
    
    !--------------------------------------------
    case( "ENMAX" )
        read(5,*) enmax(1), enmax(2), enmax(3)
        isset_ENMAX = .true.
        
        if( enmax(1) .gt. 0.d0 ) then       ! All species have same max
          Emax_keV = enmax(1)               !   energy
          
          enmax(2) = 0.d0
          enmax(3) = 0.d0
          
          Emax_keV_per_aa = 0.d0
          pmax_cgs        = 0.d0
        
        else if( enmax(2) .gt. 0.d0 ) then  ! Max energy depends on aa
          Emax_keV_per_aa = enmax(2)
          
          enmax(1) = 0.d0
          enmax(3) = 0.d0
          
          Emax_keV = 0.d0
          pmax_cgs = 0.d0
        
        else if( enmax(3) .gt. 0.d0 ) then  ! All species have same max
          pmax_cgs = enmax(3) * xmp * ccgs  !   momentum
          
          enmax(1) = 0.d0
          enmax(2) = 0.d0
          
          Emax_keV        = 0.d0
          Emax_keV_per_aa = 0.d0
        
        else
          write(*,"(2A)") 'ERROR in "ENMAX": at least one choice must ',  &!&
             'be non-zero.'
          write(*,"(A)") 'Stopping program now.'
          stop
        
        endif
    
        write(9,"(A5,2X,3ES10.2E2,5X,3A)") keyword, enmax(1), enmax(2),   &!&
            enmax(3), 'Max particle energy as Emax[keV], Emax/nuc',       &!&
            '[keV/aa], pmax/(m_pc).  1st nonzero value used.  Also used', &!&
            ' in setting PSD'
    
    !--------------------------------------------
    case( "ENXFR" )
        read(5,*) en_xfer_frac
        isset_ENXFR = .true.
        
        if( (en_xfer_frac .lt. 0.d0) .or. (en_xfer_frac .gt. 1.d0) ) then
          write(*,"(A)") 'ERROR in "ENXFR": en_xfer_frac must be in [0,1]'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, en_xfer_frac,          &!&
            'Fraction of ion energy Xferred to electrons at 1st shock ',  &!&
            'crossing. Must be [0,1]; 0 to ignore'
    
    !--------------------------------------------
    case( "FEBDW" )
        read(5,*) febdw(1), febdw(2)
        isset_FEBDW = .true.
    
          ! Make sure that rg0 has been calculated to use in unit conversion
        if( .not. isset_BMAGZ ) then
          write(*,"(2A)") 'ERROR in "FEBDW": rg0 not yet set.  Place ',   &!&
              '"FEBDW" after "BMAGZ".'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        if( febdw(1) .gt. 0.d0 ) then
          feb_DwS = febdw(1) * rg0
          febdw(2) = 0.d0
        else if( febdw(2) .gt. 0.d0 ) then
          feb_DwS = febdw(2) * pc_to_cm
          febdw(1) = 0.d0
        else
          l_use_prp = .true.
        endif
        
        write(9,"(A5,2X,2ES11.3E2,5X,2A)") keyword, febdw(1), febdw(2),   &!&
            'DwS FEB, in [rg0] or [pc].  1st non-zero value used.  Must', &!&
            ' be positive if included; set <= 0 otherwise.'
            
    !--------------------------------------------
    case( "FEBUP" )
        read(5,*) febup(1), febup(2)
        isset_FEBUP = .true.
        
          ! Make sure that rg0 has been calculated to use in unit conversion
        if( .not. isset_BMAGZ ) then
          write(*,"(2A)") 'ERROR in "FEBUP": rg0 not yet set.  Place ',   &!&
              '"FEBUP" after "BMAGZ".'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        if( febup(1) .lt. 0.d0 ) then
          feb_UpS = febup(1) * rg0
          febup(2) = 0.d0
        else if( febup(2) .lt. 0.d0 ) then
          feb_UpS = febup(2) * pc_to_cm
          febup(1) = 0.d0
        else
          write(*,"(2A)") 'ERROR in "FEBUP": at least one choice must ',  &!&
              'be negative.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
          ! Check to ensure that grid start is at least as negative as
          !   location of UpS FEB
        if( .not. isset_XGDUP ) then
          write(*,"(2A)") 'ERROR in "FEBUP": x_grid_start not yet set.',  &!&
              '  Move "FEBUP" after "XGDUP".'
          write(*,"(A)") 'Stopping program now.'
          stop
        else if( (feb_UpS/rg0) .lt. x_grid_start_rg ) then
          write(*,"(2A)") 'ERROR in "FEBUP": UpS FEB must be within ',    &!&
              'x_grid_start'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,2ES11.3E2,5X,2A)") keyword, febup(1), febup(2),   &!&
            'UpS FEB, in [rg0] or [pc].  1st non-zero value used.  Must', &!&
            ' be negative.'
    
    !--------------------------------------------
    case( "FPUSH" )
        read(5,*) i_fast_push
        isset_FPUSH = .true.
        
        if( i_fast_push .eq. 66 ) then
          do_fast_push = .true.
        else
          i_fast_push  = 0
          do_fast_push = .false.
        endif
        
        if( isset_FPSTP .and. (x_fast_stop_rg .ge. 0.d0) .and.            &!&
            do_fast_push ) then
          write(*,"(2A)") 'ERROR IN "FPSTP": If using fast UpS ',         &!&
              'transport, must stop in UpS region (i.e. x < 0)'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_fast_push, 'Enter "66" ',  &!&
            'for fast Ups transport'
    
    !--------------------------------------------
    case( "FPSTP" )
        read(5,*) x_fast_stop_rg
        isset_FPSTP = .true.
        
        if( isset_FPUSH .and. (x_fast_stop_rg .ge. 0.d0) ) then
          write(*,"(2A)") 'ERROR IN "FPSTP": If using fast UpS ',         &!&
              'transport, must stop in UpS region (i.e. x < 0)'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, x_fast_stop_rg, 'UpS', &!&
            ' x position[rg0] where PROTON fast transport stops'
    
    !--------------------------------------------
    case( "GYFAC" )
        read(5,*) eta_mfp
        isset_GYFAC = .true.
        
        write(9,"(A5,2X,ES10.3E2,5X,3A)") keyword, eta_mfp, 'Gyrofactor,',&!&
            ' i.e. MFP = gyrofac*gyroradius.  Ignored if "NWFRG" = 66'
    
    !--------------------------------------------
    case( "INDST" )
        read(5,*) i_in_distr
        isset_INDST = .true.
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_in_distr, 'Input ',        &!&
            'distribution.  1 = thermal, 2 = delta function, 3 = other'
        
        ! Currently option 3 is not available.  Flag and stop if selected,
        !  and check for other incorrect entries
        if( (i_in_distr .ne. 1) .and. (i_in_distr .ne. 2) ) then
          write(*,"(A,I0,2A)") 'ERROR: i_in_distr = ',i_in_distr,' not ', &!&
              'permitted'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
    
    !--------------------------------------------
    ! Injection probabilities for each species.  If equal to 1, particles use
    !   typical thermal leakage model.  If equal to 0, injection is blocked
    !   completely (i.e., equivalent to "NODSA" = 66).  Low values imply
    !   test-particle shocks, while high values imply nonlinear shocks.
    case( "INJFR" )
        read(5,*) i_ions
        isset_INJFR = .true.
        
        ! Error checks to make sure number of ions (a) has been set, and
        !   (b) matches what we just read
        if( .not. isset_NIONS ) then
          write(*,"(2A)") 'ERROR in "INJFR": must come after "NIONS" has',&!&
              ' been set.  Modify order.'
          write(*,"(A)") 'Stopping program now'
          stop
        elseif( i_ions .ne. n_ions ) then
          write(*,"(3(A,I0))") 'ERROR in "INJFR": i_ions (', i_ions,      &!&
              ') does not match n_ions (', n_ions, ').  Change to match'
          write(*,"(A)") 'Stopping program now'
          stop
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_ions, '# of modified ',    &!&
            'injection rates, 1 per species'
        
        do i = 1, i_ions
          read(5,*) inj_fracs(i)
          
          write(9,"(5X,ES9.2E2,5X,2A,I0)") inj_fracs(i), 'Injection ',    &!&
              'prob. for particles of species ', i
          
          ! Error check
          if( (inj_fracs(i) .lt. 0.d0) .or. (inj_fracs(i) .gt. 1.d0) ) then
            write(*,"(2A,I0,A,ES12.3E2)") 'ERROR: injection chance for ', &!&
                'species ', i, ' must be in [0,1]: ', inj_fracs(i)
            write(*,"(A)") 'Stopping program now'
            stop
          endif
        enddo
    
    !--------------------------------------------
    case( "INJWT" )
        read(5,*) i
        isset_INJWT = .true.
        
        if( i .eq. 1 ) then
          l_inj_wt = .true.
        else if( i .eq. 2 ) then
          l_inj_wt = .false.
        else
          write(*,"(A,I0,A)") 'ERROR: "INJWT" = ',i,' not permitted value.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i, 'How to assign initial ', &!&
            'weights.  1 = equal by particle, 2 = equal by bin'
    
    !--------------------------------------------
    ! Set seed for random number generator
    case( "ISEED" )
        read(5,*) iseed
        isset_ISEED = .true.
        
        iseed_in = iseed  ! iseed_in found in module randnum
        call kisset( iseed, 2*iseed, 3*iseed, 4*iseed )
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, iseed, 'Seed for random # ', &!&
            'generator'
    
    !--------------------------------------------
    case( "JETDS" )
        read(5,*) jet_dist_kpc
        isset_JETDS = .true.
        
        ! Because cosmo_calc subroutine expects only one of distance/redshift
        !   to be non-zero, check for input validity here
        if( isset_RDSHF .and. (redshift .gt. 0.d0) .and.                  &!&
           (jet_dist_kpc .gt. 0.d0) ) then
          write(*,"(2A)") 'ERROR in "JETDS": At most one of "JETDS" and', &!&
              ' "RDSHF" may be non-zero.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, jet_dist_kpc,          &!&
            'Distance from jet to observer[kpc].  Used to calc ',         &!&
            'CMB parameters'
    
    !--------------------------------------------
    case( "JETFR" )
        read(5,*) jetfr(1), jetfr(2)
        isset_JETFR = .true.
        
        if( (jetfr(1) .gt. 0.d0) .and. (jetfr(1) .le. 1.d0) ) then
          jetfr(2) = 0.d0
          
          jet_sph_frac     = jetfr(1)
          jet_open_ang_deg = acos(1.d0 - 2.d0*jet_sph_frac)
          jet_open_ang_deg = jet_open_ang_deg * radtdg
          
        else if( (jetfr(2) .gt. 0.d0) .and. (jetfr(2) .le. 180.d0) ) then
          jetfr(1) = 0.d0
          
          jet_open_ang_deg = jetfr(2)
          jet_sph_frac     = 0.5d0 * ( 1.d0 - cos(jet_open_ang_deg*degtrd) )
        
        else
          write(*,"(A)") 'ERROR IN "JETFR": Unphysical values entered.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,2ES9.2E2,5X,2A)") keyword, jetfr(1), jetfr(2),    &!&
            'Frac of sphere producing CRs, or jet open ang[deg].  1st ',  &!&
            'non-zero value used'
    
    !--------------------------------------------
    case( "JETRD" )
        read(5,*) jet_rad_pc
        isset_JETRD = .true.
        
        write(9,"(A5,2X,ES10.3E2,5X,3A)") keyword, jet_rad_pc, 'Jet ',    &!&
            'shock radius[pc] used for normalizing dN/dp.'
    
    !--------------------------------------------
    ! Number of particle species, and their properties
    case( "NIONS" )
        read(5,*) n_ions
        isset_NIONS = .true.
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, n_ions, '# of different ion',&!&
            ' species.  1st MUST be protons.  Search "NIONS" for details.'
        
        !====================
        ! Notes on usage: 
        !  1) First species must always be protons.
        !  2) If any species is given aa < 1, this will be treated as
        !    electrons with the given atomic mass number.  This means their
        !    charge is set to 1 (overriding any other value provided.  As
        !    well, the electron temperature (read in using the keyword
        !    "TELEC") will automatically be set to 0.
        !  3) Use aa = -99 to get electrons with correct mass.
        !  4) If NO addtional electron/positron pairs are to be injected, the
        !    density for the aa < 1 species must be 0.  In this case the
        !    electron density will be the charge-neutralizing density,
        !    assuming full ionization.
        !  5) If a density greater than 0 is entered for the aa < 1 species,
        !    that will be the TOTAL number of leptons ADDED.  In other words,
        !    1/2 the given density comes from electrons, and the other half
        !    from positrons.  These leptons will be added to the density of
        !    electrons required for charge neutralization, as discussed in
        !    point (2) above.
        !  6) If no species is given aa < 1, the background electron fluid
        !    will be included using charge neutrality to set density, and
        !    the electron temperature (keyword "TELEC") to set pressure.
        !====================
        
        do i = 1, n_ions
          read(5,*) aa_ion(i), zz_ion(i), tZ_ion(i), denZ_ion(i)
          
          write(9,"(5X,ES10.3E2,ES9.1E2,2ES10.2E2,5X,2A,I0)") aa_ion(i),  &!&
              zz_ion(i), tZ_ion(i), denZ_ion(i), 'A(-99 for e-), Z, ',    &!&
              'T0[K], & den0[/cc] for species ', i
          
          if( aa_ion(i) .eq. -99.d0 ) then
            aa_ion(i) = xme/xmp
          endif
          
          if( aa_ion(i) .lt. 1.d0 ) then
            zz_ion(i) = 1.d0
          endif
          
          ! Error check for a common mistake when trying to include electrons
          if( aa_ion(i) .le. 0.d0 ) then
            write(*,"(A,I0,A,ES12.3E2)") 'ERROR: species ', i,            &!&
                ' given invalid mass number: ', aa_ion(i)
            write(*,"(A)") 'Stopping program now'
            stop
          endif
        enddo
        
        
        ! Determine electron density if electrons are explicitly included,
        !  i.e. a species with aa < 1 is entered.
        if( minval(aa_ion(1:n_ions)) .lt. 1.d0 ) then
          
          den_elec = 0.d0  ! Set electron density before loop
          
          ! Determine charge neutralizing density of electrons based on
          !  inputted ion species
          do i = 1, n_ions
            if( aa_ion(i) .ge. 1.d0 ) den_elec = den_elec                 &!&
                                                +  zz_ion(i) * denZ_ion(i)
          enddo
          
          ! Include additional pairs if directed.  denZ_ion for electron
          !  species should always be total number of leptons (both electrons
          !  and positrons) per cubic centimeter.
          do i = 1, n_ions
            if( aa_ion(i) .lt. 1.d0 ) then
              if( denZ_ion(i) .le. 0.d0 ) then
                denZ_ion(i) = den_elec
              else
                denZ_ion(i) = denZ_ion(i)  +  den_elec
              endif
            endif
          enddo
          
          ! Tell code whether electrons are included self-consistently, i.e.
          !  as an individual species.  Force tZ_elec to 0 if so.
          l_sc_elec = .true.
          tZ_elec   = 0.d0
          
        else
        
          ! Tell code whether electrons are included self-consistently, i.e.
          !  as an individual species
          l_sc_elec = .false.
        
        endif
    
    !--------------------------------------------
    case( "NITRS" )
        read(5,*) n_itrs
        isset_NITRS = .true.
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, n_itrs, 'Number of ',        &!&
            'iterations to perform'
    
    !--------------------------------------------
    case( "NODSA" )
        read(5,*) i_no_DSA
        isset_NODSA = .true.
        
        if( i_no_DSA .eq. 66 ) then
          dont_DSA = .true.
        else
          i_no_DSA = 0
          dont_DSA = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_no_DSA, 'Enter "66" to ',  &!&
            'turn off DSA.  For testing only.'
    
    !--------------------------------------------
    case( "NORAD" )
        read(5,*) i_rad_losses
        isset_NORAD = .true.
        
        if( i_rad_losses .ne. 66 ) then
          i_rad_losses  = 0
          do_rad_losses = .true.
        else
          do_rad_losses = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_rad_losses, 'Enter "66" ', &!&
            'for NO radiation losses'
    
    !--------------------------------------------
    case( "NOSCT" )
        read(5,*) i_no_scatter
        isset_NOSCT = .true.
        
        if( i_no_scatter .eq. 66 ) then
          dont_scatter = .true.
        else
          i_no_scatter = 0
          dont_scatter = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_no_scatter, 'Enter "66"',  &!&
            ' for scatter-FREE propagation.  For testing only.'
    
    !--------------------------------------------
    case( "NOSHK" )
        read(5,*) i_no_shock
        isset_NOSHK = .true.
        
        if( i_no_shock .eq. 66 ) then
          dont_shock = .true.
        else
          i_no_shock = 0
          dont_shock = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_no_shock, 'Enter "66" for',&!&
            ' no shock, i.e. force r_comp = 1.  For testing only.'
    
    !--------------------------------------------
    case( "NPTHI" )
        read(5,*) n_pts_pcut_hi, en_pcut_hi
        isset_NPTHI = .true.
        
        ! Make sure particle arrays are large enough for desired population
        if( n_pts_pcut_hi .gt. na_p ) then
          write(*,"(A)") 'ERROR in "NPTHI": Array size na_p too small.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,I7,ES9.2E2,5X,2A)") keyword, n_pts_pcut_hi,       &!&
            en_pcut_hi, 'Target number of particles at hi-E pcuts, and ', &!& 
            'cutoff kinetic energy [keV/aa] between the two'
    
    !--------------------------------------------
    case( "NPTLO" )
        read(5,*) n_pts_inj, n_pts_pcut
        isset_NPTLO = .true.
        
        ! Make sure particle arrays are large enough for desired population
        if( max(n_pts_inj,n_pts_pcut) .gt. na_p ) then
          write(*,"(A)") 'ERROR in "NPTLO": Array size na_p too small.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,2I7,5X,2A)") keyword, n_pts_inj, n_pts_pcut,      &!&
            'Number of particles to inject, and target number of ',       &!&
            'particles at low-E pcuts.'
    
    !--------------------------------------------
    case( "NSPEC" )
        read(5,*) n_xspec
        isset_NSPEC = .true.
        
          ! Make sure that rg0 has been calculated to use in unit conversion
        if( .not. isset_BMAGZ ) then
          write(*,"(2A)") 'ERROR in "NSPEC": rg0 not yet set.  Place ',   &!&
              '"NSPEC" after "BMAGZ".'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, n_xspec, 'Number of x-',     &!&
            'positions [rg0] where particle spectrum is calculated'
        
        if( n_xspec .gt. 0 ) then
          
          do i = 1, n_xspec
            read(5,*) x_spec(i)
            
            write(9,"(5X,ES10.3E2,4X,2A,I0)") x_spec(i), 'x-position',    &!&
                '[rg0] #',i
            
            x_spec(i) = x_spec(i) * rg0
          enddo
          
        endif
    
    !--------------------------------------------
    ! Allows for definition of new grid variable, epsB_grid, which sets
    !   magnetic field strength everywhere.  epsB_grid must be defined in the
    !   subroutine setup_profile, and is used to set btot_grid in 
    !   smooth_grid_par.
    ! Overrides values of "BRTBF" and "BAMPF".
    case( "NWEPB" )
        read(5,*) i_custom_epsB
        isset_NWEPB = .true.
        
        if( i_custom_epsB .eq. 66 ) then
          use_custom_epsB = .true.
        else
          i_custom_epsB = 0
          use_custom_epsB = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,3A)") keyword, i_custom_epsB, 'Enter "66" ',&!&
            'to define custom epsilon_B on grid; else = B(x) controlled ',&!&
            'by "BTRBF" and "BAMPF"'
    
    !--------------------------------------------
    case( "NWFRG" )
        read(5,*) i_custom_frg
        isset_NWFRG = .true.
        
        if( i_custom_frg .eq. 66 ) then
          use_custom_frg = .true.
        else
          i_custom_frg = 0
          use_custom_frg = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,3A)") keyword, i_custom_frg, 'Enter "66" ', &!&
            'to define custom f(r_g) in subr "scattering"; else = ',      &!&
            'eta_mfp*r_g'
    
    !--------------------------------------------
    case( "OLDDT" )
        read(5,*) n_old_skip, n_old_profs, n_old_per_prof
        isset_OLDDT = .true.
        
        write(9,"(A5,2X,I6,I3,I4,5X,2A)") keyword, n_old_skip, n_old_profs,&!&
            n_old_per_prof, 'If reading old profile, # lines to skip, # ', &!&
            'profiles to average, # lines per profile'
    
    !--------------------------------------------
    case( "OLDIN" )
        read(5,*) i_read_old_prof
        isset_OLDIN = .true.
        
        if( i_read_old_prof .eq. 66 ) then
          do_old_prof = .true.
        else
          i_read_old_prof = 0
          do_old_prof = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_read_old_prof, 'Enter ',    &!&
            '"66" to read in old profile from "mc_grid_old.dat"'
    
    !--------------------------------------------
    case( "PCUTS" )
          ! pcuts must come after max energy/momentum has been set
        if( .not. isset_ENMAX ) then
          write(*,"(2A)") 'ERROR in "PCUTS": must come after "ENMAX" has',&!&
              ' been set.  Modify order.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        
        read(5,*)     ! Advance to lines of pcut input
        isset_PCUTS = .true.
        
        write(9,"(A5,5X,2A)") keyword, 'List of momentum cutoffs ',       &!&
            '[aa*m_pc] to use during iterations.'
        
        
        n_pcuts = 0
        pcut_curr = 1.d0
        
          ! Read in data from file.  A value of -1 signals end of list
        do while ( pcut_curr .ne. -1.d0 )
          read(5,*) pcut_curr
          
          n_pcuts           = n_pcuts + 1
          pcuts_in(n_pcuts) = pcut_curr
          
          write(9,"(5X,ES10.3E2,4X,A,I0)")  pcut_curr,'pcut  ',n_pcuts
        enddo
          ! Subtract 1 b/c final entry in pcuts array is -1
        n_pcuts = n_pcuts - 1
        
        ! Check to make sure we haven't used more pcuts than allowed by na_c
        if( (n_pcuts+1) .gt. na_c ) then
          write(*,"(2A)") 'ERROR in "PCUTS": parameter na_c smaller than',&!&
              ' desired number of pcuts.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        
        ! Ensure that final pcut falls above limit on energy/momentum.
        !   Include buffer pcut just to be safe.
        !--------------------------------------------------------------------
        if( Emax_keV .gt. 0.d0 ) then  ! Limit was on total energy.  Assume
                                       !   Fe for strictest limit on momentum
                                       !   per nucleon.
          
          ! Convert from momentum[m_pc/aa] to energy[keV]
          Emax_eff = 56.d0 * pcuts_in(n_pcuts-1) * xmp*ccgs * ccgs * etkev
          
          if( Emax_keV .gt. Emax_eff ) then
            write(*,"(2A)") 'ERROR in "PCUTS": max energy exceeds',       &!&
                ' highest pcut.  Add more pcuts or lower Emax_keV.'
            write(*,"(2(A,ES12.4E2))") '  Emax_keV (assuming Fe) = ',     &!&
                Emax_keV,'; Emax_eff = ',Emax_eff
            write(*,"(A)") 'Stopping program now.'
            stop
          endif
        
        else if( Emax_keV_per_aa .gt. 0.d0 ) then   ! Limit was on energy per
                                                    !   nucleon
          
          ! Convert from momentum[m_pc/aa] to energy[keV/aa]
          Emax_eff = pcuts_in(n_pcuts-1) * xmp*ccgs * ccgs * etkev
          
          if( Emax_keV_per_aa .gt. Emax_eff ) then
            write(*,"(2A)") 'ERROR in "PCUTS": max energy per aa exceeds',&!&
                ' highest pcut.  Add more pcuts or lower Emax_keV_per_aa.'
            write(*,"(2(A,ES12.4E2))") '  Emax_keV_per_aa = ',            &!&
                Emax_keV_per_aa, '; Emax_eff/aa = ',Emax_eff
            write(*,"(A)") 'Stopping program now.'
            stop
          endif
        
        else if( pmax_cgs .gt. 0.d0 ) then   ! Limit was on total momentum.
                                             !   Assume Fe for strictest
                                             !   limit on mom/nuc.
          
          pmax_eff = 56.d0*xmp*ccgs * pcuts_in(n_pcuts-1)
          
          if( pmax_cgs .gt. pmax_eff ) then
            write(*,"(2A)") 'ERROR in "PCUTS": max momentum exceeds ',    &!&
                'highest pcut.  Add more pcuts or lower pmax.'
            write(*,"(2(A,ES11.4E2))") '  pmax[m_pc] = ',pmax_cgs,        &!&
                '; pmax_eff (for Fe) = ',pmax_eff
            write(*,"(A)") 'Stopping program now.'
            stop
          endif
        
        else   ! Something unexpected has happened
          
          write(*,"(2A)") 'Unexpected result when comparing pcut max ',   &!&
              'to en/mom max'
          write(*,"(A)") 'Stopping program now.'
          stop
        
        endif
        !--------------------------------------------------------------------
    
    !--------------------------------------------
    case( "PSDBD" )
        read(5,*) psd_bins_per_dec_mom, psd_bins_per_dec_tht
        isset_PSDBD = .true.
        
        if( (psd_bins_per_dec_mom .le. 0) .or.                            &!&
            (psd_bins_per_dec_tht .le. 0) ) then
          write(*,"(2A)") 'ERROR in "PSDBD": both values must be positive.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,2I4,5X,2A)") keyword, psd_bins_per_dec_mom,       &!&
            psd_bins_per_dec_tht, '# of PSD bins per decade in (1) ',     &!&
            'momentum, (2) theta'
    
    !--------------------------------------------
    case( "PSDTB" )
        read(5,*) psd_lin_cos_bins, psd_log_tht_decs
        isset_PSDTB = .true.
        
        if( (psd_lin_cos_bins .le. 0) .or. (psd_log_tht_decs .le. 0) ) then
          write(*,"(2A)") 'ERROR in "PSDTB": both values must be positive.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,2I4,5X,2A)") keyword, psd_lin_cos_bins,           &!&
            psd_log_tht_decs, '# of PSD ang bins: (1) # linear (cosine) ',&!&
            'bins, (2) # log (theta) decades'
    
    !--------------------------------------------
    case( "RCOMP" )
        read(5,*) r_comp
        isset_RCOMP = .true.
        
        if( (.not. isset_SKSPD) .or. (.not. isset_BMAGZ) .or.             &!&
            (.not. isset_NIONS) .or. (.not. isset_THTBZ) ) then
          write(*,"(2A)") 'ERROR in "RCOMP": "SKSPD", "BMAGZ", "NIONS",', &!&
              ' and "THTBZ" must be set already.  Modify ordering.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        call calc_rRH(beta_Z, gam_Z, n_ions, aa_ion, zz_ion, denZ_ion,    &!&
           tZ_ion, l_sc_elec, tZ_elec, rRH, l_oblique, gam_adiab_2_RH)
        
        write(9,"(A5,2X,ES11.4E2,5X,A,ES12.5E2)") keyword, r_comp,        &!&
            'Target compression ratio.  = -1 to use R-H value of ', rRH
        
        if( r_comp .eq. -1.d0 ) then
          r_comp = rRH
        endif
        
        ! Calculate the downstream conditions for a test-particle shock with
        !   the inputted compression ratio; commented out version contains
        !   arguments for eventual oblique extension of code.
        call calc_DwS(l_oblique, bmag_Z, r_comp, beta_Z, beta_2, gam_2,   &!&
           bmag_2, theta_B2, theta_u2)
!comm   call calc_DwS(l_oblique, theta_BZ, bmag_Z, r_comp, beta_Z, gam_Z, &!&
!comm      n_ions, aa_ion, zz_ion, denZ_ion, tZ_ion, l_sc_elec, tZ_elec,  &!&
!comm      beta_2, gam_2, bmag_2, theta_B2, theta_u2)
        u_2 = beta_2 * ccgs
    
    !--------------------------------------------
    case( "RDSHF" )
        read(5,*) redshift
        isset_RDSHF = .true.
        
        ! Because cosmo_calc subroutine expects only one of distance/redshift
        !   to be non-zero, check for input validity here
        if( isset_JETDS .and. (jet_dist_kpc .gt. 0.d0) .and.              &!&
           (redshift .gt. 0.d0) ) then
          write(*,"(2A)") 'ERROR in "RDSHF": At most one of "JETDS" and', &!&
              ' "RDSHF" may be non-zero.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, redshift, 'Redshift ', &!&
            'of source.  Used to calc CMB parameters'
    
    !--------------------------------------------
    case( "RETRO" )
        read(5,*) i_retro
        isset_RETRO = .true.
        !DOLATER?: force do_retro = true if age_max > 0?
        
        if( i_retro .eq. 66 ) then
          do_retro = .true.
        else
          i_retro  = 0
          do_retro = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_retro, 'Enter "66" to use',&!&
            ' retro time calc DwS.  Ignored if "AGEMX" < 0'
    
    !--------------------------------------------
    ! Upstream speed in shock frame
    case( "SKSPD" )
        read(5,*) skspd(1), skspd(2), skspd(3)
        isset_SKSPD = .true.
        
        if( (skspd(1) .gt. 0.d0) .and. (skspd(1) .lt. ccgs*1.d-5) ) then
          u_Z    = skspd(1) * 1.d5
          beta_Z = u_Z / ccgs
          gam_Z  = 1.d0 / sqrt( 1.d0 - beta_Z**2 )
          
          skspd(2) = 0.d0
          skspd(3) = 0.d0
        
        else if( skspd(2) .gt. 1.d0 ) then
          gam_Z  = skspd(2)
          beta_Z = sqrt( 1.d0 - 1.d0/gam_Z**2 )
          u_Z    = beta_Z * ccgs
          
          skspd(1) = 0.d0
          skspd(3) = 0.d0
        
        else if( (skspd(3) .gt. 0.d0) .and. (skspd(3) .lt. 1.d0) ) then
          beta_Z = skspd(3)
          u_Z    = beta_Z * ccgs
          gam_Z  = 1.d0 / sqrt( 1.d0 - beta_Z**2 )
          
          skspd(1) = 0.d0
          skspd(2) = 0.d0
        
        else
          write(*,"(2A)") 'ERROR in "SKSPD": at least one choice must ',  &!&
             'be non-zero and physically reasonable.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,3ES12.4E2,5X,2A)") keyword, skspd(1), skspd(2),   &!&
            skspd(3), 'Shock speed as: km/sec, Lorentz factor, or beta.', &!&
            '  1st nonzero value used'
    
    !--------------------------------------------
    case( "SMIWT" )
        read(5,*) prof_wt_fac
        isset_SMIWT = .true.
        
        write(9,"(A5,2X,ES9.2E2,5X,3A)") keyword, prof_wt_fac, 'Factor',  &!&
            ' to weight old profile with in shock smoothing. > 1 favors', &!&
            ' old profile over new.  = 1 to use mean of the two'
    
    !--------------------------------------------
    case( "SMMOE" )
        read(5,*) smooth_mom_en_fac
        isset_SMMOE = .true.
        
        if( (smooth_mom_en_fac .lt. 0.d0) .or.                            &!&
            (smooth_mom_en_fac .gt. 1.d0) ) then
          write(*,"(2A,ES9.2E2,A)") 'ERROR in "SMMOE": smooth_mom_en_fac',&!&
              ' = ',smooth_mom_en_fac,'; must be in [0,1]'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES9.2E2,5X,2A)") keyword, smooth_mom_en_fac,      &!&
            'For profile iteration, = 0 uses only mom. flux eq, = 1 only',&!&
            ' energy.  Must be in [0,1]'
    
    !--------------------------------------------
    case( "SMPFP" )
        read(5,*) smooth_press_flux_psd_fac
        isset_SMPFP = .true.
        
        if( (smooth_press_flux_psd_fac .lt. 0.d0) .or.                    &!&
            (smooth_press_flux_psd_fac .gt. 1.d0) ) then
          write(*,"(2A,ES9.2E2,A)") 'ERROR in "SMPFP": ',                 &!&
              'smooth_press_flux_psd_fac = ',smooth_press_flux_psd_fac,   &!&
              '; must be in [0,1]'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        !DOLATER: actually get pressure calculation working properly
        if( smooth_press_flux_psd_fac .gt. 0.d0 ) then
          write(*,"(2A)") 'ERROR in "SMPFP": code does not properly ',    &!&
              'calculate pressure from PSD.'
          write(*,"(2A)") 'Set to 0 or get this code working'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES9.2E2,5X,3A)") keyword,                         &!&
            smooth_press_flux_psd_fac, 'For profile iteration, = 0 uses', &!&
            ' only flux to find pressure, = 1 only PSD data.  Must be ',  &!&
            'in [0,1]'
    
    !--------------------------------------------
    case( "SMVWT" )
        read(5,*) i_damp_prof_wt_fac
        isset_SMVWT = .true.
        
        if( i_damp_prof_wt_fac .eq. 66 ) then
          do_prof_fac_damp = .true.
        else
          i_damp_prof_wt_fac = 0
          do_prof_fac_damp   = .false.
        endif
        
        write(9,"(A5,2X,I6,5X,2A)") keyword, i_damp_prof_wt_fac, 'Enter ',&!&
            '"66" to increase old profile weighting for later iterations'
    
    !--------------------------------------------
    case( "SMSHK" )
        read(5,*) i_smooth
        isset_SMSHK = .true.
        
        if( i_smooth .eq. 66 ) then
          do_smoothing = .false.
        else
          i_smooth     = 0
          do_smoothing = .true.
        endif
        
        write(9,"(A5,2X,I6,5X,3A)") keyword, i_smooth, 'Enter "66" to ',  &!&
            'keep velocity profile const. between iterations (i.e. NO ',  &!&
            'smoothing)'
    
    !--------------------------------------------
    case( "TCUTS" )
        read(5,*)     ! Advance to lines of pcut input
        isset_TCUTS = .true.
        
        ! This tracking must be performed with an acceleration time limit,
        !   so ensure that AGEMX has been set to a positive value
        if( .not. isset_AGEMX ) then
          write(*,"(2A)") 'ERROR: keyword "TCUTS" cannot come before',    &!&
              ' keyword "AGEMX".  Adjust input file.'
          write(*,"(A)") 'Stopping program now.'
          stop
        else if( age_max .le. 0.d0 ) then
          write(*,"(2A)") 'ERROR: tcut tracking must be used with an',    &!&
              'accel time limit.  Adjust keyword "AGEMX".'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        
        write(9,"(A5,5X,2A)") keyword, 'List of time cutoffs ',           &!&
            'to use for tracking particles.'
        
        do_tcuts = .true.
        
        n_tcuts = 0
        tcut_curr = 1.d0
        
          ! Read in data from file.  A value of -1 signals end of list
        do while ( tcut_curr .ne. -1.d0 )
          read(5,*) tcut_curr
          
          n_tcuts        = n_tcuts + 1
          tcuts(n_tcuts) = tcut_curr
          
          write(9,"(5X,ES10.3E2,4X,A,I0)")  tcut_curr,'tcut  ',n_tcuts
        enddo
          ! Subtract 1 b/c final entry in pcuts array is -1
        n_tcuts = n_tcuts - 1
        
        ! Check to make sure we haven't used more tcuts than allowed by na_c
        if( (n_tcuts+1) .gt. na_c ) then
          write(*,"(2A)") 'ERROR in "TCUTS": parameter na_c smaller than',&!&
              ' desired number of tcuts.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        ! Check to make sure final tcut is much larger than age_max so that
        !   we never have to worry about exceeding it
        if( tcuts(n_tcuts) .le. (10.d0*age_max) ) then
          write(*,"(2A)") 'ERROR in "TCUTS": final tcut must be much ',   &!&
              '(10x) larger than age_max.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
    
    !--------------------------------------------
    case( "TELEC" )
        read(5,*) tZ_elec
        isset_TELEC = .true.
        
        ! This keyword can't come before the particle species are set, and
        !  the code knows whether electrons are a separate species
        if( .not. isset_NIONS ) then
          write(*,"(2A)") 'ERROR: keyword "TELEC" cannot come before',    &!&
              ' keyword "NIONS".  Adjust input file.'
          write(*,"(A)") 'Stopping program now.'
          stop
        endif
        
        if( l_sc_elec ) then
          tZ_elec = 0.d0
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, tZ_elec, 'Far UpS ',   &!&
            'elec temp.  Is set to 0 if e- is a separate species.'
    
    !--------------------------------------------
    case( "THTBZ" )
        read(5,*) theta_BZ
        isset_THTBZ = .true.
        
        ! Code can't handle oblique shocks at the moment.  Stop running if
        !  theta_BZ > 0
        if( theta_BZ .gt. 0.d0 ) then
          
          l_oblique = .true.
          
          write(*,"(2A)") 'ERROR: code cannot currently handle oblique ', &!&
             'shocks.  Adjust "THTBZ".'
          write(*,"(A)") 'Stopping program now.'
          stop
          
        else if( theta_BZ .lt. 0.d0 ) then
          
          write(*,"(2A)") 'ERROR: unphysical value for "THTBZ". Must be ',&!&
              'at least 0.'
          write(*,"(A)") 'Stopping program now.'
          stop
          
        else
          
          l_oblique = .false.
          
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, theta_BZ, 'UpS angle', &!&
            '[deg] between B-field and shock normal'
    
    !--------------------------------------------
    case( "XGDDW" )
        read(5,*) x_grid_stop_rg
        isset_XGDDW = .true.
        
        if( x_grid_stop_rg .le. 0.d0 ) then
          write(9,"(A)") 'ERROR in "XGDDW": x_grid_stop must be positive.'
          write(9,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,2A)") keyword, x_grid_stop_rg, 'DwS ',&!&
            'limit of grid[rg0]'
    
    !--------------------------------------------
    case( "XGDUP" )
        read(5,*) x_grid_start_rg
        isset_XGDUP = .true.
        
        if( x_grid_start_rg .ge. 0.d0 ) then
          write(9,"(A)") 'ERROR in "XGDUP": x_grid_start must be negative.'
          write(9,"(A)") 'Stopping program now.'
          stop
        endif
        
        write(9,"(A5,2X,ES10.3E2,5X,A)") keyword, x_grid_start_rg,        &!&
            'Start of grid[rg0].  Can not be closer than UpS FEB.'
    
    !--------------------------------------------
    case( "XNPER" )
        read(5,*) xn_per_coarse, xn_per_fine
        isset_XNPER = .true.
        
        write(9,"(A5,2X,2ES9.2E2,3X,2A)") keyword, xn_per_coarse,         &!&
            xn_per_fine, 'Number of time steps per gyro period,',         &!&
            ' for "coarse" and "fine" scattering'
    
    !--------------------------------------------
    ! Error catcher
    case default
        write(*,"(3A)") 'ERROR: unrecognized keyword "',keyword,'" detected.'
        write(*,"(A)")  '(Keywords are case sensitive)'
        write(*,"(A)")  'Stopping program now.'
        stop
  end select


  !----------------------------------------------
  ! Read the keyword for the next quantity
  read(5,"(A5)",advance='no') keyword

enddo
!----------------------------------------------------------------------------
! Loop over input finished


close(5)


write(9,"(A)") 'ENDIN'
write(9,"(A)") ' **** End of input data ****'
write(9,*)
write(9,*)


! Set defaults or warn that key parameters have not been set
!----------------------------------------------------------------------------
if( .not. isset_AGEMX ) then
  age_max = -1.d0
  
  write(9,"(2A)") 'No value provided for "AGEMX".  Assuming -1'
endif

if( .not. isset_ARTSM ) then
  x_art_start_rg = 0.d0
  x_art_scale    = 0.d0
  
  write(9,"(2A)") 'No value provided for "ARTSM".  Assuming no ',         &!&
      'artificial smoothing'
endif

if( .not. isset_BAMPF ) then
  bfield_amp = 1.d0
  
  write(9,"(2A)") 'No value provided for "BAMPF".  Assuming 1'  
endif

if( .not. isset_BTRBF ) then
  bturb_comp_frac = 0.d0
  
  write(9,"(2A)") 'No value provided for "BTRBF".  Assuming 0'
endif

if( .not. isset_BMAGZ ) then
  write(*,"(A)") 'ERROR: "BMAGZ" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_DNDPS ) then
  do_multi_dNdps = .false.
  
  write(9,"(2A)") 'No value provided for "DNDPS".  Assuming 0'
endif

if( .not. isset_EMNFC ) then
  if( (isset_INDST .and. (i_in_distr .eq. 1)) .or. (.not. isset_INDST) ) then
    emin_therm_fac = 1.d-2
    
    write(9,"(2A)") 'No value provided for "EMNFC".  Min energy for PSD ',&!&
        'will be 0.01 * thermal peak'
  endif
endif

if( .not. isset_EMNFP ) then
  en_elec_crit_keV = -1.d0
  p_elec_crit      = -1.d0
  gam_elec_crit    = -1.d0
  
  write(9,"(2A)") 'No value provided for "EMNFP".  Assuming no minimum',  &!&
      ' electron MFP'
endif

if( .not. isset_ENINJ ) then
  if( isset_INDST .and. (i_in_distr .eq. 2) ) then
    write(*,"(2A)") 'ERROR: if "INDST" = 2, "ENINJ" must be specified ',  &!&
        'manually.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
endif

if( .not. isset_ENMAX ) then
  write(*,"(A)") 'ERROR: "ENMAX" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_ENXFR ) then
  en_xfer_frac = 0.d0
  
  write(9,"(2A)") 'No value provided for "ENXFR".  Assuming 0'
endif

if( .not. isset_FEBDW ) then
  feb_DwS = -1.d0
  
  write(9,"(2A)") 'No value provided for "FEBDW".  Assuming -1'
endif

if( .not. isset_FEBUP ) then
  if( .not. isset_XGDUP ) then
    write(*,"(2A)") 'ERROR: cannot set default "FEBUP" if "XGDUP" not set.'
    write(*,"(A)") 'Stopping program now.'
    stop
  else if( .not. isset_BMAGZ ) then
    write(*,"(2A)") 'ERROR: cannot set default "FEBUP" if "BMAGZ" not set.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  feb_UpS = x_grid_start_rg * rg0
  
  write(9,"(2A)") 'No value provided for "FEBUP".  Assuming x_grid_start'
endif

if( .not. isset_FPUSH ) then
  do_fast_push = .false.
  
  write(9,"(2A)") 'No value provided for "FPUSH".  Assuming no fast push'
endif

if( .not. isset_FPSTP ) then
  if( isset_FPUSH ) then
    write(*,"(2A)") 'ERROR: If using fast push, FPSTP must be specified', &!&
        ' manually.'
    write(*,"(A)") 'Stopping program now.'
    stop  
  endif
endif

if( .not. isset_GYFAC ) then
  eta_mfp = 1.d0
  
  write(9,"(2A)") 'No value provided for "GYFAC".  Assuming 1'
endif

if( .not. isset_INDST ) then
  i_in_distr = 1
  
  write(9,"(2A)") 'No value provided for "INDST".  Assuming 1'
endif

if( .not. isset_INJFR ) then
  write(9,"(A)") 'No values provided for "INJFR".  Assuming 1.0'
endif

if( .not. isset_INJWT ) then
  l_inj_wt = .true.
  
  write(9,"(2A)") 'No value provided for "INJWT".  Assuming 1'
endif

if( .not. isset_ISEED ) then
  write(9,"(2A)") 'No value provided for "ISEED".  Assuming default ',    &!&
      'values from module randnum'
endif

if( .not. isset_JETDS ) then
  if( .not. isset_RDSHF ) then
    jet_dist_kpc = 1.d0
    redshift     = 0.d0
    
    write(9,"(2A)") 'No value provided for "JETDS" or "RDSHF". Assuming ',&!&
        '1 for "JETDS" (0 for "RDSHF")'
  endif
endif

if( .not. isset_JETFR ) then
  jet_sph_frac = 1.d0
  jet_open_ang_deg = acos(1.d0 - 2.d0*jet_sph_frac) * radtdg
  
  write(9,"(2A)") 'No value provided for "JETFR". Assuming 1 for ',       &!&
      'jet_sph_frac (whole sphere)'
endif

if( .not. isset_JETRD ) then
  jet_rad_pc = 1.d0
  
  write(9,"(A)") 'No value provided for "JETRD". Assuming 1.0'
endif

if( .not. isset_NIONS ) then
  write(*,"(A)") 'ERROR: "NIONS" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_NITRS ) then
  write(*,"(A)") 'ERROR: "NITRS" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_NODSA ) then
  dont_DSA = .false.
  
  write(9,"(2A)") 'No value provided for "NODSA".  Assuming 0'
endif

if( .not. isset_NORAD ) then
  do_rad_losses = .true.
  
  write(9,"(2A)") 'No value provided for "NORAD".  Assuming 0'
endif

if( .not. isset_NOSCT ) then
  dont_scatter = .false.
  
  write(9,"(2A)") 'No value provided for "NOSCT".  Assuming 0'
endif

if( .not. isset_NOSHK ) then
  dont_shock = .false.
  
  write(9,"(2A)") 'No value provided for "NOSHK".  Assuming 0'
endif

if( .not. isset_NPTHI ) then
  write(*,"(A)") 'ERROR: "NPTHI" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_NPTLO ) then
  write(*,"(A)") 'ERROR: "NPTLO" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_NSPEC ) then
  n_xspec = 0
  
  write(9,"(2A)") 'No value provided for "NSPEC".  Assuming 0'
endif

if( .not. isset_NWEPB ) then
  use_custom_epsB = .false.
  
  write(9,"(2A)") 'No value provided for "NWEPB".  Assuming 0'
else
  if( isset_BTRBF .or. isset_BAMPF ) then
    write(9,"(2A)") 'NOTE: "NWEPB" and 1+ of "BTRBF"/"BAMPF" have been ', &!&
        'set. "BTRBF"/"BAMPF" will be ignored'
  endif
endif

if( .not. isset_NWFRG ) then
  use_custom_frg = .false.
  
  write(9,"(2A)") 'No value provided for "NWFRG".  Assuming 0'
endif

if( .not. isset_OLDDT ) then
  if( isset_OLDIN .and. do_old_prof) then
    write(*,"(A)") 'ERROR: If reading in old profile, "OLDDT" must be ',  &!&
        'specified manually.'
    write(*,"(A)") 'Stopping program now.'
    stop
  else
    write(9,"(2A)") 'No values provided for "OLDDT", but none needed ',   &!&
        '(no old profiles read in)'
  endif
endif

if( .not. isset_OLDIN ) then
  do_old_prof = .false.
  
  write(9,"(2A)") 'No value provided for "OLDIN".  Assuming new profile'
endif

if( .not. isset_PCUTS ) then
  n_pcuts = 67  ! This must be less than the value of parameter na_c
  
  pcuts_in(1:67) = (/      1.000E-02, 2.000E-02, 5.000E-02, 1.000E-01,    &!&
     1.500E-01, 2.000E-01, 2.500E-01, 3.000E-01, 3.200E-01, 3.500E-01,    &!&
     3.600E-01, 3.800E-01, 4.000E-01, 4.200E-01, 4.400E-01, 4.500E-01,    &!&
     4.600E-01, 4.700E-01, 4.800E-01, 4.900E-01, 5.500E-01, 6.000E-01,    &!&
     6.200E-01, 6.400E-01, 6.700E-01, 7.000E-01, 7.500E-01, 8.000E-01,    &!&
     8.500E-01, 9.000E-01, 9.500E-01, 1.050E+00, 1.300E+00, 1.600E+00,    &!&
     2.000E+00, 2.500E+00, 3.000E+00, 3.800E+00, 4.500E+00, 5.500E+00,    &!&
     7.000E+00, 8.000E+00, 9.000E+00, 1.000E+01, 1.500E+01, 2.000E+01,    &!&
     3.000E+01, 5.000E+01, 8.000E+01, 1.200E+02, 2.000E+02, 3.000E+02,    &!&
     5.000E+02, 1.000E+03, 2.000E+03, 5.000E+03, 1.000E+04, 3.000E+04,    &!&
     1.000E+05, 1.000E+06, 1.000E+07, 1.000E+08, 1.000E+09, 1.000E+10,    &!&
     1.000E+11, 1.000E+12, 1.000E+13 /)
  
  if( .not. isset_ENMAX ) then
    write(*,"(2A)") 'ERROR: "ENMAX" must be specified manually.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  ! Ensure that final pcut falls above limit on energy/momentum.  Include
  !   buffer pcut just to be safe.
  !--------------------------------------------------------------------------
  if( Emax_keV .gt. 0.d0 ) then  ! Limit was on total energy.  Assume Fe for
                                 !   strictest limit on momentum per nucleon
    
    ! Convert from momentum[m_pc/aa] to energy[keV]
    Emax_eff = 56.d0 * pcuts_in(n_pcuts-1) * xmp*ccgs * ccgs * etkev
    
    if( Emax_keV .gt. Emax_eff ) then
      write(*,"(2A)") 'ERROR in "PCUTS": max energy exceeds highest ',    &!&
          'pcut.  Add more pcuts or lower Emax_keV.'
      write(*,"(2(A,ES12.4E2))") '  Emax_keV (assuming Fe) = ', Emax_keV, &!&
          '; Emax_eff = ', Emax_eff
      write(*,"(A)") 'Stopping program now.'
      stop
    endif
  
  else if( Emax_keV_per_aa .gt. 0.d0 ) then   ! Limit was on energy per
                                              !   nucleon
    
    ! Convert from momentum[m_pc/aa] to energy[keV/aa]
    Emax_eff = pcuts_in(n_pcuts-1) * xmp*ccgs * ccgs * etkev
    
    if( Emax_keV_per_aa .gt. Emax_eff ) then
      write(*,"(2A)") 'ERROR in "PCUTS": max energy per aa exceeds ',     &!&
          ' highest pcut.  Add more pcuts or lower Emax_keV_per_aa.'
      write(*,"(2(A,ES12.4E2))") '  Emax_keV_per_aa = ', Emax_keV_per_aa, &!&
          '; Emax_eff/aa = ', Emax_eff
      write(*,"(A)") 'Stopping program now.'
      stop
    endif
  
  else if( pmax_cgs .gt. 0.d0 ) then   ! Limit was on total momentum.  Assume
                                       !   Fe for strictest limit on momentum
                                       !   per nucleon.
    
    pmax_eff = 56.d0*xmp*ccgs * pcuts_in(n_pcuts-1)
    
    if( pmax_cgs .gt. pmax_eff ) then
      write(*,"(2A)") 'ERROR in "PCUTS": max momentum exceeds highest ',  &!&
          'pcut.  Add more pcuts or lower pmax.'
      write(*,"(2(A,ES11.4E2))") '  pmax[m_pc] = ',pmax_cgs,              &!&
          '; pmax_eff (for Fe) = ',pmax_eff
      write(*,"(A)") 'Stopping program now.'
      stop
    endif
  
  else   ! Something unexpected has happened
    
    write(*,"(2A)") 'Unexpected result when comparing pcut max to ',      &!&
        'en/mom max'
    write(*,"(A)") 'Stopping program now.'
    stop
  
  endif
  !--------------------------------------------------------------------------
  
endif

if( .not. isset_PSDBD ) then
  psd_bins_per_dec_mom = 10
  psd_bins_per_dec_tht = 10
  
  write(9,"(2A)") 'No value provided for "PSDBD".  Assuming 10 for both ',&!&
      'quantities'
endif

if( .not. isset_PSDTB ) then
  psd_lin_cos_bins = 119
  psd_log_tht_decs = 4
  
  write(9,"(2A)") 'No value provided for "PSDTB".  Assuming 119 linear ',&!&
      'cosine bins, 4 log theta decades'
endif

if( .not. isset_RCOMP ) then
  write(*,"(A)") 'ERROR: "RCOMP" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_RDSHF ) then
  if( .not. isset_JETDS ) then
    jet_dist_kpc = 1.d0
    redshift     = 0.d0
    
    write(9,"(2A)") 'No value provided for "JETDS" or "RDSHF". Assuming ',&!&
        '0 for "RDSHF" (1 for "JETDS")'
  else
    write(9,"(2A)") 'No value provided for "RDSHF".  Value of "JETDS" ',  &!&
        'will be used to find redshift'
  endif
endif

if( .not. isset_RETRO ) then
  if( isset_AGEMX .and. (age_max .gt. 0.d0) ) then
    do_retro = .true.
    
    write(9,"(2A)") 'No value provided for "RETRO".  Assuming "66" b/c',  &!&
        '"AGEMX" set to > 0'
  else
    do_retro = .false.
    
    write(9,"(2A)") 'No value provided for "RETRO".  Assuming 0'
  endif
endif

if( .not. isset_SKSPD ) then
  write(*,"(A)") 'ERROR: "SKSPD" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_SMIWT ) then
  prof_wt_fac = 1.d0
  
  write(9,"(2A)") 'No value provided for "SMIWT".  Assuming 1.0'
endif

if( .not. isset_SMMOE ) then
  smooth_mom_en_fac = 0.d0
  
  write(9,"(2A)") 'No value provided for "SMMOE".  Assuming 0'
endif

if( .not. isset_SMPFP ) then
  smooth_press_flux_psd_fac = 0.d0
  
  write(9,"(2A)") 'No value provided for "SMPFP".  Assuming 0'
endif

if( .not. isset_SMVWT ) then
  do_prof_fac_damp = .false.
  
  write(9,"(2A)") 'No value provided for "SMVWT".  Assuming 0'
endif

if( .not. isset_SMSHK ) then
  write(*,"(A)") 'ERROR: "SMSHK" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_TCUTS ) then
  do_tcuts = .false.
  
  write(9,"(A)") 'No value provided for "TCUTS".  Tracking will not occur.'
endif

if( .not. isset_TELEC ) then
  if( .not. isset_NIONS ) then
    write(*,"(A)") 'ERROR: Cannot set default "TELEC" unless "NIONS" ',   &!&
        'specified manually.'
    write(*,"(A)") 'Stopping program now.'
    stop
  else
    tZ_elec = tZ_ion(1)
    
    write(9,"(2A)") 'No value provided for "TELEC".  Assuming = tZ_ion(1)'
  endif
endif

if( .not. isset_THTBZ ) then
  write(*,"(A)") 'ERROR: "THTBZ" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_XGDDW ) then
  write(*,"(A)") 'ERROR: "XGDDW" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_XGDUP ) then
  write(*,"(A)") 'ERROR: "XGDUP" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif

if( .not. isset_XNPER ) then
  write(*,"(A)") 'ERROR: "XNPER" must be specified manually.'
  write(*,"(A)") 'Stopping program now.'
  stop
endif
!----------------------------------------------------------------------------
! Input check complete


if( isset_AGEMX .and. isset_ARTSM .and. isset_BAMPF .and.                 &!&
    isset_BTRBF .and. isset_BMAGZ .and. isset_EMNFP .and.                 &!&
    isset_ENINJ .and. isset_ENMAX .and. isset_ENXFR .and.                 &!&
    isset_FEBDW .and. isset_FEBUP .and. isset_FPUSH .and.                 &!&
    isset_FPSTP .and. isset_GYFAC .and. isset_INDST .and.                 &!&
    isset_INJWT .and. isset_ISEED .and. isset_JETDS .and.                 &!&
    isset_JETFR .and. isset_JETRD .and. isset_NIONS .and.                 &!&
    isset_NITRS .and. isset_NODSA .and. isset_NORAD .and.                 &!&
    isset_NOSCT .and. isset_NOSHK .and. isset_NPTHI .and.                 &!&
    isset_NPTLO .and. isset_NSPEC .and. isset_OLDIN .and.                 &!&
    isset_OLDDT .and. isset_PCUTS .and. isset_RCOMP .and.                 &!&
    isset_RETRO .and. isset_SKSPD .and. isset_SMIWT .and.                 &!&
    isset_SMMOE .and. isset_SMPFP .and. isset_SMVWT .and.                 &!&
    isset_SMSHK .and. isset_TELEC .and. isset_THTBZ .and.                 &!&
    isset_XGDDW .and. isset_XGDUP .and. isset_XNPER .and.                 &!&
    isset_PSDBD .and. isset_PSDTB .and. isset_NWFRG .and.                 &!&
    isset_EMNFC .and. isset_RDSHF .and. isset_DNDPS .and.                 &!&
    isset_TCUTS .and. isset_INJFR .and. isset_NWEPB ) then
  write(9,"(A)") 'All keywords specified manually.  No defaults set.'
endif

write(9,*)
write(9,*)
  
return
end subroutine data_input


!****************************************************************************
!****************************************************************************
subroutine get_psd_bins(pxsk, ptsk, do_th, i_pt, j_th)

! Given parallel component and total of a particle's momentum in the shock
!  frame, determine which bin of psd particle will fall into.
!
! Binning is done by the value of ptsk in code units.  Angular binning is done
!  with cos(theta) for large angles or theta for small angles.
!
! In momentum all binning is logarithmic.  In angle the value of theta_fine
!  marks the division between linear bin spacing (above theta_fine) and
!  logarithmic spacing (from theta_fine down to theta_min). In both
!  logarithmic regions the parameter bins_per_decade_*** determines the
!  fineness of the bins.
!
! WARNING: binning in angle is actually done with the *NEGATIVE* of the 
!  particle's cosine.  This lets the most finely spaced bins correspond to
!  UpS-pointing particles rather than DwS-pointing ones.  This also has the
!  effect that angles are essentially measured from the -x axis rather than
!  the +x axis.
!
! Inputs:
!  1) pbsk: component of momentum parallel to shock normal (to B-field?),
!   in code units
!  2) ptsk: total particle momentum, in code units
!  3) do_th: logical telling subroutine whether to do theta calculations
! Outputs
!  1) i_pt: bin in momentum into which particle falls
!  2) j_th: bin in angle into which particle falls

use constants
use controls, only: psd_bins_per_dec_mom, psd_bins_per_dec_tht
use psd_vars, only: psd_mom_min, num_psd_mom_bins, num_psd_tht_bins,      &!&
      psd_cos_fine, del_cos, psd_tht_min

implicit none

  ! Input variables
real(kind=8), intent(in) :: pxsk, ptsk
logical :: do_th
  ! Output variables
integer, intent(out) :: i_pt, j_th

  ! Local variables
real(kind=8) :: p_cos, theta


!--! Bin in total momentum (i_pt)
  if ( ptsk .lt. psd_mom_min ) then
    
    ! Momentum close enough to 0
    i_pt = 0
  
  else 
    
    ! Particle falls into logarithmic spacing region
    i_pt = int( log10(ptsk/psd_mom_min) * psd_bins_per_dec_mom )  +  1
  
  endif
  
  ! Sanity check
  if( i_pt .gt. num_psd_mom_bins ) then
    write (*,*) "WARNING: particle momentum exceeded PSD's bounds!"
    write (*,"(2i4,2ES15.6E2)") i_pt, num_psd_mom_bins, ptsk/(xmp*ccgs),  &!&
       psd_mom_min * 10.d0**num_psd_mom_bins / (xmp*ccgs)
    i_pt = num_psd_mom_bins
  endif


!--! Bin in angle (j_th); note that we negate the pitch angle to provide the
   !  finest resolution (i.e. the logarithimcally-spaced angle bins rather
   !  than the linearly-spaced cosine bins) for particles that are directed
   !  upstream
  if( do_th ) then
    
    p_cos = -pxsk / ptsk
    
    if ( p_cos .lt. psd_cos_fine ) then
      
      ! Pitch angle falls within linear spacing
      j_th = num_psd_tht_bins - int( (p_cos + 1.d0) / del_cos )
      
    else
      
      ! Particle falls into logarithmic spacing region
      theta = acos( p_cos )
      
      if( theta .lt. psd_tht_min ) then
        
        ! theta is close enough to zero that it might as well be
        j_th = 0
        
      else
        
        j_th = int( log10(theta/psd_tht_min) * psd_bins_per_dec_tht )  +  1
        
      endif
      
    endif
    
    ! Check to avert floating point error
    if( j_th .gt. num_psd_tht_bins ) j_th = num_psd_tht_bins
    
  endif

return
end subroutine get_psd_bins


!****************************************************************************
!****************************************************************************
subroutine get_time(time)

! Use built-in function date_and_time to return the wall-clock time.
! Built-in function cpu_time isn't used because of potential issues if the
!  code is run on more than one CPU

! Output arguments
!  1) time: elapsed time in seconds since start of the month.  Note that
!   this will lead to weird results if the month rolls over in the middle of
!   a run

implicit none

  ! Output arguments
real(kind=8), intent(out) :: time
  ! Local variables
integer, dimension(8) :: values

  ! Don't bother with optional character arguments
  call date_and_time(VALUES=values)
  
  time =  values(3)*24.0*3600.0          &!&
        + values(5)*3600.0               &!&
        + values(6)*60.0                 &!&
        + values(7)                      &!&
        + values(8)*0.001

return
end subroutine get_time


!****************************************************************************
!****************************************************************************
subroutine init_pop(do_fast_push, i_in_distr, i_ion, aa, n_pts_use,       &!&
     xwt_in, ptpf_in, pbpf_in, x_PT_cm_in, i_grid_in, pxx_flux, pxz_flux, &!&
     en_flux)

! Initializes the particle populations that will propagate through the shock
!   structure.
! Handles fast push and associated flux-tracking & changes to the population
!
! Input arguments:
!
! Output arguments:

use constants
use parameters, only: na_p, na_g, beta_rel_fl, en_rel_pt
use controls, only: tZ_ion, en_inj, l_inj_wt, n_pts_inj, denZ_ion,        &!&
     x_grid_start, rg0, eta_mfp, x_fast_stop_rg, beta_Z, gam_Z, u_Z,      &!&
     n_ions, aa_ion, zz_ion, l_sc_elec, tZ_elec, l_oblique
use grid_vars, only: n_grid, x_grid_rg, uxsk_grid, gam_sf_grid
use iteration_vars, only: ptot_inj, wt_inj, n_pts_MB
use randnum

implicit none

  ! Input arguments
integer, intent(in) :: i_in_distr, i_ion
logical, intent(in) :: do_fast_push
real(kind=8), intent(in) :: aa
  ! Output arguments
integer, intent(out) :: n_pts_use
integer, dimension(na_p), intent(out) :: i_grid_in
real(kind=8), dimension(na_p), intent(out) :: xwt_in, ptpf_in, pbpf_in,   &!&
     x_PT_cm_in
real(kind=8), dimension(na_g), intent(out) :: pxx_flux, pxz_flux, en_flux

  ! Local variables
integer :: i_prt, i, i_stop
logical :: l_nonrel
real(kind=8) :: T_or_E, rand, den_ratio, gam_sph, press_ratio, temp_ratio,&!&
     press_Z, rho_Z, den_elec, rho_curr, press_curr, beta_curr,           &!&
     gam_beta_curr, flux_px, flux_pz, flux_en, gam_ptpf, vtpf, vmin_sf,   &!&
     vmax_sf, vxsf, vxpf


!--! If not using fast push, this procedure is quite quick.  Fill the arrays
   !   and return to the main loop
  if( .not. do_fast_push ) then
    
    if( i_in_distr .eq. 1 ) then
      T_or_E = tZ_ion(i_ion)
    elseif( i_in_distr .eq. 2 ) then
      T_or_E = en_inj
    else
      write(*,"(2A)") 'ERROR in init_pop: not set to handle i_in_distr ', &!&
          ' .gt. 2'
      write(*,"(A)") 'Stopping program now.'
      stop
    endif
    
    call set_in_dist(l_inj_wt, n_pts_inj, i_in_distr, T_or_E, aa,         &!&
       denZ_ion(i_ion), ptot_inj(:,i_ion), wt_inj(:,i_ion), n_pts_MB(i_ion))
    n_pts_use = n_pts_MB(i_ion)
    
    do i_prt = 1, n_pts_use
      xwt_in(i_prt)  = wt_inj(i_prt, i_ion)
      ptpf_in(i_prt) = ptot_inj(i_prt, i_ion)
      
      call kiss(rand)
      pbpf_in(i_prt) = ptpf_in(i_prt) * 2.d0*(rand - 0.5d0)
      
      x_PT_cm_in(i_prt) = x_grid_start - 10.d0*rg0*eta_mfp
      i_grid_in(i_prt)  = 0
    enddo
    
    return
  endif
  

!--!
!--! Everything from here on assumes fast push
!--!
  
  
!--! Fast push won't work with any distribution besides thermal (yet?), so
   !   stop if this occurs
  if( i_in_distr .gt. 1 ) then
    write(*,"(2A)") 'ERROR in init_pop: fast push will only work with ',  &!&
        'thermal input distr.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  
!--! Compute density, pressure, and temperature at termination of fast push;
   !   check to ensure that #assumecold still holds
   ! Start counting at zero so that same loop applies even if fast push
   !   *isn't* enabled
  do i = 0, n_grid
    if( x_grid_rg(i+1) .gt. x_fast_stop_rg ) then
      i_stop = i
      exit
    endif
  enddo
  
  if( beta_Z .lt. beta_rel_fl ) then
    l_nonrel = .true.
    den_ratio = u_Z / uxsk_grid(i_stop)
  else
    l_nonrel = .false.
    den_ratio = gam_Z*u_Z / (gam_sf_grid(i_stop) * uxsk_grid(i_stop))
  endif
  
  ! Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
  !   to calculate the far UpS internal energy
  ! #assumecold
  gam_sph = 5.d0 / 3.d0
  
  press_ratio = den_ratio**gam_sph
  temp_ratio  = press_ratio / den_ratio
  
  if( (xkb * tZ_ion(i_ion) * temp_ratio) .gt.                             &!&
      (4.d0 * aa*xmp*ccgs**2 * en_rel_pt) ) then
    write(*,"(2A)") 'ERROR in "init_pop": fast push cannot work b/c ',    &!&
        ' highest energy thermal'
    write(*,"(A)") '    particles become mildly relativistic'
    write(*,"(A)") 'Move fast push location UpS or disable entirely.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  
!--! Only run through the flux updates for the first particle species (i.e.
   !   protons, not that it matters here); skip thereafter
  if( i_ion .eq. 1 ) then
    
    ! Update the flux arrays as if the particles had actually crossed them
    !------------------------------------------------------------------------
    ! Obtain UpS pressure and density
    press_Z  = 0.d0
    rho_Z    = 0.d0
    den_elec = 0.d0
    do i = 1, n_ions
      press_Z  =  press_Z  +  denZ_ion(i) * xkb*tZ_ion(i)
      rho_Z    =  rho_Z    +  denZ_ion(i) * xmp*aa_ion(i)
      
      if( aa_ion(i) .ge. 1.d0 ) den_elec  =  den_elec + denZ_ion(i)*zz_ion(i)
    enddo
    
    ! If electrons were not a separate species, add them in here
    if( .not. l_sc_elec ) then
      press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
      rho_Z    =  rho_Z    +  den_elec * xme
    endif
    
    
    ! Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
    !   to calculate the far UpS pressure and internal energy
    ! #assumecold
    gam_sph = 5.d0 / 3.d0
    
    
    ! Calculate fluxes and update arrays; note that if fast push isn't
    !   enabled i_stop = 0, and this loop never executes
    do i = 1, i_stop
      
      den_ratio = (gam_Z * u_Z) / (gam_sf_grid(i) * uxsk_grid(i))
      rho_curr  = rho_Z * den_ratio
      
      press_ratio = den_ratio**gam_sph     ! Note assumption that gam_sph
                                           !   doesn't change from zone to
      press_curr  = press_Z * press_ratio  !   zone: #assumecold
      
      beta_curr     = uxsk_grid(i) / ccgs
      gam_beta_curr = gam_sf_grid(i) * uxsk_grid(i) / ccgs
      
      ! Determine fluxes while handling different possible orientations and
      !   shock speeds.
      ! For non-rel fluxes, expand out to beta^2 to allow for better matching
      !   with relativistic versions.
      ! WARNING: these fluxes do not include contributions from a strong
      !   magnetic field.  This is incorporated during the smoothing process.
      !----------------------------------------------------------------------
      if( (.not. l_oblique) .and. l_nonrel ) then
        
        flux_px = rho_curr * uxsk_grid(i)**2 * ( 1.d0 + beta_curr**2 )    &!&
                 +  press_curr * ( 1.d0  +  gam_sph/(gam_sph-1.d0)        &!&
                                           * beta_curr**2 )
        flux_pz = 0.d0
        flux_en = 0.5d0 * rho_curr * uxsk_grid(i)**3 * ( 1.d0             &!&
                                                  + 1.25d0*beta_curr**2 ) &!&
                 +  press_curr * uxsk_grid(i) * gam_sph/(gam_sph-1.d0)    &!&
                   * ( 1.d0 + beta_curr**2 )
        
      else if( (.not. l_oblique) .and. (.not. l_nonrel) ) then
        
        flux_px = gam_beta_curr**2 * ( rho_curr * ccgs**2                 &!&
                                      + gam_sph/(gam_sph-1.d0)*press_curr)&!&
                 + press_curr
        flux_pz = 0.d0
        flux_en = gam_beta_curr**2 *ccgs / (uxsk_grid(i)/ccgs)            &!&
                 * (rho_curr*ccgs**2 + gam_sph/(gam_sph-1.d0)*press_curr )
        
        ! Subtract mass-energy flux from flux_en to bring it in line with
        !   non-rel calculations
        flux_en = flux_en  -  gam_beta_curr*ccgs * rho_curr * ccgs**2
        
      else
        write(*,"(2A)") 'ERROR in "init_pop": fast push cannot handle ',  &!&
            'oblique shocks yet'
        write(*,"(A)") 'Stopping program now.'
        stop        
      endif
      !----------------------------------------------------------------------
      ! Fluxes calculated
      
      
      pxx_flux(i) = flux_px
      pxz_flux(i) = flux_pz
      en_flux(i)  = flux_en
      
    enddo  ! loop over grid location
    !------------------------------------------------------------------------
    ! Arrays updated through i_fast_stop
    
  endif  ! check on i_ion
  
  
!--! With fast push fluxes taken care of, create the particle distribution
   !   that will be injected at x_fast_stop
  T_or_E = tZ_ion(i_ion) * temp_ratio
  
  call set_in_dist(l_inj_wt, n_pts_inj, i_in_distr, T_or_E, aa,           &!&
     denZ_ion(i_ion), ptot_inj(:,i_ion), wt_inj(:,i_ion), n_pts_MB(i_ion))
  n_pts_use = n_pts_MB(i_ion)
  
  do i_prt = 1, n_pts_use
    xwt_in(i_prt)  = wt_inj(i_prt, i_ion)
    ptpf_in(i_prt) = ptot_inj(i_prt, i_ion)
    
    ! Per Vladimirov+ (2009) [PhD], particle velocities should not be
    !   isotropic in plasma frame.  Must be weighted according to shock-frame
    !   pitch angle.
    !------------------------------------------------------------------------
    call kiss(rand)

    if( l_nonrel ) then
      gam_ptpf = 1.d0
      vtpf     = ptpf_in(i_prt) / (aa*xmp)
      vmin_sf  = uxsk_grid(i_stop) - vtpf
      vmax_sf  = uxsk_grid(i_stop) + vtpf
      
      vxsf     = sqrt( (vmax_sf**2 - vmin_sf**2)*rand  +  vmin_sf**2 )
      vxpf     = vxsf - uxsk_grid(i_stop)
    else
      gam_ptpf = sqrt( 1.d0  +  (ptpf_in(i_prt)/(aa*xmp*ccgs))**2 )
      vtpf     = ptpf_in(i_prt) / (gam_ptpf * aa*xmp)
      vmin_sf  = (uxsk_grid(i_stop) - vtpf)                               &!&
                / ( 1.d0  -  uxsk_grid(i_stop)*vtpf/ccgs**2 )
      vmax_sf  = (uxsk_grid(i_stop) + vtpf)                               &!&
                / ( 1.d0  +  uxsk_grid(i_stop)*vtpf/ccgs**2 )
      
      vxsf     = sqrt( (vmax_sf**2 - vmin_sf**2)*rand  +  vmin_sf**2 )
      vxpf     = (vxsf - uxsk_grid(i_stop))                               &!&
                / ( 1.d0  -  vxsf*uxsk_grid(i_stop)/ccgs**2 )
    endif
    
    pbpf_in(i_prt) = gam_ptpf * aa*xmp * vxpf
    !------------------------------------------------------------------------
    ! Velocity-weighted pitch angles finished
    
    x_PT_cm_in(i_prt) = x_fast_stop_rg * rg0
    i_grid_in(i_prt)  = i_stop
  enddo
  
return
end subroutine init_pop


!****************************************************************************
!****************************************************************************
subroutine new_pcut(n_pts_target, n_saved, l_save, i_grid_sav, l_DwS_sav, &!&
     l_inj_sav, xwt_sav, ptpf_sav, pbpf_sav, x_PT_cm_sav, xn_per_sav,     &!&
     prp_x_cm_sav, acctime_sec_sav, phi_rad_sav, tcut_sav, n_pts_use,     &!&
     i_grid_new, l_DwS_new, l_inj_new, xwt_new, ptpf_new, pbpf_new,       &!&
     x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, phi_rad_new, &!&
     tcut_new, wt_running)

! Takes the **_sav arrays filled over the course of loop_pt and splits the
!   saved particles to form the population of the next pcut
!
! Input arguments:
!
! Output arguments:
!
! Input/output arguments:
!  1) wt_running: weight factor of each particle remaining after this pcut

use parameters, only: na_p

implicit none

  ! Input arguments
integer, intent(in) :: n_pts_target, n_saved
integer, dimension(na_p), intent(in) :: i_grid_sav, tcut_sav
logical, dimension(na_p), intent(in) :: l_save, l_DwS_sav, l_inj_sav
real(kind=8), dimension(na_p), intent(in) :: xwt_sav, ptpf_sav, pbpf_sav, &!&
     x_PT_cm_sav, xn_per_sav, prp_x_cm_sav, acctime_sec_sav, phi_rad_sav
  ! Output arguments
integer, dimension(na_p), intent(out) :: i_grid_new, tcut_new
logical, dimension(na_p), intent(out) :: l_DwS_new, l_inj_new
real(kind=8), dimension(na_p), intent(out) :: xwt_new, ptpf_new, pbpf_new,&!&
     x_PT_cm_new, xn_per_new, prp_x_cm_new, acctime_sec_new, phi_rad_new
  ! Input/output variables
integer, intent(inout) :: n_pts_use
real(kind=8), intent(inout) :: wt_running
  
  ! Local variables
integer :: i_mult, n_pts_new, i, j
real(kind=8) :: wt_fac
  
  
!--! Determine multiplicity of splitting; perhaps none needed
  i_mult = n_pts_target / n_saved
  if( i_mult .lt. 1 ) i_mult = 1   ! In case n_pts_target drops btwn pcuts
  
  
!--! Calculate effect on particle weights and the weighting factor of each
   !   remaining particle in the simulation
   !CHECKTHIS: is that last claim still true if old particles are imported
   !  into the simulation?
  wt_fac = 1.d0 / float(i_mult)
  wt_running = wt_running * wt_fac
  
  
!--! Perform the splitting
  n_pts_new = 0
  
  do j = 1, n_pts_use
    
    if( .not. l_save(j) ) cycle ! Don't multiply particles that weren't
                                  !   kept, obviously
      
    do i = 1, i_mult
      
      n_pts_new = n_pts_new + 1
      
      xwt_new(n_pts_new)         = xwt_sav(j) * wt_fac
      ptpf_new(n_pts_new)        = ptpf_sav(j)
      pbpf_new(n_pts_new)        = pbpf_sav(j)
      x_PT_cm_new(n_pts_new)     = x_PT_cm_sav(j)
      i_grid_new(n_pts_new)      = i_grid_sav(j)
      l_DwS_new(n_pts_new)       = l_DwS_sav(j)
      l_inj_new(n_pts_new)       = l_inj_sav(j)
      xn_per_new(n_pts_new)      = xn_per_sav(j)
      prp_x_cm_new(n_pts_new)    = prp_x_cm_sav(j)
      acctime_sec_new(n_pts_new) = acctime_sec_sav(j)
      phi_rad_new(n_pts_new)     = phi_rad_sav(j)
      tcut_new(n_pts_new)        = tcut_sav(j)
      
    enddo  ! loop over splits
  enddo  ! loop over saved particles
  
  n_pts_use = n_pts_new
  
  
!--! Zero out the unused portions of the *_new arrays, just as a precaution
      xwt_new(n_pts_use+1:)         = 0.d0
      ptpf_new(n_pts_use+1:)        = 0.d0
      pbpf_new(n_pts_use+1:)        = 0.d0
      x_PT_cm_new(n_pts_use+1:)     = 0.d0
      i_grid_new(n_pts_use+1:)      = 0
      l_DwS_new(n_pts_use+1:)       = .false.
      l_inj_new(n_pts_use+1:)       = .false.
      xn_per_new(n_pts_use+1:)      = 0.d0
      prp_x_cm_new(n_pts_use+1:)    = 0.d0
      acctime_sec_new(n_pts_use+1:) = 0.d0
      phi_rad_new(n_pts_use+1:)     = 0.d0
      tcut_new(n_pts_use+1:)        = 0
  
return
end subroutine new_pcut


!****************************************************************************
!****************************************************************************
subroutine particle_finish(aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad, uxsk,&!&
     uzsk, utot, gam_usf, b_cos_th, b_sin_th, i_reason, xwt,              &!&
     esc_psd_feb_DwS, esc_psd_feb_UpS, esc_flux, px_esc_feb, en_esc_feb,  &!&
     esc_en_eff, esc_num_eff)

! Handles particles that leave the system during loop_helix for any reason.
!
! Input arguments:
!  1) aa: particle atomic mass
!  2) pbpf: component of ptpf parallel to magnetic field
!  3) p_perp_b_pf: component of ptpf perpendicular to magnetic field
!  4) gam_ptpf: Lorentz factor associated with ptpf
!  5) phi_rad: phase angle of gyration; looking UpS, counts clockwise from
!    +z axis
!  6) uxsk: bulk flow speed along x axis
!  7) uzsk: bulk flow speed along z axis
!  8) utot: total bulk flow speed
!  9) gam_usf: Lorentz factor associated with utot
!  10) b_cos_th: component of magnetic field along x axis
!  11) b_sin_th: component of magnetic field along z axis
!  12) i_reason: integer reason for why particle left system
!  13) xwt: particle's weight
! Output arguments: none; modifies arrays from modules as needed

use constants
use parameters, only: na_i, na_its, psd_max, en_rel_pt
use iteration_vars, only: i_itr
use species_vars, only: i_ion

implicit none

  ! Input arguments
integer, intent(in) :: i_reason
real(kind=8), intent(in) :: aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad,     &!&
     uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th, xwt
  ! Output arguments
  ! Input/output arguments
real(kind=8), dimension(na_i), intent(inout) :: esc_flux
real(kind=8), dimension(na_i, na_its), intent(inout) :: px_esc_feb,       &!&
     en_esc_feb
real(kind=8), dimension(0:psd_max, na_i), intent(inout) :: esc_en_eff,    &!&
     esc_num_eff
real(kind=8), dimension(0:psd_max, 0:psd_max), intent(inout) ::           &!&
     esc_psd_feb_UpS, esc_psd_feb_DwS

  ! Local variables
integer :: i_pt, j_th
real(kind=8), parameter :: spike_away = 1000.d0
real(kind=8) :: ptsk, pxsk, pzsk, gam_ptsk, wtfac, en_flux_add


!--! Transform plasma frame momentum into shock frame for binning
  call Xform_p_PS(aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad, uxsk, uzsk,   &!&
     utot, gam_usf, b_cos_th, b_sin_th, ptsk, pxsk, pzsk, gam_ptsk)
  
  
!--! Get PSD bins for this particle
  call get_psd_bins(pxsk, ptsk, .true., i_pt, j_th)
   
   
!--! Now take additional action based on *how* the particle left the grid
  select case( i_reason )
    case( 1 )
    ! Particle escape: DwS, with or without scattering enabled
      
      if( ptsk .gt. abs(spike_away*pxsk) ) then
        wtfac = gam_ptsk * aa*xmp * spike_away / ptsk
      else
        wtfac = gam_ptsk * aa*xmp / abs(pxsk)
      endif
      
!$omp atomic
      esc_psd_feb_DwS(i_pt, j_th) = esc_psd_feb_DwS(i_pt, j_th)           &!&
                                   +  xwt * wtfac
      
      
    case( 2 )
    ! Particle escape: pmax, UpS FEB, transverse distance
      
      if( ptsk .gt. abs(spike_away*pxsk) ) then
        wtfac = gam_ptsk * aa*xmp * spike_away / ptsk
      else
        wtfac = gam_ptsk * aa*xmp / abs(pxsk)
      endif
      
!$omp atomic
      esc_flux(i_ion) = esc_flux(i_ion) + xwt
      
!$omp atomic
      esc_psd_feb_UpS(i_pt, j_th) = esc_psd_feb_UpS(i_pt, j_th)           &!&
                                   +  xwt * wtfac
      
      if( (gam_ptsk - 1.d0) .lt. (en_rel_pt/(aa*rm_prot)) ) then
        en_flux_add = ptsk**2 / (2.d0 * aa*xmp) * xwt
      else
        en_flux_add = (gam_ptsk - 1.d0) * aa*rm_prot * xwt
      endif
      
      ! Update escape arrays that will be averaged over consecutive
      !   iterations
!$omp atomic
      px_esc_feb(i_ion, i_itr) = px_esc_feb(i_ion, i_itr)  +  abs(pxsk) * xwt
!$omp atomic
      en_esc_feb(i_ion, i_itr) = en_esc_feb(i_ion, i_itr)  +  en_flux_add
      
      ! Update escape arrays to be printed out with spectral information
!$omp atomic
      esc_en_eff(i_pt, i_ion)  = esc_en_eff(i_pt, i_ion) + en_flux_add
!$omp atomic
      esc_num_eff(i_pt, i_ion) = esc_num_eff(i_pt, i_ion) + xwt
      
    
    case( 3 )
    ! Particle escape: age_max
      
      !DOLATER: write out particles to be read in during a later run, i.e.
      !         a pre-existing population of CRs.  Also include new keyword
      !         for this purpose
      
      
    case( 4 )
    ! Zero energy after radiative losses
      
      ! Do nothing
      
      
    case default
      write(*,"(2A,I0,A)") 'ERROR in particle finish: unknown i_reason ', &!&
          'passed:',i_reason,'.  Can only handle 1-4'
      write(*,"(A)") 'Stopping program now.'
      stop
      
  end select

return
end subroutine particle_finish


!****************************************************************************
!****************************************************************************
subroutine print_input(n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,      &!&
     num_psd_mom_bins, num_psd_tht_bins, n_xspec, n_pcuts, n_grid, rRH,   &!&
     r_comp, u_Z, beta_Z, gam_Z, u_2, beta_2, gam_2, denZ_ion, bmag_Z,    &!&
     bmag_2, theta_BZ, theta_B2, theta_u2, aa_ion, tZ_ion, tZ_elec,       &!&
     mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine, feb_UpS,        &!&
     feb_DwS, rg0, age_max, en_pcut_hi, do_fast_push, bturb_comp_frac)

! Prints a whole mess of information to the screen and to mc_out
!
! Input arguments:
!  1) n_pts_inj: target # of particles for injection distribution
!  2) n_pts_pcut: target # of particle for low-E pcuts
!  3) n_pts_pcut_hi: target # of particles for hi-E pcuts
!  2) n_ions: # of ion species run
!  3) num_psd_***_bins: # of bins in momentum/theta directions in psd
!  4) n_xspec: # of additional locations to track particle spectra
!  5) n_pcuts: # of pcuts
!  6) n_grid: # of grid zone BOUNDARIES, including x_grid_stop
!  7) rRH: Rankine-Hugoniot compression ratio for this shock
!  8) r_comp: the compression ratio being used
!  9) u/beta/gam _Z/_2: total fluid velocity[cm/s, /c] and Lorentz factor
!    for far UpS and DwS regions
!  10) 
! No output arguments (it's all to the screen/file)

use constants
use parameters, only: na_p, na_i, psd_max, na_g, na_c
use controls, only: l_sc_elec

implicit none

  ! Input arguments
integer, intent(in) :: n_pts_inj, n_pts_pcut, n_pts_pcut_hi, n_ions,      &!&
     num_psd_mom_bins, num_psd_tht_bins, n_xspec, n_pcuts, n_grid
logical, intent(in) :: do_fast_push
real(kind=8), intent(in) :: rRH, r_comp, u_Z, beta_Z, gam_Z, u_2, beta_2, &!&
     gam_2, bmag_Z, bmag_2, theta_BZ, theta_B2, theta_u2, tZ_elec,        &!&
     mach_sonic, mach_alfven, xn_per_coarse, xn_per_fine, feb_UpS,        &!&
     feb_DwS, rg0, age_max, en_pcut_hi, bturb_comp_frac
real(kind=8), dimension(na_i), intent(in) :: denZ_ion, aa_ion, tZ_ion

  ! Local variables
integer :: n_pts_max, i_elec
real(kind=8) :: temp_elec, sig_KW


!--! Print array parameters & usage
  n_pts_max = max(n_pts_inj, n_pts_pcut, n_pts_pcut_hi)
  write(*,*)
  write(*,"(2X,A)") 'Array parameters/usage:'
  write(*,"(3(A,I0))") '        na_p = ', na_p, '      na_i = ', na_i,    &!&
      '        psd_max = ',psd_max
  write(*,"(4(A,I0))") '      n_pts_max = ',n_pts_max, '    n_ions = ',   &!&
      n_ions, '   psd_mom/tht_bins = ', num_psd_mom_bins,'/',num_psd_tht_bins
  
  write(9,*)
  write(9,"(2X,A)") 'Array parameters/usage:'
  write(9,"(3(A,I0))") '        na_p = ', na_p, '      na_i = ', na_i,    &!&
      '        psd_max = ',psd_max
  write(9,"(4(A,I0))") '      n_pts_max = ',n_pts_max, '    n_ions = ',   &!&
      n_ions, '   psd_mom/tht_bins = ', num_psd_mom_bins,'/',num_psd_tht_bins
  
  write(*,*)
  write(*,"(3(A,I0))") '      na_g = ', na_g, '      na_c = ', na_c,      &!&
      '       na_g = ', na_g
  write(*,"(3(A,I0))") '     n_xpsec = ', n_xspec, '     n_pcuts = ',     &!&
      n_pcuts, '     n_grid = ', n_grid
  
  write(9,*)
  write(9,"(3(A,I0))") '      na_g = ', na_g, '      na_c = ', na_c,      &!&
      '       na_g = ', na_g
  write(9,"(3(A,I0))") '     n_xpsec = ', n_xspec, '     n_pcuts = ',     &!&
      n_pcuts, '     n_grid = ', n_grid


!--! Shock speeds, escaping fluxes, and key densities
  write(*,*)
  write(*,"(2X,2(A,ES9.3E2))") 'rRH = ', rRH, '       r_comp = ', r_comp
  
  write(9,*)
  write(9,"(2X,2(A,ES9.3E2))") 'rRH = ', rRH, '       r_comp = ', r_comp
  
  write(*,*)
  write(*,"(3(A,ES9.3E2))") '  u_Z[cm/s] = ', u_Z, '    u_2[cm/s] = ',    &!&
      u_2
  write(*,"(3(A,ES9.3E2))") '     beta_Z = ', beta_Z, '       beta_2 = ', &!&
      beta_2
  write(*,"(3(A,ES9.3E2))") '      gam_Z = ', gam_Z, '        gam_2 = ',  &!&
      gam_2
  write(*,"(3(A,ES9.3E2))") '  rho_Z[prot/cm^3] = ', denZ_ion(1),         &!&
      '     rho_2[prot/cm^3] = ', denZ_ion(1)*gam_Z*beta_Z/(gam_2*beta_2)
  
  write(9,*)
  write(9,"(3(A,ES9.3E2))") '  u_Z[cm/s] = ', u_Z, '    u_2[cm/s] = ',    &!&
      u_2
  write(9,"(3(A,ES9.3E2))") '     beta_Z = ', beta_Z, '       beta_2 = ', &!&
      beta_2
  write(9,"(3(A,ES9.3E2))") '      gam_Z = ', gam_Z, '        gam_2 = ',  &!&
      gam_2
  write(9,"(3(A,ES9.3E2))") '  rho_Z[prot/cm^3] = ', denZ_ion(1),         &!&
      '     rho_2[prot/cm^3] = ', denZ_ion(1)*gam_Z*beta_Z/(gam_2*beta_2)


!--! Relevant angles and field strengths
  write(*,*)
  write(*,"(2(A,ES9.3E2))") '      bmag_Z[G] = ', bmag_Z,                 &!&
      '             bmag_2[G] = ', bmag_2
  write(*,"(2(A,ES9.3E2))") '  theta_BZ[deg] = ', theta_BZ,               &!&
      '   theta_B2[deg](calc) = ', theta_B2
  write(*,"(A,ES9.3E2)") '  theta_u2[deg](calc) = ', theta_u2
  
  write(9,*)
  write(9,"(2(A,ES9.3E2))") '      bmag_Z[G] = ', bmag_Z,                 &!&
      '             bmag_2[G] = ', bmag_2
  write(9,"(2(A,ES9.3E2))") '  theta_BZ[deg] = ', theta_BZ,               &!&
      '   theta_B2[deg](calc) = ', theta_B2
  write(9,"(A,ES9.3E2)") '  theta_u2[deg](calc) = ', theta_u2


!--! Temperatures and Mach numbers
  if( l_sc_elec ) then
    i_elec = minloc( aa_ion(1:n_ions), 1 )
    temp_elec = tZ_ion(i_elec)
  else
    temp_elec = tZ_elec
  endif
  write(*,*)
  write(*,"(3(A,ES9.3E2))") '  Temp_Z[K](prot) = ', tZ_ion(1),            &!&
      '     Temp_Z[K](elec) = ', temp_elec
  write(*,"(3(A,ES9.3E2))") '  Mach(sonic) = ', mach_sonic,               &!&
      '     Mach(Alfven) = ', mach_alfven
  
  write(9,*)
  write(9,"(3(A,ES9.3E2))") '  Temp_Z[K](prot) = ', tZ_ion(1),            &!&
      '     Temp_Z[K](elec) = ', temp_elec
  write(9,"(3(A,ES9.3E2))") '  Mach(sonic) = ', mach_sonic,               &!&
      '     Mach(Alfven) = ', mach_alfven


!--! Divisions of gyroperiod
  write(*,*)
  write(*,"(2(A,ES9.3E2))") '  N_g(coarse) = ', xn_per_coarse,            &!&
      '    N_g(fine) = ', xn_per_fine
  write(9,*)
  write(9,"(2(A,ES9.3E2))") '  N_g(coarse) = ', xn_per_coarse,            &!&
      '    N_g(fine) = ', xn_per_fine


!--! FEB info and max age
  write(*,*)
  write(*,"(2(A,ES10.3E2))") '  UpS FEB[rg0] = ', feb_UpS/rg0,             &!&
      '   UpS FEB[pc] = ', feb_UpS / pc_to_cm
  write(*,"(3(A,ES10.3E2))") '  DwS FEB[rg0] = ', feb_DwS/rg0,             &!&
      '   DwS FEB[pc] = ', feb_DwS / pc_to_cm
  write(*,"(A,ES10.3E2)") '  Max CR age[sec] = ', age_max
  
  write(9,*)
  write(9,"(2(A,ES10.3E2))") '  UpS FEB[rg0] = ', feb_UpS/rg0,             &!&
      '   UpS FEB[pc] = ', feb_UpS / pc_to_cm
  write(9,"(3(A,ES10.3E2))") '  DwS FEB[rg0] = ', feb_DwS/rg0,             &!&
      '   DwS FEB[pc] = ', feb_DwS / pc_to_cm
  write(9,"(A,ES10.3E2)") '  Max CR age[sec] = ', age_max


!--! Test-particle index from Keshet & Waxman (2005) [2005PhRvL..94k1102K]
  sig_KW = ( 3.d0*beta_Z  -  2.d0*beta_Z*beta_2**2  +  beta_2**3 )         &!&
          / (beta_Z - beta_2)
  write(*,*)
  write(*,"(A,ES9.3E2)") '  Keshet & Waxman (2005) index = ', sig_KW
  
  write(9,*)
  write(9,"(A,ES9.3E2)") '  Keshet & Waxman (2005) index = ', sig_KW


!--! Energy to switch between low- and high-pcut particle counts
  write(*,*)
  write(*,"(A,ES9.3E2)") '  High pcut energy[keV/aa] = ',en_pcut_hi
  write(*,*)
  
  write(9,*)
  write(9,"(A,ES9.3E2)") '  High pcut energy[keV/aa] = ',en_pcut_hi
  write(9,*)


!--! Finally a warning about possible complications later
  if( do_fast_push .and. (bturb_comp_frac .gt. 0.d0) ) then
    write(*,*)
    write(*,"(2A)") ' WARNING: both fast push and amplified B-field ',    &!&
        'turbulence in use'
    write(*,"(A)") ' Check flux equations in "init_pop" for consistency'
    write(*,*)
    
    write(9,*)
    write(9,"(2A)") ' WARNING: both fast push and amplified B-field ',    &!&
        'turbulence in use'
    write(9,"(A)") ' Check flux equations in "init_pop" for consistency'
    write(9,*)
  endif

return
end subroutine print_input


!****************************************************************************
!****************************************************************************
subroutine print_plot_vals(iunit)

! Prints the long list of floats at the end of each data set that are read in
!   by the plotting program to display information about the run.
!
! Input arguments:
!  1) iunit: unit number to which the line will be written

use constants
use parameters, only:
use controls, only: n_pts_inj, n_pts_pcut, do_fast_push, i_in_distr,      &!&
     dont_DSA, u_Z, gam_Z, r_comp, rRH, theta_BZ, theta_B2, theta_u2,     &!&
     bmag_Z, feb_UpS, rg0, Emax_keV, Emax_keV_per_aa, pmax_cgs,           &!&
     xn_per_coarse, xn_per_fine, mach_sonic, mach_alfven, x_grid_start_rg,&!&
     x_grid_stop_rg, x_fast_stop_rg, eta_mfp, x_art_start_rg, x_art_scale,&!&
     feb_DwS, jet_rad_pc, jet_sph_frac, jet_dist_kpc, n_ions, aa_ion,     &!&
     zz_ion, denZ_ion, tZ_ion, smooth_mom_en_fac, en_inj,                 &!&
     smooth_press_flux_psd_fac, en_xfer_frac
use randnum, only: iseed_in

implicit none

  ! Input arguments
integer, intent(in) :: iunit

  ! Local variables
integer :: iannt, idum
real(kind=8) :: x_pts_inj, x_pts_pcut, xseed, x_fast_push, x_in_distr,    &!&
     x_DSA, x_ions


      x_pts_inj  = real(n_pts_inj)
      x_pts_pcut = real(n_pts_pcut)
      xseed      = float(iseed_in)
      if( do_fast_push ) then
        x_fast_push = 66.d0
      else
        x_fast_push = 0.d0
      endif
      x_in_distr = float(i_in_distr)
      if( dont_DSA ) then
        x_DSA = 66.d0
      else
        x_DSA = 0.d0
      endif
      x_ions = float(n_ions)

      iannt = 3333
      idum  = 333
  
  ! WARNING: these column numbers are reused in both subroutine read_old_prof
  !   and the plotting program pg_color.f90.  If they are ever changed,
  !   modify the other codes accordingly!
  write(iunit,"(2I4,100ES12.3E2)") iannt, idum, &!&
                u_Z/1.d5,                       &!& ! 1   
                gam_Z,                          &!& ! 2
                r_comp,                         &!& ! 3
                rRH,                            &!& ! 4
                theta_BZ,                       &!& ! 5
                theta_B2,                       &!& ! 6
                theta_u2,                       &!& ! 7
                bmag_Z,                         &!& ! 8
                feb_UpS/rg0,                    &!& ! 9
                Emax_keV,                       &!& ! 10
                Emax_keV_per_aa,                &!& ! 11
                pmax_cgs/(xmp*ccgs),            &!& ! 12
                x_pts_inj,                      &!& ! 13
                x_pts_pcut,                     &!& ! 14
                xn_per_coarse,                  &!& ! 15
                xn_per_fine,                    &!& ! 16
                mach_sonic,                     &!& ! 17
                mach_alfven,                    &!& ! 18
                x_grid_start_rg,                &!& ! 19
                xseed,                          &!& ! 20
                x_grid_stop_rg,                 &!& ! 21
                x_fast_push,                    &!& ! 22
                x_fast_stop_rg,                 &!& ! 23
                eta_mfp,                        &!& ! 24
                x_art_start_rg,                 &!& ! 25
                x_art_scale,                    &!& ! 26
                feb_DwS/rg0,                    &!& ! 27
                jet_rad_pc,                     &!& ! 28
                jet_sph_frac,                   &!& ! 29
                jet_dist_kpc,                   &!& ! 30
                smooth_mom_en_fac,              &!& ! 31
                x_in_distr,                     &!& ! 32
                en_inj,                         &!& ! 33
                smooth_press_flux_psd_fac,      &!& ! 34
                x_DSA,                          &!& ! 35
                en_xfer_frac,                   &!& ! 36

                x_ions,                         &!&
                aa_ion(1:n_ions),               &!&
                zz_ion(1:n_ions),               &!&
                denZ_ion(1:n_ions),             &!&
                tZ_ion(1:n_ions)

return
end subroutine print_plot_vals


!****************************************************************************
!****************************************************************************
subroutine print_progress_bar(i_prt, n_pts_use)

! Displays a counter showing current particle number and percentage of
!  current pcut complete.
! (Thanks to Chris Mauney for the original version of this.)
!
! Inputs:
!  1) i_prt: current particle number
!  2) nparts_hi: number of particles in this pcut
! Outputs:
!  Message to screen

implicit none

  ! Input variables
integer, intent(in) :: i_prt, n_pts_use
  ! Local variables
real :: percent_done

  percent_done = i_prt * real(1.0/n_pts_use)
  
  open(6)!open(unit=6,carriagecontrol='fortran')
  write(6,100,advance='no') int(percent_done*1.0e2), i_prt, n_pts_use,    &!&
      char(13)
  close(6)
  
  100 format("Progress: ",I3,"% (",I0,"/",I0," complete)",A)

return
end subroutine


!****************************************************************************
!****************************************************************************
subroutine prob_return(rad_loss_fac, B_CMBz, x_PT_old, aa, zz, gyro_denom,&!&
     i_return, x_PT_cm, prp_x_cm, ptpf, gam_ptpf, pbpf, p_perp_b_pf,      &!&
     acctime_sec, phi_rad, lose_pt, helix_count, pcut_prev, xwt, tcut_curr)

! If the particle ends its movement downstream of the shock, perform a series
!   of tests to determine whether it will be culled from the simulation.
!
! Input arguments:
!  1) rad_loss_fac: constant related to radiative losses; only used to pass
!    to retro_time when needed
!  2) B_CMBz: effective magnetic field due to CMB at redshift of source; only
!    used to pass to retro_time when needed
!  3) x_PT_old: particle position before most recent move
!  4) aa: atomic mass number of ion species
!  5) gyro_denom: denominator of gyroradius fraction, zz*qcgs*bmag
!  6) helix_count: counter for number of times through main propagation loop
!    for current particle
!  7) i_cut: current pcut; needed when electrons are undergoing radiative
!    losses
!  8) pcut_prev: momentum of previous pcut; needed when electrons are
!    undergoing radiative losses
!  9) xwt: current particle weight; passed to retro_time if called
! Output arguments:
!  1) i_return: flag for fate of particle; see top of loop_helix
!  2) lose_pt: "T" if particle hit zero energy due to radiative losses while
!    in DwS region
! Input/output arguments:
!  1) x_PT_cm: position of particle after recent motion & DwS adjustment
!  2) prp_x_cm: location of probability-of-return plane
!  3) ptpf: total plasma frame momentum of particle
!  4) gam_ptpf: Lorentz factor associated with ptpf
!  5) gyro_denom: denominator of gyroradius fraction, zz*qcgs*bmag
!  6) pbpf/p_perp_b_pf: components of ptpf parallel/perpendicular to B field
!  7) acctime_sec: total accumulated acceleration time
!  8) phi_rad: phase angle of particle's gyration
!  9) tcut_curr: current tcut for particle tracking; passed to retro_time if
!    called

use constants
use controls, only: x_grid_stop, u_2, use_custom_epsB, eta_mfp, do_retro, &!&
     bmag_2
use randnum

implicit none

  ! Input arguments
integer, intent(in) :: helix_count
real(kind=8), intent(in) :: rad_loss_fac, B_CMBz, x_PT_old, aa, zz,       &!&
     pcut_prev, xwt
  ! Ouput arguments
integer, intent(out) :: i_return
logical, intent(out) :: lose_pt
  ! Input/output arguments
integer, intent(inout) :: tcut_curr
real(kind=8), intent(inout) :: x_PT_cm, prp_x_cm, ptpf, gam_ptpf,         &!&
     gyro_denom, pbpf, p_perp_b_pf, acctime_sec, phi_rad
  
  ! Local variables
real(kind=8) :: vpt_o_u2, gyro_denom_tmp, gyro_rad_tot_cm, L_diff, vtpf,  &!&
     prob_ret, rand, prp_fac
  
  
  i_return = 2   ! Presume particle didn't enter probability of return
                 !   calculation; change later as needed
  
  
!--! Test whether particle is still UpS from PRP position; don't bother with
   !   return calculations if so
  if( x_PT_cm .lt. x_grid_stop ) then
    
    ! Do nothing
    
    
!--! Particle has just crossed end of shock region as initially defined in
   !   input file
  else if( (x_PT_old .lt. x_grid_stop).and.(x_PT_cm .ge. x_grid_stop) ) then
    
    ! The following simple equation is decidedly non-trivial, and comes from
    !   two assumptions:
    !  1) the particle's diffusion coefficient D may be described
    !   by D = 1/3 * eta_mfp * r_g * v_pt (by default; a different
    !   f(r_g) may be specified as desired in place of eta*r_g)
    !  2) the relation between the diffusion coefficient D and the
    !   diffusion length L is L = D/<u>, where <u> is the average
    !   speed of diffusion.  Assuming isotropic particles in the
    !   DwS frame, <u> = u_2 since the average thermal *velocity*
    !   of the population is 0.
    ! The calculation of gyro_denom_tmp ensures that even particles that
    !   started UpS of shock still use the DwS magnetic field for their
    !   diffusion length
    !DOLATER: include f(r_g) in place of eta*r_g to allow for
    !  arbitrary diffusion
    vpt_o_u2        = ptpf/(aa*xmp*gam_ptpf * u_2)
    if( use_custom_epsB .and. (x_PT_cm .gt. x_grid_stop) ) then
      gyro_denom_tmp  = ( x_grid_stop / x_PT_cm )**0.25  /  (qcgs * bmag_2)
    else
      gyro_denom_tmp  = 1.d0 / (qcgs * bmag_2)
    endif
    gyro_rad_tot_cm = ptpf * ccgs * gyro_denom_tmp
    L_diff          = third * eta_mfp * gyro_rad_tot_cm * vpt_o_u2
    
    ! Make absolutely sure particles will have enough distance to
    !   isotropize before encountering PRP; allow for three diffusion
    !   lengths beyond *current position*, not just beyond end of grid
    prp_x_cm = x_PT_cm  +  3.d0 * L_diff
    
    
!--! Particle has crossed PRP, and we must do more complex calculations to
   !  determine if it returns
  else if( (x_PT_old .lt. prp_x_cm) .and. (x_PT_cm .ge. prp_x_cm) ) then
    
    vtpf     = ptpf / (gam_ptpf * aa*xmp)
    prob_ret = ( (vtpf - u_2) / (vtpf + u_2) )**2
    
    
  !--! If the particle's plasma frame velocity is less than u_2, or if the
     !   probability of return calculation (see Jones & Ellison 1991
     !   [1991SSRv...58..259J]) fails, the particle will not return from
     !   the DwS region.
    call kiss(rand)
    if( (vtpf .lt. u_2) .or. (rand .gt. prob_ret) ) then
    
      i_return = 0
      
      
  !--! Particle will return from DwS region.  Either analytically determine
     !   its properties upon return or use retro-time calculation
    else
      
      i_return = 1
    
    !--! Track particle histories "explicitly" (see note in subroutine)
      if( do_retro ) then
        
        call retro_time(rad_loss_fac, B_CMBz, aa, zz, gyro_denom,         &!&
           prp_x_cm, ptpf, pbpf, p_perp_b_pf, gam_ptpf, acctime_sec,      &!&
           phi_rad, lose_pt, xwt, tcut_curr)
        
        ! If electrons somehow lost all their energy due to radiative losses,
        !   flag for removal
        if( lose_pt ) i_return = 0
        
        ! Particles return from retro_time at the location of the PRP
        x_PT_cm = prp_x_cm
        
        
    !--! Analytically place particles back at PRP
      else
        
        !DOLATER: check relations in Appendix A3 of Ellison, Baring & Jones
        !   (1993) [1996ApJ...473.1029E] to verify that they are correct in
        !   relativistic limit.
        write(*,"(2A)") 'ERROR in prob_ret: code not set up for ',        &!&
            'analytical PRP calculations.'
        write(*,"(2A)") 'Must verify that EBJ1996 relations are correct ',&!&
            'in rel. case'
        write(*,"(A)") 'Stopping program now.'
        stop
        
      endif    
    
    endif ! check on prob_ret
    
  else
    
  !--! Particle is DwS from grid end, but didn't cross it this time step.
     !   Also, it is UpS from the PRP.  However, electrons experiencing
     !   radiative losses may have a smaller L_diff during their propagation,
     !   so a shorter PRP can be used.  Test for that here.
     ! Two methods for doing that:
     !  (1) If L_diff has dropped sufficiently, move the PRP far UpS from the
     !    particle's current position.  The particle will be culled at the
     !    next time step
     !  (2) Otherwise, calculate a new PRP location based on ratio of
     !    current momentum to minimum mometum for this pcut.  The strong
     !    dependence on momentum (p^5) is so that these electrons have time
     !    to isotropize DwS, even though the bulk of their motion occurred at
     !    a much higher energy and therefore mean free path
    if( (aa .lt. 1.d0) .and. (ptpf .lt. pcut_prev) ) then
      if( mod(helix_count, 1000) .eq. 0 ) then
        vpt_o_u2 = ptpf/(aa*xmp*gam_ptpf * u_2)
        gyro_rad_tot_cm = ptpf * ccgs * gyro_denom
        L_diff = third * eta_mfp * gyro_rad_tot_cm * vpt_o_u2
        
        if( x_PT_cm .gt. (2.d3*L_diff) ) then
          prp_x_cm = 0.8 * x_PT_cm
        else
          prp_fac = (pcut_prev / ptpf)**5
          if( prp_x_cm .gt. (x_grid_stop + L_diff*prp_fac) )              &!&
              prp_x_cm = x_grid_stop  +  L_diff * prp_fac
        endif
        
      endif
    endif
    
  
  endif ! check on position vs x_grid_stop & prp_x_cm

return
end subroutine prob_return


!****************************************************************************
!****************************************************************************
subroutine q_esc_calcs(gam_adiab, q_esc_cal_px, q_esc_cal_en)

! Use the Rankine-Hugoniot relations to calculate the escaping momentum &
!   energy flux.

use constants
use parameters, only: na_i, beta_rel_fl
use controls, only: r_comp, rRH, u_Z, beta_Z, gam_Z, l_oblique, n_ions,   &!&
     denZ_ion, tZ_ion, aa_ion, zz_ion, l_sc_elec, tZ_elec, gam_2, beta_2, u_2

implicit none

  ! Input arguments
real(kind=8), intent(in) :: gam_adiab
  ! Output arguments
real(kind=8), intent(out) :: q_esc_cal_px, q_esc_cal_en

  ! Local variables
integer :: j
logical :: l_nonrel
real(kind=8) :: gambeta_2, gam_fac, press_Z, rho_Z, den_elec, rho_2,      &!&
     press_2, Qen, q_fac, F_px_UpS_fl, F_en_UpS_fl, term_1, term_2, Qpx

  
  ! Quick test first.  If r_comp = rRH, we expect no escaping flux.
  !--------------------------------------------------------------------------
  if( r_comp .eq. rRH ) then
    q_esc_cal_px = 0.d0
    q_esc_cal_en = 0.d0
    return
  endif
  !--------------------------------------------------------------------------
  ! r_comp = rRH check
  
  
  !--------------------------------------------------------------------------
  !  Four possibilities for R-H relations: nonrel/rel and parallel/oblique.
  !    Determine which of the four to use.  Cutoff for nonrel/rel is set in
  !    module 'controls'
  !--------------------------------------------------------------------------
  if( beta_Z .lt. beta_rel_fl ) then
    l_nonrel = .true.
  else
    l_nonrel = .false.
  endif
  
  
  gambeta_2 = gam_2 * beta_2
  
  gam_fac   = gam_adiab / (gam_adiab - 1.d0)
      

  !--------------------------------------------------------------------------
  !  Possibility 1: Nonrel, parallel
  !    Solution comes from Ellison (1985) [1985JGR....90...29E] (Eqs 8-10).
  !    Note assumption of zero escaping momentum flux, which is good to
  !    within a couple percent for strong nonrel shocks.
  !    #DOLATER: check how much of a difference this assumption makes
  !--------------------------------------------------------------------------
  if( l_nonrel  .and.  (.not. l_oblique) ) then
    
    ! Calculate thermal pressure of far upstream gas
    press_Z  = 0.d0
    rho_Z    = 0.d0
    den_elec = 0.d0
    do j = 1, n_ions
      press_Z  =  press_Z  +  denZ_ion(j) * xkb*tZ_ion(j)
      rho_Z    =  rho_Z    +  denZ_ion(j) * xmp*aa_ion(j)
      
      if( aa_ion(j) .ge. 1.d0 ) den_elec  =  den_elec + denZ_ion(j)*zz_ion(j)
    enddo
    
    ! If electrons were not a separate species, add them in here
    if( .not. l_sc_elec ) then
      press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
      rho_Z    =  rho_Z    +  den_elec * xme
    endif
    
    ! Calculate UpS incoming energy flux
    ! #assumecold
    F_px_UpS_fl = rho_Z * u_Z**2  +  press_Z
    F_en_UpS_fl = 0.5d0 * rho_Z*u_Z**3  +  2.5d0 * press_Z * u_Z
    
    ! Calculate far DwS density (Eq 8) and pressure (Eq 9)
    rho_2   = rho_Z * gam_Z*beta_Z / gambeta_2
    press_2 = F_px_UpS_fl  -  rho_2*u_2**2
    
    ! Calculate escaping energy flux using Eq (10)
    Qen = F_en_UpS_fl  -  0.5d0*rho_Z*u_Z*u_2**2  -  press_2 * u_2 * gam_fac
    
    ! Finally, put in units of F_enZ
    q_esc_cal_en = Qen / F_en_UpS_fl
    q_esc_cal_px = 0.d0
    
    
  !--------------------------------------------------------------------------
  !  Possibility 2: Relativistic, parallel
  !   Solution comes from Ellison+ (1990) [1991ApJ...378..214E].  Uses
  !     relativistic Rankine-Hugoniot relations.  See that paper for details
  !     of equations and associated quantities.  Briefly,
  !        R-H1:             g0 * n0 * b0  =  g2 * n2 * b2
  !        R-H2:  g0^2 * w0 * b0^2  +  P0  =  g2^2 * w2 * b2^2  +  P2  + Qpx
  !        R-H3:       g0^2 * w0 * b0 * c  =  g2^2 * w2 * b2 * c       + Qen
  !     where
  !        w    = E_rm  +  E_ke  +  P,  <--- enthalpy as total energy density
  !                                           + pressure
  !        E_rm = n * m * c^2           <--- rest mass energy density
  !        E_ke = n * m * c^2 * e(p)    <---  kinetic  energy density, with
  !           e(p) =  sqrt( 1 + (p/mc)^2 )  -  1
  !        P    = 1/3 * n * p * v       <--- pressure
  !
  !   For closure, it is assumed that the two escaping fluxes are related by
  !          Qen = sqrt[0.5*(1+beta_Z)*ccgs^2] * Qpx ,
  !     i.e. the geometric mean of the arithmetic mean of u_Z and c.  This
  !     allows the solution to smoothly join with the non-rel version.
  !   Use only fluid component of fluxes, not fluid+EM, for now.
  !--------------------------------------------------------------------------
  else if( (.not. l_nonrel)  .and.  (.not. l_oblique) ) then
    
    ! Factor relating Qen and Qpx
    q_fac = sqrt( 0.5d0 * (1.d0 + beta_Z) * ccgs**2 )
    
    ! Calculate thermal pressure of far upstream gas
    press_Z  = 0.d0
    rho_Z    = 0.d0
    den_elec = 0.d0
    do j = 1, n_ions
      press_Z  =  press_Z  +  denZ_ion(j) * xkb*tZ_ion(j)
      rho_Z    =  rho_Z    +  denZ_ion(j) * xmp*aa_ion(j)
      
      if( aa_ion(j) .ge. 1.d0 ) den_elec  =  den_elec + denZ_ion(j)*zz_ion(j)
    enddo
    
    ! If electrons were not a separate species, add them in here
    if( .not. l_sc_elec ) then
      press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
      rho_Z    =  rho_Z    +  den_elec * xme
    endif
    
    
    ! Two terms to simplify the calculation of press_2.
    ! #assumecold
    F_px_UpS_fl = (gam_Z * beta_Z)**2 * (rho_Z*ccgs**2 + 2.5d0*press_Z)   &!&
                 +  press_Z
    F_en_UpS_fl = gam_Z**2 * u_Z * (rho_Z*ccgs**2 + 2.5d0*press_Z)
    term_1 = q_fac*F_px_UpS_fl  -  F_en_UpS_fl
    term_2 = q_fac*(gam_2*beta_2)**2  -  gam_2**2*u_2
    
    
    ! Calculate far DwS density and pressure
    rho_2   = rho_Z * gam_Z*beta_Z / gambeta_2
    press_2 = (term_1 - term_2*rho_2*ccgs**2) / (q_fac + gam_fac*term_2)
    
    
    ! Calculate Qpx & Qen (the physical escaping momentum & energy fluxes)
    Qpx  =  F_px_UpS_fl  -  (gam_2*beta_2)**2 * ( rho_2 * ccgs**2         &!&
                                                 +  gam_fac * press_2 )   &!&
                         -  press_2
    Qen  =  Qpx * q_fac
    
    
    ! Convert Qen & Qpx into q_esc_cal_**, which involves dividing by the
    !   far UpS values.  Subtract off mass-energy flux from F_en_UpS_fl to
    !   bring results in line with non-rel calculation.
    q_esc_cal_px = Qpx / F_px_UpS_fl
    q_esc_cal_en = Qen / ( F_en_UpS_fl                                    &!&
                          -  gam_Z * u_Z * rho_Z*ccgs**2 )
    
    
  !--------------------------------------------------------------------------
  !  Possibility 3: Oblique at any shock speed.  Not currently supported by
  !    code, but included here in case code is extended in future.
  !--------------------------------------------------------------------------
  else
    write(*,"(2A)") 'ERROR in q_esc_calcs: not implemented for oblique ', &!&
        'shocks yet.'
    write(*,"(2A)") "If this ever changes, don't forget to update the",   &!&
        ' shock profile in subroutine "setup_profile".'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif

return
end subroutine q_esc_calcs


!****************************************************************************
!****************************************************************************
subroutine read_col(read_unit, do_advance, ncol, outval)

! Skips the first two integer columns of a Monte Carlo output line, then
!  reads the *ncol*th column.  If *do_advance* is .false., then does not
!  advance to next line while doing so.

implicit none

  ! Input arguments
integer, intent(in) :: read_unit, ncol
logical, intent(in) :: do_advance
  ! Output argument
real(kind=8), intent(out) :: outval

  ! Local variables
integer :: idum1, idum2
real(kind=8), dimension(:), allocatable :: read_array

  
  allocate( read_array(ncol) )
  
  ! Read in the data up to the desired column...
  read(read_unit,*) idum1, idum2, read_array(:)
  ! ...and go back to the beginning of the line unless requested otherwise
  if( .not. do_advance ) backspace(read_unit)
  
  
  outval = read_array(ncol)
  
  deallocate( read_array )

return
end subroutine read_col


!****************************************************************************
!****************************************************************************
subroutine read_old_prof(n_old_skip, n_old_profs, n_old_per_prof,         &!&
     x_grid_rg, x_grid_cm, uxsk_grid, uzsk_grid, utot_grid, gam_sf_grid,  &!&
     beta_ef_grid, gam_ef_grid, btot_grid, epsB_grid, theta_grid, n_grid, &!&
     u_Z, gam_Z, rg0, r_comp, rRH, beta_Z, bmag_Z, u_2, beta_2, gam_2,    &!&
     theta_u2, bmag_2, theta_BZ, theta_B2, flux_px_UpS, flux_pz_UpS,      &!&
     flux_en_UpS)

! Reads in grid data from an existing run and uses that to construct grid
!   (e.g. position, bulk flow, magnetic field) for current run.
! In addition to explicit outputs, several downstream quantities are set
!   based on grid variables.
!
! Input file MUST be named "mc_grid_old.dat"
!
! Input arguments:
!  1) n_old_skip: number of lines (NOT profiles!) to skip in mc_grid_old.dat
!  2) n_old_profs: number of profiles to average from mc_grid_old.dat
!  3) n_old_per_prof: number of lines (including output from print_plot_vals)
!    per profile in mc_grid_old.dat
! Output arguments
!  #DOLATER

use constants
use parameters, only: na_g

implicit none

  ! Input arguments
integer, intent(in) :: n_old_skip, n_old_profs, n_old_per_prof
  ! Output arguments
integer, intent(out) :: n_grid
real(kind=8), intent(out) :: u_Z, gam_Z, bmag_Z, rg0, r_comp, rRH, beta_Z,&!&
     u_2, beta_2, gam_2, theta_u2, bmag_2, theta_BZ, theta_B2,            &!&
     flux_px_UpS, flux_pz_UpS, flux_en_UpS
real(kind=8), dimension(0:na_g), intent(out) :: x_grid_rg, x_grid_cm,     &!&
     uxsk_grid, uzsk_grid, utot_grid, gam_sf_grid, beta_ef_grid,          &!&
     gam_ef_grid, btot_grid, epsB_grid, theta_grid

  ! Local variables
integer :: i_ln, j_prof
real(kind=8) :: x_rg_loc, uxsk_norm_loc, uzsk_norm_loc, gam_sf_loc,       &!&
     bmag_loc, epsB_loc, theta_deg_loc, u_Z_km_loc, gam_Z_loc, r_comp_loc,&!&
     rRH_loc, theta_BZ_loc, theta_B2_loc, theta_u2_loc, bmag_Z_loc


!--! Open file for reading, and advance to the proper line
  open(unit=5, status='old', file='./mc_grid_old.dat')
  
  do i_ln = 1, n_old_skip
    read(5,*)
  enddo
  
  
!--! Zero out the arrays that will hold the sums/averages of the grid
   !   variables.  Several variables will not be used in averaging process,
   !   as they will be calculated based on the values of other quantities
  x_grid_rg(:)    = 0.d0
  x_grid_cm(:)    = 0.d0  ! Won't be used in averaging process
  uxsk_grid(:)    = 0.d0
  uzsk_grid(:)    = 0.d0
  utot_grid(:)    = 0.d0  ! Won't be used in averaging process
  gam_sf_grid(:)  = 0.d0
  beta_ef_grid(:) = 0.d0  ! Won't be used in averaging process
  gam_ef_grid(:)  = 0.d0  ! Won't be used in averaging process
  btot_grid(:)    = 0.d0
  epsB_grid(:)    = 0.d0
  theta_grid(:)   = 0.d0
  
  u_Z            = 0.d0
  gam_Z          = 0.d0
  bmag_Z         = 0.d0
  rg0            = 0.d0  ! Won't be used in averaging process
  r_comp         = 0.d0
  rRH            = 0.d0
  beta_Z         = 0.d0  ! Won't be used in averaging process
  u_2            = 0.d0  ! Won't be used in averaging process
  beta_2         = 0.d0  ! Won't be used in averaging process
  gam_2          = 0.d0  ! Won't be used in averaging process
  theta_u2       = 0.d0  ! Won't be used in averaging process
  bmag_2         = 0.d0  ! Won't be used in averaging process
  theta_BZ       = 0.d0  ! Won't be used in averaging process
  theta_B2       = 0.d0  ! Won't be used in averaging process
  
  
!--! Read in the data from mc_grid_old.dat
   ! WARNING: these column numbers are defined in subroutine smooth_grid_par,
   !   under glyph 'gggg'.  If those columns are changed, the column numbers
   !   here must be adjusted accordingly!
  do j_prof = 1, n_old_profs
    
  !--! Grid quantities
    do i_ln = 1, n_old_per_prof - 1  ! Last line is output of print_plot_vals
      call read_col(5, .false., 1,  x_rg_loc)
      call read_col(5, .false., 11, uxsk_norm_loc)
      call read_col(5, .false., 13, uzsk_norm_loc)
      call read_col(5, .false., 18, gam_sf_loc)
      call read_col(5, .false., 15, bmag_loc)
      call read_col(5, .false., 32, epsB_loc)
      call read_col(5, .true.,  17, theta_deg_loc)
      
      ! Add to running totals
      x_grid_rg(i_ln)   = x_grid_rg(i_ln)    +  x_rg_loc
      uxsk_grid(i_ln)   = uxsk_grid(i_ln)    +  uxsk_norm_loc
      uzsk_grid(i_ln)   = uzsk_grid(i_ln)    +  uzsk_norm_loc
      gam_sf_grid(i_ln) = gam_sf_grid(i_ln)  +  gam_sf_loc
      btot_grid(i_ln)   = btot_grid(i_ln)    +  bmag_loc
      epsB_grid(i_ln)   = epsB_grid(i_ln)    +  epsB_loc
      theta_grid(i_ln)  = theta_grid(i_ln)   +  theta_deg_loc*degtrd
    enddo
    
    
  !--! Additional quantities from print_plot_vals
    call read_col(5, .false., 1, u_Z_km_loc)
    call read_col(5, .false., 2, gam_Z_loc)
    call read_col(5, .false., 3, r_comp_loc)
    call read_col(5, .false., 4, rRH_loc)
    call read_col(5, .false., 5, theta_BZ_loc)
    call read_col(5, .false., 6, theta_B2_loc)
    call read_col(5, .false., 7, theta_u2_loc)
    call read_col(5, .true.,  8, bmag_Z_loc)
    
    ! Add to running totals
    u_Z      = u_Z       +  u_Z_km_loc*1.d5
    gam_Z    = gam_Z     +  gam_Z_loc
    r_comp   = r_comp    +  r_comp_loc
    rRH      = rRH       +  rRH_loc
    theta_BZ = theta_BZ  +  theta_BZ_loc
    theta_B2 = theta_B2  +  theta_B2_loc
    theta_u2 = theta_u2  +  theta_u2_loc
    bmag_Z   = bmag_Z    +  bmag_Z_loc
    
  enddo
  
  
  n_grid = n_old_per_prof - 1
  
  
!--! Average the values
  x_grid_rg(1:n_grid)   = x_grid_rg(1:n_grid)   / float(n_old_profs)
  uxsk_grid(1:n_grid)   = uxsk_grid(1:n_grid)   / float(n_old_profs)
  uzsk_grid(1:n_grid)   = uzsk_grid(1:n_grid)   / float(n_old_profs)
  gam_sf_grid(1:n_grid) = gam_sf_grid(1:n_grid) / float(n_old_profs)
  btot_grid(1:n_grid)   = btot_grid(1:n_grid)   / float(n_old_profs)
  epsB_grid(1:n_grid)   = epsB_grid(1:n_grid)   / float(n_old_profs)
  theta_grid(1:n_grid)  = theta_grid(1:n_grid)  / float(n_old_profs)
  
  u_Z      = u_Z      / float(n_old_profs)
  gam_Z    = gam_Z    / float(n_old_profs)
  r_comp   = r_comp   / float(n_old_profs)
  rRH      = rRH      / float(n_old_profs)
  theta_BZ = theta_BZ / float(n_old_profs)
  theta_B2 = theta_B2 / float(n_old_profs)
  theta_u2 = theta_u2 / float(n_old_profs)
  bmag_Z   = bmag_Z   / float(n_old_profs)


!--! Set the extreme boundaries of the grid
  x_grid_rg(0) = -1.0d30
  x_grid_cm(0) = -1.0d30 * rg0
  
  x_grid_rg(n_grid+1) = 1.0d30
  x_grid_cm(n_grid+1) = 1.0d30 * rg0
  
  
!--! Compute additional quantities based on values read in from input file
   !-------------------------------------------------------------------------
  ! Lorentz factor of 1.5 is approximately where precision of gamma exceeds
  !   precision of u_Z for a floating point number with three decimal places.
  !   So use that to re-calculate u_Z and the rest.
  if( gam_Z .gt. 1.5d0 ) then
    beta_Z = sqrt( 1.d0 - 1.d0/gam_Z**2 )
    u_Z    = beta_Z * ccgs
  else
    beta_Z = u_Z / ccgs
    gam_Z  = 1.d0 / sqrt( 1.d0 - beta_Z**2 )
  endif
  
  
  ! Downstream velocities
  u_2    = u_Z / r_comp
  beta_2 = u_2 / ccgs
  gam_2  = 1.d0 / sqrt( 1.d0 - beta_2**2 )
  
  
  ! Code <-> cgs conversion factor
  rg0 = (gam_Z * xmp * ccgs**2 * beta_Z) / (qcgs * bmag_Z)
  
  
  ! Downstream magnetic field
  bmag_2 = btot_grid(n_grid)
  
  
  ! Grid boundaries, cgs units
  x_grid_cm(1:n_grid) = x_grid_rg(1:n_grid) * rg0
  
  
  ! Other fluid quantities
  uxsk_grid(1:n_grid) = uxsk_grid(1:n_grid) * u_Z
  do i_ln = 1, n_grid
    if( uzsk_grid(i_ln) .gt. 1.d-90 )                                     &!&
        uzsk_grid(i_ln) = uzsk_grid(i_ln) * u_Z
  enddo
  
  utot_grid(1:n_grid) = sqrt(  uxsk_grid(1:n_grid)**2                     &!&
                             + uzsk_grid(1:n_grid)**2 )
  
  beta_ef_grid(1:n_grid) = (u_Z - uxsk_grid(1:n_grid))                    &!&
                          / ( 1.d0  -  u_Z*uxsk_grid(1:n_grid)/ccgs**2 )
  beta_ef_grid(1:n_grid) = beta_ef_grid(1:n_grid) / ccgs

  gam_ef_grid(1:n_grid)  = 1.d0 / sqrt( 1.d0  -  beta_ef_grid(1:n_grid)**2 )
   !-------------------------------------------------------------------------
   ! Other quantities determined
  
  
!--! Recalculate the far upstream fluxes in case those would be changed by
   !   the new fluid quantities
  call upstream_fluxes(flux_px_UpS, flux_pz_UpS, flux_en_UpS)

return
end subroutine read_old_prof


!****************************************************************************
!****************************************************************************
subroutine retro_time(rad_loss_fac, B_CMBz, aa, zz, gyro_denom, prp_x_cm, &!&
     ptpf, pbpf, p_perp_b_pf, gam_ptpf, acctime_sec, phi_rad, lose_pt,    &!&
     xwt, tcut_curr)

! Explicitly tracks particles downstream of the probability of return plane.
!   Does so by the "retrodictive approach" explained in Jones (1978)
!   [1978ApJ...222.1097J] and also by Ellison, Jones, & Reynolds (1990)
!   [1990ApJ...360..702E].
!
! Since we know the particle will return to the PRP, we run time "backwards"
!   by placing the particle in a bulk flow opposite to its motion.  Given a
!   presumably infinite DwS region of propagation, the particle is
!   guaranteed to return to the PRP over long enough time scales, and we
!   assume that the "true" history of the particle (its PRP1->DwS->PRP2 path)
!   resembles on average its backwards history (PRP1<-DwS<-PRP2).
!
! Input arguments:
!  1) rad_loss_fac: constant related to radiative losses
!  2) B_CMBz: effective magnetic field due to CMB at redshift of source
!  3) aa: atomic mass number of ion species
!  4) zz: charge number of ion species
!  5) prp_x_cm: starting (and ending) location for retro_time motion
!  6) xwt: current particle weight
! Output arguments:
!  1) lose_pt: whether particle is lost due to energy losses despite
!    "returning" from a probabilistic standpoint
!  2) phi_rad: phase angle of particle upon return
! Input/output arguments:
!  1) ptpf: total plasma frame momentum of particle
!  2) pbpf/p_perp_b_pf: components of ptpf parallel/perpendicular to B field
!  3) gam_ptpf: Lorentz factor associated with ptpf
!  4) gyro_denom: denominator of gyroradius fraction, zz*qcgs*bmag
!  5) acctime_sec: total accumulated acceleration time
!  6) tcut_curr: current tcut for particle tracking

use constants
use controls, only: use_custom_epsB, x_grid_stop, eta_mfp, do_rad_losses, &!&
     do_tcuts, tcuts
use grid_vars, only: n_grid, uxsk_grid, gam_sf_grid, gam_ef_grid,         &!&
     theta_grid, btot_grid
use species_vars, only: o_o_mc
use randnum

implicit none

  ! Input arguments
real(kind=8), intent(in) :: rad_loss_fac, B_CMBz, aa, zz, prp_x_cm, xwt
  ! Output arguments
logical, intent(out) :: lose_pt
real(kind=8), intent(out) :: phi_rad
  ! Input/output arguments
integer, intent(inout) :: tcut_curr
real(kind=8), intent(inout) :: ptpf, pbpf, p_perp_b_pf, gam_ptpf,         &!&
     gyro_denom, acctime_sec
  
  ! Local variables
logical :: keep_looping
real(kind=8) :: xn_per, phi_step, t_step_fac, uxsk, gam_usf, o_o_gam_u,   &!&
     gam_uef, bmag, b_cos_th, b_sin_th, B_CMB_loc, B_tot_sq, x_PT, rand,  &!&
     x_PT_old, phi_rad_old, ptpf_old, cos_old_pitch, sin_old_pitch,       &!&
     gyro_rad_cm, gyro_rad_tot_cm, t_step, x_move_bpar, dp_synch
  ! Floating point error can cause sin_del_phi to fall outside [-1,1]; place
  !   a limit on allowed values equal to the largest value of sin that can be
  !   distinguished from 1.0
real(kind=8), parameter :: sin_limit = 1.d0 - epsilon(1.d0)
  ! Local variables only used for pitch-angle diffusion
!comm real(kind=8) :: vp_tg, lambda_mfp, del_tht_max, cos_max, del_cos,   &!&
!comm      del_tht_scat, del_sin, phi_scat, cos_new_pitch, sin_new_pitch, &!&
!comm      phi_p_old, sin_del_phi, phi_p_new


!--! Set constants that will be used during the loop
  xn_per     = 10.d0
  phi_step   = twopi / xn_per
  t_step_fac = twopi * aa*xmp * ccgs * gyro_denom / xn_per  ! t_step/gam_ptpf
  
  uxsk      = -1.d0 * uxsk_grid(n_grid)
  gam_usf   = gam_sf_grid(n_grid)
  o_o_gam_u = 1.d0 / gam_sf_grid(n_grid)
  gam_uef   = gam_ef_grid(n_grid)
  if( use_custom_epsB ) then
    bmag    = btot_grid(n_grid) * (x_grid_stop / prp_x_cm)**0.25d0
  else
    bmag    = btot_grid(n_grid)
  endif
  b_cos_th  = cos(theta_grid(n_grid))
  b_sin_th  = sin(theta_grid(n_grid))
  
  B_CMB_loc = B_CMBz * gam_uef
  B_tot_sq  = bmag**2 + B_CMB_loc**2
  
  lose_pt   = .false.
  
  
!--! Initialize position, and phase angle
  x_PT = prp_x_cm
  
  call kiss(rand)
  phi_rad = rand * twopi
  
  
!--! Main loop; note similarity to main loop of code, as it's doing most of
   !   the same things with smaller scope
  keep_looping = .true.
  do while( keep_looping )
    
  !--! Store old values in preparation for the loop
    x_PT_old      = x_PT
    phi_rad_old   = phi_rad
    ptpf_old      = ptpf
    cos_old_pitch = pbpf / ptpf
    sin_old_pitch = p_perp_b_pf / ptpf
    
    
  !--! Calculate true gyroradius and total gyroradius for particle; if
     !   use_custom_epsB is true, then also need to adjust magnetic field for
     !   gyroradius and radiative cooling
    if( use_custom_epsB ) then
      bmag       = btot_grid(n_grid) * (x_grid_stop / x_PT)**0.25d0
      B_tot_sq   = bmag**2 + B_CMB_loc**2
      gyro_denom = 1.d0 / (zz*qcgs * bmag)
    endif
    gyro_rad_cm     = p_perp_b_pf * ccgs * gyro_denom
    gyro_rad_tot_cm = ptpf * ccgs * gyro_denom
    
    
  !--! Update phi_rad
    phi_rad = phi_rad_old  +  phi_step
    if( phi_rad .ge. twopi ) phi_rad = phi_rad - twopi
     
     
  !--! Calculate time step and movement distance; note that x_move_bpar = 
     !   pbpf*t_step*m_pt/gam_ptpf, but t_step = t_step_fac*gam_ptpf, so the
     !   factors of gam_ptpf divide out
    t_step      = t_step_fac * gam_ptpf
    x_move_bpar = pbpf * t_step_fac * aa*xmp
    
    
  !--! Move particle and update the acceleration time; note that we don't
     !   care about y or z motion here, and that uxsk is negative per the
     !   definition above the loop
    x_PT = x_PT_old  +  gam_usf * ( x_move_bpar * b_cos_th                &!&
                                   -  gyro_rad_cm * b_sin_th              &!&
                                     * (cos(phi_rad) - cos(phi_rad_old))  &!&
                                   + uxsk * t_step )
    acctime_sec = acctime_sec  +  t_step * gam_uef
    
    
  !--! If tcut tracking is enabled, it should continue even during retro_time
    if( do_tcuts ) then
      if( acctime_sec .ge. tcuts(tcut_curr) ) then
        call tcut_track(tcut_curr, xwt, ptpf)
        
        tcut_curr = tcut_curr + 1
      endif
    endif
    
    
  !--! Large-angle scattering.  Comment out pitch-angle diffusion
     !   section below if using LAS.
    call kiss(rand)
    phi_rad = twopi * rand   ! Completely randomize phase angle
    
    call kiss(rand)
    pbpf = 2.d0*(rand - 0.5d0) * ptpf  ! Completely randomize pbpf
    p_perp_b_pf = sqrt( ptpf**2  -  pbpf**2 )
    
    
!comm   !--! Pitch-angle diffusion.  Comment out large-angle scattering
!comm      !   section above if using PAD.
!comm     ! Compute del_tht_max
!comm     vp_tg = twopi * gyro_rad_tot_cm
!comm     if( use_custom_frg ) then
!comm       write(*,"(2A)") 'ERROR in retro_time: use of custom f(r_g) ', &!&
!comm           'not yet supported.  Add functionality or use standard.'
!comm       write(*,"(A)") 'Stopping program now.'
!comm       stop
!comm     else
!comm       lambda_mfp = eta_mfp * gyro_rad_tot_cm
!comm     endif
!comm     del_tht_max = sqrt( 6.d0 * vp_tg / (xn_per*lambda_mfp) )
!comm     cos_max = cos(del_tht_max)
!comm     
!comm     ! Compute change to pitch angle and roll
!comm     call kiss(rand)
!comm     del_cos = 1.d0  -  rand*(1.d0 - cos_max)
!comm     del_tht_scat = acos( del_cos )
!comm     del_sin = sin( del_tht_scat )
!comm     call kiss(rand)
!comm     phi_scat = rand*twopi - pii
!comm     
!comm     ! New pitch angle
!comm     cos_new_pitch = cos_old_pitch * del_cos                         &!&
!comm                    +  sin_old_pitch * del_sin * cos(phi_scat)
!comm     sin_new_pitch = sqrt( 1.d0 - cos_new_pitch**2 )
!comm     
!comm     ! Adjust components of ptpf and phase angle
!comm     pbpf        = ptpf * cos_new_pitch
!comm     p_perp_b_pf = ptpf * sin_new_pitch
!comm     
!comm     phi_p_old   = phi_rad + 0.5d0*pii
!comm     if( sin_new_pitch .ne. 0.d0 ) then
!comm       sin_del_phi = sin(phi_scat) * del_sin / sin_new_pitch
!comm       if( abs(sin_del_phi) .gt. sin_limit ) then
!comm         sin_del_phi = sign( sin_limit , sin_del_phi )
!comm       endif
!comm       phi_p_new = phi_p_old  +  asin(sin_del_phi)
!comm     else
!comm       phi_p_new = phi_p_old
!comm     endif
!comm     phi_rad = phi_p_new - 0.5d0*pii
    
    
  !--! Radiative losses
    if( do_rad_losses .and. (aa .lt. 1.d0) ) then
      
      ! Note that here dp_synch is actually dp/p.  If this value is too large
      !   we will directly integrate from p_i to get p_f, since the discrete
      !   approach would result in too high a loss in a single time step
      dp_synch = rad_loss_fac * B_tot_sq * ptpf * t_step
      
      ! Correction to make sure electrons don't lose too much energy in a
      !   single time step
      if( dp_synch .gt. 1.d-2 ) then
        ptpf = ptpf / (1.d0 + dp_synch)
      else
        ptpf = ptpf * (1.d0 - dp_synch) ! Put second factor of ptpf back
                                        !   into dp_synch
      endif
      
    endif  ! check on radiative losses
    
    ! Catch electrons that have somehow lost all their energy in a single
    !   time step, and update the pitch angle of particles that remain
    if( ptpf .le. 0.d0 ) then
      ptpf     = 1.d-99
      gam_ptpf = 1.d0
      
      lose_pt      = .true.
      keep_looping = .false.
    else
      pbpf        = ptpf * cos_old_pitch
      p_perp_b_pf = ptpf * sin_old_pitch
      gam_ptpf    = sqrt( 1.d0  +  (ptpf*o_o_mc)**2 )
    endif

    
  !--! Check for return to PRP
    if( x_PT .lt. prp_x_cm ) then
      keep_looping = .false.
    endif
    
    
  enddo  ! retro time loop
  
return
end subroutine retro_time


!****************************************************************************
!****************************************************************************
subroutine round(value, i_place)

! Takes a real value and rounds it to a specified decimal place.  I have no
!   clue why I can't find an intrinsic function to do this.

implicit none

  ! Input argument
integer, intent(in) :: i_place
  ! Input/output argument
real(kind=8), intent(inout) :: value

  ! Local variables
integer :: i_tmp

!--! Calculate the shift exponent
  i_tmp = i_place - floor( log10(abs(value)) )
  
  
!--! Shift, floor, shift back.
   ! The second argument to floor ensures that we use an integer type large
   !   enough to hold the shifted value
  value = value * 10.d0**i_tmp  +  0.5d0
  value = floor( value, selected_int_kind(15) )
  value = value * 10.d0**(-i_tmp)
  
  
return
end subroutine round


!****************************************************************************
!****************************************************************************
subroutine scattering(aa, gyro_denom, ptpf, gam_ptpf, xn_per, pbpf,       &!&
     p_perp_b_pf, phi_rad, gyro_period_sec)

! This is a combination of two subroutines from the old code: prob_scat and
!   scattering.
! Randomly moves the particle's momentum vector along the surface of the
!   unit sphere.
!
! Input arguments:
!  1) aa: particle atomic mass
!  2) gyro_denom: q*B, the denominator of gyroradius formula
!  3) ptpf: total plasma frame momentum
!  4) gam_ptpf: Lorentz factor associated with ptpf
!  5) xn_per: number of time steps a gyroperiod is divided into
! Output arguments:
!  1) gyro_period_sec
! Input/output arguments:
!  1) pbpf: component of ptpf parallel to magnetic field
!  2) p_perp_b_pf: component of ptpf perpendicular to magnetic field
!  3) phi_rad: phase angle of gyration

use constants
use controls, only: use_custom_frg, p_elec_crit, gam_elec_crit, eta_mfp
use randnum

implicit none

  ! Input arguments
real(kind=8), intent(in) :: aa, gyro_denom, ptpf, gam_ptpf, xn_per
  ! Output arguments
real(kind=8), intent(out) :: gyro_period_sec
  ! Input/output arguments
real(kind=8), intent(inout) :: pbpf, p_perp_b_pf, phi_rad
  
  ! Local variables
real(kind=8) :: gyro_rad_tot_cm, vp_tg, lambda_mfp, del_tht_max, cos_max, &!&
     cos_old_pitch, sin_old_pitch, rand, cos_del, del_tht_scat, sin_del,  &!&
     phi_scat, cos_new_pitch, sin_new_pitch, phi_p_old, sin_del_phi,      &!&
     phi_p_new
  ! Floating point error can cause sin_del_phi to fall outside [-1,1]; place
  !   a limit on allowed values equal to the largest value of sin that can be
  !   distinguished from 1.0
real(kind=8), parameter :: sin_limit = 1.d0 - epsilon(1.d0)
  
  
!--! If particle is an electron and p < p_elec_crit, use a constant MFP for
   !   scattering.
   ! Note that instead addition to calculating the gyro period in seconds,
   !   the code keeps vtpf times the gyro period to find del_tht_max.
  if( (aa .lt. 1.d0) .and. (ptpf .lt. p_elec_crit) ) then
    gyro_rad_tot_cm = p_elec_crit * ccgs * gyro_denom
    
    gyro_period_sec = twopi * gam_elec_crit * aa*xmp * ccgs * gyro_denom
    
    vp_tg = twopi * gyro_rad_tot_cm
  else
    gyro_rad_tot_cm = ptpf * ccgs * gyro_denom
    
    gyro_period_sec = twopi * gam_ptpf * aa*xmp * ccgs * gyro_denom
    
    vp_tg = twopi * gyro_rad_tot_cm
  endif
  
  
!--! To calculate collision time, need to know how MFP depends on gyro
   !   radius.  Can either use eta_mfp*r_g (the default) or a user-specified
   !   custom f(r_g).
   ! Note that instead of calculating the collision time in seconds, the code
   !   determines vtpf times the collision time, i.e. the mean free path
  if( use_custom_frg ) then
    write(*,"(2A)") 'ERROR in scattering: use of custom f(r_g) not yet',  &!&
        ' supported.  Add functionality or use standard.'
    write(*,"(A)") 'Stopping program now.'
    stop
  else
    lambda_mfp = eta_mfp * gyro_rad_tot_cm
  endif
    
  
!--! Calculate the maximum allowed change in pitch angle; this formula is
   !   slightly different from that used in previous version of the code
  del_tht_max = sqrt( 6.d0 * vp_tg / (xn_per * lambda_mfp) )
  cos_max = cos(del_tht_max)
  
  
!--! Compute the actual change in pitch angle, as well as its modulation due
   !   to a randomly-selected phase angle adjustment.  See Ellison+ (1990)
   !   [1990ApJ...360..702E] for more information.
  cos_old_pitch = pbpf / ptpf
  sin_old_pitch = p_perp_b_pf / ptpf
  
  call kiss(rand)
  cos_del = 1.d0  -  rand*(1.d0 - cos_max)  ! Change of cos btwn [0, cos_max]
  
  del_tht_scat = acos( cos_del )
  sin_del = sin( del_tht_scat )
  
  call kiss(rand)
  phi_scat = rand*twopi - pii
  
    ! Spherical law of cosines
  cos_new_pitch = cos_old_pitch * cos_del                                 &!&
                 +  sin_old_pitch * sin_del * cos(phi_scat)
  sin_new_pitch = sqrt( 1.d0 - cos_new_pitch**2 )
  
  
!--! Adjust the components of ptpf and the phase angle
  pbpf        = ptpf * cos_new_pitch
  p_perp_b_pf = ptpf * sin_new_pitch
  
  phi_p_old = phi_rad  +  0.5d0*pii
  
  if( sin_new_pitch .ne. 0.d0 ) then
  
    sin_del_phi = sin(phi_scat) * sin_del / sin_new_pitch
    
    if( abs(sin_del_phi) .gt. sin_limit ) then
      sin_del_phi = sign( sin_limit, sin_del_phi )
    endif
    
    phi_p_new = phi_p_old  +  asin(sin_del_phi)
    
  else
    
      ! sin_new_pitch = 0, so no change at all
    phi_p_new = phi_p_old
    
  endif
  
  phi_rad = phi_p_new  -  0.5d0*pii
  
return
end subroutine scattering


!****************************************************************************
!****************************************************************************
subroutine set_in_dist(l_inj_wt, n_pts_inj, i_in_distr, T_or_E, aa, den_Z,&!&
     ptot_out, wt_out, n_pts_use)

! Sets the injected particle distributions for all species.  Initially,
!   particles are placed in a Maxwell-Boltzmann (thermal) distribution based
!   on supplied temperature and particle mass.  This is corrected at the end
!   if a delta-function distribution was requested (not worried about the
!   wasted computation because this subroutine runs only rarely).
!
! Input arguments (also taken from module 'controls'):
!  1) l_inj_wt: whether each particle or each bin has equal weights (T for
!    equal weight particles, F for equal weight bins)
!  2) n_pts_inj: target # of particles for distribution
!  3) i_in_distr: thermal, delta function, or some other distribution
!  4) T_or_E: if using thermal distribution, this is temperature[K]; if
!    delta function, it's injection energy[keV]
!  5) aa: mass number for this particle species
!  6) den_Z: far UpS number density for this species
! Output arguments:
!  1) ptot_out: array holding plasma frame total momenta for all particles
!    in the distribution
!  2) xwt_out: array holding particle weights
!  3) n_pts_use: number of particles in the distribution; will almost surely
!    be different from n_pts_inj if using thermal dist
!CHECKTHIS: that output distribution matches M-B, just to make sure i
!   haven't made a typo

use constants
use parameters, only: na_p, na_i, num_therm_bins, en_rel_pt

implicit none

  ! Input arguments
integer, intent(in) :: n_pts_inj, i_in_distr
logical, intent(in) :: l_inj_wt
real(kind=8), intent(in) :: T_or_E, aa, den_Z
  ! Output arguments
integer, intent(out) :: n_pts_use
real(kind=8), dimension(na_p), intent(out) :: ptot_out, wt_out

  ! Local variables
integer :: n_pts_tot, n_per_bin, i, j, n_pts_this_bin, j_curr
real(kind=8) :: pt_mass, rm_en, kT, kT_o_rm, kT_min, kT_max, kT_rel_div,  &!&
     p_min, p_max, del_p, area_tot, p1, p2, en_o_kT1, en_o_kT2, log_f1,   &!&
     log_f2, f1, f2, area_per_pt, p_ave, bin_area, area_frac, en_inj_cgs


  ! Administrative constants, e.g. total number of particles to distribute
  !--------------------------------------------------------------------------
  if( l_inj_wt ) then
    
    ! Can't say anything about total number of particles in this case, b/c
    !   haven't split them into M-B distribution yet
    
  else
    n_per_bin = n_pts_inj / num_therm_bins   ! Integer math loses excess
                                             !   particles, but that's
                                             !   intended behavior
    n_pts_tot = n_per_bin * num_therm_bins
    
    if( n_per_bin .lt. 5 ) then
      write(*,"(2A,I0,2A)") 'ERROR in subroutine "set_in_dist": too few ',&!&
          'particles per bin (',n_per_bin,'; need at least 5).  ',        &!&
          'Increase n_pts_inj.'
      write(*,"(A)") 'Stopping program now.'
      stop
    endif
    
  endif
  !--------------------------------------------------------------------------
  ! End administrative section
  
  
  ! Set a mess of constants to be referred to routinely
  !------------------------------------------------------------------------
  pt_mass = aa * xmp
  rm_en   = pt_mass * ccgs**2
  
  kT      = xkb * T_or_E  ! Working under assumption of thermal dist now
  kT_o_rm = kT / rm_en
  kT_min  = 2.d-3 * kT   ! Minimum, maximum extent of Maxwell-Boltzmann
  kT_max  = 10d0  * kT   !   distribution
  
  kT_rel_div = en_rel_pt
  
  
  ! Find min and max momenta of M-B curve
  if( kT_o_rm .lt. kT_rel_div ) then
    
    ! In non-rel case, kinetic energy approx. thermal energy
    p_min = sqrt( 2.d0 * pt_mass * kT_min )
    p_max = sqrt( 2.d0 * pt_mass * kT_max )
    
  else
    
    ! Once particles are relatvistic, rest-mass energy becomes important:
    !   E^2 = p^2*c^2 + m^2*c^4  approx  (k*T + m*c^2)^2
    p_min = sqrt( (kT_min + rm_en)**2 - rm_en**2 )  /  ccgs
    p_max = sqrt( (kT_max + rm_en)**2 - rm_en**2 )  /  ccgs
    
  endif
  
  del_p = (p_max - p_min) / num_therm_bins
  !------------------------------------------------------------------------
  ! End of constants section
  
  
  ! Generate the Maxwell-Boltzmann distribution.  The actual calculation of
  !   f1 and f2 does not need modification for relativistic momenta.  Such
  !   a modification would only affect the normalization of the curve, not
  !   the dependence on momentum for a particular value of kT.  Since we
  !   only care about the relative fractional area of each bin, overall
  !   normalization doesn't matter.
  !------------------------------------------------------------------------
  ! Find total area under M-B curve
  area_tot = 0.d0
  do i = 1, num_therm_bins
    p1 = p_min + (i-1)*del_p
    p2 = p_min +    i *del_p
    
    if( kT_o_rm .lt. kT_rel_div ) then
      en_o_kT1 = p1**2 / (2.d0 * pt_mass * kT)
      en_o_kT2 = p2**2 / (2.d0 * pt_mass * kT)
    else
      en_o_kT1 = sqrt((p1*ccgs)**2 + rm_en**2) / kT
      en_o_kT2 = sqrt((p2*ccgs)**2 + rm_en**2) / kT
    endif
    
    log_f1 = 2.d0*log(p1) - en_o_kT1 ! Start working in log space because
    log_f2 = 2.d0*log(p2) - en_o_kT2 !   of potentially huge exponents
    f1 = exp(log_f1) ! And reverse the log transform
    f2 = exp(log_f2)
    
    ! Integrate using the trapezoid rule
    area_tot = area_tot  +  del_p * 0.5d0 * (f1 + f2)
  enddo
  
  
  ! Fill the bins with particles.  Needs to be done in two separate loops,
  !   despite the similarities, because particle weights are handled
  !   differently based on value of l_inj_wt
  ! First, particles have equal weight
  if( l_inj_wt ) then
    
    area_per_pt = area_tot / float(n_pts_inj)   ! Area each particle gets
                                                !   if l_inj_wt = T
    n_pts_tot  = 0       ! Total number of particles in M-B distribution
    
    do i = 1, num_therm_bins
      
      p1 = p_min + (i-1)*del_p
      p2 = p_min +    i *del_p
      
      p_ave = sqrt(p1*p2)   ! Geometric center of bin; particles in this bin
                            !   will receive this momentum
      
      if(kT_o_rm .lt. kT_rel_div) then
        en_o_kT1 = p1**2 / (2.d0 * pt_mass * kT)
        en_o_kT2 = p2**2 / (2.d0 * pt_mass * kT)
      else
        en_o_kT1 = sqrt((p1*ccgs)**2 + rm_en**2) / kT
        en_o_kT2 = sqrt((p2*ccgs)**2 + rm_en**2) / kT
      endif
      
      log_f1 = 2.d0*log(p1) - en_o_kT1 ! Start working in log space because of
      log_f2 = 2.d0*log(p2) - en_o_kT2 !  potentially huge exponents
      f1 = exp(log_f1) ! And reverse the log transform
      f2 = exp(log_f2)
      
      ! Calculate bin area using trapezoid rule
      bin_area = del_p * 0.5d0 * (f1 + f2)
      
      area_frac = bin_area / area_per_pt
      
      n_pts_this_bin = int( area_frac + 0.5d0 ) ! Rounded particle count
                                                !   for this bin
      
      do j = 1, n_pts_this_bin
        n_pts_tot = n_pts_tot + 1
        ptot_out(n_pts_tot) = p_ave
      enddo
      
    enddo
    
    
    ! If each particle has equal weight, then the total weight of the
    !   distribution should be proportional to the density of the species
    do j = 1, n_pts_tot
      wt_out(j) = den_Z  /  float(n_pts_tot)
    enddo
  
  
  ! Bins have equal weight
  else
  
    do i = 1, num_therm_bins
      
      p1 = p_min + (i-1)*del_p
      p2 = p_min +    i *del_p
      
      p_ave = sqrt(p1*p2)   ! Geometric center of bin; particles in this bin
                            !   will receive this momentum
      
      if(kT_o_rm .lt. kT_rel_div) then
        en_o_kT1 = p1**2 / (2.d0 * pt_mass * kT)
        en_o_kT2 = p2**2 / (2.d0 * pt_mass * kT)
      else
        en_o_kT1 = sqrt((p1*ccgs)**2 + rm_en**2) / kT
        en_o_kT2 = sqrt((p2*ccgs)**2 + rm_en**2) / kT
      endif
      
      log_f1 = 2.d0*log(p1) - en_o_kT1 ! Start working in log space because of
      log_f2 = 2.d0*log(p2) - en_o_kT2 !  potentially huge exponents
      f1 = exp(log_f1) ! And reverse the log transform
      f2 = exp(log_f2)
      
      ! Calculate bin area using trapezoid rule
      bin_area = del_p * 0.5d0 * (f1 + f2)
      
      area_frac = bin_area / area_tot
      
      
      do j = 1, n_per_bin
        j_curr = (i-1)*n_per_bin + j
        
        ptot_out(j_curr) = p_ave
        
        wt_out(j_curr) = area_frac / float(n_per_bin) * den_Z
      enddo
    
    enddo  ! loop over i
    
  endif  ! test of l_inj_wt
  
  n_pts_use = n_pts_tot
  !------------------------------------------------------------------------
  ! End of Maxwell-Boltzmann section
  
  
  ! If the particles are to use a delta-function distribution, set it here
  !--------------------------------------------------------------------------
  if( i_in_distr .eq. 2 ) then
    
    n_pts_use = n_pts_inj
    
    rm_en      = aa * xmp * ccgs**2
    en_inj_cgs = T_or_E * xkevte
    
    if( (en_inj_cgs/rm_en) .lt. en_rel_pt ) then
      p1 = sqrt( 2.d0 * aa * xmp * en_inj_cgs )
    else
      p1 = sqrt( en_inj_cgs**2 - rm_en**2 )  /  ccgs
    endif
    
    do j = 1, n_pts_inj
      ptot_out(j) = p1
      wt_out(j)   = den_Z / float(n_pts_tot)
    enddo
  
  endif  ! check of i_in_distr
  !--------------------------------------------------------------------------
  ! Delta function handled
  
  
  ! Error prevention
  if( (i_in_distr .ge. 3) .or. (i_in_distr .le. 0) ) then
    write(*,"(2A,I0,3A)") 'ERROR in subroutine "set_in_dist": i_in_distr',&!&
      ' = ', i_in_distr, '.  Code can only do i_in_distr = 1 or 2.'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif

return
end subroutine set_in_dist


!****************************************************************************
!****************************************************************************
subroutine set_psd_bins(psd_mom_min, psd_mom_max, psd_bins_per_dec_mom,   &!&
     psd_bins_per_dec_tht, psd_lin_cos_bins, psd_cos_fine, psd_tht_min,   &!&
     num_psd_mom_bins, num_psd_tht_bins, del_cos, psd_mom_bounds,         &!&
     psd_tht_bounds)

! Sets the BOUNDARIES of the bins of the phase space distribution.  The bins
!   are numbered from 0 to num_psd_***_bins, each boundary denotes the lower
!   edge of that # bin; the indices thus run from 0 to num_psd_***_bins + 1.
! WARNING: for angles, the number stored in psd_tht_bounds increases at
!   first (increasing theta), then decreases because increasing the angle
!   means decreasing the cosine.
!
! For total momentum: logarithmic spacing over all decades from Emin_keV
!  to Emax_keV
! For angle: linear spacing of cosine values for angles between psd_tht_fine
!  and pi.  Below psd_tht_fine spacing is logarithmic in theta for some # of
!  decades down to psd_tht_min
! For both: values less than the minimum are equivalent to 0.d0
!
! Input arguments:
!  1) psd_mom_min: minimum momentum[cgs] to use in PSD
!  2) psd_mom_max: maximum momentum[cgs] to use in PSD
!  3) psd_bins_per_dec_***: # of bins per decade to use in logarithmically-
!    spaced regions of PSD
!  4) psd_lin_cos_bins: # of bins to divide range [-1,cos(psd_tht_fine)] into
!  5) psd_cos_fine: cutoff between lin/cos and log/tht spacing of PSD bins
!  6) psd_tht_min: minimum angle for PSD
! Output arguments
!  1) num_psd_***_bins: total number of bins along given dimension, not
!    counting bin 0
!  2) del_cos: size of each linear cosine bin
!  3) psd_***_bounds: boundaries between bins, and upper edge of final bin


use constants
use parameters, only: psd_max

implicit none

  ! Input arguments
integer, intent(in) :: psd_bins_per_dec_mom, psd_bins_per_dec_tht,        &!&
     psd_lin_cos_bins
real(kind=8), intent(in) :: psd_mom_min, psd_mom_max, psd_cos_fine,       &!&
     psd_tht_min
  ! Output arguments
integer, intent(out) :: num_psd_mom_bins, num_psd_tht_bins
real(kind=8), intent(out) :: del_cos
real(kind=8), dimension(0:psd_max), intent(out) :: psd_mom_bounds,        &!&
     psd_tht_bounds

  ! Local variables
integer :: i, psd_log_tht_bins, j
real(kind=8) :: ten_root_mom, ten_root_tht, del_cos_tht, psd_tht_fine


!--! Zero out both arrays
  psd_mom_bounds(:) = 0.d0
  psd_tht_bounds(:) = 0.d0


!--! Calculate deltas for both log- and lin-spaced regions of PSD
  ten_root_mom = 10.d0**(1.d0/psd_bins_per_dec_mom)
  ten_root_tht = 10.d0**(1.d0/psd_bins_per_dec_tht)
  del_cos_tht  = (psd_cos_fine + 1.d0) / real(psd_lin_cos_bins)


  ! Set the momentum bins
  !--------------------------------------------------------------------------
  num_psd_mom_bins = int( log10(psd_mom_max / psd_mom_min)                &!&
                         *  psd_bins_per_dec_mom )
    ! Add two extra bins just to be safe
  num_psd_mom_bins = num_psd_mom_bins + 2
  
  if( (num_psd_mom_bins+1) .gt. psd_max ) then
    write(*,"(2A)") 'ERROR: more PSD momentum bins needed than allowed ', &!&
        'by psd_max'
    write(*,"(A,I3)") 'Bins required: ',num_psd_mom_bins
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  ! Fill in the array psd_mom_bounds, remembering that the array holds LOWER
  !   boundaries of that bin
  psd_mom_bounds(0) = -99.d0
  psd_mom_bounds(1) = log10(psd_mom_min)
  
  do i = 2, num_psd_mom_bins+1
    psd_mom_bounds(i) = psd_mom_bounds(i-1) + log10(ten_root_mom)
  enddo
  !--------------------------------------------------------------------------
  ! Momentum bins finished
  
  
  ! Set angle bins
  !--------------------------------------------------------------------------
  psd_tht_fine = acos(psd_cos_fine)
  
  psd_log_tht_bins = int( log10(psd_tht_fine/psd_tht_min)                 &!&
                         *  psd_bins_per_dec_tht )
  num_psd_tht_bins = psd_log_tht_bins  +  psd_lin_cos_bins
  
  if( (num_psd_tht_bins+1) .gt. psd_max ) then
    write(*,"(2A)") 'ERROR: more PSD anglular bins needed than allowed ', &!&
        'by psd_max'
    write(*,"(A,I3)") 'Bins required: ',num_psd_tht_bins
    write(*,"(A)") 'Stopping program now.'
    stop
  endif
  
  ! Fill the logarithmic part of psd_tht_bounds using the angle (in radians),
  !   NOT its logarithm
  psd_tht_bounds(0) = 1.d-99
  psd_tht_bounds(1) = psd_tht_min
  
  do i = 2, psd_log_tht_bins
    psd_tht_bounds(i) = psd_tht_bounds(i-1) * ten_root_tht
  enddo
  
  ! Now fill in the linear part of psd_tht_bounds.  Note that the lower
  !   boundary of the first cell is psd_cos_fine
  del_cos = (psd_cos_fine + 1.d0) / real(psd_lin_cos_bins)
  
  do i = 1, psd_lin_cos_bins+1
    j = psd_log_tht_bins + i
    psd_tht_bounds(j) = psd_cos_fine  -  del_cos*real(i-1)
  enddo
  !--------------------------------------------------------------------------
  ! Angle bins finished
  
return
end subroutine set_psd_bins


!****************************************************************************
!****************************************************************************
subroutine setup_grid(x_grid_start, x_grid_stop, n_grid, x_grid_rg,       &!&
     x_grid_cm)

! Sets the locations for all points on the grid.  In the process it sets the
!   value of n_grid.
! Historically n_grid = 99.  For plotting purposes it's best to keep it
!   there; otherwise old plotting scripts won't work with the output of the
!   code.
!
! No inputs; taken from module 'controls'
! Output arguments:
!  1) x_grid_start: UpS start of grid, cgs units
!  2) x_grid_stop: DwS end of grid, cgs units; this is adjusted based on
!    whether a PRP is in use
!  3) n_grid: total number of grid BOUNDARIES
!  4) x_grid_**: rg0/cgs locations of the grid boundaries

use parameters, only: na_g
use controls, only: x_grid_start_rg, l_use_prp, feb_DwS, x_grid_stop_rg,  &!&
     rg0

implicit none

  ! Output arguments
integer, intent(out) :: n_grid
real(kind=8), intent(out) :: x_grid_start, x_grid_stop
real(kind=8), dimension(0:na_g), intent(out) :: x_grid_rg, x_grid_cm

  ! Local variables
integer :: n_log_UpS, i, i_grid_ct, n_log_DwS
real(kind=8) :: del_log, x_end_man

!--! Recall that rg0 is the gyroradius of a proton with speed u_Z in
!--!   magnetic field bmag_Z.


!--! Set the start and stop positions in units of rg0
  x_grid_start = x_grid_start_rg * rg0
  
  if( l_use_prp ) then
    x_grid_stop = x_grid_stop_rg  * rg0
  else
    x_grid_stop    = feb_DwS
    x_grid_stop_rg = x_grid_stop / rg0
    
    write(9,"(A,ES9.3E2,2A)") 'DwS FEB set at x = ',x_grid_stop_rg,        &!&
        ' rg0.  Overwriting entered value for "XGDDW".'
    write(9,*)
  endif


!--! Logarithmically-spaced grid zones run from x_grid_start_rg to -10 rg0.
   !   Set them here.
  n_log_UpS = 27
  del_log = (log10(-x_grid_start_rg) - 1.d0) / float(n_log_UpS-1)
  
  do i = 1, n_log_UpS
    x_grid_rg(i) = -10.d0**( log10(-x_grid_start_rg) - (i-1)*del_log )
  enddo


!--! Many grid zones are set manually; zones can easily be added/removed, but
   !   make sure to change the number in the log-spaced regions UpS or DwS
  i_grid_ct = n_log_UpS
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -9.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -8.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -7.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -6.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -5.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -4.5d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -4.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -3.5d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -3.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -2.5d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -2.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.8d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.6d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.4d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.2d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.90d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.80d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.70d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.60d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.50d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.40d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.30d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.20d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.15d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -0.10d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -7.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -5.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -4.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -3.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -2.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.5d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -3.0d-03
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.0d-03
  
  
  ! Extremely fine spacing right around the shock
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.0d-04
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) = -1.0d-07
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.0d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  1.0d-07
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  1.0d-04
  
  
  ! Downstream from the shock, spacing doesn't need to be quite so fine
  !   because velocity gradients aren't as extreme, if they exist at all
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  1.0d-03
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  1.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  2.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  3.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  5.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  7.0d-02
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.10d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.15d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.20d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.25d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.30d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.40d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.50d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.60d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  0.80d0
  
  i_grid_ct = i_grid_ct + 1
  x_grid_rg(i_grid_ct) =  1.0d0


!--! As seen above, the manually-set grid zones end at x = +1 rg0.  DwS from
   !   there, more log-spaced zones.
  n_log_DwS = 16
  x_end_man = x_grid_rg(i_grid_ct)
  del_log   = (log10(x_grid_stop_rg) - log10(x_end_man)) / float(n_log_DwS)
  
  do i = 1, n_log_DwS
    x_grid_rg(i_grid_ct+i) = 10.d0**(log10(x_end_man) + i*del_log)
  enddo
  
  i_grid_ct = i_grid_ct + n_log_DwS
 

!--! Set n_grid, and convert everything from rg0 units to cgs units
  n_grid = i_grid_ct
  
  do i = 1, n_grid
    x_grid_cm(i) = x_grid_rg(i) * rg0
  enddo


!--! Set the extreme boundaries of the grid
  x_grid_rg(0) = -1.0d30
  x_grid_cm(0) = -1.0d30 * rg0
  
  x_grid_rg(n_grid+1) = 1.0d30
  x_grid_cm(n_grid+1) = 1.0d30 * rg0

return
end subroutine setup_grid


!****************************************************************************
!****************************************************************************
subroutine setup_profile(uxsk_grid, uzsk_grid, utot_grid, gam_sf_grid,    &!&
     beta_ef_grid, gam_ef_grid, btot_grid, theta_grid, epsB_grid, bmag_2)

! Sets the initial values of the shock profile
!
! No inputs; taken from modules 'controls' and 'grid_vars'
! Output arguments:
!  1) uxsk_grid: bulk fluid velocity along x axis (i.e., perpendicular to
!    shock face) in shock frame
!  2) uzsk_grid: bulk fluid velocity along z axis (i.e., parallel to shock
!    face) in shock frame
!  3) utot_grid: total bulk fluid velocity in shock frame
!  4) gam_sf_grid: bulk flow Lorentz factor in shock frame
!  5) beta_ef_grid: relative x-axis speed between plasma and explosion frames
!  6) gam_ef_grid: Lorentz factor associated with beta_ef_grid
!  7) btot_grid: total magnetic field
!  8) theta_grid: angle of magnetic field[radians] relative to shock normal
!   (i.e., to x axis)
!  9) epsB_grid: user-defined function for fraction of energy density in
!    magnetic field.  Sets value of btot_grid
!  10) bmag_2: field strength in DwS region.  Initially set in calc_DwS, it
!    may be reset here depending on values of bturb_comp_frac & bfield_amp

use constants
use parameters, only: na_g
use controls, only: u_Z, beta_Z, gam_Z, bmag_Z, theta_BZ, r_comp,         &!&
     bturb_comp_frac, bfield_amp, use_custom_epsB, n_ions, aa_ion,        &!&
     denZ_ion, l_sc_elec, zz_ion, flux_px_UpS, flux_en_UpS
use grid_vars, only: n_grid, x_grid_cm, x_grid_rg

implicit none

  ! Output arguments
real(kind=8), intent(inout) :: bmag_2
real(kind=8), dimension(0:na_g), intent(out) :: uxsk_grid, uzsk_grid,     &!&
     utot_grid, gam_sf_grid, beta_ef_grid, gam_ef_grid, btot_grid,        &!&
     theta_grid, epsB_grid

  ! Local variables
integer :: i, i_ion
real(kind=8) :: z_comp, comp_fac, btot_temp, amp_fac, den_Z, epsB0,       &!&
     den_Z_elec, sigma, rg2sd, sd2rg, epsB2, end_decay_rg, en_den


  do i = 0, n_grid
    if( x_grid_cm(i) .lt. 0.d0 ) then
      uxsk_grid(i)    = u_Z
      uzsk_grid(i)    = 0.d0
      utot_grid(i)    = u_Z
      gam_sf_grid(i)  = gam_Z
      beta_ef_grid(i) = 0.d0
      gam_ef_grid(i)  = 1.d0
      btot_grid(i)    = bmag_Z
      theta_grid(i)   = theta_BZ * degtrd
    else
      uxsk_grid(i)    = u_Z / r_comp
      uzsk_grid(i)    = 0.d0
      utot_grid(i)    = u_Z / r_comp
      gam_sf_grid(i)  = 1.d0 / sqrt( 1.d0 - (utot_grid(i)/ccgs)**2 )
      
      beta_ef_grid(i) = (beta_Z - uxsk_grid(i)/ccgs)                      &!&
                       /  ( 1.d0 - beta_Z*uxsk_grid(i)/ccgs )
      gam_ef_grid(i)  = 1.d0 / sqrt( 1.d0 - beta_ef_grid(i)**2 )
      
      ! When initializing magnetic field, include necessary corrections for
      !   turbulence compression
      z_comp       = (gam_Z * u_Z) / (gam_sf_grid(i) * uxsk_grid(i))
      comp_fac     = sqrt( third  +  2.d0*third*z_comp**2 )
      comp_fac     = 1.d0  +  (comp_fac - 1.d0)*bturb_comp_frac
      btot_temp    = bmag_Z * comp_fac
      ! Also include any additional amplification specified
      amp_fac      = 1.d0  +  (btot_temp/bmag_Z - 1.d0) * bfield_amp
      btot_grid(i) = bmag_Z  *  amp_fac
      
      theta_grid(i) = theta_BZ * degtrd
    endif
  enddo
  
!--! If directed in data_input, use a custom-defined epsilon_B to set
   !   btot_grid.
   !-------------------------------------------------------------------------
  if( use_custom_epsB ) then
    
    ! Calculate eps_B0, which depends on far UpS magnetic field and mass
    !   density.  If electrons aren't a separate species, they don't
    !   contribute enough mass to be important
    den_Z    = 0.d0
    do i_ion = 1, n_ions
      den_Z = den_Z + denZ_ion(i_ion)*aa_ion(i_ion)
    enddo
    epsB0 = bmag_Z**2 / (8.d0*pii * den_Z * rm_prot)
    
    
    ! The Monte Carlo length is rg0 = gamZ * betaZ * rm_prot / (e * BmagZ).
    !   The plasma skin depth is lamSD = gamZ * rm_prot / (4pi * e^2 * denZ),
    !   where denZ refers to the upstream number density of electrons
    ! With the definition 
    !       sigma = 2*epsB0/gamZ = BmagZ^2 / (4pi * gamZ * denZ * rm_prot),
    !   where denZ here refers to the number density of *protons*, one can
    !   show that in the shock frame (where grid exists),
    !       lamSD = betaZ/sqrt(sigma*den_p/den_e) * rg0.
    den_Z_elec = 0.d0
    if( l_sc_elec ) then
      den_Z_elec = denZ_ion(n_ions)
    else
      do i_ion = 1, n_ions
        den_Z_elec = den_Z_elec  +  denZ_ion(i_ion)*zz_ion(i_ion)
      enddo
    endif
    sigma = 2.d0 * epsB0 / gam_Z
    rg2sd = beta_Z / sqrt(sigma*den_Z/den_Z_elec)
    sd2rg = 1.d0 / rg2sd
    
    
    ! Also need the final value of epsilon_B downstream, in case our DwS
    !   region is long enough that the magnetic field can decay to this value
    ! Note that the R-H relations can be rearranged to read
    !      en_den(x)  =  F_en0/u(x) - F_px0
    !   assuming flux conservation everywhere.
    en_den = (flux_en_UpS + gam_Z*u_Z*den_Z*rm_prot) / uxsk_grid(n_grid)  &!&
            -  flux_px_UpS
    epsB2 = (Bmag_Z*comp_fac)**2 / (8.d0*pii * en_den)
    ! Use this value to compute the distance downstream at which the field
    !   will have decayed to it
    end_decay_rg = (50.d-4 / epsB2) * sd2rg
    
    
    ! Now we can calculate epsB_grid...
    ! Per the Blandford-McKee solution, energy \propto 1/chi \propto
    !   1/distance DwS.  Since we do not actually modify our pressures
    !   and densities according to the BM solution, instead modify
    !   epsilon_B
    do i = 1, n_grid
      if( x_grid_rg(i)*rg2sd .lt. -50.d0 ) then
        epsB_grid(i) = 10.4d-4 / abs(x_grid_rg(i)*rg2sd)**0.6d0
        epsB_grid(i) = max( epsB_grid(i), epsB0 )
      else if( x_grid_rg(i)*rg2sd .lt. 50.d0 ) then
        epsB_grid(i) = 1.d-4
      else if( x_grid_rg(i) .lt. end_decay_rg ) then
        epsB_grid(i) = 50.d-4 / (x_grid_rg(i)*rg2sd)
      else
        epsB_grid(i) = epsB2
      endif
    enddo
    
    
    ! ...and use it to calculate btot_grid
    do i = 1, n_grid
      en_den = (flux_en_UpS + gam_Z*u_Z*den_Z*rm_prot) / uxsk_grid(i)     &!&
              -  flux_px_UpS
      
      btot_grid(i) = sqrt(8.d0*pii * epsB_grid(i) * en_den)
    enddo
    
  else
    
    ! No custom epsilon_B requested
    epsB_grid(:) = 1.d-99
    
  endif
   !-------------------------------------------------------------------------
   ! Custom epsB_grid defined if needed
  
  
  bmag_2 = btot_grid(n_grid)

return
end subroutine setup_profile


!****************************************************************************
!****************************************************************************
subroutine smooth_grid_par(i_itr, i_shock, n_grid, x_grid_rg, x_grid_cm,  &!&
     gam_adiab_grid, uzsk_grid, theta_grid, press_psd_par, press_psd_perp,&!&
     flux_px_UpS, flux_en_UpS, gam_adiab_2, q_esc_cal_px, q_esc_cal_en,   &!&
     pxx_flux, en_flux, uxsk_grid, gam_sf_grid, btot_grid, utot_grid,     &!&
     gam_ef_grid, beta_ef_grid, epsB_grid)

! Uses tracked fluxes of particles, and the Rankine-Hugoniot jump conditions,
!   to determine the smoothed shock profile for the next iteration of the 
!   code.
! Only valid in parallel case, which makes subroutine, inputs and equations
!   simpler than would be required for oblique case.
!
! Input arguments:
!DOLATER: fix this so it's accurate
!comm !  1) i_itr: current iteration number
!comm !  2) n_grid: number of grid zones
!comm !  3) x_grid_rg: locations[rg0] of grid zone boundaries
!comm !  4) uzsk_grid: z-component of bulk fluid velocity, in shock frame.  This is
!comm !    a parallel shock, so uzsk_grid = 0 identically; just used for output
!comm !  5) theta_grid: angle[rad] between mean magnetic field and shock normal.
!comm !    This is a parallel shock, so theta_grid = 0 identically; just used for
!comm !    output
!comm !  6) pxx_flux: momentum flux of particles across grid zone boundaries
!comm !  7) en_flux: energy flux of particles across grid zone boundaries
!comm !  8) gam_adiab_2: DwS adiabatic index
!comm !  9) flux_px_UpS: far UpS momentum flux, calculated in upstream_fluxes
!comm !  10) flux_en_UpS: far UpS energy flux, calculated in upstream_fluxes
! Output arguments
!  1) utot_grid: total bulk flow speed of fluid in shock frame (is a pure
!    output b/c uxsk_grid is used in the calculations instead)
!  2) gam_ef_grid: Lorentz factor of flow relative to far UpS plasma (i.e.
!    in the explosion frame)
!  3) beta_ef_grid: bulk flow speed associated with gam_ef_grid
! Input/output arguments
!  1) uxsk_grid: x-component of bulk fluid velocity, in shock frame
!  2) gam_sf_grid: Lorentz factor associated with utot_grid, but since this
!    is a parallel shock it's the Lorentz factor associated with uxsk_grid
!  3) btot_grid: magnetic field strength[G] in each grid zone, including any
!    compression of turbulence or additional amplification

use constants
use parameters, only: na_g, beta_rel_fl
use controls, only: n_ions, denZ_ion, aa_ion, zz_ion, tZ_ion, l_sc_elec,  &!&
     tZ_elec, rg0, do_prof_fac_damp, prof_wt_fac, gam_Z, u_Z, beta_Z,     &!&
     gam_2, beta_2, u_2, do_smoothing, smooth_mom_en_fac,                 &!&
     smooth_press_flux_psd_fac, bturb_comp_frac, bfield_amp, bmag_Z,      &!&
     x_art_start_rg, use_custom_epsB

implicit none

  ! Input arguments
integer, intent(in) :: i_itr, i_shock, n_grid
real(kind=8), intent(in) :: flux_px_UpS, flux_en_UpS, gam_adiab_2,        &!&
     q_esc_cal_px, q_esc_cal_en
real(kind=8), dimension(na_g), intent(in) :: press_psd_par, press_psd_perp
real(kind=8), dimension(0:na_g), intent(in) :: x_grid_rg, x_grid_cm,      &!&
     uzsk_grid, theta_grid, epsB_grid
real(kind=8), dimension(na_g,2), intent(in) :: gam_adiab_grid
  ! Output arguments
real(kind=8), dimension(0:na_g), intent(out) :: utot_grid, gam_ef_grid,   &!&
     beta_ef_grid
  ! Input/output arguments
real(kind=8), dimension(na_g), intent(inout) :: pxx_flux, en_flux
real(kind=8), dimension(0:na_g), intent(inout) :: uxsk_grid, gam_sf_grid, &!&
     btot_grid

  ! Local variables
integer, parameter :: maxitrs = 10000
integer :: i, j, i_trans
logical :: lopen
real(kind=8), parameter :: target_err = 1.d-6
real(kind=8) :: den_Z, den_elec, press_Z, sigma_Z, ux, uz, utot, beta_ux, &!&
     beta_uz, gam_usf, gam_uef, beta_uef, theta_deg, bmag, gam_adiab_pre, &!&
     gam_adiab_post, ux_norm, uz_norm, gam_sq, beta_sq, gam_beta,         &!&
     den_ratio, press_elec, press_elec_Z, press_term_1, press_term_2, B_x,&!&
     B_z, pxx_EM, en_EM, pxx_norm_log, en_norm_log, pxz_norm_log,         &!&
     px_term_1, px_term_2, press_px, en_term_1, en_term_2, press_en,      &!&
     press_aniso, den_rat_tp, press_px_tp, press_en_tp, press_loc,        &!&
     ave_DwS_ux_px, ave_DwS_ux_en, Qpx, Qen, ux_scale_fac, den_loc,       &!&
     atan_loc, z_comp, comp_fac, btot_temp, amp_fac, ux_guess, ux_guess_p,&!&
     del_ux_guess, flux_diff, flux_diff_p, flux_diff_prime, ux_guess_next,&!&
     err_curr, ux_found, gb_guess, gb_guess_p, del_gb_guess,              &!&
     gb_guess_next, gb_found, gam_ux, en_den
real(kind=8), dimension(na_g) :: x_grid_log, x_grid_log_cm, pxx_tot,      &!&
     en_tot, pxx_norm, en_norm, pxz_tot, pxz_norm, press_tot_MC,          &!&
     ux_new_px, ux_new_en, ux_new


!--! Set constants
  ! First, density in units of proton rest mass, pressure, and number density
  den_Z    = 0.d0
  den_elec = 0.d0
  press_Z  = 0.d0
  do j = 1, n_ions
    den_Z    = den_Z    +  denZ_ion(j)*aa_ion(j)
    press_Z  = press_Z  +  denZ_ion(j)*xkb*tZ_ion(j)
    
    if( aa_ion(j) .ge. 1.d0 ) den_elec = den_elec +  denZ_ion(j)*zz_ion(j)
  enddo
  if( .not. l_sc_elec ) then
    den_Z   = den_Z    +  (xme/xmp)*den_elec
    press_Z = press_Z  +  den_elec*xkb*tZ_elec
  endif
  
  
  !DEBUGLINE (for now)
!--! Calculate the far upstream magnetization -- the ratio of the energy
   !   fluxes in EM fields and particles.  This will be used to scale the
   !   DwS decay of magnetic field, linking rg0 to the ion skin depth used in
   !   PIC sims.
   !#DOLATER: this uses gam_Z for the KE, rather than gam_Z - 1.  Okay for
   !   ultra-rel shocks, but badly mistaken in transrel limit.  Does this
   !   affect the results?
  sigma_Z = bmag_Z**2 / (4.d0*pii * gam_Z * den_Z * rm_prot)
  

!--! Determine weighting factor for profile averaging
  if( do_prof_fac_damp ) then
    if( i_itr .eq. 1) then
      prof_wt_fac = 1.00d0 * prof_wt_fac
    else if( i_itr .lt. 6 ) then
      prof_wt_fac = 1.15d0 * prof_wt_fac
    else
      prof_wt_fac = 1.50d0 * prof_wt_fac
    endif
    
    if(prof_wt_fac .gt. 10.0d0) prof_wt_fac = 10.0d0
  endif
  

!--! Compute a bunch of stuff about the current shock profile and print it
   !   to file; loop 4111 in old code
   !-------------------------------------------------------------------------
  do i = 1, n_grid
    
  !--! Grid coordinates in log space
    if(x_grid_rg(i) .lt. -1.0d0) then
      x_grid_log(i) = -log10(-x_grid_rg(i))
    else if(x_grid_rg(i) .gt. 1.d0) then
      x_grid_log(i) =  log10( x_grid_rg(i))
    else
      x_grid_log(i) = 0.0d0
    endif
    
    if(x_grid_rg(i) .lt.  0.0d0) then
      x_grid_log_cm(i) = -log10(-x_grid_rg(i) * rg0)
    else if(x_grid_rg(i) .gt. 0.d0) then
      x_grid_log_cm(i) =  log10( x_grid_rg(i) * rg0)
    else
      x_grid_log_cm(i) = 0.0d0
    endif
  
  
  !--! Pull from ***_grid arrays into easier-to-use variables
    ux        = uxsk_grid(i)
    uz        = uzsk_grid(i)
    utot      = utot_grid(i)
    beta_ux   = ux / ccgs
    beta_uz   = uz / ccgs
    gam_usf   = gam_sf_grid(i)
    gam_uef   = gam_ef_grid(i)
    beta_uef  = beta_ef_grid(i)
    theta_deg = theta_grid(i) * radtdg
    bmag      = btot_grid(i)
    
    
    ! "Pre" and "Post" refer respectively to the adiabatic index calculated
    !   before the particles have propagated through the profile for this
    !   iteration, and the adiabatic index calculated using the thermal
    !   crossing and PSD info of the most recent iteration
    gam_adiab_pre  = gam_adiab_grid(i,1)
    gam_adiab_post = gam_adiab_grid(i,2)
    
    
  !--! Basic calculations using those variables
    ux_norm   = ux / uxsk_grid(1)
    uz_norm   = 1.d-99  ! parallel shock; set to 0
    
    gam_sq    = gam_usf**2
    beta_sq   = 1.0d0 - 1.0d0/gam_sq
    gam_beta  = gam_usf * beta_ux

    den_ratio = gam_Z * beta_Z / gam_beta
    
    ! Compute pressure due to electrons if they weren't included self-
    !   consistently as a separate species.  Equation used to do so is an
    !   extension of Ellison & Moebius (1987) [1987ApJ...318..474E] Eq (2)
    if( l_sc_elec ) then
      press_elec = 0.d0
    else
      
      press_elec_Z = den_elec * xkb * tZ_elec
      
      ! This is (Gam_2 - 1)/(Gam_0 - 1).  Note assumption that Gam_0 = 5/3.
      !   #assumecold
      !CHECKTHIS: which adiabaic index (pre or post) gives sensible values?
      press_term_1 = (gam_adiab_post - 1.d0) * 1.5d0
      !  Here also: #assumecold
      press_term_2 = 1.d0  +  5.d0*third * gam_uef**2*beta_uef / gam_beta
      
      press_elec = press_elec_Z * press_term_1 * press_term_2
      
    endif
    
    
  !--! Magnetic field components, and associated fluxes according to
     !   Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D]
     !DOLATER: this assumes a mean field, not turbulence.  how do the
     !   equations change when there's turbulence?
    B_x    = bmag * cos(theta_grid(i))
    B_z    = bmag * sin(theta_grid(i))
    
    pxx_EM = gam_beta**2 / (8.d0*pii) * bmag**2                           &!&
            +  gam_sq / (8.d0*pii) * (B_z**2 - B_x**2)                    &!&
            - (gam_sq - gam_usf) / twopi * (beta_uz/beta_ux) * B_x * B_z
    
    en_EM  = gam_usf**2 / (4.d0*pii) * beta_ux * B_z**2                   &!&
            - (2.d0*gam_sq - gam_usf) / (4.d0*pii) * beta_uz * B_x * B_z
    
    
  !--! Total momentum/energy fluxes, including electrons (if needed) and EM.
     !   Also normalized against far UpS values and in log space for
     !   plotting.
    pxx_tot(i) = pxx_flux(i) + press_elec + pxx_EM
    en_tot(i)  = en_flux(i)               + en_EM                         &!&
                +  gam_adiab_post/(gam_adiab_post-1.d0) * press_elec * ux
    
    pxx_norm(i) = pxx_tot(i) / flux_px_UpS
    en_norm(i)  = en_tot(i)  / flux_en_UpS
    
    if( pxx_norm(i) .gt. 1.d-99 ) then
      pxx_norm_log = log10( pxx_norm(i) )
    else
      pxx_norm_log = -99.d0
    endif
    
    if( en_norm(i) .gt. 1.d-99 ) then
      en_norm_log = log10( en_norm(i) )
    else
      en_norm_log = -99.d0
    endif
    
    ! In a parallel shock, the z-momentum flux is irrelevant.  Set it to 0
    pxz_tot(i)   = 1.d-99
    pxz_norm(i)  = 1.d-99
    pxz_norm_log = -99.d0
    
    
  !--! Calculate pressure using the relativistic equations of Double+ (2004)
     !   [2004ApJ...600..485D], Eqs (27) and (28) specifically.  Combine the
     !   resultant pressure using smooth_mom_en_fac from input file.
     ! Note that flux_en_UpS has the rest mass-energy flux subtracted off,
     !   so add it back here
     ! Note, too, that q_esc_cal_** is already in units of far UpS flux
     ! DOLATER: per original code, "there is an unresolved question as to
     !   whether or not to use the escaping fluxes in these expressions".
     !   Using the escaping fluxes sounds reasonable, esp. in the nonrel
     !   case.  Make sure it's actually reasonable
    px_term_1 = flux_px_UpS * (1.d0 - q_esc_cal_px)                       &!&
               - gam_beta**2 * den_ratio * den_Z*rm_prot
    px_term_2 = 1.d0  +  gam_beta**2 * gam_adiab_pre / (gam_adiab_pre - 1.d0)
    press_px  = px_term_1 / px_term_2
    
    en_term_1 = flux_en_UpS * (1.d0 - q_esc_cal_en)                       &!&
               +  gam_Z*beta_Z*ccgs * den_Z*rm_prot                       &!&
               -  gam_sq * ux * den_ratio * den_Z*rm_prot
    en_term_2 = gam_sq * ux * gam_adiab_pre / (gam_adiab_pre - 1.d0)
    press_en  = en_term_1 / en_term_2
    
    ! These pressures can become negative if a sharp shock with high 
    !   compression ratio results in a great deal of escaping flux.  Place a
    !   floor on them for plotting purposes
    if( press_px .lt. 1.d-99 ) press_px = 1.d-99
    if( press_en .lt. 1.d-99 ) press_en = 1.d-99

    
  !--! Use the tabulated pressures from the thermal crossings and the PSD to
     !   determine two quantities: the total pressure and the degree of
     !   anisotropy.  Note that press_aniso will return exactly 1.0 if the
     !   pressure is isotropic, as press_par should be half of press_perp
    press_tot_MC(i) = press_psd_par(i) + press_psd_perp(i)
    press_aniso     = (2.d0*press_psd_par(i)) / press_psd_perp(i)
    
    
  !--! Calculate the expected downstream pressure in the absence of DSA, i.e.
     !   in the test particle limit.  No escaping flux to worry about here,
     !   but still need to add in the rest mass-energy flux
    if( i .eq. 1 ) then
    
      den_rat_tp = gam_Z * beta_Z / (gam_2*beta_2)
      
      px_term_1  = flux_px_UpS                                            &!&
                  -  (gam_2*beta_2)**2 * den_rat_tp * den_Z*rm_prot
      px_term_2  = 1.d0  +  (gam_2*beta_2)**2                             &!&
                           * gam_adiab_2 / (gam_adiab_2 - 1.d0)
      press_px_tp  = px_term_1 / px_term_2
      
      en_term_1  = flux_en_UpS  +  gam_Z*u_Z * den_Z*rm_prot              &!&
                  -  gam_2**2 * u_2 * den_rat_tp * den_Z*rm_prot
      en_term_2  = gam_2**2 * u_2 * gam_adiab_2/(gam_adiab_2 - 1.d0)
      press_en_tp  = en_term_1 / en_term_2
      
    endif
    
  !--! Write it all to file
    inquire(file='./mc_grid.dat',opened=lopen)
    if(.not. lopen) then
      open(unit=101,status='unknown',file='./mc_grid.dat')
    endif
    
    ! WARNING: these column numbers are reused in subroutine read_old_prof.
    !   If they are ever modified, change that subroutine accordingly!
    write(101,"(2I4,40ES12.3E2)")      &!&   !  mc_grid.dat  !gggg
                i_itr, i,              &!&
                x_grid_rg(i),          &!& ! 1
                x_grid_log(i),         &!& ! 2
                x_grid_cm(i),          &!& ! 3
                x_grid_log_cm(i),      &!& ! 4
                pxx_norm(i),           &!& ! 5
                pxx_norm_log,          &!& ! 6
                pxz_norm(i),           &!& ! 7
                pxz_norm_log,          &!& ! 8
                en_norm(i),            &!& ! 9
                en_norm_log,           &!& ! 10 
                ux_norm,               &!& ! 11 
          log10(ux_norm),              &!& ! 12 
                uz_norm,               &!& ! 13
          log10(uz_norm),              &!& ! 14 
                bmag,                  &!& ! 15
          log10(bmag),                 &!& ! 16
                theta_deg,             &!& ! 17
                gam_usf,               &!& ! 18
                1.d0/den_ratio,        &!& ! 19
                den_ratio,             &!& ! 20
          log10(press_px),             &!& ! 21
          log10(press_en),             &!& ! 22
          log10(press_psd_par(i)),     &!& ! 23
          log10(press_psd_perp(i)),    &!& ! 24
          log10(press_tot_MC(i)),      &!& ! 25
                press_aniso,           &!& ! 26
          log10(press_px_tp),          &!& ! 27
          log10(press_en_tp),          &!& ! 28
          log10(press_Z),              &!& ! 29
          log10(1.d0-q_esc_cal_px),    &!& ! 30  Remaining fluxes for plot:
          log10(1.d0-q_esc_cal_en),    &!& ! 31    momentum and energy
                epsB_grid(i),          &!& ! 32
          log10(epsB_grid(i))              ! 33
    
  enddo  ! loop over grid zones
  
  call print_plot_vals(101) 
   !-------------------------------------------------------------------------
   ! Grid output completed
  
  
!--! Return if keeping constant profile
  if( .not. do_smoothing ) return
  
  
!--! Non-rel calculation of new velocity profile
   !-------------------------------------------------------------------------
  if( beta_Z .lt. beta_rel_fl ) then
    
    ave_DwS_ux_px = 0.d0
    ave_DwS_ux_en = 0.d0
    Qpx = 0.d0  ! By default for nonrel shocks
    Qen = q_esc_cal_en * en_flux(1)
           
    
    do i = 1, n_grid
      ux        = uxsk_grid(i)
      beta_ux   = ux / ccgs
      gam_usf   = gam_sf_grid(i)
      gam_sq    = gam_usf**2
      beta_sq   = 1.0d0 - 1.0d0/gam_sq
      gam_beta  = gam_usf * beta_ux
      
      gam_adiab_post = gam_adiab_grid(i,2)
      
    !--! Magnetic field components, and associated fluxes according to
       !   Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D], and
       !   assume that uz = 0 (this *is* the parallel smoothing subroutine)
       !DOLATER: this assumes a mean field, not turbulence.  how do the
       !   equations change when there's turbulence?
      bmag   = btot_grid(i)
      B_x    = bmag * cos(theta_grid(i))
      B_z    = bmag * sin(theta_grid(i))
      
      pxx_EM = gam_beta**2 / (8.d0*pii) * bmag**2                         &!&
              +  gam_sq / (8.d0*pii) * (B_z**2 - B_x**2)
      
      en_EM  = gam_usf**2 / (4.d0*pii) * beta_ux * B_z**2
      
      
    !--! Calculate the pressure using the momentum equation only, since
       !   the energy equation can give negative fluxes if fast push is used.
       ! Determining pressure relies on near cancellation of two terms, so
       !   use a form for the non-rel equations that is expanded to beta^2
       !   to allow for better joining between rel and non-rel calculations.
       ! Do not include EM flux here, since pxx_flux tracked only particle
       !   contributions to F_px
       ! Also do not include escaping flux, since we only care about the
       !   particles that remain
      press_px = (pxx_flux(i)  -  den_Z*xmp * u_Z * ux * (1+beta_ux**2))  &!&
                / ( 1.d0  +  beta_ux**2                                   &!&
                            * gam_adiab_post/(gam_adiab_post - 1.d0) )
        
      
    !--! Combine flux-based pressure and PSD-based pressure as directed by
       !   user input
      press_loc = (1.d0 - smooth_press_flux_psd_fac) * press_px           &!&
                 +        smooth_press_flux_psd_fac  * press_tot_MC(i)
      
      
    !--! Find new velocity using newly-found pressure; here use both
       !   momentum *and* energy equations.  Need to include EM flux at this
       !   stage because flux_**_UpS included it.
       ! Because we are keeping terms out to beta^2 now, the momentum and
       !   energy equations go from linear/quadratic to cubic/quartic.
       !   Instead of solving the equations analytically, use Newton's
       !   method
    !--! Newton's method, momentum
       !---------------------------------------------------------------------
      ux_guess        = u_Z * 1.d-4
      ux_guess_p      = u_Z * 1.0001d0
      del_ux_guess    = ux_guess_p - ux_guess
      
      do j = 1, maxitrs
        
        ! Calculate flux differences associated with both ux_guess and
        !   ux_guess_p
        beta_ux   = ux_guess / ccgs
        px_term_1 = den_Z*xmp * u_Z * ux_guess * (1.d0 + beta_ux**2)
        px_term_2 = (1.d0 +  beta_ux**2 * gam_adiab_post                  &!&
                            / (gam_adiab_post - 1.d0) )  *  press_loc
        
        flux_diff = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2
        
        beta_ux   = ux_guess_p / ccgs
        px_term_1 = den_Z*xmp * u_Z * ux_guess_p * (1.d0 + beta_ux**2)
        px_term_2 = (1.d0 +  beta_ux**2 * gam_adiab_post                  &!&
                            / (gam_adiab_post - 1.d0) )  *  press_loc
        
        flux_diff_p = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2
        
        
        ! Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
        flux_diff_prime = (flux_diff_p - flux_diff) / del_ux_guess
        
        
        ! Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
        ux_guess_next = ux_guess  -  flux_diff/flux_diff_prime
        
        
        ! Relative change in this step
        err_curr = (ux_guess_next - ux_guess) / ux_guess
        
        ! If the relative change is small enough, we've found our solution
        !   and can exit the loop; otherwise return for another cycle
        if( abs(err_curr) .lt. target_err ) then
          ux_found = ux_guess_next
          exit
        endif
        
        ux_guess     = (ux_guess_next + ux_guess) * 0.5d0
        ux_guess_p   = 1.001d0 * ux_guess
        del_ux_guess = ux_guess_p - ux_guess
        
      enddo
      
      ! Did we hit the maximum number of iterations without finding the
      !   flux-conserving solution?
      if( j .ge. maxitrs ) then
        write(*,"(2A,I0)") "WARNING(smooth_grid_par): no Newton's method",&!&
            ' solution for px flux, grid zone ',i
        ux_found = ux_new_px(i-1)
      endif
      
      ux_new_px(i) = ux_found
      
      
    !--! Newton's method, energy
       !---------------------------------------------------------------------
      ux_guess        = u_Z * 1.d-4
      ux_guess_p      = u_Z * 1.0001d0
      del_ux_guess    = ux_guess_p - ux_guess
      
      do j = 1, maxitrs
        
        ! Calculate flux differences associated with both ux_guess and
        !   ux_guess_p
        beta_ux   = ux_guess / ccgs
        en_term_1 = 0.5d0 * den_Z*xmp * u_Z * ux_guess**2                 &!&
                   * (1.d0 + 1.25d0*beta_ux**2)
        en_term_2 = gam_adiab_post / (gam_adiab_post - 1.d0) * press_loc  &!&
                   * ux_guess * (1.d0 + beta_ux**2)
        
        flux_diff = flux_en_UpS - Qen - en_EM - en_term_1 - en_term_2
        
        beta_ux   = ux_guess_p / ccgs
        en_term_1 = 0.5d0 * den_Z*xmp * u_Z * ux_guess_p**2               &!&
                   * (1.d0 + 1.25d0*beta_ux**2)
        en_term_2 = gam_adiab_post / (gam_adiab_post - 1.d0) * press_loc  &!&
                   * ux_guess_p * (1.d0 + beta_ux**2)
        
        flux_diff_p = flux_en_UpS - Qen - en_EM - en_term_1 - en_term_2
        
        
        ! Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
        flux_diff_prime = (flux_diff_p - flux_diff) / del_ux_guess
        
        
        ! Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
        ux_guess_next = ux_guess  -  flux_diff/flux_diff_prime
        
        
        ! Relative change in this step
        err_curr = (ux_guess_next - ux_guess) / ux_guess
        
        ! If the relative change is small enough, we've found our solution
        !   and can exit the loop; otherwise return for another cycle
        if( abs(err_curr) .lt. target_err ) then
          ux_found = ux_guess_next
          exit
        endif
        
        ux_guess     = (ux_guess_next + ux_guess) * 0.5d0
        ux_guess_p   = 1.001d0 * ux_guess
        del_ux_guess = ux_guess_p - ux_guess
        
      enddo
      
      ! Did we hit the maximum number of iterations without finding the
      !   flux-conserving solution?
      if( j .ge. maxitrs ) then
        write(*,"(2A,I0)") "WARNING(smooth_grid_par): no Newton's method",&!&
            ' solution for en flux, grid zone ',i
        ux_found = ux_new_en(i-1)
      endif
      
      ux_new_en(i) = ux_found
      
      
    !--! Find the average downstream velocity, which will be used for scaling
       !   the profiles.  Only average the last 10 grid positions
      if( i .gt. (n_grid-10) ) then
        ave_DwS_ux_px = ave_DwS_ux_px  +  ux_new_px(i)
        ave_DwS_ux_en = ave_DwS_ux_en  +  ux_new_en(i)
      endif
    
    enddo  ! loop over grid positions
      
      
  !--! Scale the velocity profile, smooth it, and average the momentum and
     !   energy curves as directed by user input
    ave_DwS_ux_px = 0.1d0 * ave_DwS_ux_px
    ave_DwS_ux_en = 0.1d0 * ave_DwS_ux_en
    
    ux_scale_fac = (u_Z - u_2)/(ux_new_px(1) - ave_DwS_ux_px)  
    do i = 1, n_grid
      ux_new_px(i) = ux_scale_fac * (ux_new_px(i) - ave_DwS_ux_px)  +  u_2
      if( x_grid_rg(i) .ge. 0.d0 ) ux_new_px(i) = u_2
    enddo
    
    ux_scale_fac = (u_Z - u_2)/(ux_new_en(1) - ave_DwS_ux_en)  
    do i = 1, n_grid
      ux_new_en(i) = ux_scale_fac * (ux_new_en(i) - ave_DwS_ux_en)  +  u_2
      if( x_grid_rg(i) .ge. 0.d0 ) ux_new_en(i) = u_2
    enddo
    
    call smooth_profile(n_grid, ux_new_px)
    call smooth_profile(n_grid, ux_new_en)
    
    do i = 1, n_grid
      ux_new(i) = (1.d0 - smooth_mom_en_fac) * ux_new_px(i)               &!&
                 +        smooth_mom_en_fac  * ux_new_en(i)
    enddo
    
  endif  ! check on shock speed
   !-------------------------------------------------------------------------
   ! Non-rel shock smoothing complete
  
  
!--! Rel calculation of new velocity profile
   ! CHECKTHIS: what happens if gam*u is used as a smoothing variable instead
   !   of just u?  At high speeds gam*u is much more variable than just u
   !-------------------------------------------------------------------------
  if( beta_Z .ge. beta_rel_fl ) then
    
    ave_DwS_ux_px = 0.d0
    ave_DwS_ux_en = 0.d0
    Qpx = q_esc_cal_px * pxx_flux(1)
    Qen = q_esc_cal_en * en_flux(1)
   
    do i = 1, n_grid
      ux        = uxsk_grid(i)
      beta_ux   = ux / ccgs
      gam_usf   = gam_sf_grid(i)
      gam_sq    = gam_usf**2
      beta_sq   = 1.0d0 - 1.0d0/gam_sq
      gam_beta  = gam_usf * beta_ux
      
      den_loc = gam_Z * beta_Z / (gam_usf*beta_ux) * den_Z
      
      gam_adiab_post = gam_adiab_grid(i,2)
      
    !--! Magnetic field components, and associated fluxes according to
       !   Eqs. (27) & (28) of Double+ (2004) [2004ApJ...600..485D], and
       !   assume that uz = 0 (this *is* the parallel smoothing subroutine)
       !DOLATER: this assumes a mean field, not turbulence.  how do the
       !   equations change when there's turbulence?
      bmag   = btot_grid(i)
      B_x    = bmag * cos(theta_grid(i))
      B_z    = bmag * sin(theta_grid(i))
      
      pxx_EM = gam_beta**2 / (8.d0*pii) * bmag**2                         &!&
              +  gam_sq / (8.d0*pii) * (B_z**2 - B_x**2)
      
      en_EM  = gam_usf**2 / (4.d0*pii) * beta_ux * B_z**2
      
      
    !--! Calculate the pressure using the momentum equation only, since
       !   the energy equation can give negative fluxes if fast push is used
       ! Do not include EM flux here, since pxx_flux tracked only particle
       !   contributions to F_px
       ! Also do not include escaping flux, since we only care about the
       !   particles that remain
      px_term_1 = pxx_flux(i)  -  (gam_beta)**2 * den_loc*rm_prot
      px_term_2 = 1.d0  +  (gam_beta)**2                                  &!&
                          * gam_adiab_post / (gam_adiab_post - 1.d0)
      press_px  = px_term_1 / px_term_2
      
      
    !--! Combine flux-based pressure and PSD-based pressure as directed by
       !   user input
      press_loc = (1.d0 - smooth_press_flux_psd_fac) * press_px           &!&
                 +        smooth_press_flux_psd_fac  * press_tot_MC(i)
      
      
    !--! Find new velocity using just the momentum flux equation.  Energy
       !   flux will be used further down.
       ! DOLATER: neither of these includes EM component of flux, which will
       !   be important once turbulence is considered
    !--! Newton's method, momentum.  Use gamma*beta since that retains the
       !   sign information of just ux, but also scales to relativistic
       !   speeds
       !---------------------------------------------------------------------
      gb_guess        = gam_Z * beta_Z * 1.d-4
      gb_guess_p      = gb_guess * 1.0001d0
      del_gb_guess    = gb_guess_p - gb_guess
      
      do j = 1, maxitrs
        
        ! Calculate flux differences associated with both gb_guess and
        !   gb_guess_p
        px_term_1 = gam_Z * beta_Z * den_Z / den_loc  * gb_guess          &!& 
                   * ( den_loc*rm_prot  +  press_loc * gam_adiab_post     &!&
                                          / (gam_adiab_post - 1.d0) )
        px_term_2 = press_loc
        flux_diff = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2
        
        px_term_1 = gam_Z * beta_Z * den_Z / den_loc  * gb_guess_p        &!& 
                   * ( den_loc*rm_prot  +  press_loc * gam_adiab_post     &!&
                                          / (gam_adiab_post - 1.d0) )
        px_term_2 = press_loc
        flux_diff_p = flux_px_UpS - Qpx - pxx_EM - px_term_1 - px_term_2
        
        
        ! Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
        flux_diff_prime = (flux_diff_p - flux_diff) / del_gb_guess
        
        
        ! Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
        gb_guess_next = gb_guess  -  flux_diff/flux_diff_prime
        
        
        ! Relative change in this step
        err_curr = (gb_guess_next - gb_guess) / gb_guess
        
        ! If the relative change is small enough, we've found our solution
        !   and can exit the loop; otherwise return for another cycle
        if( abs(err_curr) .lt. target_err ) then
          gb_found = gb_guess_next
          exit
        endif
        
        gb_guess     = (gb_guess_next + gb_guess) * 0.5d0
        gb_guess_p   = 1.001d0 * gb_guess
        del_gb_guess = gb_guess_p - gb_guess
        
      enddo
      
      ! Did we hit the maximum number of iterations without finding the
      !   flux-conserving solution?
      if( j .ge. maxitrs ) then
        write(*,"(2A,I0)") "WARNING(smooth_grid_par): no Newton's method",&!&
            ' solution for px flux, grid zone ',i
        gb_found = ux_new_px(i) / sqrt( ccgs**2 - ux_new_px(i-1)**2 )
      endif
      
      ux_new_px(i) = gb_found / sqrt( 1.d0 + gb_found**2 ) * ccgs
      
      
    !--! Newton's method, energy
       !---------------------------------------------------------------------
      gb_guess        = gam_Z * beta_Z * 1.d-4
      gb_guess_p      = gb_guess * 1.0001d0
      del_gb_guess    = gb_guess_p - gb_guess
      
      do j = 1, maxitrs
        
        ! Calculate flux differences associated with both gb_guess and
        !   gb_guess_p
        gam_ux    = sqrt( 1.d0 + gb_guess**2 )
        en_term_1 = gb_guess * gam_ux * ccgs                              &!&
                   * ( den_loc*rm_prot  +  gam_adiab_post                 &!&
                                         / (gam_adiab_post - 1.d0)        &!&
                                         * press_loc )
        flux_diff = flux_en_UpS - Qen - en_EM - en_term_1
        
        gam_ux    = sqrt( 1.d0 + gb_guess_p**2 )
        en_term_1 = gb_guess_p * gam_ux * ccgs                            &!&
                   * ( den_loc*rm_prot  +  gam_adiab_post                 &!&
                                         / (gam_adiab_post - 1.d0)        &!&
                                         * press_loc )
        flux_diff_p = flux_en_UpS - Qen - en_EM - en_term_1
        
        
        ! Calculate derivative: f'(x_n)  =  (f(x_n + dx) - f(x_n)) / dx
        flux_diff_prime = (flux_diff_p - flux_diff) / del_gb_guess
        
        
        ! Actual Newton's method step: x_n+1  =  x_n  -  f(x_n)/f'(x_n)
        gb_guess_next = gb_guess  -  flux_diff/flux_diff_prime
        
        
        ! Relative change in this step
        err_curr = (gb_guess_next - gb_guess) / gb_guess
        
        ! If the relative change is small enough, we've found our solution
        !   and can exit the loop; otherwise return for another cycle
        if( abs(err_curr) .lt. target_err ) then
          gb_found = gb_guess_next
          exit
        endif
        
        gb_guess     = (gb_guess_next + gb_guess) * 0.5d0
        gb_guess_p   = 1.001d0 * gb_guess
        del_gb_guess = gb_guess_p - gb_guess
        
      enddo
      
      ! Did we hit the maximum number of iterations without finding the
      !   flux-conserving solution?
      if( j .ge. maxitrs ) then
        write(*,"(2A,I0)") "WARNING(smooth_grid_par): no Newton's method",&!&
            ' solution for en flux, grid zone ',i
        gb_found = ux_new_en(i) / sqrt( ccgs**2 - ux_new_en(i-1)**2 )
      endif
      
      ux_new_en(i) = gb_found / sqrt( 1.d0 + gb_found**2 ) * ccgs
      
      
    !--! Find the average downstream velocity, which will be used for scaling
       !   the profiles.  Only average the last 10 grid positions
      if( i .gt. (n_grid-10) ) then
        ave_DwS_ux_px = ave_DwS_ux_px  +  ux_new_px(i)
        ave_DwS_ux_en = ave_DwS_ux_en  +  ux_new_en(i)
      endif
    
    enddo  ! loop over grid positions
    
    
  !--! Smooth the velocity profile, rescale it, and average the momentum and
     !   energy curves as directed by user input
    call smooth_profile(n_grid, ux_new_px)
    call smooth_profile(n_grid, ux_new_en)
    
    ave_DwS_ux_px = 0.1d0 * ave_DwS_ux_px
    ave_DwS_ux_en = 0.1d0 * ave_DwS_ux_en
    
    ux_scale_fac = (u_Z - u_2)/(ux_new_px(1) - ave_DwS_ux_px)  
    do i = 1, n_grid
      ux_new_px(i) = ux_scale_fac * (ux_new_px(i) - ave_DwS_ux_px)  +  u_2
      if( x_grid_rg(i) .ge. 0.d0 ) ux_new_px(i) = u_2
    enddo
    
    ux_scale_fac = (u_Z - u_2)/(ux_new_en(1) - ave_DwS_ux_en)  
    do i = 1, n_grid
      ux_new_en(i) = ux_scale_fac * (ux_new_en(i) - ave_DwS_ux_en)  +  u_2
      if( x_grid_rg(i) .ge. 0.d0 ) ux_new_en(i) = u_2
    enddo
    
    do i = 1, n_grid
      ux_new(i) = (1.d0 - smooth_mom_en_fac) * ux_new_px(i)               &!&
                 +        smooth_mom_en_fac  * ux_new_en(i)
    enddo
    
  endif  ! check on shock speed
   !-------------------------------------------------------------------------
   ! Rel shock smoothing complete
  
  
!--! Artificial smoothing if directed
  if( x_art_start_rg .lt. 0.d0 ) then
    do i = 1, n_grid
      if( x_grid_rg(i) .gt. x_art_start_rg ) then
        i_trans = i-1
        exit
      endif
    enddo
    
    ux_scale_fac = -(ux_new(i_trans) - ux_new(n_grid))                    &!&
                  /  atan(x_grid_rg(i_trans))
    
    do i = i_trans, i_shock
      atan_loc  = -atan(x_grid_rg(i))
      ux_new(i) = atan_loc * ux_scale_fac  +  ux_new(n_grid)
    enddo
    
  endif
  
  
!--! Average with previous profile
  !CHECKTHIS: what happens if gam * beta is used as an averaging variable
  !  instead of just u for a rel shock?
  do i = 1, n_grid
    ux_new(i) = (ux_new(i)  +  prof_wt_fac * uxsk_grid(i))                &!&
               /  ( 1.d0 + prof_wt_fac )
  enddo
  
  
!--! Compute output arrays based on new profile
  do i = 1, n_grid
    uxsk_grid(i)    = ux_new(i)
    gam_sf_grid(i)  = 1.d0 / sqrt( 1.d0 - (uxsk_grid(i)/ccgs)**2 )
    utot_grid(i)    = ux_new(i)
    beta_ef_grid(i) = (u_Z - uxsk_grid(i)) / ccgs                         &!&
                     /  (1.d0 - u_Z*uxsk_grid(i)/ccgs**2 )
    gam_ef_grid(i)  = 1.d0 / sqrt( 1.d0 - beta_ef_grid(i)**2 )
    
    ! Include necessary corrections for turbulence compression
    z_comp       = (gam_Z * u_Z) / (gam_sf_grid(i) * uxsk_grid(i))
    comp_fac     = sqrt( third  +  2.d0*third*z_comp**2 )
    comp_fac     = 1.d0  +  (comp_fac - 1.d0)*bturb_comp_frac
    btot_temp    = bmag_Z * comp_fac
    ! Also include any additional amplification specified
    amp_fac      = 1.d0  +  (btot_temp/bmag_Z - 1.d0) * bfield_amp
    btot_grid(i) = bmag_Z  *  amp_fac
    
    ! If a custom epsilon_B is in place, use that to calculate the magnetic
    !   field, not the preceding code.
    ! Note that the R-H relations can be rearranged to read
    !      en_den(x)  =  F_en0/u(x) - F_px0
    !   assuming flux conservation everywhere.
    if( use_custom_epsB ) then
      en_den = (flux_en_UpS + gam_Z*u_Z*den_Z*rm_prot) / uxsk_grid(i)     &!&
              -  flux_px_UpS
      
      btot_grid(i) = sqrt(8.d0*pii * epsB_grid(i) * en_den)
    endif

  enddo
  
return
end subroutine smooth_grid_par


!****************************************************************************
!****************************************************************************
subroutine smooth_profile(n_grid, y_prof)

! Takes an input velocity profile and performs two tasks: enforces
!   monotonicity by smoothing out dips/bumps, and averages nearby points to
!   smooth sharp edges
!
! Input arguments:
!  1) n_grid: number of grid zones in position and velocity arrays
! Input/output arguments:
!  1) y_prof: array holding velocity profile

use constants
use parameters, only: na_g

implicit none

  ! Input arguments
integer, intent(in) :: n_grid
  ! Input/output arguments
real(kind=8), dimension(na_g), intent(inout) :: y_prof

  ! Local variables
integer :: i
real(kind=8), dimension(na_g) :: y_prof_2

!--! Run from downstream to upstream and eliminate dips/bumps
  do i = n_grid, 2, -1
    if(y_prof(i-1) .lt. y_prof(i)) y_prof(i-1) = y_prof(i) 
  enddo
  
  
!--! Now smooth sharp edges by averaging adjacent grid locations; handle the
   !   edges separately for different weighting
  y_prof_2(2)  =  0.25d0 * (2.d0*y_prof(1) + y_prof(2) + y_prof(3))
  do i = 3, n_grid-2
     y_prof_2(i)  =  third * (y_prof(i-1) + y_prof(i) + y_prof(i+1))
  enddo
  y_prof_2(n_grid-1) = 0.25d0 * ( y_prof(n_grid-2) + y_prof(n_grid-1)     &!&
                                 +  2.d0*y_prof(n_grid) )
  
  do i = 2, n_grid-1
    y_prof(i) = y_prof_2(i)
  enddo

return
end subroutine smooth_profile


!****************************************************************************
!****************************************************************************
subroutine tcut_track(tcut_curr, xwt, ptpf)

! If directed in data_input, track particles' momenta -- and whether they are
!   still interacting with the shock -- at various times since acceleration
!   began.
!
! Input arguments:
!  1) tcut_curr: current tcut for tracking
!  2) xwt: particle weight
!  3) ptpf: plasma-frame total momentum
! Output arguments:
!  None.  All adjustments made to arrays in modules

use iteration_vars, only: wt_coupled, spec_coupled
use species_vars, only: i_ion

implicit none

  ! Input arguments
integer, intent(in) :: tcut_curr
real(kind=8), intent(in) :: xwt, ptpf
  ! Output arguments
  
  ! Local variables
integer :: i_pt, j_th

!--! Since particle is still coupled to shock (i.e. being accelerated), add
   !   its weight to appropriate bin of wt_coupled
!$omp atomic
  wt_coupled(tcut_curr,i_ion) = wt_coupled(tcut_curr,i_ion) + xwt
  
  
!--! For spectra, need to convert ptpf into a psd bin, then add it to array;
   !   note that we don't care about angular component
  call get_psd_bins(0.d0, ptpf, .false., i_pt, j_th)
!$omp atomic
  spec_coupled(i_pt, tcut_curr, i_ion) =                                  &!&
                                 spec_coupled(i_pt, tcut_curr, i_ion)  +  xwt

return
end subroutine tcut_track


!****************************************************************************
!****************************************************************************
subroutine tcut_print()

! If tcuts were tracked during the run, print out the particle counts and
!   spectra at each tcut
!
! Input arguments:
!  None.  Just printing arrays from modules
! Output arguments:
!  None.

use constants
use controls, only: n_tcuts, tcuts, n_ions, n_itrs
use psd_vars, only: num_psd_mom_bins, psd_mom_bounds
use iteration_vars, only: i_itr, wt_coupled, spec_coupled

implicit none

  ! Input arguments
  ! Output arguments
  
  ! Local variables
integer :: i_ion, i_cut, i_pt, j_plot


!--! Set a floor on the array values, and normalize spec_coupled so that
   !   each spectrum has a total weight of 1.0 (actual weight, of course, is
   !   the entry in wt_coupled)
  do i_ion = 1, n_ions
    do i_cut = 1, n_tcuts
      if( wt_coupled(i_cut,i_ion) .lt. 1.d-60 ) then
        wt_coupled(i_cut,i_ion) = 1.d-99
      endif
      
      if( sum(spec_coupled(:,i_cut,i_ion)) .gt. 1.d-99 )                  &!&
        spec_coupled(:,i_cut,i_ion) = spec_coupled(:,i_cut,i_ion)         &!&
                                     /  sum(spec_coupled(:,i_cut,i_ion))
      
      do i_pt = 0, num_psd_mom_bins
        if( spec_coupled(i_pt,i_cut,i_ion) .lt. 1.d-60 ) then
          spec_coupled(i_pt,i_cut,i_ion) = 1.d-99
        endif
      enddo
    enddo
  enddo


!--! If the desired output files haven't been opened, do that now
  if( i_itr .eq. 1 ) then
    open(unit=670, file='mc_coupled_wts.dat')
    open(unit=671, file='mc_coupled_spectra.dat')
  endif
  
  
!--! Write out the particle counts at each tcut, as a histogram
  j_plot = 0
  
  do i_cut = 1, n_tcuts
    j_plot = j_plot + 1
    
    write(670,"(2i4,10es12.3e2)") i_itr, j_plot,    &!&
                log10(tcuts(i_cut)),                &!& ! 1    tcut times
                log10(wt_coupled(i_cut,1:n_ions))       ! 2-?  wts by species
    
    j_plot = j_plot + 1
    if( i_cut .lt. n_tcuts ) then
      write(670,"(2i4,10es12.3e2)") i_itr, j_plot,  &!&
                  log10(tcuts(i_cut+1)),            &!& ! 1    tcut times
                  log10(wt_coupled(i_cut,1:n_ions))     ! 2-?  wts by species
    endif
  enddo
  call print_plot_vals(670)
  
  
!--! Write out the spectra at each tcut, as a histogram.  Split into groups
   !   by species, and give each spectrum (i.e. tcut) its own column
  do i_ion = 1, n_ions
    j_plot = 0
    
    do i_pt = 0, num_psd_mom_bins
      j_plot = j_plot + 1
      
      write(671,"(2i4,f5.1,100es12.3e2)") i_itr, j_plot,  &!&
          float(i_ion),                                   &!& ! 1
                psd_mom_bounds(i_pt),                     &!& ! 2   cgs units
                psd_mom_bounds(i_pt)-log10(xmp*ccgs),     &!& ! 3   nat units
          log10(spec_coupled(i_pt,1:n_tcuts,i_ion))           ! 4-? spectra
      
      j_plot = j_plot + 1
      if( i_pt .lt. num_psd_mom_bins ) then
        write(671,"(2i4,f5.1,100es12.3e2)") i_itr, j_plot,&!&
            float(i_ion),                                 &!& ! 1
                  psd_mom_bounds(i_pt+1),                 &!& ! 2   cgs units
                  psd_mom_bounds(i_pt+1)-log10(xmp*ccgs), &!& ! 3   nat units
            log10(spec_coupled(i_pt,1:n_tcuts,i_ion))         ! 4-? spectra
      endif
    enddo
    
    call print_plot_vals(671)
  enddo
  
  
!--! If this is the last iteration, close the files
  if( i_itr .eq. n_itrs ) then
    close(670)
    close(671)
  endif

return
end subroutine tcut_print


!****************************************************************************
!****************************************************************************
subroutine thermo_calcs(num_crossings, n_cr_count, therm_grid, therm_pxsk,&!&
     therm_ptsk, therm_xwt, nc_unit, do_therm_pt_formatted, psd, zone_pop,&!&
     press_psd_par, press_psd_perp, en_den_psd)

! Reads in the PSD -- as well as the thermal arrays & scratch file -- to
!   calculate the pressure (which may be anisotropic) everywhere on the grid. 
!
! Input arguments:
!  1) num_crossings: array containing number of thermal particle crossings
!    recorded at each grid zone
!  2) n_cr_count: total number (up to na_cr) of thermal particle crossings;
!    signals whether scratch file was used
!  3) therm_grid: array of grid zone number for thermal crossings
!  4) therm_pxsk: array of x-momentum for thermal crossings
!  5) therm_ptsk: array of total momentum for thermal crossings
!  6) therm_xwt: array of flux-adjusted particle weights for thermal particle
!    crossings
!  7) nc_unit: unit number for the scratch file holding crossing data
!  8) do_therm_pt_formatted: logical telling whether the scratch file is
!    human-readable or not
!  9) psd: 3-D phase space distribution holding shock-frame ptot and
!    cos(theta) for all recorded CR grid crossings
!  10) zone_pop: (Lorentz-invariant) number of particles in each grid zone
! Output arguments
!  1) press_psd_par: component of pressure parallel to magnetic field
!  2) press_psd_perp: component of pressure perpendicular to magnetic field
!  3) en_den_psd: local kinetic energy density of fluid in every grid zone

use constants
use parameters, only: na_g, psd_max, na_cr, en_rel_pt
use controls, only: aa_ion, denZ_ion, psd_lin_cos_bins, gam_Z, beta_Z,    &!&
     tZ_ion, zz_ion, tZ_elec, l_sc_elec
use psd_vars, only: num_psd_tht_bins, psd_tht_bounds, num_psd_mom_bins,   &!&
     psd_mom_bounds
use grid_vars, only: n_grid, uxsk_grid, gam_sf_grid
use species_vars, only: i_ion, o_o_mc

implicit none

  ! Input arguments
integer, intent(in) :: n_cr_count, nc_unit
integer, dimension(na_g), intent(in) :: num_crossings
integer, dimension(na_cr), intent(in) :: therm_grid
logical, intent(in) :: do_therm_pt_formatted
real(kind=8), dimension(na_g), intent(in) :: zone_pop
real(kind=8), dimension(na_cr), intent(in) :: therm_pxsk, therm_ptsk,     &!&
     therm_xwt
real(kind=8), dimension(0:psd_max,0:psd_max,na_g), intent(in) :: psd
  ! Output arguments
real(kind=8), dimension(na_g), intent(out) :: press_psd_par,              &!&
     press_psd_perp, en_den_psd

  ! Local variables
integer :: max_cross, ntot_crossings, i, i_grid, idum, n, k_pt, j_th,     &!&
     kpt_Xf, jth_Xf
integer, dimension(na_g) :: n_cross_fill
real(kind=8) :: rest_mass_en, pxsk, ptsk, cell_wt, norm_fac, cos_hi,      &!&
     cos_lo, gam_u, beta_u, cos_theta_sk, etot_sk, px_Xf, pt_Xf, den_loc, &!&
     press_loc, den_elec, press_fac, press_par_tmp, press_perp_tmp,       &!&
     gam_tmp, en_den_fac
real(kind=8), dimension(0:psd_max) :: cos_center, pt_cgs_center, vel_ptot,&!&
     d2N_pop
real(kind=8), dimension(0:psd_max, 0:psd_max, na_g) :: d2N_pf
real(kind=8), dimension(:,:), allocatable :: therm_px, therm_pt, therm_wt


!--! Initialize d2N_pf to "zero"
   ! Dimensions of arrays (chosen for maximum speed during loops):
   !    1 - cos(theta)
   !    2 - ptot
   !    3 - grid position
   ! The omission of "dpdcos" is not a typo; the array will only ever hold a
   !   total count of particles, not a true dN/dp
  d2N_pf(:,:,:) = 1.d-99
  
  
!--! Set constants to be used repeatedly, including calculating the center
   !   points of all bins to save time later
  rest_mass_en = aa_ion(i_ion) * rm_prot  ! Rest mass-energy of the
                                          !  current particle species
   
  do j_th = 0, num_psd_tht_bins
    ! Determine current cosines, remembering that psd_tht_bounds has both a
    !  linearly-spaced region in cosine and logarithmically-spaced region in
    !  theta.
    if( j_th .gt. (num_psd_tht_bins - psd_lin_cos_bins) ) then
      cos_hi = psd_tht_bounds(j_th)
      cos_lo = psd_tht_bounds(j_th + 1)
    else if( j_th .eq. (num_psd_tht_bins - psd_lin_cos_bins) ) then
      cos_hi = cos(psd_tht_bounds(j_th))
      cos_lo = psd_tht_bounds(j_th + 1)
    else
      cos_hi = cos(psd_tht_bounds(j_th))
      cos_lo = cos(psd_tht_bounds(j_th + 1))
    endif
    
    ! Minus sign needed because finest gradations must point UpS
    cos_center(j_th) = -0.5d0 * ( cos_lo + cos_hi )
  enddo
  
  do k_pt = 0, num_psd_mom_bins
    pt_cgs_center(k_pt) = 0.5d0 * ( psd_mom_bounds(k_pt)                  &!&
                                   +  psd_mom_bounds(k_pt+1) )
    ! Convert from log to linear space
    pt_cgs_center(k_pt) = 10.d0**pt_cgs_center(k_pt)
  enddo
  
  
! Read the scratch files associated with thermal particles.  Bin them into
!  the combined d2N/(dp-dcos) for the shock frame
!----------------------------------------------------------------------------
!--! Set up the arrays that will hold the crossing data
  max_cross = maxval(num_crossings,1)
  allocate( therm_px(max_cross, 0:n_grid+1) )
  allocate( therm_pt(max_cross, 0:n_grid+1) )
  allocate( therm_wt(max_cross, 0:n_grid+1) )


!--! Fill d2N_pf with data from thermal particles
   !-------------------------------------------------------------------------
  ntot_crossings = sum(num_crossings,1)
  n_cross_fill(:) = 0 ! n_cross_fill needs to be an array because we will be
                      !  skipping around the grid as we move through the data
                      !  rather than filling a single zone at a time
  
  rewind(nc_unit)
  
    
  ! Handle crossings stored within the crossing arrays
  do i = 1, n_cr_count
    i_grid = therm_grid(i)
    
    n_cross_fill(i_grid) = n_cross_fill(i_grid) + 1
  
    ! Note: ordering of coordinates chosen to make memory accesses in next
    !  loop faster
    therm_px(n_cross_fill(i_grid), i_grid) = therm_pxsk(i)
    therm_pt(n_cross_fill(i_grid), i_grid) = therm_ptsk(i)
    therm_wt(n_cross_fill(i_grid), i_grid) = therm_xwt(i)
  enddo
    
  if( ntot_crossings .gt. na_cr ) then
    
    ! Need to go into the scratch file for the remainder of the crossings
    do i = 1, ntot_crossings-na_cr
      if( do_therm_pt_formatted ) then
        read(nc_unit,"(I3,I5,4ES20.12E2)") i_grid, idum, pxsk, ptsk, cell_wt
      else
        read(nc_unit) i_grid, idum, pxsk, ptsk, cell_wt
      endif
      
      n_cross_fill(i_grid) = n_cross_fill(i_grid) + 1
      
      ! Note: ordering of coordinates chosen to make memory accesses in next
      !  loop faster
      therm_px(n_cross_fill(i_grid), i_grid) = pxsk
      therm_pt(n_cross_fill(i_grid), i_grid) = ptsk
      therm_wt(n_cross_fill(i_grid), i_grid) = cell_wt
    enddo
    
  endif
   !-------------------------------------------------------------------------
   ! Arrays filled
  

!--! Loop over grid locations, filling in the array d2N_pf
   !-------------------------------------------------------------------------
  do i = 1, na_g
    if( num_crossings(i) .eq. 0 ) then
      cycle     ! Ignore zones that no thermal particles crossed
    endif
    
    ! Get Lorentz factor and speed relating the shock and plasma frames
    gam_u  = gam_sf_grid(i)
    beta_u = uxsk_grid(i) / ccgs
      
    
    ! Loop over all crossings for this grid zone, binning the thermal
    !  particles.
    do n = 1, num_crossings(i)
        
        ! Transform the center of the zone into the new frame
        ptsk         = therm_pt(n,i)
        pxsk         = therm_px(n,i)
        etot_sk      = sqrt( (ptsk*ccgs)**2  +  rest_mass_en**2 )
        
        px_Xf    = gam_u * (pxsk - beta_u*etot_sk/ccgs)
        pt_Xf    = sqrt( (ptsk**2 - pxsk**2) + px_Xf**2 )
          ! In rare cases, floating point roundoff can cause pt_Xf to be
          !   smaller than abs(px_Xf)
        if( abs(px_Xf) .gt. pt_Xf ) then
          px_Xf = sign(pt_Xf, px_Xf)
        endif
        
        
        ! Get location of center in transformed d2N_pf
        call get_psd_bins(px_Xf, pt_Xf, .true., kpt_Xf, jth_Xf)
        
        
        ! And add particle count to the appropriate bin in d2N_pf
        d2N_pf(jth_Xf, kpt_Xf, i) = d2N_pf(jth_Xf, kpt_Xf, i)             &!&
                                   +  therm_wt(n,i)
    enddo
    
  enddo
   !-------------------------------------------------------------------------
   ! Plasma frame d2N filled from thermal arrays
  
  
!--! Deallocate the crossing data arrays
  deallocate( therm_px ) ;  deallocate( therm_pt ) ;  deallocate( therm_wt )
!----------------------------------------------------------------------------
! Thermal particles binned


! Now, get the CRs into d2N_pf.  Do this by transforming the bin *centers* of
!   PSD into each frame, using just the center location to rebin.
! DOLATER: consider doing an exact calculation of bin overlaps rather than
!   just the center-point rebinning that currently happens
!----------------------------------------------------------------------------
!--! Loop over grid locations, rebinning all cells with particles
  do i = 1, n_grid
    
  !--! Get Lorentz factor and speed relating the shock and current frames
    gam_u  = gam_sf_grid(i)
    beta_u = uxsk_grid(i) / ccgs
      
    
  !--! Loop over cells in d2N_sf
    do j_th = 0, num_psd_mom_bins
      do k_pt = 0, num_psd_tht_bins
        
        ! Skip empty zones
        if( psd(k_pt,j_th,i) .le. 1.d-66 ) then
          cycle
        else
          cell_wt = psd(k_pt,j_th,i)
        endif
        
        
        ! Transform the center of the zone into the new frame
        cos_theta_sk = cos_center(j_th)
        ptsk         = pt_cgs_center(k_pt)
        pxsk         = ptsk * cos_theta_sk
        etot_sk      = sqrt( (ptsk*ccgs)**2  +  rest_mass_en**2 )
        
        px_Xf    = gam_u * (pxsk - beta_u*etot_sk/ccgs)
        pt_Xf    = sqrt( ptsk**2 - pxsk**2 + px_Xf**2 )
        
        
        ! Get location of center in transformed d2N_pf
        call get_psd_bins(px_Xf, pt_Xf, .true., kpt_Xf, jth_Xf)
        
        
        ! And add particle count to the appropriate bin in d2N_pf
        d2N_pf(jth_Xf, kpt_Xf, i) = d2N_pf(jth_Xf, kpt_Xf, i)  +  cell_wt
        
      enddo ! loop over angles
    enddo ! loop over ptot
    
    
    norm_fac = sum( d2N_pf(:,:,i), MASK = (d2N_pf(:,:,i).gt.1.d-66) )
    if( (num_crossings(i) .eq. 0) .and. (norm_fac .gt. 0.d0) ) then
      norm_fac = norm_fac  +  denZ_ion(i_ion) / uxsk_grid(i)
    endif
    if( norm_fac .gt. 0.d0 ) norm_fac = zone_pop(i) / norm_fac
    
    
    do k_pt = 0, psd_max-1
      do j_th = 0, psd_max-1
        if( d2N_pf(j_th, k_pt, i) .gt. 1.d-66 ) then
          d2N_pf(j_th, k_pt, i) = d2N_pf(j_th, k_pt, i) * norm_fac
        endif
      enddo
    enddo
    
    d2N_pop(i) = sum( d2N_pf(:,:,i), MASK = (d2N_pf(:,:,i).gt.1.d-66) )
  enddo ! loop over grid locations
!----------------------------------------------------------------------------
! Transformed d2N_pf calculated
  
  
! Now sum over d2N_pf to get the pressure at all grid locations.  With
!   angular information we can get both parallel and perpendicular components
!----------------------------------------------------------------------------
!--! Find the velocity associated with each momentum bin
  do k_pt = 0, num_psd_mom_bins
    vel_ptot(k_pt) = pt_cgs_center(k_pt) * ccgs * o_o_mc                  &!&
                    /  sqrt( 1.d0  +  (pt_cgs_center(k_pt)*o_o_mc)**2 )
  enddo
  
  
!--! Loop over grid zones and find pressure components everywhere
   !-------------------------------------------------------------------------
  do i = 1, n_grid
    
    den_loc = gam_Z*beta_Z*denZ_ion(i_ion)                                &!&
             /  sqrt( gam_sf_grid(i)**2 - 1.d0 )
    
  !--! Determine the normalization factors to use during the pressure
     !   calculation
     ! Three cases: (1) no detected particles, (2) only CRs detected, and
     !   (3) thermal particles detected
     !-----------------------------------------------------------------------
    if( (maxval( d2N_pf(:,:,i) ) .lt. 1.d-66) .and.                       &!&
        (num_crossings(i) .eq. 0) ) then
      
    !--! Case (1): No particles of any kind detected at this grid location.
       !   Thermal particles must have passed through, so find their density
       !   and analytically determine components of pressure
       ! #assumecold: using gam_adiab = 5/3 in the pressure calc
      press_loc = den_loc**(5.d0*third) * xkb*tZ_ion(i_ion)
      
      if( aa_ion(i_ion) .ge. 1.d0 ) then
          den_elec = den_loc * zz_ion(i_ion)
      endif
      
      ! Add pressure due to electrons if they aren't a separate species
      if( .not. l_sc_elec ) then
        press_loc = press_loc  +  den_elec * xkb*tZ_elec
      endif
      
      ! Update components of pressure, assuming isotropy; twice as many
      !   perpendicular components as parallel (y & z vs x), thus the factor
      !   of 2 in the computation
      press_psd_par(i)  = press_psd_par(i)  +        third * press_loc
      press_psd_perp(i) = press_psd_perp(i) + 2.d0 * third * press_loc
      
      ! With only thermal particles at this grid location, the energy density
      !   is easy to determine
      ! #assumecold
      en_den_psd(i) = en_den_psd(i)  +  1.5d0 * press_loc
      
      ! Don't bother running through the d2N calculation
      cycle
      
      
    else if( num_crossings(i) .eq. 0 ) then
      
    !--! Case (2): no thermal particles, but some CRs detected.  Find
       !   pressure due to un-tracked thermal particles
       ! #assumecold: using gam_adiab = 5/3 in the pressure calc
      press_loc = den_loc**(5.d0*third) * xkb*tZ_ion(i_ion)
      
      if( aa_ion(i_ion) .gt. 1.d0 ) then
        den_elec = den_loc * zz_ion(i_ion)
      endif
      
      ! Add pressure due to electrons if they aren't a separate species
      if( .not. l_sc_elec ) then
        press_loc = press_loc  +  den_elec * xkb*tZ_elec
      endif
      
      ! Scale the contributions from the thermal particles by the fraction
      !   of the zone's population they represent, and add that to the
      !   running total
      press_loc = press_loc  *  (zone_pop(i)-d2N_pop(i)) / zone_pop(i)
      press_psd_par(i)  = press_psd_par(i)  +        third * press_loc
      press_psd_perp(i) = press_psd_perp(i) + 2.d0 * third * press_loc
      
      ! Set the normalization factor for the CRs that were tracked during
      !   the iteration; summation over CRs will total d2N_pop(i), which
      !   makes up for the amount subtracted from press_loc above
      norm_fac = den_loc / zone_pop(i)
      
      ! Finally, calculate contribution of thermal particles to energy
      !   density.  This will be updated during the following loop to include
      !   CR contribution.
      ! #assumecold
      en_den_psd(i) = en_den_psd(i)  +  1.5d0 * press_loc
      
    else
      
    !--! Case (3): thermal particles detected.  If no CRs are present that
       !   just means they didn't propagate UpS to this zone.  Whether CRs
       !   are present or not, d2N_pf represents a true count of the total
       !   population of this grid zone.  Set normalization factor
       !   accordingly
       norm_fac = den_loc / zone_pop(i)
    
    endif
     !-----------------------------------------------------------------------
     ! Appropriate normalization constant(s) found
    
    
  !--! Find the pressure using the density calculated above
    do k_pt = 0, num_psd_mom_bins
      
      ! These factors are the same in all cells of this ptot column
      press_fac  = third * pt_cgs_center(k_pt) * vel_ptot(k_pt) * norm_fac
      gam_tmp    = sqrt( 1.d0  +  (pt_cgs_center(k_pt)*o_o_mc)**2 )
      en_den_fac = (gam_tmp - 1.d0) * rest_mass_en
      
      do j_th = 0, num_psd_tht_bins
        
        ! Don't bother with empty cells
        if( d2N_pf(j_th,k_pt,i) .lt. 1.d-66 ) cycle
        
        press_par_tmp  = d2N_pf(j_th, k_pt, i) * press_fac                &!&
                        *  cos_center(j_th)**2
          ! The factor of two below is because there are two perpendicular
          !   components of pressure
        press_perp_tmp = d2N_pf(j_th, k_pt, i) * press_fac                &!&
                        *  (1.d0 - cos_center(j_th)**2)
        
        press_psd_par(i)  = press_psd_par(i)  + press_par_tmp
        press_psd_perp(i) = press_psd_perp(i) + press_perp_tmp
        
        
        ! Update the energy density also
        en_den_psd(i) = en_den_psd(i)  +  en_den_fac                      &!&
                                         * d2N_pf(j_th,k_pt,i) * norm_fac
      enddo
    enddo
    
  enddo  ! loop over i
   !-------------------------------------------------------------------------
   ! Loop over grid zones finished
  
return
end subroutine thermo_calcs


!****************************************************************************
!****************************************************************************
subroutine upstream_fluxes(flux_px_UpS, flux_pz_UpS, flux_en_UpS)

! Calculates the far upstream fluxes for the shock.
!
! Two different cases considered:
!  (1) Non-relativistic oblique shock.  Uses equations of Ellison+ (1996)
!    [1996ApJ...473.1029E]
!  (2) Relativistic shock, any obliquity.  Uses equations of Double+ (2004)
!    [2004ApJ...600..485D]
! Only oblique equations used because they reduce trivially to parallel cases
!   when theta_BZ = 0
! HOWEVER, assumes that z-component of far UpS velocity is 0 in all cases;
!   in practice oblique shocks would induce some z-velocity in the shock
!   profile even though particles initially arrive with no bulk z component.
! Also assumes isotropic initial pressure, so no off-diagonal components in
!   pressure tensor.
!
! No inputs; pulls everything from module 'controls'
! Output arguments:
!  1) flux_px_UpS: far UpS momentum flux, x component
!  2) flux_pz_UpS: far UpS momentum flux, z component
!  3) flux_en_UpS: far UpS energy flux

use constants
use parameters, only: beta_rel_fl
use controls, only: l_oblique, n_ions, denZ_ion, tZ_ion, aa_ion, zz_ion,  &!&
     l_sc_elec, tZ_elec, bmag_Z, theta_BZ, gam_Z, beta_Z, u_Z

implicit none

  ! Output arguments
real(kind=8), intent(out) :: flux_px_UpS, flux_pz_UpS, flux_en_UpS
  
  ! Local variables
integer :: i
logical :: l_nonrel
real(kind=8) :: press_Z, rho_Z, den_elec, gam_sph, e_Z, B_sq, B_x, B_z,   &!&
     F_px_fl, F_px_EM, F_pz_fl, F_pz_EM, F_en_fl, F_en_EM


!--! Determine which variant of the flux equations we should use
  if( beta_Z .lt. beta_rel_fl ) then
    l_nonrel = .true.
  else
    l_nonrel = .false.
  endif
  if( l_oblique ) then
    write(*,"(2A)") 'ERROR in "upstream_fluxes": cannot handle oblique ', &!&
        'shocks yet'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif


!--! UpS internal energy density and pressure, assuming isotropic particle
   !   distribution.  Note that this INCLUDES the mass-energy density, which
   !   is typically omitted in non-rel calculations
  press_Z  = 0.d0
  rho_Z    = 0.d0
  den_elec = 0.d0
  do i = 1, n_ions
    press_Z  =  press_Z  +  denZ_ion(i) * xkb*tZ_ion(i)
    rho_Z    =  rho_Z    +  denZ_ion(i) * xmp*aa_ion(i)
    
    if( aa_ion(i) .ge. 1.d0 ) den_elec  =  den_elec + denZ_ion(i)*zz_ion(i)
  enddo
  
  ! If electrons were not a separate species, add them in here
  if( .not. l_sc_elec ) then
    press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
    rho_Z    =  rho_Z    +  den_elec * xme
  endif
    
  ! Assume an adiabatic index of 5/3, appropriate for non-rel ideal gas,
  !   to calculate the far UpS internal energy
  ! #assumecold
  gam_sph = 5.d0 / 3.d0
  e_Z     = rho_Z * ccgs**2  +  1.d0/(gam_sph - 1.d0) * press_Z


!--! Quantities related to the UpS magnetic field.  Note that B_z is the
   !   z-component of the magnetic field, not B_0
  B_sq = bmag_Z**2
  B_x  = bmag_Z * cos(theta_BZ*degtrd)
  B_z  = bmag_Z * sin(theta_BZ*degtrd)


!--! Momentum, x-component
  ! Fluid part (Double+ Eq 23)
  F_px_fl = (gam_Z * beta_Z)**2 * (e_Z + press_Z)  +  press_Z

  ! EM part (Double+ Eq 25)
  F_px_EM = (gam_Z * beta_Z)**2 * B_sq/(8.d0*pii)                         &!&
           +  gam_Z**2/(8.d0*pii) * (B_z**2 - B_x**2)
  
  ! Total
  flux_px_UpS = F_px_fl  +  F_px_EM


!--! Momentum, z-component
  ! Fluid part (Double+ Eq 24)
  F_pz_fl = 0.d0

  ! EM part (Double+ Eq 26)
  F_pz_EM = -gam_Z/(4.d0*pii) * B_x * B_z
  
  ! Total
  flux_pz_UpS = F_pz_fl  +  F_pz_EM


!--! Energy
  ! Fluid part (Double+ Eq 20)
  F_en_fl = gam_Z**2 * beta_Z * (e_Z + press_Z)

  ! EM part (Double+ Eq 21)
  F_en_EM = gam_Z**2 * beta_Z * B_z**2/(4.d0*pii)
  
  ! Total -- convert to cgs units!
  flux_en_UpS = ccgs * (F_en_fl + F_en_EM)
  
  ! And subtract off the mass-energy flux to bring it in line with non-rel
  !   calculations and what the MC code actually tracks
  flux_en_UpS = flux_en_UpS  -  gam_Z * u_Z * rho_Z*ccgs**2


!--! Non-relativistic version.  Note that it's missing the rho*c^2 flux
   !   present in the relativistic forms above.  It is also expanded to
   !   second order in beta_Z (only in the hydro terms, for now) to allow
   !   for more precise matching with the relativistic version
  if( l_nonrel ) then
    flux_px_UpS = rho_Z * u_Z**2 * (1.d0 + beta_Z**2)                     &!&
                 +  press_Z * (1.d0 + gam_sph/(gam_sph-1.d0)*beta_Z**2)   &!&
                 +  B_z**2/(8.d0*pii)
    flux_pz_UpS = -B_x * B_z / (4.d0*pii)
    flux_en_UpS = 0.5d0 * rho_Z * u_Z**3 * (1.d0 + 1.25d0 * beta_Z**2)    &!&
                 +  press_Z * u_Z * gam_sph/(gam_sph-1.d0)                &!&
                   * (1.d0 + beta_Z**2)                                   &!&
                 +  u_Z * B_z**2/(4.d0*pii)
  endif

return
end subroutine upstream_fluxes


!****************************************************************************
!****************************************************************************
subroutine upstream_machs(mach_sonic, mach_alfven)

! Calculates the sonic and Alfven mach numbers for the shock.
!
! For speed of sound, uses Equation (13) of Fujimura & Kennel (1979)
!   [1979A%26A....79..299F]
! For Alfven wave speed, uses Equation (46) of Gedalin (1993)
!   [1993PhRvE..47.4354G]
!
! No input arguments; these come from module "controls"
! Output arguments:
!  1) mach_sonic: sonic Mach number of UpS flow
!  2) mach_alfven: Alfvenic Mach number of UpS flow

use constants
use parameters, only: beta_rel_fl
use controls, only: beta_Z, n_ions, denZ_ion, tZ_ion, aa_ion, zz_ion,     &!&
     l_sc_elec, tZ_elec, bmag_Z

implicit none

  ! Input arguments
  ! Output arguments
real(kind=8), intent(out) :: mach_sonic, mach_alfven

  ! Local variables
integer :: i
logical :: l_nonrel
real(kind=8) :: gam_adiab, press_Z, rho_Z, den_elec, R_fac, c_s, v_A,     &!&
     a_fac, enthalpy


!--! Determine whether we need to use more involved formulae
  if( beta_Z .lt. beta_rel_fl ) then
    l_nonrel = .true.
  else
    l_nonrel = .false.
  endif


!--! Assume cold UpS plasma, so that the adiabatic index is 5/3 identically
   ! #assumecold
  gam_adiab = 5.d0 / 3.d0


!--! Find FK1979's R factor, the ratio of pressure to rest mass energy
   !   density
  press_Z  = 0.d0
  rho_Z    = 0.d0
  den_elec = 0.d0
  do i = 1, n_ions
    press_Z  =  press_Z  +  denZ_ion(i) * xkb*tZ_ion(i)
    rho_Z    =  rho_Z    +  denZ_ion(i) * xmp*aa_ion(i)
    
    if( aa_ion(i) .ge. 1.d0 ) den_elec  =  den_elec + denZ_ion(i)*zz_ion(i)
  enddo
  
  ! If electrons were not a separate species, add them in here
  if( .not. l_sc_elec ) then
    press_Z  =  press_Z  +  den_elec * xkb*tZ_elec
    rho_Z    =  rho_Z    +  den_elec * xme
  endif
  
  R_fac = press_Z / (rho_Z * ccgs**2)


!--! Calculate the sound speeds differently based on whether the shock is
   !   non-rel or rel
  if( l_nonrel ) then
    
    ! Use standard forms for both speeds
    c_s = sqrt( gam_adiab * press_Z / rho_Z )
    v_A = bmag_Z / sqrt( 4.d0 * pii * rho_Z )
    
  else
    
    ! Plug everything into FK1979's Equation (13)
    a_fac = gam_adiab / (gam_adiab - 1.d0)
    c_s   = ccgs * sqrt( gam_adiab * R_fac / (a_fac*R_fac + 1.d0) )
    
    ! And into Gedalin (1993)'s Equation (46); note assumption that
    !   equation of state is    e  =  rho c^2  +  P/(Gam-1)
    enthalpy = a_fac * press_Z  +  rho_Z * ccgs**2
    v_A = ccgs  /  sqrt( 1.d0  +  4.d0*pii * enthalpy / bmag_Z**2 )
    
  endif


!--! Now calculate the Mach numbers using the derived wave speeds
  mach_sonic  = beta_Z * ccgs / c_s
  mach_alfven = beta_Z * ccgs / v_A

return
end subroutine upstream_machs


!****************************************************************************
!****************************************************************************
subroutine Xform_p_PS(aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad, uxsk,     &!&
     uzsk, utot, gam_usf, b_cos_th, b_sin_th, ptsk, pxsk, pzsk, gam_ptsk)

! Takes particle's Plasma frame momentum components and transforms them to
!   the Shock frame.
! Transformation formulas are correct for all obliquities.
!
! Notes from Glen Double (31 Oct 2002):
! Method: define flow velocity and particle momentum in component vectors
!   wrt the plasma xyz frame.  Find p_para (wrt u) component of particle
!   momentum using scalar product of p and u.  Then p_perp = p - p_para
!   (using vectors).  Next, make relativistic transformation on p_para while
!   p_perp remains constant.  Finally, recreate p_rel = p_rel_para + p_perp,
!   and find p_rel components in the original xyz frame by taking scalar
!   products along the xyz axes.
!
! Input arguments:
!  1) aa: particle atomic mass
!  2) pbpf: component of ptpf parallel to magnetic field
!  3) p_perp_b_pf: component of ptpf perpendicular to magnetic field
!  4) gam_ptpf: Lorentz factor associated with ptpf
!  5) phi_rad: phase angle of gyration; looking UpS, counts clockwise from
!    +z axis
!  6) uxsk: bulk flow speed along x axis
!  7) uzsk: bulk flow speed along z axis
!  8) utot: total bulk flow speed
!  9) gam_usf: Lorentz factor associated with utot
!  10) b_cos_th: component of magnetic field along x axis
!  11) b_sin_th: component of magnetic field along z axis
! Output arguments:
!  1) ptsk: total shock frame momentum in new grid zone
!  2) pxsk/pzsk: x & z components of ptsk
!  3) gam_ptsk: Lorentz factor associated with ptsk
! Input/output arguments:

use constants
use controls, only: l_oblique
use species_vars, only: o_o_mc

implicit none

  ! Input arguments
real(kind=8), intent(in) :: aa, pbpf, p_perp_b_pf, gam_ptpf, phi_rad,     &!&
     uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th
  ! Output arguments
real(kind=8), intent(out) :: ptsk, pxsk, gam_ptsk
  
  ! Local variables
real(kind=8) :: phi_p, p_p_cos, pxpf, pypf, pzpf, o_o_utot, pysk, pzsk,   &!&
     pbsk, p_perp_b_sk


  phi_p = phi_rad + 0.5d0*pii
  
  p_p_cos = p_perp_b_pf * cos(phi_p)
  
  
!--! xyz plasma frame components
  !CHECKTHIS: theta in shock frame may be different from theta in plasma
  !  frame because of lorentz transformation between the two.  Are the next
  !  lines correct in light of this?
  pxpf = pbpf*b_cos_th  -  p_p_cos*b_sin_th
  pypf = p_perp_b_pf * sin(phi_p)
  pzpf = pbpf*b_sin_th  +  p_p_cos*b_cos_th
  
  
!--! xyz shock frame components
  if( l_oblique ) then
    o_o_utot = 1.d0 / utot
    
    pxsk = ( (gam_usf-1.d0)*(uxsk*o_o_utot)**2 + 1.d0 ) * pxpf              &!&
          + (gam_usf-1.d0)*(uxsk*uzsk*o_o_utot**2) * pzpf                   &!&
          +  gam_usf * gam_ptpf * aa*xmp * uxsk
    
    pysk = pypf
    
    pzsk = (gam_usf-1.d0)*(uxsk*uzsk*o_o_utot**2) * pxpf                    &!&
          + ( (gam_usf-1.d0)*(uzsk*o_o_utot)**2 + 1.d0 ) * pzpf             &!&
          + gam_usf * gam_ptpf * aa*xmp * uzsk
  else
    pxsk =  gam_usf * ( pxpf +  gam_ptpf * aa*xmp * uxsk )
    
    pysk = pypf
    
    pzsk = pzpf
  endif
  
  
!--! Parallel/perpendicular (new) shock frame components
  ptsk = sqrt( pxsk**2  +  pysk**2  +  pzsk**2 )
  
  pbsk = pxsk*b_cos_th  +  pzsk*b_sin_th
  
  if( ptsk .lt. abs(pbsk) ) then
    p_perp_b_sk = 1.d-6 * ptsk
    pbsk        = sign( sqrt( ptsk**2 - p_perp_b_sk**2 ), pbsk )
    !CHECKTHIS: does this *ever* happen?!
    write(*,"(A)") 'Warning: ptsk < pbsk in Xform_p_PS'
  else
    p_perp_b_sk = sqrt( ptsk**2  -  pbsk**2 )
  endif
  
  gam_ptsk = sqrt( (ptsk*o_o_mc)**2  +  1.d0 )
  
return
end subroutine Xform_p_PS


!****************************************************************************
!****************************************************************************
subroutine Xform_p_PSP(aa, ptpf, pbpf, p_perp_b_pf, gam_ptpf, phi_rad,    &!&
     uxsk_old, uzsk_old, utot_old, gam_usf_old, b_cos_old, b_sin_old,     &!&
     uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th, ptsk, pxsk, pysk,     &!&
     pzsk, pbsk, p_perp_b_sk, gam_ptsk)

! Takes particle's old plasma frame momentum components and transforms them
!   twice: from old Plasma frame to new Shock frame, then from new Shock
!   frame to new Plasma frame.
! Only ever called if particle scattered between zones with different bulk
!   flow velocities.
! Transformation formulas are correct for all obliquities.
!
! Notes from Glen Double (31 Oct 2002):
! Method: define flow velocity and particle momentum in component vectors
!   wrt the plasma xyz frame.  Find p_para (wrt u) component of particle
!   momentum using scalar product of p and u.  Then p_perp = p - p_para
!   (using vectors).  Next, make relativistic transformation on p_para while
!   p_perp remains constant.  Finally, recreate p_rel = p_rel_para + p_perp,
!   and find p_rel components in the original xyz frame by taking scalar
!   products along the xyz axes.
!
! Input arguments:
!  1) aa: particle atomic mass
!  2) uxsk/old: current and old bulk flow speed along x axis
!  3) uzsk/old: current and old bulk flow speed along z axis
!  4) utot/old: current and old total bulk flow speed
!  5) gam_usf/old: Lorentz factor associated with utot/old
!  6) b_cos_th/old: current and old component of magnetic field along x axis
!  7) b_sin_th/old: current and old component of magnetic field along z axis
! Output arguments:
!  1) ptpf: total plasma frame momentum in new grid zone
!  2) ptsk: total shock frame momentum in new grid zone
!  3) p*sk: xyz components of ptsk
!  4) pbsk: component of ptsk parallel to magnetic field
!  5) p_perp_b_sk: component of ptsk perpendicular to magnetic field
!  6) gam_ptsk: Lorentz factor associated with ptsk
! Input/output arguments:
!  1) pbpf: component of ptpf parallel to magnetic field
!  2) p_perp_b_pf: component of ptpf perpendicular to magnetic field
!  3) gam_ptpf: Lorentz factor associated with ptpf
!  4) phi_rad: phase angle of gyration; looking UpS, counts clockwise from
!    +z axis

use constants
use species_vars, only: o_o_mc

implicit none

  ! Input arguments
real(kind=8), intent(in) :: aa, uxsk_old, uzsk_old, utot_old, gam_usf_old,&!&
     b_cos_old, b_sin_old, uxsk, uzsk, utot, gam_usf, b_cos_th, b_sin_th
  ! Output arguments
real(kind=8), intent(out) :: ptpf, ptsk, pxsk, pysk, pzsk, pbsk,          &!&
     p_perp_b_sk, gam_ptsk
  ! Input/Output arguments
real(kind=8), intent(inout) :: pbpf, p_perp_b_pf, gam_ptpf, phi_rad
  
  ! Local variables
real(kind=8) :: phi_p, p_p_cos, pxpf, pypf, pzpf, p_p_z, p_p_y


  phi_p = phi_rad + 0.5d0*pii
  
  p_p_cos = p_perp_b_pf * cos(phi_p)
  
  
!--! xyz (old) plasma frame components
  pxpf = pbpf*b_cos_old  -  p_p_cos*b_sin_old
  pypf = p_perp_b_pf * sin(phi_p)
  pzpf = pbpf*b_sin_old  +  p_p_cos*b_cos_old
  
  
!--! xyz (new) shock frame components
  pxsk = ( (gam_usf_old-1.d0)*(uxsk_old/utot_old)**2 + 1.d0 ) * pxpf      &!&
        + (gam_usf_old-1.d0)*(uxsk_old*uzsk_old/utot_old**2) * pzpf       &!&
        +  gam_usf_old * gam_ptpf * aa*xmp * uxsk_old
  
  pysk = pypf
  
  pzsk = (gam_usf_old-1.d0)*(uxsk_old*uzsk_old/utot_old**2) * pxpf        &!&
        + ( (gam_usf_old-1.d0)*(uzsk_old/utot_old)**2 + 1.d0 ) * pzpf     &!&
        + gam_usf_old * gam_ptpf * aa*xmp * uzsk_old


!--! Parallel/perpendicular (new) shock frame components
  ptsk = sqrt( pxsk**2  +  pysk**2  +  pzsk**2 )
  
  pbsk = pxsk*b_cos_th  +  pzsk*b_sin_th
  
  if( ptsk .lt. abs(pbsk) ) then
    p_perp_b_sk = 1.d-6 * ptsk
    pbsk        = sign( sqrt( ptsk**2 - p_perp_b_sk**2 ), pbsk )
    write(*,"(A)") 'Warning: ptsk < pbsk in Xform_p_PSP'
  else
    p_perp_b_sk = sqrt( ptsk**2  -  pbsk**2 )
  endif
  
  gam_ptsk = sqrt( (ptsk*o_o_mc)**2  +  1.d0 )
  
  
!--! xyz (new) plasma frame components
  pxpf = ( (gam_usf - 1.d0) * (uxsk/utot)**2  +  1.d0 ) * pxsk            &!&
        + (gam_usf - 1.d0) * (uxsk*uzsk/utot**2) * pzsk                   &!&
        -  gam_usf * gam_ptsk * aa*xmp * uxsk
  
  pypf = pysk
  
  pzpf = (gam_usf - 1.d0) * (uxsk*uzsk/utot**2) * pxsk                    &!&
        + ( (gam_usf - 1.d0) * (uzsk/utot)**2  +  1.d0 ) * pzsk           &!&
        -  gam_usf * gam_ptsk * aa*xmp * uzsk
  
  
!--! Parallel/perpendicular (new) plasma frame components, including new
   !   phase angle
  ptpf = sqrt( pxpf**2  +  pypf**2  +  pzpf**2 )
  
  pbpf = pxpf*b_cos_th  +  pzpf*b_sin_th
  
  if( ptpf .lt. abs(pbpf) ) then
    p_perp_b_pf = 1.d-6 * ptpf
    pbpf        = sign( sqrt( ptpf**2 - p_perp_b_pf**2 ), pbpf )
    write(*,"(A)") 'Warning: ptpf < pbpf in Xform_p_PSP'
  else
    p_perp_b_pf = sqrt( ptpf**2  -  pbpf**2 )
  endif
  
  gam_ptpf = sqrt( (ptpf*o_o_mc)**2  +  1.d0 )
  
  p_p_z = -pxpf*b_sin_th  +  pzpf*b_cos_th
  p_p_y =  pypf
  
    ! atan2(y,x) gives value of atan(y/x) but in correct quadrant
    ! See Figure 14 of Ellison, Baring, Jones (1996) [1996ApJ...473.1029E]
    !   for more details on phi_p
  phi_p = atan2( p_p_y, p_p_z )
  
  phi_rad = phi_p - 0.5d0*pii

return
end subroutine Xform_p_PSP


!****************************************************************************
!****************************************************************************
subroutine xn_per_shift(xn_per, x_PT_cm, gyro_rad_tot_cm)

! Shifts between coarse and fine values of xn_per
!
! Input arguments:
!  1) x_PT_cm: particle position
!  2) gyro_rad_tot_cm: particle gyroradius -- scattering mfp in Bohm limit
! Input/output arguments:
!  1) xn_per: number of time steps to divide a gyroperiod into

use controls, only: xn_per_coarse, xn_per_fine

implicit none

  ! Input arguments
real(kind=8), intent(in) :: x_PT_cm, gyro_rad_tot_cm
  ! Input/output arguments
real(kind=8), intent(inout) :: xn_per
  
  
!--! Shift between coarse and fine
  if( x_pt_cm .gt. gyro_rad_tot_cm ) then
    xn_per = xn_per_coarse
  else
    xn_per = xn_per_fine
  endif

return
end subroutine xn_per_shift


!****************************************************************************
!****************************************************************************
subroutine get_normalized_dNdp(dNdp_therm, dNdp_therm_pvals, dNdp_cr,     &!&
    nc_unit, do_therm_pt_formatted, zone_pop)

! Computes the actual number of particles in each bin of dN/dp (which is
!  divided by dp, remember).  Computes total area under curve, normalizes it
!  against number of particles upstream using plasma-frame density and
!  volume, and then uses fractional area of each bin to determine number of
!  particles in it.
! Handles non-injected (i.e. thermal) particles differently than injected
!  ones, due to very small p.f. spread in momenta.
!
! The array that would be dNdp_cr_pvals is already set, as psd_mom_bounds
!
! Inputs:
!  1) nc_unit: unit number for the scratch file holding crossing data
!  2) do_therm_pt_formatted: logical telling whether the scratch file is
!    human-readable or not
! Outputs:
!  1) dNdp_therm: 3-D array, containing 1-D array for each grid zone of
!   dN/dp for the thermal particles
!  2) dNdp_therm_pvals: array of momentum bin boundaries for each row of
!   dNdp_therm; each row handled separately to maximize resolution of what
!   may be an extremely narrow peak at radically different energy from the
!   upstream population
!  3) dNdp_cr: 3-D array, containing 1-D array for each grid zone of dN/dp
!   for the population of accelerated particles
!  4) zone_pop: (Lorentz-invariant) number of particles in each grid zone

use constants
use parameters, only: psd_max, na_g, num_therm_bins
use controls, only: l_oblique, jet_rad_pc, jet_sph_frac, denZ_ion, beta_Z,&!&
     gam_Z, n_ions, do_multi_dNdps
use psd_vars, only: num_psd_mom_bins, psd_mom_bounds
use grid_vars, only: n_grid, x_grid_cm, uxsk_grid
use iteration_vars, only: i_itr
use species_vars, only: i_ion
! DEBUGLINE
use controls, only: aa_ion
use grid_vars, only: gam_sf_grid
use debug

implicit none

  ! Input variables
integer, intent(in) :: nc_unit
logical, intent(in) :: do_therm_pt_formatted
  ! Output variables
real(kind=8), dimension(na_g), intent(out) :: zone_pop
real(kind=8), dimension(0:psd_max,na_g,3), intent(out) :: dNdp_therm,     &!&
     dNdp_therm_pvals, dNdp_cr

  ! Local variables
character(len=30) :: filename
integer :: num_hist_bins, i, i_shock, m, j, j_plot
real(kind=8) :: rad_min, jet_rad_cm, rad_max, dwell_time, flux_UpS,       &!&
     area_tot_therm, area_tot_cr, area_tot, norm_factor
real(kind=8), dimension(na_g) :: surf_area
real(kind=8), dimension(0:psd_max,3) :: therm_temp, therm_pvals_temp,     &!&
     dNdp_therm_rebin


!--! Warning (stop) statement about code
  if( l_oblique ) then
    write(*,"(2A)") 'ERROR in get_dNdp subrs: Lorentz Xforms not ',       &!&
        'written to handle oblique shocks'
    write(*,"(A)") 'Stopping program now.'
    stop
  endif


!--! Administrative constants to be used during main computation loops
   !-------------------------------------------------------------------------
  num_hist_bins = num_therm_bins / 2 ! Number of bins to use in histogram for
                                     !  ptpf & ctpf; needs to be divisor of
                                     !  num_therm_bins to minimize binning
                                     !  artifacts in this subroutine
  jet_rad_cm = jet_rad_pc * pc_to_cm ! Jet radius in cm, of course


!--! Get the non-normalized dN/dp's from the crossing arrays, scratch file
   !   (if necessary), and the phase space distribution,
  call get_dNdp_therm(num_hist_bins, dNdp_therm, dNdp_therm_pvals,        &!&
     nc_unit, do_therm_pt_formatted)
  
  call get_dNdp_cr(dNdp_cr)


!--! Now have non-normalized dN/dp for both thermal and CR populations.
   !   Determine the total number of particles in each grid zone by using
   !   shock frame flux, area, crossing time:
   !       #  =  flux * area * (distance/speed)
  do i = 1, n_grid
    if( (x_grid_cm(i) .eq. 0.d0) .or.                                     &!&
        (x_grid_cm(i+1)*x_grid_cm(i) .lt. 0.d0) ) then
      ! Either current grid zone is exactly at shock, or current and next
      !  grid zones straddle shock
      i_shock = i
      exit
    endif
  enddo
  
  ! Work upstream from shock and find volume of each grid zone
  rad_min = jet_rad_cm - x_grid_cm(i_shock) ! Inner radius, including fact
                                            !   that upstream has x < 0
  do i = i_shock-1, 1, -1
    ! outer radius is inner radius plus zone width *in ISM frame*
    rad_max = rad_min + (x_grid_cm(i+1) - x_grid_cm(i))/gam_Z
  
    ! Use rad_max and rad_min to get area of jet surface; will be same in
    !   shock frame if shock motion is purely parallel to shock normal
    surf_area(i) = 4.d0*pii  *  ( 0.5d0 * (rad_max+rad_min) )**2          &!&
                  * jet_sph_frac
    
    ! finally, set rad_min for next cycle through loop
    rad_min = rad_max
  enddo
  
  ! Work downstream from shock and find volume of each grid zone
  rad_max = jet_rad_cm - x_grid_cm(i_shock)
  do i = i_shock, n_grid
    ! inner radius is outer radius minus zone width *in ISM frame*
    rad_min = rad_max - (x_grid_cm(i+1) - x_grid_cm(i))/gam_Z
  
    ! Use rad_max and rad_min to get area of jet surface
    surf_area(i) = 4.d0*pii  *  ( 0.5d0 * (rad_max+rad_min) )**2          &!&
                  * jet_sph_frac
    
    ! finally, set rad_max for next cycle through loop
    rad_max = rad_min
  enddo
  
  
!--! Now calculate the number of particles in each zone, i.e. the total area
   !   under dN/dp
  do i = 1, n_grid
  
    ! Crossing time in the shock frame, upstream particle flux in the shock
    !  frame (conserved quantity throughout the shock structure), and
    !  finally number of particles in this grid zone
    dwell_time = (x_grid_cm(i+1) - x_grid_cm(i)) / uxsk_grid(i)
    
    flux_UpS   = gam_Z * denZ_ion(i_ion) * beta_Z*ccgs
    
    zone_pop(i) = flux_UpS * surf_area(i) * dwell_time
    
    !DEBUGLINE
    den_pf = gam_Z * uxsk_grid(1) / (gam_sf_grid(i) * uxsk_grid(i))
    zone_vol(i) = zone_pop(i) / den_pf
  enddo


!--! For each grid zone, integrate the area under dNdp_therm and dNdp_cr to
   !   find the non-normalized total area under the two curves. Then
   !   normalize each bin of the dNdp curves to get the correct number of
   !   particles in each grid zone.
   ! IMPORTANT: the differential term is dp, not p^2 dp! Angular components
   !   of momentum (the p^2 sin(theta) dphi dtheta) have already been
   !   handled.
   ! Outer loop manages each of the three possible frames of interest:
   !      1  -  Shock frame
   !      2  -  Plasma frame
   !      3  -  ISM frame
   !-------------------------------------------------------------------------
  do m = 1, 3
    do i = 1, n_grid
      
      area_tot_therm = 0.d0
      area_tot_cr    = 0.d0
      area_tot       = 0.d0
      
    !--! Calculate the total area under the two curves
      do j = 0, num_hist_bins-1
        if( dNdp_therm(j, i, m) .gt. 1.d-99 )               &!& ! Increment
          area_tot_therm = area_tot_therm                   &!& ! by
                          +  dNdp_therm(j, i, m)            &!& ! dN/dp
                           * ( dNdp_therm_pvals(j+1, i, m)  &!& ! *
                              -  dNdp_therm_pvals(j, i, m) )    ! dp
      enddo
      do j = 0, num_psd_mom_bins
        if( dNdp_cr(j, i, m) .gt. 1.d-99 )                  &!& ! Increment
          area_tot_cr = area_tot_cr                         &!& ! by
                       +  dNdp_cr(j, i, m)                  &!& ! dN/dp
                        * ( 10.d0**psd_mom_bounds(j+1)      &!& ! *
                           -  10.d0**psd_mom_bounds(j) )        ! dp
      enddo
    
    
   !--! Re-scale each curve according to the normalization factor. DO NOT
      !   include dp, for easier comparison against previous subroutines that
      !   performed similar tasks; the dp can be added back in as necessary
      !   during integration of dN/dp.
      ! When dealing with fast push, need to include thermal particle
      !   population even for zones upstream of fast push location.
      !   Fortunately, area_tot_therm in this case is approximately 
      !   den_pf/u_x, the compressed plasma-frame density divided by the
      !   local shock speed.
      !DOLATER: plot area_tot_therm against local density to see if this
      !   holds even once particles start being heated
      if( (area_tot_therm .eq. 0.d0) .and. (area_tot_cr .gt. 0.d0) ) then
        den_pf = denZ_ion(i_ion) * gam_Z * uxsk_grid(1)                   &!& 
                / (gam_sf_grid(i) * uxsk_grid(i))
        area_tot = den_pf/uxsk_grid(i) + area_tot_cr
      else
        area_tot = area_tot_therm + area_tot_cr
      endif
      
      if( area_tot .gt. 0.d0 ) then
        norm_factor = zone_pop(i) / area_tot
      else
        norm_factor = 0.d0
      endif
      
      do j = 0, num_hist_bins-1
        if( dNdp_therm(j, i, m) .gt. 1.d-99 )             &!&
          dNdp_therm(j, i, m) = dNdp_therm(j, i, m)       &!& ! dN/dp
                                    *  norm_factor            ! renormalized
      enddo
      do j = 0, num_psd_mom_bins
        if( dNdp_cr(j, i, m) .gt. 1.d-99 )                &!&
          dNdp_cr(j, i, m) = dNdp_cr(j, i, m)             &!& ! dN/dp
                                 *  norm_factor               ! renormalized
      enddo
      
      
    enddo ! loop over grid zones
  
  enddo ! loop over frames
   !-------------------------------------------------------------------------
   ! Normalized dN/dp's found


!--! Plot the dN/dps as a check-by-eye on their reasonability.  Include every
   !   grid zone and frame, which is probably more data than necessary in
   !   most cases
   !-------------------------------------------------------------------------
  if( do_multi_dNdps ) then
    if( i_ion .eq. 1 ) then
      write(filename,"(A19,I0,A4)") "mc_dNdp_grid_therm_", i_itr, ".dat"
      open(unit=507,file=trim(filename))
      
      write(filename,"(A19,I0,A4)") "mc_dNdp_grid_CR_", i_itr, ".dat"
      open(unit=517,file=trim(filename))
    endif
  else
    if(i_ion .eq. 1) open(unit=507,file="mc_dNdp_grid_therm.dat")
    if(i_ion .eq. 1) open(unit=517,file="mc_dNdp_grid_CR.dat")  
  endif
  
  do i = 1, n_grid
    
  !--! Thermal particles, all three frames
    if( maxval(dNdp_therm(:,i,:)) .gt. 1.d-66 ) then
        !DEBUGLINE
        do j = 0, num_hist_bins-1
          if(dNdp_therm(j,i,2) .gt. 1.d-66) then
            p_avg = dNdp_therm_pvals(j,i,2) + dNdp_therm_pvals(j+1,i,2)
            p_avg = 0.5d0 * p_avg
            en_avg = sqrt( (aa_ion(i_ion)*rm_prot)**2  +  (p_avg*ccgs)**2 )
            
            therm_en_den(i,i_ion) = therm_en_den(i,i_ion)                 &!&
                 + dNdp_therm(j, i, 2)                                    &!&
                  * (dNdp_therm_pvals(j+1,i,2) - dNdp_therm_pvals(j,i,2)) &!&
                  * (en_avg - aa_ion(i_ion)*rm_prot)
            en_den(i,i_ion) = en_den(i,i_ion)                 &!&
                 + dNdp_therm(j, i, 2)                                    &!&
                  * (dNdp_therm_pvals(j+1,i,2) - dNdp_therm_pvals(j,i,2)) &!&
                  * (en_avg - aa_ion(i_ion)*rm_prot)
          endif
        enddo
      j_plot  = 0
      do j = 0, num_hist_bins-1
      
        j_plot = j_plot+1
        write(507,"(2I4,F5.1,20ES13.4E2)") i, j_plot,    &!&
          float(i_ion),                                  &!& ! 1
            ! Shock frame
          log10(dNdp_therm_pvals(j,i,1)),                &!& ! 2 (cgs units)
          log10(dNdp_therm_pvals(j,i,1) / (xmp*ccgs)),   &!& ! 3 (nat. units)
          log10(dNdp_therm(j, i, 1)),                    &!& ! 4
            ! Plasma frame
          log10(dNdp_therm_pvals(j,i,2)),                &!& ! 5 (cgs units)
          log10(dNdp_therm_pvals(j,i,2) / (xmp*ccgs)),   &!& ! 6 (nat. units)
          log10(dNdp_therm(j, i, 2)),                    &!& ! 7
            ! ISM frame
          log10(dNdp_therm_pvals(j,i,3)),                &!& ! 8 (cgs units)
          log10(dNdp_therm_pvals(j,i,3) / (xmp*ccgs)),   &!& ! 9 (nat. units)
          log10(dNdp_therm(j, i, 3))                         ! 10
        
        j_plot = j_plot+1
        if(j .lt. (num_hist_bins-1)) then
          write(507,"(2I4,F5.1,20ES13.4E2)") i, j_plot,        &!&
            float(i_ion),                                      &!& ! 1
              ! Shock frame
            log10(dNdp_therm_pvals(j+1,i,1)),                  &!& ! 2 (cgs)
            log10(dNdp_therm_pvals(j+1,i,1) / (xmp*ccgs)),     &!& ! 3 (nat.)
            log10(dNdp_therm(j, i, 1)),                        &!& ! 4
              ! Plasma frame
            log10(dNdp_therm_pvals(j+1,i,2)),                  &!& ! 5 (cgs)
            log10(dNdp_therm_pvals(j+1,i,2) / (xmp*ccgs)),     &!& ! 6 (nat.)
            log10(dNdp_therm(j, i, 2)),                        &!& ! 7
              ! ISM frame
            log10(dNdp_therm_pvals(j+1,i,3)),                  &!& ! 8 (cgs)
            log10(dNdp_therm_pvals(j+1,i,3) / (xmp*ccgs)),     &!& ! 9 (nat.)
            log10(dNdp_therm(j, i, 3))                             ! 10
        endif
      enddo ! loop over num_hist_bins
      
      call print_plot_vals(507)
      
    endif  ! check on existence of thermal particles
    
    
  !--! Cosmic rays
    if( maxval(dNdp_cr(:,i,:)) .gt. 1.d-66 ) then
      
      ! Convert thermal particles, if any,  to bins used for cosmic rays
      if( maxval(dNdp_therm(:,i,:)) .gt. 1.d-66 ) then
        !DEBUGLINE
        do j = 0, num_psd_mom_bins
          if(dNdp_cr(j,i,2) .gt. 1.d-66) then
            p_avg = 10.d0**psd_mom_bounds(j+1)-10.d0**psd_mom_bounds(j)
            p_avg = sqrt( p_avg )
            en_avg = sqrt( (aa_ion(i_ion)*rm_prot)**2  +  (p_avg*ccgs)**2 )
            
            en_den(i,i_ion) = en_den(i,i_ion)                 &!&
                 + dNdp_cr(j, i, 2)                                    &!&
                  * (10.d0**psd_mom_bounds(j+1)-10.d0**psd_mom_bounds(j)) &!&
                  * (en_avg - aa_ion(i_ion)*rm_prot)
          endif
        enddo
        therm_temp(:,:)       = dNdp_therm(:,i,:)
        therm_pvals_temp(:,:) = dNdp_therm_pvals(:,i,:)
        
        call rebin_dNdp_therm(num_hist_bins, therm_temp, therm_pvals_temp,&!&
           num_psd_mom_bins, psd_mom_bounds, dNdp_therm_rebin)
      else
        dNdp_therm_rebin(:,:) = 1.d-99
      endif
      
      j_plot  = 0
      do j = 0, num_psd_mom_bins
        j_plot = j_plot+1
        write(517,"(2I4,F5.1,20ES13.4E2)") i, j_plot,    &!&
          float(i_ion),                                  &!& ! 1
                psd_mom_bounds(j),                       &!& ! 2 (cgs units)
                psd_mom_bounds(j) - log10(xmp*ccgs),     &!& ! 3 (nat. units)
            ! Shock frame
          log10(dNdp_cr(j, i, 1)),                       &!& ! 4
            ! Plasma frame
          log10(dNdp_cr(j, i, 2)),                       &!& ! 5
            ! ISM frame
          log10(dNdp_cr(j, i, 3)),                       &!& ! 6
            ! Summed therm+CR dN/dps in all three frames
          log10(dNdp_cr(j,i,1)+dNdp_therm_rebin(j,1)),   &!& ! 7
          log10(dNdp_cr(j,i,2)+dNdp_therm_rebin(j,2)),   &!& ! 8
          log10(dNdp_cr(j,i,3)+dNdp_therm_rebin(j,3))        ! 9
        
        j_plot = j_plot+1
        if(j .lt. num_psd_mom_bins) then
          write(517,"(2I4,F5.1,20ES13.4E2)") i, j_plot,       &!&
            float(i_ion),                                     &!& ! 1
                  psd_mom_bounds(j+1),                        &!& ! 2 (cgs)
                  psd_mom_bounds(j+1) - log10(xmp*ccgs),      &!& ! 3 (nat.)
              ! Shock frame
            log10(dNdp_cr(j, i, 1)),                          &!& ! 4
              ! Plasma frame
            log10(dNdp_cr(j, i, 2)),                          &!& ! 5
              ! ISM frame
            log10(dNdp_cr(j, i, 3)),                          &!& ! 6
              ! Summed therm+CR dN/dps in all three frames
            log10(dNdp_cr(j,i,1)+dNdp_therm_rebin(j,1)),      &!& ! 7
            log10(dNdp_cr(j,i,2)+dNdp_therm_rebin(j,2)),      &!& ! 8
            log10(dNdp_cr(j,i,3)+dNdp_therm_rebin(j,3))           ! 9
        endif
      enddo ! loop over num_psd_mom_bins
      
      call print_plot_vals(517)
      
    endif  ! check on existence of CRs
    
  enddo  ! loop over grid zones
  
  if( i_ion .eq. n_ions ) close(507)
  if( i_ion .eq. n_ions ) close(517)
   !-------------------------------------------------------------------------
   ! end of plotting section

return
end subroutine get_normalized_dNdp


!****************************************************************************
!****************************************************************************
subroutine get_dNdp_2D(nc_unit, do_therm_pt_formatted, zone_pop,          &!&
     d2N_dpdcos_ef)

! Generates d2N/dpdcos (or something like it) for electrons in the explosion
!   frame.
!
! First, use the scratch file and the phase space distribution to generate
!  d2N/(dp-dcos), the 2-D version of dN/dp (with extra information about the
!  angular distribution).  Combine the results from thermal and CR particles
!  into a single d2N/(dp-dcos) to save memory.
! After d2N/(dp-dcos) is found in the shock frame, transform it to the plasma
!  and ISM frames by rebinning the *center* of each shock frame bin.  Leave
!  the computationally difficult problem of overlap for a later date.
! Finally, condense each d2N/(dp-dcos) into dNdp for comparison against the
!  results of the original dNdp subroutines.
!
! WARNING: d2N/(dp-dcos) isn't technically accurate.  To get number, just
!  multiply by dp and sum over the cos(theta) column.
! NOTE: Right now, this is only used to calculate the ISM-frame distribution
!  of electrons; all other frame/particle pairs can be done with the 1-D
!  dN/dp subroutines.
!
! Input arguments:
!  1) nc_unit: unit number for the scratch file holding crossing data
!  2) do_therm_pt_formatted: logical telling whether the scratch file is
!    human-readable or not
!  3) zone_pop: (Lorentz-invariant) number of particles in each grid zone
! Output arguments:
!  1) d2N_dpdcos_ef: explosion-frame array holding dN/dp spread out across
!    angular dimension


use constants
use parameters, only: psd_max, na_g, na_cr
use controls, only: aa_ion, n_ions, denZ_ion, psd_lin_cos_bins, gam_Z,    &!&
     beta_Z
use psd_vars, only: psd, num_psd_tht_bins, psd_tht_bounds,                &!&
     num_psd_mom_bins, psd_mom_bounds
use grid_vars, only: n_grid, gam_sf_grid
use species_vars, only: i_ion, num_crossings, n_cr_count, therm_grid,     &!&
     therm_pxsk, therm_ptsk, therm_xwt

implicit none


  ! Input arguments
integer, intent(in) :: nc_unit
logical, intent(in) :: do_therm_pt_formatted
real(kind=8), dimension(na_g), intent(in) :: zone_pop
  ! Output arguments
real(kind=8), dimension(0:psd_max, 0:psd_max, na_g), intent(out) ::       &!&
     d2N_dpdcos_ef

  ! Local variables
integer :: k_pt, max_cross, ntot_crossings, i, i_grid, idum, n, j_th,     &!&
     m_max, m, kpt_Xf, jth_Xf
real(kind=8) :: rest_mass_en, p_lo_cgs, p_hi_cgs, ptsk, pxsk, cell_wt,    &!&
     norm_factor, cos_hi, cos_lo, gam_u, beta_u, cos_theta_sf, ptsf_cgs,  &!&
     pxsf_cgs, etot_sf_cgs, px_Xf, pt_Xf
integer, dimension(na_g) :: n_cross_fill
real(kind=8), dimension(na_g) :: den_tot
real(kind=8), dimension(0:psd_max) :: del_p, o_o_del_p, cos_center,       &!&
     pt_cgs_center
real(kind=8), dimension(0:psd_max, 0:psd_max, na_g) :: d2N_dpdcos_sf, d2N_pf
real(kind=8), dimension(:,:), allocatable :: therm_px, therm_pt, therm_wt


!--! Initialize d2N_dpdcos_sf to "zero"
   ! Dimensions of d2N_dpdcos arrays (chosen for maximum speed during loops):
   !    1 - cos(theta)
   !    2 - ptot
   !    3 - grid position
   ! The omission of "dpdcos" in the plasma frame array is not a typo; the
   !   array will only ever hold a total count of particles, not a true dN/dp
  d2N_dpdcos_sf(:,:,:) = 1.d-99
  d2N_pf(:,:,:)        = 1.d-99
  d2N_dpdcos_ef(:,:,:) = 1.d-99


!--! Set constants to be used repeatedly
  rest_mass_en = aa_ion(i_ion) * rm_prot  ! Rest mass-energy of the
                                          !  current particle species
  do k_pt = 0, num_psd_mom_bins
    p_lo_cgs = 10.d0**psd_mom_bounds(k_pt)
    p_hi_cgs = 10.d0**psd_mom_bounds(k_pt+1)
    del_p(k_pt)     = p_hi_cgs - p_lo_cgs
    o_o_del_p(k_pt) = 1.d0 / del_p(k_pt)
  enddo


! Read the scratch files associated with thermal particles.  Bin them into
!  the combined d2N/(dp-dcos) for the shock frame
!----------------------------------------------------------------------------
!--! Set up the arrays that will hold the crossing data
  max_cross = maxval(num_crossings,1)
  allocate( therm_px(max_cross, 0:n_grid+1) )
  allocate( therm_pt(max_cross, 0:n_grid+1) )
  allocate( therm_wt(max_cross, 0:n_grid+1) )


!--! Fill the arrays with data from the scratch file
   !-------------------------------------------------------------------------
  ntot_crossings = sum(num_crossings,1)
  n_cross_fill(:) = 0 ! n_cross_fill needs to be an array because we will be
                      !  skipping around the grid as we move through the data
                      !  rather than filling a single zone at a time
  
  rewind(nc_unit)
  
    
  ! Handle crossings stored within the crossing arrays
  do i = 1, n_cr_count
    i_grid = therm_grid(i)
    
    n_cross_fill(i_grid) = n_cross_fill(i_grid) + 1
  
    ! Note: ordering of coordinates chosen to make memory accesses in next
    !  loop faster
    therm_px(n_cross_fill(i_grid), i_grid) = therm_pxsk(i)
    therm_pt(n_cross_fill(i_grid), i_grid) = therm_ptsk(i)
    therm_wt(n_cross_fill(i_grid), i_grid) = therm_xwt(i)
  enddo
    
  if( ntot_crossings .gt. na_cr ) then
    
    ! Need to go into the scratch file for the remainder of the crossings
    do i = 1, ntot_crossings-na_cr
      if( do_therm_pt_formatted ) then
        read(nc_unit,"(I3,I5,4ES20.12E2)") i_grid, idum, pxsk, ptsk, cell_wt
      else
        read(nc_unit) i_grid, idum, pxsk, ptsk, cell_wt
      endif
      
      n_cross_fill(i_grid) = n_cross_fill(i_grid) + 1
      
      ! Note: ordering of coordinates chosen to make memory accesses in next
      !  loop faster
      therm_px(n_cross_fill(i_grid), i_grid) = pxsk
      therm_pt(n_cross_fill(i_grid), i_grid) = ptsk
      therm_wt(n_cross_fill(i_grid), i_grid) = cell_wt
    enddo
    
  endif
  
  ! With the read-in complete, the scratch file can be closed
  close(nc_unit)
   !-------------------------------------------------------------------------
   ! Arrays filled
  
  
!--! Loop over grid locations, filling in the array d2N_dpdcos_sf
   !-------------------------------------------------------------------------
  do i = 1, na_g
    if( num_crossings(i) .eq. 0 ) then
      d2N_dpdcos_sf(:,:,i) = 1.d-99
      cycle     ! Ignore zones that no thermal particles crossed
    endif
    
    ! Loop over all crossings for this grid zone, binning the thermal
    !  particles.
    do n = 1, num_crossings(i)
      call get_psd_bins(therm_px(n,i), therm_pt(n,i), .true.,  k_pt, j_th)
      
      d2N_dpdcos_sf(j_th, k_pt, i) = d2N_dpdcos_sf(j_th, k_pt, i)         &!&
                                    +  therm_wt(n,i)
    enddo
  
  enddo
   !-------------------------------------------------------------------------
   ! Shock frame d2N_dpdcos filled from thermal arrays


!--! Deallocate the crossing data arrays
  deallocate( therm_px ) ;  deallocate( therm_pt ) ;  deallocate( therm_wt )
!----------------------------------------------------------------------------
! Thermal particles binned


!--! With the thermal particles taken care of in the shock frame, move on to
   !   the cosmic ray particles.  Note that the slices of psd need to be
   !   transposed from (ptot,theta) to (theta,ptot) order.
   ! Once that is taken care of, convert d2N into dN/dp by dividing by dp
   !-------------------------------------------------------------------------
  do i = 1, n_grid+1
    do k_pt = 0, psd_max-1
  
      do j_th = 0, psd_max-1
        
        if( psd(k_pt, j_th, i) .gt. 1.d-66 ) then
          d2N_dpdcos_sf(j_th, k_pt, i) = d2N_dpdcos_sf(j_th, k_pt, i)     &!&
                                        +  psd(k_pt, j_th, i)
        endif
        
      enddo
    enddo
  enddo
   !-------------------------------------------------------------------------
   ! Cosmic rays finished


!--! Calculate the density of particles represented by d2N_dpdcos_sf
  do i = 1, n_grid
    
    den_tot(i) = sum( d2N_dpdcos_sf(:,0:psd_max-1,i),                     &!&
                      MASK = (d2N_dpdcos_sf(:,0:psd_max-1,i) .gt. 1.d-66) )
    
  enddo
  
  
!--! Can convert from dN to dN/dp now
  do k_pt = 0, psd_max-1
    do j_th = 0, psd_max
    
      if( d2N_dpdcos_sf(j_th, k_pt, i) .gt. 1.d-66 ) then
        d2N_dpdcos_sf(j_th, k_pt, i) = d2N_dpdcos_sf(j_th, k_pt, i)       &!&
                                      *  o_o_del_p(k_pt)
      endif
      
    enddo
  enddo


!--! Re-scale d2N_dpdcos_sf according to the normalization factor.
   ! When dealing with fast push, need to count thermal particle population
   !   even for zones upstream of fast push location.  Fortunately the number
   !   of particles in this case is approximately denZ_ion(i_ion).  Note that
   !   this assumes minimal shock modification upstream of the fast push zone
   !   -- if (e.g. for nonrel shocks) the shock is extensively modified this
   !   will not be correct.
  do i = 1, n_grid
    
    ! No thermal particles measured, so add them in proper amount
    if( (num_crossings(i) .eq. 0) .and. (den_tot(i) .gt. 0.d0) ) then
      den_tot(i) = den_tot(i)  +  denZ_ion(i_ion)
    endif
    
    ! Determine rescaling factor
    if( den_tot(i) .gt. 0.d0 ) then
      norm_factor = zone_pop(i) / den_tot(i)
    else
      norm_factor = 0.d0
    endif
    
    ! Renormalize all nonempty cells of d2N_dpdcos_sf
    do k_pt = 0, psd_max
      do j_th = 0, psd_max
        if( (d2N_dpdcos_sf(j_th,k_pt,i) .gt. 1.d-99) .and.                &!&
            (norm_factor .gt. 0.d0) ) then
          d2N_dpdcos_sf(j_th,k_pt,i) = d2N_dpdcos_sf(j_th,k_pt,i)         &!&
                                      *  norm_factor
        else
          d2N_dpdcos_sf(j_th,k_pt,i) = 1.d-99
        endif
      enddo
    enddo
  enddo


! Now, generate d2N_dpdcos_ef and d2n_dpdcos_pf.  Do this by transforming the
!   bin *centers* of d2N_dpdcos_sf into each frame, using just the center
!   location to rebin.
! DOLATER: consider doing an exact calculation of bin overlaps rather than
!   just the center-point rebinning that currently happens
!----------------------------------------------------------------------------
!--! Calculate the center points of all the bins to save time later
  do j_th = 0, num_psd_tht_bins
    ! Determine current cosines, remembering that psd_tht_bounds has both a
    !  linearly-spaced region in cosine and logarithmically-spaced region in
    !  theta.
    if( j_th .gt. (num_psd_tht_bins - psd_lin_cos_bins) ) then
      cos_hi = psd_tht_bounds(j_th)
      cos_lo = psd_tht_bounds(j_th + 1)
    else if( j_th .eq. (num_psd_tht_bins - psd_lin_cos_bins) ) then
      cos_hi = cos(psd_tht_bounds(j_th))
      cos_lo = psd_tht_bounds(j_th + 1)
    else
      cos_hi = cos(psd_tht_bounds(j_th))
      cos_lo = cos(psd_tht_bounds(j_th + 1))
    endif
    
    ! Minus sign needed because finest gradations actually point UpS
    cos_center(j_th) = -0.5d0 * ( cos_lo + cos_hi )
  enddo
  do k_pt = 0, num_psd_mom_bins
    pt_cgs_center(k_pt) = 0.5d0 * ( psd_mom_bounds(k_pt)                  &!&
                                   +  psd_mom_bounds(k_pt+1) )
    ! Convert from log to linear space
    pt_cgs_center(k_pt) = 10.d0**pt_cgs_center(k_pt)
  enddo


!--! Loop extents reflect desired frames:
   !   (1,1)    only plasma frame
   !   (1,2)    both plasma and explosion frame
   !   (2,2)    only explosion frame
  if(i_ion .lt. n_ions) then
     m_max = 1
  else
     m_max = 2
  endif
  do m = 1, m_max  ! Loop over plasma and/or ISM frames
    
    do i = 1, n_grid
      
    !--! Get Lorentz factor and speed relating the shock and current frames
      select case( m )
        case( 1 )   ! Plasma frame
          gam_u  = gam_sf_grid(i)
          beta_u = sqrt( 1.d0 - 1.d0/gam_u**2 )
        case( 2 )   ! ISM frame
          gam_u  = gam_Z
          beta_u = beta_Z
      end select
      
      ! Loop over cells in d2N_dpdcos_sf
      do k_pt = 0, num_psd_mom_bins
        do j_th = 0, num_psd_tht_bins
          
          ! Skip empty zones
          if( d2N_dpdcos_sf(j_th,k_pt,i) .le. 1.d-66 ) then
            cycle
          else
            ! Return dN/dp to dN by multiplying by dp
            cell_wt = d2N_dpdcos_sf(j_th, k_pt, i) * del_p(k_pt)
          endif
          
          
          ! Transform the center of the zone into the new frame
          cos_theta_sf = cos_center(j_th)
          ptsf_cgs     = pt_cgs_center(k_pt)
          pxsf_cgs     = ptsf_cgs * cos_theta_sf
          etot_sf_cgs  = sqrt( (ptsf_cgs*ccgs)**2  +  rest_mass_en**2 )
          
          px_Xf    = gam_u * (pxsf_cgs - beta_u*etot_sf_cgs/ccgs)
          pt_Xf    = sqrt( ptsf_cgs**2 - pxsf_cgs**2 + px_Xf**2 )
          
          
          ! Get location of center in transformed d2N_dpdcos
          call get_psd_bins(px_Xf, pt_Xf, .true., kpt_Xf, jth_Xf)
          
          
          ! And add shock frame number to the appropriate bin in target
          !   frame; note that d2N_pf is NOT being converted back to
          !   d2N/dp
          select case( m )
            case( 1 )   ! Plasma frame
              d2N_pf(jth_Xf, kpt_Xf, i) = d2N_pf(jth_Xf, kpt_Xf, i)   &!&
                                         +  cell_wt
            case( 2 )   ! ISM frame
              d2N_dpdcos_ef(jth_Xf, kpt_Xf, i) =                      &!&
                                   d2N_dpdcos_ef(jth_Xf, kpt_Xf, i)   &!&
                                  +  cell_wt * o_o_del_p(kpt_Xf)
          end select
          
          
        enddo ! loop over angles
      enddo ! loop over ptot
      
    enddo ! loop over grid locations
  enddo ! loop over reference frames
!----------------------------------------------------------------------------
! Transformed d2N_dpdcos calculated


!--! Create NetCDF copies of all three arrays, noting that the file name
   !  must include the particle species to avoid overwriting.
!comm   write(tmp1,"(i1)") i_ion
!comm     ! Shock frame
!comm   ncfilename = "dNdp_i" // tmp1 // "_sf"
!comm   call output_netcdf( d2N_dpdcos_sf(0:num_psd_tht_bins,             &!&
!comm                                     0:num_psd_mom_bins,             &!&
!comm                                     1:n_grid),                      &!&
!comm                       num_psd_mom_bins,                             &!&
!comm                       num_psd_tht_bins,                             &!&
!comm                       n_grid,                                       &!&
!comm                       ncfilename)
!comm     ! Plasma frame
!comm   ncfilename = "dNdp_i" // tmp1 // "_pf"
!comm   call output_netcdf( d2N_dpdcos_pf(0:num_psd_tht_bins,             &!&
!comm                                     0:num_psd_mom_bins,             &!&
!comm                                     1:n_grid),                      &!&
!comm                       num_psd_mom_bins,                             &!&
!comm                       num_psd_tht_bins,                             &!&
!comm                       n_grid,                                       &!&
!comm                       ncfilename)
!comm     ! ISM frame
!comm   ncfilename = "dNdp_i" // tmp1 // "_ef"
!comm   call output_netcdf( d2N_dpdcos_ef(0:num_psd_tht_bins,             &!&
!comm                                     0:num_psd_mom_bins,             &!&
!comm                                     1:n_grid),                      &!&
!comm                       num_psd_mom_bins,                             &!&
!comm                       num_psd_tht_bins,                             &!&
!comm                       n_grid,                                       &!&
!comm                       ncfilename)

return
end subroutine get_dNdp_2D


!****************************************************************************
!****************************************************************************
subroutine get_dNdp_cr(dNdp_cr)

! Calculate dN/dp (NOT normalized) in plasma frame for particles that *have*
!  been injected into acceleration process.
!
! Input arguments:
!  None, but does use PSD information from module psd_vars
! Output arguments:
!  1) dNdp_cr: 3-D array, containing 1-D array for each grid zone, of dN/dp
!    for particles injected into acceleration process, for three different
!    reference frames

use constants
use parameters, only: psd_max, na_g
use controls, only: aa_ion, psd_lin_cos_bins, gam_Z
use psd_vars, only: num_psd_tht_bins, psd_tht_bounds, psd, psd_mom_bounds,&!&
     num_psd_mom_bins
use grid_vars, only: n_grid, gam_sf_grid
use species_vars, only: i_ion

implicit none


! Input variables
! Output variables
real(kind=8), dimension(0:psd_max,na_g,3), intent(out) :: dNdp_cr

  ! Local variables
integer :: l, i_approx, k, j, i, m, i_ct_ptsk_min, i_ct_ptsk_max,         &!&
     i_ct_ptpf_min, i_ct_ptpf_max
real(kind=8) :: rest_mass, gam_u, o_o_gam_u, p_lo_cgs, p_hi_cgs
real(kind=8), dimension(0:psd_max) :: ct_bounds
real(kind=8), dimension(0:psd_max,0:psd_max) :: Xform_corn_pt, Xform_corn_ct
real(kind=8), allocatable, dimension(:,:) :: ct_sk_xw, ct_pf_xw, ct_ef_xw


!--! Zero out the output array before summing in loops to follow
   ! Note on dNdp_cr's third dimension:
   !   1  -  Shock frame
   !   2  -  Plasma frame
   !   3  -  ISM frame
  dNdp_cr(:,:,:) = 0.d0
  
  
!--! Rest mass of particle species, used here in binning pitch angle
  rest_mass = aa_ion(i_ion) * xmp


!--! Set values needed to bin pitch angle
  ct_bounds(:) = -2.d0
  do l = 0, num_psd_tht_bins+1
    ! Determine current cosine, remembering that psd_tht_bounds has both a
    !   linearly-spaced region in cosine and logarithmically-spaced region in
    !   theta.
    ! Also need to remember that the most finely spaced bins should occur in
    !   the UpS-pointing direction, so need to negate psd_tht_bounds to get
    !   true cosine value.
    if( l .gt. (num_psd_tht_bins - psd_lin_cos_bins) ) then
      ct_bounds(l) = -psd_tht_bounds(l)
    else
      ct_bounds(l) = -cos( psd_tht_bounds(l) )
    endif
  enddo


!--! Degree of approximation to use in distributing cell_wt if NOT in shock
   !   frame (where uniform distribution may be assumed):
   !    i_approx = 0:  assume uniform distribution
   !    i_approx = 1:  assume isosceles trianglular distribution of cell_wt
   !    i_approx = 2:  assume scalene triangular distribution of cell_wt,
   !                     with peak of triangle centered on mean of ct_hi_pt
   !                     and ct_lo_pt
   !    i_approx = 3:  exact calculation of fractional area, i.e. use no 
   !                     approximations
  i_approx = 2


!--! Loop over grid zones and cos(theta)/ptot space to build dN_cr
  do k = 1, n_grid

  !--! Handle shock frame separately, as no special treatment is needed to
     !  convert from PSD to dN(p).  Conversion to dN/dp will happen at the
     !  end of this subroutine.
     !-----------------------------------------------------------------------
    do j = 0, psd_max
      do i = 0, psd_max
        if( psd(i, j, k) .gt. 0.d0 )                     &!&
          dNdp_cr(i,k,1) = dNdp_cr(i,k,1) + psd(i,j,k)
      enddo
    enddo
     !-----------------------------------------------------------------------
     ! Shock frame dN(p) found
    
    
  !--! Now handle plasma and ISM frames.  Do this sequentially, since the
     !  process is exactly the same but requires a single input that differs
     !  between the two frames.
    do m = 2, 3
      
      ! Transform corners of PSD into correct frame for use in finding
      !  p_cell_lo and p_cell_hi in next block of code
      !----------------------------------------------------------------------
      select case( m )
        case( 2 )  !  Plasma frame
          gam_u     = gam_sf_grid(k)
          o_o_gam_u = 1.d0 / gam_u
        case( 3 )  !  ISM frame
          gam_u     = gam_Z
          o_o_gam_u = 1.d0 / gam_u
        case default
          write(*,"(A)") "ERROR with frame selection in get_dNdp_cr"
          stop
      end select
      
      call Xform_psd_corners(gam_u, Xform_corn_pt, Xform_corn_ct)
      !----------------------------------------------------------------------
      ! corners transformed
      
      
    !--! Determine minimum and maximum plasma-frame momenta particles can
       !  have.  Allocate the cos(theta) arrays based on those values, since
       !  particles will be binned according to both cos(theta) and (coarse-
       !  grained) momentum.  Then zero them out.
      i_ct_ptsk_min =   floor( psd_mom_bounds(1) )
      i_ct_ptsk_max = ceiling( psd_mom_bounds(num_psd_mom_bins+1) )
      if( allocated(ct_sk_xw) ) deallocate(ct_sk_xw)
      allocate( ct_sk_xw(0:psd_max, i_ct_ptsk_min:i_ct_ptsk_max) )
      ct_sk_xw(:,:) = 0.d0
      
      i_ct_ptpf_min =   floor( minval(                                    &!&
                                      Xform_corn_pt(1:num_psd_mom_bins+1, &!&
                                                    0:num_psd_tht_bins+1) ) )
      i_ct_ptpf_max = ceiling( maxval(                                    &!&
                                      Xform_corn_pt(1:num_psd_mom_bins+1, &!&
                                                    0:num_psd_tht_bins+1) ) )
      select case( m )
        case( 2 )
          if( allocated(ct_pf_xw) ) deallocate(ct_pf_xw)
          allocate( ct_pf_xw(0:psd_max, i_ct_ptpf_min:i_ct_ptpf_max) )
          ct_pf_xw(:,:) = 0.d0
        case( 3 )
          if( allocated(ct_ef_xw) ) deallocate(ct_ef_xw)
          allocate( ct_ef_xw(0:psd_max, i_ct_ptpf_min:i_ct_ptpf_max) )
          ct_ef_xw(:,:) = 0.d0
      end select
      
      
      ! Note that dNdp_cr as calculated here is dN(p), *not* dN/dp.  That
      !   conversion happens at the end of this subroutine
      call get_Xform_dN(psd(:,:,k), m, Xform_corn_pt, Xform_corn_ct,      &!&
         o_o_gam_u, i_approx, dNdp_cr(:,k,m) )

    enddo ! loop over frames
      
      
    ! If tracking pitch angles, and any bins were filled while doing so, plot
    !   histogram of their values
    ! Because pitch angle was broken down by momentum, loop over all columns
    !   that aren't empty in ct_**_xw.
    ! WARNING: this must be re-checked and re-tested since it was brought
    !   into new version of code
    !------------------------------------------------------------------------
!comm     ! First, histograms of individual momentum decades in the shock
!comm     !   frame
!comm     do i = i_ct_ptsk_min, i_ct_ptsk_max
!comm       
!comm       do l = 0, num_psd_tht_bins
!comm         if(ct_sk_xw(l,i) .gt. 0.d0)                                 &!&
!comm            ct_sk_xw(l,i) = ct_sk_xw(l,i)                            &!&
!comm                           / (ct_bounds(l) - ct_bounds(l+1))
!comm       enddo
!comm       sum_ct_sk_xw = sum( ct_sk_xw(:,i) )
!comm       if( sum_ct_sk_xw .eq. 0.d0 ) cycle   ! skip empty rows
!comm       
!comm       k_unit = k + 8000 + 100*(numion-1); j_plot = 0
!comm       do l = 0, num_psd_tht_bins
!comm         j_plot = j_plot + 1
!comm         
!comm         if(ct_bounds(l) .eq. 1.d0) then
!comm           theta = 1.d-99
!comm         else
!comm           theta = acos(ct_bounds(l))
!comm         endif
!comm         theta_next = acos(ct_bounds(l+1))
!comm         
!comm         write(k_unit,"(2i4,1p48e15.7)") k, j_plot,      &!& ! fort.80**
!comm                       ct_bounds(l),                     &!& ! 1
!comm                       theta,                            &!& ! 2
!comm                 log10(theta),                           &!& ! 3
!comm                  real(i),                               &!& ! 4
!comm                       ct_sk_xw(l,i)/sum_ct_sk_xw            ! 5
!comm         
!comm         j_plot = j_plot + 1
!comm         if( l .lt. num_psd_tht_bins )                   &!&
!comm           write(k_unit,"(2i4,1p48e15.7)") k, j_plot,    &!& ! fort.80**
!comm                         ct_bounds(l+1),                 &!& ! 1
!comm                         theta_next,                     &!& ! 2
!comm                   log10(theta_next),                    &!& ! 3
!comm                    real(i),                             &!& ! 4
!comm                         ct_sk_xw(l,i)/sum_ct_sk_xw          ! 5
!comm   
!comm       enddo
!comm       
!comm       call print_plot_vals(k_unit)
!comm       
!comm     enddo
!comm     
!comm     ! Next, histogram for *all* momenta in the shock frame
!comm     sum_ct_sk_xw = sum( ct_sk_xw )
!comm     if( sum_ct_sk_xw .gt. 0.d0 ) then  ! skip grid locations w/no CRs
!comm       k_unit = k + 8000 + 100*(numion-1); j_plot = 0
!comm       do l = 0, num_psd_tht_bins
!comm         j_plot = j_plot + 1
!comm         
!comm         if(ct_bounds(l) .eq. 1.d0) then
!comm           theta = 1.d-99
!comm         else
!comm           theta = acos(ct_bounds(l))
!comm         endif
!comm         theta_next = acos(ct_bounds(l+1))
!comm         
!comm         write(k_unit,"(2i4,1p48e15.7)") k, j_plot,      &!& ! fort.80**
!comm                       ct_bounds(l),                     &!& ! 1
!comm                       theta,                            &!& ! 2
!comm                 log10(theta),                           &!& ! 3
!comm                   sum(ct_sk_xw(l,:))/sum_ct_sk_xw           ! 4
!comm         
!comm         j_plot = j_plot + 1
!comm         if( l .lt. num_psd_tht_bins )                   &!&
!comm           write(k_unit,"(2i4,1p48e15.7)") k, j_plot,    &!& ! fort.80**
!comm                         ct_bounds(l+1),                 &!& ! 1
!comm                         theta_next,                     &!& ! 2
!comm                   log10(theta_next),                    &!& ! 3
!comm                     sum(ct_sk_xw(l,:))/sum_ct_sk_xw         ! 4
!comm   
!comm       enddo
!comm       
!comm       call print_plot_vals(k_unit)
!comm     endif
!comm     
!comm     
!comm     ! Now histograms of individual momentum decades, in plasma frame
!comm     do i = i_ct_ptpf_min, i_ct_ptpf_max
!comm       
!comm       do l = 0, num_psd_tht_bins
!comm         if(ct_pf_xw(l,i) .gt. 0.d0)                                 &!&
!comm            ct_pf_xw(l,i) = ct_pt_xw(l,i)                            &!&
!comm                           / (ct_bounds(l) - ct_bounds(l+1))
!comm       enddo
!comm       sum_ct_pf_xw = sum( ct_pf_xw(:,i) )
!comm       if( sum_ct_pf_xw .eq. 0.d0 ) cycle   ! skip empty rows
!comm       
!comm       k_unit = k + 9000 + 100*(numion-1); j_plot = 0
!comm       do l = 0, num_psd_tht_bins
!comm         j_plot = j_plot + 1
!comm         
!comm         if(ct_bounds(l) .eq. 1.d0) then
!comm           theta = 1.d-99
!comm         else
!comm           theta = acos(ct_bounds(l))
!comm         endif
!comm         theta_next = acos(ct_bounds(l+1))
!comm         
!comm         write(k_unit,"(2i4,1p48e15.7)") k, j_plot,      &!& ! fort.90**
!comm                       ct_bounds(l),                     &!& ! 1
!comm                       theta,                            &!& ! 2
!comm                 log10(theta),                           &!& ! 3
!comm                  real(i),                               &!& ! 4
!comm                       ct_pf_xw(l,i)/sum_ct_pf_xw            ! 5
!comm         
!comm         j_plot = j_plot + 1
!comm         if( l .lt. num_psd_tht_bins )                   &!&
!comm           write(k_unit,"(2i4,1p48e15.7)") k, j_plot,    &!& ! fort.90**
!comm                         ct_bounds(l+1),                 &!& ! 1
!comm                         theta_next,                     &!& ! 2
!comm                   log10(theta_next),                    &!& ! 3
!comm                    real(i),                             &!& ! 4
!comm                         ct_pf_xw(l,i)/sum_ct_pf_xw          ! 5
!comm   
!comm       enddo
!comm       
!comm       call print_plot_vals(k_unit)
!comm   
!comm     enddo
!comm   
!comm   
!comm     ! Finally, histogram for *all* momenta in the plasma frame
!comm     sum_ct_pf_xw = sum( ct_pf_xw )
!comm     if( sum_ct_pf_xw .gt. 0.d0 ) then  ! skip grid locations w/no CRs
!comm       k_unit = k + 9000 + 100*(numion-1); j_plot = 0
!comm       do l = 0, num_psd_tht_bins
!comm         j_plot = j_plot + 1
!comm         
!comm         if(ct_bounds(l) .eq. 1.d0) then
!comm           theta = 1.d-99
!comm         else
!comm           theta = acos(ct_bounds(l))
!comm         endif
!comm         theta_next = acos(ct_bounds(l+1))
!comm         
!comm         write(k_unit,"(2i4,1p48e15.7)") k, j_plot,      &!& ! fort.90**
!comm                       ct_bounds(l),                     &!& ! 1
!comm                       theta,                            &!& ! 2
!comm                 log10(theta),                           &!& ! 3
!comm                   sum(ct_pf_xw(l,:))/sum_ct_pf_xw           ! 4
!comm         
!comm         j_plot = j_plot + 1
!comm         if( l .lt. num_psd_tht_bins )                   &!&
!comm           write(k_unit,"(2i4,1p48e15.7)") k, j_plot,    &!& ! fort.90**
!comm                         ct_bounds(l+1),                 &!& ! 1
!comm                         theta_next,                     &!& ! 2
!comm                   log10(theta_next),                    &!& ! 3
!comm                     sum(ct_pf_xw(l,:))/sum_ct_pf_xw         ! 4
!comm   
!comm       enddo
!comm       
!comm       call print_plot_vals(k_unit)
!comm     endif
    !------------------------------------------------------------------------
    ! pitch angle histograms plotted
    
    
    ! Deallocate ct arrays in preparation for next iteration of loop
    deallocate( ct_sk_xw )
    deallocate( ct_pf_xw )
    deallocate( ct_ef_xw )
    
  enddo ! loop over grid locations
  
  
  ! Finally, convert dNdp_cr into true dN/dp by dividing by dp for each cell
  do m = 1, 3
    do k = 1, n_grid
      do l = 0, psd_max-1
      
        ! Skip empty cells in PSD
        if( dNdp_cr(l, k, m) .lt. 1.d-66 ) then
          dNdp_cr(l,k,m) = 1.d-99
          cycle
        endif
        
        p_lo_cgs = 10.d0**(psd_mom_bounds(l))
        p_hi_cgs = 10.d0**(psd_mom_bounds(l+1))
        
        dNdp_cr(l, k, m) = dNdp_cr(l, k, m)/(p_hi_cgs - p_lo_cgs)
      enddo
    enddo
  enddo
  
return
end subroutine get_dNdp_cr


!****************************************************************************
!****************************************************************************
subroutine print_dNdp_esc(esc_psd_feb_UpS, esc_psd_feb_DwS)

! Calculate and print dN/dp (normalized to 1) for all particles escaping the
!   shock.  The DwS dN/dp sets the normalization factor, which is reused for
!   the UpS dN/dp.
!
! Input arguments:
!  1) esc_psd_feb_UpS: 2-D array of particles that escaped in the UpS
!    direction from the shock
!  2) esc_psd_feb_DwS: 2-D array of particles that escaped in the DwS
!    direction from the shock
! Output arguments:
!  1) dNdp_cr: 3-D array, containing 1-D array for each grid zone, of dN/dp
!    for particles injected into acceleration process, for three different
!    reference frames

use constants
use parameters, only: psd_max
use controls, only: gam_2, gam_Z, n_ions, n_itrs
use psd_vars, only: num_psd_mom_bins, psd_mom_bounds
use iteration_vars, only: i_itr
use species_vars, only: i_ion

implicit none


  ! Input arguments
real(kind=8), dimension(0:psd_max,0:psd_max), intent(in) ::               &!&
     esc_psd_feb_UpS, esc_psd_feb_DwS
  ! Output arguments

  ! Local variables
integer :: i_approx, j, i, m, i_plot
real(kind=8) :: gam_u, o_o_gam_u, p_lo_cgs, p_hi_cgs
real(kind=8), dimension(0:psd_max,3) :: dNdp_esc_UpS, dNdp_esc_DwS
real(kind=8), dimension(0:psd_max,0:psd_max) :: Xform_corn_pt, Xform_corn_ct


!--! Zero out the output array before summing in loops to follow
   ! Note on dNdp_esc_***'s second dimension:
   !   1  -  Shock frame
   !   2  -  Plasma frame
   !   3  -  ISM frame
  dNdp_esc_UpS(:,:) = 1.d-99
  dNdp_esc_DwS(:,:) = 1.d-99
  
  
!--! Degree of approximation to use in distributing cell_wt if NOT in shock
   !   frame (where uniform distribution may be assumed):
   !    i_approx = 0:  assume uniform distribution
   !    i_approx = 1:  assume isosceles trianglular distribution of cell_wt
   !    i_approx = 2:  assume scalene triangular distribution of cell_wt,
   !                     with peak of triangle centered on mean of ct_hi_pt
   !                     and ct_lo_pt
   !    i_approx = 3:  exact calculation of fractional area, i.e. use no 
   !                     approximations
  i_approx = 2


!--! Handle shock frame separately, as no special treatment is needed to
   !  convert from PSD to dN(p).  Conversion to dN/dp will happen at the
   !  end of this subroutine.
   !-------------------------------------------------------------------------
  do j = 0, psd_max
    do i = 0, psd_max
      if( esc_psd_feb_UpS(i,j) .gt. 1.d-99 )                              &!&
                 dNdp_esc_UpS(i,1) = dNdp_esc_UpS(i,1) + esc_psd_feb_UpS(i,j)
      if( esc_psd_feb_DwS(i,j) .gt. 1.d-99 )                              &!&
                 dNdp_esc_DwS(i,1) = dNdp_esc_DwS(i,1) + esc_psd_feb_DwS(i,j)
    enddo
  enddo
   !-------------------------------------------------------------------------
   ! Shock frame dN(p)s found
  
  
!--! Now handle plasma and ISM frames.  Do this sequentially, since the
   !  process is exactly the same but requires a single input that differs
   !  between the two frames.
   !-------------------------------------------------------------------------
  do m = 2, 3
    
    ! Transform corners of PSD into correct frame for use in finding
    !  p_cell_lo and p_cell_hi in next block of code
    !----------------------------------------------------------------------
    select case( m )
      case( 2 )  !  Plasma frame
        gam_u     = gam_2
        o_o_gam_u = 1.d0 / gam_u
      case( 3 )  !  ISM frame
        gam_u     = gam_Z
        o_o_gam_u = 1.d0 / gam_u
      case default
        write(*,"(A)") "ERROR with frame selection in get_dNdp_cr"
        stop
    end select
    
    call Xform_psd_corners(gam_u, Xform_corn_pt, Xform_corn_ct)
    !----------------------------------------------------------------------
    ! corners transformed
    
    
    ! Note that dNdp_cr as calculated here is dN(p), *not* dN/dp.  That
    !   conversion happens at the end of this subroutine.
    call get_Xform_dN(esc_psd_feb_DwS, m, Xform_corn_pt, Xform_corn_ct,   &!&
       o_o_gam_u, i_approx, dNdp_esc_DwS(:,m) )
    ! Also, the UpS escaping particles escape into the ISM, so the relative
    !   Lorentz factor is identical between the plasma and ISM frames
    if( m .eq. 3) then
      call get_Xform_dN(esc_psd_feb_UpS, m, Xform_corn_pt, Xform_corn_ct, &!&
         o_o_gam_u, i_approx, dNdp_esc_UpS(:,m) )
    endif
    
    
  enddo ! loop over frames
  
  ! Fill plasma-frame UpS dN/dp
  dNdp_esc_UpS(:,2) = dNdp_esc_UpS(:,3)
   !-------------------------------------------------------------------------
   ! Plasma and ISM frame dN(p)s calculated
  
  
!--! Convert dNdp_esc into true dN/dp by dividing by dp for each cell
  do i = 0, num_psd_mom_bins-1
  
    p_lo_cgs = 10.d0**(psd_mom_bounds(i))
    p_hi_cgs = 10.d0**(psd_mom_bounds(i+1))
      
    
    ! Skip empty cells in dN(p)s, otherwise convert
    do m = 1, 3
      if( dNdp_esc_UpS(i,m) .gt. 1.d-99 )                                 &!&
          dNdp_esc_UpS(i,m) = dNdp_esc_UpS(i,m) / (p_hi_cgs - p_lo_cgs)
    
      if( dNdp_esc_DwS(i,m) .gt. 1.d-99 )                                 &!&
          dNdp_esc_DwS(i,m) = dNdp_esc_DwS(i,m) / (p_hi_cgs - p_lo_cgs)
    enddo
  enddo
  
  
!--! Normalize the dN/dps to 1, based on the DwS area
  do m = 1, 3
    where( dNdp_esc_UpS(:,m) .gt. 1.d-99 )                                &!&
           dNdp_esc_UpS(:,m) = dNdp_esc_UpS(:,m) / sum( dNdp_esc_DwS(:,m) )
    
    where( dNdp_esc_DwS(:,m) .gt. 1.d-99 )                                &!&
           dNdp_esc_DwS(:,m) = dNdp_esc_DwS(:,m) / sum( dNdp_esc_DwS(:,m) )
  enddo
  
  
!--! Finally, print the dN/dps to a file
   !-------------------------------------------------------------------------
  if( (i_ion .eq. 1) .and. (i_itr .eq. 1) ) then
    open(unit=527,file="mc_dNdp_esc.dat")
  endif
  
  i_plot = 0
  do i = 0, num_psd_mom_bins
    
    i_plot = i_plot + 1
    write(527,"(2I4,F5.1,20ES13.4E2)") i_itr, i_plot,  &!&
      float(i_ion),                                    &!& ! 1
            psd_mom_bounds(i),                         &!& ! 2 (cgs units)
            psd_mom_bounds(i) - log10(xmp*ccgs),       &!& ! 3 (nat. units)
      log10(dNdp_esc_UpS(i, 1)),                       &!& ! 4  SF, UpS
      log10(dNdp_esc_UpS(i, 2)),                       &!& ! 5  PF, UpS
      log10(dNdp_esc_UpS(i, 3)),                       &!& ! 6  IF, UpS
      log10(dNdp_esc_DwS(i, 1)),                       &!& ! 7  SF, DwS
      log10(dNdp_esc_DwS(i, 2)),                       &!& ! 8  PF, DwS
      log10(dNdp_esc_DwS(i, 3))                            ! 9  IF, DwS
      
    i_plot = i_plot+1
    if(i .lt. num_psd_mom_bins) then
      write(527,"(2I4,F5.1,20ES13.4E2)") i_itr, i_plot,  &!&
        float(i_ion),                                    &!& ! 1
              psd_mom_bounds(i+1),                       &!& ! 2 (cgs)
              psd_mom_bounds(i+1) - log10(xmp*ccgs),     &!& ! 3 (nat.)
        log10(dNdp_esc_UpS(i, 1)),                       &!& ! 4  SF, UpS
        log10(dNdp_esc_UpS(i, 2)),                       &!& ! 5  PF, UpS
        log10(dNdp_esc_UpS(i, 3)),                       &!& ! 6  IF, UpS
        log10(dNdp_esc_DwS(i, 1)),                       &!& ! 7  SF, DwS
        log10(dNdp_esc_DwS(i, 2)),                       &!& ! 8  PF, DwS
        log10(dNdp_esc_DwS(i, 3))                            ! 9  IF, DwS
    endif    
    
  enddo
  
  call print_plot_vals(527)
  
  if( (i_ion .eq. n_ions) .and. (i_itr .eq. n_itrs) ) close(527)
   !-------------------------------------------------------------------------
   ! File printed
  
return
end subroutine print_dNdp_esc


!****************************************************************************
!****************************************************************************
subroutine get_Xform_dN(psd, m, Xform_corn_pt, Xform_corn_ct, o_o_gam_u,  &!&
     i_approx, dN_out)

! Calculate dN(p) (a 1-D array) for the passed slice of PSD in the specified
!   inertial frame.
!
! Input arguments:
!  1) psd: the (2-D) slice of the larger phase space distribution
!  2) m: integer specifying frame into which we're transforming
!  3) Xform_corn_**: array holding transformed corner values, both ptot and
!    cos(theta)
!  4) o_o_gam_u: conversion factor from flux to number density
!  5) i_approx: degree of approximation to use in computing the dNdp_out
! Output arguments:
!  1) dN_out: dN(p) for the given slice of PSD once transformed into the
!    specified frame

use parameters, only: psd_max
use psd_vars, only: num_psd_mom_bins, psd_mom_bounds

implicit none

  ! Input arguments
integer, intent(in) :: m, i_approx
real(kind=8), intent(in) :: o_o_gam_u
real(kind=8), dimension(0:psd_max,0:psd_max), intent(in) :: psd,          &!&
     Xform_corn_pt, Xform_corn_ct
  ! Output arguments
real(kind=8), dimension(0:psd_max), intent(out) :: dN_out
  
  ! Local variables
integer :: j, i, pt_lo_tied, pt_hi_tied, l, l_lo, l_hi
real(kind=8) :: cell_wt, pt_lo_pt, pt_lo_ct, pt_hi_pt, pt_hi_ct, ct_lo_pt,&!&
     ct_lo_ct, ct_hi_pt, ct_hi_ct, p_cell_lo, p_cell_hi, o_o_p_length_tot,&!&
     p_bottom, frac_p_length, ct_height, p_peak, p_denom_lo, p_denom_hi,  &!&
     fractional_area, p_base, ct_rh_height, ct_lh_height, partial_area,   &!&
     missing_area


!--! Loop over cos(theta) and ptot space to re-bin input PSD slice
   !-------------------------------------------------------------------------
  do j = 0, psd_max
    
    do i = 0, psd_max
         
      
      ! Skip empty cells in PSD
      if( psd(i,j) .lt. 1.d-66 ) cycle
      
          
      ! In the below block, "cell_wt" includes n0*u1 normalization and
      !   division by |vx| to change #/(cm^2-s) to #/cm^3. Divide by gamma of
      !   flow speed as part of Lorentz transformation of phase-space density
      !   (Rybicki & Lightman, p.146)
      !
      ! Obtain cell_wt, p_cell_lo and p_cell_hi
      !----------------------------------------------------------------------
      cell_wt = psd(i, j) * o_o_gam_u
      call identify_corners(i, j, Xform_corn_pt, Xform_corn_ct, pt_lo_pt, &!&
          pt_lo_ct, pt_hi_pt, pt_hi_ct, ct_lo_pt, ct_lo_ct, ct_hi_pt,     &!&
          ct_hi_ct, m, pt_lo_tied, pt_hi_tied) 
          ! Use of m in argument list signals what frame we're in; not
          !   used except in case of error/printout

      p_cell_lo = pt_lo_pt
      p_cell_hi = pt_hi_pt
      !----------------------------------------------------------------------
      ! cell_wt, p_cell_lo, p_cell_hi found
      
      
      ! Find lower and upper boundaries of psd_mom_bounds that we'll be
      !   dealing with based on p_cell_lo and p_cell_hi
      !----------------------------------------------------------------------
      do l = 0, num_psd_mom_bins+1
        if(p_cell_lo .lt. psd_mom_bounds(l)) then
          l_lo = l - 1
          exit
        endif
      enddo
      ! Error check to make sure that Lorentz transformations give reasonable
      !   results for l_lo
      if( p_cell_lo .gt. psd_mom_bounds(num_psd_mom_bins+1) ) then
        write(*,"(2A,3I4,2ES14.5E2)") "  In get_dNdp_cr, p_cell_lo ",     &!&
           "> psd_mom_max!  ", m, j, i, p_cell_lo,                        &!&
           psd_mom_bounds(num_psd_mom_bins+1)
        l_lo = num_psd_mom_bins
      endif
      
      
      do l = l_lo, num_psd_mom_bins+1
        if(p_cell_hi .lt. psd_mom_bounds(l)) then
          l_hi = l
          exit
        endif            
      enddo
      ! Error check to make sure that Lorentz transformations give reasonable
      !   results for l_hi
      if( p_cell_hi .gt. psd_mom_bounds(num_psd_mom_bins+1) ) then
        write(*,"(2A,3I4,2ES14.5E2)") "  In get_dNdp_cr, p_cell_hi ",     &!&
           "> psd_mom_max!  ", m, j, i, p_cell_hi,                        &!&
           psd_mom_bounds(num_psd_mom_bins+1)
        l_hi = num_psd_mom_bins
      endif
      !----------------------------------------------------------------------
      ! l_lo and l_hi found
      
      
      !----------------------------------------------------------------------
      ! Distribute the xwt value of cell_wt into the appropriate bins of
      !   dN_out, assigning fractional weights where the cell of PSD covers
      !   more than one bin of psd_mom_bounds
      ! NOTE: virtually every line referring to psd_mom_bounds uses the top
      !   edge of the bin under consideration (i.e. the value of the next bin
      !   of psd_mom_bounds) because this code block originally used xph
      !----------------------------------------------------------------------
      
      ! If assuming uniform distribution of cell_wt between p_cell_lo and
      !   p_cell_hi
      !----------------------------------------------------------------------
      if( i_approx .eq. 0 ) then
        
        o_o_p_length_tot = 1.d0 / (p_cell_hi - p_cell_lo)
        p_bottom         = p_cell_lo
        
        do l = l_lo, l_hi
          ! Cell fits entirely within bin of psd_mom_bounds
          if(p_cell_hi .lt. psd_mom_bounds(l_lo+1)) then
            dN_out(l) = dN_out(l) + cell_wt
            exit
          endif
          
          ! Top of bin in psd_mom_bounds less than p_cell_hi
          if(psd_mom_bounds(l+1) .lt. p_cell_hi) then
            frac_p_length = (psd_mom_bounds(l+1) - p_bottom)              &!&
                           *  o_o_p_length_tot
            dN_out(l) = dN_out(l) + cell_wt*frac_p_length
            ! Adjust p_bottom to mark counting of current bin
            p_bottom = psd_mom_bounds(l+1)
          endif
           
          ! Top of bin in psd_mom_bounds is equal to/greater than
          !   p_cell_hi
          if(psd_mom_bounds(l+1) .ge. p_cell_hi) then
            frac_p_length = (p_cell_hi - psd_mom_bounds(l))               &!&
                           *  o_o_p_length_tot
            dN_out(l) = dN_out(l) + cell_wt*frac_p_length
            exit
          endif
        enddo
      !----------------------------------------------------------------------
      ! i_approx = 0 finished
      
      
      ! If assuming isosceles trianglular distribution of cell_wt, peak is
      !   located above geometric mean of p_cell_lo and p_cell_hi.
      ! If assuming scalene triangular distribution of cell_wt, peak of
      !   triangle is located on (geometric) mean of ct_hi_pt and ct_lo_pt,
      !   noting that both are logarithms.
      ! Only difference is location of p_peak, so handle both with same block
      !   of code.
      !----------------------------------------------------------------------
      else if( (i_approx .eq. 1) .or. (i_approx .eq. 2) ) then
        
        o_o_p_length_tot = 1.d0 / (p_cell_hi - p_cell_lo)
        ct_height    = 2.d0 * cell_wt * o_o_p_length_tot ! A = 1/2*b*h
        
        p_bottom         = p_cell_lo
        if( i_approx .eq. 1 ) then
          p_peak     = 0.5d0 * ( p_cell_lo + p_cell_hi )
        else ! i_approx .eq. 2
          p_peak     = 0.5d0 * ( ct_lo_pt + ct_hi_pt )
        endif
        p_denom_lo   = 1.d0 / (p_peak - p_cell_lo)
        p_denom_hi   = 1.d0 / (p_cell_hi - p_peak)
        
        
        fractional_area = 0.d0 ! Total amount of cell_wt accounted for
        
        
        do l = l_lo, l_hi
        
          ! Cell fits entirely within bin of psd_mom_bounds
          !------------------------------------------------------------------
          if(p_cell_hi .lt. psd_mom_bounds(l_lo+1)) then
            dN_out(l) = dN_out(l) + cell_wt
            exit
          endif
          !------------------------------------------------------------------
          ! cell fits entirely within bin of psd_mom_bounds
          
          
          ! Top of bin in psd_mom_bounds less than or equal to p_peak.
          ! Area is
          !   (1) a triangle if p_bottom = p_cell_lo, or 
          !   (2) a trapezoid if p_bottom > p_cell_lo
          !------------------------------------------------------------------
          if( psd_mom_bounds(l+1) .le. p_peak ) then
            ! Calculate current base
            p_base = psd_mom_bounds(l+1) - p_bottom
            
            ! Calculate right-hand height
            ct_rh_height = (psd_mom_bounds(l+1) - p_cell_lo)              &!&
                          * p_denom_lo * ct_height
            
            ! Calculate left-hand height, taking advantage of right triangle
            !   similarity
            if( p_bottom .eq. p_cell_lo ) then
              ct_lh_height = 0.d0
            else
              ct_lh_height = (p_bottom - p_cell_lo) * p_denom_lo * ct_height
            endif
            
            ! Calculate partial area...
            partial_area = 0.5d0 * p_base * (ct_lh_height + ct_rh_height)
            
            ! ...and add it to dN_out
            dN_out(l) = dN_out(l) + partial_area
            
            ! Adjust p_bottom to mark counting of current bin, update
            !   fractional_area...
            p_bottom = psd_mom_bounds(l+1)
            fractional_area = fractional_area + partial_area
            
            ! ...and move on to next bin
            cycle
          endif
          !------------------------------------------------------------------
          ! top of bin in psd_mom_bounds between p_cell_lo and p_peak
          
          
          ! Top of bin in psd_mom_bounds between p_peak and p_cell_hi; area
          !   is remaining area minus missing triangle at right
          !------------------------------------------------------------------
          if( psd_mom_bounds(l+1) .lt. p_cell_hi ) then
            ! Calculate missing base
            p_base = p_cell_hi - psd_mom_bounds(l+1)
            
            ! Calculate left-hand height of missing area, taking advantage
            !   of right triangle similarity
            ct_lh_height = p_base * p_denom_hi * ct_height
            
            ! Calculate missing area
            missing_area = 0.5d0 * p_base * ct_lh_height
            
            ! Calculate partial area...
            partial_area = (cell_wt - fractional_area) - missing_area
            
            ! ...and add it to dN_out
            dN_out(l) = dN_out(l) + partial_area
            
            ! Adjust p_bottom to mark counting of current bin, update
            !   fractional_area
            p_bottom = psd_mom_bounds(l+1)
            fractional_area = fractional_area + partial_area
            
            ! ...and move on to next bin
            cycle
          endif
          !------------------------------------------------------------------
          ! top of bin in psd_mom_bounds between p_peak and p_cell_hi
          
          
          ! Top of bin in psd_mom_bounds above p_cell_hi; area is remaining
          !   fraction of total
          !------------------------------------------------------------------
          if( psd_mom_bounds(l+1) .ge. p_cell_hi ) then
            ! Calculate partial area...
            partial_area = cell_wt - fractional_area
            
            ! ...and add it to dN_out
            dN_out(l) = dN_out(l) + partial_area
            
            exit
          endif
          !------------------------------------------------------------------
          ! top of bin in psd_mom_bounds above p_cell_hi
        enddo
      !----------------------------------------------------------------------
      ! i_approx = 1, i_approx = 2 finished
      
      
      ! If not assuming anything about shape, i.e. performing exact
      !  calculation of fractional areas
      !----------------------------------------------------------------------
      else if( i_approx .eq. 3 ) then
        
        write(*,"(A)") "ERROR: i_approx = 3 not currently enabled"
        write(*,"(A)") "Stopping program"
        stop
        
        ! Determine which of ct_lo and ct_hi is peak_left and peak_right
        ! Find heights of both peaks, i.e. vertical distance btwn, e.g.,
        !   ct_lo and line segment connecting ct_hi and pt_**
        ! Run through same process as in original version of subroutine,
        !   dividing shape into zones for determining partial areas. Will be
        !   easier this time, though, since bottom line is flat
      !----------------------------------------------------------------------
      ! i_approx = 3 finished
      
      
      ! Invalid selection of i_approx. Flag error and stop program
      else
      
        write(*,"(A)") "ERROR: i_approx must be 0, 1, 2, or 3"
        write(*,"(A,I3)") "Provided value of i_approx: ",i_approx
        write(*,"(A)") "Stopping program"
        stop
        
      endif
      !----------------------------------------------------------------------
      ! cell_wt distributed
      
      
      ! This block tracks the distribution of pitch angles in the shock and
      !   plasma frames.
      ! In the interest of speed, it uses only the equal-weights method
      !   (i_approx = 0 above).
      ! Because ct_bounds DECREASES as index increases (since it's derived
      !   from increasing theta), many things are *slightly* different from
      !   the i_approx = 0 case for total momentum.
      ! WARNING: this must be re-checked and re-tested since it was brought
      !   into new version of code
      !----------------------------------------------------------------------
!comm       ! Adjust cell_wt to remove the velocity-weighting applied in
!comm       !   PSD
!comm       if( i .eq. 0 ) then
!comm         ptsk_cgs = 0.d0
!comm         i_pt_sk  = i_ct_ptsk_min
!comm       else
!comm         ptsk_cgs = 10.d0**(psd_mom_bounds(i) + psd_mom_bounds(i+1))
!comm         ptsk_cgs = sqrt(ptsk_cgs) * xmp * u1sk_cm
!comm         i_pt_sk  = floor( psd_mom_bounds(i) )
!comm       endif
!comm       gam_ptsk = sqrt( 1.d0 + (ptsk_cgs/(rest_mass*ccgs))**2 )
!comm       
!comm       cell_xwt = cell_wt                                        &!&
!comm                 * ptsk_cgs / (gam_ptsk * rest_mass)             &!&
!comm                 / proton_num_den_UpS
!comm       
!comm       
!comm       ! Binning shock frame values very easy; just add directly to
!comm       !  correct bin of histogram
!comm       if( m .eq. 1 )                                            &!&
!comm         ct_sk_xw(j,i_pt_sk) = ct_sk_xw(j,i_pt_sk)  +  cell_xwt
!comm       
!comm       
!comm       ! Identify the min and max extent of the plasma frame cell,
!comm       !   as well as the correct decade of momentum for binning
!comm       ct_cell_lo = min(pt_lo_ct, pt_hi_ct, ct_lo_ct, ct_hi_ct)
!comm       ct_cell_hi = max(pt_lo_ct, pt_hi_ct, ct_lo_ct, ct_hi_ct)
!comm       i_pt_pf = floor( 0.25d0*( pt_lo_pt + pt_hi_pt             &!&
!comm                                +  ct_lo_pt + ct_hi_pt) )
!comm       
!comm       
!comm       ! Determine the spread in cos(theta)
!comm       ! Remember that ct_bounds counts DOWN from +1 to -1 !!!!
!comm       do l = 0, num_psd_tht_bins+1
!comm         if(ct_cell_hi .gt. ct_bounds(l)) then
!comm           l_hi = l - 1
!comm           exit
!comm         endif
!comm       enddo
!comm       
!comm       do l = l_hi, num_psd_tht_bins+1
!comm         if(ct_cell_lo .ge. ct_bounds(l)) then
!comm           l_lo = l
!comm           exit
!comm         endif         
!comm       enddo
!comm       
!comm       
!comm       ! Distribute cell weight among all bins crossed by cell
!comm       ct_length_tot = ct_cell_hi - ct_cell_lo
!comm       ct_bottom     = ct_cell_lo
!comm       
!comm       do l = l_lo, l_hi, -1
!comm         ! Cell fits entirely within bin of ct_bounds
!comm         if(ct_cell_hi .lt. ct_bounds(l_lo-1)) then
!comm           select case( m )
!comm             case( 2 )
!comm               ct_pf_xw(l-1, i_pt_pf) = ct_pf_xw(l-1, i_pt_pf)   &!&
!comm                                       +  cell_xwt
!comm             case( 3 )
!comm               ct_ef_xw(l-1, i_pt_pf) = ct_ef_xw(l-1, i_pt_pf)   &!&
!comm                                       +  cell_xwt
!comm           end select
!comm           exit
!comm         endif
!comm         
!comm         ! Top of bin in ct_bounds less than ct_cell_hi
!comm         if(ct_bounds(l-1) .lt. ct_cell_hi) then
!comm           frac_ct_length = (ct_bounds(l-1) - ct_bottom)         &!&
!comm                           /  ct_length_tot
!comm            
!comm           select case( m )
!comm             case( 2 )
!comm               ct_pf_xw(l-1, i_pt_pf) = ct_pf_xw(l-1, i_pt_pf)   &!&
!comm                                        +  cell_xwt*frac_ct_length
!comm             case( 3 )
!comm               ct_ef_xw(l-1, i_pt_pf) = ct_ef_xw(l-1, i_pt_pf)   &!&
!comm                                       +  cell_xwt*frac_ct_length
!comm           end select
!comm           
!comm           ! Adjust ct_bottom to mark counting of current bin
!comm           ct_bottom = ct_bounds(l-1)
!comm         endif
!comm          
!comm         ! Top of bin in ct_bounds is .ge. ct_cell_hi
!comm         if(ct_bounds(l-1) .ge. ct_cell_hi) then
!comm           frac_ct_length = (ct_cell_hi - ct_bounds(l))          &!&
!comm                           /  ct_length_tot
!comm           
!comm           select case( m )
!comm             case( 2 )
!comm               ct_pf_xw(l-1, i_pt_pf) = ct_pf_xw(l-1, i_pt_pf)   &!&
!comm                                       +  cell_xwt*frac_ct_length
!comm             case( 3 )
!comm               ct_ef_xw(l-1, i_pt_pf) = ct_ef_xw(l-1, i_pt_pf)   &!&
!comm                                       +  cell_xwt*frac_ct_length
!comm           end select
!comm           
!comm           exit
!comm         endif
!comm       enddo ! loop over bins of ct_bounds      
      !------------------------------------------------------------------
      ! pitch angles tracked
      
    enddo ! loop over momenta
  enddo ! loop over angles

return
end subroutine get_Xform_dN


!****************************************************************************
!****************************************************************************
subroutine get_dNdp_therm(num_hist_bins, dNdp_therm, dNdp_therm_pvals,    &!&
     nc_unit, do_therm_pt_formatted)

! Calculate dN/dp (NOT normalized) for particles that haven't been injected
!   into acceleration process.
!
! First create arrays to hold crossing information from this ion species, and
!   fill them with data from the scratch file.
! Next, loop over grid locations, performing three main tasks:
!   a) Transform shock frame values into plasma frame
!   b) Seek out maximum and minimum momenta in plasma frame
!   c) Bin crossings according to plasma frame total momentum
! Also perform the above process for the ISM frame, and store shock frame
!   values for comparison against earlier results.
!
! Inputs:
!  1) num_hist_bins: number of bins to use in histograms
!  2) nc_unit: unit number of scratch file for crossing data
!  3) do_therm_pt_formatted: whether scratch file is human-readable
! Outputs:
!  1) dNdp_therm: 3-D array, containing 1-D array for each grid zone of
!    dN/dp for the thermal particles
!  2) dNdp_therm_pvals: array of momentum bin boundaries for each row of
!    dNdp_therm; each row handled separately to maximize resolution of what
!    may be an extremely narrow peak at radically different energy from the
!    upstream population

use constants
use parameters, only: psd_max, na_g, na_p, na_cr
use controls, only: aa_ion, gam_Z, beta_Z
use grid_vars, only: n_grid, gam_sf_grid
use species_vars, only: i_ion, therm_grid, therm_pxsk, therm_ptsk,        &!&
     therm_xwt, num_crossings, n_cr_count

implicit none

  ! Input variables
integer, intent(in) :: num_hist_bins, nc_unit
logical, intent(in) :: do_therm_pt_formatted
  ! Output variables
real(kind=8), dimension(0:psd_max,na_g,3), intent(out) :: dNdp_therm,     &!&
     dNdp_therm_pvals

  ! Local variables
integer :: max_cross, ntot_crossings, i, i_grid, idum, j, k
integer, dimension(na_g) :: n_cross_fill
real(kind=8) :: rest_mass_en, pxsk, ptsk, cell_wt, gam_u, beta_u,         &!&
     ctsk_min, ctsk_max, ptsk_min, ptsk_max, pxsk_cgs, ptsk_cgs, ctsk,    &!&
     etot_sk_cgs, pxpf_cgs, ptpf_cgs, pxef_cgs, ptef_cgs, del_ptsk,       &!&
     ptpf_max, ptpf_min, del_ptpf, ptef_max, ptef_min, del_ptef, del_ctsk,&!&
     ptsk_lo, ptpf_lo, ptef_lo, ptsk_hi, ptpf_hi, ptef_hi, ptsk_tot,      &!&
     ptpf_tot, ptef_tot, sum_ct_sk_wt, sum_ct_pf_wt, sum_ct_ef_wt
real(kind=8), dimension(:), allocatable :: ptsk_bins, ptsk_vals,          &!&
     ptpf_bins, ptpf_vals, ptef_bins, ptef_vals, ct_sk_bins, ct_sk_wt,    &!&
     ct_pf_bins, ct_pf_wt, ct_ef_bins, ct_ef_wt, ptpf, ctpf, ptef, ctef
real(kind=8), dimension(:,:), allocatable :: therm_px, therm_pt, therm_wt
  ! Variables for Maxwell-Boltzmann fitting
!comm integer :: num_skipped, num_bins, i_unit
!comm real(kind=8) :: sum_psq, sum_psq_f, sum_psd_lnpsq, sum_pfth, sum_f, &!&
!comm      sum_lnpsq, f, psq, lnA, expfac, temp, n0, pressure, gam_ptpf,  &!&
!comm      pressure_tmp
!comm real(kind=8), dimension(:), allocatable :: mb_vals


!--! Set a constant to be used repeatedly
  rest_mass_en = aa_ion(i_ion) * rm_prot  ! Rest mass-energy of the
                                          !  current particle species


!--! Also "zero" out the two output arrays to prevent issues later
  dNdp_therm(:,:,:)       = 1.d-99
  dNdp_therm_pvals(:,:,:) = 1.d-99


!--! Allocate histogram arrays to be used during subroutine
  allocate( ptsk_bins(num_hist_bins) ) ; allocate( ptsk_vals(num_hist_bins) )
  allocate( ptpf_bins(num_hist_bins) ) ; allocate( ptpf_vals(num_hist_bins) )
  allocate( ptef_bins(num_hist_bins) ) ; allocate( ptef_vals(num_hist_bins) )
  allocate( ct_sk_bins(num_hist_bins+1) )
  allocate( ct_pf_bins(num_hist_bins+1) )
  allocate( ct_ef_bins(num_hist_bins+1) )
  allocate( ct_sk_wt(num_hist_bins) )
  allocate( ct_pf_wt(num_hist_bins) )
  allocate( ct_ef_wt(num_hist_bins) )


!--! Create the arrays to hold the crossing information, and fill them with
   !   data from the scratch file
   !-------------------------------------------------------------------------
  max_cross = maxval(num_crossings(1:n_grid),1)
  allocate( therm_px(max_cross,n_grid) )
  allocate( therm_pt(max_cross,n_grid) )
  allocate( therm_wt(max_cross,n_grid) )
  
  
  ! Fill the arrays with data from the crossing arrays and (if needed) from
  !   the scratch file
  ntot_crossings = sum(num_crossings,1)
  n_cross_fill(:) = 0 ! n_cross_fill needs to be an array because we will be
                      !  skipping around the grid as we move through the data
                      !  rather than filling a single zone at a time
  
  rewind(nc_unit)
    
  ! Handle crossings stored within the crossing arrays
  do i = 1, n_cr_count
    i_grid = therm_grid(i)
    
    n_cross_fill(i_grid) = n_cross_fill(i_grid) + 1
    
    ! Note: ordering of coordinates chosen to make memory accesses in next
    !  loop faster
    therm_px(n_cross_fill(i_grid), i_grid) = therm_pxsk(i)
    therm_pt(n_cross_fill(i_grid), i_grid) = therm_ptsk(i)
    therm_wt(n_cross_fill(i_grid), i_grid) = therm_xwt(i)
  enddo
  
  if( ntot_crossings .gt. na_cr ) then
    
    ! Need to go into the scratch file for the remainder of the crossings
    do i = 1, ntot_crossings-na_cr
      if( do_therm_pt_formatted ) then
        read(nc_unit,"(I3,I5,4ES20.12E2)") i_grid, idum, pxsk, ptsk, cell_wt
      else
        read(nc_unit) i_grid, idum, pxsk, ptsk, cell_wt
      endif
      
      n_cross_fill(i_grid) = n_cross_fill(i_grid) + 1
      
      ! Note: ordering of coordinates chosen to make memory accesses in next
      !  loop faster
      therm_px(n_cross_fill(i_grid), i_grid) = pxsk
      therm_pt(n_cross_fill(i_grid), i_grid) = ptsk
      therm_wt(n_cross_fill(i_grid), i_grid) = cell_wt
    enddo
    
  endif
  
   !-------------------------------------------------------------------------
   ! Arrays created and read in
  
  
!--! Main loop over grid locations
   !-------------------------------------------------------------------------
  do i = 1, na_g
    if( num_crossings(i) .eq. 0 ) then
      dNdp_therm(:,i,:) = 1.d-99
      cycle      ! Ignore zones that no thermal particles crossed
    endif
  
    ! (1) Transform crossing data from shock frame into plasma & ISM frames
    !------------------------------------------------------------------------
    ! Get current flow speed
    gam_u  = gam_sf_grid(i)
    beta_u = sqrt( 1.d0 - 1.d0/gam_u**2 )
  
    ! Create arrays to hold plasma frame and ISM frame values, then
    !   initialize them. Since the arrays are handled independently for each
    !   grid zone, all positions in the arrays should be filled with correct
    !   data.  However, initialzing them to non-physical values serves as an
    !   additional check against error.
    allocate( ptpf(num_crossings(i)) )
    allocate( ctpf(num_crossings(i)) )
    allocate( ptef(num_crossings(i)) )
    allocate( ctef(num_crossings(i)) )
    ptpf(:) = -1.d0
    ctpf(:) = -2.d0
    ptef(:) = -1.d0
    ctef(:) = -2.d0
    
    ! Set up min. and max. values for total momentum and cos(theta) in the
    !  shock frame
    ctsk_min = 2.d0
    ctsk_max = -2.d0
    ptsk_min = 1.d99
    ptsk_max = -1.d99
    
    ! Convert information about shock frame values into plasma and ISM
    !   frames.  Even though only total momentum is needed for calculating
    !   dN/dp, histogram of pitch angle values can provide an additional
    !   check on whether the distribution is isotropic in the other frames.
    !   So convert ptsk & pxsk into ptpf/ctpf and ptef/ctef.
    ! Also, find minimum and maximum values of ptsk and ctsk.
    do j = 1, num_crossings(i)
      pxsk_cgs = therm_px(j, i)
      ptsk_cgs = therm_pt(j, i)
      
      ctsk = pxsk_cgs / ptsk_cgs
      if( ctsk .gt. ctsk_max ) ctsk_max = ctsk
      if( ctsk .lt. ctsk_min ) ctsk_min = ctsk
      
      if( ptsk_cgs .gt. ptsk_max ) ptsk_max = ptsk_cgs
      if( ptsk_cgs .lt. ptsk_min ) ptsk_min = ptsk_cgs
      
      etot_sk_cgs = sqrt( (ptsk_cgs * ccgs)**2  +  rest_mass_en**2 )
      
      pxpf_cgs = gam_u * (pxsk_cgs - beta_u*etot_sk_cgs/ccgs)
      ptpf_cgs = sqrt(ptsk_cgs**2 - pxsk_cgs**2 + pxpf_cgs**2)
      
      pxef_cgs = gam_Z * (pxsk_cgs - beta_Z*etot_sk_cgs/ccgs)
      ptef_cgs = sqrt(ptsk_cgs**2 - pxsk_cgs**2 + pxef_cgs**2)
      
      ptpf(j) = ptpf_cgs
      ctpf(j) = pxpf_cgs / ptpf_cgs
      
      ptef(j) = ptef_cgs
      ctef(j) = pxef_cgs / ptef_cgs
    enddo
    
    ! Error checks
    if( count( ptpf .eq. -1.d0 ) .gt. 0 ) then
      write(*,"(A)") "ERROR: unfilled zones in ptpf in subr. get_dNdp_therm"
      write(*,"(A)") "Stopping program"
      stop
    endif
    if( count( ctpf .eq. -2.d0 ) .gt. 0 ) then
      write(*,"(A)") "ERROR: unfilled zones in ctpf in subr. get_dNdp_therm"
      write(*,"(A)") "Stopping program"
      stop
    endif
    if( count( ptef .eq. -1.d0 ) .gt. 0 ) then
      write(*,"(A)") "ERROR: unfilled zones in ptef in subr. get_dNdp_therm"
      write(*,"(A)") "Stopping program"
      stop
    endif
    if( count( ctef .eq. -2.d0 ) .gt. 0 ) then
      write(*,"(A)") "ERROR: unfilled zones in ctef in subr. get_dNdp_therm"
      write(*,"(A)") "Stopping program"
      stop
    endif
    !------------------------------------------------------------------------
    ! Plasma frame, ISM frame values determined; Shock frame extrema found.
    
    
    ! (2) Seek out maximum and minimum total momentum in plasma frame and ISM
    !  frame. Create histogram bins in both momentum and angle as desired.
    !------------------------------------------------------------------------
    ! Already have extrema for shock frame
    del_ptsk = (ptsk_max - ptsk_min) / real(num_hist_bins)
  
    ptpf_max = maxval(ptpf, 1)
    ptpf_min = minval(ptpf, 1)
    del_ptpf = (ptpf_max - ptpf_min) / real(num_hist_bins)
  
    ptef_max = maxval(ptef, 1)
    ptef_min = minval(ptef, 1)
    del_ptef = (ptef_max - ptef_min) / real(num_hist_bins)
    
    del_ctsk = (ctsk_max - ctsk_min) / real(num_hist_bins)
    
    do k = 1, num_hist_bins
      ptsk_bins(k) = ptsk_min + (k-1)*del_ptsk
      ptpf_bins(k) = ptpf_min + (k-1)*del_ptpf
      ptef_bins(k) = ptef_min + (k-1)*del_ptef
      
      ct_sk_bins(k) = ctsk_min + (k-1) * del_ctsk
      ct_pf_bins(k) = -1.d0  +  (k-1) * 2.d0 / real(num_hist_bins)
      ct_ef_bins(k) = -1.d0  +  (k-1) * 2.d0 / real(num_hist_bins)
    enddo
    ct_sk_bins(num_hist_bins+1) = ctsk_max
    ct_pf_bins(num_hist_bins+1) = 1.d0
    ct_ef_bins(num_hist_bins+1) = 1.d0
  
    ptsk_vals(1:num_hist_bins) = 0.d0
    ptpf_vals(1:num_hist_bins) = 0.d0
    ptef_vals(1:num_hist_bins) = 0.d0
    ct_sk_wt(1:num_hist_bins) = 0.d0
    ct_pf_wt(1:num_hist_bins) = 0.d0
    ct_ef_wt(1:num_hist_bins) = 0.d0
    !------------------------------------------------------------------------
    ! Histogram bins set
    
    
    ! (3) Bin particles according to total momentum and cos(theta) in all
    !   three frames
    ! For plasma and ISM frames, divide therm_cwt and by gamma of flow speed
    !   as part of Lorentz transformation of phase-space density (Rybicki &
    !   Lightman, p.146)
    !------------------------------------------------------------------------
    do j = 1, num_crossings(i)
      ! Determine bin in shock frame momentum; add value to correct bin
      k = int((therm_pt(j,i)-ptsk_min)/del_ptsk) + 1
      if( k .gt. num_hist_bins )  k = num_hist_bins
      ptsk_vals(k) = ptsk_vals(k)  +  therm_wt(j,i)
      
      ! Determine bin in plasma frame momentum; add value to correct bin
      k = int((ptpf(j)-ptpf_min)/del_ptpf) + 1
      if( k .gt. num_hist_bins )  k = num_hist_bins
      ptpf_vals(k) = ptpf_vals(k)  +  therm_wt(j,i)/gam_u
      
      ! Determine bin in ISM frame momentum; add value to correct bin
      k = int((ptef(j)-ptef_min)/del_ptef) + 1
      if( k .gt. num_hist_bins )  k = num_hist_bins
      ptef_vals(k) = ptef_vals(k)  +  therm_wt(j,i)/gam_Z
      
      ! Determine shock frame bin of cos(theta); add value to correct bin
      ctsk = therm_px(j,i) / therm_pt(j,i)
      k = int((ctsk - ctsk_min)/del_ctsk) + 1
      if( k .lt. 1 ) k = 1 ! odd floating point error makes some k < 1
      if( k .gt. num_hist_bins )  k = num_hist_bins
      ct_sk_wt(k) = ct_sk_wt(k)  +  therm_wt(j,i)
      
      ! Determine plasma frame bin of cos(theta); add value to correct bin
      k = int((ctpf(j) + 1.d0) * real(num_hist_bins) * 0.5d0)  +  1
      if( k .gt. num_hist_bins )  k = num_hist_bins
      ct_pf_wt(k) = ct_pf_wt(k)  +  therm_wt(j,i)/gam_u
        
      ! Determine ISM frame bin of cos(theta); add value to correct bin
      k = int((ctef(j) + 1.d0) * real(num_hist_bins) * 0.5d0)  +  1
      if( k .gt. num_hist_bins )  k = num_hist_bins
      ct_ef_wt(k) = ct_ef_wt(k)  +  therm_wt(j,i)/gam_Z
    enddo
    
    ! Finish off by converting pt**_vals into dN/dp, which requires dividing
    !  by the momentum spread of the bin
    do k = 1, num_hist_bins
      ptsk_lo = ptsk_bins(k)
      ptpf_lo = ptpf_bins(k)
      ptef_lo = ptef_bins(k)
      if( k .eq. num_hist_bins ) then
        ptsk_hi = ptsk_max
        ptpf_hi = ptpf_max
        ptef_hi = ptef_max
      else
        ptsk_hi = ptsk_bins(k+1)
        ptpf_hi = ptpf_bins(k+1)
        ptef_hi = ptef_bins(k+1)
      endif
      
      ptsk_vals(k) = ptsk_vals(k) / (ptsk_hi - ptsk_lo)
      ptpf_vals(k) = ptpf_vals(k) / (ptpf_hi - ptpf_lo)
      ptef_vals(k) = ptef_vals(k) / (ptef_hi - ptef_lo)
    enddo
    
    ptsk_tot = sum(ptsk_vals, 1)
    ptpf_tot = sum(ptpf_vals, 1)
    ptef_tot = sum(ptef_vals, 1)
    sum_ct_sk_wt = sum(ct_sk_wt)
    sum_ct_pf_wt = sum(ct_pf_wt)
    sum_ct_ef_wt = sum(ct_ef_wt)
    !------------------------------------------------------------------------
    ! Particles binned
    
    
    ! Finally, copy pt**_vals and pt**_bins into the dNdp arrays; shift down by
    !  one bin to make both thermal and cosmic ray dN/dp's start at bin 0.
    !  Also, set zero values of dNdp_therm to a small but nonzero number.
    ! Note on third dimension of dNdp_therm:
    !    1  -  Shock frame
    !    2  -  Plasma frame
    !    3  -  ISM frame
    do k = 1, num_hist_bins
      ! Shock frame
      if( ptsk_vals(k) .eq. 0.d0 ) then
        dNdp_therm(k-1, i, 1) = 1.d-99
      else
        dNdp_therm(k-1, i, 1) = ptsk_vals(k)
      endif
      
      ! Plasma frame
      if( ptpf_vals(k) .eq. 0.d0 ) then
        dNdp_therm(k-1, i, 2) = 1.d-99
      else
        dNdp_therm(k-1, i, 2) = ptpf_vals(k)
      endif
      
      ! ISM frame
      if( ptef_vals(k) .eq. 0.d0 ) then
        dNdp_therm(k-1, i, 3) = 1.d-99
      else
        dNdp_therm(k-1, i, 3) = ptef_vals(k)
      endif
      
      ! Set bin boundaries
      dNdp_therm_pvals(k-1, i, 1) = ptsk_bins(k)
      dNdp_therm_pvals(k-1, i, 2) = ptpf_bins(k)
      dNdp_therm_pvals(k-1, i, 3) = ptef_bins(k)
    enddo
    dNdp_therm_pvals(num_hist_bins, i, 1) = ptsk_max
    dNdp_therm_pvals(num_hist_bins, i, 2) = ptpf_max
    dNdp_therm_pvals(num_hist_bins, i, 3) = ptef_max
    
    
    ! The next section is a set of tests and checks on thermal population.
    ! Fit Maxwell-Boltzmann distribution to plasma frame dN/dp, and
    !   determine temperature of distribution. Also, calculate pressure using
    !   integral formulation of equation (17) from Ellison & Reynolds (1991)
    !   [1991ApJ...378..214E].
    !------------------------------------------------------------------------
!comm     allocate( mb_vals(num_hist_bins) )
!comm     num_skipped   = 0
!comm     sum_psq       = 0.d0 ;  sum_psq_f     = 0.d0
!comm     sum_psq_lnpsq = 0.d0 ;  sum_pfth      = 0.d0
!comm     sum_f         = 0.d0 ;  sum_lnpsq     = 0.d0
!comm     do k = 0, num_hist_bins-1 ! now using shifted bins
!comm       ! Here, f is the natural logarithm of dN/dp because fitting occurs
!comm       !  in log-log space
!comm       if( dNdp_therm(k, i, 2) .gt. 0.d0 ) then
!comm         f = log(dNdp_therm(k, i, 2))
!comm       else
!comm         f = 0.d0
!comm         num_skipped = num_skipped + 1
!comm         cycle
!comm       endif
!comm       psq    = dNdp_therm_pvals(k,i,2)**2
!comm       sum_psq       = sum_psq       + psq
!comm       sum_psq_f     = sum_psq_f     + psq * f
!comm       sum_psq_lnpsq = sum_psq_lnpsq + psq * log(psq)
!comm       sum_pfth      = sum_pfth      + psq**2
!comm       sum_f         = sum_f         + f
!comm       sum_lnpsq     = sum_lnpsq     + log(psq)
!comm     enddo
!comm     num_bins = num_hist_bins - num_skipped
!comm     lnA    = ( sum_psq*sum_psq_f - sum_psq*sum_psq_lnpsq      &!&
!comm               - sum_pfth*sum_f + sum_pfth*sum_lnpsq )         &!&
!comm             / (sum_psq**2 - num_bins*sum_pfth)
!comm     expfac = ( -num_bins*sum_psq_f + sum_f*sum_psq            &!&
!comm               - sum_psq*sum_lnpsq + num_bins*sum_psq_lnpsq )  &!&
!comm             / (sum_psq**2 - num_bins*sum_pfth)
!comm     temp   = -0.5d0 / (xmp*aa_ion(numion) * xkb * expfac)
!comm     n0     = exp(lnA) * (xmp*aa_ion(numion)*xkb*temp)**(1.5)        &!&
!comm             *  sqrt(0.5d0 * pii)
!comm       
!comm     
!comm     ! Generate values of fitted M-B distribution
!comm     do k = 0, num_hist_bins-1
!comm       psq        = dNdp_therm_pvals(k,i,2)**2
!comm       mb_vals(k+1) = lnA + log(psq) + expfac*psq
!comm       mb_vals(k+1) = exp(mb_vals(k+1))
!comm     enddo
!comm     
!comm     ! Calculate pressure using integral formula.
!comm     ! Plot dN/dp, fitted M-B distribution, and pressure
!comm     j_plot  = 0 ;  pressure = 0.d0
!comm     i_unit  = 700 + i
!comm     do k = 1, num_hist_bins
!comm       ptpf_cgs = sqrt( dNdp_therm_pvals(k-1,i,2)         &!&
!comm                       * dNdp_therm_pvals(k,i,2) )
!comm       del_ptpf = dNdp_therm_pvals(k,i,2)                 &!&
!comm                 - dNdp_therm_pvals(k-1,i,2)
!comm       
!comm       gam_ptpf = sqrt( 1.d0 + (ptpf_cgs*ccgs/rest_mass_en)**2 )
!comm       pressure_tmp = ptpf_cgs                                 &!&
!comm                     * ptpf_cgs/(xmp*aa_ion(i_ion)*gam_ptpf)   &!&
!comm                     * dNdp_therm(k-1, i, 2)
!comm       pressure_tmp = pressure_tmp * third * del_ptpf
!comm       
!comm       pressure = pressure + pressure_tmp
!comm       
!comm       j_plot = j_plot+1
!comm       write(i_unit,"(2i4,1p48e14.5)") i, j_plot,     &!&
!comm               log10(dNdp_therm_pvals(k-1,i,2)),      &!& ! 1
!comm               dNdp_therm(k-1,i,2)/ptpf_tot,          &!& ! 2
!comm               ct_pf_bins(k),                         &!& ! 3
!comm               ct_pf_wt(k)/sum_ct_pf_wt,              &!& ! 4
!comm               ct_sk_bins(k),                         &!& ! 5
!comm               ct_sk_wt(k)/sum_ct_sk_wt,              &!& ! 6
!comm               mb_vals(k)/sum(mb_vals),               &!& ! 7
!comm               pressure                                   ! 8
!comm       j_plot = j_plot+1
!comm       if(k .lt. num_hist_bins)                       &!&
!comm         write(i_unit,"(2i4,1p48e14.5)") i, j_plot,   &!&
!comm                 log10(dNdp_therm_pvals(k,i,2)),      &!& ! 1
!comm                 dNdp_therm(k-1,i,2)/ptpf_tot,        &!& ! 2
!comm                 ct_pf_bins(k+1),                     &!& ! 3
!comm                 ct_pf_wt(k)/sum_ct_pf_wt,            &!& ! 4
!comm                 ct_sk_bins(k+1),                     &!& ! 5
!comm                 ct_sk_wt(k)/sum_ct_sk_wt,            &!& ! 6
!comm                 mb_vals(k)/sum(mb_vals),             &!& ! 7
!comm                 pressure                                 ! 8
!comm     enddo ! loop over num_hist_bins
!comm     call print_plot_vals(i_unit)
!comm     deallocate( mb_vals )
    !------------------------------------------------------------------------
    ! End of tests/checks section
    
    
    ! Deallocate arrays in preparation for next cycle through main loop
    deallocate( ptpf )
    deallocate( ctpf )
    deallocate( ptef )
    deallocate( ctef )
  enddo ! loop over grid zones
  
  
!--! Deallocate arrays to prevent memory leaks or prepare for next call of
   !   get_dNdp_therm
  deallocate( therm_px ) ; deallocate( therm_pt ) ; deallocate( therm_wt )
  deallocate( ptsk_bins ); deallocate( ptsk_vals )
  deallocate( ptpf_bins ); deallocate( ptpf_vals )
  deallocate( ptef_bins ); deallocate( ptef_vals )
  deallocate( ct_sk_bins )
  deallocate( ct_pf_bins )
  deallocate( ct_ef_bins )
  deallocate( ct_sk_wt ) ; deallocate( ct_pf_wt ) ; deallocate( ct_ef_wt )

return
end subroutine get_dNdp_therm


!****************************************************************************
!****************************************************************************
subroutine identify_corners(i, j, Xform_corn_pt, Xform_corn_ct, pt_lo_pt, &!&
     pt_lo_ct, pt_hi_pt, pt_hi_ct, ct_lo_pt, ct_lo_ct, ct_hi_pt, ct_hi_ct,&!&
     ii_sk_pf, pt_lo_tied, pt_hi_tied)

! Given particular cell in PSD, identify the corners with lowest and highest
!  total momenta and the lower/higher of the two remaining cos(theta) values.
!
! NOTE: subroutine liberally uses "minloc", "maxloc", and "count" intrinsic
!  functions
!
! Inputs:
!  1) i: marker for location on the total momentum axis
!  2) j: marker for location on the cos(theta) axis
!  3) Xform_corn_pt: values of total momentum for the corners of PSD,
!    transformed into desired frame
!  4) Xform_corn_ct: values of cos(theta) for the corners of PSD, transformed
!    into desired frame
!  5) ii_sk_pf: integer giving frame into which corners were transformed
! Outputs:
!  1) pt_lo_pt: momentum of the corner with lowest total momentum
!  2) pt_lo_ct: cos(theta) of the corner with lowest total momentum
!  3) pt_hi_pt: momentum of the corner with highest total momentum
!  4) pt_hi_ct: cos(theta) of the corner with highest total momentum
!  5) ct_lo_pt: momentum of the remaining corner with lower cos(theta)
!  6) ct_lo_ct: cos(theta) of the remaining corner with lower cos(theta)
!  7) ct_hi_pt: momentum of the remaining corner with higher cos(theta)
!  8) ct_hi_ct: cos(theta) of the remaining corner with higher cos(theta)
!  9) pt_lo_tied: flag used if multiple cells are tied for lowest ptot
! 10) pt_hi_tied: flag used if multiple cells are tied for highest ptot

use parameters, only: psd_max

implicit none

  ! Input variables
integer, intent(in) :: i, j, ii_sk_pf
real(kind=8), dimension(0:psd_max, 0:psd_max) :: Xform_corn_pt, Xform_corn_ct
  ! Output variables
integer, intent(out) :: pt_lo_tied, pt_hi_tied
real(kind=8), intent(out) :: pt_lo_pt, pt_lo_ct, pt_hi_pt, pt_hi_ct,      &!&
     ct_lo_pt, ct_lo_ct, ct_hi_pt, ct_hi_ct

  ! Local variables
integer :: i_pt_lo, i_pt_hi, j_ct_hi, j_ct_lo
logical, dimension(4) :: mask
real(kind=8), dimension(4) :: corn_pts, corn_cts


!--! Get momentum and cos(theta) values from the grid
  corn_pts(1) = Xform_corn_pt(i,  j);   corn_cts(1) = Xform_corn_ct(i,  j)
  corn_pts(2) = Xform_corn_pt(i+1,j);   corn_cts(2) = Xform_corn_ct(i+1,j)
  corn_pts(3) = Xform_corn_pt(i,  j+1); corn_cts(3) = Xform_corn_ct(i,  j+1)
  corn_pts(4) = Xform_corn_pt(i+1,j+1); corn_cts(4) = Xform_corn_ct(i+1,j+1)
  
  
!--! Pre-set logical array
  mask = .true. ! array assignment
  
  
!--! Pre-set flags for ties
  pt_lo_tied = 0
  pt_hi_tied = 0


!--! Identify *a* corner with lowest total momentum; if there are no ties,
   !   this corner is pt_lo
   !-------------------------------------------------------------------------
  i_pt_lo  = minloc(corn_pts, 1)
  pt_lo_pt = corn_pts(i_pt_lo)
  pt_lo_ct = corn_cts(i_pt_lo)
  mask(i_pt_lo) = .false.
  
  ! Count number of corners with identified lowest momentum; if this number
  !  isn't 1, flag as needing additional attention 
  if( count( corn_pts .eq. pt_lo_pt ) .gt. 1 ) then
    ! There was a tie, so set flag for special handling later
    write(*,"(A)") "ERROR: multiple corners with lowest momentum"
    write(*,"(3I4,4ES13.4E2)") ii_sk_pf, i, j, corn_pts(1), corn_pts(2),  &!&
        corn_pts(3), corn_pts(4)
    pt_lo_tied = 1
  !  write(*,"(A)") "Stopping program now."
  !  stop
  else if( count( corn_pts .eq. pt_lo_pt ) .lt. 1 ) then
    ! Who knows what happened
    write(*,"(A)") "ERROR: no corner with lowest momentum"
    write(*,"(3I4,ES13.4E2)") ii_sk_pf, i, j, pt_lo_pt
    write(*,"(4ES13.4E2)") corn_pts(1), corn_pts(2), corn_pts(3), corn_pts(4)
    write(*,"(A)") "Stopping program"
    stop
  endif
   !-------------------------------------------------------------------------
   ! pt_lo identified
  
  
!--! Identify *a* corner with highest total momentum; if there are no ties,
   !   this corner is pt_hi
   !-------------------------------------------------------------------------
  i_pt_hi  = maxloc(corn_pts, 1)
  pt_hi_pt = corn_pts(i_pt_hi)
  pt_hi_ct = corn_cts(i_pt_hi)
  mask(i_pt_hi) = .false.
  
  ! Count number of corners with identified highest momentum; if this number
  !   isn't 1, flag as needing additional attention
  if( count( corn_pts .eq. pt_hi_pt ) .gt. 1 ) then
    ! There was a tie, so set flag for special handling later
    write(*,"(A)") "ERROR: multiple corners with highest momentum"
    write(*,"(3I4,4ES13.4E2)") ii_sk_pf, i, j, corn_pts(1), corn_pts(2),  &!&
        corn_pts(3), corn_pts(4)
    pt_hi_tied = 1
  !  write(*,"(A)") "Stopping program now."
  !  stop
  else if( count( corn_pts .eq. pt_hi_pt ) .lt. 1 ) then
    ! Who knows what happened
    write(*,"(A)") "ERROR: no corner with highest momentum"
    write(*,"(3I4,ES13.4E2)") ii_sk_pf, i, j, pt_hi_pt
    write(*,"(4ES13.4E2)") corn_pts(1), corn_pts(2), corn_pts(3), corn_pts(4)
    write(*,"(A)") "Stopping program"
    stop
  endif
   !-------------------------------------------------------------------------
   ! pt_hi identified
  
  
!--! Between two remaining corners, identify higher/lower values of 
   !   cos(theta).  If the two values aren't identical, one corner will be
   !   ct_hi and other will be ct_lo.
   ! NOTE: in event that ct_hi_ct and ct_lo_ct are equal, assign corner with
   !   lower total momentum as ct_lo
   !-------------------------------------------------------------------------
  j_ct_hi  = maxloc(corn_cts, 1, mask)
  ct_hi_pt = corn_pts(j_ct_hi)
  ct_hi_ct = corn_cts(j_ct_hi)
  mask(j_ct_hi) = .false.
  
  j_ct_lo  = minloc(corn_cts, 1, mask)
  ct_lo_pt = corn_pts(j_ct_lo)
  ct_lo_ct = corn_cts(j_ct_lo)
  
  ! Make sure the two corners aren't somehow identical
  if( ct_hi_ct .eq. ct_lo_ct ) then
  
    ! ct_hi and ct_lo have same value of cos(theta); use total momentum as
    !   the tiebreaker
    if( ct_hi_pt .gt. ct_lo_pt ) then
      ! ct_hi remains ct_hi by virtue of its higher momentum
  
    else if ( ct_hi_pt .lt. ct_lo_pt ) then
      ! Switch ct_hi and ct_lo because of momentum ordering
      ct_hi_pt = corn_pts(j_ct_lo)
      ct_hi_ct = corn_cts(j_ct_lo)
      ct_lo_pt = corn_pts(j_ct_hi)
      ct_lo_ct = corn_cts(j_ct_hi)
  
    else
      ! ct_hi and ct_lo have exactly the same values in momentum and in
      !   cos(theta)
      write(*,"(A)") "ERROR: ct_hi and ct_lo identical"
      write(*,"(2I4,4ES13.4E2)") i, j, ct_hi_pt, ct_hi_ct, ct_lo_pt, ct_lo_ct
      write(*,"(A)") "Stopping program"
      stop
    
    endif
  endif
   !-------------------------------------------------------------------------
   ! ct_hi and ct_lo identified
  
  
!--! In the case of a tie for lowest/highest momentum, determine which of
   !   ct_lo and ct_hi is the source of the tie.  Also make sure the pt_**
   !   involved in the tie has the lower cos(theta).
   !-------------------------------------------------------------------------
  if( pt_lo_tied .eq. 1 ) then
  
    ! Determine which of ct_lo and ct_hi is tied with pt_lo
    if( pt_lo_pt .eq. ct_lo_pt ) then
      pt_lo_tied = 2
    else if( pt_lo_pt .eq. ct_hi_pt ) then
      pt_lo_tied = 3
    else
      write(*,"(2A)") "ERROR: something has gone horribly wrong in ",     &!&
          "identify_corners. Code 1"
      write(*,"(A)") "Stopping program now."
      stop
    endif
    
    ! Of the two tied corners, pt_lo should have the lower cos(theta)
    if( pt_lo_tied .eq. 2 ) then
      
      ! Corner tied with pt_lo is ct_lo
      if( pt_lo_ct .gt. ct_lo_ct ) then
        ! Switch pt_lo and ct_lo
        pt_lo_pt = corn_pts(j_ct_lo)
        pt_lo_ct = corn_cts(j_ct_lo)
        ct_lo_pt = corn_pts(i_pt_lo)
        ct_lo_ct = corn_cts(i_pt_lo)
      else if( pt_lo_ct .lt. ct_lo_ct ) then
        ! Nothing to do here. Assignment is correct
      else
        write(*,"(2A)") "ERROR: pt_lo and ct_lo identical"
        write(*,"(2I4,4ES13.4E2)") i, j, pt_lo_pt, pt_lo_ct, ct_lo_pt,    &!&
            ct_lo_ct
        write(*,"(A,I0,2A)") "If ",i," = 0, reduce EMNFC in mc_in.dat"
        write(*,"(A)") "Stopping program now."
        stop
      endif
      
    else if( pt_lo_tied .eq. 3 ) then
      
      ! Corner tied with pt_lo is ct_hi
      if( pt_lo_ct .gt. ct_hi_ct ) then
        ! Switch pt_lo and ct_hi
        pt_lo_pt = corn_pts(j_ct_hi)
        pt_lo_ct = corn_cts(j_ct_hi)
        ct_hi_pt = corn_pts(i_pt_lo)
        ct_hi_ct = corn_cts(i_pt_lo)
      else if( pt_lo_ct .lt. ct_hi_ct ) then
        ! Nothing to do here. Assignment is correct
      else
        write(*,"(2A)") "ERROR: pt_lo and ct_hi identical"
        write(*,"(2I4,4ES13.4E2)") i, j, pt_lo_pt, pt_lo_ct, ct_hi_pt,    &!&
            ct_hi_ct
        write(*,"(A,I0,2A)") "If ",i," = 0, reduce EMNFC in mc_in.dat"
        write(*,"(A)") "Stopping program now."
        stop
      endif
      
    endif ! check on pt_lo having lower momentum
    
  endif
  
  if( pt_hi_tied .eq. 1 ) then
  
    ! Determine which of ct_lo and ct_hi is tied with pt_hi
    if( pt_hi_pt .eq. ct_lo_pt ) then
      pt_hi_tied = 2
    else if( pt_hi_pt .eq. ct_hi_pt ) then
      pt_hi_tied = 3
    else
      write(*,"(2A)") "ERROR: something has gone horribly wrong in ",     &!&
          "identify_corners. Code 2"
      write(*,"(A)") "Stopping program now."
      stop
    endif
    
    ! Of the two tied corners, pt_lo should have the lower cos(theta)
    if( pt_hi_tied .eq. 2 ) then
      
      ! Corner tied with pt_hi is ct_lo
      if( pt_hi_ct .gt. ct_lo_ct ) then
        ! Switch pt_hi and ct_lo
        pt_hi_pt = corn_pts(j_ct_lo)
        pt_hi_ct = corn_cts(j_ct_lo)
        ct_lo_pt = corn_pts(i_pt_hi)
        ct_lo_ct = corn_cts(i_pt_hi)
      else if( pt_hi_ct .lt. ct_lo_ct ) then
        ! Nothing to do here. Assignment is correct
      else
        write(*,"(2A)") "ERROR: pt_hi and ct_lo identical"
        write(*,"(2I4,4ES13.4E2)") i, j, pt_hi_pt, pt_hi_ct, ct_lo_pt,    &!&
            ct_lo_ct
        write(*,"(A,I0,2A)") "If ",i," = 0, reduce EMNFC in mc_in.dat"
        write(*,"(A)") "Stopping program now."
        stop
      endif
      
    else if( pt_hi_tied .eq. 3 ) then
      
      ! Corner tied with pt_hi is ct_hi
      if( pt_hi_ct .gt. ct_hi_ct ) then
        ! Switch pt_hi and ct_hi
        pt_hi_pt = corn_pts(j_ct_hi)
        pt_hi_ct = corn_cts(j_ct_hi)
        ct_hi_pt = corn_pts(i_pt_hi)
        ct_hi_ct = corn_cts(i_pt_hi)
      else if( pt_hi_ct .lt. ct_hi_ct ) then
        ! Nothing to do here. Assignment is correct
      else
        write(*,"(2A)") "ERROR: pt_hi and ct_lo identical"
        write(*,"(2I4,4ES13.4E2)") i, j, pt_hi_pt, pt_hi_ct, ct_hi_pt,    &!&
            ct_hi_ct
        write(*,"(A,I0,2A)") "If ",i," = 0, reduce EMNFC in mc_in.dat"
        write(*,"(A)") "Stopping program now."
        stop
      endif
      
    endif ! check on pt_hi having lower momentum
    
  endif
   !-------------------------------------------------------------------------
   ! Ties handled

return
end subroutine identify_corners


!****************************************************************************
!****************************************************************************
subroutine rebin_dNdp_therm(num_hist_bins, dNdp_therm, dNdp_therm_pvals,  &!&
     num_psd_mom_bins, psd_mom_bounds, dNdp_therm_rebin)

! This subroutine takes the distribution of thermal particles and rebins it
!   in the bins used for the cosmic rays.
! Note that the thermal distribution is provided (and rebinned separately) in
!   all three frames: shock, plasma, and ISM.
!
! Input arguments:
!  1) num_hist_bins: number of bins used for the thermal distribution
!  2) dNdp_therm: distribution of thermal particles in all three frames
!  3) dNdp_therm_pvals: arrays holding bin boundary values for the thermal
!    distributions; as with the distributions themselves, they change based
!    on grid location
!  4) num_psd_mom_bins: number of bins used for the cosmic ray distribution
!  5) psd_mom_bounds: boundaries of the bins for the cosmic ray distribution
!    Unlike the thermal boundaries, these are constant everywhere on the
!    grid.
! Output arguments:
!  1) dNdp_therm_rebin: rebinned histogram of the thermal distribution

use constants
use parameters, only: psd_max

implicit none

  ! Input variables
integer, intent(in) :: num_hist_bins, num_psd_mom_bins
real(kind=8), dimension(0:psd_max), intent(in) :: psd_mom_bounds
real(kind=8), dimension(0:psd_max,3), intent(in) :: dNdp_therm,           &!&
     dNdp_therm_pvals
  ! Output variables
real(kind=8), dimension(0:psd_max,3), intent(out) :: dNdp_therm_rebin

  ! Local variables
integer :: i, m, j, i_lo
real(kind=8) :: dN_remaining, p_lo, p_hi, p_bottom, p_length,             &!&
     o_o_p_length, frac_p_length
real(kind=8), dimension(0:psd_max) :: lin_bounds, dN_therm, therm_pvals,  &!&
     dN_therm_rebin


!--! Zero out the output array
  dNdp_therm_rebin(:,:) = 1.d-99


!--! Convert psd_mom_bounds from log space into linear space
  do i = 0, num_psd_mom_bins
    lin_bounds(i) = 10.d0**psd_mom_bounds(i)
  enddo


!--! Loop over the various frames.
   !   1  -  Shock
   !   2  -  Plasma
   !   3  -  ISM
  do m = 1, 3
  
  !--! Convert the arrays from 2-D down to 1-D to make referring to them
     !   easier.  Also, turn dN/dp into dN, noting that dNdp_therm_pvals is
     !   already in cgs units.
    dN_therm(:)    = dNdp_therm(:,m)
    therm_pvals(:) = dNdp_therm_pvals(:,m)
    dN_therm_rebin(:) = 1.d-99
    do i = 0, num_hist_bins - 1
      dN_therm(i) = dN_therm(i) * (therm_pvals(i+1) - therm_pvals(i))
      if( dN_therm(i) .lt. 1.d-99 ) dN_therm(i) = 1.d-99
    enddo
    
    
  !--! Find the lower bound on where the thermal distribution might fall, to
     !   save cycles later
    do i = 0, num_psd_mom_bins-1
      i_lo = i
      if( lin_bounds(i+1) .gt. therm_pvals(0) ) exit
    enddo
    
    
  !--! Now, re-bin the thermal distribution.  Loop over the bins of dN_therm,
     !   assigning each one to the CR distribution bin(s) where it belongs.
    do j = 0, num_hist_bins-1
      
      if( dN_therm(j) .le. 1.d-99 ) cycle ! skip empty bins
      dN_remaining = dN_therm(j)
     
      p_lo = therm_pvals(j)
      p_hi = therm_pvals(j+1)
      
      p_bottom     = p_lo
      p_length     = p_hi - p_lo
      o_o_p_length = 1.d0 / p_length
      
      do i = i_lo, num_psd_mom_bins-1
        
        ! Determine the fraction of p_length that the current cell
        !   represents.
        frac_p_length = (lin_bounds(i+1) - p_bottom) * o_o_p_length
        
        ! If the top of the cell is above the top of the thermal bin, the
        !   remainder of the thermal bin fits entirely in this cell of the CR
        !   dN/dp.  Exit the loop and move to the next thermal bin.
        if( lin_bounds(i+1) .gt. p_hi ) then
          dN_therm_rebin(i) = dN_therm_rebin(i) + dN_remaining
          exit
        else if( frac_p_length .lt. 0.d0 ) then
          ! The current cell of the CR histogram is *entirely* below the bin
          !  of the thermal distribution.  Skip and move on to the next cycle
          cycle
        endif
        
        ! If the code reached this point, the thermal bin extends beyond the
        !   upper edge of the current cell.  Add the fraction *within* the
        !   cell to dN_therm_rebinned and reset for the next cycle
        dN_therm_rebin(i) = dN_therm_rebin(i)  +  dN_therm(j) * frac_p_length
        dN_remaining      = dN_remaining       -  dN_therm(j) * frac_p_length
        p_bottom = lin_bounds(i+1)
        
      enddo ! loop over cells in CR histogram
      
    enddo ! loop over thermal bins
    
    ! Convert the re-binned thermal distribution from dN into dN/dp.  Then
    !   copy it into the output array.
    do i = 0, num_psd_mom_bins-1
      if( dN_therm_rebin(i) .gt. 1.d-90 )                                 &!&
               dN_therm_rebin(i) = dN_therm_rebin(i)                      &!&
                                  / (lin_bounds(i+1)-lin_bounds(i))
    enddo
    dNdp_therm_rebin(:,m) = dN_therm_rebin(:)
    
  enddo ! loop over frames

return
end subroutine rebin_dNdp_therm


!****************************************************************************
!****************************************************************************
subroutine Xform_psd_corners(gam_in, Xform_corn_pt, Xform_corn_ct)

! Given a relative Lorentz factor between two frames, transform the corners
!   of the PSD into the new frame.
! Outputs total momentum as log of cgs values
!DOLATER: remove conversion to/from log space here and in get_dndp_cr
!
! Inputs:
!  1) gam_in: relative Lorentz factor between shock and resultant frame
! Outputs:
!  1) Xform_corn_pt: total momenta at the corners
!  2) Xform_corn_ct: cos(theta) (NOT theta!!!) values at the corners

use constants
use parameters, only: psd_max
use controls, only: aa_ion, psd_lin_cos_bins
use psd_vars, only: num_psd_tht_bins, psd_tht_bounds, num_psd_mom_bins,   &!&
     psd_mom_bounds
use species_vars, only: i_ion

implicit none

  ! Input variables
real(kind=8), intent(in) :: gam_in
  ! Output variables
real(kind=8), dimension(0:psd_max, 0:psd_max), intent(out) ::             &!&
     Xform_corn_pt, Xform_corn_ct

  ! Local variables
integer :: i, j
real(kind=8) :: rest_mass_en, beta_u, costht, pt_sk_cgs, px_sk_cgs,       &!&
     etot_sk_cgs, px_Xf_cgs, pt_Xf_cgs


!--! Administrative constants
  rest_mass_en = aa_ion(i_ion) * rm_prot
  beta_u = sqrt( 1.d0 - 1.d0/gam_in**2 )
  if( gam_in .lt. 1.000001d0 ) beta_u = 0.d0 ! Prevent floating point issues
  
  
!--! Fill Xform_corn_** arrays, looping over angle outermost
  do j = 0, num_psd_tht_bins + 1
    
    ! Determine current cosine, remembering that psd_tht_bounds has both a
    !   linearly-spaced region in cosine and logarithmically-spaced region in
    !   theta.
    ! Also need to remember that the most finely spaced bins should occur in
    !   the UpS-pointing direction, so need to negate psd_tht_bounds to get
    !   true cosine value.
    if( j .gt. (num_psd_tht_bins - psd_lin_cos_bins) ) then
      costht = -psd_tht_bounds(j)
    else
      costht = -cos( psd_tht_bounds(j) )
    endif
    
    ! With angle fixed, loop over total momenta
    do i = 0, num_psd_mom_bins + 1
      
      ! psd_mom_bounds uses logarithmic spacing for its bins, so undo that
      !   before continuing the calculation
      pt_sk_cgs   = 10.d0**(psd_mom_bounds(i))
!      if( i .eq. 0 ) pt_sk_cgs = 0.d0     ! Edge case when i = 0
      px_sk_cgs   = pt_sk_cgs * costht
      etot_sk_cgs = sqrt((pt_sk_cgs * ccgs)**2 + rest_mass_en**2)
    
      px_Xf_cgs  = gam_in * (px_sk_cgs - beta_u*etot_sk_cgs/ccgs) 
      pt_Xf_cgs  = sqrt(pt_sk_cgs**2 + px_Xf_cgs**2 - px_sk_cgs**2)
      
      ! Transform to log space because get_dNdp_cr expects it
      Xform_corn_pt(i,j) = log10(pt_Xf_cgs)
      Xform_corn_ct(i,j) = px_Xf_cgs / pt_Xf_cgs
      
      
    enddo ! loop over momentum
  enddo ! loop over theta
  
return
end subroutine Xform_psd_corners


!****************************************************************************
!****************************************************************************
subroutine cosmo_calc(d_CM, z)

! Calculator to shift between redshift, lookback time, and comoving distance.
! Subroutine adapted from Hogg (1999). For more detail see
! [http://adsabs.harvard.edu/abs/1999astro.ph..5116H]
!
! Note one major difference: Equation (13) in Hogg uses Omega_r where this
! program uses Omega_k, and does not include a term for what this program
! calls Omega_r. For more information, see Wright (2006):
! [http://adsabs.harvard.edu/abs/2006PASP..118.1711W]
!
! The program calculates lookback time (in years) as well, but since this
!  result is not necessary for GRB redshift calculation it is not passed
!  back to the calling subroutine.
!
! Input values:
!  1) One of (redshift, comoving distance in Mpc)
! Output values:
!  1) Other of (redshift, comoving distance in Mpc)

implicit none

  ! Input/output arguments
real(kind=8), intent(inout) :: z, d_CM

  ! Parameters; cosmology values pulled from Planck 2013 results
                     ! Hubble constant, km/s Mpc^-1
real(kind=8), parameter :: H0 = 67.8d0
                     ! Fraction of density in neutrinos?
real(kind=8), parameter :: Omega_r = 0.4165d0 / H0**2
                     ! Fraction of density in dark energy and in matter
real(kind=8), parameter :: Omega_vac = 0.683d0 - 0.5d0*Omega_r
real(kind=8), parameter :: Omega_m   = 0.317d0 - 0.5d0*Omega_r
                     ! Assume flat Universe, i.e. Omega_k = 1 - sum(Omega_*)
real(kind=8), parameter :: Omega_k = 0.d0
                     ! Speed of light, km/s
real(kind=8), parameter :: c = 2.99792458d5
                     ! Hubble distance, Mpc
real(kind=8), parameter :: d_H = c / H0
                     ! Hubble time, years
real(kind=8), parameter :: t_H = 9.778d11 / H0
                     ! Number of steps to use in integration
integer, parameter ::  num_steps = 1000

  ! Local variables
integer :: code_stat, i
real(kind=8) :: t_look, z_save, d_save, d_old, d_new, z_old, E_old, z_new,&!&
     E_new, slope

1001  format(A, ES10.4E2, A)
2001  format(2(A, ES10.4E2))


!--! Set either redshift or comoving distance to a nonzero value. Other
   !   should be set to zero
  z_save = z
  d_save = d_CM
  code_stat = -1 ! Status of code
                 !  -1: initialization value
                 !   1: value of z for provided d_CM located. Proceed to
                 !       integration for lookback time
                 !   2: d_CM less than critical value. Skip stepping through
                 !       z integral and integration to find t_look
  
  
!--! Quick error check that at least one input is physically reasonable and
   !   that only one is positive
  if( ((z .le. 0.d0) .and. (d_CM .le. 0.d0)) .or.                         &!&
      ((z .gt. 0.d0) .and. (d_CM .gt. 0.d0)) ) then
    write(*,"(2A)") 'Invalid inputs given.  Exactly one of z and d_CM ',  &!&
                    'can/must be positive.'
    write(*,2001)   'Inputted values were z = ', z_save, ' and d_CM = ', d_CM
    write(*,"(A)") 'Stoping program now.'
    stop
  endif


!--! If d_CM is given, step through integral in z (Equation (14) in Hogg
   !   1999) until next step exceeds d_CM
   !-------------------------------------------------------------------------
  if( (z .eq. 0.d0) .and. (d_CM .gt. 0.d0)                                &!&
    ! First, check to make sure d_CM is large enough to get a reasonable
    !  redshift -- otherwise return zero
    .and. (d_CM .lt. 0.443d0) ) then
    
    z      = 0.d0
    t_look = d_CM * 3.2616d6 ! convert Mpc to years
    code_stat = 2 ! d_CM less than critical value. Skip stepping through
                    !  z integral and integration to find t_look
  
  
  ! d_CM was a reasonable value
  else if( (z .eq. 0.d0) .and. (d_CM .gt. 0.d0) ) then
  
    d_old = 0.d0
    d_new = 0.d0
    z_old = 0.d0
    E_old = 1.d0 ! E(z) = sqrt[Omega_r*(1+z)^4 + Omega_m*(1+z)^3 
                 !           + Omega_k*(1+z)^2 + Omega_vac]
                 !  Equation (13) in Hogg (1999)
  
    
    do i = 1, 100
      z_new = 1.d-4 * real(i)
      E_new = sqrt( Omega_r*(1+z_new)**4                  &!& Equation (13)
                   + Omega_m*(1+z_new)**3                 &!&  in Hogg (1999)
                   + Omega_vac)
      
      ! Use Equation (14) in Hogg (1999) to update d_new
      d_new  = d_new  +  d_H*1.d-4 * 0.5d0*(1.d0/E_old + 1.d0/E_new)
      
      ! If current step in z exceeded provided value of d_CM, perform linear
      !  interpolation between the d_old and d_new
      if( (d_old .lt. d_CM) .and. (d_new .ge. d_CM) ) then
        
        slope = (d_new - d_old) / (z_new - z_old)
        z = (d_CM - d_old) / slope  +  z_old
        code_stat = 1 ! value of z for provided d_CM located. Proceed to
                      !  integration for lookback time
        exit
        
      endif
      
      ! Set variables for next pass through loop
      d_old = d_new
      E_old = E_new
      z_old = z_new
    enddo
    
    
    ! If a value of z wasn't already found stepping from z = 0 to 0.01 with
    !   steps of size 1.d-4, start stepping up from z = 0.01 with steps of
    !   size 0.01
    i = 0
    do while (code_stat .lt. 1)
      i     = i + 1
      z_new = 1.d-2 * real(i)
      E_new = sqrt( Omega_r*(1+z_new)**4                  &!& Equation (13)
                   + Omega_m*(1+z_new)**3                 &!&  in Hogg (1999)
                   + Omega_vac)
      
      ! Use Equation (14) in Hogg (1999) to update d_new
      d_new  = d_new  +  d_H*1.d-2 * 0.5d0*(1.d0/E_old + 1.d0/E_new)
      
      ! If current step in z exceeded provided value of d_CM, perform linear
      !  interpolation between the d_old and d_new
      if( (d_old .lt. d_CM) .and. (d_new .ge. d_CM) ) then
        
        slope = (d_new - d_old) / (z_new - z_old)
        z = (d_CM - d_old) / slope  +  z_old
        code_stat = 1 ! value of z for provided d_CM located. Proceed to
                      !  integration for lookback time
        exit
        
      endif
      
      ! Set variables for next pass through loop
      d_old = d_new
      E_old = E_new
      z_old = z_new
      
      
      ! Quick and dirty error check
      if( i .gt. 1000) then
        write(*,"(A,F5.2)") "z_new exceeds reasonable bounds: ",z_new
        write(*,"(A)") "Stopping program"
        stop
      endif
      
    enddo
    
  endif
   !-------------------------------------------------------------------------
   ! End calculation of z for a provided d_CM


!--! If z is given, perform integral to find d_CM using the trapezoid rule
   !   and Equation (14) in Hogg (1999).
   ! Independently of the previous, calculate t_look using the trapezoid rule
   !   and Equations (28) in Hogg (1999).
   !-------------------------------------------------------------------------
  if( ((z .gt. 0.d0).and.(d_CM .eq. 0.d0)) .or. (code_stat .eq. 1) ) then
    
    t_look = 0.d0
    z_old = 0.d0
    E_old = 1.d0 ! E(z) = sqrt[Omega_m*(1+z)^3 + Omega_r*(1+z)^2 + Omega_vac]
                 !  Equation (13) in Hogg (1999)
  
    
    do i = 1, num_steps
      z_new = z*real(i)/num_steps
      E_new = sqrt( Omega_r*(1+z_new)**4                  &!& Equation (13)
                   + Omega_m*(1+z_new)**3                 &!&  in Hogg (1999)
                   + Omega_vac)
  
      
      ! If d_CM wasn't provided, use Equation (14) in Hogg (1999) to update
      !  running total
      if( code_stat .lt. 1 )  d_CM = d_CM  +  0.5d0*(1.d0/E_old + 1.d0/E_new)
      
      ! Independently of previous calculation, use Equation (28) in Hogg
      !   (1999) to update t_look
      t_look = t_look  +  0.5d0 * (  1.d0/((1+z_old)*E_old)               &!&
                                   + 1.d0/((1+z_new)*E_new) )
      
      ! Set variables for next pass through loop
      z_old = z_new
      E_old = E_new
    enddo
    
    ! Finally, multiply by appropriate scale factors out front and width of
    !  the trapezoids used for integration
    if( code_stat .lt. 1 )  d_CM = d_H * z/real(num_steps) * d_CM
    t_look = t_H * z/real(num_steps) * t_look
  
  endif
   !-------------------------------------------------------------------------
   ! End calculation of t_look, and possibly also d_CM


!--! Debugging output lines
!comm   write(*,"(A)") "Program finished."
!comm   write(*,2001) "Inputted values were z = ", z_save,                &!&
!comm                 " and d_CM = ", d_save
!comm   write(*,*)
!comm   write(*,1001) "Calculated values: z = ", z
!comm   write(*,1001) "                d_CM = ", d_CM, " Mpc"
!comm   write(*,1001) "              t_look = ", t_look, " years"

return
end subroutine cosmo_calc
