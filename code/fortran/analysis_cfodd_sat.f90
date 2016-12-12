program analysis_cfodd_sat

  implicit none

  real(4), parameter :: undef = -999.0, delz = 240.0, rhow = 1000.0, pi = 3.141592

  integer, parameter :: nfstr = 1, nfend = 1, nskip = 1
  integer, parameter :: nmax_dbz = 25, nmax_tau = 15, nmax_cls = 5
  integer :: n, n_temp, n_tau, n_tau_min, n_tau_max, n_dbz, n_ref, n_ret, n_rwc, n_rwc_s, n_cls, n_cls_min, n_cls_max, n_aero, n_ltss
  real(4) :: taud, taud_min, taud_max
  character(6) :: character

  real(4) :: dbz_min, dbz_max, dbz_del, dbz_bin( nmax_dbz )
  real(4) :: tau_min, tau_max, tau_del, tau_bin( nmax_tau )
  real(4) :: count_cfodd( nmax_dbz,nmax_tau,nmax_cls )

  integer :: i, j, nf

  real(4) :: time, lon, lat, htop, hbtm, halt, dbz, tau, ai
  real(4) :: elev, temp, pres, ltss, pia, precip, clwpmw, aot, alfa_ocean, alfa_land
  integer :: lsflag, ipflag, ipstat
  real(4) :: pia_itpl, pia_itpl_sigma, cldf_pia, dst_pia
  real(4) :: reff_gsfc, reff_gsfc2, reff_sgm, tauc_gsfc, tauc_sgm, diff_reff_gsfc, tc_gsfc, tctop
  integer :: icphase_gsfc, icphase_ret

  logical :: ofqrt, omedn, otqrt
  real(4) :: fqrt_dbze, medn_dbze, tqrt_dbze
  real(4) :: cum_count_cfodd

  dbz_min = -30.0
  dbz_max =  20.0
  dbz_del = ( dbz_max-dbz_min )/nmax_dbz
  do n = 1, nmax_dbz
    dbz_bin( n ) = dbz_min + dbz_del*0.5 + dbz_del*( n-1 )
  end do

  tau_min = 0.0
  tau_max = 60.0
  tau_del = ( tau_max-tau_min )/nmax_tau
  do n = 1, nmax_tau
    tau_bin( n ) = tau_min + tau_del*0.5 + tau_del*( n-1 )
  end do

  open( 11, file = '/Volumes/NO NAME/A-Train/texts/2008djf_A-Train.txt', form = 'formatted' )

  open( 21, file = 'cfodd_r21_5class_djf.txt', form = 'formatted' )

  count_cfodd( :,:,: ) = 0.0

  do nf = nfstr, nfend, nskip
  write( *,* ) ' nf = ', nf
  do
    read( 10+nf,*,end=100 ) time, lon, lat, htop, hbtm, halt, dbz, &
                            reff_gsfc, tauc_gsfc, diff_reff_gsfc, tc_gsfc, icphase_gsfc, icphase_ret, &
                            reff_sgm, tauc_sgm, &
                            lsflag, elev, temp, pres, ltss, &
                            pia, precip, ipflag, ipstat, &
                            pia_itpl, pia_itpl_sigma, cldf_pia, dst_pia, &
                            clwpmw, aot, alfa_ocean, alfa_land

    reff_gsfc2 = reff_gsfc + diff_reff_gsfc

    if ( halt == htop ) tctop = temp

    if ( lsflag == 2 .and. icphase_gsfc == 1 .and. icphase_ret == 2 .and. &
         reff_gsfc > 0.0 .and. reff_gsfc <= 90.0 .and. &
         tc_gsfc > 273.15 .and. tctop > 273.15 .and. &
         tauc_gsfc > 0.0 .and. tauc_gsfc <= 100.0 ) then

      n_cls = 0

      ! effective radius (2.1um)
      if ( reff_gsfc >= 5.0 .and. reff_gsfc < 10.0 ) then
        n_cls = 1
      else if ( reff_gsfc >= 10.0 .and. reff_gsfc < 15.0 ) then
        n_cls = 2
      else if ( reff_gsfc >= 15.0 .and. reff_gsfc < 20.0 ) then
        n_cls = 3
      else if ( reff_gsfc >= 20.0 .and. reff_gsfc < 25.0 ) then
        n_cls = 4
      else if ( reff_gsfc >= 25.0 .and. reff_gsfc < 30.0 ) then
        n_cls = 5
      end if

      if ( halt == hbtm ) then
        taud = tauc_gsfc*( 1.0-( (halt-hbtm)/(htop-hbtm) )**(5.0/3.0) )
        n_tau = int( ( taud-tau_min )/tau_del ) + 1
        n_dbz = int( ( dbz-dbz_min )/dbz_del ) + 1
        if ( n_dbz >= 1 .and. n_dbz <= nmax_dbz .and. &
             n_tau >= 1 .and. n_tau <= nmax_tau .and. &
             n_cls >= 1 .and. n_cls <= nmax_cls ) then
          count_cfodd( n_dbz,n_tau,n_cls ) = count_cfodd( n_dbz,n_tau,n_cls ) + 1.0
        end if
      else
        taud = tauc_gsfc*( 1.0-( (halt-hbtm)/(htop-hbtm) )**(5.0/3.0) )
        n_tau = int( ( taud-tau_min )/tau_del ) + 1
        n_dbz = int( ( dbz-dbz_min )/dbz_del ) + 1
        if ( n_dbz >= 1 .and. n_dbz <= nmax_dbz .and. &
             n_tau >= 1 .and. n_tau <= nmax_tau .and. &
             n_cls >= 1 .and. n_cls <= nmax_cls ) then
          count_cfodd( n_dbz,n_tau,n_cls ) = count_cfodd( n_dbz,n_tau,n_cls ) + 1.0
        end if
      end if

    end if   

  end do  

100 continue

  end do

  do n_dbz = 1, nmax_dbz
  do n_tau = 1, nmax_tau
    write( 21,* ) dbz_bin( n_dbz ), tau_bin( n_tau ), &
                  ( count_cfodd( n_dbz,n_tau,n_cls )/sum( count_cfodd( :,n_tau,n_cls ) )/dbz_del, n_cls = 1, nmax_cls )
  end do
  end do

  stop

end program analysis_cfodd_sat
