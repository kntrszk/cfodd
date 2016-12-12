program analysis_cfodd_gfdl

  implicit none
  include '/usr/local/netcdf-3.6.3/include/netcdf.inc'

  integer, parameter :: imax = 48, jmax = 48, kmax = 48, nmon = 12, ncolumn = 25, nreg = 6
  integer, parameter :: nmax_dbz = 25, nmax_tau = 15, nmax_cls = 4
  real(4), parameter :: rair = 287.04, grav = 9.801, undef = -999.0
  integer :: ncid, rhid, ntmax( nmon )
  integer :: status
  integer :: kctop( imax,jmax ), kcbtm( imax,jmax ), kcctl( imax,jmax )
  logical :: octop( imax,jmax ), ocbtm( imax,jmax )
  real(4) :: tctop( imax,jmax ), delz
  real(4), allocatable :: gridxt( :,: ), gridyt( :,: )
  real(4), allocatable :: pres( : ), phalf( : ), temp( :,:,:,: ), qliq( :,:,:,: ), qice( :,:,:,: ), deff( :,:,:,: )
  real(4), allocatable :: retop( :,: ), reft( :,:,: ), cnum( :,:,: ), tauc( :,: ), cldf( :,:,: ), dtaus( :,:,:,: ), dtauc( :,:,:,: )
  real(4), allocatable :: dbze( :,:,:,:,: ), ctyp( :,:,:,:,: ), cldfrc( :,:,:,: ), zsfc( :,: )
  real(4), allocatable :: flagls( :,: ), precip( :,:,: ), prec_ls( :,:,: ), prec_conv( :,:,: ), clwp( :,:,: )
  real(4) :: cnt( nmax_dbz,nmax_tau,nmax_cls ), taud, height( imax,jmax,kmax+1 )
  real(4) :: tau_min, tau_max, tau_del, tau_bin( nmax_tau )
  real(4) :: dbz_min, dbz_max, dbz_del, dbz_bin( nmax_dbz )
  real(4) :: sum_dbz, ave_dbze, dens, qcld, dbze_max( ncolumn )
  real(4) :: ratio, temp_half
  character(150) :: path, period( nmon ), filename( nreg ), filegrid( nreg ), var_dbze( ncolumn ), var_ctyp( ncolumn )
  integer :: i, j, k, nt, n, nc, nr, nm, n_dbz, n_tau, n_ref, n_cls

  data period / '20080101', '20080201', '20080301', '20080401', '20080501', '20080601', &
                '20080701', '20080801', '20080901', '20081001', '20081101', '20081201' /

  data var_dbze / 'dbze_1' , 'dbze_2' , 'dbze_3' , 'dbze_4' , 'dbze_5' , &
                  'dbze_6' , 'dbze_7' , 'dbze_8' , 'dbze_9' , 'dbze_10', &
                  'dbze_11', 'dbze_12', 'dbze_13', 'dbze_14', 'dbze_15', &
                  'dbze_16', 'dbze_17', 'dbze_18', 'dbze_19', 'dbze_20', &
                  'dbze_21', 'dbze_22', 'dbze_23', 'dbze_24', 'dbze_25' /

  data var_ctyp / 'cloud_type_1', 'cloud_type_2', 'cloud_type_3', 'cloud_type_4', 'cloud_type_5', &
                  'cloud_type_6', 'cloud_type_7', 'cloud_type_8', 'cloud_type_9', 'cloud_type_10', &
                  'cloud_type_11', 'cloud_type_12', 'cloud_type_13', 'cloud_type_14', 'cloud_type_15', &
                  'cloud_type_16', 'cloud_type_17', 'cloud_type_18', 'cloud_type_19', 'cloud_type_20', &
                  'cloud_type_21', 'cloud_type_22', 'cloud_type_23', 'cloud_type_24', 'cloud_type_25' /

  data ntmax / 124, 116, 124, 120, 124, 120, 124, 124, 120, 124, 120, 124 /

  !! Standard AM3
  path = '../../data/GFDL_CM3/'  !! Rcrit = 8.2micron (Default)
  ! path = '../../COSP/GFDL/c48L48_am3p10_14_cospR/'  !! Rcrit = 10.6micron
  ! path = '../../COSP/GFDL/c48L48_am3p10_11_cospR/'  !! Rcrit = 6.0micron

  open( 21, file = 'cfodd_gfdl_sg_std_jan_4class.txt', form = 'formatted' )

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

  cnt( :,:,: ) = 0.0

  do nm = 1, 1

  write( *,* ) ' nm = ', nm

  filename( 1 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.atmos_4xdaily.tile1.nc'
  filename( 2 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.atmos_4xdaily.tile2.nc'
  filename( 3 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.atmos_4xdaily.tile3.nc'
  filename( 4 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.atmos_4xdaily.tile4.nc'
  filename( 5 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.atmos_4xdaily.tile5.nc'
  filename( 6 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.atmos_4xdaily.tile6.nc'

  filegrid( 1 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.grid_spec.tile1.nc'
  filegrid( 2 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.grid_spec.tile2.nc'
  filegrid( 3 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.grid_spec.tile3.nc'
  filegrid( 4 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.grid_spec.tile4.nc'
  filegrid( 5 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.grid_spec.tile5.nc'
  filegrid( 6 ) = trim(path)//trim(period( nm ))//'.nc/'//trim(period( nm ))//'.grid_spec.tile6.nc'

  allocate( gridxt( imax,jmax ) )                         ! longitude
  allocate( gridyt( imax,jmax ) )                         ! latitude
  allocate( flagls( imax,jmax ) )                         ! land-ocean flag (fractional coverage of land)
  allocate( zsfc( imax,jmax ) )                           ! surface height [m]
  allocate( pres( kmax ),phalf( kmax+1 ) )                ! pressure    [mb]
  allocate( temp( imax,jmax,kmax,ntmax( nm ) ) )          ! temperature [K]
  allocate( qliq( imax,jmax,kmax,ntmax( nm ) ) )          ! liquid water content [kg/kg]
  allocate( qice( imax,jmax,kmax,ntmax( nm ) ) )          ! ice water content [kg/kg]
  allocate( deff( imax,jmax,kmax,ntmax( nm ) ) )          ! effective diameter for stratiform liquid clouds [micron]
  allocate( cldfrc( imax,jmax,kmax,ntmax( nm ) ) )        ! cloud fraction
  allocate( retop( imax,jmax ) )                          ! effective radius @ cloud top
  allocate( reft( imax,jmax,ntmax( nm ) ) )               ! liquid droplet effective radius @ cloud top
  allocate( cnum( imax,jmax,ntmax( nm ) ) )               ! liquid droplet number concentration @ cloud top
  allocate( tauc( imax,jmax ) )                           ! total liquid cloud optical depth
  allocate( cldf( imax,jmax,ntmax( nm ) ) )               ! cloud fraction @ cloud top
  allocate( dtaus( imax,jmax,kmax,ntmax( nm ) ) )         ! in-cloud optical depth (stratiform)
  allocate( dtauc( imax,jmax,kmax,ntmax( nm ) ) )         ! in-cloud optical depth (convective)
  allocate( dbze( imax,jmax,kmax,ntmax( nm ),ncolumn ) )  ! radar reflectivity profile
  allocate( ctyp( imax,jmax,kmax,ntmax( nm ),ncolumn ) )  ! cloud type: clear (=0), stratiform (=1), convective (=2)
  allocate( precip( imax,jmax,ntmax( nm ) ) )             ! precipitation rate [kg/m2/sec]
  allocate( prec_ls( imax,jmax,ntmax( nm ) ) )            ! precipitation rate (large-scale cond.) [kg/m2/sec]
  allocate( prec_conv( imax,jmax,ntmax( nm ) ) )          ! precipitation rate (convective)        [kg/m2/sec]
  allocate( clwp( imax,jmax,ntmax( nm ) ) )               ! liquid water path (large-scale cond.) [kg/m2]

  do nr = 1, nreg

    write( *,* ) '   nr = ', nr

    !! File opening for grid configuration
    status = nf_open(trim(filegrid( nr )), nf_nowrite, ncid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Longitude
    status = nf_inq_varid(ncid, 'grid_lont', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, gridxt)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Latitude
    status = nf_inq_varid(ncid, 'grid_latt', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, gridyt)
    if ( status .ne. nf_noerr ) call handle_err(status)

    !! File closing for grid configuration
    status = nf_close(ncid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    !! File opening for atmospheric variables
    status = nf_open(trim(filename( nr )), nf_nowrite, ncid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Land-Ocean Mask
    status = nf_inq_varid(ncid, 'land_mask', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, flagls)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Pressure (full level)
    status = nf_inq_varid(ncid, 'pfull', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, pres)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Pressure (half level)
    status = nf_inq_varid(ncid, 'phalf', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, phalf)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Surface Height
    status = nf_inq_varid(ncid, 'zsurf', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, zsfc)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Temperature
    status = nf_inq_varid(ncid, 'temp', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, temp)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Liquid Water Content
    status = nf_inq_varid(ncid, 'liq_wat', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, qliq)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Ice Water Content
    status = nf_inq_varid(ncid, 'ice_wat', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, qice)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Cloud Fraction
    status = nf_inq_varid(ncid, 'cld_amt', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, cldfrc)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Cloud Top Effective Radius
    status = nf_inq_varid(ncid, 'cldtop_reff', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, reft)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Cloud Top Droplet Number Concentration
    status = nf_inq_varid(ncid, 'cldtop_dropnum', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, cnum)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Cloud Top Area
    status = nf_inq_varid(ncid, 'cldtop_area', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, cldf)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Cloud Optical Depth (Stratiform)
    status = nf_inq_varid(ncid, 'dtau_s', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, dtaus)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Cloud Optical Depth (Convective)
    status = nf_inq_varid(ncid, 'dtau_c', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, dtauc)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Precipitation (total)
    status = nf_inq_varid(ncid, 'precip', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, precip)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Precipitation (large-scale)
    status = nf_inq_varid(ncid, 'prec_ls', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, prec_ls)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Precipitation (convective)
    status = nf_inq_varid(ncid, 'prec_conv', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, prec_conv)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Particle size (stratiform)
    status = nf_inq_varid(ncid, 'strat_size_drop', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, deff)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Liquid Water Path (large-scale)
    status = nf_inq_varid(ncid, 'LWP', rhid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    status = nf_get_var_real(ncid, rhid, clwp)
    if ( status .ne. nf_noerr ) call handle_err(status)

    ! Radar Reflectivity
    do nc = 1, ncolumn
      status = nf_inq_varid(ncid, trim(var_dbze( nc )), rhid)
      if ( status .ne. nf_noerr ) call handle_err(status)

      status = nf_get_var_real(ncid, rhid, dbze( :,:,:,:,nc ))
      if ( status .ne. nf_noerr ) call handle_err(status) 
    end do

    ! Cloud Type
    do nc = 1, ncolumn
      status = nf_inq_varid(ncid, trim(var_ctyp( nc )), rhid)
      if ( status .ne. nf_noerr ) call handle_err(status)

      status = nf_get_var_real(ncid, rhid, ctyp( :,:,:,:,nc ))
      if ( status .ne. nf_noerr ) call handle_err(status)
    end do

    !! File closing for atmospheric variables
    status = nf_close(ncid)
    if ( status .ne. nf_noerr ) call handle_err(status)

    do nt = 1, ntmax( nm )

      octop( :,: ) = .false.
      ocbtm( :,: ) = .false.
      kctop( :,: ) = 0
      kcbtm( :,: ) = 0
      kcctl( :,: ) = 0
      tctop( :,: ) = undef
      retop( :,: ) = undef
      tauc ( :,: ) = 0.0

      do j = 1, jmax
      do i = 1, imax

        ! Geopotential Height
        height( i,j,kmax+1 ) = 0.0

        delz = rair*temp( i,j,kmax,nt )/grav*log( phalf( kmax+1 )/pres( kmax ) )
        height( i,j,kmax ) = height( i,j,kmax+1 ) + delz

        do k = kmax-1, 1, -1
          ratio = log( pres( k+1 )/phalf( k+1 ) )/log( pres( k+1 )/pres( k ) )
          temp_half = temp( i,j,k,nt )*ratio + temp( i,j,k+1,nt )*( 1.0-ratio )
          delz = rair*temp_half/grav*log( pres( k+1 )/pres( k ) )
          height( i,j,k ) = height( i,j,k+1 ) + delz
        end do

        ! Determining Cloud Top Temperature
        do k = 1, kmax
          if ( temp( i,j,k,nt ) /= -1.e+10 .and. qliq( i,j,k,nt ) /= -1.e+10 .and. &
               qice( i,j,k,nt ) /= -1.e+10 ) then
            dens = pres( k )*100.0/( rair*temp( i,j,k,nt ) )   !  [kg/m3]
            qcld = dens*( qliq( i,j,k,nt )+qice( i,j,k,nt ) )
            if ( ( .not. octop( i,j ) ) .and. qcld >= 1.e-04 ) then
              octop( i,j ) = .true.
              kctop( i,j ) = k
              tctop( i,j ) = temp( i,j,k,nt )
            end if
            if ( octop( i,j ) .and. qcld >= 1.e-04 ) then
              kcctl( i,j ) = k
            end if
          end if
        end do

        ! Determining Cloud Base
        do k = kmax, 1, -1
          if ( temp( i,j,k,nt ) /= -1.e+10 .and. qliq( i,j,k,nt ) /= -1.e+10 .and. &
               qice( i,j,k,nt ) /= -1.e+10 ) then
            dens = pres( k )*100.0/( rair*temp( i,j,k,nt ) )   !  [kg/m3]
            qcld = dens*( qliq( i,j,k,nt )+qice( i,j,k,nt ) )
            if ( ( .not. ocbtm( i,j ) ) .and. qcld >= 1.e-04 ) then
              ocbtm( i,j ) = .true.
              kcbtm( i,j ) = k
            end if
          end if
        end do

        if ( flagls( i,j ) == 0.0 .and. octop( i,j ) .and. ocbtm( i,j ) .and. &
             kcbtm( i,j ) == kcctl( i,j ) .and. tctop( i,j ) >= 273.15 ) then

          if ( cldf( i,j,nt ) > 0.0 .and. cldf( i,j,nt ) <= 1.0 ) then
            reft( i,j,nt ) = reft( i,j,nt )/cldf( i,j,nt ) * 1.e+06
          else
            reft( i,j,nt ) = -999.0
          end if

          if ( reft( i,j,nt ) >= 5.0 .and. reft( i,j,nt ) < 10.0 ) then
            n_cls = 1
          else if ( reft( i,j,nt ) >= 10.0 .and. reft( i,j,nt ) < 15.0 ) then
            n_cls = 2
          else if ( reft( i,j,nt ) >= 15.0 .and. reft( i,j,nt ) < 20.0 ) then
            n_cls = 3
          else if ( reft( i,j,nt ) >= 20.0 .and. reft( i,j,nt ) < 25.0 ) then
            n_cls = 4
          else
            n_cls = 0
          end if

          taud = 0.0
          dbze_max( 1:ncolumn ) = -999.0
          do k = kctop( i,j ), kmax
           if ( height( i,j,k ) >= 1000.0 ) then
            taud = taud + dtaus( i,j,k,nt )
            ! taud = tau_sum*( 1.0 - (height( i,j,k )/delz_sum)**(5.0/3.0) )
            n_tau = int( ( taud-tau_min )/tau_del ) + 1
            do nc = 1, ncolumn
              if ( dbze( i,j,k,nt,nc ) /= 1.e+20 .and. ctyp( i,j,k,nt,nc ) == 1.0 ) then ! stratiform only
                dbze_max( nc ) = max( dbze_max( nc ), dbze( i,j,k,nt,nc ) )
                n_dbz = int( ( dbze( i,j,k,nt,nc )-dbz_min )/dbz_del ) + 1
                if ( n_cls >= 1 .and. n_cls <= nmax_cls .and. &
                     n_tau >= 1 .and. n_tau <= nmax_tau .and. &
                     n_dbz >= 1 .and. n_dbz <= nmax_dbz ) then
                  cnt( n_dbz,n_tau,n_cls ) = cnt( n_dbz,n_tau,n_cls ) + 1.0
                end if
              end if
            end do
           end if
          end do 

        end if

      end do  ! for i = 1, imax
      end do  ! for j = 1, jmax

    end do  !  for nt = 1, ntmax

  end do ! for nr = 1, nreg

  deallocate( flagls,pres,phalf,temp,qliq,qice,cldfrc,deff,retop,reft,tauc,cnum,cldf,dtaus,dtauc,dbze,ctyp )
  deallocate( precip,prec_ls,prec_conv,clwp,zsfc,gridxt,gridyt )

  end do ! for nm = 1, nmon

  do n_dbz = 1, nmax_dbz
  do n_tau = 1, nmax_tau
    write( 21,* ) dbz_bin( n_dbz ), tau_bin( n_tau ), &
                  ( cnt( n_dbz,n_tau,n_cls )/sum( cnt( :,n_tau,n_cls ) )/dbz_del, n_cls = 1, nmax_cls )
  end do
  end do

  stop

contains

  subroutine handle_err(status)

  integer :: status

  write(*,*) nf_strerror(status)
  stop

  end subroutine handle_err

end program analysis_cfodd_gfdl
