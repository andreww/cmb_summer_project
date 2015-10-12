module tomo_predict

   implicit none

   ! f2py doesn't work with this at the moment, so everything is public
!   private

   ! Global noisyness setting
   logical :: debug = .false.

   ! 1D and 3D models
   ! TODO: Make model more flexible
   integer, parameter :: nlat = 89, nlon = 180, nz = 20
   integer, parameter :: nr1dmax = 10000
   real, parameter :: rmax = 6370.
   ! 3D model layout
   ! top(20): Array containing the depth to the top of each layer in the 3D
   !          model (km)
   ! bot(20): Array containing the depth to the bottom of each layer in the 3D
   !          model (km)
   ! tomo(89,180,20): Array containing velocity deviations from the 1D model (%)
   !          on a grid evaluated every 2 degrees in lat and lon between -88 to 88
   !          latitude, -178 to 180 longitude, with 20 radial layers which have
   !          their bottom depth and top depths respectively in arrays bot and top.
   !          (NB: I'm not entirely sure this is correct!)
   real, save :: tomo(nlat,nlon,nz), top_layer(nz), bot_layer(nz), v1d(nr1dmax), &
                 z1d(nr1dmax)
   integer, save :: nr1d
   ! Once this is set to .true., the model is fixed and cannot be replaced.
   ! TODO: Add ability to change model
   logical, save :: setup_done = .false.

   public :: predict, read_taup_time_file, setup, tomo_delay

contains


subroutine predict(model1d_file, model3d_file, lat, lon, dep, dt)
! predict is the main routine called by the Python module.
   character(len=*), intent(in) :: model1d_file, model3d_file
   real, intent(in) :: lat(:), lon(:), dep(:)
   real, intent(out) :: dt

   call setup(model1d_file, model3d_file)
   call tomo_delay(lat, lon, dep, dt)
end subroutine predict


subroutine setup(model1d_file, model3d_file)
! Read in 3D and 1D model files once and for all.  This routine must be called before
! calling tomo_delay
! INPUT:
!     model1d_file: Path (absolute or relative) to file containing 1D model
!     model3d_file: Path to file containing 3D model (see format below)
!
! model1d_file format:
!     This contains lines with the integer depth, followed by the wave velocity:
!     1 5.8
!     2 5.8
!
! model3d_file format:
!     This contains, one number per line, the deviation in % of the velocity at
!     each node compared to a 1D model.  The longitude loops quickest, then then
!     latitude, then the depth.  It must have 20 layers.
!     [1] layer1_top layer1_bottom (km)
!     [2] dV at lat=-88, lon=-178
!     [3] dV at lat=-88, lon=-176
!     [4+] etc.
!     [16,021] layer2_top layer2_bottom
!     [16,022+] (etc.)
   character(len=*), intent(in) :: model1d_file, model3d_file

   if (setup_done) return
   call read_1d_model(model1d_file, nr1dmax, z1d, v1d, nr1d)
   call read_3d_model(model3d_file, nlat, nlon, nz, top_layer, bot_layer, tomo)
   setup_done = .true.
end subroutine setup


subroutine tomo_delay(lat, lon, dep, dt)
! subroutine tomo_delay gives the travel time perturbation for a predefined
! ray path within a 3D model of velocity perturbations with respect to a 1D model.
! Positive times mean the phase arrives later in time than expected in the 1D model.
! INPUT:
!     lat         : Array of latitudes (degrees)
!     lon         : Array of longitudes (degrees)
!     dep         : Array of depths (km)
! OUTPUT:
!     dt          : Traveltime perturbation (s)
   ! IO
   real, intent(in) :: lat(:), lon(:), dep(:)
   real, intent(out) :: dt
   real :: path_len, dt_segment, resid
   integer :: i, iz

   if (.not.setup_done) error stop 'tomo_delay: Must call setup() first'

   if (size(lat) /= size(lon) .or. size(lon) /= size(dep)) &
      error stop 'tomo_delay: lat, lon and dep arrays must be the same length'
   if (any(dep > rmax)) error stop 'tomo_delay: maximum depth is Earth radius (6370 km)'

   dt = 0.

   do i = 2, size(lat)
      path_len = distance(lon(i), lat(i), rmax - dep(i), &
                          lon(i-1), lat(i-1), rmax - dep(i-1))
      ! Where we assume that the 1D model has z(i) == i
      iz = min(max(nint(dep(i)), 1), nr1d)
      call get_dt(lat(i), lon(i), dep(i), v1d(iz), path_len, dt_segment, resid)
      dt = dt + dt_segment
      if (debug) write(0,'(a,6(1x,f7.2),2(1x,f9.4))') &
         'tomo_delay:', dep(i), path_len, v1d(iz), resid,  lat(i), lon(i), dt_segment, dt
   enddo
end subroutine tomo_delay


function distance(lon1, lat1, r1, lon2, lat2, r2) result(d)
   ! Calculate the 3D distance between two geographic points in a sphere
   ! INPUT:
   !     lon1, lat1, lon2, lat2 : Longitudes and latitudes of points (degrees)
   !     r1, r2                 : Radii of points (km)
   ! RETURNS:
   !     distance               : Distance between points (km)
   real, intent(in) :: lon1, lat1, r1, lon2, lat2, r2
   real :: d, x1(3), x2(3)

   if (any(abs([lat1, lat2]) > 90.)) error stop 'distance: Latitude must be in rage -90 to 90 degrees'
   x1 = lon_lat_r2vec(lon1, lat1, r1)
   x2 = lon_lat_r2vec(lon2, lat2, r2)
   d = sqrt(sum((x2 - x1)**2))
end function distance


function lon_lat_r2vec(lon, lat, r) result(x)
   real, intent(in) :: lon, lat, r
   real :: x(3)
   real, parameter :: torad = atan2(1., 1.)/45.
   real :: rlon, rlat

   rlon = lon*torad
   rlat = lat*torad
   x = r*[cos(rlon)*cos(rlat), sin(rlon)*cos(rlat), sin(rlat)]
end function lon_lat_r2vec


subroutine read_1d_model(file, nmax, z, v, n)
   ! subroutine read_1d_model returns the depths and velocities of a 1D model
   ! File must contain:
   !  [1..n] depth, velocity (km, km/s)
   ! Maximum depth is determined by the file length
   character(len=*), intent(in) :: file
   integer, intent(in) :: nmax
   real, intent(out) :: z(nmax), v(nmax)
   integer, intent(out) :: n
   integer :: ier

   open(10, file=file, iostat=ier)
   if (ier /= 0) then
      write(0,'(a)') 'read_1d_model: Cannot open file "' // trim(file) // '"'
      error stop
   endif
   n = 1
   do
      read(10, *, iostat=ier) z(n), v(n)
      if (ier > 0) error stop 'read_1d_model: Error reading 1D velocity file'
      if (ier < 0) then
         n = n - 1
         exit
      endif
      if (n > 1) then
         if (z(n) <= z(n-1)) error stop 'read_1d_model: Depth must increase in file'
      endif
      n = n + 1
   enddo
   close(10)
   if (n > nmax) error stop 'read_1d_model: Model supplied is too long for arrays'
   if (debug) write(0,'(a,i0.1,a)') 'read_1d_model: Read ', n, ' points from file "' &
      // trim(file) // '"'
end subroutine read_1d_model


subroutine read_3d_model(file, nlat, nlon, nlayers, ztop, zbottom, dv)
   ! subroutine read_3d_model reads a 3D velocity perturbation model file
   character(len=*), intent(in) :: file
   integer, intent(in) :: nlat, nlon, nlayers
   real, intent(out) :: ztop(nlayers), zbottom(nlayers), dv(nlat, nlon, nlayers)
   integer :: ier, iz, ilat, ilon

   open(10, file=file, iostat=ier)
   if (ier /= 0) then
      write(0,'(a)') 'read_3d_model: Cannot open file "' // trim(file) // '"'
      error stop
   endif
   do iz = 1, nlayers
      read(10, *, iostat=ier) ztop(iz), zbottom(iz)
      if (ier /= 0) error stop 'read_3d_model: Error reading ztop, zbottom'
      do ilat = nlat, 1, -1
         do ilon = 1, nlon
            read(10, *, iostat=ier) dv(ilat,ilon,iz)
            if (ier /= 0) error stop 'read_3d_model: Error reading dv'
         enddo
      enddo
   enddo
   close(10)
   if (debug) write(0,'(a)') 'read_3d_model: Read model file "' // trim(file) // '"'
end subroutine read_3d_model


subroutine get_dt(lat1, lon1, dep, VPREM, sddp, dt, resid)
! subroutine get_dt computes the travel time perturbation between a 1D velocity
! model and a 3D one, given the location and path length of a point within the models.
! NB: The formulation assumes:
!     - A smoothly varying model evaluated on a lat-lon-depth regular grid
!       which can be linearly interpolated between points
!     - The point supplied experiences a uniform velocity along its path length
!       (despite the first assumption!)
!     The implementation assumes:
!     - The 3D model has 180 longitude points, 89 latitude points, and 20 layers,
!       where point tomo(1,2,3) is at (lat,lon) = (2,4), and so on.
! Written by Edward J. Garnero
! Modified by Andy Nowacki, University of Leeds (a.nowacki@leeds.ac.uk)
! to avoid implicit variables and declare intents, and reformat.
! INPUT:
!     lat1   : Latitude (degrees)
!     lon1   : Longitude (degrees)
!     dep    : Depth (km)
!     VPREM  : velocity of PREM at this depth
!     sddp   : Path length along which to accrue delay (km)
! OUTPUT:
!     dt     : Difference in travel time for this point between the 1D and 3D models (s)
!
   real, intent(in) :: lat1, lon1, dep, VPREM, sddp
   real, intent(out) :: dt, resid
   integer :: izi, layer, lathi, latlo, lonhi, lonlo, ilatlo, ilathi, ilonlo, ilonhi
   real :: fact, vlatlo, vlathi, Tprem, Tanom

   ! get VDH model depth index
   do izi = 1, 20
      if (dep >= top_layer(izi) .and. dep <= bot_layer(izi)) then
         layer = izi
         exit
      endif
      if (izi == 20) layer = izi
   enddo

   ! get VDH grid pt coords surrounding lat,lon
   lathi = int(lat1/2.)*2
   if (lat1 > 0.) then
       latlo = lathi
       lathi = lathi + 2
   else
       latlo = lathi - 2
   endif
   if (latlo < -88) latlo = -88
   if (lathi > 88) lathi = 88
   lonhi = int(lon1/2.)*2
   if (lon1.gt.0 ) then
       lonlo = lonhi
       lonhi = lonhi + 2
   else
       lonlo = lonhi - 2
   endif
   lonlo = lonhi - 2
   ilatlo = 45 + latlo/2
   ilathi = ilatlo + 1
   ilonlo = 90 + lonlo/2
   ilonhi = ilonlo + 1

   ! get resids @ corners, then @ lat1,lon1 (all in %)
   fact = (lat1 - real(latlo)) / 2.
   vlatlo = tomo(ilatlo,ilonlo,layer)*(1-fact) + fact*tomo(ilathi,ilonlo,layer)
   vlathi = tomo(ilatlo,ilonhi,layer)*(1-fact) + fact*tomo(ilathi,ilonhi,layer)
   fact = (lon1 - real(lonlo)) / 2.
   resid = vlatlo + fact*(vlathi - vlatlo)
   ! get time delay assoc. with resid
   Tprem = sddp/VPREM
   Tanom = sddp/(VPREM*(1 + resid/100.))
   dt = Tanom - Tprem

   if (debug) write(0,'(a,f0.4,a)') 'get_dt: ', dt, ' s'
end subroutine get_dt


subroutine read_taup_time_file(file, lat, lon, r, n)
   ! subroutine read_taup_time_file reads the output from the Java TauP package
   ! program `taup_time'.
   ! This is in the format:
   !     > Phase_name at ... [header line]
   !     distance/deg radius/km lat/deg lon/deg [n times]
   character(len=*), intent(in) :: file
   real, intent(out), dimension(:) :: lat, lon, r
   real :: d
   integer, intent(out) :: n
   integer :: ier

   if (size(lat) /= size(lon) .or. size(lon) /= size(r)) &
      error stop 'read_taup_time_file: All output arrays must be same length'
   open(10, file=file, iostat=ier)
   if (ier /= 0) error stop 'read_taup_time_file: Cannot open input file'
   read(10,*) ! Skip header
   n = 1
   do
      read(10, *, iostat=ier) d, r(n), lat(n), lon(n)
      if (ier < 0) then
         n = n - 1
         exit
      endif
      if (ier > 0) error stop 'read_taup_time_file: Problem reading ' // &
         'distance, radius, lat, lon from file'
      n = n + 1
      if (n > size(lat)) error stop 'read_taup_time_file: Input arrays are too short for file'
   enddo
   close(10)
   if (debug) write(0,'(a,i0.1,a)') 'read_taup_time_file: Got ', n, &
      ' points from file "' // trim(file) // '"'
end subroutine read_taup_time_file


end module tomo_predict
