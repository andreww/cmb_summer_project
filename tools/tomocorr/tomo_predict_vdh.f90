program tomo_predict_vdh
   use tomo_predict
   implicit none
   character(len=250) :: model1d, model3d, infile
   real :: ztop, zbot
   ! Maximum number of ray points
   integer, parameter :: npmax = 10000
   real, dimension(npmax) :: lat, lon, r, z
   real :: dt
   integer :: n

   read(*,*) model1d
   read(*,*) model3d
   read(*,*) infile
   read(*,*) ztop, zbot ! These are ignored in the updated version

   call read_taup_time_file(infile, lat, lon, r, n)

   z = 6371. - r
   call predict(model1d, model3d, lat(1:n), lon(1:n), z(1:n), dt)
   write(*,*) dt
end program tomo_predict_vdh
