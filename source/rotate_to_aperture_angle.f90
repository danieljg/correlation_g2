subroutine rotate_to_aperture_angle( rotation_angle, length_k, polar_angle,      &
                             azimuth_angle, kx, ky, kz, nn )
implicit none
 real, intent(in) :: rotation_angle
 integer, intent(in) :: nn
 real, intent(in)  :: length_k(nn), polar_angle(nn), azimuth_angle(nn)
 real, intent(out) :: kx(nn), ky(nn), kz(nn)
 real :: x(nn), z(nn)

 x(:)=length_k(:)*sin(polar_angle(:))*cos(azimuth_angle(:))
 ky(:)=length_k(:)*sin(polar_angle(:))*sin(azimuth_angle(:))
 z(:)=length_k(:)*cos(polar_angle(:))

 kx(:)=x(:)*cos(rotation_angle)+z(:)*sin(rotation_angle)
 kz(:)=-x(:)*sin(rotation_angle)+z(:)*cos(rotation_angle)

end subroutine rotate_to_aperture_angle
