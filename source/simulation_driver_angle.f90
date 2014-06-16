module vars_and_funcs
 integer, parameter :: nn = 1e6,                                               &
                       repetitions = 50
 real,    parameter ::                                                         &
   working_temp           = 51.5,                                              &
   idler_angle            = 2.3,                                               &
   min_angle              = 1.5,                                               &
   max_angle              = 3.0,                                               &
   angle_step             = 0.025,                                             &
   axial_pm_temp          = 65.0,                                              &
   crystal_dist           = 0.5,                                               &
   pump_power             = 0.030,                                             &
   pump_wavelength        = 406.118e-9,                                        &
   low_wavelength_limit   = 805.0e-9,                                          &
   high_wavelength_limit  = 820.0e-9,                                          &
   beam_waist             = 50.0e-6,                   &!!!!GET THESE NUMBERS RIGHT
   spectral_width_nm      = 0.035e-9,                  &!!!!GET THESE NUMBERS RIGHT
   crystal_length         = 0.005,                                             &
   d_eff                  = 1.0,                       &!!!!GET THESE NUMBERS RIGHT
   pi                     = 4.0*atan(1.0),                                     &
   c                      = 3.0e8
 real, parameter :: omega_pump     = 2.0*pi*c/pump_wavelength,                 &
                    spectral_width = 2.0*pi*c*spectral_width_nm/pump_wavelength**2
 real :: poling_period
 contains
real function temp_factor(temp)
real :: temp
real, parameter :: alpha=6.7e-6, beta=11.0e-9
 temp_factor=1.0+alpha*(temp-25.0)+beta*(temp-25)*(temp-25)
end function
real function crystal_length_temp()
 crystal_length_temp=crystal_length*temp_factor(temp)
end function
real function poling_period_temp()
 poling_period_temp=poling_period*temp_factor(temp)
end function

real function ktp_index(x,temp)
real :: x, temp
real :: lx
real, parameter :: na=2.12725, nb=1.18431, nc=5.14852e-2,&
                   nd=0.6603, ne=100.00506, nf=9.68956e-3
real, parameter :: a0=9.9587e-6, a1=9.9228e-6, a2=-8.9603e-6, a3=4.1010e-6
real, parameter :: b0=-1.1882e-8, b1=10.459e-8, b2=-9.8136e-8, b3=3.1481e-8
lx = x*1.0e6
ktp_index = sqrt( na + nb/(1.0-nc/(lx**2))+ nd/(1.0-ne/(lx**2)) - nf*(lx**2.0))&
      + (a0 + a1/lx + a2/(lx**2) + a3/(lx**3))*(temp-25.0)&
      + (b0 + b1/lx + b2/(lx**2) + b3/(lx**3))*(temp-25.0)**2.0
end function
 real function ktp_index_new(lambda,temp)
 implicit none
 real :: lambda,ll,n,n1,n2,temp
 real :: A=2.12725, B=1.18431, CC=5.14852e-2, D=0.6603
 real :: E= 100.00506, F=9.68956e-3
  ll  = lambda*1.0e6
  n1  = n_one(lambda)
  n2  = n_two(lambda)
  n   = sqrt(A + B/(1.0- CC/(ll**2)) + D/(1.0-E/(ll**2)) -F*ll**2)
  ktp_index_new = n + n1*(temp-25.0) + n2*(temp-25.0)**2.0
 end function
 real function n_one(lambda)
 real ::lambda,ll
  ll  = lambda*1.0e6
  n_one = 9.9587e-6 + 9.9228e-6/(ll) - 8.9603e-6/((ll)**2) + 4.1010e-6/((ll)**3)
 end function n_one
 real function n_two(lambda)
 implicit none
 real ::lambda,ll
  ll  = lambda*1.0e6
  n_two = -1.1882e-8 + 10.459e-8/(ll) - 9.8136e-8/((ll)**2) + 3.1481e-8/((ll)**3)
 end function n_two
end module vars_and_funcs


program pair_detection_simulator_angle_dependence
use vars_and_funcs
use vsl_stream
implicit none
 integer, parameter ::                                                         &
   number_of_angles = 1+ceiling(abs(max_angle-min_angle)/angle_step)
 real :: k_a, k_b, signal_gamma, idler_gamma,                                  &
         signal_aperture_angle, idler_aperture_angle
 real :: signal_aperture, idler_aperture, k_degen
 real, dimension(nn) :: signal_k, signal_polar, signal_azimuth,                &
                        idler_k,  idler_polar,  idler_azimuth,                 &
                        signal_kx, signal_ky, signal_kz,                       &
                        idler_kx,  idler_ky,  idler_kz,                        &
                        wavefunction
 real, dimension(number_of_angles) :: k_hypervolume, average,                  &
                                     point_value, phase_mismatch

 integer(kind=4) :: i = 1, j = 1
 real    :: temp = working_temp, angle = min_angle

 call read_parameters()
 call initialize_stream()

 average(:)=0.0

 do i = 1, number_of_angles

 write(*,*) 'angle ', i, ' of ',number_of_angles

 do j = 1, repetitions

 write(*,*) 'repetition ', j, ' of ',repetitions

 call set_variables()

 call determine_k_limits_and_hypervolume(k_a, k_b, k_hypervolume(i), k_degen)

 call generate_random_photons(k_a, k_b, signal_aperture_angle,                 &
                          signal_k, signal_polar, signal_azimuth, nn, k_degen)
 call generate_random_photons(k_a, k_b, idler_aperture_angle,                  &
                          idler_k,  idler_polar,  idler_azimuth,  nn, k_degen)

 call rotate_to_aperture_angle(signal_gamma,signal_k,signal_polar,signal_azimuth,&
                  signal_kx,signal_ky,signal_kz,nn)
 call rotate_to_aperture_angle(idler_gamma, idler_k, idler_polar, idler_azimuth, &
                  idler_kx, idler_ky, idler_kz, nn)

 call evaluate_wavefunction(signal_k, signal_kx, signal_ky, signal_kz,         &
                            idler_k,  idler_kx,  idler_ky,  idler_kz,          &
                            wavefunction, temp, k_a, k_b,                      &
                            point_value(i), phase_mismatch(i))!!!

 average(i) = average(i)+sum(wavefunction)/nn; wavefunction(:)=0.

 end do

 average(i) = average(i)/repetitions
 angle = min_angle+(i)*angle_step

 end do

 call write_data_to_file(min_angle, angle_step, number_of_angles,                 &
                         average , phase_mismatch, k_hypervolume, point_value)
 call deinitialize_stream()

contains

 subroutine read_parameters()
 implicit none
 character(len=30) :: value
  if(command_argument_count().eq.0)then
   write(*,*)'working with defaults, gotta give me 3 or 4 arguments like this:'
   write(*,*)'./coherence signal_aperture[mm] idler_aperture[mm]'
   write(*,*)'[defaults]./coherence 1.0 1.0'
   signal_aperture = 0.001
   idler_aperture  = 0.001
  elseif(command_argument_count().eq.2)then
   call get_command_argument(1,value)
   read(value,*) signal_aperture; signal_aperture=signal_aperture/1e3
   call get_command_argument(2,value)
   read(value,*) idler_aperture; idler_aperture=idler_aperture/1e3
  else
   write(*,*) 'wrong number of arguments, try no arguments to see some instructions'
   stop
  endif
  poling_period = pump_wavelength/(ktp_index(pump_wavelength,axial_pm_temp)    &
                                  -ktp_index(2.0*pump_wavelength,axial_pm_temp))
  poling_period = poling_period/temp_factor(axial_pm_temp)
write(*,*)'poling_period',poling_period
 end subroutine read_parameters

 subroutine set_variables()
 implicit none
  signal_aperture_angle = &
   tan(0.5*signal_aperture*cos(angle*pi/180.0)/crystal_dist)                   &
    /ktp_index(pump_wavelength*2.0,temp)
  idler_aperture_angle  = &
   tan(0.5*idler_aperture*cos(idler_angle*pi/180.0)/crystal_dist)              &
    /ktp_index(pump_wavelength*2.0,temp)
  signal_gamma  = asin(sin(angle*(pi/180.0))                                   &
                      /ktp_index(pump_wavelength*2.0,temp))
  idler_gamma   =-asin(sin(idler_angle*(pi/180.0))                             &
                      /ktp_index(pump_wavelength*2.0,temp))
 end subroutine set_variables

 subroutine determine_k_limits_and_hypervolume(k_low,k_high,hypervolume,k_degen)
 implicit none
 real, intent(out) :: k_low, k_high, hypervolume, k_degen
  k_degen = 2.0*pi*ktp_index(2.0*pump_wavelength,temp)/(2.0*pump_wavelength)
  k_low  = min( &
           ktp_index(high_wavelength_limit,temp)*2.0*pi/high_wavelength_limit, &
           ktp_index(low_wavelength_limit,temp)*2.0*pi/low_wavelength_limit)
  k_high = max( &
           ktp_index(high_wavelength_limit,temp)*2.0*pi/high_wavelength_limit, &
           ktp_index(low_wavelength_limit,temp)*2.0*pi/low_wavelength_limit)
 hypervolume = (2.0*pi/3.0)**2*(k_b**3-k_a**3)**2                              &
             * (1.0-cos(signal_aperture_angle))*(1.0-cos(idler_aperture_angle))
 end subroutine determine_k_limits_and_hypervolume

 subroutine write_data_to_file(min_angle, angle_step, number_of_angles,        &
                            average, phase_mismatch, k_hypervolume, point_value)
 implicit none
 integer, intent(in) :: number_of_angles
 real, intent(in)    :: min_angle, angle_step,                                 &
                   average(number_of_angles), phase_mismatch(number_of_angles),&
                   k_hypervolume(number_of_angles), point_value(number_of_angles)
 integer i
 character(80) fmt_list
 open(unit=10,file='coherence.dat')
 write(10,*) "#external_angle[degrees] 2nd_order_coherence phase_mismatch[pi] point_value"
 fmt_list = '(E12.4, E12.4, E12.4, E12.4)'
 do i=1,number_of_angles
 write(10,fmt_list) min_angle+(i-1)*angle_step, average(i)*k_hypervolume(i),   &
                    phase_mismatch(i)/pi, point_value(i)
 end do
 close(10)
 end subroutine
end program pair_detection_simulator_angle_dependence
