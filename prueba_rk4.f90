! Caso de prueba 1 del documento de diseño: RK4 vs solución analítica.
! Compilar: gfortran -O2 -Wall -Wextra -std=f2008 -o prueba_rk4 caida.f90 prueba_rk4.f90
program prueba_rk4
  use caida
  implicit none

  real(dp), parameter :: BETA = 573.885_dp
  real(dp), parameter :: PASOS(2) = [0.01_dp, 0.1_dp]
  real(dp) :: W, c, v, x, v_nuevo, x_nuevo, h
  integer  :: i, k, n

  W = calcular_peso(80.0_dp, 50.0_dp)
  c = calcular_c(BETA, W)

  write (*, '(a, f0.6, a, f0.6)') 'Analítica  (t = 1 s): v = ', &
    (W / BETA) * (1.0_dp - exp(-c)), '  x = ', &
    (W / BETA) * (1.0_dp + (W / (BETA * G)) * exp(-c) - W / (BETA * G))

  do k = 1, 2
    h = PASOS(k)
    n = nint(1.0_dp / h)
    v = 0.0_dp
    x = 0.0_dp
    do i = 1, n
      call paso_rk4(v, x, c, h, v_nuevo, x_nuevo)
      v = v_nuevo
      x = x_nuevo
    end do
    write (*, '(a, f4.2, a, f0.6, a, f0.6)') 'RK4, h = ', h, &
      ' (t = 1 s): v = ', v, '  x = ', x
  end do
end program prueba_rk4
