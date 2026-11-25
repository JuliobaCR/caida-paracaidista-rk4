! Proyecto programado 1 - IC4700 Lenguajes de Programación
! Caída amortiguada de un paracaidista con Runge-Kutta de orden 4.
! Versión en Fortran 2008. Sigue los diagramas Nassi-Shneiderman del
! documento de diseño (Parte 2). Este archivo es el módulo Principal
! (Figura 1); los demás módulos están en caida.f90.
!
! Compilar: gfortran -O2 -Wall -Wextra -std=f2008 -o paracaidas caida.f90 paracaidas.f90
! Ejecutar: ./paracaidas

! ---- Figura 1: Principal ----
program paracaidas
  use caida
  implicit none

  real(dp) :: v_obj, m_paracaidista, m_impedimenta, altura, W
  real(dp) :: beta, c, t_llegada, v_llegada

  v_obj = V_OBJ_KMH / FACTOR_KMH
  m_paracaidista = leer_valor_positivo('Peso del paracaidista (kg)')
  m_impedimenta  = leer_valor_positivo('Peso de la impedimenta (kg)')
  altura         = leer_valor_positivo('Altura de salto (m)')
  W = calcular_peso(m_paracaidista, m_impedimenta)

  if (altura <= altura_minima_viable(v_obj)) then
    write (*, '(a)') 'La altura es insuficiente para alcanzar la velocidad objetivo'
  else
    beta = buscar_beta(W, H_PASO, altura, v_obj, TOL_BETA)
    c = calcular_c(beta, W)
    call simular_caida(c, H_PASO, altura, t_llegada, v_llegada)
    call mostrar_resultados(W, beta, t_llegada, v_llegada)
  end if
end program paracaidas
