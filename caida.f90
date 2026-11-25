! Proyecto programado 1 - IC4700 Lenguajes de Programación
! Caída amortiguada de un paracaidista con Runge-Kutta de orden 4.
! Versión en Fortran 2008: módulo con los módulos del diseño (Parte 2).
! Un módulo del diseño = una función o subrutina de este archivo.
! Lo usan paracaidas.f90 (programa) y prueba_rk4.f90 (caso de prueba 1).

module caida
  use, intrinsic :: iso_fortran_env, only: iostat_end, error_unit
  implicit none

  integer, parameter :: dp = kind(1.0d0)   ! real(kind=8): doble precisión

  ! ---- Constantes de lectura (Sección 9 del diseño) ----
  real(dp), parameter :: G          = 9.81_dp   ! gravedad (m/s^2)
  real(dp), parameter :: H_PASO     = 0.01_dp   ! paso h de RK4 (s)
  real(dp), parameter :: V_OBJ_KMH  = 8.0_dp    ! velocidad de llegada deseada
  real(dp), parameter :: FACTOR_KMH = 3.6_dp    ! conversión m/s <-> km/h
  real(dp), parameter :: TOL_BETA   = 1.0e-6_dp ! tolerancia de la bisección
  integer,  parameter :: MAX_ITER   = 100       ! tope de iteraciones

contains

  ! ---- Figura 2: leer_valor_positivo(mensaje) -> valor ----
  function leer_valor_positivo(mensaje) result(valor)
    character(len=*), intent(in) :: mensaje
    real(dp) :: valor
    integer  :: estado
    logical  :: valido

    ! Fortran no tiene REPETIR...HASTA: se usa do while con variable lógica.
    valido = .false.
    do while (.not. valido)
      write (*, '(a, a)', advance='no') mensaje, ': '
      read (*, *, iostat=estado) valor
      if (estado == iostat_end) then      ! sin más datos: evita ciclo infinito
        write (*, *)
        write (error_unit, '(a)') 'No hay más datos de entrada.'
        stop 1
      end if
      if (estado /= 0) then               ! no era un número
        valor = 0.0_dp
      end if
      valido = (valor > 0.0_dp)
      if (.not. valido) then
        write (*, '(a)') 'El valor debe ser un número positivo.'
      end if
    end do
  end function leer_valor_positivo

  ! ---- Figura 3: calcular_peso(m_paracaidista, m_impedimenta) -> W ----
  function calcular_peso(m_paracaidista, m_impedimenta) result(W)
    real(dp), intent(in) :: m_paracaidista, m_impedimenta
    real(dp) :: W

    W = (m_paracaidista + m_impedimenta) * G
  end function calcular_peso

  ! ---- Figura 4: altura_minima_viable(v_obj) -> altura_min ----
  function altura_minima_viable(v_obj) result(altura_min)
    real(dp), intent(in) :: v_obj
    real(dp) :: altura_min

    altura_min = (v_obj * v_obj) / (2.0_dp * G)
  end function altura_minima_viable

  ! ---- Figura 7: calcular_c(beta, W) -> c ----
  function calcular_c(beta, W) result(c)
    real(dp), intent(in) :: beta, W
    real(dp) :: c

    c = beta * G / W
  end function calcular_c

  ! ---- Figura 10: derivada_velocidad(v, c) -> dvdt ----
  function derivada_velocidad(v, c) result(dvdt)
    real(dp), intent(in) :: v, c
    real(dp) :: dvdt

    dvdt = G - c * v
  end function derivada_velocidad

  ! ---- Figura 9: paso_rk4(v, x, c, h) -> v_nuevo, x_nuevo ----
  subroutine paso_rk4(v, x, c, h, v_nuevo, x_nuevo)
    real(dp), intent(in)  :: v, x, c, h
    real(dp), intent(out) :: v_nuevo, x_nuevo
    real(dp) :: k1v, k1x, k2v, k2x, k3v, k3x, k4v, k4x

    k1v = h * derivada_velocidad(v, c)
    k1x = h * v
    k2v = h * derivada_velocidad(v + k1v / 2.0_dp, c)
    k2x = h * (v + k1v / 2.0_dp)
    k3v = h * derivada_velocidad(v + k2v / 2.0_dp, c)
    k3x = h * (v + k2v / 2.0_dp)
    k4v = h * derivada_velocidad(v + k3v, c)
    k4x = h * (v + k3v)

    v_nuevo = v + (k1v + 2.0_dp * k2v + 2.0_dp * k3v + k4v) / 6.0_dp
    x_nuevo = x + (k1x + 2.0_dp * k2x + 2.0_dp * k3x + k4x) / 6.0_dp
  end subroutine paso_rk4

  ! ---- Figura 11: fraccion_hasta_suelo(x, x_nuevo, altura) -> s ----
  function fraccion_hasta_suelo(x, x_nuevo, altura) result(s)
    real(dp), intent(in) :: x, x_nuevo, altura
    real(dp) :: s

    s = (altura - x) / (x_nuevo - x)
  end function fraccion_hasta_suelo

  ! ---- Figura 12: interpolar_lineal(a, b, s) -> valor ----
  function interpolar_lineal(a, b, s) result(valor)
    real(dp), intent(in) :: a, b, s
    real(dp) :: valor

    valor = a + s * (b - a)
  end function interpolar_lineal

  ! ---- Figura 8: simular_caida(c, h, altura) -> t_llegada, v_llegada ----
  subroutine simular_caida(c, h, altura, t_llegada, v_llegada)
    real(dp), intent(in)  :: c, h, altura
    real(dp), intent(out) :: t_llegada, v_llegada
    real(dp) :: t, v, x, v_nuevo, x_nuevo, s

    t = 0.0_dp
    v = 0.0_dp
    x = 0.0_dp
    t_llegada = 0.0_dp
    v_llegada = 0.0_dp

    do while (x < altura)
      call paso_rk4(v, x, c, h, v_nuevo, x_nuevo)
      if (x_nuevo >= altura) then
        s = fraccion_hasta_suelo(x, x_nuevo, altura)
        t_llegada = interpolar_lineal(t, t + h, s)
        v_llegada = interpolar_lineal(v, v_nuevo, s)
      end if
      t = t + h
      v = v_nuevo
      x = x_nuevo
    end do
  end subroutine simular_caida

  ! ---- Figura 6: velocidad_llegada(beta, W, h, altura) -> v_llegada ----
  function velocidad_llegada(beta, W, h, altura) result(v_llegada)
    real(dp), intent(in) :: beta, W, h, altura
    real(dp) :: v_llegada
    real(dp) :: c, t_llegada

    c = calcular_c(beta, W)
    call simular_caida(c, h, altura, t_llegada, v_llegada)
  end function velocidad_llegada

  ! ---- Figura 5: buscar_beta(W, h, altura, v_obj, tol) -> beta ----
  function buscar_beta(W, h, altura, v_obj, tol) result(beta)
    real(dp), intent(in) :: W, h, altura, v_obj, tol
    real(dp) :: beta
    real(dp) :: beta_bajo, beta_alto, beta_medio, v
    integer  :: iteraciones

    beta_bajo = 0.0_dp
    beta_alto = W
    iteraciones = 0

    do while ((beta_alto - beta_bajo > tol) .and. (iteraciones < MAX_ITER))
      beta_medio = (beta_bajo + beta_alto) / 2.0_dp
      v = velocidad_llegada(beta_medio, W, h, altura)
      if (v > v_obj) then
        beta_bajo = beta_medio      ! llega muy rápido: falta resistencia
      else
        beta_alto = beta_medio
      end if
      iteraciones = iteraciones + 1
    end do
    beta = (beta_bajo + beta_alto) / 2.0_dp
  end function buscar_beta

  ! ---- Figura 13: mostrar_resultados(W, beta, t_llegada, v_llegada) ----
  subroutine mostrar_resultados(W, beta, t_llegada, v_llegada)
    real(dp), intent(in) :: W, beta, t_llegada, v_llegada

    write (*, *)
    write (*, '(a)') '--- Resultados ---'
    write (*, '(a, f0.4)') 'Peso combinado W (N): ', W
    write (*, '(a, f0.6)') 'Constante de resistencia beta (N·s/m): ', beta
    write (*, '(a, f0.6)') 'Velocidad límite W / beta (m/s): ', W / beta
    write (*, '(a, f0.4)') 'Tiempo de caída (s): ', t_llegada
    write (*, '(a, f0.6)') 'Velocidad de llegada (m/s): ', v_llegada
    write (*, '(a, f0.6)') 'Velocidad de llegada (km/h): ', v_llegada * FACTOR_KMH
  end subroutine mostrar_resultados

end module caida
