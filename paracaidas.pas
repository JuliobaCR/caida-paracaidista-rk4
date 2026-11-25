{
  Proyecto programado 1 - IC4700 Lenguajes de Programación
  Caída amortiguada de un paracaidista con Runge-Kutta de orden 4.
  Versión en Pascal. Sigue los diagramas Nassi-Shneiderman del documento de
  diseño (Parte 2). Este archivo es el módulo Principal (Figura 1); los demás
  módulos están en la unidad caida.pas.

  Compilar: fpc -O2 paracaidas.pas      (compila también caida.pas)
  Ejecutar: ./paracaidas
}
program Paracaidas;

uses
  caida;

{ ---- Figura 1: Principal ---- }
var
  v_obj, m_paracaidista, m_impedimenta, altura, W: Double;
  beta, c, t_llegada, v_llegada: Double;
begin
  v_obj := V_OBJ_KMH / FACTOR_KMH;
  m_paracaidista := leer_valor_positivo('Peso del paracaidista (kg)');
  m_impedimenta  := leer_valor_positivo('Peso de la impedimenta (kg)');
  altura         := leer_valor_positivo('Altura de salto (m)');
  W := calcular_peso(m_paracaidista, m_impedimenta);

  if altura <= altura_minima_viable(v_obj) then
    writeln('La altura es insuficiente para alcanzar la velocidad objetivo')
  else
  begin
    beta := buscar_beta(W, H_PASO, altura, v_obj, TOL_BETA);
    c := calcular_c(beta, W);
    simular_caida(c, H_PASO, altura, t_llegada, v_llegada);
    mostrar_resultados(W, beta, t_llegada, v_llegada);
  end;
end.
