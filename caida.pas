{
  Proyecto programado 1 - IC4700 Lenguajes de Programación
  Caída amortiguada de un paracaidista con Runge-Kutta de orden 4.
  Versión en Pascal: unidad con los módulos del diseño (Parte 2).
  Un módulo del diseño = una función o procedimiento de esta unidad.
  La usan paracaidas.pas (programa) y prueba_rk4.pas (caso de prueba 1).
}
unit caida;

interface

{ ---- Constantes de lectura (Sección 9 del diseño) ---- }
const
  G          = 9.81;   { gravedad (m/s^2)                 }
  H_PASO     = 0.01;   { paso h de RK4 (s)                }
  V_OBJ_KMH  = 8.0;    { velocidad de llegada deseada     }
  FACTOR_KMH = 3.6;    { conversión m/s <-> km/h          }
  TOL_BETA   = 1e-6;   { tolerancia de la bisección       }
  MAX_ITER   = 100;    { tope de iteraciones              }

function leer_valor_positivo(mensaje: string): Double;
function calcular_peso(m_paracaidista, m_impedimenta: Double): Double;
function altura_minima_viable(v_obj: Double): Double;
function calcular_c(beta, W: Double): Double;
function derivada_velocidad(v, c: Double): Double;
procedure paso_rk4(v, x, c, h: Double; var v_nuevo, x_nuevo: Double);
function fraccion_hasta_suelo(x, x_nuevo, altura: Double): Double;
function interpolar_lineal(a, b, s: Double): Double;
procedure simular_caida(c, h, altura: Double; var t_llegada, v_llegada: Double);
function velocidad_llegada(beta, W, h, altura: Double): Double;
function buscar_beta(W, h, altura, v_obj, tol: Double): Double;
procedure mostrar_resultados(W, beta, t_llegada, v_llegada: Double);

implementation

{ ---- Figura 2: leer_valor_positivo(mensaje) -> valor ---- }
function leer_valor_positivo(mensaje: string): Double;
var
  texto: string;
  valor: Double;
  error: Integer;
  valido: Boolean;
begin
  repeat
    write(mensaje, ': ');
    if eof(input) then                  { sin más datos: evita ciclo infinito }
    begin
      writeln;
      writeln(stderr, 'No hay más datos de entrada.');
      halt(1);
    end;
    readln(texto);
    { Quitar espacios, tabuladores y retornos de carro de los extremos. }
    while (length(texto) > 0) and (texto[length(texto)] in [' ', #9, #13]) do
      delete(texto, length(texto), 1);
    while (length(texto) > 0) and (texto[1] in [' ', #9]) do
      delete(texto, 1, 1);
    val(texto, valor, error);
    if (error <> 0) or (texto = '') then   { no era un número }
      valor := 0.0;
    valido := valor > 0.0;
    if not valido then
      writeln('El valor debe ser un número positivo.');
  until valido;
  leer_valor_positivo := valor;
end;

{ ---- Figura 3: calcular_peso(m_paracaidista, m_impedimenta) -> W ---- }
function calcular_peso(m_paracaidista, m_impedimenta: Double): Double;
begin
  calcular_peso := (m_paracaidista + m_impedimenta) * G;
end;

{ ---- Figura 4: altura_minima_viable(v_obj) -> altura_min ---- }
function altura_minima_viable(v_obj: Double): Double;
begin
  altura_minima_viable := (v_obj * v_obj) / (2.0 * G);
end;

{ ---- Figura 7: calcular_c(beta, W) -> c ---- }
function calcular_c(beta, W: Double): Double;
begin
  calcular_c := beta * G / W;
end;

{ ---- Figura 10: derivada_velocidad(v, c) -> dvdt ---- }
function derivada_velocidad(v, c: Double): Double;
begin
  derivada_velocidad := G - c * v;
end;

{ ---- Figura 9: paso_rk4(v, x, c, h) -> v_nuevo, x_nuevo ---- }
procedure paso_rk4(v, x, c, h: Double; var v_nuevo, x_nuevo: Double);
var
  k1v, k1x, k2v, k2x, k3v, k3x, k4v, k4x: Double;
begin
  k1v := h * derivada_velocidad(v, c);
  k1x := h * v;
  k2v := h * derivada_velocidad(v + k1v / 2.0, c);
  k2x := h * (v + k1v / 2.0);
  k3v := h * derivada_velocidad(v + k2v / 2.0, c);
  k3x := h * (v + k2v / 2.0);
  k4v := h * derivada_velocidad(v + k3v, c);
  k4x := h * (v + k3v);

  v_nuevo := v + (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0;
  x_nuevo := x + (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
end;

{ ---- Figura 11: fraccion_hasta_suelo(x, x_nuevo, altura) -> s ---- }
function fraccion_hasta_suelo(x, x_nuevo, altura: Double): Double;
begin
  fraccion_hasta_suelo := (altura - x) / (x_nuevo - x);
end;

{ ---- Figura 12: interpolar_lineal(a, b, s) -> valor ---- }
function interpolar_lineal(a, b, s: Double): Double;
begin
  interpolar_lineal := a + s * (b - a);
end;

{ ---- Figura 8: simular_caida(c, h, altura) -> t_llegada, v_llegada ---- }
procedure simular_caida(c, h, altura: Double; var t_llegada, v_llegada: Double);
var
  t, v, x, v_nuevo, x_nuevo, s: Double;
begin
  t := 0.0;  v := 0.0;  x := 0.0;
  v_nuevo := 0.0;  x_nuevo := 0.0;
  t_llegada := 0.0;
  v_llegada := 0.0;

  while x < altura do
  begin
    paso_rk4(v, x, c, h, v_nuevo, x_nuevo);
    if x_nuevo >= altura then
    begin
      s := fraccion_hasta_suelo(x, x_nuevo, altura);
      t_llegada := interpolar_lineal(t, t + h, s);
      v_llegada := interpolar_lineal(v, v_nuevo, s);
    end;
    t := t + h;
    v := v_nuevo;
    x := x_nuevo;
  end;
end;

{ ---- Figura 6: velocidad_llegada(beta, W, h, altura) -> v_llegada ---- }
function velocidad_llegada(beta, W, h, altura: Double): Double;
var
  c, t_llegada, v_llegada: Double;
begin
  c := calcular_c(beta, W);
  simular_caida(c, h, altura, t_llegada, v_llegada);
  velocidad_llegada := v_llegada;
end;

{ ---- Figura 5: buscar_beta(W, h, altura, v_obj, tol) -> beta ---- }
function buscar_beta(W, h, altura, v_obj, tol: Double): Double;
var
  beta_bajo, beta_alto, beta_medio, v: Double;
  iteraciones: Integer;
begin
  beta_bajo := 0.0;
  beta_alto := W;
  iteraciones := 0;

  while (beta_alto - beta_bajo > tol) and (iteraciones < MAX_ITER) do
  begin
    beta_medio := (beta_bajo + beta_alto) / 2.0;
    v := velocidad_llegada(beta_medio, W, h, altura);
    if v > v_obj then
      beta_bajo := beta_medio      { llega muy rápido: falta resistencia }
    else
      beta_alto := beta_medio;
    iteraciones := iteraciones + 1;
  end;
  buscar_beta := (beta_bajo + beta_alto) / 2.0;
end;

{ ---- Figura 13: mostrar_resultados(W, beta, t_llegada, v_llegada) ---- }
procedure mostrar_resultados(W, beta, t_llegada, v_llegada: Double);
begin
  writeln;
  writeln('--- Resultados ---');
  writeln('Peso combinado W (N): ', W:0:4);
  writeln('Constante de resistencia beta (N·s/m): ', beta:0:6);
  writeln('Velocidad límite W / beta (m/s): ', (W / beta):0:6);
  writeln('Tiempo de caída (s): ', t_llegada:0:4);
  writeln('Velocidad de llegada (m/s): ', v_llegada:0:6);
  writeln('Velocidad de llegada (km/h): ', (v_llegada * FACTOR_KMH):0:6);
end;

end.
