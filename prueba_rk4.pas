{ Caso de prueba 1 del documento de diseño: RK4 vs solución analítica.
  Compilar: fpc -O2 prueba_rk4.pas      (usa la unidad caida.pas)
}
program prueba_rk4;

uses
  caida;

const
  BETA = 573.885;

var
  W, c, v, x, v_nuevo, x_nuevo, h: Double;
  pasos: array[0..1] of Double = (0.01, 0.1);
  i, k, n: Integer;
begin
  W := calcular_peso(80.0, 50.0);
  c := calcular_c(BETA, W);

  writeln('Analítica  (t = 1 s): v = ', ((W / BETA) * (1.0 - exp(-c))):0:6,
          '  x = ', ((W / BETA) * (1.0 + (W / (BETA * G)) * exp(-c)
                     - W / (BETA * G))):0:6);

  for k := 0 to 1 do
  begin
    h := pasos[k];
    n := round(1.0 / h);
    v := 0.0;
    x := 0.0;
    for i := 1 to n do
    begin
      paso_rk4(v, x, c, h, v_nuevo, x_nuevo);
      v := v_nuevo;
      x := x_nuevo;
    end;
    writeln('RK4, h = ', h:0:2, ' (t = 1 s): v = ', v:0:6, '  x = ', x:0:6);
  end;
end.
