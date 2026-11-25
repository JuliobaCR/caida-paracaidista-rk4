# Caída amortiguada de un paracaidista (RK4)

Proyecto programado 1 de IC4700 Lenguajes de Programación (II semestre 2026).

El programa resuelve el modelo (W/g)·dv/dt = W − βv con Runge-Kutta de cuarto orden. Además busca por bisección la constante β del paracaídas con la que se llega al suelo a 8 km/h. Las tres versiones (C, Pascal y Fortran) traducen los mismos diagramas Nassi-Shneiderman del documento de diseño, módulo por módulo.

## Resultado

Para 80 kg + 50 kg saltando desde 1 500 m:

```
Peso combinado W (N): 1275.3000
Constante de resistencia beta (N·s/m): 573.885000
Velocidad límite W / beta (m/s): 2.222222
Tiempo de caída (s): 675.2265
Velocidad de llegada (m/s): 2.222222
Velocidad de llegada (km/h): 8.000000
```

Las tres versiones dan exactamente la misma salida.

## Archivos

| Archivo | Contenido |
|---|---|
| `paracaidas.c` | Versión en C (los 13 módulos) |
| `caida.pas`, `paracaidas.pas` | Versión en Pascal: unidad con los módulos + programa principal |
| `caida.f90`, `paracaidas.f90` | Versión en Fortran: módulo + programa principal |
| `prueba_rk4.c`, `.pas`, `.f90` | Caso de prueba 1: RK4 contra la solución analítica |
| `docs/Proyecto1_Diseno.docx` | Documento de diseño e informe |

## Compilar y ejecutar

```sh
# C
gcc -O2 -Wall -Wextra -std=c99 -pedantic -o paracaidas_c paracaidas.c -lm
gcc -O2 -Wall -Wextra -std=c99 -pedantic -o prueba_c prueba_rk4.c -lm

# Pascal (Free Pascal compila la unidad caida.pas automáticamente)
fpc -O2 -oparacaidas_pas paracaidas.pas
fpc -O2 -oprueba_pas prueba_rk4.pas

# Fortran
gfortran -O2 -Wall -Wextra -std=f2008 -o paracaidas_f90 caida.f90 paracaidas.f90
gfortran -O2 -Wall -Wextra -std=f2008 -o prueba_f90 caida.f90 prueba_rk4.f90
```

Al ejecutar `paracaidas_*`, el programa pide el peso del paracaidista, el de la impedimenta y la altura del salto.

Se probó con gcc 11.4, Free Pascal 3.2.2 y gfortran 11.4 (Ubuntu 22.04); las tres versiones compilan sin advertencias.

## Referencias

- Spiegel, M. *Ecuaciones Diferenciales Aplicadas*. Prentice-Hall, 1983, pp. 76-78.
- Burden, R. L.; Faires, J. D. *Análisis Numérico*. Grupo Editorial Iberoamérica, 1985, p. 248.
