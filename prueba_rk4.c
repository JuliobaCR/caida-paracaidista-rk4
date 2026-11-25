/* Caso de prueba 1 del documento de diseño: RK4 vs solución analítica.
 * Compilar: gcc -O2 -Wall -Wextra -std=c99 -pedantic -o prueba_rk4 prueba_rk4.c -lm
 */
#define SIN_MAIN
#include "paracaidas.c"

int main(void)
{
    const double beta = 573.885;
    const double W = calcular_peso(80.0, 50.0);
    const double c = calcular_c(beta, W);
    const double pasos[2] = { 0.01, 0.1 };
    double v_nuevo, x_nuevo, v, x, h;
    int i, k, n;

    printf("Analítica  (t = 1 s): v = %.6f  x = %.6f\n",
           (W / beta) * (1.0 - exp(-c)),
           (W / beta) * (1.0 + (W / (beta * G)) * exp(-c) - W / (beta * G)));

    for (k = 0; k < 2; k++) {
        h = pasos[k];
        n = (int)(1.0 / h + 0.5);
        v = 0.0;
        x = 0.0;
        for (i = 0; i < n; i++) {
            paso_rk4(v, x, c, h, &v_nuevo, &x_nuevo);
            v = v_nuevo;
            x = x_nuevo;
        }
        printf("RK4, h = %.2f (t = 1 s): v = %.6f  x = %.6f\n", h, v, x);
    }
    return 0;
}
