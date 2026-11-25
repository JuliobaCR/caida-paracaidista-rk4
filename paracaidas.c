/*
 * Proyecto programado 1 - IC4700 Lenguajes de Programación
 * Caída amortiguada de un paracaidista con Runge-Kutta de orden 4.
 * Versión en C. Sigue los diagramas Nassi-Shneiderman del documento de diseño
 * (Parte 2): un módulo del diseño = una función de este archivo.
 *
 * Compilar: gcc -O2 -Wall -Wextra -std=c99 -pedantic -o paracaidas paracaidas.c
 * Ejecutar: ./paracaidas
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>   /* solo para isfinite() */

/* ---- Constantes de lectura (Sección 9 del diseño) ---- */
static const double G          = 9.81;   /* gravedad (m/s^2)                 */
static const double H_PASO     = 0.01;   /* paso h de RK4 (s)                */
static const double V_OBJ_KMH  = 8.0;    /* velocidad de llegada deseada     */
static const double FACTOR_KMH = 3.6;    /* conversión m/s <-> km/h          */
static const double TOL_BETA   = 1e-6;   /* tolerancia de la bisección       */
static const int    MAX_ITER   = 100;    /* tope de iteraciones              */

/* ---- Figura 2: leer_valor_positivo(mensaje) -> valor ---- */
double leer_valor_positivo(const char *mensaje)
{
    double valor = 0.0;
    int leidos;
    int valido;
    int c;

    do {
        printf("%s: ", mensaje);
        leidos = scanf("%lf", &valor);
        if (leidos == EOF) {                 /* sin más datos: evita ciclo infinito */
            fprintf(stderr, "\nNo hay más datos de entrada.\n");
            exit(EXIT_FAILURE);
        }
        if (leidos != 1) {                   /* no era un número: descartar la línea */
            while ((c = getchar()) != '\n' && c != EOF) { }
            valor = 0.0;
        }
        valido = (valor > 0.0 && isfinite(valor));
        if (!valido) {
            printf("El valor debe ser un número positivo.\n");
        }
    } while (!valido);

    return valor;
}

/* ---- Figura 3: calcular_peso(m_paracaidista, m_impedimenta) -> W ---- */
double calcular_peso(double m_paracaidista, double m_impedimenta)
{
    return (m_paracaidista + m_impedimenta) * G;
}

/* ---- Figura 4: altura_minima_viable(v_obj) -> altura_min ---- */
double altura_minima_viable(double v_obj)
{
    return (v_obj * v_obj) / (2.0 * G);
}

/* ---- Figura 7: calcular_c(beta, W) -> c ---- */
double calcular_c(double beta, double W)
{
    return beta * G / W;
}

/* ---- Figura 10: derivada_velocidad(v, c) -> dvdt ---- */
double derivada_velocidad(double v, double c)
{
    return G - c * v;
}

/* ---- Figura 9: paso_rk4(v, x, c, h) -> v_nuevo, x_nuevo ---- */
void paso_rk4(double v, double x, double c, double h,
              double *v_nuevo, double *x_nuevo)
{
    double k1v = h * derivada_velocidad(v, c);
    double k1x = h * v;
    double k2v = h * derivada_velocidad(v + k1v / 2.0, c);
    double k2x = h * (v + k1v / 2.0);
    double k3v = h * derivada_velocidad(v + k2v / 2.0, c);
    double k3x = h * (v + k2v / 2.0);
    double k4v = h * derivada_velocidad(v + k3v, c);
    double k4x = h * (v + k3v);

    *v_nuevo = v + (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0;
    *x_nuevo = x + (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0;
}

/* ---- Figura 11: fraccion_hasta_suelo(x, x_nuevo, altura) -> s ---- */
double fraccion_hasta_suelo(double x, double x_nuevo, double altura)
{
    return (altura - x) / (x_nuevo - x);
}

/* ---- Figura 12: interpolar_lineal(a, b, s) -> valor ---- */
double interpolar_lineal(double a, double b, double s)
{
    return a + s * (b - a);
}

/* ---- Figura 8: simular_caida(c, h, altura) -> t_llegada, v_llegada ---- */
void simular_caida(double c, double h, double altura,
                   double *t_llegada, double *v_llegada)
{
    double t = 0.0, v = 0.0, x = 0.0;
    double v_nuevo, x_nuevo, s;

    *t_llegada = 0.0;
    *v_llegada = 0.0;

    while (x < altura) {
        paso_rk4(v, x, c, h, &v_nuevo, &x_nuevo);
        if (x_nuevo >= altura) {
            s = fraccion_hasta_suelo(x, x_nuevo, altura);
            *t_llegada = interpolar_lineal(t, t + h, s);
            *v_llegada = interpolar_lineal(v, v_nuevo, s);
        }
        t = t + h;
        v = v_nuevo;
        x = x_nuevo;
    }
}

/* ---- Figura 6: velocidad_llegada(beta, W, h, altura) -> v_llegada ---- */
double velocidad_llegada(double beta, double W, double h, double altura)
{
    double c = calcular_c(beta, W);
    double t_llegada, v_llegada;

    simular_caida(c, h, altura, &t_llegada, &v_llegada);
    return v_llegada;
}

/* ---- Figura 5: buscar_beta(W, h, altura, v_obj, tol) -> beta ---- */
double buscar_beta(double W, double h, double altura, double v_obj, double tol)
{
    double beta_bajo = 0.0, beta_alto = W, beta_medio, v;
    int iteraciones = 0;

    while ((beta_alto - beta_bajo > tol) && (iteraciones < MAX_ITER)) {
        beta_medio = (beta_bajo + beta_alto) / 2.0;
        v = velocidad_llegada(beta_medio, W, h, altura);
        if (v > v_obj) {
            beta_bajo = beta_medio;      /* llega muy rápido: falta resistencia */
        } else {
            beta_alto = beta_medio;
        }
        iteraciones = iteraciones + 1;
    }
    return (beta_bajo + beta_alto) / 2.0;
}

/* ---- Figura 13: mostrar_resultados(W, beta, t_llegada, v_llegada) ---- */
void mostrar_resultados(double W, double beta, double t_llegada, double v_llegada)
{
    printf("\n--- Resultados ---\n");
    printf("Peso combinado W (N): %.4f\n", W);
    printf("Constante de resistencia beta (N·s/m): %.6f\n", beta);
    printf("Velocidad límite W / beta (m/s): %.6f\n", W / beta);
    printf("Tiempo de caída (s): %.4f\n", t_llegada);
    printf("Velocidad de llegada (m/s): %.6f\n", v_llegada);
    printf("Velocidad de llegada (km/h): %.6f\n", v_llegada * FACTOR_KMH);
}

#ifndef SIN_MAIN
/* ---- Figura 1: Principal ---- */
int main(void)
{
    double v_obj, m_paracaidista, m_impedimenta, altura, W;
    double beta, c, t_llegada, v_llegada;

    v_obj = V_OBJ_KMH / FACTOR_KMH;
    m_paracaidista = leer_valor_positivo("Peso del paracaidista (kg)");
    m_impedimenta  = leer_valor_positivo("Peso de la impedimenta (kg)");
    altura         = leer_valor_positivo("Altura de salto (m)");
    W = calcular_peso(m_paracaidista, m_impedimenta);

    if (altura <= altura_minima_viable(v_obj)) {
        printf("La altura es insuficiente para alcanzar la velocidad objetivo\n");
    } else {
        beta = buscar_beta(W, H_PASO, altura, v_obj, TOL_BETA);
        c = calcular_c(beta, W);
        simular_caida(c, H_PASO, altura, &t_llegada, &v_llegada);
        mostrar_resultados(W, beta, t_llegada, v_llegada);
    }
    return 0;
}
#endif
