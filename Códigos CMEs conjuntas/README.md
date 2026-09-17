# Sobre CME_combinadas.py

Simulación cinemática y morfológica de la interacción entre dos CMEs homólogas (CME-1 "precursora" y CME-2 "sucesora") en el medio interplanetario. El código propaga ambas eyecciones, calcula sus perfiles cinemáticos y sus campos de densidad y velocidad en coordenadas polares $(\theta, r)$, detecta el solapamiento entre sus perfiles, sintetiza el campo combinado en la zona de interacción (compresión) y, cuando el solapamiento es suficientemente denso, crea una nueva nube de puntos (producto de la fusión) que se propaga de forma independiente, con su propia velocidad y densidad. Genera 4 salidas: cinemática conjunta, propagación polar (8 paneles), series temporales multipunto y un resumen numérico.

## 1. Modelo físico

- **Aceleración de dos fases** (Gallagher et al. 2003): cada CME tiene una fase de aceleración impulsiva y una de frenado exponencial, combinadas en una sola función `aceleracion(s)`.
- **Morfología no circular**: el frente (r_ext) y la retaguardia (r_int) de cada CME no son arcos concéntricos, sino perfiles deformados aleatoriamente.
- **Densidad**: perfil radial tipo gaussiano centrado en el frente (mayor compresión ahí), modulado angularmente y decayendo con la expansión (conservación aproximada de masa).
- **Interacción/colisión**: donde los perfiles de ambas CMEs se solapan, la densidad no se suma aritméticamente (crearía masa) sino que se combina como promedio ponderado cuadrático (`densidad_solapada`), y la velocidad de la zona se calcula como **promedio ponderado por densidad** (`v_pond`), consistente con conservación de momento.
- **CME nueva (producto de fusión)**: cuando el solapamiento es persistente y denso, se "nuclea" una tercera estructura (`CMENueva`) con geometría fija heredada del clúster de solapamiento en ese instante, que luego se propaga con un modelo de arrastre tipo DBM relajando hacia la velocidad del viento solar. _Nota:_ El termino `CMENueva` hace referencia a una nueva nube de puntos en la región de solapamiento, no a la producción de una nueva CME desde el Sol; solo es nomenclatura.

## 2. Parámetros globales

| Variable | Significado físico | Notas de código |
|---|---|---|
| `R_SOL_KM` | Radio solar en km (unidad de distancia $R_\odot$) | Constante |
| `DIST_TIERRA_KM`, `DIST_TIERRA_RS` | 1 UA en km y en $R_\odot$ | Límite superior de observación |
| `DENSIDAD_FONDO` | Densidad del viento solar de fondo (protones/cm³) | Piso de densidad en todo el dominio |
| `T_HORAS` | Duración de la propagación a graficar (h) | `T_CALCULO_HORAS = max(T_HORAS, 95)` amplía el cálculo interno |
| `FACTOR_ESCALA` | Factor km → $R_\odot$ (=`R_SOL_KM`) | Usado en `radio_inter`, `CMENueva.radio` |
| `V_VIENTO_SOLAR` | Velocidad asintótica del viento solar (km/s) | Velocidad hacia la que relaja toda CME nueva |
| `semilla1`, `semilla2` | Semillas de la morfología de CME-1 y CME-2 | Reproducibilidad |
| `RETRASO_CME2` | Retardo de lanzamiento de CME-2 respecto a CME-1 (s) | |
| `FACTOR_COMPRESION` | Factor de compresión adicional en la zona de solapamiento | Multiplica la densidad combinada |
| `DMIN_OVERRIDE`, `DMAX_OVERRIDE` | Límites fijos de la escala de color $\log_{10}(\rho)$ | Si `None`, se calculan del máximo real |
| `N_PUNTOS_OBS` | Nº de puntos radiales donde se muestrean series temporales | Entre `R_OBS_MIN` y `R_OBS_MAX` |

## 3. Clases principales

### `_Morfo` (mixin de forma)
Genera el contorno deformado (frente/retaguardia) y el campo de velocidad unitario compartidos por `CME` y `CMENueva`.

| Método | Descripción |
|---|---|
| `forma(th, r_cme)` | Devuelve `(r_ext, r_int)`: radios de frente y retaguardia por ángulo `th`, combinando forma base + Fourier + filamentos + ventana angular de apertura. |
| `_vel_vec_base` | Campo vectorial de velocidad `(U, V)` en cartesianas |

### `CME` (CME-1 y CME-2, "originales")

| Atributo/parámetro | Significado físico |
|---|---|
| `tr`, `td` | Tiempos característicos de la fase de aceleración ($\tau_r$) y de frenado ($\tau_d$) |
| `ar`, `ad` | Amplitudes de aceleración residual/de frenado ($a_r$, $a_d$) |
| `v0`, `x0` | Velocidad y posición iniciales (km/s, km) |
| `R0` | Radio de referencia usado en la normalización de densidad (en $R_\odot$) |
| `t0` | Instante de lanzamiento (s); permite el retraso de CME-2 |
| `semilla`, `color` | Semilla morfológica y color de graficado |

| Método | Descripción |
|---|---|
| `aceleracion(s)` | $a(s)=\dfrac{a_r a_d}{a_d e^{-s/\tau_r}+a_r e^{s/\tau_d}}$, modelo de dos fases |
| `densidad(TH,R,r_cme,t,fi,r_ext_ov)` | Campo de densidad; `fi` es un factor de intensificación (usado en la zona de interacción) y `r_ext_ov` permite forzar el frente al radio de interacción `ri` |
| `vel_vec` | Campo de velocidad vectorial normalizado y escalado por la velocidad instantánea |

### `CMENueva` (producto de fusión/interacción)

| Atributo | Significado físico |
|---|---|
| `t_nac`, `r_nac`, `v_nac`, `d_nac` | Instante, radio, velocidad y densidad de nucleación (heredados del clúster de solapamiento) |
| `theta_centro`, `apertura_angular`, `grosor_frac` | Geometría angular fija, fijada al nacer |
| `theta_puntos`, `r_puntos` | Nube de puntos del clúster real de solapamiento |

## 4. Funciones clave de interacción

| Función | Rol físico |
|---|---|
| `v_pond(d1,d2,v1,v2)` | Velocidad de la zona combinada, promedio ponderado por densidad |
| `densidad_solapada(ci1,ci2,sol)` | Densidad en la región de solapamiento: $\dfrac{\rho_1^2+\rho_2^2}{\rho_1+\rho_2}\times$ |
| `geometria_sol(sol_mask,...)` | Extrae centro angular, apertura y grosor del clúster real de solapamiento, para nuclear una `CMENueva` con esa forma |
| `calc_campos` | Ensambla el campo de densidad total combinando zonas exclusivas de cada CME y la zona de solapamiento |
| `detectar_y_registrar` | Decide si el solapamiento es suficientemente denso y persistente como para nuclear una nueva CME, y la registra sin duplicar |

## 5. Variables de código (cinemática y grillas)

| Variable | Significado |
|---|---|
| `tiempos`, `tiempos_h` | Vector temporal de cálculo (s y h) |
| `pos1, vel1, acel1` / `pos2, vel2, acel2` | Fases de propagación de CME-1/CME-2, obtenidas integrando `aceleracion` |
| `t_centros`, `pos_centros` | Instante y posición donde se cruzan los **centros** (posiciones) de ambas CMEs |
| `t_extensiones`, `pos_extensiones` | Instante y posición donde el frente de CME-2 alcanza la retaguardia de CME-1 (primer contacto físico real) |
| `t_frames`, `idx_mostrar` | Instantes de cálculo de la propagación polar y subconjunto de 8 que se grafican como paneles |
| `cmes_nuevas` | Lista acumulativa de objetos `CMENueva` detectados durante la simulación |
| `DMIN`, `DMAX` | Rango de la escala de color en $\log_{10}(\rho)$ |
| `PUNTOS_OBS` | Lista de puntos $(r,\theta)$ tipo "satélite virtual" donde se registran series temporales de densidad/velocidad |
| `dens_ser`, `vel_ser` | Diccionarios `{punto: serie_temporal}` de densidad y velocidad en cada punto de observación |

## 6. Salidas generadas

- `cinematica_conjunta_s1_{semilla1}_s2_{semilla2}.pdf` — aceleración, velocidad y posición de ambas CMEs, con línea de tiempo de fases (inicio/aceleración/propagación) y marcadores de interacción de centros/extensiones.
- `cme_conjunta_polar_s1_{semilla1}_s2_{semilla2}.pdf` — 8 paneles polares de densidad (mapa de color) y velocidad (vectores) a lo largo de la propagación.
- `cme_conjunta_polar_panel8_s1_{semilla1}_s2_{semilla2}.pdf` — panel final aislado.
- `serie_temporal_multipunto_s1_{semilla1}_s2_{semilla2}.pdf` — evolución de densidad y velocidad radial en `N_PUNTOS_OBS` puntos entre 2 $R_\odot$ y 1 UA.

## 7. Notas de uso

- Los parámetros de `cme1`/`cme2` (`tr, td, ar, ad, v0, x0, R0`), `RETRASO_CME2` y `semilla1`/`semilla2` definen completamente el caso de estudio; cambiarlos reproduce distintos escenarios de interacción.



# Sobre cinemática_CMEs_conjuntas_interactivo.py

## ¿Para qué se usa?

Es una **herramienta interactiva de exploración**, no una simulación espacial de campos (a diferencia de `CME_combinadas.py`). Grafica en tiempo real $a(t)$, $v(t)$, $x(t)$ de dos cuerpos (dos CMEs) bajo el mismo modelo cinemático de dos fases, con controles en pantalla (`TextBox`, `Button`, `RadioButtons`) para variar sus parámetros y ver de inmediato el efecto.

Su uso principal en la tesis es:
- **Visualizar los perfiles cinemáticos** ($a$, $v$, $x$) de dos CMEs candidatas antes de fijarlas en la simulación completa, ajustando $a_r,\tau_r,a_d,\tau_d,v_0,x_0$ por ensayo-error hasta reproducir una interacción.
- **Buscar umbrales de interacción**: activando el "auto-ajuste de $T_{max}$", el propio programa extiende la ventana temporal hasta encontrar la primera intersección de posiciones $x_1(t)=x_2(t)$, marcándola con el instante, la distancia (km/AU) y el punto en la curva. Esto permite explorar rápidamente qué combinaciones de parámetros (o qué retraso `t_offset` entre ambas) producen o no un cruce de trayectorias, y en qué momento.
- Exportar/cargar configuraciones de parámetros en CSV, para guardar casos de interés.

No calcula densidad, morfología angular ni deformación de perfiles espaciales: es el laboratorio 1D (cinemática pura) que precede a la simulación 2D polar completa.

## 1. Modelo físico

Mismo modelo de aceleración de dos fases que `CME_combinadas.py`:

$$a(t) = \left[\frac{1}{a_r\, e^{t/\tau_r}} + \frac{1}{a_d\, e^{-t/\tau_d}}\right]^{-1}$$

## 2. Parámetros y constantes físicas

| Variable | Significado físico | Unidad |
|---|---|---|
| `ar`, `tr` | Amplitud y tiempo característico de la fase de **aceleración residual** ($a_r,\tau_r$) | km/s², s |
| `ad`, `td` | Amplitud y tiempo característico de la fase de **frenado/arrastre** ($a_d,\tau_d$) | km/s², s |
| `v0`, `x0` | Velocidad y posición iniciales de cada cuerpo | km/s, km |
| `t_offset` (`state['offset_s']`) | Retardo de lanzamiento aplicado a uno de los dos cuerpos (equivalente a `RETRASO_CME2`) | s |
| `T_max` | Tiempo total de integración/graficado (efectivo) | s |
| `S2H` | Conversión s→h | — |
| `KM2AU`, `KM2RSOL` | Conversión de posición a UA y a radios solares | — |
| `KM2M_S2` | Conversión de aceleración km/s² → m/s² (para graficar) | — |

## 3. Indicadores mostrados en pantalla

- $a(0)$ y $v(0)$ de cada cuerpo.
- Punto y valor del **95% de $v_{max}$** (criterio usado en la tesis para marcar el fin de la fase de aceleración).
- Marcador `★` con tiempo, posición en km y en AU de la **primera intersección de posiciones** entre los dos cuerpos (umbral de interacción cinemática).

## 4. Entradas/salidas

- **Entrada**: parámetros editados directamente en los `TextBox`, o cargados desde un CSV previamente exportado.
- **Salida**: `dos_cuerpos_{timestamp}.csv` — únicamente los parámetros (no las series), para que el caso sea reproducible al recargarlo o al trasladarlo a `CME_combinadas.py`.
