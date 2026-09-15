Este script simula la evolución individual de una eyección de masa coronal (CME) desde su iniciación hasta su propagación en el medio interplanetario. Combina un modelo cinemático de dos fases (aceleración y decaimiento) para describir la posición, velocidad y aceleración de la CME en el tiempo, con una representación morfológica bidimensional en coordenadas polares que reconstruye su forma, densidad y campo de velocidad conforme se expande y avanza. La morfología incorpora asimetría angular, ruido de Fourier y filamentos generados aleatoriamente a partir de una semilla fija, buscando una apariencia más realista que una simple envolvente geométrica idealizada. El resultado son dos visualizaciones: los perfiles cinemáticos completos y un panel de instantáneas que muestra cómo se deforma y diluye la estructura de la CME a medida que se aleja del Sol.

## Parámetros de simulación

### Constantes físicas y de escala

| Parámetro | Descripción |
|---|---|
| `R_SOL_KM` | Radio solar en km, usado para convertir posiciones a unidades de $R_\odot$ |
| `FACTOR_ESCALA` | Factor de conversión de km a $R_\odot$ para las mallas polares (= `R_SOL_KM`) |
| `DENSIDAD_FONDO` | Densidad del viento solar ambiente, en protones/cm³ |
| `T_HORAS` | Duración total de la simulación, en horas |

### Cinemática (modelo de dos fases, Gallagher et al. 2003)

| Parámetro | Descripción |
|---|---|
| `tr2` | Tiempo característico de la fase de ascenso (rise) |
| `ar2` | Aceleración característica de la fase de ascenso |
| `td2` | Tiempo característico de la fase de decaimiento (decay) |
| `ad2` | Aceleración característica de la fase de decaimiento |
| `v02` | Velocidad inicial de la CME (km/s) |
| `x02` | Posición inicial de la CME (km) |
| `R_CME_INIC` | Radio inicial de la CME en $R_\odot$, usado como referencia para el factor de expansión de la densidad |

### Morfología (forma angular, coordenadas polares)

Estos parámetros controlan la asimetría, el ruido angular y los filamentos de la envolvente de la CME. Se generan aleatoriamente a partir de `SEMILLA` mediante `np.random.default_rng`, por lo que **la misma semilla siempre produce la misma morfología** — útil para reproducir figuras exactas o para variar la forma entre distintas corridas de forma controlada.

| Parámetro | Descripción |
|---|---|
| `SEMILLA` | Semilla del generador aleatorio; fija toda la morfología (ruido, asimetría, filamentos) |
| `N_MODOS` | Número de modos de Fourier usados para el ruido angular en el borde exterior y en el grosor de la envolvente |
| `amp_ext`, `fase_ext` | Amplitudes y fases (una por modo) del ruido de Fourier aplicado al radio exterior de la CME |
| `amp_gros`, `fase_gros` | Amplitudes y fases del ruido de Fourier aplicado al grosor de la envolvente |
| `asimetria` | Factor de distorsión angular entre el lado positivo y negativo de $\theta$ (1.0 = simétrico) |
| `N_FIL` | Número de "filamentos" o protuberancias locales superpuestas al borde exterior |
| `ang_fil` | Ángulo central de cada filamento (radianes, dentro de $[-\pi/2, \pi/2]$) |
| `amp_fil` | Amplitud (altura) de cada filamento, como fracción del radio de la CME |
| `ancho_fil` | Ancho angular (dispersión gaussiana) de cada filamento |

### Escala de color y densidad

| Parámetro | Descripción |
|---|---|
| `DMIN_OVERRIDE` | Límite inferior manual de la escala de color logarítmica ($\log_{10}\rho$). Si es `None`, se calcula automáticamente a partir de `DENSIDAD_FONDO` |
| `DMAX_OVERRIDE` | Límite superior manual de la misma escala. Si es `None` o `0`, se calcula a partir del máximo global de densidad detectado en la simulación |

### Visualización polar

| Parámetro | Descripción |
|---|---|
| `ESTADOS_TIEMPO` | Número de instantáneas (frames) mostradas en el panel de propagación polar |
| `ETIQUETAS` | Letras usadas para rotular cada subpanel (a, b, c, ...) |
| `COLOR_CINE` | Color de las curvas en las gráficas cinemáticas (posición/velocidad/aceleración) |

### Notas sobre la malla espacial

- La malla polar global usa 800 puntos en $\theta \in [-\pi, \pi]$ y 400 puntos en $r \in [0, r_{\max}]$, donde $r_{\max}$ se calcula como 1.4 veces la posición final de la CME.
- Los primeros 4 de los `ESTADOS_TIEMPO` frames usan una malla local adaptativa (radio máximo = 2.7 × posición de la CME en ese instante) para dar mayor detalle cerca del Sol; los frames restantes usan la malla global fija. Este umbral (`idx < 4`) está codificado directamente en el bucle y no es un parámetro configurable en el bloque de constantes.
- La malla de vectores de velocidad es más gruesa (30 × 8 puntos) que la malla de densidad, por razones de legibilidad del gráfico de flechas (`quiver`).
