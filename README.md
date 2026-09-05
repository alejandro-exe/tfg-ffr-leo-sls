# Simulador SLS LEO-FFR

Simulador **a nivel de sistema** (SLS, *system-level simulator*) en MATLAB que evalúa la
viabilidad de la **reutilización fraccionaria de frecuencias (FFR)** —incluida una variante
**adaptativa**— como técnica de gestión de interferencia en constelaciones **LEO** bajo saturación
orbital y espectral.

Material complementario del Trabajo de Fin de Grado *«Viabilidad de FFR adaptativa en
constelaciones LEO bajo saturación»* (Ingeniería de Telecomunicación, UC3M).

---

## 1. Qué es y qué alcance tiene el modelo

El simulador construye una constelación Walker, la propaga, calcula la geometría respecto de una
**rejilla de usuarios fija en tierra**, resuelve el balance de enlace y la interferencia co-canal
—**intra**-satélite (entre haces del mismo satélite) e **inter**-satélite—, reparte el ancho de
banda según el esquema de reúso y emite un **veredicto de viabilidad**.

```
constelación → propagación (ECI) → ECEF → rejilla de usuarios → geometría (elev/az/rango/vis)
   → satélite servidor → balance C/N → interferencia inter-satélite → layout multihaz
   → contexto FFR → coloreado → política (τ, α) → asignación de banda → SINR → KPIs → VEREDICTO
```

La salida de decisión es siempre `K.viab.verdict` ∈ {`inviable`, `marginal`, `viable`}, definida
sobre el **percentil 5 de SINR de los usuarios de borde**, con el suelo de demodulación del MCS más
robusto de 5G NR (−6,7 dB) y un umbral de servicio de 0 dB con 95 % de cobertura.

**Alcance del modelo (lo que SÍ hace):**

| | |
|---|---|
| Enlace | **Descendente** (downlink), banda **Ku** (12 GHz) |
| Haces | **Earth-fixed**: el teselado de celdas está fijo en el suelo, no barre con el satélite |
| Propagación orbital | **Kepleriana circular** (sin J2 ni arrastre) |
| Canal | FSPL + patrón Bessel (3GPP TR 38.811) + atmósfera y lluvia **ITU-R P.618** |
| Terminal | Envolvente **ITU-R S.1428-1** (VSAT de 60 cm) |
| Capacidad | **Fórmula de Shannon** sobre el ancho realmente asignado |
| Esquemas | reúso-1, reúso-Δ (3 y 4), FFR estática y FFR adaptativa |
| Marco | Metodología de evaluación **3GPP TR 38.821** (NTN), calibrado contra ella |

**Lo que NO hace, y conviene saber antes de leer un resultado:** no implementa SGP4 (`propagate`
aborta con cualquier método que no sea `'kepler'`); no modela tablas MODCOD reales (la capacidad es
Shannon); el patrón de haz lleva un **suelo de lóbulos plano de −30 dB** en vez del perfil riguroso
ITU-R S.1528; los haces son circulares (sin ensanchamiento fuera de boresight ni celdas elípticas);
y la convención de potencia es **PSD constante**, no hay refuerzo de potencia en la sub-banda
interior.

---

## 2. Requisitos

### 2.1 Software

| Elemento | ¿Imprescindible? | Dónde se usa exactamente |
|---|---|---|
| **MATLAB R2022b** | **Sí** | No se ha probado en otras versiones |
| **Satellite Communications Toolbox** | **Sí** | `p618Config` + `p618PropagationLosses` en `src/radio/atm_loss_dB.m`, por donde pasa **todo** el pipeline. Además `satelliteScenario`/`satellite`/`groundStation`/`aer`/`states` en `validate_geometry_satscenario` |
| **Statistics and Machine Learning Toolbox** | **Sí** | `prctile` y `ecdf` en `src/ffr/compute_kpis.m`, `compute_beam_cir`, `ffr_policy` y casi todos los runners |
| Parallel Computing Toolbox | No | Sólo `src/infra/run_sweep_points.m` (`parfor`, `parpool`, `gcp`). **Sin ella `parfor` degrada a `for` sin error y sin cambiar ningún resultado**; la función lo detecta y lo deja escrito en `INFO.mode` |
| Aerospace, Communications, Phased Array, Antenna | No | **No se usan.** La geometría orbital, el patrón Bessel y la envolvente S.1428-1 están programados a mano |

### 2.2 Datos ITU-R (**no se distribuyen con este repositorio**)

`p618PropagationLosses` necesita los mapas digitales de la UIT: **`maps.mat`, `p836.mat`,
`p837.mat`, `p840.mat`** (~71 MB en total). **No son obra propia**: proceden del paquete de soporte
de MathWorks *ITU-R Digital Maps*, derivado de material de la UIT, y por eso no se incluyen aquí.

Para instalarlos, desde MATLAB:

```
Inicio → Complementos → Obtener complementos → buscar "ITU-R Digital Maps"
```

o, si ya tienes la Satellite Communications Toolbox instalada, ejecuta cualquier función que los
pida y MATLAB ofrecerá descargarlos. Los cuatro `.mat` deben quedar accesibles desde el path
(la raíz del proyecto sirve).

> ### ⚠ Aviso del modelo de respaldo
>
> **Sin esos ficheros el pipeline NO falla: cae a un modelo de respaldo** (`0.05/sin(el)` de gases
> y `3/sin(el)` de lluvia) y emite el aviso `atm_loss_dB:p618fallback`. **Los resultados obtenidos
> con el modelo de respaldo NO son los del TFG.** Si ves ese aviso, para y instala los mapas.

### 2.3 Hardware de referencia

Portátil AMD A10-5757M, 4 núcleos lógicos, **7,4 GB de RAM**. Todo lo publicado se reproduce ahí.
**El límite real es la RAM, no la CPU** (§6).

---

## 3. Puesta en marcha

```matlab
% 1. Abre MATLAB R2022b con la RAÍZ del proyecto como directorio de trabajo.
setup_paths                 % añade src/, experimentos/, tests/ y figuras/ al path

% 2. ¿Funciona la atmósfera? (~80 s en frío; la 2ª llamada es instantánea: hay caché)
cfg = config_default();
tic; atm_loss_dB(cfg, 45); toc
tic; atm_loss_dB(cfg, 45); toc

% 3. ¿Arrancan los tres perfiles? (~4 min: tres tablas P.618 distintas)
test_profiles_boot          % debe terminar en "TEST SUPERADO"

% 4. EL CANARIO: ¿da los números del TFG? (~90 s)
run_calibration_e0
```

El paso 4 debe imprimir exactamente:

```
CIR intra-satelite [dB]  : media=-0.53  mediana=-0.59  p5=-2.28  p50=-0.59  p95=1.28
RESULTADO  : CALIBRACION SUPERADA (dentro de 1.0 dB).
```

**`CIR p5 = −2,28 dB` es el canario de regresión del proyecto.** Ese resultado ha permanecido
bit-idéntico a lo largo de todo el desarrollo; si no sale así, algo se ha roto y no tiene sentido
seguir.

> ### ⚠ Ejecuta SIEMPRE desde la raíz
>
> Los runners guardan con **nombres desnudos** (`save('ffr_results.mat', ...)`) y sus figuras en
> **carpetas relativas** (`figs_e2/`, `figs_e5/`…), es decir **en el directorio de trabajo**.
> Ejecutar desde otro sitio dispersa los resultados y rompe los anclajes entre experimentos, que
> se buscan por nombre. `setup_paths` pone el código en el path; el directorio de trabajo lo pones
> tú, y debe ser la raíz.
>
> Consecuencia: **los `.mat` y las carpetas `figs_*/` aparecen en la raíz**, no dentro de
> `experimentos/`. Es deliberado y está contemplado en el `.gitignore`.

---

## 4. Los experimentos, uno por uno

Cada carpeta de `experimentos/` es un bloque autocontenido. Los tiempos son de la máquina de
referencia e **incluyen las tablas ITU-R P.618**, que es lo que domina (§7).

---

### `E0_calibracion` — ¿el motor reproduce un resultado publicado por el 3GPP?

**La pregunta.** Antes de defender ningún número propio hay que demostrar que el motor calcula bien.
E0 reproduce el escenario de calibración **Ka Set-1 LEO-600** de 3GPP TR 38.821 y compara la CIR
intra-satélite del haz central de un clúster de 19 haces contra los valores publicados en la
Tabla 6.1.1.2-1.

| | |
|---|---|
| **Qué barre** | Nada: es un punto único con el perfil `config_calib` (Ka, 20 GHz, h = 600 km, máscara 30°) |
| **Se lanza** | `run_calibration_e0` · luego `export_e0_figures` · validación externa `validate_geometry_satscenario` |
| **Produce** | `calibration_e0_results.mat` · 7 PNG en `figs_e0/` (los guarda `export_e0_figures`, no el runner) |
| **Tarda** | ~90 s el runner · ~20 s las figuras si van a continuación (comparten la tabla P.618 de Ka) |
| **Depende de** | Nada |

**Resultado de referencia:** CIR p5 = **−2,28 dB**, dentro de 1 dB de los valores 3GPP. Es el
canario del proyecto (§3). `validate_geometry_satscenario` contrasta además la órbita propia contra
`satelliteScenario` de la toolbox: es una validación **independiente**, toolbox frente a aportación
propia.

---

### `E1_capas_base` — ¿funciona cada capa por separado?

**La pregunta.** Validación por capas antes de componerlas: geometría, balance de enlace, Doppler e
interferencia inter-satélite, cada una con sus propias comprobaciones de cordura.

| | |
|---|---|
| **Qué barre** | `run_interference_demo` barre densidad × apuntamiento (`nadir` / `boresight`); los otros tres son puntos únicos |
| **Se lanza** | `run_geometry_demo` → `run_link_budget_demo`, `run_doppler_demo`, `run_interference_demo` |
| **Produce** | `geometry_results.mat`, `link_budget_results.mat`, `doppler_results.mat`, `interference_results.mat` · 6+4+2+3 PNG en `figs_geom/`, `figs_link/`, `figs_doppler/`, `figs_interf/` |
| **Tarda** | ~40 s + ~2 min + ~1 min + ~2 min |
| **Depende de** | `run_link_budget_demo` y `run_doppler_demo` cargan `geometry_results.mat`; si falta, lo regeneran solos |

> **`run_geometry_demo` lleva ACTIVO un override de validación** (`T = 66`, `P = 6`): produce a
> propósito una constelación reducida para que la capa geométrica se pueda inspeccionar. **Los
> experimentos no lo heredan**, cada uno fija su densidad. Para verlo con la constelación completa
> 1584/72, comenta ese bloque.

**Hallazgo del bloque:** la interferencia co-canal **inter**-satélite en reúso-1 es **despreciable**
(penalización ≤ 0,16 dB incluso con 1584 satélites), por **doble discriminación**: el patrón del
satélite interferente, que hacia nuestro usuario está 6–38 anchos de haz fuera de eje, y la
directividad del terminal VSAT. La interferencia que importa es la **intra**-satélite, y eso es lo
que motiva todo el resto del trabajo.

---

### `E2_ffr_estatica` — ¿protege la FFR el borde de celda? (hipótesis H1)

**La pregunta.** El experimento de cabecera. Sobre la misma población de usuarios, los mismos
instantes y la misma clasificación centro/borde, compara los cinco esquemas de referencia y decide
si cada uno es viable.

| | |
|---|---|
| **Qué barre** | Cinco esquemas (reúso-1, reúso-Δ 3 y 4, FFR Δ=3 y Δ=4 con α = 0,4) + un barrido de α = 0,1…0,9 |
| **Se lanza** | `run_ffr_demo` · comprobación auxiliar opcional `check_buser_ratio` |
| **Produce** | `ffr_results.mat` (~37 MB) · 6 PNG en `figs_e2/` |
| **Tarda** | ~3,7 min (133 s de experimento + 88 s de tabla P.618) |
| **Depende de** | Nada. **Pero es ancla bloqueante de `densidad_orbital`** |

Escenario: Ku, Walker **53:1584/72/1** (Starlink Shell-1 completo), rejilla local de 40 km con paso
de 3 km (553 usuarios), 2 h con `dt = 60 s`, clúster de 91 haces, τ por cuantil al 50 %,
planificación `'share'`, `timeBlock = 10`.

**Resultado:** el `R_p5` de borde pasa de **5,6 Mbps** en reúso-1 a **9,9 Mbps** con FFR(Δ=3)
(×1,77), y `SINR_edge_p5` de −4,97 a **+3,80 dB**. Reúso-1 **nunca** alcanza `viable`; la FFR sí, con
cobertura 1,00. **H1 confirmada.**

---

### `E5_adaptativa` — ¿recupera eficiencia la FFR adaptativa? (hipótesis H2)

**La pregunta.** Si α se adapta a la geometría (más banda al borde cuando la elevación es mala),
¿se gana algo frente al mejor α fijo? El runner no compara sólo contra la mejor estática: construye
la **frontera de compromiso** completa barriendo α y comprueba si el punto adaptativo cae por encima.

| | |
|---|---|
| **Qué barre** | reúso-1 / FFR estática con α\* / FFR adaptativa / FFR adaptativa **invertida**, más α = 0,05…0,95 en pasos de 0,05 |
| **Se lanza** | `run_e5_adaptive` |
| **Produce** | `e5_adaptive_results.mat` (~28 MB) · 5 PNG en `figs_e5/` |
| **Tarda** | ~3,5 min (132 s + 79 s de tabla) |
| **Depende de** | Nada; usa **el mismo escenario que E2** a propósito, para que las cifras sean directamente comparables |

**Resultado: H2 NO se sostiene.** La adaptativa reproduce el compromiso de la mejor estática, no lo
mejora (−0,063 Mbps frente a un umbral de materialidad de 0,413). La causa está identificada y es
estructural: la sub-banda interior sólo renta si **Δ·SE_centro/SE_borde > 1**, y ese cociente vale
0,78–0,86 en todos los tramos de elevación. El runner incluye tres validaciones de honestidad
—caso degenerado con error exactamente 0, contraste con el sentido invertido, y barrido de la
familia de anclas— para descartar que el resultado negativo sea un artefacto.

---

### `E3_densidad_haces` — ¿satura el sistema al meter más haces? (eje espectral)

**La pregunta.** Densificar el plan de celdas mete más interferentes co-canal sobre la misma zona.
¿Hay un punto en el que el borde deja de ser viable?

| | |
|---|---|
| **Qué barre** | `beamwidth3dB_deg` = 3,0 / 2,2 / 1,5 / 1,0 / 0,7 / 0,5° → **1,85 → 66,8 haces/1000 km²** (×36), en **dos variantes**: (A) `nRings` fijo y (B) guarda preservada |
| **Se lanza** | `run_e3_beamdensity` |
| **Produce** | `e3_results.mat` (~5,7 MB) · 4 PNG en `figs_e3/` |
| **Tarda** | ~6 min |
| **Depende de** | Nada |

**Resultado:** ×36 de densidad cuesta sólo **0,84 dB** (reúso-1) y **1,89 dB** (FFR Δ=3) de
`SINR_edge_p5`, mientras el caudal agregado crece **×12,8**. Densificar no satura: multiplica la
capacidad. La variante A sale **sesgada en 3 de 6 puntos** (la rejilla de usuarios sobresale del
clúster de haces) y su aparente "saturación" es un artefacto de truncamiento; **la variante B es la
defendible**. El runner marca los puntos sesgados con dos criterios independientes.

---

### `E3b_mascara_elevacion` — ¿por qué se hunde el reúso a elevación baja?

**La pregunta.** Con haces Earth-fixed, al bajar la elevación el clúster se ve **escorzado** desde el
satélite: los haces se solapan más y la interferencia intra empeora. Este eje aísla ese mecanismo.

| | |
|---|---|
| **Qué barre** | `cfg.geom.minElev` = 45 / 40 / 35 / 30 / 25 / 20 / 15 / 10°, en **dos densidades** (T = 1584 y T = 66), con 6 esquemas |
| **Se lanza** | `run_e3b_minelev` |
| **Produce** | `minelev_results.mat` (~8,2 MB) · 4 PNG en `figs_e3b/` |
| **Tarda** | **~18 min** — es el más caro, y la razón es la caché (§7): `minElev` está en la clave, luego paga **8 tablas P.618** |
| **Depende de** | Nada. **Es ancla de `EB_ancho_banda`** |

**Hallazgo de diseño:** a T = 1584 **la máscara no muerde** (el mejor de 1584 satélites nunca baja de
52°) y el barrido sale exactamente degenerado —recorrido 0,000 dB en las 8 máscaras—. Por eso se
ejecuta también a T = 66, donde el eje sí existe. La compresión angular queda confirmada como
mecanismo: correlación **+0,9997** entre separación efectiva y `SINR_edge_p5`.

---

### `E4_interconstelacion` — ¿cuánto molesta otro operador en la misma banda?

**La pregunta.** La coordinación intra-operador protege dentro de una constelación, pero **entre
operadores no existe**. ¿Qué pasa cuando un satélite ajeno cruza cerca de la línea de vista?

| | |
|---|---|
| **Qué barre** | 1 operador (Starlink solo) vs 2 (Starlink + **OneWeb**), × apuntamiento nominal (`nadir`) y cota superior (`boresight`) |
| **Se lanza** | `run_e4_interconstellation` |
| **Produce** | `interconstellation_results.mat` (~1,5 MB) · 4 PNG en `figs_e4/` |
| **Tarda** | ~13 min |
| **Depende de** | Nada |

Trabaja en **dos escalas de tiempo**: un nivel **fino** (`dt = 1 s`, 7201 instantes) para la
estadística de **eventos in-line** —un LEO barre ~1 °/s, luego un evento dura segundos y el `dt` del
pipeline lo perdería— y un nivel KPI con la rejilla completa para los veredictos.

**Sólo entra OneWeb**, porque el criterio es **compartir la banda de usuario Ku**. Kuiper opera su
enlace de usuario en **Ka** y no puede interferir co-canal en Ku: queda fuera, y es una corrección
física, no un recorte de alcance.

**Resultado:** severidad alta, frecuencia baja. Bajo `boresight` la SINR instantánea cae 9,4 dB
durante un evento, pero `P(ψ < 2°) = 0,222 %` está por debajo del 1 %, luego los eventos no llegan
a poblar el percentil 1 y la cola sólo se mueve 0,528 dB. **Ningún veredicto cambia.**

---

### `E6_oneweb` — ¿valen las conclusiones en otra arquitectura?

**La pregunta.** Todas las conclusiones se obtuvieron con un sistema. Repetirlas sobre una
arquitectura radicalmente distinta (1200 km en vez de 550, órbita polar en vez de delta, máscara de
55° en vez de 25°, celdas de 285 km en vez de 17) las convierte de propiedades *del sistema
simulado* en propiedades *del reúso multihaz en LEO* — o acota su dominio.

| | |
|---|---|
| **Qué barre** | Perfil OneWeb completo: máscara de elevación (6 valores), densificación de haces con la huella fija (`nRings` 2…5) y **contraste de topología** hexagonal vs retícula 4×4 |
| **Se lanza** | `run_oneweb_demo` |
| **Produce** | `oneweb_results.mat` · 4 PNG en `figs_oneweb/` |
| **Tarda** | **~20 min** (6 tablas P.618 por el barrido de máscara) |
| **Depende de** | Nada |

**Resultado: H1 y la refutación de H2 SE REPRODUCEN**, con más margen incluso (+9,12 dB de ganancia
de la FFR frente a +8,82 dB en Starlink). Y **tres resultados quedan ACOTADOS**: la invariancia de
escala de la retícula sólo vale para haz estrecho (HPBW ≲ 10°), el eje de elevación sale degenerado
cuando la máscara no muerde, y **la ventaja de Δ=3 es un constructo hexagonal** —sobre celdas
cuadradas pierde 3,2 dB y conviene Δ=4—.

---

### `EB_ancho_banda` — ¿depende la ventaja de la FFR del ancho de banda del sistema?

**La pregunta.** Con más ancho de banda caben más usuarios: ¿deja de ser crítica la FFR?

| | |
|---|---|
| **Qué barre** | `cfg.radio.B_MHz` = 62,5 / 125 / 250 / 500 / 1000 (**×16**) con el **número de usuarios constante** (M = 317, verificado con `error()`) |
| **Se lanza** | `run_eB_bandwidth` |
| **Produce** | `eB_bandwidth.mat` (~1,5 MB) · 4 PNG en `figs_eB/` |
| **Tarda** | ~4 min (una sola tabla P.618: `B_MHz` **no** está en la clave de la caché) |
| **Depende de** | `minelev_results.mat` para la regresión del punto B = 250. Sin él se ejecuta igual, pero no declara la regresión superada |

**Resultado: NULO, y exacto.** La ganancia de la FFR sobre reúso-1 es idéntica a 62,5 y a 1000 MHz
(desviación **1,07e-14 dB**). B entra por igual en C, I y N, luego se cancela en la SINR y sólo
escala la tasa. El runner incluye además la métrica derivada de **capacidad en usuarios admitidos**
y un diagnóstico de **cuantización en PRB**, que es la única vía por la que un B pequeño penalizaría
de verdad a la partición.

---

### `densidad_orbital` — ¿satura el sistema al meter más satélites?

**La pregunta.** El eje de saturación orbital propiamente dicho: ¿hay una densidad de satélites a
partir de la cual el borde deja de ser viable?

| | |
|---|---|
| **Qué barre** | **T = 66 … 4000** con **P escalado** (capa Walker realista: `T/P ∈ [22, 25]`) × 6 esquemas |
| **Se lanza** | `run_density_sweep` |
| **Produce** | `e3_density_sweep.mat` (~3,4 MB) · sin figuras |
| **Tarda** | ~11 min |
| **Depende de** | **`ffr_results.mat` (E2), y de forma BLOQUEANTE**: su punto T = 1584 / P = 72 debe reproducirlo con `max|dif| = 0` o **no guarda** |

**Resultado contraintuitivo:** subir la densidad orbital **no satura el sistema, lo mejora**
(reúso-1 pasa de −11,61 a −4,51 dB), porque más satélites significan servidor a mayor elevación y
por tanto menos compresión angular del clúster. El margen de viabilidad no tiene punto de cruce: lo
que hay es que **reúso-1 nunca alcanza `viable` en todo el rango** y la FFR cruza en T ≈ 600.

Requiere `timeBlock = 4` y ejecución secuencial: su peor punto (T = 4000) pediría **18,67 GB** sin
trocear (§6).

---

### `convergencia` — ¿están bien elegidos los parámetros numéricos?

**La pregunta.** Cuántos usuarios, cuántos instantes y qué tamaño de clúster hacen falta para que
el resultado no dependa de la discretización. Es lo que legitima el coste de todos los demás
experimentos.

| script | qué barre | produce | tarda |
|---|---|---|---|
| `run_convergence_studyA` | `step_km` = 8…2 · `Nt` = 31…241 · 5 realizaciones · T = 66…4000 · informe de memoria | `convergence_studyA.mat` · 3 PNG en `figs_faseA/` | ~10–15 min |
| `run_convergence_nt_step4` | repite el eje `Nt` con la rejilla **de trabajo** (`step_km = 4`) y 3 esquemas | `convergence_nt_step4.mat` · 1 PNG | ~5 min |
| `run_convergence_step_3sch` | repite el eje **de rejilla** con 3 esquemas | `convergence_step_3sch.mat` · 1 PNG | ~4 min |
| `run_convergence_nrings` | `nRings` = 2…7 (19→169 haces), con la interferencia descompuesta **por anillo** | `convergence_nrings.mat` · 3 PNG en `figs_nrings/` | ~4 min |

**Dependencia bloqueante:** `run_convergence_nt_step4` y `run_convergence_step_3sch` **cargan
`convergence_studyA.mat` y abortan** si no reproducen sus tablas `T2` y `T1` con `max|dif| = 0`. Es
lo que impide que la copia local de la función de evaluación diverja en silencio de la original.

**Parámetros adoptados:** `step_km = 4` (más fino que el convergido, 6, por conservadurismo),
`Nt = 61`, `nRings = 5`. El estudio de `nRings` además **descarta el tamaño de clúster como eje de
saturación**: su cola la fija el suelo de lóbulos plano, es decir un artefacto del modelo de patrón,
no interferencia física.

---

### `h2_rerun` — ¿es H2 refutable por causa, o sólo falló una regla?

**La pregunta.** E5 refuta H2 con *una* regla concreta (α interpolado con la elevación). ¿Y si el
problema era esa regla y no la hipótesis? Tres piezas encadenadas cierran H2 de fuera hacia dentro.

| script | qué hace | produce | tarda |
|---|---|---|---|
| `export_h2_criterion` | Consolida el criterio Δ·SE_c/SE_b de E5, del eje de densidad y de OneWeb en una tabla plana, con verificación bloqueante de 17 cifras publicadas. **Post-proceso puro** | `h2_rerun/results/h2_criterion.mat` + log | segundos |
| `run_h2_anchors` | ¿El resultado dependía de unas anclas mal puestas? Barrido 2D de anclas de elevación, familia de α y barrido de τ | `h2_rerun/results/h2_anchors.mat` + log | ~4 min |
| `run_h2_oracle` | **Cota oráculo**: el mejor α(t) posible sin regla alguna, sobre un modelo afín validado contra el motor | `h2_rerun/results/h2_oracle_v2.mat` + log | **~71 min** |
| `figs/regen_cap7` | Regenera las figuras del capítulo desde los `.mat`, con referencias escritas en el script | 10 PNG en `h2_rerun/figs/` | segundos |

**Dependencias:** `export_h2_criterion` lee `e5_adaptive_results.mat`, `e3_density_sweep.mat` y
`oneweb_results.mat`. `run_h2_oracle` **exige `h2_anchors.mat` al día**.

**Resultado:** liberar la regla la **empeora** (las anclas nominales ya eran las mejores de la
rejilla), y el α(t) óptimo **no sigue la elevación**: donde correlaciona, lo hace con el **signo
contrario** al que postula la regla. La refutación de H2 se sostiene sobre la **forma** de la regla,
no sobre la inexistencia de margen — que es un enunciado más fuerte y más honesto.

> **`h2_rerun/` vive en la RAÍZ, no bajo `experimentos/`, y es a propósito:** sus scripts calculan
> `ROOTDIR = fileparts(H2DIR)` a partir de `mfilename('fullpath')` y cargan los `.mat` desde ahí.
> Moverlo un nivel más abajo haría que buscasen los resultados en `experimentos/` y abortarían.

---

### `tests/` — infraestructura

| test | qué demuestra | tarda |
|---|---|---|
| `test_profiles_boot` | Los tres perfiles construyen `cfg`, pasan `assert_cfg_coherent` y completan un punto del pipeline | ~4 min |
| `test_atm_cache` | La caché P.618 no altera ningún resultado (bit a bit) y sus 5 campos de clave son los correctos | ~2 min |
| `test_parfor_invariance` | `parfor` no cambia el modelo: `nWorkers = 0` vs automático, campo a campo con CDF incluidas | ~5 min |
| `test_timeblock_invariance` | El troceado temporal es **exacto**: `timeBlock` = Inf / 20 / 10 con `Nt = 41` (el bloque no divide la ventana) | ~3 min |
| `test_ctx_blocked` | `build_ctx_blocked` es exacto por tres vías | ~2 min |
| `test_mem_estimator` | El estimador de memoria reproduce los picos medidos. **No simula** | segundos |

> **`test_ctx_blocked` y el `.mat` que no se distribuye.** Su tercera comprobación compara contra
> `ffr_results_T66.mat`, un resultado archivado que **no forma parte de este repositorio** (aquí sólo
> hay código). El test lo detecta con `exist()`, imprime `ffr_results_T66.mat no encontrado: se
> omite.` y **se pronuncia sobre las otras dos**: no falla. Para recuperar la tercera hay que
> regenerar ese `.mat` con el perfil de su sección 1, que es exactamente el que lo produjo.

---

### `figuras/scripts` — regeneración de las figuras de la memoria

Leen los `.mat` publicados y regeneran las figuras del documento **sin volver a simular**
(`regen_cap5`, `regen_cap6`, `regen_cap8`, los `make_fig_4_*`, más `verificar_fidelidad` y
`auditar_regenerabilidad`). Segundos cada uno. Requieren que los `.mat` de los experimentos
correspondientes existan en la raíz.

> **Escriben en `REDACCION/`, no en `figuras/`.** La ruta de salida está escrita en el código
> (`fullfile(PROJROOT,'REDACCION','figuras')`), así que al ejecutarlos crean una carpeta `REDACCION/`
> en la raíz. Es una inconsistencia conocida entre el nombre de la carpeta de scripts y el de la
> carpeta de salida; **no se ha corregido porque exigiría editar código**, y las salidas están
> excluidas del repositorio de todos modos.

---

## 5. Orden de ejecución desde cero

Las **anclas bloqueantes** están marcadas con ⛔: el script aborta o no guarda si su entrada no está
al día.

| # | Bloque | Scripts | Tiempo |
|---|---|---|---|
| 1 | Verificación | `test_profiles_boot`, **`run_calibration_e0`** (canario) | ~6 min |
| 2 | Capas base | `run_geometry_demo` → `run_link_budget_demo`, `run_doppler_demo`, `run_interference_demo` | ~6 min |
| 3 | Calibración | `export_e0_figures`, `validate_geometry_satscenario` | ~1 min |
| 4 | Infraestructura | `test_atm_cache`, `test_mem_estimator`, `test_timeblock_invariance`, `test_parfor_invariance`, `test_ctx_blocked` | ~12 min |
| 5 | Convergencia | `run_convergence_studyA` → ⛔ `run_convergence_nt_step4`, ⛔ `run_convergence_step_3sch`, `run_convergence_nrings` | ~26 min |
| 6 | Cabecera | **`run_ffr_demo`** (E2), `run_e5_adaptive` (E5) | ~7 min |
| 7 | Ejes de saturación | `run_e3_beamdensity`, `run_e3b_minelev` → `run_eB_bandwidth`, `run_e4_interconstellation`, `run_oneweb_demo` | ~61 min |
| 8 | Densidad orbital | ⛔ `run_density_sweep` (exige `ffr_results.mat` del paso 6) | ~11 min |
| 9 | Reejecución de H2 | `export_h2_criterion`, `run_h2_anchors` → ⛔ `run_h2_oracle`, `figs/regen_cap7` | ~75 min |
| 10 | Figuras | `figuras/scripts/regen_cap5`, `regen_cap6`, `regen_cap8`, `make_fig_4_*`, `verificar_fidelidad` | ~2 min |

**Corpus principal (pasos 1–8): ≈ 2 h 10 min. Todo: ≈ 3 h 30 min.**

Resumen de anclas:

| Script | Exige | Si falta |
|---|---|---|
| `run_convergence_nt_step4` · `run_convergence_step_3sch` | `convergence_studyA.mat` | **Abortan** |
| `run_density_sweep` | `ffr_results.mat` | **No guarda**: aparta el barrido y aborta |
| `run_h2_oracle` | `h2_rerun/results/h2_anchors.mat` | **Aborta** |
| `run_eB_bandwidth` | `minelev_results.mat` | Corre, sin declarar la regresión |
| `test_ctx_blocked` | `ffr_results_T66.mat` (no distribuido) | Omite 1 de 3 comprobaciones |
| `run_link_budget_demo` · `run_doppler_demo` | `geometry_results.mat` | Lo regeneran solos |

---

## 6. ⚠ Memoria y troceado temporal

**El pico de memoria del simulador no está en los KPIs: está en el tensor de geometría
`[M × N × Nt]`** (usuarios × satélites × instantes), dentro de `compute_interference`. El modelo
medido es **≈ 74 bytes por elemento** `M·N·Nt`, más ~40 B por elemento `M·nBeams·Nt`, y es **cota
superior** de lo observado (+30 % sobre el pico real del proceso), que es lo que se quiere para
decidir cuántos workers caben.

En la configuración de cabecera (M = 553, N = 1584, Nt = 121) eso son **7,53 GB de una sola pieza**,
**más que la RAM total de la máquina de referencia (7,45 GB)**. Sin trocear, esos experimentos
**no se pueden ejecutar**.

**La solución es el troceado temporal, `cfg.compute.timeBlock`:** procesa la ventana en bloques de
instantes, acumula las matrices de muestra `[M × Nt]` y llama a `compute_kpis` **una sola vez al
final**. Ésa es exactamente la razón de que sea **EXACTO y no una aproximación**: el pipeline no
acopla instantes distintos —geometría, servidor, balance de enlace, interferencia y ganancias de haz
son por instante—, y **no se combinan KPIs parciales**, que es donde sería fácil equivocarse
(*un percentil 5 global no es la media de los percentiles 5 de cada bloque*). Verificado con
**`max|dif| = 0`** por tres vías independientes: `test_timeblock_invariance`, `test_ctx_blocked` y
`test_parfor_invariance`.

Con `timeBlock = 10` ese mismo pico baja a **0,62 GB**.

> **Valores fijados en los propios runners, y por qué:**
> - `run_ffr_demo` (E2) y `run_e5_adaptive` (E5) traen `cfg.compute.timeBlock = 10`. **Si lo pones
>   a `[]` se quedarán sin memoria.** Bajarlo a 5 da más margen y **no cambia ningún resultado**.
> - `run_density_sweep` necesita `timeBlock = 4` y `nWorkers = 0`: su peor punto (T = 4000) pediría
>   **18,67 GB** sin trocear, y con `timeBlock = 10` aún 1,54 GB por worker.

Antes de lanzar un barrido nuevo, estima:

```matlab
MEM = estimate_sweep_memory(cfg, pts, nWorkers);   % pico/worker, RAM libre, timeBlock necesario
```

**Bajo `parfor` cada worker mantiene su PROPIA geometría**: 4 workers = 4× la memoria, no se
comparte nada. En una máquina con poca RAM, `nWorkers = 0` + `timeBlock` es mejor que quedarse sin
memoria — y da exactamente el mismo resultado, por el mismo camino de código.

---

## 7. La caché ITU-R P.618

`atm_loss_dB` construye una tabla de pérdidas frente a elevación (`minElev:2:89` grados) llamando a
`p618PropagationLosses` y la interpola. Llamar a la toolbox por cada terna (usuario, satélite,
instante) sería inviable; para una ubicación fija la función es determinista, así que la tabla se
construye **una vez** y se guarda en una variable `persistent`.

Construirla cuesta **~80 s**, y **entre el 30 % y el 50 % del tiempo total del proyecto se va en
eso, no en simular.**

**La clave de la caché tiene exactamente cinco campos:**

```
freq_GHz | ground.point(1) | ground.point(2) | geom.minElev | p618_availability
```

Consecuencias operativas, que explican los tiempos de §4 y §5:

- Cambiar cualquiera de los cinco **obliga a reconstruirla**. Por eso `run_e3b_minelev`, que barre
  `minElev`, paga **8 tablas** y tarda ~18 min; y `run_oneweb_demo`, que barre la máscara en 6
  valores, paga 6.
- **`T`, `P`, `B_MHz`, `cfg.ffr.*`, `cfg.beams.*` y `ground.step_km` NO están en la clave**: un
  barrido entero sobre cualquiera de ellos comparte **una sola** tabla. Es la razón de que
  `run_eB_bandwidth` (×16 en ancho de banda) pague sólo una.
- Los runners del perfil Ku comparten la misma clave, así que **ejecutados en la misma sesión de
  MATLAB sólo pagan la primera**. Agrupa por perfil si quieres ahorrar tablas enteras.
- Bajo `parfor` cada worker es un proceso con su propio espacio `persistent`: construye la tabla
  **una vez por worker**. `run_sweep_points` la precalienta fuera del cronómetro para que el
  speedup medido no sea un artefacto.
- Para vaciarla a mano: `clear atm_loss_dB`.

---

## 8. Cómo citar

<!-- TODO (autor): completar la referencia del TFG cuando esté depositado. -->

```
TODO: Autor, «Viabilidad de FFR adaptativa en constelaciones LEO bajo saturación»,
      Trabajo de Fin de Grado, Universidad Carlos III de Madrid, 2026.
      DOI: TODO
```

**Licencia:** ver `LICENSE`. **TODO: pendiente de elegir por el autor.**

**Datos de terceros:** los mapas digitales ITU-R que necesita el simulador **no se distribuyen aquí**
y no son obra propia (paquete de soporte de MathWorks, derivado de material de la UIT). Ver §2.2.

**Normas y fuentes utilizadas** (citadas en el código, en el punto donde se aplican):
3GPP TR 38.811 (patrón de haz), TR 38.821 (metodología NTN, escenario de calibración, cadena de
ruido, separación inter-haz), TS 38.214 (umbrales MCS), TS 38.211 (numerología PRB),
ITU-R P.618 (atmósfera y lluvia), ITU-R S.1428-1 (patrón del terminal), ITU-R M.2135-1 y M.2410-0
(definición de *cell edge* y objetivos de eficiencia espectral), y las solicitudes ante la FCC de
SpaceX y OneWeb para los parámetros de sistema.
