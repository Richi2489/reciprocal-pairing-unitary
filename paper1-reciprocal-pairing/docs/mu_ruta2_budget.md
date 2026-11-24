# Ruta 2 — el presupuesto fuera-de-diagonal y la fracción de emparejamiento forzado ρ(d)

_open_math_lab (banco Hilbert–Pólya). Fundamento analítico para la cota inferior de μ(d)._
_Trabajo de sprint autónomo 2026-07-03, en paralelo a la malla de escalamiento (Ruta 1)._
_Distinción de rigor: la cota superior es un TEOREMA; la casi-tightness y el piso de ρ son
OBSERVACIONES numéricas (n=3 exacto vía `barrer`, n=4 preliminar vía SLSQP)._

## 1. Identidad del presupuesto (exacta)

Para U∈U(n), cada fila tiene norma 1: Σ_j |U_ij|² = 1. Sumando sobre i y restando la diagonal,

    Σ_{i≠j} |U_ij|²  =  n − Σ_i |U_ii|²  =  n − Σd_i     [IGUALDAD, no cota].

Llamamos **presupuesto fuera-de-diagonal** a P(d) := n − Σd. Es una cantidad fija por la
diagonal, independiente de las fases y del resto de la estructura de U.

## 2. Cota superior de μ (teorema)

M(U) = Σ_{i<j} |U_ij·U_ji| = Σ_{i<j} √(|U_ij|²·|U_ji|²)  ≤  Σ_{i<j} ½(|U_ij|²+|U_ji|²)
     = ½ Σ_{i≠j} |U_ij|²  =  P(d)/2,      (AM-GM entrada a entrada).

Por tanto, para toda diagonal realizable d:

    **μ(d) = mín_U M(U)  ≤  (n − Σd)/2.**                                      (T1)

Consecuencia: μ **no tiene escala propia**; su escala la fija P(d). μ→0 cuando Σd→n (todas
las |U_ii|→1 ⇒ U→diagonal ⇒ masa fuera-de-diagonal→0). El corner (1,…,1) es el punto de
obstrucción NULA, no las esquinas de permutación con ceros.

## 3. La fracción de emparejamiento forzado ρ(d) (el observable correcto)

Definimos la cantidad adimensional

    **ρ(d) := 2·μ(d) / (n − Σd)  =  μ(d) / [P(d)/2]  ∈ [0, 1].**               (def)

ρ mide qué fracción del presupuesto fuera-de-diagonal la unitariedad OBLIGA a que aparezca
como masa EMPAREJADA (|U_ij| y |U_ji| ambos no nulos), en vez de masa unilateral. ρ=0 sería
"toda la masa puede ser unilateral" (sin obstrucción); ρ=1 sería "todo forzosamente emparejado".

La obstrucción de Hilbert–Pólya-por-lazo (μ(d)>0) es EQUIVALENTE a ρ(d)>0. Y la pregunta de
escalamiento se reformula limpiamente, factorizando el efecto trivial del presupuesto:

    ¿ρ(d) tiene un piso positivo universal, o ρ→0 al crecer n?

## 4. Evidencia n=3 (exacta, `barrer` sobre el símplex de diagonales)

**T1 es casi-ajustada en todo el símplex:** μ ≈ P/2 (barrido de 802 diagonales). Ejemplos:
Σd=2.910: μ=0.0447 vs cota 0.050; Σd=0.717: μ=1.124 vs 1.14. La cota superior demostrada
explica el "confounder Σd" crudo de n=4: μ~P/2 domina la variación de μ en bruto.

**Estructura de ρ (barrido dedicado de 1967 diagonales interiores, minimizando ρ, no μ):**
- **ρ_min = 0.366 en d=(0.05,0.05,0.05)** (cerca del corner (0,0,0)); ρ mediana 0.851, ρ_max 0.999.
- **El ínfimo de ρ sobre el símplex CERRADO es 0**, alcanzado en las esquinas de permutación
  degeneradas (el 3-ciclo con diagonal (0,0,0) es una permutación pura ⇒ M=0 ⇒ ρ=0). Por tanto
  **NO hay piso positivo universal de ρ sobre TODAS las diagonales.** La obstrucción es
  diagonal-dependiente. [Corrección de una versión previa de este doc que reportó "ρ≳0.76 en n=3":
  era un artefacto de muestrear los extremos de μ, no de ρ — disciplina de confounder.]
- **Mecanismo:** ρ decae cuando la diagonal tiene entradas pequeñas (te acercas al régimen
  permutación-pura, donde la masa puede ser unilateral sin coste). Camino a la esquina (1,0,0):
  μ CRECE (0.697→0.935), ρ→1 (bloque 2×2 antidiagonal forzado); pero camino a (0,0,0): ρ→0.
- Casos ancla: vm-3 (1/2,1/3,1/5) Σd=1.033 → ρ=0.762; const-3 (0.4,0.4,0.4) Σd=1.2 → ρ=0.775.

**Consecuencia para von Mangoldt:** d_i=1/p_i acumula cada vez más entradas DIMINUTAS al crecer n
(la cola 1/p→0), empujando la familia vm justo hacia el régimen de ρ pequeño. Esto DA UNA RAZÓN
ESTRUCTURAL de por qué ρ(vm) cae rápido con n (§5), en vez de un simple hecho empírico: la propia
forma de las amplitudes de von Mangoldt la lleva al régimen de obstrucción débil.

## 5. Señal de escala (preliminar, a confirmar con la malla)

| n | familia | Σd | μ | P/2 = (n−Σd)/2 | ρ = 2μ/(n−Σd) |
|---|---------|-----|-------|----------------|----------------|
| 3 | vm | 1.033 | 0.749 | 0.983 | **0.762** |
| 4 | vm | 1.176 | 0.585 | 1.412 | **0.414** |
| 3 | const(0.4) | 1.200 | 0.697 | 0.900 | **0.775** |
| 4 | const(0.4) | 1.600 | 0.660 | 1.200 | **0.550** |
| 4 | random-0 | 0.964 | 0.510 | 1.518 | **0.336** |
| 4 | random-2 | 1.466 | 0.705 | 1.267 | **0.556** |

**ρ decrece de n=3 a n=4 en las dos familias comparables (vm: 0.76→0.41; const: 0.78→0.55).**
Si esta tendencia continúa, la obstrucción NORMALIZADA se diluye con la dimensión — lo que
tendría consecuencia directa para si un grafo cuántico FINITO puede realizar las amplitudes de
ζ (un ρ→0 diría que en dimensión alta la unitariedad estorba cada vez menos). Si ρ satura en un
piso positivo, la obstrucción es universal. **La malla [4,5,6,8,10] medirá ρ(n).** Sobre controles
(corrección 2026-07-06): C2 (Σd=2.0 fijo) fija P(d)=n−2 pero sus entradas 2/n ENCOGEN — y §4
muestra que entradas pequeñas per se bajan ρ; el control complementario es const(0.4) (entradas
fijas, Σd crece). Ninguno es puro por separado; la conclusión dimensional requiere que AMBOS
decrezcan (triangulación — y ambos decrecieron: ver `mu_scaling_report.md §2`).

## 6. Implicación para la cota inferior (O1/O2 de `cota_inferior_demostracion.md`)

El objetivo de la Ruta 2 se afila y se ACOTA a la familia relevante. Como el ínfimo de ρ sobre
TODAS las diagonales es 0 (§4, esquinas de permutación), NO existe cota "ρ(d)≥c>0" universal. La
pregunta correcta es **familia-específica**:

    **¿ρ_vm(n) = 2μ(d_vm(n))/(n−Σd_vm(n)) tiene un piso positivo, o →0 cuando n→∞?**

donde d_vm(n)=(1/2,1/3,…,1/p_n). La evidencia (§5) apunta a ρ_vm(n)→0, con razón estructural
(§4: la cola 1/p→0 lleva a vm al régimen de obstrucción débil). Si se confirma, la tesis del banco
cambia: la obstrucción-por-lazo NO bloquea von Mangoldt en dimensión alta — se diluye. Eso reabre
(no cierra) la realizabilidad de amplitudes por grafos grandes, y traslada el bloqueo a OTRA
restricción (escala/anclaje C6, o la estadística GUE), no a la masa emparejada.

La identidad espectral exacta Σ_{i<j}U_ijU_ji=(trU²−Σa_i²)/2 acota la masa CON SIGNO; el salto
pendiente sigue siendo de masa-con-signo a masa-en-valor-absoluto (el emparejamiento), pero
ahora medido como fracción del presupuesto, que es la forma en que la obstrucción es invariante.

## 7. Qué decide el usuario

- Si la malla confirma ρ(n) decreciente: ¿ajustamos c(n) y perseguimos su forma (ley), o eso
  cambia la tesis del banco (la obstrucción-por-lazo se diluye en dimensión alta)?
- Si ρ(n) satura: candidato a cota universal ρ≥c → objetivo de demostración O1/O2.
- El reporte de la malla (Ruta 1) se reexpresará en términos de ρ además de μ.
