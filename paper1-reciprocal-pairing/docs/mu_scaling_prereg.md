# Pre-registro — Estudio de escalamiento de μ(d) para n ≥ 4

_open_math_lab (banco Hilbert–Pólya). Reglas fijadas ANTES de ver resultados._
_Contexto: n=2 cerrado (μ(d,d)=1−d); n=3 cerrado y EXACTO por reducción unistocástica
(`docs/mu_d_n3_verification.md`, `qgraphs/unistochastic.py::barrer`). Aquí extendemos a n≥4,
donde la reducción de 1 variable ya NO aplica y hay que optimizar sobre U(n) directamente._

## Objeto

μ(d) = mín_{U∈U(n)} M(U) sujeto a |U_ii|² = d_i, con
M(U) = M1(U) = Σ_{i<j} |U_ij·U_ji|  (masa emparejada cruda; misma de `cota_inferior_demostracion.md`).

**d NO se normaliza.** d_i = |U_ii|² = B_ii son entradas diagonales de una doblemente
estocástica, cada una en [0,1] de forma independiente; Σd_i = tr(B) es libre en [0,n].
Para von Mangoldt d_i = A_{p_i}² = 1/p_i (crudo). Normalizar a suma 1 destruiría la
conexión con ζ y contradiría el trabajo de n=3 (allí Σ(1/2,1/3,1/5)=1.033, no 1).

## Pregunta (reformulada 2026-07-03)

No sólo "¿μ(d)>0 en n=10?" sino **cómo escala y de qué depende** la obstrucción al crecer n.

## Familias de diagonales (sin normalizar)

1. **von Mangoldt** `vm`: d_i = 1/p_i, p_i el i-ésimo primo (2,3,5,7,11,…). La diagonal objetivo.
2. **constante** `const`: d_i = c fijo. Se usa c = 0.4 (mismo valor que el caso limpio de n=3)
   por defecto; el valor NO depende de n (así Σd_i = 0.4·n crece linealmente — ver confounder).
3. **lineal decreciente** `linear`: d_i ∝ (n−i+1), escalada para que max d_i = 0.5 (rango fijo,
   independiente de n, para no colapsar entradas a 0 con n grande).
4. **aleatoria** `random`: d_i ~ Uniforme(0,1) i.i.d., 3 realizaciones (seeds 0,1,2). SIN normalizar.

## Confounder Σd_i (control obligatorio)

"μ vs n" mezcla DOS efectos: la dimensión n y la masa diagonal total Σd_i, que cambia con n de
forma distinta por familia. Un ajuste μ~n^{−α} puede ser artefacto de que la diagonal encoge/crece.
**Controles pre-registrados:**
- (C1) Reportar μ contra n Y contra Σd_i (dos columnas siempre).
- (C2) Barrido con Σd_i FIJADO: familia `const` con c(n) elegido para mantener Σd_i = S₀ constante
  (S₀ = 2.0) a través de n ∈ {4,5,6,8,10}. Si μ(const, Σ=2.0) sigue variando con n, el efecto es
  dimensional puro (no de magnitud). Esta es la brújula para la Ruta 2.

## Motor de optimización

- Parametrización: U = expm(A), A anti-hermitiana (A†=−A), n² parámetros reales. exp es
  suprayectiva sobre U(n) ⇒ todo unitario es alcanzable; multi-start cubre el paisaje.
- Restricción dura: g_i = |U_ii|² − d_i = 0, i=1..n (SLSQP, restricciones de igualdad).
- Objetivo primario M1; robustez con M2(U)=Σ_{i<j}|U_ij U_ji|² en ≥1 caso.
- Multi-start: num_starts arranques desde A aleatorias (rng con seed fijo), se guarda el mejor
  FACTIBLE. num_starts = 200 en la fase de calibración.
  **AJUSTE OPERATIVO (documentado 2026-07-03, tras pasar las compuertas n=3/n=4):** en la malla
  grande num_starts = min(max(200, 10·n²), 400). La cota 10·n² daba 1000 arranques en n=10, coste
  excesivo (~overnight) cuando la validación n=4 mostró que el mínimo global se encuentra con
  MUCHOS arranques (p5=mediana=μ en `const`); el tope de 400 recorta n=8/10 sin pérdida esperada
  de fidelidad. Cambio de eficiencia, no de criterio; los umbrales ε₀/δ no se tocan.
- Método base: **SLSQP** (no hay pymanopt/geoopt/autograd instalados). Riemanniano queda como
  fallback documentado si SLSQP se degrada en n≥6 (requeriría instalar pymanopt).

## Tolerancias (realistas para SLSQP — desviación declarada)

El handoff pedía error de restricción ≤ 1e-12. SLSQP con restricciones de igualdad converge
típicamente a ~1e-8, no 1e-12; ese umbral descartaría TODO. **Umbral de aceptación fijado en
`tol_feas = 1e-6`** sobre max_i|U_ii|²−d_i, y se REPORTA el error real de cada mínimo aceptado.
Un arranque con error > 1e-6 no cuenta como factible.

## Detección de infactibilidad (obligatoria)

Como en n=3 `(0.8,0.5,0.2)` no era realizable, para n≥4 hay diagonales sin ninguna U.
- Prefiltro NECESARIO (no suficiente) de doblemente-estocástica: para todo i,
  (1−d_i) ≤ Σ_{j≠i}(1−d_j), equivalente a Σd ≤ 2·d_i + n − 2 ∀i. Reproduce la infactibilidad
  de (0.8,0.5,0.2) en n=3. (Unistocasticidad es más estricta y desconocida para n≥4; el
  optimizador es el árbitro final.)
- Regla de decisión: si el prefiltro falla → **infactible (a priori)**. Si pasa el prefiltro pero
  NINGÚN arranque logra error < tol_feas → **infactible (numérico)**. En ambos casos se reporta
  como infactible, NO como μ grande.

## Reglas de decisión (go/no-go)

- **Positividad:** para un (n, familia), si el percentil-5 de los mínimos factibles > ε₀ = 0.05
  → evidencia de μ(d) > 0.
- **Falsación:** si algún mínimo factible cae por debajo de δ = 1e-3 → μ(d) puede ser 0
  (la obstrucción NO llega a esa dimensión). Resultado negativo honesto, igual de válido.
- **Escalamiento:** ajustar μ (mediana y percentil-5) vs n a ley potencial μ~n^{−α} y a
  exponencial μ~e^{−βn}; reportar R² de ambas. NO sobreinterpretar con 5 puntos.

## Estrategia STAGED (orden estricto)

1. **Cross-check n=3:** motor SLSQP sobre vm-3 (1/2,1/3,1/5) y const-3 (0.4,0.4,0.4). Debe
   reproducir μ_uni de `barrer` (0.7490 y 0.6974) a < 1e-3. Además μ_SLSQP ≥ μ_uni − tol
   (μ_uni es el mínimo global; el optimizador no puede bajar de él salvo violando restricción).
   **Si falla, el motor no es de fiar y NO se procede.**
2. **Validación n=4:** las 4 familias, 200 arranques. Verificar error de restricción bajo y
   estabilidad de mínimos (dispersión entre arranques).
3. **Sólo si 1 y 2 pasan:** lanzar la malla [4,5,6,8,10] (+12 si el coste lo permite).

## Entregables

- `qgraphs/mu_general.py` (motor), `tests/test_mu_scaling.py` (familias + cross-check n=3).
- `scripts/run_mu_scaling.py` (orquesta, escribe CSV + reporte).
- `data/mu_scaling_results.csv`, `docs/mu_scaling_report.md`.

## Alcance honesto

La Ruta 1 da una cota SUPERIOR de μ(d) (multi-start) + evidencia de positividad; NO prueba
μ(d)>0 (un valle a 0 podría escaparse). La prueba es la Ruta 2 (O1/O2 analíticas). El valor de
este estudio es dar a la Ruta 2 la forma funcional de la conjetura (¿μ constante? ¿μ~1/n?).
