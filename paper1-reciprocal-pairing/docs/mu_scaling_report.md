# Estudio de escalamiento de μ(d) / ρ(d) para n ≥ 4 — reporte final

_open_math_lab (banco Hilbert–Pólya). Ruta 1 (numérica). Sprint autónomo 2026-07-04._
_Reglas pre-registradas: `docs/mu_scaling_prereg.md`. Marco teórico: `docs/mu_ruta2_budget.md`._
_Evidencia numérica (cota superior de μ + positividad multi-start), NO prueba._

> **VEREDICTO: ESCENARIO C (indeterminado a n > 6).** El motor SLSQP es fiable hasta n=6 y
> colapsa en n=8,10. En el rango fiable, ρ(d) decrece con n en todas las familias (incluido el
> control C2 con Σd fijo), pero 4 puntos no distinguen ρ→0 de una meseta positiva pequeña. No se
> ajusta ninguna ley sobre n>6. La siguiente iteración requiere autorización humana y una
> herramienta que alcance n≥8 con factibilidad alta (`docs/next_step_prereg.md`).

Datos: `data/mu_scaling_results.csv`. Figuras: `docs/figures/qgraphs/{mu_scaling_vs_n,rho_vs_n,mu_C2_fixedS,mu_vs_sumd}.png`.

---

## Sección 1 — Factibilidad y restricciones

μ(d) = mín M(U) sobre {U∈U(n): |U_ii|²=d_i}, con U=expm(A) anti-hermitiana y SLSQP multi-start.
Un arranque cuenta como factible si `max_i |U_ii|²−d_i| < 1e-6`. **La factibilidad colapsa con n:**

| n | num_starts | % factible (media familias) | rango err_restricción (celdas ok) |
|---|-----------|------------------------------|-----------------------------------|
| 4 | 200 | 84 % | 1e-14 … 6e-11 |
| 5 | 250 | 92 % | 1e-13 … 9e-12 |
| 6 | 360 | 69 % | 5e-12 … 1e-7 |
| 8 | 400 | **17 %** | 6e-8 … 8e-7 (al borde de 1e-6) |
| 10 | 400 | **0.2 %** | 5e-7 … 1e-6 |

En n=10, 4 de 5 familias con diagonal fija dieron **0 arranques factibles** (`infeasible_numeric`);
vm y random-0 dieron 1, random-1 dio 3. Un ínfimo estimado a partir de 1 punto es un accidente,
no un mínimo. **Causa:** arranques anti-hermitianos aleatorios casi nunca caen cerca de la
superficie `|U_ii|²=1/p_i` con valores diminutos (1/17,1/19,…) en ~100 dimensiones, y SLSQP no
recupera la factibilidad. Infactibilidad a priori detectada correctamente (prefiltro DS): random
seed1 en n=4 (Σd=2.555). Un FALSIFIED: const_fixedS n=8 (μ=5e-6 < δ), analizado en §4.

## Sección 2 — Tendencias de μ y ρ (rango fiable n ≤ 6)

ρ = 2μ/(n−Σd) = fracción del presupuesto fuera-de-diagonal forzada a emparejarse (adimensional,
∈[0,1]; ver `mu_ruta2_budget.md`). Se añade n=3 EXACTO (de `unistochastic.barrer`, fuera de la malla).

| familia | n=3 (exacto) | n=4 | n=5 | n=6 | μ_abs (n=4→6) |
|---------|--------------|-----|-----|-----|----------------|
| **vm** ρ | 0.762 | 0.414 | 0.164 | 0.053 | 0.585→0.307→0.124 |
| **const(0.4)** ρ | 0.775 | 0.550 | 0.353 | 0.244 | 0.660→0.529→0.438 |
| **const_fixedS (C2, Σd=2.0)** ρ | 0.899 | 0.657 | 0.353 | 0.168 | 0.657→0.529→0.336 |
| **linear** ρ | — | 0.454 | 0.233 | 0.112 | 0.624→0.408→0.237 |

- **ρ decrece con n en TODAS las familias.** Figuras `mu_scaling_vs_n.png`, `rho_vs_n.png`
  (puntos huecos + banda gris = régimen n>6 no fiable, excluido de toda conclusión).
- **Control C2 y su confounder residual (corrección de integridad 2026-07-06):** C2 fija Σd=2.0
  pero sus entradas d_i=2/n ENCOGEN con n (0.5→0.2), y el barrido exacto n=3 mostró que entradas
  pequeñas *per se* empujan ρ hacia abajo (`mu_ruta2_budget.md §4`). Así que C2 aislado NO es
  "dimensión pura". La conclusión dimensional sobrevive por **triangulación**: `const(0.4)` tiene
  entradas FIJAS (Σd crece) y también decrece (ρ 0.775→0.550→0.353→0.244); `const_fixedS` tiene
  Σd FIJO (entradas encogen) y decrece (**0.899**→0.657→0.353→0.168; empieza en n=3). Ningún control
  aislado es puro; que AMBOS decrezcan es lo que sostiene "efecto dimensional real". Una versión
  previa de este reporte atribuía la pureza a C2 solo — corregido.
  - **Corrección de dato (2026-07-15):** el valor de `const_fixedS` en n=3 estaba escrito como 0.775
    (copiado por error de la fila `const(0.4)`); el correcto para d=(2/3,2/3,2/3) es **ρ=0.899**
    (μ=0.4495 vía `unistochastic.barrer`). Las dos familias sí coinciden en n=5 (2/5=0.4), de ahí
    el 0.353 compartido. La conclusión (ambas decrecen) se REFUERZA. Cazado al blindar el Paper 1.
- **Confounder Σd:** en bruto μ correlaciona con Σd (fig. `mu_vs_sumd.png`), pero es el efecto de la
  cota de presupuesto μ≲(n−Σd)/2 (teorema T1, `mu_ruta2_budget.md`); ρ lo normaliza, y el control
  decisivo es la PAREJA const/const_fixedS (triangulación de arriba), no C2 aislado.

## Sección 3 — Análisis de escalamiento de ρ_vm(n)  (sólo n≤6 fiable: n=3,4,5,6)

Tres modelos ajustados a ρ_vm = [0.762, 0.414, 0.164, 0.053] en n=[3,4,5,6]:

| modelo | parámetros | R² | AIC |
|--------|-----------|-----|-----|
| potencia  ρ = a·n^(−b) | a=16.9, b=2.80 | 0.971 | −20.6 |
| log-decay ρ = a/(ln n + b) | a=0.143, b=−0.913 | 0.911 | −16.1 |
| meseta+decay ρ = c + a·n^(−b) | **c=−0.612**, a=4.47, b=1.07 | 0.998 | −28.4 |

**NO se elige un modelo como verdadero.** Con 4 puntos:
- Potencia y log-decay ajustan casi igual de bien (R²>0.9). Ambos → 0 asintóticamente.
- El único modelo que *distinguiría* dilución de meseta (meseta+decay) devuelve **c=−0.612**, una
  meseta NEGATIVA sin sentido físico (ρ≥0): 4 puntos con 3 parámetros = 1 grado de libertad =
  sobreajuste. El modelo que decidiría la pregunta es precisamente el que no es fiable.
- **Conclusión de §3: los datos fiables (n≤6) son consistentes con ρ_vm→0, pero no permiten
  distinguirlo de una meseta positiva pequeña.** Para decidir hacen falta n=8,10,12 fiables.

## Sección 4 — Robustez del optimizador

**El motor NO se estabiliza en n≥8.** Evidencia dentro de la propia malla (400 arranques/celda):
- **Dispersión intra-celda enorme en n=8:** vm mejor μ=0.014 vs mediana 0.079 (5.6×);
  const_fixedS mejor 5e-6 vs mediana 0.190. Cuando el mejor y la mediana difieren en órdenes de
  magnitud, el "mejor" es un arranque afortunado en una región casi-degenerada, no un ínfimo robusto.
- **No-monotonía diagnóstica:** ρ_vm = 0.053 (n=6) → 0.004 (n=8) → 0.043 (n=10). Una señal real
  sería monótona; el zigzag es ruido de convergencia.
- **FALSIFIED (const_fixedS n=8, μ=5e-6):** con 109/400 factibles y dispersión de 5 órdenes, no se
  puede distinguir una solución casi-permutación real de un artefacto. Con Σd=2.0 y d_i=0.25 la
  diagonal es simétrica y admite configuraciones casi-permutación; queda como pregunta abierta, NO
  como evidencia de μ=0.
- **Sobre el test de seed-varianza (paso 1 del protocolo):** se inició una re-corrida de las
  celdas sospechosas con semillas distintas, pero se DESCARTÓ por redundante y costosa: los 400
  arranques por celda de la malla YA son una muestra multi-start grande con semilla fija, y su
  factibilidad (§1) y dispersión intra-celda son evidencia de inestabilidad más fuerte que una
  simple comparación entre 2–3 semillas. Re-correr n=10 (400 arranques × varias semillas) sólo
  reproduciría el ~0.2% de factibilidad ya observado, sin añadir información. La regla del protocolo
  ("si μ no se estabiliza al añadir reinicios, marcar como no fiable") se cumple con margen: en n=10
  ni siquiera hay suficientes arranques factibles para hablar de estabilización.

- **Chequeo M2 (pre-registrado "≥1 caso"; OMITIDO en la corrida original, ejecutado post-hoc
  2026-07-06 como corrección de integridad):** minimizando M2=Σ|U_ij·U_ji|² en vm n=4
  (min=0.101, feas 198/200), const n=4 (0.087, 197/200) y vm n=5 (0.024, 250/250) — todos
  estrictamente positivos con factibilidad ~100%. El argmin de M2 tiene M1 mayor que el mínimo
  de M1 (0.700 vs 0.585 en vm n=4), como corresponde a objetivos distintos. **La positividad en
  el rango fiable es robusta a la elección del objetivo.** Datos: `data/m2_robustez.json`.

**Los valores de n=8 y n=10 NO son confiables y no se usan para conclusiones cuantitativas.**

## Sección 5 — Conclusión provisional y pre-registro de la siguiente iteración

**Escenario detectado: C (indeterminado a n>6).** Justificación frente a los umbrales congelados:
la regla de Escenario A exige ρ_vm(n=10)<0.05 con R²>0.9 sobre datos fiables; n=10 tiene 0.2% de
factibilidad, así que su ρ NO es fiable y la condición de A no puede evaluarse. La regla de B exige
una mediana fiable en n=10 sin tendencia; tampoco hay dato fiable. Quedan los problemas de
convergencia en n=8,10 → **Escenario C** por definición.

Lo aprendido, que SÍ es sólido: (i) la cota de presupuesto μ≤(n−Σd)/2 (teorema); (ii) ρ es el
observable correcto; (iii) en el rango fiable n≤6, ρ decrece dimensionalmente (control C2), lo que
*sugiere* que la obstrucción-por-lazo se diluye — pero es indistinguible de una meseta pequeña con
los datos actuales.

**Acción (pre-registrada, SIN ejecutar):** `docs/next_step_prereg.md` propone dos herramientas que
alcanzarían n≥8 con factibilidad alta — **Opción A** (optimización Riemanniana sobre la fibra
unitaria, `pymanopt`) y **Opción B** (reparametrización U=D₁·V·D₂) — ambas con la misma compuerta de
validación (reproducir barrer(n=3) y n≤6 antes de creer n≥8). La tesis del banco podría reorientarse
(si ρ_vm→0, la masa emparejada no bloquea von Mangoldt en dimensión alta y el cuello pasa a otra
restricción), pero eso **requiere primero resolver la indeterminación numérica**.

---

**Los experimentos de este sprint han finalizado. La siguiente iteración requiere autorización
humana. Consulta el archivo `docs/next_step_prereg.md`.**
