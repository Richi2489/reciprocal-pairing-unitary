# Pre-registro de la siguiente iteración — μ(d)/ρ(d) para n > 6

_open_math_lab (banco Hilbert–Pólya). Escrito al cierre del sprint autónomo 2026-07-04._
_REQUIERE AUTORIZACIÓN HUMANA antes de ejecutar. Nada de esto se ha corrido._

## Veredicto del sprint anterior: ESCENARIO C (indeterminado a n > 6)

La malla `[4,5,6,8,10]` (`docs/mu_scaling_report.md`) dio evidencia **fiable sólo hasta n=6**.
El motor SLSQP sobre U=expm(A) con restricción dura `|U_ii|²=d_i` colapsa en dimensión alta:

| n | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|----|
| % arranques factibles (media familias) | 84 | 92 | 69 | 17 | **0.2** |

En n=10, 4 de 5 familias dieron 0 arranques factibles; las demás 1–3 (un ínfimo de 1 punto es un
accidente, no un mínimo). ρ_vm es no-monótono en n=8,10 (0.053→0.004→0.043), firma de ruido de
convergencia, no de señal. **Por tanto no se ajustó ninguna ley de escala sobre n>6** (regla
pre-registrada: no forzar ajuste sobre datos que el optimizador no garantiza como ínfimo).

Lo que SÍ quedó (rango fiable n≤6, con n=3 exacto de `barrer`): ρ decrece en TODAS las familias, y
el control C2 (Σd=2.0 fijo) también (ρ: 0.657→0.353→0.168 en n=4,5,6), así que el decaimiento es
**efecto dimensional real**, no dilución de la masa diagonal. Pero 4 puntos NO distinguen ρ→0 de
una meseta positiva pequeña (potencia R²=0.947 y exponencial R²=0.983 ajustan casi igual de bien).

**Causa raíz diagnosticada:** arranques anti-hermitianos aleatorios casi nunca caen cerca de la
superficie `|U_ii|²=1/p_i` cuando los valores son diminutos (1/17, 1/19, …) en ~100 dimensiones;
SLSQP no recupera la factibilidad desde ahí. Es un problema de optimización no convexa con
restricción dura y jacobiano numérico inestable, esperable — no un bug.

## Objetivo de la siguiente iteración

Alcanzar n∈{8,10,12} con factibilidad alta (>80%) para decidir entre **dilución (ρ→0)** y
**meseta (ρ≥c>0)** en la familia von Mangoldt. La pregunta científica es la misma de
`docs/mu_ruta2_budget.md §6`; sólo falta la herramienta que llegue a esas dimensiones.

## Opciones técnicas (elige el usuario; todas requieren autorización)

### Opción C — Warm-start sobre la fibra con código YA VALIDADO del repo (RECOMENDADA)
- **La causa raíz del colapso no fue el optimizador sino los ARRANQUES**: A anti-hermitianas
  aleatorias caen lejos de la fibra `|U_ii|²=d_i` y SLSQP no la recupera en ~100 dims. Pero el repo
  YA tiene un constructor de puntos SOBRE la fibra: `qgraphs/construccion.py::unitaria_con_diagonal`
  (Thompson/Horn + proyecciones alternadas, validado en Fase 4 First Light). Pasando la diagonal
  prescrita `√d_i` (real) se obtiene U₀ con `|U₀_ii|²=d_i` exacto; semillas distintas dan warm-starts
  DIVERSOS sobre la fibra. SLSQP arrancando factible solo tiene que *mantenerse* cerca de la fibra
  mientras baja M — un problema mucho mejor condicionado.
- Coste: cero dependencias nuevas, ~20 líneas (convertir U₀ → parámetros A=logm(U₀) o restar la
  optimización directamente en U con retracción polar). Nota técnica: `expm` es suprayectiva, así
  que `A₀ = logm(U₀)` (rama principal) siempre existe; verificar `expm(A₀)≈U₀` antes de usar.
- Por qué es preferible: A (pymanopt) instala dependencia + posible retracción custom; B (D₁VD₂)
  requiere derivación nueva. C reutiliza código con tests del propio repo.
- Smoke-test pre-registrado (barato, ~min): n=8 vm con 40 warm-starts; éxito si factibilidad ≥80%
  y μ ≤ el mejor de la malla (0.014). Si pasa, correr la malla n∈{8,10,12} completa con warm-start.

### Opción A — Optimización Riemanniana sobre la fibra unitaria
- Migrar a `pymanopt` (o `geoopt`) optimizando sobre la variedad `{U∈U(n): |U_ii|²=d_i}`, que
  respeta la restricción POR CONSTRUCCIÓN (retracciones sobre la variedad), eliminando las
  violaciones de restricción y el colapso de factibilidad.
- Requiere: instalar `pymanopt` + un backend de autodiff (`autograd`/`jax`) en el venv (hoy
  ausentes); parametrizar la fibra (p.ej. Stiefel/producto con ligadura de módulos diagonales) o
  penalización geométrica sobre U(n). Riesgo: la fibra `|U_ii|²=d_i` no es una subvariedad
  estándar lista en la librería; puede requerir retracción custom o Lagrangiano aumentado sobre
  la variedad U(n) plena.
- Pre-registro de validación: reproducir el μ EXACTO de `barrer` en n=3 (< 1e-3) y los valores
  fiables de n≤6 (< 2%) ANTES de creer cualquier número en n≥8 (misma compuerta que ya usamos).

### Opción B — Reparametrización U = D₁·V·D₂
- `V` unitaria real (ortogonal) y `D₁,D₂` diagonales de fase. Como `|U_ii| = |Σ_k (D₁)_i V_ik ...|`
  ... las fases D₁,D₂ dan grados de libertad que ayudan a fijar los módulos diagonales con menos
  restricciones duras (potencialmente reduciendo n restricciones a un subconjunto), y el bloque
  ortogonal V es un espacio de optimización mejor condicionado que expm(anti-hermitiana).
- Requiere: derivar cómo `|U_ii|²` depende de (V, D₁, D₂) y qué restricciones quedan realmente
  duras; validar contra n=3 exacto igual que arriba.
- Riesgo: puede que no elimine todas las restricciones duras; a analizar antes de implementar.

## Reglas congeladas para la siguiente iteración (heredadas + nuevas)

- Compuerta de validación obligatoria: reproducir barrer(n=3) < 1e-3 y n≤6 < 2% antes de reportar n≥8.
- Umbral de factibilidad para declarar un n "fiable": ≥ 80% de arranques factibles (err < 1e-6).
- Umbrales de escenario A/B/C SIN CAMBIO: ρ_vm(n=10) < 0.05 y R²>0.9 (decaimiento) → A;
  ρ_vm(n=10) mediana > 0.05 sin tendencia → B; convergencia dudosa → C.
- No ampliar a n=16 hasta tener n=10,12 fiables (evita quemar cómputo en una herramienta no validada).

## Estado

Los experimentos de este sprint han FINALIZADO. La siguiente iteración requiere autorización
humana y una decisión Opción A vs B. Consulta `docs/mu_scaling_report.md` para el diagnóstico
completo y `docs/mu_ruta2_budget.md` para el marco teórico (ρ y la cota de presupuesto).
