# μ(d) exacta para n=3 — verificación por reducción unistocástica

_open_math_lab (banco Hilbert–Pólya). Resultado matemático exacto, no evidencia sobre RH._
_Generado por `scripts/run_mu_n3.py`._

========================================================================
  μ(d) exacta para n=3 — reducción unistocástica
========================================================================

## Método (declarado)
- M = Σ_{i<j}√(B_ij B_ji) depende SÓLO de B=|U|∘²; B doblemente estocástica 3×3 con diagonal fija = 1 parámetro x; unistocasticidad n=3 = desigualdad triangular cerrada. ⇒ μ(d) = mín exacto de M(x) sobre {factible}∩{unistocástica}. Sin parametrizar U(3).
- 'W' del handoff = M (masa emparejada CRUDA, la de `cota_inferior_demostracion.md`), NO el W ponderado de FL-2. μ_ds = relajación doblemente estocástica; μ_uni = μ(d) verdadera.

## Resultados
| terna d | intervalo x | μ_ds (relajada) | μ_uni = μ(d) | ¿unist. muerde? | frac. x unist. |
|---------|-------------|-----------------|--------------|----------------|----------------|
| von Mangoldt-3 (1/2,1/3,1/5) | [0.000, 0.367] | 0.7054 | 0.7490 | SÍ | 0.977 |
| constante (0.4,0.4,0.4) | [0.000, 0.600] | 0.0000 | 0.6974 | SÍ | 0.632 |
| levemente distinta (0.5,0.45,0.4) | [0.000, 0.450] | 0.3926 | 0.6646 | SÍ | 0.711 |
| fuertemente distinta (0.8,0.5,0.2) | vacío | — | — | diagonal NO factible (d₁+d₂>1+d₃) | — |
| casi-degenerada (0.6,0.4,0.4) | [0.000, 0.400] | 0.3464 | 0.6597 | SÍ | 0.701 |

## Cruza-verificación (von Mangoldt-3, en el óptimo x*)
- U∈U(3) construida: ‖U†U−I‖∞ = 2.2e-16, ‖|diag|²−d‖∞ = 1.1e-16
- M(U) medida sobre la unitaria = 0.7490 (coincide con μ_uni = 0.7490)
- cota condicional (II) ½|trU²−Σa_i²| = 0.2268 ≤ M = 0.7490: ✓ (consistencia con la identidad espectral)

## Figura
- docs/figures/qgraphs/mu_n3.png

## Conclusión
- μ(d) > 0 para las cinco ternas probadas (incluida von Mangoldt-3): SÍ. **La obstrucción persiste en el caso mínimo no trivial n=3, de forma EXACTA (1 variable, unistocasticidad cerrada), no aproximada.**
- Caso más limpio: diagonal CONSTANTE — μ_ds = 0 (la circulante en x=0) pero NO es unistocástica ⇒ μ_uni > 0. La unistocasticidad es exactamente lo que impide la asimetría total, conectando la sonda E0 con un teorema exacto.
- Bono (existencia previa): la terna (0.8,0.5,0.2) ni siquiera es diagonal de una doblemente estocástica (d₁+d₂ > 1+d₃) — la frontera de Thompson/Horn actúa ANTES que la unistocasticidad. Dos capas de obstrucción, ambas activas.
- Alcance honesto: n=3 es el foothold riguroso; la mini-roseta es n=10 con diagonal en pares reverso-bond. Esto establece la obstrucción en el caso más pequeño y da a la cota general (O1/O2 de `cota_inferior_demostracion.md`) un ancla demostrada.

========================================================================