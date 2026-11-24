/-
P1Bounds: two elementary inequalities about unitary matrices.

Formalizes the "Lower bound" part of Theorem `thm:closed` (n = 3) and
Proposition `prop:lower4` (n = 4) of the note `paper1_mu_note`.
-/
import Mathlib

open Matrix

namespace P1Bounds

/-! ## Algebraic identities over a commutative ring -/

section CommRing

variable {R : Type*} [CommRing R]

/-- If `M * N = 1` then `adjugate M = det M • N`. -/
theorem adjugate_eq_det_smul_of_mul_eq_one {n : Type*} [Fintype n] [DecidableEq n]
    (M N : Matrix n n R) (h : M * N = 1) : adjugate M = det M • N := by
  have h' : N * M = 1 := mul_eq_one_comm.mp h
  calc adjugate M = adjugate M * (M * N) := by rw [h, Matrix.mul_one]
    _ = (adjugate M * M) * N := by rw [Matrix.mul_assoc]
    _ = det M • N := by rw [Matrix.adjugate_mul, Matrix.smul_mul, Matrix.one_mul]

/-- Cofactor identity for a `3 × 3` matrix with a (two-sided) inverse `N`:
the principal `2 × 2` minor on `{i, j}` equals `det M * N k k`. -/
theorem cofactor_three_of_mul_eq_one (M N : Matrix (Fin 3) (Fin 3) R) (h : M * N = 1)
    {i j k : Fin 3} (hij : i ≠ j) (hik : i ≠ k) (hjk : j ≠ k) :
    M i i * M j j - M i j * M j i = det M * N k k := by
  have hadj := adjugate_eq_det_smul_of_mul_eq_one M N h
  rw [Matrix.adjugate_fin_three] at hadj
  have h0 := congrFun (congrFun hadj 0) 0
  have h1 := congrFun (congrFun hadj 1) 1
  have h2 := congrFun (congrFun hadj 2) 2
  simp at h0 h1 h2
  fin_cases i <;> fin_cases j <;> fin_cases k <;> simp at hij hik hjk <;>
    simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk] <;>
    first
      | linear_combination h0
      | linear_combination h1
      | linear_combination h2

/-- Jacobi's complementary-minor identity for a `4 × 4` matrix with inverse `N`,
for the partition `{0,1} ⊔ {2,3}`.  Proof: multiply `M` by the block-triangular
matrix `T = [[1, N₁₂], [0, N₂₂]]`; then `M * T = [[M₁₁, 0], [M₂₁, 1]]`, and take
determinants. -/
theorem jacobi_block_four (M N : Matrix (Fin 4) (Fin 4) R) (h : M * N = 1) :
    M 0 0 * M 1 1 - M 0 1 * M 1 0 = det M * (N 2 2 * N 3 3 - N 2 3 * N 3 2) := by
  let T : Matrix (Fin 4) (Fin 4) R :=
    !![1, 0, N 0 2, N 0 3;
       0, 1, N 1 2, N 1 3;
       0, 0, N 2 2, N 2 3;
       0, 0, N 3 2, N 3 3]
  let P : Matrix (Fin 4) (Fin 4) R :=
    !![M 0 0, M 0 1, 0, 0;
       M 1 0, M 1 1, 0, 0;
       M 2 0, M 2 1, 1, 0;
       M 3 0, M 3 1, 0, 1]
  have e : ∀ a b : Fin 4,
      M a 0 * N 0 b + M a 1 * N 1 b + M a 2 * N 2 b + M a 3 * N 3 b
        = (1 : Matrix (Fin 4) (Fin 4) R) a b := by
    intro a b
    rw [← h, Matrix.mul_apply, Fin.sum_univ_four]
  have hMT : M * T = P := by
    ext a b
    have ea := e a b
    fin_cases a <;> fin_cases b <;>
      simp [Matrix.mul_apply, Fin.sum_univ_four, T, P] at ea ⊢ <;>
      linear_combination ea
  have hT : det T = N 2 2 * N 3 3 - N 2 3 * N 3 2 := by
    simp [T, Matrix.det_succ_row_zero, Fin.sum_univ_succ, Fin.succAbove]
    ring
  have hP : det P = M 0 0 * M 1 1 - M 0 1 * M 1 0 := by
    simp [P, Matrix.det_succ_row_zero, Fin.sum_univ_succ, Fin.succAbove]
    ring
  have hdet := congrArg det hMT
  rw [Matrix.det_mul, hT, hP] at hdet
  exact hdet.symm

/-- Jacobi's complementary-minor identity for a `4 × 4` matrix with inverse `N`,
for an arbitrary splitting `{i, j} ⊔ {k, l}` of the index set. -/
theorem jacobi_four_of_mul_eq_one (M N : Matrix (Fin 4) (Fin 4) R) (h : M * N = 1)
    {i j k l : Fin 4} (hij : i ≠ j) (hik : i ≠ k) (hil : i ≠ l)
    (hjk : j ≠ k) (hjl : j ≠ l) (hkl : k ≠ l) :
    M i i * M j j - M i j * M j i = det M * (N k k * N l l - N k l * N l k) := by
  have hinj : Function.Injective ![i, j, k, l] := by
    intro a b hab
    fin_cases a <;> fin_cases b <;> simp_all
  have hbij : Function.Bijective ![i, j, k, l] :=
    Finite.injective_iff_bijective.mp hinj
  let σ : Fin 4 ≃ Fin 4 := Equiv.ofBijective _ hbij
  have hM : (M.submatrix σ σ) * (N.submatrix σ σ) = 1 := by
    rw [Matrix.submatrix_mul_equiv, h, Matrix.submatrix_one_equiv]
  have key := jacobi_block_four (M.submatrix σ σ) (N.submatrix σ σ) hM
  rw [Matrix.det_submatrix_equiv_self] at key
  simpa [σ] using key

end CommRing

/-! ## Unitary matrices -/

section Unitary

open scoped ComplexConjugate

variable {n : Type*} [Fintype n] [DecidableEq n]

theorem unitary_mul_star (U : Matrix.unitaryGroup n ℂ) :
    (U : Matrix n n ℂ) * star (U : Matrix n n ℂ) = 1 :=
  mul_eq_one_comm.mp (Matrix.UnitaryGroup.star_mul_self U)

/-- The determinant of a unitary matrix has modulus one. -/
theorem norm_det_unitary (U : Matrix.unitaryGroup n ℂ) :
    ‖(U : Matrix n n ℂ).det‖ = 1 := by
  have h := congrArg det (unitary_mul_star U)
  rw [Matrix.det_mul, Matrix.det_one, Matrix.star_eq_conjTranspose,
    Matrix.det_conjTranspose] at h
  have h2 : ‖(U : Matrix n n ℂ).det‖ * ‖(U : Matrix n n ℂ).det‖ = 1 := by
    have := congrArg norm h
    simpa using this
  nlinarith [norm_nonneg (U : Matrix n n ℂ).det]

/-- Reverse triangle inequality in the form used twice below. -/
theorem abs_norm_sub_norm_le_of_unit (p q d : ℂ) (hd : ‖d‖ = 1) :
    |‖p‖ - ‖q‖| ≤ ‖p - d * conj q‖ := by
  have h := abs_norm_sub_norm_le p (d * conj q)
  simpa [hd] using h

end Unitary

/-! ## Dimension three (Theorem `thm:closed`, lower bound) -/

section Three

open scoped ComplexConjugate

/-- **Cofactor identity** (paper eq. `eq:cofactor`). For `U ∈ U(3)` and
`{i, j, k} = {0, 1, 2}`:  `U_ii U_jj - U_ij U_ji = det U * conj U_kk`. -/
theorem cofactor_identity_three (U : Matrix.unitaryGroup (Fin 3) ℂ)
    {i j k : Fin 3} (hij : i ≠ j) (hik : i ≠ k) (hjk : j ≠ k) :
    U.val i i * U.val j j - U.val i j * U.val j i = U.val.det * conj (U.val k k) := by
  have h := cofactor_three_of_mul_eq_one U.val (star U.val) (unitary_mul_star U) hij hik hjk
  simpa [Matrix.star_apply] using h

/-- **Pairwise lower bound, `n = 3`** (Theorem `thm:closed`, "More precisely" part):
`|U_ij U_ji| ≥ | |U_ii| |U_jj| - |U_kk| |`. -/
theorem pair_lower_bound_three (U : Matrix.unitaryGroup (Fin 3) ℂ)
    {i j k : Fin 3} (hij : i ≠ j) (hik : i ≠ k) (hjk : j ≠ k) :
    |‖U.val i i‖ * ‖U.val j j‖ - ‖U.val k k‖| ≤ ‖U.val i j * U.val j i‖ := by
  have hc := cofactor_identity_three U hij hik hjk
  have hd := norm_det_unitary U
  have h := abs_norm_sub_norm_le_of_unit (U.val i i * U.val j j) (U.val k k) U.val.det hd
  have e : U.val i i * U.val j j - U.val.det * conj (U.val k k) = U.val i j * U.val j i := by
    linear_combination hc
  rw [e, norm_mul] at h
  exact h

private theorem ne01 : (0 : Fin 3) ≠ 1 := by decide
private theorem ne02 : (0 : Fin 3) ≠ 2 := by decide
private theorem ne12 : (1 : Fin 3) ≠ 2 := by decide

/-- **Summed lower bound, `n = 3`**: with `a = |U_00|`, `b = |U_11|`, `c = |U_22|`,
`∑_{i<j} |U_ij U_ji| ≥ |ab - c| + |ac - b| + |bc - a|`. -/
theorem sum_lower_bound_three (U : Matrix.unitaryGroup (Fin 3) ℂ) :
    |‖U.val 0 0‖ * ‖U.val 1 1‖ - ‖U.val 2 2‖| + |‖U.val 0 0‖ * ‖U.val 2 2‖ - ‖U.val 1 1‖|
        + |‖U.val 1 1‖ * ‖U.val 2 2‖ - ‖U.val 0 0‖|
      ≤ ‖U.val 0 1 * U.val 1 0‖ + ‖U.val 0 2 * U.val 2 0‖ + ‖U.val 1 2 * U.val 2 1‖ := by
  have h01 := pair_lower_bound_three U (i := 0) (j := 1) (k := 2) ne01 ne02 ne12
  have h02 := pair_lower_bound_three U (i := 0) (j := 2) (k := 1) ne02 ne01 ne12.symm
  have h12 := pair_lower_bound_three U (i := 1) (j := 2) (k := 0) ne12 ne01.symm ne02.symm
  linarith

end Three

/-! ## Dimension four (Proposition `prop:lower4`) -/

section Four

open scoped ComplexConjugate

/-- **Jacobi complementary-minor identity, `n = 4`** (paper eq. `eq:jacobi4`).
For `U ∈ U(4)` and `{i, j} ⊔ {k, l} = {0, 1, 2, 3}`:
`U_ii U_jj - U_ij U_ji = det U * conj (U_kk U_ll - U_kl U_lk)`. -/
theorem jacobi_two_by_two_four (U : Matrix.unitaryGroup (Fin 4) ℂ)
    {i j k l : Fin 4} (hij : i ≠ j) (hik : i ≠ k) (hil : i ≠ l)
    (hjk : j ≠ k) (hjl : j ≠ l) (hkl : k ≠ l) :
    U.val i i * U.val j j - U.val i j * U.val j i
      = U.val.det * conj (U.val k k * U.val l l - U.val k l * U.val l k) := by
  have h := jacobi_four_of_mul_eq_one U.val (star U.val) (unitary_mul_star U)
    hij hik hil hjk hjl hkl
  rw [h]
  simp only [star_apply, RCLike.star_def, map_sub, map_mul]
  ring

/-- **Pairwise lower bound, `n = 4`** (Proposition `prop:lower4`, first display):
`|U_ij U_ji| + |U_kl U_lk| ≥ | |U_ii| |U_jj| - |U_kk| |U_ll| |`. -/
theorem pair_lower_bound_four (U : Matrix.unitaryGroup (Fin 4) ℂ)
    {i j k l : Fin 4} (hij : i ≠ j) (hik : i ≠ k) (hil : i ≠ l)
    (hjk : j ≠ k) (hjl : j ≠ l) (hkl : k ≠ l) :
    |‖U.val i i‖ * ‖U.val j j‖ - ‖U.val k k‖ * ‖U.val l l‖|
      ≤ ‖U.val i j * U.val j i‖ + ‖U.val k l * U.val l k‖ := by
  have hc := jacobi_two_by_two_four U hij hik hil hjk hjl hkl
  have hd := norm_det_unitary U
  have h := abs_norm_sub_norm_le_of_unit (U.val i i * U.val j j) (U.val k k * U.val l l)
    U.val.det hd
  have e : U.val i i * U.val j j - U.val.det * conj (U.val k k * U.val l l)
      = U.val i j * U.val j i - U.val.det * conj (U.val k l * U.val l k) := by
    rw [map_sub] at hc
    linear_combination hc
  have t : ‖U.val i j * U.val j i - U.val.det * conj (U.val k l * U.val l k)‖
      ≤ ‖U.val i j * U.val j i‖ + ‖U.val k l * U.val l k‖ := by
    refine (norm_sub_le _ _).trans ?_
    rw [norm_mul U.val.det, hd, one_mul, RCLike.norm_conj]
  rw [e] at h
  simpa [norm_mul] using h.trans t

/-- **Summed lower bound, `n = 4`** (Proposition `prop:lower4`, `μ ≥ λ`):
with `a_i = |U_ii|`,
`∑_{i<j} |U_ij U_ji| ≥ |a₀a₁ - a₂a₃| + |a₀a₂ - a₁a₃| + |a₀a₃ - a₁a₂|`. -/
theorem sum_lower_bound_four (U : Matrix.unitaryGroup (Fin 4) ℂ) :
    |‖U.val 0 0‖ * ‖U.val 1 1‖ - ‖U.val 2 2‖ * ‖U.val 3 3‖|
        + |‖U.val 0 0‖ * ‖U.val 2 2‖ - ‖U.val 1 1‖ * ‖U.val 3 3‖|
        + |‖U.val 0 0‖ * ‖U.val 3 3‖ - ‖U.val 1 1‖ * ‖U.val 2 2‖|
      ≤ ‖U.val 0 1 * U.val 1 0‖ + ‖U.val 0 2 * U.val 2 0‖ + ‖U.val 0 3 * U.val 3 0‖
        + ‖U.val 1 2 * U.val 2 1‖ + ‖U.val 1 3 * U.val 3 1‖ + ‖U.val 2 3 * U.val 3 2‖ := by
  have h1 := pair_lower_bound_four U (i := 0) (j := 1) (k := 2) (l := 3)
    (by decide) (by decide) (by decide) (by decide) (by decide) (by decide)
  have h2 := pair_lower_bound_four U (i := 0) (j := 2) (k := 1) (l := 3)
    (by decide) (by decide) (by decide) (by decide) (by decide) (by decide)
  have h3 := pair_lower_bound_four U (i := 0) (j := 3) (k := 1) (l := 2)
    (by decide) (by decide) (by decide) (by decide) (by decide) (by decide)
  linarith

/-- The same statement with the left-hand side written as a sum over the six
pairs `i < j`. -/
theorem sum_lower_bound_four' (U : Matrix.unitaryGroup (Fin 4) ℂ) :
    |‖U.val 0 0‖ * ‖U.val 1 1‖ - ‖U.val 2 2‖ * ‖U.val 3 3‖|
        + |‖U.val 0 0‖ * ‖U.val 2 2‖ - ‖U.val 1 1‖ * ‖U.val 3 3‖|
        + |‖U.val 0 0‖ * ‖U.val 3 3‖ - ‖U.val 1 1‖ * ‖U.val 2 2‖|
      ≤ ∑ p ∈ (Finset.univ : Finset (Fin 4 × Fin 4)).filter (fun p => p.1 < p.2),
          ‖U.val p.1 p.2 * U.val p.2 p.1‖ := by
  have hsum : ∑ p ∈ (Finset.univ : Finset (Fin 4 × Fin 4)).filter (fun p => p.1 < p.2),
        ‖U.val p.1 p.2 * U.val p.2 p.1‖
      = ‖U.val 0 1 * U.val 1 0‖ + ‖U.val 0 2 * U.val 2 0‖ + ‖U.val 0 3 * U.val 3 0‖
        + ‖U.val 1 2 * U.val 2 1‖ + ‖U.val 1 3 * U.val 3 1‖ + ‖U.val 2 3 * U.val 3 2‖ := by
    have hS : (Finset.univ : Finset (Fin 4 × Fin 4)).filter (fun p => p.1 < p.2)
        = {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)} := by decide
    rw [hS]
    simp [Finset.sum_insert]
    ring
  rw [hsum]
  exact sum_lower_bound_four U

end Four

end P1Bounds

/-! ## Axiom audit -/

#print axioms P1Bounds.cofactor_identity_three
#print axioms P1Bounds.pair_lower_bound_three
#print axioms P1Bounds.sum_lower_bound_three
#print axioms P1Bounds.jacobi_two_by_two_four
#print axioms P1Bounds.pair_lower_bound_four
#print axioms P1Bounds.sum_lower_bound_four
#print axioms P1Bounds.sum_lower_bound_four'
