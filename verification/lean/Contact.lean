import Mathlib.Analysis.SpecialFunctions.Trigonometric.Inverse
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.Calculus.ParametricIntegral
import Mathlib.Analysis.Calculus.LHopital
import Mathlib.Analysis.Calculus.ContDiff.Deriv
import Mathlib.Tactic

set_option autoImplicit false
namespace MRC
open Real Set

noncomputable def circularGap (h R t : ℝ) := h + R * cos t
noncomputable def noncircularGap (h R q t : ℝ) := h + R * cos t + q * cos (2*t)
noncomputable def compression (k w g : ℝ) := k*g^2/2-w
noncomputable def adhesion (w a g : ℝ) := -w*(1-(g/a)^2)^2
noncomputable def adhesionSlope (w a g : ℝ) := 4*w*g/a^2*(1-(g/a)^2)
noncomputable def potential (k w a g : ℝ) :=
  if g ≤ 0 then compression k w g else if g < a then adhesion w a g else 0
noncomputable def slope (k w a g : ℝ) :=
  if g ≤ 0 then k*g else if g < a then adhesionSlope w a g else 0

theorem compression_derivative (k w g : ℝ) :
    HasDerivAt (compression k w) (k*g) g := by
  unfold compression
  convert (((hasDerivAt_id g).pow 2).const_mul k |>.div_const 2 |>.sub_const w) using 1 <;> (try rfl) <;> simp [id] <;> ring

theorem adhesion_derivative (w a g : ℝ) (ha : a ≠ 0) :
    HasDerivAt (adhesion w a) (adhesionSlope w a g) g := by
  unfold adhesion adhesionSlope
  have H := (((hasDerivAt_id g).div_const a).pow 2).const_sub 1 |>.pow 2 |>.const_mul (-w)
  convert H using 1 <;> (try rfl) <;> simp [id] <;> field_simp <;> ring

theorem branch_joins (k w a : ℝ) (ha : a ≠ 0) :
    compression k w 0 = adhesion w a 0 ∧ adhesion w a a = 0 ∧
    adhesionSlope w a 0 = 0 ∧ adhesionSlope w a a = 0 := by
  simp [compression, adhesion, adhesionSlope, ha]

theorem branch_curvature_at_left (w a : ℝ) (ha : a ≠ 0) :
    HasDerivAt (adhesionSlope w a) (4*w/a^2) 0 := by
  unfold adhesionSlope
  have H := ((((hasDerivAt_id (0 : ℝ)).const_mul (4*w)).div_const (a^2)).mul
    ((((hasDerivAt_id (0 : ℝ)).div_const a).pow 2).const_sub 1))
  convert H using 1 <;> (try rfl) <;> simp [id]

theorem branch_curvature_at_right (w a : ℝ) (ha : a ≠ 0) :
    HasDerivAt (adhesionSlope w a) (-8*w/a^2) a := by
  unfold adhesionSlope
  have H := ((((hasDerivAt_id a).const_mul (4*w)).div_const (a^2)).mul
    ((((hasDerivAt_id a).div_const a).pow 2).const_sub 1))
  convert H using 1 <;> (try rfl) <;> simp [id] <;> field_simp <;> ring

theorem circular_interior_root (h R B : ℝ) (hR : 0 < R)
    (hB : |(B-h)/R| < 1) :
    arccos ((B-h)/R) ∈ Ioo 0 π ∧ circularGap h R (arccos ((B-h)/R)) = B := by
  rw [abs_lt] at hB
  constructor
  · exact ⟨arccos_pos.mpr hB.2, arccos_lt_pi.mpr hB.1⟩
  · unfold circularGap
    rw [cos_arccos (le_of_lt hB.1) (le_of_lt hB.2)]
    field_simp
    <;> ring

theorem circular_interior_necessary (h R B t : ℝ) (hR : 0 < R)
    (ht : t ∈ Ioo 0 π) (he : circularGap h R t = B) : |(B-h)/R| < 1 := by
  have he' : (B-h)/R = cos t := (div_eq_iff (ne_of_gt hR)).2 (by unfold circularGap at he; nlinarith)
  rw [he', abs_lt]
  have h1 := strictAntiOn_cos (show π ∈ Icc (0:ℝ) π from ⟨pi_pos.le, le_rfl⟩)
    (show t ∈ Icc (0:ℝ) π from ⟨ht.1.le, ht.2.le⟩)
  have hl := strictAntiOn_cos (show t ∈ Icc (0:ℝ) π from ⟨ht.1.le, ht.2.le⟩)
    (show π ∈ Icc (0:ℝ) π from ⟨pi_pos.le, le_rfl⟩) ht.2
  have hu := strictAntiOn_cos (show (0:ℝ) ∈ Icc (0:ℝ) π from ⟨le_rfl, pi_pos.le⟩)
    (show t ∈ Icc (0:ℝ) π from ⟨ht.1.le, ht.2.le⟩) ht.1
  simpa using And.intro hl hu

theorem circular_unique (h R : ℝ) (hR : 0 < R) :
    StrictAntiOn (circularGap h R) (Icc 0 π) := by
  intro x hx y hy hxy
  unfold circularGap
  linarith [mul_lt_mul_of_pos_left (strictAntiOn_cos hx hy hxy) hR]

theorem noncircular_quadratic (h R q B t : ℝ) :
    noncircularGap h R q t = B ↔ 2*q*(cos t)^2+R*cos t+(h-q-B) = 0 := by
  unfold noncircularGap
  rw [cos_two_mul]
  constructor <;> intro H <;> nlinarith

theorem noncircular_monotone (h R q : ℝ) (hq : 0 ≤ q) (hR : 4*q < R) :
    StrictAntiOn (noncircularGap h R q) (Icc 0 π) := by
  intro x hx y hy hxy
  have hc := strictAntiOn_cos hx hy hxy
  have hxx := neg_one_le_cos x
  have hyy := neg_one_le_cos y
  have H : 0 < R+2*q*(cos x+cos y) := by
    nlinarith [mul_nonneg hq (show 0 ≤ cos x+cos y+2 by linarith)]
  have H' := mul_pos (sub_pos.mpr hc) H
  unfold noncircularGap
  rw [cos_two_mul, cos_two_mul]
  nlinarith

theorem noncircular_arc_factor (Rx Ry q h t : ℝ) :
    sqrt ((deriv (fun u : ℝ => Rx*sin u) t)^2 +
      (deriv (noncircularGap h Ry q) t)^2) =
    sqrt ((Rx*cos t)^2 + (Ry*sin t + 2*q*sin (2*t))^2) := by
  have hx := (hasDerivAt_sin t).const_mul Rx
  have h2 := (hasDerivAt_cos (2*t)).comp t ((hasDerivAt_id t).const_mul 2)
  have hy := ((hasDerivAt_cos t).const_mul Ry |>.const_add h).add (h2.const_mul q)
  rw [hx.deriv]
  change sqrt ((Rx*cos t)^2 + (deriv (fun u => h+Ry*cos u+q*cos (2*u)) t)^2) = _
  have hd : HasDerivAt (fun u => h+Ry*cos u+q*cos (2*u))
      (-(Ry*sin t+2*q*sin (2*t))) t := by
    convert hy using 1 <;> (try rfl) <;> ring
  rw [hd.deriv, neg_sq]

theorem inclined_basis (L β : ℝ) (hL : 0 < L) :
    let ell := sqrt (L^2+β^2)
    let a := L/ell
    let b := β/ell
    a^2+b^2=1 ∧ L*(-b)+β*a=0 := by
  dsimp
  have hp : 0 < L^2+β^2 := by positivity
  have he := sqrt_pos.2 hp
  have he2 := sq_sqrt hp.le
  constructor
  · field_simp
    nlinarith
  · ring

theorem force_scaling (m estimate τ : ℝ) (hm : 0 < m) :
    m * estimate ≤ τ ↔ estimate ≤ τ/m := by
  rw [le_div_iff₀ hm, mul_comm estimate]

theorem primitive_cos2 (t : ℝ) :
    HasDerivAt (fun s : ℝ => s/2+sin (2*s)/4) ((cos t)^2) t := by
  have H := ((hasDerivAt_id t).div_const 2).add
    (((hasDerivAt_sin (2*t)).comp t ((hasDerivAt_id t).const_mul 2)).div_const 4)
  convert H using 1 <;> (try rfl) <;> simp [cos_two_mul] <;> ring

theorem primitive_cos3 (t : ℝ) :
    HasDerivAt (fun s : ℝ => sin s-(sin s)^3/3) ((cos t)^3) t := by
  have H := (hasDerivAt_sin t).sub (((hasDerivAt_sin t).pow 3).div_const 3)
  convert H using 1 <;> (try rfl)
  simp only [Nat.cast_ofNat, Nat.reduceSub]
  nlinarith [congrArg (fun x : ℝ => x*cos t) (sin_sq_add_cos_sq t)]

theorem cubic_gap_expansion (h R c : ℝ) :
    (h+R*c)^3 = h^3+3*h^2*R*c+3*h*R^2*c^2+R^3*c^3 := by ring

theorem endpoint_quadratic (R t : ℝ) :
    R*(1-cos t) = 2*R*(sin (t/2))^2 := by
  nth_rewrite 1 [show t = 2*(t/2) by ring]
  rw [cos_two_mul]
  nlinarith [congrArg (fun x : ℝ => R*x) (sin_sq_add_cos_sq (t/2))]

theorem local_parabolic_coefficient (L R m δ s : ℝ) (h : R*s^2/2=δ) :
    2*R*L*m*(δ*s - R*s^3/6) = (4/3:ℝ)*R*L*m*δ*s := by
  rw [← h]
  ring


private theorem glueLE (f g : ℝ → ℝ) (F G c x : ℝ)
    (hf : HasDerivAt f F x) (hg : HasDerivAt g G x)
    (join : x = c → f x = g x ∧ F = G) :
    HasDerivAt (fun t => if t ≤ c then f t else g t) (if x ≤ c then F else G) x := by
  rcases lt_trichotomy x c with h | rfl | h
  · rw [if_pos h.le]
    apply hf.congr_of_eventuallyEq
    filter_upwards [gt_mem_nhds h] with y hy
    simp [hy.le]
  · rw [if_pos le_rfl]
    obtain ⟨hv, hd⟩ := join rfl
    have hl : HasDerivWithinAt (fun t => if t ≤ x then f t else g t) F (Iic x) x :=
      hf.hasDerivWithinAt.congr_of_mem (fun y hy => if_pos hy) (show x ∈ Iic x from mem_Iic.mpr le_rfl)
    have hr : HasDerivWithinAt (fun t => if t ≤ x then f t else g t) F (Ici x) x := by
      rw [hd]
      apply hg.hasDerivWithinAt.congr_of_mem _ (show x ∈ Ici x from mem_Ici.mpr le_rfl)
      intro y hy
      rcases eq_or_lt_of_le (mem_Ici.mp hy) with rfl | hh
      · simpa using hv
      · simp [not_le.mpr hh]
    exact (hl.union hr).hasDerivAt (by rw [Iic_union_Ici]; exact Filter.univ_mem)
  · rw [if_neg (not_le.mpr h)]
    apply hg.congr_of_eventuallyEq
    filter_upwards [lt_mem_nhds h] with y hy
    simp [not_le.mpr hy]

theorem potential_derivative (k w a x : ℝ) (ha : 0 < a) :
    HasDerivAt (potential k w a) (slope k w a x) x := by
  have hb := branch_joins k w a ha.ne'
  have hi := glueLE (adhesion w a) (fun _ => 0) (adhesionSlope w a x) 0 a x
    (adhesion_derivative w a x ha.ne') (hasDerivAt_const x (0:ℝ))
    (fun hx => by subst x; exact ⟨hb.2.1, hb.2.2.2⟩)
  have H := glueLE (compression k w) (fun t => if t ≤ a then adhesion w a t else 0)
    (k*x) (if x ≤ a then adhesionSlope w a x else 0) 0 x
    (compression_derivative k w x) hi
    (fun hx => by subst x; simp [ha.le, hb.1, hb.2.2.1])
  convert H using 1
  · funext y
    unfold potential
    by_cases hy : y = a
    · subst y
      simp [not_le.mpr ha, hb.2.1]
    · simp [le_iff_lt_or_eq, hy]
  · unfold slope
    by_cases hy : x = a
    · subst x
      simp [not_le.mpr ha, hb.2.2.2]
    · simp [le_iff_lt_or_eq, hy]

theorem slope_continuous (k w a : ℝ) (ha : 0 < a) :
    Continuous (slope k w a) := by
  unfold slope
  have hb := branch_joins k w a ha.ne'
  apply Continuous.if
  · intro x hx
    have hh : x = 0 := by
      have he : x ∈ frontier (Iic (0:ℝ)) := hx
      simpa [frontier_Iic] using he
    subst x
    simp [ha, hb.2.2.1]
  · fun_prop
  · apply Continuous.if
    · intro x hx
      have hh : x = a := by
        have he : x ∈ frontier (Iio a) := hx
        simpa [frontier_Iio] using he
      subst x
      exact hb.2.2.2
    · unfold adhesionSlope
      fun_prop
    · exact continuous_const

theorem assembled_slope_jump_parameters :
    (20:ℝ) ≠ 4*(2/1000)/(1/10)^2 ∧ (-8:ℝ)*(2/1000)/(1/10)^2 ≠ 0 := by
  norm_num


open MeasureTheory Filter Topology

theorem translation_energy_derivative (ψ ψ' b J : ℝ → ℝ)
    (hψ : ∀ z, HasDerivAt ψ (ψ' z) z)
    (hψ' : Continuous ψ') (hb : Continuous b) (hJ : Continuous J)
    (M height : ℝ) :
    HasDerivAt (fun h => M * ∫ t in 0..π, ψ (h+b t)*J t)
      (M * ∫ t in 0..π, ψ' (height+b t)*J t) height := by
  have hc : Continuous ψ := continuous_iff_continuousAt.2 (fun z => (hψ z).continuousAt)
  let μ := volume.restrict (Ioc 0 π)
  let F : ℝ → ℝ → ℝ := fun h t => ψ (h+b t)*J t
  let D : ℝ → ℝ → ℝ := fun h t => ψ' (h+b t)*J t
  have hF (h : ℝ) : Continuous (F h) := by dsimp [F]; fun_prop
  have hD (h : ℝ) : Continuous (D h) := by dsimp [D]; fun_prop
  have hprod : Continuous (fun z : ℝ × ℝ => D z.1 z.2) := by dsimp [D]; fun_prop
  let s := Ioo (height-1) (height+1)
  have hs : s ∈ 𝓝 height := Ioo_mem_nhds (by linarith) (by linarith)
  obtain ⟨C, hC⟩ := (isCompact_Icc.prod isCompact_Icc :
      IsCompact (Icc (height-1) (height+1) ×ˢ Icc (0:ℝ) π)).exists_bound_of_continuousOn hprod.continuousOn
  have hInt : Integrable (F height) μ := by
    exact ((hF height).continuousOn.integrableOn_compact isCompact_Icc).mono_set Ioc_subset_Icc_self
  have hBound : ∀ᵐ t ∂μ, ∀ h ∈ s, ‖D h t‖ ≤ (fun _ : ℝ => C) t := by
    filter_upwards [ae_restrict_mem (μ := volume) (s := Ioc (0:ℝ) π) measurableSet_Ioc] with t ht
    intro h hh
    exact hC (h,t) ⟨⟨hh.1.le,hh.2.le⟩, Ioc_subset_Icc_self ht⟩
  have hDiff : ∀ᵐ t ∂μ, ∀ h ∈ s, HasDerivAt (F · t) (D h t) h := by
    apply Eventually.of_forall
    intro t h _
    simpa [F, D] using (((hψ (h+b t)).comp h ((hasDerivAt_id h).add_const (b t))).mul_const (J t))
  have H := (hasDerivAt_integral_of_dominated_loc_of_deriv_le hs
    (Eventually.of_forall fun h => (hF h).aestronglyMeasurable) hInt
    (hD height).aestronglyMeasurable hBound (integrable_const C) hDiff).2.const_mul M
  simpa [intervalIntegral.integral_of_le pi_pos.le, μ, F, D] using H

theorem circular_energy_force (k w a R L height : ℝ) (ha : 0 < a) :
    -deriv (fun h => 2*R*L * ∫ t in 0..π, potential k w a (circularGap h R t)) height =
    -2*R*L * ∫ t in 0..π, slope k w a (circularGap height R t) := by
  have H := translation_energy_derivative (potential k w a) (slope k w a)
    (fun t => R*cos t) (fun _ => 1) (fun z => potential_derivative k w a z ha)
    (slope_continuous k w a ha) (by fun_prop) continuous_const (2*R*L) height
  simp only [mul_one] at H
  unfold circularGap
  rw [H.deriv]
  ring

noncomputable def noncircularMetric (Rx Ry q t : ℝ) :=
  sqrt ((Rx*cos t)^2+(Ry*sin t+2*q*sin (2*t))^2)

theorem noncircular_energy_force (k w a Rx Ry q L height : ℝ) (ha : 0 < a) :
    -deriv (fun h => 2*L * ∫ t in 0..π,
        potential k w a (noncircularGap h Ry q t) * noncircularMetric Rx Ry q t) height =
    -2*L * ∫ t in 0..π, slope k w a (noncircularGap height Ry q t) * noncircularMetric Rx Ry q t := by
  have hJ : Continuous (noncircularMetric Rx Ry q) := by unfold noncircularMetric; fun_prop
  have H := translation_energy_derivative (potential k w a) (slope k w a)
    (fun t => Ry*cos t+q*cos (2*t)) (noncircularMetric Rx Ry q)
    (fun z => potential_derivative k w a z ha)
    (slope_continuous k w a ha) (by fun_prop) hJ (2*L) height
  unfold noncircularGap
  simp only [add_assoc]
  rw [H.deriv]
  ring


theorem inclined_gap_and_area (L β R h u t : ℝ) (hL : 0 < L) (hR : 0 < R) :
    let ell := sqrt (L^2+β^2)
    let xθ := R*(β/ell)*sin t
    let yθ := -R*(L/ell)*sin t
    let zθ := R*cos t
    -- y coordinate and norm of (L, β, 0) × (xθ, yθ, zθ).
    (h+β*(u-1/2)+R*(L/ell)*cos t = h+β*(u-1/2)+(R*L/ell)*cos t) ∧
    sqrt ((β*zθ)^2+(-L*zθ)^2+(L*yθ-β*xθ)^2) = R*ell := by
  dsimp
  have hp : 0 < L^2+β^2 := by positivity
  have he := sqrt_pos.2 hp
  have he2 := sq_sqrt hp.le
  constructor
  · ring
  · have hn : L * (-R * (L / sqrt (L ^ 2 + β ^ 2)) * sin t) -
          β * (R * (β / sqrt (L ^ 2 + β ^ 2)) * sin t) =
          -R*sqrt (L^2+β^2)*sin t := by
      field_simp
      linear_combination sin t * he2
    rw [hn, ← sqrt_sq (mul_nonneg hR.le he.le)]
    congr 1
    linear_combination
      (R^2*(sin t)^2-R^2) * he2 +
      (R^2*(L^2+β^2)) * sin_sq_add_cos_sq t

noncomputable def gapPrimitive (h R t : ℝ) := h*t+R*sin t
noncomputable def cubePrimitive (h R t : ℝ) :=
  h^3*t+3*h^2*R*sin t+3*h*R^2*(t/2+sin (2*t)/4)+R^3*(sin t-(sin t)^3/3)

theorem circular_gap_primitive (h R t : ℝ) :
    HasDerivAt (gapPrimitive h R) (circularGap h R t) t := by
  have H := ((hasDerivAt_id t).const_mul h).add ((hasDerivAt_sin t).const_mul R)
  convert H using 1 <;> (try rfl) <;> simp [circularGap]

theorem circular_gap_cube_primitive (h R t : ℝ) :
    HasDerivAt (cubePrimitive h R) ((circularGap h R t)^3) t := by
  have H := ((((hasDerivAt_id t).const_mul (h^3)).add
    ((hasDerivAt_sin t).const_mul (3*h^2*R))).add
    ((primitive_cos2 t).const_mul (3*h*R^2))).add ((primitive_cos3 t).const_mul (R^3))
  convert H using 1 <;> (try rfl) <;> simp [circularGap] <;> ring

theorem circular_segment_integral (h R A B lo hi : ℝ) :
    (∫ t in lo..hi, A*circularGap h R t+B*(circularGap h R t)^3) =
    (A*gapPrimitive h R hi+B*cubePrimitive h R hi) -
    (A*gapPrimitive h R lo+B*cubePrimitive h R lo) := by
  apply intervalIntegral.integral_eq_sub_of_hasDerivAt
  · intro x _
    exact ((circular_gap_primitive h R x).const_mul A).add
      ((circular_gap_cube_primitive h R x).const_mul B)
  · apply Continuous.intervalIntegrable
    unfold circularGap
    fun_prop

theorem endpoint_angle_quadratic_limit (R : ℝ) :
    Tendsto (fun t : ℝ => R*(1-cos t)/t^2) (𝓝[>] 0) (𝓝 (R/2)) := by
  have hs := (hasDerivAt_sin (0:ℝ)).tendsto_slope_zero_right
  have hsin : Tendsto (fun t : ℝ => sin t/t) (𝓝[>] 0) (𝓝 1) := by
    simpa [div_eq_mul_inv, mul_comm] using hs
  have hhalf : Tendsto (fun t : ℝ => t/2) (𝓝[>] 0) (𝓝[>] 0) := by
    apply tendsto_nhdsWithin_iff.2
    refine ⟨?_, ?_⟩
    · have hc : Continuous (fun t : ℝ => t/2) := by fun_prop
      simpa using (hc.continuousAt (x := (0:ℝ))).continuousWithinAt (s := Ioi (0:ℝ)) |>.tendsto
    · filter_upwards [self_mem_nhdsWithin] with t ht
      exact div_pos (mem_Ioi.mp ht) (by norm_num)
  have H := ((hsin.comp hhalf).pow 2).const_mul (R/2)
  simp only [one_pow, mul_one] at H
  apply H.congr'
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0)] with t ht
  have ht' : t ≠ 0 := ne_of_gt (mem_Ioi.mp ht)
  dsimp
  rw [endpoint_quadratic]
  field_simp
  <;> ring


theorem endpoint_linear_cosine_limit :
    Tendsto (fun t : ℝ => (sin t - t*cos t)/t^3) (𝓝[>] 0) (𝓝 (1/3)) := by
  have hs := (hasDerivAt_sin (0:ℝ)).tendsto_slope_zero_right
  have hsin : Tendsto (fun t : ℝ => sin t/t) (𝓝[>] 0) (𝓝 1) := by
    simpa [div_eq_mul_inv, mul_comm] using hs
  have hf (t : ℝ) : HasDerivAt (fun x => sin x-x*cos x) (t*sin t) t := by
    have H := (hasDerivAt_sin t).sub ((hasDerivAt_id t).mul (hasDerivAt_cos t))
    convert H using 1 <;> (try rfl) <;> simp [id] <;> ring
  refine HasDerivAt.lhopital_zero_right_on_Ioo (a := 0) (b := 1) (f' := fun t => t*sin t)
    (g' := fun t => 3*t^2) zero_lt_one (fun t _ => hf t) ?_ ?_ ?_ ?_ ?_
  · intro t _
    convert ((hasDerivAt_id t).pow 3) using 1 <;> (try rfl) <;> simp [id]
  · intro t ht
    have hpos : 0 < t := ht.1
    positivity
  · simpa using (hf 0).continuousAt.continuousWithinAt (s := Ioi (0:ℝ)) |>.tendsto
  · have hc : Continuous (fun t : ℝ => t^3) := by fun_prop
    simpa using (hc.continuousAt (x := (0:ℝ))).continuousWithinAt (s := Ioi (0:ℝ)) |>.tendsto
  · have H := hsin.div_const 3
    simpa only [one_div] using (H.congr' (by
      filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0)] with t ht
      have ht' : t ≠ 0 := ne_of_gt (mem_Ioi.mp ht)
      field_simp
      <;> ring))

theorem endpoint_linear_cosine_integral (R s : ℝ) :
    (∫ t in 0..s, R*(cos t - cos s)) = R*(sin s-s*cos s) := by
  have H (t : ℝ) : HasDerivAt (fun t => R*(sin t-t*cos s)) (R*(cos t-cos s)) t := by
    have HH := ((hasDerivAt_sin t).sub ((hasDerivAt_id t).mul_const (cos s))).const_mul R
    convert HH using 1 <;> (try rfl) <;> simp [id]
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => H x)
    (Continuous.intervalIntegrable (by fun_prop) 0 s)]
  simp


theorem endpoint_width_ratio (R : ℝ) (hR : 0 < R) :
    Tendsto (fun t : ℝ => t / sqrt (R*(1-cos t))) (𝓝[>] 0) (𝓝 (1/sqrt (R/2))) := by
  have H := (tendsto_const_nhds (x := (1:ℝ)) (f := 𝓝[>] (0:ℝ))).div
    (endpoint_angle_quadratic_limit R).sqrt (ne_of_gt (sqrt_pos.2 (by positivity)))
  apply H.congr'
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0)] with t ht
  have hp : 0 < t := mem_Ioi.mp ht
  have hd : 0 ≤ R*(1-cos t) := mul_nonneg hR.le (sub_nonneg.mpr (cos_le_one t))
  change 1 / sqrt (R*(1-cos t)/t^2) = t / sqrt (R*(1-cos t))
  rw [sqrt_div hd, sqrt_sq hp.le, one_div, inv_div]

theorem endpoint_scaled_linear_limit (R L jump : ℝ) (hR : 0 < R) :
    Tendsto (fun t : ℝ => (2*R*L*jump * (R*(sin t-t*cos t))) /
       (sqrt (R*(1-cos t)))^3) (𝓝[>] 0)
       (𝓝 ((2*R*L*jump*R/3)/(sqrt (R/2))^3)) := by
  have Hd := ((endpoint_angle_quadratic_limit R).sqrt).pow 3
  have Hd0 : sqrt (R/2)^3 ≠ 0 := ne_of_gt (pow_pos (sqrt_pos.2 (by positivity)) 3)
  have H := (endpoint_linear_cosine_limit.const_mul (2*R*L*jump*R)).div Hd Hd0
  have hval : 2*R*L*jump*R*(1/3) / sqrt (R/2)^3 =
      (2*R*L*jump*R/3) / sqrt (R/2)^3 := by ring
  rw [hval] at H
  apply H.congr'
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0)] with t ht
  have hp : 0 < t := mem_Ioi.mp ht
  have hd : 0 ≤ R*(1-cos t) := mul_nonneg hR.le (sub_nonneg.mpr (cos_le_one t))
  change (2*R*L*jump*R*((sin t-t*cos t)/t^3)) / sqrt (R*(1-cos t)/t^2)^3 =
    (2*R*L*jump*(R*(sin t-t*cos t))) / sqrt (R*(1-cos t))^3
  rw [sqrt_div hd, sqrt_sq hp.le]
  by_cases hz : sqrt (R*(1-cos t)) = 0
  · simp [hz]
  · field_simp
    <;> ring

theorem compression_adhesion_mismatch (k w a v : ℝ) :
    k*v - adhesionSlope w a v = (k-4*w/a^2)*v + (4*w/a^4)*v^3 := by
  unfold adhesionSlope
  by_cases ha : a = 0
  · simp [ha]
  · field_simp
    <;> ring

theorem adhesion_separation_mismatch (w a v : ℝ) (ha : a ≠ 0) :
    adhesionSlope w a (a+v) =
      (-8*w/a^2)*v - (12*w/a^3)*v^2 - (4*w/a^4)*v^3 := by
  unfold adhesionSlope
  field_simp
  <;> ring

theorem circular_patch_remainder (R s A B : ℝ) (hR : 0 < R)
    (hs : s ∈ Icc 0 π) (hd : R*(1-cos s) ≤ 1) :
    let δ := R*(1-cos s)
    let v := fun t : ℝ => R*(cos t-cos s)
    |∫ t in 0..s, A*(v t)^2 + B*(v t)^3| ≤ ((|A|+|B|)*δ^2)*s := by
  dsimp
  rw [← Real.norm_eq_abs]
  refine (intervalIntegral.norm_integral_le_of_norm_le_const
    (C := (|A|+|B|)*(R*(1-cos s))^2) ?_).trans_eq ?_
  swap
  · simp [abs_of_nonneg hs.1]
  intro t ht
  rw [uIoc_of_le hs.1] at ht
  have hcos := strictAntiOn_cos.antitoneOn (show t ∈ Icc 0 π from ⟨ht.1.le, ht.2.trans hs.2⟩) hs ht.2
  have e0 : 0 ≤ R*(cos t-cos s) := mul_nonneg hR.le (sub_nonneg.mpr hcos)
  have e1 : R*(cos t-cos s) ≤ R*(1-cos s) :=
    mul_le_mul_of_nonneg_left (sub_le_sub_right (cos_le_one t) _) hR.le
  have ep2 := pow_le_pow_left₀ e0 e1 2
  have ep3 : (R*(cos t-cos s))^3 ≤ (R*(1-cos s))^2 := by
    calc
      (R*(cos t-cos s))^3 = (R*(cos t-cos s))^2 * (R*(cos t-cos s)) := by ring
      _ ≤ (R*(cos t-cos s))^2 * 1 := mul_le_mul_of_nonneg_left (e1.trans hd) (sq_nonneg _)
      _ ≤ _ := by simpa using ep2
  rw [Real.norm_eq_abs]
  calc
    _ ≤ |A*(R*(cos t-cos s))^2| + |B*(R*(cos t-cos s))^3| := abs_add_le _ _
    _ = |A| * (R*(cos t-cos s))^2 + |B| * (R*(cos t-cos s))^3 := by
      simp only [abs_mul A, abs_mul B, abs_of_nonneg (sq_nonneg (R * (cos t - cos s))), abs_of_nonneg (pow_nonneg e0 3)]
    _ ≤ |A| * (R*(1-cos s))^2 + |B| * (R*(1-cos s))^2 :=
      add_le_add (mul_le_mul_of_nonneg_left ep2 (abs_nonneg _))
        (mul_le_mul_of_nonneg_left ep3 (abs_nonneg _))
    _ = _ := by ring


theorem circular_patch_scaled_remainder (R A B : ℝ) (hR : 0 < R) :
    Tendsto (fun s : ℝ =>
      (∫ t in 0..s, A*(R*(cos t-cos s))^2+B*(R*(cos t-cos s))^3) /
        (sqrt (R*(1-cos s)))^3) (𝓝[>] 0) (𝓝 0) := by
  let gap := fun s : ℝ => R*(1-cos s)
  have hc : Continuous gap := by dsimp [gap]; fun_prop
  have hg : Tendsto gap (𝓝[>] 0) (𝓝 0) := by
    simpa [gap] using (hc.continuousAt (x := (0:ℝ))).continuousWithinAt (s := Ioi (0:ℝ)) |>.tendsto
  have HH := ((hg.mul (endpoint_width_ratio R hR)).const_mul (|A|+|B|))
  simp only [zero_mul, mul_zero] at HH
  apply squeeze_zero_norm' _ HH
  have hSmall : ∀ᶠ s in 𝓝[>] (0:ℝ), gap s < 1 := (tendsto_order.1 hg).2 _ zero_lt_one
  have hsπ : ∀ᶠ s in 𝓝[>] (0:ℝ), s < π :=
    eventually_nhdsWithin_of_eventually_nhds (gt_mem_nhds pi_pos)
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0), hSmall, hsπ] with s hs hsmall hsπ
  have hs0 : 0 < s := mem_Ioi.mp hs
  have hmem : s ∈ Icc 0 π := ⟨hs0.le, hsπ.le⟩
  have hcos : cos s < 1 := by
    simpa using (strictAntiOn_cos (show (0:ℝ) ∈ Icc 0 π from ⟨le_rfl,pi_pos.le⟩) hmem hs0)
  have hd : 0 < gap s := mul_pos hR (sub_pos.mpr hcos)
  have hden := pow_nonneg (sqrt_nonneg (gap s)) 3
  have H := circular_patch_remainder R s A B hR hmem hsmall.le
  change |∫ t in 0..s, A*(R*(cos t-cos s))^2+B*(R*(cos t-cos s))^3| ≤
    (|A|+|B|)*(gap s)^2*s at H
  rw [Real.norm_eq_abs, abs_div, abs_of_nonneg hden]
  refine (div_le_div_of_nonneg_right H hden).trans_eq ?_
  have hsqrt := sq_sqrt hd.le
  have hz : sqrt (gap s) ≠ 0 := (sqrt_pos.2 hd).ne'
  change ((|A|+|B|)*(gap s)^2*s)/(sqrt (gap s))^3 =
    (|A|+|B|)*(gap s*(s/sqrt (gap s)))
  have heq : (gap s)^2 = (sqrt (gap s))^2 * gap s := by rw [hsqrt]; ring
  rw [heq]
  field_simp
  <;> ring

theorem circular_patch_exact_leading (R L jump A B : ℝ) (hR : 0 < R) :
    Tendsto (fun s : ℝ => 2*R*L *
      (∫ t in 0..s, jump*(R*(cos t-cos s)) +
         (A*(R*(cos t-cos s))^2+B*(R*(cos t-cos s))^3)) /
         (sqrt (R*(1-cos s)))^3) (𝓝[>] 0)
      (𝓝 ((2*R*L*jump*R/3)/(sqrt (R/2))^3)) := by
  have Hr := (circular_patch_scaled_remainder R A B hR).const_mul (2*R*L)
  simp only [mul_zero] at Hr
  have H := (endpoint_scaled_linear_limit R L jump hR).add Hr
  simp only [add_zero] at H
  apply H.congr'
  apply Eventually.of_forall
  intro s
  dsimp
  let f := fun t : ℝ => jump*(R*(cos t-cos s))
  let g := fun t : ℝ => A*(R*(cos t-cos s))^2+B*(R*(cos t-cos s))^3
  have hf : Continuous f := by dsimp [f]; fun_prop
  have hg : Continuous g := by dsimp [g]; fun_prop
  have hadd := intervalIntegral.integral_add (hf.intervalIntegrable (μ := volume) 0 s) (hg.intervalIntegrable (μ := volume) 0 s)
  change _ = 2*R*L*(∫ t in 0..s, f t + g t)/sqrt (R*(1-cos s))^3
  rw [hadd]
  unfold f
  rw [intervalIntegral.integral_const_mul, endpoint_linear_cosine_integral]
  unfold g
  ring

theorem circular_endpoint_reflection (h R t : ℝ) :
    circularGap h R (π-t) = h-R*cos t := by
  rw [circularGap, cos_pi_sub]
  ring

theorem noncircular_reduces_to_circle (h R t : ℝ) (hR : 0 ≤ R) :
    noncircularGap h R 0 t = circularGap h R t ∧ noncircularMetric R R 0 t = R := by
  constructor
  · simp [noncircularGap, circularGap]
  · unfold noncircularMetric
    have H : (R*cos t)^2+(R*sin t+2*0*sin (2*t))^2 = R^2 := by
      linear_combination (R^2) * sin_sq_add_cos_sq t
    rw [H, sqrt_sq hR]


theorem endpoint_width_from_gap (R δ : ℝ) (hR : 0 < R)
    (hd : 0 < δ) (hu : δ < 2*R) :
    let s := arccos (1-δ/R)
    s ∈ Ioo 0 π ∧ R*(1-cos s) = δ := by
  have h1 : 0 < δ/R := div_pos hd hR
  have h2 : δ/R < 2 := (div_lt_iff₀ hR).2 (by linarith)
  have hl : -1 < 1-δ/R := by linarith
  have hr : 1-δ/R < 1 := by linarith
  dsimp
  constructor
  · exact ⟨arccos_pos.mpr hr, arccos_lt_pi.mpr hl⟩
  · rw [cos_arccos hl.le hr.le]
    field_simp
    <;> ring

theorem endpoint_width_vanishes (R : ℝ) (hR : 0 < R) :
    Tendsto (fun δ : ℝ => arccos (1-δ/R)) (𝓝[>] 0) (𝓝[>] 0) := by
  apply tendsto_nhdsWithin_iff.2
  constructor
  · have hc : Continuous (fun δ : ℝ => arccos (1-δ/R)) := by fun_prop
    simpa using (hc.continuousAt (x := (0:ℝ))).continuousWithinAt (s := Ioi (0:ℝ)) |>.tendsto
  · filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0)] with δ hd
    exact arccos_pos.mpr (sub_lt_self _ (div_pos (mem_Ioi.mp hd) hR))

theorem circular_patch_delta_limit (R L jump A B : ℝ) (hR : 0 < R) :
    Tendsto (fun δ : ℝ => 2*R*L *
      (∫ t in 0..arccos (1-δ/R), jump*(R*(cos t-cos (arccos (1-δ/R)))) +
         (A*(R*(cos t-cos (arccos (1-δ/R))))^2+B*(R*(cos t-cos (arccos (1-δ/R))))^3)) /
         (sqrt δ)^3) (𝓝[>] 0)
      (𝓝 ((2*R*L*jump*R/3)/(sqrt (R/2))^3)) := by
  have H := (circular_patch_exact_leading R L jump A B hR).comp (endpoint_width_vanishes R hR)
  apply H.congr'
  have hsmall : ∀ᶠ δ : ℝ in 𝓝[>] 0, δ < 2*R :=
    eventually_nhdsWithin_of_eventually_nhds (gt_mem_nhds (by positivity))
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0), hsmall] with δ hd hu
  have he := (endpoint_width_from_gap R δ hR (mem_Ioi.mp hd) hu).2
  dsimp
  rw [he]


theorem potential_C1 (k w a : ℝ) (ha : 0 < a) : ContDiff ℝ 1 (potential k w a) := by
  rw [contDiff_one_iff_deriv]
  refine ⟨fun x => (potential_derivative k w a x ha).differentiableAt, ?_⟩
  have hd : deriv (potential k w a) = slope k w a :=
    funext fun x => (potential_derivative k w a x ha).deriv
  rw [hd]
  exact slope_continuous k w a ha

theorem circular_closed_root_range (h R B : ℝ) (hR : 0 < R) :
    (∃ t ∈ Icc 0 π, circularGap h R t = B) ↔ h-R ≤ B ∧ B ≤ h+R := by
  constructor
  · rintro ⟨t, _, ht⟩
    unfold circularGap at ht
    constructor
    · linarith [mul_le_mul_of_nonneg_left (neg_one_le_cos t) hR.le]
    · linarith [mul_le_mul_of_nonneg_left (cos_le_one t) hR.le]
  · rintro ⟨hl, hu⟩
    have hlo : -1 ≤ (B-h)/R := (le_div_iff₀ hR).2 (by linarith)
    have hhi : (B-h)/R ≤ 1 := (div_le_iff₀ hR).2 (by linarith)
    refine ⟨arccos ((B-h)/R), ⟨arccos_nonneg _, arccos_le_pi _⟩, ?_⟩
    unfold circularGap
    rw [cos_arccos hlo hhi]
    field_simp
    <;> ring

theorem circular_endpoint_cases (h R t : ℝ) (hR : 0 < R) (ht : t ∈ Icc 0 π) :
    (circularGap h R t = h+R ↔ t = 0) ∧
    (circularGap h R t = h-R ↔ t = π) := by
  have anti := circular_unique h R hR
  have h0 : (0:ℝ) ∈ Icc 0 π := ⟨le_rfl, pi_pos.le⟩
  have hp : π ∈ Icc 0 π := ⟨pi_pos.le, le_rfl⟩
  have e0 : circularGap h R 0 = h+R := by simp [circularGap]
  have ep : circularGap h R π = h-R := by simp [circularGap, sub_eq_add_neg]
  rw [← e0, ← ep]
  exact ⟨⟨anti.injOn ht h0, congrArg _⟩, ⟨anti.injOn ht hp, congrArg _⟩⟩

theorem circular_prescribed_heights (B R e : ℝ) :
    circularGap (B-R*cos e) R e = B ∧ circularGap (B+R*cos e) R (π-e) = B := by
  simp [circularGap, cos_pi_sub]

theorem compression_closed_branch (k w a g : ℝ) (hg : g ≤ 0) :
    slope k w a g = k*g := if_pos hg

theorem adhesion_closed_branch (k w a g : ℝ) (ha : 0 < a) (h0 : 0 ≤ g) (h1 : g ≤ a) :
    slope k w a g = (4*w/a^2)*g-(4*w/a^4)*g^3 := by
  rcases eq_or_lt_of_le h0 with he | hp
  · subst g
    simp [slope]
  · rcases eq_or_lt_of_le h1 with he | hlt
    · subst g
      simp [slope, not_le.mpr ha]
      field_simp
      <;> ring
    · simp only [slope, if_neg (not_le.mpr hp), if_pos hlt]
      unfold adhesionSlope
      field_simp
      <;> ring

theorem separated_closed_branch (k w a g : ℝ) (ha : 0 < a) (hg : a ≤ g) :
    slope k w a g = 0 := by simp [slope, not_le.mpr (ha.trans_le hg), not_lt.mpr hg]

theorem actual_circular_segment (k w a h R A B lo hi : ℝ)
    (hbranch : ∀ t ∈ uIcc lo hi, slope k w a (circularGap h R t) =
      A*circularGap h R t + B*(circularGap h R t)^3) :
    (∫ t in lo..hi, slope k w a (circularGap h R t)) =
    (A*gapPrimitive h R hi+B*cubePrimitive h R hi) -
    (A*gapPrimitive h R lo+B*cubePrimitive h R lo) := by
  rw [intervalIntegral.integral_congr hbranch, circular_segment_integral]

theorem circular_segment_assembly (k w a h R : ℝ) (ha : 0 < a) (ends : ℕ → ℝ) (n : ℕ) :
    (∑ i ∈ Finset.Ico 0 n, ∫ t in ends i..ends (i+1), slope k w a (circularGap h R t)) =
    (∫ t in ends 0..ends n, slope k w a (circularGap h R t)) := by
  have hc : Continuous (fun t => slope k w a (circularGap h R t)) :=
    (slope_continuous k w a ha).comp (by unfold circularGap; fun_prop)
  exact intervalIntegral.sum_integral_adjacent_intervals_Ico (Nat.zero_le n)
    (fun i _ => hc.intervalIntegrable _ _)

theorem noncircular_numeric_monotonicity :
    (0:ℝ) ≤ 45/1000 ∧ (4:ℝ)*(45/1000) < 28209479177387814/100000000000000000 := by
  norm_num

theorem inclined_axial_crossings (h β r B u : ℝ) (hb : β ≠ 0) :
    (h+β*(u-1/2)+r = B ↔ u = 1/2+(B-h-r)/β) ∧
    (h+β*(u-1/2)-r = B ↔ u = 1/2+(B-h+r)/β) := by
  constructor <;> constructor <;> intro H
  · field_simp
    linarith
  · rw [H]
    field_simp
    <;> ring
  · field_simp
    linarith
  · rw [H]
    field_simp
    <;> ring

theorem force_absolute_error_scaling (m estimate reference : ℝ) (hm : 0 ≤ m) :
    |(-m*estimate)-(-m*reference)| = m*|estimate-reference| := by
  rw [show -m*estimate - -m*reference = -m*(estimate-reference) by ring,
      abs_mul, abs_neg, abs_of_nonneg hm]

theorem endpoint_width_delta_limit (R : ℝ) (hR : 0 < R) :
    Tendsto (fun δ : ℝ => arccos (1-δ/R) / sqrt δ) (𝓝[>] 0) (𝓝 (1/sqrt (R/2))) := by
  have H := (endpoint_width_ratio R hR).comp (endpoint_width_vanishes R hR)
  apply H.congr'
  have hsmall : ∀ᶠ δ : ℝ in 𝓝[>] 0, δ < 2*R :=
    eventually_nhdsWithin_of_eventually_nhds (gt_mem_nhds (by positivity))
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0), hsmall] with δ hd hu
  have he := (endpoint_width_from_gap R δ hR (mem_Ioi.mp hd) hu).2
  dsimp
  rw [he]

theorem endpoint_limit_nonzero (R L m : ℝ) (hR : 0 < R) (hL : 0 < L) (hm : m ≠ 0) :
    (2*R*L*m*R/3)/(sqrt (R/2))^3 ≠ 0 := by
  have hd : sqrt (R/2) ≠ 0 := (sqrt_pos.2 (by positivity)).ne'
  exact div_ne_zero (div_ne_zero (by apply mul_ne_zero <;> simp [hR.ne',hL.ne',hm]) (by norm_num))
    (pow_ne_zero _ hd)


theorem actual_endpoint_branch_mismatches (k w a v : ℝ) (ha : 0 < a)
    (hv : 0 ≤ v) (hva : v ≤ a) :
    (slope k w a (-v) - adhesionSlope w a (-v) =
       (k-4*w/a^2)*(-v)+(4*w/a^4)*(-v)^3) ∧
    (slope k w a v - k*v =
       -(k-4*w/a^2)*v-(4*w/a^4)*v^3) ∧
    (slope k w a (a-v) =
       (-8*w/a^2)*(-v)-(12*w/a^3)*(-v)^2-(4*w/a^4)*(-v)^3) ∧
    (slope k w a (a+v) - adhesionSlope w a (a+v) =
       (8*w/a^2)*v+(12*w/a^3)*v^2+(4*w/a^4)*v^3) := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [compression_closed_branch k w a (-v) (by linarith)]
    exact compression_adhesion_mismatch k w a (-v)
  · rw [adhesion_closed_branch k w a v ha hv hva]
    have H := compression_adhesion_mismatch k w a v
    unfold adhesionSlope at H
    linarith [show (4*w/a^2)*v-(4*w/a^4)*v^3 =
      4*w*v/a^2*(1-(v/a)^2) by field_simp <;> ring]
  · have hv0 : 0 ≤ a-v := by linarith
    have hval : slope k w a (a-v) = adhesionSlope w a (a-v) := by
      rw [adhesion_closed_branch k w a (a-v) ha hv0 (by linarith)]
      unfold adhesionSlope
      field_simp
      <;> ring
    rw [hval, sub_eq_add_neg]
    exact adhesion_separation_mismatch w a (-v) ha.ne'
  · rw [separated_closed_branch k w a (a+v) ha (by linarith),
        adhesion_separation_mismatch w a v ha.ne']
    ring

theorem circular_crossing_chain_factor (h R t : ℝ) (hR : 0 < R) (ht : t ∈ Ioo 0 π) :
    HasDerivAt (circularGap h R) (-R*sin t) t ∧ -R*sin t ≠ 0 := by
  constructor
  · have H := ((hasDerivAt_cos t).const_mul R).const_add h
    convert H using 1 <;> (try rfl) <;> ring
  · exact mul_ne_zero (neg_ne_zero.mpr hR.ne') (sin_pos_of_pos_of_lt_pi ht.1 ht.2).ne'

theorem circular_angular_slope_jump
    (F left right : ℝ → ℝ) (h R tB mleft mright : ℝ)
    (hR : 0 < R) (ht : tB ∈ Ioo 0 π)
    (hl : HasDerivAt left mleft (circularGap h R tB))
    (hr : HasDerivAt right mright (circularGap h R tB))
    (hL : F =ᶠ[𝓝[Iic tB] tB] (fun t => left (circularGap h R t)))
    (hU : F =ᶠ[𝓝[Ici tB] tB] (fun t => right (circularGap h R t)))
    (hjoinL : F tB = left (circularGap h R tB))
    (hjoinU : F tB = right (circularGap h R tB))
    (hm : mleft ≠ mright) :
    HasDerivWithinAt F (mleft*(-R*sin tB)) (Iic tB) tB ∧
    HasDerivWithinAt F (mright*(-R*sin tB)) (Ici tB) tB ∧
    mleft*(-R*sin tB) ≠ mright*(-R*sin tB) := by
  obtain ⟨hd, hn⟩ := circular_crossing_chain_factor h R tB hR ht
  exact ⟨(hl.comp tB hd).hasDerivWithinAt.congr_of_eventuallyEq hL hjoinL,
    (hr.comp tB hd).hasDerivWithinAt.congr_of_eventuallyEq hU hjoinU,
    fun he => hm (mul_right_cancel₀ hn he)⟩


theorem circular_closed_form_actual (k w a h R : ℝ) (ha : 0 < a) (hR : 0 < R) :
    let θa := arccos ((a-h)/R)
    let θ0 := arccos ((0-h)/R)
    let P := fun t => (4*w/a^2)*gapPrimitive h R t-(4*w/a^4)*cubePrimitive h R t
    (∫ t in 0..π, slope k w a (circularGap h R t)) =
      (P θ0-P θa) + k*(gapPrimitive h R π-gapPrimitive h R θ0) := by
  dsimp
  let θa := arccos ((a-h)/R)
  let θ0 := arccos ((0-h)/R)
  have hA : 0 ≤ θa := arccos_nonneg _
  have hO : θa ≤ θ0 := antitone_arccos (div_le_div_of_nonneg_right (by linarith) hR.le)
  have hZ : θ0 ≤ π := arccos_le_pi _
  have hc : Continuous (fun t => slope k w a (circularGap h R t)) :=
    (slope_continuous k w a ha).comp (by unfold circularGap; fun_prop)
  have before (B t : ℝ) (ht : t ∈ Icc 0 π) (hle : t < arccos ((B-h)/R)) :
      B < circularGap h R t := by
    by_contra! hn
    have hd : cos t ≤ (B-h)/R := (le_div_iff₀ hR).2 (by unfold circularGap at hn; nlinarith)
    have H := antitone_arccos hd
    rw [arccos_cos ht.1 ht.2] at H
    linarith
  have after (B t : ℝ) (ht : t ∈ Icc 0 π) (hle : arccos ((B-h)/R) < t) :
      circularGap h R t < B := by
    by_contra! hn
    have hd : (B-h)/R ≤ cos t := (div_le_iff₀ hR).2 (by unfold circularGap at hn; nlinarith)
    have H := antitone_arccos hd
    rw [arccos_cos ht.1 ht.2] at H
    linarith
  have I1 : (∫ t in 0..θa, slope k w a (circularGap h R t)) = 0 := calc
    _ = ∫ _t in 0..θa, (0:ℝ) := by
      apply intervalIntegral.integral_congr_Ioo_of_le hA
      intro t ht
      exact separated_closed_branch k w a _ ha
        (before a t ⟨ht.1.le, ht.2.le.trans (hO.trans hZ)⟩ ht.2).le
    _ = 0 := by simp
  have I2 : (∫ t in θa..θ0, slope k w a (circularGap h R t)) =
      ((4*w/a^2)*gapPrimitive h R θ0-(4*w/a^4)*cubePrimitive h R θ0) -
      ((4*w/a^2)*gapPrimitive h R θa-(4*w/a^4)*cubePrimitive h R θa) := calc
    _ = ∫ t in θa..θ0, (4*w/a^2)*circularGap h R t +
       (-(4*w/a^4))*(circularGap h R t)^3 := by
      apply intervalIntegral.integral_congr_Ioo_of_le hO
      intro t ht
      have hi : t ∈ Icc 0 π := ⟨hA.trans ht.1.le,ht.2.le.trans hZ⟩
      simpa only [sub_eq_add_neg, neg_mul] using adhesion_closed_branch k w a _ ha
        (before 0 t hi ht.2).le (after a t hi ht.1).le
    _ = _ := by rw [circular_segment_integral]; ring
  have I3 : (∫ t in θ0..π, slope k w a (circularGap h R t)) =
      k*(gapPrimitive h R π-gapPrimitive h R θ0) := calc
    _ = ∫ t in θ0..π, k*circularGap h R t + 0*(circularGap h R t)^3 := by
      apply intervalIntegral.integral_congr_Ioo_of_le hZ
      intro t ht
      have hi : t ∈ Icc 0 π := ⟨hA.trans (hO.trans ht.1.le),ht.2.le⟩
      simpa using compression_closed_branch k w a _ (after 0 t hi ht.1).le
    _ = _ := by rw [circular_segment_integral]; ring
  have H1 := intervalIntegral.integral_add_adjacent_intervals
    (hc.intervalIntegrable (μ := volume) 0 θa) (hc.intervalIntegrable θa θ0)
  have H2 := intervalIntegral.integral_add_adjacent_intervals
    (hc.intervalIntegrable (μ := volume) 0 θ0) (hc.intervalIntegrable θ0 π)
  rw [I1,I2,zero_add] at H1
  rw [← H2, ← H1, I3]

theorem noncircular_metric_is_speed (Rx Ry q h t : ℝ) :
    noncircularMetric Rx Ry q t =
      sqrt ((deriv (fun u : ℝ => Rx*sin u) t)^2 +
        (deriv (noncircularGap h Ry q) t)^2) :=
  (noncircular_arc_factor Rx Ry q h t).symm


private theorem circular_actual_jump_bridge
    (k w a h R tB B ε mhi mlo : ℝ) (high low : ℝ → ℝ)
    (hR : 0 < R) (ht : tB ∈ Ioo 0 π) (hB : circularGap h R tB = B)
    (hε : 0 < ε)
    (hhi : HasDerivAt high mhi B) (hlo : HasDerivAt low mlo B)
    (up : ∀ g, B ≤ g → g < B+ε → slope k w a g = high g)
    (down : ∀ g, g ≤ B → B-ε < g → slope k w a g = low g)
    (hm : mhi ≠ mlo) :
    let F := fun t => slope k w a (circularGap h R t)
    HasDerivWithinAt F (mhi*(-R*sin tB)) (Iic tB) tB ∧
    HasDerivWithinAt F (mlo*(-R*sin tB)) (Ici tB) tB ∧
    mhi*(-R*sin tB) ≠ mlo*(-R*sin tB) := by
  dsimp
  have hc : Continuous (circularGap h R) := by unfold circularGap; fun_prop
  have hnear : ∀ᶠ t in 𝓝 tB, circularGap h R t ∈ Ioo (B-ε) (B+ε) :=
    hc.continuousAt.eventually (show Ioo (B-ε) (B+ε) ∈ 𝓝 (circularGap h R tB) by
      rw [hB]; exact Ioo_mem_nhds (by linarith) (by linarith))
  have hdomain : ∀ᶠ t in 𝓝 tB, t ∈ Ioo (0:ℝ) π := Ioo_mem_nhds ht.1 ht.2
  have anti := (circular_unique h R hR).antitoneOn
  have hbcc : tB ∈ Icc 0 π := ⟨ht.1.le,ht.2.le⟩
  apply circular_angular_slope_jump _ high low h R tB mhi mlo hR ht
  · simpa [hB] using hhi
  · simpa [hB] using hlo
  · filter_upwards [self_mem_nhdsWithin, eventually_nhdsWithin_of_eventually_nhds hnear,
      eventually_nhdsWithin_of_eventually_nhds hdomain] with t htL hn hd
    exact up _ (hB ▸ anti ⟨hd.1.le,hd.2.le⟩ hbcc (mem_Iic.mp htL)) hn.2
  · filter_upwards [self_mem_nhdsWithin, eventually_nhdsWithin_of_eventually_nhds hnear,
      eventually_nhdsWithin_of_eventually_nhds hdomain] with t htU hn hd
    exact down _ (hB ▸ anti hbcc ⟨hd.1.le,hd.2.le⟩ (mem_Ici.mp htU)) hn.1
  · rw [hB]
    exact up B le_rfl (by linarith)
  · rw [hB]
    exact down B le_rfl (by linarith)
  · exact hm

private theorem adhesive_slope_equals (k w a g : ℝ) (ha : 0 < a)
    (hg : g ∈ Icc 0 a) : slope k w a g = adhesionSlope w a g := by
  rw [adhesion_closed_branch k w a g ha hg.1 hg.2]
  unfold adhesionSlope
  field_simp
  <;> ring

theorem actual_circular_jump_at_zero (k w a h R tB : ℝ)
    (ha : 0 < a) (hR : 0 < R) (ht : tB ∈ Ioo 0 π)
    (hB : circularGap h R tB = 0) (hm : 4*w/a^2 ≠ k) :
    let F := fun t => slope k w a (circularGap h R t)
    HasDerivWithinAt F ((4*w/a^2)*(-R*sin tB)) (Iic tB) tB ∧
    HasDerivWithinAt F (k*(-R*sin tB)) (Ici tB) tB ∧
    (4*w/a^2)*(-R*sin tB) ≠ k*(-R*sin tB) := by
  apply circular_actual_jump_bridge k w a h R tB 0 a _ _ (adhesionSlope w a)
    (fun x => k*x) hR ht hB ha (branch_curvature_at_left w a ha.ne')
  · simpa using (hasDerivAt_id (0:ℝ)).const_mul k
  · intro g h0 h1
    exact adhesive_slope_equals k w a g ha ⟨h0,by linarith⟩
  · intro g hg _
    exact compression_closed_branch k w a g hg
  · exact hm

theorem actual_circular_jump_at_range (k w a h R tB : ℝ)
    (ha : 0 < a) (hR : 0 < R) (ht : tB ∈ Ioo 0 π)
    (hB : circularGap h R tB = a) (hm : (0:ℝ) ≠ -8*w/a^2) :
    let F := fun t => slope k w a (circularGap h R t)
    HasDerivWithinAt F (0*(-R*sin tB)) (Iic tB) tB ∧
    HasDerivWithinAt F ((-8*w/a^2)*(-R*sin tB)) (Ici tB) tB ∧
    0*(-R*sin tB) ≠ (-8*w/a^2)*(-R*sin tB) := by
  apply circular_actual_jump_bridge k w a h R tB a a _ _ (fun _ => 0)
    (adhesionSlope w a) hR ht hB ha (hasDerivAt_const a (0:ℝ))
    (branch_curvature_at_right w a ha.ne')
  · intro g hg _
    exact separated_closed_branch k w a g ha hg
  · intro g h1 h0
    exact adhesive_slope_equals k w a g ha ⟨by linarith,h1⟩
  · exact hm


noncomputable def patchMismatch (k w a : ℝ) (atRange positiveGap : Bool) (v : ℝ) : ℝ :=
  if atRange then
    if positiveGap then slope k w a (a+v) - adhesionSlope w a (a+v)
    else slope k w a (a-v)
  else
    if positiveGap then slope k w a v - k*v
    else slope k w a (-v) - adhesionSlope w a (-v)

noncomputable def patchCoefficients (k w a : ℝ) (atRange positiveGap : Bool) : ℝ × ℝ × ℝ :=
  if atRange then
    if positiveGap then (8*w/a^2, 12*w/a^3, 4*w/a^4)
    else (8*w/a^2, -12*w/a^3, 4*w/a^4)
  else (-(k-4*w/a^2), 0, -(4*w/a^4))

theorem actual_circular_patch_delta_limit (k w a R L : ℝ)
    (atRange positiveGap : Bool) (ha : 0 < a) (hR : 0 < R) :
    let C := patchCoefficients k w a atRange positiveGap
    Tendsto (fun δ : ℝ => 2*R*L *
      (∫ t in 0..arccos (1-δ/R),
        patchMismatch k w a atRange positiveGap (R*(cos t-cos (arccos (1-δ/R))))) /
         (sqrt δ)^3) (𝓝[>] 0)
      (𝓝 ((2*R*L*C.1*R/3)/(sqrt (R/2))^3)) := by
  let C := patchCoefficients k w a atRange positiveGap
  dsimp only
  apply (circular_patch_delta_limit R L C.1 C.2.1 C.2.2 hR).congr'
  have ha' : ∀ᶠ δ : ℝ in 𝓝[>] 0, δ < a :=
    eventually_nhdsWithin_of_eventually_nhds (gt_mem_nhds ha)
  have hR' : ∀ᶠ δ : ℝ in 𝓝[>] 0, δ < 2*R :=
    eventually_nhdsWithin_of_eventually_nhds (gt_mem_nhds (by positivity))
  filter_upwards [self_mem_nhdsWithin (a := (0:ℝ)) (s := Ioi 0),ha',hR'] with δ hd hda hdR
  obtain ⟨hs, hδ⟩ := endpoint_width_from_gap R δ hR (mem_Ioi.mp hd) hdR
  have he := intervalIntegral.integral_congr_Ioo_of_le
    (a := 0) (b := arccos (1-δ/R)) (μ := volume) hs.1.le
    (f := fun t => C.1*(R*(cos t-cos (arccos (1-δ/R)))) +
      (C.2.1*(R*(cos t-cos (arccos (1-δ/R))))^2+C.2.2*(R*(cos t-cos (arccos (1-δ/R))))^3))
    (g := fun t => patchMismatch k w a atRange positiveGap (R*(cos t-cos (arccos (1-δ/R)))))
    ?_
  · rw [he]
  · intro t ht
    let s := arccos (1-δ/R)
    let v := R*(cos t-cos s)
    have hv0 : 0 ≤ v := by
      exact mul_nonneg hR.le (sub_nonneg.mpr (strictAntiOn_cos.antitoneOn
        ⟨ht.1.le,ht.2.le.trans hs.2.le⟩ ⟨hs.1.le,hs.2.le⟩ ht.2.le))
    have hva : v ≤ a := by
      have hle : v ≤ R*(1-cos s) :=
        mul_le_mul_of_nonneg_left (sub_le_sub_right (cos_le_one t) _) hR.le
      exact hle.trans (hδ ▸ hda.le)
    have whole := actual_endpoint_branch_mismatches k w a v ha hv0 hva
    change C.1*v+(C.2.1*v^2+C.2.2*v^3) = patchMismatch k w a atRange positiveGap v
    cases atRange <;> cases positiveGap <;> dsimp [C,patchCoefficients,patchMismatch]
    · rw [whole.1]
      ring
    · rw [whole.2.1]
      ring
    · rw [whole.2.2.1]
      ring
    · rw [whole.2.2.2]
      ring

theorem noncircular_bracketed_root (h R q B : ℝ) (hq : 0 ≤ q) (hR : 4*q < R)
    (hl : noncircularGap h R q π < B) (hu : B < noncircularGap h R q 0) :
    ∃! t, t ∈ Ioo 0 π ∧ noncircularGap h R q t = B := by
  have hctx : Continuous (noncircularGap h R q) := by unfold noncircularGap; fun_prop
  obtain ⟨t, ht, hB⟩ := intermediate_value_Icc' pi_pos.le hctx.continuousOn
    (show B ∈ Icc (noncircularGap h R q π) (noncircularGap h R q 0) from ⟨hl.le,hu.le⟩)
  have hi : t ∈ Ioo 0 π := by
    constructor
    · rcases eq_or_lt_of_le ht.1 with e | e
      · subst t; exact False.elim (hu.ne hB.symm)
      · exact e
    · rcases eq_or_lt_of_le ht.2 with e | e
      · subst t; exact False.elim (hl.ne hB)
      · exact e
  refine ⟨t, ⟨hi,hB⟩, ?_⟩
  intro u hu'
  exact (noncircular_monotone h R q hq hR).injOn
    ⟨hu'.1.1.le,hu'.1.2.le⟩ ht (hu'.2.trans hB.symm)


theorem actual_patch_geometry (R B s t : ℝ) :
    circularGap (B-R*cos s) R t = B+R*(cos t-cos s) ∧
    circularGap (B+R*cos s) R (π-t) = B-R*(cos t-cos s) := by
  rw [circular_endpoint_reflection]
  unfold circularGap
  constructor <;> ring

theorem actual_patch_end_reflection (F : ℝ → ℝ) (s : ℝ) :
    (∫ t in 0..s, F (π-t)) = ∫ t in (π-s)..π, F t := by
  simpa using (intervalIntegral.integral_comp_sub_left
    (a := 0) (b := s) F π)

noncomputable def inclinedY (L β R h u t : ℝ) : ℝ :=
  h+β*(u-1/2)+R*(L/sqrt (L^2+β^2))*cos t

theorem inclined_midpoint_and_rise (L β h : ℝ) :
    (h+β*((1/2:ℝ)-1/2)=h) ∧ (h+β*(1-1/2)-(h+β*(0-1/2))=β) := by
  constructor <;> ring

theorem inclined_physical_axial_cuts (L β R h B u : ℝ) (hb : β ≠ 0)
    (hu : u ∈ Icc 0 1) :
    let r := R*L/sqrt (L^2+β^2)
    (inclinedY L β R h u 0 = B ↔ u = 1/2+(B-h-r)/β) ∧
    (inclinedY L β R h u π = B ↔ u = 1/2+(B-h+r)/β) := by
  dsimp
  have hfac : R*(L/sqrt (L^2+β^2)) = R*L/sqrt (L^2+β^2) := by ring
  simpa [inclinedY, hfac, sub_eq_add_neg] using
    inclined_axial_crossings h β (R*L/sqrt (L^2+β^2)) B u hb


noncomputable def inclinedX (L β R u t : ℝ) : ℝ :=
  L*u-R*(β/sqrt (L^2+β^2))*cos t

noncomputable def inclinedZ (R t : ℝ) : ℝ := R*sin t

theorem inclined_surface_derivatives (L β R h u t : ℝ) :
    let ell := sqrt (L^2+β^2)
    deriv (fun z => inclinedX L β R z t) u = L ∧
    deriv (fun z => inclinedY L β R h z t) u = β ∧
    deriv (fun _ : ℝ => inclinedZ R t) u = 0 ∧
    deriv (inclinedX L β R u) t = R*(β/ell)*sin t ∧
    deriv (inclinedY L β R h u) t = -R*(L/ell)*sin t ∧
    deriv (inclinedZ R) t = R*cos t := by
  dsimp
  refine ⟨?_, ?_, deriv_const _ _, ?_, ?_, ?_⟩
  · have H := ((hasDerivAt_id u).const_mul L).sub_const
      (R*(β/sqrt (L^2+β^2))*cos t)
    have Hd : HasDerivAt (fun z => inclinedX L β R z t) L u := by
      convert H using 1 <;> (try rfl) <;> simp [id]
    exact Hd.deriv
  · have H := (((hasDerivAt_id u).sub_const (1/2)).const_mul β).const_add h |>.add_const
      (R*(L/sqrt (L^2+β^2))*cos t)
    have Hd : HasDerivAt (fun z => inclinedY L β R h z t) β u := by
      convert H using 1 <;> (try rfl) <;> simp [id]
    exact Hd.deriv
  · have H := ((hasDerivAt_cos t).const_mul (R*(β/sqrt (L^2+β^2)))).const_sub (L*u)
    have Hd : HasDerivAt (inclinedX L β R u) (R*(β/sqrt (L^2+β^2))*sin t) t := by
      convert H using 1 <;> (try rfl) <;> ring
    exact Hd.deriv
  · have H := ((hasDerivAt_cos t).const_mul (R*(L/sqrt (L^2+β^2)))).const_add
      (h+β*(u-1/2))
    have Hd : HasDerivAt (inclinedY L β R h u) (-R*(L/sqrt (L^2+β^2))*sin t) t := by
      convert H using 1 <;> (try rfl) <;> ring
    exact Hd.deriv
  · exact ((hasDerivAt_sin t).const_mul R).deriv

theorem inclined_surface_area_from_map (L β R h u t : ℝ) (hL : 0 < L) (hR : 0 < R) :
    let xu := deriv (fun z => inclinedX L β R z t) u
    let yu := deriv (fun z => inclinedY L β R h z t) u
    let xt := deriv (inclinedX L β R u) t
    let yt := deriv (inclinedY L β R h u) t
    let zt := deriv (inclinedZ R) t
    sqrt ((yu*zt)^2+(-xu*zt)^2+(xu*yt-yu*xt)^2) = R*sqrt (L^2+β^2) := by
  obtain ⟨xu, yu, _, xt, yt, zt⟩ := inclined_surface_derivatives L β R h u t
  dsimp
  rw [xu,yu,xt,yt,zt]
  exact (inclined_gap_and_area L β R h u t hL hR).2


theorem circular_force_closed_form (k w a h R L : ℝ) (ha : 0 < a) (hR : 0 < R) :
    let θa := arccos ((a-h)/R)
    let θ0 := arccos ((0-h)/R)
    let P := fun t => (4*w/a^2)*gapPrimitive h R t-(4*w/a^4)*cubePrimitive h R t
    (-deriv (fun height => 2*R*L * ∫ t in 0..π, potential k w a (circularGap height R t)) h) =
      -2*R*L*((P θ0-P θa)+k*(gapPrimitive h R π-gapPrimitive h R θ0)) := by
  dsimp
  rw [circular_energy_force k w a R L h ha, circular_closed_form_actual k w a h R ha hR]

noncomputable def neighboringDerivative (k w a : ℝ) (atRange atZeroAngle : Bool) : ℝ → ℝ :=
  if atRange then
    if atZeroAngle then adhesionSlope w a else fun _ => 0
  else
    if atZeroAngle then fun g => k*g else adhesionSlope w a

noncomputable def endpointGap (R B s : ℝ) (atZeroAngle : Bool) (t : ℝ) : ℝ :=
  if atZeroAngle then circularGap (B-R*cos s) R t
  else circularGap (B+R*cos s) R (π-t)

theorem endpoint_gap_mismatch (k w a R s t : ℝ) (atRange atZeroAngle : Bool) :
    let B := if atRange then a else 0
    let g := endpointGap R B s atZeroAngle t
    slope k w a g - neighboringDerivative k w a atRange atZeroAngle g =
      patchMismatch k w a atRange atZeroAngle (R*(cos t-cos s)) := by
  have hA := actual_patch_geometry R a s t
  have h0 := actual_patch_geometry R 0 s t
  simp only [zero_sub, zero_add] at h0
  cases atRange <;> cases atZeroAngle <;>
    simp [endpointGap, neighboringDerivative, patchMismatch, hA.1, hA.2, h0.1, h0.2,
      zero_sub, add_comm]

theorem actual_patch_force_limit (k w a R L : ℝ)
    (atRange atZeroAngle : Bool) (ha : 0 < a) (hR : 0 < R) :
    let C := patchCoefficients k w a atRange atZeroAngle
    let B := if atRange then a else 0
    Tendsto (fun δ : ℝ =>
      let s := arccos (1-δ/R)
      let gap := endpointGap R B s atZeroAngle
      (-(2*R*L * (∫ t in 0..s,
        slope k w a (gap t) - neighboringDerivative k w a atRange atZeroAngle (gap t)) /
         (sqrt δ)^3))) (𝓝[>] 0)
      (𝓝 (-((2*R*L*C.1*R/3)/(sqrt (R/2))^3))) := by
  have H := (actual_circular_patch_delta_limit k w a R L atRange atZeroAngle ha hR).neg
  exact H.congr' (Eventually.of_forall fun δ => by
    dsimp only
    simp_rw [endpoint_gap_mismatch k w a R (arccos (1-δ/R)) _ atRange atZeroAngle])

theorem actual_patch_force_coefficient_nonzero (R L : ℝ) (atRange atZeroAngle : Bool)
    (hR : 0 < R) (hL : 0 < L) :
    let C := patchCoefficients 20 (2/1000) (1/10) atRange atZeroAngle
    (-((2*R*L*C.1*R/3)/(sqrt (R/2))^3)) ≠ 0 := by
  dsimp
  apply neg_ne_zero.mpr
  apply endpoint_limit_nonzero R L _ hR hL
  cases atRange <;> cases atZeroAngle <;> norm_num [patchCoefficients]

#print axioms circular_force_closed_form
#print axioms endpoint_gap_mismatch
#print axioms actual_patch_force_limit
#print axioms actual_patch_force_coefficient_nonzero
#print axioms inclined_surface_derivatives
#print axioms inclined_surface_area_from_map
#print axioms actual_patch_geometry
#print axioms actual_patch_end_reflection
#print axioms inclined_midpoint_and_rise
#print axioms inclined_physical_axial_cuts
#print axioms actual_circular_patch_delta_limit
#print axioms noncircular_bracketed_root
#print axioms actual_circular_jump_at_zero
#print axioms actual_circular_jump_at_range
#print axioms circular_closed_form_actual
#print axioms noncircular_metric_is_speed
#print axioms actual_endpoint_branch_mismatches
#print axioms circular_crossing_chain_factor
#print axioms circular_angular_slope_jump
#print axioms potential_C1
#print axioms circular_closed_root_range
#print axioms circular_endpoint_cases
#print axioms circular_prescribed_heights
#print axioms compression_closed_branch
#print axioms adhesion_closed_branch
#print axioms separated_closed_branch
#print axioms actual_circular_segment
#print axioms circular_segment_assembly
#print axioms noncircular_numeric_monotonicity
#print axioms inclined_axial_crossings
#print axioms force_absolute_error_scaling
#print axioms endpoint_width_delta_limit
#print axioms endpoint_limit_nonzero
#print axioms endpoint_width_from_gap
#print axioms endpoint_width_vanishes
#print axioms circular_patch_delta_limit
#print axioms circular_patch_exact_leading
#print axioms circular_patch_scaled_remainder
#print axioms circular_endpoint_reflection
#print axioms noncircular_reduces_to_circle
#print axioms compression_adhesion_mismatch
#print axioms adhesion_separation_mismatch
#print axioms circular_patch_remainder
#print axioms endpoint_width_ratio
#print axioms endpoint_scaled_linear_limit
#print axioms endpoint_linear_cosine_integral
#print axioms endpoint_linear_cosine_limit
#print axioms inclined_gap_and_area
#print axioms circular_gap_primitive
#print axioms circular_gap_cube_primitive
#print axioms circular_segment_integral
#print axioms endpoint_angle_quadratic_limit
#print axioms circular_energy_force
#print axioms noncircular_energy_force
#print axioms translation_energy_derivative
#print axioms potential_derivative
#print axioms slope_continuous
#print axioms assembled_slope_jump_parameters
#print axioms compression_derivative
#print axioms adhesion_derivative
#print axioms branch_joins
#print axioms branch_curvature_at_left
#print axioms branch_curvature_at_right
#print axioms circular_interior_root
#print axioms circular_interior_necessary
#print axioms circular_unique
#print axioms noncircular_quadratic
#print axioms noncircular_monotone
#print axioms noncircular_arc_factor
#print axioms inclined_basis
#print axioms force_scaling
#print axioms primitive_cos2
#print axioms primitive_cos3
#print axioms cubic_gap_expansion
#print axioms endpoint_quadratic
#print axioms local_parabolic_coefficient
end MRC
