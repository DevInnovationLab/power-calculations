# CSE AI Chatbot Power Analysis 

**Author:** Nandita (edited using Claude) · **First review by:** Claude
ESI en Valores RCT (Dominican Republic) · chatbot component


> Note: this document lays out the proposed approach for review. Results are held until the approach is signed off, then we run and fill them in. The *Claude comment* blocks under each point are a first review. **Luiza:** please mark the response line under each point — change `[ ]` to `[x]` for Yes / Discuss / No and add comments. Your edits are the record.

---

## What's available (inputs)

**Design (from the study design doc):**

| Parameter | Value |
|---|---|
| Treatment schools | 450 (225 Regular CSE / 225 CSE+Debiasing) |
| Surveyed students per school | 90 (30 per grade × 3 grades) |
| Total students in chatbot randomization | 40,500 |
| Per chatbot arm (β₂) | 20,250 |
| Per 2×2 cell (β₃) | 10,125 |
| Chatbot randomization | Individual student, within school, stratified by CSE arm |
| Power / significance | 80% / 5%, two-sided |
| ICC — pregnancy | 0.050 |
| ICC — overconfidence score | 0.065 |
| ICC — depression / knowledge | 0.100 |
| ICC — test scores | 0.200 |
| Pregnancy base rate (1 yr / 4 yr) | 0.056 / 0.168 |
| Overconfidence mean / SD | 10.4 / 27.7 |

**Other documents for reference:**
- Existing spreadsheet — MDEs for the school-level estimands (CSE vs control; debiasing vs CSE).
- Tomás's `PowerCalcsChatbot.Rmd` — analytic MDEs (`compute_mde`, lines 185–204), the saturation information-loss block (lines 329–342), and a power simulation (`powersimu`, lines 413–448).

---

## Assumptions

| Assumption | Choice | Type / note |
|---|---|---|
| Spillover channel | Untreated students benefit linearly in their school's treated share | **Our modelling choice** — the key one |
| Spillover magnitude | Swept 0.05–0.30 SD | **Placeholder** — replace with pilot estimate when available |
| Saturation levels | 25% / 50% / 75%, equal thirds, stratified by CSE arm | Design option under consideration |
| Debiasing coding | Centered (`Deb − 0.5`) | Spec choice so β₂ is the pooled effect (see §3) |
| Estimator | School fixed effects + cluster-robust SEs (`pyfixest`) | Standard library, no hand-rolled econometrics |
| Female share (girls-only rows) | 50% | From design doc |
| Simulation replications | 1,000 (fewer for heavy grids) | Numerical precision dial |


> **Luiza — do these assumptions look right?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

---

## 1. The estimands we target

The chatbot is randomized at the individual student level *within* the 450 treatment schools, on top of the school-level CSE / CSE+debiasing split. That gives three quantities of interest:

- **β₂ — direct effect of chatbot access**, averaged across both CSE arms (the study's Estimand 1: cells B+D vs A+C).
- **β₃ — chatbot × debiasing interaction** (Estimand 2): is the chatbot more effective alongside debiasing?
- **δ — spillover / saturation effect**: does a student's outcome move with the *share* of their peers who have the chatbot? From the background note, my understanding is that the whole reason we're using saturation levels is to measure spillovers. 

> **Claude comment:** Treating δ as first-class is the right call. Two caveats. (1) The design doc also raises *nonlinear* network/equilibrium effects (norms shifting only past some saturation threshold); a single linear δ can't capture those, so it's worth running at least one nonlinear saturation form as a robustness check. (2) With ~9 outcome×subgroup rows and three estimands each, multiple-testing will matter for interpretation even though it doesn't change MDEs — pre-specifying primary outcomes now would help.

> **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

---

## 2. The model 

We simulate the full design and generate each student's outcome from an explicit model. Let `T` = chatbot access (0/1), `Deb` = debiasing school (0/1), and `s_g` = realized treated share in school *g*:

```
Y_ig = γ_g + β2·T_ig + β3·(Deb_g − 0.5)·T_ig + δ·s_g·(1 − T_ig) + ε_ig
```

- `γ_g` is a school-level term carrying the between-school variance; its size is set by each outcome's **ICC** (variance ICC), and `ε` carries the within-school variance (1 − ICC). So the total variance is 1 and effects read in SD units.
- `δ·s_g·(1 − T_ig)` is the **spillover**: untreated students benefit in proportion to how saturated their school is (information passed from treated peers). This single term is what lets us study the bias and the spillover estimand.
- The debiasing indicator is **centered** (`Deb − 0.5`) so that β₂ is the *average* chatbot effect across both CSE arms (see §3).

Saturation is set per school: the uniform design fixes `s_g = 0.5` everywhere; the saturation design assigns `s_g ∈ {0.25, 0.50, 0.75}` in equal thirds, stratified by CSE arm.

> **Claude comment:** The spillover here accrues only to *untreated* students (free-riding). That's a reasonable lead case, but it's what drives the headline "−0.5·δ bias," so flag it as a functional-form assumption and stress-test it against a *symmetric* spillover (benefiting treated and untreated alike), under which the FE estimate is **not** biased — the two cases give very different recommendations. Also: real spillover probably runs at the classroom/grade level, not the whole school; if so, the relevant clustering and the saturation variation change, and school-level modelling may overstate the between-cluster variation we have. Minor: adding the spillover term mechanically raises total outcome variance, so hold σ_total fixed when reading MDEs in "SD units."

> **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

---

## 3. Estimation — two specifications

Because the chatbot varies within school, the natural estimator absorbs **school fixed effects**, which removes the between-school variance for free. We run two specifications (with school-clustered SEs):

**(A) Direct-effect spec (school FE):**
```
Y ~ chatbot + (Deb_c · chatbot) | school
```
recovers β₂ and β₃. *Centering: With a raw 0/1 debiasing variable, the chatbot coefficient would be the effect in regular-CSE schools only, estimated off half the schools with an SE root2 larger. Centering makes it the pooled average. 

**(B) Spillover spec (no FE; saturation as a regressor):**
```
Y ~ chatbot + school_share + (Deb_c · chatbot)   [clustered by school]
```
The coefficient on `chatbot` recovers the **unbiased** direct effect; the coefficient on `school_share` recovers the spillover δ. This is only identified when `s_g` varies — i.e. under the saturation design.

> **Claude comment:** Centering is correct and standard — low risk. The thing to watch is Spec B: dropping school FE to identify δ puts the between-school variance back into the error, so it's less efficient for β₂ and leans entirely on saturation being randomly assigned (it is, so no bias — but report both specs and the gap between them, since that gap *is* the spillover story). If spillover can also reach treated students, Spec B as written (`school_share` main effect only) won't fully separate direct from spillover effects — we'd want the `school_share × chatbot` term too. And for rare binary outcomes (pregnancy at 5.6%), FE on a linear-probability model is an approximation: fine for MDEs, but the final analysis may want a GLM or a robustness check.

> **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

---

## 4. Simulation + Validation

We build a virtual study with a known planted effect, run the estimator, and record whether it's detected; the fraction detected over many repetitions will be the power we need.

### How we will build the virtual study:

Each "study" is a table of 40,500 students that we construct from scratch — and because we plant the true effects ourselves, we know the right answer to check against.

0. **Set the knobs we control:** the true effects (β₂, β₃, δ), the outcome's ICC, and the design (uniform 50/50 or saturation 25/50/75). We choose these, so we know them.
1. **Make the schools and assign debiasing:** create 450 schools; mark 225 at random as debiasing, 225 as regular CSE (a school-level decision).
2. **Give each school its chatbot share:** uniform -> every school 50%; saturation -> 25% / 50% / 75% to equal thirds of schools (150 each), done separately within each CSE arm so they stay balanced.
3. **Pick which students get the chatbot, within each school:** for a 75%-share school we treat exactly 68 of its 90 students (90 × 0.75 ≈ 68), chosen at random — exact count, no accidental imbalance.
4. **Create the noise in two layers** (this is where the ICC enters): one *school-level* draw shared by all 90 students (variance = ICC, which makes classmates resemble each other) plus one *individual* draw per student (variance = 1 − ICC). They sum to total variance 1, so the outcome is in SD units. For ICC = 0.10 that's a school SD of √0.10 ≈ 0.32 and an individual SD of √0.90 ≈ 0.95.
5. **Build each student's outcome** by plugging into the model equation (§2): school draw + chatbot effect (if treated) + interaction (in debiasing schools) + spillover (for untreated, ∝ school share) + individual noise.
6. **Result:** one complete dataset (one row per student) ready to hand to the regression.

*Worked example by claude based on approach shared* (β₂ = 0.05, δ = 0.10, ICC = 0.10, a debiasing school assigned 50% saturation that drew a school value of +0.20):
- A **treated** student with individual noise −0.30 → Y = 0.20 + 0.05 − 0.30 = **−0.05**.
- An **untreated** classmate with individual noise +0.10 → Y = 0.20 + (0.10 × 0.50) + 0.10 = **0.35** (the +0.05 is the spillover from sitting in a half-saturated school).

We then repeat steps 1 to 6 a thousand times with fresh random draws; the share of runs in which the estimator detects the effect is the power. 

For validation, we check if: (i) simulated standard errors match the analytical formula to three decimals; (ii) the interaction SE comes out exactly 2× the main-effect SE, as theory predicts; (iii) the false-positive rate is 0.05 and every planted parameter is recovered. 

> **Claude comment:** Validation is solid for the simple case — but note it currently checks the *no-spillover* analytic formula and doesn't independently verify the spillover-spec SEs, so I'd add a second benchmark for δ (e.g. a school-level aggregated regression) before trusting those power numbers. Errors are assumed Gaussian; for rare binary outcomes the finite-sample coverage should be checked separately. 

> **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

---

## NG analysis on how it compares to Tomás's approach and can build on it

**What we keep:** his insight that, because the chatbot is randomized within school, the relevant noise is the within-school SD `σ_within = σ_total·√(1 − ICC)`. We reproduce his analytic MDEs (`compute_mde`, lines 185–204).

**What we add, point by point:**

1. **An explicit spillover term.** His model (and the analytic MDEs) have no peer channel, so the bias that motivates the saturation design is invisible. Adding `δ·s_g·(1−T)` lets us quantify the `≈ 0.5·δ` attenuation and recover an unbiased direct effect.

   > **Claude comment:** This is the core contribution and I think it's right — just keep the functional-form caveat from §2 attached to it, so it doesn't read as the only possible spillover model.

   > **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

2. **Power for the spillover itself.** His `powersimu` (lines 413–448) estimates a spillover coefficient, but under different assumptions from the analytic block. We power δ within the same framework.

   > **Claude comment:** This is probably the most decision-relevant finding. If δ is underpowered for high-ICC outcomes like test scores, that's an argument for *more* or *wider* saturation levels (e.g. 0/33/67/100) — worth simulating that option rather than just reporting low power.

   > **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

3. **An explicit, matched specification.** `mde1` (line 192) is implicitly the pooled effect; we make that hold by centering debiasing, so the regression the team runs matches the power calc.

   > **Claude comment:** Clearly correct and low-risk. The only thing to confirm is that the *actual* analysis code uses the centered coding — otherwise the MDE and the estimate won't line up, which is exactly the gap this point is meant to close.

   > **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:

4. **One simulation for everything.** His analytic loss (lines 329–342) and his simulation (lines 413–448) use different assumptions (direct effect fixed at 0.5 SD, line 433; `lm` without school FE, line 437). We produce the MDEs, saturation loss, bias, and spillover power from one consistent engine, validated against his formulas.

   > **Claude comment:** Right principle; the cost is speed and transparency (one big engine is harder to eyeball than a one-line formula). Mitigate by validating the engine against Tomás's closed-form numbers at each step, so consolidation can't quietly hide a bug.

   > **Luiza — agree?**  `[ ] Yes`  ·  `[ ] Discuss`  ·  `[ ] No`  — comments:


> **Claude comment (overall):** The approach is sound and a genuine improvement on a correct foundation. The single biggest risk is the spillover *functional form* (free-riding, linear, school-level): the bias and recovery conclusions are conditional on it, so I'd make testing 2–3 alternative spillover structures part of v1 rather than a follow-up. Second priority is settling whether spillover/saturation lives at the classroom or school level, since that changes what variation we actually have.

---

## Overall: green light to proceed?

If the above sounds reasonable, Nandita takes a first pass at the code on top of Tomás's framework and shares it for Luiza and Tomás to check.

> **Luiza — overall decision:**  `[ ] Yes, go ahead`  ·  `[ ] Let's discuss first`  ·  `[ ] No`  — comments:

*— Nandita (edited using Claude)*
