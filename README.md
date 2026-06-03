# One-parameter Λ-cosmology

LaTeX source for the paper **"One-parameter Λ-cosmology"** by Oleg Ponfilenok
(independent researcher, [ponfil@gmail.com](mailto:ponfil@gmail.com)).

## Summary

The paper argues that a flat dust + cosmological-constant (ΛCDM) universe is, in its
dynamics, a strictly **one-parameter** system: its physics is set by the cosmological
constant Λ alone, while present-day observables additionally require only a single
dimensionless moment of observation τ₀.

Key points:

- The whole evolution (scale factor, Hubble parameter, densities, background
  temperature, total causal radius) is written through the single function `sinh(τ)`,
  with `τ = (3/2) H_Λ t`.
- A **total causal radius** `R = R_p + R_e` (particle + event horizon) — the maximal
  causally accessible scale over all observation time — yields the closed-form constant
  `I = Γ(1/6)Γ(1/3)/(3√π) ≈ 2.804`.
- A **vacuum postulate** tied to the Károlyházy uncertainty,
  `λ̄_Λ = δR_max = ξ (R_max ℓ_P²)^{1/3} = √(R_Λ ℓ_P) = r_s`, fixes the single
  dimensionless coefficient `ξ = (2/I)^{1/3} ≈ 0.894` analytically via a holographic
  (Hawking-temperature) condition at the matter–Λ equality epoch.
- With the CMB temperature `T₀` as the only present-epoch input (fixing τ₀) and the
  measured baryon fraction, the model reproduces `H₀ ≈ 67 km/s/Mpc`, `Ω_m,0 ≈ 0.31`,
  and the photon-to-baryon ratio `N_γ/N_b ≈ 1.6×10⁹`.
- It offers a **holographic reformulation** of the cosmological constant problem,
  expressing the observed `ρ_Λ` through the IR–UV scale `λ̄_Λ = √(R_Λ ℓ_P)` instead of
  the Planck scale.

## Files

| File | Description |
|------|-------------|
| `one_parameter_cosmology_ru.tex` | Source — Russian (original) |
| `one_parameter_cosmology_en.tex` | Source — English translation |
| `refs.bib` | Shared bibliography |
| `one_parameter_cosmology_ru.pdf` | Compiled Russian PDF |

## Building

Requires a LaTeX distribution (TeX Live / MiKTeX) with the `comfortaa`, `babel`,
`amsmath`, `hyperref`, and `xcolor` packages.

```bash
pdflatex one_parameter_cosmology_en
bibtex   one_parameter_cosmology_en
pdflatex one_parameter_cosmology_en
pdflatex one_parameter_cosmology_en
```

(Replace `_en` with `_ru` to build the Russian version.) Both documents share
`refs.bib`. The English source uses `[english]{babel}` and `[T1]{fontenc}`; the
Russian source uses `[russian]{babel}` and `[T2A]{fontenc}`.

## Keywords

cosmological constant · vacuum fluctuations · Károlyházy uncertainty ·
holographic principle · photon-to-baryon ratio

## Status

Preprint / work in progress by an independent researcher. Comments welcome.
