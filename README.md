# pulsatile-transit
Stochastic process describing fluid transit in gastrointestinal tract

Fasted-state gastrointestinal (GI) fluid transit occurs not as a first-order, deterministic process, which is almost always an averaged approximation, but rather as a discrete process involving fluid packets. The transit of these fluid packets is known to occur in pulses which are irregular in both size and time. The influence this phenomenon has on the concentration presented at the absorption site should not be underestimated, and indeed this physiological factor presents a source of inter- and intra-subject variability affecting plasma levels, thus having bioequivalence (BE) implications for immediate release as well as controlled release dosage forms. 

A novel, stochastic process was developed describing cyclical fasted-state gastrointestinal motility and how GI contents are propelled in discrete fluid pockets based on motility phase. The migrating motor complex (MMC) is described as before relative to the time of dosing using a Fourier series approximation while a non-homogenous Poisson process is used to capture the inter-pulse timings. Two random processes control fluid transit: the volume of the random fluid packet formed, and the time the fluid packet pulses out of the compartment. Since a Poisson process is used, the pulsing of fluid packets is memory-less and determined by the parameter λ(t) which evolves with the MMC. The volumetric effect is also taken into account, where a larger fluid volume induces distension and is release more rapidly than a smaller volume. Since this is a stochastic model, many samples must be drawn to glean appropriate statistics.

Files:

[gi_fluids]
  Run by calling the desired function ("optimize", "emptying", or "lambda") and subsequent arguments
  - optimize dose dose_vol ITS
    - optimize the parameters give a dose and initial fluid volume (dose_vol) using ITS number of iterations  
  - emptying dose_vol ITS
    - output gastric emptying profiles of an initial fluid volume (dose_vol) using ITS number of iterations  
  - lambda
    - output the MMC-dependent parameter $\lambda$ for the Poisson process
