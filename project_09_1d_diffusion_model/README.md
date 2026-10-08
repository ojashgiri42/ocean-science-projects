# 09. 1D Vertical Heat Diffusion Model (Numerical Methods)

**Question.** How does an initial thermocline profile evolve under constant vertical diffusion with a seasonally varying surface temperature?

**Data.** No external data. The model is fully defined in the notebook.

**Method.** Solves ∂T/∂t = κ ∂²T/∂z² with an explicit forward-time, centred-space (FTCS) finite-difference scheme.

| Item | Value (from the code) |
| :--- | :--- |
| Diffusivity κ | 1×10⁻⁴ m² s⁻¹ (constant) |
| Grid | 50 levels, Δz = 10 m (nodes at 0–490 m) |
| Time step | Δt = 3600 s; 2000 steps = 83.3 days |
| Initial condition | T(z) = 5 + 20·exp(−z / 100 m) °C |
| Surface boundary | Dirichlet: T = 20 + 5·sin(2πt / 365 d) °C (this starts at 20 °C, so there is a 5 °C step from the initial 25 °C at the first time step) |
| Bottom boundary | Dirichlet: T = 5 °C at 490 m |
| Stability | From the parameters, r = κΔt/Δz² = 0.0036, well below the explicit-scheme limit of 0.5. |

There is no mixed-layer or turbulence parameterisation. This is a numerical-methods exercise, not an ocean model.

**Finding.** After 83 days, the profile differs only slightly from the initial profile. The final profile is slightly warmer between roughly 50 and 250 m. This is consistent with the small diffusive length scale √(κt) ≈ 27 m for these parameters.

![09. 1D Vertical Heat Diffusion Model (Numerical Methods)](figures/numerical_model_profile.png)
