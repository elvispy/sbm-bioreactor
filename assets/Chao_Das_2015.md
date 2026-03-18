# Numerical simulation of coupled cell motion and nutrient transport in NASA's rotating bioreactor

Tzu-Chieh Chao, Diganta B. Das

Source files:

- PDF: `assets/Chao_Das_2015_cell_motion_nutrient_NASA.pdf`
- Original text extraction: `assets/Chao_Das_2015_cell_motion_nutrient_NASA.txt`

This is a readable, equation-focused Markdown transcription reconstructed from the PDF because the plain-text extraction corrupted much of the math. I prioritized the modelling sections and notation over full prose transcription.

## Citation

Chao, T.-C., and Das, D. B. (2015). *Numerical simulation of coupled cell motion and nutrient transport in NASA's rotating bioreactor*. *Chemical Engineering Journal*. DOI: `10.1016/j.cej.2014.08.077`

## Core Symbols

| Symbol | Meaning |
| --- | --- |
| $a$ | Cell diameter |
| $C$ | Nutrient concentration |
| $c_s$ | Mass fraction of solid cell particles |
| $D$, $D_f$ | Diffusion coefficient |
| $d$ | Cell density |
| $f_h$ | Hindered settling function |
| $\mathbf{J}_s$ | Particle flux |
| $\mathbf{J}_{sc}$ | Flux contribution from particle interactions |
| $\mathbf{J}_{s\mu}$ | Flux contribution from spatial viscosity variation |
| $k_c$ | Correlation coefficient of nutrient consumption |
| $k_e$ | Cell doubling-time coefficient |
| $k_{sc}$ | Particle-interaction flux coefficient |
| $k_\mu$ | Spatial-variation flux coefficient |
| $L$ | HARV thickness |
| $p$ | Pressure |
| $r_c$ | Nutrient consumption rate |
| $t$ | Time |
| $\mathbf{u}$ | Mean fluid velocity |
| $\mathbf{u}_s$ | Solid-particle velocity |
| $\mathbf{u}_{st}$ | Stokes settling velocity |
| $\mathbf{u}_{slip}$ | Relative velocity between solid and liquid phases |
| $\mathbf{u}_r$ | Rotational wall velocity |
| $\mu$ | Effective mixture viscosity |
| $\mu_f$ | Pure-fluid viscosity |
| $\mu_c$ | Nutrient consumption coefficient |
| $\rho$ | Mixture density |
| $\rho_f$, $\rho_s$ | Fluid and solid density |
| $\rho_f^\circ$, $\rho_s^\circ$ | Phase densities of fluid and solid |
| $\tau_w$ | Shear stress at wall |
| $\Omega$, $\omega$ | Rotational speed |
| $\Phi$, $\Phi_s$ | Particle volume fraction |
| $\Phi_{\max}$ | Maximum packing fraction |
| $\dot{\gamma}$ | Magnitude of shear-rate tensor |

## 2. Methodology

### 2.1 Cell distribution

The paper models cells as non-deformable spherical microparticles suspended in a Newtonian cell culture medium, using a two-phase suspension framework for the HARV bioreactor.

#### Momentum balance

$$
\rho \frac{\partial \mathbf{u}}{\partial t}
+ \rho (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p
- \nabla\cdot\left[\rho c_s (1-c_s)\mathbf{u}_{slip}\mathbf{u}_{slip}\right]
+ \nabla\cdot\left[\mu\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{T}\right)\right]
+ \rho \mathbf{g}
\tag{1}
$$

#### Mixture density

$$
\rho = (1-\Phi_s)\rho_f^\circ + \Phi_s \rho_s^\circ
\tag{2}
$$

#### Krieger-type viscosity law

$$
\mu = \mu_f \left(1-\frac{\Phi_s}{\Phi_{\max}}\right)^{-2.5\Phi_{\max}}
\tag{3}
$$

#### Mixture continuity

$$
(\rho_s^\circ-\rho_f^\circ)\left[\nabla\cdot\left(\Phi_s(1-c_s)\mathbf{u}_{slip}\right)\right]
- \rho_f^\circ (\nabla\cdot\mathbf{u}) = 0
\tag{4}
$$

#### Solid-phase transport

$$
\frac{\partial \Phi_s}{\partial t} + \nabla\cdot(\Phi_s \mathbf{u}_s) = 0
\tag{5}
$$

with

$$
\mathbf{u}_s = \mathbf{u} + (1-c_s)\mathbf{u}_{slip}
$$

so that

$$
\frac{\partial \Phi_s}{\partial t}
- \nabla\cdot\left[\Phi_s\left(\mathbf{u} + (1-c_s)\mathbf{u}_{slip}\right)\right] = 0
\tag{6}
$$

#### Flux form used by Rao et al

$$
\frac{\partial \Phi_s}{\partial t} + \nabla\cdot(\mathbf{u}\Phi) = -\frac{\nabla\cdot\mathbf{J}_s}{\rho_s^\circ}
\tag{7}
$$

$$
\nabla\cdot\mathbf{u}
= \frac{\rho_s^\circ-\rho_f^\circ}{\rho_s^\circ\rho_f^\circ}\,\nabla\cdot\mathbf{J}_s
\tag{8}
$$

$$
\mathbf{u}_{slip} = \frac{\mathbf{J}_s}{\Phi_s d (1-c_s)}
\tag{9}
$$

$$
\frac{\partial \Phi}{\partial t} + \mathbf{u}\cdot\nabla(\Phi)
= \frac{\rho_s^\circ-\rho_f^\circ}{\rho_s^\circ\rho_f^\circ}\,\nabla\cdot\mathbf{J}_s
= \frac{\rho}{\rho_s^\circ\rho_f^\circ}\,\nabla\cdot\mathbf{J}_s
\tag{10}
$$

#### Diffusive flux decomposition

$$
\mathbf{J}_s = \mathbf{J}_{s\mu} + \mathbf{J}_{sc}
\tag{11}
$$

$$
\mathbf{J}_{sc} = -a^2\Phi^2 k_{sc}\nabla(\dot{\gamma}\Phi)
\tag{12}
$$

$$
\mathbf{J}_{s\mu} = -a^2\Phi^2 k_\mu \nabla(\ln \mu)
\tag{13}
$$

#### Governing flux equation with hindered settling

$$
\frac{\mathbf{J}_s}{\rho_s}
= -\left[\Phi D*\Phi \nabla(\dot{\gamma}\Phi) + \Phi^2 D_\mu \dot{\gamma}\nabla(\ln \mu)\right]
- f_h \mathbf{u}_{st}\Phi
\tag{14}
$$

$$
D_\Phi = 0.41 a^2
\tag{15}
$$

$$
D_\mu = 0.62 a^2
\tag{16}
$$

$$
\mathbf{u}_{st} = \frac{2a^2(\rho_s-\rho_f)}{9\mu}\,\mathbf{g}
\tag{17}
$$

$$
f_h = \frac{\mu_f(1-\Phi_{average})}{\mu}
\tag{18}
$$

#### Shear-rate tensor

$$
\dot{\boldsymbol{\gamma}} = \nabla\mathbf{u} + \nabla\mathbf{u}^{T}
\tag{19}
$$

$$
\dot{\gamma}
= \left[\frac{1}{2}\left(\dot{\boldsymbol{\gamma}}\cdot\dot{\boldsymbol{\gamma}}\right)\right]^{1/2}
= \left[\frac{1}{2}\left(4u_x^2 + 2(u_y+v_x)^2 + 4v_y^2\right)\right]^{1/2}
\tag{20}
$$

#### Rearranged cell-volume-fraction equation

$$
\frac{\partial \Phi}{\partial t} + \mathbf{u}\cdot\nabla(\Phi)
= \frac{\rho_f^\circ-\rho_s^\circ}{\rho_s^\circ\rho_f^\circ}
\nabla\cdot \rho_s^\circ
\left\{
\left[0.41a^2\Phi\nabla(\dot{\gamma}\Phi) + 0.62a^2\Phi^2\dot{\gamma}\nabla(\ln\mu)\right]
- \Phi \frac{2\mu_f a^2(1-\Phi_{average})(\rho_s-\rho_f)}{9\mu^2}\mathbf{g}
\right\}
\tag{21}
$$

#### Wall shear and rotational velocity

$$
\tau_w = \frac{4\mu_f(\mathbf{u}-\mathbf{u}_w)}{L}
\tag{22}
$$

$$
\mathbf{u}_r = \omega (y,-x) = \frac{2\pi\cdot rpm}{60}(y,-x)
\tag{23}
$$

Boundary assumptions in the text:

- no cell penetration or dispersion at the wall
- fluid velocity equals wall rotation velocity on the boundary
- suspension velocity satisfies no-slip at all walls

### 2.2 Nutrient transport

The paper considers glucose as the only nutrient and writes the transport equation as:

$$
\frac{dC}{dt} + \mathbf{u}\cdot\nabla C = D_f \nabla^2 C + r_c
\tag{24}
$$

$$
r_c = -\mu_c \cdot d
\tag{25}
$$

### 2.3 Growth kinetics

The transient cell-growth expression shown in the manuscript is:

$$
\frac{\partial d}{\partial t}
= k \cdot k_e \cdot C \cdot d_0 \cdot e^{k_e t}
= k_c \cdot C \cdot d_0 \cdot e^{k_e t}
\tag{26}
$$

Note: later prose refers to a "modified exponential growth from equation (27)", but the displayed equation in the manuscript scan/text layer is numbered (26). I did not invent a missing equation here.

### 3.2 Relation between cell density and volume fraction

$$
\Phi = \frac{\pi}{6} d a^3
\tag{28}
$$

## Parameters Reported in Table 1

| Parameter | Symbol | Value | Unit | Ref. |
| --- | --- | --- | --- | --- |
| Diffusivity of glucose | $D_f$ | $5.4\times10^{-10}$ | m$^2$ s$^{-1}$ | [26] |
| Concentration of glucose | $C$ | 5.5 | mol m$^{-3}$ | [26] |
| Fluid density | $\rho_f$ | 1050 | kg m$^{-3}$ | - |
| CHO cell density | $\rho_s$ | 1000 | kg m$^{-3}$ | - |
| Average CHO cell diameter | $a$ | $5.0\times10^{-6}$ | m | - |
| Medium viscosity | $\mu_f$ | 0.5889 | Pa s | - |
| Rotation speed | - | 7.5 | rpm | - |

## Notes on Reconstruction

- Equations (15) and (16) were absent from the text extraction but were recovered from the PDF page image.
- The original text extraction mixes $\Phi$ and $\Phi_s$ in a few places. I preserved the notation as it appears in the manuscript, even where it looks slightly inconsistent.
- The accepted-manuscript scan has minor notation inconsistencies, especially around equations (7)-(10) and the later reference to equation (27).
