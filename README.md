## Orchard Sprayer Tower Dynamics

**ORCHARD: Orchard Sprayer Tower Dynamics** is a comprehensive Matlab tool designed to simulate and animate the nonlinear stochastic dynamics of an orchard sprayer tower moving on an irregular terrain. The terrain is emulated by a stationary Gaussian random process, providing realistic and valuable insights into the dynamic behavior of the sprayer tower. Developed with an educational approach, **ORCHARD** is intuitive and user-friendly, making it accessible for researchers and engineers in mechanical and agricultural engineering.

### Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Usage](#usage)
- [Documentation](#documentation)
- [Authors](#authors)
- [Citing ORCHARD](#citing-orchard)
- [License](#license)
- [Institutional Support](#institutional-support)
- [Funding](#funding)

### Overview
**ORCHARD** was developed to conduct numerical simulations related to the stochastic dynamics of an orchard sprayer tower moving through irregular terrain. The results have been published in several peer-reviewed works:
- **A. Cunha Jr, J. L. P. Felix, and J. M. Balthazar**, *Quantification of parametric uncertainties induced by irregular soil loading in orchard tower sprayer nonlinear dynamics*, Journal of Sound and Vibration, vol. 408, pp. 252-269, 2017. [DOI](https://doi.org/10.1016/j.jsv.2017.07.023)
- **A. Cunha Jr, J. L. P. Felix, and J. M. Balthazar**, *Exploring the nonlinear stochastic dynamics of an orchard sprayer tower moving through an irregular terrain*, In: Mohamed Belhaq. (Org.). Recent Trends in Applied Nonlinear Mechanics and Physics. 1ed.: Springer International Publishing, p. 203-213, 2018. [DOI](http://dx.doi.org/10.1007/978-3-319-63937-6_11)
- **R. N. Silva et al.**, *On a Vehicular Suspension for a Non-ideal and Nonlinear Orchard Tower Sprayer Through an Inverted Pendulum Using Reologic Magneto (MR)*, In: J. M. Balthazar, Nonlinear Vibrations Excited by Limited Power Sources, Springer Cham, 2022. [DOI](https://doi.org/10.1007/978-3-030-96603-4_10)

### Features
- Simulates nonlinear stochastic dynamics of orchard sprayer towers
- Emulates irregular terrain using stationary Gaussian random processes
- Animates system dynamics in Matlab for visual analysis
- Intuitive Matlab implementation
- Educational style code for easy understanding

### Usage
To get started with **ORCHARD**, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/americocunhajr/ORCHARD.git
   ```
2. Navigate to the code directory:
   ```bash
   cd ORCHARD/ORCHARD-1.0
   ```
3. For a deterministic simulation with sinusoidal excitation, execute:
   ```bash
   main_orchard_ipv_sinusoidal
   ```
4. For a stochastic simulation with random excitation, execute:
   ```bash
   main_orchard_ipv_kl
   ```
5. For Fourier analysis in the case of stochastic excitation, execute:
   ```bash
   main_orchard_fourier_kl
   ```
6. For Monte Carlo simulation with sinusoidal excitation, execute:
   ```bash
   main_orchard_mc_sinusoidal
   ```
7. For Monte Carlo simulation with random sinusoidal excitation, execute:
   ```bash
   main_orchard_mc_sinusoidal
   ```
8. For Monte Carlo simulation with stochastic excitation, execute:
   ```bash
   main_orchard_mc_kl
   ```
9. For a parametric study changing the stiffness, execute:
   ```bash
   main_orchard_ps_sinusoidal
   ```

### Documentation
The routines in **ORCHARD** are well-commented to explain their functionality. Each routine includes a description of its purpose, as well as inputs and outputs.

### Authors
- Americo Cunha Jr
- Jorge Felix
- José Manoel Balthazar

### Citing ORCHARD
If you use **ORCHARD** in your research, please cite the following publication:
- *A. Cunha Jr, J. L. P. Felix, and J. M. Balthazar, Quantification of parametric uncertainties induced by irregular soil loading in orchard tower sprayer nonlinear dynamics, Journal of Sound and Vibration, vol. 408, pp. 252-269, 2017 https://doi.org/10.1016/j.jsv.2017.07.023*
- *A. Cunha Jr, J. L. P. Felix, and J. M. Balthazar, Exploring the nonlinear stochastic dynamics of an orchard sprayer tower moving through an irregular terrain, In: Mohamed Belhaq. (Org.). Recent Trends in Applied Nonlinear Mechanics and Physics. 1ed.: Springer International Publishing, p. 203-213, 2018 http://dx.doi.org/10.1007/978-3-319-63937-6_11*
- *R. N. Silva et al., On a Vehicular Suspension for a Non-ideal and Nonlinear Orchard Tower Sprayer Through an Inverted Pendulum Using Reologic Magneto (MR), In: J. M. Balthazar, Nonlinear Vibrations Excited by Limited Power Sources, Springer Cham, p. 151-173, 2022 https://doi.org/10.1007/978-3-030-96603-4_10*

```
@article{cunhajr2017p252,
   author  = {A. {Cunha~Jr} and J. L. P. Felix and J. M. Balthazar},
   title   = {Quantification of parametric uncertainties induced by irregular soil loading in orchard tower sprayer nonlinear dynamics},
   journal = {Journal of Sound and Vibration},
   year    = {2017},
   volume  = {408},
   pages   = {252-269},
   doi     = {10.1016/j.jsv.2017.07.023},
}
```

```
@incollection{CunhaJr2018RTANM,
   author    = {A. {Cunha~Jr} and J. L. P. Felix and J. M. Balthazar},
   title     = {Exploring the {N}onlinear {S}tochastic {D}ynamics of an {O}rchard {S}prayer {T}ower {M}oving {T}hrough an {I}rregular {T}errain},
   editor    = {M. Belhaq},
   booktitle = {Recent Trends in Applied Nonlinear Mechanics and Physics: Selected papers from CSNDD 2016 (Springer Proceedings in Physics)},
   publisher = {Springer},
   address   = {Cham},
   year      = {2018},
   pages     = {203-213},
   doi       = {10.1007/978-3-319-63937-6_11},
}
```

```
@incollection{Silva2022NVELPS,
   author    = {R. N. Silva et al.},
   title     = {On a Vehicular Suspension for a Non-ideal and Nonlinear Orchard Tower Sprayer Through an Inverted Pendulum Using Reologic Magneto {(MR)}},
   editor    = {M. Belhaq},
   booktitle = {Nonlinear Vibrations Excited by Limited Power Sources},
   publisher = {Springer},
   address   = {Cham},
   year      = {2022},
   pages     = {151-173},
   doi       = {10.1007/978-3-030-96603-4_10},
}
```

### License

**ORCHARD** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

<img src="logo/mit_license_red.png" width="10%"> 

### Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%"> &nbsp; &nbsp; <img src="logo/logo_uffs_color.png" width="25%"> &nbsp; &nbsp; <img src="logo/logo_unesp_color.png" width="20%"> 

### Funding

<img src="logo/faperj.jpg" width="20%"> &nbsp; &nbsp; <img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%">
