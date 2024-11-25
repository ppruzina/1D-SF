# 1D-SF
Codes for a one-dimensional model of salt fingering staircases, as presented in

Pružina, P., Hughes, D., & Pegler, S. (2023). Salt fingering staircases and the three-component Phillips effect. Journal of Fluid Mechanics, 968, A16. doi:10.1017/jfm.2023.534

## List of files
Code is written in Matlab
Files include: 
1. ddcsolver - solves the model equations using the MATLAB pdepe solver
2. ddcinstability - computes the values of several instability markers
3. ddcwavenumbergr - produces wavenumber - growth rate plots for given parameters
4. ddcenergyplotter - plots the steady state energy

## The model
A one-dimensional horizontally averaged model, derived from the Boussinesq equations via an averaging process, using a mixing-length assumption for closure. The dimensionless model for temperature $T(z,t)$, salinity $S(z,t)$ and kinetic energy $e(z,t)$ takes the form

$$T_t = \left(\frac{D^2}{D+1}T_z\right)_z,$$

$$S_t = \left(\frac{D^2}{D+\tau}S_z\right)_z,$$

$$e_t = \left(\frac{D^2}{D+\sigma}e_z\right)_z + \sigma e\_{zz} - \sigma \left(\frac{D^2}{D+1}T_z - \frac{D^2}{D+\tau}S_z\right) - \epsilon \frac{D^2}{e}. $$

The model depends on the dimensionless parameters $\sigma$ (Prandtl number, 10), $\tau$ (diffusivity ratio, 0.01) and $\epsilon$ (dissipation parameter, 0.02)
The system is closed using the turbulent diffusivity $D = le^{1/2}$, where $l$ is a turbulent mixing length.
The form for the diffusivity is 
$$D=\frac{\sqrt{e^2+\delta R^2}}{R}$$,
chosen for its qualitative dependence on $e$ and $R$. $\delta$ is an additional dimensionless parameter, with suitable values $<<1$



