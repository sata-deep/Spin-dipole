# Spin-dipole
**Compute atomic spin-dipole term Tz for a non-collinear (SOC) DFT calculation**

This code (SD.py) calculates the dipole term from the output of the VASP calculation.

The dipole term is given by:
T_z=\left\langle\hat{T}_z\right\rangle=\left\langle\frac{1}{2}[\sigma-3 \hat{\mathbf{r}}(\hat{\mathbf{r}} \cdot \sigma)]_z\right\rangle

$=\frac{1}{2} \int d^3 r\left[m_z(\mathbf{r})-3 \hat{r}_z(\hat{\mathbf{r}} \cdot \mathbf{m}(\mathbf{r}))\right] .$
