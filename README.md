# Dichroic-CD-model
Dichroic CD model for fitting observed ellipticities

Both folders (Q_total and Q_double_H) should be in the same directory as main script (single_fit_dichroic_estimator.py).

Three inputs are needed:
1) N- chain lenght in residues (up to 32)
2) temperature in °C
3) molar ellipticity: in deg cm2 dmol-1 per peptide unit

Program returns propagation parameter w and fractional helicity (<fH>).

Rest of the parameters are fixed according to the global fit of CD data in article "Analysis of the peptide helicity using an ensemble spectroscopic model with re-calibrated parameters".


[θ]H∞ = -41,000

∂[θ]H∞/∂T = 100

k = 3.4

[θ]C = 2,100

∂[θ]C/∂T =-45

v = 0.07
