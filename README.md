# Stark-MBL
These codes are used to probe for Stark MBL in a periodically driven system.

The IPR code obtains IPRs averaged over several vectors in the spectral bulk, and calculates them as a function of driving frequency.
The corresponding figure indicates resonant points at specific frequencies.

The entropy code calculates the time evolution of bipartite entanglement entropies, averaged overr all initial Fock states for which the number
of domains equals L-3. The values obtained at the end of each evolution are plotted against system size L to show volume law scaling at
the resonant points.

The tau code calculates the scaling exponents of the IPRs extrapolated to infinite size by plotting as a function of 1/lnD. The extrapolated
exponents indicate exponential localisation away from the resonant points and multifractality at the resonant points.
