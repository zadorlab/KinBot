# Unreleased

Transition-state conformer population filtering uses the
`IdealGasThermo(ignore_imag_modes=...)` API introduced in ASE 3.23. The full
saddle frequency list is passed to ASE so that the imaginary reaction mode is
removed without also dropping a real vibration. This does not change the
package's existing minimum dependency of ASE 3.26.
