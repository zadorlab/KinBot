# MESS execution and rate validation

An exit code of zero establishes solver execution, not usable kinetics.
`kinbot.mess_execution` records a separate `rate_validation` object for every
completed output. It also audits older reusable completion records without
resubmitting MESS. The ANL kinetics audit labels its existing rate scope as
high pressure and includes the independent finite-pressure result.

To audit an existing native output without modifying it:

```bash
python -B -m kinbot.mess_rate_audit /path/to/mess_0000.out
```

The report distinguishes high-pressure and pressure-dependent matrices,
retains each negative or nonfinite entry, and gives its temperature, pressure,
source, destination, line, and size relative to positive entries in the same
row. Unimolecular and bimolecular rows must not be compared as though their
units were identical. Native missing entries (`***`) are not zeros. Missing
rows caused by native population reduction are not deleted chemistry.
Diagonal loss summaries are not individual transition coefficients.

`passed_sign_check` means only that the parsed matrices passed these necessary
sign/format checks. It does not establish convergence, detailed balance,
coverage of the requested condition grid, correct population reduction, or
chemical accuracy. Missing sections are `not_assessed`; they cannot borrow
success from a different section. No negative rate is clipped, no tolerance
silently legitimizes a negative transition, and no retry or QC calculation is
requested by this audit. An invalid pressure table does not invalidate the
underlying stationary-point QC or automatically invalidate the HPL table.

Regression coverage includes a hash-bound, unedited native output containing
valid HPL coefficients and negative pressure coefficients. Simulated local and
queued orchestration also tests completed-record reuse, delayed output,
submission ambiguity, preservation of live jobs and repeated invocation.
These tests do not replace native solver validation.

## Check the native numerical-library interface

The MESS source's `Lapack::int_t` must match the **linked library's** Fortran
integer ABI. A 64-bit operating system does not imply an ILP64 LAPACK library.
Do not change all installations to 32-bit integers: an ILP64 library requires
the corresponding 64-bit declarations.

`validation/native_mess_lapack_probe.cc` tests this interface with a known
positive-definite factorization, a poisoned status value, and a pivoted solve.
Build it against the same headers, definitions and libraries as MESS. For a
CMake build with a `mess_driver` executable whose link interface supplies those
headers and libraries, add the following after that target is defined:

```cmake
add_executable(lapack_probe /absolute/path/to/KinBot/validation/native_mess_lapack_probe.cc)
target_link_libraries(lapack_probe PRIVATE $<TARGET_PROPERTY:mess_driver,LINK_LIBRARIES>)
# Apply the same include directories, compile definitions and required compiler
# options as mess_driver if they are not inherited from its linked libraries.
```

Run the probe before accepting kinetics. Success includes both
`cholesky_pass=1` and `solve_pass=1`. The probe supports little-endian LP64 and
ILP64 interfaces; it explicitly declines other storage layouts. It is a
library-boundary check, not a test of a chemical model or rate convergence.

In the local 2026-11-02 validation, the downloaded MESS source declared
`int64_t` but the linked Apple Accelerate interface was LP64. The old build
returned `cholesky_info=4294967296` for a valid 2-by-2 matrix. A separate build
using `int32_t` for that library passed both probes. Earlier solver outputs
were retained as historical evidence and superseded by checks with the
correct interface. Multiprecision eigensolver tests alone did not detect this
ordinary LAPACK boundary error; both checks are required. Increasing numerical
precision cannot repair an incompatible library interface.

MESS precision/reduction settings are independent checks. In particular,
`ChemicalEigenvalueMax < 1` constrains relaxation-subspace projection; it is
not an eigenvalue ratio. Changing it can change the number and grouping of
kinetic populations. A table becoming nonnegative after a change does not
alone validate that change. Preserve inputs and logs and compare precision,
energy grains, selected populations and significant rates. See the native
[MESS manual](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf).
