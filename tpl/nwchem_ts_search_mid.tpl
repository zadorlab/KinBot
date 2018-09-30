echo

title "Mid template for TS search"

restart _name_

scratch_dir _scratch_
permanent_dir ./perm

_constraints_

charge _charge_

basis
* library _basis_
end

scf
 _scfmultiplicity_
 print low
end

mp2
 print low
end

dft
 _odft_
 mult _multiplicity_
 xc _method_
 print low
end

driver
 _convergence_ # loose, default, or tight
 maxiter _maxiter_ # default is 20, for KinBot it is 100
 clear
 inhess 1
 xyz _name_
 print low
end

freq
 animate
end

task dft optimize

