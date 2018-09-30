echo

title "IRC template for TS search"

restart _name_

scratch_dir _scratch_
permanent_dir ./perm

geometry 
symmetry c1
_geom_
end

charge _charge_

basis
* library _basis_
end

dft
 _odft_
 mult _multiplicity_
 xc _method_
 print low
end

mepgs
 maxmep 30
 maxiter 20
 inhess 2
 xyz
 evib 0.0005
 stride 0.1
 print low
 _direction_
end

task dft mepgs

