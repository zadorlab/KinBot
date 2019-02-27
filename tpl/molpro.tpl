***,{name}
memory,1600,M
geomtyp=xyz
geometry={{
{natom}
{name}
{geom}
}}
{{uhf;wf,{nelectron},1,{spin},{charge}}}

basis=cc-pvdz-f12
rhf
CCSD(T)-F12

myena(1) = energy(1)
myenb(1) = energy(2)

basis=cc-pvtz-f12
rhf
CCSD(T)-F12

myenergy(2) = energy(1)
myenb(2) = energy(2)
---

