***,{name}
memory,1600,M
orient
geomtyp=xyz
geometry={{
{natom}
{name}
{geom}
}}
{{uhf;wf,{nelectron},{symm},{spin},{charge}}}

basis=cc-pvdz-f12
rhf
CCSD(T)-F12

myenergy(1) = energy(1)
myenb(1) = energy(2)

---

basis=cc-pvtz-f12
rhf
CCSD(T)-F12

myenergy(2) = energy(1)
myenb(2) = energy(2)
---

