***,{name}
memory,1600,M
orient
geomtyp=xyz
geometry={{
{natom}
{name}
{geom}
}}
basis=cc-pvdz-f12
{{rhf;wf,{nelectron},{symm},{spin},{charge}}}
{{uccsd(t)-f12b,scale_trip=1}}

! mydza is retained as a compatibility alias; both keys are explicit F12b.
mydza = energy
mydzb = energy

basis=cc-pvtz-f12
{{rhf;wf,{nelectron},{symm},{spin},{charge}}}
{{uccsd(t)-f12b,scale_trip=1}}

! mytza is retained as a compatibility alias; both keys are explicit F12b.
mytza = energy
mytzb = energy
---
