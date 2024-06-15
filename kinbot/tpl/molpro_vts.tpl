***,{fname}
{options}

text, Molpro input generated automatically by KinBot

{basis}

{geometry_block}

{methods}

{key} = energy(1)

---

RECOMMENDED TEMPLATE
This is for CASPT2
You can also create one for CC
or if you can handle the whole thing already with your script that's fine, too, using the above blocks
but the outsome should be something like this for CASPT2
***,{label}
memory,1000,M  <<<< or do it automatically >>>>
angstrom
nosym
geom={{
    {geom}
    }}

basis={basis}

{{multi
    occ,{occ}
    closed,{closed}
    wf,{nelec},1,{mult}
    }}

put, molden, {label}_{level[hl or sl]}.mld

rs2c, shift = 0.2

{key} = energy(1)
---

