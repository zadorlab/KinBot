# Completed compact-union panel, 25 October 2026, 17:25:59 PDT

All ten jobs completed with native completion checks and zero recorded startup,
manual-startup or allocation retries. `complete_review_required` is the correct
campaign state: completion does not publish or approve the orbitals. The supplied
[report](union_character_evidence.txt) and [provenance](union_report_provenance.json)
are retained verbatim. The new WFU/XML files have not yet been transferred and
independently verified locally; the report includes their remote hashes.

## CH2 + HSS C1–S2 homolysis: select CAS(5,4)

Six probes cover 2, 2.5 and 30 Å in avdz-f12 and avtz-f12. Retain natural
**19,20,21,22** from the basis-specific **2 Å** checkpoints, record 2140.2.
The canonical source at 30 Å was **19,20,21,22**: occupied S2 valence donor 19,
HSS radical and the two CH2 frontier orbitals 20–22. Source and optimized
ordinal labels refer to different orbital sets; their identical numbers are
not evidence of unchanged individual orbitals.

At 30 Å the optimized occupations are 2/1/1/1, with the occupied sulfur donor,
HSS radical and both carbon frontier directions. At 2.5 Å, 19 is predominantly
the occupied donor, 20/22 mix carbon and sulfur, and 21 is predominantly carbon
p. At 2 Å, 19/22 contain the C–S bonding/antibonding directions, 20 is sulfur
p-rich, and 21 is carbon p-rich. This retains the carbon frontier contribution
missing from the earlier near CAS(3,3) result without including other valence
bonds or an arbitrary virtual window. AO coefficients establish qualitative
character; their squared values are not AO populations.

| Active occupations 19–22 | DZ | TZ |
| --- | --- | --- |
| 2 Å | 1.96228, 1.49558, 1.00285, 0.53930 | 1.96279, 1.49549, 1.00287, 0.53885 |
| 2.5 Å | 1.98484, 1.38372, 1.00044, 0.63100 | 1.98669, 1.38287, 1.00039, 0.63004 |
| 30 Å | 2, 1, 1, 1 | 2, 1, 1, 1 |

D1 at 2 Å is -834.66814805/-834.68328100 Eh (DZ/TZ), compared with
-834.64763644/-834.66296618 Eh in CAS(3,3). Both spaces have the same reported
30 Å energy. The additional direction therefore matters for the inner bonding
region. This does not establish CASPT2 convergence or a production barrier.
One Q and one D retain total weights 0.5 each for orbital optimization; only D1
is reactive flux. The fractional averaged occupations do not by themselves
measure correlation in D1. Production retains its orbital, state, geometry and
scan review checks.

## C(SH2) + SH outer exit: the proposed CAS(9,6) enlargement did not settle the question

Four probes cover the actual L2 product complex and 30 Å in both bases.
At separation the active SH p pair is 20/21, occupation 1.5 each. At the L2
endpoint one SH p direction is again inactive, now **16**, occupation 2.0.
Active **17** is instead an occupied C/S valence orbital (occupation
1.99840/1.99836), while active 20 is carbon p-rich. The near D1 energy is only
0.554/0.565 kcal/mol below the preceding lower CAS(7,5) solution; lower energy
alone does not establish the intended fragment/state correspondence.

| Near result | DZ | TZ |
| --- | --- | --- |
| Active 17–22 occupations | 1.99840, 1.99393, 1.96762, 1.48980, 1.31127, 0.23898 | 1.99836, 1.99386, 1.96715, 1.48968, 1.31258, 0.23836 |
| Inactive SH-p orbital | 16; major S3 p coefficients -0.75521, 0.40051, -0.36189 | 16; -0.78895, 0.42605, -0.38136 |
| D1 / Eh | -834.47705480 | -834.49157188 |
| D2 / Eh | -834.44198988 | -834.45763260 |

The previous same-geometry comparison identified a six-dimensional union of
the two optimized CAS(7,5) spans. Adding far *canonical* orbital 17 was a trial
for obtaining that union, not an identity-preserving transport of the actual
near union. The present optimization did not demonstrate that it retained the
intended six-dimensional union. This is a failed selection hypothesis, not a
Molpro execution failure. It does not justify automatically increasing CAS
size again or changing root weights.

The supplied full XML comparison now confirms that difference. Against old
near natural 17–22, the new active 17–22 has minimum overlap singular value
**0.121844 (DZ), 0.122409 (TZ)**. Replacing new active 17 with inactive SH 16
raises those minima to **0.990473/0.990240**. New 18–22 individually reproduce
old 18–22 with squared overlaps above 0.9959; new 17 is instead about 98%
old occupied 12–14 by summed squared MO overlap. The dominant CI strings begin
with a doubly occupied spectator. These are complete same-geometry MO overlaps,
with reconstruction errors below 5.6e-15, not truncated AO population estimates.

Adding far canonical 17 therefore did not initialize the identified natural
union. The next four checks retain CAS(9,6), two equal doublets and both endpoint
geometries, but start directly from the old near **natural 17–22**. The near
optimized result then seeds 30 Å separately in each basis. The
[four-job recipe and handoff](../r7_natural_union_20261025/README.md) preserve
all prior attempts and require native character review; they do not freeze the
orbitals, approve the model or assume this initialization will survive optimization.
The [comparison evidence](../r7_natural_union_20261025/source_evidence.json)
records the supplied numbers, XML hashes and independently archived source hashes.

## Publication preparation and controller refresh

The eleven supported compact choices are materialized in [selections](selections).
They include the ten decisions documented in the [previous audit](README.md#disposition-of-the-other-ten-requests)
plus the CH2/HSS CAS(5,4) selection above. The three HCSSH/H CAS(1,1) choices
retain the earlier explicit limitation: their formal single-reference minimum
is selected from two completed CAS(3,2) comparisons; no separate CAS(1,1) probe
or correlated-potential convergence test was run.

[publication_plan.json](publication_plan.json) is **held** and cannot publish
partially. C(SH2)/SH is deliberately absent from its selected models until the
remaining native review is resolved. These files do not modify live inputs.
Each prepared seed is tied to the reviewed WFU hash; earlier transferred seeds
also have the independently verified native-output hash. The helper now accepts
an explicit `source_model` title for multi-candidate campaigns, selects that
candidate's configuration and manifest consistently for preview, proposal and
batch publication, and refuses ambiguous candidates, changed reviewed hashes,
native validation errors and held plans. Unselected candidates are not silently
chosen by list position.

The final all-channel handoff must still include the controller restart to load
the VRC endpoint-discovery fix, reconcile the five older published models, and
publish the twelve new ones through the owning drivers. Preserve live jobs,
limits, exclusions, operator rejections and existing sampling campaigns. No
earlier PES search, IRC, L2, ANL or CBH rerun is indicated by the synthetic-bond
bug. Verify actual scan dispatch and sampling adoption after refresh; a receipt
alone is insufficient. No production refresh has been executed in this audit.

## Latest production snapshot

The supplied [17:25 report](production_monitor_172559.txt) retains cap23,
policy revision14, 22 allocations (CH3S12, CH3SS10), zero diagnostic allocations
and 16 applied live controls. CH3S has completed L1/L2 exploration and remains
in postprocessing: four retained VRC state scans have 48 accepted entries each,
four RotdPy campaigns are running, and CBH-ANL tasks continue. CH3SS still has
five wells exploring and five in VRC; the final network is not known.

The previously reviewed CH2SHS S–S scan advanced from 2.4 Å running to three
accepted points and 2.6 Å running. Its C–S scan is still at 2.2 Å running. The
initial CH3SS C–S scan is resolved and RotdPy is running; its S–S scan is at
6.4 Å running. RotdPy submissions have moved to later surfaces (including
Surface_75 and Surface_85 for initial CH3S); the `0/116` counter measures final
surface files, not all intermediate sampling work. There is no completed MESS
solve. Historical failed TS validation and inactive diagnostic errors do not
indicate a new retry requirement. The known endpoint-graph hold remains in the
old controller, as expected before the requested final restart.
