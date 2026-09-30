"""Attach current calculation records to synthetic conformer test inputs."""
from kinbot.conformer_records import ConformerRecord, retain


def record_conformers(species):
    records = []
    for offset, index in enumerate(species.conformer_index):
        if index < 0:
            continue
        records.append(ConformerRecord(
            member_id=f'{species.name}:conformer:{index}', index=index,
            source_job=f'conf/{species.name}_{index:04d}', status='valid',
            geometry=tuple(tuple(map(float, xyz)) for xyz in species.conformer_geom[offset]),
            zero_energy_hartree=float(species.conformer_zeroenergy[offset]),
            frequencies_cm1=tuple(map(float, species.conformer_freq[offset]))))
    retain(species, records, [record.index for record in records])
