"""MESS models for two-reactant abstraction, from saved calculation records.

The ordinary MESS solver handles an accepted entrance complex. Argonne's
abstraction program handles direct bimolecular TST without an entrance well.
Both use MESS's own translation, partition sums and tunneling calculation.
"""
from collections import Counter, defaultdict
from itertools import product
import copy
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import networkx as nx
from ase import Atoms

from kinbot import constants, frequencies, geometry
from kinbot.mess import MESS, finalize_mc_mess
from kinbot.mess_mirrors import _reflect

logger = logging.getLogger('KinBot')


def _finite(value, label):
    if value is None or not np.isfinite(float(value)):
        raise ValueError(f'Missing finite {label}.')
    return float(value)


def union(blocks):
    if not blocks:
        raise ValueError('No valid stationary-state model.')
    return blocks[0] if len(blocks) == 1 else 'Union\n' + '\n'.join(blocks) + '\nEnd ! Union\n'


class RecordRenderer(MESS):
    """Reuse shared geometry, frequency and RRHO templates without QC reads."""

    @staticmethod
    def point(data):
        return SimpleNamespace(name=data.get('name', 'saved calculation'),
            atom=data['atom'], natom=len(data['atom']), geom=np.asarray(data['geometry']),
            bond=np.asarray(data['bond']), mult=data['mult'])

    @staticmethod
    def validate_record(data):
        for key, expected in (('energy_unit', 'hartree'), ('geometry_unit', 'angstrom'),
                              ('frequency_unit', 'cm^-1')):
            if data.get(key) != expected:
                raise ValueError(f'{key} must be {expected} in the saved calculation.')

    @staticmethod
    def members(data):
        RecordRenderer.validate_record(data)
        thermo = data['thermochemistry']
        if thermo.get('conformer_representation') == 'MC-RRHO':
            records = [r for r in thermo['conformer_inventory']
                       if r['retained'] and r['status'] == 'valid']
            if not records and thermo.get('rrho_representative_counting'):
                records = [thermo['rrho_representative_counting']]
            if not records:
                raise ValueError('MC calculation has no retained conformer properties.')
            return records
        optical = thermo['optical_counting']
        weight = _finite(optical.get('remaining_multiplier'), 'optical multiplier')
        decisions = (thermo.get('frequency_projection') or {}).get('rotors', [])
        accepted = {r['rotor_index'] for r in decisions if r['projected']}
        usable = {r['index'] for r in thermo['hir']['rotors'] if r['usable']}
        raw, reduced = (thermo['raw_harmonic_frequencies_cm-1'],
                        thermo['thermochemical_frequencies_cm-1'])
        if not accepted <= usable or len(raw) - len(reduced) != len(accepted):
            raise ValueError('Saved projected frequencies and active rotors disagree; recalculate the selected model.')
        return [dict(member_id=data.get('name', 'selected'), geometry=data['geometry'],
            zero_energy_hartree=data['energy'] + data['zpe'],
            frequencies_cm1=thermo['thermochemical_frequencies_cm-1'],
            raw_frequencies=thermo['raw_harmonic_frequencies_cm-1'],
            sigma_ext=thermo['sigma_ext'], remaining_optical_weight=weight,
            optical_evidence=optical)]

    def ground(self, data):
        return min(_finite(r['zero_energy_hartree'], 'E + ZPE') for r in self.members(data))

    def rotors(self, data, *, geometry=False):
        thermo = data['thermochemistry']
        if thermo.get('conformer_representation') == 'MC-RRHO':
            return ''
        blocks = []
        projection = thermo.get('frequency_projection') or {}
        projected = {r['rotor_index']: r['projected'] for r in projection.get('rotors', [])}
        for rotor in thermo['hir']['rotors']:
            if not rotor['usable'] or not projected.get(rotor['index'], False):
                blocks.append(f"! Rotor {rotor['axis']} retained harmonically: {rotor['exclusion_reason']}")
                continue
            sigma = rotor['sigma_int']
            points = rotor['points']
            count = max(3, int(np.ceil(len(points) / sigma)))
            energies = [p['relative_energy_kcal_mol'] if p['status'] == 'successful'
                        else p['fitted_relative_energy_kcal_mol'] for p in points]
            energies = [_finite(e, 'rotor potential') for e in energies]
            free = max(energies) < self.par['free_rotor_thrs']
            if len(points) % sigma or len(points) // sigma < 3:
                coefficients = np.asarray(rotor['fourier']['coefficients'])
                n = len(coefficients) // 2
                k = np.arange(1, n + 1)
                energies = [float(np.dot(coefficients[:n], 1 - np.cos(k * a)) +
                                  np.dot(coefficients[n:], np.sin(k * a)))
                            for a in np.arange(count) * 2 * np.pi / (sigma * count)]
            values = dict(group=' '.join(str(i + 1) for i in rotor['partition'][0][1:]),
                axis=' '.join(str(i + 1) for i in rotor['axis']), rotorsymm=sigma,
                nrotorpot=count, rotorpot=' '.join(f'{e:.2f}' for e in energies[:count]),
                natom=len(data['atom']), geom=self.make_geom(data['geometry'], data['atom']))
            template = (self.freerotortpl if free else
                        self.hinderedrotorgeomtpl if geometry else self.hinderedrotortpl)
            blocks.append(template.format(**values))
        if not blocks:
            blocks.append('! Harmonic model; no usable hindered rotors in the saved calculation.')
        return '\n'.join(blocks)

    @staticmethod
    def electronic(data):
        # Preserve the shared OH fine-structure convention.
        return '2\n  0. 2\n  140. 2' if str(data.get('chemid')) == '170170000000000000002' else f"1\n  0. {data['mult']}"

    def rrho(self, data, reference, *, depths=None, as_blocks=False):
        if data.get('reflected_model'):
            original = dict(data)
            original.pop('reflected_model')
            blocks = ['! Global mirror: complete calculated model reused.\n' + _reflect(block)
                      for block in self.rrho(original, reference, depths=depths, as_blocks=True)]
            return blocks if as_blocks else union(blocks)
        if len(data['atom']) == 1:
            block = f"Atom\n Name {data['atom'][0]}\n ElectronicLevels[1/cm] {self.electronic(data)}\nEnd ! Atom\n"
            return [block] if as_blocks else block
        point = self.point(data)
        blocks = []
        for record in self.members(data):
            sigma = _finite(record['sigma_ext'], 'external symmetry number')
            weight = _finite(record['remaining_optical_weight'], 'optical weight')
            if min(sigma, weight) <= 0:
                raise ValueError('Symmetry divisors must be positive.')
            energy = _finite(record['zero_energy_hartree'], 'E + ZPE')
            modes = frequencies.thermochemical_frequencies(
                record['frequencies_cm1'], depths is not None, self.par['imagfreq_threshold'])
            written_modes = modes[int(depths is not None):]
            if (not np.all(np.isfinite(modes)) or any(f <= 0 for f in written_modes)
                    or (depths is not None and (not len(modes) or modes[0] >= 0))):
                raise ValueError(f"{record['member_id']}: invalid harmonic frequencies for this stationary point.")
            tunnel = ''
            if depths is not None:
                left, right = [(energy - endpoint) * constants.AUtoKCAL for endpoint in depths]
                if min(left, right) > 0:
                    tunnel = self.tunneltpl.format(imfreq=abs(modes[0]),
                        cutoff=min(left, right), welldepth1=left, welldepth2=right)
                else:
                    tunnel = '! Eckart omitted: TS is at or below an endpoint reference.'
            geometry_modes = record.get('raw_frequencies', record['frequencies_cm1'])
            expected_counts = {1} if point.natom == 2 else {3 * point.natom - 6, 3 * point.natom - 5}
            if len(geometry_modes) not in expected_counts:
                raise ValueError(f"{record['member_id']}: incomplete harmonic frequency list.")
            text = self.rrhotpl.format(natom=point.natom,
                geom=self.rotor_geom(point, record['geometry'], geometry_modes),
                core=self.corerrtpl.format(symm=sigma / weight,
                    rotconst=self.rotor_core_line(point, record['geometry'], geometry_modes)),
                nfreq=len(modes) - int(depths is not None),
                freq=self.make_freq(modes, 1., int(depths is not None)),
                rotors=self.rotors(data), nelec=1, mult=point.mult,
                zeroenergy=(energy - reference) * constants.AUtoKCAL,
                tunneling=tunnel, shift='saved electronic energy + harmonic ZPE')
            text = text.replace(f'ElectronicLevels[1/cm]        1\n      0.    {point.mult}',
                                'ElectronicLevels[1/cm]        ' + self.electronic(data))
            blocks.append(self._counting_comment(record.get('optical_evidence'), record['member_id']) + text)
        return blocks if as_blocks else union(blocks)

    def bimolecular(self, node, reference):
        fragments = node['fragment_data']
        if len(fragments) != 2:
            raise ValueError('H-abstraction requires two separated fragments.')
        blocks = [f'Fragment f_{i}\n' + self.rrho(data, self.ground(data))
                  for i, data in enumerate(fragments)]
        return self.bimoltpl.format(chemids=node['id'], fragments='\n'.join(blocks),
            ground_energy=(sum(self.ground(d) for d in fragments) - reference) * constants.AUtoKCAL)

    def pst(self, fragments, reference, multiplicity, *, as_blocks=False):
        """One capture surface, summed over separated-fragment MC states."""
        blocks = []
        for members in product(*(self.members(d) for d in fragments)):
            geometries, modes, rotors = [], [], []
            explicit_constants, rotation_constants, masses = False, [], []
            divisor, energy = 1., 0.
            for data, member in zip(fragments, members):
                raw_modes = member.get('raw_frequencies', member['frequencies_cm1'])
                point = self.point(data)
                rendered = self.rotor_geom(point, member['geometry'], raw_modes)
                geometries.append(f"FragmentGeometry[angstrom] {len(data['atom'])}\n" + rendered)
                forced = bool(self.rotor_core_line(point, member['geometry'], raw_modes))
                explicit_constants |= forced
                masses.append(sum(constants.exact_mass[a] for a in data['atom']))
                if len(data['atom']) > 1:
                    coords = np.array([[float(x) for x in line.split()[1:]]
                                       for line in rendered.splitlines()
                                       if len(line.split()) == 4 and line.split()[0] in data['atom']])
                    moments = geometry.get_moments_of_inertia(coords, data['atom'])[0]
                    rotation_constants.extend(([frequencies.ROTATIONAL_CONSTANT_FACTOR / moments[1]] * 2)
                        if not forced and moments[0] / moments[1] < frequencies.MESS_LINEAR_MOMENT_RATIO
                        else [frequencies.ROTATIONAL_CONSTANT_FACTOR / m for m in moments])
                divisor *= member['sigma_ext'] / member['remaining_optical_weight']
                energy += member['zero_energy_hartree']
                modes.extend(frequencies.thermochemical_frequencies(member['frequencies_cm1'], False,
                                                                    self.par['imagfreq_threshold']))
                rotors.append(self.rotors(data, geometry=True))
            if explicit_constants:
                geometries = ['FragmentMasses[amu] ' + ' '.join(map(str, masses)) +
                              '\nRotationalConstants[1/cm] ' + ' '.join(map(str, rotation_constants)), '']
            stoich = ''.join(f'{atom}{count}' for atom, count in sorted(
                Counter(sum((d['atom'] for d in fragments), [])).items()))
            blocks.append(self.pstrrhotpl.format(stoich=stoich, frag1=geometries[0], frag2=geometries[1],
                pstsymm=divisor, nfreq=len(modes), freq=self.make_freq(modes, 1., 0),
                hinderedrotor='\n'.join(rotors), mult=multiplicity,
                zeroenergy=(energy - reference) * constants.AUtoKCAL))
        return blocks if as_blocks else union(blocks)


class BimolecularMESS(RecordRenderer):
    def __init__(self, par, records):
        super().__init__(par, None)
        if records.get('schema_version') != 3 or records.get('energy_unit') != 'kcal/mol':
            raise ValueError('Bimolecular MESS requires version 3 calculation records with explicit units.')
        if 'rate_models' not in records:
            raise ValueError('These saved observations predate bimolecular rate-model selection. '
                             'Restart KinBot in the calculation directory to write current records.')
        self.records = records
        self.nodes = {n['id']: n for n in records['nodes']}
        self.models = copy.deepcopy(records['rate_models'])
        self.complexes = {n['id']: n for n in self.models['complexes']}
        reactants = self.nodes[records['reference_node']]['fragment_data']
        self.reference = sum(self.ground(d) for d in reactants)
        self.mess_jobs = []
        self.warnings = list(self.models['warnings'])
        self._reuse_specified_mirror_models()
        self._complete_specified_products()

    @staticmethod
    def _identity_pair(node, mirror=False):
        return tuple(sorted((d['chemid'], d['charge'], d['mult'],
                             d['stereochemical_identity']['mirror_id' if mirror else 'id'])
                            for d in node['fragment_data']))

    @staticmethod
    def _mirror_point(point):
        result = copy.deepcopy(point)
        result['reflected_model'] = not result.get('reflected_model', False)
        identity = result['stereochemical_identity']
        if 'id' in identity and 'mirror_id' in identity:
            identity['id'], identity['mirror_id'] = identity['mirror_id'], identity['id']
        return result

    def _reuse_specified_mirror_models(self):
        """Reuse one complete endpoint calculation, as the shared writer does."""
        if self.par.get('optical_population', 'specified') != 'specified':
            return
        populations = {self._identity_pair(n): key for key, n in self.nodes.items()
                       if n['kind'] == 'bimolecular'}
        visited = set()
        ordered = [self.records['reference_node']] + sorted(self.nodes)
        for source_id in ordered:
            source = self.nodes[source_id]
            if source_id in visited or source['kind'] != 'bimolecular':
                continue
            target_id = populations.get(self._identity_pair(source, mirror=True))
            if target_id is None or target_id == source_id:
                continue
            visited.update((source_id, target_id))
            self.nodes[target_id] = dict(self.nodes[target_id],
                fragment_data=[self._mirror_point(d) for d in source['fragment_data']],
                mirrored_from=source_id)
            for role in ('rc', 'pc'):
                original = next((c for c in self.complexes.values()
                                 if c['asymptote'] == source_id and c['role'] == role), None)
                target = next((c for c in self.complexes.values()
                               if c['asymptote'] == target_id and c['role'] == role), None)
                if original and target:
                    self.complexes[target['id']] = dict(target,
                        stationary_point=self._mirror_point(original['stationary_point']))

    def _complete_specified_products(self):
        """An achiral entrance also supplies the reflected specified product.

        A specified chiral entrance does not authorize its other enantiomer.
        Racemic endpoints have already been grouped before records are saved.
        """
        if self.par.get('optical_population', 'specified') != 'specified':
            return
        entrance = self.nodes[self.records['reference_node']]
        if self._identity_pair(entrance) != self._identity_pair(entrance, mirror=True):
            return
        populations = {self._identity_pair(n): key for key, n in self.nodes.items()
                       if n['kind'] == 'bimolecular'}
        for path in list(self.models['pathways']):
            products = self.nodes[path['product']]
            if products['kind'] != 'bimolecular':
                continue
            own, mirror = self._identity_pair(products), self._identity_pair(products, mirror=True)
            if own == mirror:
                continue
            target = populations.get(mirror)
            if target is None:
                target = path['product'] + '_mirror'
                self.nodes[target] = dict(products, id=target,
                    fragment_data=[self._mirror_point(d) for d in products['fragment_data']])
                populations[mirror] = target
            if any(p['product'] == target and p['path_id'] == path['path_id']
                   for p in self.models['pathways']):
                continue
            derived = dict(path, id=path['id'] + '_mirror', product=target,
                           stationary_point=self._mirror_point(path['stationary_point']))
            if path['product_complex']:
                source = self.complexes[path['product_complex']]
                key = source['id'] + '_mirror'
                self.complexes[key] = dict(source, id=key, asymptote=target,
                    stationary_point=self._mirror_point(source['stationary_point']))
                derived['product_complex'] = key
            self.models['pathways'].append(derived)

    def _depths(self, path):
        values = []
        for side in ('reactant', 'product'):
            complex_id = path[f'{side}_complex']
            if complex_id:
                # Keep both complex energies in the tunneling model even when
                # the two-TS calculation omits the complex wells.
                data = self.complexes[complex_id]['stationary_point']
                values.append(data['energy'] + data['zpe'])
            else:
                values.append(sum(d['energy'] + d['zpe'] for d in self.nodes[path[side]]['fragment_data']))
        return values

    def _outer_choice(self, path):
        """Keep the PST surface of the higher separated-fragment E + ZPE.

        The stable endpoint key breaks exact ties independently of direction.
        Both complex energies remain available to the inner Eckart model.
        An accepted complex on either side enables the two-TS treatment, but
        the chosen PST needs only its separated fragments, not a bound well.
        """
        if not any(path[f'{side}_complex'] for side in ('reactant', 'product')):
            return None
        choices = []
        for side in ('reactant', 'product'):
            identifier = path[f'{side}_complex']
            separated = sum(self.ground(d) for d in self.nodes[path[side]]['fragment_data'])
            depth = None
            if identifier:
                point = self.complexes[identifier]['stationary_point']
                depth = (separated - point['energy'] - point['zpe']) * constants.AUtoKCAL
            choices.append(dict(side=side, complex=identifier, endpoint=path[side],
                                asymptote_energy_kcal_mol=(separated - self.reference) * constants.AUtoKCAL,
                                binding_depth_kcal_mol=depth))
        return min(choices, key=lambda c: (-round(c['asymptote_energy_kcal_mol'], 9), c['endpoint']))

    def _header(self, reactant, mass):
        p = self.par
        return self.headertpl.format(LevelOfTheory='saved final calculations; see bimol_network.json',
            TemperatureList=' '.join(map(str, p['TemperatureList'])),
            PressureList=' '.join(map(str, p['PressureList'])),
            **{key: p[key] for key in ('EnergyStepOverTemperature', 'ExcessEnergyOverTemperature',
                'ModelEnergyLimit', 'CalculationMethod', 'ChemicalEigenvalueMax',
                'EnergyRelaxationFactor', 'EnergyRelaxationPower', 'EnergyRelaxationExponentCutoff')},
            Reactant=reactant, e_coll=constants.epsilon[p['collider']],
            s_coll=constants.sigma[p['collider']], m_coll=constants.mass[p['collider']],
            e_well=p['epsilon'], s_well=p['sigma'], m_well=mass)

    def write_input(self, qc=None):
        Path('me').mkdir(exist_ok=True)
        self.mess_jobs = []
        self.warnings = list(self.models['warnings'])
        if self.par.get('uq_n', 1) != 1:
            self.warnings.append('Bimolecular uncertainty sampling is not implemented; writing the nominal model only.')
        for key, model in self.complexes.items():
            try:
                data = model['stationary_point']
                # The no-well model needs binding and tunneling energies only.
                # It must not depend on complex frequencies or symmetry.
                if not self.par['me_skip_vdW']:
                    self.rrho(data, self.reference)
                elif data.get('energy_unit') != 'hartree':
                    raise ValueError('complex energy_unit must be hartree')
                fragments = self.nodes[model['asymptote']]['fragment_data']
                complex_energy = _finite(data['energy'], 'complex energy') + _finite(data['zpe'], 'complex ZPE')
                depth = (sum(self.ground(d) for d in fragments) - complex_energy) * constants.AUtoKCAL
                if depth <= self.par['vdW_detection']:
                    raise ValueError(f'complex binding depth {depth:.6f} kcal/mol does not exceed {self.par["vdW_detection"]}')
            except ValueError as error:
                self.warnings.append(f'{key}: optional complex omitted: {error}. Saved observations retained.')
                for path in self.models['pathways']:
                    for side in ('reactant_complex', 'product_complex'):
                        if path[side] == key:
                            path[side] = None
        direct, two_ts, coupled = defaultdict(list), defaultdict(list), []
        for path in self.models['pathways']:
            try:
                data = path['stationary_point']
                for side in ('reactant', 'product'):
                    self.bimolecular(self.nodes[path[side]], self.reference)
                    pair = self._identity_pair(self.nodes[path[side]])
                    if len(pair) == 2 and pair[0] == pair[1]:
                        raise ValueError('Identical-reactant event and molecular-loss conventions are not yet supported.')
                    complex_id = path[f'{side}_complex']
                    if (complex_id and
                            any(d['charge'] for d in self.nodes[path[side]]['fragment_data'])):
                        raise ValueError('The neutral -0.6/R^6 capture approximation cannot describe charged fragments.')
                # Validate all properties before adding a pathway to a model.
                self.rrho(data, self.reference, depths=self._depths(path))
                choice = self._outer_choice(path) if self.par['me_skip_vdW'] else None
                if choice:
                    if any(d['charge'] for d in self.nodes[choice['endpoint']]['fragment_data']):
                        raise ValueError('The neutral -0.6/R^6 capture approximation cannot describe charged fragments.')
                    path['outer_pst_selection'] = choice
                    two_ts[choice['endpoint'], data['mult']].append(path)
                elif path['reactant_complex'] or path['product_complex']:
                    coupled.append(path)
                elif self.ground(data) > max(sum(self.ground(d) for d in
                        self.nodes[path[side]]['fragment_data']) for side in ('reactant', 'product')):
                    direct[path['reactant'], path['product']].append(path)
                else:
                    raise ValueError('TS is at or below a separated endpoint without a bound-complex model; an entrance model or corrected pathway is required.')
            except ValueError as error:
                self.warnings.append(f"{path['id']}: {error}")
        for endpoints, paths in sorted(direct.items()):
            stem = f'mess_direct_{len(self.mess_jobs):04d}'
            text = (f"TemperatureList[K] {' '.join(map(str, self.par['TemperatureList']))}\n"
                    f"ModelEnergyLimit[kcal/mol] {self.par['ModelEnergyLimit']}\n"
                    'AtomDistanceMin[angstrom] 0.5\n'
                    f'LogOutput {stem}.log\nRateOutput {stem}.out\n')
            text += ''.join(self.bimolecular(dict(self.nodes[key], id=name), self.reference)
                            for name, key in zip(('R', 'P'), endpoints))
            text += 'Barrier abstraction\n' + union([self.rrho(p['stationary_point'], self.reference,
                depths=self._depths(p)) for p in paths])
            self._save(stem, text, 'direct_tst', paths)
            self.mess_jobs[-1]['endpoints'] = dict(zip(('R', 'P'), endpoints))
        for paths in two_ts.values():
            self._write_two_ts(paths)
        # Only explicit-well models reach this list. The two-TS models above
        # omit both wells. Shared separated products do not connect wells.
        graph = nx.Graph()
        for path in coupled:
            wells = [path[side] for side in ('reactant_complex', 'product_complex')
                     if path[side]]
            graph.add_nodes_from(wells)
            if len(wells) == 2:
                graph.add_edge(*wells)
        for component in nx.connected_components(graph):
            paths = [p for p in coupled if p['reactant_complex'] in component or
                     p['product_complex'] in component]
            try:
                self._write_coupled(paths)
            except ValueError as error:
                self.warnings.append(f'Complex network: {error}')
        for warning in self.warnings:
            logger.warning('Bimolecular MESS: %s', warning)
        if not self.mess_jobs:
            logger.warning('No supported bimolecular rate model is ready. No MESS job will be run.')
        Path('me/bimol_models.json').write_text(json.dumps(dict(jobs=self.mess_jobs,
            warnings=self.warnings, complete=not self.warnings,
            rate_units='cm^3 molecule^-1 s^-1 for bimolecular rates; s^-1 for well loss',
            pst_potential='-0.6/R^6 in atomic units',
            energy_convention='electronic energy + saved harmonic ZPE; no second rotor ground correction',
            tunneling_convention='each TS conformer uses its own imaginary frequency and energy; Eckart depths use selected-parent endpoint energies',
            site_convention='stored rotational, internal-rotor and optical weights; no extra raw site multiplier'),
            indent=2) + '\n')

    def _write_two_ts(self, paths):
        from kinbot.bimol_two_ts import write_model
        stem = f'mess_two_ts_{len(self.mess_jobs):04d}'
        choice = paths[0]['outer_pst_selection']
        fragments = self.nodes[choice['endpoint']]['fragment_data']
        # Set the selected asymptote exactly at zero, on MESS's state grid.
        # A fractional grid offset can give negative interpolated PST counts
        # immediately above threshold. Shift every model by the same amount.
        reference = sum(self.ground(d) for d in fragments)
        outer = self.pst(fragments, reference, paths[0]['stationary_point']['mult'], as_blocks=True)
        grouped = defaultdict(list)
        for path in paths:
            grouped[path['reactant'], path['product']].append(path)
        branches = []
        for endpoints, selected in sorted(grouped.items()):
            blocks = [block for path in selected for block in self.rrho(path['stationary_point'],
                reference, depths=self._depths(path), as_blocks=True)]
            fragments_text = ''.join(self.bimolecular(dict(self.nodes[key], id=name), reference)
                                     for name, key in zip(('R', 'P'), endpoints))
            threshold = max(sum(self.ground(d) for d in self.nodes[key]['fragment_data'])
                            for key in endpoints) - reference
            branches.append(dict(endpoints=dict(zip(('R', 'P'), endpoints)),
                pathways=[p['id'] for p in selected], blocks=blocks, fragments=fragments_text,
                threshold_cm1=threshold * constants.AUtoCM,
                outer_pst_selection=[p['outer_pst_selection'] for p in selected]))
        lower = min(min(*self._depths(p), self.ground(p['stationary_point']))
                    for p in paths) - reference
        mass = sum(Atoms(d['atom']).get_masses().sum() for d in fragments)
        note = ('VDW wells omitted. Statistical two-transition-state approximation.\n'
                f"Outer PST selected by the higher separated-fragment E + ZPE: {choice['endpoint']}; "
                f"{choice['asymptote_energy_kcal_mol']:.8f} kcal/mol relative to the reactant reference.\n"
                'All energies in this calculation use that selected asymptote as zero.\n'
                'The other outer PST surface is omitted. Both complex energies set Eckart depths.\n'
                'Forward and reverse rates use the same selected surfaces. No collisional stabilization.')
        write_model(Path('me'), stem, branches, outer, mass, lower * constants.AUtoCM,
                    self.par['ModelEnergyLimit'], self.par['TemperatureList'], note,
                    0., reference)
        self.mess_jobs.append(dict(stem=stem, kind='two_ts', status='ready',
            energy_reference_hartree=reference,
            pathways=[p['id'] for p in paths], outer_pst_endpoint=choice['endpoint'],
            approximation='higher separated-fragment E + ZPE selects one outer PST surface',
            channels=[{k: v for k, v in b.items() if k not in ('blocks', 'fragments')}
                      | {'rate_output': f'{stem}_channel_{i:04d}.out'} for i, b in enumerate(branches)]))

    def _write_coupled(self, paths):
        reactant = self.records['reference_node']
        nodes, complexes = {reactant}, set()
        barriers = defaultdict(list)
        accepted = []
        for path in paths:
            rc_id = path['reactant_complex']
            pc_id = path['product_complex']
            right = pc_id or path['product']
            floor = max(self.ground(self.complexes[c]['stationary_point'])
                        for c in (rc_id, pc_id) if c)
            if self.ground(path['stationary_point']) < floor and not self.par['correct_submerged']:
                self.warnings.append(f"{path['id']}: TS lies below a connected bound well; pathway omitted.")
                continue
            # Without an entrance well, a submerged saddle still has no capture model.
            if not rc_id and self.ground(path['stationary_point']) <= self.reference:
                self.warnings.append(f"{path['id']}: submerged TS lacks a bound entrance complex; pathway omitted.")
                continue
            complexes.update(c for c in (rc_id, pc_id) if c)
            nodes.add(path['product'])
            accepted.append(path)
            barriers[rc_id or reactant, right].append(self.rrho(path['stationary_point'], self.reference,
                                                  depths=self._depths(path)))
        if not accepted:
            return
        mass = sum(Atoms(d['atom']).get_masses().sum() for d in self.nodes[reactant]['fragment_data'])
        text = 'AtomDistanceMin[angstrom] 0.5\n' + self._header(reactant, mass)
        text += ''.join(self.bimolecular(self.nodes[n], self.reference) for n in sorted(nodes))
        for identifier in sorted(complexes):
            complex_model = self.complexes[identifier]
            data = complex_model['stationary_point']
            text += f"! {complex_model['sampling']}\nWell {identifier}\nSpecies\n"
            text += self.rrho(data, self.reference) + '\nEnd ! Species\n'
            separated = complex_model['asymptote']
            if self.ground(data) >= sum(self.ground(d) for d in self.nodes[separated]['fragment_data']):
                raise ValueError(f'{identifier}: complex is not below its final separated fragments.')
            text += f'Barrier capture_{identifier} {identifier} {separated}\n'
            text += self.pst(self.nodes[separated]['fragment_data'], self.reference, data['mult'])
        for i, ((left, right), models) in enumerate(sorted(barriers.items())):
            text += f'Barrier inner_{i} {left} {right}\n' + union(models)
        text += 'End ! Model\n'
        text = finalize_mc_mess(text, correct_submerged=self.par['correct_submerged'])
        self._save(f'mess_{sum(j["kind"] == "master_equation" for j in self.mess_jobs):04d}',
                   text, 'master_equation', accepted)

    def _save(self, stem, text, kind, paths):
        text = text.replace('MicroRateOutput micro.out', f'MicroRateOutput {stem}.micro.out')
        Path(f'me/{stem}.inp').write_text(text)
        self.mess_jobs.append(dict(stem=stem, kind=kind, status='ready',
            pathways=[p['id'] for p in paths]))

    def run(self):
        import shutil
        import shlex
        import sys
        from kinbot.mess_execution import run_mess, mess_command, MESSExecutionError
        result = 0
        batches = []
        for kind in ('master_equation', 'direct_tst', 'two_ts'):
            writer = copy.copy(self)
            writer.par = dict(self.par)
            writer.mess_jobs = [job for job in self.mess_jobs if job['kind'] == kind]
            if not writer.mess_jobs:
                continue
            if kind in ('direct_tst', 'two_ts'):
                writer.par['mess_command'] = self.par.get('mess_abstraction_command', 'mess-abstraction')
                if (self.par['run_me'] and self.par['queuing'] == 'local'
                        and shutil.which(mess_command(writer.par)[0]) is None):
                    raise MESSExecutionError('Direct bimolecular TST needs the Argonne MESS abstraction executable. '
                        'Build abstraction.cc from the MESS distribution and set mess_abstraction_command '
                        'to that executable. The ordinary mess executable handles complex networks only.')
                if kind == 'two_ts':
                    writer.par['mess_command'] = shlex.join([sys.executable, '-m', 'kinbot.bimol_two_ts',
                        '--mess-command', writer.par['mess_command'],
                        '--state-command', self.par.get('mess_command', 'mess')])
            result += run_mess(writer)
            batches.append(Path('batch_me.sub').read_text())
        Path('batch_me.sub').write_text(''.join(batches))
        return result


if __name__ == '__main__':
    import argparse
    from kinbot.parameters import Parameters
    parser = argparse.ArgumentParser(description='Write/run bimolecular MESS from saved KinBot calculations.')
    parser.add_argument('input', help='KinBot JSON input with MESS settings')
    parser.add_argument('records', nargs='?', default='bimol_network.json')
    args = parser.parse_args()
    with open(args.records) as handle:
        writer = BimolecularMESS(Parameters(args.input).par, json.load(handle))
    writer.write_input()
    if writer.par['me'] == 1:
        writer.run()
