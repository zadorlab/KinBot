"""Read a pinned release of the official Active Thermochemical Tables.

ATcT identifies chemical *states*, not just formulas. A reference is selected
by its ATcT ID or by an unambiguous gas-phase structure match. The source
HTML and its digest remain part of every derived formation enthalpy.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from html import unescape
from functools import lru_cache
from pathlib import Path
import math
import os
import re
from typing import Mapping
from urllib.request import urlopen


ATCT_URL = 'https://atct.anl.gov/Thermochemical%20Data/version%20{version}/'
_VERSION = re.compile(r'1\.\d+[a-z]?\Z')
_ROW = re.compile(r'<tr\b[^>]*\bid="[^"]*\bi\d+[^" ]*[^>]*>.*?</tr>',
                  re.IGNORECASE | re.DOTALL)
_TAG = re.compile(r'<[^>]*>')
_FORMULA_ATOM = re.compile(r'([A-Z][a-z]?)(\d*)')


def _plain(value: str) -> str:
    return ' '.join(unescape(_TAG.sub(' ', value)).split())


def _field(row: str, name: str) -> str:
    match = re.search(r'<span class="' + re.escape(name) +
                      r'">(.*?)</span>', row, re.IGNORECASE | re.DOTALL)
    return _plain(match.group(1)) if match else ''


def _number(value: str, label: str) -> float | None:
    if not value:
        return None
    try:
        result = float(value)
    except ValueError as exc:
        raise ValueError(f'Invalid ATcT {label}: {value!r}.') from exc
    if not math.isfinite(result):
        raise ValueError(f'Nonfinite ATcT {label}.')
    return result


@dataclass(frozen=True)
class ATcTRecord:
    atct_id: str
    name: str
    preferred_formula: str
    phase: str
    smiles: str
    formation_0k_kj_mol: float | None
    formation_298k_kj_mol: float | None
    uncertainty_kj_mol: float | None
    exact: bool
    version: str
    source_sha256: str


@dataclass(frozen=True)
class ATcTTable:
    version: str
    source_url: str
    source_sha256: str
    records: tuple[ATcTRecord, ...]

    def by_id(self, atct_id: str) -> ATcTRecord:
        found = [item for item in self.records if item.atct_id == atct_id]
        if len(found) != 1:
            raise LookupError(f'Expected one ATcT entry for {atct_id!r}; found {len(found)}.')
        return found[0]

    def gas_by_smiles(self, smiles: str, *, formula: Mapping[str, int] | None = None
                      ) -> ATcTRecord:
        """Return a unique gas-phase structure or require an explicit ID.

        The ATcT image SMILES can omit state information. Ambiguous structures
        are never resolved by formula or by the order of table rows.
        """
        requested = _canonical_smiles(smiles)
        if not requested:
            raise ValueError(f'Invalid reference SMILES {smiles!r}.')
        found = [item for item in self.records
                 if item.phase == 'g' and item.formation_0k_kj_mol is not None
                 and item.atct_id.endswith('*0')
                 and re.search(r'\(g\)\s*$', item.preferred_formula)
                 and (formula is None or
                      _preferred_formula_atoms(item.preferred_formula) == formula)
                 and item.smiles and _canonical_smiles(item.smiles) == requested]
        if len(found) != 1:
            raise LookupError(f'{smiles!r} has {len(found)} gas-phase ATcT matches; '
                              'select an explicit ATcT ID for this electronic state.')
        return found[0]


@lru_cache(maxsize=8192)
def _canonical_smiles(value: str) -> str | None:
    try:
        import pybel
    except ImportError:
        try:
            from openbabel import pybel
        except ImportError as exc:
            raise ImportError('Open Babel Python bindings are required for ATcT matching.') from exc
    ob = pybel.ob

    old_level = ob.obErrorLog.GetOutputLevel()
    try:
        ob.obErrorLog.SetOutputLevel(0)
        molecule = pybel.readstring('smi', value)
        return molecule.write('can').split()[0]
    except (OSError, ValueError, IndexError):
        return None
    finally:
        ob.obErrorLog.SetOutputLevel(old_level)


def _preferred_formula_atoms(value: str) -> dict[str, int] | None:
    formula = value.split('(', 1)[0].replace(' ', '')
    atoms = {}
    consumed = 0
    for match in _FORMULA_ATOM.finditer(formula):
        if match.start() != consumed:
            return None
        atoms[match.group(1)] = atoms.get(match.group(1), 0) + int(match.group(2) or 1)
        consumed = match.end()
    return atoms if consumed == len(formula) and consumed else None


def parse_atct_html(raw: bytes, version: str) -> ATcTTable:
    """Parse the released ATcT HTML, including its occasionally blank fields."""
    if not _VERSION.fullmatch(version):
        raise ValueError('ATcT version must be a pinned release such as 1.222.')
    digest = sha256(raw).hexdigest()
    html = raw.decode('latin-1')
    if not re.search(r'ATcT.*?version\s+' + re.escape(version) +
                     r'\s+of the Thermochemical Network', html,
                     re.IGNORECASE | re.DOTALL):
        raise ValueError('ATcT release header does not match the pinned version.')
    records = []
    ids = set()
    for match in _ROW.finditer(html):
        row = match.group()
        atct_id = _field(row, 'ATcTID')
        if not atct_id:
            raise ValueError('ATcT species row has no ID.')
        if atct_id in ids:
            raise ValueError(f'Duplicate ATcT ID {atct_id}.')
        ids.add(atct_id)
        formula = _field(row, 'Formula')
        phase_match = re.search(r'\((g|l|cr|aq)(?:[,\s)]|$)', formula)
        phase = phase_match.group(1) if phase_match else ''
        image = re.search(r'<img\b[^>]*\balt="([^"]*)"', row,
                          re.IGNORECASE | re.DOTALL)
        units = _field(row, 'Units')
        exact = _field(row, 'Uncert').casefold() == 'exact'
        h0 = _number(_field(row, 'DHf0'), '0 K enthalpy')
        h298 = _number(_field(row, 'DHf298'), '298.15 K enthalpy')
        uncertainty = (None if exact else
                       _number(_field(row, 'Uncert').replace('±', '').strip(),
                               'uncertainty'))
        if (h0 is not None or h298 is not None) and not exact and units != 'kJ/mol':
            raise ValueError(f'{atct_id}: unsupported ATcT unit {units!r}.')
        if exact and (h0 not in (0.0, None) or h298 not in (0.0, None)):
            raise ValueError(f'{atct_id}: nonzero value marked exact.')
        records.append(ATcTRecord(
            atct_id=atct_id, name=_field(row, 'Name'),
            preferred_formula=formula, phase=phase,
            smiles=unescape(image.group(1)) if image else '',
            formation_0k_kj_mol=h0, formation_298k_kj_mol=h298,
            uncertainty_kj_mol=uncertainty, exact=exact,
            version=version, source_sha256=digest))
    if not records:
        raise ValueError('ATcT release contains no species rows.')
    return ATcTTable(version, ATCT_URL.format(version=version), digest,
                     tuple(records))


def load_atct(version: str, cache_dir: str | Path, *, refresh: bool = False) -> ATcTTable:
    """Fetch once, then use the same pinned release until refresh is requested."""
    if not _VERSION.fullmatch(version):
        raise ValueError('ATcT version must be pinned.')
    cache = Path(cache_dir)
    path = cache / f'atct_{version}.html'
    if refresh or not path.is_file():
        with urlopen(ATCT_URL.format(version=version), timeout=30) as response:
            raw = response.read()
        table = parse_atct_html(raw, version)
        cache.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
        try:
            temporary.write_bytes(raw)
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
        return table
    return parse_atct_html(path.read_bytes(), version)
