"""Read a pinned release of the official Active Thermochemical Tables.

ATcT identifies chemical *states*, not just formulas. A reference is selected
by its ATcT ID or by an unambiguous gas-phase structure match. The source
API response and its digest remain part of every derived formation enthalpy.
The legacy HTML parser remains available for archived calculations.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
from html import unescape
from functools import lru_cache
from pathlib import Path
import json
import math
import os
import re
from typing import Iterable, Mapping
from urllib.request import urlopen


ATCT_API_URL = 'https://atct.anl.gov/api/v1/'
ATCT_API_ALL_URL = ATCT_API_URL + 'all/'
ATCT_LEGACY_URL = ('https://atct.anl.gov/Thermochemical%20Data/'
                   'version%20{version}/')
# Kept as an import-compatible alias for historical callers.
ATCT_URL = ATCT_LEGACY_URL
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
    return ATcTTable(version, ATCT_LEGACY_URL.format(version=version), digest,
                     tuple(records))


def _api_field(record: Mapping, *names: str):
    for name in names:
        if name in record:
            return record[name]
    return None


def parse_atct_api(raw: bytes, version: str | None = None, *,
                   source_url: str = ATCT_API_ALL_URL) -> ATcTTable:
    """Parse either official v1 species schema into a reproducible table."""
    digest = sha256(raw).hexdigest()
    try:
        payload = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError('ATcT API response is not valid JSON.') from exc
    if isinstance(payload, dict) and isinstance(payload.get('items'), list):
        payload = payload['items']
    if not isinstance(payload, list) or not payload:
        raise ValueError('ATcT API response contains no species records.')

    versions = {str(_api_field(item, 'ATcT_TN_Version') or '')
                for item in payload if isinstance(item, Mapping)}
    if len(versions) != 1 or '' in versions:
        raise ValueError('ATcT API response does not identify one table version.')
    found_version = versions.pop()
    if not _VERSION.fullmatch(found_version):
        raise ValueError(f'Invalid ATcT API table version {found_version!r}.')
    if version is not None and version != found_version:
        raise ValueError(f'ATcT API returned {found_version}, expected {version}.')
    version = found_version

    records = []
    ids = set()
    for item in payload:
        if not isinstance(item, Mapping):
            raise ValueError('ATcT API species record is not an object.')
        atct_id = str(_api_field(item, 'ATcT_ID') or '').strip()
        if not atct_id:
            raise ValueError('ATcT API species record has no ID.')
        if atct_id in ids:
            raise ValueError(f'Duplicate ATcT ID {atct_id}.')
        ids.add(atct_id)
        formula = ' '.join(str(_api_field(item, 'Formula') or '').split())
        phase_match = re.search(r'\((g|l|cr|aq)(?:[,\s)]|$)', formula)
        phase = phase_match.group(1) if phase_match else ''
        uncertainty_raw = _api_field(
            item, 'Delta_Hf298K_uncertainty', '∆fH_298K_uncertainty')
        exact = (str(uncertainty_raw).strip().casefold() == 'exact')
        uncertainty = (None if exact else
                       _number(str(uncertainty_raw).strip()
                               if uncertainty_raw is not None else '',
                               '298 K uncertainty'))
        h0_raw = _api_field(item, 'Delta_Hf_0K', '∆fH_0K')
        h298_raw = _api_field(item, 'Delta_Hf_298K', '∆fH_298K')
        h0 = _number(str(h0_raw).strip() if h0_raw is not None else '',
                     '0 K enthalpy')
        h298 = _number(str(h298_raw).strip() if h298_raw is not None else '',
                       '298.15 K enthalpy')
        units = str(_api_field(item, 'unit', 'units') or '').strip()
        if (h0 is not None or h298 is not None) and not exact and units != 'kJ/mol':
            raise ValueError(f'{atct_id}: unsupported ATcT unit {units!r}.')
        if exact and (h0 not in (0.0, None) or h298 not in (0.0, None)):
            raise ValueError(f'{atct_id}: nonzero value marked exact.')
        records.append(ATcTRecord(
            atct_id=atct_id,
            name=str(_api_field(item, 'Name') or '').strip(),
            preferred_formula=formula,
            phase=phase,
            smiles=str(_api_field(item, 'SMILES') or '').strip(),
            formation_0k_kj_mol=h0,
            formation_298k_kj_mol=h298,
            # API v1 publishes the conventional 298 K uncertainty.  It is
            # retained as metadata and is never substituted for a 0 K value.
            uncertainty_kj_mol=uncertainty,
            exact=exact,
            version=version,
            source_sha256=digest))
    return ATcTTable(version, source_url, digest, tuple(records))


def _write_cache(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    try:
        temporary.write_bytes(raw)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def load_atct(version: str, cache_dir: str | Path, *,
              refresh: bool = False) -> ATcTTable:
    """Cache the official bulk API response for one explicitly pinned version."""
    if not _VERSION.fullmatch(version):
        raise ValueError('ATcT version must be pinned.')
    cache = Path(cache_dir)
    path = cache / f'atct_{version}.json'
    if refresh or not path.is_file():
        with urlopen(ATCT_API_ALL_URL, timeout=30) as response:
            raw = response.read()
        table = parse_atct_api(raw, version)
        _write_cache(path, raw)
        return table
    return parse_atct_api(path.read_bytes(), version)


def _fetch_reference_records(smiles: Iterable[str],
                             reference_ids: Mapping[str, str]) -> list[dict]:
    """Use the public ``atct`` client for exact, rate-limited API queries."""
    try:
        from atct.api import get_species_by_atctid, get_species_by_smiles
    except ImportError as exc:  # pragma: no cover - dependency is declared
        raise ImportError("Install KinBot's 'atct' dependency.") from exc
    records = {}
    for value in sorted(set(smiles)):
        if value in reference_ids:
            species = get_species_by_atctid(reference_ids[value], block=True)
            candidates = [species]
        else:
            page = get_species_by_smiles(value, limit=100, offset=0, block=True)
            if page.total > len(page.items):
                raise RuntimeError(
                    f'ATcT returned more than 100 states for {value!r}; '
                    'select explicit ATcT IDs.')
            candidates = page.items
        for species in candidates:
            record = species.to_dict()
            # atct 1.0.1 models the numeric values but omits the API's
            # ``units`` field when converting Species back to a dictionary.
            # API v1 thermochemical values are returned in kJ/mol.
            record['units'] = (None if
                               record.get('∆fH_298K_uncertainty') == 'exact'
                               else 'kJ/mol')
            records[record['ATcT_ID']] = record
    return [records[key] for key in sorted(records)]


def load_atct_references(smiles: Iterable[str], version: str,
                         cache_dir: str | Path, *,
                         reference_ids: Mapping[str, str] | None = None,
                         refresh: bool = False) -> ATcTTable:
    """Resolve and cache only the ATcT states needed by one CBH problem."""
    if not _VERSION.fullmatch(version):
        raise ValueError('ATcT version must be pinned.')
    values = tuple(sorted(set(smiles)))
    if not values:
        raise ValueError('At least one ATcT reference SMILES is required.')
    reference_ids = dict(reference_ids or {})
    query = json.dumps({'smiles': values, 'ids': reference_ids},
                       sort_keys=True, separators=(',', ':')).encode()
    query_hash = sha256(query).hexdigest()[:16]
    path = Path(cache_dir) / f'atct_{version}_refs_{query_hash}.json'
    source_url = ATCT_API_URL
    if refresh or not path.is_file():
        payload = _fetch_reference_records(values, reference_ids)
        raw = json.dumps(payload, ensure_ascii=False, sort_keys=True,
                         separators=(',', ':')).encode()
        table = parse_atct_api(raw, version, source_url=source_url)
        _write_cache(path, raw)
        return table
    return parse_atct_api(path.read_bytes(), version, source_url=source_url)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description='Cache and verify a pinned ATcT v1 API snapshot.')
    parser.add_argument('version', help='expected ATcT table version')
    parser.add_argument('cache_dir', help='snapshot cache directory')
    parser.add_argument('smiles', nargs='*',
                        help='optional reference SMILES; omit for /all/')
    parser.add_argument('--refresh', action='store_true')
    args = parser.parse_args(argv)
    table = (load_atct_references(
        args.smiles, args.version, args.cache_dir, refresh=args.refresh)
             if args.smiles else
             load_atct(args.version, args.cache_dir, refresh=args.refresh))
    print(json.dumps({
        'version': table.version,
        'records': len(table.records),
        'source_url': table.source_url,
        'source_sha256': table.source_sha256,
    }, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
