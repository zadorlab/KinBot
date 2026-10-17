"""Normalize Gaussian option lists before merging or ASE route rendering."""
import re


def route_options(value):
    """Accept an optional outer pair; split only top-level commas.

    Named nested options (e.g. ``ReadAtoms(C,H)``) retain their grouping.
    Anonymous nested lists cannot be emitted as Gaussian suboptions.
    """
    if value is None or value is False or value == '':
        return []
    if isinstance(value, (list, tuple)):
        return [part for item in value for part in route_options(item)]
    text = str(value).strip()
    depth = 0
    close = None
    for index, char in enumerate(text):
        if char == '(':
            depth += 1
        elif char == ')':
            depth -= 1
            if depth < 0:
                raise ValueError(f'Unbalanced Gaussian option list: {text!r}')
            if depth == 0 and close is None:
                close = index
    if depth:
        raise ValueError(f'Unbalanced Gaussian option list: {text!r}')
    if text.startswith('(') and close == len(text) - 1:
        return route_options(text[1:-1])
    parts, start, depth = [], 0, 0
    for index, char in enumerate(text + ','):
        depth += (char == '(') - (char == ')')
        if char == ',' and depth == 0:
            part = text[start:index].strip()
            if not part or part.startswith('('):
                raise ValueError(f'Invalid Gaussian option list: {text!r}')
            parts.append(part)
            start = index + 1
    return parts


def merge_route_options(base, override):
    merged = {}
    for option in route_options(base) + route_options(override):
        key = re.split(r'[=(]', option, maxsplit=1)[0].strip().lower()
        merged[key] = option
    return ','.join(merged.values())


def merge_gaussian_kwargs(base, override):
    # These are whole calculator/Link-0 values, not option lists.
    scalar = {'method', 'basis', 'fitting_basis', 'xc', 'charge', 'mult',
              'label', 'directory', 'mem', 'chk', 'oldchk', 'nprocshared',
              'nprocs', 'cpu', 'rwf', 'command', 'basisfile', 'basis_set',
              'addsec', 'ioplist', 'extra', 'spinlist', 'output_type', 'scf_recovery'}
    result = dict(base)
    for key, value in override.items():
        existing = next((k for k in result if k.lower() == key.lower()), key)
        if key.lower() in scalar or value is None or value is False or value == '':
            result[existing] = value
        else:
            result[existing] = merge_route_options(result.get(existing), value)
    return result
