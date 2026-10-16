"""Explicit numerical controls for generated Molpro probe inputs.

These controls never change the active space, roots, weights or Hamiltonian.
Native inputs and execution records retain the exact tolerances used.
"""
from copy import deepcopy
import re

from kinbot.theory import fortran_threshold, convergence_exponent


def mr_thresholds(options):
    """Molpro THRVAR is a squared coefficient norm, not an energy tolerance."""
    energy = fortran_threshold(options.get('energy_convergence', 1.e-8))
    exponent = convergence_exponent(options.get('coefficient_convergence', 1.e-5))
    return {'thrden': energy, 'thrvar': f'1.d-{2 * exponent}'}


def validate_solver_options(options):
    allowed = {
        'scf': {'max_iterations', 'energy_convergence'},
        'casscf': {'max_iterations', 'energy_convergence', 'gradient_convergence'},
        'correlation': {'max_iterations', 'max_reference_iterations', 'max_reference_vectors',
                        'energy_convergence', 'coefficient_convergence'},
    }
    if not isinstance(options, dict) or set(options) - set(allowed):
        raise ValueError('solver_options must contain scf, casscf or correlation controls.')
    for stage, controls in options.items():
        if not isinstance(controls, dict) or set(controls) - allowed[stage]:
            raise ValueError(f'Unknown solver_options.{stage} controls.')
        for key, value in controls.items():
            if key.endswith('_convergence'):
                fortran_threshold(value, f'solver_options.{stage}.{key}')
            elif type(value) is not int or value < 1:
                raise ValueError(f'solver_options.{stage}.{key} must be a positive integer.')
    return deepcopy(options)


def _option(text, key, value):
    pattern = r'\b' + key + r'\s*=\s*[^,;\s}]+'
    return re.sub(pattern, f'{key}={value}', text, flags=re.I)


def apply_solver_options(text, options):
    """Apply explicit overrides to generated decks; idempotent and fail closed.

    SCF energy uses GTHRESH and an explicit local ENERGY option; MULTI uses ACCURACY;
    correlated energies/coefficients use local, literal THRDEN/THRVAR.
    The local THRESH/ACCU exponent syntax is deliberately not mixed into this.
    """
    options = validate_solver_options(options)
    scf = options.get('scf', {})
    global_values = {}
    if 'energy_convergence' in scf:
        global_values['energy'] = fortran_threshold(scf['energy_convergence'])
    if 'coefficient_convergence' in options.get('correlation', {}):
        global_values['coeff'] = fortran_threshold(options['correlation']['coefficient_convergence'])
    if global_values:
        pattern = r'(?im)^gthresh,[^\n]+'
        def global_energy(match):
            body = match[0]
            for key, value in global_values.items():
                if re.search(r'\b' + key + r'\s*=', body, re.I):
                    body = _option(body, key, value)
                else:
                    body += f',{key}={value}'
            return body
        text, count = re.subn(pattern, global_energy, text)
        if count != 1:
            raise ValueError('Solver override requires one generated GTHRESH card.')

    def program(match):
        name, body = match[1].lower(), match[2]
        if name in ('rhf', 'uhf'):
            if 'energy_convergence' in scf:
                value = fortran_threshold(scf['energy_convergence'])
                if re.search(r'\benergy\s*=', body, re.I):
                    body = _option(body, 'energy', value)
                else:
                    head, separator, rest = body.partition(';')
                    body = head + f',energy={value}' + separator + rest
            if 'max_iterations' in scf:
                if re.search(r'(?:[;\n])\s*maxit(?:er)?\s*,', body, re.I):
                    raise ValueError('Competing SCF MAXIT card; prepare a reviewed input.')
                if re.search(r'\bmaxit\s*=', body, re.I):
                    body = _option(body, 'maxit', scf['max_iterations'])
                else:
                    body = f',maxit={scf["max_iterations"]}' + body
        elif name == 'multi':
            controls = options.get('casscf', {})
            if re.search(r'\bcas-ci\b', body, re.I):
                return match[0]  # Never turn a frozen-orbital control into CASSCF.
            if 'max_iterations' in controls:
                body, count = re.subn(r'(?i)\bmaxit\s*,\s*\d+',
                                      f'maxit,{controls["max_iterations"]}', body)
                if count != 1:
                    raise ValueError('CASSCF override requires one generated MAXIT card.')
            values = {k.removesuffix('_convergence'): fortran_threshold(v)
                      for k, v in controls.items() if k.endswith('_convergence')}
            if values:
                # Generated ACCURACY directives contain only energy/gradient.
                def accuracy(match):
                    fields = dict(re.findall(r'(\w+)\s*=\s*([^,;\s}]+)', match[1]))
                    fields.update(values)
                    return 'accuracy,' + ','.join(f'{k}={v}' for k, v in fields.items())
                body, count = re.subn(r'(?i)accuracy,([^;\n}]+)', accuracy, body)
                if count == 0:
                    body += ';accuracy,' + ','.join(f'{k}={v}' for k, v in values.items()) + ';'
                elif count != 1:
                    raise ValueError('Competing CASSCF ACCURACY cards.')
        else:
            controls = options.get('correlation', {})
            if not controls:
                return match[0]
            if re.search(r'(?:[;\n])\s*(?:thresh|accu|maxiter|maxdav)\s*,', body, re.I):
                raise ValueError('Competing MR solver cards; prepare a reviewed input.')
            head, separator, rest = body.partition(';')
            for key, value in mr_thresholds(controls).items():
                public = {'thrden': 'energy_convergence', 'thrvar': 'coefficient_convergence'}[key]
                if public not in controls:
                    continue
                if re.search(r'\b' + key + r'\s*=', body, re.I):
                    head, rest = _option(head, key, value), _option(rest, key, value)
                else:
                    head += f',{key}={value}'
            body = head + separator + rest
            for key, native in (('max_iterations', 'maxit'),
                                ('max_reference_iterations', 'maxiti'),
                                ('max_reference_vectors', 'maxvi')):
                if key not in controls:
                    continue
                if re.search(r'\b' + native + r'\s*=', body, re.I):
                    body = _option(body, native, controls[key])
                else:
                    body += ('' if body.endswith(';') else ';') + f'option,{native}={controls[key]};'
        return '{' + match[1] + body + '}'

    return re.sub(r'\{\s*(rhf|uhf|multi|rs2(?:c|-f12)?|mrci(?:-f12)?)'
                  r'(?=[,;\s}])([^{}]*)\}', program, text, flags=re.I)
