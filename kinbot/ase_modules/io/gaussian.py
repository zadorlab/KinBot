import re
import warnings
from collections.abc import Iterable
from copy import deepcopy

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.gaussian import Gaussian
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import atomic_masses_iupac2016, chemical_symbols
from ase.units import Bohr, Hartree

try:
    from ase.io.zmatrix import parse_zmatrix
except ModuleNotFoundError:
    from kinbot.ase_modules.io.zmatrix import parse_zmatrix

_link0_keys = [
    'mem',
    'chk',
    'oldchk',
    'schk',
    'rwf',
    'oldmatrix',
    'oldrawmatrix',
    'int',
    'd2e',
    'save',
    'nosave',
    'errorsave',
    'cpu',
    'nprocshared',
    'gpucpu',
    'lindaworkers',
    'usessh',
    'ssh',
    'debuglinda',
]


_link0_special = [
    'kjob',
    'subst',
]


# Certain problematic methods do not provide well-defined potential energy
# surfaces, because these "composite" methods involve geometry optimization
# and/or vibrational frequency analysis. In addition, the "energy" calculated
# by these methods are typically ZPVE corrected and/or temperature dependent
# free energies.
_problem_methods = [
    'cbs-4m', 'cbs-qb3', 'cbs-apno',
    'g1', 'g2', 'g3', 'g4', 'g2mp2', 'g3mp2', 'g3b3', 'g3mp2b3', 'g4mp4',
    'w1', 'w1u', 'w1bd', 'w1ro',
]


_xc_to_method = dict(
    pbe='pbepbe',
    pbe0='pbe1pbe',
    hse06='hseh1pbe',
    hse03='ohse2pbe',
    lda='svwn',  # gaussian "knows about" LSDA, but maybe not LDA.
    tpss='tpsstpss',
    revtpss='revtpssrevtpss',
)

_nuclear_prop_names = ['spin', 'zeff', 'qmom', 'nmagm', 'znuc',
                       'radnuclear', 'iso']


def _get_molecule_spec(atoms, nuclear_props):
    ''' Generate the molecule specification section to write
    to the Gaussian input file, from the Atoms object and dict
    of nuclear properties'''
    molecule_spec = []
    for i, atom in enumerate(atoms):
        symbol_section = atom.symbol + '('
        # Check whether any nuclear properties of the atom have been set,
        # and if so, add them to the symbol section.
        nuclear_props_set = False
        for keyword, array in nuclear_props.items():
            if array is not None and array[i] is not None:
                string = keyword + '=' + str(array[i]) + ', '
                symbol_section += string
                nuclear_props_set = True

        # Check whether the mass of the atom has been modified,
        # and if so, add it to the symbol section:
        mass_set = False
        symbol = atom.symbol
        expected_mass = atomic_masses_iupac2016[chemical_symbols.index(
            symbol)]
        if expected_mass != atoms[i].mass:
            mass_set = True
            string = 'iso' + '=' + str(atoms[i].mass)
            symbol_section += string

        if nuclear_props_set or mass_set:
            symbol_section = symbol_section.strip(', ')
            symbol_section += ')'
        else:
            symbol_section = symbol_section.strip('(')

        # Then attach the properties appropriately
        # this formatting was chosen for backwards compatibility reasons, but
        # it would probably be better to
        # 1) Ensure proper spacing between entries with explicit spaces
        # 2) Use fewer columns for the element
        # 3) Use 'e' (scientific notation) instead of 'f' for positions
        molecule_spec.append('{:<10s}{:20.10f}{:20.10f}{:20.10f}'.format(
            symbol_section, *atom.position))

    # unit cell vectors, in case of periodic boundary conditions
    for ipbc, tv in zip(atoms.pbc, atoms.cell):
        if ipbc:
            molecule_spec.append('TV {:20.10f}{:20.10f}{:20.10f}'.format(*tv))

    molecule_spec.append('')

    return molecule_spec


def _format_output_type(output_type):
    ''' Given a letter: output_type, return
    a string formatted for a gaussian input file'''
    # Change output type to P if it has been set to T, because ASE
    # does not support reading info from 'terse' output files.
    if output_type is None or output_type == '' or 't' in output_type.lower():
        output_type = 'P'

    return '#{}'.format(output_type)


def _check_problem_methods(method):
    ''' Check method string for problem methods and warn appropriately'''
    if method.lower() in _problem_methods:
        warnings.warn(
            'The requested method, {}, is a composite method. Composite '
            'methods do not have well-defined potential energy surfaces, '
            'so the energies, forces, and other properties returned by '
            'ASE may not be meaningful, or they may correspond to a '
            'different geometry than the one provided. '
            'Please use these methods with caution.'.format(method)
        )


def _pop_link0_params(params):
    '''Takes the params dict and returns a dict with the link0 keywords
    removed, and a list containing the link0 keywords and options
    to be written to the gaussian input file'''
    params = deepcopy(params)
    out = []
    # Remove keywords from params, so they are not set again later in
    # route section
    for key in _link0_keys:
        if key not in params:
            continue
        val = params.pop(key)
        if not val or (isinstance(val, str) and key.lower() == val.lower()):
            out.append('%{}'.format(key))
        else:
            out.append('%{}={}'.format(key, val))

    # These link0 keywords have a slightly different syntax
    for key in _link0_special:
        if key not in params:
            continue
        val = params.pop(key)
        if not isinstance(val, str) and isinstance(val, Iterable):
            val = ' '.join(val)
        out.append('%{} L{}'.format(key, val))

    return params, out


def _format_method_basis(output_type, method, basis, fitting_basis):
    output_string = ""
    if basis and method and fitting_basis:
        output_string = '{} {}/{}/{} ! ASE formatted method and basis'.format(
            output_type, method, basis, fitting_basis)
    elif basis and method:
        output_string = '{} {}/{} ! ASE formatted method and basis'.format(
            output_type, method, basis)
    else:
        output_string = '{}'.format(output_type)
        for value in [method, basis]:
            if value is not None:
                output_string += ' {}'.format(value)
    return output_string


def _format_route_params(params):
    '''Get keywords and values from the params dictionary and return
    as a list of lines to add to the gaussian input file'''
    out = []
    for key, val in params.items():
        # assume bare keyword if val is falsey, i.e. '', None, False, etc.
        # also, for backwards compatibility: assume bare keyword if key and
        # val are the same
        if not val or (isinstance(val, str) and key.lower() == val.lower()):
            out.append(key)
        elif not isinstance(val, str) and isinstance(val, Iterable):
            out.append('{}({})'.format(key, ','.join(val)))
        else:
            out.append('{}({})'.format(key, val))
    return out


def _format_addsec(addsec):
    '''Format addsec string as a list of lines to be added to the gaussian
     input file.'''
    out = []
    if addsec is not None:
        #out.append('')  # modification Judit Zador
        if isinstance(addsec, str):
            out.append(addsec)
        elif isinstance(addsec, Iterable):
            out += list(addsec)
    return out


def _format_basis_set(basis, basisfile, basis_set):
    '''Format either: the basis set filename (basisfile), the basis set file
    contents (from reading basisfile), or the basis_set text as a list of
    strings to be added to the gaussian input file.'''
    out = []
    # if basis='gen', set basisfile. Either give a path to a basisfile, or
    # read in the provided file and paste it verbatim
    if basisfile is not None:
        if basisfile[0] == '@':
            out.append(basisfile)
        else:
            with open(basisfile, 'r') as fd:
                out.append(fd.read())
    elif basis_set is not None:
        out.append(basis_set)
    else:
        if basis is not None and basis.lower() == 'gen':
            raise ValueError('Please set basisfile or basis_set')
    return out


def write_gaussian_in(fd, atoms, properties=['energy'],
                      method=None, basis=None, fitting_basis=None,
                      output_type='P', basisfile=None, basis_set=None,
                      xc=None, charge=None, mult=None, extra=None,
                      ioplist=None, addsec=None, spinlist=None,
                      zefflist=None, qmomlist=None, nmagmlist=None,
                      znuclist=None, radnuclearlist=None,
                      **params):
    '''
    Generates a Gaussian input file

    Parameters
    -----------
    fd: file-like
        where the Gaussian input file will be written
    atoms: Atoms
        Structure to write to the input file
    properties: list
        Properties to calculate
    method: str
        Level of theory to use, e.g. ``hf``, ``ccsd``, ``mp2``, or ``b3lyp``.
        Overrides ``xc`` (see below).
    xc: str
        Level of theory to use. Translates several XC functionals from
        their common name (e.g. ``PBE``) to their internal Gaussian name
        (e.g. ``PBEPBE``).
    basis: str
        The basis set to use. If not provided, no basis set will be requested,
        which usually results in ``STO-3G``. Maybe omitted if basisfile is set
        (see below).
    fitting_basis: str
        The name of the fitting basis set to use.
    output_type: str
        Level of output to record in the Gaussian
        output file - this may be ``N``- normal or ``P`` -
        additional.
    basisfile: str
        The name of the basis file to use. If a value is provided, basis may
        be omitted (it will be automatically set to 'gen')
    basis_set: str
        The basis set definition to use. This is an alternative
        to basisfile, and would be the same as the contents
        of such a file.
    charge: int
        The system charge. If not provided, it will be automatically
        determined from the ``Atoms`` object’s initial_charges.
    mult: int
        The system multiplicity (``spin + 1``). If not provided, it will be
        automatically determined from the ``Atoms`` object’s
        ``initial_magnetic_moments``.
    extra: str
        Extra lines to be included in the route section verbatim.
        It should not be necessary to use this, but it is included for
        backwards compatibility.
    ioplist: list
        A collection of IOPs definitions to be included in the route line.
    addsec: str
        Text to be added after the molecular geometry specification, e.g. for
        defining masses with ``freq=ReadIso``.
    spinlist: list
        A list of nuclear spins to be added into the nuclear
        propeties section of the molecule specification.
    zefflist: list
        A list of effective charges to be added into the nuclear
        propeties section of the molecule specification.
    qmomlist: list
        A list of nuclear quadropole moments to be added into
        the nuclear propeties section of the molecule
        specification.
    nmagmlist: list
        A list of nuclear magnetic moments to be added into
        the nuclear propeties section of the molecule
        specification.
    znuclist: list
        A list of nuclear charges to be added into the nuclear
        propeties section of the molecule specification.
    radnuclearlist: list
        A list of nuclear radii to be added into the nuclear
        propeties section of the molecule specification.
    params: dict
        Contains any extra keywords and values that will be included in either
        the link0 section or route section of the gaussian input file.
        To be included in the link0 section, the keyword must be one of the
        following: ``mem``, ``chk``, ``oldchk``, ``schk``, ``rwf``,
        ``oldmatrix``, ``oldrawmatrix``, ``int``, ``d2e``, ``save``,
        ``nosave``, ``errorsave``, ``cpu``, ``nprocshared``, ``gpucpu``,
        ``lindaworkers``, ``usessh``, ``ssh``, ``debuglinda``.
        Any other keywords will be placed (along with their values) in the
        route section.
    '''

    params = deepcopy(params)

    if properties is None:
        properties = ['energy']

    output_type = _format_output_type(output_type)

    # basis can be omitted if basisfile is provided
    if basis is None:
        if basisfile is not None or basis_set is not None:
            basis = 'gen'

    # determine method from xc if it is provided
    if method is None:
        if xc is not None:
            method = _xc_to_method.get(xc.lower(), xc)

    # If the user requests a problematic method, rather than raising an error
    # or proceeding blindly, give the user a warning that the results parsed
    # by ASE may not be meaningful.
    if method is not None:
        _check_problem_methods(method)

    # determine charge from initial charges if not passed explicitly
    if charge is None:
        charge = atoms.get_initial_charges().sum()

    # determine multiplicity from initial magnetic moments
    # if not passed explicitly
    if mult is None:
        mult = atoms.get_initial_magnetic_moments().sum() + 1

    # set up link0 arguments
    out = []
    params, link0_list = _pop_link0_params(params)
    out.extend(link0_list)

    # begin route line
    # note: unlike in old calculator, each route keyword is put on its own
    # line.
    out.append(_format_method_basis(output_type, method, basis, fitting_basis))

    # If the calculator's parameter dictionary contains an isolist, we ignore
    # this - it is up to the user to attach this info as the atoms' masses
    # if they wish for it to be used:
    params.pop('isolist', None)

    # Any params left will belong in the route section of the file:
    out.extend(_format_route_params(params))

    if ioplist is not None:
        out.append('IOP(' + ', '.join(ioplist) + ')')

    # raw list of explicit keywords for backwards compatibility
    if extra is not None:
        out.append(extra)

    # Add 'force' iff the user requested forces, since Gaussian crashes when
    # 'force' is combined with certain other keywords such as opt and irc.
    if 'forces' in properties and 'force' not in params:
        out.append('force')

    # header, charge, and mult
    out += ['', 'Gaussian input prepared by ASE', '',
            '{:.0f} {:.0f}'.format(charge, mult)]

    # make dict of nuclear properties:
    nuclear_props = {'spin': spinlist, 'zeff': zefflist, 'qmom': qmomlist,
                     'nmagm': nmagmlist, 'znuc': znuclist,
                     'radnuclear': radnuclearlist}
    nuclear_props = {k: v for k, v in nuclear_props.items() if v is not None}

    # atomic positions and nuclear properties:
    molecule_spec = _get_molecule_spec(atoms, nuclear_props)
    for line in molecule_spec:
        out.append(line)

    out.extend(_format_basis_set(basis, basisfile, basis_set))

    out.extend(_format_addsec(addsec))

    out += ['', '']
    fd.write('\n'.join(out))


# Regexp for reading an input file:

_re_link0 = re.compile(r'^\s*%([^\=\)\(!]+)=?([^\=\)\(!]+)?(!.+)?')
# Link0 lines are in the format:
# '% keyword = value' or '% keyword'
# (with or without whitespaces)

_re_output_type = re.compile(r'^\s*#\s*([NPTnpt]?)\s*')
# The start of the route section begins with a '#', and then may
# be followed by the desired level of output in the output file: P, N or T.

_re_method_basis = re.compile(
    r"\s*([\w-]+)\s*\/([^/=!]+)([\/]([^!]+))?\s*(!.+)?")
# Matches method, basis and optional fitting basis in the format:
# method/basis/fitting_basis ! comment
# They will appear in this format if the Gaussian file has been generated
# by ASE using a calculator with the basis and method keywords set.

_re_chgmult = re.compile(r'^\s*[+-]?\d+(?:,\s*|\s+)[+-]?\d+\s*$')
# This is a bit more complex of a regex than we typically want, but it
# can be difficult to determine whether a line contains the charge and
# multiplicity, rather than just another route keyword. By making sure
# that the line contains exactly two *integers*, separated by either
# a comma (and possibly whitespace) or some amount of whitespace, we
# can be more confident that we've actually found the charge and multiplicity.

_re_nuclear_props = re.compile(r'\(([^\)]+)\)')
# Matches the nuclear properties, which are contained in parantheses.

# The following methods are used in GaussianConfiguration's
# parse_gaussian_input method:


def _get_link0_param(link0_match):
    '''Gets link0 keyword and option from a re.Match and returns them
    in a dictionary format'''
    value = link0_match.group(2)
    if value is not None:
        value = value.strip()
    else:
        value = ''
    return {link0_match.group(1).lower().strip(): value.lower()}


def _get_all_link0_params(link0_section):
    ''' Given a string link0_section which contains the link0
    section of a gaussian input file, returns a dictionary of
    keywords and values'''
    parameters = {}
    for line in link0_section:
        link0_match = _re_link0.match(line)
        link0_param = _get_link0_param(link0_match)
        if link0_param is not None:
            parameters.update(link0_param)
    return parameters


def _convert_to_symbol(string):
    '''Converts an input string into a format
    that can be input to the 'symbol' parameter of an
    ASE Atom object (can be a chemical symbol (str)
    or an atomic number (int)).
    This is achieved by either stripping any
    integers from the string, or converting a string
    containing an atomic number to integer type'''
    symbol = _validate_symbol_string(string)
    if symbol.isnumeric():
        atomic_number = int(symbol)
        symbol = chemical_symbols[atomic_number]
    else:
        match = re.match(r'([A-Za-z]+)', symbol)
        symbol = match.group(1)
    return symbol


def _validate_symbol_string(string):
    if "-" in string:
        raise IOError("ERROR: Could not read the Gaussian input file, as"
                         " molecule specifications for molecular mechanics "
                         "calculations are not supported.")
    return string


def _get_key_value_pairs(line):
    '''Reads a line of a gaussian input file, which contains keywords and options
    separated according to the rules of the route section.

    Parameters
    ----------
    line (string)
        A line of a gaussian input file.

    Returns
    ---------
    params (dict)
        Contains the keywords and options found in the line.
    '''
    params = {}
    line = line.strip(' #')
    line = line.split('!')[0]  # removes any comments
    # First, get the keywords and options sepatated with
    # parantheses:
    match_iterator = re.finditer(r'\(([^\)]+)\)', line)
    index_ranges = []
    for match in match_iterator:
        index_range = [match.start(0), match.end(0)]
        options = match.group(1)
        # keyword is last word in previous substring:
        keyword_string = line[:match.start(0)]
        keyword_match_iter = [k for k in re.finditer(
            r'[^\,/\s]+', keyword_string) if k.group() != '=']
        keyword = keyword_match_iter[-1].group().strip(' =')
        index_range[0] = keyword_match_iter[-1].start()
        params.update({keyword.lower(): options.lower()})
        index_ranges.append(index_range)

    # remove from the line the keywords and options that we have saved:
    index_ranges.reverse()
    for index_range in index_ranges:
        start = index_range[0]
        stop = index_range[1]
        line = line[0: start:] + line[stop + 1::]

    # Next, get the keywords and options separated with
    # an equals sign, and those without an equals sign
    # must be keywords without options:

    # remove any whitespaces around '=':
    line = re.sub(r'\s*=\s*', '=', line)
    line = [x for x in re.split(r'[\s,\/]', line) if x != '']

    for s in line:
        if '=' in s:
            s = s.split('=')
            keyword = s.pop(0)
            options = s.pop(0)
            params.update({keyword.lower(): options.lower()})
        else:
            if len(s) > 0:
                params.update({s.lower(): None})

    return params


def _get_route_params(line):
    '''Reads keywords and values from a line in
    a Gaussian input file's route section,
    and returns them as a dictionary'''
    method_basis_match = _re_method_basis.match(line)
    if method_basis_match:
        params = {}
        ase_gen_comment = '! ASE formatted method and basis'
        if method_basis_match.group(5) == ase_gen_comment:
            params['method'] = method_basis_match.group(1).strip().lower()
            params['basis'] = method_basis_match.group(2).strip().lower()
            if method_basis_match.group(4):
                params['fitting_basis'] = method_basis_match.group(
                    4).strip().lower()
            return params

    return _get_key_value_pairs(line)


def _get_all_route_params(route_section):
    ''' Given a string: route_section which contains the route
    section of a gaussian input file, returns a dictionary of
    keywords and values'''
    parameters = {}

    for line in route_section:
        output_type_match = _re_output_type.match(line)
        if not parameters.get('output_type') and output_type_match:
            line = line.strip(output_type_match.group(0))
            parameters.update(
                {'output_type': output_type_match.group(1).lower()})
        # read route section
        route_params = _get_route_params(line)
        if route_params is not None:
            parameters.update(route_params)
    return parameters


def _get_charge_mult(chgmult_section):
    '''return a dict with the charge and multiplicity from
    a list chgmult_section that contains the charge and multiplicity
    line, read from a gaussian input file'''
    chgmult_match = _re_chgmult.match(str(chgmult_section))
    try:
        chgmult = chgmult_match.group(0).split()
        return {'charge': int(chgmult[0]), 'mult': int(chgmult[1])}
    except (IndexError, AttributeError):
        raise IOError("ERROR: Could not read the charge and multiplicity "
                         "from the Gaussian input file. These must be 2 "
                         "integers separated with whitespace or a comma.")


def _get_nuclear_props(line):
    ''' Reads any info in parantheses in the line and returns
    a dictionary of the nuclear properties.'''
    nuclear_props_match = _re_nuclear_props.search(line)
    nuclear_props = {}
    if nuclear_props_match:
        nuclear_props = _get_key_value_pairs(nuclear_props_match.group(1))
        updated_nuclear_props = {}
        for key, value in nuclear_props.items():
            if value.isnumeric():
                value = int(value)
            else:
                value = float(value)
            if key not in _nuclear_prop_names:
                if "fragment" in key:
                    warnings.warn("Fragments are not "
                                  "currently supported.")
                warnings.warn("The following nuclear properties "
                              "could not be saved: {}".format(
                                  {key: value}))
            else:
                updated_nuclear_props[key] = value
        nuclear_props = updated_nuclear_props

    for k in _nuclear_prop_names:
        if k not in nuclear_props:
            nuclear_props[k] = None

    return nuclear_props


def _get_atoms_info(line):
    '''Returns the symbol and position of an atom from a line
    in the molecule specification section'''
    nuclear_props_match = _re_nuclear_props.search(line)
    if nuclear_props_match:
        line = line.replace(nuclear_props_match.group(0), '')
    tokens = line.split()
    symbol = _convert_to_symbol(tokens[0])
    pos = list(tokens[1:])

    return symbol, pos


def _get_cartesian_atom_coords(symbol, pos):
    '''Returns the coordinates: pos as a list of floats if they
    are cartesian, and not in z-matrix format'''
    if len(pos) < 3 or (pos[0] == '0' and symbol != 'TV'):
        # In this case, we have a z-matrix definition, so
        # no cartesian coords.
        return
    elif len(pos) > 3:
        raise IOError("ERROR: Gaussian input file could "
                         "not be read as freeze codes are not"
                         " supported. If using cartesian "
                         "coordinates, these must be "
                         "given as 3 numbers separated "
                         "by whitespace.")
    else:
        try:
            return list(map(float, pos))
        except ValueError:
            raise(IOError(
                "ERROR: Molecule specification in"
                "Gaussian input file could not be read"))


def _get_zmatrix_line(line):
    ''' Converts line into the format needed for it to
    be added to the z-matrix contents '''
    line_list = line.split()
    if len(line_list) == 8 and line_list[7] == '1':
        raise IOError(
            "ERROR: Could not read the Gaussian input file"
            ", as the alternative Z-matrix format using "
            "two bond angles instead of a bond angle and "
            "a dihedral angle is not supported.")
    return(line.strip() + '\n')


def _read_zmatrix(zmatrix_contents, zmatrix_vars=None):
    ''' Reads a z-matrix (zmatrix_contents) using its list of variables
    (zmatrix_vars), and returns atom positions and symbols '''
    try:
        atoms = parse_zmatrix(zmatrix_contents, defs=zmatrix_vars)
    except (ValueError, AssertionError) as e:
        raise IOError("Failed to read Z-matrix from "
                         "Gaussian input file: ", e)
    except KeyError as e:
        raise IOError("Failed to read Z-matrix from "
                         "Gaussian input file, as symbol: {}"
                         "could not be recognised. Please make "
                         "sure you use element symbols, not "
                         "atomic numbers in the element labels.".format(e))
    positions = atoms.positions
    symbols = atoms.get_chemical_symbols()
    return positions, symbols


def _get_nuclear_props_for_all_atoms(nuclear_props):
    ''' Returns the nuclear properties for all atoms as a dictionary,
    in the format needed for it to be added to the parameters dictionary.'''
    params = {k + 'list': [] for k in _nuclear_prop_names}
    for dictionary in nuclear_props:
        for key, value in dictionary.items():
            params[key + 'list'].append(value)

    for key, array in params.items():
        values_set = False
        for value in array:
            if value is not None:
                values_set = True
        if not values_set:
            params[key] = None
    return params


def _get_atoms_from_molspec(molspec_section):
    ''' Takes a string: molspec_section which contains the molecule
    specification section of a gaussian input file, and returns an atoms
    object that represents this.'''
    # These will contain info that will be attached to the Atoms object:
    symbols = []
    positions = []
    pbc = np.zeros(3, dtype=bool)
    cell = np.zeros((3, 3))
    npbc = 0

    # Will contain a dictionary of nuclear properties for each atom,
    # that will later be saved to the parameters dict:
    nuclear_props = []

    # Info relating to the z-matrix definition (if set)
    zmatrix_type = False
    zmatrix_contents = ""
    zmatrix_var_section = False
    zmatrix_vars = ""

    for line in molspec_section:
        # Remove any comments and replace '/' and ',' with whitespace,
        # as these are equivalent:
        line = line.split('!')[0].replace('/', ' ').replace(',', ' ')
        if (line.split()):
            if zmatrix_type:
                # Save any variables set when defining the z-matrix:
                if zmatrix_var_section:
                    zmatrix_vars += line.strip() + '\n'
                    continue
                elif 'variables' in line.lower():
                    zmatrix_var_section = True
                    continue
                elif 'constants' in line.lower():
                    zmatrix_var_section = True
                    warnings.warn("Constants in the optimisation are "
                                  "not currently supported. Instead "
                                  "setting constants as variables.")
                    continue

            symbol, pos = _get_atoms_info(line)
            current_nuclear_props = _get_nuclear_props(line)

            if not zmatrix_type:
                pos = _get_cartesian_atom_coords(symbol, pos)
                if pos is None:
                    zmatrix_type = True

                if symbol.upper() == 'TV' and pos is not None:
                    pbc[npbc] = True
                    cell[npbc] = pos
                    npbc += 1
                else:
                    nuclear_props.append(current_nuclear_props)
                    if not zmatrix_type:
                        symbols.append(symbol)
                        positions.append(pos)

            if zmatrix_type:
                zmatrix_contents += _get_zmatrix_line(line)

    # Now that we are past the molecule spec. section, we can read
    # the entire z-matrix (if set):
    if len(positions) == 0:
        if zmatrix_type:
            if zmatrix_vars == '':
                zmatrix_vars = None
            positions, symbols = _read_zmatrix(
                zmatrix_contents, zmatrix_vars)

    try:
        atoms = Atoms(symbols, positions, pbc=pbc, cell=cell)
    except (IndexError, ValueError, KeyError) as e:
        raise IOError("ERROR: Could not read the Gaussian input file, "
                         "due to a problem with the molecule "
                         "specification: {}".format(e))

    nuclear_props = _get_nuclear_props_for_all_atoms(nuclear_props)

    return atoms, nuclear_props


def _get_readiso_param(parameters):
    ''' Returns a tuple containing the frequency
    keyword and its options, if the frequency keyword is
    present in parameters and ReadIso is one of its options'''
    freq_options = parameters.get('freq', None)
    if freq_options:
        freq_name = 'freq'
    else:
        freq_options = parameters.get('frequency', None)
        freq_name = 'frequency'
    if freq_options is not None:
        if ('readiso' or 'readisotopes') in freq_options:
            return freq_name, freq_options
    return None, None


def _get_readiso_info(line, parameters):
    '''Reads the temperature, pressure and scale from the first line
    of a ReadIso section of a Gaussian input file. Returns these in
    a dictionary.'''
    readiso_params = {}
    if _get_readiso_param(parameters)[0] is not None:
        # when count_iso is 0 we are in the line where
        # temperature, pressure, [scale] is saved
        line = line.replace(
            '[', '').replace(']', '')
        tokens = line.strip().split()
        try:
            readiso_params['temperature'] = tokens[0]
            readiso_params['pressure'] = tokens[1]
            readiso_params['scale'] = tokens[2]
        except IndexError:
            pass
    if readiso_params != {}:

        return readiso_params


def _delete_readiso_param(parameters):
    '''Removes the readiso parameter from the parameters dict'''
    parameters = deepcopy(parameters)
    freq_name, freq_options = _get_readiso_param(parameters)
    if freq_name is not None:
        if 'readisotopes' in freq_options:
            iso_name = 'readisotopes'
        else:
            iso_name = 'readiso'
        freq_options = [v.group() for v in re.finditer(
            r'[^\,/\s]+', freq_options)]
        freq_options.remove(iso_name)
        new_freq_options = ''
        for v in freq_options:
            new_freq_options += v + ' '
        if new_freq_options == '':
            new_freq_options = None
        else:
            new_freq_options = new_freq_options.strip()
        parameters[freq_name] = new_freq_options
    return parameters


def _update_readiso_params(parameters, symbols):
    ''' Deletes the ReadIso option from the route section as we
    write out the masses in the nuclear properties section
    instead of the ReadIso section.
    Ensures the masses array is the same length as the
    symbols array. This is necessary due to the way the
    ReadIso section is defined:
    The mass of each atom is listed on a separate line, in
    order of appearance in the molecule spec. A blank line
    indicates not to modify the mass for that atom.
    But you do not have to leave blank lines equal to the
    remaining atoms after you finsihed setting masses.
    E.g. if you had 10 masses and only want to set the mass
    for the first atom, you don't have to leave 9 blank lines
    after it.
    '''
    parameters = _delete_readiso_param(parameters)
    if parameters.get('isolist') is not None:
        if len(parameters['isolist']) < len(symbols):
            for i in range(0, len(symbols) - len(parameters['isolist'])):
                parameters['isolist'].append(None)
        if all(m is None for m in parameters['isolist']):
            parameters['isolist'] = None

    return parameters


def _validate_params(parameters):
    '''Checks whether all of the required parameters exist in the
    parameters dict and whether it contains any unsupported settings
    '''
    # Check for unsupported settings
    unsupported_settings = {
        "z-matrix", "modredun", "modredundant", "addredundant", "addredun",
        "readopt", "rdopt"}

    for s in unsupported_settings:
        for v in parameters.values():
            if v is not None and s in str(v):
                raise IOError(
                    "ERROR: Could not read the Gaussian input file"
                    ", as the option: {} is currently unsupported."
                    .format(s))

    for k in list(parameters.keys()):
        if "popt" in k:
            parameters["opt"] = parameters.pop(k)
            warnings.warn("The option {} is currently unsupported. "
                          "This has been replaced with {}."
                          .format("POpt", "opt"))
            return


def _get_extra_section_params(extra_section, parameters, atoms):
    ''' Takes a list of strings: extra_section, which contains
    the 'extra' lines in a gaussian input file. Also takes the parameters
    that have been read so far, and the atoms that have been read from the
    file.
    Returns an updated version of the parameters dict, containing details from
    this extra section. This may include the basis set definition or filename,
    and/or the readiso section.'''

    new_parameters = deepcopy(parameters)
    # Will store the basis set definition (if set):
    basis_set = ""
    # This will indicate whether we have a readiso section:
    readiso = _get_readiso_param(new_parameters)[0]
    # This will indicate which line of the readiso section we are reading:
    count_iso = 0
    readiso_masses = []

    for line in extra_section:
        if line.split():
            # check that the line isn't just a comment
            if line.split()[0] == '!':
                continue
            line = line.strip().split('!')[0]

        if len(line) > 0 and line[0] == '@':
            # If the name of a basis file is specified, this line
            # begins with a '@'
            new_parameters['basisfile'] = line
        elif readiso and count_iso < len(atoms.symbols) + 1:
            # The ReadIso section has 1 more line than the number of
            # symbols
            if count_iso == 0 and line != '\n':
                # The first line in the readiso section contains the
                # temperature, pressure, scale. Here we save this:
                readiso_info = _get_readiso_info(line, new_parameters)
                if readiso_info is not None:
                    new_parameters.update(readiso_info)
                # If the atom masses were set in the nuclear properties
                # section, they will be overwritten by the ReadIso
                # section
                readiso_masses = []
                count_iso += 1
            elif count_iso > 0:
                # The remaining lines in the ReadIso section are
                # the masses of the atoms
                try:
                    readiso_masses.append(float(line))
                except ValueError:
                    readiso_masses.append(None)
                count_iso += 1
        else:
            # If the rest of the file is not the ReadIso section,
            # then it must be the definition of the basis set.
            if new_parameters.get('basis', '') == 'gen' \
                    or 'gen' in new_parameters:
                if line.strip() != '':
                    basis_set += line + '\n'

    if readiso:
        new_parameters['isolist'] = readiso_masses
        new_parameters = _update_readiso_params(new_parameters, atoms.symbols)

    # Saves the basis set definition to the parameters array if
    # it has been set:
    if basis_set != '':
        new_parameters['basis_set'] = basis_set
        new_parameters['basis'] = 'gen'
        new_parameters.pop('gen', None)

    return new_parameters


def _get_gaussian_in_sections(fd):
    ''' Reads a gaussian input file and returns
    a dictionary of the sections of the file - each dictionary
    value is a list of strings - each string is a line in that
    section. '''
    # These variables indicate which section of the
    # input file is currently being read:
    route_section = False
    atoms_section = False
    atoms_saved = False
    # Here we will store the sections of the file in a dictionary,
    # as lists of strings
    gaussian_sections = {'link0': [], 'route': [],
                         'charge_mult': [], 'mol_spec': [], 'extra': []}

    for line in fd:
        line = line.strip(' ')
        link0_match = _re_link0.match(line)
        output_type_match = _re_output_type.match(line)
        chgmult_match = _re_chgmult.match(line)
        if link0_match:
            gaussian_sections['link0'].append(line)
        # The first blank line appears at the end of the route section
        # and a blank line appears at the end of the atoms section:
        elif line == '\n' and (route_section or atoms_section):
            route_section = False
            atoms_section = False
        elif output_type_match or route_section:
            route_section = True
            gaussian_sections['route'].append(line)
        elif chgmult_match:
            gaussian_sections['charge_mult'] = line
            # After the charge and multiplicty have been set, the
            # molecule specification section of the input file begins:
            atoms_section = True
        elif atoms_section:
            gaussian_sections['mol_spec'].append(line)
            atoms_saved = True
        elif atoms_saved:
            # Next we read the other sections below the molecule spec.
            # This may include the ReadIso section and the basis set
            # definition or filename
            gaussian_sections['extra'].append(line)

    return gaussian_sections


class GaussianConfiguration:

    def __init__(self, atoms, parameters):
        self.atoms = atoms.copy()
        self.parameters = deepcopy(parameters)

    def get_atoms(self):
        return self.atoms

    def get_parameters(self):
        return self.parameters

    def get_calculator(self):
        # Explicitly set parameters that must for security reasons not be
        # taken from file:
        calc = Gaussian(atoms=self.atoms, command=None, restart=None,
                        ignore_bad_restart_file=Calculator._deprecated,
                        label='Gaussian', directory='.', **self.parameters)
        return calc

    @ staticmethod
    def parse_gaussian_input(fd):
        '''Reads a gaussian input file into an atoms object and
        parameters dictionary.

        Parameters
        ----------
        fd: file-like
            Contains the contents of a  gaussian input file

        Returns
        ---------
        GaussianConfiguration
            Contains an atoms object created using the structural
            information from the input file.
            Contains a parameters dictionary, which stores any
            keywords and options found in the link-0 and route
            sections of the input file.
        '''
        # The parameters dict to attach to the calculator
        parameters = {}

        file_sections = _get_gaussian_in_sections(fd)

        # Update the parameters dictionary with the keywords and values
        # from each section of the input file:
        parameters.update(_get_all_link0_params(file_sections['link0']))
        parameters.update(_get_all_route_params(file_sections['route']))
        parameters.update(_get_charge_mult(file_sections['charge_mult']))
        atoms, nuclear_props = _get_atoms_from_molspec(
            file_sections['mol_spec'])
        parameters.update(nuclear_props)
        parameters.update(_get_extra_section_params(
            file_sections['extra'], parameters, atoms))

        _validate_params(parameters)

        return GaussianConfiguration(atoms, parameters)


def read_gaussian_in(fd, attach_calculator=False):
    '''
    Reads a gaussian input file and returns an Atoms object.

    Parameters
    ----------
    fd: file-like
        Contains the contents of a gaussian input file

    attach_calculator: bool
        When set to ``True``, a Gaussian calculator will be
        attached to the Atoms object which is returned.
        This will mean that additional information is read
        from the input file (see below for details).

    Returns
    ----------
    Atoms
        An Atoms object containing the following information that has been
        read from the input file: symbols, positions, cell.

    Notes
    ----------
    Gaussian input files can be read where the atoms' locations are set with
    cartesian coordinates or as a z-matrix. Variables may be used in the
    z-matrix definition, but currently setting constants for constraining
    a geometry optimisation is not supported. Note that the `alternative`
    z-matrix definition where two bond angles are set instead of a bond angle
    and a dihedral angle is not currently supported.

    If the parameter ``attach_calculator`` is set to ``True``, then the Atoms
    object is returned with a Gaussian calculator attached.

    This Gaussian calculator will contain a parameters dictionary which will
    contain the Link 0 commands and the options and keywords set in the route
    section of the Gaussian input file, as well as:

    • The charge and multiplicity

    • The selected level of output

    • The method, basis set and (optional) fitting basis set.

    • Basis file name/definition if set

    • Nuclear properties

    • Masses set in the nuclear properties section or in the ReadIsotopes
      section (if ``freq=ReadIso`` is set). These will be saved in an array
      with keyword ``isolist``, in the parameters dictionary. For these to
      actually be used in calculations and/or written out to input files,
      the Atoms object's masses must be manually updated to these values
      (this is not done automatically)

    If the Gaussian input file has been generated by ASE's
    ``write_gaussian_in`` method, then the basis set, method and fitting
    basis will be saved under the ``basis``, ``method`` and ``fitting_basis``
    keywords in the parameters dictionary. Otherwise they are saved as
    keywords in the parameters dict.

    Currently this does not support reading of any other sections which would
    be found below the molecule specification that have not been mentioned
    here (such as the ModRedundant section).
    '''
    gaussian_input = GaussianConfiguration.parse_gaussian_input(fd)
    atoms = gaussian_input.get_atoms()

    if attach_calculator:
        atoms.calc = gaussian_input.get_calculator()

    return atoms


# In the interest of using the same RE for both atomic positions and forces,
# we make one of the columns optional. That's because atomic positions have
# 6 columns, while forces only has 5 columns. Otherwise they are very similar.
_re_atom = re.compile(
    r'^\s*\S+\s+(\S+)\s+(?:\S+\s+)?(\S+)\s+(\S+)\s+(\S+)\s*$'
)
_re_forceblock = re.compile(r'^\s*Center\s+Atomic\s+Forces\s+\S+\s*$')
_re_l716 = re.compile(r'^\s*\(Enter .+l716.exe\)$')


def _compare_merge_configs(configs, new):
    """Append new to configs if it contains a new geometry or new data.

    Gaussian sometimes repeats a geometry, for example at the end of an
    optimization, or when a user requests vibrational frequency
    analysis in the same calculation as a geometry optimization.

    In those cases, rather than repeating the structure in the list of
    returned structures, try to merge results if doing so doesn't change
    any previously calculated values. If that's not possible, then create
    a new "image" with the new results.
    """
    if not configs:
        configs.append(new)
        return

    old = configs[-1]

    if old != new:
        configs.append(new)
        return

    oldres = old.calc.results
    newres = new.calc.results
    common_keys = set(oldres).intersection(newres)

    for key in common_keys:
        if np.any(oldres[key] != newres[key]):
            configs.append(new)
            return
    else:
        oldres.update(newres)


def read_gaussian_out(fd, index=-1):
    configs = []
    atoms = None
    energy = None
    dipole = None
    forces = None
    convergence = False
    for line in fd:
        line = line.strip()
        if line.startswith(r'1\1\GINC'):
            # We've reached the "archive" block at the bottom, stop parsing
            break

        if (line == 'Input orientation:'
                or line == 'Z-Matrix orientation:' 
                or line == "Standard orientation:"):
            if atoms is not None:
                atoms.calc = SinglePointCalculator(
                    atoms, energy=energy, dipole=dipole, forces=forces,
                )
                _compare_merge_configs(configs, atoms)
            if convergence == False:
                atoms = None
                energy = None
                dipole = None
                forces = None

            numbers = []
            positions = []
            pbc = np.zeros(3, dtype=bool)
            cell = np.zeros((3, 3))
            npbc = 0
            # skip 4 irrelevant lines
            for _ in range(4):
                fd.readline()
            while True:
                match = _re_atom.match(fd.readline())
                if match is None:
                    break
                number = int(match.group(1))
                pos = list(map(float, match.group(2, 3, 4)))
                if number == -2:
                    pbc[npbc] = True
                    cell[npbc] = pos
                    npbc += 1
                else:
                    numbers.append(max(number, 0))
                    positions.append(pos)
            atoms = Atoms(numbers, positions, pbc=pbc, cell=cell)
        elif (line == '-- Stationary point found.'):
            # for cases where convergence was achieved in one step
            convergence = True
        elif (line.startswith('Energy=')
                or line.startswith('SCF Done:')):
            # Some semi-empirical methods (Huckel, MINDO3, etc.),
            # or SCF methods (HF, DFT, etc.)
            energy = float(line.split('=')[1].split()[0].replace('D', 'e'))
            energy *= Hartree
        elif (line.startswith('E2 =') or line.startswith('E3 =')
                or line.startswith('E4(') or line.startswith('DEMP5 =')
                or line.startswith('E2(')):
            # MP{2,3,4,5} energy
            # also some double hybrid calculations, like B2PLYP
            energy = float(line.split('=')[-1].strip().replace('D', 'e'))
            energy *= Hartree
        elif line.startswith('Wavefunction amplitudes converged. E(Corr)'):
            # "correlated method" energy, e.g. CCSD
            energy = float(line.split('=')[-1].strip().replace('D', 'e'))
            energy *= Hartree
        elif _re_l716.match(line):
            # Sometimes Gaussian will print "Rotating derivatives to
            # standard orientation" after the matched line (which looks like
            # "(Enter /opt/gaussian/g16/l716.exe)", though the exact path
            # depends on where Gaussian is installed). We *skip* the dipole
            # in this case, because it might be rotated relative to the input
            # orientation (and also it is numerically different even if the
            # standard orientation is the same as the input orientation).
            line = fd.readline().strip()
            if not line.startswith('Dipole'):
                continue
            dip = line.split('=')[1].replace('D', 'e')
            tokens = dip.split()
            dipole = []
            # dipole elements can run together, depending on what method was
            # used to calculate them. First see if there is a space between
            # values.
            if len(tokens) == 3:
                dipole = list(map(float, tokens))
            elif len(dip) % 3 == 0:
                # next, check if the number of tokens is divisible by 3
                nchars = len(dip) // 3
                for i in range(3):
                    dipole.append(float(dip[nchars * i:nchars * (i + 1)]))
            else:
                # otherwise, just give up on trying to parse it.
                dipole = None
                continue
            # this dipole moment is printed in atomic units, e-Bohr
            # ASE uses e-Angstrom for dipole moments.
            dipole = np.array(dipole) * Bohr
        elif _re_forceblock.match(line):
            # skip 2 irrelevant lines
            fd.readline()
            fd.readline()
            forces = []
            while True:
                match = _re_atom.match(fd.readline())
                if match is None:
                    break
                forces.append(list(map(float, match.group(2, 3, 4))))
            forces = np.array(forces) * Hartree / Bohr
    if atoms is not None:
        atoms.calc = SinglePointCalculator(
            atoms, energy=energy, dipole=dipole, forces=forces,
        )
        _compare_merge_configs(configs, atoms)
    return configs[index]
