"""Historical source reproduction only; never used by production drivers."""
from kinbot.anl.literature import load_reference


def source_reference_spec(name, *, reference, task_id, max_nodes=1, partition=None,
                                      mrcc_command='dmrcc'):
    """Prepare one explicitly requested UHF source-reproduction task.

    Production ANL generation remains RHF/ROHF-only. This opt-in tool requires
    the supplied reference to declare UHF for the selected task.
    """
    from kinbot.anl.validation import higher_order_validation_spec

    reference = load_reference(reference)
    benchmark = reference['benchmarks'][name]
    if benchmark['tasks'][task_id].get('request', {}).get('reference') != 'UHF':
        raise ValueError('Reference task must explicitly declare UHF.')
    spec = higher_order_validation_spec(
        benchmark['molecule'], max_nodes=max_nodes, partition=partition,
        mrcc_command=mrcc_command, task_ids=(task_id,))
    task = spec['tasks'][0]
    task['input_template'] = task['input_template'].replace(
        'scftype=ROHF\n'
        'rohftype=semicanonical\n'
        'rohfcore=semicanonical\n',
        'scftype=UHF\n')
    task['result_parser']['reference'] = 'UHF'
    task['source_reference_validation'] = f'{name}-UHF'
    spec['intent'].update({
        'claim': 'pinned-literature-source-reference-validation-only',
        'source_reference_validation': f'{name}-UHF',
        'literature_benchmark': {
            'name': name,
            'source': reference['source'],
            'worksheet': benchmark['worksheet'],
            'geometry_convention': benchmark['geometry'],
        },
    })
    return spec
