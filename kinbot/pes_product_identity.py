"""Disambiguate colliding legacy fragment IDs during native PES assembly.

QC files retain their original names. Logical product populations use the full
saved graph/charge/spin identity when a legacy name has multiple meanings.
"""
from collections import defaultdict
from pathlib import Path
import re

from ase.db import connect


class ProductIdentities:
    def __init__(self, workers):
        observations = {}
        identities = defaultdict(set)
        for worker in workers:
            database = Path(worker) / 'kinbot.db'
            if not database.is_file():
                continue
            latest = {}
            for row in connect(database).select():
                name = row.get('name', '')
                if name.startswith('stereochemistry/') and name.endswith('_well'):
                    latest[name[len('stereochemistry/'):-len('_well')]] = row
            for name, row in latest.items():
                identity = row.data.get('identity', {})
                if identity.get('status') != 'assigned':
                    continue
                # Configured stereo names already distinguish their states.
                if not name.isdigit():
                    continue
                digest = identity['id']
                observations[str(worker), name] = digest
                identities[name].add(digest)
        collisions = {key for key, values in identities.items() if len(values) > 1}
        self.aliases = {}
        self.contexts = {}
        for (worker, name), digest in observations.items():
            if name not in collisions:
                continue
            alias = name + '-s' + digest[:16]
            value = dict(source_name=name, identity=digest)
            if alias in self.aliases and self.aliases[alias] != value:
                raise ValueError('Distinct product identities have colliding shortened aliases.')
            self.aliases[alias] = value
            self.contexts.setdefault(worker, {})[name] = alias

    def product(self, worker, name):
        return self.contexts.get(str(worker), {}).get(name, name)

    def original(self, name):
        return self.aliases.get(name, {}).get('source_name', name)

    def original_key(self, key):
        return '_'.join(sorted(self.original(part) for part in key.split('_')))

    def annotate(self, text, worker):
        mapping = self.contexts.get(str(worker), {})
        if not mapping:
            return text
        def replace(match):
            words = match[0].split()
            for i in range(2, len(words)):
                words[i] = '_'.join(sorted(mapping.get(part, part) for part in words[i].split('_')))
            return ' '.join(words)
        return re.sub(r'^! kinbot_(?:population|mirror_endpoints|racemic_keys|racemic_route) .+$',
                      replace, text, flags=re.M)

    def report(self):
        return dict(aliases=self.aliases, workers=self.contexts)
