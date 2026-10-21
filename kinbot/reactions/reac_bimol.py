"""Base class for force-guided bimolecular reaction families."""

from kinbot.reac_General import GeneralReac
from kinbot.stationary_pt import StationaryPoint


class BimolecularReaction(GeneralReac):
    """Common handoff from a validated L0 saddle to L1 refinement."""

    bimol_search = True
    max_step = 1
    scan = 0
    skip = 0

    def get_bimol_search_data(self):
        """Return serializable family hooks and bond changes."""
        raise NotImplementedError

    def get_bimol_site(self, instance):
        """Return atoms whose equivalence classes define site degeneracy."""
        return tuple(instance)

    def carry_out_reaction(self, step):
        if step == 0:
            if not hasattr(self, 'guess_geom'):
                return -1
            status = self.qc.check_qc(self.instance_name)
            if status == 'error':
                return -1
            if status == 'normal':
                return self.max_step + 1
            if status == 'running':
                return 1
            ts = StationaryPoint(
                self.instance_name, self.species.charge, self.species.mult,
                atom=self.species.atom, geom=self.guess_geom, wellorts=1)
            ts.characterize()
            self.qc.qc_opt_ts(ts, self.guess_geom, high_level=0)
            return 1

        if step == 1:
            status = self.qc.check_qc(self.instance_name)
            if status == 'error':
                return -1
            if status == 'normal':
                return 2
        return step
