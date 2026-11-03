"""Inner-saddle search accompanying a proposed homolytic cleavage."""
import numpy as np

from kinbot.reac_General import GeneralReac


class DissociationSaddle(GeneralReac):
    scan = 1
    skip = 0
    family_name = 'dissociation_saddle'

    def __init__(self, species, qc, par, instance, instance_name):
        super().__init__(species, qc, par, instance, instance_name)
        self.homolysis_parent = instance_name.replace('_dissociation_saddle_', '_hom_sci_', 1)
        self.max_step = int(par['homolytic_saddle_scan_points'])
        self.start_distance = float(np.linalg.norm(
            species.geom[instance[0]] - species.geom[instance[1]]))

    def get_constraints(self, step, geom):
        change = []
        if step < self.max_step:
            # Absolute targets make a repeated admission independent of the
            # previous geometry. Start at the optimized bond, not at 2 A:
            # a short-bond inner saddle can otherwise be jumped over.
            self.set_bond(0, 1, self.start_distance +
                          step * self.par['homolytic_saddle_scan_step'], change)
        return step, [], change, []
