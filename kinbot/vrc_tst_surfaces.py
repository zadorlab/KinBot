
class VRC_TST_Surface:
    def __init__(self, fragments, pp_dist):
        self.centers = {}
        for frag in fragments:
            self.centers[frag.frag_number] = frag.pivot_points

        self.distances = pp_dist

    def __repr__(self):
        return f"Surface(pivotpoints={self.centers},\n          distances={self.distances}),".replace("[[", "np.array([[").replace("]]","]])").replace("]]),","]]),\n                    ").replace("])),","])),\n\n")
