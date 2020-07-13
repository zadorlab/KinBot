import numpy as np

class B_Mat:

    def __init__(self, species):
        self.B = np.array([])
        self.natom = species.natom
        self.geom = species.geom
        species.find_bond()
        self.bond = species.bondlist
        species.find_angle()
        self.angle = species.angle
        species.find_dihedral(findall=1)
        self.dihed = species.dihed_all
        species.find_dihedral(findall=0)
        self.dihed = species.dihed

    def bond_bij(self):
        """
        Bond stretch elements of the Wilson B matrix.
        """

        for bond in self.bond:
            brow = np.zeros(self.natom * 3)

            # a - b 
            a = bond[0]  # 1
            b = bond[1]  # 2

            # st1 = e21
            st1, r1 = self.eab(b, a)
            #st2 = e12
            st2, r1 = self.eab(a, b)

            brow[3*a: 3*a+3] = st1
            brow[3*b: 3*b+3] = st2

            self.B = np.append(self.B, brow)

    def angle_bij(self):
        """
        Bond angle elements of the Wilson B matrix.
        """

        for angle in self.angle:
            brow = np.zeros(self.natom * 3)

            # a <--- b ---> c, fi
            # 1      3      2
            a = angle[0]  # 1
            b = angle[1]  # 3 
            c = angle[2]  # 2

            e31, r31 = self.eab(a, b)
            e32, r32 = self.eab(c, b)

            fi = np.arccos(np.dot(e31, e32))

            #        cos fi e31 - e32
            # st1 = ------------------
            #         r31 sin fi

            st1 = (np.cos(fi) * e31 - e32) / (r31 * np.sin(fi))

            #       cos fi e32 - e31
            # st2 = ----------------
            #         r32 sin fi

            st2 = (np.cos(fi) * e32 - e31) / (r32 * np.sin(fi))

            #       [(r31 - r32 cos fi) e31 + (r32 - r31 cos fi) e32
            # st3 = ------------------------------------------------
            #                      r31 r32 sin fi

            st3 = ((r31 - r32 * np.cos(fi)) * e31 + (r32 - r31 * np.cos(fi)) * e32)   \
                    / (r31 * r32 * np.sin(fi))
            
            brow[3*a: 3*a+3] = st1
            brow[3*c: 3*c+3] = st2  # this is c on purpose, see ordering above
            brow[3*b: 3*b+3] = st3

            self.B = np.append(self.B, brow)

        return 0

    def dihedral_bij(self):
        """
        Dihedral angle elements of the Wilson B matrix.
        """

        for dihed in self.dihed:
            brow = np.zeros(self.natom * 3)

            # a <--- b <---> c ---> d
            # 1      2       3      4
            #        fi2     fi3

            a = dihed[0]  # 1
            b = dihed[1]  # 2 
            c = dihed[2]  # 3
            d = dihed[3]  # 4

            e12, r12 = self.eab(a, b)
            e23, r23 = self.eab(b, c)
            e43, r43 = self.eab(d, c)
            e32, r32 = self.eab(c, b)

            # angles
            fi2 = np.arccos(np.dot(e12, e32))
            fi3 = np.arccos(np.dot(e43, e23))

            #           e12 x e23
            #st1 = -  -------------
            #          r12 sin^2 fi2
            st1 = -np.cross(e12, e23) / (r12 * np.sin(fi2)**2)

            #      r23 - r12 cos fi2    e12 x e23      cos fi3     e43 x e32
            #st2 = ------------------- ----------- + ------------ -----------
            #      r23 r12 sin fi2       sin fi2     r23 sin fi3     sin fi3
            st2 = (r23 - r12 * np.cos(fi2)) / (r23 * r12 * np.sin(fi2)) * (np.cross(e12, e23)) / (np.sin(fi2))  \
                  + (np.cos(fi3)) / (r23 * np.sin(fi3)) * (np.cross(e43, e32)) / (np.sin(fi3))

            # permute (14) and (23)
            st3 = (r32 - r43 * np.cos(fi2)) / (r32 * r43 * np.sin(fi2)) * (np.cross(e43, e32)) / (np.sin(fi2))  \
                  + (np.cos(fi3)) / (r32 * np.sin(fi3)) * (np.cross(e12, e23)) / (np.sin(fi3))

            st4 = -np.cross(e43, e32) / (r43 * np.sin(fi2)**2)

            brow[3*a: 3*a+3] = st1
            brow[3*b: 3*b+3] = st2  
            brow[3*c: 3*c+3] = st3
            brow[3*d: 3*d+3] = st4

            self.B = np.append(self.B, brow)

        return 0

    def eab(self, a, b):
        e = self.geom[a] - self.geom[b]
        r = np.linalg.norm(e)
        return e / r, r

def main():
    """
    This is the Wilson B matrix. 
    """

if __name__ == "__main__":
    main()
