import os
from pathlib import Path
import warnings

import numpy as np
import torch
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from util import util_path
from util.data import data
import util.nn.pes_compNet_multifid as pes
from util.sfi import daev

os.environ['KMP_DUPLICATE_LIB_OK']='True'


class Nn_surr(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, fname=None, restart=None, label='surrogate', atoms=None, 
                 tnsr=False, **kwargs):
        Calculator.__init__(self, restart=restart, label=label, atoms=atoms, 
                            tnsr=tnsr, **kwargs)
        torch.set_default_dtype(torch.float64)
        torch.set_printoptions(precision=12)
        if isinstance(fname, list) and fname.__len__() > 1:
            self.multinn = True
        elif fname is None:
            warnings.warn('No NN model provided. Falling back to C5H5. This '
                          'might lead to incorrect results.')
            self.multinn = True
            fname = [util_path + '/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand0-29540.pt',
                     util_path + '/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand1-29260.pt',
                     util_path + '/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand2-29115.pt',
                     util_path + '/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand3-0005.pt',
                     util_path + '/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand4-2500.pt']
            # fname = [util_path + '/models/new/comp-1500_rand0.pt',
            #          util_path + '/models/new/comp-1500_rand2.pt',
            #          util_path + '/models/new/comp-1500_rand4.pt',
            #          util_path + '/models/new/comp-1500_rand6.pt',
            #          util_path + '/models/new/comp-1500_rand8.pt']
        else:
            self.multinn = False
        self.surrogate = Nnpes_calc(fname, self.multinn)
        self.tnsr = tnsr

    def calculate(self, atoms=None, properties=['energy', 'forces'], system_changes=all_changes, args=None):
        Calculator.calculate(self, atoms, properties, system_changes)
        if 'forces' in properties:
            favail = True
        else:
            favail = False
        xyzd = [[[s for s in atoms.symbols], np.array(atoms.positions)]]

        self.surrogate.nforce = self.surrogate.dpes.full_symb_data[0].__len__() * 3
        self.surrogate.dpes.aev_from_xyz(xyzd, 32, 8, 8, False, self.surrogate.myaev)
        if self.multinn:
            energy, Estd = self.surrogate.eval()
            if favail:
                force, Fstd = self.surrogate.evalforce()
        else:
            energy = self.surrogate.eval()[0][0]
            Estd = torch.tensor(0.)
            if favail:
                force = self.surrogate.evalforce()
                Fstd = torch.tensor(0.)

        if self.tnsr:
            self.results['energy'] = energy
            self.results['energy_std'] = Estd
            if favail:
                self.results['forces'] = force.view(-1, 3)
                self.results['forces_std'] = Fstd
        else:
            self.results['energy'] = float(energy)
            self.results['energy_std'] = float(Estd)
            if favail:
                self.results['forces'] = np.reshape(force.detach().numpy(), (-1, 3))
                self.results['forces_std'] = Fstd.detach().numpy()


class My_args():

    def __init__(self, nntype, model_name, fid):
        self.nntype = [nntype]
        if model_name == None:
            self.load_model = False
        else:
            if isinstance(model_name, list):
                self.load_model_name = model_name
            else:
                self.load_model_name = [model_name]
        self.fid = [fid]
        self.nw = True
        self.savepth = None
        self.savepth_pars = None
        self.device = torch.device('cpu')


class Nnpes_calc():

    def __init__(self, fname, multinn=False):
        self.dpes = data.Data_pes(['C', 'H'])
        self.myaev = self.dpes.prep_aev()
        if multinn:
            self.nmodel = fname.__len__()
            options = [My_args('Comp', fnm, 'hfonly') for fnm in fname]
            self.dpes.device = options[0].device
            self.nn_pes = [pes.prep_model(True, opts, load_opt=False) for opts in options]
        else:
            self.nmodel = 1
            options = My_args('Comp', fname, 'hfonly')
            self.dpes.device = options.device
            self.nn_pes = pes.prep_model(True, options, load_opt=False)

    def eval(self, indvout=False):
        idl = list(range(0, self.dpes.ndat))
        self.dpes.prep_data(idl)
        self.dpes.indvout = indvout
        self.dpes.ymax = np.max([y[-1] for y in self.dpes.xdat])
        self.dpes.ymin = np.min([y[-1] for y in self.dpes.xdat])
        self.dpes.xb = [[xt.requires_grad_() for xt in self.dpes.xb[b]] for b in range(self.dpes.nbt)]
        if self.nmodel == 1:
            self.E_lf, self.E_hf = self.nn_pes.eval_dl(self.dpes)
            E_pred = self.E_hf
            return E_pred
        else:
            self.E_hf = torch.empty((self.dpes.ndat, self.nmodel))
            for i in range(self.nmodel):
                E_lf, E_hf = self.nn_pes[i].eval_dl(self.dpes)
                self.E_hf[:, i] = E_hf.reshape(-1)
            E_pred = torch.mean(self.E_hf)
            Estd = torch.std(self.E_hf, 1)
            return E_pred, Estd

    def evalgrad(self):
        if self.nmodel == 1:
            dEdxyz = daev.cal_dEdxyz_dl(self.dpes, self.E_hf)[0]
            return dEdxyz
        else:
            dEdxyz = torch.empty((self.dpes.ndat, self.nmodel, self.nforce))
            for i in range(self.nmodel):
                tmp = daev.cal_dEdxyz_dl(self.dpes, self.E_hf[:, i])[0]
                dEdxyz[:, i] = tmp
            gmean = torch.mean(dEdxyz, 1).reshape(-1)
            gstd  = torch.std(dEdxyz, 1).reshape(-1)
            return gmean, gstd

    def evalforce(self):
        if self.nmodel == 1:
            gradient = self.evalgrad()
            return -gradient
        else:
            gmean, gstd = self.evalgrad()
            return -gmean, gstd

# ======================================================================
def tst_anl_vs_num_forces(surr, xyz_db, d=1.e-7):
    nmolec = len(xyz_db)
    nat = [len(xyz[0]) for xyz in xyz_db]
    errmx = 0.0
    dash = '-' * 50

    for m in range(nmolec):  # loop over molecules
        atoms = Atoms(xyz_db[m][0], xyz_db[m][1])
        surr.calculate(atoms)

        # energy at nominal conditions
        E0 = surr.results['energy'].detach().numpy()

        # analytical force at nominal conditions
        F0 = surr.results['forces'].detach().numpy().flatten()
        print("E        :",E0)
        print("E_std    :",surr.results['energy_std'].detach().numpy())
        print("F_anl    :",F0)
        print("F_anl_std:",surr.results['forces_std'].detach().numpy().flatten())

        print(dash)
        print('{:>15s}{:>15s}{:>16s}'.format('F_num  ','F_anl  ','Err   '))
        print(dash)
        for i in range(3 * nat[m]):  # loop over xyz coordinates in molecule
            s = xyz_db[m][0].copy()
            x = np.copy(xyz_db[m][1])
            x[int(i / 3)][i % 3] += d
            xyz = [s, x]
            atoms = Atoms(xyz[0], xyz[1])
            surr.calculate(atoms)
            E = surr.results['energy'].detach().numpy()

            # numerical force
            Fn = -(E - E0) / d
            err = abs(Fn - F0[i])
            print('{:>15.8f}{:>15.8f}{:>16.5e}'.format(Fn,F0[i],err))

            errmx = max(errmx, err)
        print(dash)

    return errmx


# ==============================================================================================
def main():
    # set default pytorch as double precision
    torch.set_default_dtype(torch.float64)

    torch.set_printoptions(precision=12)

    home = str(Path.home())

    nn1 = home+'/mls/util/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand0-29540.pt'
    nn2 = home+'/mls/util/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand1-29260.pt'
    nn3 = home+'/mls/util/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand2-29115.pt'
    nn4 = home+'/mls/util/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand3-0005.pt'
    nn5 = home+'/mls/util/models/AL/C5H5-2b2IRC1ts_full/run1-1/C5H5_lf_reddb_2w2IRC1ts_20k_rand4-2500.pt'

    # import glob
    # fn_nn = glob.glob('../../models/AL/C5H5-2b2IRC1ts_smallset/run0/*')
    # fn_nn = ['comp-0010.pt', 'comp-0010.pt']

    fn_nn = [nn1, nn2, nn3, nn4, nn5]


    #atoms = Atoms(['C', 'H', 'C'],
    #              [(1.734285, 1.257951, 0.130304), (3.016641, 1.515035, 0.130976), (1.617565, -0.19116, -0.177416)])

    xyz_db = [[
        ['C', 'H', 'C'],
        np.array([[1.734285, 1.257951, 0.130304], [3.016641, 1.515035, 0.130976], [1.617565, -0.19116, -0.177416]])
    ]]

    # simple test
    for nn in fn_nn:
        print("nn:",nn)
        surr  = Nn_surr(nn,tnsr=False)
        for xyz in xyz_db:
            atoms = Atoms(xyz[0], xyz[1])
            atoms.set_calculator(surr)
            surr.calculate(atoms)
            print('System:',atoms)
            print('Energy:', surr.results['energy'])
            print('Forces:\n', surr.results['forces'])
            print('Numerical forces with ASE:\n', surr.calculate_numerical_forces(atoms))

    print("avg nn:")
    surr  = Nn_surr(fn_nn,tnsr=False)
    for xyz in xyz_db:
        atoms = Atoms(xyz[0], xyz[1])
        atoms.set_calculator(surr)
        surr.calculate(atoms)
        print('System:',atoms)
        print('Energy:', surr.results['energy'])
        print('Forces:\n', surr.results['forces'])
        print('Numerical forces with ASE:\n', surr.calculate_numerical_forces(atoms))

    # evaluate error between numerical and analtyical forces
    for nn in fn_nn:
        errmx = tst_anl_vs_num_forces(Nn_surr(nn), xyz_db)
        print("errmx:",errmx)
    errmx = tst_anl_vs_num_forces(Nn_surr(fn_nn), xyz_db)
    print("errmx:",errmx)

    #sys.exit()

    a = np.array([[-0.9981252206, -0.2266407005, -0.5827299231],
                 [-0.9349253557, -0.5805118019, 0.6653986964],
                 [0.2299383635, 0.7471148235, -0.9671155914],
                 [0.5032100887, -0.4581358393, 0.8855654072],
                 [1.1736620945, 0.5548549615, 0.0099222856],
                 [-1.7749905917, -0.6584807936, -1.2744784327],
                 [-1.6602525733, -1.2656303199, 1.0401661686],
                 [2.0469469305, 1.2790772568, 0.5610993450],
                 [1.1804663131, -1.1761659084, 1.4914759564],
                 [0.5184324813, 1.3801415798, -1.9407918912]])
    xyz_db = [[
        ['C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H'], a]]
    # more complex test
    xyz_db = [[
        ['C', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H', 'H'],
        np.array([[1.734285, 1.257951, 0.130304],
                  [3.016641, 1.515035, 0.130976],
                  [1.617565, -0.19116, -0.177416],
                  [1.491163, 1.573928, 1.126331],
                  [1.41661, 1.952114, -0.624898],
                  [0.231041, 0.279619, -0.018098],
                  [1.930079, -0.491148, -1.17221],
                  [2.0027, -0.86708, 0.570075],
                  [-0.30944, 0.114, 0.922074],
                  [-0.385361, 0.506492, -0.898357]])
    ], [
        ['C', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H', 'H'],
        np.array([[-0.068734, -0.040617, 0.077356],
                  [-0.096814, -0.133157, 1.168758],
                  [1.331375, -0.156964, -0.458168],
                  [-0.725694, -0.821947, -0.341217],
                  [-0.519642, 0.92543, -0.187037],
                  [1.66586, 0.331802, -1.841867],
                  [1.99176, 0.942711, -0.721455],
                  [1.984315, -0.92082, -0.043859],
                  [2.508371, -0.09514, -2.373111],
                  [0.901003, 0.842356, -2.416362]])
    ], [
        ['C', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H'],
        np.array([[1.042343, 0.333607, -0.288188],
                  [1.403663, 1.016815, 0.495818],
                  [-0.463103, 0.16199, -0.21405],
                  [1.349027, 0.754848, -1.257327],
                  [1.548637, -0.633981, -0.160537],
                  [-1.096802, -1.015782, -0.015142],
                  [-1.055864, 1.073404, -0.332644],
                  [-2.182693, -1.084731, 0.030673],
                  [-0.545227, -1.949117, 0.108606]])
    ]]

    earray = []
    farray = []

    surr = Nn_surr(fn_nn)
    for xyz in xyz_db:
        atoms = Atoms(xyz[0], xyz[1])
        surr.calculate(atoms)
        earray.append(surr.results['energy'].detach().numpy())
        farray.append(surr.results['forces'].detach().numpy().flatten())
        print('Energy:\n', surr.results['energy'])
        print('Forces:\n', surr.results['forces'])

    earray2 = []
    farray2 = []
    for j in range(fn_nn.__len__()):
        surr = Nn_surr(fn_nn[j])
        etmp = []
        ftmp = []
        for xyz in xyz_db:
            atoms = Atoms(xyz[0], xyz[1])
            surr.calculate(atoms)
            etmp.append(surr.results['energy'].detach().numpy())
            ftmp.append(surr.results['forces'].detach().numpy().flatten())
        earray2.append(etmp)
        farray2.append(ftmp)

    emn = np.mean(np.array(earray2), axis=0)
    fmn = np.mean(np.array(farray2), axis=0)

    print(np.allclose(np.array(earray).flatten(), emn.flatten()))
    for fa,fb in zip(np.array(farray),fmn.flatten()):
        print(np.allclose(fa,fb))

    errmx = tst_anl_vs_num_forces(Nn_surr(fn_nn), xyz_db)
    print("errmx:",errmx)

# ==========================================================================================

if __name__ == "__main__":
    main()


