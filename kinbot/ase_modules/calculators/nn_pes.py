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

    def __init__(self, fname, restart=None, ignore_bad_restart_file=False, 
                 label='surrogate', atoms=None, tnsr=False, **kwargs):
        Calculator.__init__(self, restart=restart, 
                            ignore_bad_restart_file=ignore_bad_restart_file, 
                            label=label,
                            atoms=atoms, tnsr=tnsr, **kwargs)
        if isinstance(fname, list) and fname.__len__() > 1:
            self.multinn = True
        else:
            self.multinn = False
        self.surrogate = Nnpes_calc(fname, self.multinn)
        self.tnsr = tnsr
    
    def calculate(self, atoms=None, properties=['energy', 'forces'], 
                  system_changes=all_changes, loaddb=None, args=None, xid=None):
        Calculator.calculate(self, atoms, properties, system_changes)
        if 'forces' in properties:
            favail = True
        else:
            favail = False

        xyzd = [[[s for s in atoms.symbols], np.array(atoms.positions)]]
        self.surrogate.dpes.aev_from_xyz(xyzd, 32, 8, 8, [4.6,3.1], False, 
                                         self.surrogate.myaev)
        self.surrogate.nforce = self.surrogate.dpes.full_symb_data[0].__len__() * 3

        if self.multinn:
            energy, Estd, E_hf = self.surrogate.eval()
            if favail:
                force, Fstd, force_ind = self.surrogate.evalforce()
        else:
            energy = self.surrogate.eval()[0][0]
            Estd = torch.tensor(0.)
            if favail:
                force = self.surrogate.evalforce()
                Fstd = torch.tensor(0.)

        if self.tnsr:
            self.results['energy'] = energy
            self.results['energy_std'] = Estd
            if self.multinn:
                self.results['all_energies'] = E_hf
            if favail:
                self.results['forces'] = force.view(-1, 3)
                self.results['forces_std'] = Fstd
                if self.multinn:
                    self.results['all_forces'] = force_ind
        else:
            self.results['energy'] = float(energy)
            self.results['energy_std'] = Estd.detach().numpy()
            if self.multinn:
                self.results['all_energies'] = E_hf.detach().numpy()
            if favail:
                self.results['forces'] = np.reshape(force.detach().numpy(), 
                                                    (-1, 3))
                self.results['forces_std'] = Fstd.detach().numpy()
                if self.multinn:
                    self.results['all_forces'] = force_ind.detach().numpy()


class My_args():

    def __init__(self, nntype, model_name, fid):
        self.nntype = [nntype]
        if model_name is None:
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
            self.nn_pes = [pes.prep_model(True, opts, load_opt=False) 
                           for opts in options]
        else:
            self.nmodel = 1
            options = My_args('Comp', fname, 'hfonly')
            self.dpes.device = options.device
            self.nn_pes = pes.prep_model(True, options, load_opt=False)
            #print('SAE energies: ', self.nn_pes.sae_energies)

    def eval(self, indvout=False):
        idl = list(range(0, self.dpes.ndat))
        self.dpes.prep_data(idl)
        self.dpes.indvout = indvout
        self.dpes.ymax = np.max([y[-1] for y in self.dpes.xdat])
        self.dpes.ymin = np.min([y[-1] for y in self.dpes.xdat])
        self.dpes.xb = [[xt.requires_grad_() for xt in self.dpes.xb[b]] 
                        for b in range(self.dpes.nbt)]
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
            return E_pred, Estd, self.E_hf

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
            return gmean, gstd, dEdxyz

    def evalforce(self):
        if self.nmodel == 1:
            gradient = self.evalgrad()
            return -gradient
        else:
            gmean, gstd, dEdxyz = self.evalgrad()
            return -gmean, gstd, -dEdxyz


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
    #              [(1.734285, 1.257951, 0.130304), 
    #               (3.016641, 1.515035, 0.130976), 
    #               (1.617565, -0.19116, -0.177416)])

    xyz_db = [[
        ['C', 'H', 'C'],
        np.array([[1.734285, 1.257951, 0.130304], 
                  [3.016641, 1.515035, 0.130976], 
                  [1.617565, -0.19116, -0.177416]])
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
            print('Numerical forces with ASE:\n', 
                  surr.calculate_numerical_forces(atoms))

    print("avg nn:")
    surr  = Nn_surr(fn_nn,tnsr=False)
    for xyz in xyz_db:
        atoms = Atoms(xyz[0], xyz[1])
        atoms.set_calculator(surr)
        surr.calculate(atoms)
        print('System:',atoms)
        print('Energy:', surr.results['energy'])
        print('Forces:\n', surr.results['forces'])
        print('Numerical forces with ASE:\n', 
              surr.calculate_numerical_forces(atoms))


if __name__ == "__main__":
    main()
