import os
import sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.append('../../../../SpecEIS')
import firedrake as df
import rasterio as rio
from firedrake.petsc import PETSc
from speceis_dg.hybrid import CoupledModel
import matplotlib.pyplot as plt
import numpy as np

class Hubbard:

    def hubbard_bed(self):
        v_dg = df.VectorFunctionSpace(self.model.mesh,self.model.E_thk)
        X = df.interpolate(self.model.mesh.coordinates,v_dg)
        meshx = X.dat.data_ro[:,0]*1000
        meshy = X.dat.data_ro[:,1]*1000
        with rio.open("./hubbard_bedrock.tif") as src:
            self.model.B.dat.data[:] = np.array([pnt[0] for pnt in src.sample(zip(meshx, meshy))])/1000
        self.model.H0.interpolate(df.Constant(0.001))

    def __init__(self,results_dir,data_dir,conservation_test=False,init_dir=None):
        if conservation_test:
            with df.CheckpointFile(f"{init_dir}/functions.h5", 'r') as afile:
                mesh = afile.load_mesh('mesh')
        else:    
            mesh = df.Mesh(f'hubbard_mesh_firedrake.msh',name='mesh')

        config = {'solver_type': 'gmres',
                  'velocity_function_space':'MTW',
                  'sliding_law': 'Budd',
                  'vel_scale': 100.,
                  'thk_scale': 1000.,
                  'len_scale': 1000.,
                  'beta_scale': 1000.,
                  'theta': 1.0,
                  'thklim': 1e-3,
                  'alpha': 1000.0,
                  'z_sea': -8,
                  'calve': True,
                  'ssa': False}
          
        model = self.model = CoupledModel(mesh,**config)
        
        #self.interpolate_bed_from_pickle(f'{data_dir}/interpolant.pkl')
        self.hubbard_bed()

        if conservation_test:
            with df.CheckpointFile(f"{init_dir}/functions.h5", 'r') as afile:
                H_in = afile.load_function(mesh, "H0", idx=399)
                model.H0.assign(H_in)
         
        model.beta2.interpolate(df.Constant(400.0)) # higher means more friction 

        z_ela = 2

        if conservation_test:
            lapse_rate=0.0
            time_step_factor = 1.01
        else:
            lapse_rate = 0/1000
            #time_step_factor = 1.05

        # For thickness forcing
        v_dg = df.VectorFunctionSpace(self.model.mesh,self.model.E_thk)
        X = df.interpolate(self.model.mesh.coordinates,v_dg)
        x = X.dat.data_ro[:,0]
        y = X.dat.data_ro[:,1]

        xvs = 793332/1e3
        yvs = 1214481/1e3
        xhs = 808883/1e3
        yhs = 1219025/1e3

        # Valerie
        mask = np.sqrt((x-xvs)**2 + (y-yvs)**2) < 1
        model.adot.dat.data[mask] = 0#.05 * min(1, dt/50)

        # Hubbard
        mask = np.sqrt((x-xhs)**2 + (y-yhs)**2) < 2
        model.adot.dat.data[mask] = 0#.075 * min(1, dt/50)

        #model.adot.dat.data[:] = (((model.B.dat.data[:] + model.H0.dat.data[:]) 
        #                            - z_ela)*lapse_rate)

        S_file = df.File(f'{results_dir}/S.pvd')
        B_file = df.File(f'{results_dir}/B.pvd')
        Us_file = df.File(f'{results_dir}/U_s.pvd')
        H_file = df.File(f'{results_dir}/H.pvd')
        N_file = df.File(f'{results_dir}/N.pvd')
        adot_file = df.File(f'{results_dir}/adot.pvd')

        Q_cg2 = df.VectorFunctionSpace(mesh,"CG",3)
        S_out = df.Function(model.Q_thk,name='S')
        N_out = df.Function(model.Q_thk,name='N')
        U_s = df.Function(Q_cg2,name='U_s')

        S_out.interpolate(model.S)
        #N_out.interpolate(model.N)
        U_s.interpolate(model.Ubar0 - 1./4*model.Udef0)

        S_file.write(S_out,time=0.)
        H_file.write(model.H0,time=0.)
        B_file.write(model.B,time=0.)
        Us_file.write(U_s,time=0.)
        adot_file.write(model.adot,time=0.)

        t = 0.0
        t_grow = 300
        t_end = 3000
        dt = 10
        max_step = 20.0

        with df.CheckpointFile(f"{results_dir}/functions.h5", 'w') as afile:

            afile.save_mesh(mesh)

            i = 0
            while t<t_end:
                #dt = min(dt*time_step_factor,max_step)

                #model.adot.dat.data[:] = (((model.B.dat.data[:] + model.H0.dat.data[:])
                #                            - z_ela)*lapse_rate)

                # Valerie
                mask = np.sqrt((x-xvs)**2 + (y-yvs)**2) < 1
                model.adot.dat.data[mask] = .05 * min(1, t/150)

                # Hubbard
                mask = np.sqrt((x-xhs)**2 + (y-yhs)**2) < 2
                model.adot.dat.data[mask] = .075 * min(1, t/150)

                converged = model.step(t,
                                       dt,
                                       picard_tol=2e-3,
                                       momentum=0.5,
                                       max_iter=20,
                                       convergence_norm='l2')

                if not converged:
                    dt*=0.5
                    continue
                t += dt
                PETSc.Sys.Print(t,dt, round(df.assemble(model.H0*df.dx), 2))
                S_out.interpolate(model.S)
                #N_out.interpolate(model.N)
                U_s.interpolate(model.Ubar0 - 1./4*model.Udef0)

                afile.save_function(model.H0, idx=i)
                afile.save_function(S_out, idx=i)
                afile.save_function(N_out, idx=i)
                afile.save_function(U_s, idx=i)

                S_file.write(S_out,time=t)
                H_file.write(model.H0,time=t)
                B_file.write(model.B,time=t)
                Us_file.write(U_s,time=t)
                N_file.write(N_out,time=t)
                adot_file.write(model.adot,time=t)
                i += 1

if __name__=='__main__':
    bc = Hubbard('./results/','./data/',conservation_test=False,init_dir='./results/')
