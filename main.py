""" 
BIOFILTER-MODEL - MAIN FILE
===========================

    Modelled components:
    1) Air convection + disperision
    2) Contaminant diffusion in the biofilm 
    3) Monod reaction kinetics in the biofilm (with oxygen limitation)
    4) Adsorption on the solid particles
    5) Air mixing model to simulate recirculation in a confined volume (e.g., an indoor space)
    
    Implemented features:
    1) Oop implementation
    2) Parallel biofilm model execution
    3) Parallel transport model
    4) Input time series
    5) Input parameters as .yml file
    6) HDF file export

(c) Matteo M. 2022

"""


import biomod as bm
import transpmod as tm
import helpers as h
import numpy as np
import pandas as pd
from scipy import special
from scipy.integrate import solve_ivp
import multiprocessing as mp
import matplotlib.pyplot as plt
import h5py
import time
import datetime, pytz


# Directory where .yml parameter files are stored (include trailing "/")
ymldirectory = 'parameters/'


# Load parameters from file
par = h.Parameters()
ymlfiles = par.listfiles(ymldirectory)
num_yml = len(ymlfiles)


for y in range(len(ymlfiles)):

    # Read parameters file
    params = par.yamlread(ymldirectory+ymlfiles[y])

    # Air flow rate conversions (m3/s)
    params['Q'] = params['Q_m3h']/3600

    # Calculate air velocity
    params['vel'] = params['Q']/params['A']     # Air velocity (m/s)

    # Set solver params
    snes_solver_parameters = params['snes_solver_parameters']
    solver_params = params['solver_params']

    # Check num processors
    max_cpu = mp.cpu_count()
    if params['num_proc_bm'] > max_cpu:
        params['num_proc_bm'] = max_cpu

    # Build time vector
    ts1 = np.arange(0,params['t_end1']+params['dt1'],params['dt1'])
    ts2 = np.arange(params['t_end1']+params['dt2'],params['t_end2']+params['dt2'],params['dt2'])
    ts3 = np.arange(params['t_end2']+params['dt3'],params['t_end3']+params['dt3'],params['dt3'])

    params['ts'] = np.concatenate((ts1,ts2,ts3))            # Build Time vector (column)

    # Number of steps
    num_t = params['ts'].size-1  

    # Create vector for profile export
    if params['export_profiles'] or params['export_bm_profiles']:
        k_profiles = np.arange(0,num_t,params['exp_profiles_step'])
        num_profiles = k_profiles.size


    ########################################################################################
    ############################   MODEL INITIALIZATION   ##################################

    # Mixing model
    # Override params['c_L']
    # Set a very low value (to initialize ESDIRK correctly e.g. 1.0e-7 -> 0.1 ug/m3 
    if params['mixing_model']:
        params['c_L'] = 1.0e-7

    ##########################  TIME SERIES ################################

    # Input time-series ppb (biofilter input)
    if params['use_ts']:
        d_in = h.FileIO(params['csv_file'])
        d_in.mw = params['mw']                                                                  # Set molar weight for conversions
        d_in.c_name = 'tVOC'                                                                    # Column name
        data_in = d_in.readfile(smooth=params['smooth'], n_samples=params['smooth_samples'])    # Return smoothed input ts
        data_ts = d_in.interp(params['ts'])*h.Conversion.ppbtog(params['mw'])                   # Return interpolated ts

        # Validation time-series (biofilter output)
        d_out = h.FileIO(params['csv_file_out'])
        d_out.mw = params['mw'] 
        data_out = d_out.readfile(smooth=params['smooth'], n_samples=params['smooth_samples'])

    # Mixing model time-series 
    if params['mixing_model']:
        if params['use_ts_mixing_ext']:
            # Load time series
            Qext_h = h.FileIO(params['csv_file_mixing_Qext'])
            Qext_h.dateformat = False
            Qext_h.c_name = 'Qext'
            Qext = Qext_h.readfile(smooth=False, n_samples=0)
            Qext_ts = Qext_h.interp(params['ts'])/3600                      # m3/s

            Cext_h = h.FileIO(params['csv_file_mixing_Cext'])
            Cext_h.dateformat = False
            Cext_h.c_name = 'Cext'
            Cext = Cext_h.readfile(smooth=False, n_samples=0)
            Cext_ts = Cext_h.interp(params['ts'])*1.0e-6                    # ug/m3 to g/m3
        else:
            # Create constant vector
            Qext_ts = params['Qext_m3h']*np.ones(params['ts'].shape)/3600   # m3/s
            Cext_ts = params['Cext']*np.ones(params['ts'].shape)*1.0e-6     # ug/m3 to g/m3


        if params['use_ts_mixing_F']:
            # Load time series
            F_h = h.FileIO(params['csv_file_mixing_F'])
            F_h.dateformat = False
            F_h.c_name = 'F'
            F = F_h.readfile(smooth=params['smooth_F'], n_samples=params['smooth_F_samples'])
            F_ts = F_h.interp(params['ts'])/1000/3600                       # mg/h to g/s
        else:
            # Create a constant vector
            F_ts = params['F']*np.ones(params['ts'].shape)/1000/3600        # g/s

    ################## INITIALIZE TRANSPORT MODEL ##########################

    # Transport model instances
    tms = [tm.TransportModel() for i in range(2)]     # 2 species (VOC, Oxygen)

    # Model initialization
    for i in range(len(tms)):
        tms[i].set_params(params)
        tms[i].set_solver_params(solver_params)
        tms[i].model_init()

    # Biofilm model mapping
    if params['export_bm_profiles']:    
        k_dof = tms[0].find_dof(params['z_profile'])


    # Function for source term calculation
    params['q'] = tms[0].params['q']     # Get q from transport model
    def source_term(q_star,q):
        f = np.zeros(2)
        f[0] = -params['q']* ( -params['alpha']*params['As']*params['Di']*dudx_i_k + (1-params['alpha'])*params['ka_g']*(q_star-q) )
        f[1] = params['q']* ( params['alpha']*params['As']*params['Do']*dudx_o_k )
        return f

    # Number of DOFs of transport model
    num_z = tms[0].c_vect.size 

    ################## INITIALIZE BIOFILM MODEL ##########################
                
    # Create num_el_z instances of the biofilm model
    bms = [bm.BiofilmModel() for i in range(num_z)]

    # Initilize all models
    for i in range(len(bms)):
        bms[i].set_params(params)
        bms[i].set_solver_params(snes_solver_parameters)
        bms[i].model_init()

    # Number of DOFs biofilm model
    num_x = bms[0].ui_0.vector().get_local().size  


    ################## ADSORPTION MODEL ##########################

    def ads_model(t,q,q_star):
        dqdt = params['alpha'] * params['ka_l'] * (q_star - q) + (1-params['alpha'])*params['ka_g']*(q_star - q)
        return dqdt

    def ads_solve(q_init,q_star):
        q = solve_ivp(fun=lambda t, q: ads_model(t, q, q_star), t_span = t_span, y0 = [q_init], method='RK45', t_eval = [t_span[1]])   # Time integration Runge-Kutta 4-5, take only the final value of the simulation
        return float(q.y)


    ##################   MIXING MODEL ############################
    if params['mixing_model']:
        def mixing_model(t,Cmix,Cout,k):
            dCdt = 1/params['V'] * (Qext_ts[k]*Cext_ts[k] + F_ts[k] + params['Q']*Cout - (params['Q'] + Qext_ts[k])*Cmix)
            return dCdt

        def mixing_solve(Cinit,Cout,k):
            Cmix = solve_ivp(fun=lambda t, Cmix: mixing_model(t, Cmix, Cout, k), t_span = t_span, y0 = [Cinit], method='RK45', t_eval = [t_span[1]])   # Time integration Runge-Kutta 4-5, take only the final value of the simulation
            return float(Cmix.y)


    ################## SCREEN OUTPUT MANAGEMENT ##########################
    def print_start():
        dt_string = datetime.datetime.now(pytz.timezone('Europe/Rome')).strftime("%d/%m/%Y %H:%M:%S")
        print('-----------------------------------------------------------')
        print(f'Simulation {y+1}/{num_yml} started - {dt_string}')
        print(f'Parameters: {ymlfiles[y]}')
        print(f'Time-steps: {num_t}')
        print('-----------------------------------------------------------')



    #################################################################################
    ###############################   MAIN LOOP    ##################################
    start_time = time.time()
    print_start()
          

    # Initialize variables
    c_vect = np.array([np.zeros(num_z), np.zeros(num_z)])
    fvect_c = np.array([np.zeros(num_z), np.zeros(num_z)])
    q_star_vect = np.zeros(num_z)           # Adsorbed concentrations q, q*
    q_init_vect = np.zeros(num_z)   
    ui_vect = np.zeros((num_z,num_x))       # Biofilm model vectors
    uo_vect = np.zeros((num_z,num_x))
    dudx_i = np.zeros(num_z)                # Derivative of biofilm profiles at x = 0
    dudx_o = np.zeros(num_z)
    c_out = np.zeros(num_t)                 # Output concentration at z = L
    if params['mixing_model']:
        c_mix = np.zeros(num_t)
    
    if params['export_profiles']:
        profile_z = np.zeros(num_z)
        profile_t = np.zeros(num_profiles)
        profile_i = np.zeros((num_profiles,num_z))
        profile_o = np.zeros((num_profiles,num_z))

    if params['export_bm_profiles']:
        profile_x = np.zeros(num_x)
        profile_t = np.zeros(num_profiles)
        profile_ui = np.zeros((num_profiles,num_x))
        profile_uo = np.zeros((num_profiles,num_x))


    # Adsorption/Biofilm model step definition - to be solved in parallel
    def biomod_step(i):
        ############ ADSORPTION MODEL ##############
        q_init = q_init_vect[i]         # Initial condition
        q_star = q_star_vect[i]         # q*
        q = ads_solve(q_init,q_star)    # Run model

        ############ BIOFILM MODEL #################
        # Assign initial conditions
        bms[i].assign_initial(ui_vect[i],uo_vect[i])
        
        # Assign boundary condition
        bms[i].assign_boundary(c_vect[0][i]/params['mi'],c_vect[1][i]/params['mo']) 

        # Solve
        bms[i].solve()

        if not bms[i].converged:
            print(f"Num cell: {i}")

        # Calculate gradients
        (dudx_i_k,dudx_o_k) = bms[i].calc_grads(0.0,0.0)

        # Save variable status for the next step
        ui_sub = bms[i].u.sub(0, deepcopy=True)
        uo_sub = bms[i].u.sub(1, deepcopy=True)
        ui_vect_end = ui_sub.vector().get_local()
        uo_vect_end = uo_sub.vector().get_local()

        # Export profiles
        (x,ui,uo) = (0,0,0)
        if params['export_bm_profiles'] and i == k_dof:
            if k in k_profiles:
                (x,ui,uo) = bms[i].dofs_values()         

        return (q,dudx_i_k,dudx_o_k,ui_vect_end,uo_vect_end,x,ui,uo)


    # Transport model i-th species - to be solved in parallel
    def transport_i(i,t_span,k):
        # Set time step
        tms[i].T = [0.0, t_span[1]-t_span[0]]

        # Assign boundary conditions
        if i == 0: 
            if (not params['use_ts']) or params['mixing_model']:
                tms[i].assign_boudary(params['c_L'])
            else:
                tms[i].assign_boudary(data_ts[k])
            
        else:
            tms[i].assign_boudary(params['co_L'])

        # Assign initial conditions and source term
        tms[i].assign_initial(c_vect[i],fvect_c[i])

        # Solve 
        tms[i].solve()

        # Set initial conditions for the next step
        c_vect_i = tms[i].rk.u.vector().get_local()
 
        # Save variables for post-processing
        c_out_i = tms[i].rk.u(params['L'])  # Concentration at the outlet of the column

        # Algorithm termination
        term = tms[i].rk.terminateReason

        # Export column profiles
        (z,c) = (0,0)
        if params['export_profiles']:
            if k in k_profiles:
                (z,c) = tms[i].dofs_values()
        
        return (c_vect_i,c_out_i,term,z,c)

    ############ MAIN LOOP ##############
    for k in range(0,params['ts'].size-1):
        # Measure step elasped time
        step_time = time.time()

        # Current t span
        t_span = [params['ts'][k], params['ts'][k+1] ]
        

        ############ ADSORPTION + BIOFILM PARALLEL LOOP ##############

        # Parallel - Async method 
        pool = mp.Pool(processes=params['num_proc_bm'])
        output = [pool.apply_async(biomod_step, args=(i,)) for i in range(num_z)]
        results = [p.get() for p in output]
        pool.close()
        pool.join()

        # Prepare data for the next step
        for i in range(num_z):
            q_init_vect[i] = q = results[i][0]
            dudx_i_k = results[i][1]
            dudx_o_k = results[i][2]
            
            ui_vect[i][:] = results[i][3]
            uo_vect[i][:] = results[i][4]

            # Calculate source term
            q_star = q_star_vect[i]
            f_term = source_term(q_star,q)
            
            # Assign f to vectors
            fvect_c[0][i] = f_term[0]
            fvect_c[1][i] = f_term[1]

            if params['export_bm_profiles']  and i == k_dof:
                x_exp = results[i][5]
                ui_exp = results[i][6]
                uo_exp = results[i][7]

        
        ################# TRANSPORT MODEL ##############

        pool = mp.Pool(processes=params['num_proc_tm'])
        output = [pool.apply_async(transport_i, args=(i,t_span,k)) for i in range(2)]
        results_t = [p.get() for p in output]
        pool.close()
        pool.join()

    
        c_vect[0][:] = results_t[0][0]
        c_vect[1][:] = results_t[1][0]
        c_out[k] = results_t[0][1]

        q_star_vect = c_vect[0]/params['mi_s']

        # Export profiles
        if params['export_profiles']:
            if k in k_profiles:
                tmp_p = np.where(k_profiles==k)
                k_p = tmp_p[0][0]
                profile_z = results_t[0][3]
                profile_t[k_p] = params['ts'][k+1]
                profile_i[k_p] = results_t[0][4]
                profile_o[k_p] = results_t[1][4]
        
        if params['export_bm_profiles']:
            if k in k_profiles:
                tmp_p = np.where(k_profiles==k)
                k_p = tmp_p[0][0]
                profile_t[k_p] = params['ts'][k+1]
                profile_x = x_exp
                profile_ui[k_p] = ui_exp
                profile_uo[k_p] = uo_exp
        
        ############### MIXING MODEL #####################
        if params['mixing_model']:
            Cout = c_out[k]
            Cmix = mixing_solve(params['c_L'],Cout,k)    # Run mixing model
            params['c_L'] = Cmix                         # Set biofilter initial condition for the next step
            c_mix[k] = Cmix

  
        # Print coupling step info
        if ( (results_t[0][2] == "Success") and (results_t[1][2] == "Success")):
            term = True
        else:
            term = False
          
        # Calculate step time
        elapsed_step_time = (time.time() - step_time)

        print(f"Step {k}. Time-step converged: {term}.\tCPU-time: {elapsed_step_time:.3f} s.\tCompleted: {(k+1)/num_t*100:.1f}%", end = '')


    # Print total elasped time
    elapsed_time = (time.time() - start_time)/60
    print('')
    print(f'Total elapsed time: {elapsed_time:.1f} min')

    #################################################################################
    #############################   POST-PROCESSING  ################################

    ######### PLOTS ##########
    if params['plots']:
        plt.figure()
        
        if params['use_ts']:
            plt.plot(data_out['t_s']/3600,data_out['tVOC_sm'],linestyle='none',marker='o',fillstyle='none', label='Data')
        
        plt.plot(params['ts'][1:]/3600,c_out*h.Conversion.gtoppb(params['mw']),linestyle='-',marker='',fillstyle='none', label='Model')
        
        plt.legend()
        plt.xlabel('Time (h)')
        plt.ylabel('Concentration (ppb)')
        plt.show()

    ###### EXPORT ############
    if params['export_h5']:
        
        date = datetime.datetime.now(pytz.timezone('Europe/Rome')).strftime("%Y%m%d_%H%M")
        curfile = f'{ymlfiles[y][:-4]}'.replace(' ','')
        filepath = f"{params['output_folder']}{curfile}-{date}.h5"
        
        with h5py.File(filepath, 'w') as hf:    
            hf.create_group('parameters')
            hf.create_group('data/inlet')
            hf.create_group('data/outlet')
            hf.create_group('data/mixing')
            hf.create_group('results/output')
            hf.create_group('results/profiles')

            # PARAMETERS - YAML FILE
            dth5 = h5py.string_dtype()
            dset = hf.create_dataset("parameters/use_ts", data=params['use_ts'])
            dset = hf.create_dataset("parameters/export_profiles", data=params['export_profiles'])
            dset = hf.create_dataset("parameters/export_bm_profiles", data=params['export_bm_profiles'])
            dset = hf.create_dataset("parameters/mixing_model", data=params['mixing_model'])
            dset = hf.create_dataset("parameters/filename", data=ymlfiles[y], dtype=dth5)
            dset = hf.create_dataset("parameters/yaml", data=par.yamlraw(ymldirectory+ymlfiles[y]), dtype=dth5)
            dset = hf.create_dataset("parameters/mw", data=params['mw'], dtype='float32')
            dset = hf.create_dataset("parameters/use_ts_mixing_ext", data=params['use_ts_mixing_ext'])
            dset = hf.create_dataset("parameters/use_ts_mixing_F", data=params['use_ts_mixing_F'])

            # DATA
            if params['use_ts']:
                # DATA - OUTLET
                dset = hf.create_dataset("data/outlet/t", data=np.array(data_out['t_s'],dtype='int64'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 's'
                dset.attrs['Description'] = 'Time'

                if params['smooth']:
                    dset = hf.create_dataset("data/outlet/c", data=np.array(data_out['Smoothed'],dtype='float32'), compression="gzip", compression_opts=9)
                else:
                    dset = hf.create_dataset("data/outlet/c", data=np.array(data_out['tVOC'],dtype='float32'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'ppb'
                dset.attrs['Description'] = 'Column outlet concentration (ppb) time-series'

                # DATA - INLET        
                dset = hf.create_dataset("data/inlet/t", data=np.array(params['ts'],dtype='int64'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 's'
                dset.attrs['Description'] = 'Time'

                dset = hf.create_dataset("data/inlet/c", data=np.array(data_ts*h.Conversion.gtoppb(params['mw']),dtype='float32'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'ppb'
                dset.attrs['Description'] = 'Column inlet concentration (ppb) time-series (interpolated)'

            # MIXING MODEL
            if params['use_ts_mixing_ext'] or params['use_ts_mixing_F']:
                # Time
                dset = hf.create_dataset("data/mixing/t", data=np.array(params['ts'],dtype='int64'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 's'
                dset.attrs['Description'] = 'Time'

            if params['use_ts_mixing_ext']:
                # Qext
                dset = hf.create_dataset("data/mixing/Qext", data=np.array(Qext_ts*3600,dtype='float32'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'm3/h'
                dset.attrs['Description'] = 'External air flow time-series (interpolated)'
                # Cext
                dset = hf.create_dataset("data/mixing/Cext", data=np.array(Cext_ts*h.Conversion.gtoppb(params['mw']),dtype='float32'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'ppb'
                dset.attrs['Description'] = 'External concentration time-series (interpolated)'

            if params['use_ts_mixing_F']:
                # F
                dset = hf.create_dataset("data/mixing/F", data=np.array(F_ts*3600*1000,dtype='float32'), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'mg/h'
                dset.attrs['Description'] = 'Indoor emission time-series interpolated'

            # RESULTS - C OUT
            c_out1 = np.insert(c_out, 0, 0, axis=0)     # Add zero at the beginning
            dset = hf.create_dataset("results/output/c", data=c_out1*h.Conversion.gtoppb(params['mw']), compression="gzip", compression_opts=9)
            dset.attrs['Unit'] = 'ppb'
            dset.attrs['Description'] = 'VOC concentration (ppb) at the outlet of the biofilter'      
            dset = hf.create_dataset("results/output/t", data=params['ts'])
            dset.attrs['Unit'] = 's'
            dset.attrs['Description'] = 'Time'
            if params['mixing_model']:
                c_mix1 = np.insert(c_mix, 0, 0, axis=0)     # Add zero at the beginning
                dset = hf.create_dataset("results/output/c_mix", data=c_mix1*h.Conversion.gtoppb(params['mw']), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'ppb'
                dset.attrs['Description'] = 'VOC concentration (ppb) in the mixed air'  

            # PROFILES
            if params['export_profiles'] or params['export_bm_profiles']:
                dset = hf.create_dataset("results/profiles/t", data=profile_t)
                dset.attrs['Unit'] = 's'
                dset.attrs['Description'] = 'Time'

            if params['export_profiles']:
                dset = hf.create_dataset("results/profiles/c", data=profile_i*h.Conversion.gtoppb(params['mw']), compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'ppb'
                dset.attrs['Description'] = 'VOC concentration (ppb) profiles. Rows: time, columns: z'      
                dset = hf.create_dataset("results/profiles/o", data=profile_o*h.Conversion.gtoppb(32.0)/1000, compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'ppm'
                dset.attrs['Description'] = 'Oxygen concentration (ppm) profiles. Rows: time, columns: z'     
                dset = hf.create_dataset("results/profiles/z", data=profile_z)
                dset.attrs['Unit'] = 'm'
                dset.attrs['Description'] = 'Space'
            
            if params['export_bm_profiles']:
                dset = hf.create_dataset("results/profiles/ui", data=profile_ui, compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'g/m3'
                dset.attrs['Description'] = 'Biofilm VOC concentration (g/m3) profiles. Rows: time, columns: x'      
                dset = hf.create_dataset("results/profiles/uo", data=profile_uo, compression="gzip", compression_opts=9)
                dset.attrs['Unit'] = 'g/m3'
                dset.attrs['Description'] = 'Biofilm Oxygen concentration (g/m3) profiles. Rows: time, columns: x'  
                dset = hf.create_dataset("results/profiles/x", data=profile_x)
                dset.attrs['Unit'] = 'um'
                dset.attrs['Description'] = 'Space (biofilm)'
