""" 
BIOFILTER-MODEL - Biofilm model class
=====================================

Modelled components:
    1) VOC/Oxygen diffusion in the biofilm 
    2) Monod reaction kinetics in the biofilm (with oxygen limitation)
    
(c) Matteo M. 2022

"""


import fenics as fe
import numpy as np


class BiofilmModel():
    def __init__(self):

        # Initialize default model parameters
        self.params = {}

        # Number of mesh elements
        self.params['num_el_x'] = 25            # Number of mesh elements for biofilm model

        # Diffusion coefficients
        self.params['Di'] = 2.0e-10             # VOC diffusion coefficient in the biofilm (m^2 / s)
        self.params['Do'] = 4.7e-10             # Oxygen diffusion coefficient in the biofilm (m^2 /s)

        # Biofilm
        self.params['delta'] = 37.6             # Biofilm thickness (um)
        self.params['mi'] = 0.27                # Distribution coefficient of contaminant air/biofilm
        self.params['mo'] = 34.4                # Distribution coefficient oxygen air/biofilm

        # Monod kinetic parameters
        self.params['mc'] = [
            4.17e-4,            # umax,i   Saturation specific growth rate (1/s) 
            11.03,              # km,i     Reaction rate, species-i        (g m-3)
            0.26,               # ko,i     Reaction rate oxygen, specie i  (g m-3)
            78.94,              # KI,i     Inibition constant              (g m-3)
            1.0e5,              # Xv       Biofilm density (g/m3)
            0.708,              # Y        Yield coefficient of a culture on VOC (g/g)
            0.341               # Yo       Yield coefficient of a culture on oxygen (g/g)
        ]

        self.params['As'] = 133          # Biofilm surface area per unit volume of soil (m^2 / m^3)
        self.params['alpha'] = 0.3       # Percentage coverage of the particle by the biofilm (-)
        self.params['p'] = 0.3           # Soil porosity (-)
        self.params['mi_s'] = 0.02       # Distribution coeff. air/solid (-)
        self.params['ka_l'] = 0.0        # Mass transfer coeff. liquid to solid (1/s)
        self.params['ka_g'] = 3.2e-4     # Mass transfer coeff. gas to solid (1/s)

        # Initialize solver parameters
        self.solver_params = {"nonlinear_solver": "snes",
                                "snes_solver": {"sign":"nonnegative",
                                                "maximum_iterations": 600,
                                                "report": False,
                                                "error_on_nonconvergence": False}}
        # Misc variables
        self.params['log_level'] = 30        # Default log level 30
        self.params['log_active'] = True     # Activate/Deactivate logging

    # Function for location of Dirichlet boundary
    @staticmethod
    def left_boundary(x):
        return fe.near(x[0], 0.0) # x = 0

    # Set the parameter values
    def set_params(self,params_in):
        for p in self.params:
            self.params[p] = params_in[p]

    # Set the solver parameters
    def set_solver_params(self,s_params):
        self.solver_params = s_params

    # Initialize model
    def model_init(self):

        # Manage FEniCS debug messages
        '''
        CRITICAL  = 50, // Critical errors that may lead to data corruption
        ERROR     = 40, // Errors
        WARNING   = 30, // Warnings
        INFO      = 20, // Information of general interest
        PROGRESS  = 16, // What's happening (broadly)
        TRACE     = 13, // What's happening (in detail)
        DBG       = 10  // All
        '''
        fe.set_log_active(self.params['log_active'])
        fe.set_log_level(self.params['log_level']) 


        # Create mesh and define function space
        self.mesh = fe.IntervalMesh(self.params['num_el_x'],0,self.params['delta'])
        
        # For the two concentrations ui, uo, we create a mixed space with
        # functions that represent the full system (ui,uo) as a single entity
        self.P = fe.FiniteElement('P', fe.interval, 1) 
        self.element = fe.MixedElement([self.P, self.P])
        self.V = fe.FunctionSpace(self.mesh, self.element)  
        
        # Function space for signle variables u_i, u_o
        self.V1 = fe.FunctionSpace(self.mesh,self.P)

        # Dirichlet condition
        self.c_L = fe.Constant((0.0,0.0)) # Initialize values at left boundary
        self.bc = fe.DirichletBC(self.V, self.c_L, self.left_boundary)

        # Define test functions
        self.v_i, self.v_o = fe.TestFunctions(self.V)

        # Define function for concentrations
        self.u = fe.Function(self.V)

        # Split system functions to access components
        self.u_i, self.u_o = fe.split(self.u)

        # Assign the initial conditions 
        self.ui_0 = fe.Function(self.V1)
        self.ui_0_vect = self.ui_0.vector().get_local()
        self.ui_0_vect = 0.0*np.ones_like(self.ui_0_vect)
        self.uo_0 = fe.Function(self.V1)
        self.uo_0_vect = self.uo_0.vector().get_local()
        self.uo_0_vect = 0.0*np.ones_like(self.uo_0_vect)

        self.assign_initial(self.ui_0_vect,self.uo_0_vect)

        # Define variational problem
        self.mc = self.params['mc']
        self.F = self.params['Di']*1e12*fe.inner(fe.grad(self.u_i), fe.grad(self.v_i))*fe.dx + self.mc[4]/self.mc[5] * ( self.mc[0]*self.u_i / ( self.mc[1] + self.u_i + self.u_i*self.u_i/self.mc[3]  ) * self.u_o/(self.mc[2] + self.u_o) ) * self.v_i * fe.dx \
            + self.params['Do']*1e12*fe.inner(fe.grad(self.u_o), fe.grad(self.v_o))*fe.dx + self.mc[4]/self.mc[6] * ( self.mc[0]*self.u_i / ( self.mc[1] + self.u_i + self.u_i*self.u_i/self.mc[3]  ) * self.u_o/(self.mc[2] + self.u_o) ) * self.v_o * fe.dx \


        # Define NonLinearVariationalProblem
        self.du = fe.TrialFunction(self.V)
        self.J = fe.derivative(self.F, self.u, self.du)
        self.problem = fe.NonlinearVariationalProblem(self.F, self.u, self.bc, self.J)
        self.solver  = fe.NonlinearVariationalSolver(self.problem)
        self.solver.parameters.update(self.solver_params)

        # Define Vector Function Spaces to evaluate the derivatives
        self.V_vec = fe.VectorFunctionSpace(self.mesh, "P", 1) 
        self.Q = fe.FunctionSpace(self.mesh, "DG", 0) # Space for derivative

    def assign_initial(self,ui_vect_in,uo_vect_in):
        self.ui_vect = ui_vect_in
        self.uo_vect = uo_vect_in
        self.ui_0.vector().set_local(self.ui_vect[:])
        self.uo_0.vector().set_local(self.uo_vect[:])
        fe.assign(self.u.sub(0), self.ui_0)        # Assign initial conditions from the previous step
        fe.assign(self.u.sub(1), self.uo_0)        # Assign initial conditions from the previous step
        
    def assign_boundary(self,ui_bound,uo_bound):
        self.ui_bound = ui_bound
        self.uo_bound = uo_bound
        self.c_L.assign(fe.Constant((ui_bound,uo_bound)))

    def solve(self):
        # STRATEGY 1:
        # Remove "nonnegative" constraint for small concentrations (better convergence)
        if self.ui_bound*self.params['mi'] < 1e-3:
            self.solver_params["snes_solver"]["sign"] = "default" 
            self.solver.parameters.update(self.solver_params)
            (n_iter,converged) = self.solver.solve()

        else:
            self.solver_params["snes_solver"]["sign"] = "nonnegative"
            self.solver.parameters.update(self.solver_params)
            (n_iter,converged) = self.solver.solve()
        
        # Save status:
        self.converged = converged

        # STRATEGY 2:
        # snes_solver_parameters["snes_solver"]["sign"] = "default"   # Switch back to default
        # solver.parameters.update(snes_solver_parameters)

        # (n_iter,converged) = solver.solve()

        # if not converged:   # Try to bound concentrations > 0
        #     snes_solver_parameters["snes_solver"]["sign"] = "nonnegative"
        #     solver.parameters.update(snes_solver_parameters)
        #     (n_iter,converged) = solver.solve()

        if not converged:
            print('Biofilm model not converged.')
            #print(f'Source term: {f_term}')
            print(f'Boundary concentrations: {[self.ui_bound,self.uo_bound]}')
            print(self.uo_vect)

    # Calculate derivatives at points (xi,xo). Approximation (faster but less accurate)
    def calc_grads(self, xi, xo):
        gradu_i = fe.project(fe.grad(self.u_i),self.V_vec)
        gradu_o = fe.project(fe.grad(self.u_o),self.V_vec)
        g_i = gradu_i(xi)*1e6
        g_o = gradu_o(xo)*1e6
        return (g_i,g_o)

    # Calculate derivatives at x=(0,0). Approximation (faster but less accurate)
    def calc_grads_fast(self,dummy_xi,dummy_xo): 
        d_i = self.u.sub(0).compute_vertex_values()[0:2] 
        d_o = self.u.sub(1).compute_vertex_values()[0:2]
        g_i = (d_i[1]-d_i[0])/self.mesh.coordinates()[1][0]*1e6   
        g_o = (d_o[1]-d_o[0])/self.mesh.coordinates()[1][0]*1e6
        return (g_i,g_o)
    
    def calc_grads_2(self,xi,xo): 
        u_i_dx = fe.project(self.u.sub(0).dx(0), self.Q)
        u_o_dx = fe.project(self.u.sub(1).dx(0), self.Q)
        g_i = u_i_dx(xi)*1e6
        g_o = u_i_dx(xo)*1e6
        return (g_i,g_o)

    def dofs_values(self):
        x = self.mesh.coordinates()
        ui = self.u.sub(0).compute_vertex_values() 
        uo = self.u.sub(1).compute_vertex_values()
        return (x,ui,uo)
