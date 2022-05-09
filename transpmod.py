""" 
BIOFILTER-MODEL - Transport model class
=======================================

Modelled components:
    1) VOC advection
    2) VOC dispersion
    3) Coupling with biofilm model
    
(c) Matteo M. 2022

"""

import fenics as fe 
import numpy as np 
import gryphon.ESDIRK as gryphon


class TransportModel():
    def __init__(self):
        
        # Initialize default model parameters below
        # Parameters will be overridden by those in the .yml configuration file
        self.params = {}

        ################## TRANSPORT PARAMS ##########################
        self.params['L'] = 0.686                 # Column length (m)
        self.params['num_el_z'] = 100            # Number of mesh elements for advection-dispersion
        self.params['DL'] = 5e-4                 # Dispersion coefficient
        self.params['vel'] = 1.81e-3             # Air velocity (m/s)
        self.params['As'] = 133                  # Biofilm surface area per unit volume of soil (m^2 / m^3)
        self.params['alpha'] = 0.3               # Percentage coverage of the particle by the biofilm (-)
        self.params['p'] = 0.3                   # Soil porosity (-)
        self.params['mi_s'] = 0.02               # Distribution coeff. air/solid (-)
        self.params['ka_l'] = 0.0                # Mass transfer coeff. liquid to solid (1/s)
        self.params['ka_g'] = 3.2e-4             # Mass transfer coeff. gas to solid (1/s)

        # INITIAL VALUES AND BOUNDARY VALUES
        self.params['c_L'] = 2.81                # Inlet concentration (Left boundary condition) (g/m3)
        self.params['c_0']  = 0.0                # Initial concentration in the column (g/m3)
        
        ################## RK SOLVER PARAMS ##########################
        self.solver_params = {}
        self.solver_params['timestepping'] = {}
        self.solver_params['output'] = {}
        self.solver_params['verbose'] = False
        self.solver_params['drawplot'] = False 
        self.solver_params['output']['statistics'] = False          # Write Output file
        self.solver_params['output']['print_statistics'] = False    # Print statistics on screen
        self.solver_params['timestepping']['adaptive'] = True
        self.solver_params['timestepping']['absolute_tolerance'] = 1e-5
        self.solver_params['timestepping']['convergence_criterion'] = "absolute"     # "absolute", "relative"

        # Misc variables
        self.params['log_level'] = 30        # Default log level 30
        self.params['log_active'] = True     # Activate/Deactivate logging

        self.T = [0.0, 1.0]                  # Time step for RK integration


    # Function for location of Dirichlet boundary
    @staticmethod
    def left_boundary(x):
        return fe.near(x[0], 0.0) # x = 0

    # RK solver instances
    def rk_instance(self):
        self.rk = gryphon.ESDIRK(self.T, self.c, self.F_c, bcs=[self.bc_c]) 
        self.set_rk_params()


    # Model initialization
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

        # Define q = (1-p)/p
        self.params['q'] = (1-self.params['p'])/self.params['p']
        
        # Mesh
        self.mesh = fe.IntervalMesh(self.params['num_el_z'],0,self.params['L']) # Create equally spaced mesh of lenght L
        # Function space
        self.V = fe.FunctionSpace(self.mesh, 'P', 1)   

        # Default boundary condition
        self.c_L = fe.Constant(self.params['c_L'])

        # Dirichlet conditions
        self.bc_c = fe.DirichletBC(self.V, self.c_L, self.left_boundary)

        # Constant velocity field
        self.W = fe.VectorFunctionSpace(self.mesh, 'P', 1)
        self.w = fe.Function(self.W)
        self.w.vector()[:] = self.params['vel']

        # Source term
        self.f_c = fe.Function(self.V)
        self.fvect_c = np.ones_like( self.f_c.vector().get_local() ) # Create a vector from f 
        self.f_c.vector().set_local(0.0*self.fvect_c) # Set f = 0.0 everywhere, for initialization 


        # Trial and test functions
        self.c = fe.TrialFunction(self.V)
        self.v = fe.TestFunction(self.V)
 
        # Define variational problem (rhs)
        self.F_c = -fe.inner(self.w,fe.grad(self.c))*self.v*fe.dx - self.params['DL']*fe.inner(fe.grad(self.c), fe.grad(self.v))*fe.dx + self.f_c*self.v*fe.dx

        # Define initial conditions for transport model
        self.c = fe.Function(self.V)
        self.c_vect = self.params['c_0']*np.ones_like( self.c.vector().get_local() )
        self.c.vector().set_local(self.c_vect)

        # Initialize RK
        self.rk_instance()

    # Set the parameter values
    def set_params(self,params_in):
        for p in self.params:
            self.params[p] = params_in[p]

    def assign_boudary(self,c_L_in):
        self.c_L.assign(fe.Constant(c_L_in))

    def assign_initial(self,c_vect_in,fvect_c_in):
        self.c_vect = c_vect_in
        self.fvect_c = fvect_c_in
        self.c.vector().set_local(self.c_vect)
        self.f_c.vector().set_local(self.fvect_c)

    def set_solver_params(self,s_params):
        for p in self.solver_params:
            self.solver_params[p] = s_params[p]
    
    
    def set_rk_params(self):
        # Solver settings
        self.rk.parameters["timestepping"]["adaptive"] = self.solver_params["timestepping"]["adaptive"]
        self.rk.parameters["timestepping"]["absolute_tolerance"] = self.solver_params["timestepping"]["absolute_tolerance"]
        self.rk.parameters["timestepping"]["convergence_criterion"] = self.solver_params["timestepping"]["convergence_criterion"]     
        if 'dtmax' in self.solver_params["timestepping"]:
            self.rk.parameters["timestepping"]["dtmax"] = self.solver_params["timestepping"]["dtmax"]
        if 'dtmin' in self.solver_params["timestepping"]:
            self.rk.parameters["timestepping"]["dtmin"] = self.solver_params["timestepping"]["dtmin"]
        self.rk.parameters["verbose"] = self.solver_params["verbose"]
        self.rk.parameters["drawplot"] = self.solver_params["drawplot"]
        self.rk.parameters["output"]["statistics"] = self.solver_params["output"]["statistics"]
        self.rk.parameters["output"]["print_statistics"] = self.solver_params["output"]["print_statistics"]
    
    # SOLVE
    def solve(self):        
        self.rk_instance()  # RK solver instances
        self.rk.solve()    

    # DOF map
    def dof_map(self):
        z = self.mesh.coordinates()
        dofmap = fe.dof_to_vertex_map(self.V)
        return (z,dofmap)

    # Find DOF number matching z == z_search
    def find_dof(self,z_search):
        z,dofmap = self.dof_map()
        k_z = np.argmin(abs(z-z_search))     # Find closest number in z
        k_dof = np.argmin(abs(dofmap-k_z))   # Find cell number matching z == z_search
        return(k_dof)

    # Evaluate solution
    def dofs_values(self): 
        z = self.mesh.coordinates()
        c = self.c.compute_vertex_values()
        return (z,c)