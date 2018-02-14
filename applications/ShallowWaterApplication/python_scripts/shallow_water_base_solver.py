from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow       # Importing our application
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.PFEM2Application as KratosPfem2
import KratosMultiphysics.MeshingApplication as KratosMeshing

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(model_part, custom_settings):
    return ShallowWaterBaseSolver(model_part, custom_settings)

class ShallowWaterBaseSolver(object):
    def __init__(self,model_part, custom_settings):  # Constructor of the class 
        self.model_part  = model_part
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"                  : "shallow_water_base_solver",
            "model_import_settings"        : {
                "input_type"          : "mdpa",
                "input_filename"      : "unknown_name"
            },
            "echo_level"                   : 0,
            "convergence_echo_level"       : 1,
            "solver_echo_level"            : 0,
            "buffer_size"                  : 2,
            "dynamic_tau"                  : 0.005,
            "relative_tolerance"           : 1e-6,
            "absolute_tolerance"           : 1e-9,
            "maximum_iterations"           : 20,
            "compute_reactions"            : false,
            "reform_dofs_at_each_step"     : false,
            "calculate_norm_dx"            : true,
            "move_mesh_flag"               : false,
            "element_name"                 : "PrimitiveVarElement2D3N",
            "condition_name"               : "NothingCondition2D2N",
            "volume_model_part_name"       : "volume_model_part",
            "skin_parts"                   : [""],
            "no_skin_parts"                : [""],
            "linear_solver_settings"       : {
                    "solver_type"     : "SkylineLUFactorizationSolver"
            },
            "time_stepping"                : {
                "automatic_time_step" : false,
                "time_step"           : 0.01
            },
            "pfem2_settings"               : {
                "convection_scalar_variable"    : "HEIGHT",
                "convection_vector_variable"    : "VELOCITY",
                "maximum_number_of_particles"   : 16
            }
        }""")
        default_settings["pfem2_settings"]["maximum_number_of_particles"].SetInt(8*self.domain_size)

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        self.Mesher = KratosMeshing.TriGenPFEMModeler()

        self.element_name = self.settings["element_name"].GetString()
        self.condition_name = self.settings["condition_name"].GetString()

        # Initialize shallow water variables utility
        self.ShallowVariableUtils = KratosShallow.ShallowWaterVariablesUtility(self.model_part)

    def AddVariables(self):
        # Basic variables
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.HEIGHT);
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY);
        # Physic problem parameters
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.FREE_SURFACE_ELEVATION);
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.GRAVITY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.BATHYMETRY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.RAIN);
        # Auxiliar variables
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)   # To compute the normals
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)         # Slip condition needs the normal
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)  # Some utilities needs the variable, but not the application
        # Variables required by TriGenPFEM
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE) # deprecated_variables.h
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)     # deprecated_variables.h
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)        # deprecated_variables.h
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)         # variables.h
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)      # variables.h

    def AddDofs(self):
        raise Exception("Calling the base class instead of the derived one")

    def GetMinimumBufferSize(self):
        return 1

    def GetComputingModelPart(self):
        return self.model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString())

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions, check the input reading
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # Initializing the remesher
        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part, number_of_avg_elems, number_of_avg_nodes)

        self.mark_outer_nodes_process = KratosPfem2.MarkOuterNodesProcess(self.computing_model_part)

        self.fluid_neigh_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part,9,18)
        self.elem_neigh_finder  = KratosMultiphysics.FindElementalNeighboursProcess(self.computing_model_part, 2, 10)
        self.cond_neigh_finder  = KratosMultiphysics.FindConditionsNeighboursProcess(self.computing_model_part,2, 10)

        (self.neighbour_search).Execute()
        (self.fluid_neigh_finder).Execute();

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = KratosFluid.EstimateDtUtility2D(self.computing_model_part,
                                                                            self.settings["time_stepping"])

        # Creating the solution strategy for the mesh stage
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                     self.settings["absolute_tolerance"].GetDouble())
        (self.conv_criteria).SetEchoLevel(self.settings["convergence_echo_level"].GetInt())

        #~ self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.domain_size,   # DomainSize
                                                                                             self.domain_size+1) # BlockSize

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())

        (self.solver).SetEchoLevel(self.settings["solver_echo_level"].GetInt())
        (self.solver).Check()

        (self.solver).Initialize()

        print ("Mesh stage solver initialization finished")

    def Solve(self):
        # Solve equations on mesh
        (self.solver).Solve()

    def Clear(self):
        (self.solver).Clear()

    def ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        # Move particles utility needs to access delta_time
        return delta_time

    def Remesh(self):

        h_factor = 0.1;
        alpha_shape = 1.2;

        for node in (self.model_part).Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H,0,0.1)

        node_erase_process = KratosMultiphysics.NodeEraseProcess(self.model_part);

        rem_nodes = False
        add_nodes = True
        (self.Mesher).ReGenerateMesh(self.element_name, self.condition_name, self.computing_model_part, node_erase_process, rem_nodes, add_nodes, alpha_shape, h_factor)

        (self.fluid_neigh_finder).Execute()
        (self.elem_neigh_finder).Execute()
        (self.cond_neigh_finder).Execute()


    #### Specific internal functions ####

    def _ModelPartReading(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.model_part)
        else:
            raise Exception("Other input options are not yet implemented")

    def _ExecuteAfterReading(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])

    def _SetBufferSize(self):
        ## Set the buffer size
        self.model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        if(self.GetMinimumBufferSize() > self.model_part.GetBufferSize() ):
            self.model_part.SetBufferSize(self.GetMinimumBufferSize())

    def _AddPrimitiveDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.VELOCITY_X);
            node.AddDof(KratosMultiphysics.VELOCITY_Y);
            node.AddDof(KratosShallow.HEIGHT);
        print ("Primitive variables for the SWE solver added correctly")

    def _AddConservedDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.MOMENTUM_X);
            node.AddDof(KratosMultiphysics.MOMENTUM_Y);
            node.AddDof(KratosShallow.HEIGHT);
        print ("Conserved variables for the SWE solver added correctly")

