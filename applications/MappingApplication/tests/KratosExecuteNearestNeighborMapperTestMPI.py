from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *

import os
import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteNearestNeighborMapperTestMPI(KratosUnittest.TestCase):

    def __init__(self):
        # Mdpa Input files
        input_file_structure = "Mapper_Test_1/FSI_Example4Mapper_1_Structural"
        input_file_fluid     = "Mapper_Test_1/FSI_Example4Mapper_1_Fluid"

        # Create and Partition Model Parts
        variable_list = [PRESSURE, VELOCITY, PARTITION_INDEX]
        self.model_part_origin  = self.partition_and_read_model_part("ModelPartNameOrigin", input_file_fluid, 3, variable_list)
        self.model_part_destination = self.partition_and_read_model_part("ModelPartNameDestination", input_file_structure, 3, variable_list)

        self.interface_sub_model_part_origin = self.model_part_origin.GetSubModelPart("FluidNoSlipInterface3D_interface_orig_fluid")

        SetMPICommunicatorProcess(self.interface_sub_model_part_origin).Execute()
        pfc_origin = ParallelFillCommunicator(self.interface_sub_model_part_origin)
        pfc_origin.Execute()

        self.interface_sub_model_part_destination = self.model_part_destination.GetSubModelPart("StructureInterface3D_interface_dest_struct")

        SetMPICommunicatorProcess(self.interface_sub_model_part_destination).Execute()
        pfc_destination = ParallelFillCommunicator(self.interface_sub_model_part_destination)
        pfc_destination.Execute()

        # Initialize Mapper
        self.nearestNeighborMapper = NearestNeighborMapper(self.interface_sub_model_part_origin, self.interface_sub_model_part_destination)


    def TestMapConstantScalarValues(self):
        map_value = 5.123
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_origin, variable_origin, map_value)

        # Overwriting Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination)
        self.CheckValues(self.interface_sub_model_part_destination, variable_destination, map_value)

        # Adding Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination, True)
        self.CheckValues(self.interface_sub_model_part_destination, variable_destination, map_value*2)

    def TestInverseMapConstantScalarValues(self):
        map_value = -8.6647
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_destination, variable_destination, map_value)

        # Overwriting Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)
        self.CheckValues(self.interface_sub_model_part_origin, variable_origin, map_value)

        # Adding Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination, True)
        self.CheckValues(self.interface_sub_model_part_origin, variable_origin, map_value*2)

    def TestMapConstantVectorValues(self):
        map_value = [15.99, -2.88, 3.123]
        variable_origin = VELOCITY
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_X, map_value[0])
        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_Y, map_value[1])
        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_Z, map_value[2])

        # Overwriting Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination)

        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_X, map_value[0])
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Y, map_value[1])
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Z, map_value[2])

        # Adding Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination, True)

        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_X, map_value[0]*2)
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Y, map_value[1]*2)
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Z, map_value[2]*2)

    def TestInverseMapConstantVectorValues(self):
        map_value = [1.4785, -0.88, -33.123]
        variable_origin = VELOCITY
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_X, map_value[0])
        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_Y, map_value[1])
        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_Z, map_value[2])

        # Overwriting Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)

        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_X, map_value[0])
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Y, map_value[1])
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Z, map_value[2])

        # Adding Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination, True)

        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_X, map_value[0]*2)
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Y, map_value[1]*2)
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Z, map_value[2]*2)

    def TestMapNonConstantScalarValues(self):
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesId(self.interface_sub_model_part_origin, variable_origin)

        self.nearestNeighborMapper.Map(variable_origin, variable_destination)

        nodal_values = self.GetNodalValuesIdMap()
        self.CheckValuesId(self.interface_sub_model_part_destination, variable_destination, nodal_values)

    def TestInverseMapNonConstantScalarValues(self):
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesId(self.interface_sub_model_part_destination, variable_destination)

        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)
        nodal_values = self.GetNodalValuesIdInverseMap()
        self.CheckValuesId(self.interface_sub_model_part_origin, variable_origin, nodal_values)


    def partition_and_read_model_part(self, model_part_name, model_part_input_file, size_domain, variable_list, number_of_partitions = mpi.size):
        model_part = ModelPart(model_part_name)
        for variable in variable_list:
            model_part.AddNodalSolutionStepVariable(variable)

        # number of partitions is by default equal to mpi.size
        if mpi.size > 1:
            if mpi.rank == 0:
                model_part_io = ReorderConsecutiveModelPartIO(model_part_input_file)

                partitioner = MetisDivideHeterogeneousInputProcess(
                    model_part_io,
                    number_of_partitions,
                    size_domain,
                    0, # verbosity, set to 1 for more detailed output
                    True)

                partitioner.Execute()

            mpi.world.barrier()
            model_part_input_file = model_part_input_file + "_" + str(mpi.rank)

        MPICommSetup = SetMPICommunicatorProcess(model_part)
        MPICommSetup.Execute()

        model_part_io = ModelPartIO(model_part_input_file)
        model_part_io.ReadModelPart(model_part)

        model_part.ProcessInfo.SetValue(DOMAIN_SIZE, size_domain)
        model_part.SetBufferSize(1)

        return model_part

    def SetValuesOnNodes(self, model_part, variable, value):
        for node in model_part.Nodes:
            if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                node.SetSolutionStepValue(variable, value)

    def SetValuesOnNodesId(self, model_part, variable):
        for node in model_part.Nodes:
            if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                node.SetSolutionStepValue(variable, node.Id)

    def CheckValues(self, model_part, variable, map_value):
        for node in model_part.Nodes:
            if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                value = node.GetSolutionStepValue(variable)
                self.assertAlmostEqual(map_value,value)

    def CheckValuesId(self, model_part, variable, nodal_values):
        i = 0
        for node in model_part.Nodes:
            if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                value = node.GetSolutionStepValue(variable)
                self.assertAlmostEqual(nodal_values[i],value)
                i += 1
        if (len(nodal_values) != i):
            raise Exception("Array Sizes are not matching")

    def GetNodalValuesIdMap(self):
        if (mpi.size == 2):
            if (mpi.rank == 0):
                nodal_ids = [56,53,6,14,1,17,25,12,57,37,61,63]
            elif (mpi.rank == 1):
                nodal_ids = [11,39,34,40,7,16,29,20,31,35,55,64,66]

        elif (mpi.size == 4):
            if (mpi.rank == 0):
                nodal_ids = [56,53,6,12,63]
            elif (mpi.rank == 1):
                nodal_ids = [14,1,25,31,57,37,61]
            elif (mpi.rank == 2):
                nodal_ids = [34,20,35,55,64]
            elif (mpi.rank == 3):
                nodal_ids = [11,39,40,7,16,5,29,66]

        elif (mpi.size == 8):
            if (mpi.rank == 0):
                nodal_ids = [6,5,12]
            elif (mpi.rank == 1):
                nodal_ids = [63]
            elif (mpi.rank == 2):
                nodal_ids = [56,53,14,1]
            elif (mpi.rank == 3):
                nodal_ids = [29,25,57,37,61]
            elif (mpi.rank == 4):
                nodal_ids = [35,55,64]
            elif (mpi.rank == 5):
                nodal_ids = [39,40]
            elif (mpi.rank == 6):
                nodal_ids = [11,34,7,16,20,31,66]
            elif (mpi.rank == 7):
                nodal_ids = []
        else:
            raise Exception("Id-Mapping \"Map\" not implemented for " + str(mpi.size) + " processors")

        return nodal_ids

    def GetNodalValuesIdInverseMap(self):
        if (mpi.size == 3):
            if (mpi.rank == 0):
                nodal_ids = [16,1,16,18,31,3,3,2,5,2,78]
            elif (mpi.rank == 1):
                nodal_ids = [48,31,10,33,33,49,49,50,60,53,6,54,74,72]
            elif (mpi.rank == 2):
                nodal_ids = [12,12,12,18,8,13,43,10,31,7,6,77]

        elif (mpi.size == 5):
            if (mpi.rank == 0):
                nodal_ids = [16,18,3,2,5,2,78]
            elif (mpi.rank == 1):
                nodal_ids = [12,12,12,18,8,13,16,43,7,6,77]
            elif (mpi.rank == 2):
                nodal_ids = [1,31,48,3,53]
            elif (mpi.rank == 3):
                nodal_ids = [31,33,49,49,50,54,74,72]
            elif (mpi.rank == 4):
                nodal_ids = [10,10,33,31,60,6]

        elif (mpi.size == 9):
            if (mpi.rank == 0):
                nodal_ids = [12,13,1]
            elif (mpi.rank == 1):
                nodal_ids = [16,16,5,78]
            elif (mpi.rank == 2):
                nodal_ids = [8,43,7,77]
            elif (mpi.rank == 3):
                nodal_ids = [12,12,18,18]
            elif (mpi.rank == 4):
                nodal_ids = [10,10,6,6]
            elif (mpi.rank == 5):
                nodal_ids = [33,49,60,54,74]
            elif (mpi.rank == 6):
                nodal_ids = [3,3,2,2]
            elif (mpi.rank == 7):
                nodal_ids = [49,50,53,72]
            elif (mpi.rank == 8):
                nodal_ids = [31,48,31,33,31]
        else:
            raise Exception("Id-Mapping \"InverseMap\" not implemented for " + str(mpi.size) + " processors")

        return nodal_ids