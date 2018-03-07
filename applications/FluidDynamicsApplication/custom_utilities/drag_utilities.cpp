//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Application includes
#include "drag_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    array_1d<double, 3> DragUtilities::CalculateBodyFittedDrag(ModelPart& rModelPart) {

        // Get the root model part
        ModelPart& root_model_part = rModelPart.GetRootModelPart();

        // Initialize the reaction variable
        const array_1d<double, 3> zero_vect(3,0.0);
        #pragma omp parallel for firstprivate(zero_vect)
        for (int i = 0; i < static_cast<int>(root_model_part.NumberOfNodes()); ++i){
            auto it_node = root_model_part.NodesBegin() + i;
            noalias(it_node->FastGetSolutionStepValue(REACTION)) = zero_vect;
        }

        Vector RHS_Contribution;
        Matrix LHS_Contribution;
        ProcessInfo& r_current_process_info = root_model_part.GetProcessInfo();
        const unsigned int domain_size = r_current_process_info[DOMAIN_SIZE];

        // Fractional step case: set the step index to compute the momentum equation
        int current_step = 0;
        if (r_current_process_info.Has(FRACTIONAL_STEP)) {
            current_step = r_current_process_info[FRACTIONAL_STEP];
            r_current_process_info[FRACTIONAL_STEP] = 1;
        }

        #pragma omp parallel for private(RHS_Contribution, LHS_Contribution) firstprivate(domain_size)
        for (int i_elem = 0; i_elem < static_cast<int>(root_model_part.NumberOfElements()); ++i_elem){
            auto it_elem = root_model_part.ElementsBegin() + i_elem;

            // Build local system
            it_elem->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, r_current_process_info);

            // Get geometry
            Element::GeometryType& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();
            const unsigned int block_size = RHS_Contribution.size()/n_nodes;

            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node){
                r_geom[i_node].SetLock();
                array_1d<double,3>& r_reaction = r_geom[i_node].FastGetSolutionStepValue(REACTION);
                for (unsigned int d = 0; d < domain_size; ++d)
                    r_reaction[d] -= RHS_Contribution[block_size*i_node + d];

                r_geom[i_node].UnSetLock();
            }
        }

        // Assemble reaction data
        root_model_part.GetCommunicator().AssembleCurrentData(REACTION);

        // Sum the reactions in the model part of interest
        VariableUtils variable_utils;
        array_1d<double, 3> drag_force = variable_utils.SumHistoricalNodeVectorVariable(REACTION, rModelPart, 0);
        drag_force *= -1.0;

        // Fractional step case: restore the step index
        if (r_current_process_info.Has(FRACTIONAL_STEP)) {
            r_current_process_info[FRACTIONAL_STEP] = current_step;
        }

        return drag_force;

    }

    array_1d<double, 3> DragUtilities::CalculateEmbeddedDrag(ModelPart& rModelPart) {
        
        // Initialize total drag force
        array_1d<double, 3> drag_force = ZeroVector(3);
        double& drag_x = drag_force[0];
        double& drag_y = drag_force[1];
        double& drag_z = drag_force[2];

        // Iterate the model part elements to compute the drag
        array_1d<double, 3> elem_drag;

        // Auxiliary var to make the reduction
        double drag_x_red = 0.0;
        double drag_y_red = 0.0;
        double drag_z_red = 0.0;

        #pragma omp parallel for reduction(+:drag_x_red) reduction(+:drag_y_red) reduction(+:drag_z_red) private(elem_drag) schedule(dynamic)
        for(int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Calculate(DRAG_FORCE, elem_drag, rModelPart.GetProcessInfo());
            drag_x_red += elem_drag[0];
            drag_y_red += elem_drag[1];
            drag_z_red += elem_drag[2];
        }
        
        drag_x += drag_x_red;
        drag_y += drag_y_red;
        drag_z += drag_z_red;

        // Perform MPI synchronization
        rModelPart.GetCommunicator().SumAll(drag_force);

        return drag_force;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const DragUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}