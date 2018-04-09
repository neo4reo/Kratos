// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//         -        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//         -        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//                  in the documentation and/or other materials provided with the distribution.
//         -        All advertising materials mentioning features or use of this software must display the following acknowledgement:
//                         This product includes Kratos Multi-Physics technology.
//         -        Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "fem_to_dem_application.h"
#include "fem_to_dem_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{

namespace Python
{

  using namespace boost::python;



  BOOST_PYTHON_MODULE(KratosFemToDemApplication)
  {

	  class_<KratosFemToDemApplication,
			  KratosFemToDemApplication::Pointer,
			  bases<KratosApplication>, boost::noncopyable >("KratosFemToDemApplication")
			;

	AddCustomStrategiesToPython();
	AddCustomUtilitiesToPython();
	AddCustomConstitutiveLawsToPython();
	AddCustomProcessesToPython();

	//registering variables in python
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DAMAGE_EDGE1);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DAMAGE_EDGE2);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DAMAGE_EDGE3);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DAMAGE_ELEMENT);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_VECTOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YIELD_STRESS_C);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YIELD_STRESS_T);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRAC_ENERGY_T)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRAC_ENERGY_C)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(ITER);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YIELD_SURFACE)
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRAIN_VECTOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_VECTOR_INTEGRATED);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(TANGENT_CONSTITUTIVE_TENSOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(SMOOTHING);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_DAMAGED);
	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(CHARACTERISTIC_LENGTH);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(MESH_REFINED);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_DYNAMIC);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_THRESHOLD);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(INITIAL_THRESHOLD);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(INTEGRATION_COEFFICIENT);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(MAPPING_PROCEDURE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_DEM);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_RADIUS);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DEM_GENERATED);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(INACTIVE_NODE);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NUMBER_OF_ACTIVE_ELEMENTS);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_FORCE_APPLIED);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_FORCE_X);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_FORCE_Y);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_FORCE_Z);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_STRESS_VECTOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(EQUIVALENT_NODAL_STRESS);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(EQUIVALENT_NODAL_STRESS_GRADIENT);
	KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(AUXILIAR_GRADIENT);

	// 3D case
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRAIN_TENSOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_TENSOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STRESS_TENSOR_INTEGRATED);
	
	// Composite calculations
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCRETE_STRESS_TENSOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STEEL_STRESS_TENSOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCRETE_STRESS_VECTOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STEEL_STRESS_VECTOR);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YOUNG_MODULUS_STEEL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(DENSITY_STEEL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(POISSON_RATIO_STEEL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(STEEL_VOLUMETRIC_PART);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(CONCRETE_STRESS_TENSOR_INTEGRATED);
	
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YIELD_STRESS_C_STEEL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(YIELD_STRESS_T_STEEL);
	KRATOS_REGISTER_IN_PYTHON_VARIABLE(FRACTURE_ENERGY_STEEL);
  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined