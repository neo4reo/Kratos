/*
==============================================================================
KratosShallowWaterApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//   
//   Project Name:        Kratos       
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                June 28th 2017
//   Revision:            1.4
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/primitive_var_element.hpp"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//----------------------------------------------------------------------

	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			
		unsigned int element_size = TNumNodes*3;
		if(rResult.size() != element_size)
			rResult.resize(element_size,false);                         // False says not to preserve existing storage!!
		
		GeometryType& rGeom = GetGeometry();
		int counter=0;
		for (unsigned int i = 0; i < TNumNodes; i++)
		{
			rResult[counter++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[counter++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
			rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
		}
		
		KRATOS_CATCH("")
	}

//----------------------------------------------------------------------

	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
		const unsigned int element_size = TNumNodes*3;
		if(rElementalDofList.size() != element_size)
			rElementalDofList.resize(element_size);
		
		GeometryType& rGeom = GetGeometry();
		int counter=0;
		for (unsigned int i = 0; i < TNumNodes; i++)
		{
			rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
			rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
			rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
		}
		
		KRATOS_CATCH("")
	}

//----------------------------------------------------------------------

	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
        // Resize of the Left and Right Hand side
		unsigned int element_size = TNumNodes*3;
		if(rLeftHandSideMatrix.size1() != element_size)
			rLeftHandSideMatrix.resize(element_size,element_size,false); // False says not to preserve existing storage!!
		
		if(rRightHandSideVector.size() != element_size)
			rRightHandSideVector.resize(element_size,false);             // False says not to preserve existing storage!!
		
		// Getting gravity
		array_1d<double,3> v_gravity = rCurrentProcessInfo[GRAVITY];
		double gravity = 9.8; //-v_gravity[2];
		
		// Getting the time step (not fixed to allow variable time step)
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		const double dt_inv = 1.0 / delta_t;
		
		// Compute the geometry
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, 2> DN_DX;
        array_1d<double,TNumNodes> N;
        double Area;
        this-> CalculateGeometry(DN_DX,Area);
        double elem_length = this->ComputeElemSize(DN_DX);
		
        // Getting the values of shape functions on Integration Points
        boost::numeric::ublas::bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;  // In this case, number of Gauss points and number of nodes coincides
        const GeometryType& rGeom = this->GetGeometry();
        Ncontainer = rGeom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );
		
		// Get nodal values for current step and projected variables
		array_1d<double, TNumNodes*3> v_depth;
		array_1d<double, TNumNodes*3> v_unknown;
		array_1d<double, TNumNodes*3> v_proj_unknown;
		double height;
		GetNodalValues(v_depth,v_unknown,v_proj_unknown,height);
		
        // Some auxilary definitions
		boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> N_vel        = ZeroMatrix(2,TNumNodes*3);  // Shape functions matrix (for velocity unknown)
		boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> N_height     = ZeroMatrix(1,TNumNodes*3);  // Shape functions vector (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,1,TNumNodes*3> DN_DX_vel    = ZeroMatrix(1,TNumNodes*3);  // Shape functions gradients vector (for velocity unknown)
		boost::numeric::ublas::bounded_matrix<double,2,TNumNodes*3> DN_DX_height = ZeroMatrix(2,TNumNodes*3);  // Shape functions gradients matrix (for height unknown)
		//
		boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> mass_matrix  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
		boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_q_div_u  = ZeroMatrix(TNumNodes*3,TNumNodes*3);
		boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_w_grad_h = ZeroMatrix(TNumNodes*3,TNumNodes*3);
		boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_u_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
		boost::numeric::ublas::bounded_matrix<double,TNumNodes*3,TNumNodes*3> aux_h_diffus = ZeroMatrix(TNumNodes*3,TNumNodes*3);
		
		// Loop on Gauss points. In this case, number of Gauss points and number of nodes coincides
        for(unsigned int igauss = 0; igauss < TNumNodes; igauss++)
        {
			noalias(N) = row(Ncontainer, igauss);
            
            // Build shape and derivatives functions at Gauss points
            for(unsigned int nnode = 0; nnode < TNumNodes; nnode++)
            {
				// Height gradient
				DN_DX_height(0, 2+nnode*3) = DN_DX(nnode,0);
				DN_DX_height(1, 2+nnode*3) = DN_DX(nnode,1);
				// Velocity divergence
				DN_DX_vel(0,   nnode*3) = DN_DX(nnode,0);
				DN_DX_vel(0, 1+nnode*3) = DN_DX(nnode,1);
				// Height shape funtions
				N_height(0, 2+nnode*3) = N[nnode];
				// Velocity shape functions
				N_vel(0,   nnode*3) = N[nnode];
				N_vel(1, 1+nnode*3) = N[nnode];
			}
			
			noalias(mass_matrix)  += prod(trans(N_vel),N_vel);
			noalias(mass_matrix)  += prod(trans(N_height),N_height);
			
			noalias(aux_q_div_u)  += prod(trans(N_height),DN_DX_vel);
			noalias(aux_w_grad_h) += prod(trans(N_vel),DN_DX_height);
			
			noalias(aux_u_diffus) += prod(trans(DN_DX_vel),DN_DX_vel);
			noalias(aux_h_diffus) += prod(trans(DN_DX_height),DN_DX_height);
		}
		
		// Copmute stabilization parameters
		bool stabilization = true;
		double Ctau = 0.01;
		double tau_u, tau_h = 0;
		if (stabilization)
		{
			if (height > 1e-6)
				tau_u = Ctau/elem_length*pow(gravity/height,0.5);
			tau_h = Ctau/elem_length*pow(height/gravity,0.5);
		}
		// Compute discontinuity capturing parameters
		bool discontinuity_capturing = true;
		double gradient_threshold = 1e-6;
		//~ double residual;
		double height_grad_norm = norm_2(prod(DN_DX_height,v_unknown));
		double k_dc = 0;
		if(discontinuity_capturing && height_grad_norm > gradient_threshold)
		{
			k_dc = 0.5*0.4*elem_length*height_grad_norm;  // Residual formulation
		}
		
		// Build LHS
		// Cross terms
		noalias(rLeftHandSideMatrix)  = height * aux_q_div_u;           // Add <q*h*div(u)> to Mass Eq.
		noalias(rLeftHandSideMatrix) += gravity * aux_w_grad_h;         // Add <w*g*grad(h)> to Momentum Eq.
		
		// Inertia terms
		noalias(rLeftHandSideMatrix) += dt_inv * mass_matrix;           // Add <N,N> to both Eq's
		
        // Stabilization terms
		noalias(rLeftHandSideMatrix) += (k_dc + tau_h) * aux_h_diffus;  // Add art. diff. to Mass Eq.
		noalias(rLeftHandSideMatrix) +=         tau_u  * aux_u_diffus;  // Add art. diff. to Momentum Eq.
		
		// Build RHS
		// Source term (bathymetry contribution)
		noalias(rRightHandSideVector)  = -gravity * prod(aux_w_grad_h, v_depth);
		
		// Inertia terms
		noalias(rRightHandSideVector) += dt_inv * prod(mass_matrix, v_proj_unknown);
		
		// Substracting the Dirichlet term (since we use a residualbased approach)
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, v_unknown);
		
		rRightHandSideVector *= Area / static_cast<double>(TNumNodes);
		rLeftHandSideMatrix *= Area  / static_cast<double>(TNumNodes);
		
		KRATOS_CATCH("")
	}

//----------------------------------------------------------------------

	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}
	
//----------------------------------------------------------------------

	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU)
		{
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
				rValues[PointNumber] = double(this->GetValue(rVariable));
		}
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::CalculateGeometry(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
    {
        const GeometryType& rGeom = this->GetGeometry();

        // We select GI_GAUSS_1 due to we are computing at the barycenter.
        const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);
        const unsigned int NumGPoints = integration_points.size();
        rArea = rGeom.Area();
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer( NumGPoints );
        rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, GeometryData::GI_GAUSS_1);

        noalias( rDN_DX ) = DN_DXContainer[0];

    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    double PrimitiveVarElement<TNumNodes>::ComputeElemSize(boost::numeric::ublas::bounded_matrix<double, TNumNodes, 2>& rDN_DX)
    {
        double l = 0.0;

        for(unsigned int i = 0; i < TNumNodes; i++)
        {
            double l_inv = 0.0;
            for(unsigned int k = 0; k < 2; k++)
            {
                l_inv += rDN_DX(i,k) * rDN_DX(i,k);
            }
            l += 1.0 / l_inv;
        }
        l = sqrt(l) / static_cast<double>(TNumNodes);
        return l;
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void PrimitiveVarElement<TNumNodes>::GetNodalValues(array_1d<double, TNumNodes*3>& rdepth,
                                                        array_1d<double, TNumNodes*3>& runkn, 
                                                        array_1d<double, TNumNodes*3>& rproj, 
                                                        double& rheight)
    {
		double lumping_factor = 1 / double(TNumNodes);
		
		rheight = 0;
		unsigned int counter = 0;
		for (unsigned int i = 0; i < TNumNodes; i++)
		{
			rdepth[counter] = 0;
			runkn[counter]  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X);
			rproj[counter]  = GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
			counter++;

			rdepth[counter] = 0;
			runkn[counter]  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y);
			rproj[counter]  = GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
			counter++;

			rdepth[counter] = GetGeometry()[i].FastGetSolutionStepValue(BATHYMETRY);
			runkn[counter]  = GetGeometry()[i].FastGetSolutionStepValue(HEIGHT);
			rproj[counter]  = GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_HEIGHT);
			counter++;
			
			rheight += GetGeometry()[i].FastGetSolutionStepValue(HEIGHT);
		}
		
		rheight *= lumping_factor;
	}

//----------------------------------------------------------------------

template class PrimitiveVarElement<3>;
template class PrimitiveVarElement<4>;

} // namespace Kratos