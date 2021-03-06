//
//   Project Name:   
//   Last modified by:    $Author:   
//   Date:                $Date:  
//   Revision:            $Revision:  
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.hpp"

namespace Kratos
{

class ThermalLinearElastic2DPlaneStressNodal : public ThermalLinearElastic2DPlaneStrainNodal
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic2DPlaneStressNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic2DPlaneStressNodal();

    // Copy Constructor
    ThermalLinearElastic2DPlaneStressNodal (const ThermalLinearElastic2DPlaneStressNodal& rOther);

    // Destructor
    virtual ~ThermalLinearElastic2DPlaneStressNodal();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    ConstitutiveLaw::Pointer Clone() const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures);

    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, const double &rYoungModulus, const double &rPoissonCoefficient );

    void CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables & rElasticVariables, double & rTemperature);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ThermalLinearElastic2DPlaneStrainNodal)
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ThermalLinearElastic2DPlaneStrainNodal)
    }

}; // Class ThermalLinearElastic2DPlaneStressNodal
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED  defined 