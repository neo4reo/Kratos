// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_elements/SprismElement3D6N.hpp"

namespace Kratos
{
/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, EAS_IMPLICIT_EXPLICIT,              4 ); // True means implicit // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, TOTAL_UPDATED_LAGRANGIAN,           5 ); // True means total lagrangian // TODO: change this using templates!!!
KRATOS_CREATE_LOCAL_FLAG( SprismElement3D6N, QUADRATIC_ELEMENT,                  6 ); // True means quadratic in-plane behaviour // TODO: Idem

// ------------------------------------------------------------------------- //
// ------------------------------ PUBLIC ----------------------------------- //
// ------------------------------------------------------------------------- //

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N( )
        : Element( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

/************************************* CONSTRUCTOR *********************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    mFinalizedStep = true; // the creation is out of the time step, it must be true

    if( GetProperties().Has(NINT_TRANS) ) {
        if (GetProperties()[NINT_TRANS] == 2) {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_1;
        } else if (GetProperties()[NINT_TRANS] == 3) {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_2;
        } else if (GetProperties()[NINT_TRANS] == 5) {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_3;
        } else if (GetProperties()[NINT_TRANS] == 7) {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_4;
        } else if (GetProperties()[NINT_TRANS] == 11) {
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_5;
        } else {
            std::cout << "The number of integration points is not defined.  NINT_TRANS: "<< GetProperties()[NINT_TRANS] << std::endl;
            std::cout << "Options are: 2, 3, 5, 7, 11  " << std::endl;
            std::cout << "Taking default number of integration points (NINT_TRANS = 2)  " << std::endl;
            mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_2;
        }
    } else {
        mThisIntegrationMethod = GeometryData::GI_EXTENDED_GAUSS_2;
    }

    //DO NOT ADD DOFS HERE!!!
}

/*********************************** COPY CONSTRUCTOR ******************************/
/***********************************************************************************/

SprismElement3D6N::SprismElement3D6N( SprismElement3D6N const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    ,mFinalizedStep(rOther.mFinalizedStep)
    ,mTotalDomainInitialSize(rOther.mTotalDomainInitialSize)
    ,mAuxMatCont(rOther.mAuxMatCont)
    ,mAuxCont(rOther.mAuxCont)
{
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

SprismElement3D6N::~SprismElement3D6N()
{
}

/********************************** ASSIGMENT OPERATOR *****************************/
/***********************************************************************************/

SprismElement3D6N&  SprismElement3D6N::operator=(SprismElement3D6N const& rOther)
{
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size() );

    mAuxMatCont.clear();
    mAuxMatCont.resize( rOther.mAuxMatCont.size());

    for(IndexType i = 0; i < mConstitutiveLawVector.size(); i++) {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
        mAuxMatCont[i]=rOther.mAuxMatCont[i];
    }

    mTotalDomainInitialSize = rOther.mTotalDomainInitialSize;
    mAuxCont = rOther.mAuxCont;

    return *this;
}

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Element::Pointer SprismElement3D6N::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<SprismElement3D6N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/*********************************** CLONE ******************************************/
/************************************************************************************/

Element::Pointer SprismElement3D6N::Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const
{
    SprismElement3D6N new_element( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    new_element.mThisIntegrationMethod = mThisIntegrationMethod;

    const SizeType integration_point_number = mConstitutiveLawVector.size();

    if ( new_element.mConstitutiveLawVector.size() != integration_point_number)
        new_element.mConstitutiveLawVector.resize(integration_point_number);

    KRATOS_ERROR_IF( new_element.mConstitutiveLawVector.size() != new_element.GetGeometry().IntegrationPointsNumber() ) << "Constitutive law not has the correct size " << new_element.mConstitutiveLawVector.size() << std::endl;
    
    for(IndexType i = 0; i < integration_point_number; i++)
        new_element.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();

    //-----------//

    if ( new_element.mAuxMatCont.size() != mAuxMatCont.size() )
        new_element.mAuxMatCont.resize(mAuxMatCont.size());

    for(IndexType i = 0; i < mAuxMatCont.size(); i++)
        new_element.mAuxMatCont[i] = mAuxMatCont[i];

    new_element.mTotalDomainInitialSize = mTotalDomainInitialSize;
    new_element.mAuxCont = mAuxCont;

    return Kratos::make_shared<SprismElement3D6N>(new_element);
}

//******************************* GETTING METHODS *********************************//
/***********************************************************************************/
/***********************************************************************************/

SprismElement3D6N::IntegrationMethod SprismElement3D6N::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< NodeType >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);

    const IndexType NumberOfNodes = GetGeometry().size() + NumberOfActiveNeighbours(NeighbourNodes);
    const IndexType dim = NumberOfNodes * 3;

    if (rResult.size() != dim)
        rResult.resize(dim, false);

    // Nodes of the central element
    IndexType index = 0;
    for (IndexType i = 0; i < 6; i++) {
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        index += 3;
    }

    // Adding the ids of the neighbouring nodes
    for (IndexType i = 0; i < 6; i++) {
        if (HasNeighbour(i, NeighbourNodes[i])) {
            rResult[index]     = NeighbourNodes[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = NeighbourNodes[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = NeighbourNodes[i].GetDof(DISPLACEMENT_Z).EquationId();
            index += 3;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& CurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< NodeType >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    rElementalDofList.resize(0);

    // Nodes of the central element
    for (IndexType i = 0; i < GetGeometry().size(); i++) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    // Adding the dofs of the neighbouring nodes
    for (IndexType i = 0; i < 6; i++) {
        if (HasNeighbour(i, NeighbourNodes[i])) {
            rElementalDofList.push_back(NeighbourNodes[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(NeighbourNodes[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(NeighbourNodes[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("");
}

/******************************** DISPLACEMENT **************************************/
/************************************************************************************/

void SprismElement3D6N::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    const SizeType mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    IndexType index = 0;

    // Nodes of the central element
    for (IndexType i = 0; i < 6; i++) {
        const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        for (IndexType j = 0; j < 3; j++)
            rValues[index + j] = disp[j];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; i++) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            const array_1d<double, 3 > & disp = p_neighbour_nodes[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            for (IndexType j = 0; j < 3; j++)
                rValues[index + j] = disp[j];
            index += 3;
        }
    }
}

/********************************** VELOCITY ****************************************/
/************************************************************************************/

void SprismElement3D6N::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    const SizeType mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    IndexType index = 0;

    // Nodes of the central element
    for (IndexType i = 0; i < 6; i++) {
        const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        for (IndexType j = 0; j < 3; j++)
            rValues[index + j] = vel[j];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; i++) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            const array_1d<double, 3 > & vel = p_neighbour_nodes[i].FastGetSolutionStepValue(VELOCITY, Step);
            for (IndexType j = 0; j < 3; j++)
                rValues[index + j] = vel[j];
            index += 3;
        }
    }
}

/******************************** ACCELERATION **************************************/
/************************************************************************************/

void SprismElement3D6N::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    const SizeType mat_size = number_of_nodes * 3;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    IndexType index = 0;

    // Nodes of the central element
    for (IndexType i = 0; i < 6; i++) {
        const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        for (IndexType j = 0; j < 3; j++)
            rValues[index + j] = acc[j];
        index += 3;
    }

    // Neighbour nodes
    for (int i = 0; i < 6; i++) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            const array_1d<double, 3 > & acc = p_neighbour_nodes[i].FastGetSolutionStepValue(ACCELERATION, Step);
            for (IndexType j = 0; j < 3; j++)
                rValues[index + j] = acc[j];
            index += 3;
        }
    }
}

//****************************** COMPUTING METHODS ********************************//
/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR);

    MatrixType left_hand_side_matrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( left_hand_side_matrix, rRightHandSideVector, local_system.CalculationFlags );

    //Set general_variables to Local system components
    local_system.SetLeftHandSideMatrix(left_hand_side_matrix);
    local_system.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateRightHandSide(
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR);
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType left_hand_side_matrix = Matrix();

    /* Initialize sizes for the system components: */
    if( rRHSVariables.size() != rRightHandSideVectors.size() ) {
        rRightHandSideVectors.resize(rRHSVariables.size());
    }

    for( IndexType i = 0; i < rRightHandSideVectors.size(); i++ ) {
        this->InitializeSystemMatrices( left_hand_side_matrix, rRightHandSideVectors[i], local_system.CalculationFlags );
    }

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrix(left_hand_side_matrix);
    local_system.SetRightHandSideVectors(rRightHandSideVectors);

    local_system.SetRightHandSideVariables(rRHSVariables);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, true);
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR, false);

    VectorType right_hand_side_vector = Vector();

    /* Initialize sizes for the system components: */
    this->InitializeSystemMatrices( rLeftHandSideMatrix, right_hand_side_vector, local_system.CalculationFlags );

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    local_system.SetRightHandSideVector(right_hand_side_vector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags */
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, true);
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR, true);

    /* Initialize sizes for the system components: */
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, local_system.CalculationFlags );

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    local_system.SetRightHandSideVector(rRightHandSideVector);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLocalSystem(
    std::vector< MatrixType >& rLeftHandSideMatrices,
    const std::vector< Variable< MatrixType > >& rLHSVariables,
    std::vector< VectorType >& rRightHandSideVectors,
    const std::vector< Variable< VectorType > >& rRHSVariables,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create local system components */
    LocalSystemComponents local_system;

    /* Calculation flags*/
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    /* Initialize sizes for the system components: */
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() ) {
        rLeftHandSideMatrices.resize(rLHSVariables.size());
    }

    if( rRHSVariables.size() != rRightHandSideVectors.size() ) {
        rRightHandSideVectors.resize(rRHSVariables.size());
    }

    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX);
    for( IndexType i = 0; i < rLeftHandSideMatrices.size(); i++ ) {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], local_system.CalculationFlags );
    }

    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_RHS_VECTOR, true);
    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, false);

    for( IndexType i = 0; i < rRightHandSideVectors.size(); i++ ) {
        this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], local_system.CalculationFlags );
    }

    local_system.CalculationFlags.Set(SprismElement3D6N::COMPUTE_LHS_MATRIX, true);

    /* Set general_variables to Local system components */
    local_system.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    local_system.SetRightHandSideVectors(rRightHandSideVectors);

    local_system.SetLeftHandSideVariables(rLHSVariables);
    local_system.SetRightHandSideVariables(rRHSVariables);

    /* Calculate elemental system */
    CalculateElementalSystem( local_system, rCurrentProcessInfo );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    const double density = GetProperties()[DENSITY]; // TODO: Take into account the volume variation

    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);
    const SizeType mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size)
        rMassMatrix.resize(mat_size, mat_size, false);
    
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);
    
    double total_mass = GetGeometry().Volume() * density;

    const bool compute_lumped_mass_matrix =  rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX) ? rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] : false;

    // LUMPED MASS MATRIX, this one is easy because each node receives the same proportion of mass
    if (compute_lumped_mass_matrix == true) {
        Vector LumpFact;
        GetGeometry().LumpingFactors(LumpFact);
        for (IndexType i = 0; i < 3; i++) {
            double temp = LumpFact[i] * total_mass;
            for (IndexType j = 0; j < 6; j++) {
                IndexType index = i * 6 + j;
                rMassMatrix(index, index) = temp;
            }
        }
    } else { // CONSISTENT MASS MATRIX
    // Manually
        total_mass /= 72.0; // Dividing for the coefficient
        for (IndexType i = 0; i < 6; i++) { // Main nodes
            for (IndexType j = 0; j < 3; j++) { // DOF (X, Y, Z)
                const IndexType index = i * 3 + j;
                if (i == 0) {
                    // Superior band
                    rMassMatrix(index, index +  3) = 2.0 * total_mass;
                    rMassMatrix(index, index +  6) = 2.0 * total_mass;
                    rMassMatrix(index, index +  9) = 2.0 * total_mass;
                    rMassMatrix(index, index + 12) =       total_mass;
                    rMassMatrix(index, index + 15) =       total_mass;
                    // Symmetric part
                    rMassMatrix(index +  3, index) = 2.0 * total_mass;
                    rMassMatrix(index +  6, index) = 2.0 * total_mass;
                    rMassMatrix(index +  9, index) = 2.0 * total_mass;
                    rMassMatrix(index + 12, index) =       total_mass;
                    rMassMatrix(index + 15, index) =       total_mass;
                } else if (i == 1) {
                    // Superior band
                    rMassMatrix(index, index +  3) = 2.0 * total_mass;
                    rMassMatrix(index, index +  6) =       total_mass;
                    rMassMatrix(index, index +  9) = 2.0 * total_mass;
                    rMassMatrix(index, index + 12) =       total_mass;
                    // Symmetric part
                    rMassMatrix(index +  3, index) = 2.0 * total_mass;
                    rMassMatrix(index +  6, index) =       total_mass;
                    rMassMatrix(index +  9, index) = 2.0 * total_mass;
                    rMassMatrix(index + 12, index) =       total_mass;
                }  else if (i == 2) {
                    // Superior band
                    rMassMatrix(index, index + 3) =       total_mass;
                    rMassMatrix(index, index + 6) =       total_mass;
                    rMassMatrix(index, index + 9) = 2.0 * total_mass;
                    // Symmetric part
                    rMassMatrix(index + 3, index) =       total_mass;
                    rMassMatrix(index + 6, index) =       total_mass;
                    rMassMatrix(index + 9, index) = 2.0 * total_mass;
                } else if (i == 3) {
                    // Superior band
                    rMassMatrix(index, index + 3) = 2.0 * total_mass;
                    rMassMatrix(index, index + 6) = 2.0 * total_mass;
                    // Symmetric part
                    rMassMatrix(index + 3, index) = 2.0 * total_mass;
                    rMassMatrix(index + 6, index) = 2.0 * total_mass;
                } else if (i == 4) {
                    // Superior band
                    rMassMatrix(index, index + 3) = 2.0 * total_mass;
                    // Symmetric part
                    rMassMatrix(index + 3, index) = 2.0 * total_mass;
                }

                // Diagonal part
                rMassMatrix(index, index)         = 4.0 * total_mass;
            }
        }
    }
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const IndexType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    // Resizing as needed the LHS
    const IndexType mat_size = number_of_nodes * 3;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );
    
    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Calculate StiffnessMatrix:

    MatrixType stiffness_matrix  = Matrix();

    this->CalculateLeftHandSide( stiffness_matrix, rCurrentProcessInfo );

    // 2.-Calculate mass matrix:

    MatrixType mass_matrix  = Matrix();

    this->CalculateMassMatrix ( mass_matrix, rCurrentProcessInfo );

    // 3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) ) {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ) {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) ) {
        beta = GetProperties()[RAYLEIGH_BETA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ) {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    // 4.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * mass_matrix;
    noalias( rDampingMatrix ) += beta  * stiffness_matrix;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const MatrixType& rStiffnessMatrix,
        const MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        )
{
    KRATOS_TRY;

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    // 0.-Initialize the DampingMatrix:
    const SizeType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * 3;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );
    
    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) ) {
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ) {
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) ) {
        beta = GetProperties()[RAYLEIGH_BETA];
    } else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ) {
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    // 2.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * rMassMatrix;
    noalias( rDampingMatrix ) += beta  * rStiffnessMatrix;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number, false );

    if ( rVariable == VON_MISES_STRESS ) {
        /* Create and initialize element variables: */
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Create constitutive law parameters: */
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        /* Set constitutive law flags: */
        Flags &ConstitutiveLawOptions = Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < IntegrationPoints.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * IntegrationPoints[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy (Values);
            
            const Matrix& stress_tensor = MathUtils<double>::StressVectorToTensor(general_variables.StressVector); //reduced dimension stress tensor


            // In general coordinates:
            double sigma_equivalent =  (0.5)*((stress_tensor(0,0)-stress_tensor(1,1))*((stress_tensor(0,0)-stress_tensor(1,1)))+
                                            (stress_tensor(1,1)-stress_tensor(2,2))*((stress_tensor(1,1)-stress_tensor(2,2)))+
                                            (stress_tensor(2,2)-stress_tensor(0,0))*((stress_tensor(2,2)-stress_tensor(0,0)))+
                                            6*(stress_tensor(0,1)*stress_tensor(1,0)+stress_tensor(1,2)*stress_tensor(2,1)+stress_tensor(2,0)*stress_tensor(0,2)));

            if( sigma_equivalent < 0 )
                sigma_equivalent = 0;

            sigma_equivalent = std::sqrt(sigma_equivalent);

            rOutput[point_number] =  sigma_equivalent;
        }
    } else if ( rVariable == NORM_ISOCHORIC_STRESS ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < IntegrationPoints.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * IntegrationPoints[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components,point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
            {
                this->GetHistoricalVariables(general_variables,point_number);
            }

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy (Values);
            
            const Matrix& stress_tensor  = MathUtils<double>::StressVectorToTensor(general_variables.StressVector); //reduced dimension stress tensor

            double stress_norm =  ((stress_tensor(0,0)*stress_tensor(0,0))+(stress_tensor(1,1)*stress_tensor(1,1))+(stress_tensor(2,2)*stress_tensor(2,2))+
                                (stress_tensor(0,1)*stress_tensor(0,1))+(stress_tensor(0,2)*stress_tensor(0,2))+(stress_tensor(1,2)*stress_tensor(1,2))+
                                (stress_tensor(1,0)*stress_tensor(1,0))+(stress_tensor(2,0)*stress_tensor(2,0))+(stress_tensor(2,1)*stress_tensor(2,1)));

            stress_norm = std::sqrt(stress_norm);

            rOutput[point_number] = stress_norm;
        }
    } else if ( rVariable == STRAIN_ENERGY ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < IntegrationPoints.size(); point_number++ ) {
            const double ZetaGauss = 2.0 * IntegrationPoints[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components,point_number, alpha_eas, ZetaGauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            double StrainEnergy = 0.0;

            // Compute stresses and constitutive parameters
            if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN))
                mConstitutiveLawVector[point_number]->CalculateMaterialResponseKirchhoff(Values);
            else
                mConstitutiveLawVector[point_number]->CalculateMaterialResponsePK2(Values);

            mConstitutiveLawVector[point_number]->GetValue(STRAIN_ENERGY, StrainEnergy);

            rOutput[point_number] = general_variables.detJ * IntegrationPoints[point_number].Weight() * StrainEnergy;  // 1/2 * sigma * epsilon
        }
    } else {
        for ( IndexType ii = 0; ii < integration_point_number; ii++ )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }

    if ( rOutput.size() != 6 ) {
        std::vector<double> r_output_aux;
        r_output_aux = rOutput;

        rOutput.resize( 6, false );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (IndexType iii = 0; iii < 6; iii++) {
            rOutput[iii] = 0.0;

            for (IndexType i_gp = 0; i_gp < integration_point_number; i_gp++)
                rOutput[iii] += interpol(i_gp, iii) * r_output_aux[i_gp];
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &constitutive_laws_options=Values.GetOptions();

        constitutive_laws_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        constitutive_laws_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables, Values, point_number);

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR)
                general_variables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;
            else
                general_variables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(Values, general_variables.StressMeasure);

            if (rOutput[point_number].size() != general_variables.StressVector.size())
                rOutput[point_number].resize( general_variables.StressVector.size(), false);
            rOutput[point_number] = general_variables.StressVector;
        }
    } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR || rVariable == HENCKY_STRAIN_VECTOR) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, point_number, alpha_eas, zeta_gauss);

            // To take in account previous step writing
            if( mFinalizedStep )
                this->GetHistoricalVariables(general_variables,point_number);

            // Compute Green-Lagrange Strain
            if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ) // TODO: Replace CL
            {
                this->CalculateGreenLagrangeStrain( general_variables.C, general_variables.StrainVector );
            }
            else if( rVariable == ALMANSI_STRAIN_VECTOR ) // TODO: Replace CL
            {
                this->CalculateAlmansiStrain( general_variables.F, general_variables.StrainVector );
            }
            else if( rVariable == HENCKY_STRAIN_VECTOR )
            {
                this->CalculateHenckyStrain( general_variables.C, general_variables.StrainVector );
            }

            if (rOutput[point_number].size() != general_variables.StrainVector.size())
                rOutput[point_number].resize( general_variables.StrainVector.size(), false );

            rOutput[point_number] = general_variables.StrainVector;
        }
    } else {
        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable , rOutput[ii]);
    }

    if ( rOutput.size() != 6 ) {
        std::vector<Vector> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (IndexType iii = 0; iii < 6; iii++) {
            rOutput[iii] = ZeroVector(rOutput[0].size());

            for (IndexType i_gp = 0; i_gp < integration_point_number; i_gp++)
                rOutput[iii] += interpol(i_gp, iii) * rOutput_aux[i_gp];
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateOnIntegrationPoints(
    const Variable<Matrix >& rVariable,
    std::vector< Matrix >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const IndexType integration_point_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_point_number )
        rOutput.resize( integration_point_number );

    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
        std::vector<Vector> stress_vector;
        if( rVariable == CAUCHY_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

        // Loop integration points
        if ( rOutput.size() != stress_vector.size() )
            rOutput.resize( stress_vector.size() );

        for ( IndexType point_number = 0; point_number < rOutput.size(); point_number++ ) {
            if (rOutput[point_number].size2() != 3)
                rOutput[point_number].resize(3, 3, false);
            rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR || rVariable == HENCKY_STRAIN_TENSOR) {
        std::vector<Vector> StrainVector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else if ( rVariable == ALMANSI_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else if ( rVariable == HENCKY_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( HENCKY_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        // Loop integration points
        if ( rOutput.size() != StrainVector.size() )
            rOutput.resize( StrainVector.size() );

        for ( IndexType point_number = 0; point_number < rOutput.size(); point_number++ ) {
            if (rOutput[point_number].size2() != 3)
                rOutput[point_number].resize(3, 3, false);

            rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(StrainVector[point_number]);
        }
    } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &constitutive_laws_options=Values.GetOptions();
        constitutive_laws_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < IntegrationPoints.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * IntegrationPoints[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, point_number, alpha_eas, zeta_gauss);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(general_variables,Values,point_number);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->CalculateMaterialResponseCauchy(Values);

            if( rOutput[point_number].size2() != general_variables.ConstitutiveMatrix.size2() )
            {
                rOutput[point_number].resize( general_variables.ConstitutiveMatrix.size1() , general_variables.ConstitutiveMatrix.size2() , false );
            }
            rOutput[point_number] = general_variables.ConstitutiveMatrix;
        }
    } else if ( rVariable == DEFORMATION_GRADIENT ) {
        // Create and initialize element variables:
        GeneralVariables general_variables;
        this->InitializeGeneralVariables(general_variables);

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        double& alpha_eas = this->GetValue(ALPHA_EAS);

        /* Calculate the cartesian derivatives */
        CartesianDerivatives this_cartesian_derivatives;
        this->CalculateCartesianDerivatives(this_cartesian_derivatives);

        /* Calculate common components (B, C) */
        CommonComponents common_components;
        common_components.clear();
        this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < IntegrationPoints.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * IntegrationPoints[point_number].Z() - 1.0;

            // Compute element kinematics C, F ...
            this->CalculateKinematics(general_variables, common_components, point_number, alpha_eas, zeta_gauss);

            if( rOutput[point_number].size2() != general_variables.F.size2() )
                rOutput[point_number].resize( general_variables.F.size1() , general_variables.F.size2() , false );
            rOutput[point_number] = general_variables.F;
        }
    } else {
        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
    }

    if ( rOutput.size() != 6 ) {
        std::vector<Matrix> rOutput_aux;
        rOutput_aux = rOutput;

        rOutput.resize( 6 );
        Matrix interpol = StructuralMechanicsMathUtilities::InterpolPrismGiD(integration_point_number);

        for (IndexType iii = 0; iii < 6; iii++) {
            rOutput[iii] = ZeroMatrix(rOutput[0].size1(), rOutput[0].size2());

            for (IndexType Gauss_Point = 0; Gauss_Point < integration_point_number; Gauss_Point++)
                rOutput[iii] += interpol(Gauss_Point, iii) * rOutput_aux[Gauss_Point];
        }
    }

    KRATOS_CATCH( "" );
}

//**************************** ON INTEGRATION POINTS ******************************//
/******************************** SET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++ ) {
        if (rVariable == DETERMINANT_F)
            mAuxCont[point_number] = rValues[point_number];

        mConstitutiveLawVector[point_number]->SetValue( rVariable, rValues[point_number], rCurrentProcessInfo );
    }
}

/******************************** SET VECTOR VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++ ) {
        mConstitutiveLawVector[point_number]->SetValue( rVariable, rValues[point_number], rCurrentProcessInfo );
    }
}

/******************************** SET MATRIX VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); point_number++ )
    {
        mConstitutiveLawVector[point_number]->SetValue( rVariable, rValues[point_number], rCurrentProcessInfo );
    }
}

/****************************** SET CONSTITUTIVE VALUE *****************************/
/***********************************************************************************/

void SprismElement3D6N::SetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rVariable == CONSTITUTIVE_LAW) {
        if ( mConstitutiveLawVector.size() != rValues.size() ) {
            mConstitutiveLawVector.resize(rValues.size());

            KRATOS_ERROR_IF(mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod )) << "Constitutive law not has the correct size " << mConstitutiveLawVector.size() << std::endl;
        }
        for(IndexType i = 0; i < rValues.size(); i++)
            mConstitutiveLawVector[i] = rValues[i];
    }
}

/******************************** GET DOUBLE VALUE *********************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if ( rVariable == VON_MISES_STRESS || rVariable == NORM_ISOCHORIC_STRESS ) {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    } else {
        const IndexType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number, false );
        for ( IndexType ii = 0; ii < integration_points_number; ii++ ) {
            if (rVariable == DETERMINANT_F)
                rValues[ii] = mAuxCont[ii];
            else
                rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }
    }
}

/********************************** GET VECTOR VALUE *******************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const IndexType integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );

    if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR ) {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    } else if ( rVariable == PK2_STRESS_VECTOR ||  rVariable == CAUCHY_STRESS_VECTOR ) {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    } else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR ||  rVariable == HENCKY_STRAIN_TENSOR) {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    } else {
        for ( IndexType point_number = 0;  point_number < integration_points_number; point_number++ ) {
            rValues[point_number] = mConstitutiveLawVector[point_number]->GetValue( rVariable, rValues[point_number] );
        }
    }
}

/*********************************** GET MATRIX VALUE ******************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const IndexType integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );

    if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR ) {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    } else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR ||  rVariable == HENCKY_STRAIN_TENSOR) {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    } else {
        for ( IndexType point_number = 0;  point_number < integration_points_number; point_number++ )
            rValues[point_number] = mConstitutiveLawVector[point_number]->GetValue( rVariable, rValues[point_number] );
    }
}

/******************************** GET CONSTITUTIVE VALUE ***************************/
/***********************************************************************************/

void SprismElement3D6N::GetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rVariable == CONSTITUTIVE_LAW) {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());
        for(IndexType i = 0; i < rValues.size(); i++)
            rValues[i] = mConstitutiveLawVector[i];
    }
}

//********************************* CHECK VALUES **********************************//
/***********************************************************************************/
/***********************************************************************************/

int  SprismElement3D6N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    /* Check the neighbours have been calculated */
    // Neighbour elements
    WeakPointerVector< Element >& p_neighbour_elements = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(p_neighbour_elements.size() == 0) << "The neighbour elements are not calculated" << std::endl;

    // Neighbour nodes
    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    KRATOS_ERROR_IF(p_neighbour_nodes.size() == 0) << "The neighbour nodes are not calculated" << std::endl;

    // Verify that nodal variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(VON_MISES_STRESS)
    KRATOS_CHECK_VARIABLE_KEY(NORM_ISOCHORIC_STRESS)
    KRATOS_CHECK_VARIABLE_KEY(CAUCHY_STRESS_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(CAUCHY_STRESS_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(PK2_STRESS_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(PK2_STRESS_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(GREEN_LAGRANGE_STRAIN_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(GREEN_LAGRANGE_STRAIN_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(ALMANSI_STRAIN_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(ALMANSI_STRAIN_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(HENCKY_STRAIN_TENSOR)
    KRATOS_CHECK_VARIABLE_KEY(HENCKY_STRAIN_VECTOR)
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_MATRIX)
    KRATOS_CHECK_VARIABLE_KEY(DEFORMATION_GRADIENT)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < this->GetGeometry().size(); i++ ) {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() == 6) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    /* Verify compatibility with the constitutive law */
    ConstitutiveLaw::Features law_features;
    this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetLawFeatures(law_features);

    // Check strain measure
    for(IndexType i = 0; i < law_features.mStrainMeasures.size(); i++) {
        KRATOS_ERROR_IF_NOT(law_features.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient) << "Constitutive law is not compatible with the element type SprismElement3D6N" << std::endl;
    }

    // Check constitutive law
    for (IndexType i = 0; i < mConstitutiveLawVector.size(); i++)
        return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo);

    return 0;

    KRATOS_CATCH( "" );
}

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    for ( IndexType i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                rCurrentProcessInfo );

    mFinalizedStep = false;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Create and initialize element variables:
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Get constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    double& AlphaEAS = this->GetValue(ALPHA_EAS);

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    // Reading integration points
    for ( IndexType PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
    {
        const double ZetaGauss = 2.0 * IntegrationPoints[PointNumber].Z() - 1.0;

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, PointNumber, AlphaEAS, ZetaGauss);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(general_variables,Values,PointNumber);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse(Values, general_variables.StressMeasure);

        // Call the constitutive law to finalize the solution step
        mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), PointNumber ),
                rCurrentProcessInfo );

        // Call the element internal variables update
        this->FinalizeStepVariables(general_variables, PointNumber);
    }

    mFinalizedStep = true;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add something if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add something if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::Initialize()
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    /* Constitutive Law initialisation */
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    /* Implicit or explicit EAS update */
    if( GetProperties().Has(EAS_IMP) )
        mELementalFlags.Set(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT, GetProperties()[EAS_IMP]);
    else
        mELementalFlags.Set(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT, true);

    /* Total or updated lagrangian */
    if( GetProperties().Has(SPRISM_TL_UL) )
        mELementalFlags.Set(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN, GetProperties()[SPRISM_TL_UL]);
    else
        mELementalFlags.Set(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN, true);

    /* Quadratic or linear element */
    if( GetProperties().Has(QUAD_ON) )
        mELementalFlags.Set(SprismElement3D6N::QUADRATIC_ELEMENT, GetProperties()[QUAD_ON]);
    else
        mELementalFlags.Set(SprismElement3D6N::QUADRATIC_ELEMENT, true);

    // Resizing the containers
    mAuxMatCont.resize( integration_points.size() );
    mAuxCont.resize( integration_points.size(), false );

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN)) { // Jacobian inverses
        // Compute jacobian inverses and set the domain initial size:
        GeometryType::JacobiansType J0;
        J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);
        mTotalDomainInitialSize = 0.0;

        /* Calculating the inverse J0 */
        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
            // Calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[point_number], mAuxMatCont[point_number], mAuxCont[point_number] );

            // Getting informations for integration
            const double integration_weight = integration_points[point_number].Weight();

            // Calculating the total volume
            mTotalDomainInitialSize += mAuxCont[point_number] * integration_weight;
        }
    } else { // Historic deformation gradient
        mTotalDomainInitialSize = 0.0; // Just initialize, not used in UL

        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
            mAuxCont[point_number] = 1.0;
            mAuxMatCont[point_number] = IdentityMatrix(3);
        }
    }

    /* Initialize AlphaEAS */
    this->SetValue(ALPHA_EAS, 0.0);

    /* Initialize EAS parameters*/
    mEAS.clear();

    /* Material initialisation */
    InitializeMaterial();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// ------------------------------------------------------------------------- //
// ----------------------------- PROTECTED --------------------------------- //
// ------------------------------------------------------------------------- //

void SprismElement3D6N::CalculateElementalSystem(
    LocalSystemComponents& rLocalSystem,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    /* Create and initialize element variables: */
    GeneralVariables general_variables;
    this->InitializeGeneralVariables(general_variables);

    /* Create constitutive law parameters: */
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    /* Set constitutive law flags: */
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    /* Reading integration points */
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    /* Getting the alpha parameter of the EAS improvement */
    double& alpha_eas = this->GetValue(ALPHA_EAS);

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR)) { // Update just if RHS is calculated
        /* Getting the increase of displacements */
        bounded_matrix<double, 36, 1 > delta_disp;

        delta_disp = GetVectorCurrentPosition() - GetVectorPreviousPosition(); // Calculates the increase of displacements

        /* Update alpha EAS */
        if (mEAS.mStiffAlpha > std::numeric_limits<double>::epsilon()) // Avoid division by zero
            alpha_eas -= prod(mEAS.mHEAS, delta_disp)(0, 0) / mEAS.mStiffAlpha;
    }

    /* Calculate the cartesian derivatives */
    CartesianDerivatives this_cartesian_derivatives;
    this->CalculateCartesianDerivatives(this_cartesian_derivatives);

    /* Calculate common components (B, C) */
    CommonComponents common_components;
    common_components.clear();
    this->CalculateCommonComponents(common_components, this_cartesian_derivatives);

    /* Reset the integrated stress components */
    StressIntegratedComponents rIntegratedStress;
    rIntegratedStress.clear();

    /* Reset the EAS integrated components */
    mEAS.clear();

    /* Auxiliary terms: Allocating the VolumeForce*/
    Vector volume_force = ZeroVector(3);

    // Reading integration points
    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

        /* Assemble B */
        this->CalculateDeformationMatrix(general_variables.B, common_components, zeta_gauss, alpha_eas);

        // Compute element kinematics C, F ...
        this->CalculateKinematics(general_variables, common_components, point_number, alpha_eas, zeta_gauss);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(general_variables, Values, point_number);

        // Compute stresses and constitutive parameters
        mConstitutiveLawVector[point_number]->CalculateMaterialResponse(Values, general_variables.StressMeasure);

        // Calculating weights for integration on the "reference configuration"
        const double integration_weight = integration_points[point_number].Weight() * general_variables.detJ;

        /* Integrate in Zeta */
        IntegrateInZeta(general_variables, rIntegratedStress, alpha_eas, zeta_gauss, integration_weight);

        if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) ) { // Calculation of the vector is required
        /* Volume forces */
            this->CalculateVolumeForce( volume_force, general_variables, integration_weight );
        }
    }

    /* Calculate the RHS */
    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) ) { // Calculation of the vector is required
        /* Contribution to external and internal forces */
        this->CalculateAndAddRHS ( rLocalSystem, general_variables, volume_force, rIntegratedStress, common_components, alpha_eas );
    }

    if ( rLocalSystem.CalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) ) { // Calculation of the matrix is required
        /* Contribution to the tangent stiffness matrix */
        this->CalculateAndAddLHS( rLocalSystem, general_variables, Values, rIntegratedStress, common_components, this_cartesian_derivatives, alpha_eas );
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::PrintElementCalculation(
    LocalSystemComponents& rLocalSystem,
    GeneralVariables& rVariables
    )
{
    KRATOS_TRY;

    std::cout << " Element: " << this->Id() << std::endl;

    WeakPointerVector< NodeType >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    const IndexType number_of_neighbours = NumberOfActiveNeighbours(NeighbourNodes);

    for ( IndexType i = 0; i < 6; i++ ) {
        const array_1d<double, 3> &current_position  = GetGeometry()[i].Coordinates();
        const array_1d<double, 3 > & current_displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & previous_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        const array_1d<double, 3> previous_position  = current_position - (current_displacement-previous_displacement);
        std::cout << " Previous  Position  node[" << GetGeometry()[i].Id() << "]: "<<previous_position << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; i++ ) {
        const array_1d<double, 3> &current_position  = NeighbourNodes[i].Coordinates();
        const array_1d<double, 3 > & current_displacement  = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & previous_displacement = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        const array_1d<double, 3> previous_position  = current_position - (current_displacement-previous_displacement);
        std::cout << " Previous  Position  neighbour node[" << NeighbourNodes[i].Id() << "]: "<<previous_position << std::endl;
    }

    for ( IndexType i = 0; i < 6; i++ ) {
        const array_1d<double, 3> & current_position  = GetGeometry()[i].Coordinates();
        std::cout << " Current  Position  node[" << GetGeometry()[i].Id()<<"]: " << current_position << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; i++ ) {
        const array_1d<double, 3> & current_position  = NeighbourNodes[i].Coordinates();
        std::cout << " Current  Position neighbour node[" << NeighbourNodes[i].Id()<<"]: " << current_position << std::endl;
    }

    for ( IndexType i = 0; i < 6; i++ ) {
        const array_1d<double, 3 > & previous_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        std::cout << " Previous Displacement node[" << GetGeometry()[i].Id() << "]: " << previous_displacement << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; i++ ) {
        const array_1d<double, 3 > & previous_displacement = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT,1);
        std::cout << " Previous Displacement neighbour node[" << NeighbourNodes[i].Id() << "]: " << previous_displacement << std::endl;
    }

    for ( IndexType i = 0; i < 6; i++ ) {
        const array_1d<double, 3 > & current_displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        std::cout << " Current  Displacement  node[" << GetGeometry()[i].Id() << "]: " << current_displacement << std::endl;
    }

    for ( IndexType i = 0; i < number_of_neighbours; i++ ) {
        const array_1d<double, 3 > & current_displacement  = NeighbourNodes[i].FastGetSolutionStepValue(DISPLACEMENT);
        std::cout << " Current  Displacement  node[" << NeighbourNodes[i].Id() << "]: " << current_displacement << std::endl;
    }

    std::cout << " Stress " << rVariables.StressVector << std::endl;
    std::cout << " Strain " << rVariables.StrainVector << std::endl;
    std::cout << " F  " << rVariables.F<<std::endl;
    std::cout << " ConstitutiveMatrix " <<rVariables.ConstitutiveMatrix << std::endl;
    std::cout << " K " << rLocalSystem.GetLeftHandSideMatrix() << std::endl;
    std::cout << " f " << rLocalSystem.GetRightHandSideVector() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bool SprismElement3D6N::HasNeighbour(
    const IndexType Index,
    const NodeType& NeighbourNode
    )
{
    if (NeighbourNode.Id() == GetGeometry()[Index].Id()) {
        return false;
    } else {
        if ( mELementalFlags.Is(SprismElement3D6N::QUADRATIC_ELEMENT) == true )
            return true;
        else
            return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t SprismElement3D6N::NumberOfActiveNeighbours(WeakPointerVector< NodeType >& pNeighbourNodes)
{
    std::size_t active_neighbours = 0;
    for (IndexType i = 0; i < pNeighbourNodes.size(); i++) {
        if (HasNeighbour(i, pNeighbourNodes[i]))
            active_neighbours++;
    }
    return active_neighbours;
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::GetNodalCoordinates(
    bounded_matrix<double, 12, 3 > & NodesCoord,
    WeakPointerVector< NodeType >& NeighbourNodes,
    const Configuration ThisConfiguration
    )
{
     NodesCoord = ZeroMatrix(12, 3);
     const IndexType number_of_neighbours = NumberOfActiveNeighbours(NeighbourNodes);

     if (ThisConfiguration == Configuration::INITIAL) {
         /* Fill the aux matrix of coordinates */
         for (IndexType i = 0; i < 6; i++) {
             const array_1d<double, 3> &initial_position = GetGeometry()[i].GetInitialPosition().Coordinates();
             for (IndexType j = 0; j < 3; j++)
                 NodesCoord(i, j) = initial_position[j];
         }

         if (number_of_neighbours == 6) { // All the possible neighours
             for (IndexType i = 0; i < 6; i++) {
                 const array_1d<double, 3> &initial_position = NeighbourNodes[i].GetInitialPosition().Coordinates();
                 for (IndexType j = 0; j < 3; j++)
                    NodesCoord(i + 6, j) = initial_position[j];
             }
         } else {
             for (IndexType i = 0; i < 6; i++) {
                 if (HasNeighbour(i, NeighbourNodes[i])) {
                     const array_1d<double, 3> &initial_position = NeighbourNodes[i].GetInitialPosition().Coordinates();

                     for (IndexType j = 0; j < 3; j++)
                        NodesCoord(i + 6, j) = initial_position[j];

                 } else {
                     for (IndexType j = 0; j < 3; j++)
                        NodesCoord(i + 6, j) = 0.0;
                 }
             }
         }
     } else if (ThisConfiguration == Configuration::CURRENT) {
         /* Fill the aux matrix of coordinates */
         for (IndexType i = 0; i < 6; i++) {
             const array_1d<double, 3> &current_position = GetGeometry()[i].Coordinates();
             for (IndexType j = 0; j < 3; j++)
                NodesCoord(i, j) = current_position[j];
         }

         if (number_of_neighbours == 6) { // All the possible neighours
             for (IndexType i = 0; i < 6; i++) {
                 const array_1d<double, 3> &current_position = NeighbourNodes[i].Coordinates();
                 for (IndexType j = 0; j < 3; j++)
                    NodesCoord(i + 6, j) = current_position[j];
             }
         } else {
             for (IndexType i = 0; i < 6; i++) {
                 if (HasNeighbour(i, NeighbourNodes[i])) {
                     const array_1d<double, 3> &current_position = NeighbourNodes[i].Coordinates();
                     for (IndexType j = 0; j < 3; j++)
                        NodesCoord(i + 6, j) = current_position[j];
                 } else {
                     for (IndexType j = 0; j < 3; j++)
                        NodesCoord(i + 6, j) = 0.0;
                 }
             }
         }
     } else {
         const std::string& config = (ThisConfiguration == Configuration::INITIAL) ? "Initial" : "Current";
         KRATOS_ERROR << " The configuration is not possible, the posibilities are Current and Initial: " << config << std::endl;
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerivatives(CartesianDerivatives& this_cartesian_derivatives)
{
    bounded_matrix<double, 12, 3 > nodes_coord; // Coordinates of the nodes
    WeakPointerVector< NodeType >& neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true ) {
        this->GetNodalCoordinates(nodes_coord, neighbour_nodes, Configuration::INITIAL);
    } else {
        this->GetNodalCoordinates(nodes_coord, neighbour_nodes, Configuration::CURRENT);
    }

    /* Calculate local system of coordinates of the element */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    const double ang_rot = GetProperties().Has(ANG_ROT) ? GetProperties()[ANG_ROT] : 0.0; // TODO: Change to consider multiple plies
    this->CalculateLocalCoordinateSystem(2, ang_rot);

    //******************************** CENTRAL POINT ******************************
    // Calculate cartesian derivatives
    bounded_matrix<double, 2, 4 > cartesian_derivatives_center_lower;
    bounded_matrix<double, 2, 4 > cartesian_derivatives_center_upper;

    // Lower face
    CalculateCartesianDerOnCenterPlane(0, nodes_coord, cartesian_derivatives_center_lower);
    // Upperr face
    CalculateCartesianDerOnCenterPlane(3, nodes_coord, cartesian_derivatives_center_upper );

    /* Transversal derivative */
    CalculateCartesianDerOnCenterTrans(this_cartesian_derivatives, nodes_coord, 0); // Center
    CalculateCartesianDerOnCenterTrans(this_cartesian_derivatives, nodes_coord, 1); // Lower part
    CalculateCartesianDerOnCenterTrans(this_cartesian_derivatives, nodes_coord, 2); // Upper part

    //******************************** GAUSS POINTS *******************************

    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 0.5;
    local_coordinates[1] = 0.5;
    local_coordinates[2] = -1.0;

    /* Transversal derivative */
    CalculateCartesianDerOnGaussTrans(nodes_coord, this_cartesian_derivatives.TransversalCartesianDerivativesGauss[0], local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(nodes_coord, this_cartesian_derivatives.TransversalCartesianDerivativesGauss[3], local_coordinates);

    /* In-plane derivative */
    if (HasNeighbour(0, neighbour_nodes[0])) { // Assuming that if the upper element has neighbours the lower has too
        CalculateCartesianDerOnGaussPlane(0, 0, nodes_coord, this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[0]);
        CalculateCartesianDerOnGaussPlane(0, 3, nodes_coord, this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[3]);
    } else {
        noalias(this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[0]) = cartesian_derivatives_center_lower;
        noalias(this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[3]) = cartesian_derivatives_center_upper;
    }

    /* Transversal derivative */
    local_coordinates[0] = 0.0;
    local_coordinates[2] = -1.0;
    CalculateCartesianDerOnGaussTrans(nodes_coord, this_cartesian_derivatives.TransversalCartesianDerivativesGauss[1], local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(nodes_coord, this_cartesian_derivatives.TransversalCartesianDerivativesGauss[4], local_coordinates);

    /* In-plane derivative */
    if (HasNeighbour(1, neighbour_nodes[1])) { //Idem
        CalculateCartesianDerOnGaussPlane(1, 0, nodes_coord, this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[1]);
        CalculateCartesianDerOnGaussPlane(1, 3, nodes_coord, this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[4]);
    } else {
        noalias(this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[1]) = cartesian_derivatives_center_lower;
        noalias(this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[4]) = cartesian_derivatives_center_upper;
    }

    /* Transversal derivative */
    local_coordinates[0] = 0.5;
    local_coordinates[1] = 0.0;
    local_coordinates[2] = -1.0;
    CalculateCartesianDerOnGaussTrans(nodes_coord, this_cartesian_derivatives.TransversalCartesianDerivativesGauss[2], local_coordinates);
    local_coordinates[2] = 1.0;
    CalculateCartesianDerOnGaussTrans(nodes_coord, this_cartesian_derivatives.TransversalCartesianDerivativesGauss[5], local_coordinates);

    /* In-plane derivative */
    if (HasNeighbour(2, neighbour_nodes[2])) { // Idem
        CalculateCartesianDerOnGaussPlane(2, 0, nodes_coord, this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[2]);
        CalculateCartesianDerOnGaussPlane(2, 3, nodes_coord, this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[5]);
    } else {
        noalias(this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[2]) = cartesian_derivatives_center_lower;
        noalias(this_cartesian_derivatives.InPlaneCartesianDerivativesGauss[5]) = cartesian_derivatives_center_upper;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCommonComponents(
    CommonComponents& rCommonComponents,
    const CartesianDerivatives& rCartesianDerivatives
    )
{
    KRATOS_TRY;

    bounded_matrix<double, 12, 3 > NodesCoord; // Coordinates of the nodes
    WeakPointerVector< NodeType >& NeighbourNodes = this->GetValue(NEIGHBOUR_NODES);
    this->GetNodalCoordinates(NodesCoord, NeighbourNodes, Configuration::CURRENT);

    /* Declare deformation Gradient F components */
    // In plane components
    bounded_matrix<double, 3, 2 > in_plane_gradient_F_gauss;
    // Transversal components
    TransverseGradient transverse_gradient;

    //*****************************************************************************

    /* COMPUTATION OF B TANGENTS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Calculating the membrane strain-displacement matrix */
    // Lower face

    // Gauss point 1
    CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[0], NodesCoord, 0, 0);
    CalculateAndAddBMembrane(rCommonComponents.BMembraneLower, rCommonComponents.CMembraneLower, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[0], in_plane_gradient_F_gauss, 0);

    // Gauss point 2
    CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[1], NodesCoord, 1, 0);
    CalculateAndAddBMembrane(rCommonComponents.BMembraneLower, rCommonComponents.CMembraneLower, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[1], in_plane_gradient_F_gauss, 1);

    // Gauss point 3
    CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[2], NodesCoord, 2, 0);
    CalculateAndAddBMembrane(rCommonComponents.BMembraneLower, rCommonComponents.CMembraneLower, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[2], in_plane_gradient_F_gauss, 2);

    rCommonComponents.BMembraneLower *= 1.0/3.0;
    rCommonComponents.CMembraneLower *= 1.0/3.0;

    // Upper face

    // Gauss point 4
    CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[3], NodesCoord, 0, 3);
    CalculateAndAddBMembrane(rCommonComponents.BMembraneUpper, rCommonComponents.CMembraneUpper, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[3], in_plane_gradient_F_gauss, 0);

    // Gauss point 5
    CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[4], NodesCoord, 1, 3);
    CalculateAndAddBMembrane(rCommonComponents.BMembraneUpper, rCommonComponents.CMembraneUpper, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[4], in_plane_gradient_F_gauss, 1);

    // Gauss point 6
    CalculateInPlaneGradientFGauss(in_plane_gradient_F_gauss, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[5], NodesCoord, 2, 3);
    CalculateAndAddBMembrane(rCommonComponents.BMembraneUpper, rCommonComponents.CMembraneUpper, rCartesianDerivatives.InPlaneCartesianDerivativesGauss[5], in_plane_gradient_F_gauss, 2);

    rCommonComponents.BMembraneUpper *= 1.0/3.0;
    rCommonComponents.CMembraneUpper *= 1.0/3.0;

    /* SHEAR CONTRIBUTION */
    /* Calculating the shear strain-displacement matrix */

    // Declaring the isoparametric transverse gradient variables
    TransverseGradientIsoParametric transverse_gradient_isoparametric;

    // Lower face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(transverse_gradient_isoparametric, NodesCoord, 0);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesGauss[0], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F1, rCartesianDerivatives.TransversalCartesianDerivativesGauss[1], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F2, rCartesianDerivatives.TransversalCartesianDerivativesGauss[2], NodesCoord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAddBShear(rCommonComponents.BShearLower, rCommonComponents.CShearLower, rCartesianDerivatives, transverse_gradient, transverse_gradient_isoparametric, 0);

    // Upper face

    /* Calculate f components in the face */
    CalculateTransverseGradientFinP(transverse_gradient_isoparametric, NodesCoord, 3);

    /* Calculate f transverse components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesGauss[3], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F1, rCartesianDerivatives.TransversalCartesianDerivativesGauss[4], NodesCoord);
    CalculateTransverseGradientF(transverse_gradient.F2, rCartesianDerivatives.TransversalCartesianDerivativesGauss[5], NodesCoord);

    /* Shear contribution to the deformation matrix */
    CalculateAndAddBShear(rCommonComponents.BShearUpper, rCommonComponents.CShearUpper, rCartesianDerivatives, transverse_gradient, transverse_gradient_isoparametric, 9);

    /* NORMAL TRANSVERSE */
    /* Calculate f normal components */
    CalculateTransverseGradientF(transverse_gradient.F0, rCartesianDerivatives.TransversalCartesianDerivativesCenter, NodesCoord);

    /* Calculating the normal transverse strain-displacement matrix */
    CalculateAndAddBNormal(rCommonComponents.BNormal, rCommonComponents.CNormal, rCartesianDerivatives.TransversalCartesianDerivativesCenter, transverse_gradient.F0);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateLocalCoordinateSystem(
    const int choose,
    const double ang
    )
{
    KRATOS_TRY;

    /* Mid-surface vectors */ // TODO: Consider CURRENT or INITIAL configuration depending
    double norm; // TODO: Use the geometry normal when avalaible
    array_1d<double, 3 > vxe, vye;
    vxe[0] = 0.5 * ((GetGeometry()[2].X0() + GetGeometry()[5].X0()) - (GetGeometry()[1].X0() + GetGeometry()[4].X0()));
    vxe[1] = 0.5 * ((GetGeometry()[2].Y0() + GetGeometry()[5].Y0()) - (GetGeometry()[1].Y0() + GetGeometry()[4].Y0()));
    vxe[2] = 0.5 * ((GetGeometry()[2].Z0() + GetGeometry()[5].Z0()) - (GetGeometry()[1].Z0() + GetGeometry()[4].Z0()));

    vye[0] = 0.5 * ((GetGeometry()[0].X0() + GetGeometry()[3].X0()) - (GetGeometry()[2].X0() + GetGeometry()[5].X0()));
    vye[1] = 0.5 * ((GetGeometry()[0].Y0() + GetGeometry()[3].Y0()) - (GetGeometry()[2].Y0() + GetGeometry()[5].Y0()));
    vye[2] = 0.5 * ((GetGeometry()[0].Z0() + GetGeometry()[3].Z0()) - (GetGeometry()[2].Z0() + GetGeometry()[5].Z0()));

    MathUtils<double>::CrossProduct(mvze, vxe, vye);
    norm = norm_2(mvze);
    mvze /= norm;

    double threshold = 1e-5;
    double OrthoComp;

    /* Performing the calculation */
    // 0- If X is the prefered normal vector
    // 1- If Y is the prefered normal vector
    // 2- If Z is the prefered normal vector

    if (choose == 0) {
        OrthoComp = mvze[1] * mvze[1] + mvze[2] * mvze[2]; // Component in th Y-Z plane
        if (OrthoComp < threshold) { // If mvze is almost orthogonal to  Y-Z plane
            mvye[0] = - mvze[2]; // Choose mvxe orthogonal to global Y direction
            mvye[1] = 0.0;
            mvye[2] = mvze[0];

            norm = norm_2(mvxe);
            mvxe /= norm;
            MathUtils<double>::CrossProduct(mvxe, mvye, mvze);
        } else { // SELECT local y=mvxe in the global YZ plane
            mvxe[0] = 0.0;
            mvxe[1] = mvze[2];
            mvxe[2] = - mvze[1];

            norm = norm_2(mvxe);
            mvxe /= norm;

            mvye[0] = OrthoComp; // Choose mvxe orthogonal to global X direction
            mvye[1] = - mvze[0] * mvze[1];
            mvye[2] = - mvze[0] * mvze[2];

            norm = norm_2(mvye);
            mvye /= norm;
        }
    } else if (choose == 1) {
        OrthoComp = mvze[0] * mvze[0] + mvze[2] * mvze[2]; // Component in th Z-X plane
        if (OrthoComp < threshold) { // If vze is almost orthogonal to  Z-X plane
            mvye[0] =       0.0; // Choose mvxe orthogonal to global X direction
            mvye[1] =   mvze[2];
            mvye[2] = - mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
            MathUtils<double>::CrossProduct(mvxe, mvye, mvze);
        } else { // SELECT local z=mvxe in the global ZX plane
            mvxe[0] = - mvze[2]; // Choose mvxe orthogonal to global Y direction
            mvxe[1] = 0.0;
            mvxe[2] = - mvze[0];

            norm = norm_2(mvxe);
            mvxe /= norm;

            mvye[0] = - mvze[0] * mvze[1];
            mvye[1] = OrthoComp;
            mvye[2] = - mvze[2] * mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
        }
    } else if (choose == 2) {
        OrthoComp = mvze[0] * mvze[0] + mvze[1] * mvze[1]; // Component in th X-Y plane
        if (OrthoComp < threshold) { // If vze is almost orthogonal to  X-Y plane
            mvye[0] = 0.0; // Choose mvxe orthogonal to global X direction
            mvye[1] = mvze[2];
            mvye[2] = - mvze[1];

            norm = norm_2(mvye);
            mvye /= norm;
            MathUtils<double>::CrossProduct(mvxe, mvye, mvze);
        } else { // SELECT local x=mvxe in the global XY plane
            mvxe[0] = - mvze[1];
            mvxe[1] = mvze[0];
            mvxe[2] = 0.0;

            norm = norm_2(mvxe);
            mvxe /= norm;

            mvye[0] = - mvze[0] * mvze[2]; // Choose mvxe orthogonal to global Z direction
            mvye[1] = - mvze[1] * mvze[2];
            mvye[2] = OrthoComp;

            norm = norm_2(mvye);
            mvye /= norm;
        }
    } else {
        mvxe[0] = 1.0;
        mvxe[1] = 0.0;
        mvxe[2] = 0.0;

        mvye[0] = 0.0;
        mvye[1] = 1.0;
        mvye[2] = 0.0;
    }

    if (ang != 0.0) {
        // Compute angle between local system mvxe-mvye and L1
        const double cosa = std::cos(ang);
        const double sina = std::sin(ang);
        // Rotate local system mvxe-mvye to best fit L1-L2
        mvze = mvxe; // Reusing as auxiliar value
        mvxe =   cosa * mvxe + sina * mvye;
        mvye = - sina * mvze  + cosa * mvye;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateIdVector(array_1d<IndexType, 18 >& rIdVector)
{
    KRATOS_TRY;

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Compute ID vector */ // TODO: Optimze this
    IndexType index = 18;
    for (IndexType i = 0; i < 6; i++) {
        if (HasNeighbour(i, p_neighbour_nodes[i])) {
            for (IndexType j = 0; j < 3; j++)
                rIdVector[i * 3 + j] = index + j;
            index += 3;
        } else {
            for (IndexType j = 0; j < 3; j++)
                rIdVector[i * 3 + j] = 36;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ComputeLocalDerivatives(
    bounded_matrix<double, 6, 3 > & LocalDerivativePatch,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    const double L_1 = 0.5 * (1.0 - rLocalCoordinates[2]);
    const double L_2 = 0.5 * (1.0 + rLocalCoordinates[2]);

    /* Derivative in direction nu and xi */
    // Lower face
    LocalDerivativePatch(0, 0) = - L_1;
    LocalDerivativePatch(1, 0) =   L_1;
    LocalDerivativePatch(2, 0) =   0.0;

    LocalDerivativePatch(0, 1) = - L_1;
    LocalDerivativePatch(1, 1) =   0.0;
    LocalDerivativePatch(2, 1) =   L_1;

    // Upper face
    LocalDerivativePatch(3, 0) = - L_2;
    LocalDerivativePatch(4, 0) =   L_2;
    LocalDerivativePatch(5, 0) =   0.0;

    LocalDerivativePatch(3, 1) = - L_2;
    LocalDerivativePatch(4, 1) =   0.0;
    LocalDerivativePatch(5, 1) =   L_2;

    /* Derivative in direction zeta */
    LocalDerivativePatch(0, 2) = - 1.0 + rLocalCoordinates[1] + rLocalCoordinates[0];
    LocalDerivativePatch(1, 2) = - rLocalCoordinates[0];
    LocalDerivativePatch(2, 2) = - rLocalCoordinates[1];
    LocalDerivativePatch(3, 2) =   1.0 - rLocalCoordinates[1] - rLocalCoordinates[0];
    LocalDerivativePatch(4, 2) =   rLocalCoordinates[0];
    LocalDerivativePatch(5, 2) =   rLocalCoordinates[1];
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ComputeLocalDerivativesQuadratic(
    bounded_matrix<double, 4, 2 >& rLocalDerivativePatch,
    const IndexType NodeGauss
    )
{
    /* Local coordinates */
    double xi  = 0.0;
    double eta = 0.0;

    if (NodeGauss == 0) {
        xi  = 0.5;
        eta = 0.5;
    } else if (NodeGauss == 1) {
        xi  = 0.0;
        eta = 0.5;
    } else if (NodeGauss == 2) {
        xi  = 0.5;
        eta = 0.0;
    }

    /* Derivative in main nodes */
    rLocalDerivativePatch(0, 0) = - 1.0 + eta;
    rLocalDerivativePatch(0, 1) = - 1.0 + xi;
    rLocalDerivativePatch(1, 0) =   1.0 - eta;
    rLocalDerivativePatch(1, 1) =   1.0 - xi - 2.0 * eta;
    rLocalDerivativePatch(2, 0) =   1.0 - 2.0 * xi - eta;
    rLocalDerivativePatch(2, 1) =   1.0 - xi;

    /* Derivative in neighbour nodes */
    if (NodeGauss == 0) {
        rLocalDerivativePatch(3, 0) = xi + eta - 0.5;
        rLocalDerivativePatch(3, 1) = xi + eta - 0.5;
    } else if (NodeGauss == 1) {
        rLocalDerivativePatch(3, 0) = xi - 0.5;
        rLocalDerivativePatch(3, 1) = 0.0;
    } else if (NodeGauss == 2) {
        rLocalDerivativePatch(3, 0) = 0.0;
        rLocalDerivativePatch(3, 1) = eta - 0.5;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianCenterGauss(
    GeometryType::JacobiansType& J,
    std::vector< Matrix >& Jinv,
    Vector& detJ,
    const IndexType rPointNumber,
    const double ZetaGauss
    )
{
    /* Fill the aux matrix of coordinates */
    bounded_matrix<double, 3, 6 > nodes_coord;
    for (IndexType i = 0; i < 6; i++) {
        const array_1d<double, 3> &current_position  = GetGeometry()[i].Coordinates();
        for (IndexType j = 0; j < 3; j++)
            nodes_coord(j, i) = current_position[j];
    }

    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 1.0/3.0;
    local_coordinates[1] = 1.0/3.0;
    local_coordinates[2] = ZetaGauss;

    /* Local derivatives patch */
    bounded_matrix<double, 6, 3 > LocalDerivativePatch;
    ComputeLocalDerivatives(LocalDerivativePatch, local_coordinates);

    /* Compute Jacobian */
    noalias(J[rPointNumber]) = prod(nodes_coord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    MathUtils<double>::InvertMatrix( J[rPointNumber], Jinv[rPointNumber], detJ[rPointNumber] );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobian(
    double & detJ,
    bounded_matrix<double, 3, 3 > & J,
    bounded_matrix<double, 6, 3 > & LocalDerivativePatch,
    const bounded_matrix<double, 12, 3 > & NodesCoord,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Auxiliar coordinates of the nodes */
    bounded_matrix<double, 3, 6 > nodes_coord_aux; // TODO: use just trans

    for (IndexType i = 0; i < 6; i++)
        for (IndexType j = 0; j < 3; j++)
            nodes_coord_aux(j, i) = NodesCoord(i, j);

    /* Local derivatives patch */
    ComputeLocalDerivatives(LocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(J) = prod(nodes_coord_aux, LocalDerivativePatch);

    /* Compute determinant */
    detJ = MathUtils<double>::Det3(J);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianAndInv(
    bounded_matrix<double, 3, 3 >& J,
    bounded_matrix<double, 3, 3 >& Jinv,
    bounded_matrix<double, 6, 3 >& LocalDerivativePatch,
    const bounded_matrix<double, 3, 6 >& NodesCoord,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Local derivatives patch */
    ComputeLocalDerivatives(LocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    Jinv = MathUtils<double>::InvertMatrix<3>(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateJacobianAndInv(
    bounded_matrix<double, 3, 3 > & J,
    bounded_matrix<double, 3, 3 > & Jinv,
    const bounded_matrix<double, 3, 6 > & NodesCoord,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Local derivatives patch */
    bounded_matrix<double, 6, 3 > LocalDerivativePatch;
    ComputeLocalDerivatives(LocalDerivativePatch, rLocalCoordinates);

    /* Compute Jacobian */
    noalias(J) = prod(NodesCoord, LocalDerivativePatch);

    /* Compute inverse of the Jaccobian */
    double detJ;
    Jinv = MathUtils<double>::InvertMatrix<3>(J, detJ);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnCenterPlane(
    const IndexType Index,
    const bounded_matrix<double, 12, 3 > & NodesCoord,
    bounded_matrix<double, 2, 4 > & CartesianDerivativesCenter
    )
{
    double norm0, norm;
    array_1d<double, 3 > vxe, vye;
    vxe[0] = GetGeometry()[2 + Index].X0() - GetGeometry()[1 + Index].X0();
    vxe[1] = GetGeometry()[2 + Index].Y0() - GetGeometry()[1 + Index].Y0();
    vxe[2] = GetGeometry()[2 + Index].Z0() - GetGeometry()[1 + Index].Z0();

    vye[0] = GetGeometry()[0 + Index].X0() - GetGeometry()[2 + Index].X0();
    vye[1] = GetGeometry()[0 + Index].Y0() - GetGeometry()[2 + Index].Y0();
    vye[2] = GetGeometry()[0 + Index].Z0() - GetGeometry()[2 + Index].Z0();

    array_1d<double, 3 > t1g, t2g, t3g;
    MathUtils<double>::CrossProduct(t3g, vxe, vye);
    norm0 = norm_2(t3g);
    t3g /= norm0;

    MathUtils<double>::CrossProduct(t2g, t3g, mvxe);
    norm = norm_2(t2g);
    t2g /= norm;

    MathUtils<double>::CrossProduct(t1g, t2g, t3g);
    norm = norm_2(t1g);
    t1g /= norm;

    array_1d<double, 3 > a, b;

    a[0] = inner_prod(vxe, t1g)/norm0;
    a[1] = inner_prod(vye, t1g)/norm0;
    a[2] = -(a[0] + a[1]);
    b[0] = inner_prod(vxe, t2g)/norm0;
    b[1] = inner_prod(vye, t2g)/norm0;
    b[2] = -(b[0] + b[1]);

    CartesianDerivativesCenter = ZeroMatrix(2, 4);
    for (IndexType i = 0; i < 3; i++) {
       CartesianDerivativesCenter(0, i) = - b[i];
       CartesianDerivativesCenter(1, i) =   a[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnGaussPlane(
    const IndexType NodeGauss,
    const IndexType Index,
    const bounded_matrix<double, 12, 3 > & NodesCoord,
    bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss
    )
{
    /* Local derivatives patch */
    bounded_matrix<double, 4, 2 > local_derivative_patch;
    ComputeLocalDerivativesQuadratic(local_derivative_patch,NodeGauss);

    /* Auxiliar coordinates of the nodes */
    bounded_matrix<double, 3, 4 > nodes_coord_aux;

    for (IndexType i = 0; i < 3; i++)
        for (IndexType j = 0; j < 3; j++)
            nodes_coord_aux(j, i) = NodesCoord(i + Index, j);

    for (IndexType j = 0; j < 3; j++)
        nodes_coord_aux(j, 3) = NodesCoord(NodeGauss + 6 + Index, j);

    /* Compute local derivatives */
    const bounded_matrix<double, 3, 2 > Xd = prod(nodes_coord_aux, local_derivative_patch);

    /* Split local derivatives */
    array_1d<double, 3 > Xdxi, Xdeta;
    Xdxi[0]  = Xd(0, 0);
    Xdxi[1]  = Xd(1, 0);
    Xdxi[2]  = Xd(2, 0);
    Xdeta[0] = Xd(0, 1);
    Xdeta[1] = Xd(1, 1);
    Xdeta[2] = Xd(2, 1);

    /* Compute orthonormal vectors */
    array_1d<double, 3 > t1g, t2g, t3g;
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t1g, t2g, t3g, mvxe, Xdxi, Xdeta);

    /* Compute Jacobian */
    bounded_matrix<double, 2, 2 > jac;
    jac(0, 0) = inner_prod(Xdxi,  t1g);
    jac(0, 1) = inner_prod(Xdxi,  t2g);
    jac(1, 0) = inner_prod(Xdeta, t1g);
    jac(1, 1) = inner_prod(Xdeta, t2g);

    /* Compute the inverse of the Jacobian */
    double AuxDet;
    const bounded_matrix<double, 2, 2 > JinvPlane = MathUtils<double>::InvertMatrix<2>(jac, AuxDet);

    /* Compute the Cartesian derivatives */
    noalias(InPlaneCartesianDerivativesGauss) = prod(JinvPlane, trans(local_derivative_patch));
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnGaussTrans(
    const bounded_matrix<double, 12, 3 > & NodesCoord,
    bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
    const array_1d<double, 3>& rLocalCoordinates
    )
{
    /* Compute local derivatives */
    double det;
    bounded_matrix<double, 3, 3 > Xd;
    bounded_matrix<double, 6, 3 > local_derivatives_patch;
    CalculateJacobian(det, Xd, local_derivatives_patch, NodesCoord, rLocalCoordinates);

    /* Split local derivatives */
    array_1d<double, 3 > Xdxi, Xdeta;
    Xdxi[0]  = Xd(0, 0);
    Xdxi[1]  = Xd(1, 0);
    Xdxi[2]  = Xd(2, 0);
    Xdeta[0] = Xd(0, 1);
    Xdeta[1] = Xd(1, 1);
    Xdeta[2] = Xd(2, 1);

    /* Compute orthonormal vectors */
    array_1d<double, 3 > t1g, t2g, t3g;
    bounded_matrix<double, 3, 3 > t = ZeroMatrix(3, 3);
    StructuralMechanicsMathUtilities::Comp_Orthonor_Base(t, t1g, t2g, t3g, mvxe, Xdxi, Xdeta);

    /* Compute Jacobian */
    bounded_matrix<double, 3, 3 > jac;
    noalias(jac) = prod(t, Xd);

    /* Compute inverse of the Jaccobian (just third column) */
    bounded_matrix<double, 3 ,1> JinvTrans;
    JinvTrans(0, 0) =   (jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1)) / det;
    JinvTrans(1, 0) = - (jac(0, 0) * jac(1, 2) - jac(0, 2) * jac(1, 0)) / det;
    JinvTrans(2, 0) =   (jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1)) / det;

    /* Compute Cartesian derivatives */
    noalias(TransversalCartesianDerivativesGauss) = prod(local_derivatives_patch, JinvTrans);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateCartesianDerOnCenterTrans(
        CartesianDerivatives& this_cartesian_derivatives,
        const bounded_matrix<double, 12, 3 > & NodesCoord,
        const IndexType Part // TODO: Replace enum
        )
{
    array_1d<double, 3> local_coordinates;
    local_coordinates[0] = 1.0/3.0;
    local_coordinates[1] = 1.0/3.0;

    if (Part == 0)
        local_coordinates[2] =   0.0;
    else if (Part == 1)
        local_coordinates[2] = - 1.0;
    else if (Part == 2)
        local_coordinates[2] =   1.0;
    else
        KRATOS_ERROR << " This part id is not possible, just 0, 1 or 2  " << Part << std::endl;

    /* Auxiliar coordinates of the nodes */
    bounded_matrix<double, 3, 6 > nodes_coord_aux;
    for (IndexType i = 0; i < 6; i++)
        for (IndexType j = 0; j < 3; j++)
            nodes_coord_aux(j, i) = NodesCoord(i, j);

    /* Auxiliar components to calculate the Jacobian and his inverse */
    bounded_matrix<double, 3, 3 > J, Jinv;

    if (Part == 0) {
        /* Calculate the Jacobian and his inverse */
        bounded_matrix<double, 6, 3 > local_derivatives_patch;
        CalculateJacobianAndInv(J, Jinv, local_derivatives_patch, nodes_coord_aux, local_coordinates);

        // Compute cartesian (y3) derivatives of the shape functions necessary to compute f_3
        /* Compute Cartesian derivatives */
        bounded_matrix<double, 6, 3 > transverse_cartesian_derivatives_gauss_aux;
        noalias(transverse_cartesian_derivatives_gauss_aux) = prod(local_derivatives_patch, Jinv);

        for (IndexType i = 0; i < 6 ; i++) {
            this_cartesian_derivatives.TransversalCartesianDerivativesCenter(i, 0) =
                     mvze[0] * transverse_cartesian_derivatives_gauss_aux(i, 0) +
                     mvze[1] * transverse_cartesian_derivatives_gauss_aux(i, 1) +
                     mvze[2] * transverse_cartesian_derivatives_gauss_aux(i, 2);
        }
     } else {
        /* Calculate the Jacobian and his inverse */
        CalculateJacobianAndInv(J, Jinv, nodes_coord_aux, local_coordinates);

         /* Split local derivatives */
         array_1d<double, 3 > Xdxi, Xdeta;
         Xdxi[0]   = Jinv(0, 0);
         Xdxi[1]   = Jinv(0, 1);
         Xdxi[2]   = Jinv(0, 2);
         Xdeta[0]  = Jinv(1, 0);
         Xdeta[1]  = Jinv(1, 1);
         Xdeta[2]  = Jinv(1, 2);

         /* Compute inverse of the Jaccobian (just in plane components)*/
         if (Part == 1) {
             this_cartesian_derivatives.JInvPlaneLower(0, 0) = inner_prod(Xdxi,  mvxe);
             this_cartesian_derivatives.JInvPlaneLower(0, 1) = inner_prod(Xdeta, mvxe);
             this_cartesian_derivatives.JInvPlaneLower(1, 0) = inner_prod(Xdxi,  mvye);
             this_cartesian_derivatives.JInvPlaneLower(1, 1) = inner_prod(Xdeta, mvye);
         } else if (Part == 2) {
             this_cartesian_derivatives.JInvPlaneUpper(0, 0) = inner_prod(Xdxi,  mvxe);
             this_cartesian_derivatives.JInvPlaneUpper(0, 1) = inner_prod(Xdeta, mvxe);
             this_cartesian_derivatives.JInvPlaneUpper(1, 0) = inner_prod(Xdxi,  mvye);
             this_cartesian_derivatives.JInvPlaneUpper(1, 1) = inner_prod(Xdeta, mvye);
         }
     }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateInPlaneGradientFGauss(
    bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
    const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
    const bounded_matrix<double, 12, 3 > & NodesCoord,
    const IndexType NodeGauss,
    const IndexType Index
    )
{
    /* Auxiliar operators */
    bounded_matrix<double, 3, 3 > nodes_coord_aux;
    bounded_matrix<double, 3, 2 > in_plane_cartesian_derivatives_gauss_aux;

    for (IndexType i = 0; i < 3; i++) {
        for (IndexType j = 0; j < 3; j++)
            nodes_coord_aux(j, i) = NodesCoord(i + Index, j);

        in_plane_cartesian_derivatives_gauss_aux(i, 0) = InPlaneCartesianDerivativesGauss(0, i);
        in_plane_cartesian_derivatives_gauss_aux(i, 1) = InPlaneCartesianDerivativesGauss(1, i);
    }

    noalias(InPlaneGradientFGauss) = prod(nodes_coord_aux, in_plane_cartesian_derivatives_gauss_aux);

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    if (HasNeighbour(NodeGauss, p_neighbour_nodes[NodeGauss])) {
        for (IndexType j = 0; j < 3 ; j++) {
            InPlaneGradientFGauss(j, 0) += NodesCoord(NodeGauss + 6 + Index, j) * InPlaneCartesianDerivativesGauss(0, 3);
            InPlaneGradientFGauss(j, 1) += NodesCoord(NodeGauss + 6 + Index, j) * InPlaneCartesianDerivativesGauss(1, 3);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateTransverseGradientF(
    array_1d<double, 3 > & TransverseGradientF,
    const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
    const bounded_matrix<double, 12, 3 > & NodesCoord
    )
{
    noalias(TransverseGradientF) = ZeroVector(3);

    for (IndexType i = 0; i < 6; i++) {
        for (IndexType j = 0; j < 3; j++) {
            TransverseGradientF[j] += TransversalCartesianDerivativesGauss(i, 0) * NodesCoord(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateTransverseGradientFinP(
    TransverseGradientIsoParametric& TransverseGradientIsoParametric,
    const bounded_matrix<double, 12, 3 > & NodesCoord,
    const IndexType Index
    )
{
    for (IndexType i = 0; i < 3; i++) {
        TransverseGradientIsoParametric.Ft[i]   = NodesCoord(2 + Index, i) - NodesCoord(1 + Index, i);
        TransverseGradientIsoParametric.Fxi[i]  = NodesCoord(0 + Index, i) - NodesCoord(2 + Index, i);
        TransverseGradientIsoParametric.Feta[i] = NodesCoord(1 + Index, i) - NodesCoord(0 + Index, i);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddBMembrane(
        bounded_matrix<double, 3, 18 > & BMembrane,
        bounded_matrix<double, 3, 1  > & CMembrane,
        const bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
        const bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
        const IndexType NodeGauss
        )
{
    for (IndexType i = 0; i < 4; i++) {
        IndexType base = i * 3;
        if (i == 3)
            base += NodeGauss * 3;

        for (IndexType j = 0; j < 3; j++) {
            BMembrane(0, base + j) += InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 0);
            BMembrane(1, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 1);
            BMembrane(2, base + j) += InPlaneCartesianDerivativesGauss(1, i) * InPlaneGradientFGauss(j, 0)
                                   +  InPlaneCartesianDerivativesGauss(0, i) * InPlaneGradientFGauss(j, 1);
        }
    }

    /* Calculate de componets of Cauchy tensor */
    // In plane auxiliar components
    array_1d<double, 3 > aux_deformation_gradient_F1, aux_deformation_gradient_F2;

    for (IndexType i = 0; i < 3; i++) {
        aux_deformation_gradient_F1[i] = InPlaneGradientFGauss(i, 0);
        aux_deformation_gradient_F2[i] = InPlaneGradientFGauss(i, 1);
    }

    CMembrane(0, 0) += inner_prod(aux_deformation_gradient_F1, aux_deformation_gradient_F1);
    CMembrane(1, 0) += inner_prod(aux_deformation_gradient_F2, aux_deformation_gradient_F2);
    CMembrane(2, 0) += inner_prod(aux_deformation_gradient_F1, aux_deformation_gradient_F2);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddMembraneKgeometric(
    bounded_matrix<double, 36, 36 > & Kgeometricmembrane,
    const CartesianDerivatives& rCartesianDerivatives,
    const array_1d<double, 3 > & SMembrane,
    const IndexType Index
    )
{
    const IndexType auxiliar_index = Index == 9 ? 3 : 0;

    bounded_matrix<double, 6, 6 > H = ZeroMatrix(6, 6);

    IndexType ii;
    IndexType jj;
    for (IndexType i = 0; i < 4; i++) {
        for (IndexType j = 0; j < 4; j++) {
            // Gauss 1
            ii = i;
            jj = j;
            H(ii, jj) += SMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, j)
                       + SMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, j)
                       + SMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, j)
                                       + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 0](0, j));

            // Gauss 2
            if (i ==  3)
                ii = 4;
            else
                ii = i;
            if (j ==  3)
                jj = 4;
            else
                jj = j;

            H(ii, jj) += SMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, j)
                       + SMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, j)
                       + SMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, j)
                                       + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 1](0, j));

            // Gauss 3
            if (i ==  3)
                ii = 5;
            else
                ii = i;
            if (j ==  3)
                jj = 5;
            else
                jj = j;

            H(ii, jj) += SMembrane[0] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, j)
                       + SMembrane[1] *  rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, j)
                       + SMembrane[2] * (rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, j)
                                       + rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](1, i) * rCartesianDerivatives.InPlaneCartesianDerivativesGauss[auxiliar_index + 2](0, j));
        }
    }

    H *= 1.0/3.0;

    // Assembling in Kgeometricmembrane
    IndexType rowindex, colindex;
    for (IndexType i = 0; i < 6; i++) {
        if (i < 3)
            rowindex = i * 3 + Index;
        else
            rowindex = i * 3 + Index + 9;

        for (IndexType j = i; j < 6; j++) {
            if (j < 3)
                colindex = j * 3 + Index;
            else
                colindex = j * 3 + Index + 9;

            for(IndexType ii = 0; ii < 3; ii++) {
                Kgeometricmembrane(rowindex + ii,colindex + ii) += H (i, j);
                if (rowindex != colindex) { // Skip diagonal
                    Kgeometricmembrane(colindex + ii, rowindex + ii) += H (i, j); // Symmetric part
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddBShear(
    bounded_matrix<double, 2, 18 >& BShear,
    bounded_matrix<double, 2, 1 >& CShear,
    const CartesianDerivatives& rCartesianDerivatives,
    const TransverseGradient& rTransverseGradient,
    const TransverseGradientIsoParametric& rTransverseGradientIsoParametric,
    const IndexType Index
    )
{
    const IndexType auxiliar_index = Index == 9 ? 3 : 0;

    const bounded_matrix<double, 2, 2 >& JInvPlane = Index == 9 ? rCartesianDerivatives.JInvPlaneUpper : rCartesianDerivatives.JInvPlaneLower;

    // Considering the Gauss point in the middle of the element
    const double eta_p = 1.0/3.0;
    const double xi_p  = 1.0/3.0;
    bounded_matrix<double, 2, 3 > Pa;
    Pa(0, 0) = - xi_p;
    Pa(0, 1) = - xi_p;
    Pa(0, 2) = 1.0 - xi_p;
    Pa(1, 0) = eta_p;
    Pa(1, 1) = eta_p - 1.0;
    Pa(1, 2) = eta_p;

    bounded_matrix<double, 3, 18 > aux_b_shear = ZeroMatrix(3, 18);

    /* First contribution*/
    for (IndexType i = 0; i < 6; i++) {
        IndexType base = i * 3;
        for (IndexType j = 0; j < 3; j++) {
            aux_b_shear(0, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliar_index + 0](i, 0) * rTransverseGradientIsoParametric.Ft[j];
            aux_b_shear(1, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliar_index + 1](i, 0) * rTransverseGradientIsoParametric.Fxi[j];
            aux_b_shear(2, base + j) += rCartesianDerivatives.TransversalCartesianDerivativesGauss[auxiliar_index + 2](i, 0) * rTransverseGradientIsoParametric.Feta[j];
        }
    }

    /* Second contibution */
    for (IndexType i = 0; i < 3; i++) {
        /* First row */
        aux_b_shear(0, i + Index + 3) -= rTransverseGradient.F0[i];
        aux_b_shear(0, i + Index + 6) += rTransverseGradient.F0[i];

        /* Second row */
        aux_b_shear(1, i + Index)     += rTransverseGradient.F1[i];
        aux_b_shear(1, i + Index + 6) -= rTransverseGradient.F1[i];

        /* Third row */
        aux_b_shear(2, i + Index)     -= rTransverseGradient.F2[i];
        aux_b_shear(2, i + Index + 3) += rTransverseGradient.F2[i];
    }

    const bounded_matrix<double, 2, 3 > aux_prod = prod(JInvPlane, Pa);
    noalias(BShear) = prod(aux_prod, aux_b_shear);

    // Calculating the components of C
    bounded_matrix<double, 3, 1 > aux_c_shear;
    aux_c_shear(0, 0) = inner_prod(rTransverseGradientIsoParametric.Ft  , rTransverseGradient.F0);
    aux_c_shear(1, 0) = inner_prod(rTransverseGradientIsoParametric.Fxi , rTransverseGradient.F1);
    aux_c_shear(2, 0) = inner_prod(rTransverseGradientIsoParametric.Feta, rTransverseGradient.F2);

    noalias(CShear) = prod(aux_prod, aux_c_shear);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddShearKgeometric(
    bounded_matrix<double, 18, 18 > & Kgeometricshear,
    const CartesianDerivatives& rCartesianDerivatives,
    const array_1d<double, 2 > & SShear,
    const IndexType Index
    )
{
    const IndexType auxiliar_index = Index == 9 ? 3 : 0;

    const bounded_matrix<double, 2, 2 >& JInvPlane = Index == 9 ? rCartesianDerivatives.JInvPlaneUpper : rCartesianDerivatives.JInvPlaneLower;

    const double Q1 = 1.0/3.0 * (SShear[0] * JInvPlane(0, 0) + SShear[1] * JInvPlane(0, 1));
    const double Q2 = 1.0/3.0 * (SShear[0] * JInvPlane(1, 0) + SShear[1] * JInvPlane(1, 1));

//    array_1d<double, 3 > q;
//    q[0] = -Q1 + Q2;
//    q[1] = -(Q1 + 2.0 * Q2);
//    q[2] = (2.0 * Q1 + Q2);

//    int delta;
//    if (index == 9)
//        delta = 3;
//    else
//        delta = 0;

//    for (IndexType i = 0; i < 3; i++) { // For each DOF
//        /* First assembling */
//        Kgeometricshear(i + index + 3, i + index + 3) -= q[0] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliar_index](1 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) += q[0] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliar_index](2 + delta, 0);

//        /* Second assembling */
//        Kgeometricshear(i + index, i + index)         += q[1] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliar_index](0 + delta, 0);
//        Kgeometricshear(i + index + 6, i + index + 6) -= q[1] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliar_index](2 + delta, 0);
//        /* Third assembling */
//        Kgeometricshear(i + index, i + index)         -= q[2] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliar_index](0 + delta, 0);
//        Kgeometricshear(i + index + 3, i + index + 3) += q[2] * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliar_index](1 + delta, 0);
//    }

    array_1d<double, 3 > n1; // Side node with + contribution (previous DOF position)
    array_1d<double, 3 > n2; // Side node with - contribution (previous DOF position)

//    if (index == 0) {
//        n1[0] = 6;
//        n1[1] = 3;
//        n1[2] = 0;

//        n2[0] = 6;
//        n2[1] = 3;
//        n2[2] = 0;
//    } else {
//        n1[0] = 15;
//        n1[1] = 12;
//        n1[2] = 9;

//        n2[0] = 15;
//        n2[1] = 12;
//        n2[2] = 9;
//    }

    // Note: Technically this is the correct one
    if (Index == 0) {
        n1[0] = 6;
        n1[1] = 0;
        n1[2] = 3;

        n2[0] = 3;
        n2[1] = 6;
        n2[2] = 0;
    } else {
        n1[0] = 15;
        n1[1] = 9;
        n1[2] = 12;

        n2[0] = 12;
        n2[1] = 15;
        n2[2] = 9;
    }

    double value = 0.0;
    for (IndexType k = 0; k < 3; k++) {
        IndexType l = 0; // Initializes DOF associated to N_3
        for (IndexType i = 0; i < 6; i++) { //  For each node
            if (k == 0)
                value = (-Q1 + Q2) *  rCartesianDerivatives.TransversalCartesianDerivativesGauss[0 + auxiliar_index](i, 0);
            else if (k == 1)
                value = -(Q1 + 2.0 * Q2) * rCartesianDerivatives.TransversalCartesianDerivativesGauss[1 + auxiliar_index](i, 0);
            else if (k == 2)
                value = (2.0 * Q1 + Q2) * rCartesianDerivatives.TransversalCartesianDerivativesGauss[2 + auxiliar_index](i, 0);

            for (IndexType j = 0; j < 3; j++) { // For each DOF (diagonal only)
                Kgeometricshear(n1[k] + j, l + j) += value;
                Kgeometricshear(l + j, n1[k] + j) += value;
            }

            for (IndexType j = 0; j < 3; j++) { // For each DOF (diagonal only)
                Kgeometricshear(n2[k] + j, l + j) -= value;
                Kgeometricshear(l + j, n2[k] + j) -= value;
            }

            l += 3; // Increment DOF position I
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddBNormal(
        bounded_matrix<double, 1, 18 > & BNormal,
        double & CNormal,
        const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesCenter,
        const array_1d<double, 3 > & TransversalDeformationGradientF
        )
{
    for (IndexType i = 0; i < 6; i++) {
        IndexType base = i * 3;
        for (IndexType j = 0; j < 3; j++)
            BNormal(0, base + j) = TransversalCartesianDerivativesCenter(i, 0) * TransversalDeformationGradientF[j];
    }

    CNormal = inner_prod(TransversalDeformationGradientF, TransversalDeformationGradientF);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddNormalKgeometric(
    bounded_matrix<double, 18, 18 > & Kgeometricnormal,
    const bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesCenter,
    const double SNormal
    )
{
    bounded_matrix<double, 6, 6 > H = ZeroMatrix(6, 6);
    for (IndexType i = 0; i < 6; i++) {
        const double aux = SNormal * TransversalCartesianDerivativesCenter(i, 0);
        for (IndexType j = 0; j < 6; j++)
            H(i, j) =  aux * TransversalCartesianDerivativesCenter(j, 0);
    }

    noalias(H) = SNormal * prod(TransversalCartesianDerivativesCenter, trans(TransversalCartesianDerivativesCenter));

    IndexType rowindex;
    IndexType colindex;
    for (IndexType i = 0; i < 6; i++) {
        rowindex = i * 3;
        for (IndexType j = 0; j < 6; j++) {
            colindex = j * 3;
            for(IndexType ii = 0; ii < 3; ii++)
                Kgeometricnormal(rowindex + ii,colindex + ii) += H(i, j);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bounded_matrix<double, 36, 1 > SprismElement3D6N::GetVectorCurrentPosition()
{
    KRATOS_TRY;

    bounded_matrix<double, 36, 1 > vector_current_position;

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (IndexType index = 0; index < 6; index++) {
        const array_1d<double,3>& current_position = GetGeometry()[index].Coordinates();
        for (IndexType j = 0; j < 3; j++)
            vector_current_position(index * 3 + j, 0) = current_position[j];
    }

    /* Neighbour nodes */
    const SizeType number_of_neighbours = NumberOfActiveNeighbours(p_neighbour_nodes);

    if (number_of_neighbours == 6) { // All the possible neighours
        for (IndexType index = 0; index < 6; index++) {
            const array_1d<double,3>& current_position = p_neighbour_nodes[index].Coordinates();
            for (IndexType j = 0; j < 3; j++)
                vector_current_position(18 + index * 3 + j, 0) = current_position[j];
        }
    } else {
        for (IndexType index = 0; index < 6; index++) {
            if (HasNeighbour(index, p_neighbour_nodes[index])) {
                const array_1d<double,3>& current_position = p_neighbour_nodes[index].Coordinates();
                for (IndexType j = 0; j < 3; j++)
                    vector_current_position(18 + index * 3 + j, 0) = current_position[j];
            } else {
                for (IndexType j = 0; j < 3; j++)
                    vector_current_position(18 + index * 3 + j, 0) = 0.0;
            }
        }
    }

    return vector_current_position;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bounded_matrix<double, 36, 1 > SprismElement3D6N::GetVectorPreviousPosition()
{
    KRATOS_TRY;

    bounded_matrix<double, 36, 1 > vector_current_position;

    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);

    /* Element nodes */
    for (IndexType index = 0; index < 6; index++) {
        const array_1d<double,3>& previous_position = GetGeometry()[index].GetInitialPosition().Coordinates()
                                                    + GetGeometry()[index].FastGetSolutionStepValue(DISPLACEMENT, 1);
        for (IndexType j = 0; j < 3; j++)
            vector_current_position(index * 3 + j, 0) = previous_position[j];
    }

    /* Neighbour nodes */
    const SizeType number_of_neighbours = NumberOfActiveNeighbours(p_neighbour_nodes);

    if (number_of_neighbours == 6) { // All the possible neighours
        for (IndexType index = 0; index < 6; index++) {
            const array_1d<double,3>& previous_position = p_neighbour_nodes[index].GetInitialPosition().Coordinates()
                                                 + p_neighbour_nodes[index].FastGetSolutionStepValue(DISPLACEMENT, 1);

            for (IndexType j = 0; j < 3; j++)
                vector_current_position(18 + index * 3 + j, 0) = previous_position[j];
        }
    } else {
        for (IndexType index = 0; index < 6; index++) {
            if (HasNeighbour(index, p_neighbour_nodes[index])) {
                const array_1d<double,3>& previous_position = p_neighbour_nodes[index].GetInitialPosition().Coordinates()
                                                     + p_neighbour_nodes[index].FastGetSolutionStepValue(DISPLACEMENT, 1);
                for (IndexType j = 0; j < 3; j++)
                    vector_current_position(18 + index * 3 + j, 0) = previous_position[j];
            } else {
                for (IndexType j = 0; j < 3; j++)
                    vector_current_position(18 + index * 3 + j, 0) = 0.0;
            }
        }
    }

    return vector_current_position;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::IntegrateInZeta(
        GeneralVariables& rVariables,
        StressIntegratedComponents& rIntegratedStress,
        const double AlphaEAS,
        const double ZetaGauss,
        const double IntegrationWeight
        )
{
    KRATOS_TRY;
    
    const double L1 = 0.5 * (1.0 - ZetaGauss);
    const double L2 = 0.5 * (1.0 + ZetaGauss);

    const double factor_eas = std::exp(2.0 * AlphaEAS * ZetaGauss);

    /* INTEGRATE PK2 IN ZETA */
    // Integrate stresses in the reference configuration
    /* In plane stresses */
    // Lower
    rIntegratedStress.SMembraneLower(0) +=  L1 * IntegrationWeight * rVariables.StressVector[0]; // xx
    rIntegratedStress.SMembraneLower(1) +=  L1 * IntegrationWeight * rVariables.StressVector[1]; // yy
    rIntegratedStress.SMembraneLower(2) +=  L1 * IntegrationWeight * rVariables.StressVector[3]; // xy
    // Upper
    rIntegratedStress.SMembraneUpper(0) +=  L2 * IntegrationWeight * rVariables.StressVector[0]; // xx
    rIntegratedStress.SMembraneUpper(1) +=  L2 * IntegrationWeight * rVariables.StressVector[1]; // yy
    rIntegratedStress.SMembraneUpper(2) +=  L2 * IntegrationWeight * rVariables.StressVector[3]; // xy

    /* Transversal stresses */ // Note: Order according to the Voigt Notation in the Wiki
    // Lower face
    rIntegratedStress.SShearLower(0)    +=  L1 * IntegrationWeight * rVariables.StressVector[5]; // xz
    rIntegratedStress.SShearLower(1)    +=  L1 * IntegrationWeight * rVariables.StressVector[4]; // yz
    // Upper face
    rIntegratedStress.SShearUpper(0)    +=  L2 * IntegrationWeight * rVariables.StressVector[5]; // xz
    rIntegratedStress.SShearUpper(1)    +=  L2 * IntegrationWeight * rVariables.StressVector[4]; // yz

    /* Normal stress */
    rIntegratedStress.SNormal           +=  factor_eas * IntegrationWeight * rVariables.StressVector[2]; // zz

    /* INTEGRATE EAS IN ZETA */
    // Calculate EAS residual
    mEAS.mRHSAlpha += IntegrationWeight * ZetaGauss * rVariables.StressVector[2] * rVariables.C[2];

    // Calculate EAS stiffness
    mEAS.mStiffAlpha += IntegrationWeight * ZetaGauss * ZetaGauss * rVariables.C[2]
            * (rVariables.ConstitutiveMatrix(2, 2) * rVariables.C[2] + 2.0 * rVariables.StressVector[2]);

    bounded_matrix<double, 1, 36 > B3;
    bounded_matrix<double, 1,  6 > D3;

    for (IndexType i = 0; i < 6; i++)
        D3(0, i) = rVariables.ConstitutiveMatrix(2, i);
    for (IndexType i = 0; i < 36; i++)
        B3(0, i) = rVariables.B(2, i);

    // Calculate H operator
    noalias(mEAS.mHEAS) += IntegrationWeight * ZetaGauss
            * (rVariables.C[2] * prod(D3, rVariables.B) + 2.0 * rVariables.StressVector[2] * B3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddLHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const StressIntegratedComponents& rIntegratedStress,
        const CommonComponents& rCommonComponents,
        const CartesianDerivatives& rCartesianDerivatives,
        double& AlphaEAS
        )
{
    /* Contributions of the stiffness matrix calculated on the reference configuration */
    if( rLocalSystem.CalculationFlags.Is( SprismElement3D6N::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) ) {
        std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
        const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

        for( IndexType i = 0; i < rLeftHandSideVariables.size(); i++ ) {
            bool calculated = false;
            /* Calculate the Material Stiffness Matrix */
            if( rLeftHandSideVariables[i] == MATERIAL_STIFFNESS_MATRIX ) {
                /* Reading integration points */
                const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

                for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
                    const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

                    /* Assemble B */
                    this->CalculateDeformationMatrix(rVariables.B, rCommonComponents, zeta_gauss, AlphaEAS);

                    // Compute element kinematics C, F ...
                    this->CalculateKinematics(rVariables, rCommonComponents, point_number, AlphaEAS, zeta_gauss);

                    // Set general variables to constitutivelaw parameters
                    this->SetGeneralVariables(rVariables, rValues, point_number);

                    // Compute stresses and constitutive parameters
                    mConstitutiveLawVector[point_number]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

                    // Calculating weights for integration on the "reference configuration"
                    const double integration_weight = integration_points[point_number].Weight() * rVariables.detJ;

                    /* Operation performed: add Km to the LefsHandSideMatrix */
                    this->CalculateAndAddKuum( rLeftHandSideMatrices[i], rVariables, integration_weight);
                }
                calculated = true;
            }

            /* Calculate the Geometric Stiffness Matrix */
            if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX ) {
                /* Operation performed: add Kg to the LefsHandSideMatrix */
                this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rIntegratedStress, rCartesianDerivatives );
                calculated = true;
            }

            /* Implicit or explicit EAS update*/
            if ( mELementalFlags.Is(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT)) {
                /* Apply EAS stabilization */
                ApplyEASLHS(rLeftHandSideMatrices[i]);
            }

            KRATOS_ERROR_IF_NOT(calculated) << " ELEMENT can not supply the required local system variable: " << rLeftHandSideVariables[i] << std::endl;
        }
    } else {
        MatrixType& LeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

        /* Reading integration points */
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        /* Calculate the Material Stiffness Matrix */
        for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
            const double zeta_gauss = 2.0 * integration_points[point_number].Z() - 1.0;

            /* Assemble B */
            this->CalculateDeformationMatrix(rVariables.B, rCommonComponents, zeta_gauss, AlphaEAS);

            // Compute element kinematics C, F ...
            this->CalculateKinematics(rVariables, rCommonComponents, point_number, AlphaEAS, zeta_gauss);

            // Set general variables to constitutivelaw parameters
            this->SetGeneralVariables(rVariables, rValues, point_number);

            // Compute stresses and constitutive parameters
            mConstitutiveLawVector[point_number]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

            // Calculating weights for integration on the "reference configuration"
            const double integration_weight = integration_points[point_number].Weight() * rVariables.detJ;

            /* Operation performed: add Km to the LefsHandSideMatrix */
            this->CalculateAndAddKuum( LeftHandSideMatrix, rVariables, integration_weight);
        }

        /* Calculate the Geometric Stiffness Matrix */
        /* Operation performed: add Kg to the LefsHandSideMatrix */
        this->CalculateAndAddKuug( LeftHandSideMatrix, rIntegratedStress, rCartesianDerivatives );

        /* Implicit or explicit EAS update*/
        if ( mELementalFlags.Is(SprismElement3D6N::EAS_IMPLICIT_EXPLICIT)) {
            /* Apply EAS stabilization */
            ApplyEASLHS(LeftHandSideMatrix);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddRHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const StressIntegratedComponents& rIntegratedStress,
        const CommonComponents& rCommonComponents,
        double& AlphaEAS
        )
{
    /* Contribution of the internal and external forces */
    if( rLocalSystem.CalculationFlags.Is( SprismElement3D6N::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) ) {
        std::vector<VectorType>& RightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
        const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
        for( IndexType i = 0; i < rRightHandSideVariables.size(); i++ ) {
            bool calculated = false;
            if( rRightHandSideVariables[i] == EXTERNAL_FORCES_VECTOR ) {
                /* Operation performed: RightHandSideVector += ExtForce */
                this->CalculateAndAddExternalForces( RightHandSideVectors[i], rVariables, rVolumeForce );
                calculated = true;
            }

            if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR ) {
                /* Operation performed: RightHandSideVector -= IntForce */
                this->CalculateAndAddInternalForces( RightHandSideVectors[i], rIntegratedStress, rCommonComponents, AlphaEAS );
                calculated = true;
            }

            KRATOS_ERROR_IF_NOT(calculated) << " ELEMENT can not supply the required local system variable: " << rRightHandSideVariables[i] << std::endl;
        }
    } else {
        VectorType& RightHandSideVector = rLocalSystem.GetRightHandSideVector();

        /* Operation performed: RightHandSideVector += ExtForce */
        this->CalculateAndAddExternalForces( RightHandSideVector, rVariables, rVolumeForce );

        /* Operation performed: RightHandSideVector -= IntForce */
        this->CalculateAndAddInternalForces( RightHandSideVector, rIntegratedStress, rCommonComponents, AlphaEAS );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddKuum(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;
    
    /* Calculate K */
    typedef bounded_matrix<double,  6, 36 > temp_type;
    const bounded_matrix<double, 36, 36 > K = IntegrationWeight * prod(trans(rVariables.B), prod<temp_type>(rVariables.ConstitutiveMatrix, rVariables.B));

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    IndexType index_i, index_j;

    for (IndexType i = 0; i < 36; i++) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (id_vector[i] < 36)
            for (IndexType j = 0; j < 36; j++) {
                index_j = j < 18 ? j : id_vector[j - 18];
                if (id_vector[j] < 36)
                    rLeftHandSideMatrix(index_i, index_j) += K(i, j);
            }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddKuug(
        MatrixType& rLeftHandSideMatrix,
        const StressIntegratedComponents& rIntegratedStress,
        const CartesianDerivatives& rCartesianDerivatives
        )
{
    KRATOS_TRY;

    /* The stress is already integrated, we just calculate it once */

    /* Auxiliar stiffness matrix */
    bounded_matrix<double, 18, 18 > aux_K = ZeroMatrix(18, 18); // Auxiliar stiffness matrix
    bounded_matrix<double, 36, 36 >     K = ZeroMatrix(36, 36); // Stiffness matrix

    /* COMPUTATION OF GEOMETRIC STIFFNESS MATRIX */

    /* MEMBRANE CONTRIBUTION */
    /* Adding the geometric membrane stiffness */
    // Lower face
    CalculateAndAddMembraneKgeometric(K, rCartesianDerivatives, rIntegratedStress.SMembraneLower, 0);
    // Upper face
    CalculateAndAddMembraneKgeometric(K, rCartesianDerivatives, rIntegratedStress.SMembraneUpper, 9);

//    /* SHEAR CONTRIBUTION */
//    /* Adding the geometric shear stiffness */
//    // Lower face
//    CalculateAndAddShearKgeometric(aux_K, rCartesianDerivatives, rIntegratedStress.SShearLower, 0);
//    // Upper face
//    CalculateAndAddShearKgeometric(aux_K, rCartesianDerivatives, rIntegratedStress.SShearUpper, 9);

    /* NORMAL TRANSVERSE */
    /* Adding the geometric normal stiffness */
    CalculateAndAddNormalKgeometric(aux_K, rCartesianDerivatives.TransversalCartesianDerivativesCenter, rIntegratedStress.SNormal);

    // Transfering to the complete stiffness matrix
    for (IndexType i = 0; i < 18; i++)
        for (IndexType j = 0; j < 18; j++)
            K(i, j) += aux_K(i, j);

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    IndexType index_i, index_j;

    for (IndexType i = 0; i < 36; i++) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (id_vector[i] < 36)
            for (IndexType j = 0; j < 36; j++) {
                index_j = j < 18 ? j : id_vector[j - 18];
                if (id_vector[j] < 36)
                    rLeftHandSideMatrix(index_i, index_j) += K(i, j);
            }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ApplyEASLHS(MatrixType& rLeftHandSideMatrix)
{
    KRATOS_TRY;

    const bounded_matrix<double, 36, 36 > lhs_aux = - prod(trans(mEAS.mHEAS), mEAS.mHEAS) / mEAS.mStiffAlpha;

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    // Note: IntegrationWeight already considered in the integration
    IndexType index_i, index_j;

    for (IndexType i = 0; i < 36; i++) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (id_vector[i] < 36)
            for (IndexType j = 0; j < 36; j++) {
                index_j = j < 18 ? j : id_vector[j - 18];
                if (id_vector[j] < 36)
                    rLeftHandSideMatrix(index_i, index_j) += lhs_aux(i, j);
            }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ApplyEASRHS(
        bounded_matrix<double, 36, 1 >& rRHSFull,
        double& AlphaEAS
        )
{
    KRATOS_TRY;

    /* Calculate the RHS */
    noalias(rRHSFull) -= trans(mEAS.mHEAS) * mEAS.mRHSAlpha / mEAS.mStiffAlpha;

    /* Update ALPHA_EAS */
    AlphaEAS -= mEAS.mRHSAlpha / mEAS.mStiffAlpha;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddExternalForces(
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce
        )
{
    KRATOS_TRY;

    const IndexType number_of_nodes = GetGeometry().PointsNumber();

    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const IndexType index = 3 * i;
        for ( IndexType j = 0; j < 3; j++ ) {
            rRightHandSideVector[index + j] += rVolumeForce[j]/static_cast<double>(number_of_nodes);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateAndAddInternalForces(
        VectorType& rRightHandSideVector,
        const StressIntegratedComponents& rIntegratedStress,
        const CommonComponents& rCommonComponents,
        double& AlphaEAS
        )
{
    KRATOS_TRY;

    bounded_matrix<double, 36, 1 > rhs_full = ZeroMatrix(36, 1);

    IndexType aux_index = 0;
    for (IndexType i = 0; i < 18; i++) {
        if (i == 9)
            aux_index += 9;

        /* Calculate residual forces */
        /* Apply membrane stress, adding the in-plane nodal force contribution */
        /* Nodes 1-3  and 7-9 */
        rhs_full(aux_index + i, 0)     += rIntegratedStress.SMembraneLower[0] * rCommonComponents.BMembraneLower(0, i); // xx
        rhs_full(aux_index + i, 0)     += rIntegratedStress.SMembraneLower[1] * rCommonComponents.BMembraneLower(1, i); // yy
        rhs_full(aux_index + i, 0)     += rIntegratedStress.SMembraneLower[2] * rCommonComponents.BMembraneLower(2, i); // xy

        /* Nodes 4-6  and 10-12 */
        rhs_full(aux_index + i + 9, 0) += rIntegratedStress.SMembraneUpper[0] * rCommonComponents.BMembraneUpper(0, i); // xx
        rhs_full(aux_index + i + 9, 0) += rIntegratedStress.SMembraneUpper[1] * rCommonComponents.BMembraneUpper(1, i); // yy
        rhs_full(aux_index + i + 9, 0) += rIntegratedStress.SMembraneUpper[2] * rCommonComponents.BMembraneUpper(2, i); // xy

        /* Apply transversal forces */
        /* Apply shear stress, adding the transverse nodal force contribution */
        rhs_full(i, 0) += rIntegratedStress.SShearLower[0] * rCommonComponents.BShearLower(0, i); // xz
        rhs_full(i, 0) += rIntegratedStress.SShearLower[1] * rCommonComponents.BShearLower(1, i); // yz
        rhs_full(i, 0) += rIntegratedStress.SShearUpper[0] * rCommonComponents.BShearUpper(0, i); // xz
        rhs_full(i, 0) += rIntegratedStress.SShearUpper[1] * rCommonComponents.BShearUpper(1, i); // yz

        /* Apply normal transverse stress */
        rhs_full(i, 0) += rIntegratedStress.SNormal * rCommonComponents.BNormal(0, i); // zz
    }

    /* Apply EAS stabilization */
    ApplyEASRHS(rhs_full, AlphaEAS);

    // Compute vector of IDs
    array_1d<IndexType, 18> id_vector;
    CalculateIdVector(id_vector);

    IndexType index_i;

    for (IndexType i = 0; i < 36; i++) {
        index_i = i < 18 ? i : id_vector[i - 18];
        if (id_vector[i] < 36)
            rRightHandSideVector[index_i] -= rhs_full(i, 0);
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::SetGeneralVariables(
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType rPointNumber
        )
{
    if(rVariables.detF < 0) {
        std::cout<<" Element: "<<this->Id()<<std::endl;
        IndexType number_of_nodes = GetGeometry().PointsNumber();

        for (IndexType i = 0; i < number_of_nodes; i++) {
            array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
            std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
            std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
        }

        KRATOS_WATCH(rVariables.F);
        KRATOS_ERROR << " SPRISM ELEMENT INVERTED: |F| < 0  detF = " << rVariables.detF << std::endl;
    }

    // Compute total F: FT
    rVariables.detFT = rVariables.detF * rVariables.detF0;
    rVariables.FT    = prod( rVariables.F, rVariables.F0 );

    rValues.SetDeterminantF(rVariables.detFT);
    rValues.SetDeformationGradientF(rVariables.FT);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);

    // Adding the standard prism shape functions
    rValues.SetShapeFunctionsValues(rVariables.N);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeSystemMatrices(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags
        )
{
    // Resizing as needed the LHS
    WeakPointerVector< NodeType >& p_neighbour_nodes = this->GetValue(NEIGHBOUR_NODES);
    const IndexType number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(p_neighbour_nodes);
    const IndexType mat_size = number_of_nodes * 3;

    if ( rCalculationFlags.Is(SprismElement3D6N::COMPUTE_LHS_MATRIX) ) {// Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); // Resetting LHS
    }

    // Resizing as needed the RHS
    if ( rCalculationFlags.Is(SprismElement3D6N::COMPUTE_RHS_VECTOR) ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); // Resetting RHS
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeMaterial()
{
    KRATOS_TRY;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL ) {
        for ( IndexType i = 0; i < mConstitutiveLawVector.size(); i++ ) {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::ResetConstitutiveLaw()
{
    KRATOS_TRY;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL ) {
        for ( IndexType i = 0; i < mConstitutiveLawVector.size(); i++ ) {
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }

    KRATOS_CATCH( "" );
}

/******************************* COMPUTE KINEMATICS ********************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateKinematics(
    GeneralVariables& rVariables,
    const CommonComponents& rCommonComponents,
    const IndexType rPointNumber,
    const double AlphaEAS,
    const double ZetaGauss
    )
{
    KRATOS_TRY;

    const double L_1 = 0.5 * (1.0 - ZetaGauss);
    const double L_2 = 0.5 * (1.0 + ZetaGauss);

    const double factor_eas = std::exp(2.0 * AlphaEAS * ZetaGauss);  // EAS factor

    /* Assemble C */
    rVariables.C[0] = L_1 * rCommonComponents.CMembraneLower(0, 0) + L_2 * rCommonComponents.CMembraneUpper(0, 0); // xx
    rVariables.C[1] = L_1 * rCommonComponents.CMembraneLower(1, 0) + L_2 * rCommonComponents.CMembraneUpper(1, 0); // yy
    rVariables.C[2] = factor_eas * rCommonComponents.CNormal;                                            // zz
    rVariables.C[3] = L_1 * rCommonComponents.CMembraneLower(2, 0) + L_2 * rCommonComponents.CMembraneUpper(2, 0); // xy
    rVariables.C[4] = L_1 * rCommonComponents.CShearLower(1, 0)    + L_2 * rCommonComponents.CShearUpper(1, 0);    // yz
    rVariables.C[5] = L_1 * rCommonComponents.CShearLower(0, 0)    + L_2 * rCommonComponents.CShearUpper(0, 0);    // xz

    rVariables.detF = rVariables.C[0] * rVariables.C[1] * rVariables.C[2] + 2 * rVariables.C[3] * rVariables.C[4] * rVariables.C[5]
                    - rVariables.C[5] * rVariables.C[5] * rVariables.C[1] -     rVariables.C[4] * rVariables.C[4] * rVariables.C[0]
                    - rVariables.C[3] * rVariables.C[3] * rVariables.C[2];


    KRATOS_ERROR_IF(rVariables.detF < std::numeric_limits<double>::epsilon()) << "The determinant of C is zero or negative.  det(C): " << rVariables.detF << std::endl;

    rVariables.detF = std::sqrt(rVariables.detF);

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // PK2 stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;

        // Jacobian Determinant for the isoparametric and numerical integration
        rVariables.detJ = mAuxCont[rPointNumber];
    } else {
        // Cauchy stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

        //Determinant of the Deformation Gradient F0
        rVariables.detF0 = mAuxCont[rPointNumber];
        rVariables.F0    = mAuxMatCont[rPointNumber];
    }

    this->CbartoFbar(rVariables, rPointNumber);

    // Get the shape functions for the order of the integration method [N]
    const Matrix& N_container = rVariables.GetShapeFunctions();

    // Set Shape Functions Values for this integration point
    rVariables.N = row( N_container, rPointNumber);

    KRATOS_CATCH( "" );
}

/***************************** COMPUTE DELTA POSITION ******************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY;

    rDeltaPosition = ZeroMatrix( 6 , 3);

    for ( IndexType i = 0; i < 6; i++ ) {
        const array_1d<double, 3 > & current_displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 > & previous_displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, 1);

        for ( IndexType j = 0; j < 3; j++ )
            rDeltaPosition(i,j) = current_displacement[j] - previous_displacement[j];
    }

    KRATOS_CATCH( "" );
}


/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CbartoFbar(
        GeneralVariables& rVariables,
        const int rPointNumber
        )
{
    KRATOS_TRY;

    /* We perform a polar decomposition of the CBar and F(regular) to obtain F_bar */

    /* Decompose C_bar */
    bounded_matrix<double, 3, 3> eigen_vector_matrix,  eigen_values_matrix;

    // Assemble matrix C_bar
    const Matrix C_bar = MathUtils<double>::VectorToSymmetricTensor(rVariables.C);

    // Decompose matrix C_bar
    MathUtils<double>::EigenSystem<3>(C_bar, eigen_vector_matrix, eigen_values_matrix, 1e-24, 100);

    for (IndexType i = 0; i < 3; i++)
        eigen_values_matrix(i, i) = std::sqrt(eigen_values_matrix(i, i));

    const Matrix U_bar = prod( eigen_values_matrix, eigen_vector_matrix );

    /* Decompose F */
    Matrix F = ZeroMatrix(3, 3);
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN)) {
        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[rPointNumber], mAuxMatCont[rPointNumber] );
    } else {
        // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
        Matrix InvJ(3, 3);
        MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

        // Deformation Gradient F [dx_n+1/dx_n]
        noalias(F) = prod( rVariables.j[rPointNumber], InvJ );
    }

    const Matrix C = prod( trans(F), F );

    // Decompose matrix C
    MathUtils<double>::EigenSystem<3>(C, eigen_vector_matrix, eigen_values_matrix, 1e-24, 100);

    for (IndexType i = 0; i < 3; i++)
        eigen_values_matrix(i, i) = std::sqrt(eigen_values_matrix(i, i));

    const Matrix U  = prod( eigen_values_matrix, eigen_vector_matrix );

    double AuxDet;
    Matrix invU(3, 3);
    MathUtils<double>::InvertMatrix(U, invU, AuxDet);
    const Matrix R  = prod( F, invU );

    /* Calculate F_bar */
    noalias(rVariables.F) = prod(R, U_bar);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateDeformationMatrix(
        Matrix& rB,
        const CommonComponents& rCommonComponents,
        const double ZetaGauss,
        const double AlphaEAS
        )
{
    KRATOS_TRY;

    rB.clear(); // Set all components to zero

    const double L1 = 0.5 * (1.0 - ZetaGauss);
    const double L2 = 0.5 * (1.0 + ZetaGauss);

    const double factor_eas = std::exp(2.0 * AlphaEAS * ZetaGauss); // EAS factor

    for (IndexType index = 0; index < 9; index++) {
        /* Element nodes */ // Note: It's important to consider the Voigt notation order considered in Kratos
        // Lower face
        rB(0, index)      = L1 * rCommonComponents.BMembraneLower(0, index);  // xx
        rB(1, index)      = L1 * rCommonComponents.BMembraneLower(1, index);  // yy
        rB(2, index)      = factor_eas * rCommonComponents.BNormal(0, index);     // zz
        rB(3, index)      = L1 * rCommonComponents.BMembraneLower(2, index);  // xy
        rB(4, index)      = L1 * rCommonComponents.BShearLower(1, index) + L2 * rCommonComponents.BShearUpper(1, index); // yz
        rB(5, index)      = L1 * rCommonComponents.BShearLower(0, index) + L2 * rCommonComponents.BShearUpper(0, index); // xz
        // Upper face
        rB(0, index + 9)  = L2 * rCommonComponents.BMembraneUpper(0, index);  // xx
        rB(1, index + 9)  = L2 * rCommonComponents.BMembraneUpper(1, index);  // yy
        rB(2, index + 9)  = factor_eas * rCommonComponents.BNormal(0, index + 9); // zz
        rB(3, index + 9)  = L2 * rCommonComponents.BMembraneUpper(2, index);  // xy
        rB(4, index + 9)  = L1 * rCommonComponents.BShearLower(1, index + 9) + L2 * rCommonComponents.BShearUpper(1, index + 9); // yz
        rB(5, index + 9)  = L1 * rCommonComponents.BShearLower(0, index + 9) + L2 * rCommonComponents.BShearUpper(0, index + 9); // xz

        /* Neighbour nodes */
        // Lower face
        rB(0, index + 18) = L1 * rCommonComponents.BMembraneLower(0, index + 9); // xx
        rB(1, index + 18) = L1 * rCommonComponents.BMembraneLower(1, index + 9); // yy
        rB(3, index + 18) = L1 * rCommonComponents.BMembraneLower(2, index + 9); // xy
        // Upper face
        rB(0, index + 27) = L2 * rCommonComponents.BMembraneUpper(0, index + 9); // xx
        rB(1, index + 27) = L2 * rCommonComponents.BMembraneUpper(1, index + 9); // yy
        rB(3, index + 27) = L2 * rCommonComponents.BMembraneUpper(2, index + 9); // xy
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::InitializeGeneralVariables(GeneralVariables& rVariables)
{
    // StressMeasure_PK1             //stress related to reference configuration non-symmetric
    // StressMeasure_PK2             //stress related to reference configuration
    // StressMeasure_Kirchhoff       //stress related to current   configuration
    // StressMeasure_Cauchy          //stress related to current   configuration

    // StressMeasure
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN))
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_PK2;
    else
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Doubles
    rVariables.detF  = 1.0;
    rVariables.detF0 = 1.0;
    rVariables.detFT = 1.0;
    rVariables.detJ  = 1.0;

    // Vectors
    noalias(rVariables.StrainVector) = ZeroVector(6);
    noalias(rVariables.StressVector) = ZeroVector(6);
    noalias(rVariables.C) = ZeroVector(6);
    noalias(rVariables.N) = ZeroVector(6);

    // Matrices
    noalias(rVariables.F)  = IdentityMatrix(3);
    noalias(rVariables.F0) = IdentityMatrix(3);
    noalias(rVariables.FT) = IdentityMatrix(3);
    noalias(rVariables.B)  = ZeroMatrix(6, 36);

    noalias(rVariables.DN_DX) = ZeroMatrix(6, 3);
    noalias(rVariables.ConstitutiveMatrix) = ZeroMatrix(6, 6);

    // Reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    // Jacobians
    rVariables.J.resize(1, false);
    rVariables.j.resize(1, false);
    rVariables.J[0] = ZeroMatrix(1, 1);
    rVariables.j[0] = ZeroMatrix(1, 1);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == false ) {
        //Calculate Delta Position
        Matrix delta_position;
        this->CalculateDeltaPosition(delta_position);
        rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, delta_position);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::FinalizeStepVariables(
        GeneralVariables & rVariables,
        const IndexType rPointNumber
        )
{
    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == false ) {
        // Update internal (historical) variables
        mAuxCont[rPointNumber] = rVariables.detF * rVariables.detF0;
        mAuxMatCont[rPointNumber] = prod(rVariables.F, rVariables.F0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::GetHistoricalVariables(
    GeneralVariables& rVariables,
    const IndexType rPointNumber
    )
{
    /* Deformation Gradient F ( set to identity ) */
    const IndexType size =  rVariables.F.size1();

    rVariables.detF  = 1.0;
    rVariables.F     = IdentityMatrix(size);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateGreenLagrangeStrain(
    const Vector& rC,
    Vector& rStrainVector
    )
{
    KRATOS_TRY;

    //Green Lagrange Strain Calculation
    if (rStrainVector.size() != 6)
        rStrainVector.resize(6, false);

    rStrainVector[0] = 0.5 * (rC[0] - 1.00); // xx
    rStrainVector[1] = 0.5 * (rC[1] - 1.00); // yy
    rStrainVector[2] = 0.5 * (rC[2] - 1.00); // zz
    rStrainVector[3] = rC[3]; // xy
    rStrainVector[4] = rC[4]; // yz
    rStrainVector[5] = rC[5]; // xz

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateGreenLagrangeStrain(
        const Matrix& rF,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    // Right Cauchy-Green Calculation
    Matrix C ( 3, 3 );

    noalias( C ) = prod( trans( rF ), rF );

    // Green Lagrange Strain Calculation
    if ( rStrainVector.size() != 6 )
        rStrainVector.resize( 6, false );

    rStrainVector[0] = 0.5 * (C(0, 0) - 1.0); // xx
    rStrainVector[1] = 0.5 * (C(1, 1) - 1.0); // yy
    rStrainVector[2] = 0.5 * (C(2, 2) - 1.0); // zz
    rStrainVector[3] = C(0, 1); // xy
    rStrainVector[4] = C(1, 2); // yz
    rStrainVector[5] = C(0, 2); // xz

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateHenckyStrain(
        const Vector& rC,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    // Declare the different matrix
    bounded_matrix<double, 3, 3> eigen_values_matrix, eigen_vectors_matrix;

    // Assemble matrix C
    const Matrix C_matrix = MathUtils<double>::VectorToSymmetricTensor(rC);

    // Decompose matrix
    MathUtils<double>::EigenSystem<3>(C_matrix, eigen_vectors_matrix, eigen_values_matrix, 1e-24, 10);

    // Calculate the eigenvalues of the E matrix
    eigen_values_matrix(0, 0) = 0.5 * std::log(eigen_values_matrix(0, 0));
    eigen_values_matrix(1, 1) = 0.5 * std::log(eigen_values_matrix(1, 1));
    eigen_values_matrix(2, 2) = 0.5 * std::log(eigen_values_matrix(2, 2));

    // Calculate E matrix
    bounded_matrix<double, 3, 3 > E_matrix;
    noalias(E_matrix) = prod(trans(eigen_vectors_matrix), eigen_values_matrix);
    noalias(E_matrix) = prod(E_matrix, eigen_vectors_matrix);

    // Hencky Strain Calculation
    if (rStrainVector.size() != 6)
        rStrainVector.resize(6, false);

    rStrainVector[0] = E_matrix(0, 0); // xx
    rStrainVector[1] = E_matrix(1, 1); // yy
    rStrainVector[2] = E_matrix(2, 2); // zz
    rStrainVector[3] = 2.0 * E_matrix(0, 1); // xy
    rStrainVector[4] = 2.0 * E_matrix(1, 2); // yz
    rStrainVector[5] = 2.0 * E_matrix(0, 2); // xz

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************

void SprismElement3D6N::CalculateAlmansiStrain(
        const Matrix& rF,
        Vector& rStrainVector
        )
{
    KRATOS_TRY;

    // Tensor Cauchy-Green Calculation
    Matrix tensor_cauchy_green = prod(rF, trans(rF));

    // Calculating the inverse of the jacobian
    Matrix inverse_tensor_cauchy_green (3, 3);
    double det_b = 0.0;
    MathUtils<double>::InvertMatrix(tensor_cauchy_green, inverse_tensor_cauchy_green, det_b);

    // Almansi Strain Calculation
    if ( rStrainVector.size() != 6)
        rStrainVector.resize(6, false);

    rStrainVector[0] = 0.5 * (1.00 - inverse_tensor_cauchy_green(0, 0));
    rStrainVector[1] = 0.5 * (1.00 - inverse_tensor_cauchy_green(1, 1));
    rStrainVector[2] = 0.5 * (1.00 - inverse_tensor_cauchy_green(2, 2));
    rStrainVector[3] = - inverse_tensor_cauchy_green(0, 1); // xy
    rStrainVector[4] = - inverse_tensor_cauchy_green(1, 2); // yz
    rStrainVector[5] = - inverse_tensor_cauchy_green(0, 2); // xz

    KRATOS_CATCH( "" );
}

/**************************** CALCULATE VOLUME CHANGE ******************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
{
    KRATOS_TRY;

    if ( mELementalFlags.Is(SprismElement3D6N::TOTAL_UPDATED_LAGRANGIAN) == true )
        rVolumeChange = 1.0;
    else
        rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    KRATOS_CATCH( "" );
}

/************************* CALCULATE VOLUME ACCELERATION ***************************/
/***********************************************************************************/

void SprismElement3D6N::CalculateVolumeForce(
    Vector& rVolumeForce,
    GeneralVariables& rVariables,
    const double IntegrationWeight
    )
{
    KRATOS_TRY;

    array_1d<double,3> volume_acceleration = ZeroVector(3);
    if (GetProperties().Has( VOLUME_ACCELERATION ))
        volume_acceleration = GetProperties()[VOLUME_ACCELERATION];
    else if( GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        for (unsigned int i_node = 0; i_node < this->GetGeometry().size(); ++i_node)
            volume_acceleration += rVariables.N[i_node] * GetGeometry()[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    // Compute volume change
    double volume_change;
    this->CalculateVolumeChange( volume_change, rVariables );

    rVolumeForce += volume_acceleration * IntegrationWeight * volume_change * GetProperties()[DENSITY];

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("FinalizedStep",mFinalizedStep);
    rSerializer.save("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.save("AuxMatCont",mAuxMatCont);
    rSerializer.save("AuxCont",mAuxCont);
}

/***********************************************************************************/
/***********************************************************************************/

void SprismElement3D6N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("FinalizedStep",mFinalizedStep);
    rSerializer.load("mTotalDomainInitialSize",mTotalDomainInitialSize);
    rSerializer.load("AuxMatCont",mAuxMatCont);
    rSerializer.load("AuxCont",mAuxCont);
}

} // Namespace Kratos.
