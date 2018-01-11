//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_STATIC_STEP_METHOD )
#define  KRATOS_STATIC_STEP_METHOD

// System includes

// External includes

// Project includes
#include "custom_strategies/time_integration_methods/static_method.hpp"

namespace Kratos
{
  ///@addtogroup SolidMechanicsApplication
  ///@{
  
  ///@name Kratos Globals
  ///@{
  
  ///@}
  ///@name Type Definitions
  ///@{
  
  ///@}
  ///@name  Enum's
  ///@{
  
  ///@}
  ///@name  Functions
  ///@{
  
  ///@}
  ///@name Kratos Classes
  ///@{

 
  /// Short class definition.
  /** Detail class definition.     
   * This class performs predict and update of dofs variables, their time derivatives and time integrals      
   */
  template<class TVariableType, class TValueType>
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) StaticStepMethod : public StaticMethod<TVariableType,TValueType>
  {   
  public:
 
    ///@name Type Definitions
    ///@{

    /// BaseType
    typedef TimeIntegrationMethod<TVariableType,TValueType>  BaseType;

    /// BaseTypePointer
    typedef typename BaseType::Pointer                BaseTypePointer;
    
    /// NodeType
    typedef typename BaseType::NodeType                      NodeType;
    
    /// KratosVariable or KratosVariableComponent    
    typedef typename BaseType::VariablePointer        VariablePointer;

    /// DerivedType
    typedef StaticMethod<TVariableType,TValueType>        DerivedType;

    
    KRATOS_CLASS_POINTER_DEFINITION( StaticStepMethod );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    StaticStepMethod() : DerivedType()
    {
      mpStepVariable = nullptr;
    }

    /// Copy Constructor.
    StaticStepMethod(StaticStepMethod& rOther)
      :DerivedType(rOther)
      ,mpStepVariable(rOther.mpStepVariable)
    {
    }

    /// Clone.
    BaseTypePointer Clone()
    {
      return BaseTypePointer( new StaticStepMethod(*this) );
    }

    /// Destructor.
    ~StaticStepMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // has step variable
    virtual bool HasStepVariable() override
    {
      return true;
    }
    
    // set step variable (step variable)
    virtual void SetStepVariable(const TVariableType& rStepVariable) override
    {
      mpStepVariable = &rStepVariable;
    }
    
    // predict
    virtual void Predict(NodeType& rNode) override
    {
     KRATOS_TRY
     
     DerivedType::Predict(rNode);

     this->PredictStepVariable(rNode);
	
     KRATOS_CATCH( "" )
    }


    
    // update
     virtual void Update(NodeType& rNode) override
    {
     KRATOS_TRY
       
     DerivedType::Update(rNode);

     this->UpdateStepVariable(rNode);


     KRATOS_CATCH( "" )
    }

     
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{


    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "StaticStepMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StaticStepMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "StaticStepMethod Data";     
    }

    
    ///@}
    ///@name Friends
    ///@{


    ///@}
    
  protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    // method variables    
    VariablePointer mpStepVariable;
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void PredictStepVariable(NodeType& rNode)
    {
      KRATOS_TRY

      this->UpdateStepVariable(rNode);
	
      KRATOS_CATCH( "" )
    }


    virtual void UpdateStepVariable(NodeType& rNode)
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
	
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      CurrentStepVariable = CurrentVariable-PreviousVariable;
	
      KRATOS_CATCH( "" )
    }
    
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
  
    ///@}

  private:

    ///@name Static Member Variables
    ///@{
  
    ///@}
    ///@name Member Variables
    ///@{
  
    ///@}
    ///@name Private Operators
    ///@{
  
    ///@}
    ///@name Private Operations
    ///@{
  
    ///@}
    ///@name Private  Access
    ///@{
  
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      // rSerializer.save("StepVariable", mpStepVariable);
    };

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      // rSerializer.load("StepVariable", mpStepVariable);
    };
    
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@}
  
  }; // Class StaticStepMethod
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  
  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, StaticStepMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const StaticStepMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }
  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_STATIC_STEP_METHOD defined
