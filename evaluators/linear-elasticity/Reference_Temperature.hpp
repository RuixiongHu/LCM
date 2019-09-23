//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef LCM_REFERENCE_TEMPERATURE_HPP
#define LCM_REFERENCE_TEMPERATURE_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Sacado_ParameterAccessor.hpp"
#include "Teuchos_Array.hpp"
#include "Albany_Layouts.hpp"

namespace LCM {
  template<typename EvalT, typename Traits>
  class Reference_Temperature : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {

  public:
    Reference_Temperature(Teuchos::ParameterList& p,
	 const Teuchos::RCP<Albany::Layouts>& dl);

    void 
    postRegistrationSetup(typename Traits::SetupData d,
			  PHX::FieldManager<Traits>& vm);

    void 
    evaluateFields(typename Traits::EvalData d);

  private:

    typedef typename EvalT::ScalarT ScalarT;
    typedef typename EvalT::MeshScalarT MeshScalarT;

    double constant_value_;
    double room_temperature_;
    double melting_temperature_;

    PHX::MDField<ScalarT,Cell,QuadPoint> psi1_;
    PHX::MDField<ScalarT,Cell,QuadPoint> reference_temperature_;
    PHX::MDField<ScalarT,Cell,QuadPoint> temperature_;
    PHX::MDField<MeshScalarT,Cell,QuadPoint,Dim> coord_;

    unsigned int num_qps_;
    unsigned int num_dims_;
    unsigned int num_nodes_;
    unsigned int workset_size_;

    bool enable_transient_;
    std::string psi1_Name_;
    std::string reference_temperature_Name_;
    std::string temperature_Name_;
  
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidReferenceTemperatureParameters() const; 

  };
}

#endif
