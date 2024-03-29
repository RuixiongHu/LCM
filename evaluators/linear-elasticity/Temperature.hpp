//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef LCM_TEMPERATURE_HPP
#define LCM_TEMPERATURE_HPP

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
  class Temperature : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {

  public:

    Temperature(Teuchos::ParameterList& p,
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

    PHX::MDField<ScalarT,Cell,QuadPoint> temperature_;

    unsigned int num_qps_;
    unsigned int num_dims_;
    unsigned int num_nodes_;
    unsigned int workset_size_;

    bool enable_transient_;
    std::string temperature_Name_;
  

  };
}

#endif
