//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#include <fstream>
#include "Sacado_ParameterRegistration.hpp"
#include "Albany_Utils.hpp"
 
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

#include <PCU.h>
#include <pumi.h>
#include <apf.h>
#include <pcu_util.h>
#include "apfShape.h"
#include "lionPrint.h"
#include "apfDynamicArray.h"
#include "apfMesh2.h"


#include "Albany_Application.hpp"
#include "Albany_APFMeshStruct.hpp"
#include "Albany_APFDiscretization.hpp"

namespace LCM {

  //**********************************************************************
  template<typename EvalT, typename Traits>
  Psi1<EvalT, Traits>::
  Psi1(Teuchos::ParameterList& p,
       const Teuchos::RCP<Albany::Layouts>& dl) :
    psi1_        (p.get<std::string>("Psi1 Name"),
		  dl->qp_scalar)
  {

    this->addEvaluatedField(psi1_);
 
    Teuchos::RCP<PHX::DataLayout> scalar_dl = dl->qp_scalar;
    std::vector<PHX::Device::size_type> dims;
    scalar_dl->dimensions(dims);
    workset_size_ = dims[0];
    num_qps_      = dims[1];

    Teuchos::ParameterList* cond_list =
      p.get<Teuchos::ParameterList*>("Parameter List");

    Teuchos::RCP<const Teuchos::ParameterList> reflist = 
      this->getValidPsi1Parameters(); 

    cond_list->validateParameters(*reflist, 0, 
				  Teuchos::VALIDATE_USED_ENABLED, Teuchos::VALIDATE_DEFAULTS_DISABLED); 

    constant_value_ = cond_list->get("Psi1", 0.0); 
  

    psi1_Name_ = p.get<std::string>("Psi1 Name")+"_old";
    this->setName("Psi1"+PHX::typeAsString<EvalT>());
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  void Psi1<EvalT, Traits>::
  postRegistrationSetup(typename Traits::SetupData d,
			PHX::FieldManager<Traits>& fm)
  {
    this->utils.setFieldData(psi1_,fm);
  }

  //**********************************************************************

  template<typename EvalT, typename Traits>
  void Psi1<EvalT, Traits>::
  evaluateFields(typename Traits::EvalData workset)
  {
        //std::cout << "psi1 has started\n" ; 
                //std::cout << "initialPsi1"<< constant_value_ << "\n" ;

    std::cout<<"this is in Psi1, has "<< workset.numCells << " cells\n";
    //grab old psi value
    Albany::MDArray psi1_old = (*workset.stateArrayPtr)[psi1_Name_];
    
    // current time
    const RealType t = workset.current_time;
    //****************************************************************
    // for actual simluation
    /*
    if (t==0.0){
      auto vector_space = workset.disc->getVectorSpace();
      Teuchos::RCP<Thyra_Vector> my_field = Thyra::createMember(vector_space);
      Teuchos::RCP<Albany::AbstractDiscretization> discPtr = workset.disc;
      auto APFDisc_ptr = Teuchos::rcp_dynamic_cast<Albany::APFDiscretization>(discPtr);
      apf::Mesh2* m = APFDisc_ptr -> getMesh();
      apf::Field* f = m -> findField("Psi1") ;
      if (f==NULL){
	std::cout<<"Psi1 field not found! Keep the old instead\n";
	for (std::size_t cell = 0; cell < workset.numCells; ++cell){
	  for (std::size_t qp = 0; qp < num_qps_; ++qp){
	    psi1_(cell,qp) = psi1_old(cell,qp);
          }
        }
      }
      else{
        APFDisc_ptr -> copyQPScalarFromAPF(1,"Psi1",f);
        APFDisc_ptr -> copyQPScalarFromAPF(1,"Psi1_old",f);        
        std::cout<<"Field copied correctly and destroyed\n";
        apf::destroyField( f );
      }
    }

    else {
      for (std::size_t cell = 0; cell < workset.numCells; ++cell){
	for (std::size_t qp = 0; qp < num_qps_; ++qp){
	  psi1_(cell,qp) = psi1_old(cell,qp);
	}
      }      
    }
      
  */

 
    // do this only at the equilibration step

    //original set up
    if (t == 0.0){
        // initializing psi_ values 
        for (std::size_t cell = 0; cell < workset.numCells; ++cell){
            for (std::size_t qp = 0; qp < num_qps_; ++qp){
                psi1_(cell, qp) = constant_value_;
	        }
	    }
    }
    

  
  }
  //**********************************************************************
  template<typename EvalT, typename Traits>
  Teuchos::RCP<const Teuchos::ParameterList>
  Psi1<EvalT, Traits>::
  getValidPsi1Parameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> valid_pl =
      rcp(new Teuchos::ParameterList("Valid Psi1 Params"));;

    valid_pl->set<double>("Psi1", 1.0);
    return valid_pl;
  }
  //**********************************************************************

}
