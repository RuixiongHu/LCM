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
  Temperature<EvalT, Traits>::
  Temperature(Teuchos::ParameterList& p,
       const Teuchos::RCP<Albany::Layouts>& dl) :
    temperature_        (p.get<std::string>("Temperature Name"),
		  dl->qp_scalar)
  {

    this->addEvaluatedField(temperature_);
 
    Teuchos::RCP<PHX::DataLayout> scalar_dl = dl->qp_scalar;
    std::vector<PHX::Device::size_type> dims;
    scalar_dl->dimensions(dims);
    workset_size_ = dims[0];
    num_qps_      = dims[1];

    temperature_Name_ = p.get<std::string>("Temperature Name")+"_old";
    this->setName("Temperature"+PHX::typeAsString<EvalT>());
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  void Temperature<EvalT, Traits>::
  postRegistrationSetup(typename Traits::SetupData d,
			PHX::FieldManager<Traits>& fm)
  {
    this->utils.setFieldData(temperature_,fm);
  }

  //**********************************************************************

  template<typename EvalT, typename Traits>
  void Temperature<EvalT, Traits>::
  evaluateFields(typename Traits::EvalData workset)
  {

    //grab old value
    Albany::MDArray temperature_old = (*workset.stateArrayPtr)[temperature_Name_];
    
    // current time
    const RealType t = workset.current_time;
    //***********************************************************
    // actual simulation
    /*
    if (t==0.0){
      auto vector_space = workset.disc->getVectorSpace();
      Teuchos::RCP<Thyra_Vector> my_field = Thyra::createMember(vector_space);
      Teuchos::RCP<Albany::AbstractDiscretization> discPtr = workset.disc;
      auto APFDisc_ptr = Teuchos::rcp_dynamic_cast<Albany::APFDiscretization>(discPtr);
      apf::Mesh2* m = APFDisc_ptr -> getMesh();
      apf::Field* f = m -> findField("Temperature_old") ;
      if (f==NULL){
	std::cout<<"Temperature field not found! Keep the old instead\n";
	for (std::size_t cell = 0; cell < workset.numCells; ++cell){
	  for (std::size_t qp = 0; qp < num_qps_; ++qp){
	    if (temperature_old(cell,qp)<5){  // if using the small scale, then convert it back to kelvin 
	      temperature_(cell,qp) = temperature_old(cell,qp)*1000;
	    }
	    else{
	      temperature_(cell,qp) = temperature_old(cell,qp);
	    }
	  }
        }
      }
      else{
        APFDisc_ptr -> copyQPScalarFromAPF(1,"Temperature",f);
        APFDisc_ptr -> copyQPScalarFromAPF(1,"Temperature_old",f);        
        std::cout<<"Temperature Field copied correctly and destroyed\n";
        apf::destroyField( f );
      }
    }

    else {
      for (std::size_t cell = 0; cell < workset.numCells; ++cell){
	for (std::size_t qp = 0; qp < num_qps_; ++qp){
	  temperature_(cell,qp) = temperature_old(cell,qp);
	}
      }      
    }
    */

        // defining temperature_
        for (std::size_t cell = 0; cell < workset.numCells; ++cell){
            for (std::size_t qp = 0; qp < num_qps_; ++qp){
              temperature_(cell,qp) = 300.0;
              temperature_old(cell,qp) = 300.0;              

	        }
	    }
    
  
  //std::cout << "psi1 has been finished\n" ; 
  
  }

  //**********************************************************************

}
