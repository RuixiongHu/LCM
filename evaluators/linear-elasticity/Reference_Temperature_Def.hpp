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
  Reference_Temperature<EvalT, Traits>::
  Reference_Temperature(Teuchos::ParameterList& p,
       const Teuchos::RCP<Albany::Layouts>& dl) :
    psi1_        (p.get<std::string>("Psi1 Name"),
		  dl->qp_scalar),
    reference_temperature_ (p.get<std::string>("Reference_Temperature Name"), 
			    dl->qp_scalar),
    coord_        (p.get<std::string>("QP Coordinate Vector Name"), dl->qp_vector),    
    temperature_ (p.get<std::string>("Temperature Name"), dl->qp_scalar)
    
  {

    this->addDependentField(psi1_);
    this->addDependentField(temperature_);
    this->addDependentField(coord_);
    this->addEvaluatedField(reference_temperature_);

    Teuchos::RCP<PHX::DataLayout> scalar_dl = dl->qp_scalar;
    std::vector<PHX::Device::size_type> dims;
    scalar_dl->dimensions(dims);
    workset_size_ = dims[0];
    num_qps_      = dims[1];

    Teuchos::ParameterList* cond_list =
      p.get<Teuchos::ParameterList*>("Parameter List");

    Teuchos::RCP<const Teuchos::ParameterList> reflist = 
      this->getValidReferenceTemperatureParameters(); 

    cond_list->validateParameters(*reflist, 0, 
				  Teuchos::VALIDATE_USED_ENABLED, Teuchos::VALIDATE_DEFAULTS_DISABLED); 

    constant_value_ = cond_list->get("Value", 0.0); 
    room_temperature_ = cond_list->get("Room Temperature",0.0);
    melting_temperature_ = cond_list->get("Melting Temperature",0.0);
  

    psi1_Name_ = p.get<std::string>("Psi1 Name")+"_old";
    reference_temperature_Name_ = p.get<std::string>("Reference_Temperature Name")+"_old";
    temperature_Name_ = p.get<std::string>("Temperature Name")+"_old";
    this->setName("Reference_Temperature"+PHX::typeAsString<EvalT>());
  }

  //**********************************************************************
  template<typename EvalT, typename Traits>
  void Reference_Temperature<EvalT, Traits>::
  postRegistrationSetup(typename Traits::SetupData d,
			PHX::FieldManager<Traits>& fm)
  {
    this->utils.setFieldData(psi1_,fm);
    this->utils.setFieldData(reference_temperature_, fm);
    this->utils.setFieldData(temperature_,fm);
    this->utils.setFieldData(coord_,fm);
  }

  //**********************************************************************

  template<typename EvalT, typename Traits>
  void Reference_Temperature<EvalT, Traits>::
  evaluateFields(typename Traits::EvalData workset)
  {
        //std::cout << "psi1 has started\n" ; 
                //std::cout << "initialPsi1"<< constant_value_ << "\n" ;

    //grab old value
    Albany::MDArray psi1_old = (*workset.stateArrayPtr)[psi1_Name_];
    Albany::MDArray reference_temperature_old = (*workset.stateArrayPtr)[reference_temperature_Name_];
    ScalarT reference_temperature_bar_;

    // current time
    const RealType t = workset.current_time;
    /*
    if (t==0.0){
      for (std::size_t cell = 0; cell < workset.numCells; ++cell){
	 for (std::size_t qp = 0; qp < num_qps_; ++qp){
	  if(psi1_old(cell,qp)>0.98){
	    reference_temperature_(cell,qp)=  std::max(melting_temperature_, temperature_(cell,qp));
	    reference_temperature_bar_ = std::max(melting_temperature_, temperature_(cell,qp));
	    //reference_temperature_old(cell,qp)=  reference_temperature_bar_;
	  }
	  else{
	    reference_temperature_(cell,qp)=  room_temperature_;
	    reference_temperature_old(cell,qp)=  room_temperature_;
	  }
	}
      }
    }
    else {
      for (std::size_t cell = 0; cell < workset.numCells; ++cell){
	for (std::size_t qp = 0; qp < num_qps_; ++qp){
	  reference_temperature_(cell,qp)= reference_temperature_old(cell,qp);
	}
      }
    }
    */ 
    if (t==0.0){
      for (std::size_t cell = 0; cell < workset.numCells; ++cell){
        for (std::size_t qp = 0; qp < num_qps_; ++qp){
          if(coord_(cell, qp, 2) > 250.0){
            reference_temperature_(cell,qp)=  melting_temperature_;
      //reference_temperature_old(cell,qp)=  reference_temperature_bar_;
          }
          else{
            reference_temperature_(cell,qp)=  room_temperature_;
          }
        }
      }
    }

  
  }
  //**********************************************************************
  template<typename EvalT, typename Traits>
  Teuchos::RCP<const Teuchos::ParameterList>
  Reference_Temperature<EvalT, Traits>::
  getValidReferenceTemperatureParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> valid_pl =
      rcp(new Teuchos::ParameterList("Valid Reference_Temperature Params"));;

    valid_pl->set<double>("Value", 1.0);
    valid_pl->set<double>("Room Temperature", 1.0);
    valid_pl->set<double>("Melting Temperature", 1.0);
    return valid_pl;
  }
  //**********************************************************************

}
