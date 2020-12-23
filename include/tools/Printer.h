#pragma once

#include <algorithm>
#include <string>
#include <cmath>
#include <cstdlib>
#include <utility>

#include <TCanvas.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>

#include <hydra/Parameter.h>

namespace dafne {

namespace Print {

void MinimizerStatus(const ROOT::Minuit2::FunctionMinimum& minimum)
{
	std::cout << "***** Minuit2 Minimizer status:" << std::endl;
	// maybe we should add more detiled explanation for the items below 
    std::cout << "- IsValid: " << minimum.IsValid() 
    		  << Form("\t\t\t [%s]", ((minimum.IsValid() == true) ? "OK":"WARN")) << std::endl;
    std::cout << "- HasValidParameters: " << minimum.HasValidParameters() 
    		  << Form("\t\t [%s]", ((minimum.HasValidParameters() == true) ? "OK":"WARN")) << std::endl;
    std::cout << "- HasValidCovariance: " << minimum.HasValidCovariance()  
    		  << Form("\t\t [%s]", ((minimum.HasValidCovariance() == true) ? "OK":"WARN")) << std::endl;
    std::cout << "- HasAccurateCovar: " << minimum.HasAccurateCovar() 
    		  << Form("\t\t [%s]", ((minimum.HasAccurateCovar() == true) ? "OK":"WARN")) 
    		  << "\n  (HasAccurateCovar means the covariance matrix is accurate)" << std::endl;	
    std::cout << "- HasPosDefCovar : " << minimum.HasPosDefCovar() 
    		  << Form("\t\t [%s]", ((minimum.HasPosDefCovar() == true) ? "OK":"WARN")) << std::endl; 
    std::cout << "- HasMadePosDefCovar : " << minimum.HasMadePosDefCovar() 
    		  << Form("\t [%s]", ((minimum.HasMadePosDefCovar() == false) ? "OK":"WARN")) 
    		  << "\n  (HasMadePosDefCovar means the covariance matrix is forced positive-definite)" << std::endl;	
    std::cout << "- HesseFailed: " << minimum.HesseFailed() 
    		  << Form("\t\t [%s]", ((minimum.HesseFailed() == false) ? "OK":"WARN")) << std::endl;
    std::cout << "- HasCovariance: " << minimum.HasCovariance() 
    		  << Form("\t\t [%s]", ((minimum.HasCovariance() == true) ? "OK":"WARN")) << std::endl; // need to find out the difference between HasCovariance() (IsAvalible) and HasValidCovariance() (IsValid)
    std::cout << "- IsAboveMaxEdm: " << minimum.IsAboveMaxEdm() 
    		  << Form("\t\t [%s]", ((minimum.IsAboveMaxEdm() == false) ? "OK":"WARN")) << std::endl;
    std::cout << "- HasReachedCallLimit: " << minimum.HasReachedCallLimit()  
    		  << Form("\t [%s]", ((minimum.HasReachedCallLimit() == false) ? "OK":"WARN")) << std::endl;

    // detailed definition of the items above, could be found in https://root.cern.ch/doc/v606/BasicFunctionMinimum_8h_source.html
}

void CovarianceMatrixStatus(const ROOT::Minuit2::FunctionMinimum& minimum)
{
	std::cout << "***** Covariance matrix status: " << minimum.UserState().CovarianceStatus() << std::endl;

	std::cout << "  * covariance matrix status = -1 : not available (inversion failed or Hesse failed) \n"
				 "  * covariance matrix status = 0 : available but not positive defined \n"
				 "  * covariance matrix status = 1 : covariance only approximate \n"
				 "  * covariance matrix status = 2 : full matrix but forced pos def \n"
				 "  * covariance matrix status = 3 : full accurate matrix \n" << std::endl;
	 
}

template<typename T>
void CheckFCN(T & fcn)
{
    // get vector of parameters
    auto parameters = fcn.GetParameters().GetVariables();

    // check number of free parameters
    int n_free = 0;
    for (size_t i = 0; i < parameters.size(); ++i) {
        if (!(parameters[i]->IsFixed())) n_free++;
    }

    std::cout << "Number of floating parameters in FCN: " << n_free << std::endl;
    if (n_free == 0) {
        std::cout << "All parameters in FCN are fixed. Exit!" << std::endl;
        exit(-1);
    }

    // check if errors of all free parameters are set
    int n_error_unset = 0;
    std::vector<hydra::Parameter*> parameters_error_unset;
    for (size_t i = 0; i < parameters.size(); ++i) {
        if (!(parameters[i]->IsFixed())) {
            if (std::isnan(parameters[i]->GetError()) || parameters[i]->GetError() == 0) {
                n_error_unset++; parameters_error_unset.push_back(parameters[i]);
            }
        }
    }

    std::cout << "Number of floating parameters with unset or zero error: " << n_error_unset << std::endl;
    if (n_error_unset != 0) {
        std::cout << "Following parameters have unset or zero errors:" << std::endl;
        for (size_t i = 0; i < parameters_error_unset.size(); ++i) {
            std::cout << "* " << parameters_error_unset[i]->GetName() << std::endl;
        }
        std::cout << "Exit!" << std::endl;
        exit(-1);
    }


}

void Canvas(const TCanvas &c, std::string outfilename)
{
	// remove illegal char in ROOT macro naming convention
	std::size_t found = outfilename.find_last_of("/");
	std::replace(outfilename.begin()+found, outfilename.end(), '.', '_');
	std::replace(outfilename.begin()+found, outfilename.end(), '-', '_');
	c.Print((outfilename + ".pdf").c_str());
	c.Print((outfilename + ".C").c_str());
}

} // namespace Print

} // namespace dafne

