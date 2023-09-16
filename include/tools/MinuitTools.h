#pragma once

#include <algorithm>
#include <string>
#include <cmath>
#include <cstdlib>
#include <utility>

#include <TCanvas.h>
#include <TGraph.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnContours.h>
#include <Minuit2/MnPlot.h>

#include <hydra/Parameter.h>
#include <hydra/LogLikelihoodFCN.h>


#include <tools/FunctorTools.h>
#include <physics/ThreeBodyPhaseSpace.h>



namespace dafne {

namespace MinuitTools {

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
    std::map<std::string, int> free_parameters;
    for (size_t i = 0; i < parameters.size(); ++i) {
        if (!(parameters[i]->IsFixed())) free_parameters[parameters[i]->GetName()] = 1;
    }

    int n_free = free_parameters.size();
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

hydra::Parameter * GetParameterPointer(std::vector<hydra::Parameter*> & parameters, const char* name)
{
    for (size_t i = 0; i < parameters.size(); ++i) {
        if (std::string(name) == std::string(parameters[i]->GetName())) return parameters[i];
    }

    return NULL;
}


template<typename FUNCTOR, typename MINIMUM>
void CheckMinimizerParameters(FUNCTOR functor, MINIMUM & minimum)
{
    // check if some plotting model parameters are not in the fitting result from the fitter (minimum)
    auto fitter_parameters = minimum.UserParameters().Parameters();
    std::vector<std::string> fitter_parameter_names;
    for (auto param : fitter_parameters) fitter_parameter_names.push_back(std::string(param.Name()));

    size_t n_pars_functor = functor.GetNumberOfParameters();
    for (size_t i = 0; i < n_pars_functor; ++i ) {
        if (find(fitter_parameter_names.begin(), fitter_parameter_names.end(), functor.GetParameter(i).GetName()) == fitter_parameter_names.end()) {
            std::cout << "The parameter \""<< functor.GetParameter(i) << "\" is not found in the ";
            std::cout << "fit result. Please have a check. Exit." << std::endl;
        }
    } 
}

// build a dummy fcn to easily synchronize the paramters in plotting functor and fitting functor, according
// to the parameter names. 
template<typename DUMMYDATA, typename FUNCTOR, typename MINIMUM>
auto UpdateParametersByBuildingDummyFCN(DUMMYDATA & data_dummy, FUNCTOR functor, MINIMUM & minimum)
{
    auto dummy_fcn_for_plotting = hydra::make_loglikehood_fcn( 
                                    hydra::make_pdf( functor, ConstantIntegrator<hydra::host::sys_t>(1.0)) ,
                                    data_dummy.begin(),
                                    data_dummy.end()
                                  );

    hydra::UserParameters& parameters = dummy_fcn_for_plotting.GetParameters();
    parameters.UpdateParameters(minimum);
    return dummy_fcn_for_plotting.GetPDF().GetFunctor();
}

template<typename MSq12, typename MSq13, typename MSq23, typename FUNCTOR, typename MINIMUM>
auto UpdateParameters(ThreeBodyPhaseSpace phsp, FUNCTOR functor, MINIMUM & minimum)
{
    // When calling, please make sure all parameters in the model (functor) are included in the 
    // fitting result from the fitter (minimum), otherwise this function will go into
    // segmentfault.


    // check if some plotting model parameters are not in the fitting result from the 
    // fitter (minimum). Under test.
    // CheckMinimizerParameters(functor, minimum);

    // build a dummy fcn(likelihood) to easily synchronize the paramters from the fitter to 
    // the plotting model
    // The building a likelihood needs data points, so build a dummy dataset and add an 
    // event inside the phase space
    auto data_dummy = phsp.GenerateDummyData<MSq12, MSq13, MSq23>(10); 
    auto normalized_model_for_plotting = UpdateParametersByBuildingDummyFCN(data_dummy, functor, minimum);

    return normalized_model_for_plotting;
}

template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename FUNCTOR, typename MINIMUM>
auto UpdateParameters(ThreeBodyPhaseSpaceWithTime phsp, FUNCTOR functor, MINIMUM & minimum)
{
    // When calling, please make sure all parameters in the model (functor) are included in the 
    // fitting result from the fitter (minimum), otherwise this function will go into
    // segmentfault.

    // build a dummy fcn(likelihood) to easily synchronize the paramters from the fitter to 
    // the plotting model
    // The building a likelihood needs data points, so build a dummy dataset and add an 
    // event inside the phase space
    auto data_dummy = phsp.GenerateDummyDataWithTime<MSq12, MSq13, MSq23, Time>(10); 
    auto normalized_model_for_plotting = UpdateParametersByBuildingDummyFCN(data_dummy, functor, minimum);

    return normalized_model_for_plotting;
}


template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename TimeError, typename FUNCTOR, typename MINIMUM>
auto UpdateParameters(ThreeBodyPhaseSpaceWithTimeAndTimeError phsp, FUNCTOR functor, MINIMUM & minimum)
{
    // When calling, please make sure all parameters in the model (functor) are included in the 
    // fitting result from the fitter (minimum), otherwise this function will go into
    // segmentfault.

    // build a dummy fcn(likelihood) to easily synchronize the paramters from the fitter to 
    // the plotting model
    // The building a likelihood needs data points, so build a dummy dataset and add an 
    // event inside the phase space
    auto data_dummy = phsp.GenerateDummyDataWithTimeAndTimeError<MSq12, MSq13, MSq23, Time, TimeError>(10); 
    auto normalized_model_for_plotting = UpdateParametersByBuildingDummyFCN(data_dummy, functor, minimum);

    return normalized_model_for_plotting;
}


TGraph* ContourToTGraph(std::vector<std::pair<double,double> > contourPoints, const char * name="contour", const bool percentage=0)
{
    std::vector<double> xs(contourPoints.size()+1);
    std::vector<double> ys(contourPoints.size()+1);

    for (size_t i = 0; i < contourPoints.size(); ++i) {
        xs[i] = contourPoints[i].first;
        ys[i] = contourPoints[i].second;
        if(percentage) { xs[i] *= 100; ys[i] *= 100;}
    }
    xs[contourPoints.size()] = xs[0];
    ys[contourPoints.size()] = ys[0];

    TGraph* tg = new TGraph(contourPoints.size()+1, xs.data(), ys.data());
    tg->SetName(name);
    return tg; 
}

/**
 * Fucntion to plot contour corresponding to certain deltaLogL. 
 * 
 * In hydra framework, FCN classes use -InL instead of -2InL:
 *
 *  https://github.com/MultithreadCorner/Hydra/blob/master/hydra/detail/LogLikelihoodFCN1.inl
 *  https://github.com/MultithreadCorner/Hydra/blob/master/hydra/detail/functors/LogLikelihood1.h 
 *
 * Some common examples are:
 *
 *  // for 1 sigma, TMath::ChisquareQuantile(1 - RooStats::SignificanceToPValue(1)*2, 2)/2 == 1.1478745
 *  PlotContour(fcn, minimum, 1.1478745); 
 *
 *  // for 2 sigma, TMath::ChisquareQuantile(1 - RooStats::SignificanceToPValue(2)*2, 2)/2 == 3.0900372
 *  PlotContour(fcn, minimum, 3.0900372); 
 *
 *  // for 3 sigma, TMath::ChisquareQuantile(1 - RooStats::SignificanceToPValue(3)*2, 2)/2 == 5.9145790
 *  PlotContour(fcn, minimum, 5.9145790); 
 *
 *  // for 4 sigma, TMath::ChisquareQuantile(1 - RooStats::SignificanceToPValue(4)*2, 2)/2 == 9.6669543
 *  PlotContour(fcn, minimum, 9.6669543); 
 *
 *  // for 5 sigma, TMath::ChisquareQuantile(1 - RooStats::SignificanceToPValue(5)*2, 2)/2 == 14.371851
 *  PlotContour(fcn, minimum, 14.371851); 
 *
 *  // for 68% CL, TMath::ChisquareQuantile(0.68, 2)/2 == 1.1394343
 *  PlotContour(fcn, minimum, 0.70); 
 *
 *  // for 70% CL, TMath::ChisquareQuantile(0.70, 2)/2 == 1.2039728 
 *  PlotContour(fcn, minimum, 0.70); 
 *
 *  // for 90% CL, TMath::ChisquareQuantile(0.90, 2)/2 == 2.3025851 
 *  PlotContour(fcn, minimum, 0.90);
 *
 *  // for 95% CL, TMath::ChisquareQuantile(0.95, 2)/2 == 2.9957323 
 *  PlotContour(fcn, minimum, 0.95); 
 *
 * The drawOption is mostly root TGraphPainter draw option, but with additional "silence" option to return TGraph* only, 
 * without plotting it, for example:
 *  
 *  PlotContour(fcn, minimum, 1.1478745, "silence");   
 *
 */
template <typename FCN, typename MINIMUM>
TGraph* Contour(FCN & fcn, MINIMUM & minimum, const char* graphName, const char* xname, const char* yname, const double deltaLogL, const int npoints, const bool percentage, const char* drawOption="al") 
{
    double previousUp = fcn.GetErrorDef(); // the up value before calling PlotContour()

    ROOT::Minuit2::MnContours contours(fcn, minimum);
    fcn.SetErrorDef(deltaLogL);
    std::vector<std::pair<double,double> > cont = contours(minimum.UserState().Index(xname), minimum.UserState().Index(yname), npoints);

    TGraph* tg = ContourToTGraph(cont, graphName, percentage);

    if (strstr(drawOption, "silence") == NULL) tg->Draw(drawOption);

    fcn.SetErrorDef(previousUp);

    return tg;
}

} // namespace MinuitTools

} // namespace dafne

