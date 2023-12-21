# Installation and Work Flow 

[1. Installation](#installation)  
[2. Work flow](#workflow)  
&ensp;&ensp;[2.1 Create a user code file](#createausercodefile)  
&ensp;&ensp;[2.2 Build an amplitude model](#buildanamplitudemodel)  
&ensp;&ensp;[2.3 Configure the detector effect](#configurethedetectoreffect)  
&ensp;&ensp;&ensp;&ensp;[2.3.1 Dalitz-plot efficiency](#dalitzplotefficiency)  
&ensp;&ensp;&ensp;&ensp;[2.3.2 Dalitz-plot resolution](#dalitzplotresolution)  
&ensp;&ensp;&ensp;&ensp;[2.3.3 Decay-time resolution](#decaytimeresolution)  
&ensp;&ensp;&ensp;&ensp;[2.4 Fitting](#fitting)  
&ensp;&ensp;&ensp;&ensp;[2.5 Monte Carlo samples generation](#montecarlosamplesgeneration)  
&ensp;&ensp;&ensp;&ensp;[2.6 Plotting](#plotting)  
&ensp;&ensp;&ensp;&ensp;[2.7 Configuration file](#configurationfile)  



### 1. Installation <a name="installation"></a> 

The DAFNE framework itself is a header only library, and don't need to compile and install. To compile the examples, one can follow the instruction below:

```
mkdir build
cd build
cmake -D CMAKE_C_COMPILER=`which gcc` -D CMAKE_CXX_COMPILER=`which g++` -D HYDRA_INCLUDE_DIR=<path to hydra> -D TCLAP_INCLUDE_PATH=<path to tclap>/include -D EIGEN3_INCLUDE_DIR=<path to eigen3>/include/eigen3 -D TBB_INCLUDE_DIRS=<path to tbb>/include/  -DTBB_tbb_LIBRARY_RELEASE=<path to tbb>/lib64/libtbb.so ../
make -j 10 # if you have 10 cores available for compilation
```

### 2. Work flow <a name="workflow"></a> 

In the following, we briefly introduce the key steps of using DAFNE. If one is starting a new practical work, we suggest to start from one of the example files in `include/*.inl`, and replace the decay model related content. 

#### 2.1 Create a user code file <a name="createausercodefile"></a> 
The first step is to create a file with the suffix ".inl" under the folder `include/`, say `include/test.inl`. Then copy the head files and namespace setting from one of the example `include/*.inl` files, and write a simple `main()` function, like
```
// some includes
using namespace dafne;
// Main
int main(int argc, char** argv)
{
	std::cout << "Hello dafne." << std::endl;
	return 0;
}
```
Then create a ".cpp" file with same file name under `src/`, say `src/test.cpp`. This file should contain only one line, that is to include the ".inl" file created just now. 
```
#include <test.inl>
```
Finally, add the file name to the bottom of `CMakeList.txt`
```
ADD_ANALYSIS(test)
```
After configuring and compiling with `cmake` and `make`, one  will get a preliminary executable based on the DAFNE framework. In the following steps, user only need to edit the ".inl".

#### 2.2 Build an amplitude model <a name="buildanamplitudemodel"></a> 
First, create the parameters:
```
// K*(892)- -> KS Pi-
auto c_re  = hydra::Parameter::Create("reKstm").Value(-1.16356).Error(0.1);
auto c_im  = hydra::Parameter::Create("imKstm").Value(1.19933).Error(0.1);
auto mass  = hydra::Parameter::Create("mKst").Value(0.8936060).Error(0.1);
auto width = hydra::Parameter::Create("wKst").Value(0.0463407).Error(0.1);
auto radius = hydra::Parameter::Create("rBW").Value(1.5).Error(0.1).Fixed();
```
Second, create amplitude components with these parameters: 
```
ThreeBodyPhaseSpace phsp( 1.8645, {0.49767, 0.139570, 0.139570} );

declarg(MSq12, double)
declarg(MSq13, double)

bool BaBarKSpipi = true; 
bool Helicity = false;
auto KSTM_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,BaBarKSpipi,
			    Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstm",
			    "#it{K}*(892)^{#minus}#it{#pi}^{+}");
```
We define a phase space object with the masses of mother ($ D^0 $) and daughters ($ K_S^0 $, $ \pi^+ $, $ \pi^- $). This object will handle the calculations for kinematics.   
We define the type of the two Dalitz-plot variables as `double`, and wrap them as `MSq12` and `MSq13` with `declarg()`[^declarg]. And we write `MSq12` and `MSq13` as amplitude class template parameters, indicating the type of variables to compute. By doing so, one can avoid applying a math function on wrong variables. Suppose when writing the code, one accidentally applys a function $ f(m_{12}, m_{13}) $ on a $ (m_{12},m_{23}) $ or $ (m_{12},t) $ dataset, the compiler will point out there is something wrong.   
The angular momentum of the resonance is 1, and this resonance decays into particle 1 and 3 ($ a $ and $ c $). The `BaBarKSpipi` and `Helicity` are switches for the amplitude definitions under different convention. The name of this amplitude is defined as `Kstm`, and this is also its name in the [configuration file](#configurationfile). The TLatex expression for the plotting, is defined. And one can also define here the color for the automatic plotting:

```
KSTM_amp.SetColor(kYellow+1);
```
Combine created components with `hydra::sum()` <a name="Af"></a> 
```
auto Af =  hydra::sum(KSTM_amp, RHO_amp, OMEGA_amp,K0STP_1430_amp);
```
Then configure the amplitude with a [configuration file](#configurationfile)
```
ConfigFile config(the_path_to_config_file);
config.ConfigureModel(Af);
```
So far, we have built a time-integrated Dalitz-plot amplitude. To convert decay amplitude to decay rate ($ {\rm pdf} $), just write:
```
auto model = rate(Af);
```
To build the time-dependent decay rate with the meson decay amplitude $ {A}_f $ and anti-meson decay amplitude $ \bar{A}_f $, just write:
```
auto tau = hydra::Parameter::Create("tau").Value(Tau::D0).Error(0.0001);
auto x   = hydra::Parameter::Create("x").Value(0.003).Error(0.0001);
auto y   = hydra::Parameter::Create("y").Value(0.006).Error(0.0001);
auto qop = hydra::Parameter::Create("qop").Value(1.0).Error(0.0001);
auto phi = hydra::Parameter::Create("phi").Value(0.0).Error(0.0001);
config.ConfigureMixingParameters(tau, x, y, qop, phi);
auto model_truth_dz =
 time_dependent_rate<Flavor::Positive,DecayTime>(tau,x,y,qop,phi,Af,Abarf);
auto model_truth_db = 
 time_dependent_rate<Flavor::Negative,DecayTime>(tau,x,y,qop,phi,Af,Abarf);
```
Here the decay-time, mixing parameters, CPV parameters are also set. 


#### 2.3 Configure the detector effect <a name="configurethedetectoreffect"></a> 
##### 2.3.1 Dalitz-plot efficiency <a name="dalitzplotefficiency"></a> 
To write an analytical Dalit-plot efficiency, just define following the scheme[^wrap_lambda]:
```
// create parameters
auto a0  = hydra::Parameter::Create("eff_a0").Value(0.0).Error(0.1);
auto a1  = hydra::Parameter::Create("eff_a1").Value(0.0).Error(0.1);
auto a2  = hydra::Parameter::Create("eff_a2").Value(0.0).Error(0.1);
auto a3  = hydra::Parameter::Create("eff_a3").Value(0.0).Error(0.1);
auto center  = hydra::Parameter::Create("eff_center").Value(1.9).Error(0.1);

// define the efficiency plane	
auto efficiency = hydra::wrap_lambda(
[phsp] __hydra_dual__ (unsigned int npar, const hydra::Parameter* params, 
MSq12 m2p, MSq13 m2m) 
	
	// judge whether in phase space or not
	if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
	
	double center = params[4];
	double x = m2m - center;
	double y = m2p - center;
	
	double a0 = params[0];
	double a1 = params[1];
	double a2 = params[2];
	double a3 = params[3];
	
	return 1.0 + a0 * (x+y)
		   + a1 * x * y
		   + a2 * (pow(x,2)+pow(y,2)) 
		   + a3 * (x-y);
	

}, a0, a1, a2, a3, center); 

// configure efficiency according to configuration file
config.Debug();
config.ConfigureEfficiency(efficiency);
```
Here users can define arbitrary number of parameters for the efficiency, with a self-defined expression.   
For the discrete Dalitz-plot efficiency, it is convenient to configure with the `ArbitraryBinningHistogram2D` class:
```
ArbitraryBinningHistogram2D efficiency_hist =
							 config.ConfigureEfficiencyHistogram();
auto efficiency = hydra::wrap_lambda(
[phsp, efficiency_hist] __hydra_dual__ (MSq12 m2p, MSq13 m2m) {

    // judge whether in phase space or not
    if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;

    double m2z = phsp.MSqJK(m2p, m2m);
    double helicity_z = phsp.HelicityJK<2,3>(m2p, m2m);

    double x = m2z;
    double y = helicity_z;
    return efficiency_hist.GetValue(x, y);

}); 
```
Multiply on the ${\rm pdf}$:
```
auto model = rate(Af)*efficiency;
```
##### 2.3.2 Dalitz-plot resolution <a name="dalitzplotresolution"></a> 
The Dalitz-plot resolution can be configure with `Resolution2D`. 
```
Resolution2D resolution2D;
config.ConfigureResolution2D(resolution2D);
TCanvas cresolution("cresolution", "cresolution", 800, 600);
```
This class could provide convenient function for the smearing, as shown in `include/generate-KSpipi-time-independent.inl`.

##### 2.3.3 Decay-time resolution <a name="decaytimeresolution"></a> 
For the decay-time resolution, we need to set the parameters for the decay-time pull distribution and the distribution of $ \sigma_t $, which is a Johnson SU distribution. 
```
auto b = hydra::Parameter::Create("b").Value(0.0).Error(0.0001);
auto s = hydra::Parameter::Create("s").Value(1.0).Error(0.0001);

auto johnson_delta = hydra::Parameter::Create().Name("johnson_delta" ).Value(1.65335e+00);
auto johnson_lambda = hydra::Parameter::Create().Name("johnson_lambda").Value(1.87922e-02);
auto johnson_gamma = hydra::Parameter::Create().Name("johnson_gamma" ).Value(-2.57429e+00);
auto johnson_xi = hydra::Parameter::Create().Name("johnson_xi").Value(4.27580e-02);

config.ConfigureTimeResolutionParameters({&b, &s, 
       &johnson_gamma, &johnson_delta, &johnson_xi, &johnson_lambda});

auto johnson_su = hydra::JohnsonSU<DecayTimeError>(
       johnson_gamma, johnson_delta, johnson_xi, johnson_lambda);
```
When building a time-dependent $ {\rm pdf} $ with detector effect considered,  just follow the form:
```
auto model_dz = 
	time_dependent_rate_with_time_resolution_pdf<Flavor::Positive,
	MSq12,MSq13,DecayTime,DecayTimeError>(
	tau,x,y,qop,phi,b,s,efficiency,Af,Afbar,johnson_su,
	{phsp.MSqMin<1,2>(),phsp.MSqMax<1,2>()},
	{phsp.MSqMin<1,3>(),phsp.MSqMax<1,3>()},
	phsp.TimeRange(),
	phsp.TimeErrorRange()); 
auto model_db = 
	time_dependent_rate_with_time_resolution_pdf<Flavor::Negative,
	MSq12,MSq13,DecayTime,DecayTimeError>(
	tau,x,y,qop,phi,b,s,efficiency,Abarf,Abarfbar,johnson_su,
	{phsp.MSqMin<1,2>(),phsp.MSqMax<1,2>()},
	{phsp.MSqMin<1,3>(),phsp.MSqMax<1,3>()},
	phsp.TimeRange(),
	phsp.TimeErrorRange()); 
```
The ranges of dalitz variables and decay-time variables are also set, which are needed by the normalization factor. Another thing to notice is that the `phsp` here is an object of `ThreeBodyPhaseSpaceWithTimeAndTimeError` class, a subclass of `ThreeBodyPhaseSpace`.
```
ThreeBodyPhaseSpaceWithTimeAndTimeError phsp( 1.8645, 
     {0.49767, 0.139570, 0.139570}, 
     {-2., 7.}, {0, 0.5});
```
Here we set the decay-time range as $ [-2, 7]~{\rm ps} $ and decay-time error range as $ [0, 0.5]~{\rm ps} $.

#### 2.4 Fitting <a name="fitting"></a> 
Create a container for the data and fill it. For example, fill the data from a CERN ROOT TTree.
```
hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23>, 
                   hydra::device::sys_t> data;
for (auto i=0; i<nentries; ++i) {
	ntp->GetEntry(i);
	data.push_back(hydra::make_tuple(MSq12(m2p), MSq13(m2m), MSq23(m2z)));
}
```
Similarly, users can build a dataset with additional decay-time or decay-time error dimensions. Normalize the $ {\rm pdf} $ and build a likelihood with the data. Check if the likelihood (FCN) and the parameters inside are set correctly.
```
auto pdf = hydra::make_pdf( model, phsp.Integrator(10*ncands) );
auto fcn = hydra::make_loglikehood_fcn( pdf, data.begin(), data.end() );
MinuitTools::CheckFCN(fcn);
```
To fit multiple datasets simultaneously, just write
```
auto merged_fcn = hydra::make_simultaneous_fcn(fcn1, fcn2);
// fcn1 and fcn2 are created with hydra::make_loglikehood_fcn()
MinuitTools::CheckFCN(merged_fcn);
```
Run the Minuit and print the fit results.
```
ROOT::Minuit2::MnPrint::SetLevel(2);
ROOT::Minuit2::MnStrategy strategy(2);
// starting values of parameters
auto pstart = fcn.GetParameters().GetMnState();
ROOT::Minuit2::MnMigrad migrad(fcn, pstart, strategy);
auto minimum = ROOT::Minuit2::FunctionMinimum( migrad(5000., 1.));

if ( !minimum.IsValid() || !minimum.HasAccurateCovar() ) {
	std::cout << "Fit did not converge or covariance matrix is not accurate.";
	std::cout << std::endl;
	return -2;
}

std::cout << "Fit results:\n" << minimum.UserState() << std::endl;
MinuitTools::CovarianceMatrixStatus(minimum);
```

#### 2.5 Monte Carlo samples generation <a name="montecarlosamplesgeneration"></a> 
To generate MC according to the time-integrated Dalitz-plot, just write
```
auto data = phsp.GenerateData<MSq12, MSq13, MSq23>(model,
            events_number,random_seed);
```
To generate a time-dependent Dalitz-plot MC sample, just write
```
auto data_dz = phsp.GenerateDataWithTime<MSq12, MSq13, MSq23, 
               DecayTime>(model_dz, tau(), y(), events_number, random_seed);
auto data_db = phsp.GenerateDataWithTime<MSq12, MSq13, MSq23, 
               DecayTime>(model_db, tau(), y(), events_number, random_seed);
```
Here, users need to provide also the mixing parameter $ y $, which is useful for a fast sampling procedure inside the `GenerateDataWithTime()`. The `phsp` here is an object of `ThreeBodyPhaseSpaceWithTime` class.
```
ThreeBodyPhaseSpaceWithTime phsp(1.8645, {0.49767, 0.139570, 0.139570}, {0., 7.});
```
The range of decay-time is set to $ [0,7]~{\rm ps} $. After considering the decay-time resolution, the generation of a time-dependent Dalitz-plot sample is a bit different, due to some technical consideration. Just write
```
auto data_dz = phsp.GenerateDataWithTimeAndTimeError<MSq12, MSq13, MSq23,
	DecayTime,DecayTimeError>(model_truth_dz, efficiency, tau(), y(), b(), s(), 
	johnson_su, events_number, random_seed);
auto data_db = phsp.GenerateDataWithTimeAndTimeError<MSq12, MSq13, MSq23,
    DecayTime,DecayTimeError>(model_truth_db, efficiency, tau(),y(), b(), s(),
    johnson_su, events_number, random_seed);
```
The `phsp` here is of  `ThreeBodyPhaseSpaceWithTimeAndTimeError` class. The `model_truth_dz`  should be created from the `time_dependent_rate()` rather than `time_dependent_rate_with_time_resolution_pdf()`.

#### 2.6 Plotting <a name="plotting"></a> 
To plot a time-integrated Daltz-plot, we firstly load the data and the decay model $ {\rm pdf} $, to a `DalitzPlotter`, setting the TLatex expression of three daughter particles. For example
```
auto plotter = DalitzPlotter<MSqZero,MSqMinus,MSqPlus>(phsp,"#it{K}^{#minus}",
               "#it{#pi}^{+}","#it{#pi}^{0}");
plotter.FillDataHistogram(data);
plotter.FillModelHistogram(model);
```
Load the histograms of the components. The component histograms should be loaded with the sum decay amplitudes [Af](#Af). The efficiency and signal fraction are also required, and can be set to constant and $ 100\% $ for the truth level simulation. The `FillComponentHistograms()` will not load the "removed" components in the [configuration file](#configurationfile). 
```
auto efficiency = ConstantFunctor(1);
double signal_fraction = 1; // pure signal
plotter.FillComponentHistograms(Af, efficiency, signal_fraction);
```
Plot the one dimensional distribution with a selected axis. 
```
int plottedAxis = 2; // 0 for m_ab^2, 1 for m_ac^2, 2 for m_bc^2. 
                     // The helicity axis is not implemented.
bool plotWithLegend = 1;
plotter.Plot1DProjections(plottedAxis, plotWithLegend);
```
Plot the two dimensional distribution of data with selected axes
```
int plottedXAxis = 0;  
int plottedYAxis = 1; 
plotter.Plot2DProjectionData(plottedXAxis, plottedYAxis);
```
Two dimensional distribution of the fit model can be ploted separately:
```
plotter.Plot2DProjectionModel(plottedXAxis, plottedYAxis);
```
The pull between data and model can also be plotted easily:
```
plotter.Plot1DPull(plottedAxis);
plotter.Plot2DPull(plottedXAxis, plottedYAxis);
```
The CERN ROOT macro file for the plots will also be generated for further modification. The ranges and binnings for plotting, can be set with the parameters of `FillModelHistograms()`, `FillDataHistograms()` and `FillComponentHistograms()`. Details can be found in `include/tools/DalitzPlotter.h`.

#### 2.7 Configuration file <a name="configurationfile"></a> 
Users can configure the DAFNE program by the configuration file. The status of amplitude components can be set as follows:
```
amplitudes_setting{
# name     status
#------------------
Kstm       floating
rho0       floating
omega      floating
K0stm_1430 removed
}amplitudes_setting
```
The first column is the component name, which should be identical to the component name in the program. The second column is the status of the component. "floating" means to add this component to the Dalitz-plot amplitude. "removed" means this component is not added to the Dalitz-plot amplitude. Users need to fix the parameters of "removed" components manually. Users should set the shared parameters between a "removed" component and a "floating" component to fixed or floating, according to the requirement of the fit. Before writing the configuration file, users need to write all possible components in the program. The configuration file can not add a component that is not in the program.

The parameters can be set as follows:
```
amplitude_parameters_setting{

# Coefficients of resonances
# name       = value  +-  uncertainty  limits_option
#-------------------------------------
rerho0        =  1.0      +- 0.1     C
reKstm       = -1.196090 +- 0.005755 F
reomega      = -0.019235 +- 0.000670 E 10
reK0stm_1430 = -0.393116 +- 0.037100 C
		
imrho0        =  0.0      +- 0.1     C
imKstm       =  1.256890 +- 0.006278 F
imomega      =  0.037376 +- 0.000529 E 10
imK0stm_1430 = -5.284940 +- 0.030800 C

}amplitude_parameters_setting
```
Each row configures a real parameter. The first column is the name, which should be identical to the parameter name in the program. The second and third are the mean value and error. The rightmost column(s) is the option for the limits. "C" for constant/fixed, "F" for floating with no constraint, "V x" for floating by $ \pm x\% $ of the mean value, "E x" for floating by $ \pm x\sigma $, and "L l h" for directly setting the limits to be $(l,h)$. The mixing parameters, decay time resolution and background parameters can be set similarly.

The discrete Dalitz-plot efficiency can be set as follows:
```
efficiency_histogram_setting{

XTicks: 0,0.1,0.3,0.6,1.9,
YTicks: -1,-0.84,-0.6,0.6,0.84,1,
	
Zs: 
0.120000,0.100000,0.100000,0.110000,
0.125000,0.110000,0.110000,0.115000,
0.121000,0.120000,0.119000,0.119000,
0.121000,0.120000,0.114000,0.117000,
0.111000,0.105000,0.104000,0.112000,
	
ZErrors: 
0.006000,0.001300,0.001100,0.000600,
0.005000,0.001200,0.000900,0.000600,
0.002500,0.000600,0.000400,0.000300,
0.005000,0.001200,0.000900,0.000600,
0.006000,0.001400,0.001000,0.000700,
	
}efficiency_histogram_setting
```
"XTicks"("YTicks") set the boundaries of x(y) axis bins. "Zs" set the efficiency of each bin. "ZErrors" set the efficiency errors. With these lines, user is setting a Dalitz-plot efficiency like this

<img src=".\figs\dalitz_plot_efficiency.png" width = "350" align=center>

The Dalitz-plot resolution can be set as follows:
```
resolution2D_setting{
	# xmin    xmax    ymin    ymax    rmsx    rmsy    rho
	0.000000  0.400000  -1.010000 -0.600000 0.001400  0.002000  -0.120
	0.000000  0.400000  -0.600000 -0.200000 0.001500  0.003000  -0.080
	0.000000  0.400000  -0.200000  0.200000 0.001600  0.004000  0.000
	0.000000  0.400000   0.200000  0.600000 0.001500  0.003000  0.080
	0.000000  0.400000   0.600000  1.010000 0.001400  0.002000  0.120
	0.400000  0.800000  -1.010000 -0.600000 0.002500  0.002000  -0.040
	0.400000  0.800000  -0.600000 -0.200000 0.002600  0.003000  -0.020
	0.400000  0.800000  -0.200000  0.200000 0.002700  0.004000  0.000
	0.400000  0.800000   0.200000  0.600000 0.002600  0.003000  0.020
	0.400000  0.800000   0.600000  1.010000 0.002500  0.002000  0.040
	0.800000  1.200000  -1.010000 -0.600000 0.003500  0.002000  0.040
	0.800000  1.200000  -0.600000 -0.200000 0.003600  0.003000  0.030
	0.800000  1.200000  -0.200000  0.200000 0.003700  0.004000  0.000
	0.800000  1.200000   0.200000  0.600000 0.003600  0.003000  -0.030
	0.800000  1.200000   0.600000  1.010000 0.003500  0.002000  -0.040
	1.200000  1.600000  -1.010000 -0.600000 0.003500  0.002000  0.050
	1.200000  1.600000  -0.600000 -0.200000 0.003600  0.003000  0.025
	1.200000  1.600000  -0.200000  0.200000 0.003700  0.004000  0.000
	1.200000  1.600000   0.200000  0.600000 0.003600  0.003000  -0.025
	1.200000  1.600000   0.600000  1.010000 0.003500  0.002000  -0.050
	1.600000  2.000000  -1.010000 -0.600000 0.003700  0.002000  0.300
	1.600000  2.000000  -0.600000 -0.200000 0.003800  0.003000  0.100
	1.600000  2.000000  -0.200000  0.200000 0.003900  0.004000  0.000
	1.600000  2.000000   0.200000  0.600000 0.003800  0.003000  -0.100
	1.600000  2.000000   0.600000  1.010000 0.003700  0.002000  -0.300    
}resolution2D_setting
```
The first four columns define the regions. The last three columns set the resolutions of the x axis variables and the correlation. These lines are setting Dalitz-plot resolution like this

<img src=".\figs\dalitz_plot_resolution.png" width = "350" align=center>




[^declarg]: This is a macro definition of the Hydra framework. For more technical details, please refer to `hydra/detail/FunctionArgument.h` in the Hydra framework.
[^wrap_lambda]: Here, the `hydra::wrap_lamda` creates a C++ functor. For more details, please refer to the document of the Hydra framework.
