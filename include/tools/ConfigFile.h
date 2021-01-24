#pragma once

#include <string>
#include <fstream>
#include <cstdlib>
#include <streambuf>
#include <sstream>
#include <utility>

#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/functions/Utils.h>
#include <hydra/Placeholders.h>
#include <hydra/detail/Sum.h>

#include <TPRegexp.h>

#include <physics/ThreeBodyPhaseSpace.h>
#include <physics/Resolution.h>
#include <tools/ArbitraryBinningHistogram2D.h>
#include <ctype.h>

namespace dafne {

class ConfigFile
{
private:
	std::string _configFileName;
	bool _debug;

	// regular expression tools
	TPMERegexp _lineSplitor; // split a paragraph to lines
	TPMERegexp _commentFinder; // filter comment line
	TPMERegexp _blankLineFilter; // filter blank line

	std::vector<std::string> _amplitudeLines;
	std::vector<std::string> _amplitudeParameterLines;
	std::vector<std::string> _mixingParameterLines;
	std::vector<std::string> _efficiencyParameterLines;
	std::vector<std::string> _resolution2DLines;
	std::vector<std::string> _pdfParameterLines; 
	std::vector<std::string> _efficiencyHistogramLines; 


	std::vector<std::string> _readBlocks; // to record blocks that are already read, avoiding reading twice

	void parse_block(const char* blockName, std::vector<std::string> & lines)
	{
		if (_debug) std::cout << "--- Reading setting block " << blockName << std::endl;
		
		if (std::find(_readBlocks.begin(), _readBlocks.end(), std::string(blockName)) != _readBlocks.end()) {
			std::cout << "WARNING: The setting block " << blockName << " has been already read, skipping it now." << std::endl;
			return;
		}

		// read the setting block from .cfg file
		TPMERegexp settingBlock(Form("%s{(.*\n)+}%s", blockName, blockName));

		std::ifstream configFileStream(_configFileName);
		std::string configFileStr((std::istreambuf_iterator<char>(configFileStream)), std::istreambuf_iterator<char>());
		if (settingBlock.Match(configFileStr) < 1) {
			if (_debug) std::cout << "WARNING: Cannot find setting block " << blockName << " in " << _configFileName << std::endl;
			return;
		}
			
		_lineSplitor.Split(settingBlock[0]);
		for (int i=1; i < _lineSplitor.NMatches()-1; ++i)
		{
			if (_commentFinder.Match(_lineSplitor[i]) != 0) continue; // skip comment line
			if (_blankLineFilter.Match(_lineSplitor[i]) == 0) continue; // skip blank line
			lines.push_back(_lineSplitor[i].Data());

			if (_debug) std::cout << _lineSplitor[i] << std::endl;
		}

		_readBlocks.push_back(blockName);
	}
	
	int find_line_index(std::vector<std::string> & lines, const std::string targetItemName) const
	{
		for (size_t i = 0; i < lines.size(); ++i)
		{
			std::istringstream is(lines[i]);
			std::string currentItemName;
			is >> currentItemName;
			if (currentItemName == targetItemName) return i;
		}
		return -1;
	}
	
	template<typename T>
	void set_amplitude(T & amplitude, std::string statusString)
	{
		if (statusString == "floating") {
			for (size_t i = 0; i < amplitude.GetNumberOfParameters(); ++i) {
				set_parameter(amplitude.Parameter(i),_amplitudeParameterLines);
			}
			amplitude.Update();
		} else if (statusString == "fixed") {
			for (size_t i = 0; i < amplitude.GetNumberOfParameters(); ++i) {
				set_parameter(amplitude.Parameter(i),_amplitudeParameterLines);
				// amplitude.Parameter(i).Fixed(1); // commented for the moment, because of the logic problem for multiple amplitudes' sharing parameters
			}
			amplitude.Update();
		} else if (statusString == "removed") {
			amplitude.Remove();
			// amplitude.FixParameters(); // commented for the moment, because of the logic problem for multiple amplitudes' sharing parameters
		} else {
			std::cout << "WARNING: Cannot recognize status " << statusString << " for amplitude " << amplitude.Name() << ", keep it floating" << std::endl;
		}
	}
		
	void set_parameter(hydra::Parameter & par, std::vector<std::string> & lines)
	{
		std::string parName = par.GetName();

		// search parameter in the config file
		const int i_par = find_line_index(lines, parName);
		if (i_par == -1) {
			std::cout << "WARNING: Cannot find parameter " << parName << " in " << _configFileName << ", using default value/settings" << std::endl;
			return;
		}

		// read parameter setting
		std::istringstream is_par(lines[i_par]);

		double value;
		double error;
		std::string separator;
		std::string statusString;

		is_par >> separator;
		is_par >> separator;

		// set value
		is_par >> value;
		par.SetValue(value);

		// set error
		is_par >> separator;
		if (separator != "+-") {
			std::cout << "WARNING: The \"+-\" syntax is not in the right place for parameter " << parName << ", keeping default uncertainty" << std::endl;
			return;
		}
		is_par >> error;
		par.SetError(error);

		// set status/limits
		is_par >> statusString;
		if (statusString == "C") { // constant/fixed parameter
			par.Fixed(1);
		} else if (statusString == "F") { // floating with no constraints
			par.Fixed(0);
		} else if (statusString == "V") { // vary by a relative x%
			double scaleFactor;
			is_par >> scaleFactor;
			par.Fixed(0);
			par.Limits(value*(1.-scaleFactor), value*(1.+scaleFactor));
		} else if (statusString == "E") { // vary by +-x sigma
			double scaleFactor;
			is_par >> scaleFactor;
			par.Fixed(0);
			par.Limits(value-scaleFactor*error, value+scaleFactor*error);
		} else if (statusString == "L") { // custom limits
			double lowerLimit;
			double upperLimit;
			is_par >> lowerLimit >> upperLimit;
			par.Fixed(0);
			par.Limits(lowerLimit, upperLimit);
		} else {
			std::cout << "WARNING: Cannot recognize status " << statusString << " for parameter " << parName << ", keeping default status" << std::endl;
			return;
		}
		
		if (_debug) std::cout << "Parameter " << parName << " set to " << value << " +- " << error << " " << statusString << std::endl;
	}
	
	template <typename T, int ...Is>
	void configure_amplitudes(T & model, std::integer_sequence<int, Is...>)
	{
		// static loop over all amplitudes
		((configure_single_amplitude<Is>(model)), ...);
	}

	template<unsigned int I, typename T>
	void configure_single_amplitude(T & model)
	{
		hydra::placeholders::placeholder<I> _X;
		ConfigureSingleAmplitude(model, _X);
	}
	
public:
	ConfigFile() = delete;

	ConfigFile(std::string configFileName, bool debug=false) :
		_configFileName(configFileName),
		_debug(debug),
		_lineSplitor("\n"),
		_commentFinder("^\\s*#"),
		_blankLineFilter("\\w")
	{
		if (configFileName.empty()) {
			std::cout << "No configuration file provided, exit." << std::endl;
			exit(-1);
		}
		
	}

	void Debug(bool value=true)
	{
		_debug = value;
	}

	// the ParseXXXSetting functions below could be unified to one function

	void ParseAmplitudeSetting()
	{
		if (_amplitudeLines.size() == 0) parse_block("amplitudes_setting", _amplitudeLines);
		if (_mixingParameterLines.size() == 0) parse_block("amplitude_parameters_setting", _amplitudeParameterLines);
	}

	void ParseMixingSetting()
	{
		if (_mixingParameterLines.size() == 0) parse_block("mixing_parameters_setting", _mixingParameterLines);
	}

	void ParseEfficiencySetting()
	{
		if (_efficiencyParameterLines.size() == 0) parse_block("efficiency_parameters_setting", _efficiencyParameterLines);
	}

	void ParseEfficiencyHistogramSetting()
	{
		if (_efficiencyHistogramLines.size() == 0) parse_block("efficiency_histogram_setting", _efficiencyHistogramLines);
	}

	void ParseResolution2DSetting()
	{
		if (_resolution2DLines.size() == 0) parse_block("resolution2D_setting", _resolution2DLines);
	}

	void ParsePdfSetting()
 	{
 		if (_pdfParameterLines.size() == 0) parse_block("pdf_parameters_setting", _pdfParameterLines);
 	}


	template<typename ...Fs>
	void ConfigureModel(hydra::Sum<Fs...> & model)
	{
		ParseAmplitudeSetting();
		configure_amplitudes(model, std::make_integer_sequence<int, sizeof...(Fs)>{});
	}
	
	template<typename T, unsigned int I>
	void ConfigureSingleAmplitude(T & model, hydra::placeholders::placeholder<I> const& X)
	{
		auto & amplitude = model.GetFunctor(X);
		
		if (_debug) std::cout << "--- Configuring amplitude " << amplitude.Name() << std::endl;
		
		const std::string amplitudeName = amplitude.Name();
		const int i_res = find_line_index(_amplitudeLines, amplitudeName);
		if (i_res == -1) {
			std::cout << "WARNING: Cannot find amplitude " << amplitudeName << " in " << _configFileName << ", continue while ignoring it." << std::endl;
			return;
		}

		std::istringstream is_res(_amplitudeLines[i_res]);
		std::string name;
		std::string statusString;

		is_res >> name;
		is_res >> statusString;
		
		if (_debug) std::cout << "Setting status for amplitude " << amplitudeName << " to " << statusString << std::endl;
		set_amplitude(amplitude, statusString);
	}

	void ConfigureMixingParameters(hydra::Parameter &tau, hydra::Parameter &x, hydra::Parameter &y, hydra::Parameter &qop, hydra::Parameter &phi)
	{
		ParseMixingSetting();

		if (_debug) std::cout << "--- Configuring mixing parameters" << std::endl;
		
		set_parameter(tau, _mixingParameterLines);
		set_parameter(  x, _mixingParameterLines);
		set_parameter(  y, _mixingParameterLines);
		set_parameter(qop, _mixingParameterLines);
		set_parameter(phi, _mixingParameterLines);
	}

	template<typename T>
	void ConfigureTimeDependentModel(T & model, std::vector<std::string> parnames={"tau","x","y","qop","phi"})
	{
		if (_debug) std::cout << "--- Configuring time-dependent model" << std::endl;

		for (auto name : parnames)
			set_parameter( model.Parameter(name.c_str()), _mixingParameterLines);
	}

	template<typename Efficiency>
	void  ConfigureEfficiency(Efficiency & efficiency)
	{
		ParseEfficiencySetting();
		for (size_t i = 0; i < efficiency.GetNumberOfParameters(); ++i) {
			set_parameter(efficiency.Parameter(i), _efficiencyParameterLines);
		}
		efficiency.Update();
	}

	void  ConfigureResolution2D(Resolution2D & resolution2D)
	{
		ParseResolution2DSetting();
		for (auto line : _resolution2DLines) {
			// read parameter setting
			std::istringstream is(line);

			double xmin;
			double xmax;
			double ymin;
			double ymax;
			double rmsx;
			double rmsy;
			double rho;

			is >> xmin;
			is >> xmax;
			is >> ymin;
			is >> ymax;
			is >> rmsx;
			is >> rmsy;
			is >> rho;

			resolution2D.Push(xmin, xmax, ymin, ymax, rmsx, rmsy, rho);
		}

	}

	// config the fcn (or the pdf in fcn) with hydra::UserParameters interface
	template<typename T>
	void ConfigurePdfFCN(T & fcn)
	{
		ParsePdfSetting();

		auto parameters = fcn.GetParameters().GetVariables();
		
		if (_debug) std::cout << "--- Configuring FCN " << std::endl;
		
		for (size_t i = 0; i < parameters.size(); ++i) {
				set_parameter(*(parameters[i]),_pdfParameterLines);
		}

		fcn.GetParameters().SetVariables(parameters);
	}

    int NumberOfDigits(std::string const str)
    {
      int n = 0;
      for (unsigned i = 0; i < str.size(); ++i) {
        if (std::isdigit(str[i])) ++n;
      }
      return n;
    }

	std::vector<double> ParseLineToVector(std::string const line, const int nprefix=0) 
	{
		std::istringstream is(line);

		if (_debug) std::cout << "converting line: " << line << " to istringstream ... ..." << std::endl;

		std::vector<double> vs;
		std::string str_sub;

		for (int i = 0; i < nprefix; ++i) {
			is >> str_sub;
			if (_debug) std::cout << "skipping prefix: " << str_sub << std::endl;
		}

		while(getline(is, str_sub, ',')) {
			std::istringstream is_sub(str_sub);

			if (_debug) std::cout << "getline: " << str_sub << std::endl;

			if (NumberOfDigits(str_sub) == 0) continue;

			double v;
			is_sub >> v;

			vs.push_back(v);
		}

		return vs;

	};


	ArbitraryBinningHistogram2D ConfigureEfficiencyHistogram()
	{
		ParseEfficiencyHistogramSetting();

		std::vector<double> xticks;
		std::vector<double> yticks;
		std::vector< std::vector<double> > zs;
		std::vector< std::vector<double> > zerrors;

		auto Transpose = [](std::vector< std::vector<double> > vin)->std::vector< std::vector<double> > {
			int nx = vin[0].size();
			int ny = vin.size();

			std::vector< std::vector<double> > vout;
			std::vector<double> init_line(ny, 0);	
			for (int i = 0; i < nx; ++i) vout.push_back(init_line);

			for (int ix = 0; ix < nx; ++ix)
			for (int iy = 0; iy < ny; ++iy) {
				vout[ix][iy] = vin[iy][ix];
			}

			return vout;
		};

		int nx = -99;
		int ny = -99;

		unsigned i_eff = 0;

		while (i_eff < _efficiencyHistogramLines.size()) {
			std::istringstream is_eff(_efficiencyHistogramLines[i_eff]);
			std::string indicator;

			is_eff >> indicator;

			if (_debug) std::cout << "handling: " << indicator << std::endl;

			if (indicator == "XTicks:") {

				xticks = ParseLineToVector(is_eff.str(), 1);
				nx = xticks.size()-1; // need to add some judgement here to judgement the dimensions are consistent
				i_eff++;
				continue;

			} else if (indicator == "YTicks:") {

				yticks = ParseLineToVector(is_eff.str(), 1);
				ny = yticks.size()-1;
				i_eff++;
				continue;

			} else if (indicator == "Zs:") {

				i_eff++;

				if (nx < 0 || ny < 0) {
					std::cout << "ny and nx seems not set, when reading the Zs. something must be wrong!" << std::endl;
					exit(-1);
				}

				std::vector< std::vector<double> > zs_in;
				for (int i = 0; i < ny; ++i) {
					std::istringstream is_eff_sub(_efficiencyHistogramLines[i_eff+i]);
					zs_in.push_back(ParseLineToVector(is_eff_sub.str()));
				}

				zs = Transpose(zs_in);

				i_eff = i_eff + ny;

				continue;

			} else if (indicator == "ZErrors:") {

				i_eff++;

				if (nx < 0 || ny < 0) {
					std::cout << "ny and nx seems not set, when reading the ZErrors. something must be wrong!" << std::endl;
					exit(-1);
				}


				std::vector< std::vector<double> > zerrors_in;
				for (int i = 0; i < ny; ++i) {
					std::istringstream is_eff_sub(_efficiencyHistogramLines[i_eff+i]);
					zerrors_in.push_back(ParseLineToVector(is_eff_sub.str()));
				}

				zerrors = Transpose(zerrors_in);

				i_eff = i_eff + ny;

				continue;
			} else {
				std::cout << "wrong indicator: " << indicator << "\nexit!" << std::endl;
				exit(-1);
			}

			if (_debug) std::cout << "finish handling: " << indicator << std::endl;


		}

		ArbitraryBinningHistogram2D hist(xticks, yticks, zs, zerrors);
		return hist;


	}


};

} // namespace dafne
