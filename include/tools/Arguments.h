#pragma once

#include "tclap/CmdLine.h"
#include "tclap/ArgException.h"

namespace dafne {

class BaseArguments {
public:
	std::string  name;
	std::string  version;
	std::string  input;
	std::string  config_file;
	std::string  outdir;
	bool         plot;
	unsigned int plotnbins;
	int          prlevel;
	bool         gentoy;
	unsigned int nevents;
	unsigned int seed;
	bool         interactive;
	TCLAP::CmdLine cmd;

	TCLAP::ValueArg<std::string> iArg;
	TCLAP::ValueArg<std::string> cArg;
	TCLAP::ValueArg<std::string> oArg;
	TCLAP::SwitchArg pArg;
	TCLAP::ValueArg<int> nbArg;
	TCLAP::ValueArg<int> lArg;
	TCLAP::SwitchArg gArg;
	TCLAP::ValueArg<unsigned int> nArg;
	TCLAP::ValueArg<unsigned int> sArg;
	TCLAP::SwitchArg xArg;

	BaseArguments(std::string exename, std::string exever, std::string definput="") : name(exename), 
		version(exever), input(definput), cmd(name,' ',version), 
		iArg("i", "input", "Input file.", false, input, "file" ),
		cArg("c", "config_file", "Model configuration file.", true, config_file, "file" ),
		oArg("o", "outdir", "Output file directory/pattern.", false, "./", "path/pattern" ),
		pArg("p", "plot", "Plot results?", false),
		nbArg("", "plotnbins","NBins for the plots", false, 50, "int"),
		lArg("", "printlevel","Printout level.", false, 1, "int"),
		gArg("g", "gentoy", "Generate toy data?", false),
		nArg("n", "nevents", "Number of events to be generated.", false, 100000, "unsigned int" ),
		sArg("", "seed", "Seed used in the generation.", false, 12345, "unsigned int" ),
		xArg("", "interactive", "Plot results in an interactive session.", false)
	{		
		// Add arguments to the command line
		cmd.add(sArg);
		cmd.add(nArg);
		cmd.add(gArg);
		cmd.add(lArg);
		cmd.add(xArg);
		cmd.add(pArg);
		cmd.add(oArg);
		cmd.add(cArg);
		cmd.add(iArg);
	}

	void LoadBasicArguments()
	{
		// Get the value parsed by each arg
		input   	  = iArg.getValue();
		config_file	  = cArg.getValue();
		outdir  	  = oArg.getValue();
		plot    	  = pArg.getValue();
		plotnbins  	  = nbArg.getValue();
		interactive   = xArg.getValue();
		prlevel 	  = lArg.getValue();
		gentoy  	  = gArg.getValue();
		nevents 	  = nArg.getValue();
		seed    	  = sArg.getValue();

		if (interactive) plot = true;

		auto last = outdir[outdir.size()-1];
		if ( last!='/' && last!='_' ) outdir += "_";
	}

	void PrintBasicArguments(std::ostream &out=std::cout) 
	{
		out << "Running with following arguments: " << std::endl;
		out << "\tInput file: " << input << std::endl;
		out << "\tModel configuration file: " << config_file << std::endl;
		out << "\tOutput file directory/pattern: " << outdir << std::endl;
		out << "\tPlot results? " << (plot ? "Yes" : "No") << std::endl;
		out << "\tPlot Nbins " << plotnbins << std::endl;
		if (plot) {
			out << "\tUse interactive session? " << (interactive ? "Yes" : "No") << std::endl;
		}
		out << "\tPrintout level: " << prlevel << std::endl;
		out << "\tGenerate toy data? " << (gentoy ? "Yes" : "No") << std::endl;
		if (gentoy) {
			out << "\tNumber of events to be generated: " << nevents << std::endl;
			out << "\tSeed: " << seed << std::endl;
		}
	}

	virtual void LoadExtendArguments() = 0;
	virtual void PrintExtendArguments(std::ostream &out=std::cout) = 0;

	void Read(const int& argc, char** argv) 
	{
		cmd.parse(argc, argv);
		LoadBasicArguments();
		LoadExtendArguments();
	}

	std::ostream &Print(std::ostream &out=std::cout) 
	{
		PrintBasicArguments(out);
		PrintExtendArguments(out);

		return out;
	}
};

class StandardArguments: public BaseArguments 
{
public:
	StandardArguments(std::string exename, std::string exever, std::string definput="") : BaseArguments(exename, exever, definput) {;}

	void LoadExtendArguments() {};
	void PrintExtendArguments(std::ostream &out=std::cout) {};
};

class D0ToKsPiPiArguments: public BaseArguments
{
public:
	std::string contour_option;
	TCLAP::ValueArg<std::string> *coArg;
	TCLAP::ValuesConstraint<std::string> *contourConstraints;

public:
	D0ToKsPiPiArguments(std::string exename, std::string exever, std::string definput="") : BaseArguments(exename, exever, definput) 
	{
		std::vector<std::string> allowedContourOptions;
		allowedContourOptions.push_back("");
		allowedContourOptions.push_back("test");
		allowedContourOptions.push_back("sigma1");
		allowedContourOptions.push_back("sigma3");
		allowedContourOptions.push_back("sigma5");
		allowedContourOptions.push_back("sigmafull");
		contourConstraints = new TCLAP::ValuesConstraint<std::string>( allowedContourOptions );

		coArg = new TCLAP::ValueArg<std::string>("", "contouroption", "option for contour plotting", false, "", contourConstraints);

		cmd.add(*coArg);
	}

	~D0ToKsPiPiArguments() {
		delete coArg;
		delete contourConstraints;
	}

	void LoadExtendArguments()
	{
		contour_option = coArg->getValue();
	}

	void PrintExtendArguments(std::ostream &out=std::cout) {
		out << "\tContour Option: " << contour_option << std::endl;
	};


};

} // namespace dafne
