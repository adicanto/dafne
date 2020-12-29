#pragma once

#include "tclap/CmdLine.h"
#include "tclap/ArgException.h"

namespace dafne {

class Arguments {
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
	std::string  contour_option;

	Arguments(std::string exename, std::string exever, std::string definput="") : name(exename), version(exever), input(definput) {;}

	void Read(const int& argc, char** argv) {
		TCLAP::CmdLine cmd(name,' ',version);

		// Define value arguments
		TCLAP::ValueArg<std::string> iArg("i", "input", "Input file.", false, input, "file" );
		TCLAP::ValueArg<std::string> cArg("c", "config_file", "Model configuration file.", true, config_file, "file" );
		TCLAP::ValueArg<std::string> oArg("o", "outdir", "Output file directory/pattern.", false, "./", "path/pattern" );
		TCLAP::SwitchArg pArg("p", "plot", "Plot results?", false);
		TCLAP::ValueArg<int> nbArg("", "plotnbins","NBins for the plots", false, 50, "int");
		TCLAP::ValueArg<int> lArg("", "printlevel","Printout level.", false, 1, "int");
		TCLAP::SwitchArg gArg("g", "gentoy", "Generate toy data?", false);
		TCLAP::ValueArg<unsigned int> nArg("n", "nevents", "Number of events to be generated.", false, 100000, "unsigned int" );
		TCLAP::ValueArg<unsigned int> sArg("", "seed", "Seed used in the generation.", false, 12345, "unsigned int" );
		TCLAP::SwitchArg xArg("", "interactive", "Plot results in an interactive session.", false);

		std::vector<std::string> allowedContourOptions;
		allowedContourOptions.push_back("");
		allowedContourOptions.push_back("test");
		allowedContourOptions.push_back("sigma1");
		allowedContourOptions.push_back("sigma3");
		allowedContourOptions.push_back("sigma5");
		allowedContourOptions.push_back("sigmafull");
		TCLAP::ValuesConstraint<std::string> contourConstraints( allowedContourOptions );
		TCLAP::ValueArg<std::string> coArg("", "contouroption", "option for contour plotting", false, "", &contourConstraints);



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
		cmd.add(coArg);

		// Parse the args.
		cmd.parse(argc, argv);

		// Get the value parsed by each arg
		input   	  = iArg.getValue();
		config_file= cArg.getValue();
		outdir  	  = oArg.getValue();
		plot    	  = pArg.getValue();
		plotnbins  	  = nbArg.getValue();
		interactive= xArg.getValue();
		prlevel 	  = lArg.getValue();
		gentoy  	  = gArg.getValue();
		nevents 	  = nArg.getValue();
		seed    	  = sArg.getValue();
		contour_option = coArg.getValue();
		
		auto last = outdir[outdir.size()-1];
		if ( last!='/' && last!='_' ) outdir += "_";
		
		if (interactive) plot = true;
	}

	std::ostream &Print(std::ostream &out=std::cout) {
		out << "Running with following arguments: " << std::endl;
		out << "\tInput file: " << input << std::endl;
		out << "\tModel configuration file: " << config_file << std::endl;
		out << "\tOutput file directory/pattern: " << outdir << std::endl;
		out << "\tPlot results? " << (plot ? "Yes" : "No") << std::endl;
		out << "\tContour Option: " << contour_option << std::endl;
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
		return out;
	}
};

} // namespace dafne
