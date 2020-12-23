#pragma once

// -----------------------------------------
// Particle masses (GeV/c2) as of PDG2019
// -----------------------------------------
namespace Mass
{
	const double E     = 0.0005109989461;
	const double Mu    = 0.1056583745;
	const double Tau   = 1.77686;
	const double K     = 0.493677;
	const double Pi    = 0.13957061;
	const double Pi0   = 0.1349770;
	const double D0    = 1.86483;
	const double Dstar = 2.01026;
	const double Dstar0= 2.00685;
	const double D     = 1.86965;
	const double Ds    = 1.96834;
	const double P     = 0.9382720467;
	const double N     = 0.939565413;
	const double Phi   = 1.019461;
	const double B     = 5.27933;
	const double B0    = 5.27964;
	const double Bs    = 5.36688;
	const double K0    = 0.497611;
	const double L0    = 1.115683;
	const double Lc    = 2.28646;

	const double Eta   = 0.54730;
	const double EtaPrime = 0.95777;
}

// -----------------------------------------
// Particle lifetimes (ps) as of PDG2019
// -----------------------------------------
namespace Tau
{
	const double D0  = 0.4101;
	const double D   = 1.040;
	const double Ds  = 0.504;
	const double Lc  = 0.200;
	const double B   = 1.638;
	const double B0  = 1.519;
	const double Bs  = 1.510;
	const double BsH = 1.619;
	const double BsL = 1.414;
	const double BsFS= 1.527;
	const double KS  = 89.54;
	const double KL  = 5116.;
}

// -----------------------------------------
// Particle MC Codes
// -----------------------------------------
namespace PdgCode
{
	const int Gamma  = 22;
	const int E      = 11;
	const int NuE    = 12;
	const int Mu     = 13;
	const int NuMu   = 14;
	const int Tau    = 15;
	const int NuTau  = 16;
	const int K      = 321;
	const int Pi     = 211;
	const int Pi0    = 111;
	const int D0     = 421;
	const int D      = 411;
	const int Dstar  = 413;
	const int Dstar0 = 423;
	const int Ds     = 431;
	const int Dsstar = 433;
	const int P      = 2212;
	const int N      = 2112;
	const int Phi    = 333;
	const int B      = 521;
	const int B0     = 511;
	const int Bs     = 531;
	const int Bc     = 541;
	const int K0     = 311;
	const int KS     = 310;
	const int KL     = 130;
	const int L0     = 3122;
	const int Lc     = 4122;
	const int Lb     = 5122;
}

