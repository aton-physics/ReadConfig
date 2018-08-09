#pragma once
#include <vector>
#include <string>

struct InputParameter {
	double BondAngle, density, temperature;
	int NumConfigs, N;	// note - I used global arrays, so reassigning N would be a massive pain. This struct will store N, but not use it ever. 
	std::string MeltCfg;
	InputParameter (double b, double d, double t, int nc, int n, std::string m)	// constructor
	{
		BondAngle = b;
		density = d;
		temperature = t;
		NumConfigs = nc;
		N = n;
		MeltCfg = m;
	}
	void print(std::ostream &strm)
	{
		strm << BondAngle << '\t' << density << '\t' << temperature << '\t' << NumConfigs << '\t' << N << '\t' << MeltCfg << '\n';
	}
};

struct BondConstraints {
	double L1 = 1.0, L2 = 1.0;					// optional
	double L3;
	BondConstraints(double l1, double l2, double bond_angle) {
		bond_angle *= 0.0174532925199;	// degrees to radians conversion factor
		L1 = l1;
		L2 = l2;
		L3 = L1 * sin(bond_angle / 2) + L2 * sin(bond_angle / 2);
	}
	std::vector<double> dsq() {	// calculate the squared bond contraints for RATTLE
		std::vector<double> temp;
		temp.push_back(L1 * L1);
		temp.push_back(L3 * L3);
		temp.push_back(L2 * L2);
		return temp;
	}
	void print(std::ostream &strm)
	{
		strm << L1 << '\t' << L2 << '\t' << L3 << '\n';
	}
};