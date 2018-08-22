#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm> // std::upper_bound to find the value for tau_orient 
#include <iterator> // istream_iterator
#include <iomanip> // std::setprecision
#include <numeric> // std::accumulate
#include <sys/stat.h> // mkdir
#include "headers/Point.h" // Point class, some storage structs
#include "headers/ParameterClass.h" // InputParameter class reads my input files, SimulationParameters has a few extra that ReadConfig needs, that aren't part of input
#include "headers/Numerical.h"

//***********************************8/6/2018**********************************
//This program takes in a potentially gigantic trajectory from molecular dynamics (gigabytes) and calculates number of nearest neighbors and all the correlation functions
//Meant to be called as a subprocess through a shell driver that also grabs task_id from the Sun Grid Engine. something like ./WriteTrajectoryToSTDOUT | xz > CompressTrajectory.xz, then unxz

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model);
void WriteConfiguration(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::ostream &strm);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
InputParameter ReadInput(std::string &tag);
SimulationParameters GetParams(InputParameter Input);
void NearestNeighbors(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::vector<int> &NeighborVector, double &NeighborCutOffSq, std::string tag);
std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag);
void FindRelaxationTime(SimulationParameters &model, std::string filename, std::vector<double> OrientCorrelation);
double PairCorrelationFirstMin(SimulationParameters &model, HistogramInfo &hist, std::string tag);
void PrintNeighbors(SimulationParameters &model, std::string filename, std::vector<int> &NumberNeighbors);
void CalculateCorrelationFunctionsWindow(SimulationParameters &model, int &NumberOfTrajectories, int &t_cor, std::vector<std::vector<Point>> &OmegaJ, std::vector<std::vector<Point>> &Cm_J,
	std::vector<double> &OrientCorrelation, std::vector<double> &MeanSqDisplacement);
void CalculateCorrelationFunctionsOrderN(SimulationParameters &model, int &NumberOfTrajectories, int &t_cor, int &NumberOfCoarseSteps, std::vector<std::vector<Point>> &OmegaJ,
	std::vector<std::vector<Point>> &Cm_J, std::vector<double> &OrientCorrelation, std::vector<double> &MeanSqDisplacement, std::vector<std::string> &OutputFiles, double &OrientationMagnitude, std::string tag);
//this function needs to be refactored

int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);
	std::vector<std::string> OutputFolders = { "OrientCorrelation", "msd", "Diffusion", "Tau_orient", "Neighbors" };
	std::vector<std::string> OutputFiles = FoldersToFiles(model, OutputFolders, tag);
	HistogramInfo HistRDF(0.01, int(model.boxl / 2 / 0.01));
	int t_cor = int(1000 / model.timestep / model.steps_between_cfgs); // correlation function is measured for this long (in units of model.timestep* #skipped configurations)
	assert(model.num_cfgs >= t_cor); //make sure there at least as many configurations as correlation steps
	std::vector<std::vector<Point>> OmegaJ(model.num_cfgs, std::vector<Point>(model.N)), Cm_J(model.num_cfgs, std::vector<Point>(model.N));
	double OrientationMagnitude = 0.0;
	double NeighborCutOffSq = pow(PairCorrelationFirstMin(model, HistRDF, tag), 2);	// find the first minimum of g(r).
	std::vector<int> NumberNeighbors(100);
	for (int n = 0; n < model.num_cfgs; n++) {
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);
		if (n == 0) OrientationMagnitude = Positions[0][2].dist_sq(Positions[0][1]);	//grab orientation vector's magnitude, but only once.
		for (int i = 0; i < model.N; i++) {	// Store orientations, center of mass of each molecule for every configuration to calculate C(t), MSD(t) respectively.
			OmegaJ[n][i] = Positions[i][2] - Positions[i][1];
			Cm_J[n][i] = (Positions[i][2] + Positions[i][1] + Positions[i][0]) / 3.0;
		}
		NearestNeighbors(model, Positions, NumberNeighbors, NeighborCutOffSq, tag);
	}
	PrintNeighbors(model, OutputFiles[4], NumberNeighbors);
	int NumberOfTrajectories = model.num_cfgs - t_cor + 1; // every new configuration constitutes new initial conditions when calculating correlation functions.
	int NumberOfCoarseSteps = 7900; // 1 * 2500 + 10 * 2500 + 25 * 2500 + 25 * 400 to get to the 10^5'th index = ~999.x tau = t_cor 
	std::vector<double> OrientCorrelation(NumberOfCoarseSteps), MeanSqDisplacement(NumberOfCoarseSteps);
	CalculateCorrelationFunctionsOrderN(model, NumberOfTrajectories, t_cor, NumberOfCoarseSteps, OmegaJ, Cm_J, OrientCorrelation, MeanSqDisplacement, OutputFiles, OrientationMagnitude, tag);
	FindRelaxationTime(model, OutputFiles[3], OrientCorrelation);
	numerical_differentiation(OutputFiles[1], "derivative", model.steps_between_cfgs*model.timestep, tag); // derivative/derivative is a temporary solution to see how long I should wait before taking the zero slope regression line
																										   //best way to do this is to not print anything but diffusion (maybe msd). send msd to numerical_diff to zero_slope > diffusion.data, don't ever want to look at derivative
	zero_slope_regression("derivative", OutputFiles[2], model.temperature, 30.0, tag);	// skip 10 tau before taking the regression
	return 0;
}

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model) {	// read in a configuration.
	std::vector<std::vector<Point>> PositionVector(model.N, std::vector<Point>(model.NA));
	double f, g;
	for (int i = 0; i < model.N; i++) {
		for (int a = 0; a < model.NA; a++) {
			std::cin >> f >> g;
			Point temp(f, g);
			PositionVector[i][a] = temp;
		}
	}
	return PositionVector;
}

void WriteConfiguration(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::ostream &strm) {	// write a configuration
	for (int i = 0; i < model.N; i++) {
		for (int a = 0; a < model.NA; a++) {
			PositionVector[i][a].print(strm);
		}
	}
}

std::ifstream& GotoLine(std::ifstream& file, unsigned int num) {//skip to a certain line in data file
	file.seekg(std::ios::beg);
	for (unsigned int i = 0; i < num - 1; ++i) {
		file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
	return file;
}

InputParameter ReadInput(std::string &tag) { // construct with all the input parameters
	char* taskID_string;
	taskID_string = getenv("SGE_TASK_ID");
	std::string task_id = taskID_string;
	tag = task_id;
	std::ifstream myfile("inputfiles/input.data"); // read parameters from this file
												   //tag = "1";
	GotoLine(myfile, std::stoi(tag));	//skip to tagged line
	int linenum, numcfg, NumMol;
	double bondangle, density, temperature;
	std::string MeltedCfg;
	myfile >> linenum >> bondangle >> density >> temperature >> numcfg >> NumMol >> MeltedCfg;
	InputParameter parameters(bondangle, density, temperature, numcfg, NumMol, MeltedCfg);
	return parameters;
}

SimulationParameters GetParams(InputParameter Input) { //read in input file, assign all the parameters.
	double density = Input.density, temperature = Input.temperature, ts = 0.001;
	int num_mol = Input.N, num_atom = 3, numcfg = Input.NumConfigs;
	double boxl = sqrt(num_mol / density);
	double boxinv = 1.0 / boxl;
	SimulationParameters model(density, temperature, boxl, boxinv, ts, numcfg, num_mol, num_atom);
	return model;
}

double PairCorrelationFirstMin(SimulationParameters &model, HistogramInfo &hist, std::string tag) {	// defines nearest neighbor cutoff distance
	std::ifstream rdfstream("PairCorrelation/" + tag + ".data");
	std::string dummyLine;
	getline(rdfstream, dummyLine);	// throw away header
	double a, b;
	std::vector<double> data;
	while (rdfstream >> a >> b) {
		data.push_back(b);	// grab elements of 2nd field, which is g(r).
	}
	if (data.size() < 100) {
		std::cout << "Tag is " << tag << "g(r) vector was too short \n";
		return 1.61;
	}
	for (int bin = 1; bin < hist.num_bins; bin++) {	// skip first bin to make range checking happy
		if (double(bin) * hist.bin_width < 1) continue;	// first minimum doesn't come before 1 sigma.
		if (data[bin] - data[bin - 1] < 0.0 && data[bin + 1] - data[bin] > 0.0) {	// if bin is less than the previous and less than the next, then it is a minimum.
			if (data[bin] < 1.0) {	// We want the first minimum that is also a trough. If a minimum represents anticorrelation (g(r) < 1), then it is likely a trough. Otherwise it is a result of molecular geometry.
				double binvalue = double(bin) * hist.bin_width;
				return binvalue;
			}
		}
	}
	double minimum = 1.61; // if the code doesn't find a trough, set this as default minimum
	std::cout << "Tag is " << tag << "found no minimum for g(r) \n";
	return minimum;
}

void NearestNeighbors(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::vector<int> &NeighborVector, double &NeighborCutOffSq, std::string tag) {
	for (int j = 0; j < model.N; j++) {				// loop over all distinct molecule pairs j,k - could presumably do j,k with k > j but I'm not sure of a neat way to count neighbors without going over every index individually, can't double count like we do in g(r).
		for (int a = 0; a < model.NA; a++) {
			int counter = 0;					// reset counter when examining new reference particle
			for (int k = 0; k < model.N; k++) {
				if (j == k) continue;			// skip neighbor check if loop is over the same molecule twice <=> make sure molecules are distinct
				for (int b = 0; b < model.NA; b++) {
					Point r_jakb = (PositionVector[j][a] - PositionVector[k][b]);
					r_jakb = r_jakb.pbc(model.boxl, model.invboxl);
					if (r_jakb.dot(r_jakb) <= NeighborCutOffSq) {	// check if close enough to be a neighbor
						counter += 1;
					}
				}
			}
			//assert(counter <= 100);	//6 neighbors for molecules with non-overlapping disks. 
			if (counter >= 25) std::cout << "counter was greater than 25, I am tag " << tag << " with cut off distance "
				<< NeighborCutOffSq << " here are my model parameters: " << model.num_cfgs << '\t' << 
				model.boxl << '\n';
			assert(counter <= 600);
			NeighborVector[counter] += 1; // If InRange, then register another near neighbor. 
		}
	}
}

std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag) {
	// Make the specified Folders, then turn them into their respective Folder/$task_id.data, then print a header (density, temperature) for plotting ease
	std::vector<std::string> vectorstring;
	for (auto i : stringvector) {
		mkdir(i.c_str(), ACCESSPERMS);
		i = i + "/" + tag + ".data";
		vectorstring.push_back(i);
		std::ofstream ofs(i);
		ofs << model.density << '\t' << model.temperature << '\n';
	}
	return vectorstring;
}

void FindRelaxationTime(SimulationParameters &model, std::string filename, std::vector<double> OrientCorrelation) {	// Want to temporarily transform OrientCorrelation, so pass by value.
	for (auto &i : OrientCorrelation) {
		i *= -1;	// flip the signature on the values of OrientCorrelation since it's monotonically decreasing and I want to use std::upper_bound
	}
	std::vector<double>::iterator it = std::lower_bound(OrientCorrelation.begin(), OrientCorrelation.end(), -0.367879);	// find first value not less than negative inverse euler's number
	auto index = it - OrientCorrelation.begin();
	if (index == int(OrientCorrelation.size())) index = 9999999;
	std::ofstream tau_orientfile(filename, std::ios_base::app);
	tau_orientfile << model.temperature << '\t' << index * model.timestep* model.steps_between_cfgs << '\n';
}

void PrintNeighbors(SimulationParameters &model, std::string filename, std::vector<int> &NumberNeighbors) {
	std::ofstream neighborfile(filename, std::ios_base::app);
	for (std::vector<int>::size_type i = 0; i < NumberNeighbors.size(); i++) {
		neighborfile << i << "\t" << NumberNeighbors[i] / double(model.N * model.NA) / model.num_cfgs << '\n';
	}
}

void CalculateCorrelationFunctionsOrderN(SimulationParameters &model, int &NumberOfTrajectories, int &t_cor, int &NumberOfCoarseSteps, std::vector<std::vector<Point>> &OmegaJ, 
	std::vector<std::vector<Point>> &Cm_J, std::vector<double> &OrientCorrelation, std::vector<double> &MeanSqDisplacement, std::vector<std::string> &OutputFiles, double &OrientationMagnitude, std::string tag) {
	//order-n method. Correlates from 0 to t_cor averaged over numberoftrajectories. Each time step is treated as a new origin. Coarser grained correlation functions every 2500 steps
	int blocksize = 2500; // number of words per block
	std::vector<int> CoarseIndex = { 1, 10, 25, 25 };	// calculate correlation functions up to (a maximum of) 1*2500 + 10*2500 + 25*2500 + 25*2500 steps
	std::vector<int> CoarseOffset;	// offset needed to index the data since we're skipping around
	for (std::vector<int>::size_type i = 0; i < CoarseIndex.size(); i++) {
		if (i == 0) {
			CoarseOffset.push_back(0);
		}
		else CoarseOffset.push_back(CoarseOffset[i - 1] + CoarseIndex[i - 1] * blocksize);
	}
	std::vector<std::vector<double>> Orient(4, std::vector<double>(blocksize)), MSD(4, std::vector<double>(blocksize));
	for (int j = 0; j < NumberOfTrajectories; j++) {
		if (j % 1000 == 0) {
			std::ofstream ofs("Track/" + tag + ".data");
			ofs << "Trajectory number " << j << '\n';
		}
		for (std::vector<int>::size_type p = 0; p < CoarseIndex.size(); p++) {
			for (int k = 0; k < blocksize; k++) {
				if (CoarseOffset[p] + (k*CoarseIndex[p]) >= t_cor) continue;	// if past t_cor, stop. (k >= 400)
				double total_omega = 0.0, total_cm = 0.0;
				for (int i = 0; i < model.N; i++) {	//CoarseOffset[3] = 90,000. j = 190,000. k can be at most 399 since max index of OmegaJ is 289,999
					total_omega += pow(OmegaJ[j][i].dot(OmegaJ[j + CoarseOffset[p] + (k*CoarseIndex[p])][i]), 2);
					Point cm_temp = Cm_J[j][i] - Cm_J[j + CoarseOffset[p] + (k*CoarseIndex[p])][i];
					total_cm += pow(cm_temp.x(), 2) + pow(cm_temp.y(), 2);
				}
				Orient[p][k] += total_omega;
				MSD[p][k] += total_cm;
			}
		}
	}
	for (int i = 0; i < NumberOfCoarseSteps; i++) {	// sort the correlation function's parts by block number and add them appropriately
			int block_number = static_cast<int>(i / blocksize);	// if i is less than blocksize, then it is 0th block, etc. Cast to int truncates to the proper block
			OrientCorrelation[i] = Orient[block_number][i%blocksize];
			MeanSqDisplacement[i] = MSD[block_number][i%blocksize];
	}
	std::ofstream orientfile(OutputFiles[0], std::ios_base::app);
	std::ofstream msdfile(OutputFiles[1], std::ios_base::app);
	for (int k = 0; k < NumberOfCoarseSteps; k++) {
		double counter = 0;
		OrientCorrelation[k] = 2.0 * OrientCorrelation[k] / pow(OrientationMagnitude, 2) / double(model.N) / NumberOfTrajectories - 1;	// divide here so I can use std::upper_bound later
		MeanSqDisplacement[k] /= (double(model.N) * double(NumberOfTrajectories));
		if (k < 2500) counter = k;
		else if (k < 5000) counter = 1 * 2500 + 10 * (k % 2500);
		else if (k < 7500) counter = 1 * 2500 + 10 * 2500 + 25 * (k % 2500);
		else counter = 1 * 2500 + 10 * 2500 + 25 * 2500 + 25 * (k % 2500);
		double what_time = counter * model.timestep* model.steps_between_cfgs;
		orientfile << what_time << '\t' << OrientCorrelation[k] << '\n';
		msdfile << what_time << '\t' << std::setprecision(12) << MeanSqDisplacement[k] << '\n';	//numerical imprecision in numerical differentiation (subtraction of nearly equal values)
	}
}

void CalculateCorrelationFunctionsWindow(SimulationParameters &model, int &NumberOfTrajectories, int &t_cor, std::vector<std::vector<Point>> &OmegaJ, std::vector<std::vector<Point>> &Cm_J,
	std::vector<double> &OrientCorrelation, std::vector<double> &MeanSqDisplacement) {
	//direct windowed method. Correlates from 0 to t_cor averaged over numberoftrajectories. Each time step is treated as a new origin
	//currently defunct, don't use. Haven't added functionality to print OrientCorrelation, MeanSqDisplacement (but it's quite easy to figure out since step size is uniform)
	for (int j = 0; j < NumberOfTrajectories; j++) {
		/*if (j % 1000 == 0) {
			std::ofstream ofs("Track/" + tag + ".data");
			ofs << "Trajectory number " << j << '\n';
		}*/
		for (int k = 0; k < t_cor; k++) {
			double total_omega = 0.0, total_cm = 0.0;
			for (int i = 0; i < model.N; i++) {
				total_omega += pow(OmegaJ[j][i].dot(OmegaJ[j + k][i]), 2);
				Point cm_temp = Cm_J[j][i] - Cm_J[j + k][i];
				total_cm += pow(cm_temp.x(), 2) + pow(cm_temp.y(), 2);
			}
			OrientCorrelation[k] += total_omega;
			MeanSqDisplacement[k] += total_cm;
		}
	}
}
