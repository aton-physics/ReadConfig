#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <random>	// need uniform distribution to select random wavevectors for the incoherent intermediate scattering function
#include <algorithm> // std::upper_bound to find the value for tau_orient 
#include <iterator> // istream_iterator
#include <iomanip> // std::setprecision
#include <numeric> // std::accumulate
#include <sys/stat.h> // mkdir
#include "headers/Point.h" // Point class, some storage structs
#include "headers/ParameterClass.h" // InputParameter class reads my input files, SimulationParameters has a few extra that ReadConfig needs, that aren't part of input
#include "headers/Numerical.h"

//***********************************8/6/2018**********************************
//This program takes in a potentially gigantic trajectory fromqs molecular dynamics (gigabytes) and calculates number of nearest neighbors and all the correlation functions
//Meant to be called as a subprocess through a shell driver that also grabs task_id from the Sun Grid Engine. something like ./WriteTrajectoryToSTDOUT | xz > CompressTrajectory.xz, then unxz
//****1/7/19****
//Make sure SimulationParameters are the same between programs that I am using. I should put this in a header, but I haven't.

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model);
double Cos1ToCos6(double cos1x);	//take in cos(1x), make it cos(6x)
void WriteConfiguration(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::ostream &strm);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
InputParameter ReadInput(std::string &tag);
SimulationParameters GetParams(InputParameter Input);
std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag);
void FindRelaxationTime(SimulationParameters &model, std::string filename, std::vector<double> OrientCorrelation);
std::vector<std::vector<int>> create_list(SimulationParameters &model, std::vector<Point> &Configuration, double length);
void CalculateCorrelationFunctionsOrderN(SimulationParameters &model, int &NumberOfTrajectories, int &t_cor, int &NumberOfCoarseSteps, std::vector<std::vector<Point>> &OmegaJ,
std::vector<std::vector<Point>> &Cm_J, std::vector<double> &OrientCorrelation, std::vector<double> &Orient6Correlation, std::vector<double> &CollectiveOrient6Correlation,
std::vector<double> &MeanSqDisplacement, std::vector<double> &ScatteringFunction, std::vector<std::string> &OutputFiles, std::string tag);

//things to change manually for different jobs: t_cor and steps_between_cfgs
int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);
	//std::vector<std::string> OutputFolders = { "Orientation2Correlation", "msd", "Diffusion", "Tau_orient", "Orientation6Correlation", "CollectiveOrientation6Correlation", "Scattering" };
	std::vector<std::string> OutputFolders = { "CRMSD" };
	std::vector<std::string> OutputFiles = FoldersToFiles(model, OutputFolders, tag);
	//CHANGE T_COR WHEN SUBMITTING JOBS OF DIFFERENT LENGTHS
	int t_cor = int(pow(10,2) / model.timestep / model.steps_between_cfgs); // correlation function is measured for this long (in units of model.timestep* #skipped configurations)
	assert(model.num_cfgs >= t_cor); //make sure there at least as many configurations as correlation steps
	int NumberOfCoarseSteps = 3250; //2500(1) + 750(10) to get to the 10^4th index (spacing is 10 tau)
	std::vector<std::vector<Point>> OmegaJ(model.num_cfgs, std::vector<Point>(model.N)), Cm_J(model.num_cfgs, std::vector<Point>(model.N));
	double OrientationMagnitude = 0.0;
	for (int n = 0; n < model.num_cfgs; n++) {
		//if (n % 1000 == 0) std::cout << n << '\n';
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);
		if (n == 0) OrientationMagnitude = sqrt(Positions[0][2].dist_sq(Positions[0][1]));	//grab orientation vector's magnitude, but only once.
		for (int i = 0; i < model.N; i++) {	// Store orientations, center of mass of each molecule for every configuration to calculate dynamical correlation functions.
			OmegaJ[n][i] = (Positions[i][2] - Positions[i][1]) / OrientationMagnitude;
			Cm_J[n][i] = (Positions[i][2] + Positions[i][1] + Positions[i][0]) / 3.0;
		}
	}
	int NumberOfTrajectories = model.num_cfgs - t_cor + 1; // every new configuration constitutes new initial conditions when calculating correlation functions.
	//int NumberOfCoarseSteps = 12390; // 2500(1 + 10 + 50 + 100) + 2390(250) to get to the 10^6'th index = ~9999.x tau = t_cor 
	//int NumberOfCoarseSteps = 6450; //2500(1 + 10) + 1450(50) to get to the 10^5th index
	std::vector<double> OrientCorrelation(NumberOfCoarseSteps), Orient6Correlation(NumberOfCoarseSteps), 
		CollectiveOrient6Correlation(NumberOfCoarseSteps), MeanSqDisplacement(NumberOfCoarseSteps), ScatteringFunction(NumberOfCoarseSteps);
	CalculateCorrelationFunctionsOrderN(model, NumberOfTrajectories, t_cor, NumberOfCoarseSteps, OmegaJ, Cm_J, 
		OrientCorrelation, Orient6Correlation, CollectiveOrient6Correlation, MeanSqDisplacement, ScatteringFunction,
		OutputFiles, tag);
	//if I calculate derivatives and diffusion, remember to alter OutputFiles[] and its calls
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
	//tag = "241";	// commented out for production runs. comment in when testing a trajectory.
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
	//IMPORTANT: make sure this matches with any other programs you use.
	//int num_mol = Input.N, num_atom = 3, numcfg = Input.NumConfigs, steps_between_cfgs = 10000;	
	int num_mol = Input.N, num_atom = 3, numcfg = Input.NumConfigs, steps_between_cfgs = 10;
	double boxl = sqrt(num_mol / density);
	double boxinv = 1.0 / boxl;
	SimulationParameters model(density, temperature, boxl, boxinv, ts, numcfg, num_mol, num_atom, steps_between_cfgs);
	return model;
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

std::vector<std::vector<int>> create_list(SimulationParameters &model, std::vector<Point> &Configuration, double length) {
	double lengthsq = length*length;
	std::vector<std::vector<int>> NeighborList(model.N);
	for (int i = 0; i < model.N; i++) {
		for (int j = 0; j < model.N; j++) {
			Point separation = Configuration[i] - Configuration[j];
			separation = separation.pbc(model.boxl, model.invboxl);
			double sepsq = pow(separation.x(), 2) + pow(separation.y(), 2);
			if (sepsq < lengthsq) NeighborList[i].push_back(j);
		}
	}
	return NeighborList;
}

void CalculateCorrelationFunctionsOrderN(SimulationParameters &model, int &NumberOfTrajectories, int &t_cor, int &NumberOfCoarseSteps, std::vector<std::vector<Point>> &OmegaJ,
	std::vector<std::vector<Point>> &Cm_J, std::vector<double> &OrientCorrelation, std::vector<double> &Orient6Correlation, std::vector<double> &CollectiveOrient6Correlation,
	std::vector<double> &MeanSqDisplacement, std::vector<double> &ScatteringFunction, std::vector<std::string> &OutputFiles, std::string tag) {
	//order-n method. Correlates from 0 to t_cor averaged over numberoftrajectories. Each time step is treated as a new origin. Coarser resolution correlation functions every 2500 steps
	int blocksize = 2500; // number of words per block
	std::vector<int> CoarseIndex = { 1, 10, 50, 100, 250 };	// calculate correlation functions up to (a maximum of) 2500*(1+10+50+100+250) steps
	std::vector<int> CoarseOffset;	// offset needed to index the data since we're skipping around
	for (std::vector<int>::size_type i = 0; i < CoarseIndex.size(); i++) {
		if (i == 0) {
			CoarseOffset.push_back(0);
		}
		else CoarseOffset.push_back(CoarseOffset[i - 1] + CoarseIndex[i - 1] * blocksize);
	}
	std::vector<std::vector<double>> Orient2(int(CoarseIndex.size()), std::vector<double>(blocksize)), 
		MSD(int(CoarseIndex.size()), std::vector<double>(blocksize));
	std::vector<std::vector<double>> Orient6(int(CoarseIndex.size()), std::vector<double>(blocksize)),
		CollectiveOrient6(int(CoarseIndex.size()), std::vector<double>(blocksize)), 
		Scattering(int(CoarseIndex.size()), std::vector<double>(blocksize));
	double NeighborCutOff = 2.23;	// 7 is very rare, only a few molecules per configuration count 7 neighbors for 2.23sigma
	/*int NumScatteringVectors = 10;
	double ScatteringMagnitude = 3.5;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> RealDist(-ScatteringMagnitude, ScatteringMagnitude);
	double qx[NumScatteringVectors], qy[NumScatteringVectors];
	for (int m = 0; m < NumScatteringVectors; m++) {
		qx[m] = RealDist(generator);
		qy[m] = (2 * (m % 2) - 1) * sqrt(ScatteringMagnitude*ScatteringMagnitude - qx[m]*qx[m]); //select qy based on magnitude of qx, and assign qy positive or negative values depending on the loop variable
	}*/
	for (int j = 0; j < NumberOfTrajectories; j++) {
		//if (j % 100 == 0) std::cout << "Trajectory # " << j << '\n';	//28 seconds as of now
		for (std::vector<int>::size_type p = 0; p < CoarseIndex.size(); p++) {
			for (int k = 0; k < blocksize; k++) {
				if (CoarseOffset[p] + (k*CoarseIndex[p]) >= t_cor) continue;	// t_cor is a hard limit to how far I want to calculate my correlation functions
				double total_omega2 = 0.0, total_omega6 = 0.0, total_collective_omega6 = 0.0, total_scattering = 0.0, total_cm = 0.0;
				std::vector<std::vector<int>> NeighborList = create_list(model, Cm_J[j + CoarseOffset[p] + (k*CoarseIndex[p])], NeighborCutOff);	// pass to create_list a whole configuration
				for (int i = 0; i < model.N; i++) {	//CoarseOffset[3] = 90,000. j = 190,000. k can be at most 399 since max index of OmegaJ is 289,999
					//total_omega6 += Cos1ToCos6(OmegaJ[j][i].dot(OmegaJ[j + CoarseOffset[p] + (k*CoarseIndex[p])][i]));
					//total_omega2 += pow(OmegaJ[j][i].dot(OmegaJ[j + CoarseOffset[p] + (k*CoarseIndex[p])][i]), 2);
					Point CR_sum = 0;
					for (auto m : NeighborList[i]) {
						CR_sum += Cm_J[j][m] - Cm_J[j + CoarseOffset[p] + (k*CoarseIndex[p])][m];
					}
					CR_sum = CR_sum / NeighborList[i].size();
					Point cm_temp = Cm_J[j][i] - Cm_J[j + CoarseOffset[p] + (k*CoarseIndex[p])][i];
					cm_temp -= CR_sum;
					total_cm += pow(cm_temp.x(), 2) + pow(cm_temp.y(), 2);
					/*for (int m = 0; m < NumScatteringVectors; m++) {		
						total_scattering += cos(qx[m]*cm_temp.x() - qy[m]*cm_temp.y());	// dot product of q and cm(t) - cm(0)
					}*/
				}
				/*for (int i = 0; i < model.N; i++) {	//higher order correlation - 2-time, 2-particle
					for (int q = 0; q < model.N; q++) {
						total_collective_omega6 += Cos1ToCos6(OmegaJ[j][i].dot(OmegaJ[j + CoarseOffset[p] + (k*CoarseIndex[p])][q]));
					}
				}*/
				//CollectiveOrient6[p][k] += total_collective_omega6;	//sorted by p because I need to later sort the correlation function's parts by block number to weight them correctly
				//Orient6[p][k] += total_omega6;
				//Orient2[p][k] += total_omega2;
				MSD[p][k] += total_cm;
				//Scattering[p][k] += total_scattering;
			}
		}
	}
	for (int i = 0; i < NumberOfCoarseSteps; i++) {	// sort the correlation function's parts by block number and add them appropriately
			int block_number = static_cast<int>(i / blocksize);	// if i is less than blocksize, then it is 0th block, etc. Cast to int truncates to the proper block
			//CollectiveOrient6Correlation[i] = CollectiveOrient6[block_number][i%blocksize];
			//Orient6Correlation[i] = Orient6[block_number][i%blocksize];
			//OrientCorrelation[i] = Orient2[block_number][i%blocksize];
			MeanSqDisplacement[i] = MSD[block_number][i%blocksize];
			//ScatteringFunction[i] = Scattering[block_number][i%blocksize];
	}
	//std::ofstream orientfile(OutputFiles[0], std::ios_base::app);
	//std::ofstream msdfile(OutputFiles[1], std::ios_base::app);
	//std::ofstream orient6file(OutputFiles[4], std::ios_base::app);
	//std::ofstream collectiveorient6file(OutputFiles[5], std::ios_base::app);
	//std::ofstream scatteringfile(OutputFiles[6], std::ios_base::app);
	std::ofstream msdfile(OutputFiles[0], std::ios_base::app);
	for (int k = 0; k < NumberOfCoarseSteps; k++) {
		double counter = 0;
		//CollectiveOrient6Correlation[k] /= (pow(double(model.N),2) * NumberOfTrajectories);
		//Orient6Correlation[k] /= (double(model.N)* NumberOfTrajectories);
		//OrientCorrelation[k] = 2.0 * OrientCorrelation[k] / double(model.N) / NumberOfTrajectories - 1;	// divide here so I can use std::upper_bound later
		MeanSqDisplacement[k] /= (double(model.N) * double(NumberOfTrajectories));
		//ScatteringFunction[k] /= (double(model.N) * double(NumberOfTrajectories) * NumScatteringVectors);
		if (k < 2500) counter = k;
		else if (k < 5000) counter = CoarseIndex[0] * blocksize + CoarseIndex[1] * (k % blocksize);	// counter is a linear way to index the piecewise way we are tracking time.
		else if (k < 7500) counter = (CoarseIndex[0] + CoarseIndex[1]) * blocksize + CoarseIndex[2] * (k % blocksize);
		else if (k < 10000) counter = (CoarseIndex[0] + CoarseIndex[1] + CoarseIndex[2]) * blocksize + CoarseIndex[3] * (k % blocksize);
		else counter = (CoarseIndex[0] + CoarseIndex[1] + CoarseIndex[2] + CoarseIndex[3]) * blocksize + CoarseIndex[4] * (k % blocksize);
		double what_time = counter * model.timestep* model.steps_between_cfgs;
		//orientfile << what_time << '\t' << OrientCorrelation[k] << '\n';
		msdfile << what_time << '\t' << std::setprecision(12) << MeanSqDisplacement[k] << '\n';	//subtraction of nearly equal values
		//orient6file << what_time << '\t' << Orient6Correlation[k] << '\n';
		//collectiveorient6file << what_time << '\t' << CollectiveOrient6Correlation[k] << '\n';
		//scatteringfile << what_time << '\t' << ScatteringFunction[k] << '\n';
	}
}

double Cos1ToCos6(double cos1x) {	//take in cos(1x), make it cos(6x)
	double c = cos1x;
	double c2 = c*c;
	double c4 = c2*c2;
	double cos6x = 4 * (8 * c4*c2 - 12 * c4 + 6 * c2 - 1) - 6 * c2 + 3;
	return cos6x;
}
