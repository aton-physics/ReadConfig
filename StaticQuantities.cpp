#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm> // std::upper_bound to find the value for tau_orient 
#include <iterator> // istream_iterator
#include <iomanip> // std::setprecision
#include <sys/stat.h> // mkdir
#include "headers/Point.h" // Point class, some storage structs
#include "headers/ParameterClass.h" // InputParameter class reads my input files, SimulationParameters has a few extra that ReadConfig needs, that aren't part of input

//***********************************8/6/2018**********************************
//This program takes in a potentially gigantic trajectory from molecular dynamics (gigabytes) and calculates various static quantities like the distribution of the order parameter, pair correlation function, also does video microscopy
//Meant to be called as a subprocess through a shell driver that also grabs task_id from the Sun Grid Engine. something like ./WriteTrajectoryToSTDOUT | xz > CompressTrajectory.xz, then unxz

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model);
void WriteConfiguration(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::ostream &strm);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
InputParameter ReadInput(std::string &tag);
SimulationParameters GetParams(InputParameter Input);
double CalculateS(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<Point>> &PositionVector);
void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode);
void radial_df(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf);
std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag);
void VideoMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<Point>> &PositionVector, std::string tag);

int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);
	std::vector<std::string> OutputFolders = {"OrderParameter", "PairCorrelation"};
	mkdir("Microscopy", ACCESSPERMS);
	std::vector<std::string> OutputFiles = FoldersToFiles(model, OutputFolders, tag);
	HistogramInfo HistOrderParameter(0.005, int(1 / 0.005)), HistRDF(0.01, int(model.boxl / 2 / 0.01));
	std::vector<double> OrderParameter(HistOrderParameter.num_bins), PairCorrelation(HistRDF.num_bins);
	double OrientationMagnitude = 0.0;
	for (int n = 0; n < model.num_cfgs; n++) {
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);
		WriteConfiguration(model, Positions, std::cout);
		if (n == 0) OrientationMagnitude = Positions[0][2].dist_sq(Positions[0][1]);	//grab orientation vector's magnitude, but only once.
		int bin = int(CalculateS(model, OrientationMagnitude, Positions) / HistOrderParameter.bin_width);
		OrderParameter[bin] += 1;
		radial_df(model, Positions, HistRDF, PairCorrelation);
		VideoMicroscopy(model, "Microscopy", Positions, tag);
	}
	PrintHistogram(model, HistOrderParameter, OutputFiles[0], OrderParameter, "prob_density");
	PrintHistogram(model, HistRDF, OutputFiles[1], PairCorrelation, "histogram");
	return 0;
}

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model) {	// read in a configuration.
	std::vector<std::vector<Point>> PositionVector(model.N, std::vector<Point> (model.NA));
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

double CalculateS(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<Point>> &PositionVector) {	//take a configuration and calculate the "orientational order parameter" see Zhou, Stratt 2018
	double s = 0, Q1 = 0, Q2 = 0;
	for (int i = 0; i < model.N; i++) {
		Point OmegaJ = PositionVector[i][2] - PositionVector[i][1];
		//Q1 += 2 * pow(OmegaJ.x(), 2) - 1;
		//Q2 += 2 * OmegaJ.x() * OmegaJ.y();
		Q1 += pow(OmegaJ.x(), 2);
		Q2 += OmegaJ.x() * OmegaJ.y();
	}
	Q1 = 2 * Q1 / NormalizationSq - model.N;	// each of OmegaJ.x and OmegaJ.y carry a factor of 1/Normalization, so their product has a factor of 1/NormalizationSq.
	Q2 *= 2 / NormalizationSq;
	s = sqrt(Q1 * Q1 + Q2 * Q2) / model.N;
	return s;
}

void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode) {
	std::ofstream ofs(filename, std::ios_base::app);
	for (int bin = 0; bin < hist.num_bins; bin++) {
		double binvalue = double(bin) * hist.bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / model.num_cfgs << '\n';
		else if (mode == "prob_density") ofs << binvalue << '\t' << data[bin] / model.num_cfgs / hist.bin_width << '\n'; // if optional parameter is specified and greater than 1, divide bin count by bin width. This turns a histogram into a probability density
		else std::cout << "Bad parameter for PrintHistogram()";
	}
}

void radial_df(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf) {	// atom-atom g(r). function takes in Positions and fills out a vector of doubles
	const double constant = 3.14159265 * model.density * model.NA;
	double rlower, rupper, nideal, rij;
	std::vector<double> histo(hist.num_bins);
	int bin;
	for (int j = 0; j < model.N - 1; j++) {				//loop over all distinct molecule pairs j,k
		for (int k = j + 1; k < model.N; k++) {
			for (int a = 0; a < model.NA; a++) {
				for (int b = 0; b < model.NA; b++) {
					Point r_jakb = PositionVector[j][a] - PositionVector[k][b]; //calculate pair separation
					r_jakb = r_jakb.pbc(model.boxl, model.invboxl);	// fold the pair separation to get the minimum image separation
					rij = sqrt(r_jakb.dot(r_jakb));
					bin = int(rij / hist.bin_width);
					if (bin < hist.num_bins) {	// make sure we only count particles out as far as some maximum (set to be half the box length)
						histo[bin] += 2;		// if rij falls within the range of the bin, count two particles (i and j)
					}
				}
			}
		}
	}
	for (bin = 0; bin < hist.num_bins; bin++) {
		rlower = bin * hist.bin_width;
		rupper = rlower + hist.bin_width;
		nideal = constant * (rupper * rupper - rlower * rlower);
		rdf[bin] += double(histo[bin]) / (model.N * model.NA) / nideal;
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

void VideoMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<Point>> &PositionVector, std::string tag) {	//print configuration with minimum images and molecule labels
	filename = filename + "/" + tag + ".data";
	std::ofstream ofs(filename, std::ios_base::app);
	for (int i = 0; i < model.N; i++) {
		for (int a = 0; a < model.NA; a++) {
			ofs << i << '\t';	// label each atom by the molecule number
			PositionVector[i][a].pbc(model.boxl, model.invboxl).print(ofs);
		}
	}
}
