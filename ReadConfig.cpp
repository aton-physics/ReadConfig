#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm> // std::upper_bound to find the value for tau_orient 
#include <iterator> // istream_iterator
#include "headers/Point.h"
#include "headers/Numerical.h"

//***********************************8/6/2018**********************************
//This program takes in a potentially gigantic trajectory file from molecular dynamics (gigabytes) and calculates various quantites including correlation functions, diffusion via msd, and histograms/probability densities.
//Input: The trajectory file (sequential lines of 'x' \t 'y'); and somehow grab number of molecules N, number of atoms NA, density, temperature. 
//Profiler: most of time is spent calculated pbc, not too surprising. 

SimulationParameters GetParams(std::string inputfile);
std::vector<std::vector<std::vector<Point>>> ParseTrajectoryFile(SimulationParameters &model, std::string TrajectoryFile);
double CalculateS(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, int &configuration_number);
void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode);
void radial_df(SimulationParameters &model, std::vector<std::vector<std::vector<Point>>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf, int &configuration_number);	
void NearestNeighbors(SimulationParameters &model, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<int> &NeighborVector, int &ConfigurationNumber, double &NeighborCutOffSq);
std::vector<std::vector<std::vector<Point>>> FillPositionVector(SimulationParameters &model, std::vector<double> &ParsedInput);
void PrintHeaders(SimulationParameters &model, std::vector<std::string> &stringvector);	
void FindRelaxationTime(SimulationParameters &model, std::string filename, std::vector<double> &OrientCorrelation);
void PrintNeighbors(SimulationParameters &model, std::string filename, std::vector<int> &NumberNeighbors);
void PrintMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<std::vector<Point>>> &PositionVector, int num_frames);

int main() {
	std::string tag = "5";
	SimulationParameters model = GetParams("dummystring");
	std::vector<std::string> AllOutputFiles = { "OrientCorrelation/OrientCorrelation" + tag + ".data", "msd/msd" + tag + ".data", "Diffusion/Diffusion" + tag + ".data",
		"FunctionsOfTemperature/Tau_orient/Tau" + tag + ".data", "OrderParameter/OrderParameter" + tag + ".data", "PairCorrelation/PairCorrelation" + tag + ".data",
		"Neighbors/Neighbors" + tag + ".data" };
	PrintHeaders(model, AllOutputFiles); // give all the output files density and temperature headers. does not include the derivative file.
	std::vector<std::vector<std::vector<Point>>> Positions = ParseTrajectoryFile(model, "Trajectory/Trajectory31.data"); // read Trajectory file, write to Positions[:,:,:] configuration:molecule:atom
	HistogramInfo HistOrderParameter(0.005, int(1 / 0.005)), HistRDF(0.01, int(model.boxl / 2 / 0.01));
	int t_cor = int(500 / model.timestep / model.steps_between_cfgs); // correlation function is measured for this long (in units of model.timestep* #skipped configurations)
	assert(model.num_cfgs >= t_cor); //make sure there at least as many configurations as correlation steps
	std::vector<double> OrderParameter(HistOrderParameter.num_bins), PairCorrelation(HistRDF.num_bins), OrientCorrelation(t_cor), MeanSqDisplacement(t_cor);
	std::vector<int> NumberNeighbors(10);
	std::vector<std::vector<Point>> OmegaJ(model.num_cfgs, std::vector<Point>(model.N)), Cm_J(model.num_cfgs, std::vector<Point>(model.N));
	double NeighborCutOffSq = 1.61 * 1.61;
	double OrientationMagnitude = Positions[5][0][2].dist_sq(Positions[5][0][1]);
	for (int n = 0; n < model.num_cfgs; n++) {
		int bin = int(CalculateS(model, OrientationMagnitude, Positions, n) / HistOrderParameter.bin_width);
		OrderParameter[bin] += 1;
		radial_df(model, Positions, HistRDF, PairCorrelation, n);
		NearestNeighbors(model, Positions, NumberNeighbors, n, NeighborCutOffSq);
		for (int i = 0; i < model.N; i++) {	// Store orientations, center of mass of each molecule for every configuration to calculate C(t), MSD(t) respectively.
			OmegaJ[n][i] = Positions[n][i][2] - Positions[n][i][1];
			Cm_J[n][i] = (Positions[n][i][2] + Positions[n][i][1] + Positions[n][i][0]) / 3.0 ;
		}
	}
	std::cout << "completed main calculation \n";
	int NumberOfTrajectories = model.num_cfgs - t_cor + 1; // every new configuration constitutes new initial conditions when calculating correlation functions. But we measure correlation up to CorrTimeScale, so we have that many fewer trajectories
	for (int j = 0; j < NumberOfTrajectories; j++) {
		for (int k = 0; k < t_cor; k++) {
			double total_omega = 0.0, total_cm = 0.0;
			for (int i = 0; i < model.N; i++) {
				total_omega += pow(OmegaJ[j][i].dot(OmegaJ[j + k][i]), 2);
				Point cm_temp = Cm_J[j][i] - Cm_J[j + k][i];
				total_cm += pow(cm_temp.x(),2) + pow(cm_temp.y(),2);
			}
			OrientCorrelation[k] += total_omega;
			MeanSqDisplacement[k] += total_cm;
		}
	}
	std::cout << "completed correlation calculation \n";
	std::ofstream orientfile("OrientCorrelation/OrientCorrelation" + tag + ".data");
	std::ofstream msdfile("msd/msd" + tag + ".data");
	for (int k = 0; k < t_cor; k++) {
		OrientCorrelation[k] = 2.0 * OrientCorrelation[k] / pow(OrientationMagnitude,2)/ double(model.N) / NumberOfTrajectories - 1;	// divide here so I can use std::upper_bound later
		MeanSqDisplacement[k] /= (double(model.N) * double(NumberOfTrajectories));
		double what_time = k * model.timestep* model.steps_between_cfgs;
		orientfile << what_time << '\t' << OrientCorrelation[k] << '\n';
		msdfile << what_time << '\t' << MeanSqDisplacement[k] << '\n';
	}
	numerical_differentiation("msd/msd" + tag + ".data", "derivative/derivative" + tag + ".data", model.steps_between_cfgs*model.timestep); // derivative/derivative is a temporary solution to see how long I should wait before taking the zero slope regression line
	//best way to do this is to not print anything but diffusion (maybe msd). send msd to numerical_diff to zero_slope > diffusion.data, don't ever want to look at derivative
	zero_slope_regression("derivative/derivative" + tag + ".data", "FunctionsOfTemperature/Diffusion_trans/Diffusion" + tag + ".data", model.temperature, 10.0);	// skip 10 tau before taking the regression
	FindRelaxationTime(model, "FunctionsOfTemperature/Tau_orient/Tau" + tag + ".data", OrientCorrelation);
	PrintHistogram(model, HistOrderParameter, "OrderParameter/OrderParameter" + tag + ".data", OrderParameter, "prob_density");
	PrintHistogram(model, HistRDF, "PairCorrelation/PairCorrelation" + tag + ".data", PairCorrelation, "histogram");
	PrintNeighbors(model, "Neighbors/Neighbors" + tag + ".data", NumberNeighbors);
	PrintMicroscopy(model, "Microscopy/Microscopy" + tag + ".data", Positions, 100);
	return 0;
}

SimulationParameters GetParams(std::string inputfile) { //read in parameter file, assign all the parameters.
	double density = 0.25, temperature = 0.5, ts = 0.001;
	int num_mol = 200, num_atom = 3, numcfg = 2*100000;
	double boxl = sqrt(num_mol / density);
	double boxinv = 1.0 / boxl;
	SimulationParameters model(density, temperature, boxl, boxinv, ts, numcfg, num_mol, num_atom);
	return model;
}

std::vector<std::vector<std::vector<Point>>> ParseTrajectoryFile(SimulationParameters &model, std::string TrajectoryFile) {
	std::ifstream ifs(TrajectoryFile);
	std::vector<double> parsed(std::istream_iterator<double>(ifs), {}); // initialize vector using the istream
	std::vector<std::vector<std::vector<Point>>> Positions = FillPositionVector(model, parsed);	// copy parsed input file into Positions
	return Positions;
}

double CalculateS(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, int &configuration_number) {	//take a configuration and calculate the "orientational order parameter" see Zhou, Stratt 2018
	double s = 0, Q1 = 0, Q2 = 0;
	for (int i = 0; i < model.N; i++) {
		Point OmegaJ = PositionVector[configuration_number][i][2] - PositionVector[configuration_number][i][1];
		//Q1 += 2 * pow(OmegaJ.x(), 2) - 1;
		//Q2 += 2 * OmegaJ.x() * OmegaJ.y();
		Q1 += pow(OmegaJ.x(), 2);
		Q2 += OmegaJ.x() * OmegaJ.y();
	}
	Q1 = 2*Q1 / NormalizationSq - model.N;	// each of OmegaJ.x and OmegaJ.y carry a factor of 1/Normalization, so their product has a factor of 1/NormalizationSq.
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

void radial_df(SimulationParameters &model, std::vector<std::vector<std::vector<Point>>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf, int &configuration_number) {	// atom-atom g(r). function takes in Positions and fills out a vector of doubles
	const double constant = 3.14159265 * model.density * model.NA;
	double rlower, rupper, nideal, rij;
	std::vector<double> histo(hist.num_bins);
	int bin;
	for (int j = 0; j < model.N - 1; j++) {				//loop over all distinct molecule pairs j,k
		for (int k = j + 1; k < model.N; k++) {
			for (int a = 0; a < model.NA; a++) {
				for (int b = 0; b < model.NA; b++) {
					Point r_jakb = PositionVector[configuration_number][j][a] - PositionVector[configuration_number][k][b]; //calculate pair separation
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

void NearestNeighbors(SimulationParameters &model, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<int> &NeighborVector, int &ConfigurationNumber, double &NeighborCutOffSq) {
	for (int j = 0; j < model.N; j++) {				// loop over all distinct molecule pairs j,k - could presumably do j,k with k > j but I'm not sure of a neat way to count neighbors without going over every index individually, can't double count like we do in g(r).
		for (int a = 0; a < model.NA; a++) {
			int counter = 0;					// reset counter when examining new reference particle
			for (int k = 0; k < model.N; k++) {
				if (j == k) continue;			// skip neighbor check if loop is over the same molecule twice <=> make sure molecules are distinct
				for (int b = 0; b < model.NA; b++) {
					//Point r_jakb = (PositionVector[ConfigurationNumber][j][a] - PositionVector[ConfigurationNumber][k][b]).pbc(model.boxl, model.invboxl);
					Point r_jakb = (PositionVector[ConfigurationNumber][j][a] - PositionVector[ConfigurationNumber][k][b]);
					r_jakb = r_jakb.pbc(model.boxl, model.invboxl);
					if (r_jakb.dot(r_jakb) <= NeighborCutOffSq) {	// check if close enough to be a neighbor
						counter += 1;
					}
				}
			}
			assert(counter <= 10);
			NeighborVector[counter] += 1; // If InRange, then register another near neighbor. 
		}
	}
}

std::vector<std::vector<std::vector<Point>>> FillPositionVector(SimulationParameters &model, std::vector<double> &ParsedInput) {
	std::vector<std::vector<std::vector<Point>>> Positions(model.num_cfgs, std::vector<std::vector<Point>>(model.N, std::vector<Point>(model.NA)));	//Now have positions[:,:,:] config:molecule:atom
	int k = 0;
	for (int n = 0; n < model.num_cfgs; n++) { // fill out Positions with the input file 
		for (int i = 0; i < model.N; i++) {
			for (int a = 0; a < model.NA; a++) {
				Point point(ParsedInput[k], ParsedInput[k + 1]);
				Positions[n][i][a] = point;
				k += 2;
			}
		}
	}
	return Positions;
}

void PrintHeaders(SimulationParameters &model, std::vector<std::string> &stringvector) {	// to each of the files in stringvector, give them a header specifying density and temperature so I can plot with titles
	for (auto &i : stringvector) {
		std::ofstream ofs(i);
		ofs << model.density << '\t' << model.temperature << '\n';
	}
}

void FindRelaxationTime(SimulationParameters &model, std::string filename, std::vector<double> &OrientCorrelation) {
	int index = distance(OrientCorrelation.begin(), std::lower_bound(OrientCorrelation.begin(), OrientCorrelation.end(), 2.71928)); // search the (ordered) vector for first element smaller than euler's number, retrieve the index
	if (index == int(OrientCorrelation.size())) index = 9999999;	//if 
	std::ofstream tau_orientfile(filename, std::ios_base::app);
	tau_orientfile << index * model.timestep* model.steps_between_cfgs << '\t' << model.temperature << '\n';
}

void PrintNeighbors(SimulationParameters &model, std::string filename, std::vector<int> &NumberNeighbors) {
	std::ofstream neighborfile(filename, std::ios_base::app);
	for (std::vector<int>::size_type i = 0; i < NumberNeighbors.size(); i++) {
		neighborfile << i << "\t" << NumberNeighbors[i] / double(model.N * model.NA) / model.num_cfgs << '\n';
	}
}

void PrintMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<std::vector<Point>>> &PositionVector, int num_frames) {
	std::ofstream ofs(filename);
	for (int n = 0; n < num_frames; n++) {
		for (int i = 0; i < model.N; i++) {
			for (int a = 0; a < model.NA; a++) {
				ofs << i << '\t';	// label each atom by the molecule number
				PositionVector[n][i][a].pbc(model.boxl, model.invboxl).print(ofs);
			}
		}
	}
}