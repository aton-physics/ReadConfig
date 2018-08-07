#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm> // std::upper_bound to find the value for tau_orient 
#include <iterator> // istream_iterator
#include "Point.h"
#include "Numerical.h"

const double density = 0.25;		//these numbers should be passed to the file as a shell variable or I can have this c++ file dig through the directory to find the input file directly
const double temperature = 0.50;
const int N = 200;
const int NA = 3;
const double boxl = sqrt(N / density);
const double boxinv = 1.0 / boxl;

double CalculateS(double NormalizationSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, const int N, int configuration_number);
void PrintHistogram(std::string filename, const int bin_count, const double bin_width, const int NumConfigurations, double density, double temperature, std::vector<double> &data, std::string mode = "");
void radial_df(std::vector<std::vector<std::vector<Point>>> &PositionVector, const int bin_count, const double bin_width, const int N, int configuration_number, std::vector<double> &rdf);
void NearestNeighbors(const int N, int NumberOfConfigurations, double NeighborCutOffSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<int> &NeighborVector);
void FillPositionVector(const int NumberOfConfigurations, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<double> &ParsedInput);

int main() {
	const int NumConfigs = 2 * 100000;
	double dt = 0.001;
	int StepsBetweenConfigurations = 10;
	std::string tag = "31";
	std::ifstream ifs("Trajectory" + tag + ".data");
	std::vector<double> parsed(std::istream_iterator<double>(ifs), {}); // initialize vector using the istream
	std::cout << "just parsed" << '\n';
	std::vector<std::vector<std::vector<Point>>> Positions(NumConfigs, std::vector<std::vector<Point>>(N, std::vector<Point>(NA)));	//Now have positions[:,:,:] config:molecule:atom
	std::cout << parsed.size() << '\n';
	FillPositionVector(NumConfigs, Positions, parsed);	// copy parsed input file into Positions
	std::cout << "just filled \n";
	//do P(s), g(r), etc on a configuration-by-configuration basis.
	double OrientationMagnitude = Positions[5][0][2].dist_sq(Positions[5][0][1]); //pbc doesn't affect orientations, don't need it here. Just be consistent - no need to worry about out of place particles since uncorrected positions don't include images.
	double OrderBinWidth = 0.005, GrBinWidth = 0.01;
	int TotalOrderBins = int(1/ OrderBinWidth), TotalGrBins = int(boxl / 2 / GrBinWidth);
	int t_cor = int(1000 / dt / StepsBetweenConfigurations); // correlation function is measured for this long (in units of dt * #skipped configurations)
	std::vector<double> OrderParameter(TotalOrderBins), PairCorrelation(TotalGrBins), OrientCorrelation(t_cor), MeanSqDisplacement(t_cor);
	std::vector<int> NumberNeighbors(10);
	std::vector<std::vector<Point>> OmegaJ(NumConfigs, std::vector<Point>(N)), cm_J(NumConfigs, std::vector<Point>(N));
	double NeighborCutOffSq = 1.61 * 1.61;
	for (int n = 0; n < NumConfigs; n++) {
		int bin = int(CalculateS(OrientationMagnitude, Positions, N, n) / OrderBinWidth);
		OrderParameter[bin] += 1;
		radial_df(Positions, TotalGrBins, GrBinWidth, N, n, PairCorrelation);
		NearestNeighbors(N, n, NeighborCutOffSq, Positions, NumberNeighbors);
		for (int i = 0; i < N; i++) {	// Store orientations, center of mass of each molecule for every configuration to calculate C(t), MSD(t) respectively.
			OmegaJ[n][i] = (Positions[n][i][2].sub(Positions[n][i][1]));
			cm_J[n][i] = Positions[n][i][2].add(Positions[n][i][1]).add(Positions[n][i][0]).mult(1.0/3.0);
		}
	}
	std::cout << "completed main calculation \n";
	int NumberOfTrajectories = NumConfigs - t_cor + 1; // every new configuration constitutes new initial conditions when calculating correlation functions. But we measure correlation up to CorrTimeScale, so we have that many fewer trajectories
	for (int j = 0; j < NumberOfTrajectories; j++) {
		for (int k = 0; k < t_cor; k++) {
			double total_omega = 0.0;
			double total_cm = 0.0;
			for (int i = 0; i < N; i++) {
				total_omega += pow(OmegaJ[j][i].dot(OmegaJ[j + k][i]), 2);
				Point cm_temp = cm_J[j][i].sub(cm_J[j + k][i]);
				total_cm += pow(cm_temp.x(),2) + pow(cm_temp.y(),2);
			}
			//OrientCorrelation[k] += 2.0 / double(N) * total_omega - 1;
			OrientCorrelation[k] += total_omega;
			MeanSqDisplacement[k] += total_cm;
		}
	}
	std::cout << "completed correlation calculation \n";
	std::ofstream orientfile("OrientCorrelation/OrientCorrelation" + tag + ".data");
	std::ofstream msdfile("msd/msd" + tag + ".data");
	orientfile << density << '\t' << temperature << '\n';
	msdfile << density << '\t' << temperature << '\n';
	for (int k = 0; k < t_cor; k++) {
		OrientCorrelation[k] = 2.0 * OrientCorrelation[k] / pow(OrientationMagnitude,2)/ double(N) / NumberOfTrajectories - 1;	// divide here so I can use std::upper_bound later
		MeanSqDisplacement[k] /= (double(N) * double(NumberOfTrajectories));
		double what_time = k * dt * StepsBetweenConfigurations;
		orientfile << what_time << '\t' << OrientCorrelation[k] << '\n';
		msdfile << what_time << '\t' << MeanSqDisplacement[k] << '\n';
	}
	std::ofstream diffusionstream("Diffusion/Diffusion" + tag + ".data");
	diffusionstream << density << '\t' << temperature << '\n';
	numerical_differentiation("msd/msd" + tag + ".data", "derivative/derivative" + tag + ".data", StepsBetweenConfigurations); // derivative/derivative is a temporary solution to see how long I should wait before taking the zero slope regression line
	//best way to do this is to not print anything but diffusion (maybe msd). send msd to numerical_diff to zero_slope > diffusion.data, don't ever want to look at derivative
	zero_slope_regression("derivative/derivative" + tag + ".data", "Diffusion/Diffusion" + tag + ".data", temperature, 100.0);	// skip 100 tau before taking the derivative
	int index = distance(OrientCorrelation.begin(), std::upper_bound(OrientCorrelation.begin(), OrientCorrelation.end(), 2.71928)); // search the (ordered) vector for the last element larger than euler's number, retrieve the index
	if (index == int(OrientCorrelation.size())) index = 9999999;
	std::ofstream tau_orientfile("FunctionsOfTemperature/Tau_orient/Tau" + tag + ".data");
	tau_orientfile << density << '\t' << temperature << '\n';
	tau_orientfile << double(index) * dt * StepsBetweenConfigurations << '\t' << temperature << '\n';
	PrintHistogram("OrderParameter/OrderParameter" + tag + ".data", TotalOrderBins, OrderBinWidth, NumConfigs, density, temperature, OrderParameter, "prob_density");
	PrintHistogram("PairCorrelation/PairCorrelation" + tag + ".data", TotalGrBins, GrBinWidth, NumConfigs, density, temperature, PairCorrelation, "histogram");
	std::ofstream neighborfile("Neighbors/Neighbors" + tag + ".data");
	neighborfile << density << '\t' << temperature << '\n';
	for (std::vector<int>::size_type i = 0; i < NumberNeighbors.size(); i++) {
		neighborfile << i << "\t" << NumberNeighbors[i] / double(N* NA) / NumConfigs << '\n';
	}
}

double CalculateS(double NormalizationSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, const int N, int configuration_number) {	//take a configuration and calculate the "orientational order parameter" see Zhou, Stratt 2018
	double s = 0, Q1 = 0, Q2 = 0;
	for (int i = 0; i < N; i++) {
		Point OmegaJ = PositionVector[configuration_number][i][2].sub(PositionVector[configuration_number][i][1]);
		//Q1 += 2 * pow(OmegaJ.x(), 2) - 1;
		//Q2 += 2 * OmegaJ.x() * OmegaJ.y();
		Q1 += pow(OmegaJ.x(), 2);
		Q2 += OmegaJ.x() * OmegaJ.y();
	}
	Q1 = 2*Q1 / NormalizationSq - N;	// each of OmegaJ.x and OmegaJ.y carry a factor of 1/Normalization, so their product has a factor of 1/NormalizationSq.
	Q2 *= 2 / NormalizationSq;
	s = sqrt(Q1 * Q1 + Q2 * Q2) / N;
	return s;
}

void PrintHistogram(std::string filename, const int bin_count, const double bin_width, const int NumConfigurations, double density, double temperature, std::vector<double> &data, std::string mode) {
	std::ofstream ofs(filename);
	ofs << density << '\t' << temperature << '\n';
	for (int bin = 0; bin < bin_count; bin++) {
		double binvalue = double(bin) * bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / NumConfigurations << '\n';
		else if (mode == "prob_density") ofs << binvalue << '\t' << data[bin] / NumConfigurations / bin_width << '\n'; // if optional parameter is specified and greater than 1, divide bin count by bin width. This turns a histogram into a probability density
		else std::cout << "Bad parameter for PrintHistogram()";
	}
}

void radial_df(std::vector<std::vector<std::vector<Point>>> &PositionVector, const int bin_count, const double bin_width, const int N, int configuration_number, std::vector<double> &rdf) {	// atom-atom g(r). function takes in Positions and fills out a vector of doubles
	const double constant = 3.14159265 * density * NA;
	double rlower, rupper, nideal, rij;
	std::vector<double> hist(bin_count);
	int bin;
	for (int j = 0; j < N - 1; j++) {				//loop over all distinct molecule pairs j,k
		for (int k = j + 1; k < N; k++) {
			for (int a = 0; a < NA; a++) {
				for (int b = 0; b < NA; b++) {
					Point r_jakb = PositionVector[configuration_number][j][a].sub(PositionVector[configuration_number][k][b]); //calculate pair separation
					r_jakb = r_jakb.pbc(boxl, boxinv);	// fold the pair separation to get the minimum image separation
					rij = sqrt(r_jakb.dot(r_jakb));	
					bin = int(rij / bin_width);	
					if (bin < bin_count) {	// make sure we only count particles out as far as some maximum (set to be half the box length)
						hist[bin] += 2;		// if rij falls within the range of the bin, count two particles (i and j)
					}
				}
			}
		}
	}
	for (bin = 0; bin < bin_count; bin++) {
		rlower = double(bin) * bin_width;
		rupper = rlower + bin_width;
		nideal = constant * (rupper * rupper - rlower * rlower);
		rdf[bin] += double(hist[bin]) / double(N * NA) / nideal;
	}
}

void NearestNeighbors(const int N, int ConfigurationNumber, double NeighborCutOffSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<int> &NeighborVector) {
	for (int j = 0; j < N; j++) {				// loop over all distinct molecule pairs j,k - could presumably do j,k with k > j but I'm not sure of a neat way to count neighbors without going over every index individually, can't double count like we do in g(r).
		for (int a = 0; a < NA; a++) {
			int counter = 0;					// reset counter when examining new reference particle
			for (int k = 0; k < N; k++) {
				if (j == k) continue;			// skip neighbor check if loop is over the same molecule twice <=> make sure molecules are distinct
				for (int b = 0; b < NA; b++) {
					Point r_jakb = PositionVector[ConfigurationNumber][j][a].sub(PositionVector[ConfigurationNumber][k][b]).pbc(boxl, boxinv);
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

void FillPositionVector(const int NumberOfConfigurations, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<double> &ParsedInput) {
	int k = 0;
	for (int n = 0; n < NumberOfConfigurations; n++) { // fill out Positions with the input file 
		for (int i = 0; i < N; i++) {
			for (int a = 0; a < NA; a++) {
				Point point(ParsedInput[k], ParsedInput[k + 1]);
				PositionVector[n][i][a] = point;
				k += 2;
			}
		}
	}
}