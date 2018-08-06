#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <iterator> // istream_iterator
#include "Point.h"

const double density = 0.25;		//these numbers should be passed to the file as a shell variable or I can have this c++ file dig through the directory to find the input file directly
const double temperature = 0.50;
const int N = 200;
const int NA = 3;
const double boxl = sqrt(N / density);
const double boxinv = 1.0 / boxl;

double CalculateS(double NormalizationSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, const int N, int configuration_number);
void PrintHistogram(std::string filename, const int bin_count, const double bin_width, const int NumConfigurations, double density, double temperature, std::vector<double> &data);
void radial_df(std::vector<std::vector<std::vector<Point>>> &PositionVector, const int bin_count, const double bin_width, const int N, int configuration_number, std::vector<double> &rdf);	
void NearestNeighbors(const int N, int NumberOfConfigurations, double NeighborCutOffSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, std::vector<int> &NeighborVector);

int main() {
	const int NumConfigs = 10000;
	std::ifstream ifs("input31.data");
	std::vector<double> parsed(std::istream_iterator<double>(ifs), {}); // initialize vector using the istream
	std::cout << "just parsed" << '\n';
	std::vector<std::vector<std::vector<Point>>> Positions(NumConfigs, std::vector<std::vector<Point>>(N, std::vector<Point>(NA, (0.0,0.0))));	//Now have positions[:,:,:] config:molecule:atom
	std::cout << parsed.size() << '\n';
	int k = 0;
	for (int n = 0; n < NumConfigs; n++) { // fill out Positions with the input file 
		for (int i = 0; i < N; i++) {
			for (int a = 0; a < NA; a++) {
				Point point(parsed[k], parsed[k + 1]);
				Positions[n][i][a] = point;
				k += 2;
			}
		}
	}
	//do P(s), g(r), etc on a configuration-by-configuration basis.
	double OrientationMagnitude = pow(Positions[0][0][2].dist(Positions[0][0][1]),2); //pbc doesn't affect orientations, don't need it here. Just be consistent - no need to worry about out of place particles since uncorrected positions don't include images.
	double OrderBinWidth = 0.005, GrBinWidth = 0.01;
	int TotalOrderBins = 1.0 / OrderBinWidth, TotalGrBins = boxl / 2 / GrBinWidth;
	std::cout << "calculating P(s) \n";
	std::vector<double> OrderParameter(TotalOrderBins), PairCorrelation(TotalGrBins);	// vector constructed with a size implies doubles are initialized to 0
	std::vector<int> NumberNeighbors(10);
	std::vector<std::vector<Point>> OmegaJ(NumConfigs, std::vector<Point>(N, (0.0,0.0)));
	double NeighborCutOffSq = 1.61 * 1.61;
	for (int n = 0; n < NumConfigs; n++) {	//Takes a long time to read these configurations. This computation part is relatively small, so might as well calculate everything at once. 
		int bin = CalculateS(OrientationMagnitude, Positions, N, n) / OrderBinWidth;
		OrderParameter[bin] += 1;
		radial_df(Positions, TotalGrBins, GrBinWidth, N, n, PairCorrelation);
		NearestNeighbors(N, n, NeighborCutOffSq, Positions, NumberNeighbors);
	}
	PrintHistogram("OrderParameter/OrderParameter1.data", TotalOrderBins, OrderBinWidth, NumConfigs, density, temperature, OrderParameter);
	PrintHistogram("PairCorrelation/PairCorrelation1.data", TotalGrBins, GrBinWidth, NumConfigs, density, temperature, PairCorrelation);
	std::ofstream neighborfile("Neighbors/Neighbors1.data");
	for (std::vector<int>::size_type i = 0; i < NumberNeighbors.size(); i++) {
		neighborfile << i << "\t" << NumberNeighbors[i] / double(N* NA) / NumConfigs << '\n';
	}
}

double CalculateS(double NormalizationSq, std::vector<std::vector<std::vector<Point>>> &PositionVector, const int N, int configuration_number) {	//take a configuration and calculate the "orientational order parameter" see Zhou, Stratt 2018
	// Advice from Yan: Pass things by reference, that was the major timesink. Advice from Andrew: pass by reference is like a check - you don't need to write one to pay for your lunch, but you might write one for a car. 
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

void PrintHistogram(std::string filename, const int bin_count, const double bin_width, const int NumConfigurations, double density, double temperature, std::vector<double> &data) {
	std::ofstream ofs(filename);
	ofs << density << '\t' << temperature << '\n';
	for (int bin = 0; bin < bin_count; bin++) {
		double binvalue = double(bin) * bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		ofs << data[bin] / NumConfigurations << '\t' << binvalue << '\n';
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