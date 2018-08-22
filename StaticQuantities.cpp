#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm> // std::upper_bound to find the value for tau_orient, std::max_element() for RoughHistogram
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
void VideoMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<Point>> &PositionVector, std::string tag,
	std::vector<int> &Rotators, std::vector<int> &Translators);
double TopTenPercent(std::vector<double> VectorToSort);
double BottomTenPercent(std::vector<double> VectorToSort);
void LabelRotatorsAndTranslators(SimulationParameters &model, std::vector<std::vector<Point>> &Microscopy_CenterOfMass, std::vector<std::vector<Point>> &Microscopy_Orientation,
	std::vector<int> &Rotators, std::vector<int> &Translators, double &OrientationMagnitude, int &Microscopy_Configurations, std::string tag);
void RoughHistogram(std::vector<double> &data, int &NumberOfBins, std::string filename);

int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);
	//std::vector<std::string> OutputFolders = { "OrderParameter", "PairCorrelation" };
	mkdir("Microscopy", ACCESSPERMS);
	mkdir("Microscopy/Rotators", ACCESSPERMS);
	mkdir("Microscopy/Translators", ACCESSPERMS);
	std::ofstream("Microscopy/" + tag + ".data");
	//std::vector<std::string> OutputFiles = FoldersToFiles(model, OutputFolders, tag);
	//HistogramInfo HistOrderParameter(0.005, int(1 / 0.005)), HistRDF(0.01, int(model.boxl / 2 / 0.01));
	//std::vector<double> OrderParameter(HistOrderParameter.num_bins), PairCorrelation(HistRDF.num_bins);
	std::vector<int> Rotators(200), Translators(200);
	int Microscopy_Configurations = 5001;
	std::vector<std::vector<Point>> Microscopy_CenterOfMass(Microscopy_Configurations, std::vector<Point>(200)), Microscopy_Orientation(Microscopy_Configurations, std::vector<Point>(200));
	double OrientationMagnitude = 0.0;
	for (int n = 0; n < model.num_cfgs; n++) {
		if (n % 1000 == 1) std::cout << n << '\n';
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);
		//WriteConfiguration(model, Positions, std::cout); // don't write yet
		if (n == 0) OrientationMagnitude = Positions[0][2].dist_sq(Positions[0][1]);	//grab orientation vector's magnitude, but only once.
		//int bin = int(CalculateS(model, OrientationMagnitude, Positions) / HistOrderParameter.bin_width);
		//OrderParameter[bin] += 1;
		//radial_df(model, Positions, HistRDF, PairCorrelation);
		if (n < Microscopy_Configurations) {	// for 5001 configurations, see who is rotating and who is translating. 
			for (int i = 0; i < model.N; i++) {
				Point Cm_i = (Positions[i][2] + Positions[i][1] + Positions[i][0]) / 3.0;
				Point Orientation_i = (Positions[i][2] - Positions[i][1]);	//don't need to normalize since I only care about relative orientations
				Microscopy_CenterOfMass[n][i] = Cm_i;
				Microscopy_Orientation[n][i] = Orientation_i;
			}
		}
		if (n == Microscopy_Configurations) { // From the 5001 configurations, calculate 5000 velocities and angular velocities. Then assign the largest 10% their colors/shapes.
			LabelRotatorsAndTranslators(model, Microscopy_CenterOfMass, Microscopy_Orientation, Rotators, Translators, OrientationMagnitude, Microscopy_Configurations, tag);
		}
		if (n >= 10000 && n < 20000 && n % 10 == 0) {
			VideoMicroscopy(model, "Microscopy", Positions, tag, Rotators, Translators); // get 1000 cfgs spaced 10 * time_betw_cfgs apart.
		}
		if (n > 20000) return 0;
	}
	//PrintHistogram(model, HistOrderParameter, OutputFiles[0], OrderParameter, "prob_density");
	//PrintHistogram(model, HistRDF, OutputFiles[1], PairCorrelation, "histogram");
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
					// bin should not be zero. if so, either exactly overlapping molecules or issue with indexing somewhere.
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

void VideoMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<Point>> &PositionVector, std::string tag, 
	std::vector<int> &Rotators, std::vector<int> &Translators) {	//print configuration with minimum images and molecule labels
	filename = filename + "/" + tag + ".data";
	std::ofstream ofs(filename, std::ios_base::app);
	for (int i = 0; i < model.N; i++) {
		for (int a = 0; a < model.NA; a++) {
			//ofs << i << '\t';	// label each atom by the molecule number
			ofs << Translators[i] << '\t' << Rotators[i] << '\t';	// label atoms by their molecule's rotation or translation
			PositionVector[i][a].pbc(model.boxl, model.invboxl).print(ofs);
		}
	}
}

double TopTenPercent(std::vector<double> VectorToSort) {
	std::sort(VectorToSort.begin(), VectorToSort.end());	//sort the thing
	int NinetiethPercentile = int(0.9 * VectorToSort.size());
	return VectorToSort[NinetiethPercentile];
}

double BottomTenPercent(std::vector<double> VectorToSort) {
	std::sort(VectorToSort.begin(), VectorToSort.end());	//sort the thing
	int TenthPercentile = int(0.1 * VectorToSort.size());
	return VectorToSort[TenthPercentile];
}

void LabelRotatorsAndTranslators(SimulationParameters &model, std::vector<std::vector<Point>> &Microscopy_CenterOfMass, std::vector<std::vector<Point>> &Microscopy_Orientation, 
	std::vector<int> &Rotators, std::vector<int> &Translators, double &OrientationMagnitude, int &Microscopy_Configurations, std::string tag) {
	int NumberOfConfigurations = Microscopy_Configurations - 1;
	std::vector<double> rotator_displacement(200), translator_displacement(200);
	for (int k = 0; k < NumberOfConfigurations; k++) {
		for (int i = 0; i < model.N; i++) {
			double cm_displacement = (Microscopy_CenterOfMass[k][i].dist(Microscopy_CenterOfMass[k + 1][i]));	// a squared displacement (not quite a correlation function)
			double orientation_displacement = (Microscopy_Orientation[k][i].dist(Microscopy_Orientation[k + 1][i]));	// a squared orientational displacement (not quite a correlation function)
			translator_displacement[i] += cm_displacement;
			rotator_displacement[i] += orientation_displacement;
		}
	}
	double Normalization = sqrt(OrientationMagnitude);
	for (auto &i : rotator_displacement) {
		i /= Normalization;
	}
	int HistogramNumBins = 20;
	RoughHistogram(rotator_displacement, HistogramNumBins, "Microscopy/Rotators/" + tag + ".data");
	RoughHistogram(translator_displacement, HistogramNumBins, "Microscopy/Translators/" + tag + ".data");
	double FastRotation = TopTenPercent(rotator_displacement);
	double FastTranslation = TopTenPercent(translator_displacement);
	double SlowRotation = BottomTenPercent(rotator_displacement);
	double SlowTranslation = BottomTenPercent(translator_displacement);
	for (std::vector<int>::size_type i = 0; i < rotator_displacement.size(); i++) {
		double translation = translator_displacement[i];
		double rotation = rotator_displacement[i];
		if (translation > FastTranslation) Translators[i] = 10;	// use parula, palette in gnuplot
		else if (translation < SlowTranslation) Translators[i] = 5;
		else Translators[i] = 1;
		if (rotation > FastRotation) Rotators[i] = 1;	// plot w circles starting from 0 degrees to 360 - 90*(integer) degrees. So fast molecules are pac-men, slow ones are half-moons.
		else if (rotation < SlowRotation) Rotators[i] = 2;
		else Rotators[i] = 0;
	}
}

void RoughHistogram(std::vector<double> &data, int &NumberOfBins, std::string filename) {
	std::ofstream ofs(filename);
	double max = *std::max_element(data.begin(), data.end());
	int num_bins = NumberOfBins;
	double bin_width = max / num_bins;
	std::vector<int> count(num_bins);
	for (auto &i : data) {	// bin values
		int bin = int(i / bin_width);	// data / (max / num_bins) == data/max * num_bins, so this line sorts the value into the proper bin.
		count[bin] += 1;
	}
	for (int bin = 0; bin < num_bins; bin++) {
		double binvalue = double(bin) * bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		ofs << binvalue << '\t' << count[bin] << '\n';
	}
}
