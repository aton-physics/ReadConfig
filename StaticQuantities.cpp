#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <assert.h>
#include <algorithm> // std::upper_bound to find the value for tau_orient, std::max_element() for RoughHistogram
#include <iterator> // istream_iterator
#include <iomanip> // std::setprecision
#include <cmath> // acos
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
double CalculateX6(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<Point>> &PositionVector);
void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode);
void radial_df(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf);
void specific_radial_df(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf, int which_atom);
std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag);
void VideoMicroscopy(SimulationParameters &model, std::string filename, std::vector<std::vector<Point>> &PositionVector, std::string tag,
	std::vector<int> &Rotators, std::vector<int> &Translators);
double TopTenPercent(std::vector<double> VectorToSort);
double BottomTenPercent(std::vector<double> VectorToSort);
void LabelRotatorsAndTranslators(SimulationParameters &model, std::vector<std::vector<Point>> &Microscopy_CenterOfMass, std::vector<std::vector<Point>> &Microscopy_Orientation,
	std::vector<int> &Rotators, std::vector<int> &Translators, double &OrientationMagnitude, int &Microscopy_Configurations, std::string tag);
void RoughHistogram(std::vector<double> &data, int &NumberOfBins, std::string filename);
//std::vector<Point> InitializeWaveVectors(int NumberOfWaveVectorsPerDirection);
//void compute_structure_factor(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &sq, std::vector<Point> &WaveVectors);

int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);
	mkdir("OrderParameter", ACCESSPERMS);
	std::vector<std::string> OutputFolders = { "OrderParameter/Nematic", "OrderParameter/Hexatic", "PairCorrelation/Atom"};
	mkdir("Microscopy", ACCESSPERMS);
	mkdir("Microscopy/Rotators", ACCESSPERMS);
	mkdir("Microscopy/Translators", ACCESSPERMS);
	std::ofstream("Microscopy/" + tag + ".data");
	std::vector<std::string> OutputFiles = FoldersToFiles(model, OutputFolders, tag);
	//int NumberOfWaveVectors = 26;
	//std::vector<Point> WaveVectors = InitializeWaveVectors(NumberOfWaveVectors);
	double deltaK = 2 * 3.14159265359 / model.boxl;	// minimum wave vector increment = 2*pi/boxl
	HistogramInfo HistOrderParameter(0.005, int(1 / 0.005)), HistRDF(0.01, int(model.boxl / 2 / 0.01)), HistX6(0.005, int(1 / 0.005));
		//HistStructureFactor(deltaK, int (sqrt(2 * pow(NumberOfWaveVectors - 1, 2))  / (deltaK)));	// histogram info - bin width, bin number
	std::vector<double> OrderParameter(HistOrderParameter.num_bins), PairCorrelation(HistRDF.num_bins), PairCorrelationCenter(HistRDF.num_bins), X6OrderParameter(HistX6.num_bins);
		//StructureFactor(HistStructureFactor.num_bins);	// histogram data 
	int Microscopy_Configurations = 5002;
	std::vector<int> Rotators(model.N), Translators(model.N);
	std::vector<std::vector<Point>> Microscopy_CenterOfMass(Microscopy_Configurations, std::vector<Point>(model.N)), Microscopy_Orientation(Microscopy_Configurations, std::vector<Point>(model.N));
	double OrientationMagnitude = 0.0;
	for (int n = 0; n < model.num_cfgs; n++) {
		//if (n >= 5000) continue;	//note that I am limiting # configurations for time reasons
		//if (n % 1000 == 1) std::cout << n << '\n';
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);
		if (n == 0) OrientationMagnitude = Positions[0][2].dist_sq(Positions[0][1]);	//grab orientation vector's magnitude, but only once.
		int NematicBin = int(CalculateS(model, OrientationMagnitude, Positions) / HistOrderParameter.bin_width);
		OrderParameter[NematicBin] += 1;
		radial_df(model, Positions, HistRDF, PairCorrelation);
		int HexaticBin = int(CalculateX6(model, OrientationMagnitude, Positions) / HistX6.bin_width);
		X6OrderParameter[HexaticBin] += 1;
		//compute_structure_factor(model, Positions, HistStructureFactor, StructureFactor, WaveVectors);
		if (n < Microscopy_Configurations) {	// for 5002 configurations, see who is rotating and who is translating. Need n+2 configurations to measure n velocities
			for (int i = 0; i < model.N; i++) {
				Point Cm_i = (Positions[i][2] + Positions[i][1] + Positions[i][0]) / 3.0;
				Point Orientation_i = (Positions[i][2] - Positions[i][1]);
				Orientation_i = Orientation_i / sqrt(Orientation_i.dot(Orientation_i));
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
	}
	PrintHistogram(model, HistX6, OutputFiles[1], X6OrderParameter, "prob_density");
	PrintHistogram(model, HistOrderParameter, OutputFiles[0], OrderParameter, "prob_density");
	PrintHistogram(model, HistRDF, OutputFiles[2], PairCorrelation, "histogram");
	//PrintHistogram(model, HistStructureFactor, OutputFiles[0], StructureFactor, "histogram");
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
	int num_mol = Input.N, num_atom = 3, numcfg = Input.NumConfigs, steps_between_cfgs = 100;
	double boxl = sqrt(num_mol / density);
	double boxinv = 1.0 / boxl;
	SimulationParameters model(density, temperature, boxl, boxinv, ts, numcfg, num_mol, num_atom, steps_between_cfgs);
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

double CalculateX6(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<Point>> &PositionVector) {	
	double X6 = 0;
	for (int j = 0; j < model.N; j++) {
		Point OmegaJ = (PositionVector[j][2] - PositionVector[j][1]);
		OmegaJ = OmegaJ / sqrt(OmegaJ.dot(OmegaJ));
		for (int k = 0; k < model.N; k++) {
			Point OmegaK = PositionVector[k][2] - PositionVector[k][1];
			OmegaK = OmegaK / sqrt(OmegaK.dot(OmegaK));
			double c = OmegaJ.dot(OmegaK); // ;
			double c2 = c*c;
			double c4 = c2*c2;
			double value = 4 * (8 * c4*c2 - 12 * c4 + 6 * c2 - 1) - 6 * c2 + 3;
			X6 += value;
		}
	}
	X6 /= pow(model.N, 2);
	return X6;
}


void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode) {
	std::ofstream ofs(filename, std::ios_base::app);
	for (int bin = 0; bin < hist.num_bins; bin++) {
		double binvalue = double(bin) * hist.bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / model.num_cfgs << '\n';
		//if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / 5000 << '\n';
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

void specific_radial_df(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &rdf, int which_atom) {	// atom-atom g(r). function takes in Positions and fills out a vector of doubles
	const double constant = 3.14159265 * model.density;
	double rlower, rupper, nideal, rij;
	int atom_index = which_atom;
	std::vector<double> histo(hist.num_bins);
	int bin;
	for (int j = 0; j < model.N - 1; j++) {				//loop over all distinct molecule pairs j,k
		for (int k = j + 1; k < model.N; k++) {
			Point r_jakb = PositionVector[j][atom_index] - PositionVector[k][atom_index]; //calculate pair separation between similar atoms
			r_jakb = r_jakb.pbc(model.boxl, model.invboxl);	// fold the pair separation to get the minimum image separation
			rij = sqrt(r_jakb.dot(r_jakb));
			bin = int(rij / hist.bin_width);
			// bin should not be zero. if so, either exactly overlapping molecules or issue with indexing somewhere.
			if (bin < hist.num_bins) {	// make sure we only count particles out as far as some maximum (set to be half the box length)
				histo[bin] += 2;		// if rij falls within the range of the bin, count two particles (i and j)
			}
		}
	}
	for (bin = 0; bin < hist.num_bins; bin++) {
		rlower = bin * hist.bin_width;
		rupper = rlower + hist.bin_width;
		nideal = constant * (rupper * rupper - rlower * rlower);
		rdf[bin] += double(histo[bin]) / model.N / nideal;
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
	int NumberOfConfigurations = Microscopy_Configurations - 2;
	std::vector<double> rotator_displacement(model.N), translator_displacement(model.N);
	for (int k = 0; k < NumberOfConfigurations; k++) {
		for (int i = 0; i < model.N; i++) {	//recover angular, cm velocity using finite difference (2-step)
			double cm_displacement = (Microscopy_CenterOfMass[k][i].dist(Microscopy_CenterOfMass[k + 2][i]));
			double angular_displacement_temporary = fabs(acos(Microscopy_Orientation[k][i].dot(Microscopy_Orientation[k + 2][i])));
			double orientation_displacement = (fabs(angular_displacement_temporary) < 0.9999 ? angular_displacement_temporary : 0.0);
			translator_displacement[i] += cm_displacement;
			rotator_displacement[i] += orientation_displacement;
		}
	}
	for (auto &i : rotator_displacement) {
		i /= (2 * NumberOfConfigurations * model.steps_between_cfgs * model.timestep);	// 2 since finite difference is over 2 steps
	}
	for (auto &i : translator_displacement) {
		i /= (2 * NumberOfConfigurations * model.steps_between_cfgs * model.timestep);
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
		if (translation >= FastTranslation) Translators[i] = 2;	// 2 is fast, 1 is slow, 0 is normal
		else if (translation < SlowTranslation) Translators[i] = 1;
		else Translators[i] = 0;
		if (rotation >= FastRotation) Rotators[i] = 2;
		else if (rotation < SlowRotation) Rotators[i] = 1;
		else Rotators[i] = 0;
	}
}

void RoughHistogram(std::vector<double> &data, int &NumberOfBins, std::string filename) {
	std::ofstream ofs(filename);
	double max = *std::max_element(data.begin(), data.end());
	int num_bins = NumberOfBins;
	double bin_width = double(max) / double(num_bins);
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

/*std::vector<Point> InitializeWaveVectors(int NumberOfWaveVectorsPerDirection) {	// create a list of wave vectors
	int NumberOfWaveVectors = pow(NumberOfWaveVectorsPerDirection, 2);
	std::vector<Point> WaveVectors;
	for (int i = 0; i < NumberOfWaveVectorsPerDirection; i++) {
		for (int a = 0; a < NumberOfWaveVectorsPerDirection; a++) {
			Point point(i, a);
			WaveVectors.push_back(point);	// build up a list of wavevectors consisting of integer pairs from 0 to NumberOfWaveVectorsPerDirection
		}
	}
	return WaveVectors;
}*/

/*void compute_structure_factor(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, HistogramInfo &hist, std::vector<double> &sq, std::vector<Point> &WaveVectors) {
	int NumberOfWaveVectors = int(WaveVectors.size());
	std::vector<double> NumberInOneBin(hist.num_bins);
	std::vector<double> histo(hist.num_bins);
	int bin;
	for (int k = 0; k < NumberOfWaveVectors; k++) {
		double c = 0, s = 0;				// cosine and sine
		for (int i = 0; i < model.N; i++) {	// loop over all distinct molecule pairs j,k
			/*for (int a = 0; a < model.NA; a++) {
				Point r_ia = PositionVector[i][a].pbc(model.boxl, model.invboxl);	// get folded position of atom
				double kr = WaveVectors[k].dot(r_ia);
				c += cos(kr);
				s += sin(kr);
			}*//*
			Point r_i = (PositionVector[i][0] + PositionVector[i][1] + PositionVector[i][2]) / 3.0;
			r_i = r_i.pbc(model.boxl, model.invboxl);
			double kr = WaveVectors[k].dot(r_i);
			c += cos(kr);
			s += sin(kr);
		}
		c = c*c;
		s = s*s;
		double k_magnitude = sqrt(WaveVectors[k].dot(WaveVectors[k]));	// compute magnitude of wavevector
		bin = int(k_magnitude / hist.bin_width);
		histo[bin] = c + s;
		if (bin < hist.num_bins) {	// don't bin anything past the max(domain)
			NumberInOneBin[bin] += 1;		// count how many wavevectors fall in the wavevector bin from (k, k + dk)
		}
	}
	for (int bin = 0; bin < hist.num_bins; bin++) {
		//sq[bin] += double(histo[bin]) / (model.N * model.NA);
		sq[bin] += double(histo[bin]) / model.N;
		if (NumberInOneBin[bin] > 1) {
			sq[bin] /= NumberInOneBin[bin];	// only divide if there is more than one contribution to the bin
		}
	}
}*/

