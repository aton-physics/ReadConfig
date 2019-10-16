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
//***1/7/2019****
std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model);
void WriteConfiguration(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::ostream &strm);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
InputParameter ReadInput(std::string &tag);
SimulationParameters GetParams(InputParameter Input);
double ComputeCosSixTheta(double &CosTheta);
double ComputeSinSixTheta(double &CosTheta, double &SinTheta);
std::vector<double> Fluctuations(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, double &magnitude_Phi6, double &Q1, double &Q2);
std::vector<double> CalculatePhi6(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector);
void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode);
std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag);

int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);
	//std::vector<std::string> OutputFolders = {"OrderParameter/6"};
	//std::vector<std::string> OutputFiles = FoldersToFiles(model, OutputFolders, tag);
	HistogramInfo HistOrderParameter(0.005, int(1 / 0.005)), HistX6(0.005, int(1 / 0.005));
	std::vector<double> OrderParameter(HistOrderParameter.num_bins, 0), BO6(HistX6.num_bins, 0), X6OrderParameter(HistX6.num_bins, 0);
	double Chi6 = 0.0;	//Chi6 = N <|phi_6|^2>
	double Avg_s6 = 0.0; //<s6> = <1/N \sum_{j}^N cos(6\theta_j)>
	double Chi6x = 0.0; //Chi6x = 1/N <\sum_{i,j=1}^N cos(6\theta_i)cos(6\theta_j)>
	double Chi6y = 0.0;
	int TemporaryNumCfg = 10000;
	for (int n = 0; n < TemporaryNumCfg; n++) {
		// for (int n = 0; n < model.num_cfgs; n++){
		//if (n % 100 == 0) std::cout << n << '\n';
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);	//read in coordinates at a certain time step
		std::vector<double> OrderParameters = CalculatePhi6(model,Positions);
		std::vector<double> Flucts = Fluctuations(model, Positions, OrderParameters[0], OrderParameters[1], OrderParameters[2]);
		Avg_s6 += Flucts[0];
		Chi6x += Flucts[1];
		Chi6y += Flucts[2];
		Chi6 += pow(OrderParameters[0], 2);
		int HexaticBin = int(OrderParameters[0] / HistX6.bin_width);
		X6OrderParameter[HexaticBin] += 1;
	}
	/*Chi6 *= (double(model.N) / model.num_cfgs);	// integer divided by another integer is not generally an integer. use double to maintain precision
	Avg_s6 /= (model.N * model.num_cfgs);
	Chi6x /= (model.N * model.num_cfgs);
	Chi6y /= (model.N * model.num_cfgs);
	Chi6x -= (model.N * pow(Avg_s6, 2));*/
	Chi6 *= (double(model.N) / TemporaryNumCfg);	// integer divided by another integer is not generally an integer. use double to maintain precision
	Avg_s6 /= (model.N * TemporaryNumCfg);
	Chi6x /= (model.N * TemporaryNumCfg);
	Chi6y /= (model.N * TemporaryNumCfg);
	Chi6x -= (model.N * pow(Avg_s6, 2));
	std::ofstream FitParams("FluctFitGauss" + tag + ".data");
	FitParams << Chi6 << '\t' << Avg_s6 << '\t' << Chi6x << '\t' << Chi6y << '\n';
	//PrintHistogram(model, HistX6, OutputFiles[0], X6OrderParameter, "prob_density");
	// cat Trajectory51 | ./ComputeFluctuations.exe
	return 0;
}

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model) {	// read in a configuration.
	std::vector<std::vector<Point>> PositionVector(model.N, std::vector<Point>(model.NA));
	double f, g;
	int throwaway1, throwaway2;
	for (int i = 0; i < model.N; i++) {
		for (int a = 0; a < model.NA; a++) {
			//std::cin >> throwaway1 >> throwaway2 >> f >> g;
			std::cin >> f >> g;
			//std::cout << throwaway1 << '\t' << throwaway2 << '\t' << f << '\t' << g << '\n';
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
	int num_mol = Input.N, num_atom = 3, numcfg = Input.NumConfigs, steps_between_cfgs = 100;	//IMPORTANT: make sure this matches with any other programs you use. I should really put this in a header.
	double boxl = sqrt(num_mol / density);
	double boxinv = 1.0 / boxl;
	SimulationParameters model(density, temperature, boxl, boxinv, ts, numcfg, num_mol, num_atom, steps_between_cfgs);
	return model;
}

double ComputeCosSixTheta(double &CosTheta) {
	double square = CosTheta*CosTheta;
	double fourth = square*square;
	return 4 * (8 * fourth*square - 12 * fourth + 6 * square - 1) - 6 * square + 3;
}

double ComputeSinSixTheta(double &CosTheta, double &SinTheta) {
	return 6 * SinTheta * pow(CosTheta, 5) - 20 * pow(SinTheta, 3) * pow(CosTheta, 3) + 6 * pow(SinTheta, 5)*CosTheta;
}

std::vector<double> Fluctuations(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, double &magnitude_Phi6, double &Q1, double &Q2) {
	//accumulate components of this vector externally
	//compute average single molecule order parameter <s6>, x-component of fluctuation Chi_6x for a configuration.
	std::vector<double> tempvector(3, 0.0);
	double avg_s6 = 0, Chi_6x = 0, Chi_6y = 0;
	Point eigenvector((magnitude_Phi6 + Q1), Q2);	//eigenvector of {{x,y},{y,-x}}
	eigenvector = eigenvector / sqrt(eigenvector.dot(eigenvector));
	double theta0 = acos(eigenvector.x()) / 3.0; //eigenvector is cos(3x), sin(3x). we want cos(x), sin(x) which is the director.
	Point Director = Point(cos(theta0), sin(theta0));
	for (int j = 0; j < model.N; j++) {
		Point OmegaJ = (PositionVector[j][2] - PositionVector[j][1]);
		OmegaJ = OmegaJ / sqrt(OmegaJ.dot(OmegaJ));
		double cos_director_j = OmegaJ.dot(Director);	//cos(theta_j)
		double det_j = OmegaJ.x()*Director.y() - OmegaJ.y()*Director.x();
		double sin_director_j = sin(atan2(det_j, OmegaJ.dot(Director)));
		double cos_six_theta_director_j = ComputeCosSixTheta(cos_director_j);
		double sin_six_theta_director_j = ComputeSinSixTheta(cos_director_j, sin_director_j);
		avg_s6 += cos_six_theta_director_j;
		for (int k = 0; k < model.N; k++) {
			Point OmegaK = PositionVector[k][2] - PositionVector[k][1];
			OmegaK = OmegaK / sqrt(OmegaK.dot(OmegaK));
			double cos_director_k = OmegaK.dot(Director);
			double det_k = OmegaK.x()*Director.y() - OmegaK.y()*Director.x();
			double sin_director_k = sin(atan2(det_k, OmegaK.dot(Director)));
			double cos_six_theta_director_k = ComputeCosSixTheta(cos_director_k);
			double sin_six_theta_director_k = ComputeSinSixTheta(cos_director_k, sin_director_k);
			Chi_6x += cos_six_theta_director_j*cos_six_theta_director_k;
			Chi_6y += sin_six_theta_director_j*sin_six_theta_director_k;
		}
	}
	tempvector = {avg_s6, Chi_6x, Chi_6y};
	return tempvector;
}

std::vector<double> CalculatePhi6(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector) {
	//avg1 is the average of the single molcule order parameter with respect to the director, avg2 is the x component of the susceptibility
	double phi6 = 0, Q1 = 0, Q2 = 0;
	for (int j = 0; j < model.N; j++) {
		Point OmegaJ = (PositionVector[j][2] - PositionVector[j][1]);
		OmegaJ = OmegaJ / sqrt(OmegaJ.dot(OmegaJ));
		for (int k = 0; k < model.N; k++) {
			Point OmegaK = PositionVector[k][2] - PositionVector[k][1];
			OmegaK = OmegaK / sqrt(OmegaK.dot(OmegaK));
			double c = OmegaJ.dot(OmegaK);
			double value6 = ComputeCosSixTheta(c);
			phi6 += value6;
		}
		double c = OmegaJ.x(); //cos of angle between orientation and x-axis
		double c2 = c*c;
		double c4 = c2*c2;
		double value = 4 * (8 * c4*c2 - 12 * c4 + 6 * c2 - 1) - 6 * c2 + 3;
		double s = OmegaJ.y();
		double sin_2x = 2 * c*s;
		double sin_value = 3 * sin_2x - 4 * pow(sin_2x, 3);
		Q1 += value;
		Q2 += sin_value;

	}
	phi6 /= pow(model.N, 2);
	std::vector<double> tempvector = { sqrt(phi6), Q1/model.N, Q2/model.N };
	return tempvector;	//magnitude(phi_6)
}

void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode) {
	std::ofstream ofs(filename, std::ios_base::app);
	for (int bin = 0; bin < hist.num_bins; bin++) {
		double binvalue = double(bin) * hist.bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / 1000 << '\n';
		//if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / 5000 << '\n';
		else if (mode == "prob_density") ofs << binvalue << '\t' << data[bin] / model.num_cfgs / hist.bin_width << '\n'; // if optional parameter is specified and greater than 1, divide bin count by bin width. This turns a histogram into a probability density
		else std::cout << "Bad parameter for PrintHistogram()";
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
