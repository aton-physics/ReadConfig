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
#include "voro++/src/voro++.hh" 	//g++-7 -O3 -o ./voronoi.exe VoronoiLab.cpp voro++/src/libvoro++.a
#include <ctime> // clock function
#include <regex>

using namespace voro;
//***********************************8/6/2018**********************************
//This program takes in a potentially gigantic trajectory from molecular dynamics (gigabytes) and calculates various static quantities like the distribution of the order parameter, pair correlation function, also does video microscopy
//Meant to be called as a subprocess through a shell driver that also grabs task_id from the Sun Grid Engine. something like ./WriteTrajectoryToSTDOUT | xz > CompressTrajectory.xz, then unxz
//***1/7/2019****
//Starting to move away from gigantic trajectory files. Trying to do everything without physical storage, so the program should read through the same stream but directly from the mouth of the integrator.
std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model);
void WriteConfiguration(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, std::ostream &strm);
std::ifstream& GotoLine(std::ifstream& file, unsigned int num);
InputParameter ReadInput(std::string &tag);
SimulationParameters GetParams(InputParameter Input);
double CalculateBO6(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, int n);
std::vector<double> CalculateX6(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<Point>> &PositionVector);
void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode);
std::vector<std::string> FoldersToFiles(SimulationParameters &model, std::vector<std::string> &stringvector, std::string tag);
std::vector<std::vector<int>> create_list(SimulationParameters &model, std::vector<Point> &Configuration, double length);
std::vector<std::vector<int>> voronoi_list(SimulationParameters &model, std::ifstream &ifs);
std::vector<double> GetDirectorProjection(SimulationParameters &model, std::vector<std::vector<int>> &NeighborList, std::vector<std::vector<Point>> &Positions);

int main() {
	std::string tag = "FailedtoAssignTag";
	InputParameter parameter = ReadInput(tag);
	SimulationParameters model = GetParams(parameter);

	/////
	const double x_min = -model.boxl / 2.0, x_max = model.boxl / 2.0;
	const double y_min = -model.boxl / 2.0, y_max = model.boxl / 2.0;
	const double z_min = -0.5, z_max = 0.5;
	const int n_x = 17, n_y = 17, n_z = 1;
	////
	HistogramInfo HistX6(0.005, int(1 / 0.005));
	std::vector<double> OrderParameter(HistX6.num_bins, 0), X6OrderParameter(HistX6.num_bins, 0), BO6(HistX6.num_bins, 0);
	std::vector<Point> Cm_J(model.N);
	//std::ofstream ofs("testconfig.data");
	for (int n = 0; n < 1; n++) {
		//if (n % 100 == 0) std::cout << n << '\n';
		std::vector<std::vector<Point>> Positions = GetConfiguration(model);	//read in coordinates at a certain time step
		container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, true, true, true, 8);	//initialize container with periodic x, y, z 
		for (int i = 0; i < model.N; i++) {	// Store center of mass
			Cm_J[i] = (Positions[i][2] + Positions[i][1] + Positions[i][0]) / 3.0;
			Cm_J[i] = Cm_J[i].pbc(model.boxl, model.invboxl);
			con.put(i, Cm_J[i].x(), Cm_J[i].y(), 0);
		}
		//con.draw_particles("voroparticles21.gnu");
		con.draw_cells_gnuplot("vorooutline.gnu");
		con.print_custom("%i %n", "neighboroutput.data"); // neighbor list
		std::ifstream ifs("neighboroutput.data");
		std::vector<std::vector<int>> NeighborList = voronoi_list(model, ifs);	//got the neighbors now
		std::vector<double> Director_angle_wrt_x = GetDirectorProjection(model, NeighborList, Positions);
		con.print_custom("%P", "vertexoutput.data"); // vertex list
		std::ifstream vertexstream("vertexoutput.data");
		con.print_custom("%C", "centroidoutput.data"); // centroid list
		std::ifstream indexstream("centroidoutput.data");
		std::vector<Point> Centroid(model.N);
		for (int i = 0; i < model.N; i++) {
			double a, b, c;
			indexstream >> a >> b >> c;
			Centroid[i] = Point(a, b);
		}
		std::vector<std::vector<Point>> vertices(model.N);	//stores each atom's voronoi polygon's vertices
		for (int i = 0; i < model.N; i++) {
			std::string s;
			std::getline(vertexstream, s);
			std::istringstream ss{ std::regex_replace(s, std::regex{ R"(\(|\)|,)" }, " ") };	//ignore formatting of voro++
			std::vector<double> vertex_list{ std::istream_iterator<double>{ss}, std::istream_iterator<double>{} };
			ss.clear();
			vertex_list.erase(std::remove(vertex_list.begin(), vertex_list.end(), 0.5), vertex_list.end());
			vertex_list.erase(std::remove(vertex_list.begin(), vertex_list.end(), -0.5), vertex_list.end());	//remove z components
			for (int j = 0; j < vertex_list.size() - 1; j += 2) {
				vertices[i].push_back({ vertex_list[j], vertex_list[j + 1] });
			}
			Point center = Centroid[i];
			std::sort(vertices[i].begin(), vertices[i].end(), [&center](Point a, Point b) ->  bool {	//sort in ascending order with respect to angle clockwise from relative y-axis
				double angle_a = 0.0, angle_b = 0.0;
				angle_a = atan2(a.y()-center.y(), a.x()-center.x());
				angle_b = atan2(b.y()-center.y(), b.x()-center.x());
				return angle_a < angle_b;
			});
		}
		for (int i = 0; i < model.N; i++) {
			std::string line_to_print = " ";
			for (int j = 0; j < vertices[i].size(); j++) {	// 2 doubles constitute a vertex, so skip by twos
				line_to_print.append(std::to_string(vertices[i][j].x()));
				line_to_print.append(",");
				line_to_print.append(std::to_string(vertices[i][j].y()));
				line_to_print.append(" ");
				line_to_print.append("to ");
			}
			line_to_print.erase(line_to_print.length() - 3);	//erase the last "to"
			std::cout << "set object " + std::to_string(i + 1) + " polygon from " + line_to_print << '\n';
			std::cout << "set object " + std::to_string(i + 1) + " fc palette fraction " + std::to_string(Director_angle_wrt_x[i] / (3.14159265 / 3.0)) +
				" fillstyle solid 1.0 border \n";
		}
	}
	return 0;
	//cat NoPBC21.data | ./voronoi.exe > tryingnewoutput.data
	//g++-7 -O3 -o ./voronoi.exe VoronoiLab.cpp voro++/src/libvoro++.a
}

std::vector<double> GetDirectorProjection(SimulationParameters &model, std::vector<std::vector<int>> &NeighborList, std::vector<std::vector<Point>> &Positions) {
	std::vector<double> DirectorProjection(model.N);
	double pi_3 = 3.14159265 / 3.0;
	for (int i = 0; i < model.N; i++) {
		double frame_angle = 0.0;
		//compute i's orientation
		//save angle between this and x axis mod pi/3
		Point Orientation_i = Positions[i][0] - ((Positions[i][2] + Positions[i][1]) / 2);
		Orientation_i = Orientation_i / Orientation_i.dot(Orientation_i);
		double angle_i = atan2(Orientation_i.y(), Orientation_i.x());
		if (angle_i < 0) angle_i += 2 * 3.14159265;	//restrict angle domain to 0,2pi
		frame_angle = fmod(angle_i, pi_3);
		for (auto j : NeighborList[i]) {
			//compute i's neighbors' orientations 
			//save angle mod pi/3
			Point Orientation_j = Positions[j][0] - ((Positions[j][2] + Positions[j][1]) / 2);
			Orientation_j = Orientation_j / Orientation_j.dot(Orientation_j);
			double angle_j = atan2(Orientation_i.y(), Orientation_i.x());
			if (angle_j < 0) angle_j += 2 * 3.14159265;
			frame_angle += fmod(angle_j, pi_3);
		}
		frame_angle /= (NeighborList[i].size() + 1); // divide by: self + number of neighbors
		DirectorProjection[i] = frame_angle;
	}
	return DirectorProjection;
}

std::vector<std::vector<Point>> GetConfiguration(SimulationParameters &model) {	// read in a configuration.
	std::vector<std::vector<Point>> PositionVector(model.N, std::vector<Point>(model.NA));
	double f, g;
	int throwaway1, throwaway2;
	for (int i = 0; i < model.N; i++) {
		for (int a = 0; a < model.NA; a++) {
			std::cin >> f >> g;
			//std::cin >> f >> g;
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
											 /*char* taskID_string;
											 taskID_string = getenv("SGE_TASK_ID");
											 std::string task_id = taskID_string;
											 tag = task_id;*/
	tag = "27";
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

/*double CalculateBO6(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, int n) {
	double global_order_parameter = 0;
	double global_order_parameter_sin = 0;
	double neighborcutoff = 2.62;
	double neighborcutoffsq = neighborcutoff*neighborcutoff;
	std::vector<Point> Cm_J(model.N,0.0);
	for (int j = 0; j < model.N; j++) {
		Cm_J[j] = (PositionVector[j][2] + PositionVector[j][1] + PositionVector[j][0]) / 3.0;
	}
	for (int j = 0; j < model.N; j++) {
		double neighbor_counter = 0;
		double local_order_parameter = 0;
		double local_order_parameter_sin = 0;
		for (int k = 0; k < model.N; k++) {
			if (k == j) continue;
			Point OmegaJK = Cm_J[k] - Cm_J[j];
			Point separation = OmegaJK.pbc(model.boxl, model.invboxl);
			double sepsq = separation.dot(separation);
			if (sepsq > neighborcutoffsq) continue;
			double magnitude = OmegaJK.dot(OmegaJK);
			double c = OmegaJK.x()/sqrt(magnitude); // dotted with the x axis	WARNING: NOT NORMALIZED
			double c2 = c*c;	//to avoid a square root, I normalize at this step (never invoke an odd power of cosine anyway, no need to normalize c)
			double c4 = c2*c2;
			double value = 4 * (8 * c4*c2 - 12 * c4 + 6 * c2 - 1) - 6 * c2 + 3;
			double s = OmegaJK.y() / sqrt(magnitude);
			double sin_2x = 2 * c*s;
			double sin_value = 3 * sin_2x - 4 * pow(sin_2x, 3);
			//std::cout << c4 << '\t' << c2 << '\t' << c << '\n';
			local_order_parameter += value;
			local_order_parameter_sin += sin_value;
			assert(value <= 1);
			neighbor_counter++;
		}
		if (neighbor_counter == 0) {
			//std::cout << j << "has no neighbors" << '\n';
			//Cm_J[j].print(std::cout);
			//std::cout << "configuration " << n << '\n';
		}
		if (neighbor_counter == 0) continue;
		local_order_parameter /= neighbor_counter;	//this is phi_m(r_j)
		local_order_parameter_sin /= neighbor_counter;
		global_order_parameter += local_order_parameter;
		global_order_parameter_sin += local_order_parameter_sin;
		//std::cout << global_order_parameter << '\n';
	}
	global_order_parameter /= model.N;
	global_order_parameter_sin /= model.N;
	//std::cout << global_order_parameter << '\t' << global_order_parameter_sin << '\n';
	double final_value = (pow(global_order_parameter, 2) + pow(global_order_parameter_sin, 2));
	return final_value;
}*/

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

std::vector<std::vector<int>> voronoi_list(SimulationParameters &model, std::ifstream &ifs) {
	std::vector<std::vector<int>> NeighborList(model.N);
	int enterNumber;
	for (int i = 0; i < model.N; i++) {
		std::string line;
		getline(ifs, line);
		std::istringstream iss(line);
		std::vector<int> numbers;
		while (iss >> enterNumber)
		{
			numbers.push_back(enterNumber);
		}
		int focal_atom_index = numbers[0];
		for (auto j : numbers) if (j != focal_atom_index) NeighborList[focal_atom_index].push_back(j);
	}
	return NeighborList;
}

double CalculateBO6(SimulationParameters &model, std::vector<std::vector<Point>> &PositionVector, int n) {
	double global_order_parameter = 0;
	double global_order_parameter_sin = 0;
	double neighborcutoff = 2.62;
	double neighborcutoffsq = neighborcutoff*neighborcutoff;
	std::vector<Point> Cm_J(model.N, 0.0);
	for (int j = 0; j < model.N; j++) {
		Cm_J[j] = (PositionVector[j][2] + PositionVector[j][1] + PositionVector[j][0]) / 3.0;
	}
	for (int j = 0; j < model.N; j++) {
		double neighbor_counter = 0;
		double local_order_parameter = 0;
		double local_order_parameter_sin = 0;
		for (int k = 0; k < model.N; k++) {
			if (k == j) continue;
			Point OmegaJK = Cm_J[k] - Cm_J[j];
			Point separation = OmegaJK.pbc(model.boxl, model.invboxl);
			double sepsq = separation.dot(separation);
			if (sepsq > neighborcutoffsq) continue;
			double magnitude = OmegaJK.dot(OmegaJK);
			double c = OmegaJK.x() / sqrt(magnitude); // dotted with the x axis	WARNING: NOT NORMALIZED
			double c2 = c*c;	//to avoid a square root, I normalize at this step (never invoke an odd power of cosine anyway, no need to normalize c)
			double c4 = c2*c2;
			double value = 4 * (8 * c4*c2 - 12 * c4 + 6 * c2 - 1) - 6 * c2 + 3;
			double s = OmegaJK.y() / sqrt(magnitude);
			double sin_2x = 2 * c*s;
			double sin_value = 3 * sin_2x - 4 * pow(sin_2x, 3);
			//std::cout << c4 << '\t' << c2 << '\t' << c << '\n';
			local_order_parameter += value;
			local_order_parameter_sin += sin_value;
			assert(value <= 1);
			neighbor_counter++;
		}
		if (neighbor_counter == 0) {
			//std::cout << j << "has no neighbors" << '\n';
			//Cm_J[j].print(std::cout);
			//std::cout << "configuration " << n << '\n';
		}
		if (neighbor_counter == 0) continue;
		local_order_parameter /= neighbor_counter;	//this is phi_m(r_j)
		local_order_parameter_sin /= neighbor_counter;
		global_order_parameter += local_order_parameter;
		global_order_parameter_sin += local_order_parameter_sin;
		//std::cout << global_order_parameter << '\n';
	}
	global_order_parameter /= model.N;
	global_order_parameter_sin /= model.N;
	//std::cout << global_order_parameter << '\t' << global_order_parameter_sin << '\n';
	double final_value = (pow(global_order_parameter, 2) + pow(global_order_parameter_sin, 2));
	return final_value;
}

std::vector<double> CalculateX6(SimulationParameters &model, double &NormalizationSq, std::vector<std::vector<Point>> &PositionVector) {
	std::vector<double> tempvector = { 0,0,0,0,0,0,0 };
	double X1 = 0, X2 = 0, X3 = 0, X4 = 0, X5 = 0, X6 = 0;
	for (int j = 0; j < model.N; j++) {
		Point OmegaJ = (PositionVector[j][2] - PositionVector[j][1]);
		OmegaJ = OmegaJ / sqrt(OmegaJ.dot(OmegaJ));
		for (int k = 0; k < model.N; k++) {
			Point OmegaK = PositionVector[k][2] - PositionVector[k][1];
			OmegaK = OmegaK / sqrt(OmegaK.dot(OmegaK));
			double c = OmegaJ.dot(OmegaK);
			double c2 = c*c;
			double c3 = c2*c;
			double c4 = c2*c2;
			double c5 = c4*c;
			double value6 = 4 * (8 * c4*c2 - 12 * c4 + 6 * c2 - 1) - 6 * c2 + 3;
			X6 += value6;
			X1 += c;
			X2 += (2 * c2 - 1);
			X3 += (4 * c3 - 3 * c);
			X4 += (1 - 8 * c2 + 8 * c4);
			X5 += (5 * c - 20 * c3 + 16 * c5);
		}
	}
	X6 /= pow(model.N, 2);
	X1 /= pow(model.N, 2);
	X2 /= pow(model.N, 2);
	X3 /= pow(model.N, 2);
	X4 /= pow(model.N, 2);
	X5 /= pow(model.N, 2);
	tempvector[0] = X1;
	tempvector[1] = X2;
	tempvector[2] = X3;
	tempvector[3] = X4;
	tempvector[4] = X5;
	tempvector[5] = X6;
	return tempvector;
}

void PrintHistogram(SimulationParameters &model, HistogramInfo &hist, std::string filename, std::vector<double> &data, std::string mode) {
	std::ofstream ofs(filename, std::ios_base::app);
	for (int bin = 0; bin < hist.num_bins; bin++) {
		double binvalue = double(bin) * hist.bin_width;	// when the value is calculated (elsewhere), it is counted into the bin associated with the value after truncation. So this prints the bin's count as well as the bin's value.
		if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / 1000 << '\n';
		//if (mode == "histogram") ofs << binvalue << '\t' << data[bin] / 5000 << '\n';
		else if (mode == "prob_density") ofs << binvalue << '\t' << data[bin] / 1000 / hist.bin_width << '\n'; // if optional parameter is specified and greater than 1, divide bin count by bin width. This turns a histogram into a probability density
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
