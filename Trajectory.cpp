#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>
#include <assert.h>
#include <limits> //numeric_limits
#include <sys/stat.h> // mkdir
#include "headers/FnDeclarations.h"
#include "headers/ParameterClass.h" // store all the input parameters
#include "unistd.h"		//a POSIX API - getopt

//***********************************1/27/18**********************************
//NA = 3, NB = 3 -> rigid nonlinear triatomic molecule.
//Held bonds fixed using RATTLE algorithm, solving equations of motion while using (linear approx. of) constraint forces to keep bond length constraints rigid
//This code is intended for study of diffusion in glass formers/supercooled liquids, and for the application of geodesics. constant NVE
//Temperature is related to kinetic energy via : T = 2K/3N. Average energy is kBT/2 per degree of freedom. Molecules have 2 translational and 1 rotational degree of freedom, so 3N/2kbT is avg energy.
// 8/9/18 Just refactored the program - size cut from 1100 lines to 653. Introduced a few objects, halved the number of global variables. TODO: Decouple the functions, get rid of more global variables.
//SHAKE, RATTLE have degenerate solutions for 180 degree bond angles, so linear molecules consisting of 3 or more atoms cannot be handled.

const int n = 10;									//2n^2 atoms
const int N = 2 * n * n;							//# molecules in simulation
const int NA = 3;									//# atoms per molecule
const int NB = 3;									//# bonds per molecule
const double m[NA] = { 1.0, 1.0, 1.0 };	//mass of atoms - in reduced units, sum(m[NA]) = 1.0 if the molecule is 1 reduced mass unit
const double dt = 0.001;							//PARAMETER: step size, reduced time scale: sqrt(m* sigma^2 / epsilon) = 1 time unit
double rhoBar = 0.10;								//PARAMETER: Initial reduced number density, N * sigma^2
double boxl = sqrt(N / rhoBar);						//periodic box length - this is modified by pack, so can leave it global.
double boxinv = 1.0 / boxl;
double boxl2 = boxl * 0.5;
const double tol = 1.0e-7;							//RATTLE algorithm tolerances - 1.0e-10 worked up until 126 degree bond angle (with maxit=50)
double timekeeping = 0;
double rx[N][NA], ry[N][NA], vx[N][NA], vy[N][NA];
double fx[N][NA], fy[N][NA];
double V, K, VC;							//Potential, Kinetic, Shifted Force Potential energies
std::string tag = "FailedToAssignTag";

int main(int argc, char ** argv) {
	std::ios::sync_with_stdio(false);	// speed up C++ I/O by disabling C I/O
	double bondlength1 = 1.0, bondlength2 = 1.0;
	int c = 1;
	const char optstring[] = "+ics"; //valid options
	while (-1 != (c = getopt(argc, argv, optstring))) {
		switch (c) {
		case 'i': {	//create input files
			InputGen();	// just run ./executable -i, don't need to qsub this one.
			break;
		}
		case 'c': {	//read input files, ouput a melted configuration. can submit with qsub -q all.q -t begin-end:increment -tc numberofcores ./script.sh 
			InputParameter parameter = ReadInput();
			BondConstraints constraints(bondlength1, bondlength2, parameter.BondAngle);
			std::vector<double> dsq = constraints.dsq();
			GenerateMeltedConfiguration(parameter, constraints); //create a melted configuration according to the specified parameters 
			break;
			// if you're getting a complaint about invalid strings, make sure to use qsub -t option
		}
		case 's': {	//run a trajectory while saving configurations and the average potential energy per particle
			InputParameter parameter = ReadInput();
			BondConstraints constraints(bondlength1, bondlength2, parameter.BondAngle);
			std::vector<double> dsq = constraints.dsq();
			AssignInitialConditions(parameter, 0.5);	// start hotter than intended - allow minor annealing
			TargetTemperature(parameter.temperature, dsq);	// cool to desired temperature
			int NumConfigs = parameter.NumConfigs;	//#of saved configurations. # configurations * skiptime + waiting time = total run length.
			int skiptime = 10;	//skip this many time steps between configurations
			run_for(3000, dsq);	//pass most relaxation times
			std::vector<double> PositionX, PositionY;
			double PE = 0.0;
			//std::ofstream ofs("CheckEnergy.data");
			for (int n = 0; n < NumConfigs; n++) {
				for (int i = 0; i < N; i++) {
					for (int a = 0; a < NA; a++) {	//store positions
						PositionX.push_back(rx[i][a]);
						PositionY.push_back(ry[i][a]);
					}
				}
				velocityverlet_ts(skiptime, dsq);
				PE += V;
				if ((n+1) % 10000 == 0 && n != 0) {	// every 10000 configurations, print and clear the positions. Reduces vmem usage, prevents bad_alloc.
					for (std::vector<int>::size_type i = 0; i < PositionX.size(); i++) {
						std::cout << PositionX[i] << '\t' << PositionY[i] << '\n';
					}
					PositionX.clear();
					PositionY.clear();
					PositionX.shrink_to_fit();
					PositionY.shrink_to_fit();
				}
				//printEnergies(ofs);
			}
			std::ofstream energyfile("Trajectory/" + std::to_string(int(parameter.BondAngle)) + "degrees/Energy" + tag + ".data");
			energyfile << parameter.temperature << '\t' << PE / double(N * NA) / NumConfigs << '\n';	//potential energy per particle, average over configurations
			break;
		}
		case '?': {
			std::cout << "Bad argument; exiting!";
			exit(1);
			break;
		}
		default: {
			std::cout << "Internal error; exiting!";
			exit(1);
			break;
		}
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

double AngMomentum() {	//calculate total angular momentum. Analytically, Ltotal = sum of system center of mass angular momentum about origin + sum of component angular momentums about the origin
	double TotalAngularMomentum = 0;
	double positionx = 0, positiony = 0, momentumx = 0, momentumy = 0;
	for (int i = 0; i < N; i++) {	//find center of mass, then find velocity of center of mass
		for (int a = 0; a < NA; a++) {
			TotalAngularMomentum += (rx[i][a] * vy[i][a] - ry[i][a] * vx[i][a]);	//component L
			positionx += rx[i][a];	//mass = 1, omitted.
			positiony += ry[i][a];
			momentumx += vx[i][a];	//mass = 1, omitted.
			momentumy += vy[i][a];
		}
	}
	positionx /= N * NA;	//mass = 1, omitted.
	positiony /= N * NA;	//to get the center of mass, take sum(mass * position) / system mass
	TotalAngularMomentum += (positionx*momentumy - positiony*momentumx);
	return TotalAngularMomentum;
}

double AngVelocity() { //Calculate instantaneous angular velocity using the Inertia Tensor
	//Inertia tensor in 2-d looks like (I11 I12 0) (I21 I22 0) (0 0 I33). Angular velocity (system) in 2-d looks like (0 0 w_z). For Angular Momentum (system) L_z, then w_z = L_z/I33.
	double AngVel = 0.0;
	for (int p = 0; p < 1; p++) {	//what does this loop structure even do here?
		double I33 = 0.0;
		for (int i = 0; i < N; i++) {	//find angular velocity
			for (int a = 0; a < NA; a++) {
				I33 += pow(rx[i][a], 2) + pow(ry[i][a], 2);
			}
		}
		AngVel = AngMomentum() / I33;
	}
	return AngVel;
}

void TargetTemperature(double &target, std::vector<double> &dsq_rattle) {	//hit target temperature with reasonable precision
	double AverageTemperature;
	std::vector<double> storetemperature;
	for (int i = 0; i < 20000; i++){
		velocityverlet_ts(5, dsq_rattle);
		storetemperature.push_back(2 * K / (3 * N));
	}
	AverageTemperature = std::accumulate(storetemperature.begin(), storetemperature.end(), 0.0) / storetemperature.size();
	storetemperature.clear();
	while (AverageTemperature < target - 0.003 || AverageTemperature > target + 0.003) {
		double a = 0;
		a = target / AverageTemperature;
		if (abs(AverageTemperature - target) < 0.02) {
			a = (1 + a) / 2; //prevent overshooting by reducing the effect of chill()
		}
		chill(a, 20, dsq_rattle);
		for (int i = 0; i < 20000; i++){
			velocityverlet_ts(4, dsq_rattle);
			storetemperature.push_back(2 * K / (3 * N));
		}
		AverageTemperature = std::accumulate(storetemperature.begin(), storetemperature.end(), 0.0) / storetemperature.size();
		storetemperature.clear();
	}
}

void AssignInitialConditions(InputParameter &Parameters, double hotter) {	//assign initial positions and velocities using a chosen temperature and an existing configuration
	std::ifstream inStream;
	inStream.exceptions(std::ifstream::failbit);
	try {
		inStream.open(Parameters.MeltCfg);
	}
	catch (const std::exception& e) {
		std::ostringstream msg;
		msg << "Opening file '" << Parameters.MeltCfg
			<< "' failed, it either doesn't exist or is not accessible.";
		throw std::runtime_error(msg.str());
	}
	std::ifstream myfile(Parameters.MeltCfg);
	double r, s;
	for (int i = 0; i < N; i++) {
		for (int a = 0; a < NA; a++) {
			myfile >> r >> s;
			rx[i][a] = r;
			ry[i][a] = s;
		}
	}
	rhoBar = Parameters.density;
	boxl = sqrt(N / rhoBar);		//fix box lengths, periodic boundary conditions
	boxinv = 1.0 / boxl;
	boxl2 = boxl * 0.5;
	InitialVelocities(Parameters.temperature + hotter); 
}

void InputGen() {		//create a bunch of input files in subdirectory "inputfiles", also mkdir all the directories of interest
	//linecount -> bond angle -> density -> temperature -> #configurations -> N -> path/to/MeltedConfiguration
	std::vector<int> BondAngle = { 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165};
	std::vector<double>densityIn = { 0.20 }, temperatureIn = { 0.7, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2 };
	int NumConfigs = 2.9 * 100000;
	mkdir("Trajectory", ACCESSPERMS);
	mkdir("inputfiles", ACCESSPERMS);
	std::ofstream inputstream("inputfiles/input.data");
	int linenumber = 1;	
	for (std::vector<int>::size_type i = 0; i < BondAngle.size(); i++) {
		std::string angleout_filename = "Trajectory/" + std::to_string(BondAngle[i]) + "degrees";
		mkdir(angleout_filename.c_str(), ACCESSPERMS);
		for (std::vector<int>::size_type j = 0; j < densityIn.size() ; j++) {
			for (std::vector<int>::size_type k = 0; k < temperatureIn.size(); k++) {
				std::string angle_filename = "MeltedConfiguration/" + std::to_string(BondAngle[i]) + "degrees";	// MeltedConfiguration/75degrees/
				mkdir("MeltedConfiguration", ACCESSPERMS);
				mkdir(angle_filename.c_str(), ACCESSPERMS);
				std::string N_filename = "N" + std::to_string(N);
				mkdir((angle_filename + "/" + N_filename).c_str(), ACCESSPERMS);
				std::string density_filename = std::to_string(densityIn[j]);
				mkdir((angle_filename + "/" + N_filename + "/").c_str(), ACCESSPERMS);
				inputstream << linenumber << '\t' << BondAngle[i] << '\t' << densityIn[j] << '\t' << temperatureIn[k] << '\t' << NumConfigs << '\t' << N << '\t' << angle_filename + "/" + N_filename + "/" + density_filename + ".data\n";
				linenumber++;
			}
		}
	}
}

void GenerateMeltedConfiguration(InputParameter &Parameters, BondConstraints &BondLengths) {	//initialize, pack to desired density, cool to desired T, relaxes, write melted positions a initialpositions file
	double AnnealTemperature = 3.0;				// start off at a high temperature so we can escape local minima and melt from the initial lattice configuration. Known as simulated annealing
	std::vector<double> dsq = BondLengths.dsq();
	InitialPositions(Parameters.BondAngle, BondLengths);
	InitialVelocities(AnnealTemperature);
	run_for(5, dsq);
	pack(Parameters.density, dsq);
	double NewTemperature = Parameters.temperature + 0.1;
	cool_past(NewTemperature, dsq);
	run_for(100, dsq);
	PrintPositions(Parameters.MeltCfg, "No");	// don't apply pbc while printing, pbc is only used to calculate forces or print human-readable configurations
}

InputParameter ReadInput() { // construct with all the input parameters
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

void PrintPositions(std::string filename, std::string YesOrNo) {	//Yes for minimum image (good for plotting), No for no minimum image (good for preventing numerical headache)
	std::ofstream positionstream(filename);
	if (YesOrNo == "Yes") {
		for (int i = 0; i < N; ++i) {
			for (int a = 0; a < NA; a++) {
				double rxia = rx[i][a];
				double ryia = ry[i][a];
				rxia -= static_cast<int>(rxia * boxinv + ((rxia >= 0.0) ? 0.5 : -0.5)) * boxl;	//minimum image algorithm, "efficient coding of the minimum image convention" by U K Deiters
				ryia -= static_cast<int>(ryia * boxinv + ((ryia >= 0.0) ? 0.5 : -0.5)) * boxl;
				positionstream <<  std::setw(12) << rxia << '\t' << std::setw(12) << ryia << '\n';
			}
		}
	}
	if (YesOrNo == "No") {
		for (int i = 0; i < N; ++i) {
			for (int a = 0; a < NA; a++) {
				double rxia = rx[i][a];
				double ryia = ry[i][a];
				positionstream << std::setw(12) << rxia << '\t' << std::setw(12) << ryia << '\n';
			}
		}
	}
}

void PrintVelocities(std::string filename) {					//velocity file mystream should be unique wrt the position file mystream
	std::ofstream velocitystream(filename);
	for (int i = 0; i < N; ++i) {
		for (int a = 0; a < NA; a++) {
			velocitystream << i << ' ' << a << '\t' << std::setw(12) << vx[i][a] << '\t' << std::setw(12) << vy[i][a] << '\n';
		}
	}
}

void printEnergies(std::ofstream &mystream) {//writes to mystream
	mystream << std::setw(12) << K << '\t' << std::setw(12) << V << '\t' << std::setw(12) << VC << '\t' << std::setw(12)
		<< std::setprecision(12) << K + VC << '\t' << std::setprecision(6) << std::setw(12) << 2 * K / (3 * N) << '\t'
		<< std::setw(12) << timekeeping << '\n';
}

void velocityverlet_ts(int ts, std::vector<double> &dsq_rattle) {			//call velocityverlet ts times without recording data
	for (int i = 0; i < ts; i++)
		velocityverlet(dsq_rattle);
}

void InitialPositions(double &bond_angle, BondConstraints &bondlengths) {
	//initialize positions
	double cell = boxl / n;
	double cell2 = cell * 0.5;
	double bond2 = bond_angle / 2;
	double L1 = bondlengths.L1;
	double L2 = bondlengths.L2;
	rx[0][0] = 0.0;						//sublattice A
	ry[0][0] = 0.0;
	rx[0][1] = -L1 * sin(bond2);
	ry[0][1] = -L1 * cos(bond2);
	rx[0][2] = L2 * sin(bond2);
	ry[0][2] = -L2 * cos(bond2);
	rx[1][0] = cell2;					//sublattice B
	ry[1][0] = cell2;
	rx[1][1] = cell2 - L1 * sin(bond2);
	ry[1][1] = cell2 - L1 * cos(bond2);
	rx[1][2] = cell2 + L2 * sin(bond2);
	ry[1][2] = cell2 - L2 * cos(bond2);
	int m = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int iref = 0; iref < 2; iref++) {	//repeatedly increment copies of sublattice A and B in each direction
				rx[iref + m][0] = rx[iref][0] + cell * i;
				rx[iref + m][1] = rx[iref][1] + cell * i;
				rx[iref + m][2] = rx[iref][2] + cell * i;
				ry[iref + m][0] = ry[iref][0] + cell * j;
				ry[iref + m][1] = ry[iref][1] + cell * j;
				ry[iref + m][2] = ry[iref][2] + cell * j;
			}
			m += 2;						//allows me to iterate over m points
		}
	}
	for (int i = 0; i < N; i++) {		//shift center of box to the origin
		for (int a = 0; a < NA; a++) {
			rx[i][a] -= boxl2;
			ry[i][a] -= boxl2;
		}
	}
}

void InitialVelocities(double InitialTemperature) {				//initialize velocities
	double sumx = 0.0, sumy = 0.0, ang_vel = 0.0;		//algorithm from A&T supplementary code. Zeroing angular momentum comes from A&T edition 2 supplementary code (on Allen's github).
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<> gauss(0, 1);
	double rtemp = sqrt(InitialTemperature);
	for (int i = 0; i < N; i++) {
		for (int a = 0; a < NA; a++) {
			vx[i][a] = rtemp * gauss(gen);
			vy[i][a] = rtemp * gauss(gen);
			sumx += vx[i][a];
			sumy += vy[i][a];
		}
	}
	sumx = sumx / (N * NA);
	sumy = sumy / (N * NA);
	for (int i = 0; i < N; i++) {
		for (int a = 0; a < NA; a++) {
			vx[i][a] -= sumx;
			vy[i][a] -= sumy;
		}
	}
	ang_vel = AngVelocity();
	for (int i = 0; i < N; i++) {	//(w cross r) gives a term that looks like (a b 0). Subtracting this value from velocities will zero the Angular momentum.
		for (int a = 0; a < NA; a++) {
			vx[i][a] += ang_vel * ry[i][a];
			vy[i][a] -= ang_vel * rx[i][a];
		}
	}
}

void computeForces() {								//input _ output _
	double rx_jakb, ry_jakb, r2, r6, r12, v_jakb, f_jakb, rxja, ryja;
	int ncut;										//counts how many rsq < rcutsq
	double fx_jakb, fy_jakb, fxj, fyj;
	double rsq;
	const double rcut = 2.5;							//PARAMETER: force truncation cut off radius, 2.5*sigma by convention
	const double rcutsq = rcut * rcut;
	const double rc2 = 1 / rcutsq;
	const double rc6 = rc2 * rc2 * rc2;
	const double rc12 = rc6 * rc6;
	for (int i = 0; i < N; i++) {
		for (int a = 0; a < NA; a++) {
			fx[i][a] = 0;							//zero atomic forces
			fy[i][a] = 0;
		}
	}
	V = 0.0;										//zero potential
	ncut = 0;
	for (int j = 0; j < N - 1; j++) {				//loop over all distinct atom pairs excluding those in same molecule
		for (int k = j + 1; k < N; k++) {
			for (int a = 0; a < NA; a++) {
				rxja = rx[j][a];
				ryja = ry[j][a];
				fxj = fx[j][a];
				fyj = fy[j][a];
				for (int b = 0; b < NA; b++) {
					rx_jakb = rxja - rx[k][b];	//pair separations
					ry_jakb = ryja - ry[k][b];
					rx_jakb -= static_cast<int>(rx_jakb * boxinv + ((rx_jakb >= 0.0) ? 0.5 : -0.5)) * boxl;	//minimum image algorithm, "efficient coding of the minimum image convention" by U K Deiters
					ry_jakb -= static_cast<int>(ry_jakb * boxinv + ((ry_jakb >= 0.0) ? 0.5 : -0.5)) * boxl;	//if ryij is more than half of the box length, minimum image vector then maps to an image particle **
					rsq = rx_jakb * rx_jakb + ry_jakb * ry_jakb;
					//assert(rsq != 0);
					if (rsq > rcutsq) continue;			//if minimum image sqd is larger than cutoff distance squared, skip interactions
					double r = sqrt(rsq);
					r2 = 1 / rsq;						//r^2, r^6, r^12
					r6 = r2 * r2 * r2;
					r12 = r6 * r6;
					v_jakb = r12 - r6;					
					V += v_jakb + (12 * rc12 - 6 * rc6) * (r / rcut - 1);
					f_jakb = (r12 + v_jakb) / rsq;		//main term of the force routine, (2r^-12 - r^-6)/r^2
					f_jakb -= (2 * rc12 - rc6) / (rcut * r);
					fx_jakb = f_jakb * rx_jakb;
					fy_jakb = f_jakb * ry_jakb;
					fxj += fx_jakb;						//accumulate forces
					fyj += fy_jakb;
					fx[k][b] -= fx_jakb;				//Newton's 3rd law, shown by r(2->1) = - r(1->2)
					fy[k][b] -= fy_jakb;
					ncut += 1;
				}
				fx[j][a] = fxj;
				fy[j][a] = fyj;
			}
		}
	}
	v_jakb = rc12 - rc6;
	VC = V - (double(ncut) * v_jakb);
	for (int i = 0; i < N; i++) {				//Multiply results by energy factors
		for (int a = 0; a < NA; a++) {
			fx[i][a] *= 24;
			fy[i][a] *= 24;
		}
	}
	V *= 4.0;
	VC *= 4.0;
}

void velocityverlet(std::vector<double> &dsq) {								//advance position and velocity for all particles, includes RATTLE
	const double dt2 = dt / 2.0;
	const double dtsq2 = dt * dt2;
	double tol2 = tol * 2.0;
	double rptol = tol;
	double rma, rmb;								//reciprocal of masses
	double rxi[NA], ryi[NA], pxi[NA], pyi[NA], vxi[NA], vyi[NA];
	double axia, ayia, pxab, pyab, pabsq, rabsq, diffsq, rxab, ryab, rpab, gab, dx, dy, vxab, vyab, rvab;
	bool moving[NA], moved[NA], done;
	int b, it, maxit = 500000;							//maxit = maximum allowed iterations
	computeForces();
	for (int i = 0; i < N; i++) {
		for (int a = 0; a < NA; a++) {				//advance positions full-step, velocities half-step (unconstrained)
			axia = fx[i][a] / m[a];					//
			ayia = fy[i][a] / m[a];
			rxi[a] = rx[i][a];
			ryi[a] = ry[i][a];
			pxi[a] = rx[i][a] + dt * vx[i][a] + dtsq2 * axia;
			pyi[a] = ry[i][a] + dt * vy[i][a] + dtsq2 * ayia;
			vxi[a] = vx[i][a] + dt2 * axia;
			vyi[a] = vy[i][a] + dt2 * ayia;
			moving[a] = false;
			moved[a] = true;
		}
		it = 0;
		done = false;
		while (!done && it <= maxit) {				//begin iterative loop
			done = true;
			for (int a = 0; a < NB; a++) {
				b = a + 1;
				if (b == NA) b = 0;
				if (moved[a] || moved[b]) {
					pxab = pxi[a] - pxi[b];
					pyab = pyi[a] - pyi[b];
					pabsq = pxab * pxab + pyab * pyab;
					rabsq = dsq[a];
					diffsq = rabsq - pabsq;
					if (std::abs(diffsq) > (rabsq * tol2)) {			// if uncorrected bond length is different from the bond length constraint up to the square root of twice the tolerance
						rxab = rxi[a] - rxi[b];
						ryab = ryi[a] - ryi[b];
						rpab = rxab * pxab + ryab * pyab;				//"s" per Andersen [1982]
						assert(rpab >= rabsq * rptol);					//make sure bond vector doesn't turn too much
						rma = 1.0 / m[a];
						rmb = 1.0 / m[b];
						gab = diffsq / (2.0 * (rma + rmb) * rpab);
						dx = rxab * gab;
						dy = ryab * gab;
						pxi[a] += rma * dx;								//p. 33 Andersen [1982]
						pyi[a] += rma * dy;
						pxi[b] -= rmb * dx;
						pyi[b] -= rmb * dy;
						dx /= dt;
						dy /= dt;
						vxi[a] += rma * dx;								//p. 32 Andersen [1982]
						vyi[a] += rma * dy;
						vxi[b] -= rmb * dx;
						vyi[b] -= rmb * dy;
						moving[a] = true;
						moving[b] = true;
						done = false;
					}
				}
			}
			for (int a = 0; a < NA; a++) {
				moved[a] = moving[a];
				moving[a] = false;
			}
			it += 1;
		}													//end of iterative loop
															//if (!done) std::cout << "molecule " << i << '\n';
		assert(done);										//too many constraint iterations
		//std::cout << "passed the tolerance, iteration count = " << it << '\n';
		for (int a = 0; a < NA; a++) {						//store new values
			rx[i][a] = pxi[a];
			ry[i][a] = pyi[a];
			vx[i][a] = vxi[a];
			vy[i][a] = vyi[a];
		}
	}														//end loop over molecules
	computeForces();										//begin 2nd stage of velocity verlet with constraints
	K = 0.0;
	for (int i = 0; i < N; i++) {
		for (int a = 0; a < NA; a++) {						//advance velocities of all atoms within a molecule a full step without constraints
			rxi[a] = rx[i][a];
			ryi[a] = ry[i][a];
			vxi[a] = vx[i][a] + dt2 * fx[i][a] / m[a];
			vyi[a] = vy[i][a] + dt2 * fy[i][a] / m[a];
			moving[a] = false;
			moved[a] = true;
		}
		it = 0;
		done = false;
		while (!done && it <= maxit) {
			done = true;
			for (int a = 0; a < NB; a++) {
				b = a + 1;
				if (b == NA) b = 0;
				if (moved[a] || moved[b]) {
					vxab = vxi[a] - vxi[b];
					vyab = vyi[a] - vyi[b];
					rxab = rxi[a] - rxi[b];
					ryab = ryi[a] - ryi[b];
					rvab = rxab * vxab + ryab * vyab;
					rma = 1.0 / m[a];
					rmb = 1.0 / m[b];
					gab = -rvab / ((rma + rmb) * dsq[a]);
					if (std::abs(gab) > tol) {
						dx = rxab * gab;
						dy = ryab * gab;
						vxi[a] += rma * dx;
						vyi[a] += rma * dy;
						vxi[b] -= rmb * dx;
						vyi[b] -= rmb * dy;
						moving[a] = true;
						moving[b] = true;
						done = false;
					}
				}
			}
			for (int a = 0; a < NA; a++) {
				moved[a] = moving[a];
				moving[a] = false;
			}
			it += 1;
		}													//end of iterative loop
		assert(done);
		for (int a = 0; a < NA; a++) {
			vx[i][a] = vxi[a];						
			vy[i][a] = vyi[a];
			K += m[a] * (vxi[a] * vxi[a] + vyi[a] * vyi[a]); //only use if m[a] is not identical for all sites
															 //K += (vxi[a] * vxi[a] + vyi[a] + vyi[a]);
		}
	}														//end loop over molecules for stage two
	K *= 0.5;
	timekeeping += dt;
}

void run_for(int howlong, std::vector<double> &dsq_rattle) {
	double num_calls = howlong / dt;
	for (int i = 0; i < num_calls; i++) {
		velocityverlet(dsq_rattle);
	}
}

void chill(double howchill, int EquilibrateTime, std::vector<double> &dsq_rattle) {	//lowers temperature slowly. howchill is a multiplier on the temperature.
	for (int i = 0; i < N; i++) {					//tBar = vBar^2 / 3N, hence frac * tBar = (vBar * sqrt(frac))^2 / 3N
		for (int a = 0; a < NA; a++) {
			vx[i][a] *= sqrt(howchill);
			vy[i][a] *= sqrt(howchill);
		}
	}
	run_for(EquilibrateTime, dsq_rattle);
}

void cool_past(double &howcold, std::vector<double> &dsq_rattle) {	//lower the temperature below howcold
	std::vector<double> storetemperature;
	double temperature = 9000;
	int i = 0;
	while (temperature > howcold) {
		chill(0.99, 3, dsq_rattle);
		storetemperature.push_back(2 * K / (3 * N));
		if (i > 10 && i % 5 == 0) {
			temperature = std::accumulate(storetemperature.begin(), storetemperature.end(), 0.0) / storetemperature.size(); //computes average value of storetemperature's contents
			storetemperature.clear();
		}//reassign temperature every now and then
		i++;
	}
}

void pack(double &newdensity, std::vector<double> &dsq_rattle) {	//increase density by reducing molecular separations and reducing the box volume
	double new_boxl_sq = N / newdensity;
	while (boxl*boxl > new_boxl_sq) { //Shrink (large steps), will overshoot target density by a little
		for (int i = 0; i < N; i++) {
			for (int a = 0; a < NA; a++) {
				rx[i][a] *= 0.995;
				ry[i][a] *= 0.995;
			}
		}
		boxl *= 0.995;
		boxinv = 1.0 / boxl;
		boxl2 = boxl * 0.5;
		rhoBar = N / (boxl * boxl);
		velocityverlet_ts(50, dsq_rattle);
	}
	while (boxl*boxl < new_boxl_sq) {	//overshot, so backtrack a little.
		for (int i = 0; i < N; i++) {
			for (int a = 0; a < NA; a++) {
				rx[i][a] *= 1.00001;
				ry[i][a] *= 1.00001;
			}
		}
		boxl *= 1.00001;
		boxinv = 1.0 / boxl;
		boxl2 = boxl * 0.5;
		rhoBar = N / (boxl * boxl);
		velocityverlet_ts(25, dsq_rattle);
	}
	while (boxl*boxl > new_boxl_sq) {
		for (int i = 0; i < N; i++) {
			for (int a = 0; a < NA; a++) {
				rx[i][a] *= 0.999999;
				ry[i][a] *= 0.999999;
			}
		}
		boxl *= 0.999999;
		boxinv = 1.0 / boxl;
		boxl2 = boxl * 0.5;
		rhoBar = N / (boxl * boxl);
		velocityverlet_ts(25, dsq_rattle);
	}
}
