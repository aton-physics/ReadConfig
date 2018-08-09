#pragma once
#include <vector>
#include "ParameterClass.h"

void PrintPositions(std::string filename, std::string YesOrNo); //Yes for minimum image (good for plotting), No for no minimum image (good for preventing numerical headache)
void PrintVelocities(std::string filename);					//velocity file mystream should be unique wrt the position file mystream
void printEnergies(std::ofstream &mystream);
void velocityverlet(std::vector<double> &dsq);							//advance position and velocity for all particles, includes RATTLE
void velocityverlet_ts(int ts, std::vector<double> &dsq_rattle);
void InitialPositions(double &bond_angle, BondConstraints &bondlengths);
void InitialVelocities(double InitialTemperature);			//initialize velocities
void computeForces();
void run_for(int howlong, std::vector<double> &dsq_rattle);
void chill(double howchill, int EquilibrateTime, std::vector<double> &dsq_rattle);
void cool_past(double &howcold, std::vector<double> &dsq_rattle);
void pack(double &newdensity, std::vector<double> &dsq_rattle);	//increase density by reducing molecular separations and reducing the box volume
InputParameter ReadInput();
void GenerateMeltedConfiguration(InputParameter &Parameters, BondConstraints &BondLengths);
void InputGen();
void TargetTemperature(double &target, std::vector<double> &dsq_rattle);			//hit target temperature exactly
void AssignInitialConditions(InputParameter &Parameters, double hotter);
double AngMomentum(); //calculate total angular momentum
double AngVelocity(); //Calculate instantaneous angular velocity using the Inertia Tensor
std::ifstream& GotoLine(std::ifstream& file, unsigned int num); //skip to a certain line in data file
