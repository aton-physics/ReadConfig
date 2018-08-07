#pragma once
#pragma once
#include <vector>
#include <cmath>
#include <string>
#include "Point.h"

//Algorithm adapted from http://web.archive.org/web/20140702223114/http://faculty.cs.niu.edu/~hutchins/csci230/best-fit.htm
struct Line {
	double _slope = 0, _yInt = 0;
	double getYforX(double x) {
		return _slope*x + _yInt;
	}
	// Construct line from points
	bool fitPoints(std::vector<Point> &pts) {
		int nPoints = pts.size();
		if (nPoints < 2) {
			// Fail: infinitely many lines passing through this single point (data is a vertical line i.e. same x coordinate)
			assert(nPoints < 2);
			return false;
		}
		double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
		for (int i = 0; i<nPoints; i++) {
			sumX += pts[i].x();
			sumY += pts[i].y();
			sumXY += pts[i].x() * pts[i].y();
			sumX2 += pts[i].x() * pts[i].x();
		}
		double xMean = sumX / nPoints;
		double yMean = sumY / nPoints;
		double denominator = sumX2 - sumX * xMean;
		// You can tune the eps (1e-7) below for your specific task
		if (std::fabs(denominator) < 1e-7) {
			// Fail: it seems a vertical line
			return false;
		}
		_slope = (sumXY - sumX * yMean) / denominator;
		_yInt = yMean - _slope * xMean;
		return true;
	}
	double get_sigma(std::vector<Point> &vect, Line myline) {
		double sigma2 = 0;
		int nPoints = vect.size();
		for (int i = 0; i < nPoints; i++) {
			double Y = vect[i].y();
			double X = vect[i].x();
			double Y1 = myline.getYforX(X);
			sigma2 += ((Y - Y1) * (Y - Y1));
		}
		sigma2 /= nPoints;
		return sqrt(sigma2);
	}
};

void get_regression(std::string inputfilename, std::string outputfilename, const int timescale, double wait) { //get regression from wait to timescale
	std::vector<Point> pointvector;
	Line myline;
	std::ifstream input_file(inputfilename);
	double x, y;
	std::string dummyLine;
	getline(input_file, dummyLine); // throw away header line
	while (input_file >> x >> y) {
		if (x >= wait) {			//throw away data points from asymptotic regime
			Point point = { x, y }; //collect points
			pointvector.push_back(point);
		}
	}
	myline.fitPoints(pointvector);
	std::ofstream regressionstream;
	regressionstream.open(outputfilename, std::ios_base::app);
	regressionstream  << myline._slope / 4 << '\t' <<  myline.get_sigma(pointvector, myline) << '\t' << myline._yInt;
};

void numerical_differentiation(std::string inputfilename, std::string outputfilename, double stepsize) {	//important: stepsize is however long the time is between f[i] and f[i+1]
	std::ifstream input_file(inputfilename);
	std::string dummyLine;
	getline(input_file, dummyLine); // throw away zeroes (zero time, zero displacement)
	std::vector<double> parsed(std::istream_iterator<double>(input_file), {});
	std::vector<Point> msdvector;
	for (int i = 0; i < int(parsed.size()); i += 2) {
		msdvector.push_back((parsed[i], parsed[i + 1]));
	}
	int counter = int(msdvector.size());
	//double x, y;
	/*while (input_file >> x >> y) {	//read, store all the ordered pairs
		Point point(x, y);
		msdvector.push_back(point);
		counter += 1;
	}*/
	std::cout << "counter is at " << counter << '\n';
	std::ofstream regressionstream(outputfilename);
	double derivative = 0.0;
	for (int i = 2; i < counter - 2; i++) {
		derivative = -msdvector[i + 2].y() + 8 * msdvector[i + 1].y() - 8 * msdvector[i - 1].y() + msdvector[i - 2].y();
		derivative /= (12 * stepsize);
		regressionstream <<  msdvector[i].x() << '\t' << derivative << '\n';
	}
};

void zero_slope_regression(std::string inputfilename, std::string outputfilename, double temperature, double wait) {	//b_fit = mean(y - mx) for (x,y) column vector
	std::vector<Point> diffvector;
	std::ifstream input_file(inputfilename);
	double x, y;
	double counter = 0;
	//std::string dummyLine;
	//getline(input_file, dummyLine); //only use if i want to parse and throw away some strings
	while (input_file >> x >> y) {	//read, store all the ordered pairs after wait
		if (x >= wait) {
			Point point = { x, y };
			diffvector.push_back(point);
			counter += 1;
		}
	}
	int nPoints = diffvector.size();
	double b_fit = 0;
	for (int i = 0; i < nPoints; i++) {	// accumulate the y values 
		b_fit += diffvector[i].y();
	}
	b_fit /= nPoints;			// b_fit = sum(y) / N
	double MSE = 0;
	for (int i = 0; i < nPoints; i++) {
		MSE += pow(b_fit - diffvector[i].y(), 2);	// squared residuals
	}
	MSE /= nPoints;
	MSE = sqrt(MSE);
	std::ofstream regressionstream(outputfilename, std::ios_base::app); // want to write a header to outputfile before calling this function
	regressionstream << temperature << '\t' << b_fit / 4 << '\t' << MSE / 4 << '\n';
};

/*void StdErrMean(std::string inputfilename, std::string inputfilename2, int NumFiles, int Resolution, int Index) {	//NumFiles specifies how many files I have, Resolution specifies how many points I chose to generate, Index labels what file to start from
	std::vector<double> MeanDiffusion(Resolution, 0.0);
	std::vector<double> StdErrDiffusion(Resolution, 0.0);
	std::vector<double> MeanTemperature(Resolution, 0.0); //this should be the same as the individual entries, redundancy is helpful for bugs
	std::vector<double> MSDT(Resolution, 0.0), MeanMSDT(Resolution, 0.0), TemperatureMSDT(Resolution, 0.0);
	double density;
	for (int i = Index; i < NumFiles + Index; i++) {
		double a, b, c;
		std::ifstream myfile(inputfilename + std::to_string(i + 1) + ".data");	//make sure inputfilename is "angell/angell"
		std::ifstream MSDTfile(inputfilename2 + std::to_string(i + 1) + ".data");
		int k = 0;
		while (MSDTfile >> a >> b >> c >> density) {	//msdt, time, T, density
			MSDT[k] += a;
			TemperatureMSDT[k] += c;
			k++;
		}
		int j = 0;
		while (myfile >> a >> b >> c >> density) {
			MeanDiffusion[j] += a;
			StdErrDiffusion[j] += b;
			MeanTemperature[j] += c;
			j++;
		}
	}
	std::ofstream file(inputfilename + std::to_string(Index / NumFiles + 1) + "Total.data");
	std::ofstream file2(inputfilename2 + std::to_string(Index / NumFiles + 1) + "Total.data");
	for (int i = 0; i < Resolution; i++) {
		MeanDiffusion[i] /= NumFiles;
		StdErrDiffusion[i] /= (NumFiles * sqrt(NumFiles)); //StdErr of the mean is average stdev divided by sqrt(# measurements)
		MeanTemperature[i] /= NumFiles;
		MeanMSDT[i] = MSDT[i] / NumFiles;
		file << MeanDiffusion[i] << '\t' << StdErrDiffusion[i] << '\t' << MeanTemperature[i] << '\t' << density << '\t' << MeanMSDT[i] << '\n';
	}
	std::vector<double> sigma(Resolution, 0.0);
	for (int i = 0; i < NumFiles; i++) {
		std::ifstream MSDTfile(inputfilename2 + std::to_string(i + 1) + ".data");
		double a, b, c;
		int k = 0;
		while (MSDTfile >> a >> b >> c >> density) {	//go through MSDTfile, grab MSDT values for each temperature
			MSDT[k] = a;
			k++;
		}
		for (int j = 0; j < Resolution; j++) {
			sigma[j] += pow(MSDT[j] - MeanMSDT[j], 2);		//accumulate residuals, first going down the file labeled i + 1 using indices j for each temperature
		}
	}
	for (int j = 0; j < Resolution; j++) {
		sigma[j] /= (NumFiles);							//stdev = sqrt(residual / n)
		sigma[j] = sqrt(sigma[j]);
		sigma[j] /= sqrt(NumFiles);						//stderr of mean = stdev / sqrt(n)
		TemperatureMSDT[j] /= NumFiles;
		file2 << MeanMSDT[j] << '\t' << sigma[j] << '\t' << TemperatureMSDT[j] << '\t' << density << '\n';
	}
};*/