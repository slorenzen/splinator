#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

using namespace std;

class Params
{
public:
	Params(double YC, double YS, double CC, double SC, double SS) : yc(YC),ys(YS),cc(CC),sc(SC),ss(SS){}
	const double yc, ys, cc, sc, ss;
};

class TimeSeries
{

public:

	TimeSeries();
	virtual ~TimeSeries();

	int length()    const {return x.size();}
	double lambda() const {return _lambda;}
	double A()      const {return _A;}
	double alpha()  const {return _alpha;}
	double phi()    const {return _phi;}
	double theta()  const {return _theta;}
	
private:

	void loadDummyFile();        // prepareonly x,y
	void resizeVectors();        // depends only on length of x
	void getBendParameters();    // depends only on x
	void getBendPenaltyMatrix(); // depends only on x
	void getCompleteMatrix();    // depends on x and lambda
	void getSolutionVectors();   // depends on x, alpha, phi
	Params getParams();
	Params getdParams();
	void solveSpline();
	double bendPenalty(const vector<double> &curve);
	double fitDist    (const vector<double> &curve);
	vector<double> x, y, ecos, esin, y_y, y_c, y_s;
	vector<vector<double> > bendParam;
	gsl_vector_view yView, ecosView, esinView, yyView, ycView, ysView;
	gsl_matrix *bendMat, *mat;
	gsl_permutation *per; // for matrix decomposition
	double _lambda, _A, _alpha, _phi, _theta;

};

#endif // TIMESERIES_H
