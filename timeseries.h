#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

using namespace std;

class Params
{
public:
	Params(double YY, double YC, double YS, double CC, double SC, double SS) : yy(YY), yc(YC),ys(YS),cc(CC),sc(SC),ss(SS){}
	const double yy, yc, ys, cc, sc, ss;
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
	size_t iter()   const {return _iter;}
	int status()    const {return _status;}
	
private:

	inline int solve(const vector<double> &y, vector<double> &x)
	{
		gsl_vector_const_view yVec = gsl_vector_const_view_array(y.data(), length());
		gsl_vector_view xVec = gsl_vector_view_array(x.data(), length());
		gsl_linalg_LU_solve (mat, per, &yVec.vector, &xVec.vector);
		for (int i = 0; i < length(); i++) x[i] *= lambda();
	}
	static inline double f_constLambda   (const gsl_vector *x, void *params) 
	{
		static gsl_vector *dummy=gsl_vector_alloc(2);
		return reinterpret_cast<TimeSeries*>(params)->fdf(gsl_vector_get(x,0), gsl_vector_get(x,1), dummy);
	}
	static inline void   df_constLambda  (const gsl_vector *x, void *params,            gsl_vector *df)
	{
		reinterpret_cast<TimeSeries*>(params)->fdf(gsl_vector_get(x,0), gsl_vector_get(x,1), df);
	}
	static inline void   fdf_constLambda (const gsl_vector *x, void *params, double *f, gsl_vector *df)
	{
		*f = reinterpret_cast<TimeSeries*>(params)->fdf(gsl_vector_get(x,0), gsl_vector_get(x,1), df);
	}
	
	void fit_constLambda();
	double fdf (const double Phi, const double Alpha, gsl_vector *df);
	void test();
	void loadDummyFile();        // prepareonly x,y
	void resizeVectors();        // depends only on length of x
	void getBendParameters();    // depends only on x
	void getBendPenaltyMatrix(); // depends only on x
	void getCompleteMatrix();    // depends on x and lambda
	void getSolutionVectors();   // depends on x, alpha, phi
	Params getParams();
	Params getdParams();
	void solveSpline();
	void writeFiles();
	void getPhiAlpha();
	double bendPenalty(const vector<double> &curve);
	double rmsd       (const vector<double> &curve);
	vector<double> x; // time vector
	vector<double> y; // original data
	vector<double> spline; // spline when not using oscillations
	vector<double> trend;  // spline when using oscillations
	vector<double> ecos, esin, ecost, esint; // exp(ax)*cos(x)  neeced for calculations
	vector<double> y_c, y_s, y_ct, y_st;
	vector<vector<double> > bendParam;
	//gsl_vector_view yView, ecosView, esinView, yyView, ycView, ysView;
	gsl_matrix *bendMat, *mat;
	gsl_permutation *per; // for matrix decomposition
	double _lambda, _A, _alpha, _phi, _theta;
	size_t _iter;
	int _status;
	static const double phi_start, alpha_start, stepSize, tolerance;
	static const int giveUpTime;
	static const int maxIter;

};

#endif // TIMESERIES_H
