#include "timeseries.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <math.h>
#include <gsl/gsl_linalg.h>

TimeSeries::TimeSeries()
{
	_A = 25000;
	_alpha = -0.4;
	_phi   = 6.2;
	_theta = 5;
	loadDummyFile();
	resizeVectors();
	getBendParameters();
	getBendPenaltyMatrix();

	_lambda = 1;
	for (int i = 0; i < 10; i++)
	{
		getCompleteMatrix();
		getSolutionVectors();
		solveSpline();
		_lambda *= 2;
	}
}

TimeSeries::~TimeSeries()
{
	gsl_matrix_free(bendMat);
	gsl_matrix_free(mat);
	gsl_permutation_free(per);
}

void TimeSeries::loadDummyFile()
{
	ifstream in("/home/stephan/work/chronostar/TopCount/bsp1.txt");
	if (!in) {cerr << "Could not open example file" << endl; exit(-1);}
	string buf;
	for (int i = 0; i < 2; i++) getline(in, buf, '\n'); // omit two outliers
	while (getline(in, buf, '\n'))
	{
		stringstream ss(buf);
		string tStr;
		getline(ss, tStr, '\t');
		x.push_back(atof(tStr.c_str()));
		getline(ss, tStr, '\t'); // temperature
		getline(ss, tStr, '\t'); // first data set
		y.push_back(atof(tStr.c_str()));
	}
	for (int i = 0; i < length(); i++) cout << x[i] << "\t"; cout << endl;
	for (int i = 0; i < length(); i++) cout << y[i] << "\t"; cout << endl;
}

void TimeSeries::resizeVectors()
{
	ecos.resize(length());
	esin.resize(length());
	y_y. resize(length());
	y_c. resize(length());
	y_s. resize(length());
	yView    = gsl_vector_view_array (y.   data(), length());
	ecosView = gsl_vector_view_array (ecos.data(), length());
	esinView = gsl_vector_view_array (esin.data(), length());
	yyView   = gsl_vector_view_array (y_y. data(), length());
	ycView   = gsl_vector_view_array (y_c. data(), length());
	ysView   = gsl_vector_view_array (y_s. data(), length());
	if (bendMat != NULL) gsl_matrix_free(bendMat);
	if (mat     != NULL) gsl_matrix_free(mat);
	if (per     != NULL) gsl_permutation_free(per);
	bendMat = gsl_matrix_calloc(length(), length());
	mat     = gsl_matrix_calloc(length(), length());
	per     = gsl_permutation_alloc (length());
	bendParam = vector<vector<double> >(length(), vector<double>(3));
}

void TimeSeries::getBendParameters()
{
	for (int midPoint = 1; midPoint < length()-1; midPoint++) // number of middle point in triplet
	{
		const double x12 = x[midPoint-1]-x[midPoint  ];
		const double x23 = x[midPoint  ]-x[midPoint+1];
		const double x31 = x[midPoint+1]-x[midPoint-1];
		bendParam[midPoint][0] = 1/x12/x31;
		bendParam[midPoint][1] = 1/x12/x23;
		bendParam[midPoint][2] = 1/x23/x31;
	}
}
void TimeSeries::getBendPenaltyMatrix()
{
	gsl_matrix_set_zero(bendMat);
	for (int midPoint = 1; midPoint < length()-1; midPoint++) // number of middle point in triplet
	{
		for (int i = -1; i < 2; i++) 
		{
			for (int k = -1; k < 2; k++) 
			{
				*gsl_matrix_ptr(bendMat, midPoint+i, midPoint+k) += bendParam[midPoint][i+1]*bendParam[midPoint][k+1];
			}
		}
	}
}

void TimeSeries::getCompleteMatrix()
{
	gsl_matrix_memcpy(mat, bendMat);
	gsl_matrix_scale(mat, lambda());
	for (int i = 0; i < length(); i++) *gsl_matrix_ptr(mat, i, i) += 1;
	int signum;
	gsl_linalg_LU_decomp (mat, per, &signum);
}

void TimeSeries::getSolutionVectors()
{
	for (int tNr = 0; tNr < length(); tNr++)
	{
		const double px = phi()*x[tNr];
		const double c = cos(px);
		const double s = sin(px);
		const double e = exp(alpha()*x[tNr]);
		ecos[tNr] = e*c;
		esin[tNr] = e*s;
	}
}

Params TimeSeries::getParams()
{
	double yc_rmsd = 0;
	double ys_rmsd = 0;
	double cc_rmsd = 0;
	double sc_rmsd = 0;
	double ss_rmsd = 0;
	// rmsd stuff
	for (int i = 0; i < length(); i++)
	{
		double Y = y   [i]-y_y[i];
		double S = esin[i]-y_s[i];
		double C = ecos[i]-y_c[i];
		yc_rmsd += Y*C;
		ys_rmsd += Y*S;
		cc_rmsd += C*C;
		sc_rmsd += S*C;
		ss_rmsd += S*S;
	}
	// bend stuff
	double yc_bend = 0;
	double ys_bend = 0;
	double cc_bend = 0;
	double sc_bend = 0;
	double ss_bend = 0;
	for (int i = 1; i < length()-1; i++)
	{
		double Y = bendParam[i][0]*y_y[i-1]+bendParam[i][1]*y_y[i]+bendParam[i][2]*y_y[i+1];
		double S = bendParam[i][0]*y_s[i-1]+bendParam[i][1]*y_s[i]+bendParam[i][2]*y_s[i+1];
		double C = bendParam[i][0]*y_c[i-1]+bendParam[i][1]*y_c[i]+bendParam[i][2]*y_c[i+1];
		yc_bend += Y*C;
		ys_bend += Y*S;
		cc_bend += C*C;
		sc_bend += S*C;
		ss_bend += S*S;
	}
	return Params(
		yc_rmsd + lambda()*yc_bend, 
		ys_rmsd + lambda()*ys_bend,
		cc_rmsd + lambda()*cc_bend, 
		sc_rmsd + lambda()*sc_bend,
		ss_rmsd + lambda()*ss_bend);
}

void TimeSeries::solveSpline()
{
	gsl_linalg_LU_solve (mat, per, &yView.   vector, &yyView.vector);
	gsl_linalg_LU_solve (mat, per, &ecosView.vector, &ycView.vector);
	gsl_linalg_LU_solve (mat, per, &esinView.vector, &ysView.vector);

	// y_y contains spline fit WITHOUT oscillation now
	for (int i = 0; i < length(); i++) cout << y_y[i] << "\t"; cout << endl;
	//
	// get amp, theta
	//
	Params p = getParams();
	_theta = atan2(p.sc*p.yc-p.cc*p.ys, p.sc*p.ys-p.ss*p.yc);
	cerr << "theta = " << theta() << endl;
	const double C = cos(theta());
	const double S = sin(theta());
	_A = (C*p.yc+S*p.ys) / (C*C*p.cc+2*S*C*p.sc+S*S*p.ss);
	cerr << "A = " << A() << endl;
	Params dP = getdParams();
	//
	//
	//
	vector<double> curve;
	double err;
	// output calculated spline, calculated cos
	curve = y_y;
	for (int i = 0; i < length(); i++) curve[i] -= A()*C*y_c[i]+A()*S*y_s[i];
	for (int i = 0; i < length(); i++) cout << curve[i] << "\t"; cout << endl;
	const double err2 = lambda()*bendPenalty(curve);
	for (int i = 0; i < length(); i++) curve[i] += A()*C*exp(alpha()*x[i])*cos(phi()*x[i])+A()*S*exp(alpha()*x[i])*sin(phi()*x[i]);
	for (int i = 0; i < length(); i++) cout << curve[i] << "\t"; cout << endl;
	const double err1 = fitDist(curve);
	cerr << A() << "\t" << err1 << "\t" << err2 << "\t" << err1+err2 << endl;
	exit(-1);
}

double TimeSeries::bendPenalty(const vector<double> &curve)
{
	double err = 0;
	for (int mid = 1; mid < length()-1; mid++)
	{
		const double a = curve[mid-1]*bendParam[mid][0] + curve[mid]*bendParam[mid][1] + curve[mid+1]*bendParam[mid][2];
		err += a*a;
	}
	return err;
}

double TimeSeries::fitDist(const std::vector< double >& curve)
{
	double err = 0;
	for (int i = 0; i < length(); i++) 
	{
		const double errPart = curve[i]-y[i];
		err += errPart*errPart;
	}
	return err;
}
