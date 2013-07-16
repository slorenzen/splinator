#include "timeseries.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

#define FVARS  double yy_rmsd = 0, yc_rmsd = 0, ys_rmsd = 0, cc_rmsd = 0, sc_rmsd = 0, ss_rmsd = 0; \
               double yy_bend = 0, yc_bend = 0, ys_bend = 0, cc_bend = 0, sc_bend = 0, ss_bend = 0; 
#define DFVARS double yct_rmsd = 0, yst_rmsd = 0, cct_rmsd = 0, sct_rmsd = 0, cst_rmsd = 0, sst_rmsd = 0; \
               double yct_bend = 0, yst_bend = 0, cct_bend = 0, sct_bend = 0, cst_bend = 0, sst_bend = 0;
#define FRMSD  double Y = spline[i]-y[i], S = esin[i]-y_s[i], C = ecos[i]-y_c[i]; \
               yy_rmsd += Y*Y; yc_rmsd += Y*C; ys_rmsd += Y*S; cc_rmsd += C*C; sc_rmsd += S*C; ss_rmsd += S*S;
#define DFRMSD double ST = esint[i]-y_st[i], CT = ecost[i]-y_ct[i]; \
               yct_rmsd += Y*CT; yst_rmsd += Y*ST; cct_rmsd += C*CT; sct_rmsd += S*CT; cst_rmsd += C*ST; sst_rmsd += S*ST;
#define FBEND  double Y  = bendParam[i][0]*spline[i-1]+bendParam[i][1]*spline[i]+bendParam[i][2]*spline[i+1]; \
               double S  = bendParam[i][0]*y_s[i-1]+bendParam[i][1]*y_s[i]+bendParam[i][2]*y_s[i+1]; \
               double C  = bendParam[i][0]*y_c[i-1]+bendParam[i][1]*y_c[i]+bendParam[i][2]*y_c[i+1]; \
               yy_bend += Y*Y; yc_bend -= Y*C; ys_bend -= Y*S; cc_bend += C*C; sc_bend += S*C; ss_bend += S*S; 
#define DFBEND double ST = bendParam[i][0]*y_st[i-1]+bendParam[i][1]*y_st[i]+bendParam[i][2]*y_st[i+1]; \
               double CT = bendParam[i][0]*y_ct[i-1]+bendParam[i][1]*y_ct[i]+bendParam[i][2]*y_ct[i+1]; \
               yct_bend -= Y*CT; yst_bend -= Y*ST; cct_bend += C*CT; sct_bend += S*CT; cst_bend += C*ST; sst_bend += S*ST;
#define CALCF  FVARS for (int i = 0; i < length(); i++) {FRMSD} for (int i = 1; i < length()-1; i++) {FBEND} \
               const double yy = yy_rmsd + yy_bend/lambda(), yc = yc_rmsd + yc_bend/lambda(), ys = ys_rmsd + ys_bend/lambda(); \
               const double cc = cc_rmsd + cc_bend/lambda(), sc = sc_rmsd + sc_bend/lambda(), ss = ss_rmsd + ss_bend/lambda();
#define CALCDF FVARS DFVARS for (int i = 0; i < length(); i++) {FRMSD DFRMSD} for (int i = 1; i < length()-1; i++) {FBEND DFBEND} \
               const double yy = yy_rmsd + yy_bend/lambda(), yc = yc_rmsd + yc_bend/lambda(), ys = ys_rmsd + ys_bend/lambda(); \
               const double cc = cc_rmsd + cc_bend/lambda(), sc = sc_rmsd + sc_bend/lambda(), ss = ss_rmsd + ss_bend/lambda(); \
               const double yct = yct_rmsd + yct_bend/lambda(), yst = yst_rmsd + yst_bend/lambda(); \
               const double cct = cct_rmsd + cct_bend/lambda(), sct = sct_rmsd + sct_bend/lambda(), cst = cst_rmsd + cst_bend/lambda(), sst = sst_rmsd + sst_bend/lambda();


const double TimeSeries::phi_start  = 2*M_PI;
const double TimeSeries::alpha_start=-.4;
const int    TimeSeries::giveUpTime = 10; // give up after 10 seconds
const double TimeSeries::stepSize   =  0.01;
const double TimeSeries::tolerance  =  0.01;
const int    TimeSeries::maxIter = 100;

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
	_lambda = .05;
	cout << "lambda\tphi\talpha\terr\tA\ttheta\tbendSpline\trmsSpline\tbendTrend\trmsFit" << endl;
	for (int i = 0; i < 25; i++)
	{
 		getCompleteMatrix();
 		//getSolutionVectors();
 		//solveSpline();
 		//writeFiles();
// 	solve(y, spline);    // do this here since lambda wil not be changed
// 		static gsl_vector *dummy=gsl_vector_alloc(2);
//  		for (_phi=5.5; _phi<7; _phi+=.05) fdf(phi(), _alpha, dummy);
		fit_constLambda();
		//writeFiles();
		_lambda *= 2;
	}
}

void TimeSeries::test() // at lambda about -1.98, -.51, -.17 and .42, matrix is singular!!!
{
	double a_data[] = { 0.18, 0.60, 0.57, 0.96,
	                    0.41, 0.24, 0.99, 0.58,
	                    0.14, 0.30, 0.97, 0.66,
	                    0.51, 0.13, 0.19, 0.85 };
	double b_data[] = { 1.0, 2.0, 3.0, 4.0 };
	gsl_matrix_view m = gsl_matrix_view_array (a_data, 4, 4);
	gsl_vector_view b = gsl_vector_view_array (b_data, 4);
	int s;
	gsl_matrix *A = gsl_matrix_alloc(4,4);
	gsl_vector *x = gsl_vector_alloc (4);
	gsl_permutation * p = gsl_permutation_alloc (4);
	gsl_matrix *inverse = gsl_matrix_alloc(4,4);
	gsl_matrix *inv2 = gsl_matrix_alloc(4,4);
	//for (int lambda = -190; lambda < -53; lambda++)
	//for (int lambda = -50; lambda < -18; lambda++)
	for (int lambdaI = -16; lambdaI < 41; lambdaI++)
	{
		const double lambda = .01*lambdaI;
		gsl_matrix_memcpy(A, &m.matrix);
		for (int i = 0; i < 4; i++) gsl_matrix_set(A, i, i, gsl_matrix_get(A, i, i)+lambda);
		gsl_linalg_LU_decomp (A, p, &s);
		gsl_linalg_LU_invert(A, p, inverse);
		gsl_matrix_memcpy(inv2, inverse);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, inverse, inverse, 0, inv2);
		if (true) // print inverse and d(inv)
		{
			cout << lambda;
			for (int i = 0; i < 4; i++)for (int k = 0; k < 4; k++) cout << "\t" << gsl_matrix_get(inverse, i, k);
			for (int i = 0; i < 4; i++)for (int k = 0; k < 4; k++) cout << "\t" << gsl_matrix_get(inv2, i, k);
			cout << endl;
		}
		if (false) // solve eqn system
		{
			double det = gsl_linalg_LU_det(A, s);
			gsl_blas_dgemv (CblasNoTrans, 1, inverse, &b.vector, 0, x);
			cout << lambda << "\t" << det;
			for (int i = 0; i < 4; i++) cout << "\t" << gsl_vector_get(x, i);
			gsl_linalg_LU_solve(A, p, &b.vector, x);
			for (int i = 0; i < 4; i++) cout << "\t" << gsl_vector_get(x, i); cout << endl;
		}
	}
	gsl_permutation_free (p);
	gsl_matrix_free(A);
	gsl_vector_free (x);
	gsl_matrix_free(inverse);
	gsl_matrix_free(inv2);
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
}

void TimeSeries::resizeVectors()
{
	spline.resize(length());
	trend. resize(length());
	fit.   resize(length());
	ecos.  resize(length());
	esin.  resize(length());
	ecost. resize(length());
	esint. resize(length());
	y_c.   resize(length());
	y_s.   resize(length());
	y_ct.  resize(length());
	y_st.  resize(length());
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
	for (int i = 0; i < length(); i++) *gsl_matrix_ptr(mat, i, i) += lambda();
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
		ecost[tNr]=ecos[tNr]*x[tNr];
		esint[tNr]=esin[tNr]*x[tNr];
	}
}

void TimeSeries::fit_constLambda() //resizeVectors(), getBendParameters(), getBendPenaltyMatrix() must have been done before!
{
	getCompleteMatrix(); // do this here since lambda wil not be changed
	solve(y, spline);    // do this here since lambda wil not be changed
	_iter = 0;
	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
	gsl_multimin_fdfminimizer *m = gsl_multimin_fdfminimizer_alloc (T, 2);
	gsl_vector *x= gsl_vector_alloc (2);
	gsl_vector_set (x, 0, phi_start);
	gsl_vector_set (x, 1, alpha_start); 
	gsl_multimin_function_fdf func;
	func.n = 2;  /* number of function components */
	func.f = &f_constLambda;
	func.df = &df_constLambda;
	func.fdf = &fdf_constLambda;
	func.params = (void *)this;
	gsl_multimin_fdfminimizer_set (m, &func, x, stepSize, tolerance);
	time_t t0 = time(0);
	time_t t1 = t0;
	do
	{
		_iter++;
		_status = gsl_multimin_fdfminimizer_iterate (m);
		if (status()) break;
		_status = gsl_multimin_test_gradient (m->gradient, 1);
		t1 = time(0);
	}
	while (status() == GSL_CONTINUE && (t1-t0) < giveUpTime && iter() < maxIter);
	//writeFiles();
	
	CALCF
	
	
	_theta = atan2(sc*yc-cc*ys, sc*ys-ss*yc);
	const double C = cos(theta());
	const double S = sin(theta());
	_A = -(C*yc+S*ys) / (C*C*cc+2*S*C*sc+S*S*ss);
	
	const double bendSpline = yy_bend;
	const double rmsSpline  = yy_rmsd;
	const double bendTrend  = yy_bend+A()*A()*C*C*cc_bend+A()*A()*S*S*ss_bend+2*A()*C*yc_bend+2*A()*S*ys_bend+2*A()*A()*S*C*sc_bend;
	const double rmsFit     = yy_rmsd+A()*A()*C*C*cc_rmsd+A()*A()*S*S*ss_rmsd+2*A()*C*yc_rmsd+2*A()*S*ys_rmsd+2*A()*A()*S*C*sc_rmsd;
	for (int i = 0; i < length(); i++) 
	{
		trend[i] = spline[i]-A()*(C*y_c[i] +S*y_s[i]);
		fit[i]   = trend[i] +A()*(C*ecos[i]+S*esin[i]);
	}
	cout << lambda() << "\t" << phi() << "\t" << alpha() << "\t" << m->f << "\t" << A() << "\t" << theta() << "\t" << bendSpline << "\t" << rmsSpline << "\t" << bendTrend << "\t" << rmsFit << endl;
	const double errControl1 = (rmsFit+bendTrend/lambda()) / (rmsSpline+bendSpline/lambda());
	const double errControl2 = (rmsd(fit)+bendPenalty(trend)/lambda()) / (rmsd(spline)+bendPenalty(spline)/lambda());
	cerr << bendSpline-bendPenalty(spline) << "\t" << rmsSpline-rmsd(spline) << "\t" << bendTrend-bendPenalty(trend) << "\t" << rmsFit-rmsd(fit) << "\t" << m->f-errControl1 << "\t" << m->f-errControl2 << endl;
	//
	//
	//
	// control output
// 	double err = p.yy + A()*A()*C*C*p.cc + A()*A()*S*S*p.ss + 2*A()*C*p.yc + 2*A()*S*p.ys + 2*A()*A()*S*C*p.sc;
// 	double nom = C*p.yc+S*p.ys;
// 	nom *= nom;
// 	double den = C*C*p.cc + S*S*p.ss + 2*S*C*p.sc;
// 	double err2 = p.yy - nom/den;
// 	double err3 = p.yy - (p.cc*p.ys*p.ys - 2*p.sc*p.yc*p.ys + p.ss*p.yc*p.yc) / (p.cc*p.ss-p.sc*p.sc);
// 	cout << p.yy << "\t" << err << "\t" << err2 << "\t" << err3 << endl;
	gsl_multimin_fdfminimizer_free(m);
	gsl_vector_free(x);
	writeFiles();
}


double TimeSeries::fdf(const double Phi, const double Alpha, gsl_vector* df)
{
	_phi=Phi;
	_alpha=Alpha;
	//getCompleteMatrix(); NOT NECESSARY SINCE LAMBDA CONSTANT
	getSolutionVectors();
	//solve(y, spline);    NOT NECESSARY SINCE LAMBDA CONSTANT
	solve(ecos, y_c);
	solve(esin, y_s);
	solve(ecost, y_ct);
	solve(esint, y_st);
	
	CALCDF
	
	//const double err = yy - (cc*ys*ys          - 2*sc*yc*ys                     + ss*yc*yc)          / (cc*ss       -sc*sc);
	double nom   = cc*ys*ys - 2*sc*yc*ys + ss*yc*yc;
	double den   = cc*ss       -sc*sc;
	double halfnom_p = (cc*ys*yct-cst*ys*ys) - (cct*yc*ys-sst*yc*ys-sc*yst*ys+sc*yc*yct) + (sct*yc*yc-ss*yc*yst);
	double halfden_p = (cc*sct-ss*cst)-(sc*cct-sc*sst);
	double halfnom_a = (cct*ys*ys+cc*yst*ys) - (sct*yc*ys + cst*yc*ys + sc*yct*ys + sc*yc*yst) + (sst*yc*yc+ss*yc*yct);
	double halfden_a = (cct*ss+cc*sst) - (sc*sct+sc*cst);
	double err = yy-nom/den;
	double err_a = -2*(halfnom_a*den-nom*halfden_a)/(den*den);
	double err_p = -2*(halfnom_p*den-nom*halfden_p)/(den*den);
	err /= yy;
	err_a /= yy;
	err_p /= yy;
	gsl_vector_set(df,0,err_p);
	gsl_vector_set(df,1,err_a);
	return err;
}

void TimeSeries::writeFiles()
{
	static bool firstCall = true;
	std::ofstream file;

	if (firstCall) 
	{
		remove("data.csv");
		file.open("data.csv"); 
		for (int i = 0; i < length(); i++) file << x[i] << "\t" << y[i] << endl; 
		file.close();
		remove("spline.csv");
		remove("trend.csv");
		remove("param.csv");
		remove("fit.csv");
		firstCall = false;
	}
	
	file.open("spline.csv", std::ios_base::app);
	for (int i = 0; i < length(); i++) file << spline[i] << "\t"; file << endl;
	file.close();
	
	file.open("trend.csv", std::ios_base::app);
	for (int i = 0; i < length(); i++) file << trend[i] << "\t"; file << endl;
	file.close();
	
	file.open("fit.csv", std::ios_base::app);
	for (int i = 0; i < length(); i++) file << fit[i] << "\t"; file << endl;
	file.close();
	
	
	file.open("param.csv", std::ios_base::app);
	file << lambda() << "\t" << A() << "\t" << alpha() << "\t" << phi() << "\t" << theta() << "\t" << bendPenalty(spline) << "\t" << rmsd(spline) << "\t" << bendPenalty(trend) << "\t" << rmsd(fit) << endl;
	file.close();
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

double TimeSeries::rmsd(const std::vector< double >& curve)
{
	double err = 0;
	for (int i = 0; i < length(); i++) 
	{
		const double errPart = curve[i]-y[i];
		err += errPart*errPart;
	}
	return err;
}
