#include <algorithm>
#include <string>
#include <iostream>
#include <tuple>
#include <vector>
#include <random>
#include <list>
#include <numeric>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <conio.h>
#include <dlib/optimization.h>
#include <Eigen/Sparse>
#include <cstdio>
#include <cmath>
#include <iomanip>

//#include <ilcplex/ilocplex.h>


//ILOSTLBEGIN
using namespace std; 
using dlib::find_min_single_variable;


double randnum(double a, double b)
{
	static default_random_engine generator;
	uniform_real_distribution<double> distribution(a, b);
	return distribution(generator);
}

bool wayToSort(double i, double j) { return i > j; }

double minDistance(vector<double> dist, vector<bool> sptSet, int T)
{
	int min = INT_MAX;
	int min_index;

	for (int v = 0; v < T; v++)
	{
		if (sptSet[v] == false && dist[v] <= min)
		{
			min = dist[v];
			min_index = v;
		}
	}
	return min_index;
}

double dijkstra(vector<vector<double>>G, int S, int T)
{
	vector<double> dist;
	vector<bool> sptSet; 

	for (int i = 0; i < T; i++)
	{
		dist.push_back(INT_MAX);
		sptSet.push_back(false);
	}

	dist[S - 1] = 0;

	for (int count = 0; count < T; count++)
	{
		int u = minDistance(dist, sptSet, T);

		sptSet[u] = true;

		for (int v = 0; v < T; v++)
		{
			if (!sptSet[v] && G[u][v] && dist[u] != INT_MAX && dist[u] + G[u][v] < dist[v])
			{
				dist[v] = dist[u] + G[u][v];
			}
		}
	}
	
	return dist[T - 1];
}

//double dijkstra(vector< vector<double> > G, int S, int T) {
//	int n = G.size();
//	vector<double> dist(n);
//	vector<bool> vis(n);
//
//	for (int i = 0; i < n; ++i) {
//		dist[i] = INT_MAX;
//	}
//	dist[S] = 0;
//
//	for (int i = 0; i < n; ++i) {
//		int cur = -1;
//		for (int j = 0; j < n; ++j) {
//			if (vis[j]) continue;
//			if (cur == -1 || dist[j] < dist[cur]) {
//				cur = j;
//			}
//		}
//
//		vis[cur] = true;
//		for (int j = 0; j < n; ++j) {
//			int path = dist[cur] + G[cur][j];
//			if (path < dist[j]) {
//				dist[j] = path;
//			}
//		}
//	}
//
//	return dist[T - 1];
//}

vector<double> CEcalculation(double alpha, vector<vector<double>> INFO)
{
	vector<double> ArcCE;

	for (int i = 0; i < INFO[0].size(); i++)
	{
		ArcCE.push_back(INFO[0][i] + pow(INFO[1][i], 2) / 2 / alpha);
	}

	return ArcCE;
}

vector<double> CE_subgradient_calculation(double alpha, vector<vector<double>> INFO)
{
	vector<double> ArcCE_subgradient;

	for (int i = 0; i < INFO[1].size(); i++)
	{
		ArcCE_subgradient.push_back(-(INFO[1][i]*INFO[1][i])/2/(alpha*alpha));
	}

	return ArcCE_subgradient;
}

bool check_lower_bound_or_not(double alpha, vector<vector<double>> INFO,
	vector<vector<double>> ArcList, int N, int S, int T, double beta)
{
	bool lower_bound_or_not = false;
	vector<double> ArcCE_subgradient = CE_subgradient_calculation(alpha, INFO);

	vector<vector<double>> G(N, vector<double>(N, 0.0));

	//for (int i = 0; i < ArcList.size(); i++)
	//{
	//	cout << "(" + to_string(ArcList[i][0]) + "," + to_string(ArcList[i][1]) + ")" << endl;
	//}

	for (int i = 0; i < ArcList.size(); i++)
	{
		G[ArcList[i][0]][ArcList[i][1]] = -(ArcCE_subgradient[i]);
	}

	double DIST = dijkstra(G, S, T);

	if ((-DIST) <= log(beta))
	{
		lower_bound_or_not = true;
	}

	else
	{
		lower_bound_or_not = false;
	}
	
	return lower_bound_or_not;
}

bool check_upper_bound_or_not(double alpha, vector<vector<double>> INFO,
	vector<double> ArcDistMean, double epsilon)
{
	bool upper_bound_or_not = true;

	vector<double> ArcCE = CEcalculation(alpha, INFO);

	for (int i = 0; i < ArcCE.size(); i++)
	{
		if (ArcCE[i] > (1 + epsilon) * ArcDistMean[i])
		{
			upper_bound_or_not = false;
			break;
		}
	}

	return upper_bound_or_not;
}

double find_alpha_lower_bound(vector<vector<double>> INFO,
	vector<vector<double>> ArcList, int N, int S, int T, double beta, double tolx)
{
	double alpha = 1.0;
	double alpha_max = 0.0;
	double alpha_min = 0.0;

	if (check_lower_bound_or_not(alpha, INFO, ArcList, N, S, T, beta))
	{
		alpha = alpha * 2;

		while (check_lower_bound_or_not(alpha, INFO, ArcList, N, S, T, beta))
		{
			alpha = alpha * 2;
		}

		alpha_max = alpha;
		alpha_min = (alpha / 2);
	}

	else
	{
		alpha = (alpha / 2);

		while (check_lower_bound_or_not(alpha, INFO, ArcList, N, S, T, beta) == 0)
		{
			alpha = alpha / 2;
		}

		alpha_max = alpha * 2;
		alpha_min = alpha;
	}

	while  ((alpha_max - alpha_min) > tolx)
	{
		alpha = (alpha_max + alpha_min) / 2;

		if (check_lower_bound_or_not(alpha, INFO, ArcList, N, S, T, beta))
		{
			alpha_min = alpha;
		}

		else
		{
			alpha_max = alpha;
		}
	}

	return alpha_min;
}

double find_alpha_upper_bound(vector<vector<double>> INFO, double epsilon)
{
	vector<double> ArcDistMean;
	double alpha = 1.0;
	double alpha_max = 0.0;

	for (int i = 0; i < INFO[0].size(); i++)
	{
		ArcDistMean.push_back(INFO[0][i]);
	}

	if (check_upper_bound_or_not(alpha, INFO, ArcDistMean, epsilon) == 1)
	{
		alpha = alpha / 2;

		while (check_upper_bound_or_not(alpha, INFO, ArcDistMean, epsilon) == 1)
		{
			alpha = alpha / 2;
		}

		alpha_max = alpha * 2;
	}

	else
	{
		alpha = alpha * 2;

		while (check_upper_bound_or_not(alpha, INFO, ArcDistMean, epsilon) == 0)
		{
			alpha = alpha * 2;
		}

		alpha_max = alpha;
	}

	return alpha_max;
}

class mathProblem
{
public:
	mathProblem(const double& C_alpha_min_iIn, const double& alpha_minIn, const double& alphaIn, const double& C_alpha_iIn, const double& epsilonIn, vector<vector<double>>& INFOIn, int& iIn)
	{
		C_alpha_min_i = C_alpha_min_iIn;
		alpha_min = alpha_minIn;
		alpha = alphaIn;
		C_alpha_i = C_alpha_iIn;
		epsilon = epsilonIn;
		INFO = INFOIn;
		i = iIn;

	}

	double operator() (const double& x) const
	{
		return -(C_alpha_min_i + (x - alpha_min) / (alpha - alpha_min) * (C_alpha_i - C_alpha_min_i) - (1 + epsilon) * CEcalculation(x, INFO)[i]);
	}

private:
	double C_alpha_min_i;
	double alpha_min;
	double alpha;
	double C_alpha_i;
	double epsilon;
	vector<vector<double>> INFO;
	int i;
};

double findMin(double C_alpha_min_iIn, double alpha_minIn, double alphaIn, double C_alpha_iIn, double epsilonIn, vector<vector<double>> INFOIn, int iIn, double x0)
{
	return find_min_single_variable(mathProblem(C_alpha_min_iIn, alpha_minIn, alphaIn, C_alpha_iIn, epsilonIn, INFOIn, iIn), x0, alpha_minIn, alphaIn, epsilonIn);
}

double binarysearch_objective(double alpha_min, double alpha, double epsilon, vector<vector<double>> INFO)
{

	vector<double> C_alpha_min;
	for (int i = 0; i < INFO[0].size(); i++)
	{
		C_alpha_min.push_back(INFO[0][i] + pow(INFO[1][i], 2) / 2 / alpha_min);
	}


	vector<double> C_alpha;
	for (int i = 0; i < INFO[0].size(); i++)
	{
		C_alpha.push_back(INFO[0][i] + pow(INFO[1][i], 2) / 2 / alpha);
	}



	int n = INFO[0].size();
	double fmax = -INFINITY;
	for (int i = 0; i < n; i++)
	{
		double C_alpha_min_i = C_alpha_min[i];
		double C_alpha_i = C_alpha[i];
		double fi = findMin(C_alpha_min_i, alpha_min, alpha, C_alpha_i, epsilon, INFO, i, (alpha_min + alpha) / 2);
		fmax = max(fmax, -fi);
	}
	
	return fmax;
}



double BinarysearchBP(double alpha_current, double alpha_max, double epsilon, vector<vector<double>> INFO)
{
	double fmax = binarysearch_objective(alpha_current, alpha_max, epsilon, INFO);
	double alpha_next = 0.0;

	if (fmax <= 0)
	{
		alpha_next = alpha_max;
	}

	else
	{
		double a = alpha_current;
		double b = alpha_max;
		double p = (a + b) / 2;
		double err = abs(binarysearch_objective(alpha_current, p, epsilon, INFO));

		while (err > 1e-7)
		{
 			if (binarysearch_objective(alpha_current, a, epsilon, INFO) * binarysearch_objective(alpha_current, p, epsilon, INFO) < 0)
			{
				b = p;
			}

			else
			{
				a = p;
			}

			p = (a + b) / 2;
			err = abs(binarysearch_objective(alpha_current, p, epsilon, INFO));
		}

		alpha_next = p;
	}

	return alpha_next;
}

int main(int argc, char **argv)
{
	const int N = 300;
	double XYpos[N][2];


	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 2; j++) {
			XYpos[i][j] = 0;
		}
	}

	XYpos[N - 1][0] = 1;
	XYpos[N - 1][1] = 1;


	//Genereate random points
	for (int i = 1; i < N - 1; i++)
	{
		int found = 0;

		while (found == 0)
		{
			double x = randnum(0.0, 1.0);
			double y = randnum(0.0, 1.0);

			if (~((x < 0.4 && x > 0.2 && y < 0.8) || (x < 0.8 && x > 0.6 && y > 0.2)))
			{
				found = 1;
				XYpos[i][0] = x;
				XYpos[i][1] = y;
			}
		}
	}


	double DivGap = 3.0;
	double DivMin = 0.0;

	double DistMat[N][N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				double x_diff = XYpos[i][0] - XYpos[j][0];
				double y_diff = XYpos[i][1] - XYpos[j][1];
				DistMat[i][j] = sqrt(x_diff*x_diff  + y_diff*y_diff);
			}

			else
			{
				DistMat[i][j] = 0;
			}
		}
	}

	int Arcs[N][N];

	//Determine distance between nodes
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double temp = randnum(0.0, 1.0);

			if (temp <= 0.4)
			{
				Arcs[i][j] = 1;
			}

			else
			{
				Arcs[i][j] = 0;
			}
		}
	}


	//arcs matrix
	for (int i = 0; i < N; i++)
	{
		double temp = randnum(0.0, 1.0);

		if (temp <= 0.2)
		{
			Arcs[0][i] = 1;
		}

		else
		{
			Arcs[0][i] = 0;
		}
	}

	for (int i = 0; i < N; i++)
	{
		double temp = randnum(0.0, 1.0);

		if (temp <= 0.2)
		{
			Arcs[i][N - 1] = 1;
		}
		else
		{
			Arcs[i][N - 1] = 0;
		}
	}

	double MaxArcLen = 0.13;
	double MinArcLen = 0.03;
	double MaxArcLen2 = 0.15;


	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i = j) {
				Arcs[i][j] = 0; //no self arcs
			}
		}
	}


	for (int i = 0; i < N; i++)
	{
		Arcs[i][0] = 0;
	}

	//for (int i = 0; i < N; i++)
	//{
	//	Arcs[N-1][i] = 0;
	//}

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if ((Arcs[i][j] == 1) && (i != 0) && (j != N - 1))
			{
				if ((DistMat[i][j] > MaxArcLen) || (DistMat[i][j] < MinArcLen))
				{
					Arcs[i][j] = 0;
				}
			}

			if ((Arcs[i][j] == 1) && (i == 0 || j == N - 1))
			{
				if (DistMat[i][j] > MaxArcLen2)
				{
					Arcs[i][j] = 0;
				}
			}
		}
	}

	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i+1; j < N; j++)
		{
			if ((Arcs[i][j] == 1) && (Arcs[j][i] == 1))
			{
				if (randnum(0.0, 1.0) > 0.5)
				{
					Arcs[i][j] = 0;
				}
				else {
					Arcs[j][i] = 0;
				}
			}
		}
	}

	int M = 0;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			M += Arcs[i][j];
		}
	}



	vector<vector<double>> ArcList;

	for (int i = 0; i < M; i++)
	{
		ArcList.push_back(vector<double>(2, 0.0));
	}

	vector<double> ArcDist;
	for (int i = 0; i < M; i++)
	{
		ArcDist.push_back(0.0);
	}


	vector<double> ArcDiv;
	for (int i = 0; i < M; i++)
	{
		ArcDiv.push_back(0.0);
	}

	vector<vector<double>> DivDistMat;
	for (int i = 0; i < N; i++)
	{
		DivDistMat.push_back(vector<double>(N, 0.0));
	}


	int idx = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (Arcs[i][j] == 1)
			{
				idx = idx + 1;
				ArcList[idx - 1][0] = i;
				ArcList[idx - 1][1] = j;
				ArcDist[idx - 1] = DistMat[i][j];
				DivDistMat[i][j] = ArcDiv[idx - 1];
			}
		}
	}

	//[Y,DistSrt] = sort(-ArcDist)
	vector<double> minus_ArcDist(ArcDist);
	for_each(minus_ArcDist.begin(), minus_ArcDist.end(), [](double & x)
	{
		x = -x;
	});


	vector <double *> Y;
	Y.reserve(minus_ArcDist.size());
	for (auto& value : minus_ArcDist)
	{
		Y.push_back(&value);
	}

	auto ascendingOrderSorter1 = [](double * i, double * j)
	{
		return *i < *j;
	};

	sort(Y.begin(), Y.end(), ascendingOrderSorter1);
	vector<int> DistSrt(minus_ArcDist.size());
	iota(DistSrt.begin(), DistSrt.end(), 0);
	auto comparator = [&minus_ArcDist](int a, int b) { return minus_ArcDist[a] < minus_ArcDist[b]; };
	sort(DistSrt.begin(), DistSrt.end(), comparator);


	for (int idx = 0; idx < M; idx++)
	{
		ArcDiv[idx] = randnum(0.0, 1.0) * DivGap + DivMin;
	}


	//[Y,I] = sort(-ArcDiv);
	vector<double> minus_ArcDiv(ArcDiv);
	for_each(minus_ArcDiv.begin(), minus_ArcDiv.end(), [](double & x)
	{
		x = -x;
	});
	vector<int> I(minus_ArcDiv.size());
	iota(I.begin(), I.end(), 0);
	auto comparator_2 = [&minus_ArcDiv](int a, int b) { return minus_ArcDiv[a] < minus_ArcDiv[b]; };
	sort(I.begin(), I.end(), comparator_2);


	for (auto& value : minus_ArcDiv)
	{
		Y.push_back(&value);
	}
	sort(Y.begin(), Y.end(), ascendingOrderSorter1);

	

	sort(ArcDiv.begin(), ArcDiv.end(), wayToSort);



	for (size_t i = 0; i < I.size(); ++i)
	{
		if (i < I[i])
		{
			swap(ArcDist[i], ArcDist[I[i]]);
		}
	}



	for (size_t i = 0; i < I.size(); ++i)
	{
		if (i < I[i])
		{
			swap(ArcList[i], ArcList[I[i]]);
		}
	}
	
	int S = 1;
	int T = N;

	//step 1

	vector<double> mean;
	for (int i = 0; i < M; i++)
	{
		mean.push_back((ArcDist[i] + (ArcDiv[i] / 2)));
	}
	
	vector<double> std;
	for (int i = 0; i < M; i++)
	{
		std.push_back((ArcDiv[i] / 4));
	}

	vector<vector<double>> INFO;
	INFO.push_back(mean);
	INFO.push_back(std);

	double epsilon = 1e-3;
	double beta = 0.01;
	double tolx = 1e-6;

	//find _alpha_upperbound
	vector <double> ArcDistMean;
	for (int i = 0; i < M; i++)
	{
		ArcDistMean.push_back(mean[i]);
	}

	double alpha = 1.0;

	double alpha_max = find_alpha_upper_bound(INFO, epsilon);
	double alpha_min = find_alpha_lower_bound(INFO, ArcList, N, S, T, beta, tolx);

	//step 3
	vector<double> alphabp;
	alphabp.push_back(alpha_min);


	double alpha_next = BinarysearchBP(alphabp.back(), alpha_max, epsilon, INFO);

	//while (alpha_max - alphabp.back() >= tolx)
	//{
	//	alphabp.push_back(alpha_next);
	//}


	system("PAUSE");
	return 0;
} 

