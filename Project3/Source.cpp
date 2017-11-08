//#include <algorithm>
//#include <string>
//#include <iostream>
//#include <tuple>
//#include <vector>
//#include <random>
//#include <ilcplex/ilocplex.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <time.h>
//#include <math.h>
//
//ILOSTLBEGIN
//
//using namespace std; //std::function name
//static string interpolate_string(string var, int i);
//static vector<double> computer_s_i(vector<double> sigma, double alpha);
//static vector<double> generate_sigma(int n);
//static double solve_integer_program(vector<double> s_i);
//
//int main(int argc, char **argv)
//{
//	double epsilon = 0.01;
//	double obj = -1e10;
//	double alpha = 1;
//
//	vector<double> sigma = generate_sigma(10);
//	vector<double> s_i = computer_s_i(sigma, alpha);
//
//	while (obj < log1p(epsilon))
//	{
//		obj = solve_integer_program(s_i);
//		alpha = alpha / 2;
//		s_i = computer_s_i(sigma, alpha);
//	}
//
//	system("PAUSE");
//	return 0;
//}
//
//static string interpolate_string(string var, int i)
//{
//	return var + "_{" + to_string(i) + "}";
//}
//
//static vector<double> computer_s_i(vector<double> sigma, double alpha) {
//	vector<double> result;
//
//	for (int i = 0; i < sigma.size(); i++)
//	{
//		double s_i = -pow(sigma[i], 2) / (2 * pow(alpha, 2));
//		result.push_back(s_i);
//	}
//
//	return result;
//}
//
//static vector<double> generate_sigma(int n)
//{
//	vector<double> result;
//
//	random_device rand_dev;
//	mt19937 mt(rand_dev());
//	uniform_real_distribution<double> dist(1.0, 10.0);
//
//	for (int i = 0; i < n; i++)
//	{
//		result.push_back(dist(mt));
//	}
//
//	return result;
//}
//
//static double solve_integer_program(vector<double> s_i)
//{
//	double result = 0.0;
//
//	IloEnv env;
//	try {
//		IloModel model(env);
//		IloBoolVarArray x_i(env, s_i.size());
//
//		IloObjective obj = IloAdd(model, IloMinimize(env));
//
//		for (int i = 0; i < s_i.size(); i++)
//		{
//			x_i[i] = IloBoolVar(env, interpolate_string("x", i).c_str());
//		}
//
//		IloNumExpr sum_over_i(env);
//
//		for (int i = 0; i < s_i.size(); i++)
//		{
//			sum_over_i += s_i[i] * x_i[i];
//		}
//
//		obj.setExpr(sum_over_i);
//		sum_over_i.end();
//
//		IloCplex cplex(model);
//
//		cplex.exportModel("test.lp");
//
//		if (!cplex.solve()) {
//			env.error() << "Failed to optimize LP." << endl;
//			throw(-1);
//		}
//
//		cout << "Solution value = " << cplex.getObjValue() << endl;
//	}
//
//	catch (IloException& e) {
//		cerr << "Concert exception caught: " << e << endl;
//	}
//
//	catch (...) {
//		cerr << "Unknown exception caught" << endl;
//	}
//
//	env.end();
//
//	return result;
//}