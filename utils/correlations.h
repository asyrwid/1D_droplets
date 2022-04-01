#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace itensor;
using namespace std;


double density_density_correlation_a(Boson sites, MPS state, int i_site, int j_site);
double density_density_correlation_b(Boson sites, MPS state, int i_site, int j_site);

double one_body_correlation_a(Boson sites, MPS state, int i_site, int j_site);
double one_body_correlation_b(Boson sites, MPS state, int i_site, int j_site);

double pair_correlation_ab(Boson sites, MPS state, int i_site, int j_site);

tuple<vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>> correlations(Boson sites, MPS state);

void save_correlations(Boson sites, MPS state, string path);
