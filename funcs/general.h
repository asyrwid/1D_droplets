
#include "itensor/all.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>


using namespace itensor;
using namespace std;

tuple<MPO, MPO, MPO, MPO, MPO, MPO, MPO> get_H(SiteSet& sites,  double t, double U, double U_ab);

MPS initial_state(Boson sites, int Natoms);

MPS annihilate_boson(MPS ket, int site);
MPS create_boson(Boson sites, MPS ket, int site);
MPS apply_operator(Boson sites, string OpName, MPS ket, int site);

double density_a(Boson sites, MPS state, int site);
double density_b(Boson sites, MPS state, int site);

double density_density_correlation_a(Boson sites, MPS state, int i_site, int j_site);
double density_density_correlation_b(Boson sites, MPS state, int i_site, int j_site);

double one_body_correlation_a(Boson sites, MPS state, int i_site, int j_site);
double one_body_correlation_b(Boson sites, MPS state, int i_site, int j_site);

double pair_correlation_ab(Boson sites, MPS state, int i_site, int j_site);
