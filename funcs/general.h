
#include "itensor/all.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>


using namespace itensor;
using namespace std;

tuple<MPO, MPO, MPO, MPO, MPO, MPO, MPO> get_H(SiteSet& sites,  double t, double U, double U_ab);

MPS initial_state(Boson sites, int Natoms);
