#include "itensor/all.h"
//#include "tdvp.h"
//#include "basisextension.h"
#include "funcs/general.h"
#include <math.h>
#include <iostream>
#include <string>
#include <stdlib.h> /* srand, rand*/
#include <time.h> /* time */

using std::vector;
using namespace itensor;
using namespace std;

int main(int argc, char** argv){

double t = 1;
double U = 1;
double U_ab = 1;
int Natoms = 10;
int L = 100;
int maxOccupation = Natoms;


int M_sites = 2*L;
cout<<"M_sites = "<<M_sites<<endl;
auto sites = Boson(M_sites, {"ConserveNb", true,
                             "MaxOcc=", int(maxOccupation),
                             "ConserveSz", true});
auto state = InitState(sites);


auto H_total = toMPO(AutoMPO(sites));
auto H_hop_a = toMPO(AutoMPO(sites));
auto H_hop_b = toMPO(AutoMPO(sites));
auto H_aa = toMPO(AutoMPO(sites));
auto H_bb = toMPO(AutoMPO(sites));
auto H_ab = toMPO(AutoMPO(sites));
auto H_edge = toMPO(AutoMPO(sites));

tie(H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab, H_edge) = get_H(sites,t,U,U_ab);
printfln("Maximum bond dimension of H is %d",maxLinkDim(H_total));


auto psi0 = initial_state(sites, Natoms);
printfln("Energy %d", innerC(psi0, H_total, psi0));


std::string densities_entropies = "/home/asyrwid/ITensor-3.1.6/programs/1D_droplets/kurwa.txt";

int NoOfSteps = 5;
int nosweeps = 4;
Real dt_bysweep = 0.02;
int MaxBondDim = 200;

MPS psi1 = imag_time_evol(sites, psi0, H_total, densities_entropies,
                   NoOfSteps,
                   nosweeps,
                   dt_bysweep,
                   MaxBondDim);
psi0 = psi1;

/*
auto psi1 = imag_time_evol(sites, psi0, H_total, densities_entropies,
  NoOfSteps, nosweeps, dt_bysweep, MaxBondDim);
*/

printfln("density_a %d", density_a(sites, psi0, 50));
printfln("density_b %d", density_b(sites, psi0, 10));

printfln("density_density_a %d", density_density_correlation_a(sites, psi0, 49, 50));
printfln("density_density_b %d", density_density_correlation_b(sites, psi0, 51, 55));

printfln("one_body_corr_a %d", one_body_correlation_a(sites, psi0, 49, 50));
printfln("one_body_corr_b %d", one_body_correlation_b(sites, psi0, 51, 55));

printfln("pair_corr_ab %d", pair_correlation_ab(sites, psi0, 51, 55));

printfln("entanglement_entropy %d", entanglement_entropy(psi0, 101));
printfln("entanglement_entropy %d", entanglement_entropy(psi0, 105));

vector<double> densities_a;
vector<double> densities_b;
vector<vector<double>> one_body_correlations_a;
vector<vector<double>> one_body_correlations_b;
vector<vector<double>> pair_correlations_ab;
vector<vector<double>> density_density_a;
vector<vector<double>> density_density_b;


tie(densities_a, densities_b) = particle_densities(sites, psi0);
tie(one_body_correlations_a,
    one_body_correlations_b,
    pair_correlations_ab,
    density_density_a,
    density_density_b) = correlations(sites, psi0);


for(int i = 0; i < L; i++){
  cout << "\nsite = " << i+1 << ";  dens_a = " << densities_a[i] << ";  dens_b = " << densities_b[i] << std::flush;
}










return 0;
}
