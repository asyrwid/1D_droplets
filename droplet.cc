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
double U_ab = -0.8;
int Natoms = 4;
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

tie(H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab, H_edge) = get_H(sites, t, U, U_ab);
printfln("Maximum bond dimension of H is %d",maxLinkDim(H_total));


auto psi0 = initial_state(sites, Natoms);
printfln("Energy %d", innerC(psi0, H_total, psi0));


std::string dir = "/home/asyrwid/ITensor-3.1.6/programs/1D_droplets/";

std::string densities_entropies = dir + "densities_entropies.txt";
std::string convergence_params = dir + "convergence_params.txt";
std::string sites_file = dir + "sites.txt";
std::string mps_file = dir + "mps.txt";
vector<string> column_names_dens_entrs = {"site", "entropy_a", "entropy_b", "density_a", "density_b"};
vector<string> column_names_conv_params = {"centr_entropy_a", "centr_entropy_b",
                                           "E_tot", "E_hop_a", "E_hop_b", "E_aa", "E_bb", "E_ab"};
prepare_file(column_names_dens_entrs, densities_entropies);
prepare_file(column_names_conv_params, convergence_params);


int NoOfSteps = 5;
int nosweeps = 4;
Real dt_bysweep = 0.02;
int MaxBondDim_tdvp = 200;

MPS psi1 = imag_time_evol(sites,
                          psi0,
                          H_total,
                          NoOfSteps,
                          nosweeps,
                          dt_bysweep,
                          MaxBondDim_tdvp);



vector<MPO> H_terms {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab};
int MaxBondDim_dmrg = 256;

dmrg_sequence(sites,
              psi1,
              H_terms,
              MaxBondDim_dmrg,
              densities_entropies,
              convergence_params,
              sites_file,
              mps_file);

/*
auto psi1 = imag_time_evol(sites, psi0, H_total, densities_entropies,
  NoOfSteps, nosweeps, dt_bysweep, MaxBondDim);
*/
/*
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
*/









return 0;
}
