#include "itensor/all.h"
#include "funcs/general.h"
#include <math.h>
#include <iostream>
#include <string>

using std::vector;
using namespace itensor;
using namespace std;

int main(int argc, char** argv){

int Natoms = 2;
int L = 50;
int maxOccupation = Natoms;
double t = 1;
double U = 5;
double r = 0.1;
double U_ab = -U*(1-r);
int MaxBondDim = 256;

string dir = "/home/asyrwid/ITensor-3.1.6/programs/1D_droplets/data/";
string parameters = parameters_to_filename(Natoms, L, maxOccupation, t, U, r, MaxBondDim);

int M_sites = 2*L; // at each physical bond we have 2 sites where even and odd ones
                   // correspond to different components
cout<<"M_sites = "<<M_sites<<endl;
auto sites = Boson(M_sites, {"ConserveNb", true, "MaxOcc=", int(maxOccupation), "ConserveSz", true});
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

string densities_entropies = dir + "densities_entropies_" + parameters;
string convergence_params = dir + "convergence_params_" + parameters;
string sites_file = dir + "sites_" + parameters;
string mps_file = dir + "mps_" + parameters;
string corrs = dir + "correlations_" + parameters;
vector<string> column_names_dens_entrs = {"site", "entropy_a", "entropy_b", "density_a", "density_b"};
vector<string> column_names_conv_params = {"centr_entropy_a", "centr_entropy_b",
                                           "E_tot", "E_hop_a", "E_hop_b", "E_aa", "E_bb", "E_ab"};
prepare_file(column_names_dens_entrs, densities_entropies);
prepare_file(column_names_conv_params, convergence_params);

// Warm up with a few steps of imaginary time evolution;
MPS psi1 = imag_time_evol(sites, psi0, H_total);
// then perform DMRG to obtain the ground state and
vector<MPO> H_terms {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab};
auto psi_GS =dmrg_sequence(sites, psi1, H_terms, MaxBondDim, // DMRG
                           densities_entropies, convergence_params, sites_file, mps_file // file paths
                          );

save_correlations(sites, psi_GS, corrs);

return 0;
}
