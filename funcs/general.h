#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace itensor;
using namespace std;

tuple<MPO, MPO, MPO, MPO, MPO, MPO, MPO> get_H(SiteSet& sites,  double t, double U, double U_ab);

MPS initial_state(Boson sites, int Natoms);

string parameters_to_filename(int Natoms,
                              int L,
                              int maxOccupation,
                              double t,
                              double U,
                              double r,
                              int MaxBondDim);

MPS apply_operator(Boson sites, string OpName, MPS ket, int site);

double density_a(Boson sites, MPS state, int site);
double density_b(Boson sites, MPS state, int site);

double density_density_correlation_a(Boson sites, MPS state, int i_site, int j_site);
double density_density_correlation_b(Boson sites, MPS state, int i_site, int j_site);

double one_body_correlation_a(Boson sites, MPS state, int i_site, int j_site);
double one_body_correlation_b(Boson sites, MPS state, int i_site, int j_site);

double pair_correlation_ab(Boson sites, MPS state, int i_site, int j_site);

double entanglement_entropy(MPS state, int lattice_site);

tuple<vector<double>, vector<double>> entropies(MPS state);

tuple<vector<double>, vector<double>> particle_densities(Boson sites, MPS state);

tuple<vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>> correlations(Boson sites, MPS state);

void save_correlations(Boson sites, MPS state, string path);

void prepare_file(vector<string> column_names, std::string path);

void collect_densities_entropies(Boson sites,
                                 MPS state,
                                 std::string densities_entropies);

void collect_convergence_parameters(Boson sites,
                                    MPS state,
                                    vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                                    std::string convergence_params);

MPS imag_time_evol(Boson sites,
                   MPS state,
                   MPO Hamiltonian,
                   int NoOfSteps = 5,
                   int nosweeps = 4,
                   Real dt_bysweep = 0.02,
                   int MaxBondDim = 200);

MPS dmrg_sequence(Boson sites,
                  MPS state,
                  vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                  int MaxBondDim,
                  std::string densities_entropies,
                  std::string convergence_params,
                  std::string sites_file,
                  std::string mps_file);
