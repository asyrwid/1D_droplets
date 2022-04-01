#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include <string>
#include <iostream>
#include <fstream>
//#include <complex>

using namespace itensor;
using namespace std;

tuple<MPO, MPO, MPO, MPO, MPO, MPO, MPO> get_H(SiteSet& sites,
                                               double t,
                                               double U,
                                               double U_ab){
  int M_sites = length(sites);
  int L = M_sites/2;

  // kinetic terms
  auto hopping_a = AutoMPO(sites);
  auto hopping_b = AutoMPO(sites);
  for (int j_site = 1; j_site < L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    hopping_a += -t, "Adag", j_site_a, "A", j_site_a + 2;
    hopping_a += -t, "Adag", j_site_a + 2, "A", j_site_a;
    hopping_b += -t, "Adag", j_site_b, "A", j_site_b + 2;
    hopping_b += -t, "Adag", j_site_b + 2, "A", j_site_b;
  }
  auto H_hop_a = toMPO(hopping_a);
  auto H_hop_b = toMPO(hopping_b);

  // intra-component interactions
  auto interaction_aa = AutoMPO(sites);
  auto interaction_bb = AutoMPO(sites);
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    interaction_aa +=  0.5*U, "N", j_site_a, "N", j_site_a;
    interaction_aa += -0.5*U, "N", j_site_a;
    interaction_bb +=  0.5*U, "N", j_site_b, "N", j_site_b;
    interaction_bb += -0.5*U, "N", j_site_b;
  }
  auto H_aa = toMPO(interaction_aa);
  auto H_bb = toMPO(interaction_bb);

/* REDUNDANT !!!
  if(PBC == 1){
    int j_site_boson_a = 2*L - 1;
    int j_site_boson_b = 2*L;
    tmp_H_a_t += -t, "Adag", j_site_boson_a, "A", 1;
    tmp_H_a_t += -t, "Adag", 1, "A", j_site_boson_a;

    tmp_H_b_t += -t,"Adag", j_site_boson_b, "A", 2;
    tmp_H_b_t += -t,"Adag", 2, "A", j_site_boson_b;
  }
*/

  // inter-component interactions
  auto interaction_ab = AutoMPO(sites);
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_a = 2*j_site - 1;
    int j_site_b = 2*j_site;
    interaction_ab += U_ab, "N", j_site_a, "N", j_site_b;
  }
  auto H_ab = toMPO(interaction_ab);

  // Hard wall at the edges choose value different than 0 if you want to
  // have additional potential at the first and last sites in the system.
  auto edge_potential = AutoMPO(sites);
  edge_potential += 0, "N", 1;
  edge_potential += 0, "N", M_sites-1;
  edge_potential += 0, "N", 2;
  edge_potential += 0, "N", M_sites;
  auto H_edge = toMPO(edge_potential);

  auto H_a = sum(H_hop_a, H_aa);
  auto H_b = sum(H_hop_b, H_bb);
  auto H_total = sum(sum(sum(H_a, H_b), H_ab), H_edge);

  return make_tuple(H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab, H_edge);
}


MPS initial_state(Boson sites, int Natoms){
// Natoms = number of atoms in each component assuming that we have Na = Nb
  auto state = InitState(sites);
  int M_sites = length(sites);

  for(int i = 1; i <= M_sites; ++i)
  { state.set(i,"0"); }

  int Q = int(M_sites/2.) - Natoms;
  for(int s = Q; s < Q + 2*Natoms;  s += 2){
    state.set(s-1, "1");
    state.set(s, "1");
  }
  auto psi = MPS(state);
  return psi;
}


string parameters_to_filename(int Natoms,
                              int L,
                              int maxOccupation,
                              double t,
                              double U,
                              double r,
                              int MaxBondDim){

string str_Natoms = "_Natoms=" + str(Natoms);
string str_L = "_L=" + str(L);
string str_maxOccupation = "_MaxOcc=" + str(maxOccupation);
string str_MaxBondDim = "_MaxBondDim=" + str(MaxBondDim);

double t1, t2, U1, U2;
t2 = std::modf(t, &t1);
U2 = std::modf(U, &U1);
string str_t = "_t=" + str(t1) + "_" + str( floor(fabs(t2)*1000) );
string str_U = "_U=" + str(U1) + "_" + str( floor(fabs(U2)*1000) );

string str_r = "_r=0_"; // we consider |r| < 1
if(r*1000<100){ str_r = str_r + "0" + str(r*1000); }
else{ str_r = str_r + str(r*1000); }

return str_Natoms + str_L + str_t + str_U + str_r + str_maxOccupation + str_MaxBondDim + ".txt";
}


MPS apply_operator(Boson sites, string OpName, MPS ket, int site){
  ket.position(site);
  auto newket = noPrime(ket(site) * op(sites, OpName, site));
  ket.set(site, newket);
  return ket;
}


double density_a(Boson sites, MPS state, int site){
  auto M = length(state);
  int L = M/2;
  int site_a = 2*site - 1;
  auto ket = apply_operator(sites, "N", state, site_a);
  auto density_boson_a = innerC(state, ket);
  return density_boson_a.real();
}


double density_b(Boson sites, MPS state, int site){
  auto M = length(state);
  int L = M/2;
  int site_b = 2*site;
  auto ket = apply_operator(sites, "N", state, site_b);
  auto density_boson_b = innerC(state, ket);
  return density_boson_b.real();
}


double density_density_correlation_a(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_a = 2*i_site - 1;
  int j_site_a = 2*j_site - 1;
  auto ket = apply_operator(sites, "N", state, j_site_a);
       ket = apply_operator(sites, "N", ket, i_site_a);
  auto N_a_i_N_a_j = innerC(state, ket);
  return N_a_i_N_a_j.real();
}


double density_density_correlation_b(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_b = 2*i_site;
  int j_site_b = 2*j_site;
  auto ket = apply_operator(sites, "N", state, j_site_b);
       ket = apply_operator(sites, "N", ket, i_site_b);
  auto N_b_i_N_b_j = innerC(state, ket);
  return N_b_i_N_b_j.real();
}


double one_body_correlation_a(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_a = 2*i_site - 1;
  int j_site_a = 2*j_site - 1;
  auto ket = apply_operator(sites, "A", state, j_site_a);
  auto bra = apply_operator(sites, "A", state, i_site_a);
  auto a_dag_i_a_j = innerC(bra, ket);
  return a_dag_i_a_j.real();
}


double one_body_correlation_b(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_b = 2*i_site;
  int j_site_b = 2*j_site;
  auto ket = apply_operator(sites, "A", state, j_site_b);
  auto bra = apply_operator(sites, "A", state, i_site_b);
  auto b_dag_i_b_j = innerC(bra, ket);
  return b_dag_i_b_j.real();
}


double pair_correlation_ab(Boson sites, MPS state, int i_site, int j_site){
  auto M = length(state);
  int L = M/2;
  int i_site_a = 2*i_site - 1;
  int j_site_a = 2*j_site - 1;
  int i_site_b = 2*i_site;
  int j_site_b = 2*j_site;

  auto ket_tmp = apply_operator(sites, "Adag", state, j_site_a);
  auto ket = apply_operator(sites, "Adag", ket_tmp, j_site_b);
  auto bra_tmp = apply_operator(sites, "Adag", state, i_site_b);
  auto bra = apply_operator(sites, "Adag", bra_tmp, i_site_a);
  auto a_dag_b_dag_b_a = innerC(bra, ket);
  return a_dag_b_dag_b_a.real();
}


double entanglement_entropy(MPS state, int lattice_site){
  auto psi = state;
  int bond = lattice_site;
  psi.position(bond);
  auto l = leftLinkIndex(psi, bond);
  auto s = siteIndex(psi, bond);
  auto [U,S,V] = svd(psi(bond), {l, s});
  auto u = commonIndex(U, S);
  double SvN = 0.;
  for(auto n : range1(dim(u)))
  {
    auto Sn = elt(S, n, n);
    auto p = sqr(Sn);
    if(p > 1E-12){ SvN += -p*log(p); }
  }
  // std::cout << "\n\nS at bond: " << bond << " = " << SvN;
  auto entanglement_entropy = SvN;
  return entanglement_entropy;
}


tuple<vector<double>, vector<double>> entropies(MPS state){
  auto M = length(state);
  int L = M/2;
  vector<double> entropies_a = {};
  vector<double> entropies_b = {};
  for(int site = 1; site <= L; site++){
    int site_a = 2*site - 1;
    int site_b = 2*site;
    if((site == 1) || (site == L)){
      entropies_a.push_back(0.);
      entropies_b.push_back(0.);
    }else{
      entropies_a.push_back(entanglement_entropy(state, site_a));
      entropies_b.push_back(entanglement_entropy(state, site_b));
    }
  }
  return make_tuple(entropies_a, entropies_b);
}


tuple<vector<double>, vector<double>> particle_densities(Boson sites, MPS state){
  auto M = length(state);
  int L = M/2;
  vector<double> densities_a = {};
  vector<double> densities_b = {};
  for(int site = 1; site <= L; site++){
    densities_a.push_back(density_a(sites, state, site));
    densities_b.push_back(density_b(sites, state, site));
  }
  return make_tuple(densities_a, densities_b);
}


tuple<vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>> correlations(Boson sites, MPS state){
// Compute all the correlations
  auto M = length(state);
  int L = M/2;
  vector<vector<double>> one_body_correlations_a = {};
  vector<vector<double>> one_body_correlations_b = {};
  vector<vector<double>> pair_correlations_ab = {};
  vector<vector<double>> density_density_a = {};
  vector<vector<double>> density_density_b = {};

  for(int i_site = 1; i_site <= L; i_site++){
    vector<double> row_one_corrs_a = {};
    vector<double> row_one_corrs_b = {};
    vector<double> row_pair_corrs_ab = {};
    vector<double> row_density_density_a = {};
    vector<double> row_density_density_b = {};

    for(int j_site = 1; j_site <= L; j_site++){
      row_one_corrs_a.push_back(one_body_correlation_a(sites, state, i_site, j_site));
      row_one_corrs_b.push_back(one_body_correlation_a(sites, state, i_site, j_site));
      row_pair_corrs_ab.push_back(pair_correlation_ab(sites, state, i_site, j_site));
      row_density_density_a.push_back(density_density_correlation_a(sites, state, i_site, j_site));
      row_density_density_b.push_back(density_density_correlation_b(sites, state, i_site, j_site));
    }
    one_body_correlations_a.push_back(row_one_corrs_a);
    one_body_correlations_b.push_back(row_one_corrs_b);
    pair_correlations_ab.push_back(row_pair_corrs_ab);
    density_density_a.push_back(row_density_density_a);
    density_density_b.push_back(row_density_density_b);
  }
  return make_tuple(one_body_correlations_a,
                    one_body_correlations_b,
                    pair_correlations_ab,
                    density_density_a,
                    density_density_b);
}


void prepare_file(vector<std::string> column_names, std::string path){
  // open file and fill the first row with names of variables (column_names)
  std::ofstream file (path);
  int n_variables = column_names.size();
  for(int i = 0; i < n_variables; i++){
    file << column_names[i] << " ";
  }file << "\n";
  file.close();
}


void collect_densities_entropies(Boson sites,
                                 MPS state,
                                 std::string densities_entropies){
  vector<double> densities_a;
  vector<double> densities_b;
  tie(densities_a, densities_b) = particle_densities(sites, state);
  vector<double> entropies_a;
  vector<double> entropies_b;
  tie(entropies_a, entropies_b) = entropies(state);

  // print total densities
  double N_a = 0;
  double N_b = 0;
  for(int i = 0; i <= int(length(state)/2)-1; i++){
    N_a += densities_a[i];
    N_b += densities_b[i];
  }
  printfln("Na = %.5f", N_a );
  printfln("Nb = %.5f", N_b );

  // write densities and entropies to file
  std::ofstream dens_entrs (densities_entropies, std::ios::app);
  for(int i = 0; i <= int(length(state)/2)-1; i++){
    int site_number = i+1;
    dens_entrs << site_number << " ";
    dens_entrs << densities_a[i] << " " << densities_b[i] << " ";
    dens_entrs << entropies_a[i] << " " << entropies_b[i] << "\n";
  }
  dens_entrs << "\n"; // separate results obtained in consecutive steps
  dens_entrs.close();
}


void collect_convergence_parameters(Boson sites,
                                    MPS state,
                                    vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                                    std::string convergence_params){
  int L = int(length(state)/2);
  int central_a = L;
  int central_b = L;
  if(L % 2 == 0){ central_a = L - 1; }
  else{ central_b = L + 1; }

  std::ofstream conv_params (convergence_params, std::ios::app);

  // entropy at central bond
  double central_entropy_a = entanglement_entropy(state, central_a);
  double central_entropy_b = entanglement_entropy(state, central_b);
  conv_params << central_entropy_a << " " << central_entropy_b << " ";

  //calculate energies
  double E_total = innerC(state, H_terms[0], state).real();
  double E_hop_a = innerC(state, H_terms[1], state).real();
  double E_hop_b = innerC(state, H_terms[2], state).real();
  double E_aa = innerC(state, H_terms[3], state).real();
  double E_bb = innerC(state, H_terms[4], state).real();
  double E_ab = innerC(state, H_terms[5], state).real();
  conv_params << E_total << " ";
  conv_params << E_hop_a << " ";
  conv_params << E_hop_b << " ";
  conv_params << E_aa << " ";
  conv_params << E_bb << " ";
  conv_params << E_ab << " ";
  conv_params << "\n";
  conv_params.close();
}


MPS imag_time_evol(Boson sites,
                   MPS state,
                   MPO H_total,
                   int NoOfSteps,
                   int nosweeps,
                   Real dt_bysweep,
                   int MaxBondDim)
{
  Real tstep = nosweeps*dt_bysweep; // time evolved during the single round of sweeps

  auto sweepst = Sweeps(nosweeps);
  sweepst.maxdim() = MaxBondDim;
  sweepst.cutoff() = 1E-18;
  sweepst.niter() = 50;

  int cut0 = 0;
  int cut1 = 0;
  auto psi = state;

  for(int n = 1; n <= NoOfSteps ; ++n){
    printf("\n ================= \n n = ", n);
    printf("\n time form ", (n-1)*tstep, " to ", n*tstep, "\n ================= \n");
    if(maxLinkDim(psi) < MaxBondDim){ // first, expand basis if the bond dimension does not exceed MaxBondDim
      cut1 = 0;
      std::vector<Real> epsilonK = {1E-10, 1E-10};
      // Play with the numbers if needed. Here, the evolution procedure is just for 'warming up' before performing DMRG
      addBasis(psi, H_total, epsilonK, {"Cutoff", 1E-10,
                                        "Method", "DensityMatrix",
                                        "KrylovOrd", 3,
                                        "DoNormalize", true,
                                        "Quiet", true});
      }else{cut0 = 1;} // if the bond dimension is already >= MaxBondDim then change cut0 --> 1
                       // so that the basis wont be expanded again
      // TDVP sweep
      if(cut0 == 0){ // evolve when the basis was expanded
        auto energy = tdvp(psi, H_total, -dt_bysweep, sweepst, {"Truncate", false, // and do not truncate
                                                                "DoNormalize", true,
                                                                "Quiet", true,
                                                                "NumCenter", 1});
      }else{
        if(cut1 == 0){
        auto energy = tdvp(psi, H_total, -dt_bysweep, sweepst, {"Truncate", true, // truncate in such a case
                                                                "DoNormalize", true,
                                                                "Quiet", true,
                                                                "NumCenter", 1});
        cut1 = 1; // after evolution with the truncation set cut1 --> 1
      }else{ // the case where cut1 = 1
        auto energy = tdvp(psi, H_total, -dt_bysweep, sweepst, {"Truncate", false, // no truncation
                                                                "DoNormalize", true,
                                                                "Quiet", true,
                                                                "NumCenter", 1});
           }
           }

    printfln("Maximum MPS bond dimension after time evolution is %d", maxLinkDim(psi));
  }

  return psi;
}


void dmrg_sequence(Boson sites,
                   MPS state,
                   vector<MPO> H_terms, // {H_total, H_hop_a, H_hop_b, H_aa, H_bb, H_ab}
                   int MaxBondDim,
                   std::string densities_entropies,
                   std::string convergence_params,
                   std::string sites_file,
                   std::string mps_file){

  MPO H_total = H_terms[0];

  int dim = 64;
  vector<int> BondDim_truncation = {};
  while(dim < 0.75*MaxBondDim){
    BondDim_truncation.push_back(dim);
    dim = 2*dim;
  }
  BondDim_truncation.push_back(MaxBondDim);

  // Initial DMRG sweeps
  auto psi = state;
  for(int n = 0; n < BondDim_truncation.size(); n++){
    int nswep = 4;
    auto sweeps = Sweeps(nswep);
    sweeps.maxdim() = BondDim_truncation[n];
    sweeps.cutoff() = 1E-16;
    for(int m = 0; m < 10; m++){
      std::cout << "\n ======= step " << m+1 << "   truncation = " << BondDim_truncation[n] << "\n";
      auto [energy0, psi0] = dmrg(H_total, psi, sweeps, {"Quiet",true});
      psi = psi0;
    }
    // collect densities and entropies
    collect_densities_entropies(sites, psi, densities_entropies);
    // collect convergence_parameters
    collect_convergence_parameters(sites, psi, H_terms, convergence_params);
  }

  // calculate entropies at central bonds
  int L = int(length(state)/2);
  int central_a = L;
  int central_b = L;
  if(L % 2 == 0){ central_a = L - 1; }
  else{ central_b = L + 1; }
  double central_entropy_a = entanglement_entropy(state, central_a);
  double central_entropy_b = entanglement_entropy(state, central_b);
  double postdmrg_central_entropy_a = 0.1*central_entropy_a; // initialization
  double postdmrg_central_entropy_b = 0.1*central_entropy_b; // initialization
  double relative_diff_a = fabs(central_entropy_a - postdmrg_central_entropy_a)/central_entropy_a;
  double relative_diff_b = fabs(central_entropy_b - postdmrg_central_entropy_b)/central_entropy_b;


  while((relative_diff_a > 0.001) || (relative_diff_b > 0.001)){
    int nswep = 4;
    auto sweeps = Sweeps(nswep);
    sweeps.maxdim() = MaxBondDim;
    sweeps.cutoff() = 1E-16;
    for(int m = 0; m < 10; m++){
      std::cout << "\n ======= step " << m+1 << "   max truncation = " << MaxBondDim << "\n";
      auto [energy0, psi0] = dmrg(H_total, psi, sweeps, {"Quiet",true});
      psi = psi0;
    }
    // collect densities and entropies
    collect_densities_entropies(sites, psi, densities_entropies);
    // collect convergence_parameters
    collect_convergence_parameters(sites, psi, H_terms, convergence_params);

    // save sites and mps`
    writeToFile(sites_file, sites);
    writeToFile(mps_file, psi);

    // recalculate entropy difference
    central_entropy_a = postdmrg_central_entropy_a;
    central_entropy_b = postdmrg_central_entropy_b;

    postdmrg_central_entropy_a = entanglement_entropy(state, central_a);
    postdmrg_central_entropy_b = entanglement_entropy(state, central_b);

    relative_diff_a = fabs(central_entropy_a - postdmrg_central_entropy_a)/central_entropy_a;
    relative_diff_b = fabs(central_entropy_b - postdmrg_central_entropy_b)/central_entropy_b;
  }

  printfln("Ground State Found!");


}
