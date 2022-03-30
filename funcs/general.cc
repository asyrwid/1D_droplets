//#include "general.h"
#include "itensor/all.h"
//#include <math.h>
#include <iostream>
#include <fstream>
#include "tdvp.h"
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


MPS annihilate_boson(Boson sites, MPS ket, int site){
  ket.position(site);
  auto newket = noPrime(ket(site) * op(sites, "A", site));
  ket.set(site, newket);
  return ket;
}


MPS create_boson(Boson sites, MPS ket, int site){
  ket.position(site);
  auto newket = noPrime(ket(site) * op(sites, "Adag", site));
  ket.set(site, newket);
  return ket;
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

/*
MPS imag_time_evol(Boson sites,
                   int nosweeps,
                   Real dt_bysweep,
                   int MaxBondDim,
                   int NoOfSteps,
                   MPS state,
                   MPO Hamiltonian,
                   std::string PrA,
                   std::string PrB)
{
  Real tstep = nosweeps*dt_bysweep; // time for one whole round of sweeps

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
    if(maxLinkDim(psi1) < MaxBondDim){ // first, expand basis if the bond dimension does not exceed MaxBondDim
      cut1 = 0;
      std::vector<Real> epsilonK = {1E-10, 1E-10};
      // Play with the numbers if needed. Here, the evolution procedure is just for 'warming up' before performing DMRG
      addBasis(psi, Hamiltonian, epsilonK, {"Cutoff", 1E-10,
                                            "Method", "DensityMatrix",
                                            "KrylovOrd", 3,
                                            "DoNormalize", true,
                                            "Quiet", true});
      }else{cut0 = 1;} // if the bond dimension is already >= MaxBondDim then change cut0 --> 1
                       // so that the basis wont be expanded again
      // TDVP sweep
      if(cut0 == 0){ // evolve when the basis was expanded
        auto energy = tdvp(psi, Hamiltonian, -dt_bysweep, sweepst, {"Truncate", false, // and do not truncate
                                                                    "DoNormalize", true,
                                                                    "Quiet", true,
                                                                    "NumCenter", 1});
      }else{
        if(cut1 == 0){
        auto energy = tdvp(psi, Hamiltonian, -dt_bysweep, sweepst, {"Truncate", true, // truncate in such a case
                                                                    "DoNormalize", true,
                                                                    "Quiet", true,
                                                                    "NumCenter", 1});
        cut1 = 1; // after evolution with the truncation set cut1 --> 1
      }else{ // the case where cut1 = 1
        auto energy = tdvp(psi, Hamiltonian, -dt_bysweep, sweepst, {"Truncate", false, // no truncation
                                                                    "DoNormalize", true,
                                                                    "Quiet", true,
                                                                    "NumCenter", 1});
           }
           }

    printfln("Maximum MPS bond dimension after time evolution is %d", maxLinkDim(psi1));

    vector<double> densities_a;
    vector<double> densities_b;

    tie(densities_a, densities_b) = particle_densities(sites, psi)


          auto psi11 = psi1;
          double nA = 0;
          double nB = 0;
          for(int j = 1; j <= Nsites; j++)
            {
              psi11.position(j);
              auto ket = psi11.A(j);
              auto bra = dag(prime(ket,"Site"));
              auto opN = op(sites,"N",j);
              auto dens = elt(bra*opN*ket);
              if( j%2==0 ){
                nA += dens;
              }else{
                nB += dens;
              }
            }
          printfln("nA = %.5f", nA );
          printfln("nb = %.5f", nB );

          std::ofstream probA0 (PrA, std::ios::app);
          std::ofstream probB0 (PrB, std::ios::app);

          for(int j = 1; j <= Nsites; j++)
            {
              psi11.position(j);
              auto ket = psi11.A(j);
              auto bra = dag(prime(ket,"Site"));
              auto opN = op(sites,"N",j);
              auto dens = elt(bra*opN*ket);
              if( j%2==1 ){
                probA0 << (j+1)/2 << " " << dens << "\n";
              }else{
                probB0 << j/2 << " " << dens << "\n";
              }
            }


            probA0 << " \n";
            probB0 << " \n";
            probA0.close();
            probB0.close();
        }

return psi1;
}

*/
