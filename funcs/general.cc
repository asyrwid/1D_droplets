#include "general.h"
#include "itensor/all.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>

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
  for(int s = Q; s < Q + 2*Natoms;  s+=2){
    state.set(s-1,"1");
    state.set(s,"1");
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
