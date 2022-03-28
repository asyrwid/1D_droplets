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
    int j_site_boson_a = 2*j_site - 1;
    int j_site_boson_b = 2*j_site;
    hopping_a += -t, "Adag", j_site_boson_a, "A", j_site_boson_a + 2;
    hopping_a += -t, "Adag", j_site_boson_a + 2, "A", j_site_boson_a;
    hopping_b += -t, "Adag", j_site_boson_b, "A", j_site_boson_b + 2;
    hopping_b += -t, "Adag", j_site_boson_b + 2, "A", j_site_boson_b;
  }
  auto H_hop_a = toMPO(hopping_a);
  auto H_hop_b = toMPO(hopping_b);

  // intra-component interactions
  auto interaction_aa = AutoMPO(sites);
  auto interaction_bb = AutoMPO(sites);
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_boson_a = 2*j_site - 1;
    int j_site_boson_b = 2*j_site;
    interaction_aa +=  0.5*U, "N", j_site_boson_a, "N", j_site_boson_a;
    interaction_aa += -0.5*U, "N", j_site_boson_a;
    interaction_bb +=  0.5*U, "N", j_site_boson_b, "N", j_site_boson_b;
    interaction_bb += -0.5*U, "N", j_site_boson_b;
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
    int j_site_boson_a = 2*L - 1;
    int j_site_boson_b = 2*L;
    interaction_ab += U_ab, "N", j_site_boson_a, "N", j_site_boson_b;
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

  auto state = InitState(sites);
  int M_sites = length(sites);

  for(int i = 1; i <= M_sites; ++i)
  { state.set(i,"0"); }

  int Q = 2*M_sites - Natoms;
  for(int s = Q; s < Q + 2*Natoms;  s+=2){
    state.set(s-1,"1");
    state.set(s,"1");
  }
  auto psi = MPS(state);
  return psi;
}

MPS annihilate_boson(Boson sites, MPS ket, int site){

  ket.position(site);
  auto newket = noprime(ket(site) * op(sites, "A", site));
  ket.set(site, newket);

  return ket;

}

MPS create_boson(Boson sites, MPS ket, int site){

  ket.position(site);
  auto newket = noprime(ket(site) * op(sites, "Adag", site));
  ket.set(site, newket);

  return ket;

}
