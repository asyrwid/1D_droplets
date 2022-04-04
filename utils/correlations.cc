#include "itensor/all.h"
#include <string>
#include <iostream>
#include <fstream>
#include "general.h"

using namespace itensor;
using namespace std;


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


tuple<vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>,
      vector<vector<double>>> correlations(Boson sites, MPS state){
// Compute all the correlations
  auto M = length(state);
  int L = M/2;
  vector<vector<double>> one_body_correlations_a(L, vector<double>(L));
  vector<vector<double>> one_body_correlations_b(L, vector<double>(L));
  vector<vector<double>> pair_correlations_ab(L, vector<double>(L));
  vector<vector<double>> density_density_a(L, vector<double>(L));
  vector<vector<double>> density_density_b(L, vector<double>(L));

  cout << "\ncalculate correlations:";
  for(int i_site = 1; i_site <= L; i_site++){
    for(int j_site = 1; j_site <= L; j_site++){
      double adagi_aj = one_body_correlation_a(sites, state, i_site, j_site);
      double bdagi_bj = one_body_correlation_b(sites, state, i_site, j_site);
      double ai_bi_bdagj_adagj = pair_correlation_ab(sites, state, i_site, j_site);
      double nai_naj = density_density_correlation_a(sites, state, i_site, j_site);
      double nbi_nbj = density_density_correlation_b(sites, state, i_site, j_site);
      double nai = density_a(sites, state, i_site);
      double naj = density_a(sites, state, j_site);
      double nbi = density_b(sites, state, i_site);
      double nbj = density_b(sites, state, j_site);

      one_body_correlations_a[i_site-1][j_site-1] = adagi_aj;
      one_body_correlations_b[i_site-1][j_site-1] = bdagi_bj;
      pair_correlations_ab[i_site-1][j_site-1] = ai_bi_bdagj_adagj - adagi_aj*bdagi_bj;
      density_density_a[i_site-1][j_site-1] = nai_naj - nai*naj;
      density_density_b[i_site-1][j_site-1] = nbi_nbj - nbi*nbj;
    }cout << "\n" << round(i_site*100./L) << "%" << std::flush;
  }
  return make_tuple(one_body_correlations_a,
                    one_body_correlations_b,
                    pair_correlations_ab,
                    density_density_a,
                    density_density_b);
}


void save_correlations(Boson sites, MPS state, string path){

  vector<vector<double>> one_body_correlations_a;
  vector<vector<double>> one_body_correlations_b;
  vector<vector<double>> pair_correlations_ab;
  vector<vector<double>> density_density_a;
  vector<vector<double>> density_density_b;

  tie(one_body_correlations_a,
      one_body_correlations_b,
      pair_correlations_ab,
      density_density_a,
      density_density_b) = correlations(sites, state);

  int n_rows = one_body_correlations_a.size();
  int n_cols = n_rows;
  std::ofstream file (path);
  // set labels
  file << "row" << " " << "col" << " ";
  file << "one_body_correlations_a" << " ";
  file << "one_body_correlations_b" << " ";
  file << "pair_correlations_ab" << " ";
  file << "density_density_a" << " ";
  file << "density_density_a" << "\n";

  for(int row = 0; row < n_rows; row++){
    for(int col = 0; col < n_cols; col++){
      file << row+1 << " " << col+1 << " ";
      file << one_body_correlations_a[row][col] << " ";
      file << one_body_correlations_b[row][col] << " ";
      file << pair_correlations_ab[row][col] << " ";
      file << density_density_a[row][col] << " ";
      file << density_density_b[row][col] << "\n";
    }
  }
  file.close();
}
