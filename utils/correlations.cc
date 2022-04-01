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
  vector<vector<double>> one_body_correlations_a = {};
  vector<vector<double>> one_body_correlations_b = {};
  vector<vector<double>> pair_correlations_ab = {};
  vector<vector<double>> density_density_a = {};
  vector<vector<double>> density_density_b = {};

  cout << "\ncalculate correlations:";
  for(int i_site = 1; i_site <= L; i_site++){
    vector<double> row_one_corrs_a = {};
    vector<double> row_one_corrs_b = {};
    vector<double> row_pair_corrs_ab = {};
    vector<double> row_density_density_a = {};
    vector<double> row_density_density_b = {};

    for(int j_site = 1; j_site <= L; j_site++){
      row_one_corrs_a.push_back(one_body_correlation_a(sites, state, i_site, j_site));
      row_one_corrs_b.push_back(one_body_correlation_b(sites, state, i_site, j_site));
      row_pair_corrs_ab.push_back(pair_correlation_ab(sites, state, i_site, j_site));
      row_density_density_a.push_back(density_density_correlation_a(sites, state, i_site, j_site));
      row_density_density_b.push_back(density_density_correlation_b(sites, state, i_site, j_site));
    }cout << "\n" << round(i_site*100./L) << "%" << std::flush;

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
