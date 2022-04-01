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
      row_pair_corrs_ab.push_back(pair_correlation_ab(sites, state, i_site, j_site) - one_body_correlation_a(sites, state, i_site, j_site)*one_body_correlation_b(sites, state, i_site, j_site) );
      row_density_density_a.push_back(density_density_correlation_a(sites, state, i_site, j_site) - one_body_correlation_a(sites, state, i_site, i_site)*one_body_correlation_a(sites, state, j_site, j_site));
      row_density_density_b.push_back(density_density_correlation_b(sites, state, i_site, j_site) - one_body_correlation_b(sites, state, i_site, i_site)*one_body_correlation_b(sites, state, j_site, j_site));
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



/// 2022.04.01 11:54


double get_string_parity_a(MPS psi, Boson sites){
  // a/ https://arxiv.org/pdf/0803.2851.pdf
  // b/ https://arxiv.org/pdf/2003.07556.pdf
  // c/ https://arxiv.org/pdf/2112.10386.pdf
  double pi = 3.14159265359;
  double string_parity_a = 0;

  int M_sites = length(sites);
  int L = M_sites/2; 

  auto psi2 = psi;
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_boson_a = 2*j_site - 1;
    double N_a = one_body_correlation_a(sites, psi, j_site, j_site);
    auto delta_N_hat = sites.op("N", j_site_boson_a) - N_a*sites.op("Id", j_site_boson_a);
    auto O_pair = expHermitian(delta_N_hat, Cplx_i * pi);   
    
    psi2.position(j_site_boson_a);
    auto psi2_tmp = noPrime(psi2(j_site_boson_a) * O_pair);
    psi2.set(j_site_boson_a, psi2_tmp);
    psi2.noPrime();
  }  
  string_parity_a = real(innerC(psi,psi2));


  return string_parity_a;
}

tuple<double, double, double, double> get_Fourier_transforms_and_field_operators_expectation_values(Boson sites,
                                                             MPS state,
                                                             vector<vector<double>> one_body_correlations_a,
                                                             vector<vector<double>> one_body_correlations_b,
                                                             vector<vector<double>> pair_correlations_ab,
                                                             vector<vector<double>> density_density_correlations_a,
                                                             vector<vector<double>> density_density_correlations_b,
                                                             string filename_parameters){

  auto M = length(state);
  int L = M/2;
  ///////////////////////////////////////////////////////////////////////////////////////////////
    //                                 Defined correlations
    // rho_a(i,j)    = <a_dag_j a_i>                                                one body density matrix for boson a
    // rho_b(i,j)    = <b_dag_j b_i>                                                one body density matrix for boson b
    // kappa_a(i,j)  = <N_a_i N_a_j> - <N_a_i><N_a_j>                               density-density correlations for boson a
    // kappa_b(i,j)  = <N_b_i N_b_j> - <N_b_i><N_b_j>                               density-density correlations for boson b
    // gamma_ab(i,j) = <a_dag_j b_dag_j b_i a_i>  - <a_dag_j a_i> <b_dag_j b_i>     two-body density matrix for pair of bosons (a,b)
    //
    // h_a  = 1/(L*N_a)<sum_{i!=j_0} rho_a(i, j0)       // expectation value of boson a hopping
    // h_b  = 1/(L*N_b)<sum_{i!=j_0} rho_b(i, j0)       // expectation value of boson b hopping
    // h_aa = 1/(L*N_a*N_b)<sum_{i!=j_0} gamma_ab(i,j0) // expectation value of pair hopping
    // f_SF_a = 1/(L*N_a)<sum_{i,j} rho_a(i,j)          // superfluid fraction for boson a
    // f_SF_b = 1/(L*N_b)<sum_{i,j} rho_b(i,j)          // superfluid fraction for boson b
    double h_a  = 0;
    double h_b  = 0;
    double h_ab = 0;
    double f_SF_a = 0;
    double f_SF_b = 0;
    double N_a = 0;
    double N_b = 0;
    for(int i_site = 1; i_site<=L; i_site++){
      for(int j_site = 1; j_site<=L; j_site++){
        if(j_site != L/2){
          h_a  = h_a  + one_body_correlations_a[i_site - 1][j_site - 1];
          h_b  = h_b  + one_body_correlations_b[i_site - 1][j_site - 1];
          h_ab = h_ab + pair_correlations_ab[i_site - 1][j_site - 1];
        }
          f_SF_a = f_SF_a + one_body_correlations_a[i_site - 1][j_site - 1];
          f_SF_b = f_SF_b + one_body_correlations_b[i_site - 1][j_site - 1];
          if(i_site == j_site){
            N_a = N_a + one_body_correlations_a[i_site - 1][j_site - 1];
            N_b = N_b + one_body_correlations_a[i_site - 1][j_site - 1];
          }
      }
    }
    f_SF_a = f_SF_a/(L*N_a);
    f_SF_b = f_SF_b/(L*N_b);
    h_a = h_a/(L*N_a);
    h_b = h_b/(L*N_b);
    h_ab = h_ab/L/(N_a*N_b);

    // Calculate structure factor S(k) and R(k) ( equiv to M(k) in arxiv:2112.10386) and Q(k) = FourierTransform[pair_correlations_ab]
    // Knowing that all correlations depends on |i-j|, we can calculate one-dimensional Fourier transform of matrix A_ij, as F[k] = 1/L^2\sum_{i,j} e^{i k |i-j|} A_ij
    // 1. S(k) = FourierTransform[density_dinsity_correlations_a]
    // 2. M(k) = FourierTransform[one_body_correlations_a];
    // 3. Q(k) = FourierTransform[pair_correlations_ab]; 
    double pi = 3.14159265359;
    vector<complex<double>> S(L, 0);
    vector<complex<double>> R(L, 0);
    vector<complex<double>> Q(L, 0);
    for (int k_int = 0; k_int<L; k_int++){
        double k = k_int*2.0*pi/L;
        for(int i = 0; i<L; i++){
            for(int j = 0; j<L; j++){
                S[k_int] = S[k_int] + exp(Cplx_i*k*(i-j))*density_density_correlations_a[i][j]/pow(L,2.0);
                R[k_int] = R[k_int] + exp(Cplx_i*k*(i-j))*one_body_correlations_a[i][j]/pow(L,2.0);
                Q[k_int] = Q[k_int] + exp(Cplx_i*k*(i-j))*pair_correlations_ab[i][j]/pow(L,2.0);
            }
        }
    }

  ofstream   file_data;
  string string_filename_data = "data_diagonal_" + filename_parameters + ".txt"; 
  file_data.open(string_filename_data,ios::out);
  for(int j_site = 1; j_site<=L; j_site++){
    file_data <<j_site<<" "<<one_body_correlations_a[j_site-1][j_site-1]<<" "<<sqrt(density_density_correlations_a[j_site-1][j_site - 1])<<" ";
    file_data <<j_site<<" "<<one_body_correlations_b[j_site-1][j_site-1]<<" "<<sqrt(density_density_correlations_b[j_site-1][j_site - 1])<<" ";
    file_data <<j_site<<" "<<(S[j_site-1]).real()<<" ";
    file_data <<j_site<<" "<<(R[j_site-1]).real()<<" ";
    file_data <<j_site<<" "<<(Q[j_site-1]).real()<<" "<<endl;
  }
  file_data.close();

 
  double string_parity_a = get_string_parity_a(state, sites);

  return make_tuple(f_SF_a, h_a, h_ab, string_parity_a);
}