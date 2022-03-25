#include "tdvp.h"
#include "basisextension.h"
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "itensor/mps/sites/boson.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>   
#include <boost/math/special_functions/factorials.hpp>

using namespace itensor;
using namespace std;


tuple<MPO,MPO,MPO,MPO,MPO,MPO,MPO> get_H(SiteSet& sites, int PBC, double t, double U, double U_ab){

  int M_sites = length(sites);
  int L = M_sites/2; 

  auto tmp_H_aa_U = AutoMPO(sites);
  auto tmp_H_bb_U = AutoMPO(sites);   
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_boson_a = 2*j_site - 1;
    int j_site_boson_b = 2*j_site;
    tmp_H_aa_U +=  0.5*U, "N",j_site_boson_a, "N",j_site_boson_a;   
    tmp_H_aa_U += -0.5*U, "N",j_site_boson_a;

    tmp_H_bb_U +=  0.5*U, "N",j_site_boson_b, "N",j_site_boson_b;   
    tmp_H_bb_U += -0.5*U, "N",j_site_boson_b;
  }    
  auto H_aa_U = toMPO(tmp_H_aa_U);
  auto H_bb_U = toMPO(tmp_H_bb_U);

  auto tmp_H_a_t = AutoMPO(sites);
  auto tmp_H_b_t = AutoMPO(sites);  
  for (int j_site = 1; j_site < L; j_site++){
    int j_site_boson_a = 2*j_site - 1;
    int j_site_boson_b = 2*j_site;
    tmp_H_a_t += -t, "Adag", j_site_boson_a, "A", j_site_boson_a + 2;  
    tmp_H_a_t += -t, "Adag", j_site_boson_a + 2, "A", j_site_boson_a;   

    tmp_H_b_t += -t,"Adag", j_site_boson_b    , "A", j_site_boson_b + 2; 
    tmp_H_b_t += -t,"Adag", j_site_boson_b + 2, "A", j_site_boson_b;   
  }
  if(PBC == 1){
    int j_site_boson_a = 2*L - 1;
    int j_site_boson_b = 2*L;    
    tmp_H_a_t += -t, "Adag", j_site_boson_a, "A", 1;
    tmp_H_a_t += -t, "Adag", 1, "A", j_site_boson_a;

    tmp_H_b_t += -t,"Adag", j_site_boson_b, "A", 2;
    tmp_H_b_t += -t,"Adag", 2, "A", j_site_boson_b;
  }
  auto H_a_t = toMPO(tmp_H_a_t);
  auto H_b_t = toMPO(tmp_H_b_t);  
 
  auto tmp_H_ab_U = AutoMPO(sites);    
  for (int j_site = 1; j_site <= L; j_site++){
    int j_site_boson_a = 2*L - 1;
    int j_site_boson_b = 2*L;       
    tmp_H_ab_U += U_ab,"N",j_site_boson_a,"N",j_site_boson_b;  
  }

    auto H_ab_U = toMPO(tmp_H_ab_U); 

    //Hard wall at the edges
    auto tmp_H_potential = AutoMPO(sites);
    tmp_H_potential += 0, "N", 1;
    tmp_H_potential += 0, "N", M_sites-1;
    tmp_H_potential += 0, "N", 2;
    tmp_H_potential += 0, "N", M_sites;
    auto H_potential = toMPO(tmp_H_potential);


    auto H_a= sum(H_a_t, H_aa_U);
    auto H_b = sum(H_b_t, H_bb_U);
    auto H_total     = sum(sum(sum(H_a, H_b),H_ab_U),H_potential);
  return make_tuple(H_total, H_a_t, H_aa_U, H_b_t, H_bb_U, H_ab_U, H_potential);
}

 

string get_string(double g){
	char g_str_format[40];
  sprintf(g_str_format,"%6.6f",g);
	return g_str_format;
}


 
 



double get_entanglement_entropy(MPS& psi, SiteSet& sites, int lattice_site, string filename_parameters)
{
  auto M = length(psi);
  int L = M/2;
  auto psi11 = psi;
  int b1 = lattice_site;
  psi11.position(b1);
  auto l11 = leftLinkIndex(psi11,b1);
  auto s11 = siteIndex(psi11,b1);
  auto [U1,S1,V1] = svd(psi11(b1),{l11,s11});
  auto u11 = commonIndex(U1,S1);
  double SvN11 = 0.;
  for(auto n : range1(dim(u11)))
  {
    auto Sn = elt(S1,n,n);
    auto p = sqr(Sn);
    if(p > 1E-12){ SvN11 += -p*log(p); }
  }
  std::cout << "\n\nS at bond: " << b1 << " = " << SvN11;
  auto entanglement_entropy = SvN11;
  
  ofstream file_entanglement;
  file_entanglement.open(filename_parameters, ios::out);
  for(int i_site = 1; i_site<L; i_site++){
    int i_boson_a = 2*i_site - 1;

    auto psi11 = psi;
    int b1 = i_boson_a;
    psi11.position(b1);
    auto l11 = leftLinkIndex(psi11,b1);
    auto s11 = siteIndex(psi11,b1);
    auto [U1,S1,V1] = svd(psi11(b1),{l11,s11});
    auto u11 = commonIndex(U1,S1);
    double SvN11 = 0.;
    for(auto n : range1(dim(u11)))
    {
      auto Sn = elt(S1,n,n);
      auto p = sqr(Sn);
      if(p > 1E-12){ SvN11 += -p*log(p); }
    }
    file_entanglement<<i_site<<" "<<SvN11<<endl;

  }

  file_entanglement.close();

  return entanglement_entropy;
}




tuple<double, double, double> get_observables(MPS& psi, SiteSet& sites,  string filename_parameters){
  //Calculate:
  // 1. densities
  // 2. correlations
  // 3. Superfluid fraction SF, expectation value of the field operator <a_i0>, and Pairsuperfluid value (SF_val, a_i0_val, PSF_val)
  auto M = length(psi);
  int L = M/2;
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //                                 Measure densities
  vector<complex<double>> density_boson_a(L), density_boson_b(L); 
  vector<complex<double>> std_density_boson_a(L), std_density_boson_b(L);
  auto psi_tmp = psi;
  for (int j_site = 1; j_site <= L; j_site++){
      int j_site_boson_a = 2*j_site - 1;
      int j_site_boson_b = 2*j_site;

      psi.position(j_site_boson_a);
      density_boson_a[j_site-1] = eltC(dag(prime(psi(j_site_boson_a), "Site")) * op(sites, "N", j_site_boson_a) * psi(j_site_boson_a));
      complex<double> na = density_boson_a[j_site-1];
      std_density_boson_a[j_site-1] = sqrt( eltC(dag(prime(psi(j_site_boson_a), "Site")) * op(sites, "N*N", j_site_boson_a) * psi(j_site_boson_a)) - na*na );
      
      psi.position(j_site_boson_b);
      density_boson_b[j_site-1] = eltC(dag(prime(psi(j_site_boson_b), "Site")) * op(sites, "N", j_site_boson_b) * psi(j_site_boson_b));
      complex<double> nb = density_boson_b[j_site-1];
      std_density_boson_b[j_site-1] = sqrt( eltC(dag(prime(psi(j_site_boson_b), "Site")) * op(sites, "N*N", j_site_boson_b) * psi(j_site_boson_b)) - nb*nb);
  }

  std::complex<double> N_a = 0;
  std::complex<double> N_b = 0;
 
    
  for (int j = 0; j < M/2; ++j){
      N_a += density_boson_a[j];
      N_b += density_boson_a[j];
  }
 
  printfln("Total a: %.10f", N_a.real());
  printfln("Total b: %.10f", N_b.real());
  println();


  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //                                 Measure correlations

    int i_0 = L/2;
    int i0_site_boson_a = 2*i_0 - 1;
    int i0_site_boson_b = 2*i_0;


    double SF_val = 0;
    vector<vector<double>> corr_one_body_boson_a_matrix(L, vector<double> (L,0));
    vector<vector<double>> corr_density_density_boson_a_matrix(L, vector<double> (L,0));
    vector<double> p_vector(L,0);

    for( int i_site = 1; i_site <= L; i_site++ ){
        for (int j_site = 1; j_site <= L; j_site++){

          int j_site_boson_a = 2*j_site;
          int i_site_boson_a = 2*i_site;
   
          auto psi_tmp = psi;
          psi_tmp.position(i_site_boson_a);
          auto a_i = op(sites,"A",i_site_boson_a);
          auto new_a_i = a_i*psi_tmp(i_site_boson_a);
          new_a_i.noPrime();
          psi_tmp.set(i_site_boson_a,new_a_i);
          psi_tmp.noPrime();

          psi_tmp.position(j_site_boson_a);
          auto a_dagger_j = op(sites,"Adag",j_site_boson_a);
          auto new_a_dagger_j = a_dagger_j*psi_tmp(j_site_boson_a);
          new_a_dagger_j.noPrime();
          psi_tmp.set(j_site_boson_a,new_a_dagger_j);
          auto corr = innerC(psi,psi_tmp);
          corr_one_body_boson_a_matrix[i_site - 1][j_site - 1] = corr.real();
          SF_val = SF_val + corr_one_body_boson_a_matrix[i_site - 1][j_site - 1];

          // <n^a_i n^a_j> - <n^a_i><n^a_j>
          psi_tmp = psi;
          psi_tmp.position(i_site_boson_a);
          auto N_i = op(sites,"N",i_site_boson_a);
          auto new_N_i = N_i*psi_tmp(i_site_boson_a);
          new_N_i.noPrime();
          psi_tmp.set(i_site_boson_a,new_N_i);
          psi_tmp.noPrime();

          psi_tmp.position(j_site_boson_a);
          auto N_j = op(sites,"N",j_site_boson_a);
          auto new_N_j = N_j*psi_tmp(j_site_boson_a);
          new_N_j.noPrime();
          psi_tmp.set(j_site_boson_a,new_N_j);
          auto corr_density_density_boson_a = innerC(psi,psi_tmp);

          complex<double> na_i_na_j = eltC(dag(prime(psi(j_site_boson_a), "Site")) * op(sites, "N", j_site_boson_a) * psi(j_site_boson_a)) - eltC(dag(prime(psi(i_site_boson_a), "Site")) * op(sites, "N", i_site_boson_a) * psi(i_site_boson_a));

          corr_density_density_boson_a_matrix[i_site - 1][j_site - 1] = corr_density_density_boson_a.real() - na_i_na_j.real();
        }
    }    
    SF_val = SF_val/L/N_a.real();

 

    // Calculate one-body density matrix rho^{a}_{j,i0} = <a^\dagger_j a_i0> for fist component "a"
    // where i_0 = L/2
    vector<complex<double>> correlations_one_body_boson_a(L);
    for (int j_site = 1; j_site <= L; j_site++){
      auto psi_tmp = psi;
      int j_site_boson_a = 2*j_site-1;
      psi_tmp.position(i0_site_boson_a);
      auto a_dagger_i = op(sites,"Adag",i0_site_boson_a);
      auto new_a_dagger_i = a_dagger_i*psi_tmp(i0_site_boson_a);
      new_a_dagger_i.noPrime();
      psi_tmp.set(i0_site_boson_a,new_a_dagger_i);

      psi_tmp.position(j_site_boson_a);
      auto a_j = op(sites,"A",j_site_boson_a);
      auto new_a_j = a_j*psi_tmp(j_site_boson_a);
      new_a_j.noPrime();
      psi_tmp.set(j_site_boson_a,new_a_j);
      auto corr_one_body = innerC(psi,psi_tmp);
      correlations_one_body_boson_a[j_site - 1] = corr_one_body;
    }

    // Calculate two-body density matrix rho^{aa}_{j,i0} = <a^\dagger_j a^\dagger_j a_i0 a_i0>
    // where i_0 = L/2  
    vector<complex<double>> correlations_two_body_boson_a_boson_a(L); 
    for (int j_site = 1; j_site <= L; j_site++){
      int j_site_boson_a = 2*j_site-1;

      auto psi_tmp = psi;
      auto psi_tmp2 = psi;

      psi_tmp.position(i0_site_boson_a);
      auto a_dagger_i = op(sites,"Adag",i0_site_boson_a);
      auto new_a_dagger_i = a_dagger_i*psi_tmp(i0_site_boson_a);
      new_a_dagger_i.noPrime();
      psi_tmp.set(i0_site_boson_a,new_a_dagger_i);
      psi_tmp.noPrime();

      psi_tmp.position(i0_site_boson_a);
      a_dagger_i = op(sites,"Adag",i0_site_boson_a);
      new_a_dagger_i = a_dagger_i*psi_tmp(i0_site_boson_a);
      new_a_dagger_i.noPrime();
      psi_tmp.set(i0_site_boson_a,new_a_dagger_i);    
      psi_tmp.noPrime();      
      
      psi_tmp2.position(j_site_boson_a);
      auto a_dagger_j = op(sites,"Adag",j_site_boson_a);
      auto new_a_dagger_j = a_dagger_j*psi_tmp2(j_site_boson_a);
      new_a_dagger_j.noPrime();
      psi_tmp2.set(j_site_boson_a,new_a_dagger_j);
      psi_tmp2.noPrime();

      psi_tmp2.position(j_site_boson_a);
      a_dagger_j = op(sites,"Adag",j_site_boson_a);
      new_a_dagger_j = a_dagger_j*psi_tmp2(j_site_boson_a);
      new_a_dagger_j.noPrime();
      psi_tmp2.set(j_site_boson_a,new_a_dagger_j);
      psi_tmp2.noPrime();
      correlations_two_body_boson_a_boson_a[j_site-1] = innerC(psi_tmp2,psi_tmp);
    }
    /////////////////////////////////////////////////////////////////////////////////

    // Calculate two-body density matrix rho^{ab}_{j,i0} = <a^\dagger_j b^\dagger_j b_i0 a_i0> - <a^dagger_j a_i0><b^dagger_j b_i0> = 
    //                                                   = <a^\dagger_j b^\dagger_j b_i0 a_i0> -  <a^dagger_j a_i0>^2
    // where i_0 = L/2  
    vector<complex<double>> correlations_two_body_boson_a_boson_b(L); 
    for (int j_site = 1; j_site <= L; j_site++){
      int j_site_boson_a = 2*j_site-1;
      int j_site_boson_b = 2*j_site;

      auto psi_tmp = psi;
      auto psi_tmp2 = psi;

      psi_tmp.position(i0_site_boson_a);
      auto a_dagger_i = op(sites,"Adag",i0_site_boson_a);
      auto new_a_dagger_i = a_dagger_i*psi_tmp(i0_site_boson_a);
      new_a_dagger_i.noPrime();
      psi_tmp.set(i0_site_boson_a,new_a_dagger_i);
      psi_tmp.noPrime();

      psi_tmp.position(i0_site_boson_b);
      auto b_dagger_i = op(sites,"Adag",i0_site_boson_b);
      auto new_b_dagger_i = b_dagger_i*psi_tmp(i0_site_boson_b);
      new_b_dagger_i.noPrime();
      psi_tmp.set(i0_site_boson_b,new_b_dagger_i);    
      psi_tmp.noPrime();      
      
      psi_tmp2.position(j_site_boson_a);
      auto a_dagger_j = op(sites,"Adag",j_site_boson_a);
      auto new_a_dagger_j = a_dagger_j*psi_tmp2(j_site_boson_a);
      new_a_dagger_j.noPrime();
      psi_tmp2.set(j_site_boson_a,new_a_dagger_j);
      psi_tmp2.noPrime();

      psi_tmp2.position(j_site_boson_b);
      auto b_dagger_j = op(sites,"Adag",j_site_boson_b);
      auto new_b_dagger_j = b_dagger_j*psi_tmp2(j_site_boson_b);
      new_b_dagger_j.noPrime();
      psi_tmp2.set(j_site_boson_b,new_b_dagger_j);
      psi_tmp2.noPrime();
      correlations_two_body_boson_a_boson_b[j_site-1] = innerC(psi_tmp2,psi_tmp) - pow(correlations_one_body_boson_a[j_site - 1],2);
    }
    /////////////////////////////////////////////////////////////////////////////////



    




    // Calculate structure factor S(k) and M(k) arxiv:2112.10386
    // 
    int N_k = L-1;
    double pi = 3.14159265359;
    double k_max = pi;
    double k_min = -pi;
    double dk = (k_max-k_min)/(N_k-1);
    vector<complex<double>> S_vector(N_k, 0);
    vector<complex<double>> M_vector(N_k, 0);
    for (int k_int = 0; k_int<N_k; k_int++){
        //double k = k_min + k_int*dk;
        double k = k_int*2.0*pi/L;
        for(int i = 0; i<L; i++){
            for(int j = 0; j<L; j++){
                S_vector[k_int] = S_vector[k_int] + exp(Cplx_i*k*(i-j))*corr_density_density_boson_a_matrix[i][j]/L/L;
                M_vector[k_int] = M_vector[k_int] + exp(Cplx_i*k*(i-j))*corr_one_body_boson_a_matrix[i][j]/L/L;
            }
        }
    }

  complex<double> a_i0_val = 0;
  complex<double> PSF_val = 0;
  for(int j_site = 1; j_site<=L; j_site++){
    if(j_site!=L/2){
      a_i0_val  = a_i0_val  + corr_one_body_boson_a_matrix[j_site-1][L/2];
      PSF_val   = PSF_val   + correlations_two_body_boson_a_boson_b[j_site-1];
    }
  }

  //Save densities; correlations; fourier transforms to file file
  ofstream   file_data;
  string string_filename_data = "data_" + filename_parameters;

  file_data.open(string_filename_data,ios::out);
  for(int j_site = 1; j_site<=L; j_site++){
    file_data <<j_site<<" "<<density_boson_a[j_site].real()<<" "<<std_density_boson_a[j_site].real();
    file_data <<" "<<density_boson_b[j_site].real()<<" "<<std_density_boson_b[j_site].real();
    file_data <<" "<<correlations_one_body_boson_a[j_site-1].real();
    file_data <<" "<<correlations_two_body_boson_a_boson_a[j_site-1].real();
    file_data <<" "<<correlations_two_body_boson_a_boson_b[j_site-1].real();
    file_data <<" "<<corr_density_density_boson_a_matrix[j_site - 1][L/2];
    file_data <<" "<<S_vector[j_site-1]<<" "<<M_vector[j_site-1]<<endl;
  }
  file_data.close();

 
  return make_tuple(SF_val, a_i0_val.real(), PSF_val.real());
} 


string get_string_filename(double g){
	char g_str_format[40];
  sprintf(g_str_format,"%04.4f",g);
	return g_str_format;
}


 

int main(int argc, char *argv[])
{
    auto CutBondDimension = atoi(argv[1]);
    auto sweeps_tdvp_input = atoi(argv[2]);
    auto nmultMPO_cutoff_exp = atof(argv[3]);
    auto nmultMPO_cutoff = pow(10,-1.0*nmultMPO_cutoff_exp);
    auto nmultMPO_MaxDim = atoi(argv[4]);
    auto nmultMPO_method = argv[5];
    auto epsilon_K_exp_tdvp = atof(argv[6]);
    auto add_basis_method_tdvp = argv[7];

    int  PBC = atoi(argv[8]);
    auto L      = atoi(argv[9]);              // Number of lattice sites
    auto N      = atoi(argv[10]);
    auto t      = atof(argv[11]);              // Hopping amplitude
    auto U      = atof(argv[12]);              // Hubbard interaction amplitude between a-type particles 
    auto r      = atof(argv[13]);              // Hubbard interaction amplitude betwen a- b- type particles
    auto maxOccupation = atoi(argv[14]);
    auto NoOfSteps_tdvp = atoi(argv[15]);
    auto dmrg_N_max_sweeps = atoi(argv[16]);

    double U_ab = (r-1)*U;

    int droplet_position_site = L/2;
    cout<<get_string_filename(U);

    double aaa;
    //cin>>aaa;
    cout<<endl;
    cout<<"PBC = "<<PBC<<endl;
    cout<<"L =  "<<L<<endl;
    cout<<"N = "<<N<<endl;
    cout<<"t = "<<t<<endl;
    cout<<"U = "<<U<<endl;
    cout<<"U_ab = "<<U_ab<<endl;
    cout<<"r = "<<r<<endl;

    double pass;
    int N_a = int(N/2);
    int N_b = int(N/2);

    char N_buffer [50];
    sprintf (N_buffer, "%03d",N);

    char L_buffer [50];
    sprintf (L_buffer, "%03d",L);

    char t_buffer [50];
    sprintf (t_buffer, "%05.2f", t);

    char U_buffer [50];
    sprintf (U_buffer, "%05.2f", U);
    cout<<U_buffer<<endl;

    char r_buffer [50];
    sprintf (r_buffer, "%03.2f", r);
    
   
    string string_parameters =  "_PBC." + to_string(PBC)  + "_L." + L_buffer + "_N." + N_buffer + "_U." + U_buffer + "_t." + t_buffer + "_r." + r_buffer; 

    double pi = 3.14159265359;
    
 

 
    Real cutoff = 1E-16;  			    //truncation error cutoff when restoring MPS form


    int M_sites = 2*L;
    cout<<"M_sites = "<<M_sites<<endl;
    auto sites = Boson(M_sites,{"ConserveNb",true,"MaxOcc=",int(maxOccupation),"ConserveSz",true});
    auto state = InitState(sites); 
 
    cout<<"N_a = "<<N_a<<" N_b = "<<N_b<<endl;
    cout<<"L = " << L << " sites"<<endl; 

    int index_right = int(droplet_position_site + N_a/2);
    int index_left;
    int idx_tmp = N_a/2;
    if(idx_tmp%2 == 0){
      index_right = index_right;
    }
    else{
      index_right = index_right+1;
    }
    for(int i = int(droplet_position_site - N_a/2); i<index_right; i++){
        cout<<i<<endl;
        state.set(2*i-1,"1");
        state.set(2*i,"1");

    }
    auto psi0 = MPS(state);
    double SF_val, PSF_val, a_i0_val;
    tie(SF_val, a_i0_val, PSF_val) = get_observables(psi0, sites, "density___initial_state_DropletPosition." + to_string(droplet_position_site) + ".txt");
    
    cout<<"Saved densities"<<endl;
   
 
   auto quiet = "true";
  

    auto H_total = toMPO(AutoMPO(sites));
    auto H_a_t   = toMPO(AutoMPO(sites));
    auto H_aa_U  = toMPO(AutoMPO(sites));
    auto H_b_t   = toMPO(AutoMPO(sites));
    auto H_bb_U  = toMPO(AutoMPO(sites));
    auto H_ab_U  = toMPO(AutoMPO(sites));
    auto H_potential =  toMPO(AutoMPO(sites));
 
  
    
  tie(H_total, H_a_t, H_aa_U, H_b_t, H_bb_U, H_ab_U, H_potential) = get_H(sites, PBC, t, U, U_ab);

  auto psi1 = psi0;
  int maxBondDimension = maxLinkDim(psi0);

    int nsweeps = 4;
    Real dt_bysweep = 0.02 ; //time step per sweep (smaller is generally more accurate)
    auto T = -1_i*dt_bysweep;
 
    
    Real tstep = nsweeps*dt_bysweep; // total time to evolve after after one round of sweeps

    auto sweeps_tdvp = Sweeps(nsweeps);
    sweeps_tdvp.maxdim() = CutBondDimension;
    sweeps_tdvp.cutoff() = 1E-18;
    sweeps_tdvp.niter() = 50;
    int cut0=0;
    int cut1=0;
    auto energy = 0.1;
  
    for(int n = 1; n <= NoOfSteps_tdvp ; ++n){
      printf("\n ================= \n n = ", n);
      printf("\n  time form ", (n-1)*tstep, " to ", n*tstep, "\n ================= \n");
      if(maxLinkDim(psi1)<CutBondDimension){
        cut1 = 0;
        double eps_K = pow(10, -1.0*epsilon_K_exp_tdvp);
        std::vector<Real> epsilonK = {eps_K, eps_K};
        addBasis(psi1,H_total,epsilonK,{"Cutoff",eps_K,"Method",add_basis_method_tdvp,"KrylovOrd",3,"DoNormalize",true,"Quiet",true});
        maxBondDimension = maxLinkDim(psi1);
      }else{ 
        cut0 = 1; 
      }
      // TDVP sweep
      if(cut0 == 0){
        energy = tdvp(psi1,H_total,-1_i*T,sweeps_tdvp,{"Truncate",false, "DoNormalize",true,   "Quiet",true,   "NumCenter",1});
      }else{
        if(cut1 == 0){
          energy = tdvp(psi1,H_total,-1_i*T,sweeps_tdvp,{"Truncate",true, "DoNormalize",true,   "Quiet",true,   "NumCenter",1});
          cut1 = 1;
        }else{
          energy = tdvp(psi1,H_total,-1_i*T,sweeps_tdvp,{"Truncate",false, "DoNormalize",true,   "Quiet",true,   "NumCenter",1});
        }
      }

      printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi1));
    }
    
    
    ////////////////////////////////////////////////////////////////////////////

    ofstream file_monitor_convergence_energy_and_entropy_vs_sweep;

    string filename_monitor_convergence_energy_and_entropy_vs_sweep = "monitor_covnergence_energy_and_entropy_vs_sweep" + string_parameters + ".txt";
    file_monitor_convergence_energy_and_entropy_vs_sweep.open(filename_monitor_convergence_energy_and_entropy_vs_sweep,ios::out);
     
  

    char nsweeps_buffer [50];
    char BondDim_buffer [50];
    auto psi2 = psi1;
    auto psi3 = psi1;

    nsweeps = 4 ;
    int BondDim = 100;
    for(int y = 1; y<=200; y++){
 
      if(y > 10){ BondDim = 200;}
      if(y > 20){ BondDim = 300;}
      if(y > 30){ BondDim = 400;}
      if(y > 40){ BondDim = 500;}
      if(y > 50){ BondDim = 600;}
      if(y > 60){ BondDim = 700;}
      if(y > 70){ BondDim = 800;}
      if(y > 80){ BondDim = 900;}
      if(y > 90){ BondDim = 1000;}
    
      auto sweeps_dmrg = Sweeps(nsweeps);
      sweeps_dmrg.maxdim() = BondDim;
      sweeps_dmrg.cutoff() = 1E-16;
  
      auto [energy, psi3] = dmrg(H_total, psi2, sweeps_dmrg, {"Quiet", true}); //get ground state
      psi2 = psi3;
      cout<<"y = "<<y<<endl;
      cout<<"sweep = "<<y*nsweeps<<endl;
      cout<<"BondDim = "<< BondDim<<" "<<maxLinkDim(psi2)<<endl;
      printfln("Ground State Energy/N = %.10f", energy/N);
      
      if(y%5==0){

      sprintf (nsweeps_buffer, "%04d",nsweeps*y);
      sprintf (BondDim_buffer, "%04d", BondDim);

        
      auto entanglement_entropy_at_droplet_position = get_entanglement_entropy(psi2, sites, droplet_position_site, "entanglement" + string_parameters + "_sweeps." + nsweeps_buffer + "_bondDim." + BondDim_buffer + ".txt" );
      
      double SF_val, a_i0_val, PSF_val;
      tie(SF_val, a_i0_val, PSF_val) = get_observables(psi2, sites, string_parameters + "_sweeps." + nsweeps_buffer + "_bondDim." + BondDim_buffer + ".txt");


       // writeToFile("sites" + string_parameters + "_dmrg_sweeps." + nsweeps_buffer + "_bondDim." + BondDim_buffer + ".site",sites);
       // writeToFile("psi" + string_parameters + "_dmrg_sweeps." + nsweeps_buffer + "_bondDim." + BondDim_buffer + ".mps",psi3);


  
      double E_total          = innerC(psi2, H_total   , psi2).real();
     
      double E_kin_boson_a    = innerC(psi2, H_a_t     , psi2).real();
      double E_int_boson_a    = innerC(psi2, H_aa_U    , psi2).real();

      double E_kin_boson_b    = innerC(psi2, H_b_t     , psi2).real();
      double E_int_boson_b    = innerC(psi2, H_bb_U    , psi2).real();

      double E_int_boson_a_boson_b  = innerC(psi2, H_ab_U    , psi2).real();



      string output_string = to_string(PBC) + " "  + to_string(L) + " " + to_string(N_a) + " " + to_string(N_b) + "  " + get_string(U) + " " + get_string(t) + " " + get_string(r);
      

      output_string = output_string + " " + to_string(y*nsweeps) + " " + to_string(BondDim);
      output_string = output_string + "   " + get_string(E_total) + " " + get_string(E_kin_boson_a) + " " + get_string(E_int_boson_a);
      output_string = output_string + " " + get_string(E_kin_boson_b) + " " + get_string(E_int_boson_b) + " " + get_string(E_int_boson_a_boson_b);
      output_string = output_string + "    " + get_string(entanglement_entropy_at_droplet_position) + "     " + get_string(SF_val) + " " + get_string(a_i0_val) + " " + get_string(PSF_val);

      file_monitor_convergence_energy_and_entropy_vs_sweep<<output_string<<endl;

    }
      if(dmrg_N_max_sweeps < y*nsweeps){ break;};
    }
    file_monitor_convergence_energy_and_entropy_vs_sweep.close();

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //                                       Save final results

    ofstream file_final_results_vs_parameters;
    string filename_final_results = "final_results_" + string_parameters + ".txt";
    file_final_results_vs_parameters.open(filename_final_results,ios::out);
    
    auto entanglement_entropy_at_droplet_position = get_entanglement_entropy(psi3, sites, droplet_position_site, "entanglement" + string_parameters + "_sweeps." + nsweeps_buffer + "_bondDim." + BondDim_buffer + ".txt" );  
    tie(SF_val, a_i0_val, PSF_val) = get_observables(psi3, sites, "_final_" + string_parameters + "_sweeps." + nsweeps_buffer + "_bondDim." + BondDim_buffer + ".txt");

   
    double E_total          = innerC(psi2, H_total   , psi2).real();
    
    double E_kin_boson_a    = innerC(psi2, H_a_t     , psi2).real();
    double E_int_boson_a    = innerC(psi2, H_aa_U    , psi2).real();

    double E_kin_boson_b    = innerC(psi2, H_b_t     , psi2).real();
    double E_int_boson_b    = innerC(psi2, H_bb_U    , psi2).real();

    double E_int_boson_a_boson_b  = innerC(psi2, H_ab_U    , psi2).real();



    string output_string = to_string(PBC) + " "  + to_string(L) + " " + to_string(N_a) + " " + to_string(N_b) + "  " + get_string(U) + " " + get_string(t) + " " + get_string(r);
    

    output_string = output_string + "   "  + get_string(E_total) + " " + get_string(E_kin_boson_a) + " " + get_string(E_int_boson_a);
    output_string = output_string + " "    + get_string(E_kin_boson_b) + " " + get_string(E_int_boson_b) + " " + get_string(E_int_boson_a_boson_b);
    output_string = output_string + "    " + get_string(entanglement_entropy_at_droplet_position) + "     " + get_string(SF_val) + " " + get_string(a_i0_val) + " " + get_string(PSF_val);

    file_monitor_convergence_energy_and_entropy_vs_sweep<<output_string<<endl;

    file_final_results_vs_parameters<<output_string<<endl;
    file_final_results_vs_parameters.close();

    
        


    /////////////////////////////////////////////////////////////////////////////////////////////////////
  return 0;
}
