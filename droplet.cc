#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
#include "funcs/general.h"
#include <math.h>
#include <iostream>
#include <string>
#include <stdlib.h> /* srand, rand*/
#include <time.h> /* time */

using std::vector;
using namespace itensor;
using namespace std;

int main(int argc, char** argv){

double t = 1;
double U = 1;
double U_ab = 1;
int Natoms = 10;
int L = 100;
int maxOccupation = Natoms;


int M_sites = 2*L;
cout<<"M_sites = "<<M_sites<<endl;
auto sites = Boson(M_sites,{"ConserveNb",true,"MaxOcc=",int(maxOccupation),"ConserveSz",true});
auto state = InitState(sites);

auto H_total = toMPO(AutoMPO(sites));
auto H_a_t   = toMPO(AutoMPO(sites));
auto H_aa_U  = toMPO(AutoMPO(sites));
auto H_b_t   = toMPO(AutoMPO(sites));
auto H_bb_U  = toMPO(AutoMPO(sites));
auto H_ab_U  = toMPO(AutoMPO(sites));
auto H_potential =  toMPO(AutoMPO(sites));

tie(H_total, H_a_t, H_aa_U, H_b_t, H_bb_U, H_ab_U, H_potential) = get_H(sites,t,U,U_ab);
printfln("Maximum bond dimension of H is %d",maxLinkDim(H_total));

MPS psi0 = initial_state(sites, Natoms);
printfln("Energy %d", innerC(psi0, H_total, psi0));

printfln("density_a %d", density_a(sites, psi0, 50));
printfln("density_b %d", density_b(sites, psi0, 10));

printfln("density_density_a %d", density_density_correlation_a(sites, psi0, 49, 50));
printfln("density_density_b %d", density_density_correlation_b(sites, psi0, 51, 55));

printfln("one_body_corr_a %d", one_body_correlation_a(sites, psi0, 49, 50));
printfln("one_body_corr_b %d", one_body_correlation_b(sites, psi0, 51, 55));

printfln("pair_corr_ab %d", pair_correlation_ab(sites, psi0, 51, 55));








return 0;
}
