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

double t=1;
double U=1;
double U_ab=1;
int L=100;
int maxOccupation=1;


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

      



return 0;
}
