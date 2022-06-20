#include <cmath>
#include <fstream>
#include <iostream>
#include <complex>
#include <RandomNumber.hpp>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace std;
using namespace Snu::Cnrc;
const int N = 6;

typedef complex<double> dcomplex;
using Size = unsigned int;

/* site direction: <-
 Ex. |pis>=|10> , i=3
 for j=0, jNext=1, spin_j=0, spin_jNext=1
 for S+_{j}S-_{j+1}, i_flip=3+pow(3,j)-pow(3,j+1)=1=|01>
*/


class QSIR_1D {
private:
    Size endTime;
    double omega;
    double kappa;
    Size number_of_ensemble;
 
public:
    QSIR_1D(const Size, const double, const double, const Size);
    ~QSIR_1D();
    void QJMC();
    int dec2tri(int decimal_number, int digit_num);
    void normalize(vector<dcomplex>& state);
    vector<dcomplex> vectorCal(vector<dcomplex>& state, vector<dcomplex>& k, const double dt, const double num);
    double calDecayDeltaP(int site, const double dt, vector<dcomplex>& state, const double kappa);
    vector<dcomplex> emissionDecay(int site, vector<dcomplex>& state);
    vector<dcomplex> dcdt(vector<dcomplex> state, const double omega, const double kappa);
};

QSIR_1D::QSIR_1D(const Size _endTime_, const double _omega_, const double _kappa_, const Size _number_of_ensemble_)
    : endTime(_endTime_), omega(_omega_), kappa(_kappa_), number_of_ensemble(_number_of_ensemble_) {
}

QSIR_1D::~QSIR_1D() {
}


void QSIR_1D::QJMC() {
    vector<dcomplex> state(pow(3, N));
    vector<dcomplex> k1(pow(3, N));
    vector<dcomplex> k2(pow(3, N));
    vector<dcomplex> k3(pow(3, N));
    vector<dcomplex> k4(pow(3, N));
    double dt = 0.01;
    const int number_of_time_step = endTime/dt;
    
    RandomRealGenerator rnd(0.0, 1.0);
    vector<double> orderParameter(number_of_time_step);
    vector<double> hist_orderparameter(N+1);
    vector<double> deltaP(N);
    vector<vector<double>> orderParameter_state(number_of_time_step, vector<double>(3, 0));
    int siteEvolved;
    double check_infection_density = 0;
    for(int ensembleNum=0; ensembleNum<number_of_ensemble; ++ensembleNum) {
        for(int i=0; i<pow(3,N); ++i) {
            state[i] = {0, 0};
        }
        state[pow(3, N/3)] = {1, 0};
        normalize(state);

        for(int t=0; t<number_of_time_step; ++t) {
            for(int site=0; site<N; ++site) deltaP[site] = calDecayDeltaP(site, dt, state, kappa);
            partial_sum(deltaP.begin(), deltaP.end(), deltaP.begin());
            double p = rnd();
            siteEvolved = distance(deltaP.begin(), lower_bound(deltaP.begin(), deltaP.end(), p));
            if(siteEvolved != N) {
                state = emissionDecay(siteEvolved, state);
            }
            else { // no emission takes place
                k1 = dcdt(state, omega, kappa);
                k2 = dcdt(vectorCal(state, k1, dt, 2), omega, kappa);
                k3 = dcdt(vectorCal(state, k2, dt, 2), omega, kappa);
                k4 = dcdt(vectorCal(state, k3, dt, 1), omega, kappa);
                for(int i=0; i<pow(3, N); ++i) {
                    state[i] = (state[i] + dt/6.*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i]));
                }
                normalize(state);
            }
            
            for(int i=0; i<pow(3, N); ++i) {
                for(int j=0; j<N; j++){
                    int spin1 = dec2tri(i,j);
                    orderParameter_state[t][spin1] += norm(state[i]);
                    if(spin1 == 1){
                        check_infection_density += norm(state[i]);
                    }
                }
            }
            if(check_infection_density/N < 1e-8) {
                vector<double> steady_value(3);
                for(int j=0; j<N; j++){
                    for(int i=0; i<pow(3, N); ++i) {
                        int spin1 = dec2tri(i,j);
                        steady_value[spin1] += norm(state[i]);
                    }
                }
                for(int j_time=t+1; j_time<number_of_time_step; j_time++){
                    for(int spinState=0; spinState<3; ++spinState) {
                        orderParameter_state[j_time][spinState] += steady_value[spinState];
                    }
                }
                break;
            }
            check_infection_density = 0;
        }
    }
    
    for(int t=0; t<number_of_time_step; ++t) {
        cout << t*dt << '\t';
        for(int spinState=0; spinState<3; ++spinState) {
            cout << orderParameter_state[t][spinState]/number_of_ensemble/N << '\t';
        }
        cout << endl;
    }
}
    

int QSIR_1D::dec2tri(int decimal_number, int digit_num) {
	return (decimal_number / (int) pow(3, digit_num)) % 3;
}

void QSIR_1D::normalize(vector<dcomplex>& state) {
    double normalizationFactor = 0;
    for(Size i=0; i<state.size(); ++i)
        normalizationFactor += norm(state[i]);
    for(Size i=0; i<state.size(); ++i)
        state[i] /= sqrt(normalizationFactor);
}

vector<dcomplex> QSIR_1D::vectorCal(vector<dcomplex>& state, vector<dcomplex>& k, const double dt, const double num) {
    vector<dcomplex> cc(state.size());
    for(Size i=0; i<state.size(); ++i)
        cc[i] = state[i]+k[i]*dt/num;
    return cc;
}

// Decay Lindblad Operator

double QSIR_1D::calDecayDeltaP(int site, const double dt, vector<dcomplex>& state, const double kappa) {
    double deltaP = 0;
    for(int i=0; i<pow(3, N); ++i) {
        int spin = dec2tri(i, site);
        if(spin == 1) {
            deltaP += kappa*dt*(norm(state[i]));
        }
    }
    return deltaP;
}

vector<dcomplex> QSIR_1D::emissionDecay(int site, vector<dcomplex>& state) {
    vector<dcomplex> cc(pow(3, N));
    for(int i=0; i<pow(3, N); ++i) {
        int spin = dec2tri(i, site);
        if(spin == 1) {
            int i_flip = i + pow(3, site);
            cc[i_flip] = state[i];
        }
    }
    normalize(cc);
    return cc;
}

// Runge-Kutta method for the differential value of wavefunction
vector<dcomplex> QSIR_1D::dcdt(vector<dcomplex> state, const double omega, const double kappa) {
    dcomplex I(0, 1);
    vector<dcomplex> cc(pow(3,N));
    for(int j=0; j<N; ++j) {
        for(int i=0; i<pow(3, N); ++i) {
            int jNext = j == N-1 ? 0 : j+1;
            int spin_j = dec2tri(i, j);
            int spin_jNext = dec2tri(i, jNext);

            if(j != N-1){
                if(spin_j == 1 && spin_jNext != 2) { // omega * n_{j}sigx_{j+1 or jNext}
                    int i_flip = i + (1-2*spin_jNext) * pow(3, jNext);
                    cc[i_flip] -= I * omega * state[i];
                }
                if(spin_j != 2 && spin_jNext == 1) { // J * S-_{j}S+_{j+1}
                    int i_flip = i + (1-2*spin_j) * pow(3, j);
                    cc[i_flip] -= I * omega * state[i];                
                } 
            }

            if(spin_j == 1) {
                cc[i] -= kappa/2.*state[i]; // decay: kappa * (Sz)^2_{j}
            }
        }
    }
    return cc;
}

int main(int argc, char *argv[]) {
    const double endTime = stoul(argv[1]);
    const double omega = stod(argv[2]);
    const double number_of_ensemble = stoul(argv[3]);
    const double kappa = 1;

    std::ios_base::sync_with_stdio(false);
    cin.tie(nullptr); cout.tie(nullptr);

    QSIR_1D* Model = new QSIR_1D(endTime, omega, kappa, number_of_ensemble);
    
    Model -> QJMC();
    delete Model;
}
