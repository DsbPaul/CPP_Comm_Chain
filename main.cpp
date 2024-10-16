/** Include zone */
#include <iomanip>
#include <iostream>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <random>
#include <string>
#include <algorithm>

/** Constant zone */
#define M_PI       3.14159265358979323846
#define N_BITS     1024
#define MOD_ORD    4
#define N_SYM      N_BITS/log2(MOD_ORD)
#define N_UPS      8
#define N_FILT     81
#define ALPHA      0.25
#define MC         10
#define N_FRAMES   10
#define CONST_REF  "QPSK"
#define L_EXC      3

/** Namespace */
using namespace std;

/** Functions declaration */
vector<int> bit_gen(int nb_bits);
vector<int> mapping_dec(vector<int> bits, int mod_order);
vector<complex<double>> symbol_map(vector<int> symb_dec, string const_ref, int mod_order);
double sig_power(vector<complex<double>> sig_in, int surech);
vector<complex<double>> upsampling(vector<complex<double>> signal, int ups_fact);
vector<complex<double>> awgn_chan(vector<complex<double>> signal, double Eb_N0, int ups_fact);
vector<complex<double>> srrc(double alpha, int Nech, int N, double Ts);
vector<complex<double>> filtering(vector<complex<double>> filtre,vector<complex<double>> signal);
vector<complex<double>> sampling(vector<complex<double>> sig_in, int ups_fact);
vector<complex<double>> normalization(vector<complex<double>> sig_in);
vector<int> ML_detector(vector<complex<double>> sig_in,vector<complex<double>> const_ref);
vector<int> symbol_demap(vector<int> symb_dec, string const_ref, int mod_order);
vector<complex<double>> constellation(string const_ref, int mod_order);
int error_calc(vector<int> bits_tx, vector<int> bits_rx, int l_exc);

/** Main code */
int main() {
    // Parameters
    vector<double> ber, ber_th;
    vector<int> bits_tx, bits_rx, symb_dec_tx, symb_dec_rx, num_err;
    vector<complex<double>> symboles_tx, symb_ups, sig_tx, sig_awgn, sig_rx, symb_down, const_ref;
    vector<complex<double>> filter = srrc(ALPHA,N_UPS,floor(N_FILT/2),1.0);
    int Eb_N0_dB[11] = {0,1,2,3,4,5,6,7,8,9,10};

    // SNR Loop
    for (int kk=0;kk<11;kk++){
        // SNR parameters
        double SNR = Eb_N0_dB[kk];
        double Eb_N0_lin = pow(10,0.1*SNR);
        int mc = 0, err_temp = 0;

        // Monte-Carlo
        while (((err_temp<100) || (mc<N_FRAMES)) && (mc<MC)){

            // Bits generation
            bits_tx = bit_gen(N_BITS);

            // Decimal mapping
            symb_dec_tx = mapping_dec(bits_tx,MOD_ORD);

            // Reference contellation
            const_ref = constellation(CONST_REF,MOD_ORD);

            // Symbol mapping
            symboles_tx = symbol_map(symb_dec_tx,CONST_REF,MOD_ORD);

            // Upsampling
            symb_ups = upsampling(symboles_tx,N_UPS);

            // Tx filtering
            sig_tx = filtering(filter,symb_ups);

            // AWGN channel
            sig_awgn = awgn_chan(sig_tx,Eb_N0_lin,N_UPS);
            //sig_awgn = sig_tx;

            // Rx filtering
            sig_rx = filtering(filter,sig_awgn);

            // Normalization
            sig_rx = normalization(sig_rx);

            // Sampling
            symb_down = sampling(sig_rx,N_UPS);

            // ML detection
            symb_dec_rx = ML_detector(symb_down,const_ref);

            // Binary demapping
            bits_rx = symbol_demap(symb_dec_rx,CONST_REF,MOD_ORD);

            // Error calculation
            err_temp += error_calc(bits_tx,bits_rx,L_EXC);

            // Clear vectors
            bits_tx.clear();
            bits_rx.clear();
            symb_dec_tx.clear();
            symb_dec_rx.clear();
            symboles_tx.clear();
            symb_ups.clear();
            sig_tx.clear();
            sig_awgn.clear();
            sig_rx.clear();
            symb_down.clear();

            // Monte-Carlo
            mc += 1;
        }

        // BER calculation
        num_err.push_back(err_temp);
        ber.push_back((double)num_err.at(kk) / ((N_BITS-6)*mc));
        ber_th.push_back(0.5*erfc(sqrt(Eb_N0_lin)));

        // Clear variables
        err_temp = 0;
        mc = 0;

        // Print results
        cout << "Eb/N0 = " << Eb_N0_dB[kk] << " dB = " << Eb_N0_lin << " W/Hz : " << endl;
        cout << "Number of error bits = " << num_err.at(kk) << ", ";
        cout << "Simulated BER = " << ber.at(kk) << ", ";
        cout << "Theoretical BER = " << ber_th.at(kk) << "." << endl;
        cout << endl;
    }

    return 0;
}

/** Function implementations */
vector<int> bit_gen(int nb_bits){
    vector<int> bits;
    for (int i=0;i<nb_bits;i++){
        bits.push_back(rand()%2);
    }
    return bits;
}

vector<int> mapping_dec(vector<int> bits, int mod_order){
    int symb;
    vector<int> symb_dec_tx;
    for (int i=0;i<(int)bits.size();i+=log2(mod_order)){
        symb = 0;
        for (int j=0;j<log2(mod_order);j++){
            symb += bits.at(i+j) * int(pow(2,log2(mod_order)-(j+1)));
        }
        symb_dec_tx.push_back(symb);
    }
    return symb_dec_tx;
}

vector<complex<double>> symbol_map(vector<int> symb_dec, string const_ref, int mod_order){
    vector<complex<double>> symbols;
    if (const_ref.compare("QPSK") == 0){
        for (int i=0;i<(int)symb_dec.size();i++){
            if (symb_dec.at(i) == 0){
                symbols.push_back((1.0 / sqrt(2.0))*complex<double>(1.0, 1.0));
            } else if(symb_dec.at(i) == 1){
                symbols.push_back((1.0 / sqrt(2.0))*complex<double>(1.0, -1.0));
            } else if (symb_dec.at(i) == 2){
                symbols.push_back((1.0 / sqrt(2.0))*complex<double>(-1.0, 1.0));
            } else if (symb_dec.at(i) == 3){
                symbols.push_back((1.0 / sqrt(2.0))*complex<double>(-1.0, -1.0));
            }
        }
    } else {
        symbols = vector<complex<double>>{(0.0)};
    }
    return symbols;
}

double sig_power(vector<complex<double>> sig_in, int surech){
    vector<double> sig_in_2;
    for (int i=0;i<(int)sig_in.size();i++){
        sig_in_2.push_back(pow(abs(sig_in.at(i)),2));
    }
    double mean_sig = accumulate(sig_in_2.begin(),sig_in_2.end(),0.0)/(sig_in_2.size());
    return mean_sig * surech;
}

vector<complex<double>> upsampling(vector<complex<double>> signal, int ups_fact){
    vector<complex<double>> symb_ups;
    for (int i=0;i<N_SYM;i++){
        symb_ups.push_back(signal.at(i));
        for (int j=0;j<(ups_fact-1);j++){
            symb_ups.push_back(complex<double>(0.0, 0.0));
        }
    }
    return symb_ups;
}

vector<complex<double>> awgn_chan(vector<complex<double>> signal, double Eb_N0, int ups_fact){
    random_device rd{};
    mt19937 gen{rd()};
    vector<complex<double>> noise, sig_awgn;
    double p_sig = sig_power(signal,ups_fact);
    double N0 = p_sig / (log2(MOD_ORD)*Eb_N0);
    normal_distribution<double> d{0.0, sqrt(N0/2)};
    auto random_sample = [&d, &gen]{return d(gen);};
    for (int i=0;i<(int)signal.size();i++){
        noise.push_back(complex<double>(random_sample(),random_sample()));
        sig_awgn.push_back(signal.at(i) + noise.at(i));
    }
    return sig_awgn;
}

vector<complex<double>> srrc(double alpha, int Nech, int N, double Ts){
    double Te = Ts/Nech, MaxVal = 1+alpha*(4/M_PI-1), t, g, g_lim;
    vector<complex<double>> temps,filter_coef;
    for (int i=0;i<N*Ts;i++){
        t = (i+1)*Te;
        g = (((4*alpha*t)/Ts)*cos((1+alpha)*M_PI*t/Ts)+sin(M_PI*(1-alpha)*t/Ts))/((M_PI*t/Ts)*(1-pow((4*alpha*t/Ts),2)));
        g_lim = alpha/sqrt(2)*((1-2/M_PI)*cos(M_PI/(4*alpha))+(1+2/M_PI)*sin(M_PI/(4*alpha)));
        if (i == round(Nech/(4*alpha)-1)){
            temps.push_back(g_lim);
        } else {
            temps.push_back(g);
        }
    }
    for (int i=0;i<N*Ts;i++){
        filter_coef.push_back(temps.at(N*Ts-(i+1)));
    }
    filter_coef.push_back(MaxVal);
    for (int i=0;i<N*Ts;i++){
        filter_coef.push_back(temps.at(i));
    }
    return filter_coef;
}

vector<complex<double>> filtering(vector<complex<double>> filtre,vector<complex<double>> signal){
    vector<complex<double>> result;
    vector<complex<double>> result_mod;
    int taille_f = filtre.size();
    int n = filtre.size() + signal.size() - 1 ;
    for (int i = 0 ; i < n ; i++){
        double tmp_i = 0.0;
        double tmp_q = 0.0;
        int j = 0;
        for (int i1 = 0 ; i1 < (int)signal.size() ; i1++){
            j = i - i1;
            if ( (j >= 0) && (j < (int)filtre.size()) ){
                tmp_i += (signal.at(i1).real() * filtre.at(j).real() - signal.at(i1).imag() * filtre.at(j).imag());
                tmp_q += (signal.at(i1).real() * filtre.at(j).imag() + signal.at(i1).imag() * filtre.at(j).real());
            }
        }
        result.push_back(complex<double>(tmp_i,tmp_q));
    }
    result_mod.insert(result_mod.begin(),result.begin() + floor(taille_f / 2), result.end() - floor(taille_f / 2));
	return result_mod;
}

vector<complex<double>> sampling(vector<complex<double>> sig_in, int ups_fact){
    vector<complex<double>> sig_out;
    for (int i=0;i<(int)sig_in.size();i+=ups_fact){
        sig_out.push_back(sig_in.at(i));
    }
    return sig_out;
}

vector<complex<double>> normalization(vector<complex<double>> sig_in){
    vector<complex<double> > sig_out;
	double sum = 0, fact = 0, square_symb;
	for (int i=0;i<(int)sig_in.size();i++){
		square_symb = pow(sig_in.at(i).real(),2.0) + pow(sig_in.at(i).imag(),2.0);
		sum += square_symb;
	}
	fact = sqrt(sum/sig_in.size());
	for (int i=0;i<(int)sig_in.size();i++){
        sig_out.push_back(std::complex<double>(sig_in.at(i).real()/fact,sig_in.at(i).imag()/fact));
	}
    return sig_out;
}

vector<complex<double>> constellation(string const_ref, int mod_order){
    vector<complex<double>> constellation;
    if (const_ref.compare("QPSK") == 0){
        for (int i=0;i<mod_order;i++){
            switch (i) {
                case 0:
                    constellation.push_back((1.0 / sqrt(2.0))*complex<double>(1.0, 1.0));
                    break;
                case 1:
                    constellation.push_back((1.0 / sqrt(2.0))*complex<double>(1.0, -1.0));
                    break;
                case 2:
                    constellation.push_back((1.0 / sqrt(2.0))*complex<double>(-1.0, 1.0));
                    break;
                case 3:
                    constellation.push_back((1.0 / sqrt(2.0))*complex<double>(-1.0, -1.0));
                    break;
            }
        }
    } else {
        constellation = vector<complex<double>>{(0.0)};
    }
    return constellation;
}

vector<int> ML_detector(vector<complex<double>> sig_in,vector<complex<double>> const_ref){
    vector<complex<double>> sig_out;
    vector<double> dist;
    complex<double> dist_temp;
    vector<int> dist_dec;
    for (int i=0;i<(int)sig_in.size();i++){
        for (int j=0;j<(int)const_ref.size();j++){
            dist.push_back(pow(abs(const_ref.at(j)-sig_in.at(i)),2));
        }
        dist_temp = *min_element(dist.begin(),dist.end());
        for (int j=0;j<(int)const_ref.size();j++){
            if (pow(abs(const_ref.at(j)-sig_in.at(i)),2) == dist_temp){
                dist_dec.push_back(j);
            }
        }
        dist.clear();
    }
    return dist_dec;
}

vector<int> symbol_demap(vector<int> symb_dec, string const_ref, int mod_order){
    vector<int> bits_out;
    if (const_ref.compare("QPSK") == 0){
        for (int i=0;i<(int)symb_dec.size();i++){
            if (symb_dec.at(i) == 0){
                bits_out.push_back(0);
                bits_out.push_back(0);
            } else if (symb_dec.at(i) == 1) {
                bits_out.push_back(0);
                bits_out.push_back(1);
            } else if (symb_dec.at(i) == 2) {
                bits_out.push_back(1);
                bits_out.push_back(0);
            } else if (symb_dec.at(i) == 3) {
                bits_out.push_back(1);
                bits_out.push_back(1);
            }
        }
    } else {
        bits_out = vector<int>{0};
    }
    return bits_out;
}

int error_calc(vector<int> bits_tx, vector<int> bits_rx, int l_exc){
    int err = 0;
    for (int i=l_exc;i<(int)bits_tx.size()-l_exc;i++){
        err += abs(bits_rx.at(i)-bits_tx.at(i));
    }
    return err;
}
