//2022-01-16
//One PaO2 <-> SpO2-file to rule them all.

#include <Rcpp.h>
//[[Rcpp::plugins("cpp14")]]
//[[Rcpp::depends(BH)]]
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
//#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>
#include <thread>

//based on Dash et Al, Eur J Appl Physiol 2016 : 
//Simple Accurate Mathematical Models of Blood HbO and HbCO2 Dissociation Curves at Varied Physiological Conditions - Evaluation and Comparison with other Models
//note, this class is stripped down (a lot) compared to previous versions, and standard values are now hard-coded to match Dash's original
class so2f_dash_f{
private:
    //standard values
    double p50_s = 26.8;
    double ph_s = 7.4;
    double temp_s = 37;
    double co2_s = 40;
    double dpg23_s = 4.65e-3;
    //the current p50
    double p50 = 26.8;
    //the values for the nH-transformation
    double alpha = 2.82;
    double beta = 1.20;
    double gamma = 29.25;
    
    //the hill-factor as function of pao2
    inline double nH(double po2){
        return alpha - beta * std::pow(10, -1 * (po2 / gamma));
    }
    

    //calculate the p50
    inline double calc_p50(double pH = 7.4, double Temp = 37, double CO2 = 40, double DPG23 = 4.65e-3){
        double p50_dph = p50_s - 25.535 * (pH - ph_s) + 10.646 * std::pow(pH - ph_s, 2.0) - 1.764 * std::pow(pH - ph_s, 3.0);
        double p50_dco2 = p50_s + 1.273e-1 * (CO2 - co2_s) + 1.083e-4 * std::pow(CO2 - co2_s, 2.0);
        double p50_ddpg23 = p50_s + 795.63 * (DPG23 - dpg23_s) + 19660.89 * std::pow(DPG23 - dpg23_s, 2.0);
        double p50_dt = p50_s + 1.435 * (Temp - temp_s) + 4.163e-2 * std::pow(Temp-temp_s, 2.0) + 6.86e-4 * std::pow(Temp - temp_s, 3.0);
        //Rcpp::Rcout << "pH = " << pH << " Temp = " << Temp << " pCO2 = " << CO2 << " p50 = " << p50_s * (p50_dph / p50_s) * (p50_dco2 / p50_s) * (p50_ddpg23 / p50_s) * (p50_dt / p50_s) <<"\n";
        return p50_s * (p50_dph / p50_s) * (p50_dco2 / p50_s) * (p50_ddpg23 / p50_s) * (p50_dt / p50_s);
    }

public:
    so2f_dash_f(){}
    ~so2f_dash_f(){}
    
    //method for setting the p50...
    inline void calc_and_set_p50(double pH = 7.4, double Temp = 37, double CO2 = 40, double DPG23 = 4.65e-3){
        p50 = calc_p50(pH, Temp, CO2, DPG23);
    }
    
    //the method to calculate spo2 from po2
    inline double po2_to_s(const double pao2){
        double nh = nH(pao2);
        double p = std::pow(pao2/p50, nh);
        return p/(1.0+p);
    }

    //the operator() if we want to use this as a functor
    inline double operator()(const double pao2){
        return po2_to_s(pao2);
    }
};

//from Scan J Clin Lab Invest. 27:239-245, 1971.
class BE_calc{
    private:
        double pCO2;

    public: 
        BE_calc(){}
        ~BE_calc(){}
        void set_pCO2(double pco2){
            pCO2 = pco2;
        }

        inline double operator()(const double pH){
            //B.E. = 0.02786 * pCO2 * 10 (pH - 6.1) + 13.77 * pH - 124.58
            return 0.02786 * 7.5 * pCO2 * std::pow(10, pH-6.1) + 13.77 * pH - 124.58;
        }
};

//from Scan J Clin Lab Invest. 27:239-245, 1971. but solved for pCO2
class pCO2_calc{
    private:
        double BE;

    public: 
        pCO2_calc(){}
        ~pCO2_calc(){}
        void set_BE(double be){
            BE = be;
        }

        inline double operator()(const double pH){
            //B.E. = 0.02786 * pCO2 * 10 (pH - 6.1) + 13.77 * pH - 124.58
            //pCO2 = (B.E. - 13.77 * pH + 124.58) / (0.02786 * 10^(pH - 6.1))
            

            return ((BE - 13.77 * pH + 124.58) / (0.02786 * std::pow(10, pH-6.1)));
        }
};






//from Kelman 1966
//J Appl Physiol. 1966 Jul;21(4):1375-6. doi: 10.1152/jappl.1966.21.4.1375.
//Digital computer subroutine for the conversion of oxygen tension into saturation
inline double so2f_kelman_f(const double po2){
  static std::array<double,7> k{{-8.53229e3, 2.121401e3, -6.707399e1, 9.359609e5,-3.134626e4,2.396167e3, -6.710441e1}};
  return 1.0 * (po2*(po2*(po2*(po2+k[2])+k[1])+k[0])) / (po2*(po2*(po2*(po2+k[6])+k[5])+k[4])+k[3]);
}
//for use from R
//[[Rcpp::export]]
Rcpp::NumericVector cpp_kelman(Rcpp::NumericVector po2){
  Rcpp::NumericVector out(po2.size());
  std::transform(po2.begin(), po2.end(), out.begin(), so2f_kelman_f);
  return out;
}


//From Aberman 1973
//J Appl Physiol. 1973 Oct;35(4):570-1. doi: 10.1152/jappl.1973.35.4.570.
//An equation for the oxygen hemoglobin dissociation curve
inline double so2f_aberman_f(const double po2){
  static std::array<double,8> k  {{51.87074,129.8325,6.828368,-223.7881,-27.9530,258.5009,21.84175,-119.2322}};
  double sum = 0.0;
  double y = (po2-27.5) / (po2 + 27.5);
  for (auto i = 0u;i<8;++i){
      sum+=k[i]*std::pow(y, i);
  }
  return sum*0.01;
}
//for use from R
//[[Rcpp::export]]
Rcpp::NumericVector cpp_aberman(Rcpp::NumericVector po2){
  Rcpp::NumericVector out(po2.size());
  std::transform(po2.begin(), po2.end(), out.begin(), so2f_aberman_f);
  return out;
}


//From severinghaus 1979
//J Appl Physiol Respir Environ Exerc Physiol. 1979 Mar;46(3):599-602. doi: 10.1152/jappl.1979.46.3.599.
//Simple, accurate equations for human blood O2 dissociation computations
inline double so2f_severinghaus_f(const double po2){
    return std::pow(std::pow(std::pow(po2, 3.0) + 150*po2, -1.0)*23400 + 1.0,-1.0);
}

//for use from R
//[[Rcpp::export]]
Rcpp::NumericVector cpp_severinghaus(Rcpp::NumericVector po2){
  Rcpp::NumericVector out(po2.size());
  std::transform(po2.begin(), po2.end(), out.begin(), so2f_severinghaus_f);
  return out;
}

//Ellis' inversion of the Severinghaus equation, from 1989.
//This will be used instead of the numerical approximation when method is 'severinghaus'
//Letter to the editor, J. Appl. Physiol 189;67;902
inline double po2_ellis(double so2){
    double fact = 11700.0/(1.0/so2-1.0);
    double f2 = std::pow(125000.0 + std::pow(fact, 2.0),0.5);
    return std::pow(fact+ f2, 1.0/3)-std::pow(-(fact - f2), 1.0/3);
}

//for use from R
//[[Rcpp::export]]
Rcpp::NumericVector cpp_ellis(Rcpp::NumericVector so2){
    //create the out vector
    Rcpp::NumericVector out = Rcpp::clone(so2);
    auto n = so2.size();
    double fact, f2;
    for (auto i = 0u;i<n;++i){
        fact = 11700.0/(1.0/so2[i]-1.0);
        f2 = std::pow(125000.0 + std::pow(fact, 2.0),0.5);
        out[i] = std::pow(fact+ f2, 1.0/3)-std::pow(-(fact - f2), 1.0/3);
    }
    return out;
}

//for use from R
//[[Rcpp::export]]
Rcpp::NumericVector cpp_ellis_kpa(Rcpp::NumericVector so2){
    //create the out vector
    auto n = so2.size();
    Rcpp::NumericVector out(n);
    
    double fact, f2;
    double conversion_factor = 1.0/7.50062;
    for (auto i = 0u;i<n;++i){
        fact = 11700.0/(1.0/(so2[i]*0.01)-1.0);
        f2 = std::pow(125000.0 + std::pow(fact, 2.0),0.5);
        out[i] = conversion_factor * (std::pow(fact+ f2, 1.0/3)-std::pow(-(fact - f2), 1.0/3));
    }
    return out;
}

//a quick-and-dirty multithreaded version of Ellis equation
//[[Rcpp::export]]
Rcpp::NumericVector cpp_ellis_thread(Rcpp::NumericVector so2, int nThreads = 0){
    //calculate maximum possible number of threads
    int max_threads = static_cast<int>(std::thread::hardware_concurrency());
	int n = so2.size();
    if (nThreads > max_threads || nThreads <= 0){
        nThreads = max_threads;
    }

    
    //keep the number of threads from being excessively large
    nThreads = std::min(n/1000 + 1,nThreads);

    //reserve space in threads vector
    std::vector<std::thread> threads;
    threads.reserve(nThreads);

    //set up the 'out'-vector
    Rcpp::NumericVector out(n);
    //some math to calculate the number of entries per thread
    int div = n / nThreads;
    int mod = n % nThreads;
    int begin = 0;
    int end = 0;
    for (auto i = 0;i<nThreads;++i){
        //the first n % nThreads threads will take care of one extra entry in the vector
        end += div + (mod-- <= 0 ? 0 : 1);
        //this lambda contains ellis' equation
        threads.emplace_back([&so2, &out, begin, end]()->void{
            double fact, f2;
            //since the output will be in kPa, not mmHg
            const double conversion_factor = 1.0/7.50062;
            for (auto j = begin;j<end;++j){
                fact = 11700.0/(1.0/(so2[j]*0.01)-1.0);
                f2 = std::pow(125000.0 + std::pow(fact, 2.0),0.5);
                out[j] = conversion_factor * (std::pow(fact + f2, 1.0/3)-std::pow(-(fact - f2), 1.0/3));
            }
        });
        begin = end;
    }

    for (auto &t: threads){
        t.join();
    }
    return out;
}

//a class for shifting po2s back and forth based on ph and temp and so on
class po2shift{
public:
    po2shift(){};
    //correcting a measured po2 into standard po2
    double po2c(const double po2, const double temp = 37.0, const double ph = 7.4, const double be = 0.0){
        return po2 *std::pow(10, (0.024*(37.0-temp) - 0.48*(7.4-ph) - 0.0013*be));
    }
    
    //2022-01-20
    //just trying out my algebra...
    //the inverse of the above function, shifting po2 back from 'standard' to expected po2 at the current ph, temp, be
    double po2e(const double po2, const double temp = 37.0, const double ph = 7.4, const double be = 0.0){
        return po2 / std::pow(10, (0.024*(37.0-temp) - 0.48*(7.4-ph) - 0.0013*be));
    }
    
    //From severinghaus
    inline double po2t(double po2s, double temp){
        //()the formula:
        //(log(po2) - log(po2s)) / (temp - 37) = 0.058 * (1/(0.243 * (0.01*po2s)^3.88+1.0))+0.013
        //right hand side of the equation
        auto r = 0.058 * (1/(0.243 * std::pow(0.01*po2s, 3.88)+1.0))+0.013;
        //multiply by temp - 37.0 and add log po2s
        auto s = (temp - 37.0) * r + std::log(po2s);
        return std::exp(s);
    }
    //from severinghaus
    inline double po2ph(double po2s, double ph){
        //the formula
        //log(po2) - log(po2s) / (ph - 7.4) = (po2s / 26.6)^0.184 - 2.2
        //right hand side of the equation
        auto r = std::pow(po2s / 26.6, 0.184) - 2.2;
        //multiply by ph - 7.4 and add log po2s
        auto s = r*(ph-7.4) + std::log(po2s);
        return std::exp(s);
        //return po2s * std::exp(r*(ph-7.4));
    }
    //from severinghaus
    inline double po2p50(double po2s, double p50){
        return (p50/26.6) * po2s;
    }
    
    //From Severinghaus, correcting a calculated standard po2s for temperature and ph into acutal po2
    double po2s(const double po2, const double temp = 37.0, const double ph = 7.4, const double p50  = 26.6){
        return po2p50(po2ph(po2t(po2, temp), ph), p50);      
    }
};

//The main function, for calculating spo2 from pao2 using any of the 4 methods.
//This also features some convenience functions such as picking units (even though the program internally uses
//Torr for pressure and
//fraction (0-1) for saturation)
//[[Rcpp::export]]
Rcpp::NumericVector pao2_to_spo2(Rcpp::NumericVector po2, 
                                Rcpp::Nullable<Rcpp::NumericVector> Temp = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> pH = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> BE = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> pCO2 = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> DPG23 = R_NilValue,
                                Rcpp::String method = "kelman", 
                                Rcpp::String unit_S = "fraction", 
                                Rcpp::String unit_P = "Torr",
                                Rcpp::String na_handling = "impute",
                                double default_temp = 37.0,
                                double default_ph = 7.4,
                                double default_be = 0,
                                double default_co2 = 40,
                                double default_dpg23 = 4.65e-3,
                                double default_p50 = 26.8,
                                bool correct = false, 
                                bool verbose = false){
    //set up conversion factors and settings for the loop later on
    double s_conversion_factor = 1.0;
    double p_conversion_factor = 1.0;
    bool shift_p50 = false;
    double n = po2.size();
    int na_method = 0;
    
    //create default vectors
    Rcpp::NumericVector temp = Rcpp::NumericVector::create();
    Rcpp::NumericVector ph = Rcpp::NumericVector::create();
    Rcpp::NumericVector be = Rcpp::NumericVector::create();
    Rcpp::NumericVector co2 = Rcpp::NumericVector::create();
    Rcpp::NumericVector dpg23 = Rcpp::NumericVector::create();
    if (Temp.isNotNull()){
        temp = Temp;
    }
    
    if (pH.isNotNull()){
        ph = pH;
    }
    
    if (BE.isNotNull()){
        be = BE;
    }
    
    if (pCO2.isNotNull()){
        co2 = pCO2;
    }
    
    if (DPG23.isNotNull()){
        dpg23 = DPG23;
    }
    
    //setting input unit
    if (unit_P == "kPa")
    {
        p_conversion_factor = 7.50062;
        if (verbose) Rcpp::Rcout << "p_conversion_factor: " << p_conversion_factor << "\n";
    } 
    else if (unit_P == "mmHg")
    {
        p_conversion_factor = 1.000000142466321;
        
        if (verbose) Rcpp::Rcout << "p_conversion_factor: " << p_conversion_factor<< "\n";
    } 
    else if (unit_P == "Torr")
    {
        p_conversion_factor = 1.0;
        if (verbose) Rcpp::Rcout << "p_conversion_factor: " << p_conversion_factor << "\n";
    } 
    else 
    {
        Rcpp::stop("'unit_P' must be one of 'Torr', 'kPa' or 'mmHg'");
    }
    
    //setting output unit
    if (unit_S == "percent"){
        s_conversion_factor = 100;
        if (verbose) Rcpp::Rcout << "s_conversion_factor: " << s_conversion_factor << "\n";
    } else if (unit_S == "fraction") {
        s_conversion_factor = 1;
        if (verbose) Rcpp::Rcout << "s_conversion_factor: " << s_conversion_factor << "\n";
    } else {
        Rcpp::stop("'unit_S' must be one of 'percent' or 'fraction'");
    }
    
    if (na_handling == "impute")
    {
        na_method = 0;
        if (verbose) Rcpp::Rcout << "na_method: " << na_method << "\n";
    } 
    //else if (!std::strncmp(na_handling.get_cstring(), "complete_cases", 500))
    else if (na_handling == "complete_cases")
    {
        na_method = 1;
        if (verbose) Rcpp::Rcout << "na_method: " << na_method << "\n";
    } 
    else 
    {
        Rcpp::stop("'na_handling' must be one of 'impute' or 'complete_cases'");
    }
    
    //throw errors if input vectors are different
    if (ph.size() > 0 && ph.size() != po2.size()){
        Rcpp::stop("'pH' is of different length than 'po2'!");
    }
    if (temp.size() > 0 && temp.size() != po2.size()){
        Rcpp::stop("'Temp' is of different length than 'po2'!");
    }
    if (be.size() > 0 && be.size() != po2.size()){
        Rcpp::stop("'BE' is of different length than 'po2'!");
    }
    if (co2.size() > 0 && co2.size() != po2.size()){
        Rcpp::stop("'CO2' is of different length than 'po2'!");
    }
    if (dpg23.size() > 0 && dpg23.size() != po2.size()){
        Rcpp::stop("'DPG23' is of different length than 'po2'!");
    }
    
    //in order to not have to rewrite the same code too often, use a function pointer to the desired method for calculation
    std::unique_ptr<so2f_dash_f> f_dash;
    double (*f)(const double) = nullptr;
    
    if (method == "severinghaus"){
        f = &so2f_severinghaus_f;
    } else if (method == "aberman"){
        f = &so2f_aberman_f;
    } else if (method == "kelman"){
        f = &so2f_kelman_f;
    } else if (method == "dash"){
        f_dash = std::make_unique<so2f_dash_f>();
        //Dash's method = correct the p50, do not transform po2 into standard po2...
        shift_p50 = true;
        correct = false;
    } else {
        Rcpp::stop("'method' must be one of 'severinghaus', 'kelman', 'dash' or 'aberman'");
    }
    
    //create our po2-shift-object to correct the values
    po2shift p;
    
    //variables to store the current or default values
    auto this_ph = default_ph;
    auto this_temp = default_temp;
    auto this_be = default_be;
    auto this_co2 = default_co2;
    auto this_dpg23 = default_dpg23;
    
    //create the out vector
    Rcpp::NumericVector out = Rcpp::clone(po2);
    
    for (auto i = 0u;i<n;++i){
        if(correct && !Rcpp::NumericVector::is_na(out[i])){           
            if (na_method == 0){
                //impute
                this_temp = default_temp;
                this_ph = default_ph;
                this_be = default_be;
                if (temp.size() > 0 && !Rcpp::NumericVector::is_na(temp[i])){
                    this_temp = temp[i];
                }
                if (ph.size() > 0 && !Rcpp::NumericVector::is_na(ph[i])){
                    this_ph = ph[i];
                }
                if (be.size() > 0 &&!Rcpp::NumericVector::is_na(be[i])){
                    this_be = be[i];
                }
            } else if (na_method == 1){
                //only count complete cases
                if (temp.size() == 0 || Rcpp::NumericVector::is_na(temp[i])){
                    out[i] = NA_REAL;
                } else {this_temp =temp[i];}
                if (ph.size() == 0 || Rcpp::NumericVector::is_na(ph[i])){
                    out[i] = NA_REAL;
                } else {this_ph =ph[i];}
                if (be.size() == 0 || Rcpp::NumericVector::is_na(be[i])){
                    out[i] = NA_REAL;
                } else {this_be =be[i];}
            }
        }
        
        if (!Rcpp::NumericVector::is_na(out[i])){
            //the most common case is first here, i.e. when using any other method than dash
            if (!shift_p50){
                out[i] = (*f)(p.po2c(out[i]*p_conversion_factor, this_temp, this_ph, this_be));
            } else {
                //calc and set the p50 for dash. Note that it will impute default values to shift the p50
                this_ph = default_ph;
                this_temp = default_temp;
                this_co2 = default_co2;
                this_dpg23 = default_dpg23;
                if (ph.size() > 0 && !Rcpp::NumericVector::is_na(ph[i])){
                    this_ph = ph[i];
                }
                if (temp.size() > 0 && !Rcpp::NumericVector::is_na(temp[i])){
                    this_temp = temp[i];
                }
                if (co2.size() > 0 && !Rcpp::NumericVector::is_na(co2[i])){
                    this_co2 = co2[i]*p_conversion_factor;
                }
                if (dpg23.size() > 0 && !Rcpp::NumericVector::is_na(dpg23[i])){
                    this_dpg23 = dpg23[i];
                }
                //shift the p50 and calculate
                f_dash->calc_and_set_p50(this_ph, this_temp, this_co2, this_dpg23);
                out[i] = f_dash->po2_to_s(out[i]*p_conversion_factor);
                
            }
            //shift into output units
            out[i]*=s_conversion_factor;
        } 
    }
    return out;
}

//The inverse of the pao2-so2-function, for calculating pao2 from spo2 using any of the 4 methods.
//This also features some convenience functions such as picking units (even though the program internally uses
//Torr for pressure and
//fraction (0-1) for saturation)
//[[Rcpp::export]]
Rcpp::NumericVector spo2_to_pao2(Rcpp::NumericVector so2, 
                                 Rcpp::Nullable<Rcpp::NumericVector> Temp = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> pH = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> BE = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> pCO2 = R_NilValue,
                                 Rcpp::Nullable<Rcpp::NumericVector> DPG23 = R_NilValue,
                                 Rcpp::String method = "kelman", 
                                 Rcpp::String unit_S = "fraction", 
                                 Rcpp::String unit_P = "Torr",
                                 Rcpp::String na_handling = "impute",
                                 double default_temp = 37.0,
                                 double default_ph = 7.4,
                                 double default_be = 0,
                                 double default_co2 = 40,
                                 double default_dpg23 = 4.65e-3,
                                 double default_p50 = 26.8,
                                 double accuracy = 0.95,
                                 std::size_t max_iterations = 30,
                                 bool correct = false,
                                 bool use_ellis = false,
                                 bool verbose = false){
    //set up conversion factors and settings for the loop later on
    double s_conversion_factor = 1.0;
    double p_conversion_factor = 1.0;
    bool shift_p50 = false;
    std::size_t n = so2.size();
    int na_method = 0;
    
    //create default vectors
    Rcpp::NumericVector temp = Rcpp::NumericVector::create();
    Rcpp::NumericVector ph = Rcpp::NumericVector::create();
    Rcpp::NumericVector be = Rcpp::NumericVector::create();
    Rcpp::NumericVector co2 = Rcpp::NumericVector::create();
    Rcpp::NumericVector dpg23 = Rcpp::NumericVector::create();
    if (Temp.isNotNull()){
        temp = Temp;
    }
    
    if (pH.isNotNull()){
        ph = pH;
    }
    
    if (BE.isNotNull()){
        be = BE;
    }
    
    if (pCO2.isNotNull()){
        co2 = pCO2;
    }
    
    if (DPG23.isNotNull()){
        dpg23 = DPG23;
    }
    
    //setting p factor
    if (unit_P == "kPa")
    {
        p_conversion_factor = 1.0/7.50062;
        if (verbose) Rcpp::Rcout << "p_conversion_factor: " << p_conversion_factor << "\n";
    } 
    else if (unit_P == "mmHg")
    {
        p_conversion_factor = 1.0/1.000000142466321;
        
        if (verbose) Rcpp::Rcout << "p_conversion_factor: " << p_conversion_factor<< "\n";
    } 
    else if (unit_P == "Torr")
    {
        p_conversion_factor = 1.0;
        if (verbose) Rcpp::Rcout << "p_conversion_factor: " << p_conversion_factor << "\n";
    } 
    else 
    {
        Rcpp::stop("'unit_P' must be one of 'Torr', 'kPa' or 'mmHg'");
    }
    
    //setting s factor
    if (unit_S == "percent"){
        s_conversion_factor = 0.01;
        if (verbose) Rcpp::Rcout << "s_conversion_factor: " << s_conversion_factor << "\n";
    } else if (unit_S == "fraction") {
        s_conversion_factor = 1;
        if (verbose) Rcpp::Rcout << "s_conversion_factor: " << s_conversion_factor << "\n";
    } else {
        Rcpp::stop("'unit_S' must be one of 'percent' or 'fraction'");
    }
    
    if (na_handling == "impute")
    {
        na_method = 0;
        if (verbose) Rcpp::Rcout << "na_method: " << na_method << "\n";
    } 
    //else if (!std::strncmp(na_handling.get_cstring(), "complete_cases", 500))
    else if (na_handling == "complete_cases")
    {
        na_method = 1;
        if (verbose) Rcpp::Rcout << "na_method: " << na_method << "\n";
    } 
    else 
    {
        Rcpp::stop("'na_handling' must be one of 'impute' or 'complete_cases'");
    }
    
    //throw errors if input vectors are different
    if (ph.size() > 0 && ph.size() != so2.size()){
        Rcpp::stop("'pH' is of different length than 'so2'!");
    }
    if (temp.size() > 0 && temp.size() != so2.size()){
        Rcpp::stop("'Temp' is of different length than 'so2'!");
    }
    if (be.size() > 0 && be.size() != so2.size()){
        Rcpp::stop("'BE' is of different length than 'so2'!");
    }
    if (co2.size() > 0 && co2.size() != so2.size()){
        Rcpp::stop("'CO2' is of different length than 'so2'!");
    }
    if (dpg23.size() > 0 && dpg23.size() != so2.size()){
        Rcpp::stop("'DPG23' is of different length than 'so2'!");
    }
    
    //in order to not have to rewrite the same code too often, use a function pointer to the desired method for calculation
    std::unique_ptr<so2f_dash_f> f_dash;
    double (*f)(const double) = nullptr;
    
    if (method == "severinghaus"){
        f = &so2f_severinghaus_f;
    } else if (method == "aberman"){
        f = &so2f_aberman_f;
    } else if (method == "kelman"){
        f = &so2f_kelman_f;
    } else if (method == "dash"){
        f_dash = std::make_unique<so2f_dash_f>();
        //Dash's method = correct the p50, do not transform po2 into standard po2...
        shift_p50 = true;
        correct = false;
    } else {
        Rcpp::stop("'method' must be one of 'severinghaus', 'kelman', 'dash' or 'aberman'");
    }
    
    //create our po2-shift-object to correct the values
    po2shift p;
    
    //variables to store the current or default values
    auto this_ph = default_ph;
    auto this_temp = default_temp;
    auto this_be = default_be;
    auto this_co2 = default_co2 / p_conversion_factor;
    auto this_dpg23 = default_dpg23;
    
    //create the out vector
    Rcpp::NumericVector out = Rcpp::clone(so2);
    
    //create the required parameters for the boost bracket_and_solve
    //set up the helper objects for the bracket-and-solve-method
    if (accuracy < 0.01 || accuracy > 0.99){
        Rcpp::stop("'accuracy' must be between 0.01 and 0.99");
    }
    int digits = std::numeric_limits<double>::digits * accuracy; 
    boost::math::tools::eps_tolerance<double> tol(digits);
    double guess = 50;
    double factor = 2;
    bool rising = true;
    std::pair<double, double> res;
    std::size_t max_iter = max_iterations;
    
    
    for (auto i = 0u;i<n;++i){
        if(correct && !Rcpp::NumericVector::is_na(out[i])){           
            if (na_method == 0){
                //impute
                this_temp = default_temp;
                this_ph = default_ph;
                this_be = default_be;
                if (temp.size() > 0 && !Rcpp::NumericVector::is_na(temp[i])){
                    this_temp = temp[i];
                }
                if (ph.size() > 0 && !Rcpp::NumericVector::is_na(ph[i])){
                    this_ph = ph[i];
                }
                if (be.size() > 0 &&!Rcpp::NumericVector::is_na(be[i])){
                    this_be = be[i];
                }
            } else if (na_method == 1){
                //only count complete cases
                //only count complete cases
                if (temp.size() == 0 || Rcpp::NumericVector::is_na(temp[i])){
                    out[i] = NA_REAL;
                } else {this_temp = temp[i];}
                if (ph.size() == 0 || Rcpp::NumericVector::is_na(ph[i])){
                    out[i] = NA_REAL;
                } else {this_ph = ph[i];}
                if (be.size() == 0 || Rcpp::NumericVector::is_na(be[i])){
                    out[i] = NA_REAL;
                } else {this_be = be[i];}
            }
        }
        
        if (!Rcpp::NumericVector::is_na(out[i])){
             //the most common case is first here, i.e. when using any other method than dash
            if (!shift_p50){
                //use boost to bracket the root of the error function f(guess) - so2
                max_iter = max_iterations;
                res = boost::math::tools::bracket_and_solve_root([&](double x){return (*f)(p.po2c(x, this_temp, this_ph, this_be)) - so2[i]*s_conversion_factor;}, guess, factor, rising, tol, max_iter);
            } else {
                //calc and set the p50 for dash. Note that it will impute default values to shift the p50
                this_ph = default_ph;
                this_temp = default_temp;
                this_co2 = default_co2;
                this_dpg23 = default_dpg23;
                if (ph.size() > 0 && !Rcpp::NumericVector::is_na(ph[i])){
                    this_ph = ph[i];
                }
                if (temp.size() > 0 && !Rcpp::NumericVector::is_na(temp[i])){
                    this_temp = temp[i];
                }
                if (co2.size() > 0 && !Rcpp::NumericVector::is_na(co2[i])){
                    this_co2 = co2[i]*p_conversion_factor;
                }
                if (dpg23.size() > 0 && !Rcpp::NumericVector::is_na(dpg23[i])){
                    this_dpg23 = dpg23[i];
                }
                //shift the p50 and calculate
                f_dash->calc_and_set_p50(this_ph, this_temp, this_co2, this_dpg23);
                //use boost to bracket the root of the error function f(guess) - so2
                max_iter = max_iterations;
                res = boost::math::tools::bracket_and_solve_root([&](double x){return f_dash->po2_to_s(x) - so2[i]*s_conversion_factor;}, guess, factor, rising, tol, max_iter);
            }
            
            //midway between the brackets!
            out[i] = (res.first + 0.5 * (res.second - res.first));
            //shift into output units
            out[i]*=p_conversion_factor;
        }
        
    }
    
    return out;
}

//calcualte BE from pH and pCO2
//[[Rcpp::export]]
Rcpp::NumericVector calc_be(Rcpp::NumericVector pH,
                                 Rcpp::NumericVector pCO2){
    
    std::size_t n = pH.size();
    Rcpp::NumericVector out = Rcpp::clone(pH);

    
    BE_calc f;
    for (auto i = 0u;i<n;++i){
        f.set_pCO2(pCO2[i]);
        out[i] = f(pH[i]);
    }

    return out;
}

//calcualte pCO2 from BE and pH
//[[Rcpp::export]]
Rcpp::NumericVector calc_pco2(Rcpp::NumericVector pH,
                            Rcpp::NumericVector BE){
    
    std::size_t n = pH.size();
    Rcpp::NumericVector out = Rcpp::clone(pH);

    
    pCO2_calc f;
    for (auto i = 0u;i<n;++i){
        f.set_BE(BE[i]);
        out[i] = f(pH[i]);
    }

    return out;
}

//calcualte pH from BE and pCO2, using bracket-and-solve for the equation by Siggaard-Andersen
//[[Rcpp::export]]
Rcpp::NumericVector calc_ph(Rcpp::NumericVector BE,
                            Rcpp::NumericVector pCO2,
                            double accuracy = 0.95,
                            std::size_t max_iterations = 30){
    std::size_t n = BE.size();
    Rcpp::NumericVector out = Rcpp::clone(BE);
    BE_calc f;

    //create the required parameters for the boost bracket_and_solve
    //set up the helper objects for the bracket-and-solve-method
    if (accuracy < 0.01 || accuracy > 0.99){
        Rcpp::stop("'accuracy' must be between 0.01 and 0.99");
    }
    int digits = std::numeric_limits<double>::digits * accuracy; 
    boost::math::tools::eps_tolerance<double> tol(digits);
    double guess = 7.4;
    double factor = 2;
    bool rising = true;
    std::pair<double, double> res;
    std::size_t max_iter = max_iterations;

    for (auto i = 0u;i<n;++i){
        f.set_pCO2(pCO2[i]);
        max_iter = max_iterations;
        res = boost::math::tools::bracket_and_solve_root([&](double x){return f(x) - BE[i];}, guess, factor, rising, tol, max_iter);

        //midway between the brackets!
        out[i] = (res.first + 0.5 * (res.second - res.first));
    }

    return out;

}