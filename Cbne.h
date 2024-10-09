#pragma once

#include "Complex.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <cmath>
#include "Config.h"
#include "Statistics.h"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/tools/polynomial.hpp>

#define choose boost::math::binomial_coefficient<double>

enum class PolyType {
    APERS,
    MUSCO,
    MUSCO_COMPRESSED
};

// helper function
std::string print_face(const Complex::Face& face)
{
    std::string res = "[";

    for (auto it = face.begin(); it != face.end(); ++it) {
        res += std::to_string(*it);
        if(it != std::prev(face.end())){
            res += ", ";
        }
    }

    res += "]";
    return res;
}


// Converts a Boost graph in adjacency list format to adjacency matrix 
// format. The Boost::copy_graph function fails silently
void convert_to_matrix(const GraphL& g1, Graph& g2)
{
    for(auto edge_it = edges(g1).first; edge_it != edges(g1).second; ++edge_it ) {
        int a = source(*edge_it, g1);
        int b = target(*edge_it, g1);
        add_edge(a, b, g2);    
    }
}

/* Samples z steps from Markov chain starting at face @param j0
 * @param n is the number of vertices in @param complex
 */
double sample_markov_chain(Complex::Face& j0, const Complex& complex, unsigned n, unsigned z, Statistics& stats) 
{

    Complex::Face face = j0;
    double product = 1;

    for(unsigned i = 0; i < z; i++) {
        unsigned up_degree = complex.number_of_up_neighbours(face);
        double  diagonal  = up_degree + face.size();
        diagonal = 1 - (diagonal / n);
        std::vector<std::pair<Complex::Face,int>> neighbours = complex.get_neighbours(face, n);
        double column_norm = (double)neighbours.size() / n + diagonal;

        std::random_device rd;  // Seed for the random number generator
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::uniform_real_distribution<double> dis(0.0, 1.0);


        // Generate a random number between 0 and 1
        double random_num = dis(gen);

        int sign = 1;

        if(random_num > (diagonal / column_norm)){
            std::uniform_int_distribution<int> dis(0, neighbours.size() - 1);
            int index = dis(gen);

            sign = neighbours[index].second;
            face = neighbours[index].first;
        }

        product = product * sign * column_norm;
    }

    if(face == j0){
        stats.incr_number_of_cycles();
        return product;
    } else
        return 0;
}

// Returns the coefficents of the @param n th 
// Chebychev polynomial of the first kind
std::vector<double> get_chebyshev_coefficients(int n) {
    using boost::math::tools::polynomial;
    
    // Base cases
    if (n == 0) return {1};
    if (n == 1) return {0, 1};
    
    // Initial polynomials
    polynomial<double> T0({1});          // T0(x) = 1
    polynomial<double> T1({0, 1});       // T1(x) = x
    
    polynomial<double> Tn;
    
    // Recursive generation
    for (unsigned i = 2; i <= n; ++i) {
        Tn = 2 * polynomial<double>({0, 1}) * T1 - T0;
        T0 = T1;
        T1 = Tn;
    }
    
    // Return coefficients of Tn
    std::vector<double> coefficients = Tn.data();

    return coefficients;
}


std::vector<double> get_coefficients_apers(int power, int degree) {
    std::vector<double>  coeffs(degree + 1);
    double denom = std::pow(2.0, power);
    for(int j = 0; j <= degree; j++){
        // check same parity
        if((j & 1) == (power & 1)){
            int multiplier = j == 0 ? 1 : 2;
            coeffs[j] = multiplier * (choose(power, (power - j) / 2) /  denom);
        } else {
            coeffs[j] = 0;
        }
    }
    return coeffs;
}

std::vector<double> get_coefficients_musco(const std::vector<double>& cheby_coeffs, double gamma, double d)
{  
    std::vector<double> res((int)(d + 1));

    double factor = 1 / (1 - gamma);

    // evaluate the chebychev polynomial at 1 / 1 - \gamma
    //double sum = cheby_coeffs[(int)d];
    //for(unsigned j = 1; j < d + 1; j++){
    //    sum += std::pow(factor, j) * cheby_coeffs[j * (int)(d + 1) + (int)d];
     // }

    double outer_factor = 1 / boost::math::chebyshev_t((unsigned)d, factor);

    for(unsigned j = 0; j < d + 1; j++){
        if(j == 0){
            res[j] = cheby_coeffs[0] * outer_factor;
        } else {
            res[j] = cheby_coeffs[j] * outer_factor * std::pow(factor, j);
        }
    }

    return res;
}

std::pair<std::vector<double>, double> get_polynomial_apers(double gamma, double epsilon, const Config& config)
{
    double d = std::ceil( (std::sqrt(2 / gamma)) * log(6 / epsilon));
    double r = std::ceil(1 / gamma * log(3 / epsilon));

    // TODO do something with config here
    // Need to check what "r" is, AYN.

    std::vector<std::vector<double>> coeffs_cheby;
    for(unsigned j = 0; j < d + 1; j++){
        coeffs_cheby.push_back(get_chebyshev_coefficients(j));
    }

    auto coeffs = get_coefficients_apers(r, d);

    assert(coeffs.size() == coeffs_cheby.size());

    std::vector<double> coeffs_new(coeffs.size());

    // multiply out the coefficients
    for(unsigned j = 0; j < coeffs.size(); j++){
        double coeff = coeffs[j];
        for(unsigned k = 0; k < coeffs_cheby[j].size(); k++){
            coeffs_cheby[j][k] = coeffs_cheby[j][k] * coeff;
        }
    }

    // collect like terms
    for(unsigned j = 0; j < coeffs.size(); j++){
        double sum = 0;
        for(unsigned k = j; k < coeffs.size(); k++){
            sum += coeffs_cheby[j][k];
        }
        coeffs_new[j] = sum;
    }

    return std::make_pair(coeffs_new, d);
}

std::pair<std::vector<double>, double> get_polynomial_musco(double gamma, double epsilon, const Config& config)
{
    double d = std::ceil( (std::sqrt(1 / gamma)) * log(4 / epsilon));

    if(config.get_deg_limit() != -1) {
        d = config.get_deg_limit();
    }

    auto coeffs_cheby =  get_chebyshev_coefficients(d);
    auto coeffs_new   = get_coefficients_musco(coeffs_cheby, gamma, d);
    return std::make_pair(coeffs_new, d);
}

std::vector<double> get_coefficients_musco_compressed(const std::vector<double>& cheby_coeffs, double gamma, double d)
{
    std::vector<double> res((int)(d + 1));

    auto cheby_coeff = [](double d, unsigned k, double gamma){
        if(((int)d & 1) != (k & 1)){
            return 0.0;
        }
        double m = (d - k) / 2;

        auto term1 = std::pow(-1, m);
        auto term2 = std::pow(2, -2 * m + d -1);
        auto term3 = choose(-m + d - 1, d - 2 * m) + choose(d - m, m);
        auto term4 = std::pow(1 - gamma, k);

        return (term1 * term2 * term3) / term4;
    };

    double factor = (1  + gamma) / (1 - gamma); 

    // evaluate the chebychev at factor and store in outer_factor
    double outer_factor = 1 / boost::math::chebyshev_t((unsigned)d, factor);

    for(unsigned j = 0; j < d + 1; j++){
        double sum = 0;
        for(unsigned k = j; k < d + 1; k++){
            double term1 = choose(k, j);
            double term2 = std::pow(-0.5 * (1 - gamma), k - j);
            double term3 = cheby_coeff(d, k, (gamma + 1) / 2);
            sum += term1 * term2 * term3;
        }
        res[j] = outer_factor * sum;
    }

    return res;
}

std::pair<std::vector<double>, double> get_polynomial_musco_compressed(double gamma, double epsilon, const Config& config)
{
    std::vector<double> coeffs_cheby, coeffs_new;

    // TODO is this the right d?
    double d = std::ceil( (std::sqrt(1 / gamma)) * log(2 / epsilon));

    if(config.get_deg_limit() != -1) {
        d = config.get_deg_limit();
    }

    while(true){
        coeffs_cheby = get_chebyshev_coefficients(d);
        coeffs_new   = get_coefficients_musco_compressed(coeffs_cheby, gamma, d);
        if(std::abs(coeffs_new[0]) < epsilon / 2){
            break;
        }
        d = d + 1;
    }

    return std::make_pair(coeffs_new, d);
}

std::optional<double> cbne_chebychev(int n, int k, const Complex& complex, double gamma, PolyType ptype, double one_norm, const Config& config, Statistics& stats)
{
    double epsilon = config.get_epsilon();

    double d;
    std::vector<double> poly;
    if(ptype == PolyType::APERS){
        //std::tie(poly, d) = get_polynomial_apers(gamma, epsilon, config);
    } else if (ptype == PolyType::MUSCO){
        std::tie(poly, d) = get_polynomial_musco(gamma, epsilon, config);
    } else {
        std::tie(poly, d) = get_polynomial_musco_compressed(gamma, epsilon, config);
    }

    auto calculate_shot_count = [&](unsigned j){
        int error_factor = ptype == PolyType::APERS ? 3 : 2;
        double delta =  epsilon / (error_factor * std::ceil(d / 2) * std::abs(poly[j]));
        //double delta =  epsilon / (error_factor * ((d + 1) / 2) * std::abs(poly[j]));

        double eta = 1 - std::pow(0.90, 1 / std::ceil(d / 2) );
        //double eta = 1 - std::pow(0.90, 2 / (d + 1)) ;
        double base = one_norm != 0 && config.use_one_norm() ? one_norm : 2;
        return std::ceil((4 * pow(base, j * 2) * log(2 / eta))  / (2 * delta * delta));
    };

    double total = 0;
    double total_step = 0;
    std::vector<double> shot_counts((unsigned)d + 1, 0);
    std::vector<double> shot_counts_completed((unsigned)d + 1, 0);

    for(unsigned j = 0; j < d + 1; j++){
        if(poly[j] != 0){
            double p = calculate_shot_count(j);
            total += p;
            total_step += p * j;
            shot_counts[j] = p;
        }
    }

    double shot_count_limit = total;
    if(config.get_iter_limit() != -1) {
        shot_count_limit = config.get_iter_limit();
    }
    if(shot_count_limit != total){
        // In this case user has asked for a lower shot count
        double frac = shot_count_limit / total;
        total = 0;
        total_step = 0;
        for(unsigned j = 0; j < d + 1; j++){
            if(poly[j] != 0){
                shot_counts[j] = std::ceil(shot_counts[j] * frac);
                total += shot_counts[j];
                total_step += shot_counts[j] * j;
            }
        }
    }

    if(config.output_count() || config.output_count_step()){

        if(config.output_count()){
            std::cout << "Shot count: " << total  << std::endl;
        }

        if(config.output_count_step()){
            std::cout << "Step count: " << total_step << std::endl;
        }

        return {};
    }

    stats.set_sample_count(total);
    stats.set_walk_length(d);

    int limit = 1;
#if OUTPUT_INTERMEDIATE
    std::ofstream out(config.get_out_path());
    limit = config.get_num_data_points();
#endif

    std::vector<double> estimates(poly.size(), 0);
    std::vector<int> shots(poly.size(), 0);
    double res = 0;

    for(unsigned count = 0; count < limit; count++){

        for(unsigned j = 0; j < d + 1; j++){
            if(poly[j] != 0){

                double sum = 0;
                int i = 0;
    
                double p; // = std::ceil(shot_counts[j]);
                if(count != limit - 1){

                    auto expr = [&](unsigned datapoint){
                        return std::floor(datapoint * (shot_counts[j] / (double)limit));
                    };

                    p = expr(count + 1) - expr(count);

                    //p = std::round(shot_counts[j] / (double)limit);
                    shot_counts_completed[j] = shot_counts_completed[j] + p;
                } else {
                    // use up the remaining budget
                    p = std::ceil(shot_counts[j] - shot_counts_completed[j]);
                }

                if(p == 0) continue; 

                while(i < p) {
                    Complex::Face j0 = complex.sample_from_complex(k);
                    double Y_z = sample_markov_chain(j0, complex, n, j, stats);
                    sum += Y_z;
                    i++;
                }

                estimates[j] += sum;
                shots[j] += (int)p;
            }
        }


        res = 0;
        for(unsigned j = 0; j < d + 1; j++){
            if(shots[j] != 0){
                res += poly[j] * (estimates[j] / shots[j]);
            }
        }

#if OUTPUT_INTERMEDIATE        
        out << res;
        if(count != limit - 1){
            out << ", ";
            out.flush();
        }
#endif
    }

#if OUTPUT_INTERMEDIATE
    out.close();
#endif
    return res;
}

std::optional<double> cbne(int n, int k, const Complex& complex, double gamma, double one_norm, const Config& config, Statistics& stats)
{
    // TODO consider using GMP library to achieve better precision
    double epsilon = config.get_epsilon();
    double delta = epsilon / 2;
    double z     = std::ceil((1 / gamma) * log(2 / epsilon));

    if(config.get_deg_limit() != -1) {
        z = config.get_deg_limit();
    }

    double base  = one_norm != 0  && config.use_one_norm() ? one_norm : 2;
    double p     = (4 * std::pow(base, 2 * z) * log(2/0.1)) / (2 * delta * delta) ;

    if(config.output_count()){
        std::cout << "Shot count: " << p << std::endl;
    }

    if(config.output_count_step()){
        std::cout << "Step count: " << p * z << std::endl;
    }

    if(config.output_count() || config.output_count_step()){
        return {};
    }

    double sum = 0;
    int i = 0;

    if(config.get_iter_limit() != -1) {
        p = config.get_iter_limit();
    }

    p = std::ceil(p);

    stats.set_sample_count(p);
    stats.set_walk_length(z);

#if USE_TIME_LIMIT
    std::chrono::duration<double> duration(0);
    auto start = std::chrono::system_clock::now();
#endif

#if OUTPUT_INTERMEDIATE
    std::ofstream out(config.get_out_path());
    int number_of_data_points = config.get_num_data_points();
    double frac = std::ceil(p / number_of_data_points);
#endif

    while(i < p 
#if USE_TIME_LIMIT
        && duration.count() < config.get_time_limit()
#endif
        ) {
        Complex::Face j0 = complex.sample_from_complex(k);
        double Y_z = sample_markov_chain(j0, complex, n, ceil(z), stats);
        sum += Y_z;
        i++;
#if USE_TIME_LIMIT
        // TODO benchmark penalty for calculating time every round of the loop
        //auto curr = std::chrono::system_clock::now();
        //duration = curr - start;
#endif

#if OUTPUT_INTERMEDIATE
        if((i % (int)frac == 0) || (i == p) || i == 10){
            out << (sum / i);
            if(i != p){
                out << ", ";
                out.flush();
            }
        }
#endif
    }

#if OUTPUT_INTERMEDIATE
    out.close();
#endif

    return sum / i;
}