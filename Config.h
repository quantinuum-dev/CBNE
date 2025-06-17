#pragma once

#include <unordered_map>
#include <string>
#include <tclap/CmdLine.h>
#include <stdexcept>

// set to 1 to output intermediate values to file
#define OUTPUT_INTERMEDIATE 1
#define USE_TIME_LIMIT 0

class Config {
public:

    Config(const TCLAP::CmdLine& cmd){
        auto args = cmd.getArgList();
        for(auto it = args.begin(); it != args.end(); it++)
        {
            std::string name = (*it)->getName();
            if(name == "cbne_version"){
                auto casted_arg = static_cast<TCLAP::ValueArg<std::string>*>(*it);
                cbne_version = casted_arg->getValue();

                if(cbne_version != "cbne" && cbne_version != "cbneCheby" && cbne_version != "cbneMusco" && cbne_version != "cbneCompressed"){
                    throw std::invalid_argument("Value " + cbne_version + " is not valid for option -a. Please set a valid algorithm version.");
                }

            } else if (name == "output_shot_count") {
                auto casted_arg = static_cast<TCLAP::SwitchArg*>(*it);
                output_shot_count = casted_arg->getValue();

            } else if (name == "output_step_count") {
                auto casted_arg = static_cast<TCLAP::SwitchArg*>(*it);
                output_step_count = casted_arg->getValue();

            } else if (name == "out") {
                auto casted_arg = static_cast<TCLAP::ValueArg<std::string>*>(*it);
                out = casted_arg->getValue();
            
            } else if (name == "num_data_points") {
                auto casted_arg = static_cast<TCLAP::ValueArg<int>*>(*it);
                number_of_data_points = casted_arg->getValue();

                if(number_of_data_points < 1){
                    throw std::invalid_argument("Value " + std::to_string(number_of_data_points) + " is not valid for option -n. Please choose a number of data points >= 1.");                
                }

            } else if (name == "time_limit"){
                auto casted_arg = static_cast<TCLAP::ValueArg<int>*>(*it);
                time_limit = casted_arg->getValue();

                if(time_limit < -1){
                    throw std::invalid_argument("Value " + std::to_string(time_limit) + " is not valid for option -t. Please set a time limit >= -1.");                
                }

            } else if (name == "path"){
                auto casted_arg = static_cast<TCLAP::ValueArg<std::string>*>(*it);
                path = casted_arg->getValue();

            } else if (name == "iter_limit") {
                auto casted_arg = static_cast<TCLAP::ValueArg<int>*>(*it);
                iter_limit = casted_arg->getValue();

                if(iter_limit < -1){
                    throw std::invalid_argument("Value " + std::to_string(iter_limit) + " is not valid for option -t. Please set an iteration limit >= -1.");                
                }

            } else if( name == "deg_limit"){
                auto casted_arg = static_cast<TCLAP::ValueArg<int>*>(*it);
                deg_limit = casted_arg->getValue();

                if(deg_limit < -1){
                    throw std::invalid_argument("Value " + std::to_string(deg_limit) + " is not valid for option -d. Please set a time limit >= -1.");                
                }

            } else if (name == "use_one_norm") {
                auto casted_arg = static_cast<TCLAP::SwitchArg*>(*it);
                one_norm = casted_arg->getValue();

            } else if (name == "epsilon"){
                auto casted_arg = static_cast<TCLAP::ValueArg<double>*>(*it);
                epsilon = casted_arg->getValue();

            } else {

                // TODO deal with unexpected arguments that are not help, version or ignore_rest
            }
        }

    }

    std::string get_path() const { 
        return path;
    }

    int get_iter_limit() const {
        return iter_limit;
    }

    int get_deg_limit() const {
        return deg_limit;
    }

    std::string get_out_path() const {
        return out;
    }

    int get_num_data_points() const {
        return number_of_data_points;
    }

    bool output_count() const {
        return output_shot_count;
    }

    bool output_count_step() const {
        return output_step_count;
    }

    bool use_one_norm() const {
        return one_norm;
    }

    double get_epsilon() const {
        return epsilon;
    }
    
    int get_time_limit() const {
        return time_limit;
    }

    std::string get_cbne_version() const {
        return cbne_version;
    }

private:
    int iter_limit;
    int deg_limit;
    std::string path;
    std::string out;
    int number_of_data_points;
    bool output_shot_count;
    bool output_step_count;
    bool one_norm;
    double epsilon;
    int time_limit;  
    std::string cbne_version;
};