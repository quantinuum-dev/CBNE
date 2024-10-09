#pragma once

#include <string>

class Statistics {

public:

    Statistics() : walk_length(0), sample_count(0), number_of_cycles(0) {}

    std::string to_string(){
        std::string str = "Statistics: \n";
        str += "Max length of walk " + std::to_string(walk_length) + "\n";
        str += "Number of samples: " + std::to_string(sample_count) + "\n";
        str += "Number of cycles: " + std::to_string(number_of_cycles) + "\n";
        str += "Number of 0s: " + std::to_string(sample_count - number_of_cycles) + "\n";
        return str;
    }

    void set_sample_count(double p){
        sample_count = p;
    }

    void set_walk_length(double z){
        walk_length = z;
    }

    void incr_number_of_cycles(){
        number_of_cycles++;
    }

private:
    unsigned walk_length;
    unsigned sample_count;
    int number_of_cycles;
};