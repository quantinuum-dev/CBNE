#include <iostream>
#include <fstream>
//#include <gmpxx.h>
#include "Complex.h"
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include "Cbne.h"
#include "Config.h"
#include "Statistics.h"
#include <boost/graph/copy.hpp>

int main(int argc, char** argv){

    try {

	TCLAP::CmdLine cmd("Command description message", ' ', "0.9");

        TCLAP::ValueArg<std::string> path_to_graph("p","path","Path to .graphml input file",false,".","string");
#if OUTPUT_INTERMEDIATE
        TCLAP::ValueArg<std::string> path_to_output("o","out","Path to output .csv file",false,"","string");
        TCLAP::ValueArg<int>         num_data_points("n","num_data_points","Number of data points to output to .csv file",false,1,"int");
#endif
        TCLAP::ValueArg<double>      epsilon("e","epsilon","Epsilon value",false,0.1,"double");
#if USE_TIME_LIMIT
        TCLAP::ValueArg<int>         time_limit("t","time_limit","Time limit for which to run algorithm",false,-1,"int");
#endif
        TCLAP::ValueArg<int>         iter_limit("i","iter_limit","Number of samples to use",false,-1,"int");
        TCLAP::ValueArg<int>         deg_limit("d","deg_limit","Degree to use",false,-1,"int");
        TCLAP::SwitchArg             output_shot_count("s","output_shot_count","Output the shot count and exit", false);
        TCLAP::SwitchArg             output_step_count("c","output_step_count","Output the step count and exit", false);
        TCLAP::SwitchArg             use_one_norm("u","use_one_norm","Use the one norm of H if it is present in the input file", true);
        TCLAP::ValueArg<std::string> cbne_version("a","cbne_version","There are various different versions of CBNE algorithm. Use this option to select amongst them. Valid values are in [cbne, cbneCheby, cbneMusco]",false,"cbne","string");


        cmd.add( path_to_graph );
#if OUTPUT_INTERMEDIATE
        cmd.add( path_to_output );
        cmd.add( num_data_points );
#endif
	cmd.add( epsilon );
#if USE_TIME_LIMIT
        cmd.add( time_limit );
#endif
        cmd.add( iter_limit );
        cmd.add( deg_limit );
        cmd.add( output_shot_count );
        cmd.add( output_step_count );
        cmd.add( use_one_norm );
        cmd.add( cbne_version );

        cmd.parse( argc, argv );

        Config config(cmd);

        std::ifstream graph_file;
        graph_file.open(config.get_path(), std::ifstream::in);

        if(!graph_file){
            std::cout << "failed to open input file " << config.get_path() << std::endl;
            return -1;
        }

        GraphL g;
        boost::dynamic_properties dp(boost::ignore_other_properties);

        auto  betti = get(&VertexProps::betti, g);
        auto  dimk  = get(&VertexProps::dimk, g);
        auto  gap   = get(&VertexProps::gap, g);
        auto  norm  = get(&VertexProps::one_norm, g);

        dp.property("betti",betti);
        dp.property("dimk",dimk);
        dp.property("gap",gap);
        dp.property("norm", norm);

        boost::read_graphml(graph_file, g, dp);

        double spectral_gap  = g[vertex(0, g)].gap;
        int    dimension     = g[vertex(0, g)].dimk;
        double one_norm      = g[vertex(0, g)].one_norm;
        double betti_est     = g[vertex(0, g)].betti;
        int    n             = num_vertices(g);

        Graph gM(num_vertices(g));
        // boost::copy_graph does not copy between adjacency lists and matrices (silently fails)
        // See https://stackoverflow.com/questions/66977885/how-can-i-read-in-a-graph-to-an-adjacency-matrix-in-the-boost-graph-library
        convert_to_matrix(g, gM);
        Complex complex(gM);


        Statistics stats;

        //double estimate = apre(n, dimension, complex, spectral_gap, eps, t_limit, outp, i_limit);
        std::optional<double> estimate;
        if(config.get_cbne_version() == "cbne"){
           estimate = cbne(n, dimension, complex, spectral_gap, one_norm, config, stats);
        } else if(config.get_cbne_version()  == "cbneCheby"){
           estimate = cbne_chebychev(n, dimension, complex, spectral_gap, PolyType::APERS, one_norm, config, stats);
        } else if(config.get_cbne_version()  == "cbneMusco"){
           estimate = cbne_chebychev(n, dimension, complex, spectral_gap, PolyType::MUSCO, one_norm, config, stats);
        } else {
           estimate = cbne_chebychev(n, dimension, complex, spectral_gap, PolyType::MUSCO_COMPRESSED, one_norm, config, stats);
        }

        if(estimate){
           std::cout << "\nBetti estimate: " << *estimate << "\n" << std::endl;
           std::cout << stats.to_string() << std::endl;
        }
        // TODO consider warning if we are outside of epsilon bound on ground truth?
        // Difficulty is that we don't always have the ground truth
    } 
    catch(const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    catch (TCLAP::ArgException &e)  // catch exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
        return 1;
    }
}
