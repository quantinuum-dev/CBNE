#include <catch2/catch_test_macros.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include "../Apres.h"
#include <fstream>

std::tuple<double, int, int> get_data(std::string path, GraphL& g)
{
    std::ifstream graph_file;
    boost::dynamic_properties dp(boost::ignore_other_properties);

    auto  betti = get(&VertexProps::betti, g);
    auto  dimk  = get(&VertexProps::dimk, g);
    auto  gap   = get(&VertexProps::gap, g);

    dp.property("betti",betti);
    dp.property("dimk",dimk);
    dp.property("gap",gap);

    graph_file.open(path, std::ifstream::in);
    boost::read_graphml(graph_file, g, dp);

    double spectral_gap  = g[vertex(0, g)].gap;
    int    dimension     = g[vertex(0, g)].dimk;
    int    n             = num_vertices(g);


    return std::tuple(spectral_gap, dimension, n);
}

TEST_CASE( "Testing sample graph 1" ) {

    GraphL g;
    auto [spectral_gap, dimension, n] = get_data("../../Apres-Python/sample-graphs/more_exampleHoles_1.graphml", g);
    Graph gM(n);

    Complex complex(gM);

    double estimate = apre(n, dimension, complex, spectral_gap, 0.1, 30);

    REQUIRE( 0.24 <= estimate );
    REQUIRE( estimate <= 0.26 );
}

TEST_CASE( "Testing sample graph 2" ) {

    GraphL g;
    auto [spectral_gap, dimension, n] = get_data("../../Apres-Python/sample-graphs/more_exampleHoles_2.graphml", g);
    Graph gM(n);

    Complex complex(gM);

    double estimate = apre(n, dimension, complex, spectral_gap, 0.1, 30);

    REQUIRE( 0.42 <= estimate );
    REQUIRE( estimate <= 0.46 );
}


TEST_CASE( "Testing sample graph 3", "[3]" ) {

    GraphL g;
    auto [spectral_gap, dimension, n] = get_data("../../Apres-Python/sample-graphs/more_exampleHoles_3.graphml", g);
    Graph gM(n);

    Complex complex(gM);

    double estimate = apre(n, dimension, complex, spectral_gap, 0.01, 130);

    REQUIRE( 0.09 <= estimate );
    REQUIRE( estimate <= 0.11 );
}
