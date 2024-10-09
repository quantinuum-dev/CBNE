#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include "../Apres.h"
#include <fstream>

// TODO share code with tests
// TODO I'm not convinced by the results I am receiving from Catch2 Benchmarking
// perhaps I may be doing something wrong, but the results don't seem to be consistent

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

TEST_CASE("benchmarking graph definition") {

    BENCHMARK_ADVANCED("benchmarking graph definition")(Catch::Benchmark::Chronometer meter) {

        GraphL g;
        auto [spectral_gap, dimension, n] = get_data("../../Apres-Python/sample-graphs/more_exampleHoles_1.graphml", g);
        Graph gM(n);

        Complex complex(gM);    

        meter.measure([spectral_gap, dimension, n, &complex] { return apre(n, dimension, complex, spectral_gap, 0.01, 30, "", 200000); });
    };
}