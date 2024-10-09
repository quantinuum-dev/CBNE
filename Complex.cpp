#include "Complex.h"
#include <numeric>
#include <random>
#include <iostream>

Complex::~Complex()
{
}

bool Complex::edge_in_graph(unsigned v1, unsigned v2) const 
{
    boost::graph_traits<Graph>::vertex_descriptor u, v;
    u = vertex(v1, graph);
    v = vertex(v2, graph);
    return edge(u,v,graph).second;
}

bool Complex::includes_face(const Face& face) const
{

    for (auto it1 = face.begin(); it1 != face.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != face.end(); ++it2) {
            if(!edge_in_graph(*it1, *it2)){
                return false;
            }
        }
    }

    return true;
}

bool Complex::connected_to_face(const Face& face, unsigned vertex, std::optional<unsigned> exclude) const
{
    for (auto it = face.begin(); it != face.end(); ++it) {
        if(!exclude || *exclude != *it){
            if(!edge_in_graph(*it, vertex)){
                return false;
            }
        }
    }
    return true;
}

int Complex::get_sign(const Face &face, Face::iterator it, unsigned new_value) const
{
    assert(*it != new_value);
    // TODO what about if this is face.end()??
    Face::iterator it_new_value = face.lower_bound(new_value);

    unsigned index1 = std::distance(face.begin(), it);
    unsigned index2 = std::distance(face.begin(), it_new_value);

    // if the value we are removing from the face is less than 
    // the value we are inserting, then we need to reduce the index by 1
    if(*it < new_value) index2--;

    return -1 * (pow(-1, index1) * (pow(-1, index2)));
}

unsigned Complex::number_of_up_neighbours(const Face& face) const
{
    // TODO make sure assertion is not called when in release mode
    assert(includes_face(face));

    unsigned num_of_up_neighbours = 0;

    unsigned n = num_vertices(graph);
    for(unsigned i = 0; i < n; i++){
        if(auto res = face.find(i); res == face.end()){
            // vertex i is not in face
            if(connected_to_face(face, i, {})){
                num_of_up_neighbours++;
            }
        }
    }

    return num_of_up_neighbours;
}

std::vector<std::pair<Complex::Face,int>>  Complex::get_neighbours(const Face &face, unsigned n) const
{
    std::vector<std::pair<Face,int>> neighbours;

    for(unsigned vertex = 0; vertex < n; vertex++){
        // vertex is not in the face
        if(auto res = face.find(vertex); res == face.end()){

            for (auto it = face.begin(); it != face.end(); ++it) {
                if(connected_to_face(face, vertex, std::optional(*it) )){
                    // found a down up neighbour
                    // now check that it is not an up down neighbour
                    if(!edge_in_graph(vertex, *it)){
                        // create neighbour by excluding *it and including vertex
                        int sign = get_sign(face, it, vertex);
                        neighbours.push_back(std::make_pair(copy_with_replacement(face, it, vertex), sign));
                    }
                }
            }
        }
    }

    return neighbours;
}

Complex::Face Complex::sample_from_complex(unsigned k) const
{
    Face face;
    
    unsigned n = num_vertices(graph);

    std::vector<unsigned> vec(n);
    std::iota(vec.begin(), vec.end(), 0);

    do {
        face.clear();
        std::sample(vec.begin(), vec.end(), std::inserter(face, face.begin()), k + 1, std::mt19937 {std::random_device{}()});
    } while(!includes_face(face));
    

    return face;
}
