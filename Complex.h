#pragma once

#include <utility>
#include <set>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

// TODO namespaces? Or overkill for a small project like this?

using namespace boost;

struct VertexProps {
    long dimk = 0;
    double betti = 0;
    double gap = 0;
    double one_norm = 0;
};

// Benchmarking suggests a slight performance benefit for using adjacency matrix over adjacency list
// This is unsurprising since we are only accessing edges which is O(1) in matrix representation
// We have to read into a list due to vagaries of the Boost library which does not let us read
// into a matrix...
typedef adjacency_matrix<undirectedS, VertexProps> Graph;
typedef adjacency_list<vecS, vecS, undirectedS, VertexProps> GraphL;

class Complex
{
public:
    // TODO check performance of changing to std::vector instead
    // given that faces are relatively small, the worse asymptotic performance
    // may be offset by having contiguous memory and reduced need allocate
    using Face = std::set<unsigned>;

private:

    const Graph& graph;
    bool edge_in_graph(unsigned v1, unsigned v2) const;
    // TODO think of a better name for the function below
    bool connected_to_face(const Face& face, unsigned vertex, std::optional<unsigned> exclude) const;
    int get_sign(const Face& face, Face::iterator it, unsigned new_value) const;

    Face copy_with_replacement(const Face& orig, Face::iterator it, unsigned new_value) const {
        Face new_face;
        
        // Copy elements from the original set to the new set, excluding the element pointed to by the iterator
        for (auto iter = orig.begin(); iter != orig.end(); ++iter) {
            if (iter != it) {
                new_face.insert(*iter);
            }
        }
        
        // Insert the new value into the new set
        new_face.insert(new_value);
        return new_face;
    }


public:
    Complex(const Graph& g) : graph(g) {}
    ~Complex();
 
    bool includes_face(const Face& face) const;
    unsigned number_of_up_neighbours(const Face& face) const;
    std::vector<std::pair<Face,int>> get_neighbours(const Face& face, unsigned n) const;
    Face sample_from_complex(unsigned k) const;
};