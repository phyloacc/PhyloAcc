#include "bpp_init.h"

namespace phyloacc
{
void InitializeTreeArrays(const PhyloTree& tree,
                          int node_count,
                          int parent_default,
                          int (*&children)[2],
                          int*& parent,
                          double*& distances)
{
    children = new int[node_count][2];
    parent = new int[node_count];
    distances = new double[node_count];

    for (int node = 0; node < node_count; ++node)
    {
        distances[node] = tree.distances[node];
        children[node][0] = -1;
        children[node][1] = -1;
        parent[node] = parent_default;
    }

    for (int node = 0; node < node_count; ++node)
    {
        int child_slot = -1;
        for (int child = 0; child < node_count; ++child)
        {
            if (tree.dag[node][child])
            {
                ++child_slot;
                children[node][child_slot] = child;
                parent[child] = node;
            }
        }
    }
}
}
