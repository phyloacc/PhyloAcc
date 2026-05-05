#ifndef PHYLOACC_BPP_INIT_H
#define PHYLOACC_BPP_INIT_H

#include "newick.h"

namespace phyloacc
{
void InitializeTreeArrays(const PhyloTree& tree,
                          int node_count,
                          int parent_default,
                          int (*&children)[2],
                          int*& parent,
                          double*& distances);
}

#endif
