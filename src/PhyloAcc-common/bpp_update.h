#ifndef PHYLOACC_COMMON_BPP_UPDATE_H
#define PHYLOACC_COMMON_BPP_UPDATE_H

#include <armadillo>

#include <vector>

namespace phyloacc
{
void MarkChangedZAncestors(const std::vector<int>& changedZ,
                          int node_count,
                          int root_sentinel,
                          const int* parent,
                          bool start_from_parent,
                          bool clear_lambda,
                          std::vector<bool>& visited,
                          std::vector<std::vector<arma::vec> >* lambda);
}

#endif
