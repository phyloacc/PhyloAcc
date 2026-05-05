#include "bpp_update.h"

#include <cassert>

namespace phyloacc
{
void MarkChangedZAncestors(const std::vector<int>& changedZ,
                          int node_count,
                          int root_sentinel,
                          const int* parent,
                          bool start_from_parent,
                          bool clear_lambda,
                          std::vector<bool>& visited,
                          std::vector<std::vector<arma::vec> >* lambda)
{
    for (int i = 0; i < node_count; ++i)
        visited[i] = true;

    if (changedZ.empty())
        return;

    for (std::vector<int>::const_reverse_iterator it = changedZ.rbegin(); it != changedZ.rend(); ++it)
    {
        int j = *it;
        if (start_from_parent)
            j = parent[j];

        while (j != root_sentinel)
        {
            if (!visited[j])
                break;
            visited[j] = false;
            if (clear_lambda)
            {
                assert(lambda != nullptr);
                for (std::size_t g = 0; g < lambda->size(); ++g)
                    (*lambda)[g][j].zeros();
            }

            j = parent[j];
            assert(j != -1);
        }
    }
}
}
