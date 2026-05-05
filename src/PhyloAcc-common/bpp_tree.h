#ifndef PHYLOACC_BPP_TREE_H
#define PHYLOACC_BPP_TREE_H

#include <set>
#include <vector>

namespace phyloacc
{
void CollectUpperTreeNodes(int root,
                           const std::vector<int>& child_nodes,
                           int* parent,
                           std::set<int>& visited_nodes);

void CollectSubtreeNodes(int root,
                         int (*children)[2],
                         std::vector<int>& visited_nodes);

void CollectSubtreeNodesUntilChildren(int root,
                                      int (*children)[2],
                                      std::set<int>& stop_children,
                                      std::vector<int>& visited_nodes);
}

#endif
