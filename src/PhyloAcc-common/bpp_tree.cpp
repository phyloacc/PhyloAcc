#include "bpp_tree.h"

namespace phyloacc
{
void CollectUpperTreeNodes(int root,
                           const std::vector<int>& child_nodes,
                           int* parent,
                           std::set<int>& visited_nodes)
{
    for (std::vector<int>::const_iterator it = child_nodes.begin(); it != child_nodes.end(); ++it)
    {
        int node = *it;
        while (node != root)
        {
            visited_nodes.insert(node);
            node = parent[node];
        }
    }

    visited_nodes.insert(root);
}

void CollectSubtreeNodes(int root,
                         int (*children)[2],
                         std::vector<int>& visited_nodes)
{
    int node = root;
    if (children[node][0] != -1)
    {
        for (int child_index = 0; child_index < 2; child_index++)
            CollectSubtreeNodes(children[node][child_index], children, visited_nodes);
    }
    visited_nodes.push_back(node);
}

void CollectSubtreeNodesUntilChildren(int root,
                                      int (*children)[2],
                                      std::set<int>& stop_children,
                                      std::vector<int>& visited_nodes)
{
    int node = root;
    if (stop_children.find(node) != stop_children.end())
    {
        visited_nodes.push_back(node);
        return;
    }

    if (children[node][0] != -1)
    {
        for (int child_index = 0; child_index < 2; child_index++)
            CollectSubtreeNodesUntilChildren(children[node][child_index], children, stop_children, visited_nodes);
    }
    visited_nodes.push_back(node);
}
}
