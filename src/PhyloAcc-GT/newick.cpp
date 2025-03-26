// This code, downloaded from Internet, is used to
// parse the phylogenetic tree in Newick format.

#include "newick.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "utils.h"

#define APPEND_LEN	256

using namespace std;
using namespace arma;

typedef struct newick_child
{
    struct newick_node *node;
    struct newick_child *next;
} newick_child;

typedef struct newick_node
{
    int id;
    char *taxon;
    char *seq;
    double dist;
    double theta = 0; // 2\mu* N
    int childNum;
    struct newick_child *child;
    struct newick_node *parent;
} newick_node;

newick_node* parseTree(char *str); // with theta
newick_node* parseTree0(char *str); // no theta

typedef struct seqMem
{
    void *pos;
    struct seqMem *next;
} seqMem;

seqMem *start;
seqMem *current;

void seqMemInit()
{
    start = NULL;
}

void* seqMalloc(int size)
{
    if (start == NULL)
    {
        start = (seqMem*)malloc(sizeof(seqMem));
        memset(start, '\0', sizeof(seqMem));
        current = start;
    }
    else
    {
        current->next = (seqMem*)malloc(sizeof(seqMem));
        memset(current->next, '\0', sizeof(seqMem));
        current = current->next;
    }
    current->pos = malloc(size);
    memset(current->pos, '\0', size);
    return(current->pos);
}

void seqFreeAll()
{
    while (start != NULL)
    {
        current = start->next;
        free(start->pos);
        free(start);
        start = current;
    }

    start = NULL;
}

void seqFree(void* pos)
{
    seqMem *node, *prenode;
    node = start;
    prenode = start;
    while (node != NULL)
    {
        if (node->pos == pos)
        {
            free(node->pos);
            if (node == start)
            {
                start = node->next;
            }
            else if (node->next == NULL)
            {
                current = prenode;
                prenode->next = NULL;
            }
            else
            {
                prenode->next = node->next;
            }
            free(node);
            break;
        }

        prenode = node;
        node = node->next;
    }
}

void inputString(char *input, char **ppcStr, int *iLen, int *iMaxLen)
{
    int inputLen;
    char *temp;
    inputLen = strlen(input);
    if (inputLen == 0)
    {
        return;
    }
    while (*iMaxLen < (*iLen + inputLen) + 1)
    {
        *iMaxLen = *iMaxLen + APPEND_LEN;
    }
    temp = (char*)seqMalloc(*iMaxLen);
    if (*ppcStr == NULL)
    {
        memcpy(temp, input, inputLen);
    }
    else
    {
        memcpy(temp, *ppcStr, *iLen);
        strcat(temp, input);
    }
    *iLen = *iLen + inputLen;
    if (*ppcStr != NULL)
    {
        seqFree(*ppcStr);
    }
    *ppcStr = temp;
}


newick_node* parseTree0(char *str)
{
    newick_node *node;
    newick_child *child;
    char *pcCurrent;
    char *pcStart;
    char *pcColon = NULL;
    char cTemp;
    int iCount;
    
    pcStart = str;
    
    if (*pcStart != '(')
    {
        // Leaf node. Separate taxon name from distance. If distance not exist then take care of taxon name only
        pcCurrent = str;
        while (*pcCurrent != '\0')
        {
            if (*pcCurrent == ':')
            {
                pcColon = pcCurrent;
            }
            pcCurrent++;
        }
        node = (newick_node*)seqMalloc(sizeof(newick_node));
        if (pcColon == NULL)
        {
            // Taxon only
            node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
        }
        else
        {
            // Taxon
            *pcColon = '\0';
            node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
            *pcColon = ':';
            // Distance
            pcColon++;
            //node->dist = (float)atof(pcColon);
            node->dist = (double)atof(pcColon);
        }
        node->childNum = 0;
    }
    else
    {
        // Create node
        node = (newick_node*)seqMalloc(sizeof(newick_node));
        child = NULL;
        // Search for all child nodes
        // Find all ',' until corresponding ')' is encountered
        iCount = 0;
        pcStart++;
        pcCurrent = pcStart;
        while (iCount >= 0)
        {
            switch (*pcCurrent)
            {
                case '(':
                    // Find corresponding ')' by counting
                    pcStart = pcCurrent;
                    pcCurrent++;
                    iCount++;
                    while (iCount > 0)
                    {
                        if (*pcCurrent == '(')
                        {
                            iCount++;
                        }
                        else if (*pcCurrent == ')')
                        {
                            iCount--;
                        }
                        pcCurrent++;
                    }
                    while (*pcCurrent != ',' && *pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    cTemp = *pcCurrent;
                    *pcCurrent = '\0';
                    // Create a child node
                    if (child == NULL)
                    {
                        node->child = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum = 1;
                        child = node->child;
                    }
                    else
                    {
                        child->next = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum++;
                        child = child->next;
                    }
                    child->node = parseTree(pcStart);
                    *pcCurrent = cTemp;
                    if (*pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    break;
                    
                case ')':
                    // End of tihs tree. Go to next part to retrieve distance
                    iCount--;
                    break;
                    
                case ',':
                    // Impossible separation since according to the algorithm, this symbol will never encountered.
                    // Currently don't handle this and don't create any node
                    break;
                    
                default:
                    // leaf node encountered
                    pcStart = pcCurrent;
                    while (*pcCurrent != ',' && *pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    cTemp = *pcCurrent;
                    *pcCurrent = '\0';
                    // Create a child node
                    if (child == NULL)
                    {
                        node->child = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum = 1;
                        child = node->child;
                    }
                    else
                    {
                        child->next = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum++;
                        child = child->next;
                    }
                    child->node = parseTree(pcStart);
                    *pcCurrent = cTemp;
                    if (*pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    break;
            }
        }
        
        // If start at ':', then the internal node has no name.
        pcCurrent++;
        if (*pcCurrent == ':')
        {
            pcStart = pcCurrent + 1;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            //node->dist = (float)atof(pcStart);
            node->dist = (double)atof(pcStart);
            *pcCurrent = cTemp;
        }
        else if (*pcCurrent != ';' && *pcCurrent != '\0')
        {
            // Find ':' to retrieve distance, if any.
            // At this time *pcCurrent should equal to ')'
            pcStart = pcCurrent;
            while (*pcCurrent != ':' && *pcCurrent!=';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
            *pcCurrent = cTemp;
            pcCurrent++;
            pcStart = pcCurrent;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            //node->dist = (float)atof(pcStart);
            node->dist = (double)atof(pcStart);
            *pcCurrent = cTemp;
        }
    }
    
    return node;
}
newick_node* parseTree(char *str)
{
    newick_node *node;
    newick_child *child;
    char *pcCurrent;
    char *pcStart;
    char *pcColon = NULL;
    char cTemp;
    int iCount;

    pcStart = str;

    if (*pcStart != '(')
    {
        // Leaf node. Separate taxon name from distance.
        pcCurrent = str;
        while (*pcCurrent != '\0')
        {
            if (*pcCurrent == ':')
            {
                pcColon = pcCurrent;
            }
            pcCurrent++;
        }
        node = (newick_node*)seqMalloc(sizeof(newick_node));
        if (pcColon == NULL)
        {
            // Taxon only
            node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
        }
        else
        {
            // Taxon and distance
            *pcColon = '\0';
            node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
            *pcColon = ':';
            pcColon++;
            node->dist = (double)atof(pcColon);
        }
        node->childNum = 0;
    }
    else
    {
        // Internal node.
        node = (newick_node*)seqMalloc(sizeof(newick_node));
        child = NULL;
        iCount = 0;
        pcStart++; // Skip the '('
        pcCurrent = pcStart;
        while (iCount >= 0)
        {
            switch (*pcCurrent)
            {
                case '(':
                    // Found a subtree, so count parentheses.
                    pcStart = pcCurrent;
                    pcCurrent++;
                    iCount++;
                    while (iCount > 0)
                    {
                        if (*pcCurrent == '(')
                        {
                            iCount++;
                        }
                        else if (*pcCurrent == ')')
                        {
                            iCount--;
                        }
                        pcCurrent++;
                    }
                    // Move until ',' or ')' encountered.
                    while (*pcCurrent != ',' && *pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    cTemp = *pcCurrent;
                    *pcCurrent = '\0';
                    // Create a child node.
                    if (child == NULL)
                    {
                        node->child = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum = 1;
                        child = node->child;
                    }
                    else
                    {
                        child->next = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum++;
                        child = child->next;
                    }
                    child->node = parseTree(pcStart);
                    *pcCurrent = cTemp;
                    if (*pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    break;

                case ')':
                    // End of this subtree.
                    iCount--;
                    break;

                case ',':
                    // According to the algorithm, this case should not occur.
                    break;

                default:
                    // Leaf node encountered among the children.
                    pcStart = pcCurrent;
                    while (*pcCurrent != ',' && *pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    cTemp = *pcCurrent;
                    *pcCurrent = '\0';
                    // Create a child node.
                    if (child == NULL)
                    {
                        node->child = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum = 1;
                        child = node->child;
                    }
                    else
                    {
                        child->next = (newick_child*)seqMalloc(sizeof(newick_child));
                        node->childNum++;
                        child = child->next;
                    }
                    child->node = parseTree(pcStart);
                    *pcCurrent = cTemp;
                    if (*pcCurrent != ')')
                    {
                        pcCurrent++;
                    }
                    break;
            }
        }

        // Now we've finished processing children.
        // Move past the closing ')'
        pcCurrent++;

        // ***** Modification: Check for the tree terminator.
        if (*pcCurrent == ';')
            return node;

        // Next, check if the internal node has additional attributes.
        if (*pcCurrent == ':')
        {
            // Parse branch length.
            pcStart = pcCurrent + 1;
            while (*pcCurrent != '#' && *pcCurrent != ';' && *pcCurrent != '\0')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->dist = (double)atof(pcStart);
            *pcCurrent = cTemp;
            pcCurrent++;
            
            // After branch length, check again if we've reached the end.
            if (*pcCurrent == ';')
                return node;
            
            // Parse theta (if present).
            pcStart = pcCurrent;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->theta = (double)atof(pcStart);
            *pcCurrent = cTemp;
        }
        else if (*pcCurrent != ';' && *pcCurrent != '\0')
        {
            // Parse taxon name then branch length and theta.
            pcStart = pcCurrent;
            while (*pcCurrent != ':' && *pcCurrent != '#' && *pcCurrent != ';' && *pcCurrent != '\0')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
            *pcCurrent = cTemp;
            pcCurrent++;
            
            pcStart = pcCurrent;
            while (*pcCurrent != '#' && *pcCurrent != ';' && *pcCurrent != '\0')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->dist = (double)atof(pcStart);
            *pcCurrent = cTemp;
            pcCurrent++;
            
            pcStart = pcCurrent;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->theta = (double)atof(pcStart);
            *pcCurrent = cTemp;
        }
    }

    return node;
}


void printTree(newick_node *root, int isRoot)
{
    newick_child *child;
    if (root->childNum == 0)
    {
        // Leaf: Always print taxon and branch length.
        printf("%s:%0.6f", root->taxon, root->dist);
    }
    else
    {
        // Internal node: print children in parentheses.
        child = root->child;
        printf("(");
        while (child != NULL)
        {
            // For recursive calls, we're not the root.
            printTree(child->node, 0);
            if (child->next != NULL)
            {
                printf(",");
            }
            child = child->next;
        }
        printf(")");
        // Print taxon if present.
        if (root->taxon != NULL)
        {
            printf("%s", root->taxon);
        }
        // Only print branch length if this isn’t the root.
        if (!isRoot)
        {
            printf(":%0.6f", root->dist);
        }
    }
}



void TravelTree1(newick_node *root, int &S)
{
    newick_child *child = root->child;

    while (child != NULL)
    {
        TravelTree1(child->node, S);
        child = child->next;
    }

// if leaf
    if (root->child == NULL)
        S = S + 1;
}

int cur_leaf_id = -1;
int cur_branch_id = -1;

void TravelTree2(newick_node *root, PhyloTree &phylo_tree, map<string, int> & speciesname)
{
    newick_child *child = root->child;

    while (child != NULL)
    {
        TravelTree2(child->node, phylo_tree, speciesname);
        child = child->next;
    }

    if (root->child == NULL)
    {
        if(speciesname.size() == 0)
        {
            root->id = ++cur_leaf_id;
        }else{
            root->id = speciesname[root->taxon];
        }
    }
    else
    {
        root->id = ++cur_branch_id + phylo_tree.S;
    }
}

void TravelTree3(newick_node *root, PhyloTree &phylo_tree)
{
    newick_child *child = root->child;

    while (child != NULL)
    {
        phylo_tree.dag[root->id][child->node->id] = true;

        TravelTree3(child->node, phylo_tree);
        child = child->next;
    }

    phylo_tree.distances[root->id] = root->dist;

    if (root == nullptr || root->taxon == nullptr) {
        cout << "\n(Error. Some nodes in the tree aren't labeled. Please label all nodes and try again!)" << endl;
        exit(1);
    }
    // if the current node does not have a label, print an error and exit

    if( root->taxon!=NULL) phylo_tree.nodes_names[root->id] = root->taxon;
    phylo_tree.thetas[root->id] = root->theta;
    
    if (root->id < phylo_tree.S)
    {
        phylo_tree.species_names[root->id] = root->taxon;
        
    }
}

// load the phylogenetic tree
PhyloTree LoadPhyloTree(string phylo_tree_path)
{
   cout << "Loading phylogenetic tree from " << phylo_tree_path << "......" << endl;
   
    
    ifstream in_prof(phylo_tree_path.c_str());
    
    if (!in_prof)
    {
        cerr << "(Error. Cannot open the phylogenetic tree input file: " << phylo_tree_path << ")" << endl;
        exit(1);
    }
    
    // count the num of species, base pairs and load the profiles
    PhyloTree phylo_tree;
    phylo_tree.pi = zeros<vec>(4);
    phylo_tree.subs_rate = zeros<mat>(4,4);
    
    newick_node *root = NULL;
    string linestr;
    
    while(!in_prof.eof())
    {
        std::getline(in_prof, linestr);
        if(linestr=="") break;
        linestr = strutils::trim(linestr);
        vector<string> line_splits = strutils::split(linestr, ' ');
        if(!strcmp(line_splits[0].c_str(), "BACKGROUND:")){
            
            int ind=0;
            for(std::size_t g=1; g<line_splits.size(); g++)
            {
                if(line_splits[g]=="") continue;
                phylo_tree.pi[ind] = atof(line_splits[g].c_str());
                ind++;
            }
            
        }
        else if(!strcmp(line_splits[0].c_str(), "RATE_MAT:")) //Han: read in Q
        {
            for(int i = 0;i<4;i++){
                std::getline(in_prof, linestr);
                vector<string> tmp = strutils::split(strutils::trim(linestr),' ');
                int ind=0;
                for(std::size_t g=0; g < tmp.size(); g++)
                {
                    if(tmp[g]=="") continue;
                    phylo_tree.subs_rate(i,ind) = atof(tmp[g].c_str());
                    ind++;
                }
            }
            
        }
        else if (!strcmp(line_splits[0].c_str(), "TREE:")){
            root = parseTree((char*)strutils::trim(line_splits[1]).c_str());
        }
    }

    printTree(root, 1);
    printf(";\n");

    // count the total num of species and label each node and
    int S = 0;
    TravelTree1(root, S); //update S to be total num of leaves.
    phylo_tree.S = S;
    int N = 2 * phylo_tree.S - 1;
    map<string, int> speciesname = map<string, int>();
    TravelTree2(root, phylo_tree, speciesname);

    // load the tree structure and species names
    phylo_tree.species_names = vector<string>(S);
    phylo_tree.nodes_names = vector<string>(N);
    phylo_tree.distances = vector<double>(N);
    phylo_tree.thetas = vector<double>(N);
    phylo_tree.dag = vector< vector<bool> >(N, vector<bool>(N, false));
    TravelTree3(root, phylo_tree);
    return phylo_tree;
}


// load the phylogenetic tree
PhyloTree LoadPhyloTree_text(string text, map<string, int>& speciesname)
{
    PhyloTree phylo_tree;
    
    newick_node *root = NULL;
    string linestr;
    
    root = parseTree0((char*)text.c_str());
     
    // label each node (root->id)
    int S = speciesname.size();
    //TravelTree1(root, S);
    phylo_tree.S = S;
    int N = 2 * phylo_tree.S - 1;
    cur_leaf_id = -1;
    cur_branch_id = -1;
    
    TravelTree2(root, phylo_tree, speciesname);
    
    // load the tree structure and species names
    phylo_tree.species_names = vector<string>(S);
    phylo_tree.nodes_names = vector<string>(N);
    phylo_tree.distances = vector<double>(N);
    phylo_tree.thetas = vector<double>(N);
    phylo_tree.dag = vector< vector<bool> >(N, vector<bool>(N, false));
    TravelTree3(root, phylo_tree);
    
    return phylo_tree;
}
