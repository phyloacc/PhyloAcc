// This code, downloaded from Internet, is used to
// parse the phylogenetic tree in Newick format.

#include "newick2.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
//#include "newick.h"

#include "utils.h"

#define APPEND_LEN	256

using namespace std;
using namespace arma;

typedef struct newick2_child
{
    struct newick2_node *node;
    struct newick2_child *next;
} newick2_child;

typedef struct newick2_node
{
    int id;
    char *taxon;
    char *seq;
    double dist;
    double theta = 0; // 2\mu* N
    int childNum;
    struct newick2_child *child;
    struct newick2_node *parent;
} newick2_node;

newick2_node* parseTree2(char *str); // with theta
newick2_node* parseTree02(char *str); // no theta


typedef struct seqMem2
{
    void *pos2;
    struct seqMem2 *next;
} seqMem2;

seqMem2 *start2;
seqMem2 *current2;

void seqMemInit2()
{
    start2 = NULL;
}

void* seqMalloc2(int size)
{
    if (start2 == NULL)
    {
        start2 = (seqMem2*)malloc(sizeof(seqMem2));
        memset(start2, '\0', sizeof(seqMem2));
        current2 = start2;
    }
    else
    {
        current2->next = (seqMem2*)malloc(sizeof(seqMem2));
        memset(current2->next, '\0', sizeof(seqMem2));
        current2 = current2->next;
    }
    current2->pos2 = malloc(size);
    memset(current2->pos2, '\0', size);
    return(current2->pos2);
}

void seqFreeAll2()
{
    while (start2 != NULL)
    {
        current2 = start2->next;
        free(start2->pos2);
        free(start2);
        start2 = current2;
    }

    start2 = NULL;
}

void seqFree2(void* pos2)
{
    seqMem2 *node, *prenode;
    node = start2;
    prenode = start2;
    while (node != NULL)
    {
        if (node->pos2 == pos2)
        {
            free(node->pos2);
            if (node == start2)
            {
                start2 = node->next;
            }
            else if (node->next == NULL)
            {
                current2 = prenode;
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

void inputString2(char *input, char **ppcStr, int *iLen, int *iMaxLen)
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
    temp = (char*)seqMalloc2(*iMaxLen);
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
        seqFree2(*ppcStr);
    }
    *ppcStr = temp;
}


newick2_node* parseTree02(char *str)
{
    newick2_node *node;
    newick2_child *child;
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
        node = (newick2_node*)seqMalloc2(sizeof(newick2_node));
        if (pcColon == NULL)
        {
            // Taxon only
            node->taxon = (char*)seqMalloc2(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
        }
        else
        {
            // Taxon
            *pcColon = '\0';
            node->taxon = (char*)seqMalloc2(strlen(pcStart) + 1);
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
        node = (newick2_node*)seqMalloc2(sizeof(newick2_node));
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
                        node->child = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                        node->childNum = 1;
                        child = node->child;
                    }
                    else
                    {
                        child->next = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                        node->childNum++;
                        child = child->next;
                    }
                    child->node = parseTree2(pcStart);
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
                        node->child = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                        node->childNum = 1;
                        child = node->child;
                    }
                    else
                    {
                        child->next = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                        node->childNum++;
                        child = child->next;
                    }
                    child->node = parseTree2(pcStart);
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
            node->taxon = (char*)seqMalloc2(strlen(pcStart) + 1);
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
newick2_node* parseTree2(char *str)
{
    newick2_node *node;
    newick2_child *child;
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
        node = (newick2_node*)seqMalloc2(sizeof(newick2_node));
        if (pcColon == NULL)
        {
            // Taxon only
            node->taxon = (char*)seqMalloc2(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
        }
        else
        {
            // Taxon
            *pcColon = '\0';
            node->taxon = (char*)seqMalloc2(strlen(pcStart) + 1);
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
        node = (newick2_node*)seqMalloc2(sizeof(newick2_node));
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
                    node->child = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                    node->childNum = 1;
                    child = node->child;
                }
                else
                {
                    child->next = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                    node->childNum++;
                    child = child->next;
                }
                child->node = parseTree2(pcStart);
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
                    node->child = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                    node->childNum = 1;
                    child = node->child;
                }
                else
                {
                    child->next = (newick2_child*)seqMalloc2(sizeof(newick2_child));
                    node->childNum++;
                    child = child->next;
                }
                child->node = parseTree2(pcStart);
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
            //while (*pcCurrent != '#' && *pcCurrent != ';')
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            //node->dist = (float)atof(pcStart);
            node->dist = (double)atof(pcStart);
            *pcCurrent = cTemp;
            pcCurrent++;
            
            /*
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            //node->theta = (float)atof(pcStart);
            node->theta = (double)atof(pcStart);
            *pcCurrent = cTemp;
            */
        }
        else if (*pcCurrent != ';' && *pcCurrent != '\0')
        {
            // Find ':' to retrieve distance, if any.
            // At this time *pcCurrent should equal to ')'
            pcStart = pcCurrent;
            while (*pcCurrent != ':' && *pcCurrent != '#' && *pcCurrent != '\0')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            node->taxon = (char*)seqMalloc2(strlen(pcStart) + 1);
            memcpy(node->taxon, pcStart, strlen(pcStart));
            node->taxon[strlen(pcStart)] = '\0';
            *pcCurrent = cTemp;
            pcCurrent++;
            
            pcStart = pcCurrent;
            //while (*pcCurrent != '#' && *pcCurrent != ';')
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            //node->dist = (float)atof(pcStart);
            node->dist = (double)atof(pcStart);
            *pcCurrent = cTemp;
            pcCurrent++;
            
            /*
            pcStart = pcCurrent;
            while (*pcCurrent != '\0' && *pcCurrent != ';')
            {
                pcCurrent++;
            }
            cTemp = *pcCurrent;
            *pcCurrent = '\0';
            //node->theta = (float)atof(pcStart);
            node->theta = (double)atof(pcStart);
            *pcCurrent = cTemp;
            */
        }
    }
    return node;
}

void printTree(newick2_node *root)
{
	newick2_child *child;
	if (root->childNum == 0)
	{
		printf("%s:%0.6f", root->taxon, root->dist);
	}
	else
	{
		child = root->child;
		printf("(");
		while (child != NULL)
		{
			printTree(child->node);
			if (child->next != NULL)
			{
				printf(",");
			}
			child = child->next;
		}
		if (root->taxon != NULL)
		{
			printf(")%s:%0.6f", root->taxon, root->dist);
		}
		else
		{
			printf("):%0.6f", root->dist);
		}
	}
}

void TravelTree1(newick2_node *root, int &S)
{
    newick2_child *child = root->child;

    while (child != NULL)
    {
        TravelTree1(child->node, S);
        child = child->next;
    }

// if leaf
    if (root->child == NULL)
        S = S + 1;
}

int cur2_leaf_id = -1;
int cur2_branch_id = -1;

void TravelTree2_theta(newick2_node *root, PhyloTree_theta &phylo_tree, map<string, int> & speciesname)
{
    newick2_child *child = root->child;

    while (child != NULL)
    {
        TravelTree2_theta(child->node, phylo_tree, speciesname);
        child = child->next;
    }

    if (root->child == NULL)
    {
        if(speciesname.size() == 0)
        {
            root->id = ++cur2_leaf_id;
        }else{
            root->id = speciesname[root->taxon];
        }
    }
    else
    {
        root->id = ++cur2_branch_id + phylo_tree.S;
    }
}

void TravelTree3_theta(newick2_node *root2, PhyloTree_theta &phylo_tree)
{
    newick2_child *child = root2->child;
    while (child != NULL)
    {   
        // printf("Assigning dag[%d][%d]\n", root2->id, child->node->id);
        phylo_tree.dag[root2->id][child->node->id] = true;

        TravelTree3_theta(child->node, phylo_tree);
        child = child->next;
    }

    phylo_tree.distances[root2->id] = root2->dist;
    if( root2->taxon!=NULL) phylo_tree.nodes_names[root2->id] = root2->taxon;
    if (root2->id < phylo_tree.S)
    {
        phylo_tree.species_names[root2->id] = root2->taxon;
    }

}


PhyloTree_theta LoadPhyloTree_theta(string tree_path)
{
    cout <<"Loading phylogenetic tree in coalescent unit from "<<tree_path<<"... "<<endl;
    ifstream in_prof(tree_path.c_str());
    if(!in_prof){
        cerr<<"(Error. Cannot open tree2 input file: " << tree_path << ")" << endl;
        exit(1);
    }

    cur2_leaf_id = -1;
    cur2_branch_id = -1;

    PhyloTree_theta tree2;
    newick2_node *root2 = NULL;
    string linestr;
    std::getline(in_prof, linestr);
    linestr = strutils::trim(linestr);
    vector<string> line_splits = strutils::split(linestr, ' ');
    root2 = parseTree2((char*)strutils::trim(line_splits[0]).c_str());

    int S2=0;
    TravelTree1(root2,S2);
    tree2.S=S2;
    int N2 = 2*tree2.S-1;
    map<string, int> speciesname = map<string, int>();
    TravelTree2_theta(root2, tree2, speciesname);
    // load the tree structure and species names
    tree2.species_names = vector<string>(S2);
    tree2.nodes_names = vector<string>(N2);
    tree2.distances = vector<double>(N2);
    tree2.dag = vector< vector<bool> >(N2, vector<bool>(N2, false));
    TravelTree3_theta(root2, tree2);

    return tree2;
}

