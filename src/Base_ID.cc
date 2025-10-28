#include "Base_ID.h"
#include <TObjArray.h>
#include <TLeaf.h>
#include <iostream>

Base_ID::Base_ID(TTree *tree, const char *className) : fTree(tree) {}

void Base_ID::PopulateLeafMap()
{
    if (!fTree)
    {
        std::cerr << "fTree is null in PopulateLeafMap." << std::endl;
        return;
    }

    TObjArray *leaves = fTree->GetListOfLeaves();
    if (!leaves)
    {
        std::cerr << "No leaves found in the tree." << std::endl;
        return;
    }

    for (int i = 0; i < leaves->GetEntries(); ++i)
    {
        TLeaf *leaf = (TLeaf *)leaves->At(i);
        leafMap[leaf->GetName()] = leaf;
    }
}
