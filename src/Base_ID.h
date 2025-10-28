#ifndef BASE_ID_H
#define BASE_ID_H

#include <TTree.h>
#include <map>
#include <string>
#include <TLeaf.h>

class Base_ID
{
public:
    Base_ID(TTree *tree, const char *className);
    virtual ~Base_ID() {}
    void PopulateLeafMap();

protected:
    TTree *fTree;
    std::map<std::string, TLeaf *> leafMap;
};

#endif // BASE_ID_H
