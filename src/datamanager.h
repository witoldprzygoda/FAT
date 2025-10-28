#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "hntuple.h"
#include "manager.h"


class DataManager : public Manager {
public:
    DataManager(const char* name = "output.root"); // name of output file + init
    virtual ~DataManager(); // finalize objects 

    virtual void initData() override;
    virtual void finalizeData() override;

};


#endif // DATAMANAGER_H

