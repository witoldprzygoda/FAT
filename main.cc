#include "PPip_ID.h"
#include <thread>
#include "data.h"
#include "src/datamanager.h"



int main() {
    DataManager manager("output10.root");
    
    PPip_ID reader("PPip_ID");
    reader.ProcessEntries();
}

