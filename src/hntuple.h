/**
 * @file hntuple.h
 * @brief Header file for the HNtuple class.
 * 
 * This file contains the declaration of the HNtuple class, which is used to improve the behavior of TNtuple.
 * The class provides a convenient way to assign values to variables without the need to create float arrays.
 * It inherits from TObject and stores the output as the original NTuple type.
 * 
 * @author: Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date Last modified: 1 March 2009
 * @date Last modified: 20 November 2023
 */


#ifndef HNTUPLE_H
#define HNTUPLE_H

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <string>
#include <utility>
#include <iterator>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

// ****************************************************************************
class HNtuple: public TObject
{
// Basic features:
// You do not need to create any float array of variables you want to fill with.
// You simply write something like: myHtuple["variable_name"] = variable_value
// and you can do it in any order you want (no more remembering of position of
// a variable in an array)

public:

   HNtuple(); // dummy constructor
   HNtuple(const char* name, const char* title = nullptr, Int_t bufsize = 32000); // late and probably most important ntuple c-tor
   HNtuple(const char* name, const char* title, const char* varlist, Int_t bufsize = 32000); // basic constructor
   virtual ~HNtuple(); // destructor closing the ntuple file

   Int_t Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0); // virtual function for a polymorphic NTulple Write call
   Float_t& operator[](const std::string& key); // the way of assigning values for variables
   const Float_t& operator[](const std::string& key) const; // the way of variable value reading
   Int_t fill(); // fills tuple, if not defined, constructs it first
   void SetDirectory(TDirectory* dir); 
   const char* getName() const { return ptrNt->GetName(); }
   const char* getTitle() const { return ptrNt->GetTitle(); }
   void setFile(TFile *ptrF) { outFile = ptrF; }

   // Query methods for ntuple state
   Bool_t isFrozen() const { return isNtuple; }
   Int_t getNVariables() const { return varArrayN; }
   std::vector<std::string> getVariableNames() const;
   Bool_t hasVariable(const std::string& key) const;

   // Diagnostics and debugging
   void printStructure(std::ostream& os = std::cout) const;
   std::string getStructureString() const;

protected:

   TFile *outFile {nullptr};

   HNtuple& operator=(const HNtuple& src) const; 
   void setMap(const std::string& varList, Bool_t& kPair); // creates a map from variable string

   const char* cname; //!
   const char* ctitle; //!
   Int_t cbufsize; //!
   std::unique_ptr<TNtuple> ptrNt; //!

   Bool_t isNtuple {kFALSE}; //! kTRUE if ntuple is defined
   Int_t varArrayN {0}; //! number of ntuple variables
   Long64_t fillCount {0}; //! number of times fill() has been called
   std::vector<Float_t> varArray; //! table of values for ntuple to be filled with
   std::map<std::string, Float_t> vKeyValue; //! pair of a variable name and a value
   std::map<std::string, Int_t> vKeyOrder; //! pair of a variable name and its position in ntuple

ClassDef(HNtuple, 0)
};
// ****************************************************************************


#endif /*!HNTUPLE_H*/

