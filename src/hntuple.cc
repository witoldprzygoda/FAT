/**
 * @file hntuple.cc
 * @brief Implementation file for the HNtuple class.
 * 
 * This file contains the implementation of the HNtuple class, which is used to create and fill ROOT NTuples.
 * The class provides constructors, member functions, and operators for creating and manipulating NTuples.
 * The implementation includes lazy construction of the NTuple, variable mapping, and filling of the NTuple.
 * 
 * @author: Witold Przygoda (witold.przygoda@uj.edu.pl)
 * @date Last modified on: 1 March 2009
 * @date Last modified on: 20 November 2023
 */


#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "hntuple.h"


ClassImp(HNtuple)
// ---------------------------------------------------------------------------------
// Dummy constructor
// ---------------------------------------------------------------------------------
HNtuple::HNtuple() : cname(nullptr), ctitle(nullptr), cbufsize(0)
{
}

// ---------------------------------------------------------------------------------
// This is ntuple constructor with "lazy"-delayed construction
// ntuple is created only after first fill attempt based on 
// variables which have been set 
// ---------------------------------------------------------------------------------
HNtuple::HNtuple(const char* name, const char* title, Int_t bufsize) : 
	cname(name), ctitle(title == nullptr ? name : title), cbufsize(bufsize)
{
}

// ---------------------------------------------------------------------------------
// Basic ntuple constructor with exactly the same parameters as in the case of ROOT TNtuple
// ---------------------------------------------------------------------------------
HNtuple::HNtuple(const char* name, const char* title, const char* varlist, Int_t bufsize) :
	cname(name), ctitle(title), cbufsize(bufsize)
{
   ptrNt = std::make_unique<TNtuple>(name, title, varlist, bufsize);
   setMap(varlist, isNtuple);
}

// ---------------------------------------------------------------------------------
void HNtuple::setMap(const std::string& varList, Bool_t& kPair)
{
   Int_t i = 0;
   std::string sName; // empty string
   for (const auto& c : varList)
   {
      if (c != ':')
      {
         sName += c;
      }
      else
      {
         // create a new pair ntuple variable - value
         vKeyOrder.insert(make_pair(sName, i++));
         if (kPair == kFALSE)
            vKeyValue.insert(make_pair(sName, 0.));
         sName.clear();
      }
   }
   // create the last new pair ntuple variable - value
   vKeyOrder.insert(make_pair(sName, i++));
   if (kPair == kFALSE)
      vKeyValue.insert(make_pair(sName, 0.));
   sName.clear();

   // create Float_t array based on variables string and number of ':' separators
   varArrayN = 1 + count_if(varList.begin(), varList.end(), [](char c){ return c == ':'; });
   varArray.resize(varArrayN);
   kPair = kTRUE;
}

// ---------------------------------------------------------------------------------
HNtuple::~HNtuple() 
{
   if (outFile)
   {
      outFile->cd();
      outFile->Close();
   }
}

// ---------------------------------------------------------------------------------
Int_t HNtuple::Write(const char* name, Int_t option, Int_t bufsize) 
{
   if (isNtuple==kTRUE)
   {
      return ptrNt->Write(name, option, bufsize); 
   }
   return -1; // error if NTuple was not booked
}

// ---------------------------------------------------------------------------------
void HNtuple::SetDirectory(TDirectory* dir)
{
   return ptrNt->SetDirectory(dir);
}

// ---------------------------------------------------------------------------------
Float_t& HNtuple::operator[](const std::string& key)
{
   if (isNtuple)
   {
      auto mIter = vKeyValue.find(key);
      if (mIter != vKeyValue.end())
      {
         return mIter->second;
      }
      throw std::invalid_argument("A variable: \"" + key + "\" tried to be assigned after HNtuple was booked, that is, e.g after the first \"fill()\" call.");
   }

return vKeyValue[key]; 
}

// ---------------------------------------------------------------------------------
const Float_t &HNtuple::operator[](const std::string &key) const
{
   const auto mcIter = vKeyValue.find(key);
   if (mcIter != vKeyValue.end())
   {
      return mcIter->second;
   }
   throw std::invalid_argument("A variable: \"" + key + "\" tried to be assigned after HNtuple was booked, that is, e.g after the first \"fill()\" call.");
}

// ---------------------------------------------------------------------------------
Int_t HNtuple::fill()
{
// This function is similar to Fill(...) from NTuple class besides
// the fact that is has no parameter and is small "f" :)

   if (isNtuple==kTRUE)
   {
      std::fill(varArray.begin(), varArray.end(), 0.0);
      for (auto& pair : vKeyValue) {
         varArray[vKeyOrder[pair.first]] = pair.second;
         // reset of map array
         pair.second = 0.0;
      }
   }
   else
   {
      // ntuple not booked yet, we create it here based on variables 
      // set with the function setVal or operator[]
      
      std::string vList;
      for (const auto &pair : vKeyValue)
      {
         vList += pair.first + ":";
      }
      vList.erase(vList.find_last_of(":"),1);

      //-------- here a part of NTuple Ctor
      if (outFile)
      {
         outFile->cd();
      }
      else
      {
         throw std::runtime_error("NTuple booked but not attached to any file. Forgot to call: void setFile(TFile *ptrF) for this ntuple?");
      }
      ptrNt = std::make_unique<TNtuple>(cname, ctitle, vList.c_str(), cbufsize);
      isNtuple = kTRUE;
      setMap(vList, isNtuple);
	 //-------- fill
      std::fill(varArray.begin(), varArray.end(), 0.0f);
      for (auto& pair : vKeyValue) {
         varArray[vKeyOrder[pair.first]] = pair.second;
         // reset of map array
         pair.second = 0.;
      }
   }

   // filling the ROOT ntuple
   return ptrNt->Fill(varArray.data());
}
// ****************************************************************************



