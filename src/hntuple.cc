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
	cname_str(name), ctitle_str(title == nullptr ? name : title), cbufsize(bufsize)
{
   cname = cname_str.c_str();
   ctitle = ctitle_str.c_str();
}

// ---------------------------------------------------------------------------------
// Basic ntuple constructor with exactly the same parameters as in the case of ROOT TNtuple
// ---------------------------------------------------------------------------------
HNtuple::HNtuple(const char* name, const char* title, const char* varlist, Int_t bufsize) :
	cname_str(name), ctitle_str(title), cbufsize(bufsize)
{
   cname = cname_str.c_str();
   ctitle = ctitle_str.c_str();
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
      // FROZEN: Build informative error message
      std::ostringstream oss;
      oss << "\n╔════════════════════════════════════════════════════════════════╗\n";
      oss << "║  HNtuple ERROR: Cannot add new variable after freeze          ║\n";
      oss << "╠════════════════════════════════════════════════════════════════╣\n";
      oss << "║ Attempted to add: \"" << key << "\"\n";
      oss << "║ NTuple name:      \"" << cname << "\"\n";
      oss << "║ Fill count:       " << fillCount << " (frozen after fill #1)\n";
      oss << "║\n";
      oss << "║ The NTuple structure is FROZEN after the first fill() call.\n";
      oss << "║ All variables must be defined BEFORE the first fill().\n";
      oss << "║\n";
      oss << "║ Current NTuple structure (" << varArrayN << " variables):\n";
      oss << "║ ┌────────────────────────────────────────────────────────────┐\n";

      std::vector<std::pair<std::string, Int_t>> sorted_vars(vKeyOrder.begin(), vKeyOrder.end());
      std::sort(sorted_vars.begin(), sorted_vars.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });

      for (const auto& var : sorted_vars) {
         oss << "║ │ [" << var.second << "] " << var.first << "\n";
      }
      oss << "║ └────────────────────────────────────────────────────────────┘\n";
      oss << "║\n";
      oss << "║ SOLUTION:\n";
      oss << "║   Add 'r[\"" << key << "\"] = value;' BEFORE the first fill() call,\n";
      oss << "║   or check for typos in the variable name.\n";
      oss << "╚════════════════════════════════════════════════════════════════╝\n";

      throw std::runtime_error(oss.str());
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

   // Variable not found - build informative error
   std::ostringstream oss;
   oss << "\nHNtuple ERROR: Variable \"" << key << "\" not found in ntuple \"" << cname << "\".\n";
   if (isNtuple) {
      oss << "NTuple is frozen (fill count: " << fillCount << "). Available variables:\n";
      std::vector<std::pair<std::string, Int_t>> sorted_vars(vKeyOrder.begin(), vKeyOrder.end());
      std::sort(sorted_vars.begin(), sorted_vars.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
      for (const auto& var : sorted_vars) {
         oss << "  [" << var.second << "] " << var.first << "\n";
      }
   }
   throw std::runtime_error(oss.str());
}

// ---------------------------------------------------------------------------------
Int_t HNtuple::fill()
{
// This function is similar to Fill(...) from NTuple class besides
// the fact that is has no parameter and is small "f" :)

   if (isNtuple==kTRUE)
   {
      // Already frozen - just fill with current values
      std::fill(varArray.begin(), varArray.end(), 0.0);
      for (auto& pair : vKeyValue) {
         varArray[vKeyOrder[pair.first]] = pair.second;
         // reset of map array
         pair.second = 0.0;
      }
   }
   else
   {
      // First fill() - create and freeze the ntuple structure

      // Validation: check if any variables were set
      if (vKeyValue.empty()) {
         throw std::runtime_error("HNtuple ERROR: Attempting to fill() without setting any variables!\n"
                                  "NTuple \"" + std::string(cname) + "\" has no variables defined.\n"
                                  "Use: r[\"variable_name\"] = value; before calling fill().");
      }

      // Build variable list
      std::string vList;
      for (const auto &pair : vKeyValue)
      {
         vList += pair.first + ":";
      }
      vList.erase(vList.find_last_of(":"),1);

      //-------- Create TNtuple
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

      // Print structure info on first fill
      std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
      std::cout << "║  HNtuple FROZEN: Structure locked after first fill()           ║\n";
      std::cout << "╠════════════════════════════════════════════════════════════════╣\n";
      
      // Name line: Total 66 - "║" (1) - " NTuple name: \"" (15) - name - "\"" (1) - " " - "║" (1)
      std::string name_str(cname);
      int name_padding = 66 - 1 - 15 - name_str.length() - 1 - 1;
      std::cout << "║ NTuple name: \"" << cname << "\"" 
                << std::string(name_padding, ' ') << "║\n";
      
      // Variables line: Total 66 - "║" (1) - " Variables:   " (14) - number - " " - "║" (1)
      std::string var_num_str = std::to_string(varArrayN);
      int var_padding = 66 - 1 - 14 - var_num_str.length() - 1;
      std::cout << "║ Variables:   " << varArrayN 
                << std::string(var_padding, ' ') << "║\n";
      
      // Inner box top: "║ ┌" + 60 dashes + "┐ ║" = 66
      std::cout << "║ ┌────────────────────────────────────────────────────────────┐ ║\n";

      std::vector<std::pair<std::string, Int_t>> sorted_vars(vKeyOrder.begin(), vKeyOrder.end());
      std::sort(sorted_vars.begin(), sorted_vars.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });

      for (const auto& var : sorted_vars) {
         // Each line: "║ │ " + content + padding + " │ ║" = 66
         // "║ │ " = 4 chars, " │ ║" = 4 chars, content area = 66 - 8 = 58 chars
         std::string var_line = "[" + std::to_string(var.second) + "] " + var.first;
         int padding = 58 - var_line.length();
         std::cout << "║ │ " << var_line << std::string(padding, ' ') << " │ ║\n";
      }
      std::cout << "║ └────────────────────────────────────────────────────────────┘ ║\n";
      std::cout << "║                                                                ║\n";
      std::cout << "║ This structure is now FROZEN. No new variables can be added.   ║\n";
      std::cout << "╚════════════════════════════════════════════════════════════════╝\n";

      //-------- Fill first entry
      std::fill(varArray.begin(), varArray.end(), 0.0f);
      for (auto& pair : vKeyValue) {
         varArray[vKeyOrder[pair.first]] = pair.second;
         // reset of map array
         pair.second = 0.;
      }
   }

   // Increment fill counter
   fillCount++;

   // filling the ROOT ntuple
   return ptrNt->Fill(varArray.data());
}

// ---------------------------------------------------------------------------------
std::vector<std::string> HNtuple::getVariableNames() const
{
   std::vector<std::pair<std::string, Int_t>> sorted_vars(vKeyOrder.begin(), vKeyOrder.end());
   std::sort(sorted_vars.begin(), sorted_vars.end(),
             [](const auto& a, const auto& b) { return a.second < b.second; });

   std::vector<std::string> names;
   names.reserve(sorted_vars.size());
   for (const auto& var : sorted_vars) {
      names.push_back(var.first);
   }
   return names;
}

// ---------------------------------------------------------------------------------
Bool_t HNtuple::hasVariable(const std::string& key) const
{
   return vKeyOrder.find(key) != vKeyOrder.end();
}

// ---------------------------------------------------------------------------------
void HNtuple::printStructure(std::ostream& os) const
{
   os << "╔════════════════════════════════════════════════════════════════╗\n";
   os << "║  HNtuple Structure                                             ║\n";
   os << "╠════════════════════════════════════════════════════════════════╣\n";
   os << "║ Name:        \"" << cname << "\"\n";
   os << "║ Title:       \"" << ctitle << "\"\n";
   os << "║ Status:      " << (isNtuple ? "FROZEN" : "UNFROZEN (can add variables)") << "\n";
   os << "║ Fill count:  " << fillCount << "\n";
   os << "║ Variables:   " << varArrayN << "\n";

   if (!vKeyOrder.empty()) {
      os << "║ ┌────────────────────────────────────────────────────────────┐\n";

      std::vector<std::pair<std::string, Int_t>> sorted_vars(vKeyOrder.begin(), vKeyOrder.end());
      std::sort(sorted_vars.begin(), sorted_vars.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });

      for (const auto& var : sorted_vars) {
         os << "║ │ [" << var.second << "] " << var.first;
         if (isNtuple) {
            auto it = vKeyValue.find(var.first);
            if (it != vKeyValue.end()) {
               os << " = " << it->second;
            }
         }
         os << "\n";
      }
      os << "║ └────────────────────────────────────────────────────────────┘\n";
   }
   os << "╚════════════════════════════════════════════════════════════════╝\n";
}

// ---------------------------------------------------------------------------------
std::string HNtuple::getStructureString() const
{
   std::ostringstream oss;
   printStructure(oss);
   return oss.str();
}

// ****************************************************************************
