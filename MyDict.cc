// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MyDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "src/hntuple.h"
#include "src/manager.h"
#include "src/histogram_registry.h"
#include "src/histogram_factory.h"
#include "src/histogram_builder.h"
#include "src/pparticle.h"
#include "src/boost_frame.h"
#include "src/ntuple_reader.h"
#include "src/cut_manager.h"
#include "src/analysis_config.h"
#include "src/progressbar.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_HNtuple(void *p = nullptr);
   static void *newArray_HNtuple(Long_t size, void *p);
   static void delete_HNtuple(void *p);
   static void deleteArray_HNtuple(void *p);
   static void destruct_HNtuple(void *p);
   static void streamer_HNtuple(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HNtuple*)
   {
      ::HNtuple *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HNtuple >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("HNtuple", ::HNtuple::Class_Version(), "src/hntuple.h", 30,
                  typeid(::HNtuple), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HNtuple::Dictionary, isa_proxy, 16,
                  sizeof(::HNtuple) );
      instance.SetNew(&new_HNtuple);
      instance.SetNewArray(&newArray_HNtuple);
      instance.SetDelete(&delete_HNtuple);
      instance.SetDeleteArray(&deleteArray_HNtuple);
      instance.SetDestructor(&destruct_HNtuple);
      instance.SetStreamerFunc(&streamer_HNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HNtuple*)
   {
      return GenerateInitInstanceLocal(static_cast<::HNtuple*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::HNtuple*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HNtuple::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *HNtuple::Class_Name()
{
   return "HNtuple";
}

//______________________________________________________________________________
const char *HNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int HNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HNtuple::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HNtuple*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class HNtuple.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HNtuple(void *p) {
      return  p ? new(p) ::HNtuple : new ::HNtuple;
   }
   static void *newArray_HNtuple(Long_t nElements, void *p) {
      return p ? new(p) ::HNtuple[nElements] : new ::HNtuple[nElements];
   }
   // Wrapper around operator delete
   static void delete_HNtuple(void *p) {
      delete (static_cast<::HNtuple*>(p));
   }
   static void deleteArray_HNtuple(void *p) {
      delete [] (static_cast<::HNtuple*>(p));
   }
   static void destruct_HNtuple(void *p) {
      typedef ::HNtuple current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_HNtuple(TBuffer &buf, void *obj) {
      ((::HNtuple*)obj)->::HNtuple::Streamer(buf);
   }
} // end of namespace ROOT for class ::HNtuple

namespace {
  void TriggerDictionaryInitialization_MyDict_Impl() {
    static const char* headers[] = {
"src/hntuple.h",
"src/manager.h",
"src/histogram_registry.h",
"src/histogram_factory.h",
"src/histogram_builder.h",
"src/pparticle.h",
"src/boost_frame.h",
"src/ntuple_reader.h",
"src/cut_manager.h",
"src/analysis_config.h",
"src/progressbar.h",
nullptr
    };
    static const char* includePaths[] = {
"/usr/local/root/include/",
"/home/przygoda/HADES/FAT/FAT/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MyDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$src/hntuple.h")))  HNtuple;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MyDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/hntuple.h"
#include "src/manager.h"
#include "src/histogram_registry.h"
#include "src/histogram_factory.h"
#include "src/histogram_builder.h"
#include "src/pparticle.h"
#include "src/boost_frame.h"
#include "src/ntuple_reader.h"
#include "src/cut_manager.h"
#include "src/analysis_config.h"
#include "src/progressbar.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"HNtuple", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MyDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MyDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MyDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MyDict() {
  TriggerDictionaryInitialization_MyDict_Impl();
}
