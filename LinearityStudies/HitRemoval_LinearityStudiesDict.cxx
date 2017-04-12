// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HitRemoval_LinearityStudiesDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "ClusterLinearityStudy.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLClusterLinearityStudy_Dictionary();
   static void larlitecLcLClusterLinearityStudy_TClassManip(TClass*);
   static void *new_larlitecLcLClusterLinearityStudy(void *p = 0);
   static void *newArray_larlitecLcLClusterLinearityStudy(Long_t size, void *p);
   static void delete_larlitecLcLClusterLinearityStudy(void *p);
   static void deleteArray_larlitecLcLClusterLinearityStudy(void *p);
   static void destruct_larlitecLcLClusterLinearityStudy(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::ClusterLinearityStudy*)
   {
      ::larlite::ClusterLinearityStudy *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::ClusterLinearityStudy));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::ClusterLinearityStudy", "ClusterLinearityStudy.h", 27,
                  typeid(::larlite::ClusterLinearityStudy), DefineBehavior(ptr, ptr),
                  &larlitecLcLClusterLinearityStudy_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::ClusterLinearityStudy) );
      instance.SetNew(&new_larlitecLcLClusterLinearityStudy);
      instance.SetNewArray(&newArray_larlitecLcLClusterLinearityStudy);
      instance.SetDelete(&delete_larlitecLcLClusterLinearityStudy);
      instance.SetDeleteArray(&deleteArray_larlitecLcLClusterLinearityStudy);
      instance.SetDestructor(&destruct_larlitecLcLClusterLinearityStudy);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::ClusterLinearityStudy*)
   {
      return GenerateInitInstanceLocal((::larlite::ClusterLinearityStudy*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::ClusterLinearityStudy*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLClusterLinearityStudy_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::ClusterLinearityStudy*)0x0)->GetClass();
      larlitecLcLClusterLinearityStudy_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLClusterLinearityStudy_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLClusterLinearityStudy(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::ClusterLinearityStudy : new ::larlite::ClusterLinearityStudy;
   }
   static void *newArray_larlitecLcLClusterLinearityStudy(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::ClusterLinearityStudy[nElements] : new ::larlite::ClusterLinearityStudy[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLClusterLinearityStudy(void *p) {
      delete ((::larlite::ClusterLinearityStudy*)p);
   }
   static void deleteArray_larlitecLcLClusterLinearityStudy(void *p) {
      delete [] ((::larlite::ClusterLinearityStudy*)p);
   }
   static void destruct_larlitecLcLClusterLinearityStudy(void *p) {
      typedef ::larlite::ClusterLinearityStudy current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::ClusterLinearityStudy

namespace {
  void TriggerDictionaryInitialization_libHitRemoval_LinearityStudies_Impl() {
    static const char* headers[] = {
"ClusterLinearityStudy.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/uboone/larlite/UserDev/BasicTool",
"/home/david/SOFTWARE/ROOT_v60410/include",
"/home/david/uboone/larlite/UserDev/HitRemoval/LinearityStudies/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$ClusterLinearityStudy.h")))  ClusterLinearityStudy;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "ClusterLinearityStudy.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::ClusterLinearityStudy", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libHitRemoval_LinearityStudies",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libHitRemoval_LinearityStudies_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libHitRemoval_LinearityStudies_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libHitRemoval_LinearityStudies() {
  TriggerDictionaryInitialization_libHitRemoval_LinearityStudies_Impl();
}
