// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME HitRemoval_NuTrackRemovalCint

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
#include "PandoraLinearRemoval.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLPandoraLinearRemoval_Dictionary();
   static void larlitecLcLPandoraLinearRemoval_TClassManip(TClass*);
   static void *new_larlitecLcLPandoraLinearRemoval(void *p = 0);
   static void *newArray_larlitecLcLPandoraLinearRemoval(Long_t size, void *p);
   static void delete_larlitecLcLPandoraLinearRemoval(void *p);
   static void deleteArray_larlitecLcLPandoraLinearRemoval(void *p);
   static void destruct_larlitecLcLPandoraLinearRemoval(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::PandoraLinearRemoval*)
   {
      ::larlite::PandoraLinearRemoval *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::PandoraLinearRemoval));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::PandoraLinearRemoval", "PandoraLinearRemoval.h", 30,
                  typeid(::larlite::PandoraLinearRemoval), DefineBehavior(ptr, ptr),
                  &larlitecLcLPandoraLinearRemoval_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::PandoraLinearRemoval) );
      instance.SetNew(&new_larlitecLcLPandoraLinearRemoval);
      instance.SetNewArray(&newArray_larlitecLcLPandoraLinearRemoval);
      instance.SetDelete(&delete_larlitecLcLPandoraLinearRemoval);
      instance.SetDeleteArray(&deleteArray_larlitecLcLPandoraLinearRemoval);
      instance.SetDestructor(&destruct_larlitecLcLPandoraLinearRemoval);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::PandoraLinearRemoval*)
   {
      return GenerateInitInstanceLocal((::larlite::PandoraLinearRemoval*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::PandoraLinearRemoval*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLPandoraLinearRemoval_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::PandoraLinearRemoval*)0x0)->GetClass();
      larlitecLcLPandoraLinearRemoval_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLPandoraLinearRemoval_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLPandoraLinearRemoval(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::PandoraLinearRemoval : new ::larlite::PandoraLinearRemoval;
   }
   static void *newArray_larlitecLcLPandoraLinearRemoval(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::PandoraLinearRemoval[nElements] : new ::larlite::PandoraLinearRemoval[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLPandoraLinearRemoval(void *p) {
      delete ((::larlite::PandoraLinearRemoval*)p);
   }
   static void deleteArray_larlitecLcLPandoraLinearRemoval(void *p) {
      delete [] ((::larlite::PandoraLinearRemoval*)p);
   }
   static void destruct_larlitecLcLPandoraLinearRemoval(void *p) {
      typedef ::larlite::PandoraLinearRemoval current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::PandoraLinearRemoval

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 214,
                  typeid(vector<int>), DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 214,
                  typeid(vector<double>), DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_libHitRemoval_NuTrackRemoval_Impl() {
    static const char* headers[] = {
"PandoraLinearRemoval.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/uboone/larlite/UserDev/BasicTool",
"/home/david/uboone/larlite/UserDev",
"/home/david/SOFTWARE/ROOT_v60410/include",
"/home/david/uboone/larlite/UserDev/HitRemoval/NuTrackRemoval/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$PandoraLinearRemoval.h")))  PandoraLinearRemoval;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "PandoraLinearRemoval.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::PandoraLinearRemoval", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libHitRemoval_NuTrackRemoval",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libHitRemoval_NuTrackRemoval_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libHitRemoval_NuTrackRemoval_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libHitRemoval_NuTrackRemoval() {
  TriggerDictionaryInitialization_libHitRemoval_NuTrackRemoval_Impl();
}
