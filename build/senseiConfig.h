#ifndef senseiConfig_h
#define senseiConfig_h

#define SENSEI_VERSION "v0.0.0"
#define SENSEI_VERSION_MAJOR 0
#define SENSEI_VERSION_MINOR 0
#define SENSEI_VERSION_PATCH 0

#define ENABLE_SENSEI
/* #undef ENABLE_PYTHON */
/* #undef ENABLE_CATALYST */
/* #undef ENABLE_CATALYST_PYTHON */
/* #undef ENABLE_LIBSIM */
#define ENABLE_VISTLE
/* #undef ENABLE_ADIOS1 */
/* #undef ENABLE_HDF5 */
/* #undef ENABLE_CONDUIT */
/* #undef ENABLE_ASCENT */
/* #undef ENABLE_VTK_GENERIC_ARRAYS */
/* #undef ENABLE_VTK_MPI */
/* #undef ENABLE_VTK_IO */
#define ENABLE_VTK_RENDERING
/* #undef ENABLE_VTK_ACCELERATORS */
/* #undef ENABLE_VTK_FILTERS */
/* #undef ENABLE_VTKM */
/* #undef ENABLE_PROFILER */

/* #undef SENSEI_PYTHON_VERSION */

#ifdef __cplusplus
// hide some differences betweem VisIt's VTK and more modern versions
#include <vtkSetGet.h>
#if defined(ENABLE_LIBSIM)
#include <vtkVersionMacros.h>
#if VTK_MAJOR_VERSION < 7
#define senseiBaseTypeMacro(a1, a2) vtkTypeMacro(a1, a2)
#else
#define senseiBaseTypeMacro(a1, a2) vtkBaseTypeMacro(a1, a2)
#endif
#define senseiTypeMacro(a1, a2) vtkTypeMacro(a1, a2)
#define senseiNewMacro(a1) \
a1 *a1::New() \
{ \
 vtkObjectBase* ret = vtkObjectFactory::CreateInstance(#a1); \
   if(ret) \
     { \
     return static_cast<a1*>(ret); \
     } \
   return new a1; \
}
#else
#define senseiBaseTypeMacro(a1, a2) vtkBaseTypeMacro(a1, a2)
#define senseiTypeMacro(a1, a2) vtkTypeMacro(a1, a2)
#define senseiNewMacro(a1) vtkStandardNewMacro(a1)
#endif
#endif

#endif
