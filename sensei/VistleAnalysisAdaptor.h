#ifndef sensei_VistleAnalysisAdaptor_h
#define sensei_VistleAnalysisAdaptor_h
#include "AnalysisAdaptor.h"
#include <string>
#include <mpi.h>

namespace sensei
{
class VistleAnalysisAdaptor : public AnalysisAdaptor
{
   public:
    static VistleAnalysisAdaptor* New();
    senseiTypeMacro(VistleAnalysisAdaptor, AnalysisAdaptor);


    // Let the caller explicitly initialize.
    void Initialize();

    bool Execute(DataAdaptor* data) override;

    int Finalize() override;

  protected:
    VistleAnalysisAdaptor();
    ~VistleAnalysisAdaptor();

  private:
    VistleAnalysisAdaptor(const LibsimAnalysisAdaptor&); // Not implemented.
    void operator=(const VistleAnalysisAdaptor&); // Not implemented.

    class PrivateData;
    PrivateData *internals;
    bool initialized = false;
};


}




#endif
