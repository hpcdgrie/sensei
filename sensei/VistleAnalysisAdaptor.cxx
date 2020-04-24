#include "VistleAnalysisAdaptor.h"

#include "MeshMetadata.h"
#include "VTKUtils.h"
#include "STLUtils.h"
#include "MPIUtils.h"
#include "Profiler.h"
#include "Error.h"
#include "BinaryStream.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkAMRBox.h>
#include <vtkAMRInformation.h>
#include <vtkOverlappingAMR.h>
#include <vtkUniformGrid.h>
#include <vtkDataObject.h>
#include <vtkXMLUniformGridAMRWriter.h>
#include <vtkUniformGridAMRDataIterator.h>

#include <sstream>
#include <algorithm>
#include <map>

#include <mpi.h>

#include <sensei.h>
#include <dataAdaptor.h>
using namespace sensei;

/ --------------------------------------------------------------------------
bool
VistleAnalysisAdaptor::Initialize()
{
    if (initialized)
        return true;

     TimeEvent<128> mark("vistle::initialize");
#ifdef VISTLE_DEBUG_LOG
    VisItDebug5("SENSEI: VistleAnalysisAdaptor::Initialize\n");
#endif

     int rank = 0, size = 1;
     MPI_Comm_rank(Comm, &rank);
     MPI_Comm_size(Comm, &size);


     // Install Comm into VisIt.
     VisItSetMPICommunicator((void *)&Comm);

     // Set up the environment.
     char *env = nullptr;
     if(rank == 0)
         env = VisItGetEnvironment();
     VisItSetupEnvironment2(env);
     if(env != nullptr)
         free(env);

     bool i0 = mode == "interactive";
     bool i1 = mode == "interactive,paused";
     if(i0 || i1)
     {
         // We can start paused if desired.
         this->paused = i1;

         // Write out .sim file that VisIt uses to connect.
         if (rank == 0)
         {
             VisItInitializeSocketAndDumpSimFile(
                 "sensei",
                 "Connected via SENSEI",
                 "/path/to/where/sim/was/started",
                 NULL, NULL, "sensei.sim2");
         }

         initialized = true;
     }
     else
     {
        // Try and initialize the runtime.
        if(VisItInitializeRuntime() == VISIT_ERROR)
        {
            SENSEI_ERROR("Could not initialize the VisIt runtime library.")
            return false;
        }
        else
        {
            // Register Libsim callbacks.
            VisItSetSlaveProcessCallback2(SlaveProcessCallback, (void*)this); // needed in batch?
            VisItSetGetMetaData(GetMetaData, (void*)this);
            VisItSetGetMesh(GetMesh, (void*)this);
            VisItSetGetVariable(GetVariable, (void*)this);
            VisItSetGetDomainList(GetDomainList, (void*)this);
            if (this->ComputeNesting)
                VisItSetGetDomainNesting(GetDomainNesting, (void*)this);

            initialized = true;
        }
    }

    return initialized;
}




// --------------------------------------------------------------------------
bool
VistleAnalysisAdaptor::Execute(sensei::DataAdaptor *DataAdaptor)
{
#ifdef VISIT_DEBUG_LOG
    VisItDebug5("SENSEI: LibsimAnalysisAdaptor::PrivateData::Execute\n");
#endif

    if (!initialized)
    {
        SENSEI_ERROR("Execute before initialization")
        return false;
    }

    // Keep a pointer to the data adaptor so the callbacks can access it.
    Adaptor = DataAdaptor;

    int rank = 0;
    MPI_Comm_rank(Comm, &rank);

    // Let visit get new metadata.
    VisItTimeStepChanged();

    // Execute
    bool retval = true;
    if (mode.substr(0, 11) == "interactive")
        retval = Execute_Interactive(rank);
    else
        retval = Execute_Batch(rank);

    // during execution data and metadata are cached due to
    // the differnece between how sensei presents data and
    // visit consumes it. you must clear the cache after each
    // execute.
    ClearCache();

    return retval;
}

// --------------------------------------------------------------------------
