#include "VistleAnalysisAdaptor.h"

#include "MeshMetadata.h"
#include "VTKUtils.h"
#include "STLUtils.h"
#include "MPIUtils.h"
#include "Profiler.h"
#include "Error.h"
#include "BinaryStream.h"
#include "DataAdaptor.h"

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
#include <memory>

#include <mpi.h>

#include <vistle/insitu/sensei/sensei.h>
#include <vistle/vtk/vtktovistle.h>

#define CERR std::cerr
using std::endl;
using vtkDataObjectPtr = vtkSmartPointer<vtkDataObject>;
// --------------------------------------------------------------------------
using namespace vistle::insitu::sensei;
namespace sensei
{

    class VistleAnalysisAdaptor::PrivateData
    {
    public:
        bool Initialize(DataAdaptor *data);
        bool Execute(DataAdaptor *data);
        int Finalize();

        // Set some Libsim startup options.
        void SetTraceFile(const std::string &traceFile);
        void SetOptions(const std::string &options);
        void SetMode(const std::string &mode);
        void SetFrequency(int f);

    private:
        std::string m_mode;
        std::string m_options;
        std::string m_tracefile;
        int m_frequency = 1;
        std::unique_ptr<SenseiAdapter> m_vistleAdaptor = nullptr;
        MPI_Comm m_comm = MPI_COMM_WORLD;
        DataAdaptor *m_senseiAdaptor = nullptr;
        bool m_initialized = false;
        std::vector<MeshMetadataPtr> m_mehsMetaData;
        MeshMetadataPtr findMeshMeta(const std::string &name);
        ObjectRetriever::PortAssignedObjectList getData(const MetaData &meta);
        struct VtkAndVistleMesh
        {
            vtkCompositeDataSetPtr vtkMesh = nullptr;
            std::vector<vistle::Object::ptr> vistleMeshes;
            operator bool() const;
        };
        std::map<std::string, VtkAndVistleMesh> m_cachedMeshes;
        std::map<std::string, MeshMetadataPtr> m_senseiMetaData;

        bool chacheMetaData(DataAdaptor *data);
        MetaData generateVistleMetaData();
        bool getMesh(const std::string &meshName, const MeshMetadataPtr meta, VtkAndVistleMesh &mesh);
        bool getVariable(ObjectRetriever::PortAssignedObjectList &output, const std::string &varName, const std::string &meshNAme, const VtkAndVistleMesh &mesh, const MeshMetadataPtr meshMeta);
        size_t getBlockIndex(vtkCompositeDataIterator *iter, vtkCompositeDataSet *mesh);
        VtkAndVistleMesh getMeshFromSim(const std::string &name, const MeshMetadataPtr &meshMeta);
        vtkDataObjectPtr getMeshFromSim(const std::string &name);

        bool addGhostCells(const std::string &meshName, const MeshMetadataPtr meshMeta, vtkDataObjectPtr vtkObject);
        VtkAndVistleMesh createVistleMesh(const vtkDataObjectPtr vtkObject);
    };

    VistleAnalysisAdaptor::PrivateData::VtkAndVistleMesh::operator bool() const
    {
        return vtkMesh && vistleMeshes.size() > 0;
    }

    bool
    sensei::VistleAnalysisAdaptor::PrivateData::Initialize(DataAdaptor *data)
    {
        if (m_initialized)
            return 0;

        TimeEvent<128> mark("vistle::initialize");
#ifdef VISTLE_DEBUG_LOG
        CERR << "SENSEI: VistleAnalysisAdaptor::Initialize" << endl;
#endif

        if (chacheMetaData(data))
        {
            MetaData vistleMeta = generateVistleMetaData();

            auto getDataFunc = std::bind(&VistleAnalysisAdaptor::PrivateData::getData, this, std::placeholders::_1);
            bool pause = m_mode == "interactive,paused";
            m_vistleAdaptor = std::unique_ptr<SenseiAdapter>(new SenseiAdapter(pause, m_comm, std::move(vistleMeta), ObjectRetriever{getDataFunc}));

            m_initialized = true;
            return 0;
        }

        return 1;
    }

    bool sensei::VistleAnalysisAdaptor::PrivateData::Execute(sensei::DataAdaptor *DataAdaptor)
    {
        if (!m_initialized)
        {
            Initialize(DataAdaptor);
        }
        // Keep a pointer to the data adaptor so the callbacks can access it.
        m_senseiAdaptor = DataAdaptor;
        return m_vistleAdaptor->Execute(DataAdaptor->GetDataTimeStep());
    }

    int sensei::VistleAnalysisAdaptor::PrivateData::Finalize()
    {
        if (m_vistleAdaptor)
        {
            return !m_vistleAdaptor->Finalize();
        }
        return false;
    }

    void sensei::VistleAnalysisAdaptor::PrivateData::SetTraceFile(const std::string &traceFile)
    {
        m_tracefile = traceFile;
    }

    void sensei::VistleAnalysisAdaptor::PrivateData::SetOptions(const std::string &options)
    {
        m_options = options;
    }

    void sensei::VistleAnalysisAdaptor::PrivateData::SetMode(const std::string &mode)
    {
        m_mode = mode;
    }

    void sensei::VistleAnalysisAdaptor::PrivateData::SetFrequency(int f)
    {
        m_frequency = f;
    }
    // private -------------------------------------------------------------------------

    bool sensei::VistleAnalysisAdaptor::PrivateData::chacheMetaData(DataAdaptor *data)
    {
        unsigned int numMeshes = 0;
        if (data->GetNumberOfMeshes(numMeshes))
        {
            SENSEI_ERROR("Failed to get number of meshes")
            return false;
        }
        for (unsigned int i = 0; i < numMeshes; i++)
        {
            MeshMetadataPtr meshMeta = MeshMetadata::New();
            if (data->GetMeshMetadata(i, meshMeta))
            {
                SENSEI_ERROR("Failed to get mesh meta data " + std::to_string(i))
                return false;
            }
            m_senseiMetaData[meshMeta->MeshName] = meshMeta;
        }
        return true;
    }

    MetaData sensei::VistleAnalysisAdaptor::PrivateData::generateVistleMetaData()
    {
        MetaData metaData;
        for (const auto meshMeta : m_senseiMetaData)
        {
            MetaMesh vistleMeshMeta(meshMeta.first);
            for (auto var : meshMeta.second->ArrayName)
            {
                vistleMeshMeta.addVar(var);
            }
            metaData.addMesh(vistleMeshMeta);
        }
        return metaData;
    }

    sensei::MeshMetadataPtr sensei::VistleAnalysisAdaptor::PrivateData::findMeshMeta(const std::string &name)
    {
        auto ptr = std::find_if(m_mehsMetaData.begin(), m_mehsMetaData.end(), [name](const sensei::MeshMetadataPtr ptr) {
            return ptr->MeshName == name;
        });
        if (ptr == m_mehsMetaData.end())
            return nullptr;
        return *ptr;
    }

    ObjectRetriever::PortAssignedObjectList sensei::VistleAnalysisAdaptor::PrivateData::getData(const MetaData &meta)
    {
        std::vector<ObjectRetriever::PortAssignedObject> outputData;
        for (const auto &meshIter : meta)
        {
            const std::string &meshName = meshIter.name();

            // get the metadata, it should already be available
            auto mdit = this->m_senseiMetaData.find(meshName);
            if (mdit == this->m_senseiMetaData.end())
            {
                SENSEI_ERROR("No metadata for mesh \"" << meshName << "\"")
                continue;
            }
            MeshMetadataPtr meshMeta = mdit->second;
            VtkAndVistleMesh mesh;
            if (!getMesh(meshName, meshMeta, mesh))
                continue;

            for (const auto &subMesh : mesh.vistleMeshes)
            {
                outputData.push_back(ObjectRetriever::PortAssignedObject{meshName, subMesh});
            }

            for (const auto &varName : meshIter)
            {
                getVariable(outputData, varName, meshName, mesh, meshMeta);
            }
        }
        return outputData;
    }

    VistleAnalysisAdaptor::PrivateData::VtkAndVistleMesh VistleAnalysisAdaptor::PrivateData::getMeshFromSim(const std::string &name, const MeshMetadataPtr &meshMeta)
    {
        auto vtkObject = getMeshFromSim(name);
        if (!vtkObject || addGhostCells(name, meshMeta, vtkObject))
            return VtkAndVistleMesh();

        auto vtkAndVistleMesh = createVistleMesh(vtkObject);
        if (meshMeta->StaticMesh)
        {
            m_cachedMeshes[name] = vtkAndVistleMesh;
        }
        return vtkAndVistleMesh;
    }

    vtkDataObjectPtr VistleAnalysisAdaptor::PrivateData::getMeshFromSim(const std::string &name)
    {
        vtkDataObject *dobj = nullptr;
        if (this->m_senseiAdaptor->GetMesh(name, false, dobj))
        {
            SENSEI_ERROR("Failed to get mesh \"" << name << "\"")
            return nullptr;
        }
        vtkDataObjectPtr vtkObject = nullptr;
        vtkObject.TakeReference(dobj);
        return vtkObject;
    }

    bool VistleAnalysisAdaptor::PrivateData::addGhostCells(const std::string &meshName, const MeshMetadataPtr meshMeta, vtkDataObjectPtr vtkObject)
    {
        if (meshMeta->NumGhostCells)
        {
            if (m_senseiAdaptor->AddGhostCellsArray(vtkObject, meshName))
            {
                SENSEI_ERROR("Failed to add ghost cells to mesh \"" << meshName << "\"")
                return -1;
            }
            if (m_senseiAdaptor->AddGhostNodesArray(vtkObject, meshName))
            {
                SENSEI_ERROR("Failed to add ghost nodes to mesh \"" << meshName << "\"")
                return -1;
            }
            return 0;
        }
    }

    VistleAnalysisAdaptor::PrivateData::VtkAndVistleMesh VistleAnalysisAdaptor::PrivateData::createVistleMesh(const vtkDataObjectPtr vtkObject)
    {
        VtkAndVistleMesh vtkAndVistleMesh;
        vtkAndVistleMesh.vtkMesh = VTKUtils::AsCompositeData(this->m_comm, vtkObject.GetPointer(), false);
        vtkCompositeDataIterator *dataSetIter = vtkAndVistleMesh.vtkMesh->NewIterator();

        for (dataSetIter->InitTraversal(); !dataSetIter->IsDoneWithTraversal(); dataSetIter->GoToNextItem())
        {
            long blockId = getBlockIndex(dataSetIter, vtkAndVistleMesh.vtkMesh.Get());
            vtkDataObject *vtkMesh = dataSetIter->GetCurrentDataObject();
            auto vistleMesh = vistle::vtk::toGrid(*m_vistleAdaptor, vtkMesh);
            vistleMesh->setBlock(blockId);

            // VTK's iterators for AMR datasets behave differently than for multiblock datasets.
            vtkUniformGridAMRDataIterator *amrIt = dynamic_cast<vtkUniformGridAMRDataIterator *>(dataSetIter);
            if (amrIt)
            {
                vtkOverlappingAMR *amrMesh = dynamic_cast<vtkOverlappingAMR *>(vtkAndVistleMesh.vtkMesh.Get());
                //As long as there is no Vistle amr structure: set more refinded partitions in front
                vistleMesh->addAttribute("POLYGON_OFFSET", std::to_string((amrMesh->GetNumberOfLevels() - amrIt->GetCurrentLevel())));
            }

            vtkAndVistleMesh.vistleMeshes.push_back(vistleMesh);
        }
        dataSetIter->Delete();
        return vtkAndVistleMesh;
    }

    bool VistleAnalysisAdaptor::PrivateData::getMesh(const std::string &meshName, const MeshMetadataPtr meta, VtkAndVistleMesh &mesh)
    {
        auto cachedMeshPair = m_cachedMeshes.find(meshName);
        if (cachedMeshPair != m_cachedMeshes.end())
        {
            mesh = cachedMeshPair->second;
        }
        else
        {
            mesh = getMeshFromSim(meshName, meta);
        }
        if (!mesh)
        {
            SENSEI_ERROR("Failed to get mesh \"" << meshName << "\"from sim")
            return false;
        }
        return true;
    }

    bool VistleAnalysisAdaptor::PrivateData::getVariable(std::vector<ObjectRetriever::PortAssignedObject> &output, const std::string &varName, const std::string &meshName, const VtkAndVistleMesh &mesh, const MeshMetadataPtr meshMeta)
    {
        auto centeringPosIt = std::find(meshMeta->ArrayName.begin(), meshMeta->ArrayName.end(), varName);
        if (centeringPosIt == meshMeta->ArrayName.end())
        {
            SENSEI_ERROR("Failed to get metadata for array \"" << varName << "\"")
            return false;
        }
        size_t centeringPos = centeringPosIt - meshMeta->ArrayName.begin();
        auto centering = meshMeta->ArrayCentering[centeringPos];

        vtkCompositeDataIterator *dataSetIter = mesh.vtkMesh->NewIterator();

        // this rank has no local data
        if (dataSetIter->IsDoneWithTraversal())
        {
            dataSetIter->Delete();
            return false;
        }

        // read the array if we have not yet
        if (!dataSetIter->GetCurrentDataObject()->GetAttributes(centering)->GetArray(varName.c_str()))
        {
            if (this->m_senseiAdaptor->AddArray(mesh.vtkMesh.GetPointer(), meshName, centering, varName))
            {
                SENSEI_ERROR("Failed to add " << VTKUtils::GetAttributesName(centering)
                                              << " data array \"" << varName << "\"")
                dataSetIter->Delete();
                return false;
            }
        }

        // VTK's iterators for AMR datasets behave differently than for multiblock
        // datasets.  we are going to have to handle AMR data as a special case for
        // now.

        size_t index = 0;
        for (dataSetIter->InitTraversal(); !dataSetIter->IsDoneWithTraversal(); dataSetIter->GoToNextItem())
        {
            long blockId = getBlockIndex(dataSetIter, mesh.vtkMesh.Get());
            auto currentObj = dataSetIter->GetCurrentDataObject();
            auto currentAttributes = currentObj->GetAttributes(centering);
            vtkDataArray *array = currentAttributes->GetArray(varName.c_str());
            if (!array)
            {
                SENSEI_ERROR("Failed to get array \"" << varName << "\"");
                ++index;
                continue;
            }
            auto vistleArray = vistle::vtk::vtkData2Vistle(*m_vistleAdaptor, array, mesh.vistleMeshes[index]);
            if (!vistleArray)
            {
                CERR << "Failed to convert vtk data array of type " << array->GetDataType() << " to vistle" << endl;
                continue;
            }

            vistleArray->setBlock(blockId);
            vistleArray->addAttribute("_species", varName);
            output.push_back(ObjectRetriever::PortAssignedObject{meshName, varName, vistleArray});
            ++index;
        }

        dataSetIter->Delete();
        return true;
    }

    size_t VistleAnalysisAdaptor::PrivateData::getBlockIndex(vtkCompositeDataIterator *iter, vtkCompositeDataSet *mesh)
    {

        vtkUniformGridAMRDataIterator *amrIt = dynamic_cast<vtkUniformGridAMRDataIterator *>(iter);
        if (amrIt)
        {
            //return 0; //GetAMRBlockSourceIndex causes segmentation fault if the sim did not call SetAMRBlockSourceIndex
            // special case for AMR
            int level = amrIt->GetCurrentLevel();
            int index = amrIt->GetCurrentIndex();
            vtkOverlappingAMR *amrMesh = dynamic_cast<vtkOverlappingAMR *>(mesh);

            return amrMesh->GetAMRBlockSourceIndex(level, index);
        }
        else
        {
            // other composite data
            return iter->GetCurrentFlatIndex() - 1;
        }
    }

    //-----------------------------------------------------------------------------
    // LibsimAnalysisAdaptor PUBLIC INTERFACE
    //-----------------------------------------------------------------------------
    senseiNewMacro(VistleAnalysisAdaptor);

    sensei::VistleAnalysisAdaptor::VistleAnalysisAdaptor()
    {
        m_internals = new VistleAnalysisAdaptor::PrivateData{};
    }

    sensei::VistleAnalysisAdaptor::~VistleAnalysisAdaptor()
    {
        delete m_internals;
    }

    bool
    sensei::VistleAnalysisAdaptor::Execute(sensei::DataAdaptor *DataAdaptor)
    {
        return m_internals->Execute(DataAdaptor);
    }

    int sensei::VistleAnalysisAdaptor::Finalize()
    {
        return m_internals->Finalize();
    }

    void sensei::VistleAnalysisAdaptor::SetTraceFile(const std::string &traceFile)
    {
        m_internals->SetTraceFile(traceFile);
    }

    void sensei::VistleAnalysisAdaptor::SetOptions(const std::string &options)
    {
        m_internals->SetOptions(options);
    }

    void sensei::VistleAnalysisAdaptor::SetMode(const std::string &mode)
    {
        m_internals->SetMode(mode);
    }

    void sensei::VistleAnalysisAdaptor::SetFrequency(int f)
    {
        m_internals->SetFrequency(f);
    }

    // --------------------------------------------------------------------------

} // namespace sensei