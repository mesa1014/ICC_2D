
#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
// #include "../src/imtiaz_2002d_noTstart_COR.hpp"
#include "../src/CellICCBioPhy.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "../src/DummyDerivedCa.hpp"
#include "Debug.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "BidomainProblem.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractElement.hpp"

// for mechanics
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "QuadraticMesh.hpp"
#include "NonlinearElasticityTools.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"

// cell factories
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"


using namespace std;

class ICCCellFactory : public AbstractCardiacCellFactory<2>
{

public:
    ICCCellFactory():AbstractCardiacCellFactory<2>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        // Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(mpSolver, mpZeroStimulus);
        // cell->SetParameter("eta", 0.04);
        CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
        cell->SetParameter("V_excitation", -55); // excitation value (threshold) default: -55mV
        cell->SetParameter("live_time", 10000); // time of resting potential
        cell->SetParameter("ode_time_step", 0.1); // Set the same as defined in HeartHandler
        cell->SetParameter("IP3Par", 0.00069); //
        cell->SetParameter("t_start", 600000); // Set larger than total simulation time
        // cell->SetParameter("ICC_Membrane__Ca_i", 0.0002); // Set larger than total simulation time

        //if(ellipseRegion.DoesContain(myPoint))
        //    cell->SetParameter("V_excitation", -68);

        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        ChastePoint<2> centre(0.06,0.3);
        ChastePoint<2> radii (0.01,0.01);
        ChasteEllipsoid<2> ellipseRegion(centre, radii);
        ChastePoint<2> myPoint(x, y);

        if(ellipseRegion.DoesContain(myPoint))
        {
          cell->SetParameter("t_start", 0); // Set larger than total simulation time
        }
        return cell;
    }
};

class Test2DMonodomain : public CxxTest::TestSuite
{
public:
    void TestSimulation() //throw(Exception)
    {
        ///// Read electric mesh file
        TetrahedralMesh<2,2> mesh;
        std::string myFile = "MeshNetwork-2D.1";
        std::string meshFile = "/home/chaste/projects/mesh/" + myFile;
        TrianglesMeshReader<2,2> mesh_reader(meshFile.c_str());
        mesh.ConstructFromMeshReader(mesh_reader);

        ///// Read mechanics mesh file
        QuadraticMesh<2> mechanics_mesh;
        std::string myFile_m = "MeshNetwork-2D.2";
        std::string meshFile_m = "/home/chaste/projects/mesh/" + myFile_m;
        TrianglesMeshReader<2,2> mesh_reader_m(meshFile_m.c_str(),2);
        mechanics_mesh.ConstructFromMeshReader(mesh_reader_m);

        ///// fixed nodes are all x = 0.0 nodes
        std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 0, 0.0);


        ///// Simulation settings
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(2.5);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12));
        HeartConfig::Instance()->SetSimulationDuration(3000);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,100);

        ///// Output file/folder
        std::string mod = myFile + "-Imtiaz";
        HeartConfig::Instance()->SetOutputDirectory(mod.c_str());
	      HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        ///// Output visualization options
        HeartConfig::Instance()->SetVisualizeWithCmgui(true);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
        HeartConfig::Instance()->SetVisualizeWithVtk(false);

        ///// Cell factory
        ICCCellFactory cell_factory;
        // PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 2> cell_factory(-5000*1000);

        // ///// MonodomainProblem
        // MonodomainProblem<2> monodomain_problem(&cell_factory);
        // monodomain_problem.SetMesh(&mesh);
        // monodomain_problem.Initialise();
        //
        // ///// Solve
        // monodomain_problem.Solve();


        ///// Cardiac ElectroMechanics Problem definition and solve
        ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
        problem_defn.SetContractionModel(NASH2004,0.01/*contraction model ODE timestep*/);
        problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
        problem_defn.SetZeroDisplacementNodes(fixed_nodes);
        problem_defn.SetMechanicsSolveTimestep(10.0);

        CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
                                                    MONODOMAIN,
                                                    &mesh,
                                                    &mechanics_mesh,
                                                    &cell_factory,
                                                    &problem_defn,
                                                    "icc2d");

        problem.Solve();


    }
};
