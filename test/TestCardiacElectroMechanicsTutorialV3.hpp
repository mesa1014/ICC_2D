
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

// #include "FileComparison.hpp"
// #include "FileFinder.hpp"

// for mechanics
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "QuadraticMesh.hpp"
#include "NonlinearElasticityTools.hpp"

#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"

#include "ExponentialMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
// cell factories
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"

// #include <GaussianQuadratureRule.hpp>
// #include <QuadraturePointsGroup.hpp>



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
    double x = pNode->rGetLocation()[0];
    double y = pNode->rGetLocation()[1];
    ChastePoint<2> centre(0.1,0.45);
    ChastePoint<2> radii (0.01,0.01);
    ChasteEllipsoid<2> ellipseRegion(centre, radii);
    ChastePoint<2> myPoint(x, y);
    // double ca = cell->GetIntracellularCalciumConcentration();
    // std::cout << "Ca: " << ca << std::endl;
    if(ellipseRegion.DoesContain(myPoint))
    {
      cell->SetParameter("t_start", 0);
      std::cout<<"I'm inside region!" << std::endl;
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
    // TetrahedralMesh<2,2> mesh;
    // std::string myFile = "MeshNetwork-2D.1";
    // std::string meshFile = "/home/chaste/projects/mesh/" + myFile;
    // TrianglesMeshReader<2,2> mesh_reader(meshFile.c_str());
    // mesh.ConstructFromMeshReader(mesh_reader);
    //
    // ///// Read mechanics mesh f1ile
    // QuadraticMesh<2> mechanics_mesh;
    // std::string myFile_m = "MeshNetwork-2D.2";
    // std::string meshFile_m = "/home/chaste/projects/mesh/" + myFile_m;
    // TrianglesMeshReader<2,2> mesh_reader_m(meshFile_m.c_str(),2);
    // mechanics_mesh.ConstructFromMeshReader(mesh_reader_m);

    TetrahedralMesh<2,2> mesh;
    mesh.ConstructRegularSlabMesh(0.01/*stepsize*/, 0.2/*length*/, 0.5/*width*/, 0.1/*depth*/);

    QuadraticMesh<2> mechanics_mesh;
    mechanics_mesh.ConstructRegularSlabMesh(0.01, 0.2, 0.5, 0.1 /*as above with a different stepsize*/);


    ///// fixed nodes are all nodes on top and bottom
    std::vector<unsigned> fixed_nodes;
    std::vector<unsigned> fixed_nodes_bottom  = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 1, 0.0);
    std::vector<unsigned> fixed_nodes_top     = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 1, 0.5);

    for(unsigned int j = 0; j < fixed_nodes_top.size(); j++)
    {
      fixed_nodes.push_back(fixed_nodes_top[j]);
      fixed_nodes.push_back(fixed_nodes_bottom[j]);
    }

    // dev test
    std::cout << "fixed node size: " << fixed_nodes.size() << std::endl;
    for(unsigned int k = 0; k < fixed_nodes.size(); k++)
    std::cout << "fixed_nodes:" << fixed_nodes[k] << ' ' << std::endl;

    // std::vector<unsigned> fixed_nodes2  = {516, 552, 588, 624, 660, 696, 732, 768, 804, 840, 876, 912, 1815, 1838, 1916, 2150, 2174, 2252, 3138, 3153, 3232, 3457, 3478, 3528, 3548, 0, 37, 74, 111, 148, 185, 222, 259, 296, 333, 370, 407, 444, 916, 935, 997, 1255, 1289, 1338, 2271, 2306, 2354, 2612, 2647, 2695};
    // std::vector<unsigned> fixed_nodes2  = {0, 516, 1338, 1838, 2695, 3548};
    // std::vector<unsigned> fixed_nodes2  = {516, 552, 588, 624, 660, 696, 732, 768, 804, 840, 876, 912,3138, 3153, 3232, 3457, 3478, 3528, 3548}; //148, 185, 222, 259, 296, 333,

    // dev test
    for (unsigned i=1; i<mechanics_mesh.GetNumNodes(); i++)
    {
      double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
      double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
      std::cout <<"node id: " << i << ", X: " << X << " , Y: " << Y << std::endl;

    }

    // Simulation settings
    HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
    HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
    HeartConfig::Instance()->SetCapacitance(2.5);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12));
    HeartConfig::Instance()->SetSimulationDuration(65000);  //ms.
    HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,100); // doesn't work here!

    // Output visualization options
    HeartConfig::Instance()->SetVisualizeWithCmgui(true);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
    HeartConfig::Instance()->SetVisualizeWithVtk(false);
    // HeartConfig::Instance()->SetOutputVariables(rY[1]);
    // Cell factory
    ICCCellFactory cell_factory;

    // Material law
    MooneyRivlinMaterialLaw<2> law(7.0);
    ExponentialMaterialLaw<2> law2(7.0, 1.5); // First parameter is 'a', second 'b', in W=a*exp(b(I1-3))

    // Fibre directions
    // coming soon!

    // Cardiac ElectroMechanics problem definition
    ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
    problem_defn.SetContractionModel(NASH2004,0.1/*contraction model ODE timestep*/);
    // problem_defn.SetSolverType(IMPLICIT);
    // problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, false /*deformation affects cell models*/);
    // problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
    problem_defn.SetMaterialLaw(INCOMPRESSIBLE, &law2);
    problem_defn.SetZeroDisplacementNodes(fixed_nodes);
    problem_defn.SetMechanicsSolveTimestep(100.0);
    // problem_defn.SetVariableFibreSheetDirectionsFile(finder, true);
    problem_defn.SetSolveUsingSnes();

    CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
      MONODOMAIN,
      &mesh,
      &mechanics_mesh,
      &cell_factory,
      &problem_defn,
      "icc2d");

      c_vector<double,2> node_to_watch;
      node_to_watch(0) = 0.1;
      node_to_watch(1) = 0.4;
      problem.SetWatchedPosition(node_to_watch);
      problem.SetOutputDeformationGradientsAndStress(100.0);

      // solve
      problem.Initialise();
      problem.Solve();


    }
  };
