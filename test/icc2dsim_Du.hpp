
#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "../src/imtiaz_2002d_noTstart_COR.hpp"
#include "../src/CellICCBioPhy.hpp"
#include "../src/Du2013_neural_sens.hpp"
#include "../src/DummyDerivedCa.hpp"
#include "../src/ICCFactory.hpp"

#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include "Debug.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "BidomainProblem.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractElement.hpp"

#include "FileComparison.hpp"
#include "FileFinder.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"


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
#include "FungMaterialLaw.hpp"
#include "SchmidCostaExponentialLaw2d.hpp"

// cell factories
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"

// #include <GaussianQuadratureRule.hpp>
// #include <QuadraturePointsGroup.hpp>



using namespace std;
class Test2DMonodomainV2 : public CxxTest::TestSuite
{
public:
  void TestSimulationV2() //throw(Exception)
  {
    ///// Read electric mesh file
    TetrahedralMesh<2,2> mesh;
    std::string myFile = "mesh2by3_s.1";
    std::string meshFile = "../projects/mesh/2by3/structured/" + myFile;
    TrianglesMeshReader<2,2> mesh_reader(meshFile.c_str());
    mesh.ConstructFromMeshReader(mesh_reader);

    ///// Read mechanics mesh f1ile
    QuadraticMesh<2> mechanics_mesh;
    std::string myFile_m = "mesh2by3_s.2";
    std::string meshFile_m = "../projects/mesh/2by3/structured/" + myFile_m;
    TrianglesMeshReader<2,2> mesh_reader_m(meshFile_m.c_str(),2);
    mechanics_mesh.ConstructFromMeshReader(mesh_reader_m);

    // TetrahedralMesh<2,2> mesh;
    // mesh.ConstructRegularSlabMesh(0.1/*stepsize*/, 2.0/*length*/, 3.0/*width*/, 0.2/*depth*/);
    //
    // QuadraticMesh<2> mechanics_mesh;Items in trash are deleted forever after 30 days
    // mechanics_mesh.ConstructRegularSlabMesh(0.2, 2.0, 3.0, 0.2 /*as above with a different stepsize*/);


    ///// fixed nodes are all nodes on top and bottom
    ///// Don't forget to change top fixed nodes, if you changed the mesh size
    std::vector<unsigned> fixed_nodes;
    std::vector<unsigned> fixed_nodes_bottom  = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 1, 0.0);
    std::vector<unsigned> fixed_nodes_top     = NonlinearElasticityTools<2>::GetNodesByComponentValue(mechanics_mesh, 1, 3.0);

    for(unsigned int j = 0; j < fixed_nodes_top.size(); j++)
    {
      fixed_nodes.push_back(fixed_nodes_top[j]);
      fixed_nodes.push_back(fixed_nodes_bottom[j]);
    }

    // dev test
    std::cout << "fixed node size: " << fixed_nodes.size() << std::endl;
    for(unsigned int k = 0; k < fixed_nodes.size(); k++)
    std::cout << "fixed_nodes:" << fixed_nodes[k] << ' ' << std::endl;

    // dev test
    // for (unsigned i=1; i<mechanics_mesh.GetNumNodes(); i++)
    // {
    //   double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
    //   double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
    //   std::cout <<"node id: " << i << ", X: " << X << " , Y: " << Y << std::endl;
    //
    // }

    // Simulation settings
    HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
    HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
    HeartConfig::Instance()->SetCapacitance(2.5);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12));
    HeartConfig::Instance()->SetSimulationDuration(60000);  //ms.
    // HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,100); // doesn't work here!
    HeartConfig::Instance()->SetPrintingTimeStep(100.0);

    // unsigned a = (unsigned) HeartConfig::Instance()->GetPrintingTimeStep();  double y = pNode->GetPoint()[1];

    // cout <<"Unsigned print timestep: " << a << endl;

    cout << "Print time step: " << HeartConfig::Instance()->GetPrintingTimeStep() << endl;
    HeartConfig::Instance()->SetPdeTimeStep(0.1);
    cout << "pde time step: " << HeartConfig::Instance()->GetPdeTimeStep() << endl;

    // Output visualization options
    HeartConfig::Instance()->SetVisualizeWithCmgui(true);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
    HeartConfig::Instance()->SetVisualizeWithVtk(false);
    // HeartConfig::Instance()->SetOutputVariables(rY[1]);

    // Cell factory
    std::set<unsigned> iccNodes;
    for (unsigned i=0; i < mesh.GetNumAllNodes() ; ++i) iccNodes.insert(i);
    ICCFactory<2> cell_factory(iccNodes);

    // Material law
    MooneyRivlinMaterialLaw<2> law(3.0); // 0.1 or 1.0
    // ExponentialMaterialLaw<2> law2(1000.0, 2.0); // First parameter is 'a', second 'b', in W=a*exp(b(I1-3))
    // FungMaterialLaw<2> law2(1400.0, 39.0, 72.0, 0.4);
    SchmidCostaExponentialLaw2d law2;

    // Fibre directions
    OutputFileHandler handler("");
    FileFinder finder = handler.FindFile("2by3_fibre_s.ortho");

    // Cardiac ElectroMechanics problem definition
    ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
    problem_defn.SetContractionModel(NASH, 0.1/*contraction model ODE timestep*/); //KERCHOFFS2003
    problem_defn.SetVariableFibreSheetDirectionsFile(finder, false);
    // problem_defn.SetSolverType(IMPLICIT);
    // problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, false /*deformation affects cell models*/);
    // problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
    problem_defn.SetMaterialLaw(INCOMPRESSIBLE, &law);
    problem_defn.SetZeroDisplacementNodes(fixed_nodes);
    problem_defn.SetMechanicsSolveTimestep(100.0);
    // problem_defn.SetVariableFibreSheetDirectionsFile(finder, true);
    problem_defn.SetSolveUsingSnes();
    problem_defn.SetVerboseDuringSolve(true);

    CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
      MONODOMAIN,
      &mesh,
      &mechanics_mesh,
      &cell_factory,
      &problem_defn,
      "icc2d_Du");

      c_vector<double,2> node_to_watch;
      node_to_watch(0) = 1.0;
      node_to_watch(1) = 1.0;
      problem.SetWatchedPosition(node_to_watch);
      problem.SetOutputDeformationGradientsAndStress(100.0);

      // solve
      // problem.SetNoElectricsOutput();
      problem.Initialise();
      problem.Solve();

      // FileFinder test_output_folder("icc2d/electrics", RelativeTo::ChasteTestOutput);
      // Hdf5ToMeshalyzerConverter<2,2> converter(test_output_folder, "voltage",
      //                                               &mesh, false,
      //                                               HeartConfig::Instance()->GetVisualizerOutputPrecision());

    }
  };
