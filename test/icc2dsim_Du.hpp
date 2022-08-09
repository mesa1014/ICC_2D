/*
This code solves electroMechanics for a 2D ICC block:
- electrics solver generates Ca and voltage based on cell model (from CellML)
- contraction model generates active tension based on Ca or voltage
- mechanics solver generates deformation based on the tension and the material law

written by: Mehrdad Sangi (msan223@aucklanduni.ac.nz)
Auckland Bioengineering Institute
The University of Auckland
August 2022
*/
#include <cxxtest/TestSuite.h>
#include <cassert>
#include <set>
#include "Debug.hpp"

//// cell factories
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "../src/imtiaz_2002d_noTstart_COR.hpp"
#include "../src/CellICCBioPhy.hpp"
#include "../src/Du2013_neural_sens.hpp"
#include "../src/DummyDerivedCa.hpp"
#include "../src/ICCFactory.hpp"

//// mesh
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

//// problem & general headers
#include "BidomainProblem.hpp"
#include "MonodomainProblem.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ChasteEllipsoid.hpp"

//// file handler
#include "FileComparison.hpp"
#include "FileFinder.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

//// for mechanics
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "QuadraticMesh.hpp"
#include "NonlinearElasticityTools.hpp"

//// solver
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"

//// material laws
#include "ExponentialMaterialLaw.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "FungMaterialLaw.hpp"
#include "SchmidCostaExponentialLaw2d.hpp"

//// not necessary
// #include <GaussianQuadratureRule.hpp>
// #include <QuadraturePointsGroup.hpp>


class Test2DMonodomainV2 : public CxxTest::TestSuite
{
public:
  void TestSimulationV2() //throw(Exception)
  {
    ///// read electric mesh from file
    //// 2cm x 3cm mesh
    TetrahedralMesh<2,2> mesh;
    std::string myFile = "mesh2by3_s.1";
    std::string meshFile = "../projects/icc2d/mesh/structured/" + myFile;
    TrianglesMeshReader<2,2> mesh_reader(meshFile.c_str());
    mesh.ConstructFromMeshReader(mesh_reader);

    ///// read mechanics mesh f1ile
    QuadraticMesh<2> mechanics_mesh;
    std::string myFile_m = "mesh2by3_s.2";
    std::string meshFile_m = "../projects/icc2d/mesh/structured/" + myFile_m;
    TrianglesMeshReader<2,2> mesh_reader_m(meshFile_m.c_str(),2);
    mechanics_mesh.ConstructFromMeshReader(mesh_reader_m);

    //// create mesh inside Chaste
    // TetrahedralMesh<2,2> mesh;
    // mesh.ConstructRegularSlabMesh(0.1/*stepsize*/, 2.0/*length*/, 3.0/*width*/, 0.2/*depth*/);
    //
    // QuadraticMesh<2> mechanics_mesh;
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

    //// dev test
    std::cout << "fixed node size: " << fixed_nodes.size() << std::endl;
    for(unsigned int k = 0; k < fixed_nodes.size(); k++)
    std::cout << "fixed_nodes:" << fixed_nodes[k] << ' ' << std::endl;

    //// dev test
    // for (unsigned i=1; i<mechanics_mesh.GetNumNodes(); i++)
    // {
    //   double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
    //   double Y = mechanics_mesh.GetNode(i)->rGetLocation()[1];
    //   std::cout <<"node id: " << i << ", X: " << X << " , Y: " << Y << std::endl;
    //
    // }

    //// cell parameters
    HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
    HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
    HeartConfig::Instance()->SetCapacitance(2.5);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12));

    //// simulation settings
    HeartConfig::Instance()->SetSimulationDuration(60000);  //ms
    //// printing time doesn't work in mechanics!
    // HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,100);
    HeartConfig::Instance()->SetPrintingTimeStep(100.0);

    //// dev test
    // cout << "Print time step: " << HeartConfig::Instance()->GetPrintingTimeStep() << endl;
    // HeartConfig::Instance()->SetPdeTimeStep(0.1);
    // cout << "pde time step: " << HeartConfig::Instance()->GetPdeTimeStep() << endl;

    //// Output visualization options (doesn't work in mechanics!)
    HeartConfig::Instance()->SetVisualizeWithCmgui(true);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);
    HeartConfig::Instance()->SetVisualizeWithVtk(false);

    //// Cell factory
    std::set<unsigned> iccNodes;
    for (unsigned i=0; i < mesh.GetNumAllNodes() ; ++i) iccNodes.insert(i);
    ICCFactory<2> cell_factory(iccNodes);

    //// material laws
    MooneyRivlinMaterialLaw<2> law(3.0); // need to check the literature for this constant
    // ExponentialMaterialLaw<2> law2(1000.0, 2.0); // First parameter is 'a', second 'b', in W=a*exp(b(I1-3))
    // FungMaterialLaw<2> law2(1400.0, 39.0, 72.0, 0.4);
    SchmidCostaExponentialLaw2d law2; //// not used here!

    //// fibre directions
    OutputFileHandler handler("");
    FileFinder finder = handler.FindFile("2by3_fibre_s.ortho");

    //// cardiac ElectroMechanics problem definition
    ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
    problem_defn.SetContractionModel(NASH, 0.1/*contraction model ODE timestep*/); //KERCHOFFS2003
    // problem_defn.SetSolverType(IMPLICIT);
    // problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, false /*deformation affects cell models*/);
    // problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
    problem_defn.SetMaterialLaw(INCOMPRESSIBLE, &law);
    problem_defn.SetZeroDisplacementNodes(fixed_nodes);
    problem_defn.SetMechanicsSolveTimestep(100.0); //// this is the same as printing time!
    problem_defn.SetVariableFibreSheetDirectionsFile(finder, false); //bool definedPerQuadPoint -> should be false
    problem_defn.SetSolveUsingSnes(); //// a more stable solver according to Chaste docs
    problem_defn.SetVerboseDuringSolve(true); //// to see more info in solver steps

    CardiacElectroMechanicsProblem<2,1> problem(INCOMPRESSIBLE,
      MONODOMAIN,
      &mesh,
      &mechanics_mesh,
      &cell_factory,
      &problem_defn,
      "icc2d_Du");

      //// set a node to watch electrics and mechanics
      c_vector<double,2> node_to_watch;
      node_to_watch(0) = 1.0;
      node_to_watch(1) = 1.0;
      problem.SetWatchedPosition(node_to_watch);
      problem.SetOutputDeformationGradientsAndStress(100.0); //// timestep in ms

      //// solve
      // problem.SetNoElectricsOutput();
      problem.Initialise();
      problem.Solve();

      //// dev test
      // FileFinder test_output_folder("icc2d/electrics", RelativeTo::ChasteTestOutput);
      // Hdf5ToMeshalyzerConverter<2,2> converter(test_output_folder, "voltage",
      //                                               &mesh, false,
      //                                               HeartConfig::Instance()->GetVisualizerOutputPrecision());

    }
  };
