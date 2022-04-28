/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
* Neither the name of the University of Oxford nor the names of its
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
*
*  Chaste tutorial - this page gets automatically changed to a wiki page
*  DO NOT remove the comments below, and if the code has to be changed in
*  order to run, please check the comments are still accurate
*
*
*/
#ifndef TESTCARDIACELECTROMECHANICSTUTORIALV3_HPP_
#define TESTCARDIACELECTROMECHANICSTUTORIALV3_HPP_

/*
* = Cardiac Electro-mechanical Problems =
*
* == Introduction ==
*
* The tutorial explains how electro-mechanics problems can be solved in Chaste. The reader should certainly read
* the electro-physiological tutorials before this tutorial, and really they should have also had a look at
* the tutorial(s) on solving general solid mechanics problems.
*
* The equations of cardiac electro-mechanics are written down in Section 4.2 of the PDF on equations and
* finite element implementations in ChasteGuides -> Miscellaneous information. '''Note:''' By default we do
* not solve these full equations: the mechanics information is not coupled back to electrics, ie by default
* the conductivities do not depend on deformation, and cell models do not get affected by stretch.
* This has to be switched on if required, as will be described further below.
*
* Before going to the code, we list the sub-models/parameters that need to be set, or can be varied,
* in electro-mechanical problems. The last five of the below are mechanics-specific.
*  * The geometry (see note 1 below)
*  * The region electrically stimulated
*  * The cell model
*  * Electro-physiological parameters (conductivity, capacitance, surface-area-to-volume ratio)
*  * Electro-physiological timesteps: ode and pde (but not printing timestep) (see note 2 below)
*  * Fibre directions (and maybe sheet/normal directions) (see note 3 below)
*  * The part of the boundary that has displacement boundary conditions
*  * Any pressure or traction boundary conditions
*  * The contraction model [the model which takes in electrical variables (voltage or calcium typically), and
*  returns cellular active tension]
*  * Whether the tissue should be treated as compressible or incompressible. (Although likely technically
*  incompressible at appropriate scales, cardiac tissue is often treated as compressible due to blood
*  squeezed out of the coronary vessels during contraction).
*  * The material law [the strain-energy function]
*  * Mechanics timesteps: mechanics update timestep, contraction model ode timestep. (see note 4 below)
*
* Notes:
*  * ''Meshes:'' Two meshes for the geometry are required, one for the electrics solve and one for the mechanics.
* The mechanics mesh would ideally be coarser but any two meshes are technically possible. The meshes should
* ideally both cover exactly the same geometry (ie either mesh being contained in the other), but the meshes
* not completely overlapping is allowed - some extrapolation of quantities will then occur.
*  * ''The electro-physiology printing timestep:'' This is not used in electro-mechanics problems; output is
* instead written after every mechanics solve, so effectively the mechanics update timestep is equal to
* the printing timestep.
*  * ''Fibres:'' In electro-physiological simulations the fibre direction is in the X-direction
* by default, but if isotropic conductivities are used the fibre direction won't be used. In mechanics
* solves, the fibres will always be used as it determines the direction of contraction. It defaults to the
* X-direction, so this is the direction the tissue will contract, unless a fibre file is given.
* If the material law is transversely isotropic, the problem is independent of sheet & normal directions.
* If the material law is anisotropic, the problem is dependent of sheet & normal directions.
*  * ''Timesteps:'' Should-divide rules are: (a) ode_timestep should-divide pde_timestep should-divide
*  mechanics_update_timestep and (b) contraction_model_ode_timestep should-divide mechanics_update_timestep.
*
* '''Another important note:''' mechanics problems are not currently implemented to scale in parallel yet. This
* is work in progress.
*
* The basic includes are */
#include <cxxtest/TestSuite.h>
#include "PlaneStimulusCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CardiacElectroMechProbRegularGeom.hpp"
#include "CardiacElectroMechanicsProblem.hpp"
#include "LuoRudy1991.hpp"
/* Some other includes that are used */
#include "NonlinearElasticityTools.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "NobleVargheseKohlNoble1998WithSac.hpp"
#include "ZeroStimulusCellFactory.hpp"
#include "FileComparison.hpp"
#include "FileFinder.hpp"

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"

#include <iostream>
#include <string>
#include <filesystem>
#include <unistd.h>

/*
* HOW_TO_TAG Cardiac/Electro-mechanics
* Run basic electro-mechanics simulations; specify different models, boundary conditions, fibres
*/

/*
* == IMPORTANT: using HYPRE ==
*class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
PointStimulus2dCellFactory()
: AbstractCardiacCellFactory<2>(),
mpStimulus(new SimpleStimulus(-5e5, 0.5))
{
}

AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
{
double x = pNode->rGetLocation()[0];
double y = pNode->rGetLocation()[1];
if (fabs(x)<0.02+1e-6 && y<-0.4+1e-6) // stimulating small region
{
return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
}
else
{
return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
}
}
};
* Mechanics solves being nonlinear are expensive, so it is recommended you also use a `Release` build type for `cmake`
* (or `build=GccOpt_ndebug` when running on the old `scons` build system)
* on larger problems. Also:
*
* Mechanics solves involve solving a nonlinear system, which is broken down into a sequence of linear solves.
* When running '''incompressible''' problems '''in 3D, or with more elements than in the first test below''',
* it is vital to change the linear solver to use HYPRE, an algebraic multigrid solver.
* Without HYRPE, the linear solve (i) may become very very slow; or
* (ii) may not converge, in which case the nonlinear solve will (probably) not converge. See the comments on using
* HYPRE in the first solid mechanics tutorial.
*
* == Simple 2d test ==
*
* This test shows how to use the `CardiacElectroMechProbRegularGeom` class, which
* inherits from a more general class, `CardiacElectroMechanicsProblem`, and
* sets up a square or cubic geometry for you. Using
* `CardiacElectroMechProbRegularGeom` is not really recommended, as the functionality
* it allows is very limited - it is better to use `CardiacElectroMechanicsProblem`, which
* is shown in the following tests. We use `CardiacElectroMechProbRegularGeom`
* in this first tutorial just to illustrate a simulation with a few lines (four!) of code.
*/
class PointStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
  boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
  PointStimulus2dCellFactory()
  : AbstractCardiacCellFactory<2>(),
  mpStimulus(new SimpleStimulus(-5e5, 0.5))
  {
  }

  AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
  {
    double x = pNode->rGetLocation()[0];
    double y = pNode->rGetLocation()[1];
    if (fabs(x)<0.02+1e-6 && y<-0.4+1e-6) // stimulating small region
    {
      return new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
    }
    else
    {
      return new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
    }
  }
};

class TestCardiacElectroMechanicsTutorialV3 : public CxxTest::TestSuite
{
public:
  void TestCardiacElectroMechanicsExampleAgainV2()
  {
    ZeroStimulusCellFactory<CML_noble_varghese_kohl_noble_1998_basic_with_sac, 2> cell_factory;

    TetrahedralMesh<2,2> electrics_mesh;
    QuadraticMesh<2> mechanics_mesh;

    std::cout << "Before read mesh ..." << std::endl;
    // could (should?) use finer electrics mesh, but keeping electrics simulation time down
    TrianglesMeshReader<2,2> reader1("/home/chaste/projects/mesh/MeshNetwork-2D.1");
    electrics_mesh.ConstructFromMeshReader(reader1);
    std::cout << "After read electrics mesh ..." << std::endl;

    TrianglesMeshReader<2,2> reader2("/home/chaste/projects/mesh/MeshNetwork-2D.2",2 /*quadratic elements*/);
    mechanics_mesh.ConstructFromMeshReader(reader2);
    std::cout << "After read mechanics mesh ..." << std::endl;

    std::vector<unsigned> fixed_nodes;
    std::vector<c_vector<double,2> > fixed_node_locations;

    fixed_nodes.push_back(0);
    fixed_node_locations.push_back(zero_vector<double>(2));

    std::cout << "Number of Nodes:" << mechanics_mesh.GetNumNodes() << std::endl;

    for (unsigned i=1; i<mechanics_mesh.GetNumNodes(); i++)
    {
      double X = mechanics_mesh.GetNode(i)->rGetLocation()[0];
      std::cout << "Node x:" << X << std::endl;
      if (fabs(X) < 1e-6) // ie, if X==0
      {
        c_vector<double,2> new_position;
        new_position(0) = 0.0;
        new_position(1) = ElectroMechanicsProblemDefinition<2>::FREE;

        fixed_nodes.push_back(i);
        fixed_node_locations.push_back(new_position);
      }
    }
    std::cout << "First loop passed!" << std::endl;
    std::vector<BoundaryElement<1,2>*> boundary_elems;
    std::vector<c_vector<double,2> > tractions;

    c_vector<double,2> traction;

    for (TetrahedralMesh<2,2>::BoundaryElementIterator iter = mechanics_mesh.GetBoundaryElementIteratorBegin();
    iter != mechanics_mesh.GetBoundaryElementIteratorEnd();
    ++iter)
    {
      if (fabs((*iter)->CalculateCentroid()[1] - 0.0) < 1e-6) // if Y=0
      {
        BoundaryElement<1,2>* p_element = *iter;
        boundary_elems.push_back(p_element);

        traction(0) =  0.0; // kPa, since the contraction model and material law use kPa for stiffnesses
        traction(1) =  1.5; // kPa, since the contraction model and material law use kPa for stiffnesses
        tractions.push_back(traction);
      }
      if (fabs((*iter)->CalculateCentroid()[1] - 0.1) < 1e-6) // if Y=0.1
      {
        BoundaryElement<1,2>* p_element = *iter;
        boundary_elems.push_back(p_element);

        traction(0) =  0.0;
        traction(1) = -1.5;
        tractions.push_back(traction);
      }
    }
    std::cout << "Second loop passed!" << std::endl;

    ElectroMechanicsProblemDefinition<2> problem_defn(mechanics_mesh);
    problem_defn.SetContractionModel(NHS,0.01/*contraction model ODE timestep*/);
    problem_defn.SetUseDefaultCardiacMaterialLaw(INCOMPRESSIBLE);
    problem_defn.SetMechanicsSolveTimestep(1.0);
    problem_defn.SetFixedNodes(fixed_nodes, fixed_node_locations);
    problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

    problem_defn.SetDeformationAffectsElectrophysiology(false /*deformation affects conductivity*/, true /*deformation affects cell models*/);

    HeartConfig::Instance()->SetSimulationDuration(40.0);

    CardiacElectroMechanicsProblem<2,2> problem(INCOMPRESSIBLE,
                                                BIDOMAIN,
                                                &electrics_mesh,
                                                &mechanics_mesh,
                                                &cell_factory,
                                                &problem_defn,
                                                "TestCardiacElectroMechanicsWithMefV2");
      // problem_defn.Validate();
      problem.Solve();

      // FileFinder test_output_folder("TestCardiacElectroMechanicsWithMef/electrics", RelativeTo::ChasteTestOutput);
      // Hdf5ToMeshalyzerConverter<2,2> converter(test_output_folder, "voltage",
      // &electrics_mesh, false,
      // HeartConfig::Instance()->GetVisualizerOutputPrecision());
      //
      // Hdf5DataReader reader("TestCardiacElectroMechanicsWithMef/electrics", "voltage");
      // unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
      // Vec voltage = PetscTools::CreateVec(electrics_mesh.GetNumNodes());
      // reader.GetVariableOverNodes(voltage, "V", num_timesteps-1);
      // ReplicatableVector voltage_repl(voltage);
      // for (unsigned i=0; i<voltage_repl.GetSize(); i++)
      // {
      //   TS_ASSERT_DELTA(voltage_repl[i], -81.9080, 1e-3);
      // }
      // PetscTools::Destroy(voltage);
    }
  };

  #endif /*TESTCARDIACELECTROMECHANICSTUTORIAL_HPP_*/
