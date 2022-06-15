#include "ICCFactory.hpp"
#include <ChasteEllipsoid.hpp>

template<unsigned DIM>
AbstractCardiacCell* ICCFactory<DIM>::CreateCardiacCellForTissueNode(Node<DIM>* pNode)
{
  // CellDu2013_neural_sensFromCellML* cell = new CellDu2013_neural_sensFromCellML(this->mpSolver, this->mpZeroStimulus);
  // double x = pNode->rGetLocation()[0];
  // double y = pNode->rGetLocation()[1];
  // ChastePoint<2> centre(1.0,2.0);
  // ChastePoint<2> radii (0.33,0.33);
  // ChasteEllipsoid<2> ellipseRegion(centre, radii);
  // ChastePoint<2> myPoint(x, y);
  // // double ca = cell->GetIntracellularCalciumConcentration();
  // // std::cout << "Ca: " << ca << std::endl;
  // if(ellipseRegion.DoesContain(myPoint))
  // {
  //   cell->SetParameter("E_K", -55.0);
  //   std::cout<<"I'm inside region!" << std::endl;
  // }
  // return cell;

  unsigned index = pNode->GetIndex();
  double y = pNode->GetPoint()[1];
  if(setICCNode.find(index) != setICCNode.end())
  {
    CellICCNeuralCalib_SMC_noCasesFromCellML* cell = new CellICCNeuralCalib_SMC_noCasesFromCellML(this->mpSolver, this->mpZeroStimulus);

    cell->SetParameter("E_K_ICC", -70.0-4.0*y);
    std::cout <<"I'm inside region!" <<std::endl;
    return cell;
  }

  return new DummyDerivedCa(this->mpSolver, this->mpZeroStimulus);

}

// Explicit instantiation
template class ICCFactory<1>;
template class ICCFactory<2>;
template class ICCFactory<3>;
