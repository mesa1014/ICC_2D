#include "ICCFactory.hpp"
#include <ChasteEllipsoid.hpp>

template<unsigned DIM>
AbstractCardiacCell* ICCFactory<DIM>::CreateCardiacCellForTissueNode(Node<DIM>* pNode)
{
  Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(this->mpSolver, this->mpZeroStimulus);
  double x = pNode->rGetLocation()[0];
  double y = pNode->rGetLocation()[1];
  ChastePoint<2> centre(1.0,2.8);
  ChastePoint<2> radii (1.5,1.5); // (1.5, 1.5)
  ChasteEllipsoid<2> ellipseRegion(centre, radii);
  ChastePoint<2> myPoint(x, y);
  // double ca = cell->GetIntracellularCalciumConcentration();
  // std::cout << "Ca: " << ca << std::endl;
  cell->SetParameter("eta", 0.045);
  if(ellipseRegion.DoesContain(myPoint))
  {
    cell->SetParameter("eta", 0.037);
    std::cout<<"I'm inside region!" << std::endl;
  }
  return cell;

// unsigned index = pNode->GetIndex();
// double y = pNode->GetPoint()[1];
// if(setICCNode.find(index) != setICCNode.end())
// {
//   CellICCNeuralCalib_SMC_noCasesFromCellML* cell = new CellICCNeuralCalib_SMC_noCasesFromCellML(this->mpSolver, this->mpZeroStimulus);
//
//   cell->SetParameter("E_K_ICC", -70.0); //-4.0*y);
//   std::cout <<"I'm inside region!" <<std::endl;
//   return cell;
// }

// return new DummyDerivedCa(this->mpSolver, this->mpZeroStimulus);

}

// Explicit instantiation
template class ICCFactory<1>;
template class ICCFactory<2>;
template class ICCFactory<3>;
