/*

Copyright (c) 2005-2022, University of Oxford.
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


#include "NashContractionModel.hpp"
#include <math.h>
#include <iostream>
// #include <cmath>

const double NashContractionModel::T0 = 24981000.0; //kPa
const double NashContractionModel::beta = 1.45; // dimensionless
const double NashContractionModel::Ca_50 = 1.95; // millimolar
const double NashContractionModel::h = 2.5; // dimensionless
const double NashContractionModel::Ca_max = 0.00055535; // millimolar
const double NashContractionModel::lambda_val = 1.0; // millimolar


NashContractionModel::NashContractionModel()
: AbstractAlgebraicContractionModel()
{
  // assert(option>=1);
  // assert(option<=3);

  // mOption = option;
  mStretch = 1.0;
}


void NashContractionModel::SetInputParameters(ContractionModelInputParameters& rInputParameters)
{
  assert(rInputParameters.intracellularCalciumConcentration != DOUBLE_UNSET);
  assert(rInputParameters.intracellularCalciumConcentration > 0.0);
  mCalciumI = rInputParameters.intracellularCalciumConcentration;
}

void NashContractionModel::SetIntracellularCalciumConcentration(double calciumConcentration)
{
  assert(calciumConcentration > 0.0);
  mCalciumI = calciumConcentration;
}
void NashContractionModel::SetStretchAndStretchRate(double stretch, double stretchRate)
{
    mStretch = stretch;
}

double NashContractionModel::GetActiveTension()
{
  double tension = (Ca_max * T0 * (1 + beta * (lambda_val - 1)) * pow(mCalciumI, h)) / (1 * (pow(mCalciumI, h) + pow(Ca_50, h)));
  std::cout << "Ca_i is: " << mCalciumI <<std::endl;
  std::cout << "tension is: " << tension <<std::endl;
  return tension;
}
