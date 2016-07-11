/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:00:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/StaticIntegrator.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/StaticIntegrator.C
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for StaticIntegrator.
// StaticIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// What: "@(#) StaticIntegrator.C, revA"

#include <StaticIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

StaticIntegrator::StaticIntegrator(int clasTag)
:IncrementalIntegrator(clasTag)
{
   
    // for subclasses
}

StaticIntegrator::~StaticIntegrator()
{
}

int
StaticIntegrator::formEleTangent(FE_Element *theEle)
{
   
  // opserr<<" this is form ELE tangent function "<<endln;//Abbas
    // only elements stiffness matrix needed

  if (statusFlag == CURRENT_TANGENT) {
  //   opserr<<" Am I in the current tangent"<<endln;//Abbas
    theEle->zeroTangent();
    theEle->addKtToTang();
    
  } else if (statusFlag == INITIAL_TANGENT) {
    // opserr<<" Am I in the initial tangent"<<endln;//Abbas
    theEle->zeroTangent();
    theEle->addKiToTang();
  } 

    return 0;
}    

int
StaticIntegrator::formEleResidual(FE_Element *theEle)
{
    // only elements residual needed
    theEle->zeroResidual();
    theEle->addRtoResidual();
         return 0;
}    

int
StaticIntegrator::formNodTangent(DOF_Group *theDof)
{
    // should never be called
    opserr << "StaticIntegrator::formNodTangent() -";
    opserr << " this method should never have been called!\n";
    return -1;
}    

int
StaticIntegrator::formNodUnbalance(DOF_Group *theDof)
{
    // only nodes unbalance need be added
    theDof->zeroUnbalance();
    theDof->addPtoUnbalance();
    return 0;
}    

//////////////////////////////////////////Abbas////////////////////

   int
StaticIntegrator::formEleTangentSensitivity(FE_Element *theEle,int gradNumber)
{
 
    if (statusFlag == CURRENT_TANGENT) {
   //  opserr<<" Am I in the current tangent sensitivity"<<endln;//Abbas
    theEle->zeroTangent();
    theEle->addKtToTangSensitivity(gradNumber,1.0);
   
  } else if (statusFlag == INITIAL_TANGENT) {
    theEle->zeroTangent();
    theEle->addKiToTang();
  } 

    return 0;
}    
////////////////////////////////////////////Abbas//////////////////

