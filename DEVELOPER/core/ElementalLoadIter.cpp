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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/ElementalLoadIter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/loadcase/ElementalLoadIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// ElementalLoadIter. ElementalLoadIter is a class for iterating through the 
// elemental loads of a load case.

#include "ElementalLoadIter.h"

#include <ElementalLoad.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// ElementalLoadIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

ElementalLoadIter::ElementalLoadIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}


ElementalLoadIter::~ElementalLoadIter()
{
}    

void
ElementalLoadIter::reset(void)
{
    myIter.reset();
}    


ElementalLoad *
ElementalLoadIter::operator()(void)
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	ElementalLoad *result = (ElementalLoad *)theComponent;
	return result;
    }
}

