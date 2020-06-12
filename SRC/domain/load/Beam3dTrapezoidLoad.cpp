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
**   WilsonWu (wilsonwu@wwforu.com)				                      **
**                                                                    **
** ****************************************************************** */                                                                       
                                                                        
// Written: WilsonWu 

// Purpose: This file contains the class implementation Beam3dGenranlPartialLoad.

#include <Beam3dTrapezoidLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dTrapezoidLoad::data(4);

Beam3dTrapezoidLoad::Beam3dTrapezoidLoad(int tag, double py, double pz,
	double aoverl, double boverl, int theElementTag)
	:ElementalLoad(tag, LOAD_TAG_Beam3dTrapezoidLoad, theElementTag),
	Py(py), Pz(pz), aOverL(aoverl), bOverL(boverl)
{

}


Beam3dTrapezoidLoad::Beam3dTrapezoidLoad()
	:ElementalLoad(LOAD_TAG_Beam3dTrapezoidLoad),
	Py(0), Pz(0), aOverL(0), bOverL(0)
{

}


Beam3dTrapezoidLoad::~Beam3dTrapezoidLoad()
{

}


const Vector & Beam3dTrapezoidLoad::getData(int &type, double loadFactor)
{
	type = LOAD_TAG_Beam3dTrapezoidLoad;
	data(0) = Py;
	data(1) = Pz;

	data(2) = aOverL;
	data(3) = bOverL;

	return data;
}


int Beam3dTrapezoidLoad::sendSelf(int commitTag, Channel &theChannel)
{
	int dbTag = this->getDbTag();

	static Vector vectData(6);
	vectData(0) = Py;
	vectData(1) = Pz;

	vectData(2) = aOverL;
	vectData(3) = bOverL;

	vectData(4) = eleTag;
	vectData(5) = this->getTag();

	int result = theChannel.sendVector(dbTag, commitTag, vectData);
	if (result < 0) {
		opserr << "Beam3dTrapezoidLoad::sendSelf - failed to send data\n";
		return result;
	}

	return 0;
}


int Beam3dTrapezoidLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int dbTag = this->getDbTag();

	static Vector vectData(6);

	int result = theChannel.recvVector(dbTag, commitTag, vectData);
	if (result < 0) {
		opserr << "Beam3dTrapezoidLoad::recvSelf - failed to recv data\n";
		return result;

	}
	this->setTag(vectData(5));
	Py = vectData(0);;
	Pz = vectData(1);;

	aOverL = vectData(2);
	bOverL = vectData(3);

	eleTag = (int)vectData(4);

	return 0;
}


void Beam3dTrapezoidLoad::Print(OPS_Stream &s, int flag /*= 0*/)
{
	s << "Beam3dTrapezoidLoad - Reference load" << endln;
	s << "  Transverse (y): " << Py << endln;
	s << "  Transverse (z): " << Pz << endln;
	s << "  Element: " << eleTag << endln;;
}
