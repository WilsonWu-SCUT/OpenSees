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

// Purpose: This file contains the class implementation Beam3dTriangularLoad.

#include <Beam3dTriangularLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dTriangularLoad::data(3);

Beam3dTriangularLoad::Beam3dTriangularLoad(int tag, double py, double pz, double aoverl, int theElementTag)
	:ElementalLoad(tag, LOAD_TAG_Beam3dTriangularLoad, theElementTag),
	Py(py), Pz(pz), aOverL(aoverl)
{

}

Beam3dTriangularLoad::Beam3dTriangularLoad()
	:ElementalLoad(LOAD_TAG_Beam3dTriangularLoad),
	Py(0), Pz(0), aOverL(0)
{

}

Beam3dTriangularLoad::~Beam3dTriangularLoad()
{

}

const Vector & Beam3dTriangularLoad::getData(int &type, double loadFactor)
{
	type = LOAD_TAG_Beam3dTriangularLoad;
	data(0) = Py;
	data(1) = Pz;

	data(3) = aOverL;

	return data;
}

int Beam3dTriangularLoad::sendSelf(int commitTag, Channel &theChannel)
{
	int dbTag = this->getDbTag();

	static Vector vectData(5);
	vectData(0) = Py;
	vectData(1) = Pz;

	vectData(2) = aOverL;

	vectData(3) = eleTag;
	vectData(4) = this->getTag();

	int result = theChannel.sendVector(dbTag, commitTag, vectData);
	if (result < 0) {
		opserr << "Beam3dTriangularLoad::sendSelf - failed to send data\n";
		return result;
	}

	return 0;
}

int Beam3dTriangularLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int dbTag = this->getDbTag();

	static Vector vectData(5);

	int result = theChannel.recvVector(dbTag, commitTag, vectData);
	if (result < 0) {
		opserr << "Beam3dTriangularLoad::recvSelf - failed to recv data\n";
		return result;

	}
	this->setTag(vectData(4));
	Py = vectData(0);
	Pz = vectData(1);

	aOverL = vectData(2);

	eleTag = (int)vectData(3);

	return 0;
}

void Beam3dTriangularLoad::Print(OPS_Stream &s, int flag /*= 0*/)
{
	s << "Beam3dTriangularLoad - Reference load" << endln;
	s << "  Transverse (y): " << Py << endln;
	s << "  Transverse (z): " << Pz << endln;
	s << "  Element: " << eleTag << endln;;
}
