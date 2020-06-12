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

#include <Beam3dGenranlPartialLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dGenranlPartialLoad::data(8);

Beam3dGenranlPartialLoad::Beam3dGenranlPartialLoad(int tag, 
	double pyi, double pzi, double ni, 
	double pyj, double pzj, double nj, 
	double aoverl, double boverl, int theElementTag)
	:ElementalLoad(tag, LOAD_TAG_Beam3dGenranlPartialLoad, theElementTag),
	Pyi(pyi), Pzi(pzi), Ni(ni),
	Pyj(pyj), Pzj(pzj), Nj(nj),
	aOverL(aoverl), bOverL(boverl)
{

}

Beam3dGenranlPartialLoad::Beam3dGenranlPartialLoad()
	:ElementalLoad(LOAD_TAG_Beam3dGenranlPartialLoad),
	Pyi(0), Pzi(0), Ni(0),
	Pyj(0), Pzj(0), Nj(0),
	aOverL(0), bOverL(0)
{

}

Beam3dGenranlPartialLoad::~Beam3dGenranlPartialLoad()
{

}


const Vector & Beam3dGenranlPartialLoad::getData(int &type, double loadFactor)
{
	type = LOAD_TAG_Beam3dGenranlPartialLoad;
	data(0) = Pyi;
	data(1) = Pzi;
	data(2) = Ni;

	data(3) = Pyj;
	data(4) = Pzj;
	data(5) = Nj;

	data(6) = aOverL;
	data(7) = bOverL;

	return data;
}


int Beam3dGenranlPartialLoad::sendSelf(int commitTag, Channel &theChannel)
{
	int dbTag = this->getDbTag();

	static Vector vectData(10);
	vectData(0) = Pyi;
	vectData(1) = Pzi;
	vectData(2) = Ni;

	vectData(3) = Pyj;
	vectData(4) = Pzj;
	vectData(5) = Nj;

	vectData(6) = aOverL;
	vectData(7) = bOverL;

	vectData(8) = eleTag;
	vectData(9) = this->getTag();

	int result = theChannel.sendVector(dbTag, commitTag, vectData);
	if (result < 0) {
		opserr << "Beam3dGenranlPartialLoad::sendSelf - failed to send data\n";
		return result;
	}

	return 0;
}


int Beam3dGenranlPartialLoad::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int dbTag = this->getDbTag();

	static Vector vectData(10);

	int result = theChannel.recvVector(dbTag, commitTag, vectData);
	if (result < 0) {
		opserr << "Beam3dGenranlPartialLoad::recvSelf - failed to recv data\n";
		return result;

	}
	this->setTag(vectData(9));
	Pyi = vectData(0);;
	Pzi = vectData(1);;
	Ni = vectData(2);;

	Pyj = vectData(3);
	Pzj = vectData(4);
	Nj = vectData(5);

	aOverL = vectData(6);
	bOverL = vectData(7);


	eleTag = (int)vectData(8);

	return 0;
}


void Beam3dGenranlPartialLoad::Print(OPS_Stream &s, int flag /*= 0*/)
{
	s << "Beam3dGenranlPartialLoad - Reference load" << endln;
	s << "  Transverse (yi): " << Pyi << endln;
	s << "  Transverse (zi): " << Pzi << endln;
	s << "  Axial (xi):      " << Ni << endln;
	s << "  Element: " << eleTag << endln;;
}
