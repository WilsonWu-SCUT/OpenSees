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
                                                                        

                                                                        
#ifndef Beam3dGenranlPartialLoad_h
#define Beam3dGenranlPartialLoad_h

// Written: WilsonWu 

// Purpose: This file contains the class definition for Beam3dGenranlPartialLoad.

#include <ElementalLoad.h>

class Beam3dGenranlPartialLoad : public ElementalLoad
{
public:
	Beam3dGenranlPartialLoad(int tag, double pyi, double pzi, double ni, double pyj, double pzj, double nj,
		double aoverl, double boverl, int theElementTag);
	Beam3dGenranlPartialLoad();
	~Beam3dGenranlPartialLoad();

	const Vector &getData(int &type, double loadFactor);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	double Pyi;	// magnitude of the transverse load at I-side
	double Pzi;	// magnitude of the transverse load at I-side
	double Ni;	// magnitude of the axial load at I-side

	double Pyj;	// magnitude of the transverse load at J-side
	double Pzj;	// magnitude of the transverse load at J-side
	double Nj;	// magnitude of the axial load at J-side

	double aOverL; //relative distance (x/L) along length from end 1 of element
	double bOverL; //relative distance (x/L) along length from THE LOAD START POINT

	static Vector data;
};

#endif

