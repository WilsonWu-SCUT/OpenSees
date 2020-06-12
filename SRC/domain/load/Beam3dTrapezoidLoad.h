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
                                                                        

                                                                        
#ifndef Beam3dTrapezoidLoad_h
#define Beam3dTrapezoidLoad_h

// Written: WilsonWu 

// Purpose: This file contains the class definition for Beam3dTrapezoidLoad.

#include <ElementalLoad.h>

class Beam3dTrapezoidLoad : public ElementalLoad
{
public:
	Beam3dTrapezoidLoad(int tag, double py, double pz,
		double aoverl, double boverl, int theElementTag);
	Beam3dTrapezoidLoad();
	~Beam3dTrapezoidLoad();

	const Vector &getData(int &type, double loadFactor);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	double Py;	// magnitude of the transverse load at I-side
	double Pz;	// magnitude of the transverse load at I-side

	double aOverL; //relative distance (x/L) along length from end 1 of element
	double bOverL; //relative distance (x/L) along length from end 1 of element

	static Vector data;
};

#endif

