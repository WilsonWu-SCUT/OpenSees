#ifndef _OPS_Globals_h
#define _OPS_Globals_h

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
                                                                        
// Written: fmk
// Created: 11/99
//
// Description: This file contains global variables used in OpenSees files.
// if you change some of the variables, you must recompile ALL the code.


#define OPS_VERSION "3.2.1"

#ifdef _WIN32
#ifndef _WIN64
#define WIN_ARCH "32-Bit"
#else
#define WIN_ARCH "64-Bit"
#endif
#else
#define WIN_ARCH "64-Bit"
#endif

#define _USING_OpenSees_STREAMS
#include <OPS_Stream.h>
#include <sstream>
#include <iostream>

//extern OPS_Stream &opserr;
extern OPS_Stream* opserrPtr;

#define _SAP

#ifdef _SAP
//#pragma warning(disable: 4101)
//#pragma warning(disable: 4244)
//#pragma warning(disable: 4267)
static std::ostringstream osbuff;
#define opserr osbuff
#define endln "\n"
#else

#define opserr (*opserrPtr)
#define endln "\n"
#endif // SAP

#include <string.h>
#include <stdlib.h>

enum NodeResponseType { Disp = 1, Vel = 2, Accel =3, IncrDisp =4, IncrDeltaDisp =5, Reaction =6, Unbalance =7, RayleighForces =8};

#ifdef _TCL85
#define TCL_Char const char
#elif _TCL84
#define TCL_Char const char
#else
#define TCL_Char char
#endif

class Domain;
class Element;

#define MAX_FILENAMELENGTH 50

extern double   ops_Dt;                // current delta T for current domain doing an update
// extern double  *ops_Gravity;        // gravity factors for current domain undergoing an update
extern int ops_Creep;
extern Domain  *ops_TheActiveDomain;   // current domain undergoing an update
extern Element *ops_TheActiveElement;  // current element undergoing an update

// global variable for initial state analysis
// added: Chris McGann, University of Washington
extern bool  ops_InitialStateAnalysis;

#define OPS_DISPLAYMODE_MATERIAL_TAG 2
#define OPS_DISPLAYMODE_ELEMENT_CLASS 3
#define OPS_DISPLAYMODE_STRESS 5
#define OPS_DISPLAYMODE_STRAIN 7
#define OPS_DISPLAYMODE_AXIAL 11
#define OPS_PRINT_CURRENTSTATE 0
#define OPS_PRINT_PRINTMODEL_SECTION  1
#define OPS_PRINT_PRINTMODEL_MATERIAL 2
#define OPS_PRINT_PRINTMODEL_JSON   25000

#endif
