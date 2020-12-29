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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
                                                                        
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam3d.
// ElasticBeam3d is a plane frame member.

#ifndef AIDMBeamColumn_h
#define AIDMBeamColumn_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <vector>

class Channel;
class Information;
class CrdTransf;
class Response;
class Renderer;
class AIDMMaterial;
class UniaxialMaterial;



class AIDMBeamColumn : public Element
{
  public:
    AIDMBeamColumn();
    AIDMBeamColumn(int tag, double A, double E, double G,
        double Jx, double Iy, double Iz,
        int Nd1, int Nd2, CrdTransf& theTransf, int numAidms, UniaxialMaterial** aidms);
    AIDMBeamColumn(int tag, double A, double E, double G,
        double Jx, double Iy, double Iz,
        int Nd1, int Nd2, CrdTransf& theTransf, const std::vector<int> tag_vec);

    AIDMBeamColumn(int tag, double A, double E, double G,
		  double Jx, double Iy, double Iz,
          int Nd1, int Nd2, CrdTransf &theTransf,
          double rho = 0.0, int cMass = 0,
          int sectionTag = 0);
    AIDMBeamColumn(int tag, double A, double E, double G,
		double Jx, double Iy, double Iz,
		int Nd1, int Nd2, CrdTransf& theTransf, 
        int Release, int monitorPointNum,
		double rho = 0.0, int cMass = 0,
		int sectionTag = 0);

    ~AIDMBeamColumn();

    const char *getClassType(void) const {return "AIDMBeamColumn";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);
    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);
    
    int update(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector& getLocalResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);    
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes = 0, int numModes = 0);

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

private:
    void initialAIDMBeamColumn(int Nd1, int Nd2, CrdTransf& theTransf);
	void setStiffMatrix(const double& length);
	void setBasicForce(const double& length, const Vector& v);
	int addPointLoad(const double& N, const double& Py, const double& Pz,
		const double& aOverL, const double& length);
	void addGeneralPartialLoad(const double& Ni, const double& Nj, const double& Pyi,
		const double& Pyj, const double& Pzi, const double& Pzj,
		const double& aOverL, const double& bOverL, const double& length);
    bool isAvailabelAIDM() const;

  private:
      int iAIDMTag;
      int jAIDMTag;
      bool isGravityConst;
      double TLoadFactor;
      double CLoadFactor;

      bool isInitial;

  private:
    double A,E,G,Jx,Iy,Iz;

    AIDMMaterial** AIDMs;
    int numAIDMs;

    static Matrix K;
    static Vector P;
    Vector Q;
    
    static Matrix kb;

    Vector q;
    double q0[5];  // Fixed end forces in basic system (no torsion)
    double p0[5];  // Reactions in basic system (no torsion)
 
    Node *theNodes[2];

    ID  connectedExternalNodes;    

    CrdTransf *theCoordTransf;
};

#endif
