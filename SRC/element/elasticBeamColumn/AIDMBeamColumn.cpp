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


// File: ~/model/ElasticBeam3d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam3d.
// ElasticBeam3d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <AIDMBeamColumn.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Renderer.h>
#include <AIDMMaterial.h>
#include <ID.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <elementAPI.h>

//#include <sstream>


void AIDMBeamColumn::setStiffMatrix(const double& L)
{
	double oneOverL = 1.0 / L;
	double EoverL = E * oneOverL;
	double EAoverL = A * EoverL;			// EA/L

	double EIzoverL2 = 2.0 * Iz * EoverL;		// 2EIz/L
	double EIzoverL3 = 3.0 * Iz * EoverL;		// 3EIz/L
	double EIzoverL4 = 2.0 * EIzoverL2;		// 4EIz/L

	double EIyoverL2 = 2.0 * Iy * EoverL;		// 2EIy/L
	double EIyoverL3 = 3.0 * Iy * EoverL;		// 3EIy/L
	double EIyoverL4 = 2.0 * EIyoverL2;		// 4EIy/L

	double GJoverL = G * Jx * oneOverL;         // GJ/L

	kb.Zero();

	kb(0, 0) = EAoverL;
	kb(5, 5) = GJoverL;
	kb(3, 3) = AIDMs[0]->getTangent();
	kb(4, 4) = AIDMs[1]->getTangent();

	if (iAIDMTag == -1)
	{
		//all release
		if (jAIDMTag == -1)
		{
			// not thing to do
		}
		//i release j elastic
		else if (jAIDMTag == -2)
		{
			kb(2, 2) = EIzoverL3;
			kb(4, 4) = EIyoverL3;
		}
		//i release j aidm
		else
		{
			kb(2, 2) = EIzoverL3;
		}
	}
	else if (iAIDMTag == -2)
	{
		// i elastic j release
		if (jAIDMTag == -1)
		{
			kb(1, 1) = EIzoverL3;
		}
		//ij elastic
		else if (jAIDMTag == -2)
		{
			kb(1, 1) = kb(2, 2) = EIzoverL3 * 2;
		}
		//i elastic j aidm
		else
		{
			kb(1, 1) = kb(2, 2) = EIzoverL3 * 2;
		}
	}
	else
	{
		//Release
		if (jAIDMTag == -1)
		{
			kb(1, 1) = EIzoverL3;
		}
		//Elastic
		else if (jAIDMTag == -2)
		{
			kb(1, 1) = kb(2, 2) = EIzoverL3 * 2;
		}
		else
		{
			kb(1, 1) = kb(2, 2) = EIzoverL3 * 2;
		}
	}

}

void AIDMBeamColumn::setBasicForce(const double& L, const Vector& v)
{
	double oneOverL = 1.0 / L;
	double EoverL = E * oneOverL;
	double EAoverL = A * EoverL;			// EA/L
	double EIzoverL2 = 2.0 * Iz * EoverL;		// 2EIz/L
	double EIzoverL3 = 3.0 * Iz * EoverL;		// 3EIz/L
	double EIzoverL4 = 2.0 * EIzoverL2;		// 4EIz/L

	double EIyoverL2 = 2.0 * Iy * EoverL;		// 2EIy/L
	double EIyoverL3 = 3.0 * Iy * EoverL;		// 3EIy/L
	double EIyoverL4 = 2.0 * EIyoverL2;		// 4EIy/L

	double GJoverL = G * Jx * oneOverL;         // GJ/L

	// Zero for integration
	q.Zero();

	q(0) = EAoverL * v(0);
	q(5) = GJoverL * v(5);
	q(3) = AIDMs[0]->getStress();
	q(4) = -AIDMs[1]->getStress();

	if (iAIDMTag == -1)
	{
		q(1) = 0; q(3) = 0;
		//all release
		if (jAIDMTag == -1)
		{
			q(2) = 0; q(4) = 0;
		}
		//i release j elastic
		else if (jAIDMTag == -2)
		{
			q(2) = EIzoverL3 * v(2);
			q(4) = EIyoverL3 * v(4);
		}
		//i release j aidm
		else
		{
			q(2) = EIzoverL3 * v(2);
		}
	}
	else if (iAIDMTag == -2)
	{
		// i elastic j release
		if (jAIDMTag == -1)
		{
			q(2) = 0; q(4) = 0;
			q(1) = EIzoverL3 * v(1);
		}
		//ij elastic
		else if (jAIDMTag == -2)
		{
			q(1) = EIzoverL3 * 2 * v(1);
			q(2) = EIzoverL3 * 2 * v(2);
		}
		//i elastic j aidm
		else
		{
			q(1) = EIzoverL3 * 2 * v(1);
			q(2) = EIzoverL3 * 2 * v(2);
		}
	}
	else
	{
		//i stiffness
		//Release
		if (jAIDMTag == -1)
		{
			q(2) = 0; q(4) = 0;
			q(1) = EIzoverL3 * v(1);
		}
		//Elastic
		else if (jAIDMTag == -2)
		{
			q(1) = EIzoverL3 * 2 * v(1);
			q(2) = EIzoverL3 * 2 * v(2);
		}
		else
		{
			q(1) = EIzoverL3 * 2 * v(1);
			q(2) = EIzoverL3 * 2 * v(2);
		}
	}

	q(0) += q0[0];
	q(1) += q0[1];
	q(2) += q0[2];
	q(3) += q0[3];
	q(4) += q0[4];
}

int AIDMBeamColumn::addPointLoad(const double& N, const double& Py, const double& Pz, const double& aOverL, const double& L)
{
	if (aOverL < 0.0 || aOverL > 1.0)
		return 0;

	double a = aOverL * L;
	double b = L - a;

	// Reactions in basic system
	p0[0] -= N;
	double V1, V2;
	V1 = Py * (1.0 - aOverL);
	V2 = Py * aOverL;
	p0[1] -= V1;
	p0[2] -= V2;
	V1 = Pz * (1.0 - aOverL);
	V2 = Pz * aOverL;
	p0[3] -= V1;
	p0[4] -= V2;

	double L2 = 1.0 / (L * L);
	double a2 = a * a;
	double b2 = b * b;

	// Fixed end forces in basic system
	q0[0] -= N * aOverL;
	double M1 = 0;
	double M2 = 0;

	if (iAIDMTag != -1 && jAIDMTag != -1) {
		M1 = -a * b2 * Py * L2;
		M2 = a2 * b * Py * L2;
		q0[1] += M1;
		q0[2] += M2;
		M1 = -a * b2 * Pz * L2;
		M2 = a2 * b * Pz * L2;
		q0[3] -= M1;
		q0[4] -= M2;
	}
	else if (jAIDMTag != -1) {
		M2 = 0.5 * Py * a * b * L2 * (a + L);
		q0[2] += M2;
		M2 = 0.5 * Pz * a * b * L2 * (a + L);
		q0[4] -= M2;
	}
	else if (iAIDMTag != -1) {
		M1 = -0.5 * Py * a * b * L2 * (b + L);
		q0[1] += M1;
		M1 = -0.5 * Pz * a * b * L2 * (b + L);
		q0[3] -= M1;
	}
	return 0;
}

void AIDMBeamColumn::addGeneralPartialLoad(const double& Ni, const double& Nj, const double& Pyi, const double& Pyj,
	const double& Pzi, const double& Pzj, const double& aOverL, const double& bOverL, const double& L)
{
	//拆分基本长度
	double deltaLengthFactor = 1E-6;
	double lengthFactor = deltaLengthFactor * L;
	//调换位置
	double minOverL = aOverL <= bOverL ? aOverL : bOverL;
	double maxOverL = aOverL <= bOverL ? bOverL : aOverL;
	//距离
	double distOverL = maxOverL - minOverL;
	//从左到右开始遍历
	for (double overL = deltaLengthFactor / 2; overL <= distOverL; overL += deltaLengthFactor)
	{
		double Py = ((Pyj - Pyi) / distOverL * overL + Pyi) * lengthFactor;
		double Pz = ((Pzj - Pzi) / distOverL * overL + Pzi) * lengthFactor;
		double N = ((Nj - Ni) / distOverL * overL + Ni) * lengthFactor;
		double targetX = minOverL + overL;
		//添加节点荷载
		addPointLoad(N, Py, Pz, targetX, L);
	}
}

Matrix AIDMBeamColumn::K(12, 12);
Vector AIDMBeamColumn::P(12);
Matrix AIDMBeamColumn::kb(6, 6);

void* OPS_AIDMBeamColumn(void)
{
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 10 && numArgs != 5) {
		opserr << "insufficient arguments:eleTag,iNode,jNode,A,E,G,J,Iy,Iz,transfTag\n";
		return 0;
	}

	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if (ndm != 3 || ndf != 6) {
		opserr << "ndm must be 3 and ndf must be 6\n";
		return 0;
	}

	// inputs: 
	int iData[3];
	int numData = 3;
	if (OPS_GetIntInput(&numData, &iData[0]) < 0) return 0;

	SectionForceDeformation* theSection = 0;
	CrdTransf* theTrans = 0;
	double data[6];
	double rigidEnd[2];
	int transfTag, secTag;

	//Section stiffness
	numData = 6;
	if (OPS_GetDoubleInput(&numData, &data[0]) < 0) return 0;
	//TransfTag
	numData = 1;
	if (OPS_GetIntInput(&numData, &transfTag) < 0) return 0;
	theTrans = OPS_getCrdTransf(transfTag);
	if (theTrans == 0) {
		opserr << "no CrdTransf is found\n";
		return 0;
	}
	//AIDMS
	int numAIDMS = 2;
	std::vector<int> vec = {};
	int aidmTags[2];
	if (OPS_GetIntInput(&numAIDMS, &aidmTags[0]) < 0) return 0;
	for (int i = 0; i < numAIDMS; i++)
	{
		vec.push_back(aidmTags[i]);
	}
	//RigidEnds
	numData = 2;
	if (OPS_GetDoubleInput(&numData, &rigidEnd[0]) < 0)
	{
		rigidEnd[0] = 0;
		rigidEnd[1] = 0;
	}
	return new AIDMBeamColumn(iData[0], data[0], data[1], data[2], data[3], data[4],
		data[5], iData[1], iData[2], *theTrans, vec, rigidEnd[0], rigidEnd[1]);
}

AIDMBeamColumn::AIDMBeamColumn(int tag, double A, double E, double G,
	double Jx, double Iy, double Iz,
	int Nd1, int Nd2, CrdTransf& theTransf, int numAidms, UniaxialMaterial** aidms)
	:Element(tag, ELE_TAG_AIDMBeamColumn),
	A(A), E(E), G(G), Jx(Jx), Iy(Iy), Iz(Iz),
	Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0),
	numAIDMs(numAidms), iAIDMTag(0), jAIDMTag(0), isGravityConst(false), CLoadFactor(0.0), TLoadFactor(0.0), isInitial(false)
{
	// allocate memory for numMaterials1d uniaxial material models
	AIDMs = new AIDMMaterial * [numAidms];

	this->initialAIDMBeamColumn(Nd1, Nd2, theTransf, 0 , 0);

	// get a copy of the material objects and check we obtained a valid copy
	for (int i = 0; i < numAIDMs; i++) {
		AIDMs[i] = dynamic_cast<AIDMMaterial*>(aidms[i]->getCopy());
		//AIDMs[i] = (aidms[i]->getCopy());
		if (AIDMs[i] == 0) {
			opserr << "AIDMBeamColumn::AIDMBeamColumn - failed to get a copy of material " << aidms[i]->getTag() << endln;
			exit(-1);
		}
	}
}

AIDMBeamColumn::AIDMBeamColumn(int tag, double A, double E, double G,
	double Jx, double Iy, double Iz,
	int Nd1, int Nd2, CrdTransf& theTransf, const std::vector<int> tag_vec)
	:Element(tag, ELE_TAG_AIDMBeamColumn),
	A(A), E(E), G(G), Jx(Jx), Iy(Iy), Iz(Iz),
	Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0),
	numAIDMs(tag_vec.size()), isGravityConst(false), CLoadFactor(0.0), TLoadFactor(0.0), isInitial(false)
{
	// allocate memory for numMaterials1d uniaxial material models
	AIDMs = new AIDMMaterial * [numAIDMs];
	iAIDMTag = tag_vec[0]; jAIDMTag = tag_vec[1];
//	iAIDMTag = -2; jAIDMTag = -2;

	this->initialAIDMBeamColumn(Nd1, Nd2, theTransf, 0, 0);

	// get a copy of the material objects and check we obtained a valid copy
	for (int i = 0; i < numAIDMs; i++) {
		if (tag_vec[i] < 0)
		{
			AIDMs[i] = new AIDMMaterial();
		}
		else
		{
			AIDMs[i] = dynamic_cast<AIDMMaterial*>(OPS_getUniaxialMaterial(tag_vec[i])->getCopy());
			if (AIDMs[i] == 0)
			{
				opserr << "WARNING no material " << tag_vec[i] <<
					"exitsts - element AIDMBeamColumn\n";
				exit(-1);
			}
		}
	}
	
}

AIDMBeamColumn::AIDMBeamColumn(int tag, double A, double E, double G,
	double Jx, double Iy, double Iz,
	int Nd1, int Nd2, CrdTransf& theTransf, const std::vector<int> tag_vec, const double& rigidILength, const double& rigidJLength)
	:Element(tag, ELE_TAG_AIDMBeamColumn),
	A(A), E(E), G(G), Jx(Jx), Iy(Iy), Iz(Iz),
	Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0),
	numAIDMs(tag_vec.size()), isGravityConst(false), CLoadFactor(0.0), TLoadFactor(0.0), isInitial(false)
{
	// allocate memory for numMaterials1d uniaxial material models
	AIDMs = new AIDMMaterial * [numAIDMs];
	iAIDMTag = tag_vec[0]; jAIDMTag = tag_vec[1];

	this->initialAIDMBeamColumn(Nd1, Nd2, theTransf, rigidILength, rigidJLength);

	// get a copy of the material objects and check we obtained a valid copy
	for (int i = 0; i < numAIDMs; i++) {
		if (tag_vec[i] < 0)
		{
			AIDMs[i] = new AIDMMaterial();
		}
		else
		{
			AIDMs[i] = dynamic_cast<AIDMMaterial*>(OPS_getUniaxialMaterial(tag_vec[i])->getCopy());
			if (AIDMs[i] == 0)
			{
				opserr << "WARNING no material " << tag_vec[i] <<
					"exitsts - element AIDMBeamColumn\n";
				exit(-1);
			}
		}
	}

}

AIDMBeamColumn::AIDMBeamColumn()
	:Element(0, ELE_TAG_AIDMBeamColumn),
	A(0.0), E(0.0), G(0.0), Jx(0.0), Iy(0.0), Iz(0.0),
	Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0), initialLength(0)
{
	// does nothing
	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;
	q0[3] = 0.0;
	q0[4] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;
	p0[3] = 0.0;
	p0[4] = 0.0;

	// set node pointers to NULL
	for (int i = 0; i < 2; i++)
		theNodes[i] = 0;
}

AIDMBeamColumn::~AIDMBeamColumn()
{
	if (theCoordTransf)
		delete theCoordTransf;
}

void AIDMBeamColumn::initialAIDMBeamColumn(int Nd1, int Nd2, CrdTransf& theTransf, const double& rigidILength, const double& rigidJLength)
{
	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;

	theCoordTransf = theTransf.getCopy3d();

	if (!theCoordTransf) {
		opserr << "AIDMBeamColumn::AIDMBeamColumn -- failed to get copy of coordinate transformation\n";
		exit(-1);
	}

	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;
	q0[3] = 0.0;
	q0[4] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;
	p0[3] = 0.0;
	p0[4] = 0.0;

	// set node pointers to NULL
	for (int i = 0; i < 2; i++)
		theNodes[i] = 0;

	//rigid End not Exist
	if (rigidILength <= 0 && rigidJLength <= 0)
		return;

	int nodeID[1];
	int size[1];
	double nodeCrdI[3];
	double nodeCrdJ[3];
	double direction[3];
	size[0] = 3;
	nodeID[0] = Nd1;
	OPS_GetNodeCrd(nodeID, size, nodeCrdI);
	nodeID[0] = Nd2;
	OPS_GetNodeCrd(nodeID, size, nodeCrdJ);
	//Calculate initial Length
	this->initialLength = sqrt(pow(nodeCrdJ[0] - nodeCrdI[0], 2) + pow(nodeCrdJ[1] - nodeCrdI[1], 2) +
		pow(nodeCrdJ[2] - nodeCrdI[2], 2));
	//is End beyond
	if (rigidILength + rigidJLength >= this->initialLength * 0.9)
	{
		opserr << "AIDMBeamColumn::AIDMBeamColumn -- failed to create rigid end: length is out of range.\n";
		return;
	}
	//Calculate vector direction
	for(int i = 0; i < 3 ; i++)
		direction[i] = nodeCrdJ[i] - nodeCrdI[i];
	//I end 
	if (rigidILength > 0)
	{
		auto factor = rigidILength / this->initialLength;
		Vector jntOffset(3);
		for (int i = 0; i < 3; i++)
			jntOffset(i) = direction[i] * factor;
		this->theCoordTransf->setnodeOffset(jntOffset, true);
	}
	if (rigidJLength > 0)
	{
		auto factor = -rigidJLength / this->initialLength;
		Vector jntOffset(3);
		for (int i = 0; i < 3; i++)
			jntOffset(i) = direction[i] * factor;
		this->theCoordTransf->setnodeOffset(jntOffset, false);
	}
}

void AIDMBeamColumn::setRigidEnd(const double& rigidLength, bool isI, double* direction)
{
	auto factor = (isI ? rigidLength : -rigidLength) / this->initialLength;
	Vector jntOffset(3);
	for (int i = 0; i < 3; i++)
		jntOffset(i) = direction[i] * factor;
}

int
AIDMBeamColumn::getNumExternalNodes(void) const
{
	return 2;
}

const ID&
AIDMBeamColumn::getExternalNodes(void)
{
	return connectedExternalNodes;
}

Node**
AIDMBeamColumn::getNodePtrs(void)
{
	return theNodes;
}

int
AIDMBeamColumn::getNumDOF(void)
{
	return 12;
}

void
AIDMBeamColumn::setDomain(Domain* theDomain)
{
	if (theDomain == 0) {
		opserr << "AIDMMaterial::setDomain -- Domain is null\n";
		exit(-1);
	}

	theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
	theNodes[1] = theDomain->getNode(connectedExternalNodes(1));


	if (theNodes[0] == 0) {
		opserr << "AIDMMaterial::setDomain  tag: " << this->getTag() << " -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
		exit(-1);
	}

	if (theNodes[1] == 0) {
		opserr << "AIDMMaterial::setDomain  tag: " << this->getTag() << " -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
		exit(-1);
	}

	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	if (dofNd1 != 6) {
		opserr << "AIDMMaterial::setDomain  tag: " << this->getTag() << " -- Node 1: " << connectedExternalNodes(0)
			<< " has incorrect number of DOF\n";
		exit(-1);
	}

	if (dofNd2 != 6) {
		opserr << "AIDMMaterial::setDomain  tag: " << this->getTag() << " -- Node 2: " << connectedExternalNodes(1)
			<< " has incorrect number of DOF\n";
		exit(-1);
	}

	this->DomainComponent::setDomain(theDomain);

	if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
		opserr << "AIDMMaterial::setDomain  tag: " << this->getTag() << " -- Error initializing coordinate transformation\n";
		exit(-1);
	}

	double L = theCoordTransf->getInitialLength();

	if (L == 0.0) {
		opserr << "AIDMMaterial::setDomain  tag: " << this->getTag() << " -- Element has zero length\n";
		exit(-1);
	}
}

int
AIDMBeamColumn::commitState()
{
	//Ger local resist force
	auto& F = this->getLocalResistingForce();
	//Calculate Shear Span
	this->AIDMs[0]->setLammda(F(4) / F(2));
	this->AIDMs[1]->setLammda(F(10) / F(8));
	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "ElasticBeam3d::commitState () - failed in base class";
	}
	//AIDM CommitState
	for (int i = 0; i < numAIDMs; i++)
	{
		retVal += AIDMs[i]->commitState();
	}
	//is Kill
	if (this->AIDMs[0]->iskill()) iAIDMTag = -1;
	if (this->AIDMs[1]->iskill()) jAIDMTag = -1;

	retVal += theCoordTransf->commitState();
	//is load already const
	if(this->isGravityConst)
		return retVal;
	//is not the same load factor
	if (this->CLoadFactor != this->TLoadFactor)
	{
		this->CLoadFactor = this->TLoadFactor;
		return retVal;
	}
	if(this->CLoadFactor == 0)
		return retVal;
	this->AIDMs[0]->checkCapacity(F(4), this->getTag());
	this->AIDMs[1]->checkCapacity(-F(10), this->getTag());
	this->isGravityConst = true;
}

int
AIDMBeamColumn::revertToLastCommit()
{
	int retVal = theCoordTransf->revertToLastCommit();
	//AIDM CommitState
	for (int i = 0; i < numAIDMs; i++)
	{
		retVal += AIDMs[i]->revertToLastCommit();
	}
	return retVal;
}

int
AIDMBeamColumn::revertToStart()
{
	return theCoordTransf->revertToStart();
}

int
AIDMBeamColumn::update(void)
{
	// Update the transformation
	int retVal = theCoordTransf->update();
	// Get basic deformations
	const Vector& v = theCoordTransf->getBasicTrialDisp();
	//initial stiffness
	if (!this->isInitial)
	{
		double L = theCoordTransf->getInitialLength();
		// Update the transformation
		double oneOverL = 1.0 / L;
		double EoverL = E * oneOverL;
		double EIyoverL3 = 3.0 * Iy * EoverL;		// 3EIy/L
		//is elatic 
		if(iAIDMTag == -1) 
			AIDMs[0]->setInitialK(0, true);
		else
		{
			bool isElastic = iAIDMTag == -2;
			if (jAIDMTag == -1)
				AIDMs[0]->setInitialK(EIyoverL3, isElastic);
			else AIDMs[0]->setInitialK(EIyoverL3 * 2, isElastic);
		}
		if (jAIDMTag == -1) AIDMs[1]->setInitialK(0, true);
		else 
		{
			bool isElastic = jAIDMTag == -2;
			if (jAIDMTag == -1)
				AIDMs[1]->setInitialK(EIyoverL3, isElastic);
			else AIDMs[1]->setInitialK(EIyoverL3 * 2, isElastic);
		}
		this->isInitial = true;
	}
	//IS initial 
	retVal += AIDMs[0]->setTrialStrain(v(3));
	retVal += AIDMs[1]->setTrialStrain(-v(4));

	return retVal < 0? -1: 0;
}

const Matrix&
AIDMBeamColumn::getTangentStiff(void)
{
	const Vector& v = theCoordTransf->getBasicTrialDisp();

	double L = theCoordTransf->getInitialLength();

	//设定刚度
	setStiffMatrix(L);
	//设定局部力
	setBasicForce(L, v);

	return theCoordTransf->getGlobalStiffMatrix(kb, q);
}


const Matrix&
AIDMBeamColumn::getInitialStiff(void)
{
	//  const Vector &v = theCoordTransf->getBasicTrialDisp();

	double L = theCoordTransf->getInitialLength();
	//设定刚度
	setStiffMatrix(L);
	return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix&
AIDMBeamColumn::getMass(void)
{
	K.Zero();
	return K;
}

void
AIDMBeamColumn::zeroLoad(void)
{
	Q.Zero();

	q0[0] = 0.0;
	q0[1] = 0.0;
	q0[2] = 0.0;
	q0[3] = 0.0;
	q0[4] = 0.0;

	p0[0] = 0.0;
	p0[1] = 0.0;
	p0[2] = 0.0;
	p0[3] = 0.0;
	p0[4] = 0.0;

	return;
}

int 
AIDMBeamColumn::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();
  TLoadFactor = loadFactor;

  if (type == LOAD_TAG_Beam3dUniformLoad) 
  {
	  double wy = data(0) * loadFactor;  // Transverse
	  double wz = data(1) * loadFactor;  // Transverse
	  double wx = data(2) * loadFactor;  // Axial (+ve from node I to J)

	  double Vy = 0.5 * wy * L;
	  double Mz = Vy * L / 6.0; // wy*L*L/12
	  double Vz = 0.5 * wz * L;
	  double My = Vz * L / 6.0; // wz*L*L/12
	  double P = wx * L;

	  // Reactions in basic system
	  p0[0] -= P;
	  p0[1] -= Vy;
	  p0[2] -= Vy;
	  p0[3] -= Vz;
	  p0[4] -= Vz;

	  // Fixed end forces in basic system
	  q0[0] -= 0.5 * P;

	  if (iAIDMTag != -1 && jAIDMTag != -1)
      {
		  q0[1] -= Mz;
		  q0[2] += Mz;
		  q0[3] += My;
		  q0[4] -= My;
	  }
	  else if (jAIDMTag != -1)
      {
		  q0[2] += wy * L * L / 8;
		  q0[4] -= wz * L * L / 8;
	  }
	  else if (iAIDMTag != -1)
      {
		  q0[1] -= wy * L * L / 8;
		  q0[3] += wz * L * L / 8;
	  }
  }
  else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
	  double wa = data(2) * loadFactor;  // Axial
	  double wy = data(0) * loadFactor;  // Transverse
	  double wz = data(1) * loadFactor;  // Transverse
	  double a = data(3) * L;
	  double b = data(4) * L;
	  double c = 0.5 * (b + a);
	  double cOverL = c / L;

	  double P = wa * (b - a);
	  double Fy = wy * (b - a);
	  double Fz = wz * (b - a);

	  // Reactions in basic system
	  p0[0] -= P;
	  double V1, V2;
	  V1 = Fy * (1.0 - cOverL);
	  V2 = Fy * cOverL;
	  p0[1] -= V1;
	  p0[2] -= V2;
	  V1 = Fz * (1.0 - cOverL);
	  V2 = Fz * cOverL;
	  p0[3] -= V1;
	  p0[4] -= V2;

	  // Fixed end forces in basic system
	  q0[0] -= P * cOverL;
	  double M1, M2;
	  double beta2 = (1 - cOverL) * (1 - cOverL);
	  double alfa2 = (cOverL) * (cOverL);
	  double gamma2 = (b - a) / L;
	  gamma2 *= gamma2;

	  M1 = -wy * (b - a) * (c * beta2 + gamma2 / 12.0 * (L - 3 * (L - c)));
	  M2 = wy * (b - a) * ((L - c) * alfa2 + gamma2 / 12.0 * (L - 3 * c));
	  q0[1] += M1;
	  q0[2] += M2;
	  M1 = -wz * (b - a) * (c * beta2 + gamma2 / 12.0 * (L - 3 * (L - c)));
	  M2 = wz * (b - a) * ((L - c) * alfa2 + gamma2 / 12.0 * (L - 3 * c));
	  q0[3] -= M1;
	  q0[4] -= M2;
  }
  else if (type == LOAD_TAG_Beam3dPointLoad) 
  {
      //General Point Load
      addPointLoad(data(2) * loadFactor, data(0) * loadFactor, data(1) * loadFactor, data(3), L);
  }
  //多态荷载模式
  //data[0]: Pyi
  //data[1]: Pzi
  //data[2]: Ni
  //data[3]: Pyj
  //data[4]: Pzj
  //data[5]: Nj
  //data[6]: aOverL
  //data[7]: bOverL
  else if (type == LOAD_TAG_Beam3dGenranlPartialLoad)
  {
	  double Pyi = data(0) * loadFactor;
	  double Pzi = data(1) * loadFactor;
	  double Ni = data(2) * loadFactor;

	  double Pyj = data(3) * loadFactor;
	  double Pzj = data(4) * loadFactor;
	  double Nj = data(5) * loadFactor;

	  double aOverL = data(6);
	  double bOverL = data(7);

	  addGeneralPartialLoad(Ni, Nj, Pyi, Pyj, Pzi, Pzj, aOverL, bOverL, L);
  }

  //三角形荷载
  //data[0]: Py
  //data[1]: Pz
  //data[2]: aOverL
  else if (type == LOAD_TAG_Beam3dTriangularLoad)
  {
	  double Py = data(0) * loadFactor;
	  double Pz = data(1) * loadFactor;

	  double aOverL = data(2);
	  //左侧三角形
	  addGeneralPartialLoad(0, 0, 0, Py, 0, Pz, 0, aOverL, L);
	  //右侧三角形
	  addGeneralPartialLoad(0, 0, Py, 0, Pz, 0, aOverL, 1, L);
  }

  //梯形荷载
  //data[0]: Py
  //data[1]: Pz
  //data[2]: aOverL
  //data[3]: bOverL
  else if (type == LOAD_TAG_Beam3dTrapezoidLoad)
  {
	  double Py = data(0) * loadFactor;
	  double Pz = data(1) * loadFactor;

	  double aOverL = data(2);
	  double bOverL = data(3);
	  //左侧三角形
	  addGeneralPartialLoad(0, 0, 0, Py, 0, Pz, 0, aOverL, L);
	  //中央三角形
	  addGeneralPartialLoad(0, 0, Py, Py, Pz, Pz, aOverL, bOverL, L);
	  //右侧三角形
	  addGeneralPartialLoad(0, 0, Py, 0, Pz, 0, bOverL, 1, L);
  }
  else 
  {
    opserr << "ElasticBeam3d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
AIDMBeamColumn::addInertiaLoadToUnbalance(const Vector& accel)
{
	return 0;
}

const Vector&
AIDMBeamColumn::getResistingForceIncInertia()
{
	P = this->getResistingForce();

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

	return P;
}

const Vector&
AIDMBeamColumn::getResistingForce()
{
	// Get basic deformations
	const Vector& v = theCoordTransf->getBasicTrialDisp();
	// Length of the element
	double L = theCoordTransf->getInitialLength();
	//BasicForce
	this->setBasicForce(L, v);
	// Vector for reactions in basic system
	Vector p0Vec(p0, 5);
	return theCoordTransf->getGlobalResistingForce(q, p0Vec);
}


const Vector& AIDMBeamColumn::getLocalResistingForce(void)
{
	double N, V, M1, M2, T;
	double L = theCoordTransf->getInitialLength();
	double oneOverL = 1.0 / L;

	// Axial
	N = q(0);
	P(6) = N;
	P(0) = -N + p0[0];

	// Torsion
	T = q(5);
	P(9) = T;
	P(3) = -T;

	// Moments about z and shears along y
	M1 = q(1);
	M2 = q(2);
	P(5) = M1;
	P(11) = M2;
	V = (M1 + M2) * oneOverL;
	P(1) = V + p0[1];
	P(7) = -V + p0[2];

	// Moments about y and shears along z
	M1 = q(3);
	M2 = q(4);
	P(4) = M1;
	P(10) = M2;
	V = (M1 + M2) * oneOverL;
	P(2) = -V + p0[3];
	P(8) = V + p0[4];

	return P;
}

int
AIDMBeamColumn::sendSelf(int cTag, Channel& theChannel)
{
	return -1;
}

int
AIDMBeamColumn::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	return -1;
}

void
AIDMBeamColumn::Print(OPS_Stream& s, int flag)
{
	this->getResistingForce();

	if (flag == -1) {
		int eleTag = this->getTag();
		s << "EL_BEAM\t" << eleTag << "\t";
		s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
		s << "\t0\t0.0000000\n";
	}

	else if (flag < -1) {
		int counter = (flag + 1) * -1;
		int eleTag = this->getTag();
		const Vector& force = this->getResistingForce();

		double P, MZ1, MZ2, VY, MY1, MY2, VZ, T;
		double L = theCoordTransf->getInitialLength();
		double oneOverL = 1.0 / L;

		P = q(0);
		MZ1 = q(1);
		MZ2 = q(2);
		VY = (MZ1 + MZ2) * oneOverL;
		MY1 = q(3);
		MY2 = q(4);
		VZ = (MY1 + MY2) * oneOverL;
		T = q(5);

		s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
		s << "\t" << -P + p0[0] << "\t" << VY + p0[1] << "\t" << -VZ + p0[3] << endln;
		s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
		s << "\t" << P << ' ' << -VY + p0[2] << ' ' << VZ + p0[4] << endln;
		s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
		s << "\t" << -T << "\t" << MY1 << "\t" << MZ1 << endln;
		s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
		s << "\t" << T << ' ' << MY2 << ' ' << MZ2 << endln;
	}

	else if (flag == 2) {
		this->getResistingForce(); // in case linear algo

		static Vector xAxis(3);
		static Vector yAxis(3);
		static Vector zAxis(3);

		theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis);

		s << "#ElasticBeamColumn3D\n";
		s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
		s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
		s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

		const Vector& node1Crd = theNodes[0]->getCrds();
		const Vector& node2Crd = theNodes[1]->getCrds();
		const Vector& node1Disp = theNodes[0]->getDisp();
		const Vector& node2Disp = theNodes[1]->getDisp();

		s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
			<< " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
			<< " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << endln;

		s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
			<< " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
			<< " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << endln;

		double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
		double L = theCoordTransf->getInitialLength();
		double oneOverL = 1.0 / L;

		N = q(0);
		Mz1 = q(1);
		Mz2 = q(2);
		Vy = (Mz1 + Mz2) * oneOverL;
		My1 = q(3);
		My2 = q(4);
		Vz = -(My1 + My2) * oneOverL;
		T = q(5);

		s << "#END_FORCES " << -N + p0[0] << ' ' << Vy + p0[1] << ' ' << Vz + p0[3] << ' '
			<< -T << ' ' << My1 << ' ' << Mz1 << endln;
		s << "#END_FORCES " << N << ' ' << -Vy + p0[2] << ' ' << -Vz + p0[4] << ' '
			<< T << ' ' << My2 << ' ' << Mz2 << endln;
	}

	if (flag == OPS_PRINT_CURRENTSTATE) {

		this->getResistingForce(); // in case linear algo

		s << "\nElasticBeam3d: " << this->getTag() << endln;
		s << "\tConnected Nodes: " << connectedExternalNodes;
		s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;

		double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
		double L = theCoordTransf->getInitialLength();
		double oneOverL = 1.0 / L;

		N = q(0);
		Mz1 = q(1);
		Mz2 = q(2);
		Vy = (Mz1 + Mz2) * oneOverL;
		My1 = q(3);
		My2 = q(4);
		Vz = -(My1 + My2) * oneOverL;
		T = q(5);

		s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
			<< -N + p0[0] << ' ' << Mz1 << ' ' << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << endln;
		s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
			<< N << ' ' << Mz2 << ' ' << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << endln;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"ElasticBeam3d\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		s << "\"E\": " << E << ", ";
		s << "\"G\": " << G << ", ";
		s << "\"A\": " << A << ", ";
		s << "\"Jx\": " << Jx << ", ";
		s << "\"Iy\": " << Iy << ", ";
		s << "\"Iz\": " << Iz << ", ";
		s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
	}
}

int
AIDMBeamColumn::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	static Vector v1(3);
	static Vector v2(3);

	theNodes[0]->getDisplayCrds(v1, fact);
	theNodes[1]->getDisplayCrds(v2, fact);
	float d1 = 0.0;
	float d2 = 0.0;
	int res = 0;

	if (displayMode > 0) {

		res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), 0);

	}
	else if (displayMode < 0) {

		theNodes[0]->getDisplayCrds(v1, 0.);
		theNodes[1]->getDisplayCrds(v2, 0.);

		// add eigenvector values
		int mode = displayMode * -1;
		const Matrix& eigen1 = theNodes[0]->getEigenvectors();
		const Matrix& eigen2 = theNodes[1]->getEigenvectors();
		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 3; i++) {
				v1(i) += eigen1(i, mode - 1) * fact;
				v2(i) += eigen2(i, mode - 1) * fact;
			}
		}
		return theViewer.drawLine(v1, v2, 0.0, 0.0, this->getTag(), 0);
	}

	if (numMode > 0) {
		// calculate q for potential need below
		this->getResistingForce();
	}

	for (int i = 0; i < numMode; i++) {

		const char* theMode = modes[i];
		if (strcmp(theMode, "axialForce") == 0) {
			d1 = q(0);
			d2 = q(0);;

			res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

		}
		else if (strcmp(theMode, "endMoments") == 0) {
			d1 = q(1);
			d2 = q(2);
			static Vector delta(3); delta = v2 - v1; delta /= 10;
			res += theViewer.drawPoint(v1 + delta, d1, this->getTag(), i);
			res += theViewer.drawPoint(v2 - delta, d2, this->getTag(), i);

		}
	}

	return res;
}

Response*
AIDMBeamColumn::setResponse(const char** argv, int argc, OPS_Stream& output)
{

	Response* theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "ElasticBeam3d");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// global forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {


		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Pz_1");
		output.tag("ResponseType", "Mx_1");
		output.tag("ResponseType", "My_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Pz_2");
		output.tag("ResponseType", "Mx_2");
		output.tag("ResponseType", "My_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 2, P);

		// local forces
	}
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

		output.tag("ResponseType", "N_1");
		output.tag("ResponseType", "Vy_1");
		output.tag("ResponseType", "Vz_1");
		output.tag("ResponseType", "T_1");
		output.tag("ResponseType", "My_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "N_2");
		output.tag("ResponseType", "Vy_2");
		output.tag("ResponseType", "Vz_2");
		output.tag("ResponseType", "T_2");
		output.tag("ResponseType", "My_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 3, P);

		// basic forces
	}
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {

		output.tag("ResponseType", "N");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Mz_2");
		output.tag("ResponseType", "My_1");
		output.tag("ResponseType", "My_2");
		output.tag("ResponseType", "T");

		theResponse = new ElementResponse(this, 4, Vector(6));

	}
	else if (strcmp(argv[0], "deformations") == 0 || strcmp(argv[0], "deformation") == 0 ||
		strcmp(argv[0], "basicDeformations") == 0) {

		output.tag("ResponseType", "eps");
		output.tag("ResponseType", "theta11");
		output.tag("ResponseType", "theta12");
		output.tag("ResponseType", "theta21");
		output.tag("ResponseType", "theta22");
		output.tag("ResponseType", "phi");
		theResponse = new ElementResponse(this, 5, Vector(6));
	}
	else if (strcmp(argv[0], "strainStress") == 0 || strcmp(argv[0], "strainstress") == 0 ||
		strcmp(argv[0], "StrainStress") == 0) 
	{
		for (int i = 0; i < this->numAIDMs; i++)
		{
			std::stringstream ss_strain;
			ss_strain << "strain" << i;
			output.tag("ResponseType", ss_strain.str().c_str());
		}
		for (int i = 0; i < this->numAIDMs; i++)
		{
			std::stringstream ss_stress;
			ss_stress << "stress" << i;
			output.tag("ResponseType", ss_stress.str().c_str());
		}
		for (int i = 0; i < this->numAIDMs; i++)
		{
			std::stringstream ss_stress;
			ss_stress << "lammda" << i;
			output.tag("ResponseType", ss_stress.str().c_str());
		}
		theResponse = new ElementResponse(this, 6, Vector(this->numAIDMs * 3));
	}
	output.endTag(); // ElementOutput

	return theResponse;
}

int
AIDMBeamColumn::getResponse(int responseID, Information& eleInfo)
{
	static Vector Res(12);
	Res = this->getResistingForce();
	static Vector strainStress(this->numAIDMs * 3);
	for (int i = 0; i < this->numAIDMs; i++)
	{
		strainStress(i) = this->AIDMs[i]->getStrain();
	}
	for (int i = 0; i < this->numAIDMs; i++)
	{
		strainStress(i + this->numAIDMs) = this->AIDMs[i]->getStress();
	}
	for (int i = 0; i < this->numAIDMs; i++)
	{
		strainStress(i + this->numAIDMs * 2) = this->AIDMs[i]->getLammda();
	}
	switch (responseID) {
	case 1: // stiffness
		return eleInfo.setMatrix(this->getTangentStiff());

	case 2: // global forces
		return eleInfo.setVector(Res);

	case 3: // local forces
		return eleInfo.setVector(this->getLocalResistingForce());

	case 4: // basic forces
		return eleInfo.setVector(q);

	case 5:
		return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());

	case 6:
		return eleInfo.setVector(strainStress);

	default:
		return -1;
	}
}


int
AIDMBeamColumn::setParameter(const char** argv, int argc, Parameter& param)
{
	return -1;
}

int
AIDMBeamColumn::updateParameter(int parameterID, Information& info)
{
	return -1;
}

