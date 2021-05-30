
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
class PMMSection;

class AIDMBeamColumn : public Element
{
  public:
    AIDMBeamColumn();
    AIDMBeamColumn(int tag, int Nd1, int Nd2, CrdTransf& theTransf, 
        const std::vector<int> tag_vec, const double& rigidILength, const double& rigidJLength);
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
    /*AIDMBeamColumn初始化*/
    void initialAIDMBeamColumn(int Nd1, int Nd2, CrdTransf& theTransf, const double& rigidILength, const double& rigidJLength);
    /*设定刚度矩阵*/
	void setStiffMatrix(const double& length, bool is_initial);
	void setBasicForce(const double& length, const Vector& v);
    /*添加节点荷载*/
	int addPointLoad(const double& N, const double& Py, const double& Pz,
		const double& aOverL, const double& length);
    /*添加任意荷载*/
	void addGeneralPartialLoad(const double& Ni, const double& Nj, const double& Pyi,
		const double& Pyj, const double& Pzi, const double& Pzj,
		const double& aOverL, const double& bOverL, const double& length);
    /*设定剪跨比*/
    void setLammda(const Vector& force_vec);

private:
    //是否硬化初始刚度
    bool ensureIniK = true;

  private:
      /*通过恒在判断是否转化为弹性*/
      bool isGravityConst;
      double TLoadFactor;
      double CLoadFactor;
      /*是否初始化*/
      bool isInitial;

private:
    PMMSection* sectionI_ptr;
    PMMSection* sectionJ_ptr;
    int sectionI_tag_;
    int sectionJ_tag_;

  private:
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
