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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterial.h,v $
                                                                        
                                                                        
#ifndef AIDMMaterial_h
#define AIDMMaterial_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMaterial. ElasticMaterial provides the abstraction
// of an viscoelastic uniaxial material,
// i.e. stress = E*strain + eta*strainrate


#include <UniaxialMaterial.h>
#include "KerasModelExport.h"

class AIDMMaterial : public UniaxialMaterial
{
  private:
      enum class AIDMParamEnum
      {
          Lammda = 0,
          LammdaS = 1,
          LammdaSV = 2,
          LammdaT = 3,

          SecantK = 10,

          StressFactor = 101,
          StrainC = 102,
          m = 103,
          n = 104,

          Afa = 201,
          Beta = 202,
          Gamma = 203,
          Eta = 204,
      };

  public:
      AIDMMaterial(int tag, double height, double width, double lammdaS, double lammdaSV, double lammdaT_pos, double Msa_pos, double Msa_neg);
      AIDMMaterial(int tag, double lammda, double lammdaS, double lammdaSV, double lammdaT_pos, double Msa_pos, double Msa_neg);
      AIDMMaterial();
    ~AIDMMaterial();

 public:
     static std::shared_ptr<KerasModel> keras_SANN_sp;
     static std::shared_ptr<KerasModel> keras_HANN_sp;
     static std::vector<double> lammda_vec;
     static std::vector<double> lammdaS_vec;
     static std::vector<double> lammdaSV_vec;
     static std::vector<double> lammdaT_vec;
     static std::vector<double> strainC_vec;
     static std::vector<double> stressFactor_vec;
     static std::vector<double> m_vec;
     static std::vector<double> n_vec;
     static std::vector<double> secantK_vec;
     static std::vector<double> afa_vec;
     static std::vector<double> beta_vec;
     static std::vector<double> gamma_vec;
     static std::vector<double> eta_vec;

public:
    const char *getClassType(void) const {return "AIDMMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return TStrain;};
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);
    
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    public:
    void setLammda(const double& shearSpan);

  protected:

  private:
      double getStressOnBackbone(const double& drift);
      void setTangentOnBackbone(const double& strain, bool loading_direct_pos);
      int setUnloadingTangent(const double& strain);
      void setReloadingTangent(const double& strain, bool loading_direct_pos);

  private:
      std::vector<double>& getRegularizedValueVector(const AIDMParamEnum& type);
      float getRegularizedValue(const double& value, const AIDMParamEnum& type);
      float getNormalValue(const double& value, const AIDMParamEnum& type);
      std::vector<float> getComponentParamsVec(bool is_pos);
      void updateSkeletonParams();
      void updateHystereticParams(bool is_pos);

  private:
      double backbone_ortStrainFactor = 1.2;
      double backbone_inidStrainFactor = 0.05;
      double initialLammda = 2;

  private:
      double lammda;
      double lammdaS;
      double lammdaSV;
      double lammdaT_pos;

  private:
      //shape paramters
      double m_pos; 
      double m_neg;
      double n_pos;
      double n_neg;
      //caping drift positive value
      double strainC_pos;
      double strainC_neg;
      double stressFactor_pos;
      double stressFactor_neg;
      //Caping Capacity based on section Analysis positve value
      double stressSA_pos;
      double stressSA_neg;
      //Hysteretic parameters
      double afa_pos;
      double afa_neg;
      double beta_pos;
      double beta_neg;
      double gamma_pos;
      double gamma_neg;
      double eta_pos;
      double eta_neg;
      //MaxDeformation
      double CstrainMax;
      double CstrainMin;
      double CstressMaxCor;
      double CStressMinCor;
      //Stiffness
      double K;
      //Need to update hysteretic
      bool needUpdateHANN_pos;
      bool needUpdateHANN_neg;
      bool needUpdateBANN;
      //Section Height
      double sectionHeight;

  private:
    double TStrain;
    double CStress;
    double CStrain;
    double TStress;
    double CK;
    bool loading_direct_pos;
    double Clammda;

};


#endif

