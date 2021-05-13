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
      AIDMMaterial(int tag, double height, double width, double lammdaS, double lammdaSV, double lammdaT_pos, double Msa_pos, double Msa_neg, bool ensureIniK = true);
      AIDMMaterial(int tag, double lammda, double lammdaS, double lammdaSV, double lammdaT_pos, double Msa_pos, double Msa_neg, bool ensureIniK = true);
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
     static int predict_num;

public:
    const char *getClassType(void) const {return "AIDMMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return TStrain;};
    bool iskill(void) { return this->isKill; };
    double getLammda(void) { return this->Clammda; }
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
    //�趨�����
    void setLammda(const double& shearSpan);
    //�жϳ������Ƿ�Խ��
    bool checkCapacity(const double& Moment, const int& eleTag);
    //�趨��ʼ�ն�
    void setInitialK(double K, bool isToElastic);

  protected:

  private:
      //�ӹǼ��ϼ���Ӧ��
      double getStressOnBackbone(const double& drift);
      //�ӹǼ��ϼ������߸ն�
      void setTangentOnBackbone(const double& strain, bool loading_direct_pos);
      //����ж�ظն�
      int setUnloadingTangent(const double& strain);
      //�����ؼ��ظն�
      void setReloadingTangent(const double& strain, bool loading_direct_pos);

  private:
      //��ô���ΪANN��������Ļ�������
      std::vector<double>& getRegularizedValueVector(const AIDMParamEnum& type);
      //����תΪANN�������
      float getRegularizedValue(const double& value, const AIDMParamEnum& type);
      //ANN�������תΪ��ͨ����
      float getNormalValue(const double& value, const AIDMParamEnum& type);
      //�����Ĺ�����������
      std::vector<float> getComponentParamsVec(bool is_pos);
      //���¹Ǽܲ���
      void updateSkeletonParams();
      //�����ͻز���
      void updateHystereticParams(bool is_pos);
      //�Ƿ�Ϊ��Ч��AIDM����
      bool isAvailabelAIDM() const;

  private:
      //�Ǽ�����ָ���Ӧ��ϵ��
      double backbone_ortStrainFactor = 1.2;
      //�����ʼ�նȲ��õķ�ֵӦ��ϵ��
      double backbone_inidStrainFactor = 0.1;
      //��ʼ�����
      double initialLammda = 4;
      //ɱ����Ԫ��Ӧ��ϵ��
      double killStressFactor = 0.2;
      int mnFactor = 1;

  private:
      //�����
      double lammda;
      //�ݽ��������ֵ
      double lammdaS;
      //����乿����ֵ
      double lammdaSV;

      //�����ݽ�����ѹ�ݽ��
      //ͨ��PMMSection����
      double lammdaT_pos;

      //��ʼ�ն�
      double initialK;

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
      double CstressMaxFactor;
      double CstressMinFactor;
      //Stiffness
      double K;
      //Need to update hysteretic
      bool needUpdateHANN_pos;
      bool needUpdateHANN_neg;
      bool needUpdateBANN;
      //Section Height
      //ͨ��PMM�����ȡ
      double sectionHeight;

  private:
    double TStrain;
    double CStress;
    double CStrain;
    double TStress;
    double CK;
    double Clammda;
    bool CLoadingDirectPos;
    bool TLoadingDirectPos;
    bool ensureIniK;
    //��Ԫ�Ƿ�ɱ��
    bool isKill;
    //��Ԫ�Ƿ���
    bool isToElastic;
};


#endif

