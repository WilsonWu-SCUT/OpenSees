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
     static int predict_num;

public:
    const char *getClassType(void) const {return "AIDMMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return TStrain;};
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
     void setLammda(const double& Lammda);
     //�趨������
     inline void setCapcacity(const double& m_pos, const double& m_neg)
     {
         this->stressSA_pos = m_pos;
         this->stressSA_neg = m_neg;
     };
    //�жϳ������Ƿ�Խ��
    bool checkCapacity(const double& Moment, const int& eleTag);
    //�趨��ʼ�ն�
    void setInitialK(const double K);
    //�趨�ݽ�����
    inline void setlammdaT_pos(const double& lammdaT)
    {
        this->lammdaT_pos = lammdaT;
    }
    //�Ƿ�Ϊ��Ч��AIDM����
    bool isAvailabelAIDM() const;

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

  private:
      //�Ǽ�����ָ���Ӧ��ϵ��
      double backbone_ortStrainFactor = 1.2;
      //�����ʼ�նȲ��õķ�ֵӦ��ϵ��
      double backbone_inidStrainFactor = 0.1;
      //��ʼ�����
      double initialLammda = 4;
      //ɱ����Ԫ��Ӧ��ϵ��
      double killStressFactor = 0.2;
     

  private: /*��ѧ����*/
      //����ȣ���Ҫrevert��
      double lammda;
      double Clammda;
      //Caping Capacity based on section Analysis positve value
      //��ֵ������
      double stressSA_pos;
      double stressSA_neg;
      //��״����ϵ��
      int mnFactor;

  private: /*��������*/
      //�ݽ��������ֵ
      double lammdaS;
      //����乿����ֵ
      double lammdaSV;

      //�����ݽ�����ѹ�ݽ��
      //ͨ��PMMSection����
      double lammdaT_pos;

    #pragma region AIDM����
  private:   /*�Ƿ���±�������*/
      bool needUpdateHANN_pos;
      bool needUpdateHANN_neg;
      bool needUpdateBANN;

  private: /*�Ǽܲ���*/
      //shape paramters
      double m_pos;
      double m_neg;
      double n_pos;
      double n_neg;
      //caping drift positive value
      double strainC_pos;
      double strainC_neg;
      //������ϵ��
      double stressFactor_pos;
      double stressFactor_neg;

  private: /*�ͻز���*/
      double afa_pos;
      double afa_neg;
      double beta_pos;
      double beta_neg;
      double gamma_pos;
      double gamma_neg;
      double eta_pos;
      double eta_neg;

    #pragma endregion

   private:
      //Stiffness �ն�
      double K;
      double CK;

  private:
      //MaxDeformation ��ʷ������
      double CstrainMax;
      double CstrainMin;
      //��ʷ�����ζ�Ӧ�ĳ�����ϵ��
      double CstressMaxFactor;
      double CstressMinFactor;

  private:
    //����������������Ӧ��Ӧ��
    double TStrain;
    double CStress;
    double CStrain;
    double TStress;
    /*���ط���*/
    bool CLoadingDirectPos;
    bool TLoadingDirectPos;

  private:
    //�ն��Ƿ�תΪ����
    bool isConstant;
};


#endif

