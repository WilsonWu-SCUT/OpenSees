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
    //设定剪跨比
     void setLammda(const double& Lammda);
     //设定承载力
     inline void setCapcacity(const double& m_pos, const double& m_neg)
     {
         this->stressSA_pos = m_pos;
         this->stressSA_neg = m_neg;
     };
    //判断承载力是否越界
    bool checkCapacity(const double& Moment, const int& eleTag);
    //设定初始刚度
    void setInitialK(const double K);
    //设定纵筋配筋比
    inline void setlammdaT_pos(const double& lammdaT)
    {
        this->lammdaT_pos = lammdaT;
    }
    //是否为有效的AIDM对象
    bool isAvailabelAIDM() const;

  private:
      //从骨架上计算应力
      double getStressOnBackbone(const double& drift);
      //从骨架上计算切线刚度
      void setTangentOnBackbone(const double& strain, bool loading_direct_pos);
      //计算卸载刚度
      int setUnloadingTangent(const double& strain);
      //计算重加载刚度
      void setReloadingTangent(const double& strain, bool loading_direct_pos);

  private:
      //获得处理为ANN输出参数的基本变量
      std::vector<double>& getRegularizedValueVector(const AIDMParamEnum& type);
      //参数转为ANN输入参数
      float getRegularizedValue(const double& value, const AIDMParamEnum& type);
      //ANN输出参数转为普通参数
      float getNormalValue(const double& value, const AIDMParamEnum& type);
      //输入层的构件特征参数
      std::vector<float> getComponentParamsVec(bool is_pos);
      //更新骨架参数
      void updateSkeletonParams();
      //更新滞回参数
      void updateHystereticParams(bool is_pos);

  private:
      //骨架曲线指向点应变系数
      double backbone_ortStrainFactor = 1.2;
      //计算初始刚度采用的峰值应变系数
      double backbone_inidStrainFactor = 0.1;
      //初始剪跨比
      double initialLammda = 4;
      //杀死单元的应力系数
      double killStressFactor = 0.2;
     

  private: /*力学参数*/
      //剪跨比（需要revert）
      double lammda;
      double Clammda;
      //Caping Capacity based on section Analysis positve value
      //峰值承载力
      double stressSA_pos;
      double stressSA_neg;
      //形状调整系数
      int mnFactor;

  private: /*基本参数*/
      //纵筋配筋特征值
      double lammdaS;
      //面积配箍特征值
      double lammdaSV;

      //受拉纵筋与受压纵筋比
      //通过PMMSection计算
      double lammdaT_pos;

    #pragma region AIDM参数
  private:   /*是否更新本构参数*/
      bool needUpdateHANN_pos;
      bool needUpdateHANN_neg;
      bool needUpdateBANN;

  private: /*骨架参数*/
      //shape paramters
      double m_pos;
      double m_neg;
      double n_pos;
      double n_neg;
      //caping drift positive value
      double strainC_pos;
      double strainC_neg;
      //承载力系数
      double stressFactor_pos;
      double stressFactor_neg;

  private: /*滞回参数*/
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
      //Stiffness 刚度
      double K;
      double CK;

  private:
      //MaxDeformation 历史最大变形
      double CstrainMax;
      double CstrainMin;
      //历史最大变形对应的承载力系数
      double CstressMaxFactor;
      double CstressMinFactor;

  private:
    //迭代步及分析步的应力应变
    double TStrain;
    double CStress;
    double CStrain;
    double TStress;
    /*加载方向*/
    bool CLoadingDirectPos;
    bool TLoadingDirectPos;

  private:
    //刚度是否转为弹性
    bool isConstant;
};


#endif

