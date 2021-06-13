                                                                                                                                
#ifndef AIDMMaterial_h
#define AIDMMaterial_h

#include "KerasModelExport.h"
#include "DALRMaterial.h"

class AIDMMaterial
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
    AIDMMaterial(const double& lammdaSV, const double& lammdaS, const double& lammdaT_pos);
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
    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain(void) {return TStrain;};
    double getLammda(void) { return this->lammda; }
    double getAxialRatio(void) { return this->axialRatio; }
    int getLoadingType(void) { return this->TLoadingTag; }
    double getInitialStrain(bool isPos) 
    {
        return this->getCstrain(isPos) * this->backbone_inidStrainFactor;
    }
    double getCstrain(bool isPos)
    {
        return (isPos ? this->strainC_pos : this->strainC_neg);
    }
    double getCALR(void) {
        return this->dalr_sp->GetCALR();
    }

    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    AIDMMaterial* getCopy(void); 

  public:
      //设定斜率
      void setARK(const double& ark, const double& factor)
      {
          this->dalr_sp->setARK(ark);
          this->dalr_sp->setFactor(factor);
      }
    //设定剪跨比
     void setLammda(const double& Lammda);
     //设定承载力
     void setCapcacity(const double& m_pos, const double& m_neg, const double& axialRatio);
    //判断承载力是否越界
    bool checkCapacity(const double& Moment);
    //设定初始刚度
    void setInitialK(const double K, bool ensureIniK);
    //是否为有效的AIDM对象
    bool isAvailabelAIDM() const;
    //单元是否已被杀死
    bool isKill() const;
    void Kill();
    //重设刚度
    void ResetTrialStrain(const double& strain, const int& loadingTag);
    //计算骨架承载力退化系数
    double GetSoftenFactor(bool isPos);

  private:
      
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
      //更新滞回参数
      void updateHystereticParams(bool is_pos);
      //是否需要更新剪跨或轴压系数
      bool isValidUpdate();


      public:
          //更新骨架参数
          void updateSkeletonParams();
          //从骨架上计算应力
          double getStressOnBackbone(const double& drift, const int& mnFactor);
          //获得附加约束轴力
          double GetALR();

  private:
      //骨架曲线指向点应变系数
      double backbone_ortStrainFactor = 1.2;
      //弹性段（超过该段才开始更新滞回参数）
      double backbone_inidStrainFactor = 0.1;
      //初始刚度
      double backbone_iniKFactor = 0.05;
      //初始剪跨比
      double initialLammda = 4;
      //杀死单元的应力系数
      double killStressFactor = 0.1;
      //更新剪跨比 轴压系数的变形边界
      double updateDeformationFactor = 0.5;
      //更新剪跨比 轴压系数的承载力边界
      double updateForceFactor = 0.5;

  private: /*力学参数*/
      //剪跨比（需要revert）
      double lammda;
      //轴压比
      double axialRatio;
      //Caping Capacity based on section Analysis positve value
      //峰值承载力
      double stressSA_pos;
      double stressSA_neg;
      //形状调整系数
      int mnFactor_;

  private: /*基本参数*/
      //面积配箍特征值
      double lammdaSV;
      //受拉纵筋与受压纵筋比（通过PMMSection计算）
      double lammdaT_pos; 
      //纵筋配筋特征值（通过PMMSection计算）
      double lammdaS;

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
    //加载模式 1: 骨架 2：重加载 3：卸载
    int TLoadingTag;

    private:
        //轴向约束恢复力模型
        std::shared_ptr<DALRMaterial> dalr_sp;
};


#endif

