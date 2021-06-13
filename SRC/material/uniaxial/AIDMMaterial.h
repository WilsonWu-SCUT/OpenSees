                                                                                                                                
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
      //�趨б��
      void setARK(const double& ark, const double& factor)
      {
          this->dalr_sp->setARK(ark);
          this->dalr_sp->setFactor(factor);
      }
    //�趨�����
     void setLammda(const double& Lammda);
     //�趨������
     void setCapcacity(const double& m_pos, const double& m_neg, const double& axialRatio);
    //�жϳ������Ƿ�Խ��
    bool checkCapacity(const double& Moment);
    //�趨��ʼ�ն�
    void setInitialK(const double K, bool ensureIniK);
    //�Ƿ�Ϊ��Ч��AIDM����
    bool isAvailabelAIDM() const;
    //��Ԫ�Ƿ��ѱ�ɱ��
    bool isKill() const;
    void Kill();
    //����ն�
    void ResetTrialStrain(const double& strain, const int& loadingTag);
    //����Ǽܳ������˻�ϵ��
    double GetSoftenFactor(bool isPos);

  private:
      
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
      //�����ͻز���
      void updateHystereticParams(bool is_pos);
      //�Ƿ���Ҫ���¼������ѹϵ��
      bool isValidUpdate();


      public:
          //���¹Ǽܲ���
          void updateSkeletonParams();
          //�ӹǼ��ϼ���Ӧ��
          double getStressOnBackbone(const double& drift, const int& mnFactor);
          //��ø���Լ������
          double GetALR();

  private:
      //�Ǽ�����ָ���Ӧ��ϵ��
      double backbone_ortStrainFactor = 1.2;
      //���ԶΣ������öβſ�ʼ�����ͻز�����
      double backbone_inidStrainFactor = 0.1;
      //��ʼ�ն�
      double backbone_iniKFactor = 0.05;
      //��ʼ�����
      double initialLammda = 4;
      //ɱ����Ԫ��Ӧ��ϵ��
      double killStressFactor = 0.1;
      //���¼���� ��ѹϵ���ı��α߽�
      double updateDeformationFactor = 0.5;
      //���¼���� ��ѹϵ���ĳ������߽�
      double updateForceFactor = 0.5;

  private: /*��ѧ����*/
      //����ȣ���Ҫrevert��
      double lammda;
      //��ѹ��
      double axialRatio;
      //Caping Capacity based on section Analysis positve value
      //��ֵ������
      double stressSA_pos;
      double stressSA_neg;
      //��״����ϵ��
      int mnFactor_;

  private: /*��������*/
      //����乿����ֵ
      double lammdaSV;
      //�����ݽ�����ѹ�ݽ�ȣ�ͨ��PMMSection���㣩
      double lammdaT_pos; 
      //�ݽ��������ֵ��ͨ��PMMSection���㣩
      double lammdaS;

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
    //����ģʽ 1: �Ǽ� 2���ؼ��� 3��ж��
    int TLoadingTag;

    private:
        //����Լ���ָ���ģ��
        std::shared_ptr<DALRMaterial> dalr_sp;
};


#endif

