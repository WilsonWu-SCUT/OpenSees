                                                                                                                                
#ifndef DALRMaterial_h
#define DALRMaterial_h

class DALRMaterial
{

public:
    DALRMaterial();
    ~DALRMaterial();

public:
    //�����ѹϵ��
    double GetALROnBackbone(const double& TStrain, const double& dc, 
        const double& lammda, const double& CStressFactor);
    double GetUnloadALR(const double& TStrain);
    double GetReloadALR(const double& TStrain, const double& dc, const double& strainMax,
        const double& lammda, const double& CStressFactor);
    double GetCALR(void) 
    {
        return -this->CALR;
    }
    //�趨��ʼ�ն�
    void setARK(const double& ark) { this->ARK = ark; }
    double getARK(void) { return this->ARK; }

private:
    //��ʼ�ն�
    double GetInitialK(const double& lammda);
    //��ֵ��ѹϵ��
    double GetALRMax(const double& lammda, const double& dc);

private:
    double Cstrain;
    //��ѹϵ��
    double CALR;
    //��ʼ�ն�
    double ARK;
};
#endif

