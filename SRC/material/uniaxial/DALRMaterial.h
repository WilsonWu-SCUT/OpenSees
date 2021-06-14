                                                                                                                                
#ifndef DALRMaterial_h
#define DALRMaterial_h

class DALRMaterial
{

public:
    DALRMaterial();
    ~DALRMaterial();

public:
    //获得轴压系数
    double GetALROnBackbone(const double& TStrain, const double& dc, 
        const double& lammda, const double& CStressFactor);
    double GetUnloadALR(const double& TStrain);
    double GetReloadALR(const double& TStrain, const double& dc, const double& strainMax,
        const double& lammda, const double& CStressFactor);
    double GetCALR(void) 
    {
        return -this->CALR * this->ReductFactor;
    }
    //设定初始刚度
    void setARK(const int& ark) { this->ARK = ark; }
    //设定折减系数
    void setFactor(const double& factor) { this->ReductFactor = factor; }
    int getARK(void) { return this->ARK; }
    double getFactor(void) { return this->ReductFactor; }

private:
    //初始刚度
    double GetInitialK(const double& lammda);
    //峰值轴压系数
    double GetALRMax(const double& lammda, const double& dc);

private:
    double Cstrain;
    //轴压系数
    double CALR;
    //初始刚度
    int ARK;
    //轴压折减系数 考虑翼缘
    double ReductFactor;
};
#endif

