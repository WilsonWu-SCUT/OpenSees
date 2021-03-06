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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterial.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) ElasticMaterial.C, revA"

#include <AIDMMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_AIDMMaterial(void)
{
    opserr << "ERROR AIDMMaterial NOT AVAILABLE! \n";
    return nullptr;
}

std::shared_ptr<KerasModel> AIDMMaterial::keras_SANN_sp(new KerasModel("AIDMBB.model"));
std::shared_ptr<KerasModel> AIDMMaterial::keras_HANN_sp(new KerasModel("AIDMHY.model"));

std::vector<double> AIDMMaterial::lammda_vec = { 0, 5, 1 };
std::vector<double> AIDMMaterial::lammdaS_vec = { 0, 2, 0.6 };
std::vector<double> AIDMMaterial::lammdaSV_vec = { 0, 0.3, 0.6 };
std::vector<double> AIDMMaterial::lammdaT_vec = { 0, 6, 0.4 };
std::vector<double> AIDMMaterial::strainC_vec = { 0, 0.06, 1 };
std::vector<double> AIDMMaterial::stressFactor_vec = { 0.4, 1.3, 1 };
std::vector<double> AIDMMaterial::m_vec = { 2, 9, 1 };
std::vector<double> AIDMMaterial::n_vec = { 1, 8, 1 };
std::vector<double> AIDMMaterial::secantK_vec = { 0, 5, 0.6 };
std::vector<double> AIDMMaterial::afa_vec = { 0, 5, 1 };
std::vector<double> AIDMMaterial::beta_vec = { 0, 3, 0.4 };
std::vector<double> AIDMMaterial::gamma_vec = { 0, 1, 1 };
std::vector<double> AIDMMaterial::eta_vec = { 0, 1, 6 };
int AIDMMaterial::predict_num = 0;

AIDMMaterial::AIDMMaterial(const double& lammdaSV, const double& lammdaS, const double& lammdaT_pos)
:TStrain(0.0), CStress(0.0), TStress(0.0), CStrain(0.0),
K(0), CK(0.0),
CstrainMax(0), CstrainMin(0), CstressMaxFactor(0), CstressMinFactor(0),
lammda(initialLammda), lammdaS(lammdaS), lammdaSV(lammdaSV), lammdaT_pos(lammdaT_pos),
needUpdateHANN_pos(false), needUpdateHANN_neg(false), needUpdateBANN(true),
afa_pos(0.0), afa_neg(0.0), isConstant(false),TLoadingDirectPos(true), mnFactor_(1),
stressSA_pos(0), stressSA_neg(0),axialRatio(0), eta_pos(0), eta_neg(0), 
dalr_sp(new DALRMaterial()), TLoadingTag(0)
{
}


AIDMMaterial::AIDMMaterial()
: TStrain(0.0), CStress(0.0), TStress(0.0), CStrain(0.0),
K(0), CK(0.0),
CstrainMax(0), CstrainMin(0), CstressMaxFactor(0), CstressMinFactor(0),
lammda(initialLammda), lammdaS(0.0), lammdaSV(0.0), lammdaT_pos(0.0),
needUpdateHANN_pos(false), needUpdateHANN_neg(false), needUpdateBANN(true),
afa_pos(0.0), afa_neg(0.0), isConstant(false), TLoadingDirectPos(true), mnFactor_(1),
stressSA_pos(0), stressSA_neg(0), axialRatio(0), eta_pos(0), eta_neg(0), 
dalr_sp(new DALRMaterial()), TLoadingTag(0)
{

}


AIDMMaterial::~AIDMMaterial()
{
    if (AIDMMaterial::predict_num != 0)
    {
        opserr << "Neural networks were utilizied for " << AIDMMaterial::predict_num << " times\n";
        AIDMMaterial::predict_num = 0;
    }
  // does nothing
    

}

int 
AIDMMaterial::setTrialStrain(double strain, double strainRate)
{
    //is elastic 单元为弹性
    if (this->isConstant || !this->isAvailabelAIDM())
    {
        TStrain = strain;
        TStress = TStrain * K;
        return 0;
    }
    //割线法
    //if (this->K != this->CK)
    //{
    //    TStrain = strain;
    //    TStress = this->CStress + K * (strain - this->CStrain);
    //}


    //判断加载方向 产生非常小的变形才判断加载方向
    if (abs(CStrain - strain) > 1E-16)  TLoadingDirectPos = CStrain < strain;
    //未发生任何加载 且应变都不为0（非初始状态）
    else if (strain != 0 && CStrain != 0)  return 0;
    //history maximum strain 获得历史最大变形
    auto max_strain = TLoadingDirectPos ? CstrainMax : CstrainMin;
    //strain larger than history maximum value -> go on backbone  
    //超过历史最大变形 指向骨架曲线
    if ((TLoadingDirectPos ? strain >= max_strain : strain <= max_strain))
    {
        this->setTangentOnBackbone(strain, TLoadingDirectPos);
        return 0;
    }
    //Unloading 如果应变呈现卸载趋势（同向，绝对值小） 获得卸载刚度
    else if (abs(strain) < abs(CStrain) && CStress * CStrain > 0)
    {
        auto info  = this->setUnloadingTangent(strain);
        if(info == 0) return 0;
    }
    //Reload path not created 如果没有历史最大变形 则指向骨架
    if (max_strain == 0)
    {
        this->setTangentOnBackbone(strain, TLoadingDirectPos);
    }
    else
    {
        this->setReloadingTangent(strain, TLoadingDirectPos);
    }

    return 0;
}

void AIDMMaterial::ResetTrialStrain(const double& strain, const int& loadingTag)
{
    //is elastic 单元为弹性
    if (this->isConstant || !this->isAvailabelAIDM())
    {
        TStrain = strain;
        TStress = TStrain * K;
    }
    //骨架
    else if (loadingTag == 1)
        this->setTangentOnBackbone(strain, TLoadingDirectPos);
    //重加载
    else if (loadingTag == 2)
        this->setReloadingTangent(strain, TLoadingDirectPos);
    //卸载
    else if (loadingTag == 3)
        this->setUnloadingTangent(strain);
}

double AIDMMaterial::GetSoftenFactor(bool isPos)
{
    auto strainMax = isPos ? this->CstrainMax : this->CstrainMin;
    auto Mcstrain = isPos ? this->strainC_pos : this->strainC_neg;
    auto stressMax = isPos ? this->stressSA_pos : this->stressSA_neg;
    //判断是否软化
    if (std::abs(strainMax) <= std::abs(Mcstrain)) return 1;
    //计算承载力
    auto stress = this->getStressOnBackbone(strainMax, this->mnFactor_);
    return std::abs(stress) / stressMax;
}

int 
AIDMMaterial::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
    this->setTrialStrain(strain, strainRate);
    stress = TStress;
    tangent = K;
    return 0;
}

void AIDMMaterial::setTangentOnBackbone(const double& strain, bool loading_direct_pos)
{
    //Update params 更新骨架参数
    this->updateSkeletonParams();
    //Oriented Strain Stress 指向的应变（初值）
    auto oriented_strain = backbone_inidStrainFactor * (loading_direct_pos ? strainC_pos : -strainC_neg);
    //In same side 当应变较大时，根据上次提交的应变确定指向应变
    if ((loading_direct_pos ? strain >= 0 : strain <= 0))
    {
        //通过strain确定刚度 否则刚度不发生变化
        oriented_strain = abs(oriented_strain) > abs(this->CStrain * backbone_ortStrainFactor) ?
            oriented_strain : this->CStrain * backbone_ortStrainFactor;
        auto strainmax = loading_direct_pos ? this->CstrainMax : this->CstrainMin;
        //判断指向的位移（离最大变形值太远 则仍走重加载路线）
        if (std::abs(strain) < std::abs(strainmax) && std::abs(oriented_strain) < std::abs(strainmax))
            return;
    }
    //更新加载模式
    this->TLoadingTag = 1;
    // strain out of oriented 如果指向应变小于目标应变  则指向应变改为目标应变
    if((loading_direct_pos? oriented_strain < strain: oriented_strain > strain))
        oriented_strain = strain;
    //通过应变获得骨架上的应力
    auto oriented_stress = this->getStressOnBackbone(oriented_strain, this->mnFactor_);

    /*智能使用Cstrain*/
    auto strain0 = this->CStrain;
    auto stress0 = this->CStress;
    //Updating stiffness 如果没有应变增量 则取上一步采用的刚度值
    if (oriented_strain - strain0 == 0) return;
    //Updating stiffness  更新刚度
    K = (oriented_stress - stress0) / (oriented_strain - strain0);
    // Updating TStrain TStress 更新该迭代步下的应力应变
    TStrain = strain;
    TStress = stress0 + K * (strain - strain0);
    //is kill the element 变形超过峰值变形（有上升下降）
    bool isPost = abs(strain) > (loading_direct_pos ? strainC_pos : strainC_neg);
    //峰值承载力
    auto backboneMaxStress = loading_direct_pos ? stressFactor_pos * stressSA_pos :
        stressFactor_neg * stressSA_neg;
    //残余应力系数是否越界
    bool isStressBelow = abs(this->TStress) / abs(backboneMaxStress) < this->killStressFactor;
    //软化段、残余应力越界、负刚度 则认为构件可以被杀死
    if (K < 0 && isPost && isStressBelow)
    {
        if (this->needToKill)
        {
            K = 0;
            TStress = 0;
        }
        else
        {
            TStress = (TStress < 0 ? -1 : 1) * backboneMaxStress * this->killStressFactor;
            K = (TStress - stress0) / (TStrain - strain0);
        }
    }
    
}

int AIDMMaterial::setUnloadingTangent(const double& strain)
{
    //更新加载模式
    this->TLoadingTag = 3;
    //Update params 更新滞回参数
    this->updateHystereticParams(strain > 0);
    //Maximum unloading stiffness 最大卸载刚度 指向零点
    auto min_unloading_K = CStress / CStrain;
    //Oriented Strain Stress 最大卸载刚度也不应大于初始刚度（针对y轴附近的震荡）
    auto oriented_strain = backbone_inidStrainFactor * (strain > 0 ? strainC_pos : -strainC_neg);
    auto max_unloading_K = this->getStressOnBackbone(oriented_strain, this->mnFactor_) / oriented_strain;
    //Caping secant stiffness 峰值点的割线刚度系数
    auto secant_Kc = CStress > 0 ? stressFactor_pos * stressSA_pos / strainC_pos :
        stressFactor_neg * stressSA_neg / strainC_neg;
    //Unloading stiffness 计算卸载刚度
    auto unloading_K = CStress > 0 ? secant_Kc * afa_pos : secant_Kc * afa_neg;
    //尚未达到更新滞回参数的条件
    if (unloading_K == 0 && min_unloading_K < secant_Kc)
    {
        /*采用前一分析步的K值*/
    }
    //Unloading stiffness is to rigid 卸载刚度不可越界
    else
    {
        K = unloading_K < min_unloading_K ? min_unloading_K : unloading_K;
        K = K > max_unloading_K ? max_unloading_K : K;
    }
    //Updating TStrain TStress 更新应力应变
    auto stress = CStress + K * (strain - CStrain);
    //is balance 如果超过卸载边界 获得重加载刚度
    if (stress * strain <= 0)
        return -1;
    TStrain = strain;
    TStress = stress;
    return 0;
}

void AIDMMaterial::setReloadingTangent(const double& strain, bool loading_direct_pos)
{
    //更新加载模式
    this->TLoadingTag = 2;
    //Update params 更新滞回参数
    this->updateHystereticParams(loading_direct_pos);
    //Maximum deformation 最大变形
    auto max_strain = loading_direct_pos ? CstrainMax : CstrainMin;
    //最大变形对应的应力
    auto max_stress = loading_direct_pos ? CstressMaxFactor * stressFactor_pos * stressSA_pos :
        CstressMinFactor * stressFactor_neg * stressSA_neg;
    //prk point capacity 拐点应力
    auto brk_prt_stress = max_stress * (loading_direct_pos ? gamma_pos : gamma_neg);
    //oriented point capacity 指向点应力
    auto oreinted_prt_stress = max_stress * (loading_direct_pos ? eta_pos : eta_neg);
    //Caping secant stiffness 割线刚度
    auto secant_Kc = loading_direct_pos ? stressFactor_pos * stressSA_pos / strainC_pos :
        stressFactor_neg * stressSA_neg / strainC_neg;
    //break proint maxstress with sign 二段重加载刚度不可大于初始刚度
    auto initial_strain_pos = loading_direct_pos? this->backbone_inidStrainFactor * strainC_pos: -this->backbone_inidStrainFactor * strainC_neg;
    auto iniK = this->getStressOnBackbone(initial_strain_pos, this->mnFactor_) / initial_strain_pos;
    auto brk_prt_strain_max = max_strain - (oreinted_prt_stress - brk_prt_stress) / iniK;
    //towards to oriented prt: 是否已经过拐点应力
    bool isLargerthanBrkprtStress = loading_direct_pos ? CStress > brk_prt_stress: CStress < brk_prt_stress;
    //是否已经经过 拐点最大应变
    bool isbeyondBrkprtStrainMax = loading_direct_pos ? CStrain > brk_prt_strain_max : CStrain < brk_prt_strain_max;
    //拐点应力是否超过指向点应力
    bool isBrkprtBeyondOrientedPrt = abs(brk_prt_stress) > abs(oreinted_prt_stress);
    //Oriented stiffness 计算二段重加载刚度
    auto reloadingB_K = (oreinted_prt_stress - CStress) / (max_strain - CStrain);
    //如果经过 拐点最大应变 或 拐点应力 或 拐点应力超过指向点应力
    //二段重加载
    if (isLargerthanBrkprtStress || isbeyondBrkprtStrainMax || isBrkprtBeyondOrientedPrt)
    {
        //应变无发生任何变化 ？ 
        if (max_strain - CStrain == 0) return;
        K = reloadingB_K;
    }
    //maybe towards to break prt 指向拐点
    else
    {
        //First reloading stiffness 计算第一段重加载刚度
        auto reloadingA_K = loading_direct_pos ? secant_Kc * beta_pos : secant_Kc * beta_neg;
        //Reloading target strain 拐点真实应变
        auto reloadingA_strain = CStrain + (brk_prt_stress - CStress) / reloadingA_K;
        //拐点真实应变太靠后（2加载段刚度太大不允许）  当前目标应变要小于拐点真实应变
        bool isStrainBeyond = loading_direct_pos ? (reloadingA_strain > brk_prt_strain_max || reloadingA_strain < strain):
           ( reloadingA_strain < brk_prt_strain_max || reloadingA_strain > strain);
        //超过采用二段应变
        K = isStrainBeyond ? reloadingB_K : reloadingA_K;
    }
    // Updating TStrain TStress
    TStrain = strain;
    TStress = CStress + K * (strain - CStrain);
    return;
}

#pragma region 获得输入层和输出层参数

std::vector<double>& AIDMMaterial::getRegularizedValueVector(const AIDMParamEnum& type)
{
    switch (type)
    {
    case AIDMParamEnum::Lammda: return AIDMMaterial::lammda_vec;
    case AIDMParamEnum::LammdaS: return AIDMMaterial::lammdaS_vec;
    case AIDMParamEnum::LammdaSV: return AIDMMaterial::lammdaSV_vec;
    case AIDMParamEnum::LammdaT: return AIDMMaterial::lammdaT_vec;

    case AIDMParamEnum::StressFactor: return AIDMMaterial::stressFactor_vec;
    case AIDMParamEnum::StrainC: return AIDMMaterial::strainC_vec;
    case AIDMParamEnum::m: return AIDMMaterial::m_vec;
    case AIDMParamEnum::n: return AIDMMaterial::n_vec;

    case AIDMParamEnum::SecantK: return AIDMMaterial::secantK_vec;

    case AIDMParamEnum::Afa: return AIDMMaterial::afa_vec;
    case AIDMParamEnum::Beta: return AIDMMaterial::beta_vec;
    case AIDMParamEnum::Gamma: return AIDMMaterial::gamma_vec;
    case AIDMParamEnum::Eta: return AIDMMaterial::eta_vec;

    default: return AIDMMaterial::secantK_vec;
    }
}

float AIDMMaterial::getRegularizedValue(const double& value, const AIDMParamEnum& type)
{
    auto& vec = this->getRegularizedValueVector(type);
    return pow((value - vec[0]) / (vec[1] - vec[0]), vec[2]);
}

float AIDMMaterial::getNormalValue(const double& value, const AIDMParamEnum& type)
{
    if (type == AIDMParamEnum::Afa || type == AIDMParamEnum::Beta ||
        type == AIDMParamEnum::Gamma || type == AIDMParamEnum::Eta)
    {
        if (value < 0)
            return 0.001;
    }
    auto& vec = getRegularizedValueVector(type);
    auto powValue = pow(value, 1 / vec[2]);
    auto normalValue = powValue * (vec[1] - vec[0]) + vec[0];



    if (type == AIDMParamEnum::m || type == AIDMParamEnum::n)
    {
        normalValue = normalValue > 1 ? normalValue : 1.5;
    }
    else  if (type == AIDMParamEnum::StressFactor)
    {
        normalValue = normalValue > 1.2 ? 1.2 : normalValue;
        normalValue = normalValue < 0.5 ? 0.5 : normalValue;
    }
    else if (type == AIDMParamEnum::StrainC)
    {
        normalValue = normalValue < 0.005 ? 0.005 : normalValue;
    }

    else if (type == AIDMParamEnum::Afa || type == AIDMParamEnum::Beta)
    {
        if (value < 0)
            return 0.001;
    }

    else if (type == AIDMParamEnum::Gamma || type == AIDMParamEnum::Eta)
    {
        if (value < 0)
            return 0.001;
        else if (value > 1)
            return 1;
    }

    return normalValue;
}

std::vector<float> AIDMMaterial::getComponentParamsVec(bool is_pos)
{
    std::vector<float> vec = {};
    auto lammda_reg = this->getRegularizedValue(this->lammda, AIDMParamEnum::Lammda);
    auto lammdas_reg = this->getRegularizedValue(lammdaS, AIDMParamEnum::LammdaS);
    auto lammdasv_reg = this->getRegularizedValue(lammdaSV, AIDMParamEnum::LammdaSV);
    auto lammdat_reg = this->getRegularizedValue(is_pos ? lammdaT_pos : 1 / lammdaT_pos, AIDMParamEnum::LammdaT);
    return { lammda_reg, lammdas_reg,  lammdasv_reg, lammdat_reg };
}

#pragma endregion

void AIDMMaterial::updateSkeletonParams()
{
    //Upadte Backbone curves
    if (!this->isAvailabelAIDM() || !needUpdateBANN) return;
    //Input params
    auto& componentParamvec_pos = this->getComponentParamsVec(true);
    auto& componentParamvec_neg = this->getComponentParamsVec(false);
    //Predict
    auto& outputParams_pos = AIDMMaterial::keras_SANN_sp->predict(componentParamvec_pos);
    auto& outputParams_neg = AIDMMaterial::keras_SANN_sp->predict(componentParamvec_neg);
    //Getvalue
    m_pos = this->getNormalValue(outputParams_pos[0], AIDMParamEnum::m);
    n_pos = this->getNormalValue(outputParams_pos[1], AIDMParamEnum::n);
    strainC_pos = this->getNormalValue(outputParams_pos[2], AIDMParamEnum::StrainC);
    stressFactor_pos = this->getNormalValue(outputParams_pos[3], AIDMParamEnum::StressFactor);
    m_neg = this->getNormalValue(outputParams_neg[0], AIDMParamEnum::m);
    n_neg = this->getNormalValue(outputParams_neg[1], AIDMParamEnum::n);
    strainC_neg = this->getNormalValue(outputParams_neg[2], AIDMParamEnum::StrainC);
    stressFactor_neg = this->getNormalValue(outputParams_neg[3], AIDMParamEnum::StressFactor);
    needUpdateBANN = false;
    AIDMMaterial::predict_num += 2;
}

void AIDMMaterial::updateHystereticParams(bool is_pos)
{
    //Try to update skeleton
    if (!this->isAvailabelAIDM()) return;
    this->updateSkeletonParams();
    //is need to update
    if ((is_pos ? !needUpdateHANN_pos : !needUpdateHANN_neg))
        return;
    //Maximum deformation
    auto max_strain = is_pos ? CstrainMax : CstrainMin;
    auto max_stress = is_pos ? CstressMaxFactor * stressFactor_pos * stressSA_pos :
        CstressMinFactor * stressFactor_neg * stressSA_neg;
    auto skeleton_stress = this->getStressOnBackbone(max_strain, this->mnFactor_);
    //Caping secant stiffness
    auto secant_Kc = is_pos ? stressFactor_pos * stressSA_pos / strainC_pos :
        stressFactor_neg * stressSA_neg / strainC_neg;
    //Input params
    auto secantK = max_stress / max_strain;
    //is Pre capping
    bool isPreCapping = (is_pos ? strainC_pos : strainC_neg) > abs(max_strain);
    //Not Reach Max deformation Prt
    if (max_strain == 0)
    {
        //Oriented Strain Stress
        auto oriented_strain = backbone_inidStrainFactor * (is_pos ? strainC_pos : -strainC_neg);
        secantK = this->getStressOnBackbone(oriented_strain, this->mnFactor_) / (oriented_strain);
    }
    //post capping ignore the stiffness dagradation
    else if (isPreCapping && abs(max_stress) / abs(skeleton_stress) < 0.8)
        secantK = this->getStressOnBackbone(max_strain, this->mnFactor_) / (max_strain);
    auto secantKFactor = this->getRegularizedValue(secantK / secant_Kc, AIDMParamEnum::SecantK);
    auto& componentParamvec = this->getComponentParamsVec(is_pos);
    componentParamvec.push_back(secantKFactor);
    //Predict
    auto& outputParams = AIDMMaterial::keras_HANN_sp->predict(componentParamvec);
    //Getvalue
    if (is_pos)
    {
        afa_pos = this->getNormalValue(outputParams[0], AIDMParamEnum::Afa);
        beta_pos = this->getNormalValue(outputParams[1], AIDMParamEnum::Beta);
        gamma_pos = this->getNormalValue(outputParams[2], AIDMParamEnum::Gamma);
        eta_pos = this->getNormalValue(outputParams[3], AIDMParamEnum::Eta);
        needUpdateHANN_pos = false;
    }
    else
    {
        afa_neg = this->getNormalValue(outputParams[0], AIDMParamEnum::Afa);
        beta_neg = this->getNormalValue(outputParams[1], AIDMParamEnum::Beta);
        gamma_neg = this->getNormalValue(outputParams[2], AIDMParamEnum::Gamma);
        eta_neg = this->getNormalValue(outputParams[3], AIDMParamEnum::Eta);
        needUpdateHANN_neg = false;
    }
    AIDMMaterial::predict_num += 1;
}

double 
AIDMMaterial::getStress(void)
{
    return TStress;
}


double 
AIDMMaterial::getTangent(void)
{
    return K;
}


double 
AIDMMaterial::getInitialTangent(void)
{
    if (!this->isAvailabelAIDM()) return this->K;
    auto pos_dStrain = this->backbone_iniKFactor * strainC_pos;
    auto neg_dStrain = this->backbone_iniKFactor * strainC_neg;
    auto K0_pos = this->getStressOnBackbone(pos_dStrain, this->mnFactor_) / pos_dStrain;
    auto K0_neg = this->getStressOnBackbone(-neg_dStrain, this->mnFactor_) / (-neg_dStrain);
    return K0_pos > K0_neg ? K0_pos : K0_neg;
}


int 
AIDMMaterial::commitState(void)
{
    //是否软化状态
    bool isSoften = abs(CStress) > abs(TStress);
    //提交应力应变、刚度、加载方向
    CStrain = TStrain;
    CStress = TStress;
    CK = K;
    CLoadingDirectPos = TLoadingDirectPos;
    //Need to kill 判断是否为无效单元或单元已被杀死
    if (this->isConstant || !this->isAvailabelAIDM()) return 0;

    //Record Maxmum deformation
    if (CStrain > 0)
    {
        // not the maximum deformation 未达到历史最大变形
        if (abs(CstrainMax) > abs(CStrain) || abs(CStrain) < backbone_inidStrainFactor * strainC_pos)
            return 0;
        // is the new maximum deformation 达到历史最大变形
        //获得新的历史最大变形 获得最大变形点的应力系数 滞回参数需要被更新
        CstrainMax = CStrain;
        CstressMaxFactor = CStress / (stressFactor_pos * stressSA_pos);
        needUpdateHANN_pos = true;
        //处于软化段
        bool isPost = abs(CstrainMax) > abs(strainC_pos);
        //承载力低于限值要求
        bool isCapacitybelow = abs(CstressMaxFactor) < this->killStressFactor;
        // 处于软化段 承载力低于限值要求 承载力在退化过程 （为何不通过负刚度做判断）
        if (isSoften && isPost && isCapacitybelow) 
            this->Kill();
    }
    else
    {
        // not the maximum deformation  未达到历史最大变形
        if (abs(CstrainMin) > abs(CStrain) || abs(CStrain) < backbone_inidStrainFactor * strainC_neg)
            return 0;
        // is the new maximum deformation 达到历史最大变形
        //获得新的历史最大变形 获得最大变形点的应力系数 滞回参数需要被更新
        CstrainMin = CStrain;
        CstressMinFactor = CStress / (stressFactor_neg * stressSA_neg);
        needUpdateHANN_neg = true;
        // is post 处于软化段
        bool isPost = abs(CstrainMin) > abs(strainC_neg);
        //承载力低于限值要求
        bool isCapacitybelow = abs(CstressMinFactor) < this->killStressFactor;
        // 处于软化段 承载力低于限值要求 承载力在退化过程 （为何不通过负刚度做判断）
        if (isSoften && isPost && isCapacitybelow) 
            this->Kill();
    }
    return 0;
}


int 
AIDMMaterial::revertToLastCommit(void)
{
    /*迭代不收敛时调用的方法  返回前一个状态*/
    this->TStrain = CStrain;
    this->TStress = CStress;
    this->K = CK;
    this->TLoadingDirectPos = CLoadingDirectPos;
    this->needUpdateBANN = true;
    this->needUpdateHANN_neg = true;
    this->needUpdateHANN_pos = true;
    return 0;
}


int 
AIDMMaterial::revertToStart(void)
{
    this->TStrain = 0.0;
    this->TStress = 0.0;
    this->K = 0.0;

    this->lammda = this->initialLammda;
    this->axialRatio = 0;

    this->needUpdateBANN = true;
    this->needUpdateHANN_neg = false;
    this->needUpdateHANN_pos = false;
    this->isConstant = false;
    this->TLoadingDirectPos = true;
  
    this->mnFactor_ = 1;
    this->CstrainMax = 0;   this->CstrainMin = 0;
    this->CstressMaxFactor = 0;   this->CstressMinFactor = 0;

    return 0;
}


AIDMMaterial*
AIDMMaterial::getCopy(void)
{
    AIDMMaterial* mat = new AIDMMaterial(this->lammdaSV, this->lammdaS, this->lammdaT_pos);
    mat->setARK(this->dalr_sp->getARK(), this->dalr_sp->getFactor());
    return mat;
}

void AIDMMaterial::setLammda(const double& Lammda)
{
    if (this->lammda != this->initialLammda)
    {
        if (!this->isValidUpdate()) return;
    }
    auto target_lammda = Lammda < 0.5 ? 0.5 : Lammda;
    target_lammda = target_lammda > 5 ? 5 : target_lammda;
    if (abs(this->lammda - target_lammda) > 0.25)
    {
        this->lammda = target_lammda;
        needUpdateBANN = true;
        if(CstrainMax != 0)
            needUpdateHANN_pos = true;
        if(CstrainMin != 0)
            needUpdateHANN_neg = true;
    }
}

void AIDMMaterial::setCapcacity(const double& m_pos, const double& m_neg, const double& axial_ratio)
{
    if (this->stressSA_pos != 0 && this->stressSA_neg != 0)
    {
        if (!this->isValidUpdate()) return;
    }
    auto target_AR = axial_ratio < -1 ? -1 : axial_ratio;
    target_AR = target_AR > 1 ? 1 : target_AR;
    if (abs(this->axialRatio - target_AR) > 0.05)
    {
        this->axialRatio = target_AR;
        this->stressSA_pos = m_pos;
        this->stressSA_neg = m_neg;
    }
    else if (this->stressSA_pos == 0 || this->stressSA_neg == 0)
    {
        this->axialRatio = target_AR;
        this->stressSA_pos = m_pos;
        this->stressSA_neg = m_neg;
    }
    else
    {
        if (std::abs(m_pos - this->stressSA_pos) / this->stressSA_pos > 0.05)
            this->stressSA_pos = m_pos;
        if (std::abs(m_neg - this->stressSA_neg) / this->stressSA_neg > 0.05)
            this->stressSA_neg = m_neg;
    }
}

bool AIDMMaterial::checkCapacity(const double& Moment)
{
    /*非有效的AIDM对象*/
    if (!this->isAvailabelAIDM()) return true;
    /*承载力限值系数  超过转为弹性*/
    auto limitFactor = 0.8;
    /*计算构件承载力*/
    auto capacity = Moment < 0 ? this->stressSA_neg * this->stressFactor_neg : 
         this->stressSA_pos * this->stressFactor_pos;
    /*承载力满足需求*/
     if (abs(Moment) < capacity * limitFactor)  return true;
     //刚度转为弹性 刚度为割线刚度 
     this->isConstant = true;
     K = this->CStress / this->CStrain;
     return false;
}

double AIDMMaterial::getStressOnBackbone(const double& drift, const int& mnFactor)
{
    auto m = drift > 0 ? m_pos : m_neg;
    auto n = drift > 0 ? n_pos : n_neg;
    auto dc = drift > 0 ? strainC_pos : strainC_neg;
    //interpolation for m value
    if (abs(drift) < dc)
    {
        auto mMax = m * mnFactor;
        auto nMin = n / mnFactor < 1 ? 1.5 : n / mnFactor;
        auto coe_BM = mMax;
        auto coe_AM = (m - coe_BM) / pow(dc, 1);
        auto coe_BN = nMin;
        auto coe_AN = (n - coe_BN) / pow(dc, 1);
        m = coe_AM * abs(pow(drift, 1)) + coe_BM;
        n = coe_AN * abs(pow(drift, 1)) + coe_BN;
    }
    auto capacity = drift > 0 ? stressFactor_pos * stressSA_pos : stressFactor_neg * stressSA_neg;
    auto x = abs(drift / dc);
    auto y = (m * x) / (1 + (m - (n / (n - 1))) * x + pow(x, n) / (n - 1));
    return drift > 0? y * capacity : -y * capacity;
}

double AIDMMaterial::GetALR()
{
    //初始状态
    if (this->TLoadingTag == 0) return 0;
    //峰值承载力及其对应的变形
    auto dc = this->TLoadingDirectPos ? this->strainC_pos : -this->strainC_neg;
    auto stressMax = this->TLoadingDirectPos ? (stressFactor_pos * stressSA_pos) :
        (stressFactor_neg * stressSA_neg);
    //骨架
    if (this->TLoadingTag == 1)
    {
        //承载力退化
        auto CStressFactor = this->CStress / stressMax;
        return -this->dalr_sp->GetALROnBackbone(this->TStrain, dc, this->lammda, CStressFactor);
    }
    //重加载
    else if (this->TLoadingTag == 2)
    {
        //最大变形点
        auto strainMax = this->TLoadingDirectPos ? this->CstrainMax : this->CstrainMin;
        //承载力退化
        auto CStressFactor = this->TLoadingDirectPos ? (this->CstressMaxFactor * this->eta_pos) :
            (this->CstressMinFactor * this->eta_neg);
        return -this->dalr_sp->GetReloadALR(this->TStrain, dc, strainMax, this->lammda, CStressFactor);
    }
    else if (this->TLoadingTag == 3)
    {
        return -this->dalr_sp->GetUnloadALR(this->TStrain);
    }
    else
    {
        opserr << "Unknow TLoadingTag For AIDMMaterial" << endln;
        return 0;
    }
}

void AIDMMaterial::setInitialK(const double iniK, bool ensureIniK)
{
    //是否为有效单元
    if (!this->isAvailabelAIDM())
    {
        this->isConstant = true;
        this->K = iniK; return;
    }
    //不需要硬化单元
    if (!ensureIniK) return;
    /*初始化调整系数*/
    this->mnFactor_ = 0;
    /*用于计算初始刚度的应变*/
    auto initial_strain_pos = this->backbone_iniKFactor * strainC_pos;
    auto initial_strain_neg = this->backbone_iniKFactor * -strainC_neg;
    /*刚度系数*/
    auto factorpos = 0.0;
    auto factorneg = 0.0;
    while (factorpos < 1.0 && factorneg < 1.0)
    {
        /*第一次进入 等于1*/
        this->mnFactor_ += 1;
        /*计算初始刚度*/
        auto initial_strain_K_pos = this->getStressOnBackbone(initial_strain_pos, this->mnFactor_) / initial_strain_pos;
        auto initial_strain_K_neg = this->getStressOnBackbone(initial_strain_neg, this->mnFactor_) / initial_strain_neg;
        /*更新刚度系数*/
        factorpos = initial_strain_K_pos / iniK;
        factorneg = initial_strain_K_neg / iniK;
        /*设定上限值 不然将导致曲线畸形*/
        if (this->mnFactor_ == 5) return;
    }
}

bool AIDMMaterial::isAvailabelAIDM() const
{
    if (this->lammdaS == 0 || this->lammdaSV == 0 || this->lammdaT_pos == 0)
        return false;
    return true;
}

bool AIDMMaterial::isKill() const
{
    //是否为有效的AIDM
    if (!this->isAvailabelAIDM()) return false;
    //常刚度 且刚度为零
    return (this->isConstant && this->K == 0);
}

void AIDMMaterial::Kill()
{
    //是否为有效的AIDM
    if (!this->isAvailabelAIDM()) return;
    //不需要被杀死
    if (!this->needToKill) return;
    //常刚度 且刚度为零
    this->isConstant = true;
    this->K = 0;
}

bool AIDMMaterial::isValidUpdate()
{
    //如果还未超过最大变形
    if (this->TLoadingDirectPos ? (CstrainMax == 0) : (CstrainMin == 0))
        return false;
    //变形是否合法
    auto max_strain = this->TLoadingDirectPos ? CstrainMax : CstrainMin;
    auto limit_strain = max_strain * this->updateDeformationFactor;
    //变形是否超过限值
    bool isStrainOut = this->TLoadingDirectPos ? (this->TStrain > limit_strain) : (this->TStrain < limit_strain);
    //最大变形对应的应力
    auto max_stress = this->TLoadingDirectPos ? CstressMaxFactor * stressFactor_pos * stressSA_pos :
        CstressMinFactor * stressFactor_neg * stressSA_neg;
    //oriented point capacity 指向点应力
    auto oreinted_prt_stress = max_stress * (this->TLoadingDirectPos ? eta_pos : eta_neg);
    auto limit_stress = oreinted_prt_stress * this->updateForceFactor;
    //判断承载力是否超过更新限值
    bool isForceOut = this->TLoadingDirectPos ? (this->TStress > limit_stress) : (this->TStress < limit_stress);
    return isStrainOut || isForceOut;
}
