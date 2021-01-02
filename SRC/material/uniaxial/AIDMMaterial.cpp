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

#ifdef _SAP
    return nullptr;
#else
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? <eta?> <Eneg?> ... " << endln;
		return 0;
	}

	int iData[1];
	double dData[7];
	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
		return 0;
	}

	numData = OPS_GetNumRemainingInputArgs();

	if (numData >= 7) {
		numData = 7;
		if (OPS_GetDoubleInput(&numData, dData) != 0) {
			opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
			return 0;
		}
	}
	else {
        numData = 6;
        if (OPS_GetDoubleInput(&numData, dData) != 0) {
            opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
            return 0;
        }
	}
    if(numData == 7)
        theMaterial = new AIDMMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);
	// Parsing was successful, allocate the material
	else theMaterial = new AIDMMaterial(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type ElasticMaterial\n";
		return 0;
	}
	return theMaterial;
#endif // SAP
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

AIDMMaterial::AIDMMaterial(int tag, double height, double width, double lammdaS, double lammdaSV, double lammdaT_pos, 
    double Msa_pos, double Msa_neg, bool ensureIniK)
:UniaxialMaterial(tag, MAT_TAG_AIDMMaterial),
TStrain(0.0), CStress(0.0), TStress(0.0), CStrain(0.0),
K(0), CK(0.0),
stressSA_pos(Msa_pos), stressSA_neg(Msa_neg),
CstrainMax(0), CstrainMin(0),
lammda(initialLammda), lammdaS(lammdaS), lammdaSV(lammdaSV), lammdaT_pos(lammdaT_pos),
needUpdateHANN_pos(false), needUpdateHANN_neg(false), needUpdateBANN(true), sectionHeight(height),
afa_pos(0.0), afa_neg(0.0), isKill(false), isToElastic(false), initialK(0.0), TLoadingDirectPos(true), ensureIniK(ensureIniK)
{
    this->updateSkeletonParams();
}


AIDMMaterial::AIDMMaterial(int tag, double lammda, double lammdaS, double lammdaSV, double lammdaT_pos, 
    double Msa_pos, double Msa_neg, bool ensureIniK)
    :UniaxialMaterial(tag, MAT_TAG_AIDMMaterial),
    TStrain(0.0), CStress(0.0), TStress(0.0), CStrain(0.0),
    K(0), CK(0.0),
    stressSA_pos(Msa_pos), stressSA_neg(Msa_neg),
    CstrainMax(0), CstrainMin(0),
    lammda(lammda), lammdaS(lammdaS), lammdaSV(lammdaSV), lammdaT_pos(lammdaT_pos),
    needUpdateHANN_pos(false), needUpdateHANN_neg(false), needUpdateBANN(true), sectionHeight(0.0),
    afa_pos(0.0), afa_neg(0.0), isKill(false), isToElastic(false), initialK(0.0), TLoadingDirectPos(true), ensureIniK(ensureIniK)
{
    this->updateSkeletonParams();
}


AIDMMaterial::AIDMMaterial()
:UniaxialMaterial(0, MAT_TAG_AIDMMaterial),
TStrain(0.0), CStress(0.0), TStress(0.0), CStrain(0.0),
K(0), CK(0.0),
stressSA_pos(0.0), stressSA_neg(0.0),
CstrainMax(0), CstrainMin(0),
lammda(initialLammda), lammdaS(0.0), lammdaSV(0.0), lammdaT_pos(0.0),
needUpdateHANN_pos(false), needUpdateHANN_neg(false), needUpdateBANN(true), sectionHeight(0.0),
afa_pos(0.0), afa_neg(0.0), isKill(false), isToElastic(false), initialK(0.0), TLoadingDirectPos(true), ensureIniK(false)
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
    //is elastic
    if (this->isToElastic)
    {
        TStrain = strain;
        TStress = TStrain * K;
        return 0;
    }
    //Ensure is null or not
    if (!this->isAvailabelAIDM() || this->isKill)
    {
        TStrain = strain;
        return 0;
    }
    if (abs(CStrain - strain) > 1E-16)
        TLoadingDirectPos = CStrain < strain;
    else if (strain != 0 && CStrain != 0)
        return 0;
    //history maximum strain
    auto max_strain = TLoadingDirectPos ? CstrainMax : CstrainMin;
    //strain larger than history maximum value -> go on backbone
    if ((TLoadingDirectPos ? strain >= max_strain : strain <= max_strain))
    {
        this->setTangentOnBackbone(strain, TLoadingDirectPos);
        return 0;
    }
    //Unloading
    else if (abs(strain) < abs(CStrain) && CStress * CStrain > 0)
    {
        auto info  = this->setUnloadingTangent(strain);
        if(info == 0) return 0;
    }
    //Reload path not created
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
    //Update params
    this->updateSkeletonParams();
    //Oriented Strain Stress
    auto oriented_strain = backbone_inidStrainFactor * (loading_direct_pos ? strainC_pos : -strainC_neg);
    //In same side
    if ((loading_direct_pos ? strain >= 0 : strain <= 0))
    {
        oriented_strain = abs(oriented_strain) > abs(CStrain * backbone_ortStrainFactor) ?
            oriented_strain : CStrain * backbone_ortStrainFactor;
    }
    // strain out of oriented
    if((loading_direct_pos? oriented_strain < strain: oriented_strain > strain))
        oriented_strain = strain;
    auto oriented_stress = this->getStressOnBackbone(oriented_strain);
    //Updating stiffness
    if (oriented_strain - CStrain == 0) return;
    K = (oriented_stress - CStress) / (oriented_strain - CStrain);
    // Updating TStrain TStress
    TStrain = strain;
    TStress = CStress + K * (strain - CStrain);
    //is kill the element
    bool isPost = abs(strain) > (loading_direct_pos ? strainC_pos : strainC_neg);
    auto backboneMaxStress = loading_direct_pos ? stressFactor_pos * stressSA_pos :
        stressFactor_neg * stressSA_neg;
    bool isStressBelow = abs(TStress) / abs(backboneMaxStress) < this->killStressFactor;
    if (K < 0 && isPost && isStressBelow)
    {
        K = 0;
        TStress = 0;
    }
}

int AIDMMaterial::setUnloadingTangent(const double& strain)
{
    //Update params
    this->updateHystereticParams(strain > 0);
    //Maximum unloading stiffness
    auto min_unloading_K = CStress / CStrain;
    //Oriented Strain Stress
    auto oriented_strain = backbone_inidStrainFactor * (strain > 0 ? strainC_pos : -strainC_neg);
    auto max_unloading_K = this->getStressOnBackbone(oriented_strain) / oriented_strain;
    //Caping secant stiffness
    auto secant_Kc = CStress > 0 ? stressFactor_pos * stressSA_pos / strainC_pos :
        stressFactor_neg * stressSA_neg / strainC_neg;
    //Unloading stiffness
    auto unloading_K = CStress > 0 ? secant_Kc * afa_pos : secant_Kc * afa_neg;
    //Unloading stiffness is to rigid
    K = unloading_K < min_unloading_K ? min_unloading_K : unloading_K;
    K = K > max_unloading_K ? max_unloading_K : K;
    //Updating TStrain TStress
    auto stress = CStress + K * (strain - CStrain);
    //is balance
    if (stress * strain <= 0)
        return -1;
    TStrain = strain;
    TStress = stress;
    return 0;
}

void AIDMMaterial::setReloadingTangent(const double& strain, bool loading_direct_pos)
{
    //Update params
    this->updateHystereticParams(loading_direct_pos);
    //Maximum deformation
    auto max_strain = loading_direct_pos ? CstrainMax : CstrainMin;
    auto max_stress = loading_direct_pos ? CstressMaxFactor * stressFactor_pos * stressSA_pos :
        CstressMinFactor * stressFactor_neg * stressSA_neg;
    //prk point capacity
    auto brk_prt_stress = max_stress * (loading_direct_pos ? gamma_pos : gamma_neg);
    //oriented point capacity
    auto oreinted_prt_stress = max_stress * (loading_direct_pos ? eta_pos : eta_neg);
    //Caping secant stiffness
    auto secant_Kc = loading_direct_pos ? stressFactor_pos * stressSA_pos / strainC_pos :
        stressFactor_neg * stressSA_neg / strainC_neg;
    //break proint maxstress with sign
    auto brk_prt_strain_max = max_strain - (oreinted_prt_stress - brk_prt_stress) / this->getInitialTangent();
    //towards to oriented prt:
    bool isLargerthanBrkprtStress = loading_direct_pos ? CStress > brk_prt_stress: CStress < brk_prt_stress;
    bool isbeyondBrkprtStrainMax = loading_direct_pos ? CStrain > brk_prt_strain_max : CStrain < brk_prt_strain_max;
    bool isBrkprtBeyondOrientedPrt = abs(brk_prt_stress) > abs(oreinted_prt_stress);
    //Oriented stiffness
    auto reloadingB_K = (oreinted_prt_stress - CStress) / (max_strain - CStrain);
    if (isLargerthanBrkprtStress || isbeyondBrkprtStrainMax || isBrkprtBeyondOrientedPrt)
    {
        if (max_strain - CStrain == 0) 
            return;
        K = reloadingB_K;
    }
    //maybe towards to break prt
    else
    {
        //First reloading stiffness
        auto reloadingA_K = loading_direct_pos ? secant_Kc * beta_pos : secant_Kc * beta_neg;
        //Reloading target strain
        auto reloadingA_strain = CStrain + (brk_prt_stress - CStress) / reloadingA_K;
        bool isStrainBeyond = loading_direct_pos ? reloadingA_strain > brk_prt_strain_max: reloadingA_strain < brk_prt_strain_max;
        K = isStrainBeyond ? reloadingB_K : reloadingA_K;
    }
    // Updating TStrain TStress
    TStrain = strain;
    TStress = CStress + K * (strain - CStrain);
    return;
}

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
    auto normalValue =  powValue * (vec[1] - vec[0]) + vec[0];



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
    auto lammdat_reg = this->getRegularizedValue(is_pos? lammdaT_pos : 1/ lammdaT_pos, AIDMParamEnum::LammdaT);
    return { lammda_reg, lammdas_reg,  lammdasv_reg, lammdat_reg };
}

void AIDMMaterial::updateSkeletonParams()
{
    //Upadte Backbone curves
    if (!needUpdateBANN) return;
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
    this->updateSkeletonParams();
    //is need to update
    if ((is_pos ? !needUpdateHANN_pos : !needUpdateHANN_neg))
        return;
    //Maximum deformation
    auto max_strain = is_pos ? CstrainMax : CstrainMin;
    auto max_stress = is_pos ? CstressMaxFactor * stressFactor_pos * stressSA_pos :
        CstressMinFactor * stressFactor_neg * stressSA_neg;
    auto skeleton_stress = this->getStressOnBackbone(max_strain);
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
        secantK = this->getStressOnBackbone(oriented_strain) / (oriented_strain);
    }
    //post capping ignore the stiffness dagradation
    else if (isPreCapping && abs(max_stress) / abs(skeleton_stress) < 0.8)
        secantK = this->getStressOnBackbone(max_strain) / (max_strain);
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
    auto pos_dStrain = this->backbone_inidStrainFactor * strainC_pos;
    auto neg_dStrain = this->backbone_inidStrainFactor * strainC_neg;
    auto K0_pos = this->getStressOnBackbone(pos_dStrain) / pos_dStrain;
    auto K0_neg = this->getStressOnBackbone(-neg_dStrain) / (-neg_dStrain);
    return K0_pos > K0_neg ? K0_pos : K0_neg;
}


int 
AIDMMaterial::commitState(void)
{
    bool isSoften = abs(CStress) > abs(TStress);
    CStrain = TStrain;
    CStress = TStress;
    CK = K;
    Clammda = lammda;
    CLoadingDirectPos = TLoadingDirectPos;
    //Need to kill
    if (this->isKill || !this->isAvailabelAIDM())
        return 0;
    //Record Maxmum deformation
    if (CStrain > 0)
    {
        /*if (this->initialK != 0 && abs(CStrain) > abs(strainC_pos * backbone_iniKStrainFactor))
        {
            this->initialK = 0;
        }*/
        // not the maximum deformation
        if (abs(CstrainMax) > abs(CStrain) || abs(CStrain) < backbone_inidStrainFactor * strainC_pos)
            return 0;
        // is the new maximum deformation
        CstrainMax = CStrain;
        CstressMaxFactor = CStress / (stressFactor_pos * stressSA_pos);
        needUpdateHANN_pos = true;
        // is post
        bool isPost = abs(CstrainMax) > abs(strainC_pos);
        bool isCapacitybelow = abs(CstressMaxFactor) < this->killStressFactor;
        // need to kill
        if (isSoften && isPost && isCapacitybelow)
        {
            isKill = true;
            K = 0;
            TStress = 0;
        }
    }
    else
    {
        /*if (this->initialK != 0 && abs(CStrain) > abs(strainC_neg * backbone_iniKStrainFactor))
        {
            this->initialK = 0;
        }*/
        // not the maximum deformation
        if (abs(CstrainMin) > abs(CStrain) || abs(CStrain) < backbone_inidStrainFactor * strainC_neg)
            return 0;
        // is the new maximum deformation
        CstrainMin = CStrain;
        CstressMinFactor = CStress / (stressFactor_neg * stressSA_neg);
        needUpdateHANN_neg = true;
        // is post
        bool isPost = abs(CstrainMin) > abs(strainC_neg);
        bool isCapacitybelow = abs(CstressMinFactor) < this->killStressFactor;
        // need to kill
        if (isSoften && isPost && isCapacitybelow)
        {
            isKill = true;
            K = 0;
            TStress = 0;
        }
    }
    return 0;
}


int 
AIDMMaterial::revertToLastCommit(void)
{
    TStrain = CStrain;
    TStress = CStress;
    K = CK;
    lammda = Clammda;
    TLoadingDirectPos = CLoadingDirectPos;
    needUpdateBANN = true;
    needUpdateHANN_neg = true;
    needUpdateHANN_pos = true;
    return 0;
}


int 
AIDMMaterial::revertToStart(void)
{
    TStrain = 0.0;
    TStress = 0.0;
    K = 0.0;
    lammda = initialLammda;
    return 0;
}


UniaxialMaterial *
AIDMMaterial::getCopy(void)
{
    AIDMMaterial*theCopy = new AIDMMaterial(this->getTag(), 
        lammda, lammdaS, lammdaSV, lammdaT_pos, stressSA_pos, stressSA_neg, this->ensureIniK);
    theCopy->TStrain = TStrain;
    theCopy->CStress = CStress;
    theCopy->CStrain = CStrain;
    theCopy->TStress = TStress;
    theCopy->CK = CK;
    theCopy->Clammda = Clammda;
    theCopy->sectionHeight = sectionHeight;
    return theCopy;
}


int 
AIDMMaterial::sendSelf(int cTag, Channel &theChannel)
{
    opserr << "AIDMMaterial::sendSelf() - failed to send data\n";
    return -1;
}


int 
AIDMMaterial::recvSelf(int cTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker)
{
    opserr << "AIDMMaterial::recvSelf() - failed to receive data\n";
    return -1;
}


void 
AIDMMaterial::Print(OPS_Stream &s, int flag)
{
	opserr << "AIDMMaterial::Print() - failed to print\n";
}


int
AIDMMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return -1;
}

int 
AIDMMaterial::updateParameter(int parameterID, Information &info)
{
    return -1;
}

void AIDMMaterial::setLammda(const double& shearSpan)
{
    //return;
    if (this->sectionHeight == 0)
        return;
    auto Lammda = abs(shearSpan) / this->sectionHeight;
    Lammda = Lammda < 0.5 ? 0.5 : Lammda;
    Lammda = Lammda > 5 ? 5 : Lammda;
    if (abs(this->lammda - Lammda) > 0.2)
    {
        this->lammda = Lammda;
        needUpdateBANN = true;
        if(CstrainMax != 0)
            needUpdateHANN_pos = true;
        if(CstrainMin != 0)
            needUpdateHANN_neg = true;
    }
}

bool AIDMMaterial::checkCapacity(const double& Moment, const int& eleTag)
{
    if (!isAvailabelAIDM())
        return true;
    auto limitFactor = 0.8;
    auto capacity = Moment < 0 ? this->stressSA_neg * this->stressFactor_neg : 
         this->stressSA_pos * this->stressFactor_pos;
     if (abs(Moment) < capacity * limitFactor)
         return true;
     /*if (Moment > 0)
         this->stressSA_pos = abs(Moment) * ( 1 / limitFactor);
     else this->stressSA_neg = abs(Moment)* (1 / limitFactor);*/
     opserr << "Warinning: the capacity of AIDMBeamColumn " << eleTag << " with AIDMMaterial " << this->getTag() <<
         " is too weak: convert to elastic automatically."<< endln;
     this->isToElastic = true;
     K = this->CStress / this->CStrain;
     return false;
}

double AIDMMaterial::getStressOnBackbone(const double& drift)
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
        m = coe_AM * pow(drift, 1) + coe_BM;
        n = coe_AN * pow(drift, 1) + coe_BN;
    }
    auto capacity = drift > 0 ? stressFactor_pos * stressSA_pos : stressFactor_neg * stressSA_neg;
    auto x = abs(drift / dc);
    auto y = (m * x) / (1 + (m - (n / (n - 1))) * x + pow(x, n) / (n - 1));
    return drift > 0? y * capacity : -y * capacity;
}

void AIDMMaterial::setInitialK(double iniK, bool isToElastic)
{
    //Initial Boundary
    this->initialK = iniK;
    mnFactor = 1;
    //Elastic
    this->isToElastic = isToElastic;
    if (this->isToElastic)
    {
        K = this->initialK;
        return;
    }
    else if (!this->ensureIniK) return;
    auto initial_strain_pos = backbone_inidStrainFactor * strainC_pos;
    auto initial_strain_neg = backbone_inidStrainFactor * -strainC_neg;
    auto factorpos = 0.0;
    auto factorneg = 0.0;
    while (factorpos < 0.8 && factorneg < 0.8)
    {
        mnFactor += 1;
        auto initial_strain_K_pos = this->getStressOnBackbone(initial_strain_pos) / initial_strain_pos;
        auto initial_strain_K_neg = this->getStressOnBackbone(initial_strain_neg) / initial_strain_neg;
        factorpos = initial_strain_K_pos / this->initialK;
        factorneg = initial_strain_K_neg / this->initialK;
        if (mnFactor == 5)
            return;
    }
}

bool AIDMMaterial::isAvailabelAIDM() const
{
    if (this->lammdaS == 0 || this->lammdaSV == 0 || this->lammdaT_pos == 0)
        return false;
    return true;
}
