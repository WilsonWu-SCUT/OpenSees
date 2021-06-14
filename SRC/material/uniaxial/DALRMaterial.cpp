#include <DALRMaterial.h>
#include <algorithm>
#include <cmath>

DALRMaterial::DALRMaterial()
	:Cstrain(0), CALR(0), ARK(0), ReductFactor(1.0)
{

}

DALRMaterial::~DALRMaterial()
{

}

double DALRMaterial::GetALROnBackbone(const double& TStrain, const double& dc,
	const double& lammda, const double& CStressFactor)
{
		//判断变形与峰值变形反号
	if(this->ARK == 0)  this->CALR = 0;
	//第一次反向加载可能出现
	else if (TStrain * dc <= 0)
	{
		//不做事
	}
	//未超过峰值承载力对应变形
	else if (std::abs(TStrain) < std::abs(dc))
	{
		auto oriented_strain = dc;
		auto oriented_ALR = this->GetALRMax(lammda, dc);
		auto slope = (oriented_ALR - this->CALR) / (oriented_strain - this->Cstrain);
		this->CALR = slope * (TStrain - this->Cstrain) + this->CALR;
	}
	//超过峰值承载力对应变形
	else
	{
		auto ALRMax = this->GetALRMax(lammda, dc);
		this->CALR = ALRMax * std::abs(CStressFactor);
	}
	//更新变形
	this->Cstrain = TStrain;
	//返回约束周亚系数
	return this->CALR * this->ReductFactor;
}

double DALRMaterial::GetUnloadALR(const double& TStrain)
{
	//判断变形与Cstrain反号
	if (TStrain * this->Cstrain <= 0 || this->ARK == 0)
	{
		this->Cstrain = TStrain;
		this->CALR = 0;
	}
	else
	{
		this->CALR = this->CALR / this->Cstrain * TStrain;
		this->Cstrain = TStrain;
	}
	return this->CALR * this->ReductFactor;
}

double DALRMaterial::GetReloadALR(const double& TStrain, const double& dc, const double& strainMax,
	const double& lammda, const double& CStressFactor)
{
	//是否和最大变形点同向
	if (TStrain * strainMax <= 0 || this->ARK == 0)
	{
		this->Cstrain = TStrain;
	}
	//是否超过峰值承载力变形
	else if (std::abs(dc) > std::abs(strainMax))
	{
		auto oriented_strain = dc;
		auto oriented_ALR = this->GetALRMax(lammda, dc);
		auto slope = (oriented_ALR - this->CALR) / (oriented_strain - this->Cstrain);
		this->CALR = slope * (TStrain - this->Cstrain) + this->CALR;
		//计算线性段
		//auto calr = std::abs(TStrain) * this->GetInitialK(lammda);
		//判断是否超过
		//if (this->CALR < calr) this->CALR = calr;
		this->Cstrain = TStrain;
	}
	//超过峰值
	else
	{
		auto oriented_strain = strainMax;
		auto oriented_ALR = this->GetALRMax(lammda, dc) * std::abs(CStressFactor);
		auto slope = (oriented_ALR - this->CALR) / (oriented_strain - this->Cstrain);
		this->CALR = slope * (TStrain - this->Cstrain) + this->CALR;
		this->Cstrain = TStrain;
	}
	return this->CALR * this->ReductFactor;
}

double DALRMaterial::GetInitialK(const double& lammda)
{
	auto a = -lammda * 2E-5 - 1E-5;
	auto b = lammda >= 1 ? 0.023 * std::log(lammda) + 0.017: 0.017;
	auto xmax = -0.5 * b / a;
	auto x = this->ARK <= 0 || this->ARK > xmax ? xmax : this->ARK;
	auto k = a * std::pow(x, 2) + b * x;
	return k;
}

double DALRMaterial::GetALRMax(const double& lammda, const double& dc)
{
	return this->GetInitialK(lammda) * std::abs(dc);
}