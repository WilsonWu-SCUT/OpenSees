#include <DALRMaterial.h>
#include <algorithm>
#include <cmath>

DALRMaterial::DALRMaterial()
	:Cstrain(0), CALR(0), ARK(0)
{

}

DALRMaterial::~DALRMaterial()
{

}

double DALRMaterial::GetALROnBackbone(const double& TStrain, const double& dc,
	const double& lammda, const double& CStressFactor)
{
	//���±���
	this->Cstrain = TStrain;
	//�жϱ������ֵ���η���
	if (TStrain * dc <= 0 || this->ARK < 0)
	{
		this->CALR = 0;
	}
	//δ������ֵ��������Ӧ����
	else if (std::abs(TStrain) < std::abs(dc))
	{
		this->CALR = std::abs(TStrain) * this->GetInitialK(lammda);
		
	}
	//������ֵ��������Ӧ����
	else
	{
		auto ALRMax = this->GetALRMax(lammda, dc);
		this->CALR = ALRMax * std::abs(CStressFactor);
	}
	return this->CALR;
}

double DALRMaterial::GetUnloadALR(const double& TStrain)
{
	//�жϱ�����Cstrain����
	if (TStrain * this->Cstrain <= 0 || this->ARK < 0)
	{
		this->Cstrain = TStrain;
		this->CALR = 0;
	}
	else
	{
		this->CALR = this->CALR / this->Cstrain * TStrain;
		this->Cstrain = TStrain;
	}
	return this->CALR;
}

double DALRMaterial::GetReloadALR(const double& TStrain, const double& dc, const double& strainMax,
	const double& lammda, const double& CStressFactor)
{
	//�Ƿ�������ε�ͬ��
	if (TStrain * strainMax <= 0 || this->ARK < 0)
	{
		this->Cstrain = TStrain;
	}
	//�Ƿ񳬹���ֵ����������
	else if (std::abs(dc) > std::abs(strainMax))
	{
		auto oriented_strain = dc;
		auto oriented_ALR = this->GetALRMax(lammda, dc);
		auto slope = (oriented_ALR - this->CALR) / (oriented_strain - this->Cstrain);
		this->CALR = slope * (TStrain - this->Cstrain) + this->CALR;
		this->Cstrain = TStrain;
	}
	else
	{
		auto oriented_strain = strainMax;
		auto oriented_ALR = this->GetALRMax(lammda, dc) * std::abs(CStressFactor);
		auto slope = (oriented_ALR - this->CALR) / (oriented_strain - this->Cstrain);
		this->CALR = slope * (TStrain - this->Cstrain) + this->CALR;
		this->Cstrain = TStrain;
	}
	return this->CALR;
}

double DALRMaterial::GetInitialK(const double& lammda)
{
	auto a = -lammda * 2E-5 - 1E-5;
	auto b = lammda >= 1 ? 0.023 * std::log(lammda) + 0.017: 0.017;
	auto xmax = this->ARK <= 0 ? -0.5 * b / a : this->ARK;
	auto k = a * std::pow(xmax, 2) + b * xmax;
	return k;
}

double DALRMaterial::GetALRMax(const double& lammda, const double& dc)
{
	return this->GetInitialK(lammda) * std::abs(dc);
}