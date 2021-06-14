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
		//�жϱ������ֵ���η���
	if(this->ARK == 0)  this->CALR = 0;
	//��һ�η�����ؿ��ܳ���
	else if (TStrain * dc <= 0)
	{
		//������
	}
	//δ������ֵ��������Ӧ����
	else if (std::abs(TStrain) < std::abs(dc))
	{
		auto oriented_strain = dc;
		auto oriented_ALR = this->GetALRMax(lammda, dc);
		auto slope = (oriented_ALR - this->CALR) / (oriented_strain - this->Cstrain);
		this->CALR = slope * (TStrain - this->Cstrain) + this->CALR;
	}
	//������ֵ��������Ӧ����
	else
	{
		auto ALRMax = this->GetALRMax(lammda, dc);
		this->CALR = ALRMax * std::abs(CStressFactor);
	}
	//���±���
	this->Cstrain = TStrain;
	//����Լ������ϵ��
	return this->CALR * this->ReductFactor;
}

double DALRMaterial::GetUnloadALR(const double& TStrain)
{
	//�жϱ�����Cstrain����
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
	//�Ƿ�������ε�ͬ��
	if (TStrain * strainMax <= 0 || this->ARK == 0)
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
		//�������Զ�
		//auto calr = std::abs(TStrain) * this->GetInitialK(lammda);
		//�ж��Ƿ񳬹�
		//if (this->CALR < calr) this->CALR = calr;
		this->Cstrain = TStrain;
	}
	//������ֵ
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