#include "PMMSection.h"
#include "AutoMesh.h"
#include "SectionAnalysisExport.h"
#include <elementAPI.h>
#include "AIDMMaterial.h"

using namespace AutoMesh;



void* OPS_PMMSectionRectBeam()
{
	//�����ǲ�������
	if (OPS_GetNumRemainingInputArgs() != 9 &&
		OPS_GetNumRemainingInputArgs() != 3) {
		opserr << "insufficient arguments for PMMSectionRectBeam\n";
		return 0;
	}
	//�ж��Ƿ�δ����
	bool isPlastic = OPS_GetNumRemainingInputArgs() > 3;

	// get tag
	int tag;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) < 0) return 0;

	//��������
	std::vector<double> dimension_vec;
	std::vector<int> As_vec = {0};
	double fck;
	double fy;
	double lammdaSV_y, lammdaSV_z;

	//�ߴ�
	numData = 2;
	double double_data[2];
	if (OPS_GetDoubleInput(&numData, &double_data[0]) < 0) return 0;
	for (auto info : double_data)  dimension_vec.push_back(info);

	if (isPlastic)
	{
		//���
		numData = 2;
		int int_data[2];
		if (OPS_GetIntInput(&numData, &int_data[0]) < 0) return 0;
		for (auto info : int_data)  As_vec.push_back(info);
		//����ǿ��
		numData = 1;
		if (OPS_GetDoubleInput(&numData, &fck) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &fy) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &lammdaSV_y) < 0) return 0;
		if (OPS_GetDoubleInput(&numData, &lammdaSV_z) < 0) return 0;
		return new PMMSection(tag, 1, dimension_vec, As_vec, true,
			fck, fy, 345, AutoMesh::SectionType::RC, lammdaSV_y, lammdaSV_z);
	}
	else return new PMMSection(tag, 1, dimension_vec, AutoMesh::SectionType::RC);	
}

PMMSection::PMMSection()
	:SectionForceDeformation(0, SEC_TAG_PMMSection), 
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(0), steel_fy_(0)
{
	//��ʼ��ָ��
	this->AIDM_2_ptr = new AIDMMaterial();
	this->AIDM_3_ptr = new AIDMMaterial();
}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, const std::vector<int> As_vec,
	bool is_beam, const double& fck, const double& bar_fy, const double& steel_fy,
	AutoMesh::SectionType matType, const double& lammdaSV_y, const double& lammdaSV_z)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(fck), steel_fy_(steel_fy)
{
	//��ʼ������
	this->as_vec_ = As_vec;
	if (!this->iniSection(section_type, dimension_vec, matType, bar_fy, is_beam)) return;
	//����������ѹ����
	auto lammdaT_2_pos = 1;
	auto lammdaT_3_pos = this->FRAMSection_sp_->isBeam() ? this->as_vec_[2] / this->as_vec_[1] : 1;
	//�������
	this->section_sp_->analysis(4);
	//���������
	this->capacity_3pos_ini_ = this->section_sp_->get_moment(0, 180) * 1E6;
	this->capacity_3neg_ini_ = this->section_sp_->get_moment(0, 0) * 1E6;
	this->capacity_2pos_ini_ = this->section_sp_->get_moment(0, 270) * 1E6;
	this->capacity_2neg_ini_ = this->section_sp_->get_moment(0, 90) * 1E6;
	//�ݽ��������ֵ
	auto lammdaS = this->section_sp_->A(AutoMesh::MatType::ReinforceBar) / this->section_sp_->A()
		* bar_fy / this->fck_;
	//AIDMTAG
	int tag2 = tag * 100 + 2;
	int tag3 = tag * 100 + 3;
	//��������ָ��
	this->AIDM_2_ptr = lammdaSV_y <= 0 ? new AIDMMaterial() : new AIDMMaterial(tag2, lammdaSV_y, lammdaS, lammdaT_2_pos);
	this->AIDM_3_ptr = lammdaSV_z <= 0 ? new AIDMMaterial() : new AIDMMaterial(tag3, lammdaSV_z, lammdaS, lammdaT_3_pos);
	//�趨����
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	this->AIDM_3_ptr->updateSkeletonParams();
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	this->AIDM_2_ptr->updateSkeletonParams();
	
}

PMMSection::PMMSection(const int& tag, const int& section_type,
	const std::vector<double> dimension_vec, AutoMesh::SectionType matType)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(30), steel_fy_(345)
{
	//��ʼ������
	if (!this->iniSection(section_type, dimension_vec, matType, 400, true)) return;
	//��ʼ��ָ��
	this->AIDM_2_ptr = new AIDMMaterial();
	this->AIDM_3_ptr = new AIDMMaterial();
}

PMMSection::PMMSection(const int& tag)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection),
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	sectionHeight_3_(0), sectionHeight_2_(0),
	capacity_3pos_ini_(0), capacity_3neg_ini_(0),
	capacity_2pos_ini_(0), capacity_2neg_ini_(0),
	fck_(0), steel_fy_(0)
{
	//��ʼ��ָ��
	this->AIDM_2_ptr = new AIDMMaterial();
	this->AIDM_3_ptr = new AIDMMaterial();
}

PMMSection::~PMMSection()
{
	if (this->AIDM_3_ptr) delete this->AIDM_3_ptr;
	if (this->AIDM_2_ptr) delete this->AIDM_2_ptr;
}

bool PMMSection::iniSection(const int& section_type, const std::vector<double> dimension_vec,
	AutoMesh::SectionType& matType, const double& bar_fy,bool isbeam)
{
	//�������
	this->FRAMSection_sp_.reset(new AutoMesh::FRAMSection(section_type, matType));
	//�趨�������������������Ĭ��20��
	if (!this->FRAMSection_sp_->set_dimension(dimension_vec, 20))
	{
		opserr << "PMMSection::PMMSection Error to set_dimension\n";
		return false;
	}
	//�������߶�
	this->sectionHeight_2_ = this->FRAMSection_sp_->get_B();
	this->sectionHeight_3_ = this->FRAMSection_sp_->get_H();
	//�ж��Ƿ�Ϊ��Ч���
	if (this->as_vec_.size() != 0)
	{
		if (!this->FRAMSection_sp_->set_As(this->as_vec_, isbeam))
		{
			opserr << "PMMSection::PMMSection Error to set_As\n";
			return false;
		}
		
	}
	//�ʷֽ���
	auto mesh_ptr = this->FRAMSection_sp_->get_mesh_sp(50, AutoMesh::MeshMethod::LoopingTri, false);
	if (!mesh_ptr->mesh())
	{
		opserr << "PMMSection::PMMSection Error to mesh\n";
		return false;
	}
	//�趨�������
	this->section_sp_.reset(new SectionAnalysis::Section(mesh_ptr, this->fck_, 0.1 * this->fck_, this->steel_fy_));
	//�����������
	if (this->as_vec_.size() != 0)
	{
		auto as_pos_vec = this->FRAMSection_sp_->get_as_pos_vec();
		auto as_A_vec = this->FRAMSection_sp_->get_as_A_vec();
		//��������
		for (int i = 0; i < std::min(as_pos_vec.size(), as_A_vec.size()); i++)
			this->section_sp_->add_reinforced_bar(as_pos_vec[i]->get_x(),
				as_pos_vec[i]->get_y(), as_A_vec[i], bar_fy);
	}
	return true;
}

void PMMSection::setInitialK(const double& L, const int& factor)
{
	/*�趨��ʼ�ն�*/
	auto initial_K2 = this->E() * this->Iz() * factor / L;
	auto initial_K3 = this->E() * this->Iy() * factor / L;
	//������ʼ�ն�
	this->AIDM_2_ptr->setInitialK(initial_K2);
	this->AIDM_3_ptr->setInitialK(initial_K3);
}

void PMMSection::setLammda(const Vector& force_vec, bool is_I)
{
	//��ֹ�Ƿ�����
	if (force_vec(2) == 0 || force_vec(8) == 0 || force_vec(1) == 0 || force_vec(7) == 0 || 
		this->sectionHeight_3_ == 0 || this->sectionHeight_2_ == 0)
		return;
	//�������
	auto shearSpan_3 = is_I ?
		force_vec(4) / force_vec(2) : force_vec(10) / force_vec(8);
	auto shearSpan_2 = is_I ?
		force_vec(5) / force_vec(1) : force_vec(11) / force_vec(7);
	//��������
	auto shearspanRatio_3 = std::abs(shearSpan_3) / this->sectionHeight_3_;
	auto shearspanRatio_2 = std::abs(shearSpan_2) / this->sectionHeight_2_;
	//�趨�����
	this->AIDM_2_ptr->setLammda(shearspanRatio_2);
	this->AIDM_3_ptr->setLammda(shearspanRatio_3);
}

void PMMSection::setCapacity(const Vector& force_vec, bool is_I)
{
	//AIDM������
	if (this->as_vec_.size() == 0) return;
	//��������
	auto my = (is_I ? force_vec(4) : force_vec(10) * -1) / 1E6;
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1) / 1E6;
	auto P = (is_I ? -1 * force_vec(0) : force_vec(6)) / 1E3;
	//��ʼ��������
	double myca_pos, myca_neg, mzca_pos, mzca_neg;
	//��������� ��ƫѹ
	if (!this->consider_double_bending)
	{
		myca_pos = this->section_sp_->get_moment(P, 180);
		myca_neg = this->section_sp_->get_moment(P, 0);
		mzca_pos = this->section_sp_->get_moment(P, 270);
		mzca_neg = this->section_sp_->get_moment(P, 90);
	}
	//����˫ƫѹ
	else this->section_sp_->get_moment(my, mz, P, myca_pos, myca_neg, mzca_pos, mzca_neg);
	//���µ�λ
	myca_pos *= 1E6; myca_neg *= 1E6; mzca_pos *= 1E6; mzca_neg *= 1E6;
	//��ֹ��������С
	myca_pos = myca_pos < this->capacity_3pos_ini_* this->min_capacity_factor ? 
		this->capacity_3pos_ini_ * this->min_capacity_factor : myca_pos;

	myca_neg = myca_neg < this->capacity_3neg_ini_* this->min_capacity_factor ?
		this->capacity_3neg_ini_ * this->min_capacity_factor : myca_neg;

	mzca_pos = mzca_pos < this->capacity_2pos_ini_* this->min_capacity_factor ?
		this->capacity_2pos_ini_ * this->min_capacity_factor : mzca_pos;

	mzca_neg = mzca_neg < this->capacity_2neg_ini_* this->min_capacity_factor ?
		this->capacity_2neg_ini_ * this->min_capacity_factor : mzca_neg;
	//������ѹϵ��
	auto Ac = this->section_sp_->A(AutoMesh::MatType::Concrete) +
		this->section_sp_->A(AutoMesh::MatType::CoverConcrete);
	auto fA = this->steel_fy_ * this->section_sp_->A(AutoMesh::MatType::Steel) + Ac * 
		(P <= 0 ? this->fck_: 0.1 * this->fck_);
	auto axialRatio = fA == 0 ? 0 : P / fA * 1E3;
	//�趨������
	this->AIDM_3_ptr->setCapcacity(myca_pos, myca_neg, axialRatio);
	this->AIDM_2_ptr->setCapcacity(mzca_pos, mzca_neg, axialRatio);
}

bool PMMSection::checkCapacity(const Vector& force_vec, const int& eleTag, bool is_I)
{
	//��������
	auto my = (is_I ? force_vec(4) : force_vec(10) * -1);
	auto mz = (is_I ? force_vec(5) : force_vec(11) * -1);
	//�ж������Ƿ���������
	bool is2Bool = this->AIDM_2_ptr->checkCapacity(mz);
	bool is3Bool = this->AIDM_3_ptr->checkCapacity(my);
	if (!is2Bool || !is3Bool)
		opserr << "Warinning: the capacity of AIDMBeamColumn " << eleTag << " with PMMSection " << this->getTag() <<
		" is too weak: convert to elastic automatically." << endln;
	return (is2Bool && is3Bool);
}

int PMMSection::setTrialDeformation(const Vector& deformation_vec, bool is_I)
{
	int retVal = this->AIDM_3_ptr->setTrialStrain(
		is_I ? deformation_vec(3): -deformation_vec(4));
	retVal += this->AIDM_2_ptr->setTrialStrain(
		is_I ? deformation_vec(1) : -deformation_vec(2));
	return retVal;
}

std::vector<std::string> PMMSection::getResponseStrVec(bool is_I)
{
	std::string descp = is_I ? "I" : "J";
	return
	{
		"strain2" + descp,
		"stress2" + descp,
		"lammda2" + descp,
		"strain3" + descp,
		"stress3" + descp,
		"lammda3" + descp,
	};
}

std::vector<double> PMMSection::getResponseVec()
{
	return
	{
		this->AIDM_2_ptr->getStrain(),
		this->AIDM_2_ptr->getStress(),
		this->AIDM_2_ptr->getLammda(),
		this->AIDM_3_ptr->getStrain(),
		this->AIDM_3_ptr->getStress(),
		this->AIDM_3_ptr->getLammda(),
	};
}

SectionForceDeformation* PMMSection::getCopy(void)
{
	//��ʼ��
	auto section = new PMMSection(this->getTag());
	section->section_sp_ = this->section_sp_;
	section->FRAMSection_sp_ = this->FRAMSection_sp_;
	//��������
	section->sectionHeight_3_ = this->sectionHeight_3_;
	section->sectionHeight_2_ = this->sectionHeight_2_;
	//��ʼ������
	section->capacity_3pos_ini_ = this->capacity_3pos_ini_;
	section->capacity_3neg_ini_ = this->capacity_3neg_ini_;
	section->capacity_2pos_ini_ = this->capacity_2pos_ini_;
	section->capacity_2neg_ini_ = this->capacity_2neg_ini_;
	//���
	section->as_vec_ = this->as_vec_;
	section->fck_ = this->fck_;
	section->steel_fy_ = this->steel_fy_;

	//�ع�AIDMָ��
	section->AIDM_2_ptr = (AIDMMaterial*)this->AIDM_2_ptr->getCopy();
	section->AIDM_3_ptr = (AIDMMaterial*)this->AIDM_3_ptr->getCopy();
	section->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	section->AIDM_3_ptr->updateSkeletonParams();
	section->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	section->AIDM_2_ptr->updateSkeletonParams();

	return section;
}

#pragma region Section basic information

double PMMSection::A(void) const
{
	return this->section_sp_->A();
}

double PMMSection::Iy(void) const
{
	return this->section_sp_->Iy();
}

double PMMSection::Iz(void) const
{
	return this->section_sp_->Iz();
}

double PMMSection::E(void) const
{
	return this->section_sp_->E();
}

double PMMSection::Jx(void) const
{
	return this->Iy() + this->Iz();
}

double PMMSection::G(void) const
{
	return  0.4 * this->E();
}

#pragma endregion

#pragma region unavailable

const Vector& PMMSection::getSectionDeformation(void)
{
	//��ʼ��
	s.Zero();
	//��ó�ʼ�ն�
	s(0) = this->AIDM_2_ptr->getStrain();
	s(1) = this->AIDM_3_ptr->getStrain();
	return s;
}

const Vector& PMMSection::getStressResultant(void)
{
	//��ʼ��
	s.Zero();
	//��ó�ʼ�ն�
	s(0) = this->AIDM_2_ptr->getStress();
	s(1) = this->AIDM_3_ptr->getStress();
	return s;
}

const Matrix& PMMSection::getSectionTangent(void)
{
	//��ʼ��
	ks.Zero();
	ks(0, 0) = this->AIDM_2_ptr->getTangent();
	ks(1, 1) = this->AIDM_3_ptr->getTangent();	
	return ks;
}

const Matrix& PMMSection::getInitialTangent(void)
{
	//��ʼ��
	ks.Zero();
	//��ó�ʼ�ն�
	ks(0, 0) = this->AIDM_2_ptr->getInitialTangent();
	ks(1, 1) = this->AIDM_3_ptr->getInitialTangent();
	return ks;
}

int PMMSection::commitState(void)
{
	//AIDM CommitState
	int retVal = this->AIDM_2_ptr->commitState();
	retVal += this->AIDM_3_ptr->commitState();
	//�Ƿ�ɱ����Ԫ
	if (this->AIDM_2_ptr->isKill() || this->AIDM_3_ptr->isKill())
	{
		this->AIDM_2_ptr->Kill();
		this->AIDM_3_ptr->Kill();
	}
	//�Ƿ�ɹ�
	return retVal;
}

int PMMSection::revertToLastCommit(void)
{
	int retVal = this->AIDM_2_ptr->revertToLastCommit();
	retVal += this->AIDM_3_ptr->revertToLastCommit();
	return retVal;
}

int PMMSection::revertToStart(void)
{
	this->AIDM_2_ptr->revertToStart();
	this->AIDM_3_ptr->setCapcacity(this->capacity_3pos_ini_, this->capacity_3neg_ini_, 0);
	this->AIDM_3_ptr->revertToStart();
	this->AIDM_2_ptr->setCapcacity(this->capacity_2pos_ini_, this->capacity_2neg_ini_, 0);
	return 0;
}

const ID& PMMSection::getType(void)
{
	return NULL;
}

int PMMSection::getOrder(void) const
{
	return 2;
}

int PMMSection::sendSelf(int commitTag, Channel& theChannel)
{
	opserr << "PMMSection::sendSelf() - failed to send data\n";
	return -1;
}

int PMMSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	opserr << "PMMSection::recvSelf() - failed to receive data\n";
	return -1;
}

void PMMSection::Print(OPS_Stream& s, int flag)
{
	opserr << "PMMSection::Print() - failed to print\n";
}

Vector PMMSection::s(2);

Matrix PMMSection::ks(2, 2);

#pragma endregion