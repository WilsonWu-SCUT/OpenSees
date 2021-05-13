#include "PMMSection.h"
#include "SectionAnalysisExport.h"
#include "AutoMesh.h"

using namespace AutoMesh;

#include <elementAPI.h>
#include "AIDMMaterial.h"

void* OPS_PMMSectionRectBeam()
{
	if (OPS_GetNumRemainingInputArgs() < 7) {
		opserr << "insufficient arguments for PMMSectionRectBeam\n";
		return 0;
	}

	// get tag
	int tag, AIDMtag_3, AIDMtag_2;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) < 0) return 0;
	if (OPS_GetIntInput(&numData, &AIDMtag_3) < 0) return 0;
	if (OPS_GetIntInput(&numData, &AIDMtag_2) < 0) return 0;

	//参数容器
	std::vector<double> dimension_vec;
	std::vector<int> As_vec = {0};
	int fcu;
	double fy;

	// get diemention
	numData = 2;
	double double_data[2];
	if (OPS_GetDoubleInput(&numData, &double_data[0]) < 0) return 0;
	for (auto info : double_data)  dimension_vec.push_back(info);

	// get As
	numData = 2;
	int int_data[2];
	if (OPS_GetIntInput(&numData, &int_data[0]) < 0) return 0;
	for (auto info : int_data)  As_vec.push_back(info);

	numData = 1;
	if (OPS_GetIntInput(&numData, &fcu) < 0) return 0;
	if (OPS_GetDoubleInput(&numData, &fy) < 0) return 0;


	return new PMMSection(tag, AIDMtag_3, AIDMtag_2, 1, dimension_vec, As_vec, true, fcu, fy, 345);
}

PMMSection::PMMSection()
	:SectionForceDeformation(0, SEC_TAG_PMMSection), 
	section_sp_(), FRAMSection_sp_(),
	AIDM_3_ptr(), AIDM_2_ptr(),
	AIDM_3_Tag_(-9), AIDM_2_Tag_(-9)
{

}

PMMSection::PMMSection(const int& tag, const int& AIDMTag_3, const int& AIDMTag_2,
	const int& section_type,
	const std::vector<double> dimension_vec,
	const std::vector<int> As_vec, bool is_beam,
	const int& fcu, const double& bar_fy, const double& steel_fy)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection), AIDM_3_Tag_(AIDMTag_3), AIDM_2_Tag_(AIDMTag_2)
{
	//构造截面
	this->FRAMSection_sp_.reset(new AutoMesh::FRAMSection(section_type,
		AutoMesh::SectionType::RC));
	//设定截面基本参数（保护层默认20）
	auto isSucess = this->FRAMSection_sp_->set_dimension(dimension_vec, 20);
	//失败直接返回
	if (!isSucess)
	{
		opserr << "PMMSection::PMMSection Error to set_dimension\n";
		return;
	}
	//设定配筋
	isSucess = this->FRAMSection_sp_->set_As(As_vec, is_beam);
	//失败直接返回
	if (!isSucess)
	{
		opserr << "PMMSection::PMMSection Error to set_As\n";
		return;
	}
	//剖分截面
	auto mesh_ptr = this->FRAMSection_sp_->get_mesh_sp(50, AutoMesh::MeshMethod::LoopingTri, true);
	isSucess = mesh_ptr->mesh();	
	//失败直接返回
	if (!isSucess)
	{
		opserr << "PMMSection::PMMSection Error to mesh\n";
		return;
	}
	//设定截面参数
	this->section_sp_.reset(new SectionAnalysis::Section(mesh_ptr, fcu, steel_fy));
	//获得配筋及配筋面积
	auto as_pos_vec = this->FRAMSection_sp_->get_as_pos_vec();
	auto as_A_vec = this->FRAMSection_sp_->get_as_A_vec();
	//遍历向量
	for (int i = 0; i < std::min(as_pos_vec.size(), as_A_vec.size()); i++)
		this->section_sp_->add_reinforced_bar(as_pos_vec[i]->get_x(), 
			as_pos_vec[i]->get_y(), as_A_vec[i], bar_fy);
	//截面分析
	this->section_sp_->analysis(4);
	//获得AIDM指针
	if (this->AIDM_3_Tag_ < 0) this->AIDM_3_ptr = new AIDMMaterial();
	else this->AIDM_3_ptr = dynamic_cast<AIDMMaterial*>(OPS_getUniaxialMaterial(this->AIDM_3_Tag_)->getCopy());
	if (this->AIDM_2_Tag_ < 0) this->AIDM_2_ptr = new AIDMMaterial();
	else this->AIDM_2_ptr = dynamic_cast<AIDMMaterial*>(OPS_getUniaxialMaterial(this->AIDM_2_Tag_)->getCopy());
	//判断是否读取成功
	if (this->AIDM_3_ptr == 0)
	{
		opserr << "WARNING no AIDMmaterial " << this->AIDM_3_Tag_ <<
			"exitsts - section PMMSection\n";
		this->AIDM_3_Tag_ = -2;
	}
	if (this->AIDM_2_ptr == 0)
	{
		opserr << "WARNING no AIDMmaterial " << this->AIDM_2_Tag_ <<
			"exitsts - section PMMSection\n";
		this->AIDM_2_Tag_ = -2;
	}
}

PMMSection::PMMSection(const int& tag, std::shared_ptr<SectionAnalysis::Section> section_sp)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection)
{
	this->section_sp_ = section_sp;
}

PMMSection::~PMMSection()
{

}

int PMMSection::get_moment(const double& My, const double& Mz, const double& axial_load, bool isI)
{
	auto axial_load_kn = axial_load / 1E3;
	auto moment = this->section_sp_->get_moment(My, Mz, axial_load_kn, isI);
	return moment * 1E6;
}

void PMMSection::setLammda(const Vector& force_vec, bool is_I)
{
	this->AIDM_3_ptr->setLammda(is_I?
		force_vec(4) / force_vec(2): force_vec(10) / force_vec(8));
	this->AIDM_2_ptr->setLammda(is_I ?
		force_vec(5) / force_vec(1) : force_vec(11) / force_vec(7));
}

int PMMSection::setTrialDeformation(const Vector& deformation_vec, bool is_I)
{
	int retVal = this->AIDM_3_ptr->setTrialStrain(
		is_I ? deformation_vec(3): -deformation_vec(4));
	retVal += this->AIDM_2_ptr->setTrialStrain(
		is_I ? deformation_vec(1) : -deformation_vec(2));
	return retVal;
}

SectionForceDeformation* PMMSection::getCopy(void)
{
	auto section = new PMMSection(this->getTag(), this->section_sp_);
	section->AIDM_2_ptr = dynamic_cast<AIDMMaterial*>(this->AIDM_2_ptr->getCopy());
	section->AIDM_3_ptr = dynamic_cast<AIDMMaterial*>(this->AIDM_3_ptr->getCopy());
	return section;
}

#pragma region unavailable

const Vector& PMMSection::getSectionDeformation(void)
{
	return s;
}

const Vector& PMMSection::getStressResultant(void)
{
	return s;
}

const Matrix& PMMSection::getSectionTangent(void)
{
	return ks;
}

const Matrix& PMMSection::getInitialTangent(void)
{
	return ks;
}

int PMMSection::commitState(void)
{
	//AIDM CommitState
	int retVal = this->AIDM_2_ptr->commitState();
	retVal += this->AIDM_3_ptr->commitState();
	//is kill
	if (this->AIDM_2_ptr->iskill()) this->AIDM_2_Tag_ = -1;
	if (this->AIDM_3_ptr->iskill()) this->AIDM_3_Tag_ = -1;
	//是否成功
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
	this->AIDM_3_ptr->revertToStart();
	return 0;
}

const ID& PMMSection::getType(void)
{
	return NULL;
}

int PMMSection::getOrder(void) const
{
	return 0;
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


Vector PMMSection::s(0);

Matrix PMMSection::ks(2, 2);

#pragma endregion


