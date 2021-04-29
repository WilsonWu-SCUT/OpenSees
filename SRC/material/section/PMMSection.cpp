#include "PMMSection.h"
#include "AutoMesh.h"

using namespace AutoMesh;

#include <elementAPI.h>



void* OPS_PMMSectionRectBeam()
{
	if (OPS_GetNumRemainingInputArgs() < 7) {
		opserr << "insufficient arguments for PMMSectionRectBeam\n";
		return 0;
	}

	// get tag
	int tag;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) < 0) return 0;

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


	return new PMMSection(tag, 1, dimension_vec, As_vec, true, fcu, fy, 345);
}

PMMSection::PMMSection()
	:SectionForceDeformation(0, SEC_TAG_PMMSection)
{

}

PMMSection::PMMSection(const int& tag,
	const int& section_type,
	const std::vector<double> dimension_vec,
	const std::vector<int> As_vec, bool is_beam,
	const int& fcu, const double& bar_fy, const double& steel_fy)
	: SectionForceDeformation(tag, SEC_TAG_PMMSection)
{
	//构造截面
	std::shared_ptr<AutoMesh::FRAMSection> section(new AutoMesh::FRAMSection(section_type, 
		AutoMesh::SectionType::RC));
	//设定截面基本参数（保护层默认20）
	auto isSucess = section->set_dimension(dimension_vec, 20);
	//失败直接返回
	if (!isSucess)
	{
		opserr << "PMMSection::PMMSection Error to set_dimension\n";
		return;
	}
	//设定配筋
	isSucess = section->set_As(As_vec, is_beam);
	//失败直接返回
	if (!isSucess)
	{
		opserr << "PMMSection::PMMSection Error to set_As\n";
		return;
	}
	//剖分截面
	auto mesh_ptr = section->get_mesh_sp(50, AutoMesh::MeshMethod::LoopingTri, true);
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
	auto as_pos_vec = section->get_as_pos_vec();
	auto as_A_vec = section->get_as_A_vec();
	//遍历向量
	for (int i = 0; i < std::min(as_pos_vec.size(), as_A_vec.size()); i++)
		this->section_sp_->add_reinforced_bar(as_pos_vec[i]->get_x(), 
			as_pos_vec[i]->get_y(), as_A_vec[i], bar_fy);
	//截面分析
	this->section_sp_->analysis(4);
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
	return 0;
}

int PMMSection::revertToLastCommit(void)
{
	return 0;
}
int PMMSection::revertToStart(void)
{
	return 0;
}

SectionForceDeformation* PMMSection::getCopy(void)
{
	auto section = new PMMSection(this->getTag(), this->section_sp_);
	return section;
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
	return -1;
}

int PMMSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	return -1;
}

void PMMSection::Print(OPS_Stream& s, int flag)
{

}


Vector PMMSection::s(0);
Matrix PMMSection::ks(0, 0);

#pragma endregion


