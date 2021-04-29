#pragma once

namespace AutoMesh
{
	//�͸�����
	enum class SectionBasicShapeType
	{
		//���ָ�
		IShaped = 0,

		//����
		BoxShaped = 1,

		//����
		RingShaped = 2,

		//ʮ�ָ�
		XShaped = 3,

		//����
		RectangularShaped = 4,

		//Բ��
		CircleShaped = 5,
	};

	//��������
	enum class MatType
	{
		//������
		Concrete = 0,

		//�ݽ�
		ReinforceBar = 1,

		//�ֲ�
		Steel = 2,

		//�߽������
		CoverConcrete = 3,
	};

	//��������
	enum class SectionType
	{
		RC = 6,

		Steel = 5,
	};

	//�ʷַ���
	enum class MeshMethod
	{
		//�ı����ʷ�
		LoopingQuad = 0,
		SubmappingQuad = 1,
		TriToQuad = 2,
		PavingQuad = 3,
		MappingQuad = 4,
		GridQuadTri = 5,
		DelaunayTri = 6,
		//�������ʷ�
		LoopingTri = 7,
		AFMTri = 8,
		CADTri = 9,
	};
}
