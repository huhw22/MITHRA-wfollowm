/********************************************************************************************************
* fdtd.h 实现了真正的fdtd时间解代码
********************************************************************************************************/

#ifndef FDTD_HH_
#define FDTD_HH_

#include "classes.h"
#include "solver.h"

namespace MITHRA
{

  /*该类函数用于求解时域中的时域有限差分域。*/
  class FdTd : public Solver
  {
  public:
    FdTd (Mesh& 				mesh,
	  Bunch& 				bunch,
	  Seed& 				seed,
	  std::vector<Undulator>&		undulator,
	  std::vector<ExtField>& 		extField,
	  std::vector<FreeElectronLaser>& 	FEL);

    /*将电流复位为零。*/
    void currentReset ();

    /*为字段更新更新单元格点上的电流。*/
    void currentUpdate ();

    /*在不同的处理器之间传输电流。*/
    void currentCommunicate ();

    /*为一个时间步更新字段*/
    void fieldUpdate ();

    /*从电位中计算第m个像素的场。*/
    void fieldShift ();

    /*从电位中计算第m个像素的场。*/
    void fieldEvaluate (long int m);

    /*采样字段并将其保存到给定的文件中。*/
    void fieldSample ();

    /*将字段可视化为整个域上的vtk文件，并将它们保存到给定名称的文件中。*/
    void fieldVisualizeAllDomain 	(unsigned int ivtk);

    /*将字段可视化为平面中的vtk文件，并将它们保存到具有给定名称的文件中。*/
    void fieldVisualizeInPlane 		(unsigned int ivtk);

    /*将字段可视化为垂直于x轴的平面上的vtk文件，并将它们保存到文件中
    *名字。*/
    void fieldVisualizeInPlaneXNormal 	(unsigned int ivtk);

    /*将字段可视化为垂直于y轴的平面上的vtk文件，并将它们保存到文件中
    *名字。*/
    void fieldVisualizeInPlaneYNormal 	(unsigned int ivtk);

    /*将字段可视化为垂直于z轴的平面上的vtk文件，并将它们保存到文件中
    *名字。*/
    void fieldVisualizeInPlaneZNormal 	(unsigned int ivtk);

    /*将字段的总概要文件写入给定的文件名中。*/
    void fieldProfile ();


  };
}
#endif
