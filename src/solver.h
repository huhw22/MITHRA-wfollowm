/********************************************************************************************************
* solve.h: darius代码的求解器类的实现
********************************************************************************************************/

#ifndef SOLVER_H_
#define SOLVER_H_

#include <iomanip>
#include <list>
#include <vector>

#include "classes.h"
#include "database.h"
#include "fieldvector.h"
#include "boostframe.h"

namespace MITHRA
{

  /*该类函数用于求解时域中的时域有限差分域。*/
  class Solver
  {

  public:

    Solver( Mesh& 				mesh,
	    Bunch& 				bunch,
	    Seed& 				seed,
	    std::vector<Undulator>&		undulator,
	    std::vector<ExtField>& 		extField,
	    std::vector<FreeElectronLaser>& 	FEL);

    /*使用解析后的数据，设置仿真所需的参数。*/
    void 		setSimulationParameters 	();

    /*将网格推进电子静止框架。*/
    void 		lorentzBoostMesh		();

    /*推动粒子进入电子静止框架。*/
    void 		lorentzBoostBunch 		();

    /*根据粒子的纵向坐标，在它们各自的处理器中分布粒子。*/
    void 		distributeParticles 		(std::list<Charge>& chargeVector);

    /*回收颗粒去除不再属于处理器的颗粒。*/
    void		recycleParticles		();

    /*得到从文件中读取的一堆的平均伽马和平均方向。*/
    void 		computeFileGamma 		(BunchInitialize & bunchInit);

    /*为字段值和坐标初始化矩阵。*/
    void 		initialize			();

    /*初始化问题的时间和空间网格。*/
    void 		initializeMesh			();

    /*初始化数据以更新fdtd算法中的字段。*/
    void 		initializeField			();

    /*初始化采样种子或计算域中的总字段所需的数据。*/
    void 		initializeSeedSampling		();

    /*初始化可视化和保存字段所需的数据。*/
    void 		initializeSeedVTK		();

    /*初始化分析和保存字段所需的数据。*/
    void 		initializeSeedProfile		();

    /*初始化更新堆所需的数据。*/
    void 		initializeBunchUpdate		();

    /*初始化包含用户给出的束的电荷向量。*/
    void 		initializeBunch			();

    /*用于求解时域域的函数。*/
    void 		solve					();    

    /*为一个时间步更新字段。*/
    void 		bunchUpdate			();

    /*采样这些数据并将其保存到给定的文件中。*/
    void 		bunchSample			();

    /*将这些文件可视化为vtk文件，并将它们保存到给定名称的文件中。*/
    void 		bunchVisualize			();

    /*将字段的总概要文件写入给定的文件名中。*/
    void 		bunchProfile			();

    /*计算出波动器的磁场，并将其与种子磁场相加。*/
    void 		undulatorField 			(UpdateBunchParallel& ubp, FieldVector<Double>& r);

    /*计算外场的场并将其添加到种子的场中。*/
    void 		externalField			(UpdateBunchParallel& ubp, FieldVector<Double>& r);

    /*初始化采样所需的数据，并在给定位置保存辐射功率。*/
    void 		initializePowerSample		();

    /*在给定位置取样辐射功率并保存到文件中。*/
    void 		powerSample			();

    /*初始化在给定位置显示辐射功率所需的数据。*/
    void 		initializePowerVisualize	();

    /*可视化给定位置的辐射功率并将其保存到文件中。*/
    void 		powerVisualize			();

    /*初始化采样所需的数据，并在给定位置保存辐射能量。*/
    void 		initializeEnergySample		();

    /*在给定位置采集辐射能量并保存到文件中。*/
    void 		energySample			();
    
    /*初始化在给定lab系探测平面探测的数据*/
    void initializeDetector();

    /*把探测结果保存到文件中*/
    void detectorSample();

    /*初始化在给定位置存储击中屏幕的粒子所需的数据。*/
    void 		initializeScreenProfile		();

    /*存储粒子在给定位置撞击屏幕的束状轮廓，并将其保存到
    *文件。*/
    void 		screenProfile			();

    /*完成字段计算。*/
    void 		finalize			();

    /*定义一个布尔函数来比较波动量的起始值。*/
    static bool 	undulatorCompare 		(Undulator i, Undulator j);

    /*定义线性插值函数。*/
    Double	 	interp				(Double x0, Double x1, Double y0, Double y1, Double x);

    /*定义检查粒子是否属于处理器的函数。*/
    bool 		particleInProcessor		(const Double& z );

    /*根据给定的时移，移动束和波动器的时间。*/
    void    		shiftBackInTime			();

    /*从网格索引值返回实验室帧中的实际坐标。*/
    FieldVector<Double> rc				(const long int& i);

    /*将电流复位为零。*/
    virtual void currentReset () = 0;

    /*为字段更新更新单元格点上的电流。*/
    virtual void currentUpdate () = 0;

    /*在不同的处理器之间传输电流。*/
    virtual void currentCommunicate () = 0;

    /*为一个时间步更新字段*/
    virtual void fieldUpdate () = 0;

    /*从电位中计算第m个像素的场。*/
    virtual void fieldShift () = 0;

    /*从电位中计算第m个像素的场。*/
    virtual void fieldEvaluate (long int m) = 0;

    /*采样字段并将其保存到给定的文件中。*/
    virtual void fieldSample () = 0;

    /*将字段可视化为整个域上的vtk文件，并将它们保存到给定名称的文件中。*/
    virtual void fieldVisualizeAllDomain 	(unsigned int ivtk) = 0;

    /*将字段可视化为平面中的vtk文件，并将它们保存到具有给定名称的文件中。*/
    virtual void fieldVisualizeInPlane 		(unsigned int ivtk) = 0;

    /*将字段可视化为垂直于x轴的平面上的vtk文件，并将它们保存到文件中
    *名字。*/
    virtual void fieldVisualizeInPlaneXNormal 	(unsigned int ivtk) = 0;

    /*将字段可视化为垂直于y轴的平面上的vtk文件，并将它们保存到文件中
    *名字。*/
    virtual void fieldVisualizeInPlaneYNormal 	(unsigned int ivtk) = 0;

    /*将字段可视化为垂直于z轴的平面上的vtk文件，并将它们保存到文件中
    *名字。*/
    virtual void fieldVisualizeInPlaneZNormal 	(unsigned int ivtk) = 0;

    /*将字段的总概要文件写入给定的文件名中。*/
    virtual void fieldProfile () = 0;

    /*计算静态波动器的场。*/
    void 	staticUndulator			(UpdateBunchParallel& ubp, typename std::vector<Undulator>::iterator& iter);

    /*计算平面波的场。*/
    template<class T>
    void 	planeWave			(UpdateBunchParallel& ubp, T& s);

    /*计算截断平面波的场。*/
    template<class T>
    void 	planeWaveTruncated		(UpdateBunchParallel& ubp, T& s);

    /*计算高斯光束的场。*/
    template<class T>
    void 	gaussianBeam			(UpdateBunchParallel& ubp, T& s);

    /*计算超高斯光束的场。*/
    template<class T>
    void 	superGaussianBeam		(UpdateBunchParallel& ubp, T& s);

    /*计算驻平面波的场。*/
    template<class T>
    void 	standingPlaneWave		(UpdateBunchParallel& ubp, T& s);

    /*计算截断驻平面波的场。*/
    template<class T>
    void 	standingPlaneWaveTruncated	(UpdateBunchParallel& ubp, T& s);

    /*计算固定高斯光束的场。*/
    template<class T>
    void 	standingGaussianBeam		(UpdateBunchParallel& ubp, T& s);

    /*计算固定高斯光束的场。*/
    template<class T>
    void 	standingSuperGaussianBeam	(UpdateBunchParallel& ubp, T& s);


    /****************************************************************************************************
    * FdTd代码中所需参数的列表。
    ****************************************************************************************************/

    /*仿真的解析参数分为四类：mesh， bunch, seed, and
    *und。*/
    Mesh& 								mesh_;
    Bunch&								bunch_;
    Seed&								seed_;
    std::vector<Undulator>&						undulator_;
    std::vector<ExtField>&                                              extField_;
    std::vector<FreeElectronLaser>&					FEL_;

    /*计算网格中节点在三个不同时间点的矢量势。的
    *矢量anp1_也用于存储场移位后的电流。*/
    std::vector<FieldVector<Double> >* 					anp1_;
    std::vector<FieldVector<Double> >* 					an_;
    std::vector<FieldVector<Double> >* 					anm1_;

    /*计算网格中节点在三个不同时间点的静态电位。*/
    std::vector<Double>* 						fnp1_;
    std::vector<Double>* 						fn_;
    std::vector<Double>* 						fnm1_;

    /*布尔向量决定粒子的包含。*/
    std::vector<bool>                                                   pic_;

    /*计算网格中节点在三个不同时间点的矢量势。*/
    std::vector<FieldVector<float> > 					en_;
    std::vector<FieldVector<float> > 					bn_;

    /*每个方向上的节点数。*/
    int									N0_, N1_, N2_, N1N0_;

    /*模拟中的电荷总数。*/
    unsigned int							Nc_;

    /*特定处理器中沿z的节点数和第一列的索引。*/
    int									np_, k0_;

    /*网格划分中关键点的Z坐标。*/
    Double								zp_ [2];

    /*计算网格的边界。*/
    Double                                                              xmin_, xmax_;
    Double                                                              ymin_, ymax_;
    Double                                                              zmin_, zmax_;

    /*字段计算的时间和时间步长数。*/
    Double								timep1_;
    Double								time_;
    Double								timem1_;
    unsigned int 							nTime_;

    /*电荷结构的矢量，包含电荷点的位置和动量。*/
    std::list<Charge>							chargeVectorn_;

    /*时间和时间步长数的束计算。*/
    Double								timeBunch_;
    unsigned int 							nTimeBunch_;

    /*每个字段更新中的束更新数。*/
    Double 								nUpdateBunch_;

    /*由第一个波动参数导出的运动帧的gamma， beta和dt因子。*/
    Double								gamma_;
    Double								beta_;
    Double								dt_;

    /*用于通信收费的MPI数据类型。*/
    MPI_Datatype 							MPI_CHARGE;

    /*定义一个包含更新值所需参数的结构。这些参数是
    *在类中定义一次，以避免每次字段更新时都声明它们。*/
    UpdateField								uf_;

    /*定义一个结构，其中包含采样字段和值所需的参数。这些
    *参数在类中定义一次，以避免每次字段更新时都声明它们。*/
    SampleField								sf_;

    /*定义一个结构，其中包含可视化字段配置文件所需的参数。这些
    *参数在类中定义一次，以避免每次字段更新时都声明它们。*/
    std::vector<VisualizeField>						vf_;

    /*定义一个包含保存字段概要文件所需参数的结构。这些
    *参数在类中定义一次，以避免每次字段更新时都声明它们。*/
    ProfileField							pf_;

    /*定义一个包含更新束分布所需参数的结构。这些
    *参数在类中定义一次，以避免每次字段更新时都声明它们。*/
    UpdateBunch						                ub_;

    /*定义一个结构，其中包含采样组值所需的参数。这些参数
    *在类中只定义一次，以避免每次字段更新时都声明它们。*/
    SampleBunch								sb_;

    /*定义一个包含可视化束分布所需参数的结构。这些
    *参数在类中定义一次，以避免每次字段更新时都声明它们。*/
    VisualizeBunch							vb_;

    /*定义一个包含保存堆配置文件所需参数的结构。这些参数
    *在类中只定义一次，以避免每次字段更新时都声明它们。*/
    ProfileBunch							pb_;

    /*定义一个结构，其中包含更新当前对象所需的参数。这些参数
    *在类中只定义一次，以避免每次字段更新时都声明它们。*/
    UpdateCurrent						        uc_;

    /*定义一个包含采样FEL辐射功率值所需参数的结构。
    *这些参数在类中定义一次，以避免每次字段都声明它们
    *更新。*/
    std::vector<SampleRadiationPower>				        rp_;

    /*定义一个包含采样FEL辐射功率值所需参数的结构。
    *这些参数在类中定义一次，以避免每次字段都声明它们
    *更新。*/
    std::vector<SampleRadiationEnergy>				        re_;

    /*定义一个结构，保存所有lab系探测器值的结构*/
    std ::vector<SampleRadiationDetector>   rd_;

    /*定义一个结构，其中包含存储击中屏幕的粒子所需的参数
    *这些参数在类中定义一次，以避免每次都声明它们
    *更新。*/
    std::vector<SampleScreenProfile>				        scrp_;

    /*定义MPI变量的值。*/
    int									rank_, size_;
    int									rankB_, rankF_;

    /*光速在给定的长度尺度和时间尺度下的值。*/
    Double								c0_, m0_, e0_;

    /*相对论变换框架*/
    BoostFrameTransform boostFrame_;
  };

}

#endif
