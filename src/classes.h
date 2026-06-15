/********************************************************************************************************
*classes.hh。mithra中使用的类的实现
********************************************************************************************************/

#ifndef CLASSES_H_
#define CLASSES_H_

#include <list>
#include <stdlib.h>
#include <time.h>

#include "database.h"
#include "fieldvector.h"
#include "stdinclude.h"

namespace MITHRA
{

  /*结构，包含所有已解析的网格参数。*/
  struct Mesh
  {
    /*与网格相关的解析数据（ls）。*/
    Double		        lengthScale_;

    /*计算域（x0, y0, z0）的中心坐标。*/
    FieldVector<Double>		meshCenter_;

    /*网格盒在三个维度（lx, ly, lz）中的长度。*/
    FieldVector<Double>		meshLength_;

    /*三维（dx, dy, dz）的网格分辨率。*/
    FieldVector<Double>		meshResolution_;

    /*解析后的数据与时间推进方案有关。*/
    Double    		        timeScale_;
    Double		        timeStep_;
    Double		        totalTime_;
    Double		        totalDist_;

    /*有限差分网格的截断顺序。可以是1也可以是2。*/
    unsigned int	        truncationOrder_;

    /*返回空间电荷假设状态的布尔标志。*/
    bool			spaceCharge_;

    /*布尔值，决定是否在初始设置中移动数据集。*/
    bool			optimizePosition_;

    /*相对于标准初始条件的初始时间位移。*/
    Double			timeShift_;

    /*将用于更新字段值的求解器类型。*/
    SolverType			solver_;

    /*用户给出的运动网格的洛伦兹因子。*/
    Double			gamma_;

    /*显示存储的网格值。*/
    void show ();

    /*在网格初始化器中初始化参数。*/
    void initialize ();
  };

  /*类中包含用于Bunch属性的数据和用于初始化类的函数
  *束，也评估束的性质。*/
  class Bunch
  {

  public:
    typedef std::list <Charge> ChargeVector;

    /*构造函数清除并初始化内部数据结构。*/
    Bunch ();

    /*用手动类型初始化一堆。这一束产生一个电荷，等于cloudCharge_。*/
    void initializeManual (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

    /*在实验室框架中用椭球类型初始化一束。这群人产生了很多
    *电荷等于numberOfParticles_与总电荷等于cloudCharge_这是
    *分布在一个椭球体中，尺寸由sigmposition_给出，中心由
    *位置向量。粒子具有以初始能量为中心的均匀能量分布
    *方差由sigmaGammaBeta_确定。*/
    void initializeEllipsoid (BunchInitialize bunchInit, ChargeVector & chargeVector, int rank, int size, int ia);

    /*用3d水晶类型初始化一束。这一堆电荷相等
    *到numberOfParticles_与总电荷等于cloudCharge_排列在a
    * 3D晶体。每个方向上的粒子数用数字_表示。因此,
    * numberOfParticles_应该是这三个数的乘积的倍数。这个比率给出
    *每个晶体点的粒子数。每个粒子的位置由
    *晶格常数，晶体以位置向量为中心。在每一点，电荷
    *在晶体周围有一个小的高斯分布。*/
    void initialize3DCrystal (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

    /*用一个文件类型初始化一个束。这个集合产生从给定文件读取的许多电荷。
    *初始化的次数等于文本文件中表的垂直长度。的
    *文件格式应包含电荷值、3个位置坐标和3个动量坐标
    *电荷分布。*/
    void initializeFile (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

    /****************************************************************************************************/

    /*包含簇初始化参数的数据。*/
    std::vector<BunchInitialize>	bunchInit_;

    /*为整个项目解析的目录。*/
    std::string				directory_;

    /*书写电子加速度分析输出的基本名称。*/
    std::string				basename_;

    /*布尔变量，决定是否应该进行堆抽样。*/
    bool				sampling_;

    /*存储更新电子运动的时间步长。*/
    Double				timeStep_;

    /*在输出文件中写入串宏观值的节奏。*/
    Double         			rhythm_;

    /*布尔参数，确定是否应该进行vtk可视化。*/
    bool				bunchVTK_;

    /*应该保存一堆vtk文件的目录。*/
    std::string				bunchVTKDirectory_;

    /*保存vtk可视化文件的文件名。*/
    std::string				bunchVTKBasename_;

    /*生成vtk文件的节奏。它应该比时间步长大两倍。*/
    Double				bunchVTKRhythm_;

    /*布尔参数，用于确定是否应该保存束配置文件。*/
    bool				bunchProfile_;

    /*应该保存束配置文件的目录。*/
    std::string				bunchProfileDirectory_;

    /*应该保存束配置文件的文件名。*/
    std::string				bunchProfileBasename_;

    /*应该保存束剖面的时间点矢量。*/
    std::vector<Double>      		bunchProfileTime_;

    /*保存束配置文件的节奏。它应该是一个比时间步长大的双精度值。*/
    Double				bunchProfileRhythm_;

    /*波动器的位置从串初始化实例开始。*/
    Double				zu_;

    /*运动框架在静止实验室框架中的Beta。*/
    Double				beta_;

    /*显示集群的存储值。*/
    void show ();
  };

  /*定义主信号类。*/
  class Signal
  {
  public:

    /*初始化参数值。*/
    Signal ();

    /*带有信号类型、时间偏移、方差、频率和载波包络相位的初始化器。*/
    void initialize (std::string type, Double l0, Double s, Double l, Double cep, unsigned int nR, std::vector<Double> sigmaInvG);

  public:

    /*信号的存储类型。*/
    SignalType				signalType_;

    /*信号的时间偏移量。*/
    Double     				t0_;

    /*信号的方差。方差根据点的强度来定义
    *信号，即信号的平方是最大值的一半。*/
    Double     				s_;

    /*调制频率。*/
    Double     				f0_;

    /*平顶脉冲的上升周期。*/
    unsigned int			nR_;

    /*调制的载波包络相位。*/
    Double     				cep_;

    /*反高斯时间剖面中的σ值。*/
    std::vector<Double>			sigmaInvG_;

    /*提供时刻t的信号。*/
    Double self (Double& t, Double& phase);

    /*显示此信号的存储值。*/
    void show ();
  };

  /*定义主种子类。*/
  class Seed
  {
  public:
    Seed ();

    void initialize (std::string        	type,
		     std::vector<Double>    	position,
		     std::vector<Double>    	direction,
		     std::vector<Double>    	polarization,
		     Double                 	amplitude,
		     std::vector<Double>	radius,
		     std::vector<int>		order,
		     Signal                   	signal);

    /*存储全域辐射场可视化所需的数据。*/
    struct vtk
    {
      bool 				sample_;
      std::vector <FieldType> 		field_;
      std::string 			directory_;
      SamplingType 			type_;
      std::string 			basename_;
      Double 				rhythm_;
      PlaneType 			plane_;
      FieldVector <Double> 		position_;

      /*初始化数据库以实现字段可视化。*/
      vtk ();
    };

  public:

    /*种子的存储类型。*/
    SeedType           			seedType_;

    /*光速在给定的长度尺度和时间尺度下的值。*/
    Double				c0_;

    /*种子的存储参考位置。对于波，它是一个参考位置，对于
    赫兹偶极子，它是偶极子的位置。*/
    FieldVector<Double>      		position_;

    /*储存种子的方向。对于波，它是传播方向，对于a
    赫兹偶极子，是偶极子的方向。*/
    FieldVector<Double>       		direction_;

    /*把波的偏振储存在种子里。*/
    FieldVector<Double>       		polarization_;

    /*储存种子的振幅。*/
    Double                  		amplitude_;

    /*存储种子的归一化振幅。*/
    Double                  		a0_;

    /*在平行和垂直方向上存储高斯光束的瑞利半径。*/
    std::vector<Double>      		radius_;

    /*存储此种子的信号类。*/
    Signal                   		signal_;

    /*超高斯光束的阶数。*/
    std::vector<int>			order_;

    /*光波动器的波长。*/
    Double				l_;

    /*光波动器激发的瑞利长度。*/
    std::vector<Double>                 zR_;

    /*洛伦兹变换的参数。*/
    Double				beta_;
    Double				gamma_;
    Double				dt_;

    /*在计算中存储所有必需的变量。*/
    Double                    		gamma;
    Double				tsignal;
    Double				d, l, zRp, wrp, zRs, wrs, x, y, x0, y0, z, p, t;
    FieldVector<Double>			rv, yv, ax, az;
    FieldVector<Double>			rl;
    Double				tl;

  public:

    /*返回任意位置和时间的电位。*/
    void fields (const FieldVector<Double>& aufpunkt, const Double& time, FieldVector<Double>& a);

    /*将采样辐射场所需的数据存储在一个点上。*/
    bool                                sampling_;
    SamplingType                        samplingType_;
    std::vector<FieldType>              samplingField_;
    std::string                         samplingDirectory_;
    std::string                         samplingBasename_;
    Double                        	samplingRhythm_;
    std::vector<FieldVector<Double> >   samplingPosition_;
    FieldVector<Double>                 samplingLineBegin_;
    FieldVector<Double>                 samplingLineEnd_;
    FieldVector<Double>                 samplingSurfaceBegin_;
    FieldVector<Double>                 samplingSurfaceEnd_;
    unsigned int                        samplingRes_;

    std::vector <vtk> vtk_;

    /*存储在all-domain中编写字段配置文件所需的数据。*/
    bool                                profile_;
    std::vector<FieldType>              profileField_;
    std::string                         profileDirectory_;
    std::string                         profileBasename_;
    std::vector<Double>                 profileTime_;
    Double                        	profileRhythm_;

    /*设置种子的采样类型。*/
    SamplingType 	samplingType 	(std::string samplingType);

    /*设置种子的采样类型。*/
    SamplingType 	vtkType 	(std::string vtkType);

    /*在平面可视化中设置vtk的平面类型。*/
    PlaneType 		planeType 	(std::string planeType);

    /*设置种子的田间采样类型。*/
    FieldType 		fieldType 	(std::string fieldType);

    /*显示此信号的存储值。*/
    void show ();
  };

  /*定义包含波动器主要参数的结构。*/
  class Undulator
  {
  public:
    Undulator ();

    /*波动器的磁场。*/
    Double				k_;

    /*静态波动的周期。*/
    Double				lu_;

    /*波动器的起始位置。*/
    Double				rb_;

    /*波荡器的长度。*/
    unsigned int			length_;

    /*二级铁磁场强度*/
    Double bd_;

    /*二级铁z方向物理长度*/
    Double ld_;

    /*束头和波动器之间的初始距离开始。*/
    Double				dist_;
    
    /*归一化速度和波动运动的等效伽马。*/
    Double				beta_;
    Double				gamma_;
    Double				dt_;

    /*波动偏振相对于x轴的角度。*/
    Double				theta_;

    /*类型的波动器，它可以是一个光学或静态波动器。*/
    UndulatorType			type_;

    /*种子的存储类型。*/
    SeedType           			seedType_;

    /*光速在给定的长度尺度和时间尺度下的值。*/
    Double				c0_;

    /*种子的存储参考位置。对于波，它是一个参考位置，对于
    赫兹偶极子，它是偶极子的位置。*/
    FieldVector<Double>      		position_;

    /*储存种子的方向。对于波，它是传播方向，对于a
    赫兹偶极子，是偶极子的方向。*/
    FieldVector<Double>       		direction_;

    /*把波的偏振储存在种子里。*/
    FieldVector<Double>       		polarization_;

    /*储存种子的振幅。*/
    Double                  		amplitude_;

    /*存储种子的归一化振幅。*/
    Double				a0_;

    /*在平行和垂直方向上存储高斯光束的瑞利半径。*/
    std::vector<Double>      		radius_;

    /*存储此种子的信号类。*/
    Signal                   		signal_;

    /*超高斯光束的阶数。*/
    std::vector<int>			order_;

    /*光波动器的波长。*/
    Double				l_;

    /*光波动器激发的瑞利长度。*/
    std::vector<Double>                 zR_;

    /*设置波动器的类型。*/
    UndulatorType undulatorType (std::string undulatorType);

    /*根据输入参数初始化波动器的数据。*/
    void initialize (std::string        	type,
		     std::vector<Double>    	position,
		     std::vector<Double>    	direction,
		     std::vector<Double>    	polarization,
		     Double                 	amplitude,
		     std::vector<Double>	radius,
		     Double			wavelength,
		     std::vector<int>		order,
		     Signal                   	signal);

    /*显示波动器的存储值。*/
    void show ();
  };

  class ExtField
  {
  public:

    ExtField();

    /*外场的类型，它可以是一个电磁场或一个腔场。*/
    ExtFieldType                        type_;

    /*种子的存储类型。*/
    SeedType                            seedType_;

    /*光速在给定的长度尺度和时间尺度下的值。*/
    Double                              c0_;

    /*种子的存储参考位置。对于波，它是一个参考位置，对于赫兹，它是一个参考位置
    *偶极子，它是偶极子的位置。*/
    FieldVector<Double>                 position_;

    /*储存种子的方向。对于波，它是传播方向，对于赫兹，它是传播方向
    *偶极子，它是偶极子的方向。*/
    FieldVector<Double>                 direction_;

    /*把波的偏振储存在种子里。*/
    FieldVector<Double>                 polarization_;

    /*储存种子的振幅。*/
    Double                              amplitude_;

    /*存储种子的归一化振幅。*/
    Double				a0_;

    /*在平行和垂直方向上存储高斯光束的瑞利半径。*/
    std::vector<Double>                 radius_;

    /*存储此种子的信号类。*/
    Signal                              signal_;

    /*超高斯光束的阶数。*/
    std::vector<int>			order_;

    /*外场的波长。*/
    Double				l_;

    /*外激励的瑞利长度。*/
    std::vector<Double>                 zR_;

    /*根据输入参数初始化波动器的数据。*/
    void initialize (std::string                type,
		     std::vector<Double>        position,
		     std::vector<Double>        direction,
		     std::vector<Double>        polarization,
		     Double                     amplitude,
		     std::vector<Double>        radius,
		     Double                     wavelength,
		     std::vector<int>        	order,
		     Signal                     signal);

    /*显示外部字段的存储值。*/
    void show ();
  };

  struct FreeElectronLaser
  {
    /*结构，包含解析后的辐射功率参数。*/
    struct RadiationSampling
    {
      /*飞机到这群人的距离来测量辐射功率。*/
      std::vector<Double>			z_;

      /*存储辐射功率采样所需的数据。*/
      bool					sampling_;

      /*开启场原始值写入的操作*/
      bool writeField_;
      bool fieldUseFloat32_;

      /*存储保存辐射数据的目录。*/
      std::string				directory_;

      /*存储保存辐射数据的文件的基本名称。*/
      std::string				basename_;

      /*线路的开始和结束以及保存辐射数据的分辨率。*/
      Double                 			lineBegin_;
      Double                 			lineEnd_;
      unsigned int                              res_;

      /*存储辐射功率的采样类型。*/
      SamplingType                        	samplingType_;

      /*谐波的波长，它的功率应该被绘制出来。*/
      std::vector<Double>			lambda_;

      /*功率计算的波长扫描数据。*/
      Double					lambdaMin_;
      Double					lambdaMax_;
      unsigned int				lambdaRes_;

      /*设置辐射功率的采样类型。*/
      void samplingType(std::string samplingType);

      /*初始化用于初始化辐射功率的值。*/
      RadiationSampling();
    };

    /*lab 系固定探测面的采样配置*/
    struct RadiationDetector 
    {
      bool sampling_;
      /*lab 系 plane-position*/
      Double zLab_;                 
      std::string directory_;
      std::string basename_;

      /*写功率开关，第一版只有这个*/
      bool writePowerLine_;         
      /*写场开关，第一版暂无*/
      bool writeField_;     
      
      /* detector 使用的单频归一化波长，第一版默认取 1.0 */
      Double lambda_;

      RadiationDetector();
    };

    /*结构，其中包含用于电力或能源可视化的已解析参数。*/
    struct RadiationVisualization
    {
      /*平面距离采样束的辐射功率或能量。*/
      Double					z_;

      /*存储采样辐射功率或能量所需的数据。*/
      bool					sampling_;

      /*存储保存辐射数据的目录。*/
      std::string				directory_;

      /*存储保存辐射数据的文件的基本名称。*/
      std::string				basename_;

      /*存储节奏，以便计算和节省辐射功率或能量。*/
      Double					rhythm_;

      /*谐波的波长，其功率或能量应绘制出来。*/
      Double					lambda_;

      /*初始化用于初始化辐射功率或能量的值。*/
      RadiationVisualization ();
    };

    /*定义包含功率采样和可视化所需结构的变量。*/
    RadiationSampling 				radiationPower_;
    RadiationVisualization 			vtkPower_;

    /*定义包含能量采样和可视化所需结构的变量。*/
    RadiationSampling 				radiationEnergy_;
    RadiationVisualization 			vtkEnergy_;

    /*定义得到探测面属性的变量*/
    RadiationDetector radiationDetector_;
    
    /*记录束配置文件的屏幕的解析参数。这将产生束配置文件
    *在实验室框架。*/
    struct ScreenProfile
    {
      /*激活在屏幕上存储采样所需的数据的标志。*/
      bool					sampling_;

      /*存储实验室框架中束配置文件保存的目录。*/
      std::string				directory_;

      /*存储保存此束配置文件的文件的基本名称。*/
      std::string				basename_;

      /*将位置存储在保存束配置文件的波动器中。*/
      std::vector<Double>      			pos_;

      /*将节奏存储在保存束配置文件的位置。*/
      Double					rhythm_;

      ScreenProfile ();
    };

    /*定义包含分析实验室框架中的束所需结构的变量。*/
    ScreenProfile 				screenProfile_;
  };
}
#endif
