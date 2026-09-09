/********************************************************************************************************
* stdinclude.hh，代码中使用的不同标准数据库集。
********************************************************************************************************/

#ifndef STDINCLUDE_HH_
#define STDINCLUDE_HH_

#include <mpi.h>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "fieldvector.h"

namespace MITHRA
{
  /*定义时域信号的类型。*/
  enum SignalType       	{NEUMANN, GAUSSIAN, SECANT, FLATTOP, INVGAUSSIAN};

  /*定义激励的类型。*/
  enum SeedType   		{PLANEWAVE, PLANEWAVETRUNCATED, GAUSSIANBEAM, SUPERGAUSSIANBEAM,
				 STANDINGPLANEWAVE, STANDINGPLANEWAVETRUNCATED, STANDINGGAUSSIANBEAM, STANDINGSUPERGAUSSIANBEAM};

  /*定义外部字段的类型。*/
  enum ExtFieldType             {EMWAVE};

  /*定义激励的类型。*/
  enum SamplingType     	{ATPOINT, OVERLINE, INPLANE, ALLDOMAIN};

  /*定义垂直于轴的平面类型。*/
  enum PlaneType     		{XNORMAL, YNORMAL, ZNORMAL};

  /*定义要采样的字段的类型。*/
  enum FieldType		{Ex, Ey, Ez, Bx, By, Bz, Ax, Ay, Az, F};

  /*定义代码支持的波动器的类型。*/
  enum UndulatorType		{STATIC, OPTICAL, DIPOLE};

  /*定义用于FEL相互作用的求解器的类型。*/
  enum SolverType		{FD, NSFD};

  /*科学常数。*/
  const Double PI           	= 3.1415926535;
  const Double TPI		= 2.0 * PI;
  const Double EPSILON_ZERO  	= 8.85418782e-12;
  const Double MU_ZERO        	= 4.0 * PI * 1.0e-7;
  const Double C0             	= 1.0 / sqrt(EPSILON_ZERO * MU_ZERO);
  const Double Z0             	= sqrt(MU_ZERO / EPSILON_ZERO);
  const Double EC             	= 1.602e-19;
  const Double EM             	= 9.109e-31;
  const Double HB              	= 1.054e-34;
  const Double KB             	= 1.381e-23;

  /*文件后缀。*/
  const std::string VTS_FILE_SUFFIX = ".vts";
  const std::string PTS_FILE_SUFFIX = ".pvts";
  const std::string PTU_FILE_SUFFIX = ".pvtu";
  const std::string TXT_FILE_SUFFIX = ".txt";
  const std::string VTU_FILE_SUFFIX = ".vtu";

  /*单位虚数。*/
  const Complex I = Complex (0.0,1.0);

  /*检查路径是否为绝对路径（以‘/’开头）。*/
  inline bool isabsolute(std::string filename)
  {
    return (filename.compare(0,1,"/") == 0);
  }

  /*返回参数的符号。*/
  template <typename T> inline int signof(T x)
  {
    return ( (x > 0) ? 1 : ( (x < 0) ? -1 : 0 ) );
  }

  /*在终端窗口上打印一条消息。*/
  inline void printmessage(std::string filename, unsigned int linenumber, std::string message)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank == 0)
      {
	/*设置要为消息创建的流。*/
	std::stringstream printedMessage;
	printedMessage.clear();

	/*记下时间。*/
	time_t rawtime; time(&rawtime);
	std::string timeStr = ctime(&rawtime);

	/*将时间打印到字符串中。*/
	printedMessage << timeStr.substr(0,timeStr.size()-1) << " ::: ";

	/*将缩短的文件名和linenumber打印到字符串中。*/
	const char * elem = filename.c_str();
	const char * shortfn = elem;
	while ( *elem != '\0' ){
	    if ( *elem == '/' )
	      shortfn = elem + 1;
	    elem = elem + 1;
	}
	printedMessage  << shortfn << ":" << linenumber << " ::: \t \t " << message;

	/*如果流中写入了内容，则在终端上打印。*/
	if (printedMessage.str().length() > 0 )     std::cout << printedMessage.str() << std::endl;
      }
  }

  /*将任意数字转换为字符串。*/
  template<typename numbertype>
  inline std::string stringify(numbertype value)
  {
    std::ostringstream oStream;
    try
    {
	oStream << value;
    }
    catch(std::exception& error)
    {
	printmessage(std::string(__FILE__), __LINE__, std::string("Cannot convert this variable to a string."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Error:") + error.what());
	exit(1);
    }

    return oStream.str();
  }

  /*定义每个充电点的结构。*/
  struct Charge
  {

    Double			q;		/*点的电荷以电子电荷为单位。*/
    FieldVector<Double>		rnp, rnm;	/*电荷的位置矢量。*/
    FieldVector<Double>		gb;		/*电荷的标准化速度矢量。*/

    /*确定粒子是否正在通过波动器的入口点的双标志。这个标志
    *可以用来更好地促进束的运动框架。我们需要考虑它是两倍的，
    *因为这个标志需要在群集更新期间进行通信。*/
    Double			e;

    Double w;   // 软删除权重，1->0
    Double wm;  // 上一步权重

    /*输入文件中的粒子编号（从1开始）；非文件粒子为0。*/
    unsigned long long id;

    Charge();

  };

  /*检查保存数据的目录是否存在。*/
  bool pathExist (std::string const & s);

  /*将文件名拆分为两个字符串，包括其路径和文件名。*/
  void splitFilename (std::string const & str, std::string & path, std::string & file);

  /*检查文件名引用的目录是否存在。如果没有，创建目录。*/
  void createDirectory (std::string filename, unsigned int rank);

  /*函数创建一组霍尔顿序列，用于随机粒子的生成。*/
  Double halton (unsigned int i, unsigned int j);

  /*初始化密特拉时的Hello消息。*/
  void helloMessage();

  /*两个双精度值之间的正余数。*/
  Double pmod ( const Double& a, const Double& b);

}
#endif
