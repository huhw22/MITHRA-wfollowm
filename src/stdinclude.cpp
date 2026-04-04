/********************************************************************************************************
* stinclude.cpp，代码中使用的不同标准函数集。
********************************************************************************************************/

#include "stdinclude.h"

namespace MITHRA
{

  /*检查保存数据的目录是否存在。*/
  bool pathExist (std::string const & s)
  {
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
  }

  /*将文件名拆分为两个字符串，包括其路径和文件名。*/
  void splitFilename (const std::string& str, std::string& path, std::string& file)
  {
    unsigned found = str.find_last_of("/");
    path = str.substr(0,found+1);
    file = str.substr(found+1);
  }

  /*检查文件名引用的目录是否存在。如果没有，创建目录。*/
  void createDirectory(std::string filename, unsigned int rank)
  {
    std::string path, file;
    splitFilename(filename, path, file);
    if ( !(pathExist(path)) && rank == 0 )
      if ( mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1 )
	{
	  std::cout << "Could not create the directory " << path << ". Probably the given address does not exist." << std::endl;
	  exit(1);
	}
  }

  Charge::Charge ()
  {
    q = 0.0; rnp = rnm = 0.0; gb = 0.0;
    e = 0.0;
    w  = 1.0;
    wm = 1.0;
  }

  /*函数创建一组霍尔顿序列，用于随机粒子的生成。*/
  Double halton (unsigned int i, unsigned int j)
  {
    if (i > 20)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" dimension can not be larger than 20. ") );
	exit(1);
      }

    unsigned int prime [20] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    int p0, p, k, k0, a;
    Double x = 0.0;

    k0 = j + 1;

    p = prime[i];

    p0 = p;
    k  = k0;
    x  = 0.0;
    while (k > 0)
      {
	a   = k % p;
	x  += a / (double) p0;
	k   = int (k/p);
	p0 *= p;
      }

    return 1.0 - x;
  }

  /*初始化密特拉时的Hello消息。*/
  void helloMessage()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::") );
    printmessage(std::string(__FILE__), __LINE__, std::string("MITHRA-2.0: Completely Numerical Calculation of Free Electron Laser Radiation)") );
    printmessage(std::string(__FILE__), __LINE__, std::string("Version 2.0, Copyright 2019, Arya Fallahi") );
    printmessage(std::string(__FILE__), __LINE__, std::string("Code developers: ") );
    printmessage(std::string(__FILE__), __LINE__, std::string("- Arya Fallahi ( IT'IS Foundation, Zurich, Switzerland )") );
    printmessage(std::string(__FILE__), __LINE__, std::string("- Arnau Alba   ( Paul Scherrer Institut (PSI), Villigen, Switzerland )") );
    printmessage(std::string(__FILE__), __LINE__, std::string(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::") );
  }

  /*两个双精度值之间的正余数。*/
  Double pmod ( const Double& a, const Double& b)
  {
    Double x = fmod(a, b);
    x += ( x < 0.0 ) ? b : 0.0;
    return (x);
  }
}
