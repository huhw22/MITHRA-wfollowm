/********************************************************************************************************
*class.hh：mithra中使用的类相关功能的实现
********************************************************************************************************/

#include <fstream>

#include "classes.h"

namespace MITHRA
{

  /** *网类  ***************************************************************************************/

  /*显示存储的网格值。*/
  void Mesh::show ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(" Length scale = ") + stringify(lengthScale_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Total mesh length vector = ") + stringify(meshLength_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Meshing resolution vector = ") + stringify(meshResolution_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Time scale = ") + stringify(timeScale_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Total simulation time = ") + stringify(totalTime_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Mesh truncation order = ") + stringify(truncationOrder_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Initial shift back in time = ") + stringify(timeShift_));
    if 	( spaceCharge_ == true )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Space-charge = true "));
    else if 	( spaceCharge_ == false )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Space-charge = false "));
    if 	( optimizePosition_ == true )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Optimize initial bunch position = true "));
    if 	( solver_ == NSFD )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Solver = Non-standard finite-difference "));
    else if 	( solver_ == FD )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Solver = finite-difference "));

  }

  /*在网格初始化器中初始化参数。*/
  void Mesh::initialize ()
  {
    spaceCharge_ 	=  false;
    optimizePosition_	=  false;
    solver_		=  NSFD;
    totalDist_		=  0.0;
    timeShift_		=  0.0;
    gamma_ 		= -1.0;
  }

  /** *群类  **************************************************************************************/

  /*构造函数清除并初始化内部数据结构。*/
  Bunch::Bunch ()
  {
    /*初始化bunchInit向量。*/
    bunchInit_.clear();

    /*将堆的参数初始化为一些初始值。*/
    timeStep_			= 0.0;

    /**初始化堆的参数为一些初始值。*/
    sampling_			= false;
    directory_			= "./";
    basename_             	= "";
    rhythm_			= 0.0;

    /**初始化参数为vtk可视化的一些第一个值。*/
    bunchVTK_             	= false;
    bunchVTKDirectory_    	= "./";
    bunchVTKBasename_    	= "";
    bunchVTKRhythm_      	= 0.0;

    /**初始化保存束配置文件的参数为一些初始值。*/
    bunchProfile_            	= false;
    bunchProfileDirectory_    	= "./";
    bunchProfileBasename_     	= "";
    bunchProfileTime_.clear();
    bunchProfileRhythm_		= 0.0;

    /*初始化zu和beta参数。*/
    zu_ 			= 0.0;
    beta_			= 0.0;
  }

  /*用手动类型初始化一堆。这一束产生一个电荷，等于cloudCharge_。*/
  void Bunch::initializeManual (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {
    /*声明初始化电荷矢量所需的参数。*/
    Charge        charge;

    /*确定每个电荷点的性质，并将它们添加到电荷矢量中。*/
    charge.q  	= bunchInit.cloudCharge_;
    charge.rnp  = bunchInit.position_[ia];
    charge.gb.mv( bunchInit.initialGamma_, bunchInit.betaVector_ );

    /*当且仅当该费用位于处理器的部分时，将该费用插入到费用列表中。*/
    if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
      chargeVector.push_back(charge);
  }

  /*用椭球类型初始化束。这一束产生的电荷数等于
  * numberOfParticles_与总电荷等于cloudCharge_分布在一个
  *椭球体，尺寸由sigmapposition_给出，中心由位置向量给出。的
  *粒子具有以initialEnergy_为中心的均匀能量分布，其方差由
  * sigmaGammaBeta_。*/
  void Bunch::initializeEllipsoid (BunchInitialize bunchInit, ChargeVector & chargeVector, int rank, int size, int ia)
  {
    /*如果不是4的倍数，则纠正粒子数。*/
    if ( bunchInit.numberOfParticles_ % 4 != 0 )
      {
	unsigned int n = bunchInit.numberOfParticles_ % 4;
	bunchInit.numberOfParticles_ += 4 - n;
	printmessage(std::string(__FILE__), __LINE__, std::string("Warning: The number of particles in the bunch is not a multiple of four. ") +
		     std::string("It is corrected to ") +  stringify(bunchInit.numberOfParticles_) );
      }

    /*保存最初给定的粒子数。*/
    unsigned int	Np = bunchInit.numberOfParticles_, i, Np0 = chargeVector.size();

    /*声明初始化电荷矢量所需的参数。*/
    Charge            	charge; charge.q  = bunchInit.cloudCharge_ / Np;
    FieldVector<Double> gb (0.0); gb.mv( bunchInit.initialGamma_, bunchInit.betaVector_ );
    FieldVector<Double> r  (0.0);
    FieldVector<Double> t  (0.0);
    Double            	t0, g;
    Double		zmin = 1e100;
    Double		Ne, bF, bFi;
    unsigned int	bmi;
    std::vector<Double>	randomNumbers;

    /*四粒子群的初始化只有在存在波动时才能进行
    *互动。*/
    unsigned int	ng = ( bunchInit.lambda_ == 0.0 ) ? 1 : 4;

    /*检查聚束系数。*/
    if ( bunchInit.bF_ > 2.0 || bunchInit.bF_ < 0.0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The bunching factor can not be larger than one or a negative value !!!") );
	exit(1);
      }

    /*如果生成器是随机的，我们应该确保不同的处理器不会产生相同的结果
    *随机数。*/
    if 	( bunchInit.generator_ == "random" )
      {
	/*初始化随机数生成器。*/
	srand ( time(NULL) );
	/*Np / ng * 20为最大粒子数。*/
	randomNumbers.resize( Np / ng * 20, 0.0);
	for ( unsigned int ri = 0; ri < Np / ng * 20; ri++)
	  randomNumbers[ri] = ( (double) rand() ) / RAND_MAX;
      }

    /*根据输入声明生成器函数。*/
    auto generate = [&] (unsigned int n, unsigned int m) {
      if 	( bunchInit.generator_ == "random" )
	return  ( randomNumbers[ n * 2 * Np/ng + m ] );
      else
	return 	( halton(n,m) );
    };

    /*声明注入枪弹噪声的函数。*/
    auto insertCharge = [&] (Charge q) {

      for ( unsigned int ii = 0; ii < ng; ii++ )
	{
	  /*随机调制是根据被激活的弹噪声引入的。*/
	  if ( bunchInit.shotNoise_ )
	    {
	      /*获取光束数。*/
	      bmi = int( ( charge.rnp[2] - zmin ) / bunchInit.lambda_ );

	      /*得到调制的相位和幅度。*/
	      bFi = bF * sqrt( - 2.0 * log( generate( 8 , bmi ) ) );

	      q.rnp[2]  = charge.rnp[2] - bunchInit.lambda_ / 4 * ii;

	      q.rnp[2] -= bunchInit.lambda_ / PI * bFi * sin( 2.0 * PI / bunchInit.lambda_ * q.rnp[2] + 2.0 * PI * generate( 9 , bmi ) );
	    }
	  else if ( bunchInit.lambda_ != 0.0)
	    {
	      q.rnp[2]  = charge.rnp[2] - bunchInit.lambda_ / 4 * ii;

	      q.rnp[2] -= bunchInit.lambda_ / PI * bunchInit.bF_ * sin( 2.0 * PI / bunchInit.lambda_ * q.rnp[2] + bunchInit.bFP_ * PI / 180.0 );
	    }

	  /*将这个电荷设置为电荷矢量。*/
	  chargeVector.push_back(q);
	}
    };

    /*如果有射击噪声，我们需要z轴坐标的最小值
    *计算FEL桶数。*/
    if ( bunchInit.shotNoise_ )
      {
	for (i = 0; i < Np / ng; i++)
	  {
	    if ( bunchInit.distribution_ == "uniform" )
	      zmin = std::min(   ( 2.0 * generate(2, i + Np0) - 1.0 ) * bunchInit.sigmaPosition_[2] , zmin );
	    else if ( bunchInit.distribution_ == "gaussian" )
	      zmin = std::min(   bunchInit.sigmaPosition_[2] * sqrt( - 2.0 * log( generate(2, i + Np0) ) ) * sin( 2.0 * PI * generate(3, i + Np0) ) , zmin );
	    else
	      {
		printmessage(std::string(__FILE__), __LINE__, std::string("The longitudinal type is not correctly given to the code !!!") );
		exit(1);
	      }
	  }

	if ( bunchInit.distribution_ == "uniform" )
	  for ( ; i < unsigned( Np / ng * ( 1.0 + 2.0 * bunchInit.lambda_ * sqrt( 2.0 * PI ) / ( 2.0 * bunchInit.sigmaPosition_[2] ) ) ); i++)
	    {
	      t0  = 2.0 * bunchInit.lambda_ * sqrt( - 2.0 * log( generate( 2, i + Np0 ) ) ) * sin( 2.0 * PI * generate( 3, i + Np0 ) );
	      t0 += ( t0 < 0.0 ) ? ( - bunchInit.sigmaPosition_[2] ) : ( bunchInit.sigmaPosition_[2] );

	      zmin = std::min(   t0 , zmin );
	    }

	zmin = zmin + bunchInit.position_[ia][2];

	/*得到每个自由电子束流的平均电子数。*/
	Ne = bunchInit.cloudCharge_ * bunchInit.lambda_ / ( 2.0 * bunchInit.sigmaPosition_[2] );

	/*根据给定的值设置镜头噪声的聚束系数级别。*/
	bF = ( bunchInit.bF_ == 0.0 ) ? 1.0 / sqrt(Ne) : bunchInit.bF_;

	printmessage(std::string(__FILE__), __LINE__, std::string("The standard deviation of the bunching factor for the shot noise implementation is set to ") + stringify(bF) );
      }

    /*确定每个电荷点的性质，并将它们添加到电荷矢量中。*/
    for (i = rank; i < Np / ng; i += size)
      {
	/*确定横向坐标。*/
	r[0] = bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( generate(0, i + Np0) ) ) * cos( 2.0 * PI * generate(1, i + Np0) );
	r[1] = bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( generate(0, i + Np0) ) ) * sin( 2.0 * PI * generate(1, i + Np0) );

	/*确定纵坐标。*/
	if ( bunchInit.distribution_ == "uniform" )
	  r[2] = ( 2.0 * generate(2, i + Np0) - 1.0 ) * bunchInit.sigmaPosition_[2];
	else if ( bunchInit.distribution_ == "gaussian" )
	  r[2] = bunchInit.sigmaPosition_[2] * sqrt( - 2.0 * log( generate(2, i + Np0) ) ) * sin( 2.0 * PI * generate(3, i + Np0) );
	else
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The longitudinal type is not correctly given to the code !!!") );
	    exit(1);
	  }

	/*确定横向动量。*/
	t[0] = bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( generate(4, i + Np0) ) ) * cos( 2.0 * PI * generate(5, i + Np0) );
	t[1] = bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( generate(4, i + Np0) ) ) * sin( 2.0 * PI * generate(5, i + Np0) );
	t[2] = bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( generate(6, i + Np0) ) ) * cos( 2.0 * PI * generate(7, i + Np0) );

	if ( fabs(r[0]) < bunchInit.tranTrun_ && fabs(r[1]) < bunchInit.tranTrun_ && fabs(r[2]) < bunchInit.longTrun_)
	  {
	    /*将产生的电荷移动到中心位置和动量空间。*/
	    charge.rnp    = bunchInit.position_[ia];
	    charge.rnp   += r;

	    charge.gb   = gb;
	    charge.gb  += t;

	    /*将这个电荷和镜像电荷插入到电荷矢量中。*/
	    insertCharge(charge);
	  }
      }

    /*如果束的纵向类型是均匀的，则需要添加一个锥形部分以去除
    * CSE来自一群人的尾部。*/
    if ( bunchInit.distribution_ == "uniform" )
      for ( ; i < unsigned( Np / ng * ( 1.0 + 2.0 * bunchInit.lambda_ * sqrt( 2.0 * PI ) / ( 2.0 * bunchInit.sigmaPosition_[2] ) ) ); i += size)
	{
	  r[0] = bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( generate(0, i + Np0) ) ) * cos( 2.0 * PI * generate(1, i + Np0) );
	  r[1] = bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( generate(0, i + Np0) ) ) * sin( 2.0 * PI * generate(1, i + Np0) );

	  /*确定纵坐标。*/
	  r[2] = 2.0 * bunchInit.lambda_ * sqrt( - 2.0 * log( generate(2, i + Np0) ) ) * sin( 2.0 * PI * generate(3, i + Np0) );
	  r[2] += ( r[2] < 0.0 ) ? ( - bunchInit.sigmaPosition_[2] ) : ( bunchInit.sigmaPosition_[2] );

	  /*确定横向动量。*/
	  t[0] = bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( generate(4, i + Np0) ) ) * cos( 2.0 * PI * generate(5, i + Np0) );
	  t[1] = bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( generate(4, i + Np0) ) ) * sin( 2.0 * PI * generate(5, i + Np0) );
	  t[2] = bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( generate(6, i + Np0) ) ) * cos( 2.0 * PI * generate(7, i + Np0) );

	  if ( fabs(r[0]) < bunchInit.tranTrun_ && fabs(r[1]) < bunchInit.tranTrun_ && fabs(r[2]) < bunchInit.longTrun_)
	    {
	      /*将产生的电荷移动到中心位置和动量空间。*/
	      charge.rnp   = bunchInit.position_[ia];
	      charge.rnp  += r;

	      charge.gb  = gb;
	      charge.gb += t;

	      /*将这个电荷和镜像电荷插入到电荷矢量中。*/
	      insertCharge(charge);
	    }
	}

    /*根据安装的粒子数重新设置粒子数变量的值
    *宏粒子，并执行相应的更改。*/
    bunchInit.numberOfParticles_ = chargeVector.size();
  }

  /*用3d水晶类型初始化一束。这一束产生的电荷数等于
  * numberOfParticles_与总电荷等于cloudCharge_排列在一个3D晶体。
  *每个方向的粒子数由numbers_给出。因此，numberOfParticles_应该
  *是这三个数的乘积的倍数。这个比率给出了粒子的数量
  *每个水晶点。每个粒子的位置由晶格常数和晶体决定
  *以position_向量为中心。在每一点上，电荷都有一个小的高斯分布
  *在水晶周围。*/
  void Bunch::initialize3DCrystal (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {
    /*检查numberOfParticles_是否是numbers_中值的乘积的倍数。*/
    if ( bunchInit.numberOfParticles_ % (bunchInit.numbers_[0] * bunchInit.numbers_[1] * bunchInit.numbers_[2]) != 0 )
      {
	printmessage(std::string(__FILE__), __LINE__,
		     std::string("The number of the particles and their lattice numbers do not match !!!") );
	exit(1);
      }

    /*声明初始化电荷矢量所需的参数。*/
    Charge            	charge;
    FieldVector<Double> 	gb (0.0);
    gb.mv( bunchInit.initialGamma_, bunchInit.betaVector_ );
    unsigned int np = bunchInit.numberOfParticles_ / (bunchInit.numbers_[0] * bunchInit.numbers_[1] * bunchInit.numbers_[2]);

    /*清除电荷矢量以添加电荷。*/
    chargeVector.clear();

    /*确定每个电荷点的性质，并将它们添加到电荷矢量中。*/
    for (unsigned int i = 0; i < bunchInit.numbers_[0]; i++)
      {
	for (unsigned int j = 0; j < bunchInit.numbers_[1]; j++)
	  {
	    for (unsigned int k = 0; k < bunchInit.numbers_[2]; k++)
	      {
		for (unsigned int l = 0; l < np; l++)
		  {
		    charge.q  = bunchInit.cloudCharge_ / bunchInit.numberOfParticles_;

		    charge.rnp[0]  = bunchInit.position_[ia][0] + ( i + 1.0 - 0.5 * bunchInit.numbers_[0] ) * bunchInit.latticeConstants_[0];
		    charge.rnp[1]  = bunchInit.position_[ia][1] + ( j + 1.0 - 0.5 * bunchInit.numbers_[1] ) * bunchInit.latticeConstants_[1];
		    charge.rnp[2]  = bunchInit.position_[ia][2] + ( k + 1.0 - 0.5 * bunchInit.numbers_[2] ) * bunchInit.latticeConstants_[2];
		    charge.rnp[0] += 0.5 * bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( halton(0,i) ) ) * sin( 2.0 * PI * halton(1,i) );
		    charge.rnp[1] += 0.5 * bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( halton(2,i) ) ) * sin( 2.0 * PI * halton(3,i) );
		    charge.rnp[2] += 0.5 * bunchInit.sigmaPosition_[2] * sqrt( - 2.0 * log( halton(4,i) ) ) * sin( 2.0 * PI * halton(5,i) );

		    charge.gb     = gb;
		    charge.gb[0] += bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( halton(6,i) ) ) * sin( 2.0 * PI * halton(7,i) );
		    charge.gb[1] += bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( halton(8,i) ) ) * sin( 2.0 * PI * halton(9,i) );
		    charge.gb[2] += bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(10,i)) ) * sin( 2.0 * PI * halton(11,i));

		    /*当且仅当该费用存在于处理器中时，将该费用插入到费用列表中
        *部分。*/
		    if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
		      chargeVector.push_back(charge);
		  }
	      }
	  }
      }
  }

  /*用一个文件类型初始化一个束。这个集合产生从给定文件读取的许多电荷。
  *初始化的次数等于文本文件中表的垂直长度。的
  *文件格式应包含电荷值、3个位置坐标和3个动量坐标
  *电荷分布。*/
  void Bunch::initializeFile (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {

    /*声明初始化电荷矢量所需的参数。*/
    Charge                    	charge;
    int 			saveRank = 0;
    bool 			flag = false;
    unsigned long long          particleId = 0;

    /*清除电荷矢量以添加电荷。*/
    chargeVector.clear();

    /*读取文件，并根据保存的值填充位置和电场矢量。*/
    std::ifstream myfile ( bunchInit.fileName_.c_str() );

    charge.q  = bunchInit.cloudCharge_ / bunchInit.numberOfParticles_;

    while (myfile >> charge.rnp[0]
                  >> charge.rnp[1]
                  >> charge.rnp[2]
                  >> charge.gb[0]
                  >> charge.gb[1]
                  >> charge.gb[2])
      {
	/* 使用输入文件中的有效数据行顺序作为稳定的粒子编号。 */
	charge.id = ++particleId;

	charge.rnp += bunchInit.position_[ia];

	/*将此费用仅插入到一个处理器中的费用列表中。*/
	if (saveRank == rank)
	  chargeVector.push_back(charge);
	saveRank = ( saveRank == size - 1 ) ? 0 : saveRank + 1;

	if ( bunchInit.tranTrun_ > 0.0 )
	  if ( fabs(charge.rnp[0]) > bunchInit.tranTrun_ || fabs(charge.rnp[1]) > bunchInit.tranTrun_ )
	    flag = true;
      }

    /*编写关于截断参数的警告。*/
    if ( flag )
      printmessage(std::string(__FILE__), __LINE__, std::string("Warning: Some particle coordinates are out of the transverse truncation length for the bunch. "
    	      "The results may be inaccurate !!!") );

    /*计算安装粒子的总量。*/
    unsigned int NqL = chargeVector.size(), NqG = 0;
    MPI_Reduce(&NqL,&NqG,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

    /*检查电荷矢量的大小与粒子的数量。*/
    if ( bunchInit.numberOfParticles_ != NqG && rank == 0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The number of the particles and the file size do not match !!! The file contains " + stringify(NqG) + " particles.") );
	exit(1);
      }
  }

  /*显示集群的存储值。*/
  void Bunch::show ()
  {
    for (unsigned int i = 0; i < bunchInit_.size(); i++)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Type of the bunch initialization = ") + bunchInit_[i].bunchType_ );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Type of the current profile = ") + bunchInit_[i].distribution_ );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Number of macro-particles = ") + stringify(bunchInit_[i].numberOfParticles_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Total number of electrons in the cloud = ") + stringify(bunchInit_[i].cloudCharge_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Total charge of the cloud [Coulombs] = ") + stringify(-bunchInit_[i].cloudCharge_ * EC) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial mean gamma of the bunch = ") + stringify(bunchInit_[i].initialGamma_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial mean speed of the bunch = ") + stringify(bunchInit_[i].initialBeta_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial direction of the bunch speed = ") + stringify(bunchInit_[i].initialDirection_) );
	for ( unsigned int ia = 0; ia < bunchInit_[i].position_.size(); ia++)
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Initial position of the bunch = ") + stringify(bunchInit_[i].position_[ia]));
	printmessage(std::string(__FILE__), __LINE__, std::string(" Position spread of the bunch = ") + stringify(bunchInit_[i].sigmaPosition_));
	printmessage(std::string(__FILE__), __LINE__, std::string(" Momentum spread of the bunch = ") + stringify(bunchInit_[i].sigmaGammaBeta_));
	if ( bunchInit_[i].bunchType_ == "3D-crystal" )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Number of crystal points = ") + stringify(bunchInit_[i].numbers_) );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Lattice constant of the crystal = ") + stringify(bunchInit_[i].latticeConstants_) );
	  }
	if ( bunchInit_[i].bunchType_ == "file" )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" File name of the bunch = ") + bunchInit_[i].fileName_ );
	  }
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial bunching factor in the bunch = ") + stringify(bunchInit_[i].bF_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Time step for the bunch calculation = ") + stringify(timeStep_) );
      }
    if ( sampling_ )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Sampling of the bunch data is enabled. ") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Directory for the bunch sampling = ") + stringify(directory_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Base name for the bunch sampling = ") + stringify(basename_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Rhythm the bunch sampling = ") + stringify(rhythm_) );
      }
    if ( bunchVTK_ )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" vtk visualization of the bunch is enabled. ") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Directory for the bunch visualization = ") + stringify(bunchVTKDirectory_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Base name for the bunch visualization = ") + stringify(bunchVTKBasename_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Rhythm the bunch visualization = ") + stringify(bunchVTKRhythm_) );
      }
    if ( bunchProfile_ )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Save the bunch profile is enabled. ") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Directory for the bunch profile = ") + stringify(bunchProfileDirectory_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Base name for the bunch profile = ") + stringify(bunchProfileBasename_) );
	for (unsigned int i = 0; i < bunchProfileTime_.size(); i++)
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Bunch profile will be saved at = ") + stringify(bunchProfileTime_[i]) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Rhythm the bunch profiling = ") + stringify(bunchProfileRhythm_) );
      }
  }

  /** 信号类  *************************************************************************************/

  /*初始化参数值。*/
  Signal::Signal ()
  {
    t0_ 		= 0.0;
    s_  		= 0.0;
    f0_ 		= 1.0;
    cep_		= 0.0;
    signalType_		= GAUSSIAN;
    nR_			= 1;
    sigmaInvG_.resize(2,0.0);
  }

  /*带有信号类型、时间偏移、方差、频率和载波包络相位的初始化器。*/
  void Signal::initialize (std::string type, Double l0, Double s, Double l, Double cep, unsigned int nR, std::vector<Double> sigmaInvG)
  {
    /*初始化信号类型。*/
    if      ( type.compare("neumann") == 0 )            signalType_ = NEUMANN;
    else if ( type.compare("gaussian") == 0 )           signalType_ = GAUSSIAN;
    else if ( type.compare("secant-hyperbolic") == 0 )  signalType_ = SECANT;
    else if ( type.compare("flat-top") == 0 )   	signalType_ = FLATTOP;
    else if ( type.compare("inverse-gaussian") == 0 )   signalType_ = INVGAUSSIAN;
    else { std::cout << type << " is an unknown signal type for the given set of parameters." << std::endl; exit(1); }

    /*初始化载波的时延、方差和频率。*/
    t0_ = l0;
    s_  = s;
    f0_ = 1 / l;

    /*初始化载波信封阶段。*/
    cep_ = cep * PI / 180;

    /*初始化平顶脉冲的上升周期数。*/
    nR_  = nR;

    /*初始化反高斯信号中sigma的值。*/
    sigmaInvG_ = sigmaInvG;

    /*检查方差是否不等于零。*/
    if (s_ == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Variance of signal is set to zero. "));
	printmessage(std::string(__FILE__), __LINE__, std::string(" This is not allowed because we divide through the variance. "));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /*检查反高斯信号的σ值是否不等于零。*/
    if ( signalType_ == INVGAUSSIAN )
      {
	if ( sigmaInvG_[0] * sigmaInvG_[1] == 0.0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" sigma of the inverse-gaussian signal is set to zero. "));
	    printmessage(std::string(__FILE__), __LINE__, std::string(" This is not allowed because we divide through the sigma value. "));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	    exit(1);
	  }
      }
  }

  Double Signal::self (Double & t, Double & phase)
  {
    /*如果信号在中心周围的20*s范围之外，则认为它是零。*/
    if ( fabs(t - t0_) > 10.0 * s_ )	return ( 0.0 );
    else
      {
	if  ( signalType_ == NEUMANN )
	  return ( - cos( 2 * PI * f0_ * (t-t0_) + cep_ + phase) * 2.7724 * (t-t0_) / (s_*s_) * exp( -1.3863 * (t-t0_) * (t-t0_) / (s_*s_) ) );

	else if  ( signalType_ == GAUSSIAN )
	  return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase) * exp( -1.3863 * pow( (t-t0_) / s_ , 2 ) ) );

	else if  ( signalType_ == SECANT )
	  return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase) / cosh( (t-t0_) / s_ ) );

	else if  ( signalType_ == FLATTOP )
	  {
	    if      ( t - t0_ <= - s_ / 2.0 )
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) * exp( - pow( (t-t0_+s_/2.0)*f0_/nR_, 2) ) );
	    else if ( t - t0_ <= s_ / 2.0   )
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) );
	    else
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) * exp( - pow( (t-t0_-s_/2.0)*f0_/nR_, 2) ) );
	  }

	else if  ( signalType_ == INVGAUSSIAN )
	  {
	    if      ( t - t0_ <= - s_ / 2.0 )
	      return ( cos( 2 * PI * f0_ * ( t - t0_ ) + cep_ + phase ) *
		       pow( ( 1.0 + pow( ( t - t0_ ) / sigmaInvG_[0] , 2 ) ) * ( 1.0 + pow( ( t - t0_ ) / sigmaInvG_[1] , 2 ) ) , 0.25 ) *
		       exp( - pow( ( t - t0_ + s_/2.0 ) * f0_ / nR_, 2 ) ) );
	    else if ( t - t0_ <= s_ / 2.0   )
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) *
		       pow( ( 1.0 + pow( ( t - t0_ ) / sigmaInvG_[0] , 2 ) ) * ( 1.0 + pow( ( t - t0_ ) / sigmaInvG_[1] , 2 ) ) , 0.25 ) );
	    else
	      return ( cos( 2 * PI * f0_ * ( t - t0_ ) + cep_ + phase ) *
		       pow( ( 1.0 + pow( ( t - t0_ ) / sigmaInvG_[0] , 2 ) ) * ( 1.0 + pow( ( t - t0_ ) / sigmaInvG_[1] , 2 ) ) , 0.25 ) *
		       exp( - pow( ( t - t0_ - s_/2.0 ) * f0_ / nR_, 2 ) ) );
	  }
      }
    return (0.0);
  }

  /*显示此信号的存储值。*/
  void Signal::show ()
  {
    if      (signalType_ == NEUMANN)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = neumann pulse"));
    else if (signalType_ == GAUSSIAN)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = gaussian pulse"));
    else if (signalType_ == SECANT)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = secant hyperbolic pulse"));
    else if (signalType_ == FLATTOP)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = flat-top pulse"));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal offset         = ") + stringify(t0_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal variance       = ") + stringify(s_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal frequency  = ") + stringify(f0_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal carrier envelope phase = ") + stringify(cep_));
  }

  /** *种子类  ***************************************************************************************/

  Seed::Seed ()
  {
    seedType_ 			= PLANEWAVE;
    c0_				= 0.0;
    position_			= 0.0;
    direction_			= 0.0;
    polarization_		= 0.0;
    amplitude_			= 0.0;
    a0_				= 0.0;
    radius_.resize(2,0.0);
    order_.resize(2,0);

    l_				= 0.0;

    zR_.resize(2,0.0);

    beta_ 			= 0.0;
    gamma_			= 1.0;
    dt_				= 0.0;

    gamma 			= 1.0;
    tsignal			= 0.0;
    d = l = zRp = wrp = zRs = wrs = x = y = x0 = y0 = z = p = t = 0.0;
    rv = yv = ax = az 		= 0.0;
    rl 				= 0.0;
    tl 				= 0.0;

    sampling_			= false;
    samplingType_		= ATPOINT;
    samplingField_.clear();
    samplingDirectory_		= "";
    samplingBasename_		= "";
    samplingRhythm_		= 0.0;
    samplingPosition_.clear();
    samplingLineBegin_		= 0.0;
    samplingLineEnd_		= 0.0;
    samplingSurfaceBegin_	= 0.0;
    samplingSurfaceEnd_		= 0.0;
    samplingRes_		= 0.0;

    vtk_.clear();

    profile_			= false;
    profileField_.clear();
    profileDirectory_		= "";
    profileBasename_		= "";
    profileTime_.clear();
    profileRhythm_		= 0.0;
  }

  void Seed::initialize (std::string        	type,
			 std::vector<Double>    position,
			 std::vector<Double>    direction,
			 std::vector<Double>    polarization,
			 Double                 a0,
			 std::vector<Double>	radius,
			 std::vector<int>	order,
			 Signal                 signal)
  {
    /*根据seedType返回的字符串设置种子类型。*/
    if      ( type.compare("plane-wave"         	  ) == 0 ) seedType_ = PLANEWAVE;
    else if ( type.compare("truncated-plane-wave"	  ) == 0 ) seedType_ = PLANEWAVETRUNCATED;
    else if ( type.compare("gaussian-beam"      	  ) == 0 ) seedType_ = GAUSSIANBEAM;
    else if ( type.compare("super-gaussian-beam"      	  ) == 0 ) seedType_ = SUPERGAUSSIANBEAM;
    else if ( type.compare("standing-plane-wave"	  ) == 0 ) seedType_ = STANDINGPLANEWAVE;
    else if ( type.compare("standing-truncated-plane-wave") == 0 ) seedType_ = STANDINGPLANEWAVETRUNCATED;
    else if ( type.compare("standing-gaussian-beam"       ) == 0 ) seedType_ = STANDINGGAUSSIANBEAM;
    else if ( type.compare("standing-super-gaussian-beam" ) == 0 ) seedType_ = STANDINGSUPERGAUSSIANBEAM;
    else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

    /*设置种子类的矢量位置、方向和极化。*/
    position_ 		= position;
    polarization_ 	= polarization;
    direction_ 		= direction;

    /*检查方向向量的长度是否为零，并对向量进行归一化。*/
    if ( direction_.norm2() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The direction vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( direction_.norm2() );
	direction_ /= vl;
      }

    /*检查偏振矢量的长度是否为零，并对矢量进行归一化。*/
    if ( polarization_.norm2() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The polarization vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = polarization_.norm();
	polarization_ /= vl;
      }

    /*检查极化和方向是否互为法向。*/
    if ( fabs( polarization_ * direction_ ) > 1.0e-50 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Polarization is not normal to the direction."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /*初始化种子的振幅。*/
    a0_ 		= a0;

    /*初始化高斯光束的瑞利半径。*/
    radius_ 		= radius;

    /*检查偏振矢量的长度是否为零，并对矢量进行归一化。*/
    if ( seedType_ == GAUSSIANBEAM || seedType_ == STANDINGGAUSSIANBEAM || seedType_ == SUPERGAUSSIANBEAM || seedType_ == STANDINGSUPERGAUSSIANBEAM )
      {
	if ( radius_[0] * radius_[1] == 0.0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	    exit(1);
	  }
      }

    /*初始化种子的信号。*/
    signal_           = signal;

    /*初始化超高斯光束的值。*/
    if ( seedType_ == SUPERGAUSSIANBEAM || seedType_ == STANDINGSUPERGAUSSIANBEAM )
      {
	order_ = order;

	Double d1 = 0.0, d2 = 0.0;
	for ( int i = -order_[0]; i <= order_[0]; i++ ) d1 += exp(-i*i);
	for ( int i = -order_[1]; i <= order_[1]; i++ ) d2 += exp(-i*i);
	radius_[0] /= order_[0] + sqrt( 1.0 - log(d1) );
	radius_[0] /= order_[0] + sqrt( 1.0 - log(d2) );
	a0_ /= d1  * d2;
      }
  }

  /*返回任意位置和时间的电位。*/
  void Seed::fields (const FieldVector<Double>& aufpunkt, const Double& time, FieldVector<Double>& a)
  {
    /*将坐标系从束静止坐标系转移到实验坐标系。*/
    rl[0] = aufpunkt[0]; rl[1] = aufpunkt[1];
    rl[2] = gamma_ * ( aufpunkt[2] + beta_ * c0_ * ( time + dt_ ) );
    tl    = gamma_ * ( time + dt_  + beta_ / c0_ * aufpunkt[2]    );

    /*计算沿传播方向到参考位置的距离。*/
    rv = rl; rv -= position_;
    z  = rv * direction_ ;

    /*计算传播延迟并从时间中减去它。*/
    tl -= z / c0_;

    /*重置脉冲的载波包络相位。*/
    p = 0.0;

    /*现在根据给定的特定种子来操作电场矢量。*/
    if ( seedType_ == PLANEWAVE )
      {
	/*在校正时间检索信号值。*/
	tsignal = signal_.self(tl, p);

	/*仅当信号值大于限制时才计算该字段。*/
	if ( fabs(tsignal) < 1.0e-6 )
	  a = 0.0;
	else
	  a.mv( amplitude_ * tsignal , polarization_ );
      }
    else if ( seedType_ == PLANEWAVETRUNCATED )
      {
	/*在校正时间检索信号值。*/
	tsignal = signal_.self(tl, p);

	/*计算到中心线的横向距离。*/
	x  = rv * polarization_;
	yv = cross(direction_, polarization_);
	y  = rv * yv;

	/*仅当信号值大于限制时才计算该字段。*/
	if ( fabs(tsignal) < 1.0e-6  || fabs(x) > radius_[0] || fabs(y) > radius_[1] )
	  a = 0.0;
	else
	  a.mv( amplitude_ * tsignal , polarization_ );
      }
    else if ( seedType_ == GAUSSIANBEAM )
      {
	/*在校正时间检索信号值。*/
	tsignal = signal_.self(tl, p);

	if ( fabs(tsignal) < 1.0e-6 ) a = 0.0;
	else
	  {
	    /*计算到中心线的横向距离。*/
	    x  = rv * polarization_;
	    yv = cross(direction_, polarization_);
	    y  = rv * yv;

	    /*计算与给定中心频率相对应的波长。*/
	    l = c0_ / signal_.f0_;

	    /*计算光束的瑞利长度和相对半径。*/
	    zRp = PI * radius_[0] * radius_[0] / l;
	    wrp = sqrt(1.0 + z * z / ( zRp * zRp ));
	    zRs = PI * radius_[1] * radius_[1] / l;
	    wrs = sqrt(1.0 + z * z / ( zRs * zRs ));

	    /*计算点和参考点之间的横向矢量。*/
	    p         = 0.5 * ( atan(z/zRp) + atan(z/zRs) - PI ) - PI*z/l * ( pow(x/(zRp*wrp),2) + pow(y/(zRs*wrs),2) );
	    tsignal   = signal_.self(tl, p);
	    t         = exp( - pow(x/(radius_[0]*wrp),2) - pow(y/(radius_[1]*wrs),2) ) / sqrt(wrs*wrp) * amplitude_;
	    a.mv( t * tsignal, polarization_);
	  }
      }
    else if ( seedType_ == SUPERGAUSSIANBEAM )
      {
	/*在校正时间检索信号值。*/
	tsignal = signal_.self(tl, p);

	if ( fabs(tsignal) < 1.0e-6 ) a = 0.0;
	else
	  {
	    /*计算到中心线的横向距离。*/
	    x  = rv * polarization_;
	    yv = cross(direction_, polarization_);
	    y  = rv * yv;

	    /*计算与给定中心频率相对应的波长。*/
	    l = c0_ / signal_.f0_;

	    /*计算光束的瑞利长度和相对半径。*/
	    zRp = PI * radius_[0] * radius_[0] / l;
	    wrp = sqrt(1.0 + z * z / ( zRp * zRp ));
	    zRs = PI * radius_[1] * radius_[1] / l;
	    wrs = sqrt(1.0 + z * z / ( zRs * zRs ));

	    /*循环超高斯光束的元素并添加它们的场。*/
	    for ( int i = - order_[0]; i <= order_[0]; i++ )
	      for ( int j = - order_[1]; j <= order_[1]; j++ )
		{
		  /*计算点和参考点之间的横向矢量。*/
		  x0 = ( x - i * radius_[0] ) / wrp;
		  y0 = ( y - j * radius_[1] ) / wrs;

		  /*计算点和参考点之间的横向矢量。*/
		  p         = 0.5 * ( atan(z/zRp) + atan(z/zRs) - PI ) - PI*z/l * ( pow(x/(zRp*wrp),2) + pow(y/(zRs*wrs),2) );
		  tsignal   = signal_.self(tl, p);
		  t         = exp( - pow(x/(radius_[0]*wrp),2) - pow(y/(radius_[1]*wrs),2) ) / sqrt(wrs*wrp) * amplitude_;
		  a.pmv( t * tsignal, polarization_);
		}
	  }
      }

    /*现在将计算得到的磁矢量势转移到束静止坐标系中。*/
    a[2] *= gamma_;
  }

  /*初始化数据库以实现字段可视化。*/
  Seed::vtk::vtk ()
  {
    sample_			= false;
    field_.clear();
    directory_			= "";
    basename_			= "";
    rhythm_			= 0.0;
    type_			= ALLDOMAIN;
    plane_			= ZNORMAL;
    position_			= 0.0;
  }

  /*设置种子的采样类型。*/
  SamplingType Seed::samplingType (std::string samplingType)
  {
    if      ( samplingType.compare("at-point")   == 0 )	return( ATPOINT  );
    else if ( samplingType.compare("over-line")  == 0 )	return( OVERLINE );
    else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /*设置种子的采样类型。*/
  SamplingType Seed::vtkType (std::string vtkType)
  {
    if      ( vtkType.compare("in-plane") == 0 )	return( INPLANE   );
    else if ( vtkType.compare("all-domain") == 0 )	return( ALLDOMAIN );
    else { std::cout << vtkType << " is an unknown vtk type." << std::endl; exit(1); }
  }

  /*在平面可视化中设置vtk的平面类型。*/
  PlaneType Seed::planeType (std::string planeType)
  {
    if      ( planeType.compare("yz")     == 0 )	return( XNORMAL );
    else if ( planeType.compare("xz")     == 0 )	return( YNORMAL );
    else if ( planeType.compare("xy")     == 0 )	return( ZNORMAL );
    else { std::cout << planeType << " is an unknown vtk plane type." << std::endl; exit(1); }
  }

  /*设置种子的田间采样类型。*/
  FieldType Seed::fieldType (std::string fieldType)
  {
    if      ( fieldType.compare("Ex") == 0 )  return(Ex);
    else if ( fieldType.compare("Ey") == 0 )  return(Ey);
    else if ( fieldType.compare("Ez") == 0 )  return(Ez);
    else if ( fieldType.compare("Bx") == 0 )  return(Bx);
    else if ( fieldType.compare("By") == 0 )  return(By);
    else if ( fieldType.compare("Bz") == 0 )  return(Bz);
    else if ( fieldType.compare("Ax") == 0 )  return(Ax);
    else if ( fieldType.compare("Ay") == 0 )  return(Ay);
    else if ( fieldType.compare("Az") == 0 )  return(Az);
    else if ( fieldType.compare("F") == 0 )   return(F);
    else { std::cout << fieldType << " is an unknown sampling field." << std::endl; exit(1); }
  }

  /*显示此信号的存储值。*/
  void Seed::show ()
  {
    if        (seedType_ == PLANEWAVE)      		printmessage(std::string(__FILE__), __LINE__, std::string("Seed type = plane-wave"));
    else if   (seedType_ == PLANEWAVETRUNCATED)   	printmessage(std::string(__FILE__), __LINE__, std::string("Seed type = truncated-plane-wave"));
    else if   (seedType_ == GAUSSIANBEAM)   		printmessage(std::string(__FILE__), __LINE__, std::string("Seed type = gaussian-beam"));
    printmessage(std::string(__FILE__), __LINE__, std::string("Seed position [")
    + stringify(position_[0]) + std::string("; ")
    + stringify(position_[1]) + std::string("; ")
    + stringify(position_[2]) + std::string("]"));
    printmessage(std::string(__FILE__), __LINE__, std::string("Seed direction [")
    + stringify(direction_[0]) + std::string("; ")
    + stringify(direction_[1]) + std::string("; ")
    + stringify(direction_[2]) + std::string("]"));
    printmessage(std::string(__FILE__), __LINE__, std::string("Seed normalized amplitude = ") + stringify(a0_));
    if      ( seedType_ == PLANEWAVE || seedType_ == GAUSSIANBEAM )
      printmessage(std::string(__FILE__), __LINE__, std::string("Seed polarization [")
    + stringify(polarization_[0]) + std::string("; ")
    + stringify(polarization_[1]) + std::string("; ")
    + stringify(polarization_[2]) + std::string("]"));
    if      ( seedType_ == GAUSSIANBEAM )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed parallel Rayleigh-radius = ") + stringify(radius_[0]));
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed perpendicular Rayleigh-radius = ") + stringify(radius_[1]));
      }
    if (sampling_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling at a position is enabled.") );
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling directory = ") +  samplingDirectory_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling base-name = ") +  samplingBasename_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling rhythm = ") +  stringify(samplingRhythm_) );
	for (unsigned int i=0; i < samplingPosition_.size(); i++)
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling position = ") +  stringify(samplingPosition_[i]) );
	  }
      }
    for (unsigned int i = 0; i < vtk_.size(); i++)
      {
	if (vtk_[i].sample_)
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization in all domain is enabled.") );
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization directory = ") +  vtk_[i].directory_);
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization base-name = ") +  vtk_[i].basename_);
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization rhythm = ") +  stringify(vtk_[i].rhythm_) );
	  }
      }
    if (profile_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile saving at a time is enabled.") );
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile directory = ") +  profileDirectory_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile base-name = ") +  profileBasename_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile rhythm = ") +  stringify(profileRhythm_) );
	for (unsigned int i=0; i < profileTime_.size(); i++)
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile time = ") +  stringify(profileTime_[i]) );
	  }
      }
    signal_.show();
  }

  /** 波荡器类  **********************************************************************************/

  Undulator::Undulator ()
  {
    k_			= 0.0;
    lu_			= 0.0;
    rb_			= 0.0;
    length_		= 0;
    dist_		= 0.0;
    bd_ = 0.0;
    ld_ = 0.0;
    beta_		= 0.0;
    gamma_		= 1.0;
    dt_			= 0.0;
    theta_		= 0.0;
    type_		= STATIC;
    seedType_ 		= PLANEWAVE;
    c0_			= 0.0;
    position_		= 0.0;
    direction_		= 0.0;
    polarization_	= 0.0;
    amplitude_		= 0.0;
    a0_			= 0.0;
    radius_.resize(2,0.0);
    order_.resize(2,0);
    l_			= 0.0;
    zR_.resize(2,0.0);
  }

  /*设置波动器的类型。*/
  UndulatorType Undulator::undulatorType (std::string undulatorType)
  {
    if      ( undulatorType.compare("static")   == 0 )   return( STATIC  );
    else if ( undulatorType.compare("optical")  == 0 )   return( OPTICAL );
    else if ( undulatorType.compare("dipole") == 0 ) return( DIPOLE );
    else { std::cout << undulatorType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /*根据输入参数初始化波动器的数据。*/
  void Undulator::initialize (std::string        	type,
			      std::vector<Double>    	position,
			      std::vector<Double>    	direction,
			      std::vector<Double>    	polarization,
			      Double                 	a0,
			      std::vector<Double>	radius,
			      Double			wavelength,
			      std::vector<int>		order,
			      Signal                   	signal)
  {
    /*根据seedType返回的字符串设置种子类型。*/
    if      ( type.compare("plane-wave"         	  ) == 0 ) seedType_ = PLANEWAVE;
    else if ( type.compare("truncated-plane-wave"	  ) == 0 ) seedType_ = PLANEWAVETRUNCATED;
    else if ( type.compare("gaussian-beam"      	  ) == 0 ) seedType_ = GAUSSIANBEAM;
    else if ( type.compare("super-gaussian-beam"      	  ) == 0 ) seedType_ = SUPERGAUSSIANBEAM;
    else if ( type.compare("standing-plane-wave"	  ) == 0 ) seedType_ = STANDINGPLANEWAVE;
    else if ( type.compare("standing-truncated-plane-wave") == 0 ) seedType_ = STANDINGPLANEWAVETRUNCATED;
    else if ( type.compare("standing-gaussian-beam"       ) == 0 ) seedType_ = STANDINGGAUSSIANBEAM;
    else if ( type.compare("standing-super-gaussian-beam" ) == 0 ) seedType_ = STANDINGSUPERGAUSSIANBEAM;
    else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

    /*设置种子类的矢量位置、方向和极化。*/
    position_ 		= position;
    polarization_ 	= polarization;
    direction_ 		= direction;

    /*检查方向向量的长度是否为零，并对向量进行归一化。*/
    if ( direction_.norm2() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The direction vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( direction_.norm2() );
	direction_ /= vl;
      }

    /*检查偏振矢量的长度是否为零，并对矢量进行归一化。*/
    if ( polarization_.norm2() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The polarization vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = polarization_.norm();
	polarization_ /= vl;
      }

    /*检查极化和方向是否互为法向。*/
    if ( fabs( polarization_ * direction_ ) > 1.0e-50 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Polarization is not normal to the direction."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /*初始化种子和波动器的振幅。*/
    a0_ 		= a0;

    /*初始化高斯光束的瑞利半径。*/
    radius_ 		= radius;

    /*根据信号的给定波长初始化波动周期。*/
    lu_			= wavelength;

    /*检查偏振矢量的长度是否为零，并对矢量进行归一化。*/
    if ( seedType_ == GAUSSIANBEAM || seedType_ == STANDINGGAUSSIANBEAM || seedType_ == SUPERGAUSSIANBEAM || seedType_ == STANDINGSUPERGAUSSIANBEAM )
      {
	if ( radius_[0] * radius_[1] == 0.0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	    exit(1);
	  }
      }

    /*初始化种子的信号。*/
    signal_           = signal;

    /*初始化超高斯光束的值。*/
    if ( seedType_ == SUPERGAUSSIANBEAM || seedType_ == STANDINGSUPERGAUSSIANBEAM )
      {
	order_ = order;

	Double d1 = 0.0, d2 = 0.0;
	for ( int i = -order_[0]; i <= order_[0]; i++ ) d1 += exp(-i*i);
	for ( int i = -order_[1]; i <= order_[1]; i++ ) d2 += exp(-i*i);
	radius_[0] /= order_[0] + sqrt( 1.0 - log(d1) );
	radius_[0] /= order_[0] + sqrt( 1.0 - log(d2) );
	a0_ /= d1  * d2;
      }
  }

  /*显示波动器的存储值。*/
  void Undulator::show ()
  {
    if ( type_ == STATIC )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator type is static.") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator parameter = ") + stringify(k_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator period = ") + stringify(lu_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator begin = ") + stringify(rb_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator length = ") + stringify(length_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Magnetic field angle with respect to x = ") + stringify(theta_ * 180 / PI));
      }
    else if ( type_ == OPTICAL )
      {
	if 		(seedType_ == PLANEWAVE)      printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = plane-wave"));
	else if 	(seedType_ == GAUSSIANBEAM)   printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = gaussian-beam"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam position [")
	+ stringify(position_[0]) + std::string("; ")
	+ stringify(position_[1]) + std::string("; ")
	+ stringify(position_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam direction [")
	+ stringify(direction_[0]) + std::string("; ")
	+ stringify(direction_[1]) + std::string("; ")
	+ stringify(direction_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam normalized amplitude = ")
	+ stringify(a0_));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam polarization [")
	+ stringify(polarization_[0]) + std::string("; ")
	+ stringify(polarization_[1]) + std::string("; ")
	+ stringify(polarization_[2]) + std::string("]"));
	if      ( seedType_ == GAUSSIANBEAM )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam parallel Rayleigh-radius = ") + stringify(radius_[0]));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam perpendicular Rayleigh-radius = ") + stringify(radius_[1]));
	  }
	signal_.show();
      }
  }

  /*** ExtField 类 ***********************************************************************************/

  ExtField::ExtField ()
  {
    type_ 		= EMWAVE;
    seedType_		= PLANEWAVE;
    c0_			= 0.0;
    position_		= 0.0;
    direction_		= 0.0;
    polarization_	= 0.0;
    amplitude_ 		= 0.0;
    a0_			= 0.0;
    radius_.resize(2,0.0);
    order_.resize(2,0);
    l_			= 0.0;
    zR_.resize(2,0.0);
  }

  /*根据输入参数初始化波动器的数据。*/
  void ExtField::initialize (std::string                type,
			     std::vector<Double>        position,
			     std::vector<Double>        direction,
			     std::vector<Double>        polarization,
			     Double                     a0,
			     std::vector<Double>        radius,
			     Double                     wavelength,
			     std::vector<int>		order,
			     Signal                     signal)
  {
    /*根据seedType返回的字符串设置种子类型。*/
    if      ( type.compare("plane-wave"         	  ) == 0 ) seedType_ = PLANEWAVE;
    else if ( type.compare("truncated-plane-wave"	  ) == 0 ) seedType_ = PLANEWAVETRUNCATED;
    else if ( type.compare("gaussian-beam"      	  ) == 0 ) seedType_ = GAUSSIANBEAM;
    else if ( type.compare("super-gaussian-beam"      	  ) == 0 ) seedType_ = SUPERGAUSSIANBEAM;
    else if ( type.compare("standing-plane-wave"	  ) == 0 ) seedType_ = STANDINGPLANEWAVE;
    else if ( type.compare("standing-truncated-plane-wave") == 0 ) seedType_ = STANDINGPLANEWAVETRUNCATED;
    else if ( type.compare("standing-gaussian-beam"       ) == 0 ) seedType_ = STANDINGGAUSSIANBEAM;
    else if ( type.compare("standing-super-gaussian-beam" ) == 0 ) seedType_ = STANDINGSUPERGAUSSIANBEAM;
    else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

    /*设置种子类的矢量位置、方向和极化。*/
    position_         = position;
    polarization_     = polarization;
    direction_        = direction;

    /*检查方向向量的长度是否为零，并对向量进行归一化。*/
    if ( direction_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The direction vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( direction_.norm() );
	direction_ /= vl;
      }

    /*检查偏振矢量的长度是否为零，并对矢量进行归一化。*/
    if ( polarization_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The polarization vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( polarization_.norm() );
	polarization_ /= vl;
      }

    /*检查极化和方向是否互为法向。*/
    if ( fabs( polarization_ * direction_ ) > 1.0e-50 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Polarization is not normal to the direction."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /*初始化种子的振幅。*/
    a0_        		= a0;

    /*初始化高斯光束的瑞利半径。*/
    radius_           	= radius;

    /*检查偏振矢量的长度是否为零，并对矢量进行归一化。*/
    if ( seedType_ == GAUSSIANBEAM || seedType_ == STANDINGGAUSSIANBEAM || seedType_ == SUPERGAUSSIANBEAM || seedType_ == STANDINGSUPERGAUSSIANBEAM )
      {
	if ( radius_[0] * radius_[1] == 0.0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	    exit(1);
	  }
      }

    /*初始化种子的信号。*/
    signal_           = signal;

    /*初始化超高斯光束的值。*/
    if ( seedType_ == SUPERGAUSSIANBEAM || seedType_ == STANDINGSUPERGAUSSIANBEAM )
      {
	order_ = order;

	Double d1 = 0.0, d2 = 0.0;
	for ( int i = -order_[0]; i <= order_[0]; i++ ) d1 += exp(-i*i);
	for ( int i = -order_[1]; i <= order_[1]; i++ ) d2 += exp(-i*i);
	radius_[0] /= order_[0] + sqrt( 1.0 - log(d1) );
	radius_[0] /= order_[0] + sqrt( 1.0 - log(d2) );
	a0_ /= d1  * d2;
      }
  }

  /*显示外部字段的存储值。*/
  void ExtField::show ()
  {
    if ( type_ == EMWAVE )
      {
	if            (seedType_ == PLANEWAVE)      printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = plane-wave"));
	else if       (seedType_ == GAUSSIANBEAM)   printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = gaussian-beam"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam position [")
	+ stringify(position_[0]) + std::string("; ")
	+ stringify(position_[1]) + std::string("; ")
	+ stringify(position_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam direction [")
	+ stringify(direction_[0]) + std::string("; ")
	+ stringify(direction_[1]) + std::string("; ")
	+ stringify(direction_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam normalized amplitude = ") + stringify(a0_));
	if      ( seedType_ == PLANEWAVE || seedType_ == GAUSSIANBEAM )
	  printmessage(std::string(__FILE__), __LINE__, std::string("Beam polarization [")
	+ stringify(polarization_[0]) + std::string("; ")
	+ stringify(polarization_[1]) + std::string("; ")
	+ stringify(polarization_[2]) + std::string("]"));
	if      ( seedType_ == GAUSSIANBEAM )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam parallel Rayleigh-radius = ") + stringify(radius_[0]));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam perpendicular Rayleigh-radius = ") + stringify(radius_[1]));
	  }
	signal_.show();
      }
  }

  /*** FreeElectronLaser **************************************************************************/

  /*设置辐射功率的采样类型。*/
  void FreeElectronLaser::RadiationSampling::samplingType (std::string samplingType)
  {
    if      ( samplingType.compare("at-point")   == 0 )   samplingType_ = ATPOINT;
    else if ( samplingType.compare("over-line")  == 0 )   samplingType_ = OVERLINE;
    else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /*初始化用于初始化辐射功率的值。*/
  FreeElectronLaser::RadiationSampling::RadiationSampling ()
  {
    z_.clear();
    sampling_		= false;
    directory_		= "";
    basename_		= "";
    lineBegin_		= 0.0;
    lineEnd_		= 0.0;
    res_		= 0.0;
    samplingType_	= ATPOINT;
    lambda_.clear();
    lambdaMin_		= 0.0;
    lambdaMax_		= 0.0;
    lambdaRes_		= 0.0;
  }

  /*初始化探测模块初始值*/
  FreeElectronLaser::RadiationDetector::RadiationDetector()
  {
    sampling_ = false;
    zLab_ = 0.0;
    directory_ = ".";
    basename_ = "detector";
    writePowerLine_ = true;
    writeField_ = false;
    lambda_ = 1.0;
  }

  /*初始化用于初始化辐射功率的值。*/
  FreeElectronLaser::RadiationVisualization::RadiationVisualization ()
  {
    z_			= 0.0;
    sampling_		= false;
    directory_		= "";
    basename_		= "";
    rhythm_		= 0.0;
  }

  /*初始化保存撞击屏幕的粒子的值。*/
  FreeElectronLaser::ScreenProfile::ScreenProfile ()
  {
    sampling_          	= false;
    directory_    	= "./";
    basename_     	= "";
    rhythm_		= 0.0;
    pos_.clear();
  }
}
