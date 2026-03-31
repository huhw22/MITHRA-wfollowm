/********************************************************************************************************
* solver.cpp：在mithra代码中实现求解器类的函数
********************************************************************************************************/

#include <algorithm>
#include <list>
#include <sys/time.h>

#include "solver.h"
#include "beam.cc"

namespace MITHRA
{
  Solver::Solver (Mesh& 				mesh,
		  Bunch& 				bunch,
		  Seed& 				seed,
		  std::vector<Undulator>&		undulator,
		  std::vector<ExtField>& 		extField,
		  std::vector<FreeElectronLaser>& 	FEL)
  : mesh_ 	( mesh ),
    bunch_ 	( bunch ),
    seed_ 	( seed ),
    undulator_ 	( undulator ),
    extField_ 	( extField ),
    FEL_ 	( FEL )
  {
    /*清除FdTd类中的向量。*/
    anp1_ = new std::vector<FieldVector<Double> > ();
    an_   = new std::vector<FieldVector<Double> > ();
    anm1_ = new std::vector<FieldVector<Double> > ();

    fnp1_ = new std::vector<Double> ();
    fn_   = new std::vector<Double> ();
    fnm1_ = new std::vector<Double> ();

    chargeVectorn_.clear();

    /*将节点数重置为零。*/
    N0_ = N1_ = N2_ = 0;

    /*重置FdTd类的时间和时间号。*/
    timep1_	   =  0.0;
    time_  	   =  0.0;
    timem1_	   =  0.0;
    nTime_ 	   =  0;
    nTimeBunch_    =  0;

    /*初始化MPI变量值。*/
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&size_);
    rankB_ = ( rank_ == 0 ) ? size_ - 1 : rank_ - 1;
    rankF_ = ( rank_ == size_ - 1 ) ? 0 : rank_ + 1;

    /*初始化收费的MPI数据类型。*/
    MPI_Type_contiguous(11, MPI_DOUBLE, &MPI_CHARGE);
    MPI_Type_commit(&MPI_CHARGE);

    /*根据给定的长度尺度和时间尺度初始化光速值。*/
    c0_ = C0 / mesh_.lengthScale_ * mesh_.timeScale_;
    m0_ = MU_ZERO / mesh_.lengthScale_;
    e0_ = 1.0 / ( c0_ * c0_ * m0_ );
  }

  /******************************************************************************************************
  *使用解析后的数据，设置仿真所需的参数。
  ******************************************************************************************************/

  void Solver::setSimulationParameters ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Boosting the given mesh parameters into the electron rest frame ") );

    /****************************************************************************************************/

    /*在求解器的其他部分初始化相应的长度和时间尺度。*/
    seed_.c0_ 			 = c0_;
    seed_.signal_.t0_ 		/= c0_;
    seed_.signal_.f0_ 		*= c0_;
    seed_.signal_.s_          	/= c0_;

    seed_.l_ 		 	 = c0_ / seed_.signal_.f0_;
    seed_.zR_.resize(2,0.0);

    seed_.zR_[0] 		 = PI * seed_.radius_[0] * seed_.radius_[0] / seed_.l_;
    seed_.zR_[1] 		 = PI * seed_.radius_[1] * seed_.radius_[1] / seed_.l_;

    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      {
	iter->c0_ 	       	 = c0_;
	iter->signal_.t0_     	/= c0_;
	iter->signal_.f0_ 	*= c0_;
	iter->signal_.s_ 	/= c0_;

	iter->l_ 		 = c0_ / iter->signal_.f0_;
	iter->zR_.resize(2,0.0);

	iter->zR_[0] 		 = PI * iter->radius_[0] * iter->radius_[0] / iter->l_;
	iter->zR_[1] 		 = PI * iter->radius_[1] * iter->radius_[1] / iter->l_;
      }

    for (std::vector<ExtField>::iterator iter = extField_.begin(); iter != extField_.end(); iter++)
      {
	iter->c0_ 	         = c0_;
	iter->signal_.t0_     	/= c0_;
	iter->signal_.f0_ 	*= c0_;
	iter->signal_.s_ 	/= c0_;

	iter->l_ 		 = c0_ / iter->signal_.f0_;
	iter->zR_.resize(2,0.0);

	iter->zR_[0] 		 = PI * iter->radius_[0] * iter->radius_[0] / iter->l_;
	iter->zR_[1] 		 = PI * iter->radius_[1] * iter->radius_[1] / iter->l_;
      }

    /****************************************************************************************************/

	/*根据设定的长度和时间尺度参数校正种子的振幅
	*和波动。注意，种子振幅是矢量势的振幅，而
	*波动振幅是场的振幅。*/
    seed_.amplitude_ 	        = seed_.a0_ * EM * c0_ / EC;
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      iter->amplitude_ 		= iter->a0_ * EM * c0_ * 2 * PI * iter->signal_.f0_ / EC;
    for (std::vector<ExtField>::iterator iter = extField_.begin(); iter != extField_.end(); iter++)
      iter->amplitude_        	= iter->a0_ * EM * c0_ * 2 * PI * iter->signal_.f0_ / EC;

    /****************************************************************************************************/

    /*计算输入束的平均值。*/
    Double gamma = 0.0;
    for (unsigned int i = 0; i < bunch_.bunchInit_.size(); i++)
      {
	if ( bunch_.bunchInit_[i].bunchType_ == "file" )
	  computeFileGamma(bunch_.bunchInit_[i]);
	if ( bunch_.bunchInit_[i].bunchType_ == "other" )
	  printmessage(std::string(__FILE__), __LINE__, std::string("Bunch mean gamma and direction are given by an external program. " ) );
	gamma += bunch_.bunchInit_[i].initialGamma_ / bunch_.bunchInit_.size();
      }

    /*现在，根据波动器的类型确定这群的最大值和最小值
	*穿过波动器。*/
    Double gmin = 1.0e100, gmax = -1.0e100, g = 0.0;
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
	{
		if ( iter->type_ == STATIC )
		{
		g  	 = gamma / sqrt( 1.0 + iter->k_ * iter->k_ / 2.0 );
		gmin = ( gmin < g ) ? gmin : g;
		gmax = ( gmax > g ) ? gmax : g;
		}
		else if ( iter->signal_.signalType_ == FLATTOP )
		{
		g  = gamma / sqrt( 1.0 + iter->a0_ * iter->a0_ / 2.0 );
		gmin = ( gmin < g ) ? gmin : g;
		gmax = ( gmax > g ) ? gmax : g;
		}
		else
		{
		/*对于光波动器，当脉冲不是平顶的脉冲格式时，的伽马值
		电子在相互作用过程中发生变化。*/
		g  = gamma / sqrt( 1.0 + iter->a0_ * iter->a0_ / 2.0 );
		gmin = ( gmin < g ) ? gmin : g;
		g  = gamma;
		gmax = ( gmax > g ) ? gmax : g;
		}

	}

    /*将移动帧的伽马值设置为最大值和最小值的平均值。*/
    if ( mesh_.gamma_ == -1.0 )
      gamma_ = ( undulator_.size() == 0 ) ? gamma : ( gmin + gmax ) / 2.0;
    else
      gamma_ = mesh_.gamma_;

    /****************************************************************************************************/

    /*升压和初始化波动器相关参数。*/
    beta_ = sqrt( 1.0 - 1.0 / ( gamma_ * gamma_ ) );
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      if ( iter->type_ == OPTICAL ) iter->lu_ /= ( 1 + beta_ );

    /*根据起始点对波动器进行排序。*/
    std::sort(undulator_.begin(), undulator_.end(), undulatorCompare);

    /*现在移动所有的波动器模块，使第一个模块从0开始。*/
    for (std::vector<Undulator>::reverse_iterator iter = undulator_.rbegin(); iter != undulator_.rend(); iter++)
      iter->rb_ -= undulator_[0].rb_;

    /****************************************************************************************************/

    /*使用束的给定伽马来设置每个束的调制波长。*/
    for (unsigned int i = 0; i < bunch_.bunchInit_.size(); i++)
	{
		/*首先确定这个群的向量。*/
		bunch_.bunchInit_[i].initialBeta_	= sqrt( 1.0 - 1.0 / pow( bunch_.bunchInit_[i].initialGamma_ , 2 ) );
		bunch_.bunchInit_[i].betaVector_.mv( bunch_.bunchInit_[i].initialBeta_, bunch_.bunchInit_[i].initialDirection_);

		/*计算束的调制波长，作为初始聚束因子或射束
		*噪声实现。*/
		if ( undulator_.size() > 0 )
			bunch_.bunchInit_[i].lambda_		= undulator_[0].lu_ / ( 2.0 * gamma_ * gamma_ ) * bunch_.bunchInit_[i].betaVector_[2] / beta_;
		else
			bunch_.bunchInit_[i].lambda_		= 0.0;

		printmessage(std::string(__FILE__), __LINE__, std::string("Modulation wavelength of the bunch outside the undulator is set to " + stringify( bunch_.bunchInit_[i].lambda_ ) ) );
	}

    printmessage(std::string(__FILE__), __LINE__, std::string("The given mesh parameters are boosted into the electron rest frame :::") );

    /****************************************************************************************************/

    /*洛伦兹增强参数也应该被转移到种子类，以便正确地
	*计算计算域内的字段。*/
    seed_.beta_    	= beta_;
    seed_.gamma_   	= gamma_;
  }

  /******************************************************************************************************
  *将网片推入电子休息架。
  ******************************************************************************************************/

  void Solver::lorentzBoostMesh ()
  {
    /****************************************************************************************************/

    /*将网格数据提升到电子静止框架中。*/
    mesh_.meshLength_[2] 	*= gamma_;
    mesh_.meshResolution_[2] 	*= gamma_;
    mesh_.meshCenter_[2] 	*= gamma_;
    mesh_.totalTime_		/= gamma_;
    mesh_.timeShift_		/= gamma_;

    /****************************************************************************************************/

    /*根据给定的仿真设置网格参数。*/
    if ( mesh_.solver_ == NSFD )
      {
	/*调整横向网格分辨率以匹配稳定性准则。*/
	Double t = 1.0 / sqrt( pow( mesh_.meshResolution_[2] / mesh_.meshResolution_[0], 2.0 ) + pow( mesh_.meshResolution_[2] / mesh_.meshResolution_[1], 2.0 ) );
	if ( t < 1.02 )
	  {
	    mesh_.meshResolution_[0] *= 1.02 / t;
	    mesh_.meshResolution_[1] *= 1.02 / t;

	    printmessage(std::string(__FILE__), __LINE__, std::string("Transverse discretization along x is set to " + stringify(mesh_.meshResolution_[0]) ) );
	    printmessage(std::string(__FILE__), __LINE__, std::string("Transverse discretization along y is set to " + stringify(mesh_.meshResolution_[1]) ) );
	  }

	/*根据色散条件，可以得到场时间步长。*/
	mesh_.timeStep_		 = mesh_.meshResolution_[2] / c0_;
	printmessage(std::string(__FILE__), __LINE__, std::string("Time step for the field update is set to " + stringify(mesh_.timeStep_ * gamma_) ) );
      }
    else if ( mesh_.solver_ == FD )
      {
	/*根据色散条件，可以得到场时间步长。*/
	mesh_.timeStep_		 = 0.98 / ( c0_ * sqrt( 1.0 / pow(mesh_.meshResolution_[0], 2.0) + 1.0 / pow(mesh_.meshResolution_[1], 2.0) + 1.0 / pow(mesh_.meshResolution_[2], 2.0) ) );
	printmessage(std::string(__FILE__), __LINE__, std::string("Time step for the field update is set to " + stringify(mesh_.timeStep_ * gamma_) ) );
      }
  }

  /******************************************************************************************************
	*推动粒子进入电子静止框架。
	******************************************************************************************************/

  void Solver::lorentzBoostBunch ()
  {
    /****************************************************************************************************/

    /*如果给定了簇更新时间步长，则设置它，否则根据MITHRA规则设置它。*/
    bunch_.timeStep_			/= gamma_;

    /*根据给定的场时间步长调整给定的束时间步长。*/
    if (bunch_.timeStep_ == 0)
      bunch_.timeStep_ 			 = mesh_.timeStep_ ;
    else
      bunch_.timeStep_ 			 = mesh_.timeStep_ / ceil(mesh_.timeStep_ / bunch_.timeStep_);
    nUpdateBunch_    			 = mesh_.timeStep_ / bunch_.timeStep_;
    printmessage(std::string(__FILE__), __LINE__, std::string("Time step for the bunch update is set to " + stringify(bunch_.timeStep_ * gamma_) ) );

    /****************************************************************************************************/

    /*将束状采样参数提升到电子静止系。*/
    bunch_.rhythm_			/= gamma_;
    bunch_.bunchVTKRhythm_		/= gamma_;
    for (unsigned int i = 0; i < bunch_.bunchProfileTime_.size(); i++)
      bunch_.bunchProfileTime_[i] 	/= gamma_;
    bunch_.bunchProfileRhythm_		/= gamma_;

    /****************************************************************************************************/

    /*将宏观粒子的坐标提升到束静止坐标系中。在助推过程中
	* z坐标的最大值很重要，因为它需要在移动中使用
	被介绍给一群人。*/

    Double zmaxL = -1.0e100, zmaxG;
    for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
      {
	Double g  	= std::sqrt(1.0 + iterQ->gb.norm2());
	Double bz 	= iterQ->gb[2] / g;
	iterQ->rnp[2]  *= gamma_;
	iterQ->gb[2] 	= gamma_ * g * ( bz - beta_ );

	zmaxL 		= std::max( zmaxL , iterQ->rnp[2] );
      }
    MPI_Allreduce(&zmaxL, &zmaxG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    /****************************************************************************************************/

    /*在这里，我们定义了在时间上的位移，使得束端在边缘场的开始处
	*部分。对于光波动器，在输入参数中应考虑这种分离
	*其中给出了偏移量。*/

    /* 这一时间偏移操作旨在确保：在 t=0 时刻，束团内 z 坐标的最大值恰好位于距离波荡器起始端 undulator[0].dist_ 的位置处。
     * （其中，默认距离设定为：对于静态波荡器取 2 个波荡器周期，对于光学波荡器取 10 个波荡器周期；
     * 这一差异是由于这两种情况下所采用的边缘场模型格式不同所致。）												*/

    if ( undulator_.size() > 0 )
      {
	Double nl = ( undulator_[0].type_ == STATIC ) ? 2.0 : 5.0 * undulator_[0].signal_.nR_;
	if (undulator_[0].dist_ == 0.0)
	  undulator_[0].dist_ = nl * undulator_[0].lu_;
	else if (undulator_[0].dist_ < nl * undulator_[0].lu_)
	  printmessage(std::string(__FILE__), __LINE__, std::string("Warning: the undulator is set very close to the bunch, the results may be inaccurate.") );
	dt_ 		= - 1.0 / ( beta_ * undulator_[0].c0_ ) * ( zmaxG + undulator_[0].dist_ / gamma_ );
	printmessage(std::string(__FILE__), __LINE__, std::string("Initial distance from bunch head to undulator is ") + stringify(undulator_[0].dist_) );
      }

    /*对种子田也应进行同样的时间转移。*/
    seed_.dt_		= dt_;

    /* 基于上述定义，在初始时刻——即粒子进入波荡器入口（实验室坐标系下 z = 0 处）之时——该位置对应于束团静止坐标系下的 z = zmax + undulator_[0].dist_ / gamma_（此时束团时间 timeBunch = 0.0）。
     * 在 MITHRA 模拟中，我们假定粒子在抵达起始点之前是沿直线运动的。
     * 这一假设要求束团的各项属性必须与波荡器起始点处的设定值相吻合。此处所谓的“起始点”，具体是指波荡器入口位置减去波荡器距离后的那个点。				*/
    bunch_.zu_ 		= zmaxG;
    bunch_.beta_ 	= beta_;

    /****************************************************************************************************/

    for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
      {
	Double g	= std::sqrt(1.0 + iterQ->gb.norm2());
	iterQ->rnp[0]  += iterQ->gb[0] / g * ( iterQ->rnp[2] - bunch_.zu_ ) * beta_;
	iterQ->rnp[1]  += iterQ->gb[1] / g * ( iterQ->rnp[2] - bunch_.zu_ ) * beta_;
	iterQ->rnp[2]  += iterQ->gb[2] / g * ( iterQ->rnp[2] - bunch_.zu_ ) * beta_;
      }


    /****************************************************************************************************/

    /* 束团需要进行平移，以确保其在进入波荡器时位于计算域的中心。为此，需要获取束团的平均 z 坐标和 beta_z 值。*/
    if ( mesh_.optimizePosition_ && ( undulator_.size() > 0 ))
      {
	Double zL  = 0.0, zG;
	Double bzL = 0.0, bzG;
	for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
	  {
	    zL += iterQ->rnp[2];
	    bzL += iterQ->gb[2] / std::sqrt( 1 + iterQ->gb.norm2() );
	  }
	MPI_Allreduce(&zL, &zG, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&bzL, &bzG, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	unsigned int NqL = chargeVectorn_.size(), NqG = 0;
	MPI_Allreduce(&NqL, &NqG, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	zG /= NqG;
	bzG /= NqG;

	Double shift 	 = bzG * (zmaxG + undulator_[0].dist_ / gamma_ - zG) / (bzG + beta_) + zG;
	zmaxG 		-= shift;
	bunch_.zu_ 	 = zmaxG;
	dt_ 		 = - 1.0 / ( beta_ * undulator_[0].c0_ ) * ( zmaxG + undulator_[0].dist_ / gamma_ );
	seed_.dt_ 	 = dt_;

	for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
	  iterQ->rnp[2] -= shift;

	printmessage(std::string(__FILE__), __LINE__, std::string("The bunch center is shifted back by ") + stringify(shift) + std::string(" .") );
      }

    /****************************************************************************************************/

    /* 根据粒子的纵向坐标，将其分配至相应的处理器*/
    distributeParticles(chargeVectorn_);

    /* 向用户打印宏粒子的总数。	*/
    unsigned int NqL = chargeVectorn_.size(), NqG = 0;
    MPI_Reduce(&NqL,&NqG,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    printmessage(std::string(__FILE__), __LINE__, std::string("The total number of macro-particles is equal to ") + stringify(NqG) + std::string(" .") );

    /* 初始化总电荷数。 */
    Nc_ = chargeVectorn_.size();

    /****************************************************************************************************/

    /* 如果仿真的总行驶距离值非零，则修正仿真中的总时间因子。 */
    if ( mesh_.totalDist_ > 0.0 )
      {
	/* 定义必要的变量，并获取 zmin 和平均 beta_z。*/
	double Lu = 0.0;
	for (auto und = undulator_.begin(); und != undulator_.end(); und++)
	  Lu += und->lu_ * und->length_ / gamma_;
	double zEnd = mesh_.totalDist_ / gamma_;
	double zMin = 1e100;
	double bz = 0;
	for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
	  {
	    zMin = std::min(zMin, iter->rnp[2]);
	    bz += iter->gb[2] / std::sqrt(1 + iter->gb.norm2());
	  }
	MPI_Allreduce(MPI_IN_PLACE, &zMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &bz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	unsigned int Nq = chargeVectorn_.size();
	MPI_Allreduce(MPI_IN_PLACE, &Nq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	bz /= Nq;

	mesh_.totalTime_ = 1 / (c0_ * (bz + beta_)) * (zEnd - beta_ * c0_ * dt_ - zMin + bz / beta_* Lu);

	printmessage(std::string(__FILE__), __LINE__, std::string("The total time to simulate has been set to ") + stringify(mesh_.totalTime_ * gamma_) + std::string(" .") );

      }
  }

  /******************************************************************************************************
   根据纵向坐标，将粒子分配至各自的处理器。
   ******************************************************************************************************/

  void Solver::distributeParticles (std::list<Charge>& chargeVector)
  {
    /* 请注意，此函数仅重新分配 q、rnp 和 gbnp，但不包括 rnm 和 gbnm。*/
    std::vector<Double> sendCV;
    std::list<Charge>::iterator it = chargeVector.begin();
    while(it != chargeVector.end())
      {
	if ( particleInProcessor( it->rnp[2] ) )
	  it++;
	else
	  {
	    sendCV.push_back(it->q);
	    sendCV.push_back(it->rnp[0]);
	    sendCV.push_back(it->rnp[1]);
	    sendCV.push_back(it->rnp[2]);
	    sendCV.push_back(it->gb[0]);
	    sendCV.push_back(it->gb[1]);
	    sendCV.push_back(it->gb[2]);
	    it = chargeVector.erase(it);
	  }
      }

    /* 遍历处理器，并将电荷放置到相应的处理器中。 */
    for ( unsigned int ip = 0; ip < size_; ip++ )
      {
	/* 获取第 i 个处理器待发送数据的大小。 */
	int sizeSend = sendCV.size();

	/* 将大小广播给所有其他处理器。 */
	MPI_Bcast(&sizeSend, 1, MPI_INT, ip, MPI_COMM_WORLD);

	/* 初始化接收缓冲区。 */
	std::vector<Double> recvCV (sizeSend);
	if ( ip == rank_ ) recvCV = sendCV;

	/* 现在，将数据从第 i 个处理器广播给所有其他处理器。 */
	MPI_Bcast(&recvCV[0], sizeSend, MPI_DOUBLE, ip, MPI_COMM_WORLD);

	/* 现在，将所有电荷放入相应处理器的电荷向量中。 */
	unsigned int i = 0;
	Charge charge;
	while (i < recvCV.size() )
	  {
	    if ( particleInProcessor(recvCV[i+3]) )
	      {
		charge.q 	= recvCV[i++];
		charge.rnp[0] 	= recvCV[i++];
		charge.rnp[1] 	= recvCV[i++];
		charge.rnp[2] 	= recvCV[i++];
		charge.gb[0] 	= recvCV[i++];
		charge.gb[1] 	= recvCV[i++];
		charge.gb[2] 	= recvCV[i++];
		chargeVector.push_back(charge);
	      }
	    else
	      i += 7;
	  }
      }
  }

  /******************************************************************************************************
   “回收粒子”功能移除那些不再属于该处理器的粒子。
   ******************************************************************************************************/

  void Solver::recycleParticles ()
  {
    std::list<Charge>::iterator it = chargeVectorn_.begin();
    while(it != chargeVectorn_.end())
      {
	if ( particleInProcessor( it->rnp[2] ) )
	  it++;
	else
	  it = chargeVectorn_.erase(it);
      }
  }

  /******************************************************************************************************
   获取从文件中读取的一组数据的平均伽马值和平均方向。
   ******************************************************************************************************/
  void Solver::computeFileGamma 		(BunchInitialize & bunchInit)
  {
    /* 声明用于保存平均值所需的参数。 */
    Double ignore;
    FieldVector<Double> gb (0.0);
    bunchInit.initialGamma_ = 0.0;
    bunchInit.initialDirection_ = 0.0;

    /* 读取文件，对动量 gbnp 进行求和，以计算其平均值。 */
    std::ifstream myfile ( bunchInit.fileName_.c_str() );

    while (myfile.good())
      {
	/* 忽略每行前三个属于粒子位置的数值。 */
	myfile >> ignore;
	myfile >> ignore;
	myfile >> ignore;

	myfile >> gb[0];
	myfile >> gb[1];
	myfile >> gb[2];

	bunchInit.initialGamma_ += std::sqrt( 1 + gb.norm2() );
	bunchInit.initialDirection_ += gb;
      }

    bunchInit.initialGamma_ /= bunchInit.numberOfParticles_;
    bunchInit.initialDirection_ /= bunchInit.initialDirection_.norm();

    printmessage(std::string(__FILE__), __LINE__, std::string("Computed average gamma from file is " + stringify( bunchInit.initialGamma_ ) ) );
    printmessage(std::string(__FILE__), __LINE__, std::string("Computed average direction from file is " + stringify( bunchInit.initialDirection_ ) ) );


  }

  /******************************************************************************************************
   初始化用于存储场值和坐标的矩阵。
   ******************************************************************************************************/

  void Solver::initialize ()
  {
  	/* 利用已解析的数据，设置仿真的所需参数。 */
    setSimulationParameters();

    /* 作为第一步，根据 datainput 文件中给定的参数，对给定的束团进行初始化。*/
    initializeBunch();
    timeBunch_ = time_;

    /* 将所有物理量转换至电子静止系。 */
    lorentzBoostMesh();

    /* 初始化问题的空间和时间网格。 */
    initializeMesh();

    /* 将束团转换至实验室参考系，并针对 FEL 模拟对束团进行适配——即：添加散粒噪声、分布镜像宏粒子，并添加束团尾部。*/
    lorentzBoostBunch();

	/* 初始化 lab <-> box 参考系变换：
	* 约束为 tLab=0,zLab=0 <-> tBox=0,zBox=zRef */
	boostFrame_.set(gamma_, beta_, c0_, - beta_ * c0_ * dt_);
	printmessage(std::string(__FILE__), __LINE__,
                       std::string("::: Set boostFrame:  gamma :"+ stringify(gamma_ ) + " beta: " 
					+ stringify(beta_) + " c0: " + stringify(c0_) + " zRef: " + stringify(- beta_ * c0_ * dt_)) );

    /* 初始化字段的更新数据。 */
    initializeField();

    /* 如果启用了采样功能，则初始化对该字段进行采样及保存所需的数据。 */
    if (seed_.sampling_)			initializeSeedSampling();

    /* 初始化用于场可视化与保存所需的数据。 */
    initializeSeedVTK();

    /* 如果已启用性能分析，则初始化对该字段进行分析及保存所需的数据。 */
    if (seed_.profile_)				initializeSeedProfile();

    /* 初始化更新束所需的数据。 */
    initializeBunchUpdate();

	

    /* 如果启用了辐射功率的采样或可视化功能，则初始化计算及保存辐射功率所需的各项数据。*/
    initializePowerSample(); initializePowerVisualize();

	/*如果启用了固定探测面，则初始化探测面计算的各项数据*/
	initializeDetector();

    /* 如果启用了辐射能量采样，则初始化用于计算及保存辐射能量所需的各项数据。*/
    initializeEnergySample();

    /* 初始化用于保存撞击屏幕粒子的所需数据。 */
    initializeScreenProfile();

    /* 根据给定的时间偏移量，对束团和波荡器的时间进行偏移。 */
    shiftBackInTime();
  }

  /******************************************************************************************************
   初始化问题的时空网格。
   ******************************************************************************************************/

  void Solver::initializeMesh ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the temporal and spatial mesh of the problem.") );

    /* 声明计算所需的变量，以避免冗余的数据声明。 */
    unsigned int 	m = 0;

    /* 首先，需要调整网格长度，使其在各个方向上均为网格分辨率的整数倍。 */
    N0_ = (int) ( mesh_.meshLength_[0] / mesh_.meshResolution_[0] ) + 2;
    N1_ = (int) ( mesh_.meshLength_[1] / mesh_.meshResolution_[1] ) + 2;
    N2_ = (int) ( mesh_.meshLength_[2] / mesh_.meshResolution_[2] ) + 2;
    N1N0_ = N1_*N0_;
    mesh_.meshLength_[0] = ( N0_ - 1 ) * mesh_.meshResolution_[0];
    mesh_.meshLength_[1] = ( N1_ - 1 ) * mesh_.meshResolution_[1];
    mesh_.meshLength_[2] = ( N2_ - 1 ) * mesh_.meshResolution_[2];

    /* 计算每个处理器中的节点数量。 */
    if ( size_ > 1 )
      {
	if ( rank_ == 0 )
	  {
	    np_ = N2_ / size_ + 1;
	    k0_ = 0;
	  }
	else if ( rank_ == size_ - 1 )
	  {
	    np_ = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 1;
	    k0_ = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
	  }
	else
	  {
	    np_ = N2_ / size_ + 2;
	    k0_ = rank_ * ( N2_ / size_ ) - 1;
	  }
      }
    else
      {
	np_ = N2_;
	k0_ = 0;
      }

    /* 接下来，根据节点数量，初始化各场量矩阵和坐标矩阵。 */
    FieldVector<Double> ZERO_VECTOR 		(0.0);
    FieldVector<float>  ZERO_VECTOR_FLOAT 	(0.0);
    anp1_ = new std::vector<FieldVector<Double> > (N1N0_*np_, ZERO_VECTOR);
    an_   = new std::vector<FieldVector<Double> > (N1N0_*np_, ZERO_VECTOR);
    anm1_ = new std::vector<FieldVector<Double> > (N1N0_*np_, ZERO_VECTOR);
    en_  .resize(N1N0_*np_, ZERO_VECTOR_FLOAT);
    bn_  .resize(N1N0_*np_, ZERO_VECTOR_FLOAT);
    pic_ .resize(N1N0_*np_, false);

    if ( mesh_.spaceCharge_ )
      {
	fnp1_ = new std::vector<Double> (N1N0_*np_, 0.0 );
	fn_   = new std::vector<Double> (N1N0_*np_, 0.0 );
	fnm1_ = new std::vector<Double> (N1N0_*np_, 0.0 );
      }

    /* 设置计算网格的边界。 */
    xmin_ = mesh_.meshCenter_[0] - mesh_.meshLength_[0] / 2.0;
    xmax_ = mesh_.meshCenter_[0] + mesh_.meshLength_[0] / 2.0;
    ymin_ = mesh_.meshCenter_[1] - mesh_.meshLength_[1] / 2.0;
    ymax_ = mesh_.meshCenter_[1] + mesh_.meshLength_[1] / 2.0;
    zmin_ = mesh_.meshCenter_[2] - mesh_.meshLength_[2] / 2.0;
    zmax_ = mesh_.meshCenter_[2] + mesh_.meshLength_[2] / 2.0;

    /* 现在，根据节点数量和网格长度设置节点的坐标。 */
    FieldVector<Double> r (0.0);
    for (int i = 0; i < N0_; i++)
      for (int j = 0; j < N1_; j++)
	for (int k = 0; k < np_; k++)
	  {
	    m = N1N0_*k+N1_*i+j;
	    r = rc(m);

	    if 		( k == 0 ) 		
	      zp_[0] = r[2];
	    else if 	( k == np_ - ( ( rank_ == size_ - 1 ) ? 1 : 2 ) ) 	
	      zp_[1] = r[2];
	  }

    /* 初始化用于字段更新的时间值。 */
    timep1_	   =  mesh_.timeStep_;
    time_  	   =  0.0;
    timem1_	   = -mesh_.timeStep_;

    printmessage(std::string(__FILE__), __LINE__, std::string("The temporal and spatial mesh of the problem is initialized. :::") );
  }

  /******************************************************************************************************
   初始化用于FDTD算法中场更新的数据。
   ******************************************************************************************************/

  void Solver::initializeField ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(" ::: Initializing the field update data") );

    /* 声明计算所需的变量，以避免冗余的数据声明。 */
    unsigned int      m = 0;

    /* 初始化用于更新电流的数据。 */
    uc_.dx = mesh_.meshResolution_[0];
    uc_.dy = mesh_.meshResolution_[1];
    uc_.dz = mesh_.meshResolution_[2];
    uc_.dv = - m0_ * EC / mesh_.timeStep_ /   ( uc_.dx * uc_.dy * uc_.dz );
    uc_.rc = - EC / e0_ /	              ( uc_.dx * uc_.dy * uc_.dz );
    FieldVector<Double> ZERO_VECTOR (0.0);
    uc_.jt.resize(N1N0_,ZERO_VECTOR);
    if ( mesh_.spaceCharge_ ) uc_.rt.resize(N1N0_,0.0);

    /* 现在计算随时间更新场量所需的所有系数。*/

    uf_.dt	 = mesh_.timeStep_;

    uf_.dx	 = mesh_.meshResolution_[0];
    uf_.dy	 = mesh_.meshResolution_[1];
    uf_.dz	 = mesh_.meshResolution_[2];

    uf_.dx2    = 2.0 * uf_.dx;
    uf_.dy2    = 2.0 * uf_.dy;
    uf_.dz2    = 2.0 * uf_.dz;

    uf_.N0m1   = N0_ - 1;
    uf_.N1m1   = N1_ - 1;
    uf_.npm1   = np_ - 1;

    /* 非标准有限差分法的系数 */
    Double beta   = ( 1.0 + 0.02 / ( pow(uf_.dz/uf_.dx,2.0) + pow(uf_.dz/uf_.dy,2.0) ) ) / 4.0;
    Double alpha  = 1.0 - 2.0 * beta;
    uf_.af.alpha_ = alpha;
    uf_.af.beta_  = beta / alpha;

    uf_.a[0]	 = 2.0 * ( 1.0 - alpha * pow(c0_*uf_.dt/uf_.dx,2) - alpha * pow(c0_*uf_.dt/uf_.dy,2) - pow(c0_*uf_.dt/uf_.dz,2) );
    uf_.a[1]	 = pow(c0_*uf_.dt/uf_.dx,2.0);
    uf_.a[2]	 = pow(c0_*uf_.dt/uf_.dy,2.0);
    uf_.a[3]	 = pow(c0_*uf_.dt/uf_.dz,2.0) - 2.0 * ( beta * pow(c0_*uf_.dt/uf_.dx,2) + beta * pow(c0_*uf_.dt/uf_.dy,2) );
    uf_.a[4]	 = pow(c0_*uf_.dt,2.0) * uc_.dv;
    uf_.a[5]   = pow(c0_*uf_.dt,2.0) * uc_.rc;

    uf_.af.ufa_ = &uf_.a[0];

    /* x = 0 和 x = h 处的边界系数。 */

    Double alpha1 = 0.0;
    Double alpha2 = 0.0;
    Double p      = ( 1.0 + cos(alpha1) * cos(alpha2) ) / ( cos(alpha1) + cos(alpha2) );
    Double q      = - 1.0 / ( cos(alpha1) + cos(alpha2) );

    Double d   = 1.0 / ( 2.0 * uf_.dt * uf_.dx ) + p / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    uf_.bB[0]	= (   1.0 / ( 2.0 * uf_.dt * uf_.dx ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.bB[1]	= ( - 1.0 / ( 2.0 * uf_.dt * uf_.dx ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.bB[2]	= (   p   / ( c0_ * uf_.dt * uf_.dt ) + q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( uf_.dy * uf_.dy ) + c0_ / ( uf_.dz * uf_.dz ) ) ) / d;
    uf_.bB[3]	= - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dy * uf_.dy ) ) / d ;
    uf_.bB[4]	= - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dz * uf_.dz ) ) / d ;

    /* y = 0 和 y = h 处的边界系数。 */
    d  	 = 1.0 / ( 2.0 * uf_.dt * uf_.dy ) + p / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    uf_.cB[0]	= (   1.0 / ( 2.0 * uf_.dt * uf_.dy ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.cB[1]	= ( - 1.0 / ( 2.0 * uf_.dt * uf_.dy ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.cB[2]	= (   p / ( c0_ * uf_.dt * uf_.dt ) + q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( uf_.dx * uf_.dx ) + c0_ / ( uf_.dz * uf_.dz ) ) ) / d;
    uf_.cB[3]	= - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dx * uf_.dx ) ) / d ;
    uf_.cB[4]	= - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dz * uf_.dz ) ) / d ;

    /* z = 0 和 z = h 处的边界系数。 */

    alpha1 = 0.0;
    alpha2 = 0.0;
    p      = ( 1.0 + cos(alpha1) * cos(alpha2) ) / ( cos(alpha1) + cos(alpha2) );
    q      = - 1.0 / ( cos(alpha1) + cos(alpha2) );

    d  	 = 1.0 / ( 2.0 * uf_.dt * uf_.dz ) + p / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    uf_.dB[0]	= (   1.0 / ( 2.0 * uf_.dt * uf_.dz ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.dB[1]	= ( - 1.0 / ( 2.0 * uf_.dt * uf_.dz ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.dB[2]  	= (   p / ( c0_ * uf_.dt * uf_.dt ) + q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( uf_.dx * uf_.dx ) + c0_ / ( uf_.dy * uf_.dy ) ) ) / d;
    uf_.dB[3]  	= - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dx * uf_.dx ) ) / d ;
    uf_.dB[4]  	= - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dy * uf_.dy ) ) / d ;

    /* 沿 z 轴边的边系数 */
    d  	= ( 1.0 / uf_.dy + 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) + 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt );
    uf_.eE[0] = ( - ( 1.0 / uf_.dy - 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.eE[1] = (   ( 1.0 / uf_.dy - 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.eE[2] = (   ( 1.0 / uf_.dy + 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.eE[3] = ( 3.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt ) - c0_ / ( 4.0 * uf_.dz * uf_.dz ) ) / d;
    uf_.eE[4] = c0_ / ( 8.0 * uf_.dz * uf_.dz ) / d;

    /* 沿 x 轴边的边系数 */
    d  	= ( 1.0 / uf_.dz + 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) + 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt );
    uf_.fE[0] = ( - ( 1.0 / uf_.dz - 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.fE[1] = (   ( 1.0 / uf_.dz - 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.fE[2] = (   ( 1.0 / uf_.dz + 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.fE[3] = ( 3.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt ) - c0_ / ( 4.0 * uf_.dx * uf_.dx ) ) / d;
    uf_.fE[4] = c0_ / ( 8.0 * uf_.dx * uf_.dx ) / d;

    /* 沿 y 轴边的边系数 */
    d  	= ( 1.0 / uf_.dx + 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) + 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt );
    uf_.gE[0] = ( - ( 1.0 / uf_.dx - 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.gE[1] = (   ( 1.0 / uf_.dx - 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.gE[2] = (   ( 1.0 / uf_.dx + 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.gE[3] = ( 3.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt ) - c0_ / ( 4.0 * uf_.dy * uf_.dy ) ) / d;
    uf_.gE[4] = c0_ / ( 8.0 * uf_.dy * uf_.dy ) / d;

    /* 角系数 */
    uf_.hC[0]  =   ( - 1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[1]  =   (   1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[2]  =   ( - 1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[3]  =   ( - 1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[4]  =   (   1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[5]  =   (   1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[6]  =   ( - 1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[7]  =   (   1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[8]  = - ( - 1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[9]  = - (   1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[10] = - ( - 1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[11] = - ( - 1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[12] = - (   1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[13] = - (   1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[14] = - ( - 1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[15] = - (   1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[16] = 1.0 / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    /* 种子场也应在计算网格中进行初始化。此操作应针对 TF 域内的网格点执行。*/
    if (seed_.amplitude_ > 1.0e-50)
      {
	FieldVector<Double> r (0.0);
	for (int i = 2; i < N0_ - 2; i++)
	  for (int j = 2; j < N1_ - 2; j++)
	    for (int k = 2 * abs(signof(rank_) - 1 ) ; k < np_ - 2 * abs(signof(rank_ - size_ + 1) + 1 ); k++)
	      {
		m = N1N0_ * k + N1_ * i + j; r = rc(m);
		seed_.fields(r, time_,   (*an_)[m]   );
		seed_.fields(r, timem1_, (*anm1_)[m] );
	      }
      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The field update data are initialized. :::") );
  }

  /******************************************************************************************************
   初始化计算域内采样种子或整个场所需的数据。
   ******************************************************************************************************/

  void Solver::initializeSeedSampling ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(" ::: Initializing the field sampling data") );

    /* 如果种子采样节奏仍为零，则返回错误。 */
    if ( seed_.samplingRhythm_ == 0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The sampling rhythm of the field is zero although sampling is activated !!!") );
	exit(1);
      }

    /* 将种子采样节奏进行洛伦兹提升，转换至电子静止系。 */
    seed_.samplingRhythm_		/= gamma_;

    /* 对采样数据执行洛伦兹变换。 */
    for (unsigned i = 0; i < seed_.samplingPosition_.size(); i++)
      seed_.samplingPosition_[i][2] 	*= gamma_;
    seed_.samplingLineBegin_[2]	*= gamma_;
    seed_.samplingLineEnd_  [2]	*= gamma_;

    /* 如果采样类型为“沿线绘制”，则根据线的起点和终点初始化位置。 */
    if ( seed_.samplingType_ == OVERLINE )
      {
	FieldVector<Double> l = seed_.samplingLineEnd_;
	l       -= seed_.samplingLineBegin_;
	Double n = seed_.samplingRes_ ;
	l       /= n;

	FieldVector<Double> position, r;

	for (unsigned i = 0; i < n; i++)
	  {
	    position[0] = seed_.samplingLineBegin_[0] + i * l[0];
	    position[1] = seed_.samplingLineBegin_[1] + i * l[1];
	    position[2] = seed_.samplingLineBegin_[2] + i * l[2];
	    seed_.samplingPosition_.push_back(position);
	  }

      }

    /* 如果采样点与对应的处理器不匹配，则将其从采样位置中移除。*/
    std::vector<FieldVector<Double> >	samplingPosition; samplingPosition.clear();
    for (unsigned int n = 0; n < seed_.samplingPosition_.size(); ++n)
      {
	/* 获取采样点的位置。 */
	sf_.position 	= seed_.samplingPosition_[n];

	/* 检查采样点是否位于计算网格内。 */
	if ( sf_.position[0] < xmax_ - ub_.dx && sf_.position[0] > xmin_ + ub_.dx &&
	    sf_.position[1] < ymax_ - ub_.dy && sf_.position[1] > ymin_ + ub_.dy &&
	    sf_.position[2] < zmax_ - ub_.dz && sf_.position[2] > zmin_ + ub_.dz )
	  {
	    if ( sf_.position[2] < zp_[1] && sf_.position[2] >= zp_[0] )
	      samplingPosition.push_back(seed_.samplingPosition_[n]);
	  }
	else
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The sampling point does not reside in the grid. No data is saved.") );
	  }
      }
    seed_.samplingPosition_ = samplingPosition;
    sf_.N  = seed_.samplingPosition_.size();

    /* 创建用于保存种子采样数据的文件流及目录。 */
    if ( sf_.N > 0 )
      {
	std::string baseFilename = "";
	if (!(isabsolute(seed_.samplingBasename_))) baseFilename = seed_.samplingDirectory_;
	baseFilename += seed_.samplingBasename_ + "-" + stringify(rank_) + TXT_FILE_SUFFIX;

	createDirectory(baseFilename, 0);
	sf_.file = new std::ofstream(baseFilename.c_str(),std::ios::trunc);
      }

    /* 设置常量以更改输出文件中的单位。 */
    sf_.Ce = mesh_.lengthScale_ / pow( mesh_.timeScale_ , 2 );
    sf_.Cb = 1.0 / ( mesh_.lengthScale_ * mesh_.timeScale_  );
    sf_.Ca = 1.0 / mesh_.timeScale_;
    sf_.Cj = 1.0 / ( pow( mesh_.lengthScale_ , 2 ) * mesh_.timeScale_ );
    sf_.Cf = pow( mesh_.lengthScale_ / mesh_.timeScale_ , 2 );

    printmessage(std::string(__FILE__), __LINE__, std::string(" The field sampling data are initialized. :::") );
  }

  /******************************************************************************************************
   初始化用于场可视化与保存所需的数据。
   ******************************************************************************************************/

  void Solver::initializeSeedVTK ()
  {
    /* 根据已定义的可视化数量，调整可视化字段数据库的大小。 */
    vf_.resize( seed_.vtk_.size() );

    /* 对可视化数据执行洛伦兹变换。 */
    for (unsigned int i = 0; i < seed_.vtk_.size(); i++)
      {
	/* 如果此 vtk 文件的采样选项未启用，则继续循环。 */
	if ( !seed_.vtk_[i].sample_ ) continue;

	/* 如果种子可视化节奏仍为零，则返回错误。 */
	if ( seed_.vtk_[i].rhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The visualization rhythm of the field is zero although visualization is activated !!!") );
	    exit(1);
	  }

	/* 将种子可视化节奏进行洛伦兹提升，转换至电子静止系。 */
	seed_.vtk_[i].rhythm_			/= gamma_;

	/* 将平面位置进行洛伦兹提升，转换至电子静止系。 */
	seed_.vtk_[i].position_[2]		*= gamma_;

	if (!(isabsolute(seed_.vtk_[i].basename_)))
	  seed_.vtk_[i].basename_ = seed_.vtk_[i].directory_ + seed_.vtk_[i].basename_;
	splitFilename(seed_.vtk_[i].basename_, vf_[i].path, vf_[i].name);

	/* 如果 baseFilename 所在的目录不存在，则创建该目录。*/
	createDirectory(seed_.vtk_[i].basename_, rank_);

	/* 初始化用于存储数据的向量，并将其写入 VTK 文件。 */
	std::vector<Double> ZERO_VECTOR ( (seed_.vtk_[i].field_).size(), 0.0);
	if 	  ( seed_.vtk_[i].type_ == ALLDOMAIN )
	  vf_[i].v.resize( N1N0_*np_, ZERO_VECTOR);
	else if ( seed_.vtk_[i].type_ == INPLANE )
	  {
	    if      ( seed_.vtk_[i].plane_ == XNORMAL )
	      {
		/* 检查给定平面的位置是否位于计算域内。*/
		if ( seed_.vtk_[i].position_[0] > xmax_ - ub_.dx || seed_.vtk_[i].position_[0] < xmin_ + ub_.dx )
		  {
		    printmessage(std::string(__FILE__), __LINE__, std::string("The plane does not reside in the grid. No data is saved.") );
		    seed_.vtk_.erase( seed_.vtk_.begin() + i );
		    vf_.erase( vf_.begin() + i );
		  }
		else
		  vf_[i].v.resize( N1_*np_, ZERO_VECTOR);
	      }
	    else if ( seed_.vtk_[i].plane_ == YNORMAL )
	      {
		/* 检查给定平面的位置是否位于计算域内。*/
		if ( seed_.vtk_[i].position_[1] > ymax_ - ub_.dy || seed_.vtk_[i].position_[1] < ymin_ + ub_.dy )
		  {
		    printmessage(std::string(__FILE__), __LINE__, std::string("The plane does not reside in the grid. No data is saved.") );
		    seed_.vtk_.erase( seed_.vtk_.begin() + i );
		    vf_.erase( vf_.begin() + i );
		  }
		else
		  vf_[i].v.resize( N0_*np_, ZERO_VECTOR);
	      }
	    else if ( seed_.vtk_[i].plane_ == ZNORMAL )
	      {
		/* 检查给定平面的位置是否位于计算域内。*/
		if ( seed_.vtk_[i].position_[2] > zmax_ - ub_.dz || seed_.vtk_[i].position_[2] < zmin_ + ub_.dz )
		  {
		    printmessage(std::string(__FILE__), __LINE__, std::string("The plane does not reside in the grid. No data is saved.") );
		    seed_.vtk_.erase( seed_.vtk_.begin() + i );
		    vf_.erase( vf_.begin() + i );
		  }
		else
		  {
		    if ( seed_.vtk_[i].position_[2] < zp_[1] && seed_.vtk_[i].position_[2] >= zp_[0] )
		      vf_[i].v.resize( N1N0_,   ZERO_VECTOR);
		  }
	      }
	  }
      }
  }

  /******************************************************************************************************
   初始化用于性能分析及字段保存的所需数据。
   ******************************************************************************************************/

  void Solver::initializeSeedProfile ()
  {
    /* 如果批量分析的节奏仍为零且未设置时间，则返回错误。 */
    if ( seed_.profileRhythm_ == 0 && seed_.profileTime_.size() == 0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The profiling rhythm of the field is zero and no time is set although profiling of the field is activated !!!") );
	exit(1);
      }

    /* 将种子剖面节奏进行洛伦兹提升，转换至电子静止系。 */
    seed_.profileRhythm_	/= gamma_;

    /* 对剖析数据执行洛伦兹变换。 */
    for (unsigned i = 0; i < seed_.profileTime_.size(); i++)
      seed_.profileTime_[i] 	/= gamma_;

    if (!(isabsolute(seed_.profileBasename_))) seed_.profileBasename_ = seed_.profileDirectory_ + seed_.profileBasename_;

    /* 如果 baseFilename 所在的目录不存在，则创建该目录。*/
    createDirectory(seed_.profileBasename_, rank_);

    pf_.dt = 2.0 * mesh_.timeStep_;
  }

  /******************************************************************************************************
   初始化更新批次所需的必要数据。
   ******************************************************************************************************/

  void Solver::initializeBunchUpdate ()
  {
    /* 计算时间步长和网格分辨率的绝对值。 */
    ub_.dt	= mesh_.timeStep_;
    ub_.dtb	= c0_ * bunch_.timeStep_;
    ub_.dx	= mesh_.meshResolution_[0];
    ub_.dy	= mesh_.meshResolution_[1];
    ub_.dz	= mesh_.meshResolution_[2];
    ub_.r1    	= - EC / ( EM * c0_ ) * bunch_.timeStep_ / 2.0;
    ub_.r2    	= - EC / EM * bunch_.timeStep_ / 2.0;

    /* 如果启用了成组采样，则初始化用于采样及保存成组数据所需的各项数据。 */
    if (bunch_.sampling_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the bunch sampling data.") );

	std::string baseFilename = "";
	if (!(isabsolute(bunch_.basename_))) baseFilename = bunch_.directory_;
	baseFilename += bunch_.basename_ + TXT_FILE_SUFFIX;

	/* 如果 baseFilename 所在的目录不存在，则创建该目录。*/
	createDirectory(baseFilename, rank_);

	sb_.file = new std::ofstream(baseFilename.c_str(),std::ios::trunc);

	/* 如果批量采样节奏仍为零，则返回错误。 */
	if ( bunch_.rhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The sampling rhythm of the bunch is zero although sampling is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The sampling data are initialized. :::") );
      }

    /* 如果启用了束可视化功能，则初始化用于可视化该束所需的必要数据。*/
    if (bunch_.bunchVTK_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the bunch visualization data.") );

	if (!(isabsolute(bunch_.bunchVTKBasename_))) bunch_.bunchVTKBasename_ = bunch_.bunchVTKDirectory_ + bunch_.bunchVTKBasename_;

	/* 如果 baseFilename 所在的目录不存在，则创建该目录。*/
	createDirectory(bunch_.bunchVTKBasename_, rank_);

	/* 如果簇可视化节奏仍为零，则返回错误。 */
	if ( bunch_.bunchVTKRhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The visualization rhythm of the bunch is zero although visualization is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The bunch visualization data are initialized. :::") );
      }

    /* 如果启用了束剖面写入功能，则初始化用于束剖面分析所需的各项数据。 */
    if (bunch_.bunchProfile_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the bunch profiling data.") );

	if (!(isabsolute(bunch_.bunchProfileBasename_))) bunch_.bunchProfileBasename_ = bunch_.bunchProfileDirectory_ + bunch_.bunchProfileBasename_;

	/* 如果 baseFilename 所在的目录不存在，则创建该目录。*/
	createDirectory(bunch_.bunchProfileBasename_, rank_);

	/* 如果批量分析的节奏仍为零且未设置时间，则返回错误。 */
	if ( bunch_.bunchProfileRhythm_ == 0 && bunch_.bunchProfileTime_.size() == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The profiling rhythm of the bunch is zero and no time is set although profiling of the bunch is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The bunch profiling data are initialized. :::") );
      }
  }

  /******************************************************************************************************
   初始化包含用户指定粒子束的电荷矢量。
   ******************************************************************************************************/

  void Solver::initializeBunch ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("[[[ Initializing the bunch and prepare the charge vector ") );

    /* 定义一个临时电荷向量。 */
    std::list<Charge> qv;

    for (unsigned int i = 0; i < bunch_.bunchInit_.size(); i++)
      {
	/* 清除临时电荷向量。 */
	qv.clear();

	if ( bunch_.bunchInit_[i].position_.size() == 0 )
	  {
	    bunch_.bunchInit_[i].position_.push_back( FieldVector<Double>(0.0) );
	    printmessage(std::string(__FILE__), __LINE__, std::string("No bunch position was given, using default " + stringify(bunch_.bunchInit_[i].position_[0]) ) );
	  }      


	/* 在代码中初始化该组对象。 */
	if 	  ( bunch_.bunchInit_[i].bunchType_ == "manual" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initializeManual(	bunch_.bunchInit_[i], qv, zp_, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "ellipsoid" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initializeEllipsoid( bunch_.bunchInit_[i], qv, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "3D-crystal" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initialize3DCrystal(	bunch_.bunchInit_[i], qv, zp_, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "file" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initializeFile(		bunch_.bunchInit_[i], qv, zp_, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "other" )
	  printmessage(std::string(__FILE__), __LINE__, std::string("The charge vector has been filled in by an external program. ") );

	/* 将束团分布添加到全局电荷矢量中。 */
	chargeVectorn_.splice(chargeVectorn_.end(),qv);
      }

    printmessage(std::string(__FILE__), __LINE__, std::string("The bunch is initialized and the charge vector is prepared. ]]]") );
  }

  /******************************************************************************************************
   数据初始化完成后，若时间偏移量非零，则将束团和波荡器在时间上向后平移。
   ******************************************************************************************************/
  void Solver::shiftBackInTime()
  {
    /* 首先检查时间偏移是否非零。 */
    if ( mesh_.timeShift_ != 0.0 )
      {
	/* 若要将波荡器在时间上向后推移，只需调整 dt_ 因子即可。 */
	timem1_ 	-= mesh_.timeShift_;
	time_   	-= mesh_.timeShift_;
	timep1_ 	-= mesh_.timeShift_;
	timeBunch_ 	-= mesh_.timeShift_;

	/* 若要将束团在时间上向后平移，应调整 rnm 和 rnp 的数值。 */
	for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
	  {
	    Double t  	= c0_ * mesh_.timeShift_ / std::sqrt(1.0 + iterQ->gb.norm2());
	    iterQ->rnp.mmv( t , iterQ->gb );
	  }

	/* 根据粒子的纵向坐标，将其分配至各自的处理器。 */
	distributeParticles(chargeVectorn_);
      }
  }

  /******************************************************************************************************
   用于求解时域场量的函数。
   ******************************************************************************************************/

  void Solver::solve ()
  {
    /* 声明计算所需的变量，以避免冗余的数据声明。 */
    timeval           			simulationStart, simulationEnd;
    Double 				deltaTime, p = 0.0;
    std::stringstream 			printedMessage;
    std::list<Charge>::iterator 	iter;

    /* 在开始进行仿真之前，必须对用于存储场值及坐标的矩阵进行初始化。
     * 这一过程将依据网格结构中给定的网格长度及网格分辨率来完成。											*/
    initialize();

    /* 获取开始执行电气更新部分的时间 */
    gettimeofday(&simulationStart, NULL);

    /* 现在，从零时刻开始，对整个仿真时段内的场进行更新。 */
    printmessage(std::string(__FILE__), __LINE__, std::string("-> Run the time domain inital particle simulation ...") );

    /* 首先运行求解器，计算粒子运动直至初始时刻。 */
    while (time_ < 0.0)
      {
	/* 更新位置和速度参数。 */
	for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
	  iter->rnm  = iter->rnp;

	/* 更新粒子束，直至其属性时间达到场的时间点。 */
	for (Double t = 0.0; t < nUpdateBunch_; t += 1.0)
	  {
	    bunchUpdate();
	    timeBunch_ += bunch_.timeStep_;
	    ++nTimeBunch_;
	  }

	/* 记录已穿过诊断屏的粒子。 */
	screenProfile();

	/* 如果已启用对束的采样，且已达到采样的节奏，则对该束进行采样，
	 * 并将其保存至文件中。													*/
	if ( bunch_.sampling_ && fmod(time_ + mesh_.timeShift_ , bunch_.rhythm_) < mesh_.timeStep_ && ( time_ + mesh_.timeShift_ > 0.0 ) ) bunchSample();

	/* 如果已启用束团可视化功能，且已达到可视化所需的节奏，
	 * 则对束团进行可视化，并将 VTK 数据保存至指定的文件名中。				*/
	if ( bunch_.bunchVTK_ && fmod(time_ + mesh_.timeShift_ , bunch_.bunchVTKRhythm_) < mesh_.timeStep_ && ( time_ + mesh_.timeShift_ > 0.0 ) ) bunchVisualize();

	/* 如果已启用批次性能分析功能，且已达到分析时间点，则写入该批次的性能分析数据，
	 * 并将其保存至指定的文件中。						*/
	if (bunch_.bunchProfile_)
	  {
	    for (unsigned int i = 0; i < (bunch_.bunchProfileTime_).size(); i++)
	      if ( time_ - bunch_.bunchProfileTime_[i] < mesh_.timeStep_ && time_ > bunch_.bunchProfileTime_[i] )
		bunchProfile();
	    if ( fmod(time_ + mesh_.timeShift_ , bunch_.bunchProfileRhythm_) < mesh_.timeStep_ && ( time_ + mesh_.timeShift_ > 0.0 ) && ( bunch_.bunchProfileRhythm_ != 0.0 ) )
	      bunchProfile();
	  }

	/* 在当前沉积过程结束后，需要对电荷进行回收，以便将已离开处理器域的电荷从处理器中移除。*/
	recycleParticles();

	timem1_ += mesh_.timeStep_;
	time_   += mesh_.timeStep_;
	timep1_ += mesh_.timeStep_;
	++nTime_;

	gettimeofday(&simulationEnd, NULL);
	deltaTime  = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
	deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );

	if ( rank_ == 0 && ( int(time_/mesh_.timeShift_ * 100.0) !=  int(timem1_/mesh_.timeShift_ * 100.0) ) )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Percentage of the initial simulation completed (%)      = ") +
			 stringify( ( 1.0 - fabs( time_ / mesh_.timeShift_ ) ) * 100.0 ) );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Average calculation time for each time step (s) = ") +
			 stringify(deltaTime/(double)(nTime_))     );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Estimated remaining time of the initial simulation (min)                  = ") +
			 stringify( fabs( time_ / mesh_.timeShift_ ) / ( 1.0 - fabs( time_ / mesh_.timeShift_ ) ) * deltaTime / 60 ) );
	  }
      }

    /* 获取开始执行电气更新部分的时间 */
    gettimeofday(&simulationStart, NULL);

    /* 现在，从零时刻开始，对整个仿真时段内的场进行更新。 */
    printmessage(std::string(__FILE__), __LINE__, std::string("-> Run the time domain field simulation ...") );

    /* 现在，执行完整的辐射计算，直至最终时刻。 */
    while (time_ < mesh_.totalTime_)
      {
	/* 利用 FDTD 算法更新一个时间步长的场量。 */
	fieldUpdate();

	/* 为了确保 PIC 模型实现中电荷守恒的正确性，
	 * 我们首先需要利用内存中保存的初始位置来更新电荷的运动状态。
	 * 随后，电流及电荷的更新操作必须完全基于电荷的“初始位置”与“最终位置”来进行。
	 * 这一点至关重要，在未来的版本中绝不应予以更改。		*/

	/* 更新位置和速度参数。 */
	for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
	  iter->rnm  = iter->rnp;

	/* 更新粒子束，直至其属性时间达到场的时间点。 */
	for (Double t = 0.0; t < nUpdateBunch_; t += 1.0)
	  {
	    bunchUpdate();
	    timeBunch_ += bunch_.timeStep_;
	    ++nTimeBunch_;
	  }

	/* 在当前沉积过程结束后，需要对电荷进行回收，以便将已离开处理器域的电荷从处理器中移除。*/
	recycleParticles();

	/* 如果已启用该字段的采样功能，且已达到预设的采样节奏，则在指定位置对该字段进行采样，
	 * 并将采样数据保存至文件中。													*/
	if ( seed_.sampling_ && fmod(time_, seed_.samplingRhythm_) < mesh_.timeStep_ && time_ > 0.0 ) fieldSample();

	/* 如果已启用场的可视化功能，且已达到可视化所需的节奏，
	 * 则对场进行可视化，并将 VTK 数据保存至指定的文件名。				*/
	for (unsigned int i = 0; i < seed_.vtk_.size(); i++)
	  {
	    if ( seed_.vtk_[i].sample_ && fmod(time_, seed_.vtk_[i].rhythm_) < mesh_.timeStep_ && time_ > 0.0 )
	      {
		if 		( seed_.vtk_[i].type_ == ALLDOMAIN ) fieldVisualizeAllDomain(i);
		else if 	( seed_.vtk_[i].type_ == INPLANE   ) fieldVisualizeInPlane(i);
	      }
	  }

	/* 如果已启用该字段的性能分析功能，且已达到分析时间，则写入字段分析数据，
	 * 并将其保存至指定的文件中。													*/
	if (seed_.profile_)
	  {
	    for (unsigned int i = 0; i < seed_.profileTime_.size(); i++)
	      if ( time_ - seed_.profileTime_[i] < mesh_.timeStep_ && time_ > seed_.profileTime_[i] )
		fieldProfile();

	    if ( fmod(time_, seed_.profileRhythm_) < mesh_.timeStep_ && time_ > 0.0 && seed_.profileRhythm_ != 0 )
	      fieldProfile();
	  }

	/* 如果已启用对束的采样，且已达到采样的节奏，则对该束进行采样，
	 * 并将其保存至文件中。													*/
	if ( bunch_.sampling_ && fmod(time_ + mesh_.timeShift_ , bunch_.rhythm_) < mesh_.timeStep_ && ( time_ + mesh_.timeShift_ > 0.0 ) ) bunchSample();

	/* 如果已启用束团可视化功能，且已达到可视化所需的节奏，
	 * 则对束团进行可视化，并将 VTK 数据保存至指定的文件名中。				*/
	if ( bunch_.bunchVTK_ && fmod(time_ + mesh_.timeShift_ , bunch_.bunchVTKRhythm_) < mesh_.timeStep_ && ( time_ + mesh_.timeShift_ > 0.0 ) ) bunchVisualize();

	/* 如果已启用批次性能分析功能，且已达到分析时间点，则写入该批次的性能分析数据，
	 * 并将其保存至指定的文件中。						*/
	if (bunch_.bunchProfile_)
	  {
	    for (unsigned int i = 0; i < (bunch_.bunchProfileTime_).size(); i++)
	      if ( time_ - bunch_.bunchProfileTime_[i] < mesh_.timeStep_ && time_ > bunch_.bunchProfileTime_[i] )
		bunchProfile();
	    if ( fmod(time_ + mesh_.timeShift_ , bunch_.bunchProfileRhythm_) < mesh_.timeStep_ && ( time_ + mesh_.timeShift_ > 0.0 ) && ( bunch_.bunchProfileRhythm_ != 0.0 ) )
	      bunchProfile();
	  }

	/* 记录已穿过诊断屏的粒子。 */
	screenProfile();

	/* 若 FEL 输出的辐射功率监测已启用，且采样周期条件已满足：
	 * 在指定位置对辐射功率进行采样，并将其保存至文件中。		*/
	powerSample(); powerVisualize();

	/* 若 FEL 针对探测面的监测已启用，且到达探测面范围
	 * 针对探测面进行取出，然后把数据保存在文件里 */
	detectorSample();

	/* 若 FEL 输出的辐射能量监测已启用，且已达到采样节拍，
	 * 则在指定位置对辐射能量进行采样，并将其保存至文件中。		*/
	energySample();

	/* 移动计算字段及其对应的时间点。 */
	fieldShift();

	/* 将电荷和电流值重置为零。 */
	currentReset();

	/* 更新当前对象的值。 */
	currentUpdate();

	/* 在处理器之间传递电流。 */
	currentCommunicate();

	timem1_ += mesh_.timeStep_;
	time_   += mesh_.timeStep_;
	timep1_ += mesh_.timeStep_;
	++nTime_;

	gettimeofday(&simulationEnd, NULL);
	deltaTime  = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
	deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );

	if ( rank_ == 0 && ( int(time_/mesh_.totalTime_ * 1000.0) !=  int(timem1_/mesh_.totalTime_ * 1000.0) ) )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Percentage of the total simulation completed (%)      = ") +
			 stringify( ( time_ ) / ( mesh_.totalTime_ ) * 100.0 ) );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Average calculation time for each time step (s) = ") +
			 stringify(deltaTime/(double)(nTime_))     );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Estimated remaining time (min)                  = ") +
			 stringify( ( mesh_.totalTime_ / time_ - 1) * deltaTime / 60 ) );
	  }
      }

    /* 完成计算及数据保存。 */
    finalize();
  }

  /******************************************************************************************************
   更新单个时间步的场
   ******************************************************************************************************/

  void Solver::bunchUpdate ()
  {
    /* 首先定义一个用于指定处理器编号的参数。 */
    UpdateBunchParallel	        ubp;
    MPI_Status                  status;
    int                         msgtag1 = 1, msgtag2 = 2;

    /* 将位于束团横向域之外的粒子计数设为零。 */
    ubp.nt = 0;

    /* 遍历束团中的电荷点，提取各点处的种子场实数值，
     * 与波荡器场进行叠加，并最终在场内对粒子进行加速。	*/

    for ( auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++ )
      {
	/* 如果粒子不属于当前处理器，则继续粒子的循环 */
	ubp.zr = pmod( iter->rnp[2] - zmin_ , mesh_.meshLength_[2] ) + zmin_;
	if ( ! ( ( ubp.zr >= zp_[0] ) && ( ubp.zr < zp_[1] ) ) ) continue;

	/* 获取布尔标志，用于判断粒子是否位于计算域内。 */
	ubp.b1x = ( iter->rnp[0] < xmax_ - ub_.dx && iter->rnp[0] > xmin_ + ub_.dx );
	ubp.b1y = ( iter->rnp[1] < ymax_ - ub_.dy && iter->rnp[1] > ymin_ + ub_.dy );
	ubp.b1z = ( iter->rnp[2] < zp_[1]         && iter->rnp[2] >= zp_[0] );

	/* 初始化波荡器场。 */
	ubp.bt = 0.0;
	ubp.et = 0.0;

	/* 计算粒子位置处的波荡器场。 */
	undulatorField(ubp, iter->rnp);

	/* 计算粒子位置处的外部场，并将其叠加到波荡器场上。 */
	externalField(ubp, iter->rnp);

	/* 若粒子已穿过波荡器入口区，则计算场量。 */
	if ( iter->e == 1.0 )
	  {

	    if ( ubp.b1x && ubp.b1y && ubp.b1z )
	      {
		ubp.dxr = modf( ( iter->rnp[0] - xmin_ ) / ub_.dx , &ubp.d1 ); ubp.i = (int) ubp.d1;
		ubp.dyr = modf( ( iter->rnp[1] - ymin_ ) / ub_.dy , &ubp.d1 ); ubp.j = (int) ubp.d1;
		ubp.dzr = modf( ( iter->rnp[2] - zmin_ ) / ub_.dz , &ubp.d1 ); ubp.k = (int) ubp.d1;

		/* 获取单元格的索引。 */
		ubp.m   = ( ubp.k - k0_) * N1N0_ + ubp.i * N1_ + ubp.j;

		if (!pic_[ubp.m            ])     fieldEvaluate(ubp.m            );
		if (!pic_[ubp.m+N1_        ])     fieldEvaluate(ubp.m+N1_        );
		if (!pic_[ubp.m+1          ])     fieldEvaluate(ubp.m+1          );
		if (!pic_[ubp.m+N1_+1      ])     fieldEvaluate(ubp.m+N1_+1      );
		if (!pic_[ubp.m+N1N0_      ])     fieldEvaluate(ubp.m+N1N0_      );
		if (!pic_[ubp.m+N1N0_+N1_  ])     fieldEvaluate(ubp.m+N1N0_+N1_  );
		if (!pic_[ubp.m+N1N0_+1    ])     fieldEvaluate(ubp.m+N1N0_+1    );
		if (!pic_[ubp.m+N1N0_+N1_+1])     fieldEvaluate(ubp.m+N1N0_+N1_+1);

		/* 计算并插值电场，以求得束团点处的数值。 */
		ubp.et.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , en_[ubp.m		]);
		ubp.et.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , en_[ubp.m+N1_	  	]);
		ubp.et.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   * ( 1.0 - ubp.dzr) , en_[ubp.m+1	  	]);
		ubp.et.pmv(         ubp.dxr   *         ubp.dyr   * ( 1.0 - ubp.dzr) , en_[ubp.m+N1_+1	  	]);
		ubp.et.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) *         ubp.dzr  , en_[ubp.m+N1N0_	  	]);
		ubp.et.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) *         ubp.dzr  , en_[ubp.m+N1N0_+N1_  	]);
		ubp.et.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   *         ubp.dzr  , en_[ubp.m+N1N0_+1	]);
		ubp.et.pmv(         ubp.dxr   *         ubp.dyr   *         ubp.dzr  , en_[ubp.m+N1N0_+N1_+1	]);

		/* 计算并插值磁场，以求得束团点处的数值。*/
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , bn_[ubp.m		]);
		ubp.bt.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , bn_[ubp.m+N1_		]);
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   * ( 1.0 - ubp.dzr) , bn_[ubp.m+1		]);
		ubp.bt.pmv(         ubp.dxr   *         ubp.dyr   * ( 1.0 - ubp.dzr) , bn_[ubp.m+N1_+1		]);
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) *         ubp.dzr  , bn_[ubp.m+N1N0_		]);
		ubp.bt.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) *         ubp.dzr  , bn_[ubp.m+N1N0_+N1_	]);
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   *         ubp.dzr  , bn_[ubp.m+N1N0_+1	]);
		ubp.bt.pmv(         ubp.dxr   *         ubp.dyr   *         ubp.dzr  , bn_[ubp.m+N1N0_+N1_+1	]);
	      }
	    else if ( !(ubp.b1x) && !(ubp.b1y) && ubp.b1z )
	      {
		/* 将位于计算域横向尺寸范围之外的粒子数量加一。 */
		ubp.nt++;
	      }
	  }
	else if ( undulator_.size() > 0 )
	  {
	    /* 根据粒子在实验室坐标系中的位置，更新发射矢量标志。 */
	    ubp.lz 	= gamma_ * ( iter->rnp[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) );
	    iter->e 	= ( ubp.lz > - undulator_[0].dist_ ) ? 1.0 : 0.0;
	  }
	else
	  iter->e	= 1.0;

	/* 根据计算出的电场和磁场，更新粒子的速度。 */

	/* 首先，计算 (gamma*beta)- 的值。 */
	ubp.gbm = iter->gb;
	ubp.gbm.pmv( ub_.r1 , ubp.et );

	/* 其次，计算 (gamma*beta)' 的值。 */
	ubp.gbp   = cross( ubp.gbm , ubp.bt );
	ubp.d1    = sqrt( 1.0 + ubp.gbm.norm2() );
	ubp.gbp.mv( ub_.r2 / ubp.d1, ubp.gbp);
	ubp.gbp  += ubp.gbm;

	/* 第三步，计算 (gamma*beta)+。 */
	ubp.gbpl  = cross( ubp.gbp , ubp.bt );
	ubp.gbpl.mv( 2.0 / ( ubp.d1 / ub_.r2 + ub_.r2 / ubp.d1 * ubp.bt.norm2() ), ubp.gbpl);
	ubp.gbpl += ubp.gbm;

	/* 第四步，更新 (gamma*beta) 向量。 */
	iter->gb = ubp.gbpl;
	iter->gb.pmv( ub_.r1 , ubp.et );

	/* 确定粒子的运动。 */
	ubp.dr.mv( ub_.dtb / sqrt (1.0 + iter->gb.norm2()) , iter->gb );

	/* 确定粒子的最终位置。 */
	iter->rnp += ubp.dr;

	/* 计算用于处理器关联的相对坐标。 */
	ubp.zr += ubp.dr[2];

	/* 若粒子进入相邻计算域，将其存入通信缓冲区。 */
	if 	( ubp.zr <  zp_[0] )	ubp.qSB.push_back( *iter );
	else if ( ubp.zr >= zp_[1] )	ubp.qSF.push_back( *iter );
      }

    /* 现在，将那些跨越边界传播的电荷发送给其他处理器。 */
    MPI_Send(&ubp.qSB[0],ubp.qSB.size(),MPI_CHARGE,rankB_,msgtag1,MPI_COMM_WORLD);

    MPI_Probe(rankF_,msgtag1,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_CHARGE,&ub_.nL);
    ubp.qRF.resize(ub_.nL);
    MPI_Recv(&ubp.qRF[0],ub_.nL,MPI_CHARGE,rankF_,msgtag1,MPI_COMM_WORLD,&status);

    MPI_Send(&ubp.qSF[0],ubp.qSF.size(),MPI_CHARGE,rankF_,msgtag2,MPI_COMM_WORLD);

    MPI_Probe(rankB_,msgtag2,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_CHARGE,&ub_.nL);
    ubp.qRB.resize(ub_.nL);
    MPI_Recv(&ubp.qRB[0],ub_.nL,MPI_CHARGE,rankB_,msgtag2,MPI_COMM_WORLD,&status);

    /* 现在，将进入当前处理器的粒子插入到粒子列表中。 */
    std::copy( ubp.qRF.begin(), ubp.qRF.end(), std::back_inserter(chargeVectorn_) );
    std::copy( ubp.qRB.begin(), ubp.qRB.end(), std::back_inserter(chargeVectorn_) );

    /* 添加来自不同处理器、且位于域横向范围之外的电荷点数量。 */
    MPI_Reduce(&ubp.nt, &ub_.nt,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	if ( rank_ == 0 )
	{
		static int last_reported_nt = 0;

	if ( ub_.nt > last_reported_nt )
	{
		printmessage(std::string(__FILE__), __LINE__,
					std::string("Warning: " + stringify(ub_.nt) +
					" particles have transverse dimensions that are larger than the computational domain size "
					"(+" + stringify(ub_.nt - last_reported_nt) + " newly added).") );

		last_reported_nt = ub_.nt;
	}
	}
  }

  /******************************************************************************************************
   对批次数据进行采样，并将其保存到指定文件中。
   ******************************************************************************************************/

  void Solver::bunchSample ()
  {
    /* 首先，我们需要评估电荷云的属性；为此，应初始化相应的数据。 */
    sb_.q   = 0.0;
    sb_.r   = 0.0;
    sb_.r2  = 0.0;
    sb_.gb  = 0.0;
    sb_.gb2 = 0.0;

    /* 现在，我们对所有电荷进行求和。由于这些点电荷大小相等，因此无需进行加权求和。 */
    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( particleInProcessor(iter->rnp[2]) )
	  {
	    sb_.q	 += iter->q;
	    sb_.r .pmv(   iter->q, iter->rnp );
	    sb_.gb.pmv(   iter->q, iter->gb);

	    for (int l = 0; l < 3; l++)
	      {
		sb_.r2 [l] += iter->rnp[l]  * iter->rnp[l]  * iter->q;
		sb_.gb2[l] += iter->gb[l] * iter->gb[l] * iter->q;
	      }
	  }
      }

    /* 将每个处理器的贡献累加到第一个元素上。 */
    MPI_Reduce(&sb_.q    , &sb_.qT    , 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.r[0] , &sb_.rT[0] , 3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.r2[0], &sb_.r2T[0], 3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.gb[0], &sb_.gbT[0], 3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.gb2[0],&sb_.gb2T[0],3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if ( rank_ == 0 )
      {
	/* 将所得数值除以粒子数，以计算平均值。 */
	sb_.rT   /= sb_.qT;
	sb_.r2T  /= sb_.qT;
	sb_.gbT  /= sb_.qT;
	sb_.gb2T /= sb_.qT;

	(*sb_.file).setf(std::ios::scientific);
	(*sb_.file).precision(4);

	/** 将时间写入第一列。 */
	*sb_.file << timeBunch_ << "\t";

	/** 现在，写入本行电荷分布的计算值。 */
	*sb_.file << sb_.rT[0]   << "\t" << sb_.rT[1]   << "\t" << sb_.rT[2]   << "\t";
	*sb_.file << sb_.gbT[0]  << "\t" << sb_.gbT[1]  << "\t" << sb_.gbT[2]  << "\t";
	*sb_.file << sqrt( sb_.r2T[0]  - sb_.rT[0]  * sb_.rT[0]  ) << "\t" ;
	*sb_.file << sqrt( sb_.r2T[1]  - sb_.rT[1]  * sb_.rT[1]  ) << "\t" ;
	*sb_.file << sqrt( sb_.r2T[2]  - sb_.rT[2]  * sb_.rT[2]  ) << "\t" ;
	*sb_.file << sqrt( sb_.gb2T[0] - sb_.gbT[0] * sb_.gbT[0] ) << "\t" ;
	*sb_.file << sqrt( sb_.gb2T[1] - sb_.gbT[1] * sb_.gbT[1] ) << "\t" ;
	*sb_.file << sqrt( sb_.gb2T[2] - sb_.gbT[2] * sb_.gbT[2] ) << std::endl ;
      }
  }

  /******************************************************************************************************
   将粒子束可视化为 VTK 文件，并以指定名称保存到文件中。
   ******************************************************************************************************/

  void Solver::bunchVisualize ()
  {
    /* 首先定义参数。 */
    Double                            gamma, beta;

    /* 如果旧文件存在，应将其删除。 */
    vb_.fileName = bunch_.bunchVTKBasename_ + "-p" + stringify(rank_) + "-" + stringify(nTimeBunch_) + VTU_FILE_SUFFIX;
    vb_.file = new std::ofstream(vb_.fileName.c_str(),std::ios::trunc);

    vb_.file->setf(std::ios::scientific);
    vb_.file->precision(4);

    /* 存储模拟中的粒子数量。 */
    vb_.N = 0;
    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      if ( particleInProcessor(iter->rnp[2]) )
	vb_.N++;

    /* 写入 VTK 文件的初始数据。 */
    *vb_.file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
	<< std::endl;
    *vb_.file << "<UnstructuredGrid>"                                                     	<< std::endl;
    *vb_.file << "<Piece NumberOfPoints=\"" << vb_.N + 1 << "\" NumberOfCells=\"" << 1 << "\">"
	<< std::endl;

    /* 插入充电桩网格的坐标。 */
    *vb_.file << "<Points>"                                                             	<< std::endl;
    *vb_.file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" 	<< std::endl;

    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( particleInProcessor(iter->rnp[2]) )
	  *vb_.file << iter->rnp[0] << " " << iter->rnp[1] << " " << iter->rnp[2] 		<< std::endl;
      }

    *vb_.file << xmin_ << " " << ymin_ << " " << zmin_        					<< std::endl;
    *vb_.file << "</DataArray>"                                                          	<< std::endl;
    *vb_.file << "</Points>"                                                              	<< std::endl;

    /* 将每个单元的顶点数写入 VTK 文件。 */
    *vb_.file << "<Cells>"                                                               	<< std::endl;
    *vb_.file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"       	<< std::endl;
    for (unsigned i = 0; i < vb_.N + 1 ; ++i) *vb_.file << i << " ";
    *vb_.file                                                                         		<< std::endl;
    *vb_.file << "</DataArray>"                                                               	<< std::endl;
    *vb_.file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"         	<< std::endl;
    *vb_.file << vb_.N + 1                                                    			<< std::endl;
    *vb_.file << "</DataArray>"                                                            	<< std::endl;
    *vb_.file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"              	<< std::endl;
    *vb_.file << 2                                                                       	<< std::endl;
    *vb_.file << "</DataArray>"                                                           	<< std::endl;
    *vb_.file << "</Cells>"                                                                	<< std::endl;

    *vb_.file << "<PointData Vectors = \"charge\">"                                    		<< std::endl;
    *vb_.file << "<DataArray type=\"Float64\" Name=\"charge\" NumberOfComponents=\"3\" format=\"ascii\">"
	<< std::endl;
    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( particleInProcessor(iter->rnp[2]) )
	  {
	    gamma = sqrt( 1.0 + iter->gb.norm2() );
	    beta  = iter->gb[2] / gamma;
	    *vb_.file << iter->q << " " <<  gamma * gamma_ * ( 1.0 + beta_ * beta )
            												    << " " << gamma * gamma_ * ( 1.0 + beta_ * beta ) * 0.512   << std::endl;
	  }
      }
    *vb_.file << 0.0 << " " << 0.0 << " " << 0.0						<< std::endl;
    *vb_.file << 0.0 << " " << 0.0 << " " << 0.0						<< std::endl;
    *vb_.file << "</DataArray>"                                                             	<< std::endl;
    *vb_.file << "</PointData>"                                                             	<< std::endl;
    *vb_.file << "</Piece>"                                                                 	<< std::endl;
    *vb_.file << "</UnstructuredGrid>"                                                      	<< std::endl;
    *vb_.file << "</VTKFile>"                                                              	<< std::endl;

    /* 关闭文件。 */
    (*vb_.file).close();

    /* 在根处理器上连接 VTK 文件。*/
    if ( rank_ == 0 )
      {
	/* 编写用于合并文件的并行 VTK 文件。 */
	vb_.fileName = bunch_.bunchVTKBasename_ + "-" + stringify(nTimeBunch_) + PTU_FILE_SUFFIX;
	vb_.file = new std::ofstream(vb_.fileName.c_str(),std::ios::trunc);

	/* 获取文件的路径和名称。 */
	splitFilename(bunch_.bunchVTKBasename_, vb_.path, vb_.name);

	*vb_.file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<< std::endl;
	*vb_.file << "<PUnstructuredGrid> GhostLevel = \"0\""                                        	<< std::endl;

	/* 插入电荷云的网格坐标。 */
	*vb_.file << "<PPoints>"                                                                    	<< std::endl;
	*vb_.file << "<PDataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"
	    << std::endl;
	*vb_.file << "</PPoints>"                                                                  	<< std::endl;

	*vb_.file << "<PPointData>"                                                                	<< std::endl;
	*vb_.file << "<PDataArray type=\"Float64\" Name=\"charge\" NumberOfComponents=\"3\" format=\"ascii\" />"
	    << std::endl;
	*vb_.file << "</PPointData>"                                                                	<< std::endl;
	for (int i = 0; i < size_; ++i )
	  {
	    vb_.fileName = vb_.name + "-p" + stringify(i) + "-" + stringify(nTimeBunch_) + VTU_FILE_SUFFIX;
	    *vb_.file << "<Piece  Source=\"" << vb_.fileName << "\"/>"                         		<< std::endl;
	  }
	*vb_.file << "</PUnstructuredGrid>"                                                         	<< std::endl;
	*vb_.file << "</VTKFile>"                                                                  	<< std::endl;

	(*vb_.file).close();
      }
  }

  /******************************************************************************************************
   将该字段的完整概况写入指定的文件中。
   ******************************************************************************************************/

  void Solver::bunchProfile ()
  {
    /* 如果旧文件存在，应将其删除。 */
    pb_.fileName = bunch_.bunchProfileBasename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + TXT_FILE_SUFFIX;
    pb_.file = new std::ofstream(pb_.fileName.c_str(),std::ios::trunc);

    (*pb_.file).setf(std::ios::scientific);
    (*pb_.file).precision(15);
    (*pb_.file).width(40);

    *pb_.file << time_ * gamma_ << std::endl;

    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( particleInProcessor(iter->rnp[2]) )
	  {
	    /* 遍历粒子，并将每个粒子的数据写入文件中。 */
	    *pb_.file << iter->q  	<< "\t";
	    *pb_.file << iter->rnp[0]  	<< "\t";
	    *pb_.file << iter->rnp[1]  	<< "\t";
	    *pb_.file << iter->rnp[2]  	<< "\t";
	    *pb_.file << iter->gb[0] 	<< "\t";
	    *pb_.file << iter->gb[1] 	<< "\t";
	    *pb_.file << iter->gb[2] 	<< std::endl;
	  }
      }

    /* 关闭文件。											*/
    (*pb_.file).close();
  }

  /******************************************************************************************************
   计算波荡器的磁场，并将其叠加到种子场的磁场上。
   ******************************************************************************************************/

  void Solver::undulatorField (UpdateBunchParallel & ubp, FieldVector <Double> & r)
  {
    /* 对于用户指定的每一个外部场，将其添加到波荡器场中。 */
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      {

	/* 计算波荡器磁场。 */
	ub_.b0	= (iter->lu_ != 0.0 ) ? EM * c0_ * 2 * PI / iter->lu_ * iter->k_ / EC : 0.0;

	/* 计算波荡器波数。 */
	ub_.ku 	= (iter->lu_ != 0.0 ) ? 2 * PI / iter->lu_ : 0.0;

	/* 计算波荡器角度的正弦和余弦函数。 */
	ub_.ct	= cos( iter->theta_ );
	ub_.st	= sin( iter->theta_ );

	if ( iter->type_ == STATIC )
	  {
	    /* 首先，确定相对于波荡器起点的相对位置。下方的方程
	     * 假定在 t=0 时刻，粒子束团位于距离第一个
	     * 波荡器 gamma*rb_ 的位置处。										*/
	    ubp.lz = gamma_ * ( r[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) ) - iter->rb_;
	    ubp.ly = r[0] * ub_.ct + r[1] * ub_.st;

	    /* 现在，根据获取的位置计算波荡器场。 */
	    this->staticUndulator(ubp, iter);
	  }
	else if ( iter->type_ == OPTICAL )
	  {
	    /* 将坐标从束团静止系转换至实验室系。 */
	    ubp.rl[0] = r[0]; ubp.rl[1] = r[1];
	    ubp.rl[2] = gamma_ * ( r[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) );
	    ubp.t0    = gamma_ * ( timeBunch_ + dt_  + beta_ / c0_ * r[2] );

	    /* 沿传播方向计算至参考位置的距离。 */
	    ubp.rv = ubp.rl; ubp.rv -= iter->position_;
	    ubp.z  = ubp.rv * iter->direction_ ;

	    /* 计算传播延迟，并将其从时间中减去。 */
	    ubp.tl = ubp.t0 - ubp.z / c0_;

	    /* 重置脉冲的载波包络相位。 */
	    ubp.p0 = 0.0;

	    /* 现在，根据给定的特定种子，对电场矢量进行操作。 */
	    switch ( iter->seedType_ )
	    {
	      case PLANEWAVE:
		this->planeWave(ubp, *iter);			break;

	      case PLANEWAVETRUNCATED:
		this->planeWaveTruncated(ubp, *iter);		break;

	      case GAUSSIANBEAM:
		this->gaussianBeam(ubp, *iter);			break;

	      case SUPERGAUSSIANBEAM:
		this->superGaussianBeam(ubp, *iter);		break;

	      case STANDINGPLANEWAVE:
		this->standingPlaneWave(ubp, *iter);		break;

	      case STANDINGPLANEWAVETRUNCATED:
		this->standingPlaneWaveTruncated(ubp, *iter);	break;

	      case STANDINGGAUSSIANBEAM:
		this->standingGaussianBeam(ubp, *iter);		break;

	      case STANDINGSUPERGAUSSIANBEAM:
		this->standingSuperGaussianBeam(ubp, *iter);	break;
	    }

	    /* 现在，将计算所得的磁矢量势变换至束团静止系。 */
	    ubp.bt[0] += gamma_ * ( ubp.bT[0] + beta_ / c0_ * ubp.eT[1] );
	    ubp.bt[1] += gamma_ * ( ubp.bT[1] - beta_ / c0_ * ubp.eT[0] );
	    ubp.bt[2] += ubp.bT[2];

	    ubp.et[0] += gamma_ * ( ubp.eT[0] - beta_ * c0_ * ubp.bT[1] );
	    ubp.et[1] += gamma_ * ( ubp.eT[1] + beta_ * c0_ * ubp.bT[0] );
	    ubp.et[2] += ubp.eT[2];
	  }
      }
  }

  /******************************************************************************************************
   计算外场场强，并将其叠加到种子的场强上。
   ******************************************************************************************************/

  void Solver::externalField (UpdateBunchParallel & ubp, FieldVector <Double> & r)
  {

    /* 将坐标从束团静止系转换至实验室系。 */
    ubp.rl[0] = r[0];
    ubp.rl[1] = r[1];
    ubp.rl[2] = gamma_ * ( r[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) );
    ubp.t0    = gamma_ * ( timeBunch_ + dt_  + beta_ / c0_ * r[2] );

    /* 对于用户指定的每一个外部场，将其添加到波荡器场中。 */
    for (std::vector<ExtField>::iterator iter = extField_.begin(); iter != extField_.end(); iter++)
      {
	/* 沿传播方向计算至参考位置的距离。 */
	ubp.rv  = ubp.rl; ubp.rv -= iter->position_;
	ubp.z   = ubp.rv * iter->direction_ ;

	/* 计算传播延迟，并将其从时间中减去。 */
	ubp.tl  = ubp.t0 - ubp.z / c0_;
	ubp.tlm = ubp.t0 + ubp.z / c0_;

	/* 重置脉冲的载波包络相位。 */
	ubp.p0  = 0.0;

	/* 现在，根据给定的特定种子，对电场矢量进行操作。 */
	switch ( iter->seedType_ )
	{
	  case PLANEWAVE:
	    this->planeWave(ubp, *iter);			break;

	  case PLANEWAVETRUNCATED:
	    this->planeWaveTruncated(ubp, *iter);		break;

	  case GAUSSIANBEAM:
	    this->gaussianBeam(ubp, *iter);			break;

	  case SUPERGAUSSIANBEAM:
	    this->superGaussianBeam(ubp, *iter);		break;

	  case STANDINGPLANEWAVE:
	    this->standingPlaneWave(ubp, *iter);		break;

	  case STANDINGPLANEWAVETRUNCATED:
	    this->standingPlaneWaveTruncated(ubp, *iter);	break;

	  case STANDINGGAUSSIANBEAM:
	    this->standingGaussianBeam(ubp, *iter);		break;

	  case STANDINGSUPERGAUSSIANBEAM:
	    this->standingSuperGaussianBeam(ubp, *iter);	break;
	}

	/* 现在，将计算所得的磁矢量势变换至束团静止系。 */
	ubp.bt[0] += gamma_ * ( ubp.bT[0] + beta_ / c0_ * ubp.eT[1] );
	ubp.bt[1] += gamma_ * ( ubp.bT[1] - beta_ / c0_ * ubp.eT[0] );
	ubp.bt[2] += ubp.bT[2];

	ubp.et[0] += gamma_ * ( ubp.eT[0] - beta_ * c0_ * ubp.bT[1] );
	ubp.et[1] += gamma_ * ( ubp.eT[1] + beta_ * c0_ * ubp.bT[0] );
	ubp.et[2] += ubp.eT[2];
      }

  }

  /******************************************************************************************************
   初始化用于在给定位置采样并保存辐射能量所需的数据。
   ******************************************************************************************************/

  void Solver::initializeEnergySample ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the FEL radiation energy data.") );

    re_.clear(); re_.resize(FEL_.size());

    /* 遍历不同的 FEL 输出参数，若已启用能量计算功能，则计算辐射能量。 */
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* 仅当电源采样功能启用时进行初始化。 */
	if (!FEL_[jf].radiationEnergy_.sampling_) continue;

	/* 对采样数据执行洛伦兹变换。 */
	for (unsigned int i = 0; i < FEL_[jf].radiationEnergy_.z_.size(); i++)
	  FEL_[jf].radiationEnergy_.z_[i] = FEL_[jf].radiationEnergy_.z_[i] * gamma_;
	FEL_[jf].radiationEnergy_.lineBegin_ = FEL_[jf].radiationEnergy_.lineBegin_ * gamma_;
	FEL_[jf].radiationEnergy_.lineEnd_   = FEL_[jf].radiationEnergy_.lineEnd_   * gamma_;

	/* 如果采样类型为“沿线绘制”，则根据线的起点和终点初始化位置。 */
	if ( FEL_[jf].radiationEnergy_.samplingType_ == OVERLINE )
	  {
	    Double l = 0.0;
	    FieldVector<Double> position;
	    while ( fabs(l) < fabs(FEL_[jf].radiationEnergy_.lineEnd_ - FEL_[jf].radiationEnergy_.lineBegin_) )
	      {
		FEL_[jf].radiationEnergy_.z_.push_back( FEL_[jf].radiationEnergy_.lineBegin_ + l);
		l += FEL_[jf].radiationEnergy_.res_ * gamma_;
	      }
	  }

	re_[jf].N  = FEL_[jf].radiationEnergy_.z_.size();

	/* 将归一化波长扫描添加到波长向量中。*/
	for (Double rw = FEL_[jf].radiationEnergy_.lambdaMin_; rw < FEL_[jf].radiationEnergy_.lambdaMax_; rw += FEL_[jf].radiationEnergy_.lambdaRes_)
	  FEL_[jf].radiationEnergy_.lambda_.push_back(rw);
	re_[jf].Nl = FEL_[jf].radiationEnergy_.lambda_.size();

	/* 根据待计算的点数及线程大小，分配用于存储幂值的内存。 */
	re_[jf].pL.resize(re_[jf].Nl * re_[jf].N, 0.0);
	re_[jf].pG.resize(re_[jf].Nl * re_[jf].N, 0.0);

	re_[jf].Nf = 0;

	/* 初始化用于保存数据的文件流。 */
	re_[jf].file.resize(re_[jf].Nl);
	re_[jf].w.resize(re_[jf].Nl);
	for (unsigned int i = 0; i < re_[jf].Nl; i++)
	  {
	    std::string baseFilename = "";
	    if (!(isabsolute(FEL_[jf].radiationEnergy_.basename_))) baseFilename = FEL_[jf].radiationEnergy_.directory_;
	    baseFilename += FEL_[jf].radiationEnergy_.basename_ + "-" + stringify(i) + TXT_FILE_SUFFIX;

	    /* 如果 baseFilename 所在的目录不存在，则创建该目录。 */
	    createDirectory(baseFilename, rank_);

	    re_[jf].file[i] = new std::ofstream(baseFilename.c_str(),std::ios::trunc);
	    ( *(re_[jf].file[i]) ).setf(std::ios::scientific);
	    ( *(re_[jf].file[i]) ).precision(15);
	    ( *(re_[jf].file[i]) ).width(40);

	    /* 确定计算各辐射谐波幅值所需的时刻点数量。 */
	    Double dt = undulator_[0].lu_ / FEL_[jf].radiationEnergy_.lambda_[i] / ( gamma_ * c0_ );
	    re_[jf].Nf = ( unsigned( dt / mesh_.timeStep_ ) > re_[jf].Nf ) ? unsigned( dt/mesh_.timeStep_ ) : re_[jf].Nf;

	    /* 计算每个波长的角频率。 */
	    re_[jf].w[i] = 2 * PI / dt;
	  }

	/* 根据获取到的 Nf，调整用于存储时域数据的向量大小。*/
	re_[jf].fdt.resize(re_[jf].Nf, std::vector<std::vector<Double> > (re_[jf].N * N1_ * N0_, std::vector<Double> (4,0.0) ) );

	re_[jf].dt = mesh_.timeStep_;
	re_[jf].dx = mesh_.meshResolution_[0];
	re_[jf].dy = mesh_.meshResolution_[1];
	re_[jf].dz = mesh_.meshResolution_[2];

	re_[jf].pc  = 2.0 * re_[jf].dx * re_[jf].dy / ( m0_ * re_[jf].Nf * re_[jf].Nf ) * pow(mesh_.lengthScale_,2) / pow(mesh_.timeScale_,3);

      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The FEL radiation energy data is initialized. :::") );
  }

  /******************************************************************************************************
   在给定位置对辐射能量进行采样，并将其保存至文件。
   ******************************************************************************************************/

  void Solver::energySample ()
  {
    /* 声明计算辐射功率所需的临时参数。 */
    unsigned int              mi, ni;
    FieldVector<Double>       et, bt;
    Complex                   ew1, bw1, ew2, bw2, ex;

    /* 遍历不同的 FEL 输出参数，若已启用能量计算功能，则计算辐射能量。 */
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* 仅当电源采样功能启用时进行初始化。 */
	if (!FEL_[jf].radiationEnergy_.sampling_) continue;

	/* 首先重置所有先前计算的幂。 */
	for (unsigned k = 0; k < re_[jf].N; ++k)
	  for (unsigned l = 0; l < re_[jf].Nl; ++l)
	    re_[jf].pL[k * re_[jf].Nl + l] = 0.0;

	/* 遍历采样位置、横向离散点及频率，以计算特定点和频率处的辐射功率。 */
	for (unsigned k = 0; k < re_[jf].N; ++k)
	  {
	    /* 获取包含该点的单元格的 z-index。 */
	    re_[jf].dzr = modf( ( FEL_[jf].radiationEnergy_.z_[k] - ( mesh_.meshCenter_[2] - mesh_.meshLength_[2] / 2.0 ) ) / mesh_.meshResolution_[2] , &re_[jf].c);
	    re_[jf].k   = (int) re_[jf].c;

	    /* 如果处理器不支持此索引，则不继续执行。 */
	    if ( re_[jf].k < k0_ || re_[jf].k > k0_ + np_ - 1) continue;

	    /* 遍历横向指标。 */
	    for (int i = 2; i < N0_ - 2; i += 1)
	      for (int j = 2; j < N1_ - 2; j += 1)
		{
		  /* 获取计算网格以及场存储网格中的索引。 */
		  mi = ( re_[jf].k - k0_) * N1_* N0_ + i*N1_ + j;
		  ni = k * N1_* N0_ + i*N1_ + j;

		  /* 计算并插值电场，以求得束团点处的数值。 */
		  et[0] = ( 1.0 - re_[jf].dzr ) * en_[mi][0] + re_[jf].dzr * en_[mi+N1N0_][0];
		  et[1] = ( 1.0 - re_[jf].dzr ) * en_[mi][1] + re_[jf].dzr * en_[mi+N1N0_][1];

		  /* 计算并插值磁场，以求得其在束团点处的值。*/
		  bt[0] = ( 1.0 - re_[jf].dzr ) * bn_[mi][0] + re_[jf].dzr * bn_[mi+N1N0_][0];
		  bt[1] = ( 1.0 - re_[jf].dzr ) * bn_[mi][1] + re_[jf].dzr * bn_[mi+N1N0_][1];

		  /* 将场量变换至实验室坐标系。 */
		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );

		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );

		  /* 将该场的贡献计入辐射功率。 */
		  for (unsigned l = 0; l < re_[jf].Nl; l++)
		    {
		      ew1 = Complex (0.0, 0.0);
		      bw1 = Complex (0.0, 0.0);
		      for (unsigned m = 0; m < re_[jf].Nf; m++)
			{
			  ex  = exp( I * ( re_[jf].w[l] * m * mesh_.timeStep_ ) );
			  ew1 += re_[jf].fdt[m][ni][0] * ex;
			  bw1 += re_[jf].fdt[m][ni][3] / ex;
			}

		      ew2 = Complex (0.0, 0.0);
		      bw2 = Complex (0.0, 0.0);
		      for (unsigned m = 0; m < re_[jf].Nf; m++)
			{
			  ex  = exp( I * ( re_[jf].w[l] * m * mesh_.timeStep_ ) );
			  ew2 += re_[jf].fdt[m][ni][1] * ex;
			  bw2 += re_[jf].fdt[m][ni][2] / ex;
			}

		      re_[jf].pL[k * re_[jf].Nl + l] += re_[jf].pc * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
		    }
		}
	  }

	/* 在根处理器上将各处理器的数据汇总。 */
	MPI_Allreduce(&re_[jf].pL[0],&re_[jf].pG[0],re_[jf].N*re_[jf].Nl,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	/* 如果处理器的秩（rank）为零——即为主处理器——则将字段保存到指定文件中。 */
	for (unsigned l = 0; l < re_[jf].Nl; l++)
	  {
	    if ( rank_ == int( l % size_ ) )
	      {
		for (unsigned k = 0; k < re_[jf].N; ++k)
		  *(re_[jf].file[l]) << gamma_ * ( FEL_[jf].radiationEnergy_.z_[k] + beta_ * c0_ * timeBunch_ )
		  - pow( gamma_ * beta_ , 2 ) * undulator_[0].rb_ << "\t" << re_[jf].pG[k * re_[jf].Nl + l] << "\t";
		*(re_[jf].file[l]) << std::endl;
	      }
	  }
      }
  }

  /******************************************************************************************************
   初始化用于存储在给定位置撞击屏幕的粒子所需的数据。
   ******************************************************************************************************/

  void Solver::initializeScreenProfile ()
  {
    scrp_.clear(); scrp_.resize(FEL_.size());
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the screen to measure bunch profile.") );
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* 仅当屏幕已启用时进行初始化。 */
	if (!FEL_[jf].screenProfile_.sampling_) continue;

	if (!(isabsolute(FEL_[jf].screenProfile_.basename_))) FEL_[jf].screenProfile_.basename_ = FEL_[jf].screenProfile_.directory_ + FEL_[jf].screenProfile_.basename_;

	/* 根据给定的节奏，向位置向量添加元素。 */
	if ( FEL_[jf].screenProfile_.rhythm_ > 0.0 )
	  {
	    Double z = 0.0;
	    while ( z < ( undulator_.end()->rb_ + undulator_.end()->length_ * undulator_.end()->lu_ ) )
	      {
		FEL_[jf].screenProfile_.pos_.push_back(z);
		z += FEL_[jf].screenProfile_.rhythm_;
	      }
	  }

	/*  调整存储文件名的向量大小  */
	scrp_[jf].fileNames.resize((FEL_[jf].screenProfile_.pos_).size());
	scrp_[jf].files.resize((FEL_[jf].screenProfile_.pos_).size());

	/* 如果 baseFilename 所在的目录不存在，则创建该目录。*/
	createDirectory(FEL_[jf].screenProfile_.basename_, rank_);

	/* 订单界面 */
	std::sort(FEL_[jf].screenProfile_.pos_.begin(), FEL_[jf].screenProfile_.pos_.end());

	MPI_Barrier(MPI_COMM_WORLD);

	/* 为每个屏幕打开文件，并在第一行写入其位置。 */
	for ( unsigned int i = 0; i < (FEL_[jf].screenProfile_.pos_).size(); i++ )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Screen ") + stringify(i) + std::string(" is at distance ") + stringify(FEL_[jf].screenProfile_.pos_[i]) + std::string(" from the undulator beginning.") );
	    scrp_[jf].fileNames[i] = FEL_[jf].screenProfile_.basename_ + "-p" + stringify(rank_) + "-screen" + stringify(i) + TXT_FILE_SUFFIX;
	    scrp_[jf].files[i] = new std::ofstream(scrp_[jf].fileNames[i].c_str(),std::ios::trunc);
	    (*scrp_[jf].files[i]).setf(std::ios::scientific);
	    (*scrp_[jf].files[i]).precision(15);
	    (*scrp_[jf].files[i]).width(40);
	  }

	/* 如果屏幕采样已启用但未指定屏幕，则返回错误。 */
	if ( FEL_[jf].screenProfile_.pos_.size() == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("No position is set for the screen although the screen sampling is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The screen profiling data is initialized. :::") );
      }
  }

  /******************************************************************************************************
   存储撞击指定位置屏幕的粒子束剖面，并将其保存至文件。
   ******************************************************************************************************/

  void Solver::screenProfile ()
  {
    for ( unsigned jf = 0; jf < FEL_.size(); jf++)
      {
	/* 如果未启用屏幕配置文件采样，则跳过此函数。 */
	if (!FEL_[jf].screenProfile_.sampling_) continue;

	/* 遍历所有屏幕。 */
	for (unsigned i = 0; i < (FEL_[jf].screenProfile_.pos_).size(); i++ )
	  {
	    Double lzScreen = FEL_[jf].screenProfile_.pos_[i];
	    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
	      {
		/* 仅考虑属于当前处理器域的粒子。 */
		if ( particleInProcessor(iter->rnp[2]) )
		  {
		    Double lzm = gamma_ * ( iter->rnm[2] + beta_ * c0_ * ( timeBunch_- mesh_.timeStep_ + dt_ ) );
		    if ( lzm >= lzScreen )	continue;

		    Double lzp = gamma_ * ( iter->rnp[2] + beta_ * c0_ * ( timeBunch_ + dt_) );
		    if ( lzp <  lzScreen )	continue;

		    /* 对量进行插值，并将其写入文件。 */
		    // *scrp_[jf].files[i] << iter->q  	        << "\t";
		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->rnm[0], iter->rnp[0], lzScreen ) << "\t";
		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->rnm[1], iter->rnp[1], lzScreen ) << "\t";
		    Double tm  = gamma_ * ( timeBunch_ + dt_ - mesh_.timeStep_	+ beta_ / c0_ * iter->rnm[2] );
		    Double tp  = gamma_ * ( timeBunch_ + dt_	  		+ beta_ / c0_ * iter->rnp[2] );
		    *scrp_[jf].files[i] << interp( lzm, lzp, tm, tp, lzScreen ) << "\t";

		    //		    /* Use different positions since momenta are found at (time - 1/2*bunchTimeStep).	*/
		    //		    Double gm = sqrt( 1 + pow(iter->gbnm[0], 2) + pow(iter->gbnm[1], 2) + pow(iter->gbnm[2], 2) );
		    //		    lzm -= gamma_ * c0_ *  .5 * bunch_.timeStep_ * ( iter->gbnm[2] / gm + beta_ );
		    //		    Double gp = sqrt( 1 + pow(iter->gbnp[0], 2) + pow(iter->gbnp[1], 2) + pow(iter->gbnp[2], 2) );
		    //		    lzp -= gamma_ * c0_ *  .5 * bunch_.timeStep_ * ( iter->gbnp[2] / gp + beta_ );
		    //		    Double gbzm = gamma_ * ( iter->gbnm[2] + beta_ * gm );
		    //		    Double gbzp = gamma_ * ( iter->gbnp[2] + beta_ * gp );
		    //
		    //		    /* Write the momentum data into the file.						*/
		    //		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->gbnm[0], iter->gbnp[0], lzScreen ) << "\t";
		    //		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->gbnm[1], iter->gbnp[1], lzScreen ) << "\t";
		    //		    *scrp_[jf].files[i] << interp( lzm, lzp, gbzm, gbzp, lzScreen ) << std::endl;

		    /* Here, we use the stair-case approximation and set the momentum during the
		     * whole time step constant.							*/
		    *scrp_[jf].files[i] << iter->gb[0] << "\t";
		    *scrp_[jf].files[i] << iter->gb[1] << "\t";
		    *scrp_[jf].files[i] << gamma_ * (iter->gb[2] + beta_ * std::sqrt(1.0 + iter->gb.norm2())) << std::endl;
		  }
	      }
	  }
      }
  }

  /******************************************************************************************************
   完成字段计算。
   ******************************************************************************************************/

  void Solver::finalize ()
  {
    /* 如果采样处于活动状态，则关闭文件，因为仿真现已结束。 */
    if (seed_.sampling_ && sf_.N > 0) 	(*sf_.file).close();

    /* 如果启用了成组采样，则关闭文件，因为仿真现已结束。 */
    if (bunch_.sampling_) 	(*sb_.file).close();
  }

  /******************************************************************************************************
   
   *****************************************************************定义用于比较波荡器起始位置的布尔函数。*************************************/

  bool Solver::undulatorCompare (Undulator i, Undulator j)
  { return ( i.rb_ < j.rb_ ); }

  /****************************************************************************************************
   定义线性插值函数。
   ****************************************************************************************************/

  Double Solver::interp( Double x0, Double x1, Double y0, Double y1, Double x )
  {
    return y0 + ( x - x0 ) / ( x1 - x0 ) * ( y1 - y0 );
  }

  /****************************************************************************************************
   定义用于检查粒子是否属于处理器的函数。
   ****************************************************************************************************/

  bool Solver::particleInProcessor( const Double& z )
  {
    Double zr = pmod( z - zmin_ , mesh_.meshLength_[2] ) + zmin_;
    return ( ( zr >= zp_[0] ) && ( zr < zp_[1] ) );
  }

  /****************************************************************************************************
   定义一个函数，用于根据单元格索引返回实际坐标。
   ****************************************************************************************************/

  FieldVector<Double> Solver::rc( const long int& m )
  {
    unsigned int i, j, k;
    k = m / N1N0_;
    i = ( m % N1N0_ ) / N1_;
    j = m - ( N1N0_*k+N1_*i );

    FieldVector<Double> v;
    v[0] = xmin_ + i         * mesh_.meshResolution_[0];
    v[1] = ymin_ + j         * mesh_.meshResolution_[1];
    v[2] = zmin_ + (k + k0_) * mesh_.meshResolution_[2];

    return (v);
  }

}
