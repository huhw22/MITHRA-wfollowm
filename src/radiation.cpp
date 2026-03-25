/********************************************************************************************************
* radiation.cpp：实现辐射功率的计算功能
********************************************************************************************************/

#include <string>

#include "fieldvector.h"
#include "solver.h"
#include "stdinclude.h"

namespace MITHRA
{

  /******************************************************************************************************
  *初始化采样所需的数据，并在给定位置保存辐射功率。
  ******************************************************************************************************/

  void Solver::initializePowerSample()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the data for FEL radiation power sampling.") );
    rp_.clear(); rp_.resize(FEL_.size());

    /*对不同的FEL输出参数进行循环，如果是，则初始化功率计算
	*激活。*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/*当且仅当电源采样启用时初始化。*/
	if (!FEL_[jf].radiationPower_.sampling_) continue;

	/*对采样数据执行洛伦兹升压。*/
	for (unsigned int i = 0; i < FEL_[jf].radiationPower_.z_.size(); i++)
	  FEL_[jf].radiationPower_.z_[i] 	*= gamma_;
	FEL_[jf].radiationPower_.lineBegin_ 	*= gamma_;
	FEL_[jf].radiationPower_.lineEnd_   	*= gamma_;

	Double dl = fabs(FEL_[jf].radiationPower_.lineEnd_ - FEL_[jf].radiationPower_.lineBegin_) / FEL_[jf].radiationPower_.res_;

	/*如果采样类型为行上绘图，则根据行开始和初始化位置
	*行结束。*/
	if ( FEL_[jf].radiationPower_.samplingType_ == OVERLINE )
	  {
	    Double l = 0.0;
	    FieldVector<Double> position;
	    while ( fabs(l) < fabs(FEL_[jf].radiationPower_.lineEnd_ - FEL_[jf].radiationPower_.lineBegin_) )
	      {
		FEL_[jf].radiationPower_.z_.push_back( FEL_[jf].radiationPower_.lineBegin_ + l);
		l += dl;
	      }
	  }

	/*设置每个处理器的采样点数。*/
	rp_[jf].N  = FEL_[jf].radiationPower_.z_.size();
	rp_[jf].Nz = 0;
	for (unsigned int i = 0; i < rp_[jf].N; i++)
	  if ( FEL_[jf].radiationPower_.z_[i] < zp_[1] && FEL_[jf].radiationPower_.z_[i] >= zp_[0] )
	    ++rp_[jf].Nz;

	/*将归一化波长扫描添加到波长矢量中。*/
	dl = ( FEL_[jf].radiationPower_.lambdaMax_ - FEL_[jf].radiationPower_.lambdaMin_ ) / FEL_[jf].radiationPower_.lambdaRes_;
	for (Double rw = FEL_[jf].radiationPower_.lambdaMin_; rw < FEL_[jf].radiationPower_.lambdaMax_; rw += dl)
	  FEL_[jf].radiationPower_.lambda_.push_back(rw);
	rp_[jf].Nl = FEL_[jf].radiationPower_.lambda_.size();

	/*根据要计算的点数和线程的大小为其分配内存
	*拯救力量。*/
	rp_[jf].pL.resize(rp_[jf].Nl * rp_[jf].N, 0.0);
	rp_[jf].pG.resize(rp_[jf].Nl * rp_[jf].N, 0.0);

	rp_[jf].Nf = 0;

	/*初始化文件流以保存数据。*/
	rp_[jf].file.resize(rp_[jf].Nl);
	rp_[jf].w.resize(rp_[jf].Nl);
	for (unsigned int i = 0; i < rp_[jf].Nl; i++)
	  {
	    std::string baseFilename = "";
	    if (!(isabsolute(FEL_[jf].radiationPower_.basename_))) baseFilename = FEL_[jf].radiationPower_.directory_;
	    baseFilename += FEL_[jf].radiationPower_.basename_ + "-" + stringify(i) + TXT_FILE_SUFFIX;

	    /*如果baseFilename所在目录不存在，则创建该目录。*/
	    createDirectory(baseFilename, rank_);

	    rp_[jf].file[i] = new std::ofstream(baseFilename.c_str(),std::ios::trunc);
	    ( *(rp_[jf].file[i]) ).setf(std::ios::scientific);
	    ( *(rp_[jf].file[i]) ).precision(15);
	    ( *(rp_[jf].file[i]) ).width(40);

	    /*确定计算每次辐射振幅所需的时间点数目
		*谐波。在这里，我们使用三个辐射周期内的功率来计算
		*瞬时功率在选定的谐波。*/
	    Double dt = undulator_[0].lu_ / FEL_[jf].radiationPower_.lambda_[i] / ( gamma_ * c0_ );
	    rp_[jf].Nf = ( unsigned( 3.0 * dt / mesh_.timeStep_ ) > rp_[jf].Nf ) ? unsigned( 3.0 * dt / mesh_.timeStep_ ) : rp_[jf].Nf;

	    /*计算每个波长的角频率。*/
	    rp_[jf].w[i] = 2 * PI / dt;
	  }

	/*根据得到的Nf，调整矢量大小以保存时域数据。*/
	rp_[jf].fdt.resize(rp_[jf].Nf, std::vector<std::vector<Double> > (rp_[jf].Nz * N1_ * N0_, std::vector<Double> (4,0.0) ) );

	rp_[jf].dt  = mesh_.timeStep_;
	rp_[jf].dx  = mesh_.meshResolution_[0];
	rp_[jf].dy  = mesh_.meshResolution_[1];
	rp_[jf].dz  = mesh_.meshResolution_[2];

	rp_[jf].pc  = 2.0 * rp_[jf].dx * rp_[jf].dy / ( m0_ * rp_[jf].Nf * rp_[jf].Nf ) * pow(mesh_.lengthScale_,2) / pow(mesh_.timeScale_,3);

	/*从获得的波长矢量设置傅里叶系数的大小并初始化
	*数据。*/
	rp_[jf].ep.resize(rp_[jf].Nl, std::vector<Complex> (rp_[jf].Nf, Complex (0.0, 0.0) ) );
	rp_[jf].em.resize(rp_[jf].Nl, std::vector<Complex> (rp_[jf].Nf, Complex (0.0, 0.0) ) );
	for (unsigned int i = 0; i < rp_[jf].Nl; i++)
	  for (unsigned int j = 0; j < rp_[jf].Nf; j++)
	    {
	      rp_[jf].ep[i][j] = cos( rp_[jf].w[i] * j * mesh_.timeStep_ ) + I * sin( rp_[jf].w[i] * j * mesh_.timeStep_ );
	      rp_[jf].em[i][j] = cos( rp_[jf].w[i] * j * mesh_.timeStep_ ) - I * sin( rp_[jf].w[i] * j * mesh_.timeStep_ );
	    }
      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The data for sampling FEL radiation power is initialized. :::") );
  }

  /******************************************************************************************************
  *在给定位置取样辐射功率并保存到文件中。
  ******************************************************************************************************/

  void Solver::powerSample()
  {
    /*声明计算辐射功率所需的临时参数。*/
    unsigned int              	kz;
    long int                  	mi, ni;
    FieldVector<Double>       	et, bt;
    Complex                   	ew1, bw1, ew2, bw2;

    /*在不同的FEL输出参数上进行循环，计算出辐射能量
	*计算被激活。*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/*当且仅当电源采样启用时初始化。*/
	if (!FEL_[jf].radiationPower_.sampling_) continue;

	/*首先重置之前计算的所有能量。*/
	for (unsigned k = 0; k < rp_[jf].N; ++k)
	  for (unsigned l = 0; l < rp_[jf].Nl; ++l)
	    rp_[jf].pL[k * rp_[jf].Nl + l] = 0.0;

	/*设置采样点的索引为0。*/
	kz = 0;

	/*循环遍历采样位置，横向定向，和频率来计算
	*在特定点和频率处的辐射功率。*/

	/*K个索引在采样位置上循环。*/
	for (unsigned k = 0; k < rp_[jf].N; ++k)
	  {
	    /*如果处理器不支持此索引，请不要继续。*/
	    if ( !( FEL_[jf].radiationPower_.z_[k] < zp_[1] && FEL_[jf].radiationPower_.z_[k] >= zp_[0] ) ) continue;

	    /*获得包含该点的单元格的z索引。*/
	    rp_[jf].dzr = modf( ( FEL_[jf].radiationPower_.z_[k] - zmin_ ) / mesh_.meshResolution_[2] , &rp_[jf].c);
	    rp_[jf].k   = (int) rp_[jf].c;

	    /*获取用于功率计算的时间序列中的索引。*/
	    rp_[jf].m	= nTime_ % rp_[jf].Nf;

	    /*在横向指标上绕圈。*/
	    for (int i = 2; i < N0_ - 2; i += 1)
	      for (int j = 2; j < N1_ - 2; j += 1)
		{
		  /*获取计算网格和字段存储网格中的索引。*/
		  mi = ( rp_[jf].k - k0_) * N1_* N0_ + i * N1_ + j;
		  ni = kz * N1_* N0_ + i * N1_ + j;

		  /*计算相应像素的字段。*/
		  if (!pic_[mi      ]) fieldEvaluate(mi      );
		  if (!pic_[mi+N1N0_]) fieldEvaluate(mi+N1N0_);

		  /*计算并插值电场，求出束点处的值。*/
		  et[0] = ( 1.0 - rp_[jf].dzr ) * en_[mi][0] + rp_[jf].dzr * en_[mi+N1N0_][0];
		  et[1] = ( 1.0 - rp_[jf].dzr ) * en_[mi][1] + rp_[jf].dzr * en_[mi+N1N0_][1];

		  /*计算并插值磁场以求其在束点处的值。*/
		  bt[0] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][0] + rp_[jf].dzr * bn_[mi+N1N0_][0];
		  bt[1] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][1] + rp_[jf].dzr * bn_[mi+N1N0_][1];

		  /*把力场转换成实验室框架。*/
		  rp_[jf].fdt[rp_[jf].m][ni][0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
		  rp_[jf].fdt[rp_[jf].m][ni][1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );

		  rp_[jf].fdt[rp_[jf].m][ni][2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
		  rp_[jf].fdt[rp_[jf].m][ni][3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );

		  /*把这个场的贡献加到辐射功率上。*/

		  /*L指数环在给定波长的功率采样。*/
		  for ( unsigned l = 0; l < rp_[jf].Nl; l++)
		    {
		      ew1 = Complex (0.0, 0.0); bw1 = ew1; ew2 = ew1; bw2 = ew1;

		      for ( unsigned m = 0; m < rp_[jf].Nf; m++)
			{
			  ew1 += rp_[jf].fdt[m][ni][0] * rp_[jf].ep[l][m];
			  bw1 += rp_[jf].fdt[m][ni][3] * rp_[jf].em[l][m];
			  ew2 += rp_[jf].fdt[m][ni][1] * rp_[jf].ep[l][m];
			  bw2 += rp_[jf].fdt[m][ni][2] * rp_[jf].em[l][m];
			}

		      /*把贡献加到幂级数上。*/
		      rp_[jf].pL[k * rp_[jf].Nl + l] += rp_[jf].pc * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
		    }
		}

	    /*将采样点的迭代器加1。*/
	    kz += 1;
	  }

	/*将来自每个处理器的数据一起添加到根处理器。*/
	MPI_Allreduce(&rp_[jf].pL[0],&rp_[jf].pG[0],rp_[jf].N*rp_[jf].Nl,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	/*如果处理器的秩等于零，即根处理器将字段保存到
	*给定文件。*/
	for ( unsigned l = 0; l < rp_[jf].Nl; l++)
	  {
	    if ( rank_ == int( l % size_ ) )
	      {
		for (unsigned k = 0; k < rp_[jf].N; ++k)
		  *(rp_[jf].file[l]) << gamma_ * ( FEL_[jf].radiationPower_.z_[k] + beta_ * c0_ * ( timeBunch_ + dt_ ) ) << "\t" << rp_[jf].pG[k * rp_[jf].Nl + l] << "\t";
		*(rp_[jf].file[l]) << std::endl;
	      }
	  }
      }
  }


  void Solver::initializeDetector ()
  {
	printmessage(std::string(__FILE__), __LINE__,
				std::string("::: Initializing the data for FEL radiation detector sampling.") );

	rd_.clear();
	rd_.resize(FEL_.size());

	for (unsigned int jf = 0; jf < FEL_.size(); jf++)
	{
	  /* 先把运行态清掉 */
	  rd_[jf].active     = false;
	  rd_[jf].started    = false;
	  rd_[jf].finished   = false;
	  rd_[jf].file       = NULL;

	  rd_[jf].zLab       = FEL_[jf].radiationDetector_.zLab_;
	  rd_[jf].zBox       = 0.0;
	  rd_[jf].zLabCheck  = 0.0;

	  rd_[jf].tLab0      = 0.0;
	  rd_[jf].tLab       = 0.0;

	  rd_[jf].ownerRank  = -1;

	  rd_[jf].w          = 0.0;
	  rd_[jf].Nf         = 0;
	  rd_[jf].m          = 0;

	  rd_[jf].k          = -1;

	  rd_[jf].dt         = 0.0;
	  rd_[jf].dx         = 0.0;
	  rd_[jf].dy         = 0.0;
	  rd_[jf].dz         = 0.0;
	  rd_[jf].dzr        = 0.0;
	  rd_[jf].c          = 0.0;

	  rd_[jf].pc         = 0.0;

	  rd_[jf].pL         = 0.0;
	  rd_[jf].pG         = 0.0;

	  rd_[jf].planeL.clear();
	  rd_[jf].planeG.clear();
	  rd_[jf].fdt.clear();
	  rd_[jf].ep.clear();
	  rd_[jf].em.clear();

	  /* 没开 detector 就跳过 */
	  if (!FEL_[jf].radiationDetector_.sampling_) continue;

	  /* 第一版只处理 power-line */
	  if (!FEL_[jf].radiationDetector_.writePowerLine_) continue;

	  /* -------- 以下完全仿照 initializePowerSample 的频域初始化逻辑 -------- */

	  /* 单频 detector 的等效周期 */
	  Double dt = undulator_[0].lu_ /
				FEL_[jf].radiationDetector_.lambda_ /
				( gamma_ * c0_ );

	  /* 和 radiation-power 完全同一规则：用三个辐射周期决定时间窗长度 */
	  rd_[jf].Nf = unsigned( 3.0 * dt / mesh_.timeStep_ );

	  /* 单频角频率 */
	  rd_[jf].w  = 2.0 * PI / dt;

	  /* 时间/空间步长 */
	  rd_[jf].dt = mesh_.timeStep_;
	  rd_[jf].dx = mesh_.meshResolution_[0];
	  rd_[jf].dy = mesh_.meshResolution_[1];
	  rd_[jf].dz = mesh_.meshResolution_[2];

	  /* 功率归一化系数：直接照抄 radiation-power 的 pc 定义 */
	  rd_[jf].pc = 2.0 * rd_[jf].dx * rd_[jf].dy
				/ ( m0_ * rd_[jf].Nf * rd_[jf].Nf )
				* pow(mesh_.lengthScale_, 2)
				/ pow(mesh_.timeScale_, 3);

	  /* 当前整张 detector 面的局部/全局场快照
	  * 每个 xy 点存 4 个量：Ex, Ey, Bx, By（lab 系） */
	  rd_[jf].planeL.resize(N1_ * N0_ * 4, 0.0);
	  rd_[jf].planeG.resize(N1_ * N0_ * 4, 0.0);

	  /* 时间环形缓冲：
	  * [Nf][N1_*N0_][4] */
	  rd_[jf].fdt.resize(
		rd_[jf].Nf,
		std::vector< std::vector<Double> >(N1_ * N0_, std::vector<Double>(4, 0.0))
	  );

	  /* 单频 Fourier 系数：
	  * detector 这里只有一个频点，所以 ep/em 只是一维 */
	  rd_[jf].ep.resize(rd_[jf].Nf, Complex(0.0, 0.0));
	  rd_[jf].em.resize(rd_[jf].Nf, Complex(0.0, 0.0));

	  for (unsigned int j = 0; j < rd_[jf].Nf; j++)
	  {
		rd_[jf].ep[j] = cos( rd_[jf].w * j * mesh_.timeStep_ )
						+ I * sin( rd_[jf].w * j * mesh_.timeStep_ );

		rd_[jf].em[j] = cos( rd_[jf].w * j * mesh_.timeStep_ )
						- I * sin( rd_[jf].w * j * mesh_.timeStep_ );
	  }
	}

	printmessage(std::string(__FILE__), __LINE__,
				std::string(" The data for sampling FEL radiation detector is initialized. :::") );
  }


  void Solver::detectorSample ()
{
  /* 临时变量 */
  long int                  mi, ni;
  FieldVector<Double>       et, bt;
  Complex                   ew1, bw1, ew2, bw2;

  int                       ownerLocal, ownerGlobal;
  bool                      activeNow;

  for (unsigned int jf = 0; jf < FEL_.size(); jf++)
    {
      if (!FEL_[jf].radiationDetector_.sampling_) continue;
      if (!FEL_[jf].radiationDetector_.writePowerLine_) continue;

      /* 每一步先清当前瞬时状态 */
      rd_[jf].active     = false;
      rd_[jf].ownerRank  = -1;
      rd_[jf].pL         = 0.0;
      rd_[jf].pG         = 0.0;

      std::fill(rd_[jf].planeL.begin(), rd_[jf].planeL.end(), 0.0);
      std::fill(rd_[jf].planeG.begin(), rd_[jf].planeG.end(), 0.0);

      ownerLocal  = -1;
      ownerGlobal = -1;
      activeNow   = false;

      /* 与 powerSample 保持同一位置映射时刻定义 */
	  rd_[jf].zBox      = boostFrame_.boxZFromLabZAndBoxT(rd_[jf].zLab, timeBunch_);
	  rd_[jf].zLabCheck = boostFrame_.labZFromBoxZT(rd_[jf].zBox, timeBunch_);
	  rd_[jf].tLab      = boostFrame_.labTFromBoxZT(rd_[jf].zBox, timeBunch_);

      /* ----------- 1) 先判断当前 detector 面是否在全局有效区间 ----------- */
      if ( rd_[jf].zBox >= zmin_ + mesh_.meshResolution_[2] &&
           rd_[jf].zBox <= zmax_ - 2.0 * mesh_.meshResolution_[2] )
        activeNow = true;

      /* ----------- 2) 如果已经开始过，但这一步已经不在有效区，则结束记录 ----------- */
      if (!activeNow && rd_[jf].started && !rd_[jf].finished)
        {
          if (rank_ == 0 && rd_[jf].file != NULL)
            {
              rd_[jf].file->flush();
              rd_[jf].file->close();
              delete rd_[jf].file;
              rd_[jf].file = NULL;
            }

          rd_[jf].finished = true;

          printmessage(std::string(__FILE__), __LINE__,
                       std::string("::: Detector recording ends.") );
        }

      /* 当前无效，或者已经结束，后面都不再处理 */
      if (!activeNow || rd_[jf].finished)
        continue;

      rd_[jf].active = true;

      /* ----------- 3) 第一次进入有效区：打开文件并打印提示 ----------- */
      if (!rd_[jf].started)
        {
          rd_[jf].tLab0    = rd_[jf].tLab;
          rd_[jf].started  = true;
          rd_[jf].finished = false;

          if (rank_ == 0)
            {
              std::string baseFilename = "";
              if (!(isabsolute(FEL_[jf].radiationDetector_.basename_)))
                baseFilename = FEL_[jf].radiationDetector_.directory_;

              baseFilename += FEL_[jf].radiationDetector_.basename_
                           + "-power-line" + TXT_FILE_SUFFIX;

              createDirectory(baseFilename, rank_);

              rd_[jf].file = new std::ofstream(baseFilename.c_str(), std::ios::trunc);
              (*(rd_[jf].file)).setf(std::ios::scientific);
              (*(rd_[jf].file)).precision(15);
              (*(rd_[jf].file)).width(40);

              (*(rd_[jf].file))
                << "# time_rel\t"
                << "z_box_abs\t"
                << "owner_rank\t"
                << "power"
                << std::endl;
            }

          printmessage(std::string(__FILE__), __LINE__,
                       std::string("::: Detector recording starts.") );
        }

      /* ----------- 4) 判断当前 detector 面属于哪个 rank ----------- */
      if ( rd_[jf].zBox >= zp_[0] &&
           rd_[jf].zBox <  zp_[1] - mesh_.meshResolution_[2] )
        {
          ownerLocal = rank_;

          /* 与 powerSample 一样，按 z 方向线性插值 */
          rd_[jf].dzr = modf( ( rd_[jf].zBox - zmin_ ) / mesh_.meshResolution_[2], &rd_[jf].c );
          rd_[jf].k   = (long int) rd_[jf].c;

          /* 只由 owner rank 生成当前整张 detector 面的局部场快照 planeL
           * 每个 xy 点存四个量：Ex, Ey, Bx, By（lab 系） */
          for (int i = 2; i < N0_ - 2; i += 1)
            for (int j = 2; j < N1_ - 2; j += 1)
              {
                mi = ( rd_[jf].k - k0_ ) * N1_ * N0_ + i * N1_ + j;
                ni = i * N1_ + j;

                if (!pic_[mi        ]) fieldEvaluate(mi        );
                if (!pic_[mi+N1N0_  ]) fieldEvaluate(mi+N1N0_  );

                et[0] = ( 1.0 - rd_[jf].dzr ) * en_[mi][0] + rd_[jf].dzr * en_[mi+N1N0_][0];
                et[1] = ( 1.0 - rd_[jf].dzr ) * en_[mi][1] + rd_[jf].dzr * en_[mi+N1N0_][1];

                bt[0] = ( 1.0 - rd_[jf].dzr ) * bn_[mi][0] + rd_[jf].dzr * bn_[mi+N1N0_][0];
                bt[1] = ( 1.0 - rd_[jf].dzr ) * bn_[mi][1] + rd_[jf].dzr * bn_[mi+N1N0_][1];

                /* 与 powerSample 完全相同：先变到 lab 系，再存入时间缓冲 */
                rd_[jf].planeL[4 * ni + 0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
                rd_[jf].planeL[4 * ni + 1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );
                rd_[jf].planeL[4 * ni + 2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
                rd_[jf].planeL[4 * ni + 3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );
              }
        }

      /* ----------- 5) 先把 owner rank 算出来，再聚合整张 detector 面 ----------- */
      MPI_Allreduce(&ownerLocal, &ownerGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      rd_[jf].ownerRank = ownerGlobal;

      /* 如果当前没有任何 rank 持有这个面，则认为本步无效并结束 */
      if (rd_[jf].ownerRank < 0)
        {
          rd_[jf].active = false;

          if (rd_[jf].started && !rd_[jf].finished)
            {
              if (rank_ == 0 && rd_[jf].file != NULL)
                {
                  rd_[jf].file->flush();
                  rd_[jf].file->close();
                  delete rd_[jf].file;
                  rd_[jf].file = NULL;
                }

              rd_[jf].finished = true;

              printmessage(std::string(__FILE__), __LINE__,
                           std::string("::: Detector recording ends.") );
            }

          continue;
        }

      MPI_Allreduce(&rd_[jf].planeL[0], &rd_[jf].planeG[0],
                    int(rd_[jf].planeL.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      /* ----------- 6) 与 powerSample 一样，把当前 plane 写入时间环形缓冲 ----------- */
      rd_[jf].m = nTime_ % rd_[jf].Nf;

      for (int i = 0; i < N0_; i += 1)
        for (int j = 0; j < N1_; j += 1)
          {
            ni = i * N1_ + j;

            rd_[jf].fdt[rd_[jf].m][ni][0] = rd_[jf].planeG[4 * ni + 0];
            rd_[jf].fdt[rd_[jf].m][ni][1] = rd_[jf].planeG[4 * ni + 1];
            rd_[jf].fdt[rd_[jf].m][ni][2] = rd_[jf].planeG[4 * ni + 2];
            rd_[jf].fdt[rd_[jf].m][ni][3] = rd_[jf].planeG[4 * ni + 3];
          }

      /* ----------- 7) 与 powerSample 完全相同的 Fourier 功率计算 ----------- */
      if (rank_ == rd_[jf].ownerRank)
        {
          rd_[jf].pL = 0.0;

          for (int i = 2; i < N0_ - 2; i += 1)
            for (int j = 2; j < N1_ - 2; j += 1)
              {
                ni = i * N1_ + j;

                ew1 = Complex(0.0, 0.0);
                bw1 = Complex(0.0, 0.0);
                ew2 = Complex(0.0, 0.0);
                bw2 = Complex(0.0, 0.0);

                for (unsigned int m = 0; m < rd_[jf].Nf; m += 1)
                  {
                    ew1 += rd_[jf].fdt[m][ni][0] * rd_[jf].ep[m];
                    bw1 += rd_[jf].fdt[m][ni][3] * rd_[jf].em[m];
                    ew2 += rd_[jf].fdt[m][ni][1] * rd_[jf].ep[m];
                    bw2 += rd_[jf].fdt[m][ni][2] * rd_[jf].em[m];
                  }

                rd_[jf].pL += rd_[jf].pc
                            * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
              }
        }

      MPI_Allreduce(&rd_[jf].pL, &rd_[jf].pG, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      /* ----------- 8) 只由 rank 0 写 power-line ----------- */
      if (rank_ == 0 && rd_[jf].file != NULL)
        {
          (*(rd_[jf].file))
            << rd_[jf].tLab - rd_[jf].tLab0 << "\t"
            << rd_[jf].zBox                  << "\t"
            << rd_[jf].ownerRank             << "\t"
            << rd_[jf].pG
            << std::endl;
        }

      /* ----------- 9) 如果下一步就到仿真终点，则做最终收尾 ----------- */
      if ( time_ + mesh_.timeStep_ >= mesh_.totalTime_ &&
           rd_[jf].started && !rd_[jf].finished )
        {
          if (rank_ == 0 && rd_[jf].file != NULL)
            {
              rd_[jf].file->flush();
              rd_[jf].file->close();
              delete rd_[jf].file;
              rd_[jf].file = NULL;
            }

          rd_[jf].finished = true;

          printmessage(std::string(__FILE__), __LINE__,
                       std::string("::: Detector recording ends at simulation stop.") );
        }
    }
}
  /******************************************************************************************************
  *初始化在给定位置显示辐射功率所需的数据。
  ******************************************************************************************************/

  void Solver::initializePowerVisualize()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the data for FEL radiation power visualization.") );

    /*对不同的FEL输出参数进行循环，如果是，则初始化功率计算
	*激活。*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/*当且仅当电源采样启用时初始化。*/
	if (!FEL_[jf].vtkPower_.sampling_) continue;

	/*如果功率显示节奏仍然为零，则返回错误。*/
	if ( FEL_[jf].vtkPower_.rhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The power visualization rhythm of the field is zero although power visualization is activated !!!") );
	    exit(1);
	  }

	/*洛伦兹将功率可视化采样节奏提高到电子静止框架。*/
	FEL_[jf].vtkPower_.rhythm_	/= gamma_;

	/*对采样数据执行洛伦兹升压。*/
	FEL_[jf].vtkPower_.z_ 		*= gamma_;

	/*创建用于保存可视化数据的文件名。*/
	if (!(isabsolute(FEL_[jf].vtkPower_.basename_)))
	  FEL_[jf].vtkPower_.basename_ = FEL_[jf].vtkPower_.directory_ + FEL_[jf].vtkPower_.basename_;

	/*如果baseFilename所在目录不存在，则创建该目录。*/
	createDirectory(FEL_[jf].vtkPower_.basename_, rank_);

	/*设置每个处理器的采样点数。*/
	rp_[jf].N  = 1;
	rp_[jf].Nz = ( FEL_[jf].vtkPower_.z_ < zp_[1] && FEL_[jf].vtkPower_.z_ >= zp_[0] ) ? 1 : 0;
	rp_[jf].Nl = 1;

	/*如果Nz不等于1，则不要继续循环。*/
	if ( rp_[jf].Nz == 0 ) continue;

	/*根据要计算的点数和线程的大小为其分配内存
	*拯救力量。*/
	rp_[jf].pL.resize(N1_*N0_, 0.0);
	rp_[jf].pG.clear();

	/*初始化文件流以保存数据。*/
	rp_[jf].file.resize(rp_[jf].Nz);
	rp_[jf].w.resize(rp_[jf].Nz);

	/*确定计算每次辐射振幅所需的时间点数目
	*谐波。在这里，我们使用三个辐射周期内的功率来计算瞬时
	*所选谐波处的功率。*/
	Double dt = undulator_[0].lu_ / FEL_[jf].vtkPower_.lambda_ / ( gamma_ * c0_ );
	rp_[jf].Nf = unsigned( 3.0 * dt / mesh_.timeStep_ );

	/*计算每个波长的角频率。*/
	rp_[jf].w[0] = 2 * PI / dt;

	/*根据得到的Nf，调整矢量大小以保存时域数据。*/
	rp_[jf].fdt.resize(rp_[jf].Nf, std::vector<std::vector<Double> > (N1_ * N0_, std::vector<Double> (4, 0.0) ) );

	rp_[jf].dt  = mesh_.timeStep_;
	rp_[jf].dx  = mesh_.meshResolution_[0];
	rp_[jf].dy  = mesh_.meshResolution_[1];
	rp_[jf].dz  = mesh_.meshResolution_[2];

	rp_[jf].pc  = 2.0 * rp_[jf].dx * rp_[jf].dy / ( m0_ * rp_[jf].Nf * rp_[jf].Nf ) * pow(mesh_.lengthScale_,2) / pow(mesh_.timeScale_,3);

	/*从获得的波长矢量设置傅里叶系数的大小并初始化
	*数据。*/
	rp_[jf].ep.resize(rp_[jf].Nl, std::vector<Complex> (rp_[jf].Nf, Complex (0.0, 0.0) ) );
	rp_[jf].em.resize(rp_[jf].Nl, std::vector<Complex> (rp_[jf].Nf, Complex (0.0, 0.0) ) );
	for (unsigned int i = 0; i < rp_[jf].Nl; i++)
	  for (unsigned int j = 0; j < rp_[jf].Nf; j++)
	    {
	      rp_[jf].ep[i][j] = cos( rp_[jf].w[i] * j * mesh_.timeStep_ ) + I * sin( rp_[jf].w[i] * j * mesh_.timeStep_ );
	      rp_[jf].em[i][j] = cos( rp_[jf].w[i] * j * mesh_.timeStep_ ) - I * sin( rp_[jf].w[i] * j * mesh_.timeStep_ );
	    }
      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The data for FEL radiation power visualization is initialized. :::") );
  }

  /******************************************************************************************************
  *将给定位置的辐射功率可视化并保存到文件中。
  ******************************************************************************************************/

  void Solver::powerVisualize()
  {

    /*声明计算辐射功率所需的临时参数。*/
    unsigned int              	m;
    long int                  	mi, ni;
    FieldVector<Double>       	et, bt;
    Complex                   	ew1, bw1, ew2, bw2;

    /*在不同的FEL输出参数上进行循环，并将辐射功率可视化
	*被激活。*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {

	/*当且仅当电源采样启用时初始化。*/
	if ( FEL_[jf].vtkPower_.sampling_ && rp_[jf].Nz == 1 )
	  {

	    /*计算平面所在单元格的索引。*/
	    rp_[jf].dzr 	= modf( ( FEL_[jf].vtkPower_.z_ - zmin_ ) / mesh_.meshResolution_[2] , &rp_[jf].c);
	    rp_[jf].k   	= (int) rp_[jf].c - k0_;

	    /*获取用于功率计算的时间序列中的索引。*/
	    rp_[jf].m	= nTime_ % rp_[jf].Nf;

	    /*在横向指标上绕圈。*/
	    for (int i = 1; i < N0_ - 1; i++)
	      for (int j = 1; j < N1_ - 1; j++)
		{
		  /*获取计算网格和字段存储网格中的索引。*/
		  mi = rp_[jf].k * N1_* N0_ + i * N1_ + j;
		  ni = i * N1_ + j;

		  /*计算相应像素的字段。*/
		  if (!pic_[mi      ]) fieldEvaluate(mi      );
		  if (!pic_[mi+N1N0_]) fieldEvaluate(mi+N1N0_);

		  /*计算并插值电场，求出束点处的值。*/
		  et[0] = ( 1.0 - rp_[jf].dzr ) * en_[mi][0] + rp_[jf].dzr * en_[mi+N1N0_][0];
		  et[1] = ( 1.0 - rp_[jf].dzr ) * en_[mi][1] + rp_[jf].dzr * en_[mi+N1N0_][1];

		  /*计算并插值磁场以求其在束点处的值。*/
		  bt[0] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][0] + rp_[jf].dzr * bn_[mi+N1N0_][0];
		  bt[1] = ( 1.0 - rp_[jf].dzr ) * bn_[mi][1] + rp_[jf].dzr * bn_[mi+N1N0_][1];

		  /*把力场转换成实验室框架。*/
		  rp_[jf].fdt[rp_[jf].m][ni][0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
		  rp_[jf].fdt[rp_[jf].m][ni][1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );

		  rp_[jf].fdt[rp_[jf].m][ni][2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
		  rp_[jf].fdt[rp_[jf].m][ni][3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );

		  /*L指数环在给定波长的功率采样。*/
		  for ( unsigned l = 0; l < rp_[jf].Nl; l++)
		    {
		      ew1 = Complex (0.0, 0.0); bw1 = ew1; ew2 = ew1; bw2 = ew1;

		      for (unsigned m = 0; m < rp_[jf].Nf; m++)
			{
			  ew1 += rp_[jf].fdt[m][ni][0] * rp_[jf].ep[l][m];
			  bw1 += rp_[jf].fdt[m][ni][3] * rp_[jf].em[l][m];
			  ew2 += rp_[jf].fdt[m][ni][1] * rp_[jf].ep[l][m];
			  bw2 += rp_[jf].fdt[m][ni][2] * rp_[jf].em[l][m];
			}

		      rp_[jf].pL[ni] = rp_[jf].pc * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
		    }
		}

	    if ( fmod(time_, FEL_[jf].vtkPower_.rhythm_) < mesh_.timeStep_ )
	      {

		/*旧文件如果存在，应该删除。*/
		std::string baseFilename = FEL_[jf].vtkPower_.basename_ + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
		rp_[jf].file[0] = new std::ofstream(baseFilename.c_str(),std::ios::trunc);

		rp_[jf].file[0]->setf(std::ios::scientific);
		rp_[jf].file[0]->precision(4);

		/*为vtk文件写入初始数据。*/
		*rp_[jf].file[0] << "<?xml version=\"1.0\"?>"						<< std::endl;
		*rp_[jf].file[0] << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
		    "compressor=\"vtkZLibDataCompressor\">" 						<< std::endl;
		*rp_[jf].file[0] << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
		    0 << " " << 0 << "\">"								<< std::endl;
		*rp_[jf].file[0] << "<Piece Extent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
		    0 << " " << 0 << "\">"								<< std::endl;

		/*插入充电点的网格坐标。*/
		*rp_[jf].file[0] << "<Points>"                                                        	<< std::endl;
		*rp_[jf].file[0] << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<< std::endl;
		FieldVector<Double> r1 (0.0), r2 (0.0);
		for (int j = 0; j < N1_; j++)
		  for (int i = 0; i < N0_; i++)
		    {
		      m = rp_[jf].k * N1_ * N0_ + i * N1_ + j;
		      r1 = rc(m); r2 = rc(m + N1N0_);
		      *rp_[jf].file[0] << r1[0] * ( 1.0 - rp_[jf].dzr ) + r2[0] * rp_[jf].dzr << " "
			  << r1[1] << " " << r1[2] 							<< std::endl;
		    }
		*rp_[jf].file[0] << "</DataArray>"                                                    	<< std::endl;
		*rp_[jf].file[0] << "</Points>"                                                       	<< std::endl;

		/*将每个单元格数据插入到vtk文件中。*/
		*rp_[jf].file[0] << "<CellData>"                                                       	<< std::endl;
		*rp_[jf].file[0] << "</CellData>"                                                      	<< std::endl;

		/*根据计算的电场插入点数据。*/
		*rp_[jf].file[0] << "<PointData Vectors = \"power\">"                                 	<< std::endl;
		*rp_[jf].file[0] << "<DataArray type=\"Float64\" Name=\"power\" NumberOfComponents=\"" << 1 << "\" format=\"ascii\">"
		    << std::endl;
		for (int j = 0; j < N1_; j++)
		  for (int i = 0; i < N0_; i++)
		    *rp_[jf].file[0] << rp_[jf].pL[i * N1_ + j] 					<< std::endl;

		*rp_[jf].file[0] << "</DataArray>"                                                   	<< std::endl;
		*rp_[jf].file[0] << "</PointData>"                                                    	<< std::endl;
		*rp_[jf].file[0] << "</Piece>"                                                        	<< std::endl;
		*rp_[jf].file[0] << "</StructuredGrid>"                                               	<< std::endl;
		*rp_[jf].file[0] << "</VTKFile>"                                                      	<< std::endl;

		/*关闭文件。*/
		(*rp_[jf].file[0]).close();
	      }
	  }
      }
  }

}       /*命名空间Darius结束。*/
