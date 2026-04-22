/********************************************************************************************************
* fdtd.cpp：实现函数的实时时域解
*代码
********************************************************************************************************/

#include "fdtd.h"

namespace MITHRA
{
  FdTd::FdTd( Mesh& 				mesh,
	      Bunch& 				bunch,
	      Seed& 				seed,
	      std::vector<Undulator>&		undulator,
	      std::vector<ExtField>& 		extField,
	      std::vector<FreeElectronLaser>& 	FEL )
  : Solver ( mesh, bunch, seed, undulator, extField, FEL )
  {};

  /******************************************************************************************************
  *将电流复位为零。
  ******************************************************************************************************/

  void FdTd::currentReset ()
  {
    /*为了减少对内存消耗的要求，我们使用向量anp1来存储电流。*/

    Double*  jn = &(*anp1_)[0][0];
    Double*  je = &(*anp1_)[(long)N1N0_*np_-1][2];
    while ( jn != je )
      *(jn++) = 0.0;
    *je = 0.0;
  }

  /******************************************************************************************************
  *更新单元格点的电流以进行字段更新。
  ******************************************************************************************************/

  void FdTd::currentUpdate ()
  {
    /*为了减少对内存消耗的要求，我们使用向量anp1来存储电流。*/

    FieldVector<Double>*	      jn   = &(*anp1_)[0];
    bool                              bp, bm;
    std::list<Charge>::iterator       it = chargeVectorn_.begin();

	long long nUsedLocal = 0;         // 本步真正参与沉积的粒子数
	long long nSuppressedLocal = 0;   // qeff 比原始 q 明显变小的粒子数
	Double qAbsLocal = 0.0;           // sum |q|
	Double qEffAbsLocal = 0.0;        // sum |qeff|

    /*现在，应该在电荷上执行一个循环，并且应该更新电流。*/
    for (it = chargeVectorn_.begin(); it != chargeVectorn_.end(); it++)
      {
	/*求出粒子的电流。*/
	uc_.rp  = it->rnp;
	uc_.rm  = it->rnm;

	/*保存检测粒子在处理器域中的标志。*/
	bp = (
	    uc_.rp[0] < xmax_ - uc_.dx && uc_.rp[0] > xmin_ + uc_.dx &&
	    uc_.rp[1] < ymax_ - uc_.dy && uc_.rp[1] > ymin_ + uc_.dy &&
	    uc_.rp[2] < zp_[1]         && uc_.rp[2] >= zp_[0] );
	bm = (
	    uc_.rm[0] < xmax_ - uc_.dx && uc_.rm[0] > xmin_ + uc_.dx &&
	    uc_.rm[1] < ymax_ - uc_.dy && uc_.rm[1] > ymin_ + uc_.dy &&
	    uc_.rm[2] < zp_[1]         && uc_.rm[2] >= zp_[0] );

	/*如果以上条件都不满足，则继续循环。*/
	if ( ! (bp || bm) ) continue;

	/*得到粒子的电荷。*/
	uc_.q = it->q * 0.5 * (it->wm + it->w);
	++nUsedLocal;
	qAbsLocal    += std::abs(it->q);
	qEffAbsLocal += std::abs(uc_.q);

	if (std::abs(uc_.q) < 0.999999 * std::abs(it->q))
		++nSuppressedLocal;

	if (std::abs(uc_.q) < 1.0e-20)
		continue;

	/*得到下一个时间步长的宏观粒子的指数。*/
	uc_.ip  = (int) floor( ( uc_.rp[0] - xmin_ ) / uc_.dx );
	uc_.jp  = (int) floor( ( uc_.rp[1] - ymin_ ) / uc_.dy );
	uc_.kp  = (int) floor( ( uc_.rp[2] - zmin_ ) / uc_.dz );

	/*得到上一个时间步长的宏观粒子的指数。*/
	uc_.im  = (int) floor( ( uc_.rm[0] - xmin_ ) / uc_.dx );
	uc_.jm  = (int) floor( ( uc_.rm[1] - ymin_ ) / uc_.dy );
	uc_.km  = (int) floor( ( uc_.rm[2] - zmin_ ) / uc_.dz );

	/*计算中继点。*/
	uc_.r[0]  = std::min( std::min( uc_.im, uc_.ip ) * uc_.dx + uc_.dx + xmin_,
			      std::max( std::max( uc_.im, uc_.ip ) * uc_.dx + xmin_, 0.5 * (uc_.rm[0] + uc_.rp[0]) ) );
	uc_.r[1]  = std::min( std::min( uc_.jm, uc_.jp ) * uc_.dy + uc_.dy + ymin_,
			      std::max( std::max( uc_.jm, uc_.jp ) * uc_.dy + ymin_, 0.5 * (uc_.rm[1] + uc_.rp[1]) ) );
	uc_.r[2]  = std::min( std::min( uc_.km, uc_.kp ) * uc_.dz + uc_.dz + zmin_,
			      std::max( std::max( uc_.km, uc_.kp ) * uc_.dz + zmin_, 0.5 * (uc_.rm[2] + uc_.rp[2]) ) );

	/*计算电荷通量。*/
	uc_.jcm  = uc_.r;
	uc_.jcm -= uc_.rm;
	uc_.jcp  = uc_.rp;
	uc_.jcp -= uc_.r;

	/*如果电荷在计算域之外，则停止模拟。*/
	if ( bp )
	  {
	    /*找出要考虑其电流的节点的索引。*/
	    uc_.m   = N1N0_ * ( uc_.kp - k0_ ) + N1_ * uc_.ip + uc_.jp;

	    uc_.dxp = modf( ( 0.5 * ( uc_.rp[0] + uc_.r[0] ) - xmin_ ) / uc_.dx , &uc_.c );
	    uc_.dyp = modf( ( 0.5 * ( uc_.rp[1] + uc_.r[1] ) - ymin_ ) / uc_.dy , &uc_.c );
	    uc_.dzp = modf( ( 0.5 * ( uc_.rp[2] + uc_.r[2] ) - zmin_ ) / uc_.dz , &uc_.c );

	    /*计算每个顶点对电流的贡献。*/
	    uc_.x1  = 1.0 - uc_.dxp;
	    uc_.x2  = uc_.dxp;
	    uc_.y1  = 1.0 - uc_.dyp;
	    uc_.y2  = uc_.dyp;
	    uc_.z1  = 1.0 - uc_.dzp;
	    uc_.z2  = uc_.dzp;

	    (*(jn+uc_.m)            )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+N1_)        )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+1  )        )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+N1_+1)      )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_)      )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_+N1_)  )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_+1)    )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcp[0];
	    (*(jn+uc_.m+N1N0_+N1_+1))[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcp[0];

	    (*(jn+uc_.m)            )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+N1_)        )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+1  )        )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+N1_+1)      )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_)      )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_+N1_)  )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_+1)    )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcp[1];
	    (*(jn+uc_.m+N1N0_+N1_+1))[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcp[1];

	    (*(jn+uc_.m)            )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+N1_)        )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+1  )        )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcp[2];
	    (*(jn+uc_.m+N1_+1)      )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_)      )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_+N1_)  )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_+1)    )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcp[2];
	    (*(jn+uc_.m+N1N0_+N1_+1))[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcp[2];
	  }

	/*如果电荷在计算域之外，则停止模拟。*/
	if ( bm )
	  {
	    /*找出要考虑其电流的节点的索引。*/
	    uc_.m   = N1N0_ * ( uc_.km - k0_ ) + N1_ * uc_.im + uc_.jm;

	    uc_.dxm = modf( ( 0.5 * ( uc_.rm[0] + uc_.r[0] ) - xmin_ ) / uc_.dx , &uc_.c );
	    uc_.dym = modf( ( 0.5 * ( uc_.rm[1] + uc_.r[1] ) - ymin_ ) / uc_.dy , &uc_.c );
	    uc_.dzm = modf( ( 0.5 * ( uc_.rm[2] + uc_.r[2] ) - zmin_ ) / uc_.dz , &uc_.c );

	    /*计算每个顶点对电流的贡献。*/
	    uc_.x1  = 1.0 - uc_.dxm;
	    uc_.x2  = uc_.dxm;
	    uc_.y1  = 1.0 - uc_.dym;
	    uc_.y2  = uc_.dym;
	    uc_.z1  = 1.0 - uc_.dzm;
	    uc_.z2  = uc_.dzm;

	    (*(jn+uc_.m)            )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+N1_)        )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+1  )        )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+N1_+1)      )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z1 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_)      )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_+N1_)  )[0] += uc_.q * 0.5 * uc_.y1 * uc_.z2 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_+1)    )[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcm[0];
	    (*(jn+uc_.m+N1N0_+N1_+1))[0] += uc_.q * 0.5 * uc_.y2 * uc_.z2 * uc_.jcm[0];

	    (*(jn+uc_.m)            )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+N1_)        )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+1  )        )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+N1_+1)      )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z1 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_)      )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_+N1_)  )[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_+1)    )[1] += uc_.q * 0.5 * uc_.x1 * uc_.z2 * uc_.jcm[1];
	    (*(jn+uc_.m+N1N0_+N1_+1))[1] += uc_.q * 0.5 * uc_.x2 * uc_.z2 * uc_.jcm[1];

	    (*(jn+uc_.m)            )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+N1_)        )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+1  )        )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcm[2];
	    (*(jn+uc_.m+N1_+1)      )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_)      )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_+N1_)  )[2] += uc_.q * 0.5 * uc_.x2 * uc_.y1 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_+1)    )[2] += uc_.q * 0.5 * uc_.x1 * uc_.y2 * uc_.jcm[2];
	    (*(jn+uc_.m+N1N0_+N1_+1))[2] += uc_.q * 0.5 * uc_.x2 * uc_.y2 * uc_.jcm[2];
	  }
    }
	// long long nUsedGlobal = 0;
	// long long nSuppressedGlobal = 0;
	// Double qAbsGlobal = 0.0;
	// Double qEffAbsGlobal = 0.0;

	// MPI_Allreduce(&nUsedLocal, &nUsedGlobal, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	// MPI_Allreduce(&nSuppressedLocal, &nSuppressedGlobal, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	// MPI_Allreduce(&qAbsLocal, &qAbsGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// MPI_Allreduce(&qEffAbsLocal, &qEffAbsGlobal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// if (rank_ == 0 && (nTime_ % 50 == 0))
	// {
	// 	printmessage(std::string(__FILE__), __LINE__,
	// 		std::string("[softkill fdtd currentUpdate] nUsed=") + stringify(nUsedGlobal) +
	// 		" nSuppressed=" + stringify(nSuppressedGlobal) +
	// 		" sum|q|=" + stringify(qAbsGlobal) +
	// 		" sum|qeff|=" + stringify(qEffAbsGlobal) +
	// 		" ratio=" + stringify(qEffAbsGlobal / (qAbsGlobal + 1.0e-300)));
	// }
  }

  /******************************************************************************************************
  *在不同的处理器之间传输电流。
  ******************************************************************************************************/

  void FdTd::currentCommunicate ()
  {
    /*为了减少对内存消耗的要求，我们使用向量anp1来存储电流。*/

    int                       	msgtag9 = 9;
    MPI_Status                	status;
    std::list<Charge>::iterator	it = chargeVectorn_.begin();
    Double			zr;

    /*将每个处理器对边界处电荷和电流密度的贡献相加。*/
    if (rank_ != 0)
      MPI_Send(&(*anp1_)[0][0],                                    3*N1N0_,MPI_DOUBLE,rank_-1,msgtag9, MPI_COMM_WORLD);

    if (rank_ != size_ - 1)
      {
	MPI_Recv(&uc_.jt[0][0],                               3*N1N0_,MPI_DOUBLE,rank_+1,msgtag9, MPI_COMM_WORLD,&status);

	for (int i = 0; i < N1N0_; i++)
	  (*anp1_)[(np_-2)*N1N0_+i] += uc_.jt[i];
      }

    /*现在电荷和电流密度已经沉积，从
	*列表。*/
    it = chargeVectorn_.begin();
    while ( it != chargeVectorn_.end() )
      {
	zr = pmod( it->rnp[2] - zmin_ , mesh_.meshLength_[2] ) + zmin_;
	if      ( zr <  zp_[0] )
	  it = chargeVectorn_.erase(it);
	else if ( zr >= zp_[1] )
	  it = chargeVectorn_.erase(it);
	else
	  ++it;
      }
  }

  /******************************************************************************************************
  *更新一个时间步的字段
  ******************************************************************************************************/

  void FdTd::fieldUpdate ()
  {
    /*定义更新字段时临时需要的值。*/
    unsigned int 	i, j, k;
    long int            m, l;
    MPI_Status 		status;
    int 		msgtag1 = 1, msgtag3 = 3;
    int                 msgtag5 = 5, msgtag6 = 6, msgtag7 = 7, msgtag8 = 8;
    FieldVector<Double>	atemp; atemp = 0.0;

    uf_.anp1 = &(*anp1_)[0][0];
    uf_.an   = &(*an_)  [0][0];
    uf_.anm1 = &(*anm1_)[0][0];
    uf_.jn   = &(*anp1_)[0][0];
    uf_.en   = &en_[0][0];
    uf_.bn   = &bn_[0][0];

    const long int L0  = 3*N1_;
    const long int L1  = 3*N1N0_;
    const long int L2  = 3*(N1_+N1N0_);
    const long int L3  = 3*(N1_-N1N0_);
    const long int L4  = 3*(1+N1N0_);
    const long int L5  = 3*(1-N1N0_);
    const long int L6  = 3*(1+N1_);
    const long int L7  = 3*(1-N1_);
    const long int L8  = 3*(1+N1_+N1N0_);
    const long int L9  = 3*(1+N1_-N1N0_);
    const long int L10 = 3*(1-N1_+N1N0_);
    const long int L11 = 3*(1-N1_-N1N0_);

    /*首先将所有单元格中的粒子标记设置为零。*/
    std::vector<bool>::iterator  pn = pic_.begin();
    std::vector<bool>::iterator  pe = pic_.end();
    while ( pn != pe ) *(pn++) = false;

    /*循环遍历网格中的点并按一个时间步更新字段。自，时间更新
	*对于计算域中的点与边界上的点不同，我们将
	*在内部点上先做一个循环。在所有这些更新之后，边界上的点
	*将相应更新。*/
    if ( mesh_.solver_ == NSFD )
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    for (unsigned k = 1; k < uf_.npm1; k++)
	    {
			m = N1N0_ * k + N1_ * i + j;
			l = 3 * m;

			const bool useCPML = cpmlXY_.enabled && !inPhysicalXY(i, j);

			if (useCPML)
			{
				uf_.af.advanceMagneticPotentialNSFD_CPMLXY(
					uf_.anp1+l,    uf_.anm1+l,    uf_.an+l,
					uf_.an  +l+L0, uf_.an  +l+L2, uf_.an+l+L3,
					uf_.an  +l-L0, uf_.an  +l-L3, uf_.an+l-L2,
					uf_.an  +l+3 , uf_.an  +l+L4, uf_.an+l+L5,
					uf_.an  +l-3 , uf_.an  +l-L5, uf_.an+l-L4,
					uf_.an  +l+L1, uf_.an  +l-L1, uf_.jn+l,

					&cpmlXY_.psiXnp1_[m][0],
					&cpmlXY_.psiXn_[m][0],
					&cpmlXY_.psiXn_[m + N1_][0],
					&cpmlXY_.psiXn_[m - N1_][0],

					&cpmlXY_.psiYnp1_[m][0],
					&cpmlXY_.psiYn_[m][0],
					&cpmlXY_.psiYn_[m + 1][0],
					&cpmlXY_.psiYn_[m - 1][0],

					cpmlXY_.bx[i], cpmlXY_.cx[i], cpmlXY_.kappaX[i],
					cpmlXY_.by[j], cpmlXY_.cy[j], cpmlXY_.kappaY[j],
					uf_.dx, uf_.dy);
			}
			else
			{
				uf_.af.advanceMagneticPotentialNSFD(
					uf_.anp1+l,    uf_.anm1+l,    uf_.an+l,
					uf_.an  +l+L0, uf_.an  +l+L2, uf_.an+l+L3,
					uf_.an  +l-L0, uf_.an  +l-L3, uf_.an+l-L2,
					uf_.an  +l+3 , uf_.an  +l+L4, uf_.an+l+L5,
					uf_.an  +l-3 , uf_.an  +l-L5, uf_.an+l-L4,
					uf_.an  +l+L1, uf_.an  +l-L1, uf_.jn+l);
			}
	    }
      }
    else if ( mesh_.solver_ == FD )
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    for (unsigned k = 1; k < uf_.npm1; k++)
	    {
			m = N1N0_ * k + N1_ * i + j;
			l = 3 * m;

			const bool useCPML = cpmlXY_.enabled && !inPhysicalXY(i, j);

			if (useCPML)
			{
				uf_.af.advanceMagneticPotentialFD_CPMLXY(
					uf_.anp1+l,    uf_.anm1+l,    uf_.an+l,
					uf_.an  +l+L0, uf_.an  +l+L2, uf_.an+l+L3,
					uf_.an  +l-L0, uf_.an  +l-L3, uf_.an+l-L2,
					uf_.an  +l+3 , uf_.an  +l+L4, uf_.an+l+L5,
					uf_.an  +l-3 , uf_.an  +l-L5, uf_.an+l-L4,
					uf_.an  +l+L1, uf_.an  +l-L1, uf_.jn+l,

					&cpmlXY_.psiXnp1_[m][0],
					&cpmlXY_.psiXn_[m][0],
					&cpmlXY_.psiXn_[m + N1_][0],
					&cpmlXY_.psiXn_[m - N1_][0],

					&cpmlXY_.psiYnp1_[m][0],
					&cpmlXY_.psiYn_[m][0],
					&cpmlXY_.psiYn_[m + 1][0],
					&cpmlXY_.psiYn_[m - 1][0],

					cpmlXY_.bx[i], cpmlXY_.cx[i], cpmlXY_.kappaX[i],
					cpmlXY_.by[j], cpmlXY_.cy[j], cpmlXY_.kappaY[j],
					uf_.dx, uf_.dy);
			}
			else
			{
				uf_.af.advanceMagneticPotentialFD(
					uf_.anp1+l,    uf_.anm1+l,    uf_.an+l,
					uf_.an  +l+L0, uf_.an  +l+L2, uf_.an+l+L3,
					uf_.an  +l-L0, uf_.an  +l-L3, uf_.an+l-L2,
					uf_.an  +l+3 , uf_.an  +l+L4, uf_.an+l+L5,
					uf_.an  +l-3 , uf_.an  +l-L5, uf_.an+l-L4,
					uf_.an  +l+L1, uf_.an  +l-L1, uf_.jn+l);
			}
	    }
      }

    /*如果种子的振幅超过一定的极限，则将种子注入计算域
	*通过TF/SF边界。*/
    if ( fabs(seed_.amplitude_) > 1.0e-50 )
      {
	unsigned int KI = ( rank_ == 0         ) ? 2       : 1;
	unsigned int KF = ( rank_ == size_ - 1 ) ? np_ - 2 : np_ - 1;

	for ( int j = 2; j < N1_-2; j++)
	  for ( unsigned k = KI; k < KF; k++)
	    {
	      i = 1;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m+N1_),	time_, atemp); (*anp1_)[m].mmv(uf_.a[1], atemp);
	      i = 2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m-N1_),	time_, atemp); (*anp1_)[m].pmv(uf_.a[1], atemp);
	      i = N0_-2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m-N1_),	time_, atemp); (*anp1_)[m].mmv(uf_.a[1], atemp);
	      i = N0_-3;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m+N1_),	time_, atemp); (*anp1_)[m].pmv(uf_.a[1], atemp);
	    }

	for ( int i = 2; i < N0_-2; i++)
	  for ( unsigned k = KI; k < KF; k++)
	    {
	      j = 1;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m+1), time_, atemp); (*anp1_)[m].mmv(uf_.a[2], atemp);
	      j = 2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m-1), time_, atemp); (*anp1_)[m].pmv(uf_.a[2], atemp);
	      j = N1_-2;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m-1), time_, atemp); (*anp1_)[m].mmv(uf_.a[2], atemp);
	      j = N1_-3;
	      m = N1N0_ * k + N1_ * i + j;
	      seed_.fields(rc(m+1),		time_, atemp); (*anp1_)[m].pmv(uf_.a[2], atemp);
	    }

	if ( rank_ == 0 )
	  {
	    for (int i = 2; i < N0_-2; i++)
	      for (int j = 2; j < N1_-2; j++)
		{
		  k = 1;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(rc(m+N1N0_), time_, atemp); (*anp1_)[m].mmv(uf_.a[3], atemp);
		  k = 2;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(rc(m-N1N0_), time_, atemp); (*anp1_)[m].pmv(uf_.a[3], atemp);
		}
	  }

	if ( rank_ == size_-1 )
	  {
	    for (int i = 2; i < N0_-2; i++)
	      for (int j = 2; j < N1_-2; j++)
		{
		  k = np_-2;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(rc(m-N1N0_),	time_, atemp); (*anp1_)[m].mmv(uf_.a[3], atemp);
		  k = np_-3;
		  m = N1N0_ * k + N1_ * i + j;
		  seed_.fields(rc(m+N1N0_),	time_, atemp); (*anp1_)[m].pmv(uf_.a[3], atemp);
		}
	  }
      }

    /*在x = xmin边界上循环网格中的点，并使用第一个更新字段
	*阶吸收边界条件。*/
    uf_.af.ufB_ = &uf_.bB[0];
    for (unsigned j = 1; j < uf_.N1m1; j++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  l = 3 * ( N1N0_ * k + j );

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l+L0,	uf_.an  +l+L0, 	uf_.anp1+l+L0,
	      uf_.an  +l+L6,  	uf_.an  +l-L7,	uf_.an  +l+L2,
	      uf_.an  +l+L3, 	uf_.an  +l+3,   uf_.an  +l-3,
	      uf_.an  +l+L1,  	uf_.an  +l-L1);
	}

    /*在x = xmax边界上循环网格中的点，并使用第一个更新字段
	*阶吸收边界条件。*/
    for (unsigned j = 1; j < uf_.N1m1; j++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  l = 3 * ( N1N0_ * k + N1N0_ - N1_ + j );

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l-L0,	uf_.an  +l-L0,	uf_.anp1+l-L0,
	      uf_.an  +l+L7,	uf_.an  +l-L6,	uf_.an  +l-L3,
	      uf_.an  +l-L2, 	uf_.an  +l+3,	uf_.an  +l-3,
	      uf_.an  +l+L1,  	uf_.an  +l-L1 );
	}

    /*循环网格中y = ymin边界上的点，并使用第一个更新字段
	*阶吸收边界条件。*/
    uf_.af.ufB_ = &uf_.cB[0];
    for (unsigned i = 1; i < uf_.N0m1; i++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  l = 3 * ( N1N0_ * k + N1_* i );

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l+3,   	uf_.an  +l+3,   uf_.anp1+l+3,
	      uf_.an  +l+L6,  	uf_.an  +l+L7, 	uf_.an  +l+L4,
	      uf_.an  +l+L5, 	uf_.an  +l+L0,  uf_.an  +l-L0,
	      uf_.an  +l+L1,  	uf_.an  +l-L1);
	}

    /*循环网格中y = ymax边界上的点，并使用第一个更新字段
	*阶吸收边界条件。*/
    for (unsigned i = 1; i < uf_.N0m1; i++)
      for (unsigned k = 1; k < uf_.npm1; k++)
	{
	  l = 3 * ( N1N0_ * k + N1_* i + N1_ - 1 );

	  uf_.af.advanceBoundaryF(
	      uf_.anp1+l,	uf_.anm1+l,	uf_.an  +l,
	      uf_.anm1+l-3, 	uf_.an  +l-3,  	uf_.anp1+l-3,
	      uf_.an  +l-L7,  	uf_.an  +l-L6, 	uf_.an  +l-L5,
	      uf_.an  +l-L4,	uf_.an  +l+L0,  uf_.an  +l-L0,
	      uf_.an  +l+L1,  	uf_.an  +l-L1);
	}

    /*在z = zmin边界上循环网格中的点，并使用第一个更新字段
	*阶吸收边界条件。*/
    uf_.af.ufB_ = &uf_.dB[0];
    if ( rank_ == 0 )
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    {
	      l = 3 * ( N1_ * i + j );

	      uf_.af.advanceBoundaryF(
		  uf_.anp1+l,		uf_.anm1+l,		uf_.an  +l,
		  uf_.anm1+l+L1,     	uf_.an  +l+L1,   	uf_.anp1+l+L1,
		  uf_.an  +l+L2,   	uf_.an  +l-L3, 		uf_.an  +l+L4,
		  uf_.an  +l-L5, 	uf_.an  +l+L0,     	uf_.an  +l-L0,
		  uf_.an  +l+3,     	uf_.an  +l-3);
	    }
      }

    /*在z = zmax边界上循环网格中的点，并使用第一个更新字段
	*阶吸收边界条件。*/
    if (rank_ == size_ - 1)
      {
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  for (unsigned j = 1; j < uf_.N1m1; j++)
	    {
	      l = 3 * ( N1N0_ * uf_.npm1 + N1_ * i + j );

	      uf_.af.advanceBoundaryF(
		  uf_.anp1+l,		uf_.anm1+l,		uf_.an  +l,
		  uf_.anm1+l-L1,     	uf_.an  +l-L1,   	uf_.anp1+l-L1,
		  uf_.an  +l+L3,   	uf_.an  +l-L2, 		uf_.an  +l+L5,
		  uf_.an  +l-L4, 	uf_.an  +l+L0,     	uf_.an  +l-L0,
		  uf_.an  +l+3,    	uf_.an  +l-3);
	    }
      }

    if ( mesh_.truncationOrder_ == 2 )
      {

	/*循环网格中x = (xmin,xmax)和y = (ymin,ymax)边界上的边缘点
	*并使用一阶吸收边界条件更新字段。理解
	*下面几行代码，最好在旁边列出Li的值。*/
	uf_.af.ufB_ = &uf_.eE[0];
	for (unsigned k = 1; k < uf_.npm1; k++)
	  {
	    l = 3 * ( N1N0_ * k );

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l+L0,		uf_.an+l+L0,	uf_.anm1+l+L0,
		uf_.anp1+l+3,		uf_.an+l+3,	uf_.anm1+l+3,
		uf_.anp1+l+L6,		uf_.an+l+L6,	uf_.anm1+l+L6,
		uf_.an  +l-L1,		uf_.an+l+L3,	uf_.an  +l+L5,	uf_.an+l+L9,
		uf_.an  +l+L1,		uf_.an+l+L2,	uf_.an  +l+L4,	uf_.an+l+L8);

	    l = 3 * ( N1N0_ * k + N1_ * uf_.N0m1 );

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l-L0,		uf_.an+l-L0,	uf_.anm1+l-L0,
		uf_.anp1+l+3,		uf_.an+l+3,	uf_.anm1+l+3,
		uf_.anp1+l+L7,		uf_.an+l+L7,	uf_.anm1+l+L7,
		uf_.an  +l-L1,		uf_.an+l-L2,	uf_.an  +l+L5,	uf_.an+l+L11,
		uf_.an  +l+L1,		uf_.an+l-L3,	uf_.an  +l+L4,	uf_.an+l+L10);

	    l = 3 * ( N1N0_ * k + uf_.N1m1 );

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l+L0,		uf_.an+l+L0,	uf_.anm1+l+L0,
		uf_.anp1+l-3,		uf_.an+l-3,	uf_.anm1+l-3,
		uf_.anp1+l-L7,		uf_.an+l-L7,	uf_.anm1+l-L7,
		uf_.an  +l-L1,		uf_.an+l+L3,	uf_.an  +l-L4,	uf_.an+l-L10,
		uf_.an  +l+L1,		uf_.an+l+L2,	uf_.an  +l-L5,	uf_.an+l-L11);

	    l = 3 * ( N1N0_ * k + N1_ * uf_.N0m1 + uf_.N1m1 );

	    uf_.af.advanceEdgeF(
		uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		uf_.anp1+l-L0,		uf_.an+l-L0,	uf_.anm1+l-L0,
		uf_.anp1+l-3,		uf_.an+l-3,	uf_.anm1+l-3,
		uf_.anp1+l-L6,		uf_.an+l-L6,	uf_.anm1+l-L6,
		uf_.an  +l-L1,		uf_.an+l-L2,	uf_.an  +l-L4,	uf_.an+l-L8,
		uf_.an  +l+L1,		uf_.an+l-L3,	uf_.an  +l-L5,	uf_.an+l-L9);
	  }

	/*循环网格中z = (zmin,zmax)和y = (ymin,ymax)边界上的边缘点和
	*使用一阶吸收边界条件更新场。要了解
	*在代码行之后，最好在旁边列出Li的值。*/
	uf_.af.ufB_ = &uf_.fE[0];
	for (unsigned i = 1; i < uf_.N0m1; i++)
	  {
	    if ( rank_ == 0 )
	      {
		l = 3 * ( N1_ * i );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+3,	uf_.an+l+3,	uf_.anm1+l+3,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l+L4,	uf_.an+l+L4,	uf_.anm1+l+L4,
		    uf_.an  +l-L0,	uf_.an+l+L7,	uf_.an  +l-L3,	uf_.an+l+L10,
		    uf_.an  +l+L0,	uf_.an+l+L6,	uf_.an  +l+L2,	uf_.an+l+L8);

		l = 3 * ( N1_ * i + uf_.N1m1 );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-3,	uf_.an+l-3,	uf_.anm1+l-3,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l-L5,	uf_.an+l-L5,	uf_.anm1+l-L5,
		    uf_.an  +l-L0,	uf_.an+l-L6,	uf_.an  +l-L3,	uf_.an+l-L9,
		    uf_.an  +l+L0,	uf_.an+l-L7,	uf_.an  +l+L2,	uf_.an+l-L11);
	      }

	    if ( rank_ == size_ - 1 )
	      {
		l = 3 * ( N1_ * i + N1N0_ * uf_.npm1 );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+3,	uf_.an+l+3,	uf_.anm1+l+3,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l+L5,	uf_.an+l+L5,	uf_.anm1+l+L5,
		    uf_.an  +l-L0,	uf_.an+l+L7,	uf_.an  +l-L2,	uf_.an+l+L11,
		    uf_.an  +l+L0,	uf_.an+l+L6,	uf_.an  +l+L3,	uf_.an+l+L9);

		l = 3 * ( N1_ * i + N1N0_ * uf_.npm1 + uf_.N1m1 );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-3,	uf_.an+l-3,	uf_.anm1+l-3,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l-L4,	uf_.an+l-L4,	uf_.anm1+l-L4,
		    uf_.an  +l-L0,	uf_.an+l-L6,	uf_.an  +l-L2,	uf_.an+l-L8,
		    uf_.an  +l+L0,	uf_.an+l-L7,	uf_.an  +l+L3,	uf_.an+l-L10);
	      }
	  }

	/*循环网格中z = (zmin,zmax)和x = (xmin,xmax)边界上的边缘点和
	*使用一阶吸收边界条件更新场。要了解
	*在代码行之后，最好在旁边列出Li的值。*/
	uf_.af.ufB_ = &uf_.gE[0];
	for (unsigned j = 1; j < uf_.N1m1; j++)
	  {
	    if ( rank_ == 0 )
	      {
		l = 3 * j;

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l+L0,	uf_.an+l+L0,	uf_.anm1+l+L0,
		    uf_.anp1+l+L2,	uf_.an+l+L2,	uf_.anm1+l+L2,
		    uf_.an+l-3,		uf_.an+l-L5,	uf_.an+l-L7,	uf_.an+l-L11,
		    uf_.an+l+3,		uf_.an+l+L4,	uf_.an+l+L6,	uf_.an+l+L8);

		l = 3 * ( N1N0_ - N1_ + j );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l+L1,	uf_.an+l+L1,	uf_.anm1+l+L1,
		    uf_.anp1+l-L0,	uf_.an+l-L0,	uf_.anm1+l-L0,
		    uf_.anp1+l-L3,	uf_.an+l-L3,	uf_.anm1+l-L3,
		    uf_.an+l-3,		uf_.an+l-L5,	uf_.an+l-L6,	uf_.an+l-L9,
		    uf_.an+l+3,		uf_.an+l+L4,	uf_.an+l+L7,	uf_.an+l+L10);
	      }

	    if ( rank_ == size_ - 1 )
	      {
		l = 3 * ( N1N0_ * uf_.npm1 + j );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l+L0,	uf_.an+l+L0,	uf_.anm1+l+L0,
		    uf_.anp1+l+L3,	uf_.an+l+L3,	uf_.anm1+l+L3,
		    uf_.an+l-3,		uf_.an+l-L4,	uf_.an+l-L7,	uf_.an+l-L10,
		    uf_.an+l+3,		uf_.an+l+L5,	uf_.an+l+L6,	uf_.an+l+L9);

		l = 3 * ( N1N0_ * uf_.npm1 + N1_ * uf_.N0m1 + j );

		uf_.af.advanceEdgeF(
		    uf_.anp1+l,		uf_.an+l,	uf_.anm1+l,
		    uf_.anp1+l-L1,	uf_.an+l-L1,	uf_.anm1+l-L1,
		    uf_.anp1+l-L0,	uf_.an+l-L0,	uf_.anm1+l-L0,
		    uf_.anp1+l-L2,	uf_.an+l-L2,	uf_.anm1+l-L2,
		    uf_.an+l-3,		uf_.an+l-L4,	uf_.an+l-L6,	uf_.an+l-L8,
		    uf_.an+l+3,		uf_.an+l+L5,	uf_.an+l+L7,	uf_.an+l+L11);
	      }
	  }

	/*现在更新计算域中八个角的字段。*/
	uf_.af.ufB_ = &uf_.hC[0];
	if ( rank_ == 0 )
	  {
	    m = 0;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m+N1_+1),		uf_.an+3*(m+N1_+1),		uf_.anm1+3*(m+N1_+1),
		uf_.anp1+3*(m+N1N0_+N1_),	uf_.an+3*(m+N1N0_+N1_),		uf_.anm1+3*(m+N1N0_+N1_),
		uf_.anp1+3*(m+N1N0_+1),		uf_.an+3*(m+N1N0_+1),		uf_.anm1+3*(m+N1N0_+1),
		uf_.anp1+3*(m+N1N0_+N1_+1),	uf_.an+3*(m+N1N0_+N1_+1),	uf_.anm1+3*(m+N1N0_+N1_+1));

	    m = N1N0_ - N1_;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m-N1_+1),		uf_.an+3*(m-N1_+1),		uf_.anm1+3*(m-N1_+1),
		uf_.anp1+3*(m+N1N0_-N1_),	uf_.an+3*(m+N1N0_-N1_),		uf_.anm1+3*(m+N1N0_-N1_),
		uf_.anp1+3*(m+N1N0_+1),		uf_.an+3*(m+N1N0_+1),		uf_.anm1+3*(m+N1N0_+1),
		uf_.anp1+3*(m+N1N0_-N1_+1),	uf_.an+3*(m+N1N0_-N1_+1),	uf_.anm1+3*(m+N1N0_-N1_+1));

	    m = uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m+N1_-1),		uf_.an+3*(m+N1_-1),		uf_.anm1+3*(m+N1_-1),
		uf_.anp1+3*(m+N1N0_+N1_),	uf_.an+3*(m+N1N0_+N1_),		uf_.anm1+3*(m+N1N0_+N1_),
		uf_.anp1+3*(m+N1N0_-1),		uf_.an+3*(m+N1N0_-1),		uf_.anm1+3*(m+N1N0_-1),
		uf_.anp1+3*(m+N1N0_+N1_-1),	uf_.an+3*(m+N1N0_+N1_-1),	uf_.anm1+3*(m+N1N0_+N1_-1));

	    m = N1N0_ - N1_ + uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m+N1N0_),		uf_.an+3*(m+N1N0_),		uf_.anm1+3*(m+N1N0_),
		uf_.anp1+3*(m-N1_-1),		uf_.an+3*(m-N1_-1),		uf_.anm1+3*(m-N1_-1),
		uf_.anp1+3*(m+N1N0_-N1_),	uf_.an+3*(m+N1N0_-N1_),		uf_.anm1+3*(m+N1N0_-N1_),
		uf_.anp1+3*(m+N1N0_-1),		uf_.an+3*(m+N1N0_-1),		uf_.anm1+3*(m+N1N0_-1),
		uf_.anp1+3*(m+N1N0_-N1_-1),	uf_.an+3*(m+N1N0_-N1_-1),	uf_.anm1+3*(m+N1N0_-N1_-1));
	  }

	if ( rank_ == size_ - 1 )
	  {
	    m = N1N0_ * uf_.npm1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m+N1_+1),		uf_.an+3*(m+N1_+1),		uf_.anm1+3*(m+N1_+1),
		uf_.anp1+3*(m-N1N0_+N1_),	uf_.an+3*(m-N1N0_+N1_),		uf_.anm1+3*(m-N1N0_+N1_),
		uf_.anp1+3*(m-N1N0_+1),		uf_.an+3*(m-N1N0_+1),		uf_.anm1+3*(m-N1N0_+1),
		uf_.anp1+3*(m-N1N0_+N1_+1),	uf_.an+3*(m-N1N0_+N1_+1),	uf_.anm1+3*(m-N1N0_+N1_+1));

	    m = N1N0_ * uf_.npm1 + N1N0_ - N1_;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m+1),		uf_.an+3*(m+1),			uf_.anm1+3*(m+1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m-N1_+1),		uf_.an+3*(m-N1_+1),		uf_.anm1+3*(m-N1_+1),
		uf_.anp1+3*(m-N1N0_-N1_),	uf_.an+3*(m-N1N0_-N1_),		uf_.anm1+3*(m-N1N0_-N1_),
		uf_.anp1+3*(m-N1N0_+1),		uf_.an+3*(m-N1N0_+1),		uf_.anm1+3*(m-N1N0_+1),
		uf_.anp1+3*(m-N1N0_-N1_+1),	uf_.an+3*(m-N1N0_-N1_+1),	uf_.anm1+3*(m-N1N0_-N1_+1));

	    m = N1N0_ * uf_.npm1 + uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m+N1_),		uf_.an+3*(m+N1_),		uf_.anm1+3*(m+N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m+N1_-1),		uf_.an+3*(m+N1_-1),		uf_.anm1+3*(m+N1_-1),
		uf_.anp1+3*(m-N1N0_+N1_),	uf_.an+3*(m-N1N0_+N1_),		uf_.anm1+3*(m-N1N0_+N1_),
		uf_.anp1+3*(m-N1N0_-1),		uf_.an+3*(m-N1N0_-1),		uf_.anm1+3*(m-N1N0_-1),
		uf_.anp1+3*(m-N1N0_+N1_-1),	uf_.an+3*(m-N1N0_+N1_-1),	uf_.anm1+3*(m-N1N0_+N1_-1));

	    m = N1N0_ * uf_.npm1 + N1N0_ - N1_ + uf_.N1m1;
	    uf_.af.advanceCornerF(
		uf_.anp1+3*m,			uf_.an+3*m,			uf_.anm1+3*m,
		uf_.anp1+3*(m-N1_),		uf_.an+3*(m-N1_),		uf_.anm1+3*(m-N1_),
		uf_.anp1+3*(m-1),		uf_.an+3*(m-1),			uf_.anm1+3*(m-1),
		uf_.anp1+3*(m-N1N0_),		uf_.an+3*(m-N1N0_),		uf_.anm1+3*(m-N1N0_),
		uf_.anp1+3*(m-N1_-1),		uf_.an+3*(m-N1_-1),		uf_.anm1+3*(m-N1_-1),
		uf_.anp1+3*(m-N1N0_-N1_),	uf_.an+3*(m-N1N0_-N1_),		uf_.anm1+3*(m-N1N0_-N1_),
		uf_.anp1+3*(m-N1N0_-1),		uf_.an+3*(m-N1N0_-1),		uf_.anm1+3*(m-N1N0_-1),
		uf_.anp1+3*(m-N1N0_-N1_-1),	uf_.an+3*(m-N1N0_-N1_-1),	uf_.anm1+3*(m-N1N0_-N1_-1));
	  }
      }

    /*在整个处理器中通信计算字段。*/
    if (rank_ != size_ - 1)
      MPI_Send(uf_.anp1+3*(np_-2)*N1N0_, 	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag1,MPI_COMM_WORLD);

    if (rank_ != 0)
      MPI_Recv(uf_.anp1,			3*N1N0_,MPI_DOUBLE,rank_-1,msgtag1,MPI_COMM_WORLD,&status);

    if (rank_ != 0)
      MPI_Send(uf_.anp1+3*N1N0_,	 	3*N1N0_,MPI_DOUBLE,rank_-1,msgtag3,MPI_COMM_WORLD);

    if (rank_ != size_ - 1)
      MPI_Recv(uf_.anp1+3*(np_-1)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag3,MPI_COMM_WORLD,&status);

    /*现在更新了A和phi量，计算边界网格点的E和B以供稍后使用
	*加速度和功率测量。*/
    for (unsigned i = 1; i < uf_.N0m1; i++)
      for (unsigned j = 1; j < uf_.N1m1; j++)
	{
	  m = N1N0_ + N1_ * i + j;

	  /*计算这个像素处的场。*/
	  fieldEvaluate(m);

	  /*为这个像素和它之前的像素设置布尔标志为true。*/
	  pic_[m-N1N0_] = true;

	  /*对于左边界（z=zmin），只需将场设置为下一个z平面。*/
	  if ( rank_ == 0 )
	    {
	      en_[m-N1N0_] = en_[m];
	      bn_[m-N1N0_] = bn_[m];
	    }

	  m = N1N0_ * ( np_ - 2 ) + N1_ * i + j;

	  /*计算这个像素处的场。*/
	  fieldEvaluate(m);

	  /*为这个像素和它前面的像素设置布尔标志为true。*/
	  pic_[m+N1N0_] = true;

	  /*对于右边界（z=zmax），只需将字段设置为与之前的z平面相等。*/
	  if ( rank_ == size_ - 1 )
	    {
	      en_[m+N1N0_] = en_[m];
	      bn_[m+N1N0_] = bn_[m];
	    }
	}

    /*在整个处理器中通信计算字段。*/
    if (rank_ != size_ - 1)
      {
	MPI_Send(uf_.en+3*(np_-2)*N1N0_, 	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag5,MPI_COMM_WORLD);
	MPI_Send(uf_.bn+3*(np_-2)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag6,MPI_COMM_WORLD);
      }

    if (rank_ != 0)
      {
	MPI_Recv(uf_.en,			3*N1N0_,MPI_DOUBLE,rank_-1,msgtag5,MPI_COMM_WORLD,&status);
	MPI_Recv(uf_.bn,		  	3*N1N0_,MPI_DOUBLE,rank_-1,msgtag6,MPI_COMM_WORLD,&status);
      }

    if (rank_ != 0)
      {
	MPI_Send(uf_.en+3*N1N0_,	 	3*N1N0_,MPI_DOUBLE,rank_-1,msgtag7,MPI_COMM_WORLD);
	MPI_Send(uf_.bn+3*N1N0_,		3*N1N0_,MPI_DOUBLE,rank_-1,msgtag8,MPI_COMM_WORLD);
      }

    if (rank_ != size_ - 1)
      {
	MPI_Recv(uf_.en+3*(np_-1)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag7,MPI_COMM_WORLD,&status);
	MPI_Recv(uf_.bn+3*(np_-1)*N1N0_,	3*N1N0_,MPI_DOUBLE,rank_+1,msgtag8,MPI_COMM_WORLD,&status);
      }

	/* 整个场步完成后，再把 CPML 记忆变量从 n+1 滚到 n */
	if (cpmlXY_.enabled) shiftCPMLXY();
  }

  /******************************************************************************************************
  *移动计算字段和字段的时间点。
  ******************************************************************************************************/

  void FdTd::fieldShift ()
  {
    std::vector<FieldVector<Double> >* at = anm1_;
    anm1_ = an_;
    an_   = anp1_;
    anp1_ = at;
  }

  /******************************************************************************************************
  *从电位中评估第m个像素的场。
  ******************************************************************************************************/

  void FdTd::fieldEvaluate (long int m)
  {
    /*计算电场。*/
    en_[m].dv ( - uf_.dt, (*anp1_)[m] );
    en_[m].mdv( - uf_.dt, (*an_)  [m] );

    /*计算磁场。*/
    bn_[m][0] = 0.5 * (
	( *(uf_.an  +3*(m+1    )+2 ) - *(uf_.an  +3*(m-1    )+2  ) ) / uf_.dy2 -
	( *(uf_.an  +3*(m+N1N0_)+1 ) - *(uf_.an  +3*(m-N1N0_)+1  ) ) / uf_.dz2 +
	( *(uf_.anp1+3*(m+1    )+2 ) - *(uf_.anp1+3*(m-1    )+2  ) ) / uf_.dy2 -
	( *(uf_.anp1+3*(m+N1N0_)+1 ) - *(uf_.anp1+3*(m-N1N0_)+1  ) ) / uf_.dz2 );

    bn_[m][1] = 0.5 * (
	( *(uf_.an  +3*(m+N1N0_)   ) - *(uf_.an  +3*(m-N1N0_)    ) ) / uf_.dz2 -
	( *(uf_.an  +3*(m+N1_  )+2 ) - *(uf_.an  +3*(m-N1_  )+2  ) ) / uf_.dx2 +
	( *(uf_.anp1+3*(m+N1N0_)   ) - *(uf_.anp1+3*(m-N1N0_)    ) ) / uf_.dz2 -
	( *(uf_.anp1+3*(m+N1_  )+2 ) - *(uf_.anp1+3*(m-N1_  )+2  ) ) / uf_.dx2 );

    bn_[m][2] = 0.5 * (
	( *(uf_.an  +3*(m+N1_  )+1 ) - *(uf_.an  +3*(m-N1_  )+1  ) ) / uf_.dx2 -
	( *(uf_.an  +3*(m+1    )   ) - *(uf_.an  +3*(m-1    )    ) ) / uf_.dy2 +
	( *(uf_.anp1+3*(m+N1_  )+1 ) - *(uf_.anp1+3*(m-N1_  )+1  ) ) / uf_.dx2 -
	( *(uf_.anp1+3*(m+1    )   ) - *(uf_.anp1+3*(m-1    )    ) ) / uf_.dy2 );

    /*将该像素的布尔标志设置为true。*/
    pic_[m] = true;
  }

  /******************************************************************************************************
  *采样字段并将其保存到给定的文件。
  ******************************************************************************************************/

  void FdTd::fieldSample ()
  {

    /*当且仅当采样点位于此处理器范围内时，执行字段采样。*/
    if ( sf_.N > 0 )
      {
	( *(sf_.file) ).setf(std::ios::scientific);
	( *(sf_.file) ).precision(4);

	/*把时间写在第一列。*/
	*(sf_.file) << time_ * gamma_ << "\t";

	for (unsigned int n = 0; n < sf_.N; ++n)
	  {
	    /*得到采样点的位置。*/
	    sf_.position 	= seed_.samplingPosition_[n];

	    /*得到采样点的指数。*/
	    sf_.dxr = modf( ( sf_.position[0] - xmin_ ) / mesh_.meshResolution_[0] , &sf_.c1);
	    sf_.i   = (int) sf_.c1;
	    sf_.dyr = modf( ( sf_.position[1] - ymin_ ) / mesh_.meshResolution_[1] , &sf_.c1);
	    sf_.j   = (int) sf_.c1;
	    sf_.dzr = modf( ( sf_.position[2] - zmin_ ) / mesh_.meshResolution_[2] , &sf_.c1);
	    sf_.k   = (int) sf_.c1;
	    sf_.m   = ( sf_.k - k0_ ) * N1N0_ + sf_.i * N1_ + sf_.j;

	    /*计算字段以找到采样点处的值。*/
	    if (!pic_[sf_.m            ])     fieldEvaluate(sf_.m            );
	    if (!pic_[sf_.m+N1_        ])     fieldEvaluate(sf_.m+N1_        );
	    if (!pic_[sf_.m+1          ])     fieldEvaluate(sf_.m+1          );
	    if (!pic_[sf_.m+N1_+1      ])     fieldEvaluate(sf_.m+N1_+1      );
	    if (!pic_[sf_.m+N1N0_      ])     fieldEvaluate(sf_.m+N1N0_      );
	    if (!pic_[sf_.m+N1N0_+N1_  ])     fieldEvaluate(sf_.m+N1N0_+N1_  );
	    if (!pic_[sf_.m+N1N0_+1    ])     fieldEvaluate(sf_.m+N1N0_+1    );
	    if (!pic_[sf_.m+N1N0_+N1_+1])     fieldEvaluate(sf_.m+N1N0_+N1_+1);

	    /*计算并插值电场，求出采样点处的值。*/
	    sf_.et.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), en_[sf_.m]);
	    sf_.et.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), en_[sf_.m+N1_]);
	    sf_.et.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), en_[sf_.m+1]);
	    sf_.et.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), en_[sf_.m+N1_+1]);
	    sf_.et.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         en_[sf_.m+N1N0_]);
	    sf_.et.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         en_[sf_.m+N1N0_+N1_]);
	    sf_.et.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         en_[sf_.m+N1N0_+1]);
	    sf_.et.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         en_[sf_.m+N1N0_+N1_+1]);

	    /*计算并插值磁场，求其在采样点处的值。*/
	    sf_.bt.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), bn_[sf_.m]);
	    sf_.bt.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), bn_[sf_.m+N1_]);
	    sf_.bt.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), bn_[sf_.m+1]);
	    sf_.bt.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), bn_[sf_.m+N1_+1]);
	    sf_.bt.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         bn_[sf_.m+N1N0_]);
	    sf_.bt.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         bn_[sf_.m+N1N0_+N1_]);
	    sf_.bt.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         bn_[sf_.m+N1N0_+1]);
	    sf_.bt.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         bn_[sf_.m+N1N0_+N1_+1]);

	    sf_.at.mv ((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), (*an_)[sf_.m]);
	    sf_.at.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * (1.0 - sf_.dzr), (*an_)[sf_.m+N1_]);
	    sf_.at.pmv((1.0 - sf_.dxr) * sf_.dyr           * (1.0 - sf_.dzr), (*an_)[sf_.m+1]);
	    sf_.at.pmv(sf_.dxr         * sf_.dyr           * (1.0 - sf_.dzr), (*an_)[sf_.m+N1_+1]);
	    sf_.at.pmv((1.0 - sf_.dxr) * (1.0 - sf_.dyr)   * sf_.dzr,         (*an_)[sf_.m+N1N0_]);
	    sf_.at.pmv(sf_.dxr         * (1.0 - sf_.dyr)   * sf_.dzr,         (*an_)[sf_.m+N1N0_+N1_]);
	    sf_.at.pmv((1.0 - sf_.dxr) * sf_.dyr           * sf_.dzr,         (*an_)[sf_.m+N1N0_+1]);
	    sf_.at.pmv(sf_.dxr         * sf_.dyr           * sf_.dzr,         (*an_)[sf_.m+N1N0_+N1_+1]);

	    /*写下下一列的坐标。*/
	    *(sf_.file) << sf_.position[0] << "\t";
	    *(sf_.file) << sf_.position[1] << "\t";
	    *(sf_.file) << sf_.position[2] << "\t";

	    /*写下下一列中的字段。*/
	    for (unsigned int i = 0; i < seed_.samplingField_.size(); i++)
	      {
		if 		( seed_.samplingField_[i] == Ex )
		  *(sf_.file) << ( gamma_ * sf_.et[0] + c0_ * sqrt( pow(gamma_, 2) - 1 ) * sf_.bt[1] ) * sf_.Ce << "\t";
		else if 	( seed_.samplingField_[i] == Ey )
		  *(sf_.file) << ( gamma_ * sf_.et[1] - c0_ * sqrt( pow(gamma_, 2) - 1 ) * sf_.bt[0] ) * sf_.Ce << "\t";
		else if	( seed_.samplingField_[i] == Ez )
		  *(sf_.file) << sf_.et[2] * sf_.Ce << "\t";

		else if	( seed_.samplingField_[i] == Bx )
		  *(sf_.file) << ( gamma_ * sf_.bt[0] - sqrt( pow(gamma_, 2) - 1 ) / c0_ * sf_.et[1] ) * sf_.Cb << "\t";
		else if	( seed_.samplingField_[i] == By )
		  *(sf_.file) << ( gamma_ * sf_.bt[1] + sqrt( pow(gamma_, 2) - 1 ) / c0_ * sf_.et[0] ) * sf_.Cb << "\t";
		else if	( seed_.samplingField_[i] == Bz )
		  *(sf_.file) << sf_.bt[2] * sf_.Cb << "\t";

		else if	( seed_.samplingField_[i] == Ax )
		  *(sf_.file) << sf_.at[0] * sf_.Ca << "\t";
		else if	( seed_.samplingField_[i] == Ay )
		  *(sf_.file) << sf_.at[1] * sf_.Ca << "\t";
		else if	( seed_.samplingField_[i] == Az )
		  *(sf_.file) << sf_.at[2] * sf_.Ca << "\t";
	      }
	  }

	/*将文件光标移到下一行。*/
	*(sf_.file) << std::endl;
      }
  }

  /******************************************************************************************************
  *可视化的领域作为vtk文件在整个领域，并将其保存到文件与给定的名称。
  ******************************************************************************************************/

  void FdTd::fieldVisualizeAllDomain (unsigned int ivtk)
  {
    long int			m;

    /*旧文件如果存在，应该删除。*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /*计算要在vtk文件中可视化的字段。*/
    for (int k = 0; k < np_; k++ )
      for (int j = 1; j < N1_-1; j++)
	for (int i = 1; i < N0_-1; i++)
	  {
	    m = k * N1_ * N0_ + i * N1_ + j;

	    if (!pic_[m]) fieldEvaluate(m);

	    for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	      {
		if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[m][l] = en_[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[m][l] = en_[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[m][l] = en_[m][2];
		else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[m][l] = bn_[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[m][l] = bn_[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[m][l] = bn_[m][2];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[m][l] = (*an_)[m][0];
		else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[m][l] = (*an_)[m][1];
		else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[m][l] = (*an_)[m][2];
	      }
	  }

    /*为vtk文件写入初始数据。*/
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;

    /*插入充电点的网格坐标。*/
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    FieldVector<Double> r (0.0);
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++)
      for (int j = 0; j < N1_; j++)
	for (int i = 0; i < N0_; i++)
	  {
	    m = k * N1_ * N0_ + i * N1_ + j;
	    r = rc(m);
	    *vf_[ivtk].file << r[0] << " " << r[1] << " " << r[2] 					<< std::endl;
	  }
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /*将每个单元格数据插入到vtk文件中。*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          	<< std::endl;

    /*根据计算的电场插入点数据。*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    		<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" format=\"ascii\">"
	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ) ; k++ )
      for (int j = 0; j < N1_; j++)
	for (int i = 0; i < N0_; i++)
	  {
	    m = k * N1_ * N0_ + i * N1_ + j;
	    *vf_[ivtk].file << vf_[ivtk].v[m][0];
	    for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[m][l];
	    *vf_[ivtk].file << std::endl;
	  }
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    	<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /*关闭文件。*/
    (*vf_[ivtk].file).close();

    /*写入连接并行文件的文件。*/

    if ( rank_ == 0)
      {
	/*在文件名后加上vtk后缀和vtk文件编号。*/
	vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-" + stringify(nTime_) + PTS_FILE_SUFFIX;
	unsigned int k0, np;

	/*假设在运行代码之前删除了目录中的所有文件。*/
	vf_[ivtk].file = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

	*vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
	*vf_[ivtk].file << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" >"			<< std::endl;
	*vf_[ivtk].file << "<PStructuredGrid WholeExtent=\"0 " << N0_-1 << " 0 " << N1_-1 << " 0 "  <<
	    N2_-1 << "\" GhostLevel = \"0\" >"                                    			<< std::endl;

	/*插入电荷云的网格坐标。*/
	*vf_[ivtk].file << "<PPoints>"                                                          	<< std::endl;
	*vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"	<< std::endl;
	*vf_[ivtk].file << "</PPoints>"                                                             	<< std::endl;

	*vf_[ivtk].file << "<PPointData>"                                                           	<< std::endl;
	*vf_[ivtk].file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" Name=\"field\" format=\"ascii\" />"
	    << std::endl;
	*vf_[ivtk].file << "</PPointData>"                                                          	<< std::endl;

	for (int i = 0; i < size_; ++i)
	  {
	    /*评估每个处理器中的节点数量。*/
	    if ( size_ > 1 )
	      {
		if ( i == 0 )
		  {
		    np = N2_ / size_ + 1;
		    k0 = 0;
		  }
		else if ( i == size_ - 1 )
		  {
		    np = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 1;
		    k0 = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
		  }
		else
		  {
		    np = N2_ / size_ + 2;
		    k0 = i * ( N2_ / size_ ) - 1;
		  }
	      }
	    else
	      {
		np = N2_;
		k0 = 0;
	      }

	    vf_[ivtk].fileName = vf_[ivtk].name + "-p" + stringify(i) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
	    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_-1 << " 0 " << N1_-1 << " " <<
		k0 << " " << k0 + np - 2 + ( (i == size_ - 1) ? 1 : 0 ) << "\""
		<< " Source=\"" << vf_[ivtk].fileName << "\" />"                       		<< std::endl;
	  }
	*vf_[ivtk].file << "</PStructuredGrid>"                                                     	<< std::endl;
	*vf_[ivtk].file << "</VTKFile>"                                                             	<< std::endl;

	/*关闭文件。*/
	(*vf_[ivtk].file).close();
      }
  }

  /******************************************************************************************************
  *可视化领域的vtk文件在平面上，并将其保存到文件与给定的名称。
  ******************************************************************************************************/

  void FdTd::fieldVisualizeInPlane (unsigned int ivtk)
  {
    if	( seed_.vtk_[ivtk].plane_ == XNORMAL )  fieldVisualizeInPlaneXNormal(ivtk);
    else if   ( seed_.vtk_[ivtk].plane_ == YNORMAL )  fieldVisualizeInPlaneYNormal(ivtk);
    else if   ( seed_.vtk_[ivtk].plane_ == ZNORMAL )
      {
	if ( seed_.vtk_[ivtk].position_[2] < zp_[1] && seed_.vtk_[ivtk].position_[2] >= zp_[0] )
	  fieldVisualizeInPlaneZNormal(ivtk);
      }
  }

  /******************************************************************************************************
  *将字段可视化为vtk文件，在平面上垂直于x轴，并将它们保存到文件中
  *名字。
  ******************************************************************************************************/

  void FdTd::fieldVisualizeInPlaneXNormal (unsigned int ivtk)
  {
    unsigned int		i, n;
    long int			m;
    Double			dxr, c;

    /*旧文件如果存在，应该删除。*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /*计算平面所在单元格的索引。*/
    dxr = modf( ( seed_.vtk_[ivtk].position_[0] - xmin_ ) / mesh_.meshResolution_[0] , &c);
    i   = (int) c;

    /*计算要在vtk文件中可视化的字段。*/
    for (int k = 0; k < np_; k++ )
      for (int j = 1; j < N1_-1; j++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  n = k * N1_ + j;

	  if (!pic_[m]) 		fieldEvaluate(m);
	  if (!pic_[m + N1_]) 	fieldEvaluate(m + N1_);

	  for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	    {
	      if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[n][l] = en_[m][0] * ( 1.0 - dxr ) + en_[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[n][l] = en_[m][1] * ( 1.0 - dxr ) + en_[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[n][l] = en_[m][2] * ( 1.0 - dxr ) + en_[m + N1_][2] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[n][l] = bn_[m][0] * ( 1.0 - dxr ) + bn_[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[n][l] = bn_[m][1] * ( 1.0 - dxr ) + bn_[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[n][l] = bn_[m][2] * ( 1.0 - dxr ) + bn_[m + N1_][2] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[n][l] = (*an_)[m][0] * ( 1.0 - dxr ) + (*an_)[m + N1_][0] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[n][l] = (*an_)[m][1] * ( 1.0 - dxr ) + (*an_)[m + N1_][1] * dxr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[n][l] = (*an_)[m][2] * ( 1.0 - dxr ) + (*an_)[m + N1_][2] * dxr;
	    }
	}

    /*为vtk文件写入初始数据。*/
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0  0  0 " << N1_ - 1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\" 0 0 0 " << N1_-1 << " " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;

    /*插入充电点的网格坐标。*/
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    FieldVector<Double> r1 (0.0), r2 (0.0);
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++)
      for (int j = 0; j < N1_; j++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  r1 = rc(m); r2 = rc(m + N1_);
	  *vf_[ivtk].file << r1[0] * ( 1.0 - dxr ) + r2[0] * dxr << " " << r1[1] << " " << r1[2] 	<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /*将每个单元格数据插入到vtk文件中。*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          	<< std::endl;

    /*根据计算的电场插入点数据。*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    	<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" format=\"ascii\">"
	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ) ; k++ )
      for (int j = 0; j < N1_; j++)
	{
	  n = k * N1_ + j;
	  *vf_[ivtk].file << vf_[ivtk].v[n][0];
	  for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[n][l];
	  *vf_[ivtk].file << std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    	<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /*关闭文件。*/
    (*vf_[ivtk].file).close();

    /*写入连接并行文件的文件。*/

    if ( rank_ == 0)
      {
	/*在文件名后加上vtk后缀和vtk文件编号。*/
	vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-" + stringify(nTime_) + PTS_FILE_SUFFIX;
	unsigned int k0, np;

	/*假设在运行代码之前删除了目录中的所有文件。*/
	vf_[ivtk].file = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

	*vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
	*vf_[ivtk].file << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" >"			<< std::endl;
	*vf_[ivtk].file << "<PStructuredGrid WholeExtent=\" 0 0 0 " << N1_-1 << " 0 "  <<
	    N2_-1 << "\" GhostLevel = \"0\" >"                                    			<< std::endl;

	/*插入电荷云的网格坐标。*/
	*vf_[ivtk].file << "<PPoints>"                                                          	<< std::endl;
	*vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"	<< std::endl;
	*vf_[ivtk].file << "</PPoints>"                                                             	<< std::endl;

	*vf_[ivtk].file << "<PPointData>"                                                           	<< std::endl;
	*vf_[ivtk].file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" Name=\"field\" format=\"ascii\" />"
	    << std::endl;
	*vf_[ivtk].file << "</PPointData>"                                                          	<< std::endl;

	for (int i = 0; i < size_; ++i)
	  {
	    /*评估每个处理器中的节点数量。*/
	    if ( size_ > 1 )
	      {
		if ( i == 0 )
		  {
		    np = N2_ / size_ + 1;
		    k0 = 0;
		  }
		else if ( i == size_ - 1 )
		  {
		    np = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 1;
		    k0 = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
		  }
		else
		  {
		    np = N2_ / size_ + 2;
		    k0 = i * ( N2_ / size_ ) - 1;
		  }
	      }
	    else
	      {
		np = N2_;
		k0 = 0;
	      }

	    vf_[ivtk].fileName = vf_[ivtk].name + "-p" + stringify(i) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
	    *vf_[ivtk].file << "<Piece Extent=\"0 0 0 " << N1_-1 << " " <<
		k0 << " " << k0 + np - 2 + ( (i == size_ - 1) ? 1 : 0 ) << "\""
		<< " Source=\"" << vf_[ivtk].fileName << "\" />"                       		<< std::endl;
	  }
	*vf_[ivtk].file << "</PStructuredGrid>"                                                     	<< std::endl;
	*vf_[ivtk].file << "</VTKFile>"                                                             	<< std::endl;

	/*关闭文件。*/
	(*vf_[ivtk].file).close();
      }
  }

  /******************************************************************************************************
  *将字段可视化为垂直于y轴的平面中的vtk文件，并将它们保存到文件中
  *名字。
  ******************************************************************************************************/

  void FdTd::fieldVisualizeInPlaneYNormal (unsigned int ivtk)
  {
    unsigned int		j, n;
    long int			m;
    Double			dyr, c;

    /*旧文件如果存在，应该删除。*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /*计算平面所在单元格的索引。*/
    dyr = modf( ( seed_.vtk_[ivtk].position_[1] - ymin_ ) / mesh_.meshResolution_[1] , &c);
    j   = (int) c;

    /*计算要在vtk文件中可视化的字段。*/
    for (int k = 0; k < np_; k++ )
      for (int i = 1; i < N0_-1; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  n = k * N0_ + i;

	  if (!pic_[m]) 	fieldEvaluate(m);
	  if (!pic_[m + 1]) 	fieldEvaluate(m + 1);

	  for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	    {
	      if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[n][l] = en_[m][0] * ( 1.0 - dyr ) + en_[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[n][l] = en_[m][1] * ( 1.0 - dyr ) + en_[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[n][l] = en_[m][2] * ( 1.0 - dyr ) + en_[m + 1][2] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[n][l] = bn_[m][0] * ( 1.0 - dyr ) + bn_[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[n][l] = bn_[m][1] * ( 1.0 - dyr ) + bn_[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[n][l] = bn_[m][2] * ( 1.0 - dyr ) + bn_[m + 1][2] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[n][l] = (*an_)[m][0] * ( 1.0 - dyr ) + (*an_)[m + 1][0] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[n][l] = (*an_)[m][1] * ( 1.0 - dyr ) + (*an_)[m + 1][1] * dyr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[n][l] = (*an_)[m][2] * ( 1.0 - dyr ) + (*an_)[m + 1][2] * dyr;
	    }
	}

    /*为vtk文件写入初始数据。*/
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 0 " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 0 " <<
	k0_ << " " << k0_ + np_ - 2 + ( (rank_ == size_ - 1) ? 1 : 0 )
	<< "\">"											<< std::endl;

    /*插入充电点的网格坐标。*/
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    FieldVector<Double> r1 (0.0), r2 (0.0);
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++)
      for (int i = 0; i < N0_; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  r1 = rc(m); r2 = rc(m + 1);
	  *vf_[ivtk].file << r1[0] * ( 1.0 - dyr ) + r2[0] * dyr << " " << r1[1] << " " << r1[2] 	<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /*将每个单元格数据插入到vtk文件中。*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          	<< std::endl;

    /*根据计算的电场插入点数据。*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    	<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" format=\"ascii\">"
	<< std::endl;
    for (int k = 0 ; k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ) ; k++ )
      for (int i = 0; i < N0_; i++)
	{
	  n = k * N0_ + i;
	  *vf_[ivtk].file << vf_[ivtk].v[n][0];
	  for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[n][l];
	  *vf_[ivtk].file << std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    	<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /*关闭文件。*/
    (*vf_[ivtk].file).close();

    /*写入连接并行文件的文件。*/

    if ( rank_ == 0)
      {
	/*在文件名后加上vtk后缀和vtk文件编号。*/
	vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-" + stringify(nTime_) + PTS_FILE_SUFFIX;
	unsigned int k0, np;

	/*假设在运行代码之前删除了目录中的所有文件。*/
	vf_[ivtk].file = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

	*vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
	*vf_[ivtk].file << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" >"			<< std::endl;
	*vf_[ivtk].file << "<PStructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 0 0 "  <<
	    N2_-1 << "\" GhostLevel = \"0\" >"                                    			<< std::endl;

	/*插入电荷云的网格坐标。*/
	*vf_[ivtk].file << "<PPoints>"                                                          	<< std::endl;
	*vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"	<< std::endl;
	*vf_[ivtk].file << "</PPoints>"                                                             	<< std::endl;

	*vf_[ivtk].file << "<PPointData>"                                                           	<< std::endl;
	*vf_[ivtk].file << "<DataArray type=\"Float64\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size() << "\" Name=\"field\" format=\"ascii\" />"
	    << std::endl;
	*vf_[ivtk].file << "</PPointData>"                                                          	<< std::endl;

	for (int i = 0; i < size_; ++i)
	  {
	    /*评估每个处理器中的节点数量。*/
	    if ( size_ > 1 )
	      {
		if ( i == 0 )
		  {
		    np = N2_ / size_ + 1;
		    k0 = 0;
		  }
		else if ( i == size_ - 1 )
		  {
		    np = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 1;
		    k0 = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
		  }
		else
		  {
		    np = N2_ / size_ + 2;
		    k0 = i * ( N2_ / size_ ) - 1;
		  }
	      }
	    else
	      {
		np = N2_;
		k0 = 0;
	      }

	    vf_[ivtk].fileName = vf_[ivtk].name + "-p" + stringify(i) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
	    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 0 " <<
		k0 << " " << k0 + np - 2 + ( (i == size_ - 1) ? 1 : 0 ) << "\""
		<< " Source=\"" << vf_[ivtk].fileName << "\" />"                       			<< std::endl;
	  }
	*vf_[ivtk].file << "</PStructuredGrid>"                                                     	<< std::endl;
	*vf_[ivtk].file << "</VTKFile>"                                                             	<< std::endl;

	/*关闭文件。*/
	(*vf_[ivtk].file).close();
      }
  }

  /******************************************************************************************************
  *将字段可视化为垂直于z轴的平面上的vtk文件，并将它们保存到文件中
  *名字。
  ******************************************************************************************************/

  void FdTd::fieldVisualizeInPlaneZNormal (unsigned int ivtk)
  {
    unsigned int	        k, n;
    long int			m;
    Double		        dzr, c;

    /*旧文件如果存在，应该删除。*/
    vf_[ivtk].fileName = seed_.vtk_[ivtk].basename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + VTS_FILE_SUFFIX;
    (vf_[ivtk].file) = new std::ofstream(vf_[ivtk].fileName.c_str(),std::ios::trunc);

    vf_[ivtk].file->setf(std::ios::scientific);
    vf_[ivtk].file->precision(4);

    /*计算平面所在单元格的索引。*/
    dzr = modf( ( seed_.vtk_[ivtk].position_[2] - zmin_ ) / mesh_.meshResolution_[2] , &c);
    k   = (int) c - k0_;

    /*计算要在vtk文件中可视化的字段。*/
    for (int j = 1; j < N1_-1; j++)
      for (int i = 1; i < N0_-1; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  n = i * N1_ + j;

	  if (!pic_[m]) 		fieldEvaluate(m);
	  if (!pic_[m + N1N0_]) 	fieldEvaluate(m + N1N0_);

	  for (unsigned l = 0; l < seed_.vtk_[ivtk].field_.size(); l++ )
	    {
	      if 		( seed_.vtk_[ivtk].field_[l] == Ex )	vf_[ivtk].v[n][l] = en_[m][0] * ( 1.0 - dzr ) + en_[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ey )	vf_[ivtk].v[n][l] = en_[m][1] * ( 1.0 - dzr ) + en_[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ez )	vf_[ivtk].v[n][l] = en_[m][2] * ( 1.0 - dzr ) + en_[m + N1N0_][2] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bx )	vf_[ivtk].v[n][l] = bn_[m][0] * ( 1.0 - dzr ) + bn_[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == By )	vf_[ivtk].v[n][l] = bn_[m][1] * ( 1.0 - dzr ) + bn_[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Bz )	vf_[ivtk].v[n][l] = bn_[m][2] * ( 1.0 - dzr ) + bn_[m + N1N0_][2] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ax )	vf_[ivtk].v[n][l] = (*an_)[m][0] * ( 1.0 - dzr ) + (*an_)[m + N1N0_][0] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Ay )	vf_[ivtk].v[n][l] = (*an_)[m][1] * ( 1.0 - dzr ) + (*an_)[m + N1N0_][1] * dzr;
	      else if 	( seed_.vtk_[ivtk].field_[l] == Az )	vf_[ivtk].v[n][l] = (*an_)[m][2] * ( 1.0 - dzr ) + (*an_)[m + N1N0_][2] * dzr;
	    }
	}

    /*为vtk文件写入初始数据。*/
    *vf_[ivtk].file << "<?xml version=\"1.0\"?>"							<< std::endl;
    *vf_[ivtk].file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
	"compressor=\"vtkZLibDataCompressor\">" 							<< std::endl;
    *vf_[ivtk].file << "<StructuredGrid WholeExtent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	0 << " " << 0 << "\">"									<< std::endl;
    *vf_[ivtk].file << "<Piece Extent=\"0 " << N0_ - 1 << " 0 " << N1_ - 1 << " " <<
	0 << " " << 0 << "\">"									<< std::endl;

    /*插入充电点的网格坐标。*/
    *vf_[ivtk].file << "<Points>"                                                                	<< std::endl;
    *vf_[ivtk].file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"	<< std::endl;
    FieldVector<Double> r1 (0.0), r2 (0.0);
    for (int j = 0; j < N1_; j++)
      for (int i = 0; i < N0_; i++)
	{
	  m = k * N1_ * N0_ + i * N1_ + j;
	  r1 = rc(m); r2 = rc(m + N1N0_);
	  *vf_[ivtk].file << r1[0] * ( 1.0 - dzr ) + r2[0] * dzr << " " << r1[1] << " " << r1[2] 	<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                       		<< std::endl;
    *vf_[ivtk].file << "</Points>"                                                         		<< std::endl;

    /*将每个单元格数据插入到vtk文件中。*/
    *vf_[ivtk].file << "<CellData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</CellData>"                                                          	<< std::endl;

    /*根据计算的电场插入点数据。*/
    *vf_[ivtk].file << "<PointData Vectors = \"field\">"                                    		<< std::endl;
    *vf_[ivtk].file << "<DataArray type=\"Float64\" Name=\"field\" NumberOfComponents=\"" << seed_.vtk_[ivtk].field_.size()
		      << "\" format=\"ascii\">"									<< std::endl;
    for (int j = 0; j < N1_; j++)
      for (int i = 0; i < N0_; i++)
	{
	  n = i * N1_ + j;
	  *vf_[ivtk].file << vf_[ivtk].v[n][0];
	  for (unsigned l = 1; l < seed_.vtk_[ivtk].field_.size(); l++) *vf_[ivtk].file << " " << vf_[ivtk].v[n][l];
	  *vf_[ivtk].file 										<< std::endl;
	}
    *vf_[ivtk].file << "</DataArray>"                                                   		<< std::endl;
    *vf_[ivtk].file << "</PointData>"                                                        		<< std::endl;
    *vf_[ivtk].file << "</Piece>"                                                            		<< std::endl;
    *vf_[ivtk].file << "</StructuredGrid>"                                                    	<< std::endl;
    *vf_[ivtk].file << "</VTKFile>"                                                         		<< std::endl;

    /*关闭文件。*/
    (*vf_[ivtk].file).close();
  }

  /******************************************************************************************************
  *将字段的总概要文件写入给定的文件名中。
  ******************************************************************************************************/

  void FdTd::fieldProfile ()
  {
    /*声明点上循环的迭代器。*/
    pf_.fileName = seed_.profileBasename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + TXT_FILE_SUFFIX;
    pf_.file = new std::ofstream(pf_.fileName.c_str(),std::ios::trunc);

    pf_.file->setf(std::ios::scientific);
    pf_.file->precision(4);

    /*在网格点上执行循环，并将字段数据保存到文本文件中。*/
    FieldVector<Double> r (0.0);
    for (int i = 0; i < N0_; i++ )
      for (int j = 0; j < N1_; j++ )
	for (int k = ( (rank_ == 0) ? 0 : 1 ); k < np_ - ( ( rank_ == size_  - 1 ) ? 0 : 1 ); k++ )
	  {
	    pf_.m = k * N1_ * N0_ + i * N1_ + j;
	    r = rc(pf_.m);
	    *pf_.file << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" ;

	    for (unsigned l = 0; l < seed_.profileField_.size(); l++)
	      {
		if 		( seed_.profileField_[l] == Ex )
		  *pf_.file << en_[pf_.m][0] << "\t";
		else if	( seed_.profileField_[l] == Ey )
		  *pf_.file << en_[pf_.m][1] << "\t";
		else if	( seed_.profileField_[l] == Ez )
		  *pf_.file << en_[pf_.m][2] << "\t";

		else if	( seed_.profileField_[l] == Bx )
		  *pf_.file << bn_[pf_.m][0] << "\t";
		else if	( seed_.profileField_[l] == By )
		  *pf_.file << bn_[pf_.m][1] << "\t";
		else if	( seed_.profileField_[l] == Bz )
		  *pf_.file << bn_[pf_.m][2] << "\t";

		else if 	( seed_.profileField_[l] == Ax )
		  *pf_.file << (*an_)[pf_.m][0] << "\t";
		else if 	( seed_.profileField_[l] == Ay )
		  *pf_.file << (*an_)[pf_.m][1] << "\t";
		else if 	( seed_.profileField_[l] == Az )
		  *pf_.file << (*an_)[pf_.m][2] << "\t";
	      }

	    *pf_.file << std::endl;
	  }

    /*关闭文件。*/
    (*pf_.file).close();
  }
}
