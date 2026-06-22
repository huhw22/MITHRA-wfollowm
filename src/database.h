/********************************************************************************************************
* database.h：与数据库相关的头文件的实现
********************************************************************************************************/

#ifndef DATABASE_H_
#define DATABASE_H_

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstddef>

#include "fieldvector.h"
#include "stdinclude.h"
#include "detector_field_writer.h"

namespace MITHRA
{

  /*用于初始化一组的数据结构。*/
  struct BunchInitialize
  {
    /*束的类型是手动，椭球，圆柱体，立方体和3d晶体之一。如果是的话
    *手动电荷在点的位置矢量将产生。*/
    std::string     			bunchType_;

    /*束中的分布类型（横向或纵向）。*/
    std::string     			distribution_;

    /*用于创建束分布的生成器类型。*/
    std::string				generator_;

    /*群中宏观粒子的总数。*/
    unsigned int       			numberOfParticles_;

    /*在pC中的总电荷。*/
    Double  				cloudCharge_;

    /*束的初始能量，单位是MeV。*/
    Double            			initialGamma_;

    /*串的初始归一化速度。*/
    Double            			initialBeta_;

    /*束的初始运动方向，这是一个单位矢量。*/
    FieldVector<Double>			initialDirection_;

    /*束的中心在单位长度尺度上的位置。*/
    std::vector<FieldVector<Double> >	position_;

    /*三维晶型在每个方向上的宏观粒子数。*/
    FieldVector<unsigned int>		numbers_;

    /*三维晶体类型在x， y和z方向上的晶格常数。*/
    FieldVector<Double>			latticeConstants_;

    /*在每个方向的位置上以长度为单位展开。对于3D晶体
    类型，它将是每个微束晶体在位置上的分布。*/
    FieldVector<Double>			sigmaPosition_;

    /*向各个方向扩散能量。*/
    FieldVector<Double>			sigmaGammaBeta_;

    /*存储电子产生的截断横向距离。*/
    Double				tranTrun_;

    /*存储电子产生的截断纵向距离。*/
    Double				longTrun_;

    /*读取电子分布的文件名称。*/
    std::string				fileName_;

    /*与波动器外束长相对应的辐射波长*/
    Double				lambda_;

    /*初始化束的聚束因子。*/
    Double				bF_;

    /*初始化束的聚束因子的相位。*/
    Double				bFP_;

    /*确定射击噪声激活的布尔标志。*/
    bool				shotNoise_;

    /*束的初始向量，它是和方向的乘积。*/
    FieldVector<Double>			betaVector_;

    /*将堆初始化的参数初始化为一些初始值。*/
    BunchInitialize ();
  };

  /*更新字段所需的数据结构。*/
  class AdvanceField
  {
  public:

    typedef const Double* 	P;
    typedef Double* 		Q;

    void advanceMagneticPotentialNSFD 	(Q v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceScalarPotentialNSFD 	(Q v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceMagneticPotentialFD 	(Q v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceScalarPotentialFD 	(Q v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );


    void advanceBoundaryF 		(Q v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13);

    void advanceBoundaryS 		(Q v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13);

    void advanceEdgeF 			(Q v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19);

    void advanceEdgeS 			(Q v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19);

    void advanceCornerF 		(Q v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19, P v20,
					 P v21, P v22, P v23);

    void advanceCornerS 		(Q v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19, P v20,
					 P v21, P v22, P v23);
    P 					ufa_;
    P 					ufB_;
    Double 				alpha_;
    Double 				beta_;
  };

  /*更新字段所需的数据结构。*/
  struct UpdateField
  {
    Double				jw, jwt, jb;
    Double				w,  b;
    Double				dt,  dx,  dy,  dz;
    Double				dx2, dy2, dz2;
    Double				a 	[6];

    Double				bB  	[5];
    Double				cB  	[5];
    Double				dB  	[5];

    Double				eE	[5];
    Double				fE	[5];
    Double				gE	[5];

    Double				hC	[17];

    AdvanceField			af;

    Double 				*v0, *v1, *v2;
    Double 				*v3, *v31,*v32;
    Double 				*v4, *v41,*v42;
    Double 				*v5, *v51,*v52;
    Double 				*v6, *v61,*v62;
    Double 				*v7, *v8, *v9;
    Double 				*v10,*v11,*v12,*v13;
    Double 				*v14,*v15,*v16,*v17,*v18,*v19;

    Double				*anp1, *an, *anm1;
    Double				*fnp1, *fn, *fnm1;
    Double				*jn,   *rn;
    Double				*en,   *bn;

    unsigned int			N0m1, N1m1, npm1;
  };

  /*采样字段所需的数据结构。*/
  struct SampleField
  {
    std::ofstream*			file;
    FieldVector<Double>			position;
    int					i, j, k, m;
    FieldVector<Double>			et, at, bt;
    Double				f;
    Double				dxr, dyr, dzr;
    Double				c2, jw, p;
    double				c1;
    unsigned int			N;
    Double				Ce, Cb, Ca, Cj, Cf;
  };

  /*可视化字段所需的数据结构。*/
  struct VisualizeField
  {
    std::ofstream* 			file;
    std::string				fileName, path, name;
    unsigned int 			i, j, k, m;
    std::vector<std::vector<Double> >	v;
    Double				jw;
  };

  /*保存字段配置文件所需的数据结构。*/
  struct ProfileField
  {
    unsigned int 			i, j, k, l, m;
    std::ofstream* 			file;
    std::string				fileName;
    FieldVector<Double>			e, a, b;
    Double				dt;
    Double				jw;
    int                                 kmin, kmax;
  };

  /*更新群并行操作所需的参数。*/
  struct UpdateBunchParallel
  {
    int					i, j, k, n;
    long int                            m;
    bool                                b1x, b1y, b1z;
    int nOutX;         // x 方向出横向域
    int nOutY;         // y 方向出横向域
    int nOutXY;        // x 或 y 任一方向出横向域
    int nCrossZSlab;   // 跨出当前 MPI z 子域
    int nOutBoxAny;    // 任一出 box 事件总数
    std::vector<Charge>                 qSB, qSF, qRB, qRF;
    bool				dq;

    FieldVector<Double>			et, bt;
    Double				dxr, dyr, dzr;
    Double				d2;
    double				d1;
    FieldVector<Double>			gbm, gbp, gbpl, dr;
    Double				sz, cy;

    Double                              lz, ly;
    Double				tsignal, tsignalm, tsignale, tsignalb;
    Double				d, l, wrp, wrs, p0, p1, t;
    Double				x, y, z, x0, y0, r0;
    Double				af;
    Double                              atanS, atanP;
    FieldVector<Double>			rv, yv, rl;
    FieldVector<Double>			ex, ez, eT;
    FieldVector<Double>			by, bz, bT;
    Double				tl, t0, tlm;

    Double				zr;

    UpdateBunchParallel();
  };

  /*保存更新群所需的数据结构。*/
  struct UpdateBunch
  {
    Double				dx, dy, dz, dt, dtb;
    Double				ku, b0;
    Double				r1, r2;
    Double				ct, st;
    int                                 nL, i;
    Charge                              Q;
    int nOutX;         // x 方向出横向域
    int nOutY;         // y 方向出横向域
    int nOutXY;        // x 或 y 任一方向出横向域
    int nCrossZSlab;   // 跨出当前 MPI z 子域
    int nOutBoxAny;    // 任一出 box 事件总数
  };

  /*采样群所需的数据结构。*/
  struct SampleBunch
  {
    std::ofstream*			file;
    Double				g, q, qT;
    FieldVector<Double>			r,  gb,  r2,  gb2;
    FieldVector<Double>			rT, gbT, r2T, gb2T;
  };

  /*可视化群所需的数据结构。*/
  struct VisualizeBunch
  {
    std::ofstream* 			file;
    std::string				fileName, path, name;
    unsigned int			i, N;
  };

  /*剖析群所需的数据结构。*/
  struct ProfileBunch
  {
    std::ofstream* 			file;
    std::string				fileName;
  };

  /*结构的数据需要更新当前。*/
  struct UpdateCurrent
  {
    int					ip, jp, kp, im, jm, km;
    long int				m;
    Double				dx, dy, dz;
    Double				dxp, dyp, dzp, dxm, dym, dzm;
    Double				x1, x2, y1, y2, z1, z2, q;
    FieldVector<Double>			rp, rm, r;
    FieldVector<Double>			jcp, jcm;
    Double				rc;
    Double				dv;
    Double				c;
    std::vector<FieldVector<Double> >	jt;
    std::vector<Double>			rt;
  };

    struct SamplePlaneFieldH5
  {
    DetectorFieldWriter writer;

    std::vector<double> xCoord;
    std::vector<double> yCoord;

    std::vector<double> exFrame;
    std::vector<double> eyFrame;
    std::vector<double> ezFrame;
    std::vector<double> bxFrame;
    std::vector<double> byFrame;
    std::vector<double> bzFrame;

    int nx = 0;
    int ny = 0;

    std::string fileName;
    bool useFloat32 = true;

    void resize(int nxIn, int nyIn)
    {
      nx = nxIn;
      ny = nyIn;

      const std::size_t n = static_cast<std::size_t>(nx) *
                            static_cast<std::size_t>(ny);

      xCoord.resize(nx);
      yCoord.resize(ny);

      exFrame.assign(n, 0.0);
      eyFrame.assign(n, 0.0);
      ezFrame.assign(n, 0.0);
      bxFrame.assign(n, 0.0);
      byFrame.assign(n, 0.0);
      bzFrame.assign(n, 0.0);
    }

    void clearFrame()
    {
      std::fill(exFrame.begin(), exFrame.end(), 0.0);
      std::fill(eyFrame.begin(), eyFrame.end(), 0.0);
      std::fill(ezFrame.begin(), ezFrame.end(), 0.0);
      std::fill(bxFrame.begin(), bxFrame.end(), 0.0);
      std::fill(byFrame.begin(), byFrame.end(), 0.0);
      std::fill(bzFrame.begin(), bzFrame.end(), 0.0);
    }

    void close()
    {
      writer.close();
    }
  };


  /*采样字段所需的数据结构。*/
  struct SampleRadiationPower
  {
    std::vector<std::ofstream*>				file;
    std::vector<Double>					w;
    int							k, l;
    Double						dxr, dyr, dzr;
    double						c;
    std::vector<std::vector<std::vector<Double> > >	pLi;
    std::vector<Double>					pL, pG;
    Double						dt, dx, dy, dz;
    Double						pc;
    unsigned int					N, Nl, Nf, Nz;
    std::vector<std::vector<std::vector<Double> > >	fdt;
    unsigned int					m;
    std::vector<std::vector<Complex> >			ep, em;
    std::vector<SamplePlaneFieldH5> fieldPlanes;
  };

  /*采样字段所需的数据结构。*/
  struct SampleRadiationEnergy
  {
    std::vector<std::ofstream*>				file;
    std::vector<Double>					w;
    int							k, l;
    Double						dxr, dyr, dzr;
    Double						c;
    std::vector<std::vector<std::vector<Double> > >     pLi;
    std::vector<Double>					pL, pG;
    Double						dt, dx, dy, dz;
    Double						pc;
    unsigned int					N, Nl, Nf;
    std::vector<std::vector<std::vector<Double> > >     fdt;
  };


  struct SampleRadiationDetector
  {
    /* 当前 detector 面是否处于有效采样区间 */
    bool                                      active;

    /* 是否已经第一次进入有效采样区间，用于定义 time=0 */
    bool                                      started;

    /* 是否已经完成探测 */
    bool finished;

    /* 输出文件：第一版只写 power-line */
    std::ofstream*                            file;

    /* detector 面在 lab 系的目标位置 */
    Double                                    zLab;

    /* 当前时刻 detector 面在 box 系中的位置 */
    Double                                    zBox;

    /* 将 zBox 反算回 lab 系的位置，用于诊断检查 */
    Double                                    zLabCheck;

    /* 当前输出的 lab 系时间，以及第一次有效采样时刻 */
    Double                                    tLab0;
    Double                                    tLab;

    /* 当前拥有 detector 面的 rank */
    int                                       ownerRank;

    /* 当前 detector 面对应的全局 z 索引 */
    long int                                  k;

    /* 采样间隔与插值相关参数 */
    Double                                    dzr;
    Double                                    c;

    /* 当前步功率：本地与全局 */
    Double                                    pL;
    Double                                    pG;

    /*全场输出*/
    DetectorFieldWriter fieldWriter;
    std::vector<double> xCoord;
    std::vector<double> yCoord;
    std::vector<double>  exFrame, eyFrame, ezFrame, bxFrame, byFrame, bzFrame;

    int                 fieldNx;
    int                 fieldNy;
    std::string         fieldFileName;
    bool                fieldUseFloat32;
  };

  
  /*从屏幕上获取束形文件所需的数据结构。*/
  struct SampleScreenProfile
  {
    std::vector<std::ofstream*> 	files;
    std::vector<std::string>		fileNames;
  };


  struct ScalarCPML
  {
    /*
    * 开关与厚度
    */
    bool enabled_ = false;

    int px_ = 0;
    int py_ = 0;
    int pz_ = 0;

    /*
    * 网格尺寸
    *
    * n0_ : x direction grid size
    * n1_ : y direction grid size
    * np_ : local z direction grid size on this rank
    */
    int n0_ = 0;
    int n1_ = 0;
    int np_ = 0;

    /*
    * slot arrays
    *
    * xSlot_[i] >= 0 means x-index i is inside x-PML.
    * ySlot_[j] >= 0 means y-index j is inside y-PML.
    * zSlot_[k] >= 0 means local z-index k is inside global z-PML.
    *
    * slot value is the compact index used for memory arrays.
    */
    std::vector<int> xSlot_;
    std::vector<int> ySlot_;
    std::vector<int> zSlot_;

    int nxSlot_ = 0;
    int nySlot_ = 0;
    int nzSlot_ = 0;

    /*
    * CPML coefficients
    *
    * D f = df / kappa + psi
    * psi_new = b * psi_old + a * df
    */
    std::vector<Double> kx_;
    std::vector<Double> ky_;
    std::vector<Double> kz_;

    std::vector<Double> ax_;
    std::vector<Double> ay_;
    std::vector<Double> az_;

    std::vector<Double> bx_;
    std::vector<Double> by_;
    std::vector<Double> bz_;

    /*
    * Directional NSFD-CPML scaling
    *
    * Cx = Cy = Cz = (c dt)^2.  The grid spacing is applied by the
    * first-difference operators, so for example
    *
    *   (sx / dx)^2 D_x^- D_x^+ W_z A
    *
    * is the x contribution to the dimensionless time update.  Do not
    * derive Cz from the expanded NSFD coefficient a3; a3 already includes
    * the transverse W_z correction carried by z-neighbour points.
    */
    Double Cx_ = 0.0;
    Double Cy_ = 0.0;
    Double Cz_ = 0.0;

    Double sx_ = 0.0;
    Double sy_ = 0.0;
    Double sz_ = 0.0;

    Double invDx_ = 0.0;
    Double invDy_ = 0.0;
    Double invDz_ = 0.0;

    Double sxInvDx_ = 0.0;
    Double syInvDy_ = 0.0;
    Double szInvDz_ = 0.0;

    /*
    * First-layer derivative buffers:
    *
    * gx_ = sx * D_x^+_pml W_z A
    * gy_ = sy * D_y^+_pml W_z A
    * gz_ = sz * D_z^+_pml A
    *
    * Each one is a 3-component vector field.
    */
    std::vector<FieldVector<Double>> gx_;
    std::vector<FieldVector<Double>> gy_;
    std::vector<FieldVector<Double>> gz_;

    /*
    * CPML memories for scalar second-order operator.
    *
    * xp_A : x plus  derivative acting on W_z A
    * xm_G : x minus derivative acting on gx
    *
    * yp_A : y plus  derivative acting on W_z A
    * ym_G : y minus derivative acting on gy
    *
    * zp_A : z plus  derivative acting on A
    * zm_G : z minus derivative acting on gz
    *
    * Each memory entry is FieldVector<Double>, so Ax/Ay/Az are stored
    * independently but share the same spatial slot.
    */
    std::vector<FieldVector<Double>> psi_xp_A_;
    std::vector<FieldVector<Double>> psi_xm_G_;

    std::vector<FieldVector<Double>> psi_yp_A_;
    std::vector<FieldVector<Double>> psi_ym_G_;

    std::vector<FieldVector<Double>> psi_zp_A_;
    std::vector<FieldVector<Double>> psi_zm_G_;

    /*
    * Index helpers.
    *
    * x-memory layout:
    *   [k][j][sx]
    *
    * y-memory layout:
    *   [k][i][sy]
    *
    * z-memory layout:
    *   [sz][i][j]
    */
    inline long idxX(int sx, int j, int k) const
    {
      return ((long)k * n1_ + j) * nxSlot_ + sx;
    }

    inline long idxY(int i, int sy, int k) const
    {
      return ((long)k * n0_ + i) * nySlot_ + sy;
    }

    inline long idxZ(int i, int j, int sz) const
    {
      return ((long)sz * n0_ + i) * n1_ + j;
    }

    inline long gridSize() const
    {
      return (long)n0_ * n1_ * np_;
    }

    void clear()
    {
      enabled_ = false;

      px_ = py_ = pz_ = 0;
      n0_ = n1_ = np_ = 0;

      nxSlot_ = nySlot_ = nzSlot_ = 0;

      xSlot_.clear();
      ySlot_.clear();
      zSlot_.clear();

      kx_.clear();
      ky_.clear();
      kz_.clear();

      ax_.clear();
      ay_.clear();
      az_.clear();

      bx_.clear();
      by_.clear();
      bz_.clear();

      gx_.clear();
      gy_.clear();
      gz_.clear();

      psi_xp_A_.clear();
      psi_xm_G_.clear();

      psi_yp_A_.clear();
      psi_ym_G_.clear();

      psi_zp_A_.clear();
      psi_zm_G_.clear();

      Cx_ = Cy_ = Cz_ = 0.0;
      sx_ = sy_ = sz_ = 0.0;

      invDx_ = invDy_ = invDz_ = 0.0;
      sxInvDx_ = syInvDy_ = szInvDz_ = 0.0;
    }

    void resetMemories()
    {
      std::fill(gx_.begin(), gx_.end(), FieldVector<Double>(0.0));
      std::fill(gy_.begin(), gy_.end(), FieldVector<Double>(0.0));
      std::fill(gz_.begin(), gz_.end(), FieldVector<Double>(0.0));

      std::fill(psi_xp_A_.begin(), psi_xp_A_.end(), FieldVector<Double>(0.0));
      std::fill(psi_xm_G_.begin(), psi_xm_G_.end(), FieldVector<Double>(0.0));

      std::fill(psi_yp_A_.begin(), psi_yp_A_.end(), FieldVector<Double>(0.0));
      std::fill(psi_ym_G_.begin(), psi_ym_G_.end(), FieldVector<Double>(0.0));

      std::fill(psi_zp_A_.begin(), psi_zp_A_.end(), FieldVector<Double>(0.0));
      std::fill(psi_zm_G_.begin(), psi_zm_G_.end(), FieldVector<Double>(0.0));
    }
  };
  
}

#endif
