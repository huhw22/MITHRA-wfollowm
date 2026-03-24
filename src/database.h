/********************************************************************************************************
* database.h：与数据库相关的头文件的实现
********************************************************************************************************/

#ifndef DATABASE_H_
#define DATABASE_H_

#include <fstream>
#include <string>
#include <vector>

#include "fieldvector.h"
#include "stdinclude.h"

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
    float				*en,   *bn;

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
    int					nt;
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
    int					nt;
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

  
  /*从屏幕上获取束形文件所需的数据结构。*/
  struct SampleScreenProfile
  {
    std::vector<std::ofstream*> 	files;
    std::vector<std::string>		fileNames;
  };
}

#endif
