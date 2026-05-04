/********************************************************************************************************
* database.cpp: mithra中数据库相关功能的实现
********************************************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "database.h"
#include "stdinclude.h"

namespace MITHRA
{
  /*将堆初始化的参数初始化为一些初始值。*/
  BunchInitialize::BunchInitialize ()
  {
    bunchType_			= "";
    distribution_		= "";
    numberOfParticles_ 		= 0;
    cloudCharge_		= 0.0;
    initialGamma_		= 0.0;
    initialBeta_		= 0.0;
    initialDirection_		= 0.0;
    position_.clear();
    numbers_			= 0;
    latticeConstants_		= 0.0;
    sigmaPosition_		= 0.0;
    sigmaGammaBeta_		= 0.0;
    tranTrun_			= 0.0;
    longTrun_			= 0.0;
    fileName_			= "";
    bF_				= 0.0;
    bFP_			= 0.0;
    shotNoise_			= false;
    lambda_			= 0.0;
  }

  UpdateBunchParallel::UpdateBunchParallel ()
  {
    qSB.clear(); qSF.clear(); qRB.clear(); qRF.clear();
  }

  /*利用非标准有限差分算法推进磁势。*/
  void AdvanceField::advanceMagneticPotentialNSFD (Q v0 , P v1 , P v2 ,
						   P v3 , P v31, P v32,
						   P v4 , P v41, P v42,
						   P v5 , P v51, P v52,
						   P v6 , P v61, P v62,
						   P v7 , P v8 , P v9 )
  {

    *v0 =
	*ufa_     * *v2 - *v1 + alpha_ * (
	*(ufa_+1) * ( *v3 + *v4 + beta_ * ( *v31 + *v32 + *v41 + *v42 ) )   +
	*(ufa_+2) * ( *v5 + *v6 + beta_ * ( *v51 + *v52 + *v61 + *v62 ) ) ) +
	*(ufa_+3) * ( *v7 + *v8 ) +
	*(ufa_+4) * *v9;

    *(v0+1) =
	*ufa_     * *(v2+1) - *(v1+1) + alpha_ * (
	*(ufa_+1) * ( *(v3+1) + *(v4+1) + beta_ * ( *(v31+1) + *(v32+1) + *(v41+1) + *(v42+1) ) )   +
	*(ufa_+2) * ( *(v5+1) + *(v6+1) + beta_ * ( *(v51+1) + *(v52+1) + *(v61+1) + *(v62+1) ) ) ) +
	*(ufa_+3) * ( *(v7+1) + *(v8+1) ) +
	*(ufa_+4) * *(v9+1);

    *(v0+2) =
	*ufa_     * *(v2+2) - *(v1+2) + alpha_ * (
	*(ufa_+1) * ( *(v3+2) + *(v4+2) + beta_ * ( *(v31+2) + *(v32+2) + *(v41+2) + *(v42+2) ) )   +
	*(ufa_+2) * ( *(v5+2) + *(v6+2) + beta_ * ( *(v51+2) + *(v52+2) + *(v61+2) + *(v62+2) ) ) ) +
	*(ufa_+3) * ( *(v7+2) + *(v8+2) ) +
	*(ufa_+4) * *(v9+2);
  };

  /*利用非标准有限差分算法推进标量势。*/
  void AdvanceField::advanceScalarPotentialNSFD (Q v0 , P v1 , P v2 ,
						 P v3 , P v31, P v32,
						 P v4 , P v41, P v42,
						 P v5 , P v51, P v52,
						 P v6 , P v61, P v62,
						 P v7 , P v8 , P v9 )
  {
    *v0 =
	*ufa_     * *v2 - *v1 + alpha_ * (
	*(ufa_+1) * ( *v3 + *v4 + beta_ * ( *v31 + *v32 + *v41 + *v42 ) )   +
	*(ufa_+2) * ( *v5 + *v6 + beta_ * ( *v51 + *v52 + *v61 + *v62 ) ) ) +
	*(ufa_+3) * ( *v7 + *v8 ) +
	*(ufa_+5) * *v9;
  };

  /*利用标准有限差分算法推进磁势。*/
  void AdvanceField::advanceMagneticPotentialFD (Q v0 , P v1 , P v2 ,
						 P v3 , P v31, P v32,
						 P v4 , P v41, P v42,
						 P v5 , P v51, P v52,
						 P v6 , P v61, P v62,
						 P v7 , P v8 , P v9 )
  {
    *v0 =
	*ufa_     * *v2 - *v1 +
	*(ufa_+1) * ( *v3 + *v4 ) +
	*(ufa_+2) * ( *v5 + *v6 ) +
	*(ufa_+3) * ( *v7 + *v8 ) +
	*(ufa_+4) * *v9;

    *(v0+1) =
	*ufa_     * *(v2+1) - *(v1+1) +
	*(ufa_+1) * ( *(v3+1) + *(v4+1) ) +
	*(ufa_+2) * ( *(v5+1) + *(v6+1) ) +
	*(ufa_+3) * ( *(v7+1) + *(v8+1) ) +
	*(ufa_+4) * *(v9+1);

    *(v0+2) =
	*ufa_     * *(v2+2) - *(v1+2) +
	*(ufa_+1) * ( *(v3+2) + *(v4+2) ) +
	*(ufa_+2) * ( *(v5+2) + *(v6+2) ) +
	*(ufa_+3) * ( *(v7+2) + *(v8+2) ) +
	*(ufa_+4) * *(v9+2);
  };

  /*利用标准有限差分算法推进标量势。*/
  void AdvanceField::advanceScalarPotentialFD	(Q v0 , P v1 , P v2 ,
						 P v3 , P v31, P v32,
						 P v4 , P v41, P v42,
						 P v5 , P v51, P v52,
						 P v6 , P v61, P v62,
						 P v7 , P v8 , P v9 )
  {
    *v0 =
	*ufa_     * *v2 - *v1 +
	*(ufa_+1) * ( *v3 + *v4 ) +
	*(ufa_+2) * ( *v5 + *v6 ) +
	*(ufa_+3) * ( *v7 + *v8 ) +
	*(ufa_+5) * *v9;
  }

	void AdvanceField::advanceMagneticPotentialFD_CPMLXY(
		Q v0 , P v1 , P v2 ,
		P v3 , P v31, P v32,
		P v4 , P v41, P v42,
		P v5 , P v51, P v52,
		P v6 , P v61, P v62,
		P v7 , P v8 , P v9,

		Q psiXnp1, P psiXn_c, P psiXn_p, P psiXn_m,
		Q psiYnp1, P psiYn_c, P psiYn_p, P psiYn_m,

		Double bx, Double cx, Double kappaX,
		Double by, Double cy, Double kappaY,
		Double dx, Double dy)
	{
		/* 1) 先按原始 FD 矢势核推进 */
		advanceMagneticPotentialFD(
			v0, v1, v2,
			v3, v31, v32,
			v4, v41, v42,
			v5, v51, v52,
			v6, v61, v62,
			v7, v8, v9);

	}

	void AdvanceField::advanceMagneticPotentialNSFD_CPML(
		Q v0 , P v1 , P v2 ,
		P v3 , P v31, P v32,
		P v4 , P v41, P v42,
		P v5 , P v51, P v52,
		P v6 , P v61, P v62,
		P v7 , P v8 , P v9,
		Double eta, Double dt)
	{
		(void)v9; // sponge region: source is explicitly ignored

		const Double coeffT = *ufa_;
		const Double coeffX = alpha_ * (*(ufa_ + 1));
		const Double coeffY = alpha_ * (*(ufa_ + 2));
		const Double coeffZ = *(ufa_ + 3);

		const Double damp  = 0.5 * ((eta > 0.0) ? eta : 0.0) * dt;
		const Double denom = 1.0 + damp;

		for (int c = 0; c < 3; ++c)
		{
			const Double uc   = *(v2 + c);
			const Double upre = *(v1 + c);

			const Double ux_p = *(v3 + c);
			const Double ux_m = *(v4 + c);

			const Double uy_p = *(v5 + c);
			const Double uy_m = *(v6 + c);

			const Double uz_p = *(v7 + c);
			const Double uz_m = *(v8 + c);

			const Double ufree =
				coeffT * uc - upre
				+ coeffX * (
					ux_p + ux_m
					+ beta_ * (
						*(v31 + c) + *(v32 + c)
					+ *(v41 + c) + *(v42 + c)
					)
				)
				+ coeffY * (
					uy_p + uy_m
					+ beta_ * (
						*(v51 + c) + *(v52 + c)
					+ *(v61 + c) + *(v62 + c)
					)
				)
				+ coeffZ * (uz_p + uz_m);

			*(v0 + c) = (ufree + damp * upre) / denom;
		}
	}

	void AdvanceField::advanceScalarPotentialFD_CPMLXY(
		Q v0 , P v1 , P v2 ,
		P v3 , P v31, P v32,
		P v4 , P v41, P v42,
		P v5 , P v51, P v52,
		P v6 , P v61, P v62,
		P v7 , P v8 , P v9,

		Q psiXnp1, P psiXn_c, P psiXn_p, P psiXn_m,
		Q psiYnp1, P psiYn_c, P psiYn_p, P psiYn_m,

		Double bx, Double cx, Double kappaX,
		Double by, Double cy, Double kappaY,
		Double dx, Double dy)
	{
		/* 1) 先按原始 FD 标量势核推进 */
		advanceScalarPotentialFD(
			v0, v1, v2,
			v3, v31, v32,
			v4, v41, v42,
			v5, v51, v52,
			v6, v61, v62,
			v7, v8, v9);

		/* 2) 用当前时间层场，更新 x/y 方向的 CPML 记忆变量 */
		const Double dphidx = (*v3 - *v4) / (2.0 * dx);
		const Double dphidy = (*v5 - *v6) / (2.0 * dy);

		*psiXnp1 = bx * (*psiXn_c) + cx * dphidx;
		*psiYnp1 = by * (*psiYn_c) + cy * dphidy;

		/* 3) 用旧层 psi 的邻点，计算 d(psi)/dx, d(psi)/dy */
		const Double dpsiXdx = (*psiXn_p - *psiXn_m) / (2.0 * dx);
		const Double dpsiYdy = (*psiYn_p - *psiYn_m) / (2.0 * dy);

		/* 4) 对 x/y 二阶导数做 kappa 修正
			原核已经包含了普通 x/y Laplacian，这里只加“相对旧核的修正量” */
		const Double lapX = (*v3 - 2.0 * (*v2) + *v4);
		const Double lapY = (*v5 - 2.0 * (*v2) + *v6);

		*v0 += (*(ufa_ + 1)) * (1.0 / (kappaX * kappaX) - 1.0) * lapX;
		*v0 += (*(ufa_ + 2)) * (1.0 / (kappaY * kappaY) - 1.0) * lapY;

		/* 5) 加上 psi 导数修正项
			这里先用最直接的一阶显式写法 */
		*v0 += (*(ufa_ + 1)) * dx * dpsiXdx;
		*v0 += (*(ufa_ + 2)) * dy * dpsiYdy;
	}

	void AdvanceField::advanceScalarPotentialNSFD_CPMLXY(
		Q v0 , P v1 , P v2 ,
		P v3 , P v31, P v32,
		P v4 , P v41, P v42,
		P v5 , P v51, P v52,
		P v6 , P v61, P v62,
		P v7 , P v8 , P v9,

		Q psiXnp1, P psiXn_c, P psiXn_p, P psiXn_m,
		Q psiYnp1, P psiYn_c, P psiYn_p, P psiYn_m,

		Double bx, Double cx, Double kappaX,
		Double by, Double cy, Double kappaY,
		Double dx, Double dy)
	{
		/* 1) 先按原始 NSFD 标量势核推进 */
		advanceScalarPotentialNSFD(
			v0, v1, v2,
			v3, v31, v32,
			v4, v41, v42,
			v5, v51, v52,
			v6, v61, v62,
			v7, v8, v9);

		/* 2) 用当前时间层标量势，更新 x/y 方向的 CPML 记忆变量 */
		const Double dphidx = (*v3 - *v4) / (2.0 * dx);
		const Double dphidy = (*v5 - *v6) / (2.0 * dy);

		*psiXnp1 = bx * (*psiXn_c) + cx * dphidx;
		*psiYnp1 = by * (*psiYn_c) + cy * dphidy;

		/* 3) 用旧层 psi 的邻点，计算 dpsi/dx, dpsi/dy */
		const Double dpsiXdx = (*psiXn_p - *psiXn_m) / (2.0 * dx);
		const Double dpsiYdy = (*psiYn_p - *psiYn_m) / (2.0 * dy);

		/* 4) 对 x/y 二阶导数做 kappa 修正
			NSFD 原核已经包含普通 x/y 项，这里只加“相对原核的修正量” */
		const Double lapX =
			(*v3 + *v4 + beta_ * (*v31 + *v32 + *v41 + *v42) - 2.0 * (1.0 + 2.0 * beta_) * (*v2));

		const Double lapY =
			(*v5 + *v6 + beta_ * (*v51 + *v52 + *v61 + *v62) - 2.0 * (1.0 + 2.0 * beta_) * (*v2));

		*v0 += alpha_ * (*(ufa_ + 1)) * (1.0 / (kappaX * kappaX) - 1.0) * lapX;
		*v0 += alpha_ * (*(ufa_ + 2)) * (1.0 / (kappaY * kappaY) - 1.0) * lapY;

		/* 5) 加上 psi 导数修正项 */
		*v0 += alpha_ * (*(ufa_ + 1)) * dx * dpsiXdx;
		*v0 += alpha_ * (*(ufa_ + 2)) * dy * dpsiYdy;
	}

  /*在计算域的边界处推进磁势。*/
  void AdvanceField::advanceBoundaryF 		(Q v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13)
  {
    *v0 =
	*ufB_     * ( *v1 + *v5 ) +
	*(ufB_+1) * *v3 +
	*(ufB_+2) * ( *v2 + *v4 ) +
	*(ufB_+3) * ( *v6 + *v7 + *v10 + *v11 ) +
	*(ufB_+4) * ( *v8 + *v9 + *v12 + *v13 );
    *(v0+1) =
	*ufB_     * ( *(v1+1) + *(v5+1) ) +
	*(ufB_+1) * *(v3+1) +
	*(ufB_+2) * ( *(v2+1) + *(v4+1) ) +
	*(ufB_+3) * ( *(v6+1) + *(v7+1) + *(v10+1) + *(v11+1) ) +
	*(ufB_+4) * ( *(v8+1) + *(v9+1) + *(v12+1) + *(v13+1) );
    *(v0+2) =
	*ufB_     * ( *(v1+2) + *(v5+2) ) +
	*(ufB_+1) * *(v3+2) +
	*(ufB_+2) * ( *(v2+2) + *(v4+2) ) +
	*(ufB_+3) * ( *(v6+2) + *(v7+2) + *(v10+2) + *(v11+2) ) +
	*(ufB_+4) * ( *(v8+2) + *(v9+2) + *(v12+2) + *(v13+2) );
  }

  /*在计算域的边界处推进标量势。*/
  void AdvanceField::advanceBoundaryS 		(Q v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13)
  {
    *v0 =
	*ufB_     * ( *v1 + *v5 ) +
	*(ufB_+1) * *v3 +
	*(ufB_+2) * ( *v2 + *v4 ) +
	*(ufB_+3) * ( *v6 + *v7 + *v10 + *v11 ) +
	*(ufB_+4) * ( *v8 + *v9 + *v12 + *v13 );
  }

  /*在计算域的边缘推进磁势。*/
  void AdvanceField::advanceEdgeF 		(Q v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19)
  {
    *v0     =
	*ufB_ *     ( *v3 + *v8 ) +
	*(ufB_+1) * ( *v5 + *v6 ) +
	*(ufB_+2) * ( *v2 + *v9 ) +
	*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
	*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
	*v11;

    *(v0+1) =
	*ufB_ *     ( *(v3+1) + *(v8+1) ) +
	*(ufB_+1) * ( *(v5+1) + *(v6+1) ) +
	*(ufB_+2) * ( *(v2+1) + *(v9+1) ) +
	*(ufB_+3) * ( *(v1+1) + *(v4+1) + *(v7+1) + *(v10+1) ) +
	*(ufB_+4) * ( *(v12+1)+ *(v13+1)+ *(v14+1)+ *(v15+1)+ *(v16+1)+ *(v17+1)+ *(v18+1)+ *(v19+1)) -
	*(v11+1);

    *(v0+2) =
	*ufB_ *     ( *(v3+2) + *(v8+2) ) +
	*(ufB_+1) * ( *(v5+2) + *(v6+2) ) +
	*(ufB_+2) * ( *(v2+2) + *(v9+2) ) +
	*(ufB_+3) * ( *(v1+2) + *(v4+2) + *(v7+2) + *(v10+2) ) +
	*(ufB_+4) * ( *(v12+2)+ *(v13+2)+ *(v14+2)+ *(v15+2)+ *(v16+2)+ *(v17+2)+ *(v18+2)+ *(v19+2)) -
	*(v11+2);

  }

  /*在计算域的边缘推进标量势。*/
  void AdvanceField::advanceEdgeS 		(Q v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19)
  {
    *v0     =
	*ufB_ *     ( *v3 + *v8 ) +
	*(ufB_+1) * ( *v5 + *v6 ) +
	*(ufB_+2) * ( *v2 + *v9 ) +
	*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
	*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
	*v11;
  }

  /*推进计算域边角处的磁势。*/
  void AdvanceField::advanceCornerF 		(Q v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19, P v20,
						 P v21, P v22, P v23)
  {
    *v0     = - ( *v1  * *(ufB_+16) + *v2  * *(ufB_+8) +
	*v3  * *(ufB_+1) + *v4  * *(ufB_+16) + *v5  * *(ufB_+9) +
	*v6  * *(ufB_+2) + *v7  * *(ufB_+16) + *v8  * *(ufB_+10) +
	*v9  * *(ufB_+3) + *v10 * *(ufB_+16) + *v11 * *(ufB_+11) +
	*v12 * *(ufB_+4) + *v13 * *(ufB_+16) + *v14 * *(ufB_+12) +
	*v15 * *(ufB_+5) + *v16 * *(ufB_+16) + *v17 * *(ufB_+13) +
	*v18 * *(ufB_+6) + *v19 * *(ufB_+16) + *v20 * *(ufB_+14) +
	*v21 * *(ufB_+7) + *v22 * *(ufB_+16) + *v23 * *(ufB_+15)) / *ufB_;

    *(v0+1) = - ( *(v1+1)  * *(ufB_+16) + *(v2+1)  * *(ufB_+8) +
	*(v3+1)  * *(ufB_+1) + *(v4+1)  * *(ufB_+16) + *(v5+1)  * *(ufB_+9) +
	*(v6+1)  * *(ufB_+2) + *(v7+1)  * *(ufB_+16) + *(v8+1)  * *(ufB_+10) +
	*(v9+1)  * *(ufB_+3) + *(v10+1) * *(ufB_+16) + *(v11+1) * *(ufB_+11) +
	*(v12+1) * *(ufB_+4) + *(v13+1) * *(ufB_+16) + *(v14+1) * *(ufB_+12) +
	*(v15+1) * *(ufB_+5) + *(v16+1) * *(ufB_+16) + *(v17+1) * *(ufB_+13) +
	*(v18+1) * *(ufB_+6) + *(v19+1) * *(ufB_+16) + *(v20+1) * *(ufB_+14) +
	*(v21+1) * *(ufB_+7) + *(v22+1) * *(ufB_+16) + *(v23+1) * *(ufB_+15)) / *ufB_;

    *(v0+2) = - ( *(v1+2)  * *(ufB_+16) + *(v2+2)  * *(ufB_+8) +
	*(v3+2)  * *(ufB_+1) + *(v4+2)  * *(ufB_+16) + *(v5+2)  * *(ufB_+9) +
	*(v6+2)  * *(ufB_+2) + *(v7+2)  * *(ufB_+16) + *(v8+2)  * *(ufB_+10) +
	*(v9+2)  * *(ufB_+3) + *(v10+2) * *(ufB_+16) + *(v11+2) * *(ufB_+11) +
	*(v12+2) * *(ufB_+4) + *(v13+2) * *(ufB_+16) + *(v14+2) * *(ufB_+12) +
	*(v15+2) * *(ufB_+5) + *(v16+2) * *(ufB_+16) + *(v17+2) * *(ufB_+13) +
	*(v18+2) * *(ufB_+6) + *(v19+2) * *(ufB_+16) + *(v20+2) * *(ufB_+14) +
	*(v21+2) * *(ufB_+7) + *(v22+2) * *(ufB_+16) + *(v23+2) * *(ufB_+15)) / *ufB_;

  }

  /*在计算域的边角处推进标量势。*/
  void AdvanceField::advanceCornerS 		(Q v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19, P v20,
						 P v21, P v22, P v23)
  {
    *v0     = - ( *v1  * *(ufB_+16) + *v2  * *(ufB_+8) +
	*v3  * *(ufB_+1) + *v4  * *(ufB_+16) + *v5  * *(ufB_+9) +
	*v6  * *(ufB_+2) + *v7  * *(ufB_+16) + *v8  * *(ufB_+10) +
	*v9  * *(ufB_+3) + *v10 * *(ufB_+16) + *v11 * *(ufB_+11) +
	*v12 * *(ufB_+4) + *v13 * *(ufB_+16) + *v14 * *(ufB_+12) +
	*v15 * *(ufB_+5) + *v16 * *(ufB_+16) + *v17 * *(ufB_+13) +
	*v18 * *(ufB_+6) + *v19 * *(ufB_+16) + *v20 * *(ufB_+14) +
	*v21 * *(ufB_+7) + *v22 * *(ufB_+16) + *v23 * *(ufB_+15)) / *ufB_;
  }
}
