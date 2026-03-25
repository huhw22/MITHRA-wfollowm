#include "boostframe.h"

namespace MITHRA
{
  BoostFrameTransform::BoostFrameTransform()
    : gamma_(1.0), beta_(0.0), c0_(1.0), zRef_(0.0)
  {}

  BoostFrameTransform::BoostFrameTransform(Double gamma, Double beta,
                                                 Double c0, Double zRef)
    : gamma_(gamma), beta_(beta), c0_(c0), zRef_(zRef)
  {}

  void BoostFrameTransform::set(Double gamma, Double beta,
                                   Double c0, Double zRef)
  {
    gamma_ = gamma;
    beta_  = beta;
    c0_    = c0;
    zRef_  = zRef;
  }

  Double BoostFrameTransform::gamma() const { return gamma_; }
  Double BoostFrameTransform::beta()  const { return beta_;  }
  Double BoostFrameTransform::c0()    const { return c0_;    }
  Double BoostFrameTransform::zRef()  const { return zRef_;  }

  void BoostFrameTransform::labToBox(Double tLab, Double zLab,
                                        Double& tBox, Double& zBox) const
  {
    tBox = gamma_ * ( tLab - beta_ * zLab / c0_ );
    zBox = zRef_  + gamma_ * ( zLab - beta_ * c0_ * tLab );
  }

  void BoostFrameTransform::boxToLab(Double tBox, Double zBox,
                                        Double& tLab, Double& zLab) const
  {
    Double zShift = zBox - zRef_;
    zLab = gamma_ * ( zShift + beta_ * c0_ * tBox );
    tLab = gamma_ * ( tBox   + beta_ * zShift / c0_ );
  }

  Double BoostFrameTransform::boxZFromLabZAndBoxT(Double zLab, Double tBox) const
  {
    return zRef_ + zLab / gamma_ - beta_ * c0_ * tBox;
  }

  Double BoostFrameTransform::boxTFromLabZT(Double zLab, Double tLab) const
  {
    return gamma_ * ( tLab - beta_ * zLab / c0_ );
  }

  Double BoostFrameTransform::labZFromBoxZT(Double zBox, Double tBox) const
  {
    return gamma_ * ( (zBox - zRef_) + beta_ * c0_ * tBox );
  }

  Double BoostFrameTransform::labTFromBoxZT(Double zBox, Double tBox) const
  {
    return gamma_ * ( tBox + beta_ * (zBox - zRef_) / c0_ );
  }

  Double BoostFrameTransform::boxZFromLabZT(Double zLab, Double tLab) const
  {
    return zRef_ + gamma_ * ( zLab - beta_ * c0_ * tLab );
  }

  Double BoostFrameTransform::boxTFromBoxZLabT(Double zBox, Double tLab) const
  {
    return tLab / gamma_ - beta_ * ( zBox - zRef_ ) / c0_;
  }

  Double BoostFrameTransform::labZFromBoxZLabT(Double zBox, Double tLab) const
  {
    return ( zBox - zRef_ ) / gamma_ + beta_ * c0_ * tLab;
  }

  Double BoostFrameTransform::labTFromLabZBoxT(Double zLab, Double tBox) const
  {
    return tBox / gamma_ + beta_ * zLab / c0_;
  }
}