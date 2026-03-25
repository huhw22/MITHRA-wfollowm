#ifndef BOOSTFRAME_H
#define BOOSTFRAME_H

#include "fieldvector.h"

namespace MITHRA
{
  class BoostFrameTransform
  {
    private:
      Double gamma_;
      Double beta_;
      Double c0_;
      Double zRef_;

    public:
      BoostFrameTransform();
      BoostFrameTransform(Double gamma, Double beta, Double c0, Double zRef);

      void set(Double gamma, Double beta, Double c0, Double zRef);

      Double gamma() const;
      Double beta() const;
      Double c0() const;
      Double zRef() const;

      /* 完整变换：lab -> box */
      void labToBox(Double tLab, Double zLab, Double& tBox, Double& zBox) const;

      /* 完整变换：box -> lab */
      void boxToLab(Double tBox, Double zBox, Double& tLab, Double& zLab) const;

      /* 常用便捷函数 */
      Double boxZFromLabZAndBoxT(Double zLab, Double tBox) const;
      Double boxTFromLabZT(Double zLab, Double tLab) const;
      Double labZFromBoxZT(Double zBox, Double tBox) const;
      Double labTFromBoxZT(Double zBox, Double tBox) const;

      Double boxZFromLabZT(Double zLab, Double tLab) const;
      Double boxTFromBoxZLabT(Double zBox, Double tLab) const;
      Double labZFromBoxZLabT(Double zBox, Double tLab) const;
      Double labTFromLabZBoxT(Double zLab, Double tBox) const;
  };
}

#endif