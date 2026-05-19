#ifndef DETECTOR_FIELD_WRITER_H
#define DETECTOR_FIELD_WRITER_H

#include <string>
#include <vector>
#include <cstddef>

#include "hdf5.h"

namespace MITHRA
{
  class DetectorFieldWriter
  {
  public:
    DetectorFieldWriter();
    ~DetectorFieldWriter();

    // 禁止拷贝，避免 HDF5 句柄重复释放
    DetectorFieldWriter(const DetectorFieldWriter&) = delete;
    DetectorFieldWriter& operator=(const DetectorFieldWriter&) = delete;

    // 允许移动（可选；如果你不需要，也可以删掉）
    DetectorFieldWriter(DetectorFieldWriter&& other) noexcept;
    DetectorFieldWriter& operator=(DetectorFieldWriter&& other) noexcept;

    void open(const std::string& filename,
              int nx, int ny,
              const std::vector<double>& x,
              const std::vector<double>& y,
              double zLabFixed,
              bool useFloat32 = true);

    void append(double time,
                double zBox,
                const float* ex,
                const float* ey,
                const float* ez,
                const float* bx,
                const float* by,
                const float* bz);

    void append(double time,
                double zBox,
                const double* ex,
                const double* ey,
                const double* ez,
                const double* bx,
                const double* by,
                const double* bz);

    void close();
    void flush();

    bool isOpen() const { return isOpen_; }
    std::size_t nFrames() const { return nFrames_; }

  private:
    void resetHandles_();
    void closeDataset_(hid_t& dset);
    void closeSpace_(hid_t& space);
    void closeProp_(hid_t& prop);
    void closeFile_(hid_t& file);

    hid_t createFixed1DDataset_(const char* name,
                                hid_t dtype,
                                hsize_t n,
                                const void* data);

    hid_t createExtendable1DDataset_(const char* name,
                                     hid_t dtype,
                                     hsize_t chunk0);

    hid_t createExtendable3DDataset_(const char* name,
                                     hid_t dtype,
                                     hsize_t chunk0,
                                     hsize_t chunk1,
                                     hsize_t chunk2);

    void writeScalarAttributeDouble_(hid_t obj,
                                     const char* name,
                                     double value);

    void writeScalarAttributeInt_(hid_t obj,
                                  const char* name,
                                  int value);

    void writeStringAttribute_(hid_t obj,
                               const char* name,
                               const std::string& value);

    template<typename T>
    void appendFrameImpl_(double time,
                          double zBox,
                          const T* ex,
                          const T* ey,
                          const T* ez,
                          const T* bx,
                          const T* by,
                          const T* bz,
                          hid_t memType);

  private:
    bool isOpen_;
    bool useFloat32_;

    int nx_;
    int ny_;
    std::size_t nFrames_;

    hid_t file_;

    hid_t dsetX_;
    hid_t dsetY_;
    hid_t dsetTime_;
    hid_t dsetZBox_;

    hid_t dsetEx_;
    hid_t dsetEy_;
    hid_t dsetEz_;
    hid_t dsetBx_;
    hid_t dsetBy_;
    hid_t dsetBz_;
  };
}

#endif