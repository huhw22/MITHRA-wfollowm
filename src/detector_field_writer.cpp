#include "detector_field_writer.h"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstring>

namespace MITHRA
{
  namespace
  {
    inline void throwIfNeg(herr_t status, const std::string& msg)
    {
      if (status < 0)
        throw std::runtime_error(msg);
    }

    inline void throwIfInvalid(hid_t id, const std::string& msg)
    {
      if (id < 0)
        throw std::runtime_error(msg);
    }
  }

  DetectorFieldWriter::DetectorFieldWriter()
    : isOpen_(false),
      useFloat32_(true),
      nx_(0),
      ny_(0),
      nFrames_(0),
      file_(-1),
      dsetX_(-1),
      dsetY_(-1),
      dsetTime_(-1),
      dsetZBox_(-1),
      dsetEx_(-1),
      dsetEy_(-1),
      dsetEz_(-1),
      dsetBx_(-1),
      dsetBy_(-1),
      dsetBz_(-1)
  {}

  DetectorFieldWriter::~DetectorFieldWriter()
  {
    close();
  }

  DetectorFieldWriter::DetectorFieldWriter(DetectorFieldWriter&& other) noexcept
  {
    *this = std::move(other);
  }

  DetectorFieldWriter& DetectorFieldWriter::operator=(DetectorFieldWriter&& other) noexcept
  {
    if (this != &other)
    {
      close();

      isOpen_     = other.isOpen_;
      useFloat32_ = other.useFloat32_;
      nx_         = other.nx_;
      ny_         = other.ny_;
      nFrames_    = other.nFrames_;

      file_       = other.file_;
      dsetX_      = other.dsetX_;
      dsetY_      = other.dsetY_;
      dsetTime_   = other.dsetTime_;
      dsetZBox_   = other.dsetZBox_;
      dsetEx_     = other.dsetEx_;
      dsetEy_     = other.dsetEy_;
      dsetEz_     = other.dsetEz_;
      dsetBx_     = other.dsetBx_;
      dsetBy_     = other.dsetBy_;
      dsetBz_     = other.dsetBz_;

      other.resetHandles_();
      other.isOpen_  = false;
      other.nx_      = 0;
      other.ny_      = 0;
      other.nFrames_ = 0;
    }
    return *this;
  }

  void DetectorFieldWriter::resetHandles_()
  {
    file_     = -1;
    dsetX_    = -1;
    dsetY_    = -1;
    dsetTime_ = -1;
    dsetZBox_ = -1;
    dsetEx_   = -1;
    dsetEy_   = -1;
    dsetEz_   = -1;
    dsetBx_   = -1;
    dsetBy_   = -1;
    dsetBz_   = -1;
  }

  void DetectorFieldWriter::closeDataset_(hid_t& dset)
  {
    if (dset >= 0)
    {
      H5Dclose(dset);
      dset = -1;
    }
  }

  void DetectorFieldWriter::closeSpace_(hid_t& space)
  {
    if (space >= 0)
    {
      H5Sclose(space);
      space = -1;
    }
  }

  void DetectorFieldWriter::closeProp_(hid_t& prop)
  {
    if (prop >= 0)
    {
      H5Pclose(prop);
      prop = -1;
    }
  }

  void DetectorFieldWriter::closeFile_(hid_t& file)
  {
    if (file >= 0)
    {
      H5Fclose(file);
      file = -1;
    }
  }

  hid_t DetectorFieldWriter::createFixed1DDataset_(const char* name,
                                                   hid_t dtype,
                                                   hsize_t n,
                                                   const void* data)
  {
    hsize_t dims[1] = { n };
    hid_t space = H5Screate_simple(1, dims, NULL);
    throwIfInvalid(space, std::string("Failed to create dataspace for ") + name);

    hid_t dset = H5Dcreate2(file_, name, dtype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0)
    {
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create dataset ") + name);
    }

    throwIfNeg(H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
               std::string("Failed to write dataset ") + name);

    H5Sclose(space);
    return dset;
  }

  hid_t DetectorFieldWriter::createExtendable1DDataset_(const char* name,
                                                        hid_t dtype,
                                                        hsize_t chunk0)
  {
    hsize_t dims[1]     = { 0 };
    hsize_t maxdims[1]  = { H5S_UNLIMITED };
    hsize_t chunks[1]   = { chunk0 };

    hid_t space = H5Screate_simple(1, dims, maxdims);
    throwIfInvalid(space, std::string("Failed to create dataspace for ") + name);

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0)
    {
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create dcpl for ") + name);
    }

    throwIfNeg(H5Pset_chunk(dcpl, 1, chunks),
               std::string("Failed to set chunk for ") + name);

    // 开启轻量压缩；如果你那边 HDF5 没有 zlib，可先删掉这行
    H5Pset_deflate(dcpl, 4);

    hid_t dset = H5Dcreate2(file_, name, dtype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (dset < 0)
    {
      H5Pclose(dcpl);
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create dataset ") + name);
    }

    H5Pclose(dcpl);
    H5Sclose(space);
    return dset;
  }

  hid_t DetectorFieldWriter::createExtendable3DDataset_(const char* name,
                                                        hid_t dtype,
                                                        hsize_t chunk0,
                                                        hsize_t chunk1,
                                                        hsize_t chunk2)
  {
    hsize_t dims[3]     = { 0, static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_) };
    hsize_t maxdims[3]  = { H5S_UNLIMITED, static_cast<hsize_t>(ny_), static_cast<hsize_t>(nx_) };
    hsize_t chunks[3]   = { chunk0, chunk1, chunk2 };

    hid_t space = H5Screate_simple(3, dims, maxdims);
    throwIfInvalid(space, std::string("Failed to create dataspace for ") + name);

    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl < 0)
    {
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create dcpl for ") + name);
    }

    throwIfNeg(H5Pset_chunk(dcpl, 3, chunks),
               std::string("Failed to set chunk for ") + name);

    H5Pset_deflate(dcpl, 4);

    hid_t dset = H5Dcreate2(file_, name, dtype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (dset < 0)
    {
      H5Pclose(dcpl);
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create dataset ") + name);
    }

    H5Pclose(dcpl);
    H5Sclose(space);
    return dset;
  }

  void DetectorFieldWriter::writeScalarAttributeDouble_(hid_t obj,
                                                        const char* name,
                                                        double value)
  {
    hsize_t dims[1] = { 1 };
    hid_t space = H5Screate_simple(1, dims, NULL);
    throwIfInvalid(space, std::string("Failed to create attr space for ") + name);

    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0)
    {
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create attribute ") + name);
    }

    throwIfNeg(H5Awrite(attr, H5T_NATIVE_DOUBLE, &value),
               std::string("Failed to write attribute ") + name);

    H5Aclose(attr);
    H5Sclose(space);
  }

  void DetectorFieldWriter::writeScalarAttributeInt_(hid_t obj,
                                                     const char* name,
                                                     int value)
  {
    hsize_t dims[1] = { 1 };
    hid_t space = H5Screate_simple(1, dims, NULL);
    throwIfInvalid(space, std::string("Failed to create attr space for ") + name);

    hid_t attr = H5Acreate2(obj, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0)
    {
      H5Sclose(space);
      throw std::runtime_error(std::string("Failed to create attribute ") + name);
    }

    throwIfNeg(H5Awrite(attr, H5T_NATIVE_INT, &value),
               std::string("Failed to write attribute ") + name);

    H5Aclose(attr);
    H5Sclose(space);
  }

  void DetectorFieldWriter::writeStringAttribute_(hid_t obj,
                                                  const char* name,
                                                  const std::string& value)
  {
    hid_t type = H5Tcopy(H5T_C_S1);
    throwIfInvalid(type, std::string("Failed to create string type for ") + name);

    throwIfNeg(H5Tset_size(type, value.size()),
               std::string("Failed to set string size for ") + name);

    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0)
    {
      H5Tclose(type);
      throw std::runtime_error(std::string("Failed to create attr space for ") + name);
    }

    hid_t attr = H5Acreate2(obj, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0)
    {
      H5Sclose(space);
      H5Tclose(type);
      throw std::runtime_error(std::string("Failed to create attribute ") + name);
    }

    const char* s = value.c_str();
    throwIfNeg(H5Awrite(attr, type, s),
               std::string("Failed to write attribute ") + name);

    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(type);
  }

  void DetectorFieldWriter::open(const std::string& filename,
                                 int nx, int ny,
                                 const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 double zLabFixed,
                                 bool useFloat32)
  {
    if (isOpen_)
      throw std::runtime_error("DetectorFieldWriter::open called while file is already open.");

    if (nx <= 0 || ny <= 0)
      throw std::runtime_error("DetectorFieldWriter::open got non-positive nx or ny.");

    if (static_cast<int>(x.size()) != nx)
      throw std::runtime_error("DetectorFieldWriter::open x.size() != nx.");

    if (static_cast<int>(y.size()) != ny)
      throw std::runtime_error("DetectorFieldWriter::open y.size() != ny.");

    useFloat32_ = useFloat32;
    nx_         = nx;
    ny_         = ny;
    nFrames_    = 0;

    file_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    throwIfInvalid(file_, "Failed to create HDF5 file: " + filename);

    try
    {
      dsetX_    = createFixed1DDataset_("x", H5T_NATIVE_DOUBLE, static_cast<hsize_t>(nx), x.data());
      dsetY_    = createFixed1DDataset_("y", H5T_NATIVE_DOUBLE, static_cast<hsize_t>(ny), y.data());
      dsetTime_ = createExtendable1DDataset_("time", H5T_NATIVE_DOUBLE, 256);
      dsetZBox_ = createExtendable1DDataset_("z_box", H5T_NATIVE_DOUBLE, 256);

      const hid_t fieldType = useFloat32_ ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;

      dsetEx_ = createExtendable3DDataset_("Ex", fieldType, 1, static_cast<hsize_t>(ny), static_cast<hsize_t>(nx));
      dsetEy_ = createExtendable3DDataset_("Ey", fieldType, 1, static_cast<hsize_t>(ny), static_cast<hsize_t>(nx));
      dsetEz_ = createExtendable3DDataset_("Ez", fieldType, 1, static_cast<hsize_t>(ny), static_cast<hsize_t>(nx));
      dsetBx_ = createExtendable3DDataset_("Bx", fieldType, 1, static_cast<hsize_t>(ny), static_cast<hsize_t>(nx));
      dsetBy_ = createExtendable3DDataset_("By", fieldType, 1, static_cast<hsize_t>(ny), static_cast<hsize_t>(nx));
      dsetBz_ = createExtendable3DDataset_("Bz", fieldType, 1, static_cast<hsize_t>(ny), static_cast<hsize_t>(nx));

      writeStringAttribute_(file_, "detector_type", "plane");
      writeStringAttribute_(file_, "storage_dtype", useFloat32_ ? "float32" : "float64");
      writeScalarAttributeInt_(file_, "nx", nx_);
      writeScalarAttributeInt_(file_, "ny", ny_);
      writeScalarAttributeDouble_(file_, "z_lab_fixed", zLabFixed);
    }
    catch (...)
    {
      close();
      throw;
    }

    isOpen_ = true;
  }

  template<typename T>
  void DetectorFieldWriter::appendFrameImpl_(double time,
                                             double zBox,
                                             const T* ex,
                                             const T* ey,
                                             const T* ez,
                                             const T* bx,
                                             const T* by,
                                             const T* bz,
                                             hid_t memType)
  {
    if (!isOpen_)
      throw std::runtime_error("DetectorFieldWriter::append called before open().");

    if (!ex || !ey || !ez || !bx || !by || !bz)
      throw std::runtime_error("DetectorFieldWriter::append got null field pointer.");

    const hsize_t oldNt = static_cast<hsize_t>(nFrames_);
    const hsize_t newNt = oldNt + 1;

    // -------- time --------
    {
      hsize_t newDims[1] = { newNt };
      throwIfNeg(H5Dset_extent(dsetTime_, newDims), "Failed to extend dataset time.");

      hid_t fileSpace = H5Dget_space(dsetTime_);
      throwIfInvalid(fileSpace, "Failed to get file space for time.");

      hsize_t start[1] = { oldNt };
      hsize_t count[1] = { 1 };
      throwIfNeg(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, NULL, count, NULL),
                 "Failed to select hyperslab for time.");

      hid_t memSpace = H5Screate_simple(1, count, NULL);
      throwIfInvalid(memSpace, "Failed to create mem space for time.");

      double value = time;
      throwIfNeg(H5Dwrite(dsetTime_, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, &value),
                 "Failed to write time.");

      H5Sclose(memSpace);
      H5Sclose(fileSpace);
    }

    // -------- z_box --------
    {
      hsize_t newDims[1] = { newNt };
      throwIfNeg(H5Dset_extent(dsetZBox_, newDims), "Failed to extend dataset z_box.");

      hid_t fileSpace = H5Dget_space(dsetZBox_);
      throwIfInvalid(fileSpace, "Failed to get file space for z_box.");

      hsize_t start[1] = { oldNt };
      hsize_t count[1] = { 1 };
      throwIfNeg(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, NULL, count, NULL),
                 "Failed to select hyperslab for z_box.");

      hid_t memSpace = H5Screate_simple(1, count, NULL);
      throwIfInvalid(memSpace, "Failed to create mem space for z_box.");

      double value = zBox;
      throwIfNeg(H5Dwrite(dsetZBox_, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, &value),
                 "Failed to write z_box.");

      H5Sclose(memSpace);
      H5Sclose(fileSpace);
    }

    // -------- field helper lambda --------
    auto writeFrame3D = [&](hid_t dset, const T* buf, const char* name)
    {
      hsize_t newDims[3] = { newNt,
                             static_cast<hsize_t>(ny_),
                             static_cast<hsize_t>(nx_) };
      throwIfNeg(H5Dset_extent(dset, newDims),
                 std::string("Failed to extend dataset ") + name);

      hid_t fileSpace = H5Dget_space(dset);
      throwIfInvalid(fileSpace, std::string("Failed to get file space for ") + name);

      hsize_t start[3] = { oldNt, 0, 0 };
      hsize_t count[3] = { 1,
                           static_cast<hsize_t>(ny_),
                           static_cast<hsize_t>(nx_) };

      throwIfNeg(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, NULL, count, NULL),
                 std::string("Failed to select hyperslab for ") + name);

      hid_t memSpace = H5Screate_simple(3, count, NULL);
      throwIfInvalid(memSpace, std::string("Failed to create mem space for ") + name);

      throwIfNeg(H5Dwrite(dset, memType, memSpace, fileSpace, H5P_DEFAULT, buf),
                 std::string("Failed to write field dataset ") + name);

      H5Sclose(memSpace);
      H5Sclose(fileSpace);
    };

    writeFrame3D(dsetEx_, ex, "Ex");
    writeFrame3D(dsetEy_, ey, "Ey");
    writeFrame3D(dsetEz_, ez, "Ez");
    writeFrame3D(dsetBx_, bx, "Bx");
    writeFrame3D(dsetBy_, by, "By");
    writeFrame3D(dsetBz_, bz, "Bz");

    ++nFrames_;
  }

  void DetectorFieldWriter::append(double time,
                                   double zBox,
                                   const float* ex,
                                   const float* ey,
                                   const float* ez,
                                   const float* bx,
                                   const float* by,
                                   const float* bz)
  {
    if (!useFloat32_)
      throw std::runtime_error("DetectorFieldWriter::append(float*) called but file storage is float64.");

    appendFrameImpl_(time, zBox, ex, ey, ez, bx, by, bz, H5T_NATIVE_FLOAT);
  }

  void DetectorFieldWriter::append(double time,
                                   double zBox,
                                   const double* ex,
                                   const double* ey,
                                   const double* ez,
                                   const double* bx,
                                   const double* by,
                                   const double* bz)
  {
    if (useFloat32_)
      throw std::runtime_error("DetectorFieldWriter::append(double*) called but file storage is float32.");

    appendFrameImpl_(time, zBox, ex, ey, ez, bx, by, bz, H5T_NATIVE_DOUBLE);
  }

  void DetectorFieldWriter::close()
  {
    closeDataset_(dsetBz_);
    closeDataset_(dsetBy_);
    closeDataset_(dsetBx_);
    closeDataset_(dsetEz_);
    closeDataset_(dsetEy_);
    closeDataset_(dsetEx_);
    closeDataset_(dsetZBox_);
    closeDataset_(dsetTime_);
    closeDataset_(dsetY_);
    closeDataset_(dsetX_);
    closeFile_(file_);

    isOpen_  = false;
    nx_      = 0;
    ny_      = 0;
    nFrames_ = 0;
    useFloat32_ = true;
  }
}