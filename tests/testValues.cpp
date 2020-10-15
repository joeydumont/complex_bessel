/*! ------------------------------------------------------------------------- *
 * \file testValues.cpp                                                       *
 *                                                                            *
 * \author Joey Dumont    <joey.dumont@gmail.com>                             *
 * \author Denis Gagnon   <gagnon88@gmail.com>                                *
 * \since 2014-11-26                                                          *
 * \date  2014-12-03                                                          *
 *                                                                            *
 * \brief Compares the values of Bessel functions computed with other imple-  *
 *        implementations.                                                    *
 *                                                                            *
 * This compares the values of the Bessel functions computed via the          *
 * complex_bessel library to the values computed via other implementations.   *
 * The values must be stored in a HDF5 file with the following hierarchy:     *
 *  /                                                                         *
 *    /besselJ                                                                *
 *    /besselY                                                                *
 *    /besselI                                                                *
 *    /besselK                                                                *
 *    /hankelH1                                                               *
 *    /hankelH2                                                               *
 * where the names represent different datasets, not groups. Each dataset     *
 * contains a 4-dimensional array of complex numbers. The indices refer to:   *
 *  - the order of differentiation;                                           *
 *  - the order of the Bessel function;                                       *
 *  - the real part of the argument of the Bessel function;                   *
 *  - the imaginary part of the argument of the Bessel function.              *
 * Since the three-dimensional grid we use to evaluate the Bessel functions   *
 * (order, Re(z), Im(z)) have equally spaced nodes, we only store the min/max *
 * values of each parameter as an attribute in each dataset. We store them as *
 *  - vmax: maximum order we consider, the range is actually (-vmax, vmax);   *
 *  - zrmin/zrmax: min/max values of real part of argument z.                 *
 *  - zimin/zimax: min/max values of imaginary part of argument z.            *
 * --------------------------------------------------------------------------*/

#include <armadillo>
#include <complex_bessel.h>
#include <hdf5.h>
#include <mpi.h>
#include <string>
#include <vector>

/*! Proper definition of modular arithmetic.
 *    @param[in] value Value of which we wish to take the modulo.
 *    @param[in] modulo Module of the arithmetic.
 *    @retval user_mod value mod modulo. */
inline double
user_mod(double value, double modulo)
{
  return value - modulo * floor(value / modulo);
} // user_mod()

int
main(int argc, char* argv[])
{
  // Library initilization.
  MPI_Init(&argc, &argv);
  herr_t      H5open();
  std::string filename = argv[1];

  // We determine the size and ranks of our process.
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // We distribute the indices among the processors.
  // For proper load balancing, we should have size % 6 = 0.
  std::vector<std::vector<unsigned int>> indices;
  indices.resize(size);

  for (unsigned int i = 0; i < 6; i++) {
    int idx = user_mod(i, size);
    indices[idx].push_back(i);
  }

  // Open existing file.
  hid_t plist_id;
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);

  // We store the name of the datasets we want to open.
  std::vector<std::string> dataset_names;
  dataset_names.push_back(std::string("/besselJ"));
  dataset_names.push_back(std::string("/besselY"));
  dataset_names.push_back(std::string("/besselI"));
  dataset_names.push_back(std::string("/besselK"));
  dataset_names.push_back(std::string("/hankelH1"));
  dataset_names.push_back(std::string("/hankelH2"));

  // We store the appropriate functions pointers in an array.
  std::vector<std::complex<double> (*)(double, std::complex<double>, int)> f_ptr;
  f_ptr.push_back(sp_bessel::besselJp);
  f_ptr.push_back(sp_bessel::besselYp);
  f_ptr.push_back(sp_bessel::besselIp);
  f_ptr.push_back(sp_bessel::besselKp);
  f_ptr.push_back(sp_bessel::hankelH1p);
  f_ptr.push_back(sp_bessel::hankelH2p);

  // We loop over the datasets.
  for (auto iter = indices[rank].begin(); iter != indices[rank].end(); iter++) {
    // Open dataset.
    hid_t dataset_id = H5Dopen(file_id, dataset_names[*iter].c_str(), H5P_DEFAULT);

    // Obtain the dataspace
    hid_t dspace = H5Dget_space(dataset_id);

    // We obtain the dimensions of the dataset.
    const int ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t   dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);

    // We read the dataset.
    std::complex<double> values[dims[0]][dims[1]][dims[2]][dims[3]];
    hid_t                complex_id = H5Dget_type(dataset_id);
    H5Dread(dataset_id, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);

    // We now open/read the attributes.
    double vmax, zimin, zimax, zrmin, zrmax;
    hid_t  vmax_id  = H5Aopen(dataset_id, "vmax", H5P_DEFAULT);
    hid_t  zimin_id = H5Aopen(dataset_id, "zimin", H5P_DEFAULT);
    hid_t  zimax_id = H5Aopen(dataset_id, "zrmax", H5P_DEFAULT);
    hid_t  zrmin_id = H5Aopen(dataset_id, "zrmin", H5P_DEFAULT);
    hid_t  zrmax_id = H5Aopen(dataset_id, "zrmax", H5P_DEFAULT);

    H5Aread(vmax_id, H5T_NATIVE_DOUBLE, &vmax);
    H5Aread(zimin_id, H5T_NATIVE_DOUBLE, &zimin);
    H5Aread(zimax_id, H5T_NATIVE_DOUBLE, &zimax);
    H5Aread(zrmin_id, H5T_NATIVE_DOUBLE, &zrmin);
    H5Aread(zrmax_id, H5T_NATIVE_DOUBLE, &zrmax);

    // We now evaluate the Bessel functions at the computed values.
    arma::vec orders = arma::linspace(-vmax, vmax, dims[1]);
    arma::vec realZ  = arma::linspace(zrmin, zrmax, dims[2]);
    arma::vec imagZ  = arma::linspace(zimin, zimax, dims[3]);

    unsigned int count = 0;
    for (int i = 0; i < dims[0]; i++) {
      for (int j = 0; j < dims[1]; j++) {
        for (int k = 0; k < dims[2]; k++) {
          for (int l = 0; l < dims[3]; l++) {
            double eps = std::abs(
              f_ptr[*iter](orders(j), std::complex<double>(realZ(k), imagZ(l)), i) -
              values[i][j][k][l]);

            if (eps > 1.0e-13) {
              std::cout << "Issue in " << dataset_names[*iter].c_str() << " at nu = " << orders(j)
                        << " z = " << realZ(k) << " + i" << imagZ(l) << "."
                        << " and p = " << i << std::endl;
              std::cout << "Epsilon is " << eps << std::endl;
              count++;
            }
          }
        }
      }
    }
  }

  // Library closures
  herr_t H5close();
  MPI_Finalize();

  return 0;
}