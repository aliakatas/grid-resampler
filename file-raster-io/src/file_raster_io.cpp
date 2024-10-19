#include "file_raster_io.h"
#include "common_maths_and_geometry_ops.h"

#ifdef __WIN32
#include "ogrsf_frmts/ogrsf_frmts.h"
#elif __linux__
#include "ogrsf_frmts.h"
#endif 

#include <string>
#include <vector>
#include <filesystem>

#define DEFAULT_DEFLATION_LEVEL		4
#define DEFLATE_STRING		        "DEFLATE"
#define LZW_STRING			        "LZW"
#define LZMA_STRING					"LZMA"

// These are the projections tested to work
#define GTIFF_PROJ_SIZE				4
const char* gtiff_projections[GTIFF_PROJ_SIZE] = {
	"NAD27",
	"NAD83",
	"WGS72",
	"WGS84"
};

namespace file_raster_io
{
	std::vector<std::string> messages;
}

#if __cplusplus
extern "C" {
#endif

	//#######################################################################
	DllExport void get_file_raster_io_mod_errors(char* buffer, int* nchars, int* err)
	{
		// this should be a failsafe...
		std::string temp;
		for (auto& s : file_raster_io::messages)
			temp = s + std::string("\n");

		if (buffer)
		{
			strcpy(buffer, temp.c_str());
			file_raster_io::messages.clear();
		}
		else
		{
			*nchars = temp.length() + 1;
		}
		*err = FILE_RASTER_IO_EXIT_SUCCESS;
	}

	//#######################################################################
	DllExport void initialise_file_raster_io()
	{
		GDALAllRegister();
	}

	//#######################################################################
	DllExport void write_file_raster(const char* fname,
		const int nrows, const int ncolumns,
		const double dx, const double dy, const double angle_deg,
		const double xll, const double yll,
		const double nodata,
		const int item_byte_size, const int is_float,
		void* data,
		const char* units,
		const char* compression, const int deflation_level,
		const char* projection, const int init_gdal_drivers,
		int* err) 
	{
		// Presume the worst
		*err = FILE_RASTER_IO_EXIT_ERROR;
		
		try
		{
			if (init_gdal_drivers == 1)
				// This is not ideal... Best to be called centrally!
				GDALAllRegister();
			
			std::string file_extension = std::filesystem::path(fname).extension().string();
			fm_string::convert_to_lowercase(file_extension);

			std::string pszFormat;
			char** papszOptions = NULL;
			if (file_extension.compare(".tif") == 0)
			{
				pszFormat = "GTiff";
				papszOptions = CSLSetNameValue(papszOptions, "TFW", "YES");

				if (compression)
				{
					if (strcasecmp(compression, "DEFLATE") == 0)
					{
						papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
						int def_lev = deflation_level;
						if (deflation_level < 1 || deflation_level > 12)
						{
							std::string temp = std::string("Compression level ") + std::to_string(deflation_level) + std::string(" was ignored - using 6 instead");
							file_raster_io::messages.push_back(temp.c_str());
							def_lev = 6;
						}
						papszOptions = CSLSetNameValue(papszOptions, "ZLEVEL", std::to_string(def_lev).c_str());
					}
					else if (strcasecmp(compression, "LZW") == 0)
					{
						papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
					}
					else
					{
						std::string temp = std::string("Compression option ") + std::string(compression) + std::string(" was ignored - using DEFLATE instead");
						file_raster_io::messages.push_back(temp.c_str());
						papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
						papszOptions = CSLSetNameValue(papszOptions, "ZLEVEL", "6");
					}

					papszOptions = CSLSetNameValue(papszOptions, "NUM_THREADS", "2");
				}
			}
			else if (file_extension.compare(".asc") == 0)
			{
				pszFormat = "AAIGrid";
			}
			/*else if (file_extension.compare(".hdr") == 0 || file_extension.compare(".flt") == 0)
			{

			}*/
			else
			{
				std::string temp{ "Failed to recognise output format by its extension: " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat.c_str());
			if (poDriver == NULL)
			{
				CPLFree(papszOptions);
				std::string temp{ "Failed to open driver for output raster dataset: " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			GDALDataType data_type;
			if (is_float)
			{
				if (item_byte_size == sizeof(float))
					data_type = GDT_Float32;
				else if (item_byte_size == sizeof(double))
					data_type = GDT_Float64;
				else {
					CPLFree(papszOptions);
					std::string temp{ "Floating point output data type is not supported for " };
					temp += std::string(fname);
					throw std::runtime_error(temp.c_str());
				}
			}
			else
			{
				if (item_byte_size == sizeof(int))
					data_type = GDT_Int32;
				else if (item_byte_size == sizeof(long long))
					data_type = GDT_Int64;
				else {
					CPLFree(papszOptions);
					std::string temp{ "Integer output data type is not supported for " };
					temp += std::string(fname);
					throw std::runtime_error(temp.c_str());
				}
			}

			// Let's create the dataset
			GDALDataset* poDstDS = poDriver->Create(fname, ncolumns, nrows, 1, data_type, papszOptions);
			if (poDstDS == NULL)
			{
				CPLFree(papszOptions);
				std::string temp{ "Failed to create dataset for " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			double adfGeoTransform[6];
			double xtl = xll;
			double ytl = yll + dy * nrows;
			double rotation_matrix[4];
			math_ops::calculate_rotation_matrix(math_ops::deg_to_rad(angle_deg), rotation_matrix);
			math_ops::rotate_point_about(&xtl, &ytl, rotation_matrix, xll, yll);
			math_ops::set_affine_transformation(adfGeoTransform, xtl, ytl, dx, dy, angle_deg);
			
			if (poDstDS->SetGeoTransform(adfGeoTransform) != CE_None)
			{
				GDALClose((GDALDatasetH)poDstDS);		// if there is an error here, we don't care
				CPLFree(papszOptions);
				std::string temp{ "Failed to set geotransformation for " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}
			
			/*OGRSpatialReference oSRS;
			char* pszSRS_WKT = NULL;
			oSRS.SetUTM(11, TRUE);
			oSRS.SetWellKnownGeogCS("NAD27");
			oSRS.exportToWkt(&pszSRS_WKT);
			poDstDS->SetProjection(pszSRS_WKT);
			CPLFree(pszSRS_WKT);
			*/

			GDALRasterBand* poBand = poDstDS->GetRasterBand(1);
			poBand->SetNoDataValue(nodata);
			if (poBand->RasterIO(GF_Write, 0, 0, ncolumns, nrows,
				data, ncolumns, nrows, data_type, 0, 0) != CE_None)
			{
				GDALClose((GDALDatasetH)poDstDS);		// if there is an error here, we don't care
				CPLFree(papszOptions);
				std::string temp{ "Failed to write grid data values to " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			/* Once we're done, close properly the dataset */
			if (GDALClose((GDALDatasetH)poDstDS) != CE_None)
			{
				CPLFree(papszOptions);
				std::string temp{ "Failed to close raster dataset " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}
			CPLFree(papszOptions);
		}
		catch (const std::exception& e) {
			file_raster_io::messages.push_back(e.what());
			*err = FILE_RASTER_IO_EXIT_ERROR;
			return;
		}
		*err = FILE_RASTER_IO_EXIT_SUCCESS;
	}

	//#######################################################################
	DllExport void read_file_raster(const char* fname,
		int* nrows, int* ncolumns,
		double* dx, double* dy, double* angle_deg,
		double* xll, double* yll,
		double* nodata,
		int* item_byte_size, int* is_float,
		void* data,
		char* units,
		char* projection, const int init_gdal_drivers,
		int* err)
	{
		// Presume the worst
		*err = FILE_RASTER_IO_EXIT_ERROR;

		try
		{
			if (init_gdal_drivers == 1)
				// This is not ideal... Best to be called centrally!
				GDALAllRegister();

			const GDALAccess eAccess = GA_ReadOnly;
			GDALDatasetUniquePtr poDataset = GDALDatasetUniquePtr(GDALDataset::FromHandle(GDALOpen(fname, eAccess)));
			if (!poDataset)
			{
				std::string temp{ "Failed to open raster dataset: " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			if (poDataset->GetRasterCount() != 1)
			{
				std::string temp{ "File is not supported - raster dataset must contain 1 band: " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			*nrows = poDataset->GetRasterYSize();
			*ncolumns = poDataset->GetRasterXSize();
			
			/*if (poDataset->GetProjectionRef() != NULL)
				printf("Projection is `%s'\n", poDataset->GetProjectionRef());*/
			double adfGeoTransform[6];
			if (poDataset->GetGeoTransform(adfGeoTransform) != CE_None)
			{
				std::string temp{ "Failed to get geotransformation from raster dataset: " };
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			double xtl, ytl, angle_r;
			math_ops::decrypt_affine_transformation_matrix(&xtl, &ytl, dx, dy, &angle_r, adfGeoTransform);
			*angle_deg = math_ops::rad_to_deg(angle_r);

			math_ops::apply_affine_transformation(xll, yll, *nrows + 1, 0, adfGeoTransform);

			GDALRasterBand* poBand = poDataset->GetRasterBand(1);
			*nodata = poBand->GetNoDataValue();

			GDALDataType data_type = poBand->GetRasterDataType();
			switch (data_type)
			{
			case GDT_Float32:
			case GDT_CFloat32:
				*is_float = 1;
				*item_byte_size = sizeof(float);
				break;
			case GDT_Float64:
			case GDT_CFloat64:
				*is_float = 1;
				*item_byte_size = sizeof(double);
				break;
			case GDT_Int32:
			case GDT_CInt32:
				*is_float = 0;
				*item_byte_size = sizeof(int);
				break;
			case GDT_Int64:
				*is_float = 0;
				*item_byte_size = sizeof(long long);
				break;
			default:
				std::string temp{ "Data type " };
				temp += std::string(GDALGetDataTypeName(data_type));
				temp += std::string(" not supported for raster dataset: ");
				temp += std::string(fname);
				throw std::runtime_error(temp.c_str());
			}

			if (data)
			{
				if (poBand->RasterIO(GF_Read, 0, 0, *ncolumns, *nrows, data, *ncolumns, *nrows, data_type, 0, 0) != CE_None)
				{
					std::string temp{ "Failed to get grid data values from raster dataset: " };
					temp += std::string(fname);
					throw std::runtime_error(temp.c_str());
				}
			}
			// since the dataset is unique ptr, it will be closed/destroyed here
		}
		catch (const std::exception& e) {
			file_raster_io::messages.push_back(e.what());
			*err = FILE_RASTER_IO_EXIT_ERROR;
			return;
		}

		*err = FILE_RASTER_IO_EXIT_SUCCESS;
	}

#if __cplusplus
}
#endif
