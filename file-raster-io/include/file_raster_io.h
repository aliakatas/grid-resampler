#ifndef FILE_RASTER_IO_H
#define FILE_RASTER_IO_H

#define FILE_RASTER_IO_EXIT_SUCCESS			0
#define FILE_RASTER_IO_EXIT_ERROR			-1
#define FILE_RASTER_IO_EXIT_WARNING			1

#ifdef _DLLBUILD
#define DllExport   __declspec( dllexport )
#else
#define DllExport 
#endif

#if __cplusplus
extern "C" {
#endif

	/**
	 * @brief Retrieves the messages from the File_Raster_IO module.
	 * Once the messages are returned, the message buffer is cleared.
	 *
	 * @param buffer [out] The buffer to copy the messages to. It can be NULL so that the number of character needed is returned.
	 * @param nchars [inout] The number of characters needed for the buffer (if buffer is NULL).
	 * @param err [out] Exit error code. 0 indicates success. Negative value indicates severe error - the program may not work as expected. Positive value indicates a warning.
	 */
	DllExport void get_file_raster_io_mod_errors(char* buffer, int* nchars, int* err);

	/**
	* @brief Interfacing function to request loading of all GDAL drivers.
	*/
	DllExport void initialise_file_raster_io();
	
	/**
	* @brief Creates a GeoTIFF/ASCII file and the accompanying "world file" (for GeoTIFF) 
	* for integer and floating point data (single/double precision supported).
	* 
	* The origin of the data is considered to be at the upper-left vertex of the upper-left cell.
	* That means that element (0,0) is located at the upper-left corner of the domain.
	* Any rotation is assumed to be about lower-left corner of the lower-left cell. 
	* ASCII files do NOT support rotated grids through this function.
	* 
	* @param fname [in] The full path to the file where data is written.
	* @param nrows [in] The number of rows of the array of data.
	* @param ncolumns [in] The number of columns of the array of data.
	* @param dx [in] The x-direction resolution of the data.
	* @param dy [in] The y-direction resolution of the data.
	* @param angle_deg [in] The angle of rotation about the lower-left corner in degrees. Ignored for ASCII.
	* @param xll [in] The x-coordinate of the lower-left corner of the domain.
	* @param yll [in] The y-coordinate of the lower-left corner of the domain.
	* @param nodata [in] The value representing "no data".
	* @param item_byte_size [in] The byte size of one element of data array. Eg sizeof(float), sizeof(double), sizeof(int) etc.
	* @param is_float [in] If set to 0, indicates working with integers, floating point numbers are assumed otherwise.
	* @param data [in] The data to be written.
	* @param units [in] A null-terminated string of the length unit. Ignored for ASCII.
	* @param compression [in] A null-terminated string to pass the name of the compression algorithm. DEFLATE and LZW are supported. For no compression, can be NULL. Ignored for ASCII.
	* @param deflation_level [in] If compression == DEFLATE, an integer between 1 (min compression, fast) and 9 (max compression, slower). Defaults to 6 if out-of-bounds. Ignored for ASCII.
	* @param projection [in] Should be set to NULL for now. This is placeholder to pass the null-terminated string for the projection. Ignored for ASCII.
	* @param init_gdal_drivers [in] Indicates whether to initialise GDAL drivers (1) or not (0).
	* @param err [out] Returns 0 if all operations succeeded. Negative indicates errors. Positive indicates a warning. 
	*/
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
		int* err);

	/**
	* @brief Reads data and attributes from a GeoTIFF/ASCII file. Integer and
	* floating point data (single/double precision) supported.
	*
	* @param fname [in] The full path to the file where data is written.
	* @param nrows [out] The number of rows of the array of data.
	* @param ncolumns [out] The number of columns of the array of data.
	* @param dx [out] The x-direction resolution of the data.
	* @param dy [out] The y-direction resolution of the data.
	* @param angle_deg [out] The angle of rotation about the lower-left corner in degrees. For ASCII, this is always 0.
	* @param xll [out] The x-coordinate of the lower-left corner of the domain.
	* @param yll [out] The y-coordinate of the lower-left corner of the domain.
	* @param nodata [out] The value representing "no data".
	* @param item_byte_size [out] The byte size of one element of data array. 
	* @param is_float [out] A value of 0 indicates integers, floating point numbers should be assumed otherwise.
	* @param data [out] The data array read from the file. 
	* @param units [out] A null-terminated string of the length unit. Ignored for ASCII. 
	* @param projection [out] Should be set to NULL for now. This is placeholder to retrieve the null-terminated string for the projection. Ignored for ASCII.
	* @param init_gdal_drivers [in] Indicates whether to initialise GDAL drivers (1) or not (0).
	* @param err [out] Returns 0 if all operations succeeded. Negative indicates errors. Positive indicates a warning.
	*/
	DllExport void read_file_raster(const char* fname,
		int* nrows, int* ncolumns,
		double* dx, double* dy, double* angle_deg,
		double* xll, double* yll,
		double* nodata,
		int* item_byte_size, int* is_float,
		void* data,
		char* units,
		char* projection, const int init_gdal_drivers,
		int* err);

#if __cplusplus
}
#endif

#endif // FILE_RASTER_IO_H
