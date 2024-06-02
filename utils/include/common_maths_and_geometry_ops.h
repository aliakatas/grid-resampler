#ifndef COMMON_MATHS_AND_GEOMETRY_OPS_H
#define COMMON_MATHS_AND_GEOMETRY_OPS_H

#include <cmath>
#ifdef __NVCC__
#include <device_launch_parameters.h>
#define CFUNC_DECORATION __host__ __device__
#else 
#define CFUNC_DECORATION
#endif

#define ONE_EIGHTY_DEG		    180.0
#define DEFAULT_ZERO_THRESHOLD	0.00001

//********************************************// 
namespace index_ops
{
    //##################################################################
	CFUNC_DECORATION inline int get_row_major_linear_index(const int irow, const int icol, const int ncols)
	{
		return irow * ncols + icol;
	}

	//##################################################################
	CFUNC_DECORATION inline int get_column_major_linear_index(const int irow, const int icol, const int nrows)
	{
		return icol * nrows + irow;
	}

	//##################################################################
	CFUNC_DECORATION inline void get_2D_indices_from_row_major_linear_index(int* irow, int* icol, const int ncols, const int idx)
	{
		*irow = idx / ncols;
		*icol = idx % ncols;
	}

	//##################################################################
	CFUNC_DECORATION inline void get_2D_indices_from_column_major_linear_index(int* irow, int* icol, const int nrows, const int idx)
	{
		*irow = idx % nrows;
		*icol = idx / nrows;
	}
}

//********************************************//
namespace maths_ops
{
    //##################################################################
	template <typename T>
	CFUNC_DECORATION int are_equal(const T val1, const T val2, const T tolerance = static_cast<T>(DEFAULT_ZERO_THRESHOLD))
	{
		static_assert(false, "Unsupported type for are_equal");
	}

    template <>
	CFUNC_DECORATION int are_equal(const float val1, const float val2, 
        const float tolerance)
	{
		return (std::abs(val1 - val2) < tolerance) ? 1 : 0;
	}

    template <>
	CFUNC_DECORATION int are_equal(const double val1, const double val2, 
        const double tolerance)
	{
		return (std::abs(val1 - val2) < tolerance) ? 1 : 0;
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION T deg_to_rad(const T d)
	{
		static_assert(false, "Unsupported type for deg_to_rad");
	}

    template <>
	CFUNC_DECORATION float deg_to_rad(const float d)
	{
		return static_cast<float>(M_PI) * d / static_cast<float>(ONE_EIGHTY_DEG);
	}

    template <>
	CFUNC_DECORATION double deg_to_rad(const double d)
	{
		return M_PI * d / ONE_EIGHTY_DEG;
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION T rad_to_deg(const T r)
	{
		static_assert(false, "Unsupported type for rad_to_deg");
	}

    template <>
	CFUNC_DECORATION float rad_to_deg(const float r)
	{
		return static_cast<float>(ONE_EIGHTY_DEG) * r / static_cast<float>(M_PI);
	}

    template <>
	CFUNC_DECORATION double rad_to_deg(const double r)
	{
		return ONE_EIGHTY_DEG * r / M_PI;
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION T dot_product(const T a_x, const T a_y, const T b_x, const T b_y)
	{
		static_assert(false, "Unsupported type for dot_product");
	}

    template <>
	CFUNC_DECORATION float dot_product(const float a_x, const float a_y, const float b_x, const float b_y)
	{
		return a_x * b_x + a_y * b_y;
	}

    template <>
	CFUNC_DECORATION double dot_product(const double a_x, const double a_y, const double b_x, const double b_y)
	{
		return a_x * b_x + a_y * b_y;
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION T cross_product(const T a_x, const T a_y, const T b_x, const T b_y)
	{
		static_assert(false, "Unsupported type for cross_product");
	}

    template <>
	CFUNC_DECORATION float cross_product(const float a_x, const float a_y, const float b_x, const float b_y)
	{
		return a_x * b_y - a_y * b_x;
	}

    template <>
	CFUNC_DECORATION double cross_product(const double a_x, const double a_y, const double b_x, const double b_y)
	{
		return a_x * b_y - a_y * b_x;
	}
}

//********************************************//
namespace geom_ops
{
	//##################################################################
	template <typename T>
	CFUNC_DECORATION void create_vector(T* vecx, T* vecy, 
        const T start_x, const T start_y,
        const T end_x, const T end_y)
	{
		static_assert(false, "Unsupported type for create_vector");
	}

    template <>
	CFUNC_DECORATION void create_vector(float* vecx, float* vecy, 
        const float start_x, const float start_y,
        const float end_x, const float end_y)
	{
		*vecx = end_x - start_x;
        *vecy = end_y - start_y;
	}

    template <>
	CFUNC_DECORATION void create_vector(double* vecx, double* vecy, 
        const double start_x, const double start_y,
        const double end_x, const double end_y)
	{
		*vecx = end_x - start_x;
        *vecy = end_y - start_y;
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION T vector_direction(const T x, const T y)
	{
		static_assert(false, "Unsupported type for vector_direction");
	}

    template <>
	CFUNC_DECORATION float vector_direction(const float x, const float y)
	{
		return std::atan2(y, x);
	}

    template <>
	CFUNC_DECORATION double vector_direction(const double x, const double y)
	{
		return std::atan2(y, x);
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION T vector_magnitude(const T x, const T y)
	{
		static_assert(false, "Unsupported type for vector_magnitude");
	}

    template <>
	CFUNC_DECORATION float vector_magnitude(const float x, const float y)
	{
        return std::sqrt(maths_ops::dot_product<float>(x, y, x, y));
	}

    template <>
	CFUNC_DECORATION double vector_magnitude(const double x, const double y)
	{
        return std::sqrt(maths_ops::dot_product<double>(x, y, x, y));
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION T squared_distance(const T x, const T y, const T x_origin, const T y_origin)
	{
		static_assert(false, "Unsupported type for squared_distance");
	}

    template <>
	CFUNC_DECORATION float squared_distance(const float x, const float y, const float x_origin, const float y_origin)
	{
		return (x - x_origin) * (x - x_origin) + (y - y_origin) * (y - y_origin);
	}

    template <>
	CFUNC_DECORATION double squared_distance(const double x, const double y, const double x_origin, const double y_origin)
	{
		return (x - x_origin) * (x - x_origin) + (y - y_origin) * (y - y_origin);
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION T distance(const T x, const T y, const T x_origin, const T y_origin)
	{
		static_assert(false, "Unsupported type for distance");
	}

    template <>
	CFUNC_DECORATION float distance(const float x, const float y, const float x_origin, const float y_origin)
	{
		return std::sqrt(squared_distance<float>(x, y, x_origin, y_origin));
	}

    template <>
	CFUNC_DECORATION double distance(const double x, const double y, const double x_origin, const double y_origin)
	{
		return std::sqrt(squared_distance<double>(x, y, x_origin, y_origin));
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION void get_line_coeff(T x1, T y1, T x2, T y2, T* a, T* b, T* c)
	{
		static_assert(false, "Unsupported type for get_line_coeff");
	}

    template <>
	CFUNC_DECORATION void get_line_coeff(float x1, float y1, float x2, float y2, float* a, float* b, float* c)
	{
		*a = y1 - y2;
		*b = x2 - x1;
		*c = x1 * y2 - x2 * y1;
	}

    template <>
	CFUNC_DECORATION void get_line_coeff(double x1, double y1, double x2, double y2, double* a, double* b, double* c)
	{
		*a = y1 - y2;
		*b = x2 - x1;
		*c = x1 * y2 - x2 * y1;
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION T distance_point_to_line(T x1, T y1, T x2, T y2, T xp, T yp)
	{
        static_assert(false, "Unsupported type for distance_point_to_line");
	}

    template <>
	CFUNC_DECORATION float distance_point_to_line(float x1, float y1, float x2, float y2, float xp, float yp)
	{
		float a, b, c;
		get_line_coeff<float>(x1, y1, x2, y2, &a, &b, &c);
		return std::abs(a * xp + b * yp + c) / std::sqrt(a * a + b * b);
	}

    template <>
	CFUNC_DECORATION double distance_point_to_line(double x1, double y1, double x2, double y2, double xp, double yp)
	{
		double a, b, c;
		get_line_coeff<double>(x1, y1, x2, y2, &a, &b, &c);
		return std::abs(a * xp + b * yp + c) / std::sqrt(a * a + b * b);
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION int is_point_to_left_of_line_segment(T xstart, T ystart, T xend, T yend, T xp, T yp)
	{
		static_assert(false, "Unsupported type for is_point_to_left_of_line_segment");
	}

    template <>
	CFUNC_DECORATION int is_point_to_left_of_line_segment(float xstart, float ystart, float xend, float yend, float xp, float yp)
	{
		return (xend - xstart) * (yp - ystart) - (yend - ystart) * (xp - xstart) > 0 ? 1 : 0;
	}

    template <>
	CFUNC_DECORATION int is_point_to_left_of_line_segment(double xstart, double ystart, double xend, double yend, double xp, double yp)
	{
		return (xend - xstart) * (yp - ystart) - (yend - ystart) * (xp - xstart) > 0 ? 1 : 0;
	}
    
    //##################################################################
	template <typename T>
	CFUNC_DECORATION T shoelace_term(T x1, T y1, T x2, T y2)
	{
		static_assert(false, "Unsupported type for shoelace_term");
	}

    template <>
	CFUNC_DECORATION float shoelace_term(float x1, float y1, float x2, float y2)
	{
		return maths_ops::cross_product<float>(x1, y1, x2, y2);
	}

    template <>
	CFUNC_DECORATION double shoelace_term(double x1, double y1, double x2, double y2)
	{
		return maths_ops::cross_product<double>(x1, y1, x2, y2);
	}

    //##################################################################
	template <typename T, typename RM>
	CFUNC_DECORATION void calculate_rotation_matrix(T angle_rad, RM rotation_matrix[4])
	{
        static_assert(false, "Unsupported type for calculate_rotation_matrix");
	}

    template <>
	CFUNC_DECORATION void calculate_rotation_matrix(float angle_rad, float rotation_matrix[4])
	{
		rotation_matrix[0] = std::cos(angle_rad);
		rotation_matrix[1] = -std::sin(angle_rad);
		rotation_matrix[2] = -rotation_matrix[1];
		rotation_matrix[3] = rotation_matrix[0];
	}

    template <>
	CFUNC_DECORATION void calculate_rotation_matrix(double angle_rad, double rotation_matrix[4])
	{
		rotation_matrix[0] = std::cos(angle_rad);
		rotation_matrix[1] = -std::sin(angle_rad);
		rotation_matrix[2] = -rotation_matrix[1];
		rotation_matrix[3] = rotation_matrix[0];
	}

    template <>
	CFUNC_DECORATION void calculate_rotation_matrix(float angle_rad, double rotation_matrix[4])
	{
		rotation_matrix[0] = std::cos(angle_rad);
		rotation_matrix[1] = -std::sin(angle_rad);
		rotation_matrix[2] = -rotation_matrix[1];
		rotation_matrix[3] = rotation_matrix[0];
	}

    template <>
	CFUNC_DECORATION void calculate_rotation_matrix(double angle_rad, float rotation_matrix[4])
	{
		rotation_matrix[0] = std::cos(angle_rad);
		rotation_matrix[1] = -std::sin(angle_rad);
		rotation_matrix[2] = -rotation_matrix[1];
		rotation_matrix[3] = rotation_matrix[0];
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION void translate_point(T* x, T* y, const T dx, const T dy)
	{
        static_assert(false, "Unsupported type for translate_point");
	}

    template <>
	CFUNC_DECORATION void translate_point(float* x, float* y, const float dx, const float dy)
	{
		*x += dx;
		*y += dy;
	}

    template <>
	CFUNC_DECORATION void translate_point(double* x, double* y, const double dx, const double dy)
	{
		*x += dx;
		*y += dy;
	}

    //##################################################################
	template <typename T, typename RM>
	CFUNC_DECORATION void rotate_point(T* x, T* y, const RM* rotation_matrix)
	{
		static_assert(false, "Unsupported type for rotate_point");
	}

    template <>
	CFUNC_DECORATION void rotate_point(float* x, float* y, const float rotation_matrix[4])
	{
		float xtemp = *x;
		float ytemp = *y;
		*x = rotation_matrix[0] * xtemp + rotation_matrix[1] * ytemp;
		*y = rotation_matrix[2] * xtemp + rotation_matrix[3] * ytemp;
	}

    template <>
	CFUNC_DECORATION void rotate_point(double* x, double* y, const double rotation_matrix[4])
	{
		double xtemp = *x;
		double ytemp = *y;
		*x = rotation_matrix[0] * xtemp + rotation_matrix[1] * ytemp;
		*y = rotation_matrix[2] * xtemp + rotation_matrix[3] * ytemp;
	}

    template <>
	CFUNC_DECORATION void rotate_point(float* x, float* y, const double rotation_matrix[4])
	{
		float xtemp = *x;
		float ytemp = *y;
		*x = rotation_matrix[0] * xtemp + rotation_matrix[1] * ytemp;
		*y = rotation_matrix[2] * xtemp + rotation_matrix[3] * ytemp;
	}

    template <>
	CFUNC_DECORATION void rotate_point(double* x, double* y, const float rotation_matrix[4])
	{
		double xtemp = *x;
		double ytemp = *y;
		*x = rotation_matrix[0] * xtemp + rotation_matrix[1] * ytemp;
		*y = rotation_matrix[2] * xtemp + rotation_matrix[3] * ytemp;
	}

    //##################################################################
	template <typename T, typename RM>
	CFUNC_DECORATION void rotate_point_about(T* x, T* y, const RM* rot_mat, const T xref, const T yref)
	{
		static_assert(false, "Unsupported type for rotate_point_about");
	}

    template <>
	CFUNC_DECORATION void rotate_point_about(float* x, float* y, const float rot_mat[4], const float xref, const float yref)
	{
		*x -= xref;
		*y -= yref;
		float xnew = *x * rot_mat[0] + *y * rot_mat[1];
		float ynew = *x * rot_mat[2] + *y * rot_mat[3];
		*x = xnew + xref;
		*y = ynew + yref;
	}

    template <>
	CFUNC_DECORATION void rotate_point_about(double* x, double* y, const double rot_mat[4], const double xref, const double yref)
	{
		*x -= xref;
		*y -= yref;
		double xnew = *x * rot_mat[0] + *y * rot_mat[1];
		double ynew = *x * rot_mat[2] + *y * rot_mat[3];
		*x = xnew + xref;
		*y = ynew + yref;
	}

    template <>
	CFUNC_DECORATION void rotate_point_about(float* x, float* y, const double rot_mat[4], const float xref, const float yref)
	{
		*x -= xref;
		*y -= yref;
		float xnew = *x * rot_mat[0] + *y * rot_mat[1];
		float ynew = *x * rot_mat[2] + *y * rot_mat[3];
		*x = xnew + xref;
		*y = ynew + yref;
	}

    template <>
	CFUNC_DECORATION void rotate_point_about(double* x, double* y, const float rot_mat[4], const double xref, const double yref)
	{
		*x -= xref;
		*y -= yref;
		double xnew = *x * rot_mat[0] + *y * rot_mat[1];
		double ynew = *x * rot_mat[2] + *y * rot_mat[3];
		*x = xnew + xref;
		*y = ynew + yref;
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION void unit_vector(T* unit_vec_x, T* unit_vec_y,
		const T original_vec_x, const T original_vec_y, const T zero_theshold = static_cast<T>(DEFAULT_ZERO_THRESHOLD))
	{
		static_assert(false, "Unsupported type for unit_vector");
	}

    template <>
	CFUNC_DECORATION void unit_vector(float* unit_vec_x, float* unit_vec_y,
		const float original_vec_x, const float original_vec_y, const float zero_theshold)
	{
		float mag = vector_magnitude<float>(original_vec_x, original_vec_y);
		if (mag > zero_theshold) {
			*unit_vec_x = original_vec_x / mag;
			*unit_vec_y = original_vec_y / mag;
		}
		else {
			*unit_vec_x = original_vec_x;
			*unit_vec_y = original_vec_y;
		}
	}

    template <>
	CFUNC_DECORATION void unit_vector(double* unit_vec_x, double* unit_vec_y,
		const double original_vec_x, const double original_vec_y, const double zero_theshold)
	{
		double mag = vector_magnitude<double>(original_vec_x, original_vec_y);
		if (mag > zero_theshold) {
			*unit_vec_x = original_vec_x / mag;
			*unit_vec_y = original_vec_y / mag;
		}
		else {
			*unit_vec_x = original_vec_x;
			*unit_vec_y = original_vec_y;
		}
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION T parallel_vector_component(
		const T vec_x, const T vec_y,
		const T ref_line_start_x, const T ref_line_start_y,
		const T ref_line_end_x, const T ref_line_end_y,
		const T zero_theshold = static_cast<T>(DEFAULT_ZERO_THRESHOLD))
	{
		static_assert(false, "Unsupported type for parallel_vector_component");
	}

    template <>
	CFUNC_DECORATION float parallel_vector_component(
		const float vec_x, const float vec_y,
		const float ref_line_start_x, const float ref_line_start_y,
		const float ref_line_end_x, const float ref_line_end_y,
		const float zero_theshold)
	{
		float line_vec_x = ref_line_end_x - ref_line_start_x;
		float line_vec_y = ref_line_end_y - ref_line_start_y;
		float unit_vec_x = 0.f, unit_vec_y = 0.f;
		unit_vector<float>(&unit_vec_x, &unit_vec_y, line_vec_x, line_vec_y, zero_theshold);

		return maths_ops::dot_product<float>(vec_x, vec_y, unit_vec_x, unit_vec_y);
	}

    template <>
	CFUNC_DECORATION double parallel_vector_component(
		const double vec_x, const double vec_y,
		const double ref_line_start_x, const double ref_line_start_y,
		const double ref_line_end_x, const double ref_line_end_y,
		const double zero_theshold)
	{
		double line_vec_x = ref_line_end_x - ref_line_start_x;
		double line_vec_y = ref_line_end_y - ref_line_start_y;
		double unit_vec_x = 0., unit_vec_y = 0.;
		unit_vector<double>(&unit_vec_x, &unit_vec_y, line_vec_x, line_vec_y, zero_theshold);

		return maths_ops::dot_product<double>(vec_x, vec_y, unit_vec_x, unit_vec_y);
	}

	//##################################################################
	template <typename T>
	CFUNC_DECORATION T perpendicular_vector_component(
		const T vec_x, const T vec_y,
		const T ref_line_start_x, const T ref_line_start_y,
		const T ref_line_end_x, const T ref_line_end_y,
		const T zero_theshold = static_cast<T>(DEFAULT_ZERO_THRESHOLD))
	{
		static_assert(false, "Unsupported type for perpendicular_vector_component");
	}

    template <>
	CFUNC_DECORATION float perpendicular_vector_component(
		const float vec_x, const float vec_y,
		const float ref_line_start_x, const float ref_line_start_y,
		const float ref_line_end_x, const float ref_line_end_y,
		const float zero_theshold)
	{
		float line_vec_x = ref_line_end_x - ref_line_start_x;
		float line_vec_y = ref_line_end_y - ref_line_start_y;
		float unit_vec_x = 0.f, unit_vec_y = 0.f;
		unit_vector<float>(&unit_vec_x, &unit_vec_y, line_vec_x, line_vec_y, zero_theshold);

		return maths_ops::cross_product<float>(unit_vec_x, unit_vec_y, vec_x, vec_y);
	}

    template <>
	CFUNC_DECORATION double perpendicular_vector_component(
		const double vec_x, const double vec_y,
		const double ref_line_start_x, const double ref_line_start_y,
		const double ref_line_end_x, const double ref_line_end_y,
		const double zero_theshold)
	{
		double line_vec_x = ref_line_end_x - ref_line_start_x;
		double line_vec_y = ref_line_end_y - ref_line_start_y;
		double unit_vec_x = 0., unit_vec_y = 0.;
		unit_vector<double>(&unit_vec_x, &unit_vec_y, line_vec_x, line_vec_y, zero_theshold);

		return maths_ops::cross_product<double>(unit_vec_x, unit_vec_y, vec_x, vec_y);
	}

    //##################################################################
	template <typename T>
	CFUNC_DECORATION void tangent_and_normal_vector_components(
        T* tang, T* norm, 
		const T vec_x, const T vec_y,
		const T ref_line_start_x, const T ref_line_start_y,
		const T ref_line_end_x, const T ref_line_end_y,
		const T zero_theshold = static_cast<T>(DEFAULT_ZERO_THRESHOLD))
	{
		static_assert(false, "Unsupported type for tangent_and_normal_vector_components");
	}

    template <>
	CFUNC_DECORATION void tangent_and_normal_vector_components(
        float* tang, float* norm, 
		const float vec_x, const float vec_y,
		const float ref_line_start_x, const float ref_line_start_y,
		const float ref_line_end_x, const float ref_line_end_y,
		const float zero_theshold)
	{
		float line_vec_x = ref_line_end_x - ref_line_start_x;
		float line_vec_y = ref_line_end_y - ref_line_start_y;
		float unit_vec_x = 0.f, unit_vec_y = 0.f;
		unit_vector<float>(&unit_vec_x, &unit_vec_y, line_vec_x, line_vec_y, zero_theshold);

        *tang = maths_ops::dot_product<float>(vec_x, vec_y, unit_vec_x, unit_vec_y);
        *norm = maths_ops::cross_product<float>(unit_vec_x, unit_vec_y, vec_x, vec_y);
	}

    template <>
	CFUNC_DECORATION void tangent_and_normal_vector_components(
        double* tang, double* norm, 
		const double vec_x, const double vec_y,
		const double ref_line_start_x, const double ref_line_start_y,
		const double ref_line_end_x, const double ref_line_end_y,
		const double zero_theshold)
	{
		double line_vec_x = ref_line_end_x - ref_line_start_x;
		double line_vec_y = ref_line_end_y - ref_line_start_y;
		double unit_vec_x = 0., unit_vec_y = 0.;
		unit_vector<double>(&unit_vec_x, &unit_vec_y, line_vec_x, line_vec_y, zero_theshold);

        *tang = maths_ops::dot_product<double>(vec_x, vec_y, unit_vec_x, unit_vec_y);
        *norm = maths_ops::cross_product<double>(unit_vec_x, unit_vec_y, vec_x, vec_y);
	}
}

namespace gis_ops
{
    //##################################################################
	template <typename T, typename GT>
	CFUNC_DECORATION void apply_affine_transformation(T* x, T* y, const int irow, const int icol, const GT geotransform[6])
	{
		static_assert(false, "Unsupported type for apply_affine_transformation");
	}

    template <>
	CFUNC_DECORATION void apply_affine_transformation(float* x, float* y, const int irow, const int icol, const float geotransform[6])
	{
		*x = geotransform[0] + icol * geotransform[1] + irow * geotransform[2];
		*y = geotransform[3] + icol * geotransform[4] + irow * geotransform[5];
	}

    template <>
	CFUNC_DECORATION void apply_affine_transformation(double* x, double* y, const int irow, const int icol, const double geotransform[6])
	{
		*x = geotransform[0] + icol * geotransform[1] + irow * geotransform[2];
		*y = geotransform[3] + icol * geotransform[4] + irow * geotransform[5];
	}

    template <>
	CFUNC_DECORATION void apply_affine_transformation(float* x, float* y, const int irow, const int icol, const double geotransform[6])
	{
		*x = geotransform[0] + icol * geotransform[1] + irow * geotransform[2];
		*y = geotransform[3] + icol * geotransform[4] + irow * geotransform[5];
	}

    template <>
	CFUNC_DECORATION void apply_affine_transformation(double* x, double* y, const int irow, const int icol, const float geotransform[6])
	{
		*x = geotransform[0] + icol * geotransform[1] + irow * geotransform[2];
		*y = geotransform[3] + icol * geotransform[4] + irow * geotransform[5];
	}

    //##################################################################
	template <typename T, typename GT>
	CFUNC_DECORATION void apply_inverse_affine_transformation(int* irow, int* icol, const T x, const T y, const GT transform_matrix[6])
	{
		static_assert(false, "Unsupported type for apply_inverse_affine_transformation");
	}

    template <>
	CFUNC_DECORATION void apply_inverse_affine_transformation(int* irow, int* icol, const float x, const float y, const float transform_matrix[6])
	{
		float Dx = x - transform_matrix[0];
		float Dy = y - transform_matrix[3];

		float det = transform_matrix[1] * transform_matrix[5] - transform_matrix[2] * transform_matrix[4];
		float det_icol = Dx * transform_matrix[5] - transform_matrix[2] * Dy;
		float det_irow = transform_matrix[1] * Dy - Dx * transform_matrix[4];

		*irow = static_cast<int>(std::round(det_irow / det));
		*icol = static_cast<int>(std::round(det_icol / det));
	}

    template <>
	CFUNC_DECORATION void apply_inverse_affine_transformation(int* irow, int* icol, const double x, const double y, const double transform_matrix[6])
	{
		double Dx = x - transform_matrix[0];
		double Dy = y - transform_matrix[3];

		double det = transform_matrix[1] * transform_matrix[5] - transform_matrix[2] * transform_matrix[4];
		double det_icol = Dx * transform_matrix[5] - transform_matrix[2] * Dy;
		double det_irow = transform_matrix[1] * Dy - Dx * transform_matrix[4];

		*irow = static_cast<int>(std::round(det_irow / det));
		*icol = static_cast<int>(std::round(det_icol / det));
	}

    template <>
	CFUNC_DECORATION void apply_inverse_affine_transformation(int* irow, int* icol, const float x, const float y, const double transform_matrix[6])
	{
		float Dx = x - transform_matrix[0];
		float Dy = y - transform_matrix[3];

		float det = transform_matrix[1] * transform_matrix[5] - transform_matrix[2] * transform_matrix[4];
		float det_icol = Dx * transform_matrix[5] - transform_matrix[2] * Dy;
		float det_irow = transform_matrix[1] * Dy - Dx * transform_matrix[4];

		*irow = static_cast<int>(std::round(det_irow / det));
		*icol = static_cast<int>(std::round(det_icol / det));
	}

    template <>
	CFUNC_DECORATION void apply_inverse_affine_transformation(int* irow, int* icol, const double x, const double y, const float transform_matrix[6])
	{
		double Dx = x - transform_matrix[0];
		double Dy = y - transform_matrix[3];

		double det = transform_matrix[1] * transform_matrix[5] - transform_matrix[2] * transform_matrix[4];
		double det_icol = Dx * transform_matrix[5] - transform_matrix[2] * Dy;
		double det_irow = transform_matrix[1] * Dy - Dx * transform_matrix[4];

		*irow = static_cast<int>(std::round(det_irow / det));
		*icol = static_cast<int>(std::round(det_icol / det));
	}

    //##################################################################
	template <typename T, typename GT>
	CFUNC_DECORATION void set_affine_transformation(GT geotransform[6], 
        const T x_top_left, const T y_top_left, const T dx, const T dy, const T angle_deg)
	{
		static_assert(false, "Unsupported type for set_affine_transformation");
	}

    template <>
	CFUNC_DECORATION void set_affine_transformation(float geotransform[6], 
        const float x_top_left, const float y_top_left, const float dx, const float dy, const float angle_deg)
	{
		float angle_rad = maths_ops::deg_to_rad<float>(angle_deg);
		geotransform[0] = x_top_left;
		geotransform[1] = dx * std::cos(angle_rad);
		geotransform[2] = dy * std::sin(angle_rad);
		geotransform[3] = y_top_left;
		geotransform[4] = dx * std::sin(angle_rad);
		geotransform[5] = -dy * std::cos(angle_rad);
	}

    template <>
	CFUNC_DECORATION void set_affine_transformation(double geotransform[6], 
        const double x_top_left, const double y_top_left, const double dx, const double dy, const double angle_deg)
	{
		double angle_rad = maths_ops::deg_to_rad<double>(angle_deg);
		geotransform[0] = x_top_left;
		geotransform[1] = dx * std::cos(angle_rad);
		geotransform[2] = dy * std::sin(angle_rad);
		geotransform[3] = y_top_left;
		geotransform[4] = dx * std::sin(angle_rad);
		geotransform[5] = -dy * std::cos(angle_rad);
	}

    template <>
	CFUNC_DECORATION void set_affine_transformation(float geotransform[6], 
        const double x_top_left, const double y_top_left, const double dx, const double dy, const double angle_deg)
	{
		double angle_rad = maths_ops::deg_to_rad<double>(angle_deg);
		geotransform[0] = x_top_left;
		geotransform[1] = dx * std::cos(angle_rad);
		geotransform[2] = dy * std::sin(angle_rad);
		geotransform[3] = y_top_left;
		geotransform[4] = dx * std::sin(angle_rad);
		geotransform[5] = -dy * std::cos(angle_rad);
	}

    template <>
	CFUNC_DECORATION void set_affine_transformation(double geotransform[6], 
        const float x_top_left, const float y_top_left, const float dx, const float dy, const float angle_deg)
	{
		float angle_rad = maths_ops::deg_to_rad<float>(angle_deg);
		geotransform[0] = x_top_left;
		geotransform[1] = dx * std::cos(angle_rad);
		geotransform[2] = dy * std::sin(angle_rad);
		geotransform[3] = y_top_left;
		geotransform[4] = dx * std::sin(angle_rad);
		geotransform[5] = -dy * std::cos(angle_rad);
	}

	//##################################################################
	template <typename T, typename GT>
	CFUNC_DECORATION void decrypt_affine_transformation_matrix(T* top_left_x, T* top_left_y,
		T* dx, T* dy, T* angle_rad, const GT transform_matrix[6])
	{
		static_assert(false, "Unsupported type for decrypt_affine_transformation_matrix");
	}

    template <>
	CFUNC_DECORATION void decrypt_affine_transformation_matrix(float* top_left_x, float* top_left_y,
		float* dx, float* dy, float* angle_rad, const float transform_matrix[6])
	{
		// Get the origin
		*top_left_x = transform_matrix[0];
		*top_left_y = transform_matrix[3];

		// Get the angle as well as the resolution
		*angle_rad = std::atan(transform_matrix[4] / transform_matrix[1]);
		if (std::abs(*angle_rad) > 1.e-6) {
			*dx = transform_matrix[1] / std::cos(*angle_rad);
			*dy = -transform_matrix[5] / std::cos(*angle_rad);
		}
		else {
			*dx = transform_matrix[1];
			*dy = -transform_matrix[5];
		}
	}

    template <>
	CFUNC_DECORATION void decrypt_affine_transformation_matrix(double* top_left_x, double* top_left_y,
		double* dx, double* dy, double* angle_rad, const double transform_matrix[6])
	{
		// Get the origin
		*top_left_x = transform_matrix[0];
		*top_left_y = transform_matrix[3];

		// Get the angle as well as the resolution
		*angle_rad = std::atan(transform_matrix[4] / transform_matrix[1]);
		if (std::abs(*angle_rad) > 1.e-6) {
			*dx = transform_matrix[1] / std::cos(*angle_rad);
			*dy = -transform_matrix[5] / std::cos(*angle_rad);
		}
		else {
			*dx = transform_matrix[1];
			*dy = -transform_matrix[5];
		}
	}

    template <>
	CFUNC_DECORATION void decrypt_affine_transformation_matrix(float* top_left_x, float* top_left_y,
		float* dx, float* dy, float* angle_rad, const double transform_matrix[6])
	{
		// Get the origin
		*top_left_x = transform_matrix[0];
		*top_left_y = transform_matrix[3];

		// Get the angle as well as the resolution
		*angle_rad = std::atan(transform_matrix[4] / transform_matrix[1]);
		if (std::abs(*angle_rad) > 1.e-6) {
			*dx = transform_matrix[1] / std::cos(*angle_rad);
			*dy = -transform_matrix[5] / std::cos(*angle_rad);
		}
		else {
			*dx = transform_matrix[1];
			*dy = -transform_matrix[5];
		}
	}

    template <>
	CFUNC_DECORATION void decrypt_affine_transformation_matrix(double* top_left_x, double* top_left_y,
		double* dx, double* dy, double* angle_rad, const float transform_matrix[6])
	{
		// Get the origin
		*top_left_x = transform_matrix[0];
		*top_left_y = transform_matrix[3];

		// Get the angle as well as the resolution
		*angle_rad = std::atan(transform_matrix[4] / transform_matrix[1]);
		if (std::abs(*angle_rad) > 1.e-6) {
			*dx = transform_matrix[1] / std::cos(*angle_rad);
			*dy = -transform_matrix[5] / std::cos(*angle_rad);
		}
		else {
			*dx = transform_matrix[1];
			*dy = -transform_matrix[5];
		}
	}
}

#endif // COMMON_MATHS_AND_GEOMETRY_OPS_H