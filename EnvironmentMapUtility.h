#ifndef ENVIRONMENT_MAP_UTILITY_H
#define ENVIRONMENT_MAP_UTILITY_H

#include <opencv2/opencv.hpp>

#include <map>
#include <string>

namespace EnvironmentMap
{
	const float PI = 3.14159265f;
	const float PI2 = 6.2831853f;
	const float PI4 = 12.5663706f;

	enum Interpolation
	{
		Nearest = 0,
		Bilinear = 1,
		Bicubic = 2,
	};

	class EnvironmentMapUtility
	{
	private:
		/// <summary>
		/// cube map 6个面的名称
		/// </summary>
		const static std::string cubeMapNames[6];

		/// <summary>
		/// 相机面朝cube map 6个面时方向的基向量
		/// </summary>
		const static cv::Vec3f dirBasis[6][3];

	private:
		/// <summary>
		/// 把v截断在区间 [minv, maxv]
		/// </summary>
		/// <param name="v"></param>
		/// <param name="minv"></param>
		/// <param name="maxv"></param>
		/// <returns></returns>
		static float clamp(float v, float minv, float maxv);

		/// <summary>
		/// 返回正浮点数x所在的最小整数区间 [x0, x1]
		/// </summary>
		/// <param name="x"></param>
		/// <param name="x0"></param>
		/// <param name="x1"></param>
		static void getIntegerInterval(float x, int& x0, int& x1);

		/// <summary>
		/// 双三次插值权值
		/// </summary>
		static float weightBicubic(float x, float a = -0.5f);

		/// <summary>
		/// 从环境贴图上读取像素 (row, col), 可能为负值
		/// </summary>
		static cv::Vec3f samplePixel(const cv::Mat& image, int row, int col);

		/// <summary>
		/// 最邻近插值采样非整数坐标 (row, col)
		/// </summary>
		static cv::Vec3f sampleNearest(const cv::Mat& image, float row, float col);

		/// <summary>
		/// 双线性插值采样非整数坐标 (row, col)
		/// </summary>
		static cv::Vec3f sampleBilinear(const cv::Mat& image, float row, float col);

		/// <summary>
		/// 双三次插值采样非整数坐标
		/// </summary>
		static cv::Vec3f sampleBicubic(const cv::Mat& image, float row, float col);

		/// <summary>
		/// 从环境贴图插值采样纹理坐标(col, row) = (u, v)
		/// </summary>
		static cv::Vec3f sampleUV(const cv::Mat& image, float u, float v);

		/// <summary>
		/// 从环境贴图采样方向(phi, theta): phi 取值 [0, 2pi], theta 取值 [0, pi] 
		/// </summary>
		static cv::Vec3f sampleDirection(const cv::Mat& image, float phi, float theta);

	public:
		/// <summary>
		/// 判断图像是否为有效的环境贴图
		/// </summary>
		/// <param name="image"></param>
		/// <returns></returns>
		static bool isValidEnvironmentMap(const cv::Mat& image);

		/// <summary>
		/// 将环境贴图转换为6个面的cube map
		/// </summary>
		/// <param name="image"></param>
		/// <returns></returns>
		static std::map<std::string, cv::Mat> convertToCubeMap(const cv::Mat& image);
	};
}

#endif