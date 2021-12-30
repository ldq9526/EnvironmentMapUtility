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
		/// cube map 6���������
		/// </summary>
		const static std::string cubeMapNames[6];

		/// <summary>
		/// ����泯cube map 6����ʱ����Ļ�����
		/// </summary>
		const static cv::Vec3f dirBasis[6][3];

	private:
		/// <summary>
		/// ��v�ض������� [minv, maxv]
		/// </summary>
		/// <param name="v"></param>
		/// <param name="minv"></param>
		/// <param name="maxv"></param>
		/// <returns></returns>
		static float clamp(float v, float minv, float maxv);

		/// <summary>
		/// ������������x���ڵ���С�������� [x0, x1]
		/// </summary>
		/// <param name="x"></param>
		/// <param name="x0"></param>
		/// <param name="x1"></param>
		static void getIntegerInterval(float x, int& x0, int& x1);

		/// <summary>
		/// ˫���β�ֵȨֵ
		/// </summary>
		static float weightBicubic(float x, float a = -0.5f);

		/// <summary>
		/// �ӻ�����ͼ�϶�ȡ���� (row, col), ����Ϊ��ֵ
		/// </summary>
		static cv::Vec3f samplePixel(const cv::Mat& image, int row, int col);

		/// <summary>
		/// ���ڽ���ֵ�������������� (row, col)
		/// </summary>
		static cv::Vec3f sampleNearest(const cv::Mat& image, float row, float col);

		/// <summary>
		/// ˫���Բ�ֵ�������������� (row, col)
		/// </summary>
		static cv::Vec3f sampleBilinear(const cv::Mat& image, float row, float col);

		/// <summary>
		/// ˫���β�ֵ��������������
		/// </summary>
		static cv::Vec3f sampleBicubic(const cv::Mat& image, float row, float col);

		/// <summary>
		/// �ӻ�����ͼ��ֵ������������(col, row) = (u, v)
		/// </summary>
		static cv::Vec3f sampleUV(const cv::Mat& image, float u, float v);

		/// <summary>
		/// �ӻ�����ͼ��������(phi, theta): phi ȡֵ [0, 2pi], theta ȡֵ [0, pi] 
		/// </summary>
		static cv::Vec3f sampleDirection(const cv::Mat& image, float phi, float theta);

	public:
		/// <summary>
		/// �ж�ͼ���Ƿ�Ϊ��Ч�Ļ�����ͼ
		/// </summary>
		/// <param name="image"></param>
		/// <returns></returns>
		static bool isValidEnvironmentMap(const cv::Mat& image);

		/// <summary>
		/// ��������ͼת��Ϊ6�����cube map
		/// </summary>
		/// <param name="image"></param>
		/// <returns></returns>
		static std::map<std::string, cv::Mat> convertToCubeMap(const cv::Mat& image);
	};
}

#endif