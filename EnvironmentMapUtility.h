#ifndef ENVIRONMENT_MAP_UTILITY_H
#define ENVIRONMENT_MAP_UTILITY_H

#include <opencv2/opencv.hpp>

#include <array>
#include <functional>
#include <map>
#include <memory>
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
		Lanczos = 3,
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

		/// <summary>
		/// ������ֵ�㷨
		/// </summary>
		static Interpolation interpolation;

		/// <summary>
		/// ǰ3����г������
		/// </summary>
		static std::function<double(const cv::Vec3d&)> SHB[25];

	private:
		/// <summary>
		/// ��v�ض������� [minv, maxv]
		/// </summary>
		template <typename T>
		static T clamp(T v, T minv, T maxv);

		/// <summary>
		/// ������������x���ڵ���С�������� [x0, x1]
		/// </summary>
		static void getIntegerInterval(float x, int& x0, int& x1);

		/// <summary>
		/// ˫���β�ֵȨֵ
		/// </summary>
		static float weightBicubic(float x, float a = -0.5f);

		/// <summary>
		/// Lanczos��ֵȨֵ
		/// </summary>
		static float weightLanczos(float x, float a = 2.f);

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
		/// Lanczos��ֵ��������������
		/// </summary>
		static cv::Vec3f sampleLanczos(const cv::Mat& image, float row, float col);

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
		static bool isValidEnvironmentMap(const cv::Mat& image);

		/// <summary>
		/// ��������ͼת��Ϊ6�����cube map
		/// </summary>
		static std::map<std::string, cv::Mat> convertToCubeMap(const cv::Mat& image, Interpolation interpolation = Interpolation::Bilinear);

		/// <summary>
		/// ��ȡ diffuse ���ն���ͼ, pixel(row, col) = (1/pi) * (4*pi/N) * sum[ I * max(dot(n, d), 0.0) ]
		/// </summary>
		static cv::Mat getDiffuseIrradianceMap(const cv::Mat& image, int sampleCount = 10000, Interpolation interpolation = Interpolation::Nearest);

		/// <summary>
		/// ͨ������ÿ�����ص�����ǻ�ȡ��гϵ��
		/// </summary>
		static std::unique_ptr<std::vector<cv::Vec3d> > getCoefficientsSH(const cv::Mat& image, int order);

		/// <summary>
		/// ͨ��������Ȳ�����ȡ��гϵ��
		/// </summary>
		static std::unique_ptr<std::vector<cv::Vec3d> > uniformSampleCoefficientsSH(const cv::Mat& image, int order, int sampleCount = 10000, Interpolation interpolation = Interpolation::Nearest);

		/// <summary>
		/// ����гϵ�����ӻ�������ͼ
		/// </summary>
		static cv::Mat getEnvironmentMapFromSH(const std::unique_ptr<std::vector<cv::Vec3d> >& coeffSH, int order);
	};
}

#endif