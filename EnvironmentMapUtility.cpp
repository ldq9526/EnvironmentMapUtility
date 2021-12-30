#include "EnvironmentMapUtility.h"

#include <cmath>

namespace EnvironmentMap
{
	const std::string EnvironmentMapUtility::cubeMapNames[6] = { "px", "py", "pz", "nx", "ny", "nz" };

	const cv::Vec3f EnvironmentMapUtility::dirBasis[6][3] = {
		{cv::Vec3f(0,0,-1), cv::Vec3f(0,-1,0), cv::Vec3f(1,0,0)},// +x
		{cv::Vec3f(1,0,0), cv::Vec3f(0,0,1), cv::Vec3f(0,1,0)},// +y
		{cv::Vec3f(1,0,0), cv::Vec3f(0,-1,0), cv::Vec3f(0,0,1)},// +z
		{cv::Vec3f(0,0,1), cv::Vec3f(0,-1,0), cv::Vec3f(-1,0,0)},// -x
		{cv::Vec3f(1,0,0), cv::Vec3f(0,0,-1), cv::Vec3f(0,-1,0)},// -y
		{cv::Vec3f(-1,0,0), cv::Vec3f(0,-1,0), cv::Vec3f(0,0,-1)},// -z
	};

	Interpolation EnvironmentMapUtility::interpolation = Interpolation::Bilinear;

	float EnvironmentMapUtility::clamp(float v, float minv, float maxv)
	{
		if (minv > maxv)
			std::swap(minv, maxv);
		if (v < minv)
			return minv;
		if (v > maxv)
			return maxv;
		return v;
	}

	void EnvironmentMapUtility::getIntegerInterval(float x, int& x0, int& x1)
	{
		int px = int(std::floor(x));
		float dx = x - px - 0.5f;
		if (dx < 0.f)
		{
			x0 = px - 1;
			x1 = px;
		}
		else if (dx > 0.f)
		{
			x0 = px;
			x1 = px + 1;
		}
		else
		{
			x0 = px;
			x1 = px;
		}
	}

	float EnvironmentMapUtility::weightBicubic(float x, float a)
	{
		x = std::abs(x);
		float x2 = x * x;
		float x3 = x2 * x;
		if (x < 1.f)
			return (a + 2) * x3 - (a + 3) * x2 + 1;
		else if (x < 2.f)
			return a * x3 - 5 * a * x2 + 8 * a * x - 4 * a;
		else
			return 0.f;
	}

	float EnvironmentMapUtility::weightLanczos(float x, float a)
	{
		if (std::fabs(x) < 1e-4f)
			return 1.f;

		float xpi = PI * x;
		return a * std::sin(xpi) * std::sin(xpi / a) / (xpi * xpi);
	}

	cv::Vec3f EnvironmentMapUtility::samplePixel(const cv::Mat& image, int row, int col)
	{
		if (col < 0)
			col += image.cols;
		if (col >= image.cols)
			col %= image.cols;

		if (row < 0)
		{
			row = -1 - row;
			col = (col + image.cols / 2) % image.cols;
		}
		if (row >= image.rows)
		{
			row = 2 * image.rows - row - 1;
			col = (col + image.cols / 2) % image.cols;
		}

		return image.at<cv::Vec3f>(row, col);
	}

	cv::Vec3f EnvironmentMapUtility::sampleNearest(const cv::Mat& image, float row, float col)
	{
		return image.at<cv::Vec3f>(int(row), int(col));
	}

	cv::Vec3f EnvironmentMapUtility::sampleBilinear(const cv::Mat& image, float row, float col)
	{
		int x0, x1;//双线性插值所用4点的行范围
		int y0, y1;//双线性插值所用4点的列范围
		getIntegerInterval(col, x0, x1);
		getIntegerInterval(row, y0, y1);

		float wx = x1 + 0.5f - col;
		float wy = y1 + 0.5f - row;
		return
			wx * wy * samplePixel(image, y0, x0) +
			(1.f - wx) * wy * samplePixel(image, y0, x1) +
			wx * (1.f - wy) * samplePixel(image, y1, x0) + 
			(1.f - wx) * (1.f - wy) * samplePixel(image, y1, x1);
	}

	cv::Vec3f EnvironmentMapUtility::sampleBicubic(const cv::Mat& image, float row, float col)
	{
		int x0 = int(std::floor(col)), y0 = int(std::floor(row));

		cv::Vec3f result(0, 0, 0);
		for (int dx = -1; dx <= 2; dx++)
		{
			int x = x0 + dx;
			float wx = weightBicubic(x - col);
			for (int dy = -1; dy <= 2; dy++)
			{
				int y = y0 + dy;
				float wy = weightBicubic(y - row);
				result += (wx * wy * samplePixel(image, y, x));
			}
		}
		return result;
	}

	cv::Vec3f EnvironmentMapUtility::sampleLanczos(const cv::Mat& image, float row, float col)
	{
		int x0 = int(std::floor(col)), y0 = int(std::floor(row));

		cv::Vec3f result(0, 0, 0);
		for (int dx = -4; dx <= 5; dx++)
		{
			int x = x0 + dx;
			float wx = weightLanczos(x - col, 5);
			for (int dy = -4; dy <= 5; dy++)
			{
				int y = y0 + dy;
				float wy = weightLanczos(y - row, 5);
				result += (wx * wy * samplePixel(image, y, x));
			}
		}
		return result;
	}

	cv::Vec3f EnvironmentMapUtility::sampleUV(const cv::Mat& image, float u, float v)
	{
		float col = u * image.cols;
		float row = v * image.rows;

		switch (interpolation)
		{
		case Interpolation::Nearest:
			return sampleNearest(image, row, col);
		case Interpolation::Bicubic:
			return sampleBicubic(image, row, col);
		case Interpolation::Lanczos:
			return sampleLanczos(image, row, col);
		default:
			return sampleBilinear(image, row, col);
		}
	}

	cv::Vec3f EnvironmentMapUtility::sampleDirection(const cv::Mat& image, float phi, float theta)
	{
		float u = clamp(phi / PI2, 0.f, 1.f);
		float v = clamp(theta / PI, 0.f, 1.f);
		return sampleUV(image, u, v);
	}

	bool EnvironmentMapUtility::isValidEnvironmentMap(const cv::Mat& image)
	{
		return (!image.empty()) && (image.cols == image.rows * 2);
	}

	std::map<std::string, cv::Mat> EnvironmentMapUtility::convertToCubeMap(const cv::Mat& image, Interpolation interpolation)
	{
		EnvironmentMapUtility::interpolation = interpolation;

		std::map<std::string, cv::Mat> result;
		if (!isValidEnvironmentMap(image))
			return result;

		cv::Mat environmentMap;
		image.copyTo(environmentMap);
		int imageType = environmentMap.type();
		if (imageType == CV_8UC3)
			environmentMap.convertTo(environmentMap, CV_32FC3, 1.0 / 255);

		int size = environmentMap.rows / 2;
		for (int i = 0; i < 6; i++)
		{
			cv::Mat face = cv::Mat_<cv::Vec3f>(size, size);
			for (int row = 0; row < size; row++)
			{
				for (int col = 0; col < size; col++)
				{
					float u = (2 * col + 1.f) / size - 1.f;
					float v = (2 * row + 1.f) / size - 1.f;
					cv::Vec3f dir = cv::normalize(u * dirBasis[i][0] + v * dirBasis[i][1] + dirBasis[i][2]);
					float theta = std::acosf(dir[1]);
					float phi = std::atan2(-dir[0], -dir[2]);
					if (phi < 0.f)
						phi += PI2;
					face.at<cv::Vec3f>(row, col) = sampleDirection(environmentMap, phi, theta);
				}
			}

			if (imageType == CV_8UC3)
				face.convertTo(face, CV_8UC3, 255);
			result[cubeMapNames[i]] = face;
		}

		return result;
	}
}