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

	double basisSH00(const cv::Vec3d& d) {
		// 0.5 * sqrt(1/pi)
		return 0.282095;
	}

	double basisSH1n1(const cv::Vec3d& d) {
		// -sqrt(3/(4pi)) * y
		return -0.488603 * d[1];
	}

	double basisSH10(const cv::Vec3d& d) {
		// sqrt(3/(4pi)) * z
		return 0.488603 * d[2];
	}

	double basisSH1p1(const cv::Vec3d& d) {
		// -sqrt(3/(4pi)) * x
		return -0.488603 * d[0];
	}

	double basisSH2n2(const cv::Vec3d& d) {
		// 0.5 * sqrt(15/pi) * x * y
		return 1.092548 * d[0] * d[1];
	}

	double basisSH2n1(const cv::Vec3d& d) {
		// -0.5 * sqrt(15/pi) * y * z
		return -1.092548 * d[1] * d[2];
	}

	double basisSH20(const cv::Vec3d& d) {
		// 0.25 * sqrt(5/pi) * (-x^2-y^2+2z^2)
		return 0.315392 * (-d[0] * d[0] - d[1] * d[1] + 2.0 * d[2] * d[2]);
	}

	double basisSH2p1(const cv::Vec3d& d) {
		// -0.5 * sqrt(15/pi) * x * z
		return -1.092548 * d[0] * d[2];
	}

	double basisSH2p2(const cv::Vec3d& d) {
		// 0.25 * sqrt(15/pi) * (x^2 - y^2)
		return 0.546274 * (d[0] * d[0] - d[1] * d[1]);
	}

	double basisSH3n3(const cv::Vec3d& d) {
		// -0.25 * sqrt(35/(2pi)) * y * (3x^2 - y^2)
		return -0.590044 * d[1] * (3.0 * d[0] * d[0] - d[1] * d[1]);
	}

	double basisSH3n2(const cv::Vec3d& d) {
		// 0.5 * sqrt(105/pi) * x * y * z
		return 2.890611 * d[0] * d[1] * d[2];
	}

	double basisSH3n1(const cv::Vec3d& d) {
		// -0.25 * sqrt(21/(2pi)) * y * (4z^2-x^2-y^2)
		return -0.457046 * d[1] * (4.0 * d[2] * d[2] - d[0] * d[0]
			- d[1] * d[1]);
	}

	double basisSH30(const cv::Vec3d& d) {
		// 0.25 * sqrt(7/pi) * z * (2z^2 - 3x^2 - 3y^2)
		return 0.373176 * d[2] * (2.0 * d[2] * d[2] - 3.0 * d[0] * d[0]
			- 3.0 * d[1] * d[1]);
	}

	double basisSH3p1(const cv::Vec3d& d) {
		// -0.25 * sqrt(21/(2pi)) * x * (4z^2-x^2-y^2)
		return -0.457046 * d[0] * (4.0 * d[2] * d[2] - d[0] * d[0]
			- d[1] * d[1]);
	}

	double basisSH3p2(const cv::Vec3d& d) {
		// 0.25 * sqrt(105/pi) * z * (x^2 - y^2)
		return 1.445306 * d[2] * (d[0] * d[0] - d[1] * d[1]);
	}

	double basisSH3p3(const cv::Vec3d& d) {
		// -0.25 * sqrt(35/(2pi)) * x * (x^2-3y^2)
		return -0.590044 * d[0] * (d[0] * d[0] - 3.0 * d[1] * d[1]);
	}

	double basisSH4n4(const cv::Vec3d& d) {
		// 0.75 * sqrt(35/pi) * x * y * (x^2-y^2)
		return 2.503343 * d[0] * d[1] * (d[0] * d[0] - d[1] * d[1]);
	}

	double basisSH4n3(const cv::Vec3d& d) {
		// -0.75 * sqrt(35/(2pi)) * y * z * (3x^2-y^2)
		return -1.770131 * d[1] * d[2] * (3.0 * d[0] * d[0] - d[1] * d[1]);
	}

	double basisSH4n2(const cv::Vec3d& d) {
		// 0.75 * sqrt(5/pi) * x * y * (7z^2-1)
		return 0.946175 * d[0] * d[1] * (7.0 * d[2] * d[2] - 1.0);
	}

	double basisSH4n1(const cv::Vec3d& d) {
		// -0.75 * sqrt(5/(2pi)) * y * z * (7z^2-3)
		return -0.669047 * d[1] * d[2] * (7.0 * d[2] * d[2] - 3.0);
	}

	double basisSH40(const cv::Vec3d& d) {
		// 3/16 * sqrt(1/pi) * (35z^4-30z^2+3)
		double z2 = d[2] * d[2];
		return 0.105786 * (35.0 * z2 * z2 - 30.0 * z2 + 3.0);
	}

	double basisSH4p1(const cv::Vec3d& d) {
		// -0.75 * sqrt(5/(2pi)) * x * z * (7z^2-3)
		return -0.669047 * d[0] * d[2] * (7.0 * d[2] * d[2] - 3.0);
	}

	double basisSH4p2(const cv::Vec3d& d) {
		// 3/8 * sqrt(5/pi) * (x^2 - y^2) * (7z^2 - 1)
		return 0.473087 * (d[0] * d[0] - d[1] * d[1])
			* (7.0 * d[2] * d[2] - 1.0);
	}

	double basisSH4p3(const cv::Vec3d& d) {
		// -0.75 * sqrt(35/(2pi)) * x * z * (x^2 - 3y^2)
		return -1.770131 * d[0] * d[2] * (d[0] * d[0] - 3.0 * d[1] * d[1]);
	}

	double basisSH4p4(const cv::Vec3d& d) {
		// 3/16*sqrt(35/pi) * (x^2 * (x^2 - 3y^2) - y^2 * (3x^2 - y^2))
		double x2 = d[0] * d[0];
		double y2 = d[1] * d[1];
		return 0.625836 * (x2 * (x2 - 3.0 * y2) - y2 * (3.0 * x2 - y2));
	}

	std::function<double(const cv::Vec3d&)> EnvironmentMapUtility::SHB[25] = {
		basisSH00,
		basisSH1n1,basisSH10,basisSH1p1,
		basisSH2n2,basisSH2n1,basisSH20,basisSH2p1,basisSH2p2,
		basisSH3n3,basisSH3n2,basisSH3n1,basisSH30,basisSH3p1,basisSH3p2,basisSH3p3,
		basisSH4n4,basisSH4n3,basisSH4n2,basisSH4n1,basisSH40,basisSH4p1,basisSH4p2,basisSH4p3,basisSH4p4
	};

	template <typename T>
	T EnvironmentMapUtility::clamp(T v, T minv, T maxv)
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
		return image.at<cv::Vec3f>(clamp(int(row), 0, image.rows - 1), clamp(int(col), 0, image.cols - 1));
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

	cv::Mat EnvironmentMapUtility::getDiffuseIrradianceMap(const cv::Mat& image, int sampleCount, Interpolation interpolation)
	{
		EnvironmentMapUtility::interpolation = interpolation;

		cv::Mat result = cv::Mat_<cv::Vec3f>(128, 256);
		if (!isValidEnvironmentMap(image))
			return result;

		cv::Mat environmentMap;
		image.copyTo(environmentMap);
		int imageType = environmentMap.type();
		if (imageType == CV_8UC3)
		{
			environmentMap.convertTo(environmentMap, CV_32FC3, 1.0 / 255);
			cv::pow(environmentMap, 2.2, environmentMap);
		}

		std::vector<float> samplesPhi(sampleCount, 0.f), samplesTheta(sampleCount, 0.f);
		cv::RNG rng(0);
		for (int i = 0; i < sampleCount; i++)
		{
			float u = rng.uniform(0.f, 1.f), v = rng.uniform(0.f, 1.f);
			samplesPhi[i] = PI2 * u;
			samplesTheta[i] = std::acos(2 * v - 1.f);
		}

		for (int x = 0; x < result.cols; x++)
		{
			float phi = PI2 * (x + 0.5f) / result.cols;
			for (int y = 0; y < result.rows; y++)
			{
				float theta = PI * (y + 0.5f) / result.rows;
				cv::Vec3f normal(-std::sin(theta) * std::sin(phi), std::cos(theta), std::sin(theta) * std::cos(phi));

				cv::Vec3f value(0, 0, 0);
				for (int i = 0; i < sampleCount; i++)
				{
					cv::Vec3f dir(-std::sin(samplesTheta[i]) * std::sin(samplesPhi[i]), std::cos(samplesTheta[i]), std::sin(samplesTheta[i]) * std::cos(samplesPhi[i]));
					float H = normal.dot(dir);
					cv::Vec3f sample(0, 0, 0);
					if (H > 0.f)
					{
						sample = sampleDirection(environmentMap, samplesPhi[i], samplesTheta[i]);
					}
					else if (H < 0.f)
					{
						H = 0.f - H;
						if(samplesPhi[i] < PI)
							sample = sampleDirection(environmentMap, samplesPhi[i] + PI, PI - samplesTheta[i]);
						else if(samplesPhi[i] < PI)
							sample = sampleDirection(environmentMap, samplesPhi[i] - PI, PI - samplesTheta[i]);
						else
							sample = sampleDirection(environmentMap, 0.f, PI - samplesTheta[i]);
					}
					value += (H * sample);
				}
				result.at<cv::Vec3f>(y, x) = value;
			}
		}

		result *= (2.f / sampleCount);
		if (imageType == CV_8UC3)
		{
			cv::pow(result, 1.0 / 2.2, result);
			result.convertTo(result, CV_8UC3, 255.0);
		}

		return result;
	}


	std::unique_ptr<std::vector<cv::Vec3d> > EnvironmentMapUtility::getCoefficientsSH(const cv::Mat& image, int order)
	{
		if (order < 0 || order > 4)
			order = 2;

		int coeffCount = (order + 1) * (order + 1);
		cv::Vec3d zero(0, 0, 0);
		std::unique_ptr<std::vector<cv::Vec3d> > coeffSH(new std::vector<cv::Vec3d>(coeffCount, zero));

		if (!isValidEnvironmentMap(image))
			return coeffSH;

		cv::Mat environmentMap;
		image.copyTo(environmentMap);
		int imageType = environmentMap.type();
		if (imageType == CV_8UC3)
			environmentMap.convertTo(environmentMap, CV_32FC3, 1.0 / 255);

		float pixelArea = (PI * PI2 / (image.cols * image.rows));
		for (int y = 0; y < image.rows; y++)
		{
			float theta = PI * (y + 0.5f) / image.rows;
			double r = std::sin(theta);
			double weight = r * pixelArea;
			for (int x = 0; x < image.cols; x++)
			{
				float phi = PI2 * (x + 0.5f) / image.cols;
				cv::Vec3f color = environmentMap.at<cv::Vec3f>(y, x);

				cv::Vec3d dir(-r * std::sin(phi), std::cos(theta), r * cos(phi));
				for (int k = 0; k < coeffCount; k++)
					(*coeffSH)[k] += (weight * SHB[k](dir) * color);
			}
		}

		return coeffSH;
	}

	std::unique_ptr<std::vector<cv::Vec3d> > EnvironmentMapUtility::uniformSampleCoefficientsSH(const cv::Mat& image, int order, int sampleCount, Interpolation interpolation)
	{
		EnvironmentMapUtility::interpolation = interpolation;

		if (order < 0 || order > 4)
			order = 2;

		int coeffCount = (order + 1) * (order + 1);
		cv::Vec3d zero(0, 0, 0);
		std::unique_ptr<std::vector<cv::Vec3d> > coeffSH(new std::vector<cv::Vec3d>(coeffCount, zero));

		if (!isValidEnvironmentMap(image))
			return coeffSH;

		cv::Mat environmentMap;
		image.copyTo(environmentMap);
		int imageType = environmentMap.type();
		if (imageType == CV_8UC3)
			environmentMap.convertTo(environmentMap, CV_32FC3, 1.0 / 255);

		cv::RNG rng(0);
		for (int i = 0; i < sampleCount; i++)
		{
			float phi = PI2 * rng.uniform(0.f, 1.f);
			float theta = std::acos(2 * rng.uniform(0.f, 1.f) - 1);
			double r = std::sin(theta);
			cv::Vec3d dir(-r * std::sin(phi), std::cos(theta), r * cos(phi));

			for (int k = 0; k < coeffCount; k++)
				(*coeffSH)[k] += (SHB[k](dir) * sampleDirection(environmentMap, phi, theta));
		}

		double weight = PI4 / sampleCount;
		for (int k = 0; k < coeffCount; k++)
			(*coeffSH)[k] *= weight;

		return coeffSH;
	}

	cv::Mat EnvironmentMapUtility::getEnvironmentMapFromSH(const std::unique_ptr<std::vector<cv::Vec3d> >& coeffSH, int order)
	{
		cv::Mat result = cv::Mat_<cv::Vec3d>::zeros(256, 512);

		if (coeffSH == nullptr)
			return result;

		if (order < 0 || order > 4)
			order = 2;

		int coeffCount = std::min(int(coeffSH->size()), (order + 1) * (order + 1));
		for (int x = 0; x < result.cols; x++)
		{
			double phi = PI2 * (x + 0.5) / result.cols;
			for (int y = 0; y < result.rows; y++)
			{
				double theta = PI * (y + 0.5) / result.rows;
				double r = std::sin(theta);
				cv::Vec3d dir(-r * std::sin(phi), std::cos(theta), r * cos(phi));

				for (int k = 0; k < coeffCount; k++)
					result.at<cv::Vec3d>(y, x) += (SHB[k](dir) * (*coeffSH)[k]);
			}
		}

		return result;
	}
}