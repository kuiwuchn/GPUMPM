#if 1
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <ctime>

#include "cySampleElim.h"
#include "cyPoint.h"

class SampleGenerator
{
public:
	SampleGenerator();
	~SampleGenerator();

	void test();
	void LoadSDF(std::string filename, float& pDx, float& minx, float& miny, float& minz, int& ni, int& nj, int& nk);
	void GeneratePoissonSamples(float samplesPerVol, std::vector<float>& outputSamples, int inputScale = 5);

	void GeneratePoissonSamples(int length, int width, int height, float scale, int outputSamplesNumber, std::vector<float>& outputSamples, int inputScale = 5);

	int GenerateUniformSamples(float samplesPerVol, std::vector<float>& outputSamples);

protected:
	inline float fetchGrid(int i, int j, int k) { return m_phiGrid[i + m_ni*(j + m_nj*k)];}
	inline float fetchGridTrilinear(float x, float y, float z)
	{
		float dx = x - floor(x);
		float dy = y - floor(y);
		float dz = z - floor(z);

		float c000 = fetchGrid((int)floor(x), (int)floor(y), (int)floor(z));
		float c001 = fetchGrid((int)floor(x), (int)floor(y), (int)ceil(z));
		float c010 = fetchGrid((int)floor(x), (int)ceil(y), (int)floor(z));
		float c011 = fetchGrid((int)floor(x), (int)ceil(y), (int)ceil(z));
		float c100 = fetchGrid((int)ceil(x), (int)floor(y), (int)floor(z));
		float c101 = fetchGrid((int)ceil(x), (int)floor(y), (int)ceil(z));
		float c110 = fetchGrid((int)ceil(x), (int)ceil(y), (int)floor(z));
		float c111 = fetchGrid((int)ceil(x), (int)ceil(y), (int)ceil(z));

		float c00 = c000 * (1 - dx) + c100 * dx;
		float c01 = c001 * (1 - dx) + c101 * dx;
		float c10 = c010 * (1 - dx) + c110 * dx;
		float c11 = c011 * (1 - dx) + c111 * dx;

		float c0 = c00 * (1 - dy) + c10 * dy;
		float c1 = c01 * (1 - dy) + c11 * dy;

		return c0 * (1 - dz) + c1 * dz;
	}
private:
	int					m_ni, m_nj, m_nk;
	float				m_dx;
	int					m_padding;
	cyPoint3f			m_minBox;
	std::vector<float>	m_phiGrid;
	cy::WeightedSampleElimination<cy::Point3f, float, 3, int> m_wse;
};

SampleGenerator::SampleGenerator()
{
}

SampleGenerator::~SampleGenerator()
{
}

void SampleGenerator::LoadSDF(std::string filename, float& pDx, float& minx, float& miny, float& minz, int& ni, int& nj, int& nk)
{
	std::ifstream infile(filename);
	infile >> m_ni >> m_nj >> m_nk;
	infile >> m_minBox[0] >> m_minBox[1] >> m_minBox[2];
	infile >> m_dx;

	//std::cout << m_dx << std::endl;

	pDx = m_dx;

	ni = m_ni;
	nj = m_nj;
	nk = m_nk;
	minx = m_minBox[0];
	miny = m_minBox[1];
	minz = m_minBox[2];

	std::cout << "load grid size " << m_ni << ", " << m_nj << ", " << m_nk << std::endl;

	int gridSize = m_ni * m_nj * m_nk;
	m_phiGrid.resize(gridSize);
	for (int i = 0; i < gridSize; ++i) {
		infile >> m_phiGrid[i];
	}
	infile.close();
}

void SampleGenerator::GeneratePoissonSamples(int length, int width, int height, float scale, int outputSamplesNumber, std::vector<float>& outputSamples, int inputScale)
{
	int numSamples = outputSamplesNumber;

	std::cout << "Generate input " << numSamples * inputScale << " samples...";
	std::vector<cy::Point3f> inputPoints(numSamples * inputScale);
	for (size_t i = 0; i < inputPoints.size(); i++) {
		inputPoints[i].x = (float)rand() / RAND_MAX * length;
		inputPoints[i].y = (float)rand() / RAND_MAX * width;
		inputPoints[i].z = (float)rand() / RAND_MAX * height;
	}
	std::cout << "done\n";

	std::cout << "Eliminate samples...";
	std::vector<cy::Point3f> tmpPoints;
	tmpPoints.resize(numSamples);
	m_wse.Eliminate(inputPoints.data(), inputPoints.size(), tmpPoints.data(), tmpPoints.size());
	std::cout << "done\n";

	for (int i = 0; i < (int)tmpPoints.size(); i++)
	{
		outputSamples.push_back(tmpPoints[i].x * scale);
		outputSamples.push_back(tmpPoints[i].y * scale);
		outputSamples.push_back(tmpPoints[i].z * scale);
	}
}

void SampleGenerator::GeneratePoissonSamples(float samplesPerVol, std::vector<float>& outputSamples, int inputScale)
{
	int numSamples = m_phiGrid.size() * samplesPerVol;

	std::cout << "Generate input " << numSamples * inputScale << " samples...";
	std::vector<cy::Point3f> inputPoints(numSamples * inputScale);
	for (size_t i = 0; i < inputPoints.size(); i++) {
		inputPoints[i].x = (float)rand() / RAND_MAX * (m_ni - 1);
		inputPoints[i].y = (float)rand() / RAND_MAX * (m_nj - 1);
		inputPoints[i].z = (float)rand() / RAND_MAX * (m_nk - 1);
	}
	std::cout << "done\n";

	std::cout << "Eliminate samples...";
	std::vector<cy::Point3f> tmpPoints;
	tmpPoints.resize(numSamples);
	m_wse.Eliminate(inputPoints.data(), inputPoints.size(), tmpPoints.data(), tmpPoints.size());
	std::cout << "done\n";

	outputSamples.clear();
	for (int i = 0; i < (int)tmpPoints.size(); i++)
	{
		if (fetchGridTrilinear(tmpPoints[i].x, tmpPoints[i].y, tmpPoints[i].z) < 0)
		{
			outputSamples.push_back(tmpPoints[i].x);
			outputSamples.push_back(tmpPoints[i].y);
			outputSamples.push_back(tmpPoints[i].z);
		}
	}
}

#if 1
int SampleGenerator::GenerateUniformSamples(float samplesPerCell, std::vector<float>& outputSamples)
{
	// get total sample number
	int validCellNum = 0;
	for (int i = 0; i < m_ni - 1; i++)
	{
		for (int j = 0; j < m_nj - 1; j++)
		{
			for (int k = 0; k < m_nk - 1; k++)
			{
				if (fetchGrid(i, j, k) < 0 ||
					fetchGrid(i, j, k + 1) < 0 ||
					fetchGrid(i, j + 1, k) < 0 ||
					fetchGrid(i, j + 1, k + 1) < 0 ||
					fetchGrid(i + 1, j, k) < 0 ||
					fetchGrid(i + 1, j, k + 1) < 0 ||
					fetchGrid(i + 1, j + 1, k) < 0 ||
					fetchGrid(i + 1, j + 1, k + 1) < 0)
					validCellNum++;
			}
		}
	}

	int sampleNum = validCellNum * samplesPerCell;

	for (int i = 0; i < sampleNum; i++)
	{
		cyPoint3f tmpPoint;
		do 
		{
			tmpPoint.x = (float)rand() / RAND_MAX * (m_ni - 1);
			tmpPoint.y = (float)rand() / RAND_MAX * (m_nj - 1);
			tmpPoint.z = (float)rand() / RAND_MAX * (m_nk - 1);
		} while (fetchGridTrilinear(tmpPoint.x, tmpPoint.y, tmpPoint.z) >= 0);

		outputSamples.push_back(tmpPoint.x);
		outputSamples.push_back(tmpPoint.y);
		outputSamples.push_back(tmpPoint.z);
	}

	return sampleNum;
}
#endif 

void SampleGenerator::test()
{
       /* std::vector< cy::Point3f > inputPoints(160000);*/
	//for (size_t i = 0; i < inputPoints.size(); i++) {
		//inputPoints[i].x = (float)rand() / RAND_MAX;
		//inputPoints[i].y = (float)rand() / RAND_MAX;
		//inputPoints[i].z = (float)rand() / RAND_MAX;
	//}
	
	//std::vector< cy::Point3f > outputPoints(32000);
	//std::clock_t start;
	//double duration;
	//start = std::clock();
	//m_dx = m_wse.Eliminate(inputPoints.data(), inputPoints.size(),
	//outputPoints.data(), outputPoints.size());

	//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//std::cout << "time = " << duration << '\n';

	//float minRadius = 255;
	//for (int i = 0; i < 12000; i++)
		//for (int j = 0; j < 12000; j++)
		//{
			//if (i == j) continue;
			//float tmp = (outputPoints[i] - outputPoints[j]).Length();
			//minRadius = fminf(minRadius, tmp);
		//}
	//std::cout << "min radius = " << minRadius << std::endl;
	//std::cout << m_dx << std::endl;
	/*std::cout << minRadius / m_dx << std::endl;*/

	//int gNi = 10, gNj = 10, gNk = 10;
	//int samplesPerSubcell = 12;
	//int inputScale = 5;

	//int sampleVolume = gNi * gNj * gNk;
	//int inputSampleCount = sampleVolume * samplesPerSubcell * inputScale;
	//int outputSampleCount = sampleVolume * samplesPerSubcell;

	//std::vector<cy::Point3f> inputSamples;
	//inputSamples.resize(inputSampleCount);

	//int c = 0;
	//for (int i = 0; i < gNi; i++) for (int j = 0; j < gNj; j++) for (int k = 0; k < gNk; k++)
	//	for (int si = 0; si < samplesPerSubcell * inputScale; si++)
	//	{
	//		inputSamples[c][0] = (float)rand() / RAND_MAX;
	//		inputSamples[c][1] = (float)rand() / RAND_MAX;
	//		inputSamples[c][2] = (float)rand() / RAND_MAX;
	//		c++;
	//	}

	//std::vector<cy::Point3f> tmpPoints;
	//tmpPoints.resize(outputSampleCount);
	//m_wse.Eliminate(inputSamples.data(), inputSamples.size(), tmpPoints.data(), tmpPoints.size());

	//float minRadius = 255;
	//for (int i = 0; i < outputSampleCount; i++)
	//	for (int j = 0; j < outputSampleCount; j++)
	//	{
	//		if (i == j) continue;
	//		float tmp = (tmpPoints[i] - tmpPoints[j]).Length();
	//		minRadius = fminf(minRadius, tmp);
	//	}
	//std::cout << "min radius = " << minRadius << std::endl;
}
#else
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <ctime>

#include "cySampleElim.h"
#include "cyPoint.h"

class SampleGenerator
{
public:
	SampleGenerator();
	~SampleGenerator();

	void test();
	void LoadSDF(std::string filename);
	void GeneratePoissonSamples(int samplesPerVol, std::vector<float>& outputSamples, int inputScale = 5);

	int GenerateUniformSamples(int samplesPerVol, std::vector<float>& outputSamples);

protected:
	inline float fetchGrid(int i, int j, int k) { return m_phiGrid[i + m_ni*(j + m_nj*k)];}
	inline float fetchGridTrilinear(float x, float y, float z)
	{
		float dx = x - floor(x);
		float dy = y - floor(y);
		float dz = z - floor(z);

		float c000 = fetchGrid((int)floor(x), (int)floor(y), (int)floor(z));
		float c001 = fetchGrid((int)floor(x), (int)floor(y), (int)ceil(z));
		float c010 = fetchGrid((int)floor(x), (int)ceil(y), (int)floor(z));
		float c011 = fetchGrid((int)floor(x), (int)ceil(y), (int)ceil(z));
		float c100 = fetchGrid((int)ceil(x), (int)floor(y), (int)floor(z));
		float c101 = fetchGrid((int)ceil(x), (int)floor(y), (int)ceil(z));
		float c110 = fetchGrid((int)ceil(x), (int)ceil(y), (int)floor(z));
		float c111 = fetchGrid((int)ceil(x), (int)ceil(y), (int)ceil(z));

		float c00 = c000 * (1 - dx) + c100 * dx;
		float c01 = c001 * (1 - dx) + c101 * dx;
		float c10 = c010 * (1 - dx) + c110 * dx;
		float c11 = c011 * (1 - dx) + c111 * dx;

		float c0 = c00 * (1 - dy) + c10 * dy;
		float c1 = c01 * (1 - dy) + c11 * dy;

		return c0 * (1 - dz) + c1 * dz;
	}
private:
	int					m_ni, m_nj, m_nk;
	float				m_dx;
	int					m_padding;
	cyPoint3f			m_minBox;
	std::vector<float>	m_phiGrid;
	cy::WeightedSampleElimination<cy::Point3f, float, 3, int> m_wse;
};

SampleGenerator::SampleGenerator()
{
}

SampleGenerator::~SampleGenerator()
{
}

void SampleGenerator::LoadSDF(std::string filename)
{
	std::ifstream infile(filename);
	infile >> m_ni >> m_nj >> m_nk;
	infile >> m_minBox[0] >> m_minBox[1] >> m_minBox[2];
	infile >> m_dx;

	std::cout << "load grid size " << m_ni << ", " << m_nj << ", " << m_nk << std::endl;

	int gridSize = m_ni * m_nj * m_nk;
	m_phiGrid.resize(gridSize);
	for (int i = 0; i < gridSize; ++i) {
		infile >> m_phiGrid[i];
	}
	infile.close();
}

void SampleGenerator::GeneratePoissonSamples(int samplesPerVol, std::vector<float>& outputSamples, int inputScale)
{
	int numSamples = m_phiGrid.size() * samplesPerVol;

	std::cout << "Generate input " << numSamples * inputScale << " samples...";
	std::vector<cy::Point3f> inputPoints(numSamples * inputScale);
	for (size_t i = 0; i < inputPoints.size(); i++) {
		inputPoints[i].x = (float)rand() / RAND_MAX * (m_ni - 1);
		inputPoints[i].y = (float)rand() / RAND_MAX * (m_nj - 1);
		inputPoints[i].z = (float)rand() / RAND_MAX * (m_nk - 1);
	}
	std::cout << "done\n";

	std::cout << "Eliminate samples...";
	std::vector<cy::Point3f> tmpPoints;
	tmpPoints.resize(numSamples);
	m_wse.Eliminate(inputPoints.data(), inputPoints.size(), tmpPoints.data(), tmpPoints.size());
	std::cout << "done\n";

	outputSamples.clear();
	for (int i = 0; i < (int)tmpPoints.size(); i++)
	{
		if (fetchGridTrilinear(tmpPoints[i].x, tmpPoints[i].y, tmpPoints[i].z) < 0)
		{
			outputSamples.push_back(tmpPoints[i].x);
			outputSamples.push_back(tmpPoints[i].y);
			outputSamples.push_back(tmpPoints[i].z);
		}
	}
}

int SampleGenerator::GenerateUniformSamples(int samplesPerCell, std::vector<float>& outputSamples)
{
	// get total sample number
	int validCellNum = 0;
	for (int i = 0; i < m_ni - 1; i++)
	{
		for (int j = 0; j < m_nj - 1; j++)
		{
			for (int k = 0; k < m_nk - 1; k++)
			{
				if (fetchGrid(i, j, k) < 0 ||
					fetchGrid(i, j, k + 1) < 0 ||
					fetchGrid(i, j + 1, k) < 0 ||
					fetchGrid(i, j + 1, k + 1) < 0 ||
					fetchGrid(i + 1, j, k) < 0 ||
					fetchGrid(i + 1, j, k + 1) < 0 ||
					fetchGrid(i + 1, j + 1, k) < 0 ||
					fetchGrid(i + 1, j + 1, k + 1) < 0)
					validCellNum++;
			}
		}
	}

	int sampleNum = validCellNum * samplesPerCell;

	for (int i = 0; i < sampleNum; i++)
	{
		cyPoint3f tmpPoint;
		do 
		{
			tmpPoint.x = (float)rand() / RAND_MAX * (m_ni - 1);
			tmpPoint.y = (float)rand() / RAND_MAX * (m_nj - 1);
			tmpPoint.z = (float)rand() / RAND_MAX * (m_nk - 1);
		} while (fetchGridTrilinear(tmpPoint.x, tmpPoint.y, tmpPoint.z) >= 0);

		outputSamples.push_back(tmpPoint.x);
		outputSamples.push_back(tmpPoint.y);
		outputSamples.push_back(tmpPoint.z);
	}

	return sampleNum;
}

void SampleGenerator::test()
{
	std::vector< cy::Point3f > inputPoints(160000);
	for (size_t i = 0; i < inputPoints.size(); i++) {
		inputPoints[i].x = (float)rand() / RAND_MAX;
		inputPoints[i].y = (float)rand() / RAND_MAX;
		inputPoints[i].z = (float)rand() / RAND_MAX;
	}
	
	std::vector< cy::Point3f > outputPoints(32000);
	std::clock_t start;
	double duration;
	start = std::clock();
	m_wse.Eliminate(inputPoints.data(), inputPoints.size(),
	outputPoints.data(), outputPoints.size());

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "time = " << duration << '\n';

	float minRadius = 255;
	for (int i = 0; i < 12000; i++)
		for (int j = 0; j < 12000; j++)
		{
			if (i == j) continue;
			float tmp = (outputPoints[i] - outputPoints[j]).Length();
			minRadius = fminf(minRadius, tmp);
		}
	std::cout << "min radius = " << minRadius << std::endl;

	//int gNi = 10, gNj = 10, gNk = 10;
	//int samplesPerSubcell = 12;
	//int inputScale = 5;

	//int sampleVolume = gNi * gNj * gNk;
	//int inputSampleCount = sampleVolume * samplesPerSubcell * inputScale;
	//int outputSampleCount = sampleVolume * samplesPerSubcell;

	//std::vector<cy::Point3f> inputSamples;
	//inputSamples.resize(inputSampleCount);

	//int c = 0;
	//for (int i = 0; i < gNi; i++) for (int j = 0; j < gNj; j++) for (int k = 0; k < gNk; k++)
	//	for (int si = 0; si < samplesPerSubcell * inputScale; si++)
	//	{
	//		inputSamples[c][0] = (float)rand() / RAND_MAX;
	//		inputSamples[c][1] = (float)rand() / RAND_MAX;
	//		inputSamples[c][2] = (float)rand() / RAND_MAX;
	//		c++;
	//	}

	//std::vector<cy::Point3f> tmpPoints;
	//tmpPoints.resize(outputSampleCount);
	//m_wse.Eliminate(inputSamples.data(), inputSamples.size(), tmpPoints.data(), tmpPoints.size());

	//float minRadius = 255;
	//for (int i = 0; i < outputSampleCount; i++)
	//	for (int j = 0; j < outputSampleCount; j++)
	//	{
	//		if (i == j) continue;
	//		float tmp = (tmpPoints[i] - tmpPoints[j]).Length();
	//		minRadius = fminf(minRadius, tmp);
	//	}
	//std::cout << "min radius = " << minRadius << std::endl;
}
#endif
