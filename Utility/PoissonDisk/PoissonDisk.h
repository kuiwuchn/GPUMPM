#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "cySampleElim.h"
#include "cyPoint.h"

class PoissonDisk
{
public:
	PoissonDisk();
	~PoissonDisk();

	void LoadSDF(std::string filename);
	void GenerateSamples(int samplesPerVol, std::vector<float>& outputSamples, int inputScale = 5);

protected:
	inline float fetchGrid(int i, int j, int k) { return m_phiGrid[i + m_ni*(j + m_nj*k)];}
	inline float fetchGridTrilinear(float x, float y, float z)
	{
		float dx = x - floor(x);
		float dy = y - floor(y);
		float dz = z - floor(z);

		float c000 = fetchGrid(floor(x), floor(y), floor(z));
		float c001 = fetchGrid(floor(x), floor(y), ceil(z));
		float c010 = fetchGrid(floor(x), ceil(y), floor(z));
		float c011 = fetchGrid(floor(x), ceil(y), ceil(z));
		float c100 = fetchGrid(ceil(x), floor(y), floor(z));
		float c101 = fetchGrid(ceil(x), floor(y), ceil(z));
		float c110 = fetchGrid(ceil(x), ceil(y), floor(z));
		float c111 = fetchGrid(ceil(x), ceil(y), ceil(z));

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

PoissonDisk::PoissonDisk()
{
}

PoissonDisk::~PoissonDisk()
{
}

void PoissonDisk::LoadSDF(std::string filename)
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

void PoissonDisk::GenerateSamples(int samplesPerVol, std::vector<float>& outputSamples, int inputScale)
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
