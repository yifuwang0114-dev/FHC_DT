#include "dt.h"
//-------------------------mesh improvement-------------------------

// Smooth for Volume-based smoothing
int DT::smooth_volume(int iNod, bool equalAngle)
{
	auto t1 = getTime_now();
	int i, m, k, iElem, iSecond;
	double newpos[3], oldEnergy, newEnergy = 1e10, minq = 1e10;
	double *verts[4];
	std::vector<int> sph;
	std::vector<double> area, oldVolume, newVolume;
	double oldpos[3] = {Nodes[iNod].pt[0], Nodes[iNod].pt[1], Nodes[iNod].pt[2]};

	findSphere(iNod, sph);
	int n = sph.size();

	if (!improve_step)
	{
		for (i = 0; i < n; i++)
		{
			if (ishulltet(sph[i]))
			{
				return 0;
			}

			Elems[sph[i]].q = calVolume(sph[i]);
		}
	}

	// get worse quality
	if (!equalAngle)
	{
		area.resize(0);
	}
	else
	{
		area.resize(n);
	}
	newVolume.resize(n);
	oldVolume.resize(n);

	for (i = 0; i < n; i++)
	{
		iElem = sph[i];
		oldVolume[i] = Elems[iElem].q;

		if (equalAngle)
		{
			for (m = 0; m < 4; m++)
			{
				iSecond = Elems[iElem].form[m];
				if (iSecond == iNod)
				{
					area[i] = calArea(Nodes[Elems[iElem].form[(m + 1) % 4]].pt,
									  Nodes[Elems[iElem].form[(m + 2) % 4]].pt, Nodes[Elems[iElem].form[(m + 3) % 4]].pt);
				}
			}
		}
	}
	oldEnergy = getVolEnergy(oldVolume, area);

	int descentNum = 0;
	while (1)
	{
		descentNum++; // Number of descents
		double VolGrad_initial[3] = {0};
		double VolGrad[3] = {0};
		double Hessian[9] = {0};
		double HessianT[9] = {0};
		getVolGrad(iNod, sph, VolGrad_initial, area);
		getHessian(iNod, sph, Hessian, area);
		if (!inverseM(Hessian, HessianT))
		{
			return 0;
		}
		vecTimesMatrix13_33(VolGrad_initial, HessianT, VolGrad);

		double alpha = 1;
		double beta = 0.8;
		double gamma = 1e-4; // 0.01;
		while (1)
		{
			// new position
			for (i = 0; i < 3; i++)
				newpos[i] = Nodes[iNod].pt[i] - VolGrad[i] * alpha;
			for (i = 0; i < n; i++)
			{
				for (int m = 0; m <= 3; m++)
				{
					int iElemNd = Elems[sph[i]].form[m];
					verts[m] = (iElemNd != iNod) ? Nodes[iElemNd].pt : newpos;
				}
				newVolume[i] = tetquality(verts[0], verts[1], verts[2], verts[3], improve_Metric);
				if (newVolume[i] < 0)
					break;
			}

			if (i == n)
			{
				// all tet positive now
				newEnergy = getVolEnergy(newVolume, area);
				if (newEnergy <= oldEnergy + gamma * alpha * dot(VolGrad, VolGrad_initial))
				{
					break;
				}
			}

			alpha *= beta;

			if (alpha < /*1e-8*/ 1e-10)
				break;
		}
		if (newEnergy < oldEnergy)
		{
			// can't improvement more
			for (i = 0; i < n; i++)
			{
				Elems[sph[i]].q = newVolume[i];
			}

			for (int i = 0; i < 3; i++)
				Nodes[iNod].pt[i] = newpos[i];
			// if (newEnergy / oldEnergy > 0.9999) {
			if (fabs((newEnergy - oldEnergy) / oldEnergy) < 1e-5)
			{
				auto t2 = getTime_now();
				volumesmoothtime += getTime(t1, t2);
				return 1;
			}
			else
			{
				oldEnergy = newEnergy;
			}
		}
		else
		{
			if (descentNum > 1)
			{
				auto t2 = getTime_now();
				volumesmoothtime += getTime(t1, t2);
				return 1;
			}
			else
			{
				auto t2 = getTime_now();
				volumesmoothtime += getTime(t1, t2);
				return 0;
			}
		}
	}
	auto t2 = getTime_now();
	volumesmoothtime += getTime(t1, t2);
	return 0;
}

int DT::getVolGrad(int iNod, std::vector<int> sph, double *VolGrad, std::vector<double> area)
{
	VolGrad[0] = VolGrad[1] = VolGrad[2] = 0;
	for (int i = 0; i < sph.size(); i++)
	{
		int j = 0, a, b, c, d;
		for (; j < 4; j++)
		{
			if (Elems[sph[i]].form[j] == iNod)
				break;
		}
		DFC(j, d, a, b, c);
		double *pa = Nodes[Elems[sph[i]].form[a]].pt;
		double *pb = Nodes[Elems[sph[i]].form[b]].pt;
		double *pc = Nodes[Elems[sph[i]].form[c]].pt;
		double papb[3] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
		double papc[3] = {pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
		double ans[3] = {0};
		cross(papb, papc, ans);
		if (area.size() == 0)
		{
			VolGrad[0] += Elems[sph[i]].q * ans[0] / 3.0;
			VolGrad[1] += Elems[sph[i]].q * ans[1] / 3.0;
			VolGrad[2] += Elems[sph[i]].q * ans[2] / 3.0;
		}
		else
		{
			VolGrad[0] += Elems[sph[i]].q * ans[0] / area[i] / 3.0;
			VolGrad[1] += Elems[sph[i]].q * ans[1] / area[i] / 3.0;
			VolGrad[2] += Elems[sph[i]].q * ans[2] / area[i] / 3.0;
		}
	}
	return 0;
}

int DT::getHessian(int iNod, std::vector<int> sph, double *Hessian, std::vector<double> area)
{
	Hessian[0] = Hessian[1] = Hessian[2] = Hessian[3] = Hessian[4] = Hessian[5] = Hessian[6] = Hessian[7] = Hessian[8] = 0;
	for (int i = 0; i < sph.size(); i++)
	{
		int j = 0, a, b, c, d;
		for (; j < 4; j++)
		{
			if (Elems[sph[i]].form[j] == iNod)
				break;
		}
		DFC(j, d, a, b, c);
		double *pa = Nodes[Elems[sph[i]].form[a]].pt;
		double *pb = Nodes[Elems[sph[i]].form[b]].pt;
		double *pc = Nodes[Elems[sph[i]].form[c]].pt;
		double papb[3] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
		double papc[3] = {pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
		double ans[3] = {0};
		double hes[9] = {0};
		cross(papb, papc, ans);
		tensorproduct33(ans, ans, hes);
		if (area.size() == 0)
		{
			for (j = 0; j < 9; j++)
			{
				Hessian[j] += hes[j] / 18.0;
			}
		}
		else
		{
			for (j = 0; j < 9; j++)
			{
				Hessian[j] += hes[j] / area[i] / 18.0;
			}
		}
	}

	for (int r = 0; r < 3; ++r)
	{
		for (int c = r + 1; c < 3; ++c)
		{
			double sym = 0.5 * (Hessian[r * 3 + c] + Hessian[c * 3 + r]);
			Hessian[r * 3 + c] = Hessian[c * 3 + r] = sym;
		}
	}
	return 0;
}

double DT::getVolEnergy(std::vector<double> volume, std::vector<double> area)
{
	double energy = 0;
	if (area.size() == 0)
	{
		for (int i = 0; i < volume.size(); i++)
			energy += volume[i] * volume[i];
	}
	else
	{
		for (int i = 0; i < volume.size(); i++)
			energy += volume[i] * volume[i] / area[i];
	}
	return energy;
}
