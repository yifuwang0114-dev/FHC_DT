#include "dt.h"

using namespace std;

DT::DT() {};
DT::~DT() {};

/* Main process*/
int DT::tetrahedralize(Mesh &mesh, Args &args)
{
	auto t_1 = getTime_now();
	dt_init(args);
	if (!BndPntInst(mesh, args))
		return 1;
	BoudaryRecover(mesh, args);
	ColorVirtualTet(args);
	MeshRefine(args);

	auto t_2 = getTime_now();

	if (infolevel > 0)
	{
		printf("Mesh Generation cost  : %.3f s\n", getTime(t_1, t_2));
		printf("Mesh Generation create: %d tet\n", (int)Elems.size());
		printf("Mesh Generation Speed : %.3f W/s\n", Elems.size() / 10000.0 / getTime(t_1, t_2));
	}
	auto t_3 = getTime_now();

	MeshImprove(args);
	auto t_4 = getTime_now();
	RemoveTet(args);
	outMesh(mesh, args);

	if (infolevel > 0)
	{
		printf("Mesh improve cost     : %.3f s\n", getTime(t_3, t_4));
		printf("Mesh improve Speed    : %.3f W/s\n", mesh.T.size() / 10000.0 / getTime(t_3, t_4));
		printf("Total cost            : %.3f s\n", getTime(t_1, t_4));
		printf("Final tet num         : %zu tet\n", mesh.T.size());
		printf("DT speed              : %.3f W/s\n", mesh.T.size() / 10000.0 / getTime(t_1, t_4));
		printf("Memory cost           : %.3f MB\n", getPeakMegabytesUsed());
		printf("Tetrahedralize success!\n");
	}
	return 1;
}

/*Boundary Recover process*/
int DT::BoudaryRecover(Mesh &mesh, Args &args)
{
	auto t_begin = getTime_now();
	if (infolevel > 0)
		spdlog::info("Recovering boundaries.");

	auto t1 = getTime_now();
	// Read input Faces
	buildBndInfo(mesh, args);

	/************************** Recover Edges **************************/
	AutorecoverEdges(args);

	auto t2 = getTime_now();
	/************************** Recover Faces **************************/
	recoverFacesPass(args);

	auto t3 = getTime_now();

	removeStPass();

	auto t4 = getTime_now();

	auto t_end = getTime_now();
	if (infolevel > 0)
	{
		spdlog::info("Recovering Edge: {:.6f}", getTime(t1, t2));
		spdlog::info("Recovering Face: {:.6f}", getTime(t2, t3));
		spdlog::info("Remove Boundary st: {:.6f}", getTime(t3, t4));
		spdlog::info("Add steiner point: {}", addst);
		spdlog::info("Add Boundary point: {}", addstbnd);
		spdlog::info("Remove steiner point: {}", rmvst);
		spdlog::info("Add inner st: {:.6f}", addinnersttime);
		spdlog::info("Remove inner st: {:.6f}", Removeinnersttime);
		spdlog::info("Volume Smooth: {:.6f}", volumesmoothtime);
		spdlog::info("Recovering boundaries cost time: {:.3f}", getTime(t_begin, t_end));
	}

	return 0;
}

/*Segment Recover Process*/
int DT::AutorecoverEdges(Args &args)
{
	std::queue<int> lost;
	std::map<int, int> N_lostE_pre;

	// find all lost edges
	for (int i = 0; i < Elems.size(); i++)
	{
		if (isDelEle(i))
			continue;
		for (int j = 0; j < 6; j++)
		{
			int *edgeid = BndEdg.find(Elems[i].form[Egid[j][0]], Elems[i].form[Egid[j][1]]);
			if (edgeid)
			{
				SurEdgs[*edgeid].info = 1;
			}
		}
	}

	for (int i = 0; i < SurEdgs.size(); i++)
	{
		if (!isRecBndEdg(i))
		{
			lost.push(i);
			N_lostE_pre[SurEdgs[i].iStart]++;
			N_lostE_pre[SurEdgs[i].iEnd]++;
		}
	}

	if (infolevel > 0)
		spdlog::info("Lost Edges: {} / {}", (int)lost.size(), (int)SurEdgs.size());

	/**************************** Start Recover Edge *******************************/
	int tempfliplevel = 1;
	// int BRlevel = 0;
	while (lost.size() != 0)
	{
		int nlost = lost.size();
		int success = 0;

		for (int i = 0; i < nlost; i++)
		{
			int te = lost.front();
			lost.pop();
			if (isDelSurEdg(te))
				continue;

			int fullsearch = 0;
			int stflag = 0;
			fliplevel = tempfliplevel - SurEdgs[te].info * 10;

			if (SurEdgs[te].info <= -3)
			{
				fliplevel = std::max(1000, fliplevel);
				fullsearch = 1;
			}

			if (SurEdgs[te].info <= -4)
				stflag = 1;
			if (SurEdgs[te].info <= -5)
				stflag = 2;

			int ret = recoverEdge(te, fullsearch, stflag);

			if (ret == 0)
			{ // recover edge fail
				lost.push(te);
			}
			else if (ret > 1)
			{ // recover edge fail,and split it
				for (int j = ret; j < SurEdgs.size(); j++)
				{
					if (recoverEdge(j, 1, 0) == 0)
					{
						SurEdgs[j].info = -3;
						lost.push(j); // add new edge wait for recover
					}
				}
				success++;
			}
			else if (ret == 1)
			{
				success++;
			}
		}

		for (int iloop = 0; iloop < 3; iloop++)
		{
			int nrmv = 0;
			for (int i = nSurNodes; i < Nodes.size(); i++)
			{
				if (isDelNod(i) || isbndpnt(i) || i == ghost)
				{
					continue;
				}
				if (!removePnt(i))
					smooth_volume(i, true);
			}
			if (nrmv == 0)
				break;
		}

		updateFliptype(N_lostE_pre, lost);

		if (infolevel > 1)
			spdlog::info("level:{} lost edges:{}", tempfliplevel, (int)lost.size());
		tempfliplevel++;
		if (tempfliplevel > 1000)
		{
			if (spdlog::get_level() != spdlog::level::off)
				printf("Boundary recovery failed!\n");
			spdlog::error("Boundary recovery failed!");
			throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
		}
	}

	for (int iloop = 0; iloop < 3; iloop++)
	{
		int nrmv = 0;
		for (int i = nSurNodes; i < Nodes.size(); i++)
		{
			if (isDelNod(i) || isbndpnt(i) || i == ghost)
			{
				continue;
			}
			if (!removePnt(i))
				smooth_volume(i, true);
			else
				nrmv++;
		}
		if (nrmv == 0)
			break;
	}

	seg[0] = seg[1] = -1;

	if (lost.empty() && infolevel > 0)
		spdlog::info("Finish Recover Edges.");

	return 0;
}

/*Segment Recover Type Update*/
void DT::updateFliptype(std::map<int, int> &N_lostE_pre, std::queue<int> &lost)
{
	std::queue<int> tempQ = lost;
	std::map<int, int> N_lostE_now;
	while (!tempQ.empty())
	{
		int x = tempQ.front();
		tempQ.pop();
		N_lostE_now[SurEdgs[x].iStart]++;
		N_lostE_now[SurEdgs[x].iEnd]++;
	}

	tempQ = lost;
	while (!tempQ.empty())
	{
		int x = tempQ.front();
		tempQ.pop();
		int pa = SurEdgs[x].iStart;
		int pb = SurEdgs[x].iEnd;
		int olda = N_lostE_pre[pa];
		int oldb = N_lostE_pre[pb];
		int newa = N_lostE_now[pa];
		int newb = N_lostE_now[pb];

		std::swap(SurEdgs[x].iStart, SurEdgs[x].iEnd);

		if (olda <= newa && oldb <= newb)
		{
			SurEdgs[x].info = std::max(SurEdgs[x].info - 1, -90);
			if (newa + newb > 100)
			{
				SurEdgs[x].info = std::min(SurEdgs[x].info, -4);
			}
		}
		else
		{
			SurEdgs[x].info = std::min(SurEdgs[x].info + 1, 0); // 0
		}
	}

	N_lostE_pre.clear();
	for (auto it : N_lostE_now)
	{
		N_lostE_pre[it.first] = it.second;
	}

	return;
}

/* Recover a segment process*/
int DT::recoverEdge(const int targetE, int fullsearch, int info)
{
	if (isDelSurEdg(targetE) || isRecBndEdg(targetE))
		return 1;

	SurEdg *lostE = &SurEdgs[targetE];
	seg[0] = lostE->iStart;
	seg[1] = lostE->iEnd;
	int filpdepth = fliplevel, ret = 0;
	std::vector<int> newN;

	ret = recoverEdgebyFlip(targetE, 0, filpdepth << 1 | 0);

	if (ret == 1)
	{
		lostE->info = 1; // set this edge is recovered
		return 1;
	}
	else if (ret > 1)
		return ret;

	ret = recoverEdgebyFlip(targetE, 1, filpdepth << 1 | 0);

	if (ret == 1)
	{
		lostE->info = 1; // set this edge is recovered
		return 1;
	}
	else if (ret > 1)
		return ret;

	// fullsearch
	if (fullsearch == 1)
	{
		ret = recoverEdgebyFlip(targetE, 0, filpdepth << 1 | 1);

		if (ret == 1)
		{
			lostE->info = 1; // set this edge is recovered
			return 1;
		}
		else if (ret > 1)
			return ret;
	}

	if (info > 0)
	{
		ret = FHCSteinerInsert(targetE, newN, 0);
		if (ret == 0 && info > 1)
			ret = splitBndEdge(targetE, 1); // intersec point plus

		for (auto iNod : newN)
		{
			if (isDelNod(iNod))
				continue;
			if (!removePnt(iNod, 20))
			{
				smooth_volume(iNod, true);
			}
		}
	}

	return ret;
}

/* Recover a segment*/
int DT::recoverEdgebyFlip(const int targetE, int dirflag, int info)
{
	int dir, srctet, ia, ib, ic, id, pa, pb, pc, pd, i, tried = 0;
	SurEdg *lostE = &SurEdgs[targetE];
	int p1 = lostE->iStart;
	int p2 = lostE->iEnd;
	if (dirflag == 1)
		std::swap(p1, p2);

	while (1)
	{
		if (isDelSurEdg(targetE))
			return 2;
		// find search direction
		dir = finddirection(p1, p2, srctet);

		if (dir == -20 || tried++ > 1000)
		{
			dir = finddirection(p1, p2, srctet);
			// some error happen
			return 0;
		}
		else if (0 <= dir && dir <= 3)
		{ // Across Vertex
			int IntersectPnt = Elems[srctet].form[dir];
			if (IntersectPnt == p2)
				return 1; // Success recover
			else
			{
				if (!isbndpnt(IntersectPnt))
				{
					if (removePnt(IntersectPnt))
					{
						// success
						continue;
					}
					if (disturbPnt(IntersectPnt))
					{
						continue;
					}
					else
					{
						return splitBndEdge(targetE, -IntersectPnt);
					}
				}
				if (spdlog::get_level() != spdlog::level::off)
					printf("Colline happen, %d located in %d,%d\n", IntersectPnt, p1, p2);
				spdlog::error("Colline happen,{} located in {} {}", IntersectPnt, p1, p2);
				// outTempMesh("./temp.vtk");
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
		}
		else if (4 <= dir && dir <= 7)
		{ // Across Face
			dir -= 4;
			DNC(dir, id, ia, ib, ic);
			pa = Elems[srctet].form[ia];
			pb = Elems[srctet].form[ib];
			pc = Elems[srctet].form[ic];
			if (isBndTri(pa, pb, pc))
			{ // chcek if it is bnd Tri
				spdlog::error("Intersection Face: {}:{} {} {} | Edge:{} {}", BndTri.get(pa, pb, pc), pa, pb, pc, p1, p2);
				// outTempMesh("./temp.vtk");
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}
			std::vector<int> rF = {srctet}; // store init and new tet

			if (removeface(rF, dir, info >> 1) == 1)
			{
				continue; // success remove a face,go to next loop
			}
		}
		else if (-14 <= dir && dir <= -1)
		{ // Across Edge
			ia = ((-dir) >> 2) & 3;
			ib = (-dir) & 3;
			pa = Elems[srctet].form[ia];
			pb = Elems[srctet].form[ib];
			if (isBndEdg(pa, pb))
			{ // chcek if it is bnd Tri

				spdlog::error("Intersection Edge: {} {} | {} {}", pa, pb, p1, p2);
				throw EXCEPTIONSTRING(std::string("error exit in") + std::string(__FILE__) + std::to_string(__LINE__));
			}

			std::vector<int> rE = {srctet}; // store init and new tet
			if (removeEdge(rE, ia, ib, info >> 1))
			{
				continue;
			}
		}
		else
		{
			spdlog::error("Error happen when recover Edge:{} {}", p1, p2);
		}

		/*if info&1==1 ,found all face and edge intersect and remove them*/
		if ((info & 1) == 1)
		{
			std::vector<std::array<int, 3>> Vid;
			int ret = FindIntersect(targetE, Vid);
			int sus = 0;
			std::vector<int> mp;

			auto seen = [&](int x)
			{
				return std::find(mp.begin(), mp.end(), x) != mp.end();
			};

			for (const auto &it : Vid)
			{
				if (seen(it[0]) || seen(it[1]) || seen(it[2]))
					continue;

				if (it[2] != -1)
				{
					// try to remove face
					if (isMeshFace(it[0], it[1], it[2], &srctet))
					{
						for (int t = 0; t < 4; t++)
						{
							if (Elems[srctet].form[t] != it[0] && Elems[srctet].form[t] != it[1] && Elems[srctet].form[t] != it[2])
							{
								dir = t;
								break;
							}
						}
						std::vector<int> rf = {srctet};
						if (removeface(rf, dir, /* info >> 1*/ 1) == 1)
							sus++;
						else
						{
							mp.push_back(it[0]);
							mp.push_back(it[1]);
							mp.push_back(it[2]);
						}
					}
				}
				else
				{
					// try to remove edge
					if (isMeshEdge(it[0], it[1], &srctet))
					{
						for (int t = 0; t < 4; t++)
						{
							if (Elems[srctet].form[t] == it[0])
							{
								ia = t;
							}
							else if (Elems[srctet].form[t] == it[1])
							{
								ib = t;
							}
						}
						std::vector<int> re = {srctet};
						ret = removeEdge(re, ia, ib, /*info >> 1*/ 1);

						if (ret == 1)
							sus++;
						else
						{
							mp.push_back(it[0]);
							mp.push_back(it[1]);
						}
					}
				}
			}
			if (sus > 0)
				return recoverEdgebyFlip(targetE, dirflag, info);
		}
		break;
	}
	return 0;
}

/* FHC- based steiner point insertion*/
int DT::FHCSteinerInsert(const int lostE, std::vector<int> &newN, int info)
{
	std::vector<std::array<int, 3>> Vid;
	if (FindIntersect(lostE, Vid))
	{ // have been recovered
		return 1;
	}

	if ((Vid.size() >= info && info != 0) || Vid.size() == 0)
		return 0;

	int pa = SurEdgs[lostE].iStart;
	int pb = SurEdgs[lostE].iEnd;
	double *p1 = Nodes[pa].pt;
	double *p2 = Nodes[pb].pt;
	double LenEdg = distance(p1, p2);
	double space = (Nodes[pa].space + Nodes[pb].space) / 2.0;
	double linep[2][3], facept[3][3], intPnt[3], pnt[3];
	int intTyp, intCod, srctet, nVid = Vid.size();

	for (int k = 0; k < 3; k++)
	{
		linep[0][k] = p1[k];
		linep[1][k] = p2[k];
	}

	int nume = 0, numf = 0;
	for (auto it : Vid)
	{
		if (it[2] == -1)
			nume++;
		else
			numf++;
	}

	// printf("lostE:%d  Vid:%d  info:%d  e:%d  f:%d\n",lostE, Vid.size(),info, nume,numf);

	for (auto it : Vid)
	{
		if (it[2] != -1)
		{ // Face
			if (isMeshFace(it[0], it[1], it[2], &srctet))
			{
				std::vector<int> oldtet = {srctet};
				int dir = isNod_in_Tet(Elems[srctet].form[0] + Elems[srctet].form[1] + Elems[srctet].form[2] + Elems[srctet].form[3] - it[0] - it[1] - it[2], srctet);
				int a, b, c, d;
				DFC(dir, a, b, c, d);
				b = Elems[srctet].form[b];
				c = Elems[srctet].form[c];
				d = Elems[srctet].form[d];
				int adddir = removeface(oldtet, dir, 10, -1);

				if (adddir < 0)
				{

					int rmvEdg = -adddir - 1, pc = -1, pd = -1;
					if (rmvEdg == 0)
					{
						pc = b, pd = d;
					}
					if (rmvEdg == 1)
					{
						pc = d, pd = c;
					}
					if (rmvEdg == 2)
					{
						pc = c, pd = b;
					}

					// try to remove face
					for (int k = 0; k < 3; k++)
					{
						facept[0][k] = Nodes[it[0]].pt[k];
						facept[1][k] = Nodes[it[1]].pt[k];
						facept[2][k] = Nodes[it[2]].pt[k];
					}

					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);

					int newp = addNode();
					Nodes[newp].space = space;

					for (int j = 0; j < 3; j++)
						Nodes[newp].pt[j] = (Nodes[pc].pt[j] + Nodes[pd].pt[j] + intPnt[j]) / 3.0;

					isMeshFace(it[0], it[1], it[2], &srctet);
					dir = isNod_in_Tet(Elems[srctet].form[0] + Elems[srctet].form[1] + Elems[srctet].form[2] + Elems[srctet].form[3] - it[0] - it[1] - it[2], srctet);

					int BW_tet = getP2T(pa);
					int loc = locate_pnt(newp, BW_tet);
					if (loc < 1)
						continue;

					std::vector<int> BW_vec = {srctet, getNeig(srctet, dir)};

					double ori1 = calVolume(BW_vec[0]);
					double ori2 = calVolume(BW_vec[1]);
					double d1 = distance(Nodes[pc].pt, Nodes[newp].pt);
					double d2 = distance(Nodes[pd].pt, Nodes[newp].pt);
					double d3 = distance(intPnt, Nodes[newp].pt);

					int ret = BW_insert_vertex(newp, BW_vec, 3);

					if (ret <= 0)
						DelNod(newp);
					else
						newN.push_back(newp);
				}

				if (recoverEdge(lostE, 0, 0) == 1)
				{
					return 1;
				}
			}
		}
		else
		{ // Edge
			if (isMeshEdge(it[0], it[1], &srctet))
			{
				std::vector<int> shell, shellp;
				findShell(srctet, isNod_in_Tet(it[0], srctet), isNod_in_Tet(it[1], srctet), shell, shellp);

				bool hullflag = false;
				for (int i = 0; i < shellp.size(); i++)
				{
					if (shellp[i] == ghost)
					{
						hullflag = true;
						break;
					}
				}

				for (int i = 0; i < shellp.size(); i++)
				{
					if (shellp[i] == ghost || shellp[i] == pa || shellp[i] == pb)
						continue;
					for (int j = 0; j < 3; j++)
					{
						facept[0][j] = Nodes[it[0]].pt[j];
						facept[1][j] = Nodes[it[1]].pt[j];
						facept[2][j] = Nodes[shellp[i]].pt[j];
					}
					dt::GEOM_FUNC::lin_tri_intersect3d(linep, facept, &intTyp, &intCod, intPnt);

					if (intCod != 0)
						continue;

					int px = it[1];
					if (distance2(intPnt, Nodes[it[0]].pt) > distance2(intPnt, Nodes[it[1]].pt))
						px = it[0];

					int newp = addNode();
					Nodes[newp].space = space;

					for (int k = 0; k < 3; k++)
						Nodes[newp].pt[k] = (intPnt[k] + Nodes[px].pt[k]) / 2.0;

					int BW_tet = getP2T(pa);
					int loc = locate_pnt(newp, BW_tet);

					if (loc < 1 || loc == 100 || ishulltet(BW_tet))
					{
						DelNod(newp);
						continue;
					}

					if (BW_insert_vertex(newp, shell, 3) <= 0)
					{
						DelNod(newp);
						continue;
					}
					else
						newN.push_back(newp);

					// direction
					double norm[3] = {0};
					calnormal(pa, pb, newp, norm);
					double normLen = lenvec(norm);
					for (int j = 0; j < 3; j++)
						norm[j] /= normLen;

					double SA[3] = {Nodes[shellp[i]].pt[0] - p1[0], Nodes[shellp[i]].pt[1] - p1[1], Nodes[shellp[i]].pt[2] - p1[2]};
					double nrom_sa = dot(norm, SA);
					if (nrom_sa > 0.0)
					{
						norm[0] = -norm[0];
						norm[1] = -norm[1];
						norm[2] = -norm[2];
					}

					double area = calArea(p1, p2, Nodes[newp].pt);
					double d = area / LenEdg * 2.0;
					double newpos[3] = {0}, alpha = 0.5;
					double *verts[4];

					std::vector<int> sph;
					findSphere(newp, sph);

					int iter = 0;
					while (iter++ < 16)
					{
						bool moveflag = true;
						for (int j = 0; j < 3; j++)
							newpos[j] = Nodes[newp].pt[j] + d * alpha * norm[j];
						alpha /= 2.0;
						for (int j = 0; j < sph.size(); j++)
						{
							if (ishulltet(sph[j]))
								continue;
							for (int m = 0; m < 4; m++)
							{
								int iElemNd = Elems[sph[j]].form[m];
								verts[m] = (iElemNd != newp) ? Nodes[iElemNd].pt : newpos;
							}
							double ori = dt::GEOM_FUNC::orient3d(verts[0], verts[1], verts[2], verts[3]);
							if (ori >= 0)
							{
								moveflag = false;
								break; // This tet becomes invalid.
							}
						}
						if (moveflag)
						{
							for (int j = 0; j < 3; j++)
								Nodes[newp].pt[j] = newpos[j];
							break;
						}
					} // while (iter < 16)
					break;
				} // for (int i = 0; i < shellp.size(); i++)
				if (recoverEdge(lostE, 0, 0) == 1)
				{
					// printf(" success:e\n");
					return 1;
				}
			} // if (isMeshEdge(it[0], it[1], &srctet))
		}
	} // for(auto it:Vid

	return FHCSteinerInsert(lostE, newN, Vid.size());
}

/* Boundary Steiner point removal*/
int DT::removeEdgStiner(const int idx, int level)
{
	// if (level > 3) return 0;
	int iNod = EdgSteiner[idx].first; // steiner point
	if (iNod == -1)
		return 1;
	int lostedg = EdgSteiner[idx].second; // lost edges
	int p1 = SurEdgs[lostedg].iStart;
	int p2 = SurEdgs[lostedg].iEnd;
	int subedg0 = SurEdgs[lostedg].info;	 // sub edge 0
	int subedg1 = SurEdgs[lostedg].info + 1; // sub edge 1
	int mainfold = SurEdgs[lostedg].face.size();
	int a, b, c, d, i, j, k, m;
	double edgelen = distance(Nodes[p1].pt, Nodes[p2].pt);
	std::vector<int> newN;
	std::vector<double> volvec;

	newN.push_back(iNod);
	// If this edge has been recovered due to accuracy issues
	int tempE, trytmvst = 0;
	if (isMeshEdge(p1, p2, &tempE))
	{
		std::vector<int> oldtet = {tempE};
		if (removeEdge(oldtet, isNod_in_Tet(p1, tempE), isNod_in_Tet(p2, tempE), 1000, -1) == 0)
		{
			trytmvst = 1;
		}
	}

	if (trytmvst == 1)
	{
		clearbndpnt(iNod);
		smooth_volume(iNod, true); // smooth_diff(iNod, 0);
	}
	else
	{
		if (mainfold == 1)
		{
			if (infolevel > 2)
				spdlog::info("Can't remove single edge steiner now!");
			return 2;
		}
		/************************** Build subTri Hash ************************/
		TriHasher<int> subTri;
		for (i = 0; i < SurEdgs[subedg0].face.size(); i++)
		{
			SurTri *s = &SurTris[SurEdgs[subedg0].face[i]];
			subTri.add(s->form[0], s->form[1], s->form[2], SurEdgs[subedg0].face[i]);
			if (!isMeshFace(s->form[0], s->form[1], s->form[2]))
			{
				int ret = recoverFace(SurEdgs[subedg0].face[i], 0);
			}
		}
		for (i = 0; i < SurEdgs[subedg1].face.size(); i++)
		{
			SurTri *s = &SurTris[SurEdgs[subedg1].face[i]];
			subTri.add(s->form[0], s->form[1], s->form[2], SurEdgs[subedg1].face[i]);
			if (!isMeshFace(s->form[0], s->form[1], s->form[2]))
			{
				int ret = recoverFace(SurEdgs[subedg1].face[i], 0);
			}
		}
		/**************************** color sphere **************************/
		std::vector<int> sph;
		findSphere(iNod, sph);

		std::vector<int> sphColor(sph.size(), -1);
		std::map<int, int> mp2sph;
		for (i = 0; i < sph.size(); i++)
			mp2sph[sph[i]] = i;
		int color = 0;
		for (i = 0; i < sph.size(); i++)
		{
			if (sphColor[i] != -1)
				continue;
			std::queue<int> q; // for BFS
			q.push(i);
			sphColor[i] = color;
			while (!q.empty())
			{
				int src = sph[q.front()];
				q.pop();
				for (m = 0; m < 4; m++)
				{
					int neig = getNeig(src, m);
					if (mp2sph.find(neig) == mp2sph.end())
						continue; // neig don't in sph
					if (sphColor[mp2sph[neig]] != -1)
						continue;
					DNC(m, a, b, c, d);
					b = Elems[src].form[b];
					c = Elems[src].form[c];
					d = Elems[src].form[d];
					if (subTri.find(b, c, d))
						continue;
					q.push(mp2sph[neig]);
					sphColor[mp2sph[neig]] = color;
				}
			}
			color++;
		}
		if (color != mainfold)
		{
			if (infolevel > 0)
				spdlog::info("Remove Edge steiner point {} error,classify sph fail!", idx);
			return 0;
		}

		// if (idx == 9)
		//	printSph_VTK(sph,"./temp.vtk");
		/*********************** determine norm and new pt ***********************/
		std::vector<std::vector<int>> sphcla(color);
		std::vector<std::vector<double>> normal(color, std::vector<double>(3));
		std::vector<std::vector<double>> newPcoord(color, std::vector<double>(3));
		std::vector<std::vector<int>> bndElem(mainfold * 2, std::vector<int>(3, -1));
		for (i = 0; i < mainfold; i++)
		{
			// continue hull tet,ignore volume
			bool hullflag = false;
			//------------- Classifiy sphere ------------
			for (j = 0; j < sph.size(); j++)
			{
				if (sphColor[j] == i)
				{
					sphcla[i].push_back(sph[j]);
					if (ishulltet(sph[j]))
						hullflag = true;
				}
			}

			//--------------- get sub Tri and bnd ---------------
			int nsubsph = sphcla[i].size();
			std::vector<std::vector<int>> bndpt(nsubsph + 2, std::vector<int>(3, -1));
			std::map<int, int> PTV; // parent bnd tri have been visited?
			int SubTriord = 0;
			for (k = 0; k < 3; k++)
				normal[i][k] = 0;
			for (j = 0; j < nsubsph; j++)
			{
				int t0 = sphcla[i][j];
				double vol = calVolume(t0);
				// printf("%d %e | ", j, vol);
				for (m = 0; m < 4; m++)
				{
					int neig = getNeig(t0, m);
					DFC(m, a, b, c, d);
					b = Elems[t0].form[b];
					c = Elems[t0].form[c];
					d = Elems[t0].form[d];

					if (mp2sph.find(neig) == mp2sph.end())
					{
						bndpt[j][0] = b; // store bnd point
						bndpt[j][1] = c;
						bndpt[j][2] = d;
						continue; // neig don't in sph
					}
					// In it's sphere
					if (sphColor[mp2sph[neig]] == i)
						continue;

					// determine ord of last two tri
					int tt = subTri.get(b, c, d);
					int Tid = SurTris[tt].parent; // parent tri's idx
					if (std::find(SurEdgs[lostedg].face.begin(), SurEdgs[lostedg].face.end(), Tid) == SurEdgs[lostedg].face.end())
					{
						continue;
					}
					if (PTV.find(Tid) == PTV.end())
					{
						PTV[Tid] = SubTriord++;
						// first bnd
						bndpt[nsubsph + PTV[Tid]][0] = b;
						bndpt[nsubsph + PTV[Tid]][1] = c;
						bndpt[nsubsph + PTV[Tid]][2] = d;
						for (k = 0; k < 3; k++)
						{
							if (bndpt[nsubsph + PTV[Tid]][k] == iNod)
							{
								for (int k1 = 0; k1 < 3; k1++)
								{
									if (SurTris[Tid].form[k1] != b && SurTris[Tid].form[k1] != c && SurTris[Tid].form[k1] != d)
									{
										bndpt[nsubsph + PTV[Tid]][k] = SurTris[Tid].form[k1];
										break;
									}
								}
								break;
							}
						}
						// normal [b,c,d]->a
						double f[3];
						calnormal(bndpt[nsubsph + PTV[Tid]][0], bndpt[nsubsph + PTV[Tid]][1], bndpt[nsubsph + PTV[Tid]][2], f);
						double flen = lenvec(f);
						for (k = 0; k < 3; k++)
							f[k] /= flen;
						for (k = 0; k < 3; k++)
							normal[i][k] += f[k];
					}
				}
			}
			// printf("\n");
			if (PTV.size() != 2)
			{
				return 0;
			}

			// store bnd infomation
			for (j = 0; j < 2; j++)
				for (k = 0; k < 3; k++)
					bndElem[i * 2 + j][k] = bndpt[nsubsph + j][k];

			// Normalization
			for (k = 0; k < 3; k++)
				normal[i][k] /= 2.0;

			//----------------- determine the new point ----------------
			double oldp[3] = {Nodes[iNod].pt[0], Nodes[iNod].pt[1], Nodes[iNod].pt[2]};
			double newp[3] = {oldp[0] + edgelen * normal[i][0], oldp[1] + edgelen * normal[i][1], oldp[2] + edgelen * normal[i][2]};
			double ori, len = edgelen / 10.0;
			double ShortestDistance = 1e-16;

			for (j = bndpt.size() - 1; j >= bndpt.size() - 2; j--)
			{
				ori = dt::GEOM_FUNC::orient3d(Nodes[bndpt[j][0]].pt, Nodes[bndpt[j][1]].pt, newp, Nodes[bndpt[j][2]].pt);
				if (ori <= 0)
				{
					std::swap(bndpt[j][1], bndpt[j][2]);
				}
			}

			// try to find a position,let all tet's volume is positive
			while (len > ShortestDistance)
			{
				bool allpositive = true;
				double newp[3] = {oldp[0] + len * normal[i][0], oldp[1] + len * normal[i][1], oldp[2] + len * normal[i][2]};
				for (j = bndpt.size() - 1; j >= 0; j--)
				{
					if (bndpt[j][0] == ghost || bndpt[j][1] == ghost || bndpt[j][2] == ghost)
						continue;
					ori = dt::GEOM_FUNC::orient3d(Nodes[bndpt[j][0]].pt, Nodes[bndpt[j][1]].pt, newp, Nodes[bndpt[j][2]].pt);
					if (ori <= 0)
					{
						// double vol = 0;
						// if (j < bndpt.size() - 2)
						//	vol = calVolume(sphcla[i][j]);
						// else {
						//	printf("0");
						// }
						// printf("%d %e | ", j, vol);

						// if (j < bndpt.size() - 2) {
						//	double G_pnt[3] = {0};
						//	for (int k = 0; k < 3; k++)
						//		for (int t = 0; t < 4; t++)
						//			G_pnt[k] += Nodes[Elems[sphcla[i][j]].form[t]].pt[k];
						//	double tempnormal[3] = { 0 };
						//	for (int k = 0; k < 3; k++)
						//		tempnormal[k] = G_pnt[k] / 4.0 - oldp[k];
						//	double nlen = lenvec(tempnormal);
						//	printf("%e %e %e ->", normal[i][0], normal[i][1], normal[i][2]);
						//	for (k = 0; k < 3; k++)
						//		normal[i][k] = tempnormal[k] / nlen;
						//	printf("%e %e %e\n", normal[i][0], normal[i][1], normal[i][2]);
						// }
						// if (j == bndpt.size() - 1 || j == bndpt.size() - 2) {
						//	len = ShortestDistance;
						// }
						len *= 0.5;
						allpositive = false;
						break;
					}
				}
				if (allpositive)
				{
					for (k = 0; k < 3; k++)
						newPcoord[i][k] = newp[k]; // store new position in normal for easily
					break;
				}
			}
			// printf("\n");
			// can't find a position let all tet's volume is positive
			if (len < ShortestDistance)
			{
				if (level < 3 && !hullflag)
				{
					int ret = -1;
					for (auto it : sph)
					{
						if (ishulltet(it) || isDelEle(it) || isNod_in_Tet(iNod, it) == -1)
							continue;
						double vol = calVolume(it);
						improve_step = true;
						improve_Metric = 3;
						ret = removebadtet(it, /* vol < 1e-10*/ 0, -1);
						improve_step = false;
					}
					return removeEdgStiner(idx, level + 1);
				}

				for (k = 0; k < 3; k++)
					newPcoord[i][k] = oldp[k]; // store new position in normal for easily
			}
		}

		/************************ creat new tet and neig info **************************/
		std::vector<int> newElemIdx(2 * mainfold);
		for (i = 0; i < mainfold; i++)
		{
			if (i == 0)
			{
				for (k = 0; k < 3; k++)
					Nodes[newN[i]].pt[k] = newPcoord[i][k]; // obtain new position
			} // iNod used for first class
			else
			{
				newN.push_back(addNode(newPcoord[i][0], newPcoord[i][1], newPcoord[i][2], Nodes[iNod].space));
				//---------- update elems connect to iNod -----------
				for (j = 0; j < sphcla[i].size(); j++)
				{
					int t0 = sphcla[i][j];
					for (m = 0; m < 4; m++)
					{
						if (Elems[t0].form[m] == iNod)
						{
							Elems[t0].form[m] = newN[i];
							break;
						}
					}
				}
			}
			clearbndpnt(newN[i]);
			// reset point to tet
			setP2T(newN[i], sphcla[i][0]);
			// create new tet
			for (j = 0; j < 2; j++)
			{
				newElemIdx[i * 2 + j] = addElem(bndElem[i * 2 + j][0], bndElem[i * 2 + j][1], bndElem[i * 2 + j][2], newN[i]);
			}

			// connect 2 neig information
			for (j = 0; j < 2; j++)
			{
				int t0 = newElemIdx[i * 2 + j];
				for (k = 0; k < 3; k++)
				{ // neig[0] connect to another class
					if (getNeig(t0, k) != -1)
						continue;
					DFC(k, a, b, c, d);
					b = Elems[t0].form[b];
					c = Elems[t0].form[c];
					d = Elems[t0].form[d];

					for (m = 0; m < sphcla[i].size(); m++)
					{
						int t1 = sphcla[i][m];
						int *form = Elems[t1].form;
						int ia = -1, ib = -1, ic = -1, id = -1;
						for (int k1 = 0; k1 < 4; k1++)
						{
							if (form[k1] == b)
								ib = k1;
							else if (form[k1] == c)
								ic = k1;
							else if (form[k1] == d)
								id = k1;
							else
								ia = k1;
						}
						if (ib == -1 || ic == -1 || id == -1)
							continue;
						// find neig
						bond(t0, k, t1, ia);
						break;
					}
				}
			}
		}
		// conect last one neig
		for (i = 0; i < mainfold * 2; i++)
		{
			int t0 = newElemIdx[i];
			for (k = 0; k < 4; k++)
			{ // neig[0] connect to another class
				if (getNeig(t0, k) != -1)
					continue;
				DNC(k, a, b, c, d);
				b = Elems[t0].form[b];
				c = Elems[t0].form[c];
				d = Elems[t0].form[d];
				for (j = i + 1; j < mainfold * 2; j++)
				{
					int t1 = newElemIdx[j];
					int *form = Elems[t1].form;
					int ia = -1, ib = -1, ic = -1, id = -1;
					for (int k1 = 0; k1 < 4; k1++)
					{
						if (form[k1] == b)
							ib = k1;
						else if (form[k1] == c)
							ic = k1;
						else if (form[k1] == d)
							id = k1;
						else
							ia = k1;
					}
					if (ib == -1 || ic == -1 || id == -1)
						continue;
					// find neig
					bond(t0, k, t1, ia);
					break;
				}
			}
		}
	}
	/************************ update SurTris and SurEdgs **************************/
	// SurEdgs
	setDelSurEdg(subedg0);	   // sub edge 0
	setDelSurEdg(subedg0 + 1); // sub edge 1
	setDelSurEdg(subedg0 + 2); // sub edge 2
	setDelSurEdg(subedg0 + 3); // sub edge 3
	SurEdgs[lostedg].info = 1; // recovered
	BndEdg.add(p1, p2, lostedg);

	BndEdg.erase(SurEdgs[subedg0 + 3].iStart, SurEdgs[subedg0 + 3].iEnd); // update BndEdg Hash
	BndEdg.erase(SurEdgs[subedg0 + 2].iStart, SurEdgs[subedg0 + 2].iEnd); // update BndEdg Hash
	BndEdg.erase(SurEdgs[subedg0 + 1].iStart, SurEdgs[subedg0 + 1].iEnd); // update BndEdg Hash
	BndEdg.erase(SurEdgs[subedg0].iStart, SurEdgs[subedg0].iEnd);		  // update BndEdg Hash

	for (i = 0; i < SurEdgs[subedg0].face.size(); i++)
	{
		SurTri *s = &SurTris[SurEdgs[subedg0].face[i]]; // sub Tri
		setDelSurTri(SurEdgs[subedg0].face[i]);
		// SurTris[s->parent].info = 1;                      //set parent bnd Tri is recovered, Comment for findTriParent
		BndTri.erase(s->form[0], s->form[1], s->form[2]); // update BndTri Hash
	}
	for (i = 0; i < SurEdgs[subedg1].face.size(); i++)
	{
		SurTri *s = &SurTris[SurEdgs[subedg1].face[i]]; // sub Tri
		setDelSurTri(SurEdgs[subedg1].face[i]);
		// SurTris[s->parent].info = 1;                      //set parent bnd Tri is recovered
		BndTri.erase(s->form[0], s->form[1], s->form[2]); // update BndTri Hash
	}

	for (i = 0; i < newN.size(); i++)
	{
		if (!removePnt(newN[i]))
		{
			smooth_volume(newN[i], true);
		}
	}

	// Recover BndTri
	for (i = 0; i < mainfold; i++)
	{
		int Fidx = SurEdgs[lostedg].face[i];
		BndTri.add(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2], Fidx);
		if (!isMeshFace(SurTris[Fidx].form[0], SurTris[Fidx].form[1], SurTris[Fidx].form[2]))
		{
			SurTris[Fidx].info = -1;
			int ret = recoverFace(Fidx, 0);
			if (ret < 0)
				spdlog::warn("May Error happen when remove steiner point in edge {}", idx);
		}
		SurTris[Fidx].info = 1;
		for (int j = 0; j < 3; j++)
		{
			int tempp1 = SurTris[Fidx].form[j];
			int tempp2 = SurTris[Fidx].form[(j + 1) % 3];
			int tempEdg = BndEdg.get(tempp1, tempp2);
			if (tempEdg == -1)
				continue;
			for (int k = 0; k < SurEdgs[tempEdg].face.size(); k++)
			{
				int subsubsubf = SurEdgs[tempEdg].face[k];
				for (int t = 0; t < 3; t++)
				{
					if (SurTris[subsubsubf].form[t] == iNod)
					{
						SurEdgs[tempEdg].face[k] = Fidx;
						break;
					}
				}
			}
		}
	}

	return 1;
}
