#include "hierarchy.h"
#include "positions.h"
#include "timer.h"
#include"quadratic.h"
//3D===========================================================================================================//
std::vector<std::vector<uint32_t>> mTs;

std::vector<tuple_E> mpEs;
std::vector<std::vector<uint32_t>> mpFvs, mpFes, mpPs;
std::vector<bool>mpF_boundary_flag; std::vector<std::vector<uint32_t>> PV_npfs, PE_npfs, PF_npps;

std::vector<short> mpV_flag, mpE_flag, mpF_flag, mpP_flag;

std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan> Es_red;
std::vector<uint32_t> pV_map;
std::vector<std::vector<uint32_t>> Reverse_pV_map;
std::vector<uint32_t> pE_map;
std::vector<std::vector<uint32_t>> Reverse_pE_map;
std::vector<uint32_t> pF_map;
std::vector<std::vector<uint32_t>> Reverse_pF_map;


std::vector<std::vector<uint32_t>> PV_npvs;
std::vector<std::vector<uint32_t>> PV_npes_sudo;
std::vector<std::vector<uint32_t>> PV_npfs_sudo;
std::vector<std::vector<uint32_t>> PE_npfs_sudo;
std::vector<std::vector<uint32_t>> PF_npps_sudo;

//re-coloring
MatrixXf mQ_copy, mO_copy, mN_copy;
vector<Quadric> Quadric_copy, newQu3D;
MatrixXf newQ, newN3D, newV3D, newC3D;

//timeing
long long topo_check_time = 0, decomposition_time = 0, total_time = 0;
//===========================================================================================================//

void construct_Es_TetEs_Fs_TetFs_FEs()
{
	mpEs.clear(); mpFvs.clear(); mpFes.clear(); mpPs.clear(); mpF_boundary_flag.clear();
	//mpFvs, mpPs, mpF_boundary_flag
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> tempF;
	tempF.reserve(mTs.size() * 4);
	mpPs.resize(mTs.size());
	for (uint32_t t = 0; t < mTs.size(); ++t) {
		for (uint32_t f = 0; f < 4; ++f) {
			uint32_t v0 = mTs[t][tet_faces[f][0]], v1 = mTs[t][tet_faces[f][1]], v2 = mTs[t][tet_faces[f][2]];
			if (v0 > v1) std::swap(v0, v1);
			if (v1 > v2) std::swap(v2, v1);
			if (v0 > v1) std::swap(v0, v1);
			tempF.push_back(std::make_tuple(v0, v1, v2, t, f));
		}
		std::vector<uint32_t> fs(4);
		mpPs[t] = fs;
	}
	std::sort(tempF.begin(), tempF.end());
	mpFvs.reserve(tempF.size() / 3); mpF_boundary_flag.reserve(tempF.size() / 3);
	int F_num = -1;
	std::vector<uint32_t> fi(3);
	for (uint32_t i = 0; i < tempF.size(); ++i) {
		if (i == 0 || (i != 0 &&
			(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) ||
				std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
				std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1])))) {
			F_num++;
			fi[0] = std::get<0>(tempF[i]); fi[1] = std::get<1>(tempF[i]); fi[2] = std::get<2>(tempF[i]);
			mpFvs.push_back(fi); mpF_boundary_flag.push_back(true);
		}
		else if (i != 0 && (std::get<0>(tempF[i]) == std::get<0>(tempF[i - 1]) &&
			std::get<1>(tempF[i]) == std::get<1>(tempF[i - 1]) &&
			std::get<2>(tempF[i]) == std::get<2>(tempF[i - 1])))
			mpF_boundary_flag[F_num] = false;

		mpPs[std::get<3>(tempF[i])][std::get<4>(tempF[i])] = F_num;
	}
	//mpFes, mEs
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> temp;
	temp.reserve(mpFvs.size() * 3);
	mpFes.resize(mpFvs.size()); std::vector<uint32_t> fes(3);
	for (uint32_t i = 0; i < mpFvs.size(); ++i) {
		for (uint32_t e = 0; e < 3; ++e) {
			uint32_t v0 = mpFvs[i][e], v1 = mpFvs[i][(e + 1) % 3];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, i, e));
		}
		mpFes[i] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mpEs.reserve(temp.size() / 3);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mpEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), false, 0, Edge_tag::B, E_num, -1, 0));
		}
		mpFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
	//
	for (uint32_t i = 0; i < mpFes.size(); ++i)
			if (mpF_boundary_flag[i]) for (uint32_t e = 0; e < 3; ++e) std::get<2>(mpEs[mpFes[i][e]]) = true;


	PE_npfs.clear(); PE_npfs.resize(mpEs.size());
	PF_npps.clear(); PF_npps.resize(mpFvs.size());
	PV_npvs.clear(); PV_npvs.resize(mO_copy.cols());
	for (uint32_t i = 0; i < mpEs.size(); i++) {
		uint32_t v0 = get<0>(mpEs[i]), v1 = get<1>(mpEs[i]);
		PV_npvs[v0].push_back(v1);
		PV_npvs[v1].push_back(v0);
	}
	for (uint32_t i = 0; i < mpFes.size(); i++) for (auto eid:mpFes[i]) { PE_npfs[eid].push_back(i); }
	for (uint32_t i = 0; i < mpPs.size(); i++) for (auto fid:mpPs[i]) PF_npps[fid].push_back(i);
}

bool simple_polygon_3D(std::vector<std::vector<uint32_t>> &fvs, std::vector<std::vector<uint32_t>> &fes, std::vector<uint32_t> &pvs,
	std::vector<uint32_t> &pes, std::vector<uint32_t> &vs_disgard, std::vector<uint32_t> &es_disgard, bool v_involve)
{
	es_disgard.clear();//es_disgard is in the interior
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				es_disgard.push_back(fes[i][j]);
				mpE_flag[fes[i][j]] = false;
			}
			else mpE_flag[fes[i][j]] = true;
	}
	short which_polygon = 0;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				if (!pes.size())
					which_polygon = i;

				pes.push_back(fes[i][j]); mpE_flag[fes[i][j]] = false;
			}
	}
	//test nvs for each v
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = pV_map[std::get<0>(mpEs[pes[i]])], v1 = pV_map[std::get<1>(mpEs[pes[i]])];
		mpV_flag[v0]++; mpV_flag[v1]++;
		if (mpV_flag[v0] > 2 || mpV_flag[v1] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				uint32_t v0_ = pV_map[std::get<0>(mpEs[pes[j]])], v1_ = pV_map[std::get<1>(mpEs[pes[j]])];
				mpV_flag[v0_] = mpV_flag[v1_] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = pV_map[std::get<0>(mpEs[pes[i]])], v1 = pV_map[std::get<1>(mpEs[pes[i]])];
		mpV_flag[v0] = mpV_flag[v1] = 0;
	}
	//extract the polygon	
	if (!pes.size()) return false;
	pvs.clear();
	pvs.reserve(pes.size());
	std::vector<bool> e_flag(pes.size(), false);
	pvs.push_back(pV_map[std::get<0>(mpEs[pes[0]])]);
	pvs.push_back(pV_map[std::get<1>(mpEs[pes[0]])]);
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {
				if (pV_map[std::get<0>(mpEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(pV_map[std::get<1>(mpEs[pes[j]])]);
					start_v = pV_map[std::get<1>(mpEs[pes[j]])];
					break;
				}
				else if (pV_map[std::get<1>(mpEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(pV_map[std::get<0>(mpEs[pes[j]])]);
					start_v = pV_map[std::get<0>(mpEs[pes[j]])];
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size())
		return false;
	if (!v_involve) return true;
	//judge direction
	bool correct = true;
	for (int i = 0; i < fvs[which_polygon].size(); i++)
		if (fvs[which_polygon][i] == pV_map[std::get<0>(mpEs[pes[0]])])
			if (fvs[which_polygon][(i + 1) % fvs[which_polygon].size()] != pV_map[std::get<1>(mpEs[pes[0]])])
				correct = false;
	if (!correct)
		std::reverse(pvs.begin(), pvs.end());
	//vs_disgard
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			mpV_flag[fvs[i][j]] = 1;
	for (uint32_t i = 0; i < pvs.size(); ++i) mpV_flag[pvs[i]] = 0;
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			if (mpV_flag[fvs[i][j]] == 1) {
				vs_disgard.push_back(fvs[i][j]); mpV_flag[fvs[i][j]] = 0;
			}
	return true;
}
bool simple_polygon_3D_v2(vector<vector<uint32_t>> &fvs, vector<vector<uint32_t>> &fes, vector<uint32_t> &pvs,
	vector<uint32_t> &pes, vector<uint32_t> &vs_disgard, vector<uint32_t> &es_disgard, bool v_involve){
	es_disgard.clear();//es_disgard is in the interior
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				es_disgard.push_back(fes[i][j]);
				mpE_flag[fes[i][j]] = false;
			}
			else mpE_flag[fes[i][j]] = true;
	}
	short which_polygon = 0;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mpE_flag[fes[i][j]]) {
				if (!pes.size())
					which_polygon = i;

				pes.push_back(fes[i][j]); mpE_flag[fes[i][j]] = false;
			}
	}
	//test nvs for each v
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = get<0>(mpEs[pes[i]]), v1 = get<1>(mpEs[pes[i]]);
		mpV_flag[v0]++; mpV_flag[v1]++;
		if (mpV_flag[v0] > 2 || mpV_flag[v1] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				uint32_t v0_ = get<0>(mpEs[pes[j]]), v1_ = get<1>(mpEs[pes[j]]);
				mpV_flag[v0_] = mpV_flag[v1_] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = get<0>(mpEs[pes[i]]), v1 = get<1>(mpEs[pes[i]]);
		mpV_flag[v0] = mpV_flag[v1] = 0;
	}
	//extract the polygon	
	if (!pes.size()) return false;
	pvs.clear();
	pvs.reserve(pes.size());
	std::vector<bool> e_flag(pes.size(), false);
	pvs.push_back(get<0>(mpEs[pes[0]]));
	pvs.push_back(get<1>(mpEs[pes[0]]));
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {
				if (get<0>(mpEs[pes[j]]) == start_v) {
					e_flag[j] = true;
					pvs.push_back(get<1>(mpEs[pes[j]]));
					start_v = get<1>(mpEs[pes[j]]);
					break;
				}
				else if (get<1>(mpEs[pes[j]]) == start_v) {
					e_flag[j] = true;
					pvs.push_back(get<0>(mpEs[pes[j]]));
					start_v = get<0>(mpEs[pes[j]]);
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size())
		return false;
	if (!v_involve) return true;
	//judge direction
	bool correct = true;
	for (int i = 0; i < fvs[which_polygon].size(); i++)
		if (fvs[which_polygon][i] == get<0>(mpEs[pes[0]]))
			if (fvs[which_polygon][(i + 1) % fvs[which_polygon].size()] != get<1>(mpEs[pes[0]]))
				correct = false;
	if (!correct)
		std::reverse(pvs.begin(), pvs.end());
	//vs_disgard
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			mpV_flag[fvs[i][j]] = 1;
	for (uint32_t i = 0; i < pvs.size(); ++i) mpV_flag[pvs[i]] = 0;
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			if (mpV_flag[fvs[i][j]] == 1) {
				vs_disgard.push_back(fvs[i][j]); mpV_flag[fvs[i][j]] = 0;
			}
	return true;
}
bool simple_polyhedral(std::vector<std::vector<uint32_t>> &pfs, std::vector<uint32_t> &pf,
	std::vector<uint32_t> &vs_disgard, std::vector<uint32_t> &es_disgard, std::vector<uint32_t> &fs_disgard)
{
	fs_disgard.clear();//fs_disgard is in the interior
	for (int i = 0; i < pfs.size(); i++) {
		for (int j = 0; j < pfs[i].size(); j++)
			if (mpF_flag[pfs[i][j]]) {
				fs_disgard.push_back(pfs[i][j]);
				mpF_flag[pfs[i][j]] = false;
			}
			else mpF_flag[pfs[i][j]] = true;
	}
	for (int i = 0; i < pfs.size(); i++) {
		for (int j = 0; j < pfs[i].size(); j++)
			if (mpF_flag[pfs[i][j]]) {
				pf.push_back(pfs[i][j]); mpF_flag[pfs[i][j]] = false;
			}
	}
	//test each e whether non-manifold
	bool non_simple = false;
	for (uint32_t i = 0; i < pf.size(); ++i) {
		for (uint32_t j = 0; j < mpFes[pf[i]].size(); ++j) {
			mpE_flag[mpFes[pf[i]][j]]++;
			if (mpE_flag[mpFes[pf[i]][j]] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pf.size(); ++k) for (uint32_t j = 0; j < mpFes[pf[k]].size(); ++j) mpE_flag[mpFes[pf[k]][j]] = false;
			return false;
		}
	}
	for (uint32_t k = 0; k < pf.size(); ++k) for (uint32_t j = 0; j < mpFes[pf[k]].size(); ++j) if (mpE_flag[mpFes[pf[k]][j]] != 2) non_simple = true;
	for (uint32_t i = 0; i < pf.size(); ++i) for (uint32_t j = 0; j < mpFes[pf[i]].size(); ++j) mpE_flag[mpFes[pf[i]][j]] = false;
	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pf.size(); ++i) {
		for (uint32_t j = 0; j <mpFes[pf[i]].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[mpFes[pf[i]][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[mpFes[pf[i]][j]])];
			PV_npfs_sudo[v0].push_back(i);
			PV_npfs_sudo[v1].push_back(i);
			if (!mpV_flag[v0]) { vs_set.push_back(v0); mpV_flag[v0] = true; }
			if (!mpV_flag[v1]) { vs_set.push_back(v1); mpV_flag[v1] = true; }
		}
	}
	for (uint32_t i = 0; i < pf.size(); ++i) {
		for (uint32_t j = 0; j <mpFes[pf[i]].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[mpFes[pf[i]][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[mpFes[pf[i]][j]])];
			mpV_flag[v0] = mpV_flag[v1] = false;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++) {
		std::sort(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end());
		PV_npfs_sudo[vs_set[m]].erase(std::unique(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()), PV_npfs_sudo[vs_set[m]].end());
		if (PV_npfs_sudo[vs_set[m]].size() == pf.size()) continue;

		std::vector<std::vector<uint32_t>> fes(PV_npfs_sudo[vs_set[m]].size()), fvs;
		for (uint32_t k = 0; k < fes.size(); k++) fes[k] = mpFes[pf[PV_npfs_sudo[vs_set[m]][k]]];
		std::vector<uint32_t> fv, fe, vs_dis, es_dis;
		if (!simple_polygon_3D(fvs, fes, fv, fe, vs_dis, es_dis, false))
		{
			non_simple = true; break;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++)
		PV_npfs_sudo[vs_set[m]].clear();
	if (non_simple) return false;
	//test single layer of polyhedral, non non-manifold v and multi-layers of polyhedral
	std::vector<uint32_t> pf_temp, pf_sudo; pf_sudo.reserve(pf.size()); pf_temp.reserve(pf.size());
	for (uint32_t i = 0; i < pf.size(); ++i) mpF_flag[pf[i]] = true;
	pf_temp.push_back(pf[0]); pf_sudo = pf_temp; mpF_flag[pf_temp[0]] = false;
	while (pf_temp.size()) {
		std::vector<uint32_t> pf_;
		for (int i = 0; i < pf_temp.size(); i++)
			for (int j = 0; j < mpFes[pf_temp[i]].size(); j++) {
				uint32_t eid = mpFes[pf_temp[i]][j];
				for (int k = 0; k<PE_npfs[eid].size(); k++)
					if (mpF_flag[PE_npfs[eid][k]]) {
						pf_.push_back(PE_npfs[eid][k]);
						mpF_flag[PE_npfs[eid][k]] = false;
					}
			}
		if (pf_.size()) {
			pf_temp = pf_;
			pf_sudo.insert(pf_sudo.end(), pf_temp.begin(), pf_temp.end());
		}
		else break;
	}

	if (pf_sudo.size() != pf.size()) {
		for (uint32_t i = 0; i < pf.size(); ++i) mpF_flag[pf[i]] = false;
		return false;
	}
	//es_disgard, vs_disgard
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			for (uint32_t k = 0; k < mpFes[pfs[i][j]].size(); k++)
				mpE_flag[mpFes[pfs[i][j]][k]] = true;
			for (uint32_t k = 0; k < mpFvs[pfs[i][j]].size(); k++)
				mpV_flag[mpFvs[pfs[i][j]][k]] = true;
		}
	}
	for (uint32_t j = 0; j < pf.size(); ++j) {
		for (uint32_t k = 0; k < mpFes[pf[j]].size(); k++)
			mpE_flag[mpFes[pf[j]][k]] = false;
		for (uint32_t k = 0; k < mpFvs[pf[j]].size(); k++)
			mpV_flag[mpFvs[pf[j]][k]] = false;
	}
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			for (uint32_t k = 0; k < mpFes[pfs[i][j]].size(); k++)
				if (mpE_flag[mpFes[pfs[i][j]][k]]) {
					mpE_flag[mpFes[pfs[i][j]][k]] = false;
					es_disgard.push_back(mpFes[pfs[i][j]][k]);
				}
			for (uint32_t k = 0; k < mpFvs[pfs[i][j]].size(); k++)
				if (mpV_flag[mpFvs[pfs[i][j]][k]]) {
					mpV_flag[mpFvs[pfs[i][j]][k]] = false;
					vs_disgard.push_back(mpFvs[pfs[i][j]][k]);
				}
		}
	}

	return true;
}
bool simple_polyhedral_v2(std::vector<std::vector<uint32_t>> &pfs)
{
	//test each e whether non-manifold
	bool non_simple = false;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			mpE_flag[pfs[i][j]]++;
			if (mpE_flag[pfs[i][j]] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pfs.size(); ++k) for (uint32_t j = 0; j < pfs[k].size(); ++j) mpE_flag[pfs[k][j]] = false;
			return false;
		}
	}
	for (uint32_t k = 0; k < pfs.size(); ++k) 
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			if (mpE_flag[pfs[k][j]] != 2) non_simple = true;
	for (uint32_t k = 0; k < pfs.size(); ++k)
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			mpE_flag[pfs[k][j]] = false;

	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[pfs[i][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[pfs[i][j]])];
			PV_npfs_sudo[v0].push_back(i);
			PV_npfs_sudo[v1].push_back(i);
			if (!mpV_flag[v0]) { vs_set.push_back(v0); mpV_flag[v0] = true; }
			if (!mpV_flag[v1]) { vs_set.push_back(v1); mpV_flag[v1] = true; }
		}
	}
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = pV_map[std::get<0>(mpEs[pfs[i][j]])];
			uint32_t v1 = pV_map[std::get<1>(mpEs[pfs[i][j]])];
			mpV_flag[v0] = mpV_flag[v1] = false;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++){
		std::sort(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()); 
		PV_npfs_sudo[vs_set[m]].erase(std::unique(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()), PV_npfs_sudo[vs_set[m]].end());
		if (PV_npfs_sudo[vs_set[m]].size() == pfs.size()) continue;

		std::vector<std::vector<uint32_t>> fes(PV_npfs_sudo[vs_set[m]].size()), fvs;
		for (uint32_t k = 0; k < fes.size(); k++) fes[k] = pfs[PV_npfs_sudo[vs_set[m]][k]];	
		std::vector<uint32_t> fv, fe, vs_dis, es_dis;
		if (!simple_polygon_3D(fvs, fes, fv, fe, vs_dis, es_dis, false))
		{
			non_simple = true; break;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++)
		PV_npfs_sudo[vs_set[m]].clear();
	if (non_simple) return false;

	return true;
}
bool simple_polyhedral_v3(vector<vector<uint32_t>> &pfs)
{
	//test each e whether non-manifold
	bool non_simple = false;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			mpE_flag[pfs[i][j]]++;
			if (mpE_flag[pfs[i][j]] > 2) non_simple = true;
		}
		if (non_simple) {
			for (uint32_t k = 0; k < pfs.size(); ++k) for (uint32_t j = 0; j < pfs[k].size(); ++j) mpE_flag[pfs[k][j]] = false;
			return false;
		}
	}
	for (uint32_t k = 0; k < pfs.size(); ++k)
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			if (mpE_flag[pfs[k][j]] != 2) non_simple = true;
	for (uint32_t k = 0; k < pfs.size(); ++k)
		for (uint32_t j = 0; j < pfs[k].size(); ++j)
			mpE_flag[pfs[k][j]] = false;

	if (non_simple) return false;
	//test each v whether non-manifold
	std::vector<uint32_t> vs_set;
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = get<0>(mpEs[pfs[i][j]]);
			uint32_t v1 = get<1>(mpEs[pfs[i][j]]);
			PV_npfs_sudo[v0].push_back(i);
			PV_npfs_sudo[v1].push_back(i);
			if (!mpV_flag[v0]) { vs_set.push_back(v0); mpV_flag[v0] = true; }
			if (!mpV_flag[v1]) { vs_set.push_back(v1); mpV_flag[v1] = true; }
		}
	}
	for (uint32_t i = 0; i < pfs.size(); ++i) {
		for (uint32_t j = 0; j < pfs[i].size(); ++j) {
			uint32_t v0 = get<0>(mpEs[pfs[i][j]]);
			uint32_t v1 = get<1>(mpEs[pfs[i][j]]);
			mpV_flag[v0] = mpV_flag[v1] = false;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++) {
		std::sort(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end());
		PV_npfs_sudo[vs_set[m]].erase(std::unique(PV_npfs_sudo[vs_set[m]].begin(), PV_npfs_sudo[vs_set[m]].end()), PV_npfs_sudo[vs_set[m]].end());
		if (PV_npfs_sudo[vs_set[m]].size() == pfs.size()) continue;

		std::vector<std::vector<uint32_t>> fes(PV_npfs_sudo[vs_set[m]].size()), fvs;
		for (uint32_t k = 0; k < fes.size(); k++) fes[k] = pfs[PV_npfs_sudo[vs_set[m]][k]];
		std::vector<uint32_t> fv, fe, vs_dis, es_dis;
		if (!simple_polygon_3D_v2(fvs, fes, fv, fe, vs_dis, es_dis, false))
		{
			non_simple = true; break;
		}
	}
	for (uint32_t m = 0; m<vs_set.size(); m++)
		PV_npfs_sudo[vs_set[m]].clear();
	if (non_simple) return false;

	return true;
}
void cut_a_polyhedral(std::vector<uint32_t> &ps, std::vector<std::vector<uint32_t>> &pfs, std::vector<uint32_t> &e_circle, 
	std::vector<uint32_t> &ps0, std::vector<uint32_t> &ps1)
{
	for (uint32_t i = 0; i < pfs.size(); ++i) for (uint32_t j = 0; j < pfs[i].size(); ++j) PE_npfs_sudo[pfs[i][j]].push_back(i);
	for (uint32_t i = 0; i < ps.size(); ++i) mpF_flag[i] = true;

	std::vector<uint32_t> pf_temp; pf_temp.push_back(0); ps0 = pf_temp; mpF_flag[0] = false;
	while (pf_temp.size()) {
		std::vector<uint32_t> pf_;
		for (int i = 0; i < pf_temp.size(); i++)
			for (int j = 0; j < pfs[pf_temp[i]].size(); j++) {
				uint32_t eid = pfs[pf_temp[i]][j];
				if (std::find(e_circle.begin(), e_circle.end(), eid) != e_circle.end()) continue;
				for (int k = 0; k<PE_npfs_sudo[eid].size(); k++)
					if (mpF_flag[PE_npfs_sudo[eid][k]]) {
						pf_.push_back(PE_npfs_sudo[eid][k]);
						mpF_flag[PE_npfs_sudo[eid][k]] = false;
					}
			}
		if (pf_.size()) {
			pf_temp = pf_;
			ps0.insert(ps0.end(), pf_temp.begin(), pf_temp.end());
		}
		else break;
	}
	for (uint32_t i = 0; i < ps.size(); ++i) mpF_flag[i] = false;
	for (uint32_t i = 0; i < pfs.size(); ++i) for (uint32_t j = 0; j < pfs[i].size(); ++j) if(PE_npfs_sudo[pfs[i][j]].size()) PE_npfs_sudo[pfs[i][j]].clear();
	for (uint32_t i = 0; i < ps0.size(); ++i) ps0[i] = ps[ps0[i]];
	std::set<uint32_t> s_model(ps.begin(), ps.end());
	std::set<uint32_t> s_pattern(ps0.begin(), ps0.end());
	std::set_difference(s_model.begin(), s_model.end(), s_pattern.begin(), s_pattern.end(), std::back_inserter(ps1));
}
void reindex_3D(MatrixXf &HV, std::vector<std::vector<uint32_t>> &HFv, std::vector<std::vector<uint32_t>> &HPf){
	//re-index F
	std::vector<int32_t> F_flag(HFv.size(), -1);
	for (auto pfs : HPf) for (auto fid : pfs) F_flag[fid] = 0;
	uint32_t f_num = 0;
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) F_flag[i] = f_num++;
	std::vector<std::vector<uint32_t>> mFs_local_(f_num);
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) mFs_local_[F_flag[i]] = HFv[i];
	mFs_local_.swap(HFv);
	for (auto &pfs : HPf)for (uint32_t i = 0; i < pfs.size(); i++)pfs[i] = F_flag[pfs[i]];
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs:HFv) for (auto vid:fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) mV_local_.col(V_flag[i]) = HV.col(i);
	for (auto &fvs:HFv) for (uint32_t j = 0; j < fvs.size(); j++)
			fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
}
void reindex_3D(MatrixXf &HV, MatrixXf &HQ, std::vector<std::vector<uint32_t>> &HFv, std::vector<std::vector<uint32_t>> &HPf) {
	//re-index F
	std::vector<int32_t> F_flag(HFv.size(), -1);
	for (auto pfs : HPf) for (auto fid : pfs) F_flag[fid] = 0;
	uint32_t f_num = 0;
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) F_flag[i] = f_num++;
	std::vector<std::vector<uint32_t>> mFs_local_(f_num);
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) mFs_local_[F_flag[i]] = HFv[i];
	mFs_local_.swap(HFv);
	for (auto &pfs : HPf)for (uint32_t i = 0; i < pfs.size(); i++)pfs[i] = F_flag[pfs[i]];
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs : HFv) for (auto vid : fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num), mQ_local_(4, v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) {
			mV_local_.col(V_flag[i]) = HV.col(i);
			mQ_local_.col(V_flag[i]) = HQ.col(i);
		}
	for (auto &fvs : HFv) for (uint32_t j = 0; j < fvs.size(); j++)
		fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
	mQ_local_.swap(HQ);
}
void reindex_3D(MatrixXf &HV, MatrixXf &HQ, vector<Quadric> &HQU, std::vector<std::vector<uint32_t>> &HFv, std::vector<std::vector<uint32_t>> &HPf) {
	//re-index F
	std::vector<int32_t> F_flag(HFv.size(), -1);
	for (auto pfs : HPf) for (auto fid : pfs) F_flag[fid] = 0;
	uint32_t f_num = 0;
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) F_flag[i] = f_num++;
	std::vector<std::vector<uint32_t>> mFs_local_(f_num);
	for (uint32_t i = 0; i < F_flag.size(); i++)
		if (F_flag[i] != -1) mFs_local_[F_flag[i]] = HFv[i];
	mFs_local_.swap(HFv);
	for (auto &pfs : HPf)for (uint32_t i = 0; i < pfs.size(); i++)pfs[i] = F_flag[pfs[i]];
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs : HFv) for (auto vid : fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num), mQ_local_(4, v_num);
	vector<Quadric> HQU_local_(v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) {
			mV_local_.col(V_flag[i]) = HV.col(i);
			mQ_local_.col(V_flag[i]) = HQ.col(i);
			HQU_local_[V_flag[i]] = HQU[i];
		}
	for (auto &fvs : HFv) for (uint32_t j = 0; j < fvs.size(); j++)
		fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
	mQ_local_.swap(HQ);
	HQU_local_.swap(HQU);
}

bool MultiResolutionHierarchy::recursive_ring(vector<uint32_t> &rvs, vector<uint32_t> &pvs, vector<uint32_t> &pes, uint32_t v0, uint32_t ve[2]) {
	if (v0 == ve[0] || v0 == ve[1]) return true;

	std::vector<std::tuple<Float, uint32_t, uint32_t>> potentials;
	for (auto eid : PV_npes_sudo[v0]) {

		if (find(pes.begin(), pes.end(), eid) != pes.end()) continue;

		int32_t fid = share_the_same_fs(pes, eid);
		if (fid != -1) continue; //share boundary triangles/quads/pentagons

		uint32_t vid = get<0>(mpEs[eid]);
		if (get<0>(mpEs[eid]) == v0) vid = get<1>(mpEs[eid]);

		Float cost_ = compute_cost_face3D(rvs, vid, true);
		potentials.push_back(make_tuple(cost_, eid, vid));
	}
	sort(potentials.begin(), potentials.end());
	if (!potentials.size()) return false;
	for (auto a_tuple : potentials) {
		pvs.push_back(get<2>(a_tuple));
		pes.push_back(get<1>(a_tuple));
		rvs.push_back(get<2>(a_tuple));

		if (!recursive_ring(rvs, pvs, pes, get<2>(a_tuple), ve)) {
			pvs.pop_back();
			pes.pop_back();
			rvs.pop_back();
			continue;
		}
		return true;
	}
	return false;
};
int32_t MultiResolutionHierarchy::share_the_same_fs(vector<uint32_t> &es, const uint32_t eid) {
	for (uint32_t i = 0; i < es.size(); i++) {
		std::vector<uint32_t> sharedf;
		std::set_intersection(PE_npfs[es[i]].begin(), PE_npfs[es[i]].end(), PE_npfs[eid].begin(), PE_npfs[eid].end(), back_inserter(sharedf));
		if (sharedf.size()) {
			if (mpFes[sharedf[0]].size() <= 4) return sharedf[0]; //share boundary triangles/quads
			uint32_t v0_e0 = get<0>(mpEs[eid]), v1_e0 = get<1>(mpEs[eid]),
				v0_ei = get<0>(mpEs[es[i]]), v1_ei = get<1>(mpEs[es[i]]);
			if (!(v0_e0 == v0_ei || v1_e0 == v0_ei || v0_e0 == v1_ei || v1_e0 == v1_ei)) return sharedf[0];//isolated two edges
																										   //not in the same direction
		}
	}
	return -1;
}
Float MultiResolutionHierarchy::compute_cost_face3D(vector<uint32_t> &vs, uint32_t rv, bool single) {
	vector<uint32_t> vs0; vector<Quaternion> qs(vs.size());
	Float min_cost = numeric_limits<Float>::max();
	vector<Vector3f> dirs;
	if (single) {
		vs0.push_back(rv);
		Quaternion q0 = mQ_copy.col(rv);
		for (uint32_t i = 0; i < vs.size(); i++) qs[i] = Quaternion::applyRotation(mQ_copy.col(vs[i]), q0);

		dirs.resize(vs.size());
		for (uint32_t i = 0; i < vs.size(); i++)
			dirs[i] = exact_3dir(mO_copy.col(rv), q0, mO_copy.col(vs[i]), qs[i], mScale, mInvScale);
	}
	else {
		vs0 = vs;
		Quaternion q0 = mQ_copy.col(vs[0]); qs[0] = q0;
		for (uint32_t i = 1; i < vs.size(); i++) qs[i] = Quaternion::applyRotation(mQ_copy.col(vs[i]), q0);

		vector<Vector3f> dirs;
		for (uint32_t i = 0; i<vs0.size(); i++)for (uint32_t j = i + 1; j<vs.size(); j++)
			dirs.push_back(exact_3dir(mO_copy.col(vs0[i]), qs[i], mO_copy.col(vs[j]), qs[j], mScale, mInvScale));
	}
	for (uint32_t i = 0; i < 3; i++) {
		Float cost_i = 0;
		for (auto dir : dirs) {
			dir[i] = std::abs(dir[i]);
			cost_i += dir[i];
		}
		if (cost_i < min_cost) min_cost = cost_i;
	}
	return min_cost;
}
