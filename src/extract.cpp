#include "hierarchy.h"
#include "positions.h"
//2D&3D========================================================================================================//
std::priority_queue<tuple_E, std::vector<tuple_E>, LessThan> Es_reddash;
std::vector<uint32_t> V_map;
std::vector<std::vector<uint32_t>> Reverse_V_map;

std::vector<tuple_E> mEs;
std::vector<tuple_F> mFs;
std::vector<std::vector<uint32_t>> mFs2D;
std::vector<std::vector<uint32_t>> FEs;

vector<vector<uint32_t>> V_pvs, V_pes, E_pfs;

std::vector<short> mV_flag;
std::vector<bool> mE_flag, mF_flag;
//re-coloring
MatrixXf mQ_copy2D, mO_copy2D, mN_copy2D, mV_copy2D;
MatrixXf newQ2D, newN2D, newV2D;
bool non_manifold = false;
//2D===========================================================================================================//
void construct_Es_FEs()
{
	FEs.clear(); mEs.clear(); V_pvs.clear(); E_pfs.clear();

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(mFs2D.size() * 3);
	FEs.resize(mFs2D.size());
	for (uint32_t f = 0; f < mFs2D.size(); ++f) {
		for (uint32_t e = 0; e < mFs2D[f].size(); ++e) {
			uint32_t v0 = mFs2D[f][e], v1 = mFs2D[f][(e + 1) % mFs2D[f].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(mFs2D[f].size());
		FEs[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mEs.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(mEs[E_num]) = false;

		FEs[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}

	if (mV_copy2D.cols()) {
		V_pvs.resize(mV_copy2D.cols());
		for (uint32_t i = 0; i < mEs.size(); i++) {
			uint32_t v0 = get<0>(mEs[i]), v1 = get<1>(mEs[i]);
			V_pvs[v0].push_back(v1);
			V_pvs[v1].push_back(v0);
		}
	}
	E_pfs.resize(mEs.size());	
	for (uint32_t i = 0; i < FEs.size(); i++) for (auto eid : FEs[i]) { E_pfs[eid].push_back(i); }
}
bool simple_polygon(std::vector<std::vector<uint32_t>> &fvs, std::vector<std::vector<uint32_t>> &fes, std::vector<uint32_t> &pvs, 
	std::vector<uint32_t> &pes, std::vector<uint32_t> &vs_disgard, std::vector<uint32_t> &es_disgard)
{
	es_disgard.clear();//es_disgard is in the interior
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++) 
			if (mE_flag[fes[i][j]]) {
				es_disgard.push_back(fes[i][j]); 
				mE_flag[fes[i][j]] = false;
			}
			else mE_flag[fes[i][j]] = true;
	}
	short which_polygon = 0;
	for (int i = 0; i < fes.size(); i++) {
		for (int j = 0; j < fes[i].size(); j++)
			if (mE_flag[fes[i][j]]) { 
				if (!pes.size()) 
					which_polygon = i;

				pes.push_back(fes[i][j]); mE_flag[fes[i][j]] = false;
			}
	}
	//test nvs for each v
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = V_map[std::get<0>(mEs[pes[i]])], v1 = V_map[std::get<1>(mEs[pes[i]])];
		mV_flag[v0]++; mV_flag[v1]++;
		if (mV_flag[v0] > 2 || mV_flag[v1] > 2) {
			for (uint32_t j = 0; j < pes.size(); ++j) {
				uint32_t v0_ = V_map[std::get<0>(mEs[pes[j]])], v1_ = V_map[std::get<1>(mEs[pes[j]])];
				mV_flag[v0_] = mV_flag[v1_] = 0;
			}
			return false;
		}
	}
	for (uint32_t i = 0; i < pes.size(); ++i) {
		uint32_t v0 = V_map[std::get<0>(mEs[pes[i]])], v1 = V_map[std::get<1>(mEs[pes[i]])];
		mV_flag[v0] = mV_flag[v1] = 0;
	}
	//extract the polygon	
	pvs.clear();
	pvs.reserve(pes.size());
	std::vector<bool> e_flag(pes.size(), false);
	pvs.push_back(V_map[std::get<0>(mEs[pes[0]])]);
	pvs.push_back(V_map[std::get<1>(mEs[pes[0]])]);
	e_flag[0] = true;
	uint32_t start_v = pvs[1];
	for (uint32_t i = 2; i < pes.size(); i++) {
		for (uint32_t j = 1; j < pes.size(); j++) {
			if (!e_flag[j]) {
				if (V_map[std::get<0>(mEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(V_map[std::get<1>(mEs[pes[j]])]);
					start_v = V_map[std::get<1>(mEs[pes[j]])];
					break;
				}
				else if (V_map[std::get<1>(mEs[pes[j]])] == start_v) {
					e_flag[j] = true;
					pvs.push_back(V_map[std::get<0>(mEs[pes[j]])]);
					start_v = V_map[std::get<0>(mEs[pes[j]])];
					break;
				}
			}
		}
	}
	if (pvs.size() != pes.size())
		return false;
	//judge direction
	bool correct = true;
	for (int i = 0; i < fvs[which_polygon].size(); i++)
		if (fvs[which_polygon][i] == V_map[std::get<0>(mEs[pes[0]])])
			if (fvs[which_polygon][(i + 1) % fvs[which_polygon].size()] != V_map[std::get<1>(mEs[pes[0]])])
				correct = false;
	if (!correct)
		std::reverse(pvs.begin(), pvs.end());
	//vs_disgard
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			mV_flag[fvs[i][j]] = 1; 
	for (uint32_t i = 0; i < pvs.size(); ++i) mV_flag[pvs[i]] = 0;
	for (uint32_t i = 0; i < fvs.size(); ++i)
		for (uint32_t j = 0; j < fvs[i].size(); ++j)
			if (mV_flag[fvs[i][j]] == 1) {
				vs_disgard.push_back(fvs[i][j]); mV_flag[fvs[i][j]] = 0;
			}
	return true;
}
void reindex_2D(MatrixXf &Vs, std::vector<std::vector<uint32_t>> &F_Vs)
{
	std::vector<int> V_flag(Vs.size(), -1);
	for (uint32_t i = 0; i < F_Vs.size(); i++)
		for (uint32_t j = 0; j < F_Vs[i].size(); j++)
			V_flag[F_Vs[i][j]] = 0;

	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1)
			V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num);
	v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1)
			mV_local_.col(v_num++) = Vs.col(i);
	std::vector<std::vector<uint32_t>> mFs_local_(F_Vs.size());
	for (uint32_t i = 0; i < F_Vs.size(); i++)
		for (uint32_t j = 0; j < F_Vs[i].size(); j++)
			mFs_local_[i].push_back(V_flag[F_Vs[i][j]]);

	mV_local_.swap(Vs);
	mFs_local_.swap(F_Vs);
}
void topology_check_2D(std::vector<std::vector<uint32_t>> &F_Vs, std::vector<std::vector<uint32_t>> &F_Es, int &genus, bool &manifoldness)
{
	std::set<uint32_t> Vs_;
	for (uint32_t i = 0; i < F_Vs.size(); i++) {
		for (uint32_t j = 0; j < F_Vs[i].size(); j++)
			Vs_.insert(F_Vs[i][j]);
	}
	int v_num = Vs_.size();

	std::set<uint32_t> Es_; uint32_t max_eid = 0;
	for (uint32_t i = 0; i < F_Es.size(); i++) {
		for (uint32_t j = 0; j < F_Es[i].size(); j++) {
			Es_.insert(F_Es[i][j]);
			if (max_eid < F_Es[i][j])
				max_eid = F_Es[i][j];
		}
	}
	int e_num = Es_.size();

	std::vector<std::vector<uint32_t>> E_nts(max_eid+1);
	for (uint32_t i = 0; i < F_Es.size(); i++)
		for (uint32_t j = 0; j < F_Es[i].size(); j++)
			E_nts[F_Es[i][j]].push_back(i);

	genus = -(v_num + F_Vs.size() - e_num - 2) / 2;

	manifoldness = true;
	for (uint32_t i = 0; i < E_nts.size(); i++)
		if (E_nts[i].size() > 2) manifoldness = false;
}

void reindex_2D(MatrixXf &HV, MatrixXf &HQ, MatrixXf &HN, MatrixXf &HO, std::vector<std::vector<uint32_t>> &HFv) {
	//re-index V
	std::vector<int32_t> V_flag(HV.size(), -1);
	for (auto fvs : HFv) for (auto vid : fvs) V_flag[vid] = 0;
	uint32_t v_num = 0;
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) V_flag[i] = v_num++;
	MatrixXf mV_local_(3, v_num), mQ_local_(3, v_num), mN_local_(3, v_num), mO_local_(3, v_num);
	for (uint32_t i = 0; i < V_flag.size(); i++)
		if (V_flag[i] != -1) {
			mV_local_.col(V_flag[i]) = HV.col(i);
			mQ_local_.col(V_flag[i]) = HQ.col(i);
			mN_local_.col(V_flag[i]) = HN.col(i);
			mO_local_.col(V_flag[i]) = HO.col(i);
		}
	for (auto &fvs : HFv) for (uint32_t j = 0; j < fvs.size(); j++)
		fvs[j] = V_flag[fvs[j]];
	mV_local_.swap(HV);
	mQ_local_.swap(HQ);
	mN_local_.swap(HN);
	mO_local_.swap(HO);
}

void MultiResolutionHierarchy::tagging_singularities_T_nodes(MatrixXf &V_tagging, vector<tuple_E> &E_tagging, vector<vector<uint32_t>> &F_tagging) {
	enum  V_type
	{
		regular = 0,
		singular,
		t_node,
		boundary,
		s_t_node,
	};

	vector<std::vector<uint32_t>> Vs_nes(V_tagging.cols()), Vs_nfs(V_tagging.cols());
	vector<bool> mV_B_flag(V_tagging.cols(), false);
	V_flag.resize(V_tagging.cols());
	fill(V_flag.begin(), V_flag.end(), V_type::regular);//boundary flag

	for (auto e : E_tagging) if (std::get<2>(e) == 1) { mV_B_flag[std::get<1>(e)] = mV_B_flag[std::get<0>(e)] = true; }
	for (uint32_t i = 0; i < E_tagging.size(); i++) {
		uint32_t v0 = std::get<0>(E_tagging[i]), v1 = std::get<1>(E_tagging[i]);
		Vs_nes[v0].push_back(i);
		Vs_nes[v1].push_back(i);
	}
	for (uint32_t f = 0; f < F_tagging.size(); ++f) for (auto vid : F_tagging[f]) Vs_nfs[vid].push_back(f);

	vector<int32_t> V_tags(V_tagging.cols(), 0);//0-regular, 1-singular, 2-t_node, 3 -boundary, and 4 both singular & t_node
	for (uint32_t i = 0; i < Vs_nfs.size(); i++) {
		if (mV_B_flag[i]) {
			V_flag[i] = V_type::boundary;
			continue;
		}

		if (Vs_nes[i].size() != 4) {
			V_flag[i] = V_type::singular;
			continue;
		}
	}
	for (uint32_t i = 0; i < F_tagging.size(); i++) {
		if (F_tagging[i].size() == 5) {
			vector<int32_t> t_candidates;
			for (uint32_t j = 0; j < F_tagging[i].size(); j++) {
				//if (Vs_nes[F_tag[i][j]].size() == 3 && !mV_B_flag[F_tag[i][j]]) t_candidates.push_back(j);
				if (mV_B_flag[F_tagging[i][j]]) {
					if (Vs_nes[F_tagging[i][j]].size() == 3 || Vs_nes[F_tagging[i][j]].size() == 2) t_candidates.push_back(j);
				}
				else {
					if (Vs_nes[F_tagging[i][j]].size() == 3) t_candidates.push_back(j);
				}
			}
			if (t_candidates.size()) {

				vector<tuple<double, uint32_t>> vs_rank;
				for (auto v_id : t_candidates) {
					int32_t v0_pre = (v_id - 1 + F_tagging[i].size()) % F_tagging[i].size(), v0_aft = (v_id + 1) % F_tagging[i].size();

					MatrixXf bes_vec(3, 2);
					bes_vec.col(0) = (V_tagging.col(F_tagging[i][v_id]) - V_tagging.col(F_tagging[i][v0_pre])).normalized();
					bes_vec.col(1) = (V_tagging.col(F_tagging[i][v_id]) - V_tagging.col(F_tagging[i][v0_aft])).normalized();

					Float dot_ = bes_vec.col(0).dot(bes_vec.col(1));
					Float angle = std::acos(dot_);
					Float cost = std::abs(angle - PAI);// (std::abs(n.sum()));
					vs_rank.push_back(std::make_tuple(cost, F_tagging[i][v_id]));
				}
				sort(vs_rank.begin(), vs_rank.end());
				
				V_flag[get<1>(vs_rank[0])] = V_type::t_node;
			}
			else {
				cout << "singular & t-node pentagon " << i << endl;
			}
		}
	}
}

