#include "hierarchy.h"
#include "timer.h"
#include "quat.h"
#include "tri_tri_intersection.h"
#include <numeric>
#include "hierarchy.h"
#include "quat.h"
#include "timer.h"
#include "orientations.h"
#include <fstream>


void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V) {
	/// Vertex indices used by the OBJ format
	struct obj_vertex {
	  uint32_t p = (uint32_t)-1;
	  uint32_t n = (uint32_t)-1;
	  uint32_t uv = (uint32_t)-1;

	  inline obj_vertex() { }

	  inline obj_vertex(const std::string &string) {
		  std::vector<std::string> tokens = str_tokenize(string, '/', true);

		  if (tokens.size() < 1 || tokens.size() > 3)
			  throw std::runtime_error("Invalid vertex data: \"" + string + "\"");

		  p = str_to_uint32_t(tokens[0]);

#if 0
		  if (tokens.size() >= 2 && !tokens[1].empty())
				uv = str_to_uint32_t(tokens[1]);

			if (tokens.size() >= 3 && !tokens[2].empty())
				n = str_to_uint32_t(tokens[2]);
#endif
	  }

	  inline bool operator==(const obj_vertex &v) const {
		  return v.p == p && v.n == n && v.uv == uv;
	  }
	};

	/// Hash function for obj_vertex
	struct obj_vertexHash : std::unary_function<obj_vertex, size_t> {
	  std::size_t operator()(const obj_vertex &v) const {
		  size_t hash = std::hash<uint32_t>()(v.p);
		  hash = hash * 37 + std::hash<uint32_t>()(v.uv);
		  hash = hash * 37 + std::hash<uint32_t>()(v.n);
		  return hash;
	  }
	};

	typedef std::unordered_map<obj_vertex, uint32_t, obj_vertexHash> VertexMap;

	size_t last_slash_idx = filename.rfind('.');
	if (filename.substr(last_slash_idx) != ".OBJ" && filename.substr(last_slash_idx) != ".obj")
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");

	std::ifstream is(filename);
	if (is.fail())
		throw std::runtime_error("Unable to open OBJ file \"" + filename + "\"!");
	cout << "Loading \"" << filename << "\" .. ";
	cout.flush();
	Timer<> timer;

	std::vector<Vector3f>   positions;
	//std::vector<Vector2f>   texcoords;
	//std::vector<Vector3f>   normals;
	std::vector<uint32_t>   indices;
	std::vector<obj_vertex> vertices;
	VertexMap vertexMap;

	std::string line_str;
	while (std::getline(is, line_str)) {
		std::istringstream line(line_str);

		std::string prefix;
		line >> prefix;

		if (prefix == "v") {
			Vector3f p;
			line >> p.x() >> p.y() >> p.z();
			positions.push_back(p);
		}
		else if (prefix == "vt") {
			/*
			Vector2f tc;
			line >> tc.x() >> tc.y();
			texcoords.push_back(tc);
			*/
		}
		else if (prefix == "vn") {
			/*
			Vector3f n;
			line >> n.x() >> n.y() >> n.z();
			normals.push_back(n);
			*/
		}
		else if (prefix == "f") {
			std::string v1, v2, v3, v4;
			line >> v1 >> v2 >> v3 >> v4;
			obj_vertex tri[6];
			int nVertices = 3;

			tri[0] = obj_vertex(v1);
			tri[1] = obj_vertex(v2);
			tri[2] = obj_vertex(v3);

			if (!v4.empty()) {
				/* This is a quad, split into two triangles */
				tri[3] = obj_vertex(v4);
				tri[4] = tri[0];
				tri[5] = tri[2];
				nVertices = 6;
			}
			/* Convert to an indexed vertex list */
			for (int i = 0; i<nVertices; ++i) {
				const obj_vertex &v = tri[i];
				VertexMap::const_iterator it = vertexMap.find(v);
				if (it == vertexMap.end()) {
					vertexMap[v] = (uint32_t)vertices.size();
					indices.push_back((uint32_t)vertices.size());
					vertices.push_back(v);
				}
				else {
					indices.push_back(it->second);
				}
			}
		}
	}
	F.resize(3, indices.size() / 3);
	memcpy(F.data(), indices.data(), sizeof(uint32_t)*indices.size());
	V.resize(3, vertices.size());
	for (uint32_t i = 0; i<vertices.size(); ++i)
		V.col(i) = positions.at(vertices[i].p - 1);
}

MultiResolutionHierarchy::MultiResolutionHierarchy() {
    mV = { MatrixXf::Zero(3, 1) };
    mN = { MatrixXf::Zero(3, 1) };
    mO = { MatrixXf::Zero(3, 1) };
    mQ = { MatrixXf::Zero(4, 1) };
	ratio_scale = 3.0;
    splitting = triangles = false;
}

bool MultiResolutionHierarchy::load(const MatrixXf& V, const MatrixXu& F) {
	mV.resize(1);
	mV[0] = V;
	mF = F;

	mAABB = AABB(
			mV[0].rowwise().minCoeff(),
			mV[0].rowwise().maxCoeff()
	);

	diagonalLen = 3 * (mAABB.max - mAABB.min).norm() / 100;
	//ratio_scale = ms.mAverageEdgeLength * 3.5 / diagonalLen;

	return true;
}

bool MultiResolutionHierarchy::load(const std::string &filename) {

    mV.resize(1);
	mV[0] = MatrixXf::Zero(3, 1);
	mF = MatrixXu::Zero(3, 1);
	
	try {
		load_obj(filename, mF, mV[0]);
	}
	catch (const std::exception &e) {
		std::cout << "failed loading obj file." << std::endl;
		return false;
	}

//	mV.resize(1);
	mAABB = AABB(
		mV[0].rowwise().minCoeff(),
		mV[0].rowwise().maxCoeff()
	);

	diagonalLen = 3 * (mAABB.max - mAABB.min).norm() / 100;
	//ratio_scale = ms.mAverageEdgeLength * 3.5 / diagonalLen;

    return true;
}


void MultiResolutionHierarchy::build() {
	Timer<> timer;
	mV.resize(1);

	timer.beginStage("Computing face and vertex normals");
	mN.resize(1);
	mC.resize(1);
	nV_boundary_flag.resize(1);
	mN[0].setZero(3, mV[0].cols());
	mC[0].setZero(3, mV[0].cols());
	mNF.resize(3, mF.cols()); mCF.resize(3, mF.cols());
	VectorXi count(mV[0].cols());
	count.setZero();
	for (uint32_t i = 0; i<mF.cols(); ++i) {
		uint32_t i0 = mF(0, i), i1 = mF(1, i), i2 = mF(2, i);
		Vector3f v0 = mV[0].col(i0), v1 = mV[0].col(i1), v2 = mV[0].col(i2);
		Vector3f n = (v1 - v0).cross(v2 - v0).normalized();
		mNF.col(i) = n;
		mCF.col(i) += (v0 + v1 + v2) / 3;
		mN[0].col(i0) += n; mN[0].col(i1) += n; mN[0].col(i2) += n;
		count[i0]++; count[i1]++; count[i2]++;
	}

	for (uint32_t i = 0; i<mN[0].cols(); ++i) {
		if (mN[0].col(i) != Vector3f::Zero()) {
			Vector3f d1 = mN[0].col(i) / count[i],
				d2 = mN[0].col(i).normalized();
			mN[0].col(i) = d2;
			if (d2 != Vector3f::Zero())
				mC[0].col(i) = mV[0].col(i);
		}
	}

	vnfs.clear();
	vnfs.resize(mV[0].cols());
	for (uint32_t i = 0; i < mF.cols(); ++i) for (uint32_t j = 0; j < 3; j++) vnfs[mF(j, i)].push_back(i);

	timer.endStage();

	timer.beginStage("Computing adjacency data structure");
	mL.clear(); mP.clear();
	nV_boundary_flag[0].clear(); nV_boundary_flag[0].resize(mV[0].cols(), false);


		construct_tEs_tFEs(mF, nFes, nEs);
		//nV_nes, tag boundary V
		nV_nes.clear(); nV_nes.resize(mV[0].cols());
		for (uint32_t i = 0; i < nEs.size(); i++) {
			uint32_t v0 = std::get<0>(nEs[i]);
			uint32_t v1 = std::get<1>(nEs[i]);
			nV_nes[v0].push_back(i);
			nV_nes[v1].push_back(i);
			if (std::get<2>(nEs[i])) {
				nV_boundary_flag[0][v0] = nV_boundary_flag[0][v1] = true;
			}
		}

		std::vector<std::pair<uint32_t, uint32_t>> adj;
		adj.reserve(mF.cols() * 6);
		for (uint32_t f = 0; f < mF.cols(); ++f) {
			for (int i = 0; i < 3; ++i) {
				uint32_t v0 = mF(i, f);
				uint32_t v1 = mF((i + 1) % 3, f);
				adj.push_back(std::make_pair(v0, v1));
				adj.push_back(std::make_pair(v1, v0));
			}
		}
		std::sort(adj.begin(), adj.end());
		adj.erase(std::unique(adj.begin(), adj.end()), adj.end());

		std::vector<Triplet> triplets;
		for (auto item : adj)
			triplets.push_back(Triplet(item.first, item.second, 1.f));
		mL.resize(1);
		mL[0].resize(mV[0].cols(), mV[0].cols());
		mL[0].setFromTriplets(triplets.begin(), triplets.end());
	for (uint32_t i = 0; i < (uint32_t)mL[0].rows(); ++i) {
		Float sum = 1 / mL[0].row(i).sum();
		mL[0].row(i) *= sum;
		mL[0].coeffRef(i, i) = -sum;
	}
	mL[0].makeCompressed();
	timer.endStage();

	struct WeightedEdge {
		WeightedEdge(uint32_t _i0, uint32_t _i1, Float weight)
			: weight(weight), i0(_i0), i1(_i1) {
			if (i0 > i1)
				std::swap(i0, i1);
		}

		bool operator<(const WeightedEdge &e) const {
			return std::tie(weight, i0, i1) < std::tie(e.weight, e.i0, e.i1);
		}

		Float weight;
		uint32_t i0, i1;
	};

	timer.beginStage("Building hierarchy");
	while (mL[mL.size() - 1].cols() > 1) {
		const MatrixXf &V = mV[mV.size() - 1];
		const MatrixXf &N = mN[mN.size() - 1];
		const MatrixXf &C = mC[mC.size() - 1];
		const vector<bool> &VB = nV_boundary_flag[nV_boundary_flag.size() - 1];
		const SMatrix &L = mL[mL.size() - 1];
		std::vector<bool> collapsed(L.cols(), false);
		std::vector<bool> visited(L.cols(), false);
		std::set<WeightedEdge> edges;

		double edgeSum = 0;
		size_t edgeCount = 0;
		for (int k = 0; k < L.outerSize(); ++k) {
			for (SMatrix::InnerIterator it(L, k); it; ++it) {
				if (it.col() == it.row())
					continue;
				Float length = (V.col(it.row()) - V.col(it.col())).norm();
				edgeSum += length;
				edgeCount += 1;
				edges.insert(WeightedEdge(it.row(), it.col(), length));
			}
		}

		std::vector<Triplet> P_triplets, R_triplets;
		std::vector<Vector3f> V_next, N_next, C_next;
		std::map<uint32_t, uint32_t> vertex_map;

		uint32_t nVertices = 0; vector<bool> vb_flag(V.cols(), false);
		for (auto const &e : edges) {
			visited[e.i0] = visited[e.i1] = true;
			if (collapsed[e.i0] || collapsed[e.i1])
				continue;
			collapsed[e.i0] = true;
			collapsed[e.i1] = true;
			P_triplets.push_back(Triplet(e.i0, nVertices, 1.0f));
			P_triplets.push_back(Triplet(e.i1, nVertices, 1.0f));
			R_triplets.push_back(Triplet(nVertices, e.i0, 0.5f));
			R_triplets.push_back(Triplet(nVertices, e.i1, 0.5f));
			V_next.push_back(0.5f * (V.col(e.i0) + V.col(e.i1)));

			if (VB[e.i0] || VB[e.i1]) vb_flag[nVertices] = true;

			Vector3f n = N.col(e.i0) + N.col(e.i1);
			Vector3f c = C.col(e.i0) + C.col(e.i1);
			if (N.col(e.i0) != Vector3f::Zero() &&
				N.col(e.i1) != Vector3f::Zero()) {
					n.normalize();
					c *= 0.5f;
			}

			N_next.push_back(n);
			C_next.push_back(c);

			vertex_map[e.i0] = nVertices;
			vertex_map[e.i1] = nVertices;
			nVertices++;
		}

		for (uint32_t i = 0; i<V.cols(); ++i) {
			if (collapsed[i] || !visited[i])
				continue;
			P_triplets.push_back(Triplet(i, nVertices, 1.0f));
			R_triplets.push_back(Triplet(nVertices, i, 1.0f));
			V_next.push_back(V.col(i));
			N_next.push_back(N.col(i));
			C_next.push_back(C.col(i));
			vertex_map[i] = nVertices;

			if (VB[i]) vb_flag[nVertices] = true;

			nVertices++;
		}
		vb_flag.resize(nVertices);

		if (mL.size() != 1)
			std::cout << ", ";
		std::cout << nVertices;
		std::cout.flush();

		SMatrix P(V.cols(), nVertices), R(nVertices, V.cols());

		P.setFromTriplets(P_triplets.begin(), P_triplets.end());
		R.setFromTriplets(R_triplets.begin(), R_triplets.end());

		SMatrix L2 = R*L*P;
		MatrixXf V2(3, nVertices), N2(3, nVertices), C2(3, nVertices), Q2(4, nVertices);
		for (uint32_t i = 0; i<nVertices; ++i) {
			V2.col(i) = V_next[i];
			N2.col(i) = N_next[i];
			C2.col(i) = C_next[i];
		}

		nV_boundary_flag.push_back(vb_flag);
		mP.push_back(std::move(P));
		mN.push_back(std::move(N2));
		mV.push_back(std::move(V2));
		mC.push_back(std::move(C2));
		mL.push_back(L2);
	}
	std::cout << " ";
	timer.endStage();

	mQ.resize(mL.size());
	mO.resize(mL.size());

	pcg32 rng;
		for (uint32_t i = 0; i < mL.size(); ++i) {
			mQ[i].resize(3, mV[i].cols());
			mO[i].resize(3, mV[i].cols());

			for (uint32_t j = 0; j < mV[i].cols(); ++j) {
				if (i == 0 && nV_boundary_flag[i][j]) {
					std::vector<uint32_t> vs;
					for (auto eid : nV_nes[j]) {
						if (!std::get<2>(nEs[eid])) continue;

						uint32_t v0 = std::get<0>(nEs[eid]);
						uint32_t v1 = std::get<1>(nEs[eid]);
						if (v0 == j) vs.push_back(v1);
						else vs.push_back(v0);
					}
					Vector3f direct0, direct1; direct1.setZero();
					direct0 = (mV[0].col(j) - mV[0].col(vs[0])).normalized();
					for (uint32_t k = 1; k < vs.size(); k++) {
						Vector3f direct_ = (mV[0].col(vs[k]) - mV[0].col(j)).normalized();
						direct1 += direct_;
					}

					auto newdir = direct0;
//					if (std::abs(direct0.dot(direct1)) < 0.5)
						newdir = direct0;
//					else
//						newdir = (direct0 + direct1).normalized();
//					if(newdir.hasNaN()) {
//						std::cout<<newdir<<std::endl;
//						std::cout<<"La";
//					}
                    mQ[i].col(j) = newdir;
					continue;
				}

				Vector3f n = mN[i].col(j), v = mV[i].col(j);
				Vector3f s, t;
				coordinate_system(n, s, t);
				float angle = rng.nextFloat() * 2 * M_PI;
				mQ[i].col(j) = s * std::cos(angle) + t * std::sin(angle);
			}


			for (uint32_t j = 0; j < mV[i].cols(); ++j) {
				Vector3f n = mN[i].col(j), v = mV[i].col(j);
				rng.nextFloat();
				Vector3f o = aabbRand(mAABB, rng);
				o -= n.dot(o - v) * n;
				mO[i].col(j) = o;
			}
		}
		//propagate up
		for (uint32_t i = 1; i < mL.size(); ++i) {
			for (int k = 0; k < mP[i - 1].outerSize(); ++k) {
				SMatrix::InnerIterator it(mP[i - 1], k);
				for (; it; ++it) {
					if (nV_boundary_flag[i - 1][it.row()])
						mQ[i].col(it.col()) = mQ[i - 1].col(it.row());
				}
			}
		}
	mOrientationIterations = 0;


	sta.tN = mF.cols();
	sta.tetN = mT.cols();
}
void MultiResolutionHierarchy::construct_tEs_tFEs(MatrixXu & F, std::vector<std::vector<uint32_t>> &mtFes, std::vector<tuple_E> &mtEs) {
	mtFes.clear(); mtEs.clear();

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(F.cols() * 3);
	mtFes.resize(F.cols());
	for (uint32_t f = 0; f < F.cols(); ++f) {
		for (uint32_t e = 0; e < 3; ++e) {
			uint32_t v0 = F(e, f), v1 = F((e + 1) % 3, f);
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(3);
		mtFes[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mtEs.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) || std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mtEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(mtEs[E_num]) = false;

		mtFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
}


void MultiResolutionHierarchy::smoothOrientationsTri(uint32_t l, bool alignment, bool randomization, bool extrinsic) {
	const MatrixXf &N = mN[l];
	const SMatrix &L = mL[l];
	MatrixXf &Q = mQ[l];

	Timer<> timer;
	double error = 0;
	int nLinks = 0;
	MatrixXf Q_new(Q.rows(), Q.cols());

#if PARALLEL
	tbb::spin_mutex mutex;
	tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t) L.outerSize(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
            std::vector<std::pair<uint32_t, Float>> neighbors;
            double errorLocal = 0;
            int nLinksLocal = 0;
            for (uint32_t k = range.begin(); k != range.end(); ++k) {
#else
	std::vector<std::pair<uint32_t, Float>> neighbors;
	double errorLocal = 0;
	int nLinksLocal = 0;
	for (uint32_t k = 0; k < (uint32_t) L.outerSize(); k++) {
#endif
		SMatrix::InnerIterator it(L, k);

		uint32_t i = it.row();
		Vector3f q_i = Vector3f::Zero();
		Vector3f n_i = N.col(i);

		if (nV_boundary_flag[l][i]) {
			Q_new.col(i) = Q.col(i);
			continue;
		}

		neighbors.clear();
		for (; it; ++it) {
			uint32_t j = it.col();
			if (i == j)
				continue;
			neighbors.push_back(std::make_pair(j, it.value()));
		}

		if (randomization && neighbors.size() > 0)
			pcg32(mOrientationIterations, k)
					.shuffle(neighbors.begin(), neighbors.end());

		for (auto n : neighbors) {
			uint32_t j = n.first;
			Float value = n.second;
			Float dp;

			Vector3f q_j = Q.col(j), n_j = N.col(j);
			if (extrinsic) {
				q_j = applyRotationExtrinsic((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
				dp = Q.col(i).dot(applyRotation(Q.col(i), N.col(i), Q.col(j), N.col(j)));
			} else {
				q_j = applyRotation((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);
				dp = Q.col(i).dot(applyRotation(Q.col(i), N.col(i), Q.col(j), N.col(j)));
			}

			errorLocal += std::abs(std::acos(std::min(dp, (Float) 1)));
			++nLinksLocal;

			q_i += q_j * value;
		}

		if (q_i != Vector3f::Zero())
			Q_new.col(i) = (q_i - n_i.dot(q_i) * n_i).normalized();
	}
	error += errorLocal;
	nLinks += nLinksLocal;
#if PARALLEL
	tbb::spin_mutex::scoped_lock guard(mutex);
	}
    );
#else
#endif
	mOrientationIterations++;
	Q = std::move(Q_new);
}

void MultiResolutionHierarchy::prolongOrientations(int level) {
	const SMatrix &P = mP[level];
	const MatrixXf &N = mN[level];
	for (int k = 0; k < P.outerSize(); ++k) {
		SMatrix::InnerIterator it(P, k);
		for (; it; ++it) {
			Quaternion q_j = mQ[level + 1].col(it.col());
			Vector3f n_i = N.col(it.row());
			if (n_i != Vector3f::Zero()) {
				Float magnitude = q_j.norm();
				q_j = Quaternion(q_j / magnitude).align(n_i) * magnitude;
			}
			mQ[level].col(it.row()) = q_j;
		}
	}
}

