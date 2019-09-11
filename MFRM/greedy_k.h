#pragma once

#include <string> // std::string
#include <fstream> // std::ifstream, std::ofstream
#include <chrono> // std::chrono::high_resolution_clock
#include <limits> // FLT_MAX
#include <algorithm> // std::set_difference, std::set_intersection
#include <boost/unordered_set.hpp> // for consistency reasons, we use boost::unordered_set instead of std::unordered_set

#include "geo_defines.h" // R*-tree, EARTH_RADIUS, PI
#include "typed_directed_graph.h"

//================================================================================
// facility-replacement pair associated with the expected change of total distance (ED)

struct facility_replacement
{
	int f_id; // facility index number
	int c_id; // candidate index number
	float ED; // the expected change of total distance (ED)

	facility_replacement()
		: f_id(-1)
		, c_id(-1)
		, ED(0.0f) {}

	facility_replacement(int f_id_in, int c_id_in, float ED_in)
		: f_id(f_id_in)
		, c_id(c_id_in)
		, ED(ED_in) {}

	// remark: the same facility-replacement pair validation is not implemented
	bool operator<(const facility_replacement &rhs) const
	{
		return ED < rhs.ED; // only compare ED
	};

	// boost::hash is implemented by calling the function hash_value
	friend std::size_t hash_value(const facility_replacement &v)
	{
		std::size_t seed = 0;
		boost::hash_combine(seed, v.f_id);
		boost::hash_combine(seed, v.c_id);
		return seed;
	}
};

//--------------------------------------------------------------------------------
// a mutable fibonacci max-heap for ordering facility_replacement by ED;

typedef boost::heap::fibonacci_heap<facility_replacement> fac_rep_max_fibonacci_heap;


//================================================================================
// REMARK:
//	if a reference location has no nearest or sub-nearest facility (namely isolate sub-grahp), we ignore it and regard its contribute on ED as 0


//================================================================================
// testing

//#define ALGO_TESTING
//#define ALGO_TESTING_TRACE_GRAPH
//#define ALGO_TESTING_NO_EARLY_STOPPING
//#define ALGO_TESTING_GREEDY

static std::ofstream algo_out_file("E:\\Experiment\\MFRM\\datasets\\algo_testing.txt", std::ofstream::out | std::ofstream::trunc);


//================================================================================
// Local Network Table (LNT) for LNB algorithm;
// an assistant hash map for local network edges, which consists of shortest distances from a reference location to a vertex, and that distance plus out-edge weight (distance)
// an assistant hash map for information of candidates
// an assistant hash map for reference locations that overlap some vertices
// Remark: "refloc_hash_map" is not discussed in the paper

typedef boost::unordered_map < std::pair < int, int >, // ids of starting and ending vertices of a directed edge, whose format is <vs_id, ve_id>
	std::pair < float, // dist(r, vs)
	float > > // dist(r, vs) + dist(vs, ve)
loc_edges_hash_map; // assistant hash map

typedef boost::tuple < int, // refloc_id, REFLOC
	float, // Dd+ (Sigma upper bound), DDU
	float, // Dd- (Sigma lower bound), DDL
	float > // offset, OFFSET
ref_loc_entry;

enum ref_loc_entry_element { REFLOC = 0, DDU, DDL, OFFSET };

typedef boost::unordered_map < std::pair < int, int >, // ids of starting and ending vertices of a directed edge, whose format is <vs_id, ve_id>
	std::pair < float, // ED+ (ED upper bound on a directed edge)
	std::vector<ref_loc_entry> > >
LNT_hash_map;

typedef boost::tuple < int, // refloc id, use REFLOC of "ref_loc_entry_element"
	float, // Dd+ (Sigma upper bound), DDU_2
	float, // Dd- (Sigma lower bound), DDL_2
	float > // offset, OFFSET_2
ref_loc_2_entry;

enum ref_loc_2_entry_element { DDU_2 = 1, DDL_2, OFFSET_2 };

typedef boost::unordered_map < std::pair < int, int >, // ids of starting and ending vertices of a directed edge, whose format is <vs_id, ve_id>
	std::vector<ref_loc_2_entry> >
L2NT_hash_map;

typedef boost::unordered_map < int, // candidate id
	boost::tuple < int, // vs_id, C_VS
	int, // ve_id, C_VE
	float, // dist_vs
	float > > // dist_ve
cand_info_hash_map; // assistant hash map

enum cand_info_hash_map_element { C_VS = 0, C_VE }; // the value of VS_DIST, VE_DIST are defined and consistent with "fac_hash_map_element"

typedef boost::unordered_map < int, // vertex id
	float > // accumulated dnn
vertex_dnn_hash_map; // assistant hash map

typedef boost::unordered_map < int, // facility id
	std::pair < float, // Delta(f), aggregate (dnn - d2nn) * Pr(r)
	std::set < int > > > // rnn(F,f)
delta_f_hash_map; // <f, \Delta(f), rnn(F,f)>

typedef boost::unordered_map < int, // refloc id
	boost::tuple < float, // prob, PROB
	int, // nn facility id, NN
	int, // 2nn facility id, SNN
	float > > // (d2nn - dnn) * prob, DIFF
refloc_hash_map; // not discussed in the paper

enum refloc_element { PROB = 0, NN, SNN, DIFF };

typedef boost::unordered_map < int, // refloc id
	std::set < std::pair <int, int> > > // <vs_id, ve_id>, which is inside LN(r) or L2N(r)
refloc_LN_or_L2N_hash_map;

#pragma region TOP_VISTOR
//--------------------------------------------------------------------------------
// top event visitor for LNB algorithm
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//			LNT and vertex_nnd hash maps will never record any local network of the reference location, and we view its ERD as 0

class LNB_top_visitor
	: public top_event_visitor
{
public:
	LNB_top_visitor() { ptr_obsolete_facs = NULL; } // default must be NULL

	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map)
	{
		if (v.type == 'f')
		{
			if (ptr_obsolete_facs != NULL && ptr_obsolete_facs->find(v.id) != ptr_obsolete_facs->end()) // the v.id facility has been obsolete
				return false; // a flat to continue graph traversal

			if (nn_id == -1) // nearest facility is found
			{
				nn_id = v.id; // nearest facility
				dnn = v.dist; // distance from a reference location to its nearest facility
			}
			else // sub-nearest facility is found
			{
				if (ptr_dnn_hash_map != NULL // vertex_dnn_hash_map is not set
					&& r_v_id != -1) // the reference location overlaps a vertex (vertex id is r_v_id)
				{
#ifdef ALGO_TESTING
					//algo_out_file << "r overlaps v" << r_v_id << ", dnn = " << v.dist << '\n';
#endif
					vertex_dnn_hash_map::iterator iter_dnn = ptr_dnn_hash_map->find(r_v_id);
					if (iter_dnn != ptr_dnn_hash_map->end()) // the vertex has been overlapped by some reference location
						iter_dnn->second += dnn; // accumulate dnn
					else // the vertex is overlapped by this reference location first time
						(*ptr_dnn_hash_map)[r_v_id] = dnn; // initial dnn
				}

				// consturct reflocs hash map, ie., <refloc_id, prob, nn, snn, (d2nn - dnn) * Pr(r)>
				float diff = (dnn - v.dist) * prob; // (dnn - d2nn) * Pr(r)
				(*ptr_refloc_hash_map)[refloc_id] = boost::make_tuple(prob, nn_id, v.id,
					-diff); // (d2nn - dnn) * Pr(r)

				// construct <f, rnn(F,r), \Delta(f)> tuples
				delta_f_hash_map::iterator iter_delta = ptr_delta_f_hash_map->find(nn_id);
				if (iter_delta != ptr_delta_f_hash_map->end()) // the facility has been in the tuples
				{
					iter_delta->second.first += diff; // (dnn - d2nn) * prob, aggregate \Delta(f)
					iter_delta->second.second.insert(refloc_id); // r \in rnn(F,f)
				}
				else // the first time presents the facility
				{
					std::set<int> rnn_set;
					rnn_set.insert(refloc_id);
					(*ptr_delta_f_hash_map)[nn_id] = std::make_pair(diff, rnn_set);
				}

				// deal with each edge in the local network/sub-network of a reference location, ie., edges in LN(r)/L2N(r)
				for (loc_edges_hash_map::const_iterator iter_edge = ptr_loc_edges_hash_map->cbegin(); iter_edge != ptr_loc_edges_hash_map->cend(); ++iter_edge)
				{
					// for clarity, we define a variable for edge id
					std::pair<int, int> vs_ve = iter_edge->first; // the edge vs->ve, <vs_id, ve_id>

					// determine whether <vs_id, ve_id> is in the local network
					bool is_in_LNT = true; // indicate that <vs_id, ve_id> is in the local network
					if (iter_edge->second.first >= dnn) // d(r, vs) ">" for a facility spliting an edge; "==" for a facility overlapping an ordinary vertex
						is_in_LNT = false; // not located in LNT but in L2NT

					bool is_in_L2NT = true; // indicate that <vs_id, ve_id> is in the local sub-network
											// 1) no definite shortest path for ve (ie, d(r, ve) is not definite), namely frontier edge
											// 2) definite d(r, ve) and d(r, ve) > dnn
	
					//notice:when Dd- + d(c, vj) >= 0,"<vi, vj>(ie, <vs, ve>)"is divided into two ends "<vi,nn>" and "<nn,ve>" by "nn" to,c is in <vi,nn>,
					//       as is in LN,is not in L2NT,so this value is not calculated when calculating L2NT; 
					//       when Dd- + d(c, vj) < dnn - d2nn, c is out of 2nn,as c is not in L2N,so C doesn't affect f;
					// compute offset
					float offset = 0.0f;
					definite_hash_map::const_iterator iter_definite = hash_map.find(vertex(vs_ve.second, 'v'));
					if (iter_definite != hash_map.cend()) // ve has definite shortest path distance d(r, ve)
					{
						if (iter_definite->second < iter_edge->second.second) // d(r, ve) < d(r, vs) + d(vs, ve)
							offset = (iter_edge->second.second - iter_definite->second) / 2.0f; // (d(r, vs) + d(vs, ve) - d(r, ve)) / 2

						// determine whether <vs_id, ve_id> is in local network and/or local sub-network
						if (iter_edge->second.second <= dnn) // d(r, vs) + d(vs, ve) <= dnn
							is_in_L2NT = false; // not located in L2NT
					}

					// create new entry in LNT
					float Dd_u = dnn - iter_edge->second.first; // Dd+ = dnn - dist(r, vs)
					if (is_in_LNT)
					{
						// handle LNT
						LNT_hash_map::iterator iter_LNT = ptr_LNT_hash_map->find(vs_ve);
						if (iter_LNT != ptr_LNT_hash_map->end()) // edge <vs_id, ve_id> has already been in local network of other reference locations
						{
							iter_LNT->second.first += (Dd_u * prob); // new_ED+ = old_ED+ + (Dd+ * prob)
							iter_LNT->second.second.push_back(boost::make_tuple(refloc_id,
								Dd_u, // Dd+
								dnn - iter_edge->second.second, // Dd- = dnn - (dist(r, vs) + dist(vs, ve))
								offset)); // offset
						}
						else // <vs_id, ve_id> first time presents
						{
							(*ptr_LNT_hash_map)[vs_ve] = std::make_pair(Dd_u * prob, // Dd+ as initial ED+
								std::vector<ref_loc_entry>());
							(*ptr_LNT_hash_map)[vs_ve].second.push_back(boost::make_tuple(refloc_id,
								Dd_u, // Dd+
								dnn - iter_edge->second.second, // Dd- = dnn - (dist(r, vs) + dist(vs, ve))
								offset)); // offset
						}

						// handle refloc_LN_hash_map
						(*ptr_refloc_LN_hash_map)[refloc_id].insert(vs_ve);
					}

					// create new entry in L2NT
					if (is_in_L2NT)
					{
						// handle L2NT
						L2NT_hash_map::iterator iter_L2NT = ptr_L2NT_hash_map->find(vs_ve);
						if (iter_L2NT != ptr_L2NT_hash_map->end()) // edge <vs_id, ve_id> has already been in local sub-network of other reference locations
						{
							iter_L2NT->second.push_back(boost::make_tuple(refloc_id,
								Dd_u, // Dd+
								dnn - iter_edge->second.second, // Dd- = dnn - (dist(r, vs) + dist(vs, ve))
								offset)); // offset
						}
						else // <vs_id, ve_id> first time presents
						{
							(*ptr_L2NT_hash_map)[vs_ve] = std::vector<ref_loc_2_entry>();
							(*ptr_L2NT_hash_map)[vs_ve].push_back(boost::make_tuple(refloc_id,
								Dd_u, // Dd+
								dnn - iter_edge->second.second, // Dd- = dnn - (dist(r, vs) + dist(vs, ve))
								offset)); // offset
						}

						// handle refloc_L2N_hash_map
						(*ptr_refloc_L2N_hash_map)[refloc_id].insert(vs_ve);
					}
				}
				return true; // a flag to terminate graph traversal
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(LNT_hash_map *LNT_hash_map_in, L2NT_hash_map *L2NT_hash_map_in, delta_f_hash_map *delta_f_hash_map_in,
		refloc_hash_map *refloc_hash_map_in, loc_edges_hash_map *loc_edges_hash_map_in, vertex_dnn_hash_map *dnn_hash_map,
		refloc_LN_or_L2N_hash_map *refloc_LN_hash_map_in, refloc_LN_or_L2N_hash_map *refloc_L2N_hash_map_in)
	{
		ptr_LNT_hash_map = LNT_hash_map_in;
		ptr_L2NT_hash_map = L2NT_hash_map_in;
		ptr_delta_f_hash_map = delta_f_hash_map_in;
		ptr_refloc_hash_map = refloc_hash_map_in;
		ptr_loc_edges_hash_map = loc_edges_hash_map_in;
		ptr_dnn_hash_map = dnn_hash_map;
		ptr_refloc_LN_hash_map = refloc_LN_hash_map_in;
		ptr_refloc_L2N_hash_map = refloc_L2N_hash_map_in;
	}
	void set_overlap(int r_v_id_in) { r_v_id = r_v_id_in; };
	void set_refloc(int refloc_id_in, float prob_in) {
		refloc_id = refloc_id_in; prob = prob_in; // set the reference location and its present probability
		nn_id = -1; dnn = 0.0f;	} // reset nearest facility
	void set_obsolete_facs(std::set<int> *obsolete_facs_in) { ptr_obsolete_facs = obsolete_facs_in; };

private:
	LNT_hash_map *ptr_LNT_hash_map; // a pointer to a hash map as Local Network Table
	L2NT_hash_map *ptr_L2NT_hash_map; // a pointer to a hash map as Local Sub-Network Table
	delta_f_hash_map *ptr_delta_f_hash_map; // a pointer to a hash map for <f, \Delta(f)> pairs
	refloc_hash_map *ptr_refloc_hash_map; // a pointer to a hash map of reference locations
	loc_edges_hash_map *ptr_loc_edges_hash_map; // a pointer to an assistant hash map for local network edges
	vertex_dnn_hash_map *ptr_dnn_hash_map; // a pointer to an assistant hash map for reference locations that overlap some vertices
	refloc_LN_or_L2N_hash_map *ptr_refloc_LN_hash_map; // a pointer to hash map for LN edges of reference location
	refloc_LN_or_L2N_hash_map *ptr_refloc_L2N_hash_map; // a pointer to hash map for L2N edges of reference location
	std::set<int> *ptr_obsolete_facs; // a set for obsolete facilities
	int r_v_id; // vertex id, which indicates whether this reference location overlaps a vertex (v_id) or not (-1)
	int refloc_id; // id of the reference location which is the source vertex of the current traversal
	float prob; // present probability of the reference location
	int nn_id; // id of the nearest facility of a reference location
	float dnn; // the distance from a reference location to its nearest facility
};
#pragma endregion TOP_VISTOR

#pragma region EXT_VISTOR
//--------------------------------------------------------------------------------
// extend event visitor for LNB algorithm

class LNB_extend_visitor
	: public extend_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const out_edge &e)
	{
		if (v.type == 'v')
		{
			if (e.ve.type != 'f') // must be 'v' (e.dist is impossible to be zero), then v->e.ve
			{
				(*ptr_loc_edges_hash_map)[std::make_pair(v.id, e.ve.id)] = std::make_pair(v.dist, v.dist + e.dist);
#ifdef ALGO_TESTING
				//algo_out_file << v.id << "-" << e.ve.id << ", d(r,vs) = " << v.dist << ", d(r,vs)+d(vs,ve) = " << v.dist + e.dist << '\n';
#endif
				return false; // a flag to continue graph traversal
			}
			else // e.ve.type == 'f'
			{
				if (e.dist == 0) // f overlaps v, then no need to consider this virtual edge
					return false; // a flag to continue graph traversal
				else // e.dist != 0, e.ve is a facility that must split an edge
				{
					fac_hash_map::const_iterator iter_fac = ptr_fac_hash_map->find(e.ve.id); // the facility must exist
					if (iter_fac == ptr_fac_hash_map->cend()) // fault tolerance, because of obsolete facilities
						return false; // a flag to continue graph traversal

					if (v.id == iter_fac->second.get<VS_ID>()) // f.vs(v)->f.ve
					{
						(*ptr_loc_edges_hash_map)[std::make_pair(v.id, iter_fac->second.get<VE_ID>())] = std::make_pair(v.dist, v.dist + iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>());
#ifdef ALGO_TESTING
						//algo_out_file << v.id << "-" << iter_fac->second.get<VE_ID>() << ", d(r,vs) = " << v.dist << ", d(r,vs)+d(vs,ve) = " << v.dist + iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>() << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
					else // v.id == iter_fac->second.get<VE_ID>(), f.ve(v)->f.vs
					{
						(*ptr_loc_edges_hash_map)[std::make_pair(v.id, iter_fac->second.get<VS_ID>())] = std::make_pair(v.dist, v.dist + iter_fac->second.get<VE_DIST>() + iter_fac->second.get<VS_DIST>());
#ifdef ALGO_TESTING
						//algo_out_file << v.id << "-" << iter_fac->second.get<VS_ID>() << ", d(r,vs) = " << v.dist << ", d(r,vs)+d(vs,ve) = " << v.dist + iter_fac->second.get<VE_DIST>() + iter_fac->second.get<VS_DIST>() << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
				}
			}
		}
		else if (v.type == 'r')
		{
			if (e.dist == 0) // v(r) overlaps e.ve, whatever e.ve.type is 'v' or 'f', no need to consider the virtual edge
				return false; // a flag to continue graph traversal
			else // e.dist != 0, v(r) must split an edge
			{
				if (e.ve.type != 'f') // e.ve is ordinary vertex
				{
					if (e.ve.id == ve_id) // r.vs->r.ve
					{
						(*ptr_loc_edges_hash_map)[std::make_pair(vs_id, ve_id)] = std::make_pair(0.0f, e.dist);
#ifdef ALGO_TESTING
						//algo_out_file << vs_id << "-" << ve_id << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << e.dist << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
					else // e.ve.id == vs_id, r.ve->r.vs
					{
						(*ptr_loc_edges_hash_map)[std::make_pair(ve_id, vs_id)] = std::make_pair(0.0f, e.dist);
#ifdef ALGO_TESTING
						//algo_out_file << ve_id << "-" << vs_id << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << e.dist << '\n';
#endif
						return false; // a flag to continue graph traversal
					}
				}
				else // e.ve.type == 'f', r and f locate on and split the same edge
				{
					fac_hash_map::const_iterator iter_fac = ptr_fac_hash_map->find(e.ve.id); // the facility must exist
					if (iter_fac == ptr_fac_hash_map->cend()) // fault tolerance, because of obsolete facilities
						return false; // a flag to continue graph traversal

					if (iter_fac->second.get<VS_ID>() == vs_id)
					{
						if (iter_fac->second.get<VS_DIST>() < dist_vs) // v(r).ve->v(r).vs
						{
							(*ptr_loc_edges_hash_map)[std::make_pair(ve_id, vs_id)] = std::make_pair(0.0f, dist_vs);
#ifdef ALGO_TESTING
							//algo_out_file << ve_id << "-" << vs_id << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_vs << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
						else // iter_fac->second.get<VS_DIST>() > dist_vs (== never holds), v(r).vs->v(r).ve
						{
							(*ptr_loc_edges_hash_map)[std::make_pair(vs_id, ve_id)] = std::make_pair(0.0f, dist_ve);
#ifdef ALGO_TESTING
							//algo_out_file << vs_id << "-" << ve_id << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_ve << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
					}
					else // iter_fac->second.get<VS_ID>() == ve_id
					{
						if (iter_fac->second.get<VS_DIST>() < dist_ve) // v(r).vs->v(r).ve
						{
							(*ptr_loc_edges_hash_map)[std::make_pair(vs_id, ve_id)] = std::make_pair(0.0f, dist_ve);
#ifdef ALGO_TESTING
							//algo_out_file << vs_id << "-" << ve_id << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_ve << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
						else // iter_fac->second.get<VS_DIST>() > dist_ve (== never holds), v(r).ve->v(r).vs
						{
							(*ptr_loc_edges_hash_map)[std::make_pair(ve_id, vs_id)] = std::make_pair(0.0f, dist_vs);
#ifdef ALGO_TESTING
							//algo_out_file << ve_id << "-" << vs_id << ", d(r,vs) = " << 0 << ", d(r,vs)+d(vs,ve) = " << dist_vs << '\n';
#endif
							return false; // a flag to continue graph traversal
						}
					}
				}
			}
		}
		else // if (v.type == 'f'), the edge that v (ie., nn) locates on has been recorded when encountered v
			return false; // a flag to continue graph traversal
	}

	void set_data_structure(fac_hash_map *fac_hash_map_in, loc_edges_hash_map *loc_edges_hash_map_in) {
		ptr_fac_hash_map = fac_hash_map_in; ptr_loc_edges_hash_map = loc_edges_hash_map_in; }
	void set_vs_ve(int vs_in, int ve_in, float dist_vs_in, float dist_ve_in) { vs_id = vs_in; ve_id = ve_in; dist_vs = dist_vs_in; dist_ve = dist_ve_in; }

private:
	fac_hash_map *ptr_fac_hash_map; // a pointer to a hash map of facilities associated with endpoints
	loc_edges_hash_map *ptr_loc_edges_hash_map; // a pointer to an assistant hash map for local network edges
	int vs_id, ve_id; // the original starting and ending ordinary vertices of the edge which the reference location splits
	float dist_vs, dist_ve; // the distances from reference location to vs and ve
};
#pragma endregion EXT_VISTOR

//--------------------------------------------------------------------------------
// LNB algorithm: update LNT, L2NT and \Delta(f) if <f, c> is the previously selected optimal facility-candidate pair
// remark: dnn_map (in the following two functions) is not considered in this function 
//         facs_map must be carefully addressed (but not yet now)
// [in] op: the last optimal facility_replacement
// [in/out] LNT: a hash map as Local Network Table, which needs to be updated in this function
// [in/out] L2NT: a hash map as Local Sub-Network Table, which needs to be updated in this function
// [in/out] delta_f_map: a hash map for <f, \Delta(f), rnn(F,f)>, which needs to be updated in this function
// [in/out] max_delta_f: the max \Delta(f) for all facilities (with c without f), which needs to be reset in this function
// [in/out] graph: the typed directed graph, which needs to be updated with c as a new facility in this function (fac_id is as 10000 + cand_id)
// [in/out] facs_map: a hash map for storing facilities associated with the endpoints, which needs to be constructed in this function
// [in/out] reflocs_map: a hash map for storing reference locations associated with the probabilities and (dnn - d2nn) * Pr(r), which needs to be updated in this function
// [in/out] refloc_LN_map: a hash map for storing LN edges of reference locations, which needs to be updated in this function
// [in/out] refloc_L2N_map: a hash map for storing L2N edges of reference locations, which needs to be updated in this function
// [in/out] cand_hash_map: a hash map for storing information of candidates (cand_id, vs_id, ve_id, dist_vs, dist_ve), which needs to be updated in this function
// [in/out] obsolete_facs: a set for storing obsolete facilities, which needs to be updated in this function
// [in] reflocs_file_path: file path of reference locations

void LNB_update_LNT_L2NT_Delta(const facility_replacement &op, LNT_hash_map &LNT, L2NT_hash_map &L2NT, delta_f_hash_map &delta_f_map,
	std::pair<int, float> &max_delta_f, typed_directed_graph &graph, fac_hash_map &facs_map, refloc_hash_map &reflocs_map,
	refloc_LN_or_L2N_hash_map &refloc_LN_map, refloc_LN_or_L2N_hash_map &refloc_L2N_map, cand_info_hash_map &cand_hash_map, std::set<int> &obsolete_facs,
	const char *reflocs_file_path)
{
	boost::tuple<int, int, float, float> op_c = cand_hash_map[op.c_id]; // op.c's info, C_VS = 0, C_VE, VS_DIST, VE_DIST
#ifdef ALGO_TESTING_GREEDY
	algo_out_file << op.c_id << ", " << op_c.get<C_VS>() << ", " << op_c.get<C_VE>() << '\n';
#endif
	std::pair<int, int> vs_ve = std::make_pair(op_c.get<C_VS>(), op_c.get<C_VE>()), // <vs, ve> that op.c locates on
		ve_vs = std::make_pair(op_c.get<C_VS>(), op_c.get<C_VE>()); // <ve, vs> for the inverse direction
	cand_hash_map.erase(op.c_id); // remove op.c from candidates
	int new_fac_id = op.c_id + 10000;
	delta_f_map[new_fac_id] = std::make_pair(0.0f, std::set<int>()); // op.c as another facility
#ifdef ALGO_TESTING_GREEDY
	algo_out_file << new_fac_id << ", " << delta_f_map[new_fac_id].first << ", " << delta_f_map[new_fac_id].second.size() << '\n';
#endif
	graph.insert_vertex(vertex(op_c.get<C_VS>(), 'v'), vertex(op_c.get<C_VE>(), 'v'), vertex(new_fac_id, 'f'),
		op_c.get<VS_DIST>(), op_c.get<VE_DIST>()); // update graph with op.c as a new facility
	facs_map[new_fac_id] = boost::make_tuple(op_c.get<C_VS>(), op_c.get<C_VE>(), op_c.get<VS_DIST>(), op_c.get<VE_DIST>());
	facs_map.erase(op.f_id);
	obsolete_facs.insert(op.f_id); // another obsolete facility

	// 1£©for Delta(f) we will directly rerun at present,if f=op.f we will skip not execute,
	//    when f=2nn ,puting edges in LNT and L2NT,we just need to deal with top-vistor.
	// 2) we don't conduct the situation of 2NN == op.f,replacement facilities are some ref loc'2NN,
	// future,we wil add r2nn(F,r) to Delta(f),when we deal with Delta(f),we will manage L2N(r) that may change by dealing with rnn(F,r)'method.
	// Considering that the greedy algorithm is an approximate solution at the moment, we have not yet processed this one.
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// op.c locates on LN of some reference locations
	std::set<int> refloc_LN_vs_ve; // all the reference locations whose LNs overlaps <vs, ve>, then no need to handle <ve, vs> for these reference locations
	bool is_from_vs_to_ve = true; // indicate we consider vs->ve first, which controls the "break (false)" of "do-while"
	LNT_hash_map::iterator iter_vs_ve = LNT.find(vs_ve); // vs->ve
#ifdef ALGO_TESTING_GREEDY
	algo_out_file << vs_ve.first << ", " << vs_ve.second << '\n';
#endif
	do
	{
		if (iter_vs_ve != LNT.end()) // vs->ve exists in LNT
		{
			boost::unordered_map<int, float> refloc_Rd_map; // <refloc_id, Rd> for each reference location whose LN overlaps <vs, ve>
			std::vector<ref_loc_entry>::iterator iter_r_entry = iter_vs_ve->second.second.begin(), iter_r_entry_end = iter_vs_ve->second.second.end();
			for (; iter_r_entry != iter_r_entry_end; ++iter_r_entry) // iterate each reference location whose LN overlaps <vs, ve>
			{
				int the_refloc = iter_r_entry->get<REFLOC>(); // the reference location currently handled
				if (is_from_vs_to_ve) // vs->ve, 1st time into the "do-while"
					refloc_LN_vs_ve.insert(the_refloc); // the reference location has been handled
				else // ve->vs, 2nd time into the "do-while",
					if (refloc_LN_vs_ve.find(the_refloc) != refloc_LN_vs_ve.end()) // and the reference location has been handled when vs->ve
						continue; // go to next reference location

				// first calculate Rd = Dd- + d(c, vj) = dnn - d(r, c)
				if (op_c.get<VE_DIST>() < iter_r_entry->get<OFFSET>()) // the condition that offset takes effect
				{
					// utilizing virtual candidate c', which is the mapping of c with respect to lc, to compute Dd- + d(c, vj)
					float dist_v_c_ve = 2.0f * iter_r_entry->get<OFFSET>() - op_c.get<VE_DIST>(); // d(c', ve) = 2 * offset - d(c, ve)
					refloc_Rd_map[the_refloc] = (iter_r_entry->get<DDL>() + dist_v_c_ve); // Rd = Dd- + d(c, vj)
				}
				else
				{
					float Rd = iter_r_entry->get<DDL>() + op_c.get<VE_DIST>(); // Dd- + d(c, ve)
					if (Rd <= 0) // c locates outside local network (<), maybe locate on L2N
						// meaningless (==), this case means op.c overlap a facility
						continue; // consider next reference location whose LN overlaps <vs, ve>
					else // Rd > 0
					{
						if (Rd > iter_r_entry->get<DDU>()) // > Dd+, must consider reverse direction
						{
							if (is_from_vs_to_ve) // vs->ve, 1st time into the "do-while"
								refloc_LN_vs_ve.erase(the_refloc); // reference location needs to be handled in reverse direction <ve, vs>
							continue; // the reference location will be considered for ve->vs
						}
						else
							refloc_Rd_map[the_refloc] = Rd;
					}
				}
			}

			// handle each reference location whose LN overlaps <vs, ve>
			for (boost::unordered_map<int, float>::iterator iter_r_Rd = refloc_Rd_map.begin(); iter_r_Rd != refloc_Rd_map.end(); ++iter_r_Rd)
			{
				int the_refloc = iter_r_Rd->first; // the reference location currently handled
				float Rd = iter_r_Rd->second; // Rd = dnn - d(r, c) = Dd - +d(c, vj)

				// update \Delta(f) for c (and old NN if necessary), refloc_hash_map, LNT and L2NT
				if (reflocs_map[the_refloc].get<NN>() == op.f_id) // nn(F,r) == op.f, i.e., NN is exactly the replaced facility
				{
					// no need to update \Delta(NN_old == op.f_id) value as the facility will be removed, only delete the reference location which won't be affected by old NN
					int nn_old = reflocs_map[the_refloc].get<NN>(); // old NN is obsolete
					delta_f_map[nn_old].second.erase(the_refloc); // remove the reference location from been affected by old NN

					// update op.c as a new facility in \Delta(f)
					float diff = Rd * reflocs_map[the_refloc].get<PROB>(); // (dnn - d(r, c)) * Pr(r)
					float diff_c_2nn = -(reflocs_map[the_refloc].get<DIFF>()) - diff; // -((d2nn - dnn) * Pr(r)) - (dnn - d(r, c)) * Pr(r) = (d(r, c) - d2nn) * Pr(r)
					delta_f_map[new_fac_id].first += diff_c_2nn;
					delta_f_map[new_fac_id].second.insert(the_refloc); // the reference location is affected by new facility op.c now
#ifdef ALGO_TESTING_GREEDY
					algo_out_file << new_fac_id << ", " << delta_f_map[new_fac_id].first << ", " << delta_f_map[new_fac_id].second.size() << '\n';
#endif

					// update reference location info in refloc_hash_map
					reflocs_map[the_refloc].get<NN>() = new_fac_id; // op.c as new NN
					reflocs_map[the_refloc].get<DIFF>() = -diff_c_2nn; // (d2nn - d(r, c)) * Pr(r)

					// only update L2N(r) with new Dd+ and Dd- values
					std::set<std::pair<int, int>>::iterator iter_L2N = refloc_L2N_map[the_refloc].begin(), iter_L2N_end = refloc_L2N_map[the_refloc].end();
					for (; iter_L2N != iter_L2N_end; ++iter_L2N) // iterate each <vi, vj> in L2N(r)
					{
						std::vector<ref_loc_2_entry>::iterator iter_L2N_r = L2NT[*iter_L2N].begin(), iter_L2N_r_end = L2NT[*iter_L2N].end();
						for (; iter_L2N_r != iter_L2N_r_end; ++iter_L2N_r) // search the entry of currently handled reference location in <vi, vj> edge of L2NT
						{
							if (iter_L2N_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								// update Dd+ and Dd-
								iter_L2N_r->get<DDU_2>() -= Rd;
								iter_L2N_r->get<DDL_2>() -= Rd;
								break; // go on to the next <vi, vj> in L2N(r)
							}
						}
					}

					// move some parts of LN(r) into L2N(r), update other parts with new Dd+ and Dd-
					std::set<std::pair<int, int>> need_to_delete_from_LN; // as we iterate refloc_LN_map[the_refloc], we cannot modify it during the iteration
					std::set<std::pair<int, int>>::iterator iter_LN = refloc_LN_map[the_refloc].begin(), iter_LN_end = refloc_LN_map[the_refloc].end();
					for (; iter_LN != iter_LN_end; ++iter_LN) // iterate each <vi, vj> in LN(r)
					{
						std::vector<ref_loc_entry>::iterator iter_LNT_r = LNT[*iter_LN].second.begin(), iter_LNT_r_end = LNT[*iter_LN].second.end();
						for (; iter_LNT_r != iter_LNT_r_end; ++iter_LNT_r) // search the entry of currently handled reference location in <vi, vj> edge of LNT
						{
							if (iter_LNT_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								float Dd_u = iter_LNT_r->get<DDU>() - Rd; // new Dd+
								if (Dd_u <= 0) // move the entry from LNT to L2NT
								{
									if (refloc_L2N_map[the_refloc].find(*iter_LN) == refloc_L2N_map[the_refloc].end()) // <vi, vj> is not exist in L2N(r), then add to L2NT
									{
										L2NT_hash_map::iterator iter_L2NT = L2NT.find(*iter_LN);
										if (iter_L2NT != L2NT.end()) // edge <vi_id, vj_id> has already been in L2NT, push_back the entry
										{
											iter_L2NT->second.push_back(boost::make_tuple(the_refloc,
												Dd_u, // Dd+
												iter_LNT_r->get<DDL>() - Rd, // Dd-
												iter_LNT_r->get<OFFSET>())); // offset
										}
										else // <vs_id, ve_id> first time presents in L2NT
										{
											L2NT[*iter_LN] = std::vector<ref_loc_2_entry>();
											L2NT[*iter_LN].push_back(boost::make_tuple(the_refloc,
												Dd_u, // Dd+
												iter_LNT_r->get<DDL>() - Rd, // Dd-
												iter_LNT_r->get<OFFSET>())); // offset
										}
										refloc_L2N_map[the_refloc].insert(*iter_LN); // add <vi, vj> as L2N(r)
									}

									// delete entry from LNT
									LNT[*iter_LN].first -= iter_LNT_r->get<DDU>() * reflocs_map[the_refloc].get<PROB>(); // \Delta+(<vi, vj>) -= Dd+ * Pr(r)
									LNT[*iter_LN].second.erase(iter_LNT_r);
									need_to_delete_from_LN.insert(*iter_LN); // later, remove <vi, vj> from LN(r)
								}
								else // still in LNT, only update entry values
								{
									LNT[*iter_LN].first -= diff; // \Delta+(<vi, vj>) -= (dnn - d(r, c)) * Pr(r)
									iter_LNT_r->get<DDU>() = Dd_u;
									iter_LNT_r->get<DDL>() -= Rd;
								}
								break; // go on to the next <vi, vj> in LN(r)
							}
						} //~ search the entry of currently handled reference location in <vi, vj> edge of LNT
					} //~ iterate each <vi, vj> in LN(r)

					// now, do the remove of some <vi, vj> from LN(r)
					for (std::set<std::pair<int, int>>::iterator iter_remove = need_to_delete_from_LN.begin(); iter_remove != need_to_delete_from_LN.end(); ++iter_remove)
						refloc_LN_map[the_refloc].erase(*iter_remove);
				}
				else // nn(F,r) != op.f, NN is not the replaced facility
				{
					// update \Delta(NN_old) value
					int nn_old = reflocs_map[the_refloc].get<NN>(); // old NN of reference location
					delta_f_map[nn_old].first += reflocs_map[the_refloc].get<DIFF>(); // cut the value related to the_refloc
					delta_f_map[nn_old].second.erase(the_refloc); // remove the reference location from been affected by old NN

					// update op.c as a new facility in \Delta(f)
					float diff = Rd * reflocs_map[the_refloc].get<PROB>(); // (dnn - d(r, c)) * Pr(r);
					delta_f_map[new_fac_id].first += -diff; // (d(r, c) - dnn)) * Pr(r), dnn is d2nn now
					delta_f_map[new_fac_id].second.insert(the_refloc); // the reference location is affected by new facility op.c now
#ifdef ALGO_TESTING_GREEDY
					algo_out_file << new_fac_id << ", " << delta_f_map[new_fac_id].first << ", " << delta_f_map[new_fac_id].second.size() << '\n';
#endif

					// update reference location info in refloc_hash_map
					reflocs_map[the_refloc].get<SNN>() = nn_old; // NN as new 2NN
					reflocs_map[the_refloc].get<NN>() = new_fac_id; // op.c as new NN
					reflocs_map[the_refloc].get<DIFF>() = -diff;

					// delete L2N(r) from L2NT, as edges between d(r,c) and dnn are as new L2N(r)
					std::set<std::pair<int, int>>::iterator iter_L2N = refloc_L2N_map[the_refloc].begin(), iter_L2N_end = refloc_L2N_map[the_refloc].end();
					for (; iter_L2N != iter_L2N_end; ++iter_L2N) // iterate each <vi, vj> in L2N(r)
					{
						std::vector<ref_loc_2_entry>::iterator iter_L2N_r = L2NT[*iter_L2N].begin(), iter_L2N_r_end = L2NT[*iter_L2N].end();
						for (; iter_L2N_r != iter_L2N_r_end; ++iter_L2N_r) // search the entry of currently handled reference location in <vi, vj> edge of L2NT
						{
							if (iter_L2N_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								L2NT[*iter_L2N].erase(iter_L2N_r); // delete entry from <vi, vj> in L2NT
								break; // go on to the next <vi, vj> in L2N(r)
							}
						}
					}
					refloc_L2N_map[the_refloc].clear(); // delete L2N(r)

					// update Dd+ and Dd- values in LN(r), and some become in L2N(r)
					std::set<std::pair<int, int>> need_to_delete_from_LN; // as we iterate refloc_LN_map[the_refloc], we cannot modify it during the iteration
					std::set<std::pair<int, int>>::iterator iter_LN = refloc_LN_map[the_refloc].begin(), iter_LN_end = refloc_LN_map[the_refloc].end();
					for (; iter_LN != iter_LN_end; ++iter_LN) // iterate each <vi, vj> in LN(r)
					{
						std::vector<ref_loc_entry>::iterator iter_LN_r = LNT[*iter_LN].second.begin(), iter_LN_r_end = LNT[*iter_LN].second.end();
						for (; iter_LN_r != iter_LN_r_end; ++iter_LN_r) // search the entry of currently handled reference location in <vi, vj> edge of LNT
						{
							if (iter_LN_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								float Dd_u = iter_LN_r->get<DDU>() - Rd; // new Dd+
								if (Dd_u <= 0) // move the entry from LNT to L2NT
								{
									// L2N(r) now is empty, then add directly
									L2NT_hash_map::iterator iter_L2NT = L2NT.find(*iter_LN);
									if (iter_L2NT != L2NT.end()) // edge <vi_id, vj_id> has already been in L2NT, push_back the entry
									{
										iter_L2NT->second.push_back(boost::make_tuple(the_refloc,
											Dd_u, // Dd+
											iter_LN_r->get<DDL>() - Rd, // Dd-
											iter_LN_r->get<OFFSET>())); // offset
									}
									else // <vs_id, ve_id> first time presents in L2NT
									{
										L2NT[*iter_LN] = std::vector<ref_loc_2_entry>();
										L2NT[*iter_LN].push_back(boost::make_tuple(the_refloc,
											Dd_u, // Dd+
											iter_LN_r->get<DDL>() - Rd, // Dd-
											iter_LN_r->get<OFFSET>())); // offset
									}
									refloc_L2N_map[the_refloc].insert(*iter_LN); // add <vi, vj> as L2N(r)

									// delete entry from LNT
									LNT[*iter_LN].first -= iter_LN_r->get<DDU>() * reflocs_map[the_refloc].get<PROB>(); // \Delta+(<vi, vj>) -= Dd+ * Pr(r)
									LNT[*iter_LN].second.erase(iter_LN_r);
									need_to_delete_from_LN.insert(*iter_LN); // later, remove <vi, vj> from LN(r)
								}
								else // still in LNT, only update entry values
								{
									LNT[*iter_LN].first -= diff; // \Delta+(<vi, vj>) -= (dnn - d(r, c)) * Pr(r)
									iter_LN_r->get<DDU>() = Dd_u;
									iter_LN_r->get<DDL>() -= Rd;
								}
								break; // go on to the next <vi, vj> in LN(r)
							}
						} //~ search the entry of currently handled reference location in <vi, vj> edge of LNT
					} //~ iterate each <vi, vj> in LN(r)

					// now, do the remove of some <vi, vj> from LN(r)
					for (std::set<std::pair<int, int>>::iterator iter_remove = need_to_delete_from_LN.begin(); iter_remove != need_to_delete_from_LN.end(); ++iter_remove)
						refloc_LN_map[the_refloc].erase(*iter_remove);
				} //~ nn(F,r) != op.f, NN is not the replaced facility
			} //~ iterate each reference location whose LN overlaps <vs, ve>
		} //~ vs->ve exists in LNT

		if (is_from_vs_to_ve) // 1st time into "do-while", then prepare for the reverse direction, i.e., ve->vs
		{
			is_from_vs_to_ve = false; // indicate we change to consider ve->vs now
			iter_vs_ve = LNT.find(ve_vs); // after this line, become ve->vs!
		}
		else // ve->vs, i.e., 2st time into "do-while"
			break; // out of the "do-while"
	} while (true);

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// op.c locates on L2N of some reference locations
	std::set<int> refloc_L2N_vs_ve; // all the reference locations whose L2Ns overlaps <vs, ve>, then no need to handle <ve, vs> for these reference locations
	is_from_vs_to_ve = true; // again, we consider vs->ve first, which controls the "break (false)" of "do-while"
	L2NT_hash_map::iterator iter_L2NT_vs_ve = L2NT.find(vs_ve); // vs->ve
	do
	{
		if (iter_L2NT_vs_ve != L2NT.end()) // vs->ve exists in L2NT
		{
			boost::unordered_map<int, float> refloc_Rd_map; // <refloc_id, Rd> for each reference location whose L2N overlaps <vs, ve>
			std::vector<ref_loc_2_entry>::iterator iter_r_entry = iter_L2NT_vs_ve->second.begin(), iter_r_entry_end = iter_L2NT_vs_ve->second.end();
			for (; iter_r_entry != iter_r_entry_end; ++iter_r_entry) // iterate each reference location whose L2N overlaps <vs, ve>
			{
				int the_refloc = iter_r_entry->get<REFLOC>(); // the reference location currently handled
				if (is_from_vs_to_ve) // vs->ve, 1st time into the "do-while" (L2NT)
					refloc_L2N_vs_ve.insert(the_refloc); // the reference location has been handled
				else // ve->vs, 2nd time into the "do-while" (L2NT),
					if (refloc_L2N_vs_ve.find(the_refloc) != refloc_L2N_vs_ve.end()) // and the reference location has been handled when vs->ve
						continue; // go to next reference location

				// first calculate Rd = Dd- + d(c, vj) = dnn - d(r, c)
				if (op_c.get<VE_DIST>() < iter_r_entry->get<OFFSET_2>()) // the condition that offset takes effect
				{
					// utilizing virtual candidate c', which is the mapping of c with respect to lc, to compute Dd- + d(c, vj)
					float dist_v_c_ve = 2.0f * iter_r_entry->get<OFFSET_2>() - op_c.get<VE_DIST>(); // d(c', ve) = 2 * offset - d(c, ve)
					refloc_Rd_map[the_refloc] = (iter_r_entry->get<DDL_2>() + dist_v_c_ve); // Dd- + d(c, vj)
				}
				else
				{
					float Rd = iter_r_entry->get<DDL_2>() + op_c.get<VE_DIST>(); // Dd- + d(c, ve)
					if (Rd > 0 && Rd > iter_r_entry->get<DDU_2>()) // > Dd+, must consider reverse direction
					{
						if (is_from_vs_to_ve) // vs->ve, 1st time into the "do-while" (L2NT)
							refloc_L2N_vs_ve.erase(the_refloc); // reference location needs to be handled in reverse direction <ve, vs>
						continue; // the reference location will be considered for ve->vs
					}
					else
						refloc_Rd_map[the_refloc] = Rd;
				}
			}

			// handle each reference location whose L2N overlaps <vs, ve>
			for (boost::unordered_map<int, float>::iterator iter_r_Rd = refloc_Rd_map.begin(); iter_r_Rd != refloc_Rd_map.end(); ++iter_r_Rd)
			{
				int the_refloc = iter_r_Rd->first; // the reference location currently handled
				float Rd = iter_r_Rd->second; // Rd = dnn - d(r, c) = Dd - +d(c, vj)

				// update \Delta(f) for c (and old NN if necessary), refloc_hash_map, LNT and L2NT
				if (reflocs_map[the_refloc].get<NN>() == op.f_id) // nn(F,r) == op.f, i.e., NN is exactly the replaced facility
				{
					// no need to update \Delta(NN_old == op.f_id) value as the facility will be removed, only delete the reference location which won't be affected by old NN
					int nn_old = reflocs_map[the_refloc].get<NN>(); // old NN is obsolete
					delta_f_map[nn_old].second.erase(the_refloc); // remove the reference location from been affected by old NN

					// update op.c as a new facility in \Delta(f), as old NN will be replaced
					float diff = Rd * reflocs_map[the_refloc].get<PROB>(); // (dnn - d(r, c)) * Pr(r)
					float diff_c_2nn = -(reflocs_map[the_refloc].get<DIFF>()) - diff; // -((d2nn - dnn) * Pr(r)) - (dnn - d(r, c)) * Pr(r) = (d(r, c) - d2nn) * Pr(r)
					delta_f_map[new_fac_id].first += diff_c_2nn;
					delta_f_map[new_fac_id].second.insert(the_refloc); // the reference location is affected by new facility op.c now
#ifdef ALGO_TESTING_GREEDY
					algo_out_file << new_fac_id << ", " << delta_f_map[new_fac_id].first << ", " << delta_f_map[new_fac_id].second.size() << '\n';
#endif

					// update reference location info in refloc_hash_map
					reflocs_map[the_refloc].get<NN>() = new_fac_id; // op.c as new NN
					reflocs_map[the_refloc].get<DIFF>() = -diff_c_2nn; // (d2nn - d(r, c)) * Pr(r)

					// only update LN(r) with new Dd+ and Dd- values
					std::set<std::pair<int, int>>::iterator iter_LN = refloc_LN_map[the_refloc].begin(), iter_LN_end = refloc_LN_map[the_refloc].end();
					for (; iter_LN != iter_LN_end; ++iter_LN) // iterate each <vi, vj> in LN(r)
					{
						std::vector<ref_loc_entry>::iterator iter_LN_r = LNT[*iter_LN].second.begin(), iter_LN_r_end = LNT[*iter_LN].second.end();
						for (; iter_LN_r != iter_LN_r_end; ++iter_LN_r) // search the entry of currently handled reference location in <vi, vj> edge of LNT
						{
							if (iter_LN_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								// update Dd+ and Dd-
								iter_LN_r->get<DDU>() -= Rd;
								iter_LN_r->get<DDL>() -= Rd;
								break; // go on to the next <vi, vj> in LN(r)
							}
						}
					}

					// move some parts of L2N(r) into LN(r), update other parts in L2N(r) with new Dd+ and Dd-
					std::set<std::pair<int, int>> need_to_delete_from_L2N; // as we iterate refloc_L2N_map[the_refloc], we cannot modify it during the iteration
					std::set<std::pair<int, int>>::iterator iter_L2N = refloc_L2N_map[the_refloc].begin(), iter_L2N_end = refloc_L2N_map[the_refloc].end();
					for (; iter_L2N != iter_L2N_end; ++iter_L2N) // iterate each <vi, vj> in L2N(r)
					{
						std::vector<ref_loc_2_entry>::iterator iter_L2NT_r = L2NT[*iter_L2N].begin(), iter_L2NT_r_end = L2NT[*iter_L2N].end();
						for (; iter_L2NT_r != iter_L2NT_r_end; ++iter_L2NT_r) // search the entry of currently handled reference location in <vi, vj> edge of L2NT
						{
							if (iter_L2NT_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								float Dd_u = iter_L2NT_r->get<DDU_2>() - Rd; // new Dd+, as Rd = dnn - d(r, c) < 0, then new Dd+ > old Dd+
								if (Dd_u > 0) // move the entry from L2NT to LNT
								{
									if (refloc_LN_map[the_refloc].find(*iter_L2N) == refloc_LN_map[the_refloc].end()) // <vi, vj> is not exist in LN(r), then add to LNT
									{
										LNT_hash_map::iterator iter_LNT = LNT.find(*iter_L2N);
										if (iter_LNT != LNT.end()) // edge <vi_id, vj_id> has already been in LNT, push_back the entry
										{
											iter_LNT->second.first += Dd_u * reflocs_map[the_refloc].get<PROB>();
											iter_LNT->second.second.push_back(boost::make_tuple(the_refloc,
												Dd_u, // Dd+
												iter_L2NT_r->get<DDL_2>() - Rd, // Dd-
												iter_L2NT_r->get<OFFSET_2>())); // offset
										}
										else // <vs_id, ve_id> first time presents in LNT
										{
											LNT[*iter_L2N] = std::make_pair(Dd_u * reflocs_map[the_refloc].get<PROB>(), // Dd+ as initial ED+
												std::vector<ref_loc_entry>());
											LNT[*iter_L2N].second.push_back(boost::make_tuple(the_refloc,
												Dd_u, // Dd+
												iter_L2NT_r->get<DDL_2>() - Rd, // Dd-
												iter_L2NT_r->get<OFFSET_2>())); // offset
										}
										refloc_LN_map[the_refloc].insert(*iter_L2N); // add <vi, vj> as LN(r)
									}

									// delete entry from L2NT
									L2NT[*iter_L2N].erase(iter_L2NT_r);
									need_to_delete_from_L2N.insert(*iter_L2N); // later, remove <vi, vj> from L2N(r)
								}
								else // still in L2NT, only update entry values
								{
									iter_L2NT_r->get<DDU_2>() = Dd_u;
									iter_L2NT_r->get<DDL_2>() -= Rd;
								}
								break; // go on to the next <vi, vj> in LN(r)
							}
						}
					}

					// now, do the remove of some <vi, vj> from L2N(r)
					for (std::set<std::pair<int, int>>::iterator iter_remove = need_to_delete_from_L2N.begin(); iter_remove != need_to_delete_from_L2N.end(); ++iter_remove)
						refloc_L2N_map[the_refloc].erase(*iter_remove);
				}
				else // nn(F,r) != op.f, NN is not the replaced facility
				{
					// update \Delta(NN) value
					int nn = reflocs_map[the_refloc].get<NN>(); // NN is not changed
					float diff = Rd * reflocs_map[the_refloc].get<PROB>(); // (dnn - d(r, c)) * Pr(r);
					delta_f_map[nn].first += reflocs_map[the_refloc].get<DIFF>() + diff; // (dnn - d2nn) * Pr(r) -> (dnn - d(r, c)) * Pr(r)

					// no need to add op.c (2NN) to be a new facility in \Delta(f)

					// update reference location info in refloc_hash_map
					reflocs_map[the_refloc].get<SNN>() = new_fac_id; // op.c as new 2NN
					reflocs_map[the_refloc].get<DIFF>() = -diff;

					// no need to update LN(r) (and LNT)

					// delete parts of L2N(r) which are out of range <r, op.c>
					std::set<std::pair<int, int>> need_to_delete_from_L2N; // as we iterate refloc_L2N_map[the_refloc], we cannot modify it during the iteration
					std::set<std::pair<int, int>>::iterator iter_L2N = refloc_L2N_map[the_refloc].begin(), iter_L2N_end = refloc_L2N_map[the_refloc].end();
					for (; iter_L2N != iter_L2N_end; ++iter_L2N) // iterate each <vi, vj> in L2N(r)
					{
						std::vector<ref_loc_2_entry>::iterator iter_L2NT_r = L2NT[*iter_L2N].begin(), iter_L2NT_r_end = L2NT[*iter_L2N].end();
						for (; iter_L2NT_r != iter_L2NT_r_end; ++iter_L2NT_r) // search the entry of currently handled reference location in <vi, vj> edge of L2NT
						{
							if (iter_L2NT_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
							{
								if (iter_L2NT_r->get<DDU_2>() <= Rd) // dnn - d(r, vi) <= dnn - d(r, c), i.e., d(r, vi) >= d(r, c), which means out of L2N(r)
								{
									L2NT[*iter_L2N].erase(iter_L2NT_r); // delete entry from <vi, vj> in L2NT
									need_to_delete_from_L2N.insert(*iter_L2N); // later, delete <vi, vj> from L2N(r)
								}
								break; // go on to the next <vi, vj> in L2N(r)
							}
						}
					}

					// now, do the remove of some <vi, vj> from L2N(r)
					for (std::set<std::pair<int, int>>::iterator iter_remove = need_to_delete_from_L2N.begin(); iter_remove != need_to_delete_from_L2N.end(); ++iter_remove)
						refloc_L2N_map[the_refloc].erase(*iter_remove);					
				} //~ nn(F,r) != op.f, NN is not the replaced facility
			} //~ iterate each reference location whose L2N overlaps <vs, ve>
		} //~ vs->ve exists in L2NT

		if (is_from_vs_to_ve) // 1st time into "do-while", then prepare for the reverse direction, i.e., ve->vs
		{
			is_from_vs_to_ve = false; // indicate we change to consider ve->vs now
			iter_L2NT_vs_ve = L2NT.find(ve_vs); // after this line, become ve->vs!
		}
		else // ve->vs, i.e., 2st time into "do-while"
			break; // out of the "do-while"
	} while (true);

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// clear \Delta(op.f), LN(r), L2N(r) and reflocs_map[r], where r in rnn(F, op.f)
	std::set<int> need_reconstruct; // reference locations whose LNT and L2NT entries need to be re-constructed
	std::set<int>::iterator iter_rnn_r = delta_f_map[op.f_id].second.begin(), iter_rnn_r_end = delta_f_map[op.f_id].second.end();
	for (; iter_rnn_r != iter_rnn_r_end; ++iter_rnn_r) // iterate each reference location in rnn(F, op.f)
	{
		int the_refloc = *iter_rnn_r;
		need_reconstruct.insert(the_refloc); // record reference locations which need to be re-constructed

		// delete entries in LNT and LN(r)
		std::set<std::pair<int, int>>::iterator iter_LN = refloc_LN_map[the_refloc].begin(), iter_LN_end = refloc_LN_map[the_refloc].end();
		for (; iter_LN != iter_LN_end; ++iter_LN) // iterate <vs, ve> in LN(r)
		{
			std::vector<ref_loc_entry>::iterator iter_LNT_r = LNT[*iter_LN].second.begin(), iter_LNT_r_end = LNT[*iter_LN].second.end();
			for (; iter_LNT_r != iter_LNT_r_end; ++iter_LNT_r) // search the entry of currently handled reference location in <vs, ve> edge of LNT
			{
				if (iter_LNT_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
				{
					LNT[*iter_LN].first += reflocs_map[the_refloc].get<DIFF>(); // cut the value related to the_refloc
					LNT[*iter_LN].second.erase(iter_LNT_r); // delete entry from <vs, ve> in LNT
					break; // go on to the next <vs, ve> of LN(r)
				}
			}
		}
		refloc_LN_map.erase(the_refloc); // delete <refloc_id, LN(r)>

		// delete entries in L2NT and L2N(r)
		std::set<std::pair<int, int>>::iterator iter_L2N = refloc_L2N_map[the_refloc].begin(), iter_L2N_end = refloc_L2N_map[the_refloc].end();
		for (; iter_L2N != iter_L2N_end; ++iter_L2N) // iterate <vs, ve> in L2N(r)
		{
			std::vector<ref_loc_2_entry>::iterator iter_L2NT_r = L2NT[*iter_L2N].begin(), iter_L2NT_r_end = L2NT[*iter_L2N].end();
			for (; iter_L2NT_r != iter_L2NT_r_end; ++iter_L2NT_r) // search the entry of currently handled reference location in <vs, ve> edge of L2NT
			{
				if (iter_L2NT_r->get<REFLOC>() == the_refloc) // currently handled reference location entry is found
				{
					L2NT[*iter_L2N].erase(iter_L2NT_r); // delete entry from <vs, ve> in L2NT
					break; // go on to the next <vs, ve> of L2N(r)
				}
			}
		}
		refloc_L2N_map.erase(the_refloc); // delete <refloc_id, L2N(r)>
	}

	delta_f_map.erase(op.f_id); // remove the obsolete facility
	for (std::set<int>::iterator iter_r = need_reconstruct.begin(); iter_r != need_reconstruct.end(); ++iter_r)
		reflocs_map.erase(*iter_r); // remove refloc information, and these information will be re-construct in the following

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// re-construct partly LNT, L2NT, \Delta(f) and reflocs_map
	// create an extend event visitor for LNB algorithm
	LNB_extend_visitor extend_visitor;
	loc_edges_hash_map refloc_loc_edges_hash_map; // an assistant hash map for local network edges
	extend_visitor.set_data_structure(&facs_map, &refloc_loc_edges_hash_map); // set data structures

	// create a top event visitor for LNB algorithm
	LNB_top_visitor top_visitor;
	top_visitor.set_data_structures(&LNT, &L2NT, &delta_f_map, &reflocs_map, &refloc_loc_edges_hash_map,
		NULL, // vertex_dnn_hash_map has been set when the 1st time we run a LNB_top_visitor
		&refloc_LN_map, &refloc_L2N_map); // set data structures
	if (!obsolete_facs.empty())
		top_visitor.set_obsolete_facs(&obsolete_facs); // set obsolete facilities

	// scan reference locations, and deal with those need to be re-constructed
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!need_reconstruct.empty() && !ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			
			if (need_reconstruct.find(refloc_id) == need_reconstruct.end()) // this reference location needs no re-construction
				continue;
			need_reconstruct.erase(refloc_id); // this reference location will be re-constructed

			std::string::size_type ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos;
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);

			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));

			// insert the reference location into refloc_LN_map & refloc_L2N_map
			refloc_LN_map[refloc_id] = std::set< std::pair <int, int> >();
			refloc_L2N_map[refloc_id] = std::set< std::pair <int, int> >();

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			extend_visitor.set_vs_ve(vs_id, ve_id, dist_vs, dist_ve); // the original starting and ending vertices of the reference location
			top_visitor.set_refloc(refloc_id, prob); // set the reference location and its present probability, also reset nearest facility
			refloc_loc_edges_hash_map.clear(); // must reset the assistant hash map for each new reference location
			if (vs_id == ve_id)
				top_visitor.set_overlap(vs_id); // reference location overlaps a vertex
			else
				top_visitor.set_overlap(-1); // reference location doesn't overlap a vertex
			graph.dijkstra(refloc_v, &top_visitor, &extend_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// re-evaluate the max \Delta(f)
	max_delta_f.first = -1; // facility id, null
	max_delta_f.second = -FLT_MAX; // \Delta(f), min negative float
	for (delta_f_hash_map::const_iterator iter_delta = delta_f_map.cbegin(); iter_delta != delta_f_map.cend(); ++iter_delta)
	{
		if (iter_delta->second.first > max_delta_f.second) // \Delta(f) > max_delta_f
		{
			max_delta_f.first = iter_delta->first; // f
			max_delta_f.second = iter_delta->second.first; // \Delta(f)
		}
	}
}

#pragma region CONS_LNT
//--------------------------------------------------------------------------------
// LNB algorithm: construct LNT and L2NT
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [out] LNT: a hash map as Local Network Table, which needs to be constructed in this function
// [out] L2NT: a hash map as Local Sub-Network Table, which needs to be constructed in this function
// [out] dnn_map; a hash map for reference locations that overlap some vertices, which needs to be constructed in this function
// [out] delta_f_map: a hash map for <f, \Delta(f), rnn(F,f)>, which needs to be constructed in this function
// [out] max_delta_f: the max \Delta(f) for all facilities, which needs to be set in this function
// [out] graph: the typed directed graph, which needs to be constructed in this function
// [out] facs_map: a hash map for storing facilities associated with the endpoints, which needs to be constructed in this function
// [out] reflocs_map: a hash map for storing reference locations associated with the probabilities and (dnn - d2nn) * Pr(r), which needs to be constructed in this function
// [out] refloc_LN_map: a hash map for storing LN edges of reference locations, which needs to be constructed in this function
// [out] refloc_L2N_map: a hash map for storing L2N edges of reference locations, which needs to be constructed in this function

void LNB_construct_LNT(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path,
	LNT_hash_map &LNT, L2NT_hash_map &L2NT, vertex_dnn_hash_map &dnn_map, delta_f_hash_map &delta_f_map, std::pair<int, float> &max_delta_f,
	typed_directed_graph &graph, fac_hash_map &facs_map, refloc_hash_map &reflocs_map,
	refloc_LN_or_L2N_hash_map &refloc_LN_map, refloc_LN_or_L2N_hash_map &refloc_L2N_map)
{
#ifdef ALGO_TESTING_TRACE_GRAPH
	graph.set_testing(&algo_out_file);
#endif

	// init graph, and deploy facilities
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);

	// deploy facilities and construct related data structure
	graph.deploy_facilities(facs_file_path, &facs_map);

	if (!graph.get_is_facilities_deployed()) // check graph and facilities status
		return;

	// create an extend event visitor for LNB algorithm
	LNB_extend_visitor extend_visitor;
	loc_edges_hash_map refloc_loc_edges_hash_map; // an assistant hash map for local network edges
	extend_visitor.set_data_structure(&facs_map, &refloc_loc_edges_hash_map); // set data structures

	// create a top event visitor for LNB algorithm
	LNB_top_visitor top_visitor;
	top_visitor.set_data_structures(&LNT, &L2NT, &delta_f_map, &reflocs_map, &refloc_loc_edges_hash_map, &dnn_map,
		&refloc_LN_map, &refloc_L2N_map); // set data structures

	// deal with each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);

			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			// lon and lat are useless for network traversal

			// insert the reference location into refloc_LN_map & refloc_L2N_map
			refloc_LN_map[refloc_id] = std::set< std::pair <int, int> >();
			refloc_L2N_map[refloc_id] = std::set< std::pair <int, int> >();

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			extend_visitor.set_vs_ve(vs_id, ve_id, dist_vs, dist_ve); // the original starting and ending vertices of the reference location
			top_visitor.set_refloc(refloc_id, prob); // set the reference location and its present probability, also reset nearest facility
			refloc_loc_edges_hash_map.clear(); // must reset the assistant hash map for each new reference location
			if (vs_id == ve_id)
				top_visitor.set_overlap(vs_id); // reference location overlaps a vertex
			else
				top_visitor.set_overlap(-1); // reference location doesn't overlap a vertex
			graph.dijkstra(refloc_v, &top_visitor, &extend_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}
	}

	// evaluate the max \Delta(f)
	max_delta_f.first = -1; // facility id, null
	max_delta_f.second = -FLT_MAX; // \Delta(f), min negative float
	for (delta_f_hash_map::const_iterator iter_delta = delta_f_map.cbegin(); iter_delta != delta_f_map.cend(); ++iter_delta)
	{
		if (iter_delta->second.first > max_delta_f.second) // \Delta(f) > max_delta_f
		{
			max_delta_f.first = iter_delta->first; // f
			max_delta_f.second = iter_delta->second.first; // \Delta(f)
		}
	}
}
#pragma endregion CONS_LNT

//--------------------------------------------------------------------------------
// LNB algorithm: query optimal facility-replacement pair based on LNT & L2NT
// remark: the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat";
//		   for simplicity, we use LNT, L2NT, rnn_f_tuples & reflocs_map as non-const parameters (for [] operator), hence, in this function, be careful to use them, and DON'T change their values anytime;
//		   if a reference location has no nearest or sub-nearest facility (namely isolate sub-grahp), we ignore it and regard its contribute on ED as 0
// [re] facility_replacement: the optimal facility-replacement pair; default facility-replacement object means opening candidates file fails
// [out] checked_num: the number of facility-replacement pairs actually checked
// [in] LNT: the constructed hash map as LNT
// [in] L2NT: the constructed hash map as Local Sub-Network Table
// [in] dnn_map; the constructed hash map for reference locations that overlap some vertices
// [in] delta_f_map: a hash map for <f, \Delta(f), rnn(F,f)>
// [in] max_delta_f: the max \Delta(f) for all facilities
// [in] graph: the constructed graph
// [in] facs_map: a hash map for storing facilities associated with the endpoints
// [in] reflocs_map: a hash map for storing reference locations associated with the probabilities and (dnn - d2nn) * Pr(r)
// [in] cands_file_path: file path of candidates
// [out] cands_hash_map: a hash map for storing information of candidates (cand_id, vs_id, ve_id, dist_vs, dist_ve)

facility_replacement LNB_query(int &checked_num, LNT_hash_map &LNT, L2NT_hash_map &L2NT, const vertex_dnn_hash_map &dnn_map, delta_f_hash_map &delta_f_map,
	std::pair<int, float> &max_delta_f, const typed_directed_graph &graph, const fac_hash_map &facs_map, refloc_hash_map &reflocs_map, const char *cands_file_path,
	cand_info_hash_map &cands_map)
{
#ifdef ALGO_TESTING
	algo_out_file << "LNB max_delta_f: <" << max_delta_f.first << ", " << max_delta_f.second << ">\n";
#endif

	cand_max_fibonacci_heap max_heap; // max-heap ordered by ED+ in LNT (LNH in paper)

	boost::unordered_map < int, // candidate id
		std::set < int > > // refloc ids
		Rc_u_sets; // Rc+ set for candidates, each of which overlaps with some vertex
				   // NOTE: this set variable is only for each candidate that overlaps a vertex, whose Rc+ set is intialized when constructing LNH

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// construct a max-heap (LNH) for candidates
	if (!cands_map.empty()) // then Do Not set cands_hash_map, as it has been initialized already
	{
		for (cand_info_hash_map::const_iterator iter_c = cands_map.cbegin(); iter_c != cands_map.cend(); ++iter_c)
		{
			if (iter_c->second.get<C_VS>() != iter_c->second.get<C_VE>()) // candidate splits an edge
				// NOTE: Rc+ set for this candidate c will be constuct when accumulating \Delta(f_any, c)
			{
				// here, we explicitly indicate "vs" and "ve"
				std::pair<int, int> vs_ve(iter_c->second.get<C_VS>(), iter_c->second.get<C_VE>()), // vs->ve
					ve_vs(iter_c->second.get<C_VE>(), iter_c->second.get<C_VS>()); // ve->vs, reverse direction

				// retrieve ED+ of edge
				float ED_u = 0.0f; // ED+ for ordering candidates in max-heap
				LNT_hash_map::const_iterator iter_vs_ve = LNT.find(vs_ve); // vs->ve
				LNT_hash_map::const_iterator iter_ve_vs = LNT.find(ve_vs); // ve->vs
				if (iter_vs_ve != LNT.cend() && iter_ve_vs != LNT.cend()) // vs->ve & ve->vs are both valid
					ED_u = iter_vs_ve->second.first + iter_ve_vs->second.first; // larger than upper
				else if (iter_vs_ve != LNT.cend()) // only vs->ve is valid
					ED_u = iter_vs_ve->second.first;
				else if (iter_ve_vs != LNT.cend()) // only ve->vs is valid
					ED_u = iter_ve_vs->second.first;
				else // neither are valid, ED+ should be 0.0f
					continue; // no need to consider this candidate

				max_heap.push(candidate(iter_c->first, ED_u)); // push each candidate
			}
			else // vs_id == ve_id, candidate overlaps a vertex
				// NOTE: 1) facility cannot overlap the vertex, otherwise, candidate overlaps facility, which is strictly forbidden;
				//			 hence, any network traversal cannot end at the vertex, thus considering only out-edges can avoid double counting,
				//			 where candidate is viewed as an intermediate vertex in a path, namely as vs and ve for two different edges, for computing ERd;
				//		  2) also, if traversal ends immediately before the vertex (i.e., r->f->c), candidate is outside local network and has no positive benefit;
				//		  3) moreover, if a reference location splits an edge with candidate (the vertex) as one of the two endpoints,
				//			 namely out-edge c->r->ve, then Rd of c must > nnd (Dd+), thus it has to consider the reverse direction
				//		  4) conversely, only considering in-edges is not enough, the counter example is reference location overlaps the vertex,
				//			 and for this reference location, there is only out-edge
				//		  5) in the special case, where c overlaps v and some r(s) also overlap v, we compute ED with the help of vertex_nnd_hash_map
				//		  6) as ED will be computed in this case code block, then Rc+ set for this candidate, which overlaps with the vertex, is intialized here
			{
				float ED = 0.0f; // for a candidate that overlaps a vertex, we directly compute its ED, because Ed can be obtained via Dd+ or Dd-
				Rc_u_sets[iter_c->first] = std::set<int>(); // Rc+ set for this candidate is intialized in this case code block

				adjacent_edges_hash_map::const_iterator iter_vertex = graph.adjacent_edges.find(vertex(iter_c->second.get<C_VS>(), 'v')); // the overlapped vertex

				// considering only out-edges is enough; hence, no need to consider in-edges
				std::set<out_edge>::const_iterator iter_out_edge = iter_vertex->second.first.cbegin();
				for (; iter_out_edge != iter_vertex->second.first.cend(); ++iter_out_edge) // compute ERd for each out-edge
				{
					// fault-tolerant check in case out-vertex is a facility
					int out_v_id = iter_out_edge->ve.id;
					float dist_vs_ve = iter_out_edge->dist;
					if (iter_out_edge->ve.type != 'v') // must be 'f'
					{
						fac_hash_map::const_iterator iter_fac = facs_map.find(out_v_id); // as facility exists, this iterator must not be cend()
						if (iter_fac->second.get<VS_ID>() == iter_c->second.get<C_VS>())
							out_v_id = iter_fac->second.get<VE_ID>();
						else
							out_v_id = iter_fac->second.get<VS_ID>();
						dist_vs_ve = iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>(); // edge distance
					}

					LNT_hash_map::const_iterator iter_vs_ve = LNT.find(std::make_pair(iter_c->second.get<C_VS>(), out_v_id)); // vs->ve
					if (iter_vs_ve != LNT.cend()) // vs->ve takes effect on ERd
					{
						std::vector<ref_loc_entry>::const_iterator iter_entry = iter_vs_ve->second.second.cbegin();
						for (; iter_entry != iter_vs_ve->second.second.cend(); ++iter_entry) // each reference location related to this out-edge
						{
							// 1) we view edge distance as d(c, ve), then d(c, ve) must > offset, thus offset has no effect
							// 2) c (say the overlapped vertex) has definite distance from r, hence d(r, vs) <= nnd, thus Rd is impossible < 0 (locates outside local network)
							float Rd = iter_entry->get<DDL>() + dist_vs_ve; // Dd- + d(c, ve) >= 0
							if (Rd > iter_entry->get<DDU>()) // > Dd+, must consider reverse direction
								continue;
							else // c takes effect
							{
								ED += Rd * reflocs_map[iter_entry->get<REFLOC>()].get<PROB>(); // Rd * Pr(r)
								Rc_u_sets[iter_c->first].insert(iter_entry->get<REFLOC>()); // insert refloc into Rc+ set
							}
						}
					}
				}

				// fault-tolerant check in case the candidate overlaps a vertex, which is also overlapped by some reference location(s)
				// P.S.: I have no idea why I wrote these codes below; maybe, it prevent to re-add ED of overlapped reference locations for each out-edge?
				unsigned out_edges_size = static_cast<unsigned>(iter_vertex->second.first.size());
				vertex_dnn_hash_map::const_iterator iter_dnn = dnn_map.find(iter_c->second.get<C_VS>());
				if (iter_dnn != dnn_map.cend()) // some reference location(s) overlaps the vertex
					ED -= (out_edges_size - 1) * iter_dnn->second; // ED = ED - size * accumulated dnn (minus all out-edges) + accumulated dnn (reserve one out-edge)

				max_heap.push(candidate(iter_c->first, ED)); // push candidate and actual ED
			}
		}
	}
	else
	{
		std::ifstream ifs_cands(cands_file_path); // read candidates file
		if (!ifs_cands.fail()) // candidates file exists
		{
			while (!ifs_cands.bad() && ifs_cands.good())
			{
				char buf[1024];
				ifs_cands.getline(buf, 1024);
				std::string str_buf(buf);

				// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
				std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos;
				vs_pos = str_buf.find(' ', begin_pos);
				ve_pos = str_buf.find(' ', vs_pos + 1);
				dist_vs_pos = str_buf.find(' ', ve_pos + 1);
				dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
				lon_pos = str_buf.find(' ', dist_ve_pos + 1);

				int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
				int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
				int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
				float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
				float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, lon_pos - dist_ve_pos - 1).c_str()));
				// lon and lat are useless for query

				if (vs_id != ve_id) // candidate splits an edge
					// NOTE: Rc+ set for this candidate c will be constuct when accumulating \Delta(f_any, c)
				{
					// here, we explicitly indicate "vs" and "ve"
					std::pair<int, int> vs_ve(vs_id, ve_id), // vs->ve
						ve_vs(ve_id, vs_id); // ve->vs, reverse direction

					// retrieve ED+ of edge
					float ED_u = 0.0f; // ED+ for ordering candidates in max-heap
					LNT_hash_map::const_iterator iter_vs_ve = LNT.find(vs_ve); // vs->ve
					LNT_hash_map::const_iterator iter_ve_vs = LNT.find(ve_vs); // ve->vs
					if (iter_vs_ve != LNT.cend() && iter_ve_vs != LNT.cend()) // vs->ve & ve->vs are both valid
						ED_u = iter_vs_ve->second.first + iter_ve_vs->second.first; // larger than upper
					else if (iter_vs_ve != LNT.cend()) // only vs->ve is valid
						ED_u = iter_vs_ve->second.first;
					else if (iter_ve_vs != LNT.cend()) // only ve->vs is valid
						ED_u = iter_ve_vs->second.first;
					else // neither are valid, ED+ should be 0.0f
						continue; // no need to consider this candidate

					max_heap.push(candidate(cand_id, ED_u)); // push each candidate
					cands_map[cand_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve); // record information of each candidate
				}
				else // vs_id == ve_id, candidate overlaps a vertex
					// NOTE: 1) facility cannot overlap the vertex, otherwise, candidate overlaps facility, which is strictly forbidden;
					//			 hence, any network traversal cannot end at the vertex, thus considering only out-edges can avoid double counting,
					//			 where candidate is viewed as an intermediate vertex in a path, namely as vs and ve for two different edges, for computing ERd;
					//		  2) also, if traversal ends immediately before the vertex (i.e., r->f->c), candidate is outside local network and has no positive benefit;
					//		  3) moreover, if a reference location splits an edge with candidate (the vertex) as one of the two endpoints,
					//			 namely out-edge c->r->ve, then Rd of c must > nnd (Dd+), thus it has to consider the reverse direction
					//		  4) conversely, only considering in-edges is not enough, the counter example is reference location overlaps the vertex,
					//			 and for this reference location, there is only out-edge
					//		  5) in the special case, where c overlaps v and some r(s) also overlap v, we compute ED with the help of vertex_nnd_hash_map
					//		  6) as ED will be computed in this case code block, then Rc+ set for this candidate, which overlaps with the vertex, is intialized here
				{
					float ED = 0.0f; // for a candidate that overlaps a vertex, we directly compute its ED, because Ed can be obtained via Dd+ or Dd-
					Rc_u_sets[cand_id] = std::set<int>(); // Rc+ set for this candidate is intialized in this case code block

					adjacent_edges_hash_map::const_iterator iter_vertex = graph.adjacent_edges.find(vertex(vs_id, 'v')); // the overlapped vertex

					// considering only out-edges is enough; hence, no need to consider in-edges
					std::set<out_edge>::const_iterator iter_out_edge = iter_vertex->second.first.cbegin();
					for (; iter_out_edge != iter_vertex->second.first.cend(); ++iter_out_edge) // compute ERd for each out-edge
					{
						// fault-tolerant check in case out-vertex is a facility
						int out_v_id = iter_out_edge->ve.id;
						float dist_vs_ve = iter_out_edge->dist;
						if (iter_out_edge->ve.type != 'v') // must be 'f'
						{
							fac_hash_map::const_iterator iter_fac = facs_map.find(out_v_id); // as facility exists, this iterator must not be cend()
							if (iter_fac->second.get<VS_ID>() == vs_id)
								out_v_id = iter_fac->second.get<VE_ID>();
							else
								out_v_id = iter_fac->second.get<VS_ID>();
							dist_vs_ve = iter_fac->second.get<VS_DIST>() + iter_fac->second.get<VE_DIST>(); // edge distance
						}

						LNT_hash_map::const_iterator iter_vs_ve = LNT.find(std::make_pair(vs_id, out_v_id)); // vs->ve
						if (iter_vs_ve != LNT.cend()) // vs->ve takes effect on ERd
						{
							std::vector<ref_loc_entry>::const_iterator iter_entry = iter_vs_ve->second.second.cbegin();
							for (; iter_entry != iter_vs_ve->second.second.cend(); ++iter_entry) // each reference location related to this out-edge
							{
								// 1) we view edge distance as d(c, ve), then d(c, ve) must > offset, thus offset has no effect
								// 2) c (say the overlapped vertex) has definite distance from r, hence d(r, vs) <= nnd, thus Rd is impossible < 0 (locates outside local network)
								float Rd = iter_entry->get<DDL>() + dist_vs_ve; // Dd- + d(c, ve) >= 0
								if (Rd > iter_entry->get<DDU>()) // > Dd+, must consider reverse direction
									continue;
								else // c takes effect
								{
									ED += Rd * reflocs_map[iter_entry->get<REFLOC>()].get<PROB>(); // Rd * Pr(r)
									Rc_u_sets[cand_id].insert(iter_entry->get<REFLOC>()); // insert refloc into Rc+ set
								}
							}
						}
					}

					// fault-tolerant check in case the candidate overlaps a vertex, which is also overlapped by some reference location(s)
					// P.S.: I have no idea why I wrote these codes below; maybe, it prevent to re-add ED of overlapped reference locations for each out-edge?
					unsigned out_edges_size = static_cast<unsigned>(iter_vertex->second.first.size());
					vertex_dnn_hash_map::const_iterator iter_dnn = dnn_map.find(vs_id);
					if (iter_dnn != dnn_map.cend()) // some reference location(s) overlaps the vertex
						ED -= (out_edges_size - 1) * iter_dnn->second; // ED = ED - size * accumulated dnn (minus all out-edges) + accumulated dnn (reserve one out-edge)

					max_heap.push(candidate(cand_id, ED)); // push candidate and actual ED
					cands_map[cand_id] = boost::make_tuple(vs_id, ve_id, 0.0f, 0.0f); // set overlap flags (i.e., LNT.cend() and 0.0f) for the candidate
				}
			}
		}
		else // opening candidates file fails
			return facility_replacement();
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// iterate all facility-replacements, based on the ordered candidates in LNH, max_delta_f (ie., max{\Delta(f)}) & Rc+ sets
	facility_replacement op(-1, -1, 0.0f); // initialize optimal pair
	checked_num = 0; // initialize
	while (!max_heap.empty()) // LNH
	{
		// retrieve the top candidate in max-heap
		candidate top_c = max_heap.top();
		boost::unordered_map<int, float> delta_f_offset; // <f, \Delta(f)_offset> based on the current top_c
#ifdef ALGO_TESTING
		//algo_out_file << "ED+: c" << top_c.id << ": " << top_c.ERD << '\n'; // for candidate class, "ERD" member name is still used w.r.t. its essence
#endif

#ifndef ALGO_TESTING_NO_EARLY_STOPPING
		if (top_c.ERD < op.ED) // current op pair is the optima
			break;
#endif
		
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// deal with \Delta+(c), the positive effect
		if (cands_map[top_c.id].get<VE_DIST>() != 0	&& cands_map[top_c.id].get<VS_DIST>() != 0) // c splits an edge, but not overlaps a vertex, as overlapping has definite ED+
		{
			top_c.ERD = 0.0f; // prepare to accumulate actual ED+ of this candidate
							  // and Rc+ set will be construct in this "if" case code block

			// iteratively deal with two possible edges vs->ve and ve->vs
			LNT_hash_map::const_iterator iters[2] = { LNT.find(std::make_pair(cands_map[top_c.id].get<C_VS>(), cands_map[top_c.id].get<C_VE>())), // vs->ve
				LNT.find(std::make_pair(cands_map[top_c.id].get<C_VE>(), cands_map[top_c.id].get<C_VS>())) }; // ve->vs
			float dist_c_ve[2] = { cands_map[top_c.id].get<VE_DIST>(), cands_map[top_c.id].get<VS_DIST>() }; // for vs->ve, VE_DIST; for ve->vs, VS_DIST
			for (int i = 0; i < 2; ++i)
			{
				std::vector<ref_loc_entry>::const_iterator iter_entry, iter_entry_end; // use "const", as ITER_VS_VE and ITER_VE_VS are both const_iterator
				if (iters[i] != LNT.cend()) // vs->ve or ve->vs exists
				{
					// iterate all related reference locations
					iter_entry_end = iters[i]->second.second.cend();
					for (iter_entry = iters[i]->second.second.cbegin(); iter_entry != iter_entry_end; ++iter_entry)
					{
						// the case that offset > 0
						if (dist_c_ve[i] < iter_entry->get<OFFSET>()) // the condition that offset takes effect
						{
#ifdef ALGO_TESTING
							//algo_out_file << iters[i]->first.first << "-" << iters[i]->first.second << ", d(ve) = " << dist_c_ve[i] << " < offset = " << iter_entry->get<OFFSET>() << '\n';
#endif

							// utilizing virtual candidate c', which is the mapping of c with respect to lc, to compute ERd
							float dist_v_c_ve = 2.0f * iter_entry->get<OFFSET>() - dist_c_ve[i]; // d(c', ve) = 2 * offset - d(c, ve)
							top_c.ERD += (iter_entry->get<DDL>() + dist_v_c_ve) // Dd- + d(c, vj)
								* reflocs_map[iter_entry->get<REFLOC>()].get<PROB>(); // Pr(r)

							// set the \Delta(f)_offset
							boost::unordered_map<int, float>::iterator iter_delta = delta_f_offset.find(reflocs_map[iter_entry->get<REFLOC>()].get<NN>());
							if (iter_delta != delta_f_offset.end()) // f (nn_r) has been affected for other reference location(s) in Rc+
								iter_delta->second += reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();
							else
								delta_f_offset[reflocs_map[iter_entry->get<REFLOC>()].get<NN>()] = reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();

							continue;
						}

						float Rd = iter_entry->get<DDL>() + dist_c_ve[i]; // Dd- + d(c, ve)
						if (Rd <= 0) // c locates outside local network (<), or has no benefit (==)
							continue;
						else // Rd > 0
						{
							if (Rd > iter_entry->get<DDU>()) // > Dd+, must consider reverse direction
								continue;
							else // c takes effect
							{
								top_c.ERD += Rd * reflocs_map[iter_entry->get<REFLOC>()].get<PROB>(); // Rd * Pr(r)
								
								// set the \Delta(f)_offset
								boost::unordered_map<int, float>::iterator iter_delta = delta_f_offset.find(reflocs_map[iter_entry->get<REFLOC>()].get<NN>());
								if (iter_delta != delta_f_offset.end()) // f (nn_r) has been affected for other reference location(s) in Rc+
									iter_delta->second += reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();
								else
									delta_f_offset[reflocs_map[iter_entry->get<REFLOC>()].get<NN>()] = reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();
							}
						}
					}
				}
			}
		}
		else // top_c overlaps a vertex
		{
			// assign each reference location in top_c's Rc+ set to the corresponding nearest facility, as well as the \Delta(f)_offset value
			for (std::set<int>::const_iterator iter_Rc = Rc_u_sets[top_c.id].cbegin(); iter_Rc != Rc_u_sets[top_c.id].cend(); ++iter_Rc)
			{
				int nn_r = reflocs_map[*iter_Rc].get<NN>(); // nn(F,r)
				float delta_offset = reflocs_map[*iter_Rc].get<DIFF>(); // (d2nn - dnn) * prob
				
				boost::unordered_map<int, float>::iterator iter_delta = delta_f_offset.find(nn_r);
				if (iter_delta != delta_f_offset.end()) // f (nn_r) has been affected for other reference location(s) in Rc+
					iter_delta->second += delta_offset;
				else
					delta_f_offset[nn_r] = delta_offset;
			}
		}

#ifdef ALGO_TESTING
		//algo_out_file << "EDc+: c" << top_c.id << ": " << top_c.ERD << '\n';
#endif

		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// deal with Rc- reflecting by the \Delta(f)_offset
		// iteratively deal with two possible edges vs->ve and ve->vs
		L2NT_hash_map::const_iterator iters[2] = { L2NT.find(std::make_pair(cands_map[top_c.id].get<C_VS>(), cands_map[top_c.id].get<C_VE>())), // vs->ve
			L2NT.find(std::make_pair(cands_map[top_c.id].get<C_VE>(), cands_map[top_c.id].get<C_VS>())) }; // ve->vs
		float dist_c_ve[2] = { cands_map[top_c.id].get<VE_DIST>(), cands_map[top_c.id].get<VS_DIST>() }; // for vs->ve, VE_DIST; for ve->vs, VS_DIST
		for (int i = 0; i < 2; ++i)
		{
			std::vector<ref_loc_2_entry>::const_iterator iter_entry, iter_entry_end;
			if (iters[i] != L2NT.cend()) // vs->ve or ve->vs exists
			{
				// iterate all related reference locations
				iter_entry_end = iters[i]->second.cend();
				for (iter_entry = iters[i]->second.cbegin(); iter_entry != iter_entry_end; ++iter_entry)
				{
					// the case that offset > 0
					if (dist_c_ve[i] < iter_entry->get<OFFSET_2>()) // the condition that offset takes effect
					{
						// utilizing virtual candidate c', which is the mapping of c with respect to lc, to compute Ed
						float dist_v_c_ve = 2.0f * iter_entry->get<OFFSET_2>() - dist_c_ve[i]; // d(c', ve) = 2 * offset - d(c, ve)
						float Dd_l_d_c_vj = (iter_entry->get<DDL_2>() + dist_v_c_ve) // Dd- + d(c, vj)
							* reflocs_map[iter_entry->get<REFLOC>()].get<PROB>(); // Pr(r)
						if (Dd_l_d_c_vj >= 0) // c locates on LN not L2N
							continue;

						// set the \Delta(f)_offset
						boost::unordered_map<int, float>::iterator iter_delta = delta_f_offset.find(reflocs_map[iter_entry->get<REFLOC>()].get<NN>());
						if (iter_delta != delta_f_offset.end()) // f (nn_r) has been affected for other reference location(s) in Rc+ or Rc-
							iter_delta->second += Dd_l_d_c_vj + reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>(); // (dnn - d(r, c)) - (d2nn - dnn) = d2nn - d(r, c) to \Delta(f)_offset
						else
							delta_f_offset[reflocs_map[iter_entry->get<REFLOC>()].get<NN>()] = Dd_l_d_c_vj + reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();

						continue;
					}

					float Rd = iter_entry->get<DDL_2>() + dist_c_ve[i]; // Dd- + d(c, ve)
					if (Rd >= 0) // c locates on LN not L2N
						continue;
					else // Rd < 0
					{
						float ERd = Rd * reflocs_map[iter_entry->get<REFLOC>()].get<PROB>(); // Rd * Pr(r)

						// set the \Delta(f)_offset
						boost::unordered_map<int, float>::iterator iter_delta = delta_f_offset.find(reflocs_map[iter_entry->get<REFLOC>()].get<NN>());
						if (iter_delta != delta_f_offset.end()) // f (nn_r) has been affected for other reference location(s) in Rc+
						{
							if (-ERd <= reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>()) // Dd- + d(c, vj) < dnn - d2nn, and c is outside L2N
								iter_delta->second += ERd + reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();
						}
						else
						{
							if (-ERd <= reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>()) // Dd- + d(c, vj) < dnn - d2nn, and c is outside L2N
								delta_f_offset[reflocs_map[iter_entry->get<REFLOC>()].get<NN>()] = ERd + reflocs_map[iter_entry->get<REFLOC>()].get<DIFF>();
						}
					}
				}
			}
		}
		
		//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// evaluate the max negative effect for affected facilities
		std::pair<int, float> max_neg(-1, -FLT_MAX); // <f=null, -infinity>
		for (boost::unordered_map<int, float>::const_iterator iter_f_offset = delta_f_offset.cbegin(); iter_f_offset != delta_f_offset.cend(); ++iter_f_offset)
		{
			++checked_num; // check one more <f,r> pairs
			float delta_f = delta_f_map[iter_f_offset->first].first; // \Delta(f)
			delta_f += iter_f_offset->second; // after offset
			if (delta_f > max_neg.second) // larger \Delta(f)
			{
				max_neg.first = iter_f_offset->first;
				max_neg.second = delta_f;
			}
		}
		if (max_delta_f.second > max_neg.second) // compare with max \Delta(f)
		{
			++checked_num; // check one more <f,r> pairs
			max_neg.first = max_delta_f.first;
			max_neg.second = max_delta_f.second;
		}

#ifdef ALGO_TESTING
		algo_out_file << "LNB c" << top_c.id << ", ED = " << top_c.ERD << '\n';
#endif
		if (top_c.ERD + max_neg.second > op.ED)
		{
			op.f_id = max_neg.first;
			op.c_id = top_c.id;
			op.ED = top_c.ERD + max_neg.second;
		}
		
		max_heap.pop(); // pop the top vertex
	}

	return op; // the optimal facility-replacement pair
}