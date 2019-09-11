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

	// only check f_id and c_id, while ED value is ignored
	bool operator==(const facility_replacement &rhs) const
	{
		return f_id == rhs.f_id && c_id == rhs.c_id;
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

static std::ofstream algo_out_file("D:\\Experiment\\MFRM\\datasets\\algo_testing.txt", std::ofstream::out | std::ofstream::trunc);


#pragma region ALGO_EN
//================================================================================
// Candidate Accumulation Table (CAT) & Facility Accumulation Table (FAT) for EN algorithm

typedef boost::unordered_map < int, // candidate id
	float > // ED+
cand_ED_hash_map; // CAT, <candidate id, RD+> pair

typedef boost::unordered_map < std::pair <int, // facility id
	int >, // candidate id
	float > // offset value
fac_cand_offest_hash_map; // FAT, < <facility id, candidate id>, offset> pair

//--------------------------------------------------------------------------------
// top event visitor for EN algorithm
// remark: if a reference location has no nearest or sub-nearest facility (namely isolate sub-grahp), we ignore it and regard its contribute on ED as 0;
//		   then CAT & FAT will never record any ED+ or offset information of the reference location

class EN_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'c')
		{
			if (nn_id == -1) // nearest facility has not been found yet
				vn_dists[v.id] = v.dist; // record the shortest path distance from reference location to the candidate closer than nn
			else // nearest facility has been found
				v2_dists[v.id] = v.dist; // record the shortest path distance from reference location to the candidate farther than nn
#ifdef ALGO_TESTING
			//algo_out_file << "top: c" << v.id << ", d(r,c) = " << v.dist << '\n';
#endif
		}
		else if (v.type == 'f') // nearest/sub-nearest facility is encountered
		{
			if (nn_id == -1) // v is the nearest facility
			{
				nn_id = v.id;
				dnn = v.dist;
#ifdef ALGO_TESTING
				//algo_out_file << "dnn = " << v.dist << '\n';
#endif
			}
			else // v is the sub-nearest facility
			{
				float d2nn = v.dist; // define this variable for clarity
#ifdef ALGO_TESTING
				//algo_out_file << "d2nn = " << v.dist << '\n';
#endif

				// Evaluate \Delta(f), where f is nn_id and <nn_id, -1> do is \Delta(f)
				fac_cand_offest_hash_map::iterator iter_FAT = ptr_FAT_hash_map->find(std::make_pair(nn_id, -1));
				if (iter_FAT != ptr_FAT_hash_map->end()) // <facility, null> is exist in FAT
					iter_FAT->second += (dnn - d2nn) * prob; // aggregate (dnn - d2nn) * Pr(r) to <facility, null>
				else // <facility, null> is not exist in FAT yet
					ptr_FAT_hash_map->insert(std::make_pair(std::make_pair(nn_id, -1), (dnn - d2nn) * prob));
#ifdef ALGO_TESTING
				//algo_out_file << "<" << nn_id << ",-1> offset: " << (*ptr_FAT_hash_map)[std::make_pair(nn_id, -1)] << '\n';
#endif

				// iterate each encountered candidate (in Vn) before the nearest facility
				for (boost::unordered_map<int, float>::const_iterator iter_vn = vn_dists.cbegin(); iter_vn != vn_dists.cend(); ++iter_vn)
				{
					// aggregate (dnn - d(r,vk)) * Pr(r) to candidate vk ED+ into CAT
					cand_ED_hash_map::iterator iter_CAT = ptr_CAT_hash_map->find(iter_vn->first);
					if (iter_CAT != ptr_CAT_hash_map->end()) // candidate is exist in CAT
						iter_CAT->second += (dnn - iter_vn->second) * prob;
					else // candidate is not exist in CAT yet
						ptr_CAT_hash_map->insert(std::make_pair(iter_vn->first, (dnn - iter_vn->second) * prob));

					// aggregate (d2nn - dnn) * Pr(r) to <facility, vk> offset into FAT
					iter_FAT = ptr_FAT_hash_map->find(std::make_pair(nn_id, iter_vn->first));
					if (iter_FAT != ptr_FAT_hash_map->end()) // <facility, vk> is exist in FAT
						iter_FAT->second += (d2nn - dnn) * prob; 
					else // <facility, vk> is not exist in FAT yet
						ptr_FAT_hash_map->insert(std::make_pair(std::make_pair(nn_id, iter_vn->first), (d2nn - dnn) * prob));
#ifdef ALGO_TESTING
					//algo_out_file << "c" << iter_vn->first << " ED+: " << (*ptr_CAT_hash_map)[iter_vn->first] << '\n';
					//algo_out_file << "<" << nn_id << "," << iter_vn->first << "> offset: " << (*ptr_FAT_hash_map)[std::make_pair(nn_id, iter_vn->first)] << '\n';
#endif
				}

				// iterate each encountered candidate (in V2) after the nearest facility
				for (boost::unordered_map<int, float>::const_iterator iter_v2 = v2_dists.cbegin(); iter_v2 != v2_dists.cend(); ++iter_v2)
				{
					// aggregate (d2nn - d(r,vk)) * Pr(r) to <facility, vk> offset into FAT
					iter_FAT = ptr_FAT_hash_map->find(std::make_pair(nn_id, iter_v2->first));
					if (iter_FAT != ptr_FAT_hash_map->end()) // <facility, vk> is exist in FAT
						iter_FAT->second += (d2nn - iter_v2->second) * prob;
					else // <facility, vk> is not exist in FAT yet
						ptr_FAT_hash_map->insert(std::make_pair(std::make_pair(nn_id, iter_v2->first), (d2nn - iter_v2->second) * prob));
#ifdef ALGO_TESTING
					//algo_out_file << "<" << nn_id << "," << iter_v2->first << "> offset: " << (*ptr_FAT_hash_map)[std::make_pair(nn_id, iter_v2->first)] << '\n';
#endif
				}

				return true; // a flag to terminate graph traversal
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(cand_ED_hash_map *CAT, fac_cand_offest_hash_map *FAT) {
		ptr_CAT_hash_map = CAT; ptr_FAT_hash_map = FAT; }
	void set_prob(float prob_in) { prob = prob_in; }
	void reset() {
		vn_dists.clear(); // Vn for distances to candidates before nearest facility is found
		v2_dists.clear(); // V2 for distances to candidates after nearest facility is found
		nn_id = -1; // id of the nearest facility
		dnn = 0.0f; } // distance to the nearest facility

private:
	cand_ED_hash_map *ptr_CAT_hash_map; // a pointer to a hash map for storing <candidate id, ED+> pairs of CAT
	fac_cand_offest_hash_map *ptr_FAT_hash_map; // a pointer to a hash map for storing < <facility id, candidate id>, value> pairs of FAT
	float prob; // present probability of the reference location
	int nn_id; // id of the nearest facility of a reference location
	float dnn; // the distance from a reference location to the nearest facility
	boost::unordered_map < int, // candidate id
		float > // shortest path distance from reference location to the candidate
		vn_dists, v2_dists; // hash map for temporarily storing shortest distances from a reference location to encountered candidates before/after nearest facility is found

};

//--------------------------------------------------------------------------------
// EN algorithm
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [re] facility_replacement: the optimal facility-replacement pair; default facility_replacement object (f_id and c_id are both -1) means
//							  that candidates are not deployed successfully or opening reference locations file fails
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [in] cands_file_path: file path of candidates
// [out] traversal_time: the time (ms) of network traversals for all reference locations
// [in] k: return top-k facility-replacement pairs, defualt is 1; normally k >= 1, which means a definite top-k; 0 means output all FR pairs
// [out] topk: the vector for recording top-k facility-replacement pairs, default is NULL

facility_replacement EN(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path, const char *cands_file_path,
	__int64 &traversal_time, unsigned k = 1, std::vector<facility_replacement> *topk = NULL)
{
	// init graph
	typed_directed_graph graph;
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);

	// deploy facilities and candidates
	std::vector<int> fac_ids, cand_ids; // record ids of facilities and candidates, respectively
	graph.deploy_facilities(facs_file_path, fac_ids);
	graph.deploy_candidates(cands_file_path, cand_ids);

	if (!graph.get_is_candidates_deployed()) // check graph, facilities and candidates status
		return facility_replacement();

	// create an top event visitor for EN algorithm
	EN_visitor visitor;
	cand_ED_hash_map CAT; // <candidate, ED+> pairs in CAT
	fac_cand_offest_hash_map FAT; // < <facility id, candidate id>, value> pairs in FAT
	visitor.set_data_structures(&CAT, &FAT); // set data structures

	// deal with each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		auto begin_time = std::chrono::high_resolution_clock::now(); // starting time of network traversals

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
			// lon and lat are useless for graph

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it
#ifdef ALGO_TESTING
			//algo_out_file << "Ref Loc: " << refloc_id << ", prob = " << prob << '\n';
#endif

			// traverse graph from the reference location by dijkstra algorithm
			visitor.set_prob(prob); // set the present probability of the reference location
			visitor.reset(); // reset the hash maps that temporarily store shortest distances from a reference location to encountered candidates
			// before/after nearest facility is found; and reset nearest facility and the corresponding distance
			graph.dijkstra(refloc_v, &visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

		auto end_time = std::chrono::high_resolution_clock::now(); // ending time of network traversals
		traversal_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count();

		// rank facility-replacement pairs
		fac_rep_max_fibonacci_heap max_heap;
		for (std::vector<int>::const_iterator iter_cand = cand_ids.cbegin(); iter_cand != cand_ids.cend(); ++iter_cand)
		{
			// Rc+ aggregation
			float ED_c = 0.0f;
			cand_ED_hash_map::const_iterator iter_CAT = CAT.find(*iter_cand);
			if (iter_CAT != CAT.cend()) // candidate is exist in CAT ("not exist" means candidate affects no reference location)
				ED_c += iter_CAT->second;
#ifdef ALGO_TESTING
			algo_out_file << "EN c" << *iter_cand << ", ED = " << ED_c << '\n';
#endif

			// iterate all facilities
			for (std::vector<int>::const_iterator iter_fac = fac_ids.cbegin(); iter_fac != fac_ids.cend(); ++iter_fac)
			{
				float ED = ED_c;

				// Rf aggregation
				fac_cand_offest_hash_map::const_iterator iter_FAT = FAT.find(std::make_pair(*iter_fac, -1));
				if (iter_FAT != FAT.cend()) // <facility, null> is exist in FAT ("not exist" means facility affects no reference location)
				{
					ED += iter_FAT->second;
#ifdef ALGO_TESTING
					if (*iter_cand == 0)
						algo_out_file << "EN f" << iter_FAT->first.first << ", ED = " << iter_FAT->second << '\n';
#endif
				}

				// Rc- aggregation
				iter_FAT = FAT.find(std::make_pair(*iter_fac, *iter_cand));
				if (iter_FAT != FAT.cend()) // <facility, candidate> is exist in FAT ("not exist" means candidate doesn't affect facility)
					ED += iter_FAT->second;

				max_heap.push(facility_replacement(*iter_fac, *iter_cand, ED));
#ifdef ALGO_TESTING
				//algo_out_file << "<" << *iter_fac << ", " << *iter_cand << ">, ED = " << ED << '\n';
#endif
			}
		}

		facility_replacement op = max_heap.top();
#ifdef ALGO_TESTING
		//algo_out_file << "optimal: <" << op.f_id << ", " << op.c_id << ">, ED = " << op.ED << '\n';
		//algo_out_file << std::endl;
#endif

		// return k-top facility-replacement pairs
		if (topk != NULL)
		{
			if (k > 1)
			{
				unsigned n = 0;
				while (n < k)
				{
					topk->push_back(max_heap.top());
					++n;
					max_heap.pop(); // pop the top vertex
				}
			}
			else if (k == 0)
			{
				while (!max_heap.empty())
				{
					topk->push_back(max_heap.top());
					max_heap.pop(); // pop the top vertex
				}
			}
		}

		return op; // the optimal facility-replacement pair
	}
	return facility_replacement(); // opening reference locations file fails
}
#pragma endregion ALGO_EN

#pragma region ALGO_LNB
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
	float, // Dd- (Sigma lower bound), DDL_2
	float > // offset, OFFSET_2
ref_loc_2_entry;

enum ref_loc_2_entry_element { DDL_2 = 1, OFFSET_2 };

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
	float > // Delta(f)
delta_f_hash_map; // <f, \Delta(f)>

typedef boost::unordered_map < int, // refloc id
	boost::tuple < float, // prob, PROB
	int, // nn facility id, NN
	float > > // (d2nn - dnn) * prob, DIFF
refloc_hash_map; // not discussed in the paper

enum refloc_element { PROB = 0, NN, DIFF };

//--------------------------------------------------------------------------------
// top event visitor for LNB algorithm
// remark: if a candidate locates in a local network of a reference location, which has no nearest facility (namely isolate sub-grahp),
//			LNT and vertex_nnd hash maps will never record any local network of the reference location, and we view its ERD as 0

class LNB_top_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map)
	{
		if (v.type == 'f')
		{
			if (nn_id == -1) // nearest facility is found
			{
				nn_id = v.id; // nearest facility
				dnn = v.dist; // distance from a reference location to its nearest facility
			}
			else // sub-nearest facility is found
			{
				if (r_v_id != -1) // the reference location overlaps a vertex (vertex id is r_v_id)
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

				// consturct reflocs hash map, ie., <refloc_id, prob, nn, (d2nn - dnn) * Pr(r)>
				float diff = (dnn - v.dist) * prob; // (dnn - d2nn) * Pr(r)
				(*ptr_refloc_hash_map)[refloc_id] = boost::make_tuple(prob, nn_id, -diff); // (d2nn - dnn) * Pr(r)

				// construct <f, rnn(F,r), \Delta(f)> tuples
				delta_f_hash_map::iterator iter_delta = ptr_delta_f_hash_map->find(nn_id);
				if (iter_delta != ptr_delta_f_hash_map->end()) // the facility has been in the tuples
					iter_delta->second += diff; // (dnn - d2nn) * prob, aggregate \Delta(f)
				else // the first time presents the facility
					(*ptr_delta_f_hash_map)[nn_id] = diff;

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
					// 注意：当Dd- + d(c, vj) >= 0, 此时<vi, vj>(ie, <vs, ve>)被nn分割为<vi,nn>和<nn,ve>两段，c落在<vi,nn>上，即落在LN上，而非L2N上，故在计算L2NT时不计算此值
					//		 当Dd- + d(c, vj) < dnn - d2nn, 此时c落在2nn之外，即L2N之外，故c不影响f

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
					if (is_in_LNT)
					{
						float Dd_u = dnn - iter_edge->second.first; // Dd+ = dnn - dist(r, vs)
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
					}

					// create new entry in L2NT
					if (is_in_L2NT)
					{
						L2NT_hash_map::iterator iter_L2NT = ptr_L2NT_hash_map->find(vs_ve);
						if (iter_L2NT != ptr_L2NT_hash_map->end()) // edge <vs_id, ve_id> has already been in local sub-network of other reference locations
						{
							iter_L2NT->second.push_back(boost::make_tuple(refloc_id,
								dnn - iter_edge->second.second, // Dd- = dnn - (dist(r, vs) + dist(vs, ve))
								offset)); // offset
						}
						else // <vs_id, ve_id> first time presents
						{
							(*ptr_L2NT_hash_map)[vs_ve] = std::vector<ref_loc_2_entry>();
							(*ptr_L2NT_hash_map)[vs_ve].push_back(boost::make_tuple(refloc_id,
								dnn - iter_edge->second.second, // Dd- = dnn - (dist(r, vs) + dist(vs, ve))
								offset)); // offset
						}
					}
				}
				return true; // a flag to terminate graph traversal
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(LNT_hash_map *LNT_hash_map_in, L2NT_hash_map *L2NT_hash_map_in, delta_f_hash_map *delta_f_hash_map_in,
		refloc_hash_map *refloc_hash_map_in, loc_edges_hash_map *loc_edges_hash_map_in, vertex_dnn_hash_map *dnn_hash_map)
	{
		ptr_LNT_hash_map = LNT_hash_map_in;
		ptr_L2NT_hash_map = L2NT_hash_map_in;
		ptr_delta_f_hash_map = delta_f_hash_map_in;
		ptr_refloc_hash_map = refloc_hash_map_in;
		ptr_loc_edges_hash_map = loc_edges_hash_map_in;
		ptr_dnn_hash_map = dnn_hash_map;
	}
	void set_overlap(int r_v_id_in) { r_v_id = r_v_id_in; };
	void set_refloc(int refloc_id_in, float prob_in) {
		refloc_id = refloc_id_in; prob = prob_in; // set the reference location and its present probability
		nn_id = -1; dnn = 0.0f;	} // reset nearest facility

private:
	LNT_hash_map *ptr_LNT_hash_map; // a pointer to a hash map as Local Network Table
	L2NT_hash_map *ptr_L2NT_hash_map; // a pointer to a hash map as Local Sub-Network Table
	delta_f_hash_map *ptr_delta_f_hash_map; // a pointer to a hash map for <f, \Delta(f)> pairs
	refloc_hash_map *ptr_refloc_hash_map; // a pointer to a hash map of reference locations
	loc_edges_hash_map *ptr_loc_edges_hash_map; // a pointer to an assistant hash map for local network edges
	vertex_dnn_hash_map *ptr_dnn_hash_map; // a pointer to an assistant hash map for reference locations that overlap some vertices
	int r_v_id; // vertex id, which indicates whether this reference location overlaps a vertex (v_id) or not (-1)
	int refloc_id; // id of the reference location which is the source vertex of the current traversal
	float prob; // present probability of the reference location
	int nn_id; // id of the nearest facility of a reference location
	float dnn; // the distance from a reference location to its nearest facility
};

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

//--------------------------------------------------------------------------------
// LNB algorithm: construct LNT and L2NT
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [out] LNT; a hash map as Local Network Table, which needs to be constructed in this function
// [out] L2NT: a hash map as Local Sub-Network Table, which needs to be constructed in this function
// [out] dnn_map; a hash map for reference locations that overlap some vertices, which needs to be constructed in this function
// [out] delta_f_map: a hash map for <f, \Delta(f)> pairs, which needs to be constructed in this function
// [out] max_delta_f: the max \Delta(f) for all facilities, which needs to be set in this function
// [out] graph: the typed directed graph, which needs to be constructed in this function
// [out] facs_map: a hash map for storing facilities associated with the endpoints, which needs to be constructed in this function
// [out] reflocs_map: a hash map for storing reference locations associated with the probabilities and (dnn - d2nn) * Pr(r), which needs to be constructed in this function

void LNB_construct_LNT(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path,
	LNT_hash_map &LNT, L2NT_hash_map &L2NT, vertex_dnn_hash_map &dnn_map, delta_f_hash_map &delta_f_map, std::pair<int, float> &max_delta_f,
	typed_directed_graph &graph, fac_hash_map &facs_map, refloc_hash_map &reflocs_map)
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
	top_visitor.set_data_structures(&LNT, &L2NT, &delta_f_map, &reflocs_map, &refloc_loc_edges_hash_map, &dnn_map); // set data structures

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
		if (iter_delta->second > max_delta_f.second) // \Delta(f) > max_delta_f
		{
			max_delta_f.first = iter_delta->first; // f
			max_delta_f.second = iter_delta->second; // \Delta(f)
		}
	}
}

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
// [in] delta_f_map: a hash map for <f, \Delta(f)> pairs
// [in] max_delta_f: the max \Delta(f) for all facilities
// [in] graph: the constructed graph
// [in] facs_map: a hash map for storing facilities associated with the endpoints
// [in] reflocs_map: a hash map for storing reference locations associated with the probabilities and (dnn - d2nn) * Pr(r)
// [in] cands_file_path: file path of candidates

facility_replacement LNB_query(int &checked_num, LNT_hash_map &LNT, L2NT_hash_map &L2NT, const vertex_dnn_hash_map &dnn_map, delta_f_hash_map &delta_f_map,
	std::pair<int, float> &max_delta_f, const typed_directed_graph &graph, const fac_hash_map &facs_map, refloc_hash_map &reflocs_map, const char *cands_file_path)
{
#ifdef ALGO_TESTING
	algo_out_file << "LNB max_delta_f: <" << max_delta_f.first << ", " << max_delta_f.second << ">\n";
#endif

	cand_max_fibonacci_heap max_heap; // max-heap ordered by ED+ in LNT (LNH in paper)
	cand_info_hash_map hash_map; // hash map for information of candidates

	boost::unordered_map < int, // candidate id
		std::set < int > > // refloc ids
		Rc_u_sets; // Rc+ set for candidates, each of which overlaps with some vertex
				   // NOTE: this set variable is only for each candidate that overlaps a vertex, whose Rc+ set is intialized when constructing LNH

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// construct a max-heap (LNH) for candidates
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
				hash_map[cand_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve); // record information of each candidate
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
				hash_map[cand_id] = boost::make_tuple(vs_id, ve_id, 0.0f, 0.0f); // set overlap flags (i.e., LNT.cend() and 0.0f) for the candidate
			}
		}
	}
	else // opening candidates file fails
		return facility_replacement();

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
		if (hash_map[top_c.id].get<VE_DIST>() != 0 && hash_map[top_c.id].get<VS_DIST>() != 0) // c splits an edge, but not overlaps a vertex, as overlapping has definite ED+
		{
			top_c.ERD = 0.0f; // prepare to accumulate actual ED+ of this candidate
							  // and Rc+ set will be construct in this "if" case code block

			// iteratively deal with two possible edges vs->ve and ve->vs
			LNT_hash_map::const_iterator iters[2] = { LNT.find(std::make_pair(hash_map[top_c.id].get<C_VS>(), hash_map[top_c.id].get<C_VE>())), // vs->ve
				LNT.find(std::make_pair(hash_map[top_c.id].get<C_VE>(), hash_map[top_c.id].get<C_VS>())) }; // ve->vs
			float dist_c_ve[2] = { hash_map[top_c.id].get<VE_DIST>(), hash_map[top_c.id].get<VS_DIST>() }; // for vs->ve, VE_DIST; for ve->vs, VS_DIST
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
		L2NT_hash_map::const_iterator iters[2] = { L2NT.find(std::make_pair(hash_map[top_c.id].get<C_VS>(), hash_map[top_c.id].get<C_VE>())), // vs->ve
			L2NT.find(std::make_pair(hash_map[top_c.id].get<C_VE>(), hash_map[top_c.id].get<C_VS>())) }; // ve->vs
		float dist_c_ve[2] = { hash_map[top_c.id].get<VE_DIST>(), hash_map[top_c.id].get<VS_DIST>() }; // for vs->ve, VE_DIST; for ve->vs, VS_DIST
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
			float delta_f = delta_f_map[iter_f_offset->first]; // \Delta(f)
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

#pragma endregion ALGO_LNB

#pragma region ALGO_NSJ
//================================================================================
// NNFC hash map for NSJ algorithm;
// an assistant hash map for sets of reference locations which cover some specific candidate;
// an assistant hash map for Euclidean distance dE(r,c)

typedef boost::unordered_map < int, // reference location id
	boost::tuple < float, // prob, use PROB in "refloc_element"
	int, // nn facility id, use NN in "refloc_element"
	float, // (d2nn - dnn) * Pr(r), use DIFF in "refloc_element"
	geo_point, // coordinate, COORDINATE
	float, // dnn, DNN
	float, // d2nn, D2NN
	geo_box, // NNFC MBR, NNFC_MBR
	geo_box, // 2NNFC MBR, NNFC_2_MBR
	int > > // number of covered candidates (HM.r), COVER_NUM
nnfc_hash_map; // NNFC hash map

enum nnfc_hash_map_element { //PROB = 0, use this enumeration value in "refloc_element"
	// NN = 1, use NN in "refloc_element"
	//DIFF = 2, use this enumeration value in "refloc_element"
	COORDINATE = 3, DNN, D2NN, NNFC_MBR, NNFC_2_MBR, COVER_NUM };

enum direction { NORTHWARD = 0, SOUTHWARD, EASTWARD, WESTWARD };

typedef boost::unordered_map < int, // candidate id
	std::pair < float, // \Delta+(c)
	std::set <int> > > // ids of reference locations whose 2-NNFCs cover the candidate, namely Rc'
cover_hash_map;

typedef boost::unordered_map < int, // candidate id
	std::map < int, // facility id
	float > > // offset for \Delta(f,c)
offset_c_f_hash_map;

typedef boost::unordered_map < int, // candidate id
	float > // \Delta(c)
delta_c_hash_map;

//--------------------------------------------------------------------------------
// top event visitor for calculate NNFC
// remark: if a reference location has no nearest or sub-nearest facility (namely isolate sub-grahp), we ignore it and regard its contribute on ED as 0;
//		   then NNFC hash map will never record the reference location

class NNFC_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'f')
		{
			if (nn == -1) // nearest facility is found
			{
				nn = v.id; // nearest neighbor facility
				dnn = v.dist; // nearset neighbor distance, namely the NNFC radius
#ifdef ALGO_TESTING
				//algo_out_file << "dnn(r" << refloc_id << ") = " << dnn << '\n';
#endif
			}
			else // nn != -1, sub-nearest facility is found
			{
#ifdef ALGO_TESTING
				//algo_out_file << "d2nn(r" << refloc_id << ") = " << v.dist << '\n';
#endif

				// construct <f, \Delta(f), rnn(F,f)>
				delta_f_hash_map::iterator iter_find = ptr_f_delta_f_hash_map->find(nn); // iter_find (must) != end(), as initialized in "NSJ_construct_NNFC"
				iter_find->second += (dnn - v.dist) * prob; // \Delta(f) += (dnn - d2nn) * Pr(r), where r \in rnn(F,f)
#ifdef ALGO_TESTING
				//algo_out_file << "delta(f" << v.id << ") = " << iter_find->second.first << '\n';
#endif

				// calculate offsets coordinate for NNFC & 2-NNFC MBR
				float min_lon, min_lat, max_lon, max_lat, // for NNFC
					min_lon_2, min_lat_2, max_lon_2, max_lat_2; // for 2-NNFC
				geo_offset(lon, lat, dnn, WESTWARD, min_lon);
				geo_offset(lon, lat, dnn, SOUTHWARD, min_lat);
				geo_offset(lon, lat, dnn, EASTWARD, max_lon);
				geo_offset(lon, lat, dnn, NORTHWARD, max_lat);
				geo_offset(lon, lat, v.dist, WESTWARD, min_lon_2);
				geo_offset(lon, lat, v.dist, SOUTHWARD, min_lat_2);
				geo_offset(lon, lat, v.dist, EASTWARD, max_lon_2);
				geo_offset(lon, lat, v.dist, NORTHWARD, max_lat_2);

				(*ptr_NNFC_hash_map)[refloc_id] = // reference location id
					boost::make_tuple(prob, // present probability
					nn, // nn facility id
					(v.dist - dnn) * prob, // DIFF = (d2nn - dnn) * Pr(r)
					geo_point(lon, lat), // coordinate
					dnn, // dnn
					v.dist, // d2nn
					geo_box(geo_point(min_lon, min_lat), geo_point(max_lon, max_lat)), // NNFC MBR
					geo_box(geo_point(min_lon_2, min_lat_2), geo_point(max_lon_2, max_lat_2)), // 2-NNFC MBR
					0); // candidates being covered by 2-NNFC (HM.r)

				return true; // a flag to terminate graph traversal
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structure(nnfc_hash_map *ptr_NNFC_hash_map_in, delta_f_hash_map *ptr_f_delta_f_hash_map_in) {
		ptr_NNFC_hash_map = ptr_NNFC_hash_map_in; ptr_f_delta_f_hash_map = ptr_f_delta_f_hash_map_in; }
	void set_refloc(int refloc_id_in, float lon_in, float lat_in, float prob_in) { refloc_id = refloc_id_in; lon = lon_in; lat = lat_in; prob = prob_in;
		nn = -1; dnn = 0.0f; } // also initialize dnn

public:
	// calculate offset coordinate after moving to a direction
	// [in] lon_c: longitude of source coordinate
	// [in] lat_c: latitude of source coordinate
	// [in] distance: distance to be moved
	// [in] dir: the direction
	// [out] lon_or_lat: resulting longitude or latitude
	static void geo_offset(float lon_c, float lat_c, float distance, direction dir, float &lon_or_lat)
	{
		if (dir == NORTHWARD || dir == SOUTHWARD) // the same longitude
		{
			float degree = 360.0f * distance / (2.0f * PI * EARTH_RADIUS);

			// lon_or_lat now is lat
			if (dir == NORTHWARD)
				lon_or_lat = lat_c + degree;
			else // dir == SOUTHWARD
				lon_or_lat = lat_c - degree;
		}
		else if (dir == EASTWARD || dir == WESTWARD) // the same latitude
		{
			float radian_latitude = lat_c * PI / 180.0f;
			float degree = 360.0f * distance / (2.0f * PI * EARTH_RADIUS * cos(radian_latitude));

			// lon_or_lat now is lon
			if (dir == EASTWARD)
				lon_or_lat = lon_c + degree;
			else // dir == WESTWARD
				lon_or_lat = lon_c - degree;
		}
	}

protected:
	int nn; // nn of r
	float dnn; // dnn of r
	nnfc_hash_map *ptr_NNFC_hash_map; // a pointer to a hash map for NNFCs
	delta_f_hash_map *ptr_f_delta_f_hash_map; // a pointer to a hash map for <f, \Delta(f), rnn(F,f)>
	int refloc_id; // reference location id
	float lon; // longitude of the reference location
	float lat; // latitude of the reference location
	float prob; // present probability of the reference location
};

//--------------------------------------------------------------------------------
// top event visitor for NSJ algorithm

class NSJ_visitor
	: public top_event_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'c') // c must be inside L2N of reference location
		{
			if (!nn_met) // nn hasn't been encountered, in LN
			{
				// accumulate \Delta(c)
				delta_c_hash_map::iterator iter_delta_c = ptr_delta_c->find(v.id);
				if (iter_delta_c != ptr_delta_c->end()) // candidate id exists
					iter_delta_c->second += ((*ptr_NNFCs)[refloc_id].get<DNN>() - v.dist) * (*ptr_NNFCs)[refloc_id].get<PROB>(); // \Delta(c) += (dnn(r) - d(r,c)) * Pr(r)
				else
					(*ptr_delta_c)[v.id] = ((*ptr_NNFCs)[refloc_id].get<DNN>() - v.dist) * (*ptr_NNFCs)[refloc_id].get<PROB>(); // init \Delta(c) = (dnn(r) - d(r,c)) * Pr(r)
#ifdef ALGO_TESTING
				algo_out_file << "r" << refloc_id << ", Delta(c" << v.id << ") = " << (*ptr_delta_c)[v.id] << '\n';
#endif

				// accumulate offset of \Delta(<nn(r),c>)
				offset_c_f_hash_map::iterator iter_c = ptr_offset_c_f->find(v.id);
				if (iter_c != ptr_offset_c_f->end()) // candidate id exists
				{
					int f_id = (*ptr_NNFCs)[refloc_id].get<NN>();
					std::map<int, float>::iterator iter_f = iter_c->second.find(f_id);
					if (iter_f != iter_c->second.end()) // facility id exists
						iter_f->second += (*ptr_NNFCs)[refloc_id].get<DIFF>(); // += (d2nn(r) - dnn) * Pr(r)
					else
						(iter_c->second)[f_id] = (*ptr_NNFCs)[refloc_id].get<DIFF>(); // init (d2nn(r) - dnn) * Pr(r)
				}
				else
				{
					(*ptr_offset_c_f)[v.id] = std::map<int, float>();
					(*ptr_offset_c_f)[v.id][(*ptr_NNFCs)[refloc_id].get<NN>()] = (*ptr_NNFCs)[refloc_id].get<DIFF>(); // init (d2nn(r) - dnn) * Pr(r)
				}
#ifdef ALGO_TESTING
				algo_out_file << "r" << refloc_id << ", Delta(f" << (*ptr_NNFCs)[refloc_id].get<NN>() << ", c" << v.id << ")offset = d2nn " << (*ptr_NNFCs)[refloc_id].get<DIFF>() << '\n';
#endif
			}
			else // in L2N
			{
				// accumulate offset of \Delta(<nn(r),c>)
				offset_c_f_hash_map::iterator iter_c = ptr_offset_c_f->find(v.id);
				if (iter_c != ptr_offset_c_f->end()) // candidate id exists
				{
					int f_id = (*ptr_NNFCs)[refloc_id].get<NN>();
					std::map<int, float>::iterator iter_f = iter_c->second.find(f_id);
					if (iter_f != iter_c->second.end()) // facility id exists
						iter_f->second += ((*ptr_NNFCs)[refloc_id].get<D2NN>() - v.dist) * (*ptr_NNFCs)[refloc_id].get<PROB>(); // += (d2nn(r) - d(r,c)) * Pr(r)
					else
						(iter_c->second)[f_id] = ((*ptr_NNFCs)[refloc_id].get<D2NN>() - v.dist) * (*ptr_NNFCs)[refloc_id].get<PROB>(); // init (d2nn(r) - d(r,c)) * Pr(r)
				}
				else
				{
					(*ptr_offset_c_f)[v.id] = std::map<int, float>();
					(*ptr_offset_c_f)[v.id][(*ptr_NNFCs)[refloc_id].get<NN>()] = ((*ptr_NNFCs)[refloc_id].get<D2NN>() - v.dist)
						* (*ptr_NNFCs)[refloc_id].get<PROB>(); // init (d2nn(r) - d(r,c)) * Pr(r)
				}
#ifdef ALGO_TESTING
				algo_out_file << "r" << refloc_id << ", Delta(f" << (*ptr_NNFCs)[refloc_id].get<NN>() << ", c" << v.id << ")offset = d(r,c) "
					<< ((*ptr_NNFCs)[refloc_id].get<D2NN>() - v.dist) * (*ptr_NNFCs)[refloc_id].get<PROB>() << '\n';
#endif
			}

			if (--((*ptr_NNFCs)[refloc_id].get<COVER_NUM>()) // a candidate is encountered, then the number of remnant candidates is decreases by 1
				== 0) // all d(r,c) <= dnn/d2nn, which means all candidates have benefits/negative effects, hence, no need to traverse until 2nn facility
				return true; // a flag to terminate graph traversal
		}
		else if (v.type == 'f') // nearest facility is found
		{
			if (!nn_met) // v is nn
				nn_met = true; // nn has been met
			else // v is 2nn
			{
				(*ptr_NNFCs)[refloc_id].get<COVER_NUM>() = 0; // already traversed from this reference location, then reset this value as a flag
				return true; // a flag to terminate graph traversal
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structures(nnfc_hash_map *NNFCs, offset_c_f_hash_map *offset_c_f, delta_c_hash_map *delta_c) {
		ptr_NNFCs = NNFCs; ptr_offset_c_f = offset_c_f;	ptr_delta_c = delta_c; }
	void set_refloc(int refloc_id_in) {	refloc_id = refloc_id_in; nn_met = false; }

protected:
	bool nn_met; // indicate whether nn facility has been met
	int refloc_id; // reference location id
	std::set<int> *ptr_Rc; // a pointer to actual Rc set, NOTE must be reset for each <f,c> pair
	nnfc_hash_map *ptr_NNFCs; // a pointer to a hash map for NNFCs
	offset_c_f_hash_map *ptr_offset_c_f; // a pointer to a hash map for offset of \Delta(<f,c>)
	delta_c_hash_map *ptr_delta_c; // a pointer to a hash map for \Delta(c)
};

//--------------------------------------------------------------------------------
// NSJ algorithm: construct NNFCs R-tree
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [out] f_r_hash_map: a hash map for <f, \Delta(f), rnn(F,f)> tuples
// [out] max_delta_f: the facility with max \Delta(f)
// [out] reflocs: a hash map for the network location information of reference locations, which uses a "fac_hash_map" as alternative data structure
// [out] NNFCs; a hash map for NNFCs, which needs to be constructed in this function
// [out] graph: the typed directed graph, which needs to be constructed in this function
// [reserved/out] NNFC_rtree; this parameter is now reserved; an R-tree for NNFCs, which needs to be constructed in this function

void NSJ_construct_NNFC(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path,
	delta_f_hash_map &f_r_hash_map, std::pair<int, float> &max_delta_f, fac_hash_map &reflocs, nnfc_hash_map &NNFCs, typed_directed_graph &graph)
{
	// init graph
#ifdef ALGO_TESTING_TRACE_GRAPH
	graph.set_testing(&algo_out_file);
#endif
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);
	
	// deploy facilities and initialize <f, \Delta(f), rnn(F,f)> tuples
	std::vector<int> fac_ids; // record all ids of facilities
	graph.deploy_facilities(facs_file_path, fac_ids);
	for (std::vector<int>::iterator iter_f = fac_ids.begin(); iter_f != fac_ids.end(); ++iter_f)
		f_r_hash_map[*iter_f] = 0.0f; // must initialize the hash map, as some facility may affect no reference location

	if (!graph.get_is_facilities_deployed()) // check graph and facilities status
		return;

	// create a top event visitor for NNFCs
	NNFC_visitor top_visitor;
	top_visitor.set_data_structure(&NNFCs, &f_r_hash_map); // set data structure

	// compute NNFC/2NNFC of each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// record network location information for each reference location
			reflocs[refloc_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve);

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id, lon, lat, prob); // set id, geo-information and present probability of a reference location
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}
	}

	// evaluate the max \Delta(f)
	max_delta_f.first = -1;
	max_delta_f.second = -FLT_MAX;
	for (delta_f_hash_map::iterator iter_f = f_r_hash_map.begin(); iter_f != f_r_hash_map.end(); ++iter_f)
	{
		if (iter_f->second > max_delta_f.second) // \Delta(f) > max_delta_f
		{
			max_delta_f.first = iter_f->first; // f
			max_delta_f.second = iter_f->second; // \Delta(f)
		}
	}
}


//--------------------------------------------------------------------------------
// NSJ algorithm: query optimal candidate based on NNFCs
// remark: the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
// [re] facility_replacement: the id of optimal facility replacement pair
// [out] checked_num: the number of candidates actually checked
// [out] time_costs: the time costs (ms) of querying R*-tree and constructing FRH
// [in] f_r_hash_map: a hash map for <f, \Delta(f), rnn(F,f)> tuples
// [in] max_delta_f: the facility with max \Delta(f)
// [in] reflocs: a hash map for the network location information of reference locations, which uses a "fac_hash_map" as alternative data structure
// [in/out] NNFCs; a hash map for NNFCs, COVER_NUM element will be changed
// [in/out] graph: the constructed graph, candidates will be further deployed
// [in] cands_file_path: file path of candidates

facility_replacement NSJ_query(int &checked_num, __int64 *time_costs, delta_f_hash_map &f_r_hash_map, std::pair<int, float> &max_delta_f,
	fac_hash_map &reflocs, nnfc_hash_map &NNFCs, typed_directed_graph &graph, const char *cands_file_path)
{
#ifdef ALGO_TESTING
	algo_out_file << "NSJ max_delta_f: <" << max_delta_f.first << ", " << max_delta_f.second << ">\n";
#endif

	auto begin_time = std::chrono::high_resolution_clock::now();

	// deploy candidates onto the directed graph and construct candidates R*-tree
	geo_cand_rtree cand_rtree;
	boost::unordered_map<int, geo_point> cand_geos;
	graph.deploy_candidates_rtree(cands_file_path, cand_rtree, cand_geos);

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join (alternative by R*-tree range query): construct hash maps for <c, <\Delta+(c), Rc'> >
	cover_hash_map Rc_map; // sets <c, <\Delta+(c), Rc'> >
	for (nnfc_hash_map::iterator iter_NNFC = NNFCs.begin(); iter_NNFC != NNFCs.end(); ++iter_NNFC)
	{
		// range query for candidates which are covered by a 2-NNFC
		std::vector<geo_cand> covered_cands; // covered by 2-NNFC
		cand_rtree.query(boost::geometry::index::within(iter_NNFC->second.get<NNFC_2_MBR>()), std::back_inserter(covered_cands));
		iter_NNFC->second.get<COVER_NUM>() = static_cast<int>(covered_cands.size()); // record the maximum number of (potentially) covered (by 2-NNFC) candidates

		// record the set of ids of reference locations whose 2-NNFCs cover the candidate, namely Rc'
		for (std::vector<geo_cand>::iterator iter_cand = covered_cands.begin(); iter_cand != covered_cands.end(); ++iter_cand)
		{
			// Rc+': inside NNFC
			if (boost::geometry::within(iter_cand->first, iter_NNFC->second.get<NNFC_MBR>()))
			{
				// compute Ed+ for pair <r,c>
				float distE = boost::geometry::distance(iter_NNFC->second.get<COORDINATE>(), iter_cand->first) * EARTH_RADIUS; // dE(r,c)
				float Ed_u = (iter_NNFC->second.get<DNN>() - distE) * iter_NNFC->second.get<PROB>(); // Ed+r(c) = (dnn - dE(r,c)) * Pr(r)
				if (Ed_u > 0.0f) // possible benefit
					// >: candidate locates inside NNFC
					// <: candidate locates inside the 4 corners between NNFC and its MBR
					// ==: exactly on the NNFC bound, maybe < for road network distance
				{
					// \Delta+(c)'
					cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
					if (iter_Rc_find != Rc_map.end()) // candidate has already been created
					{
						iter_Rc_find->second.first += Ed_u; // \Delta+(c) += (dnn - dE(r,c)) * Pr(r)
						iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					else // need to create a set for the candidate
					{
						Rc_map[iter_cand->second] = std::make_pair(Ed_u, std::set<int>()); // Ed_u as the initial \Delta+(c)' value
						Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					continue; // next candidate
				}
			}
			
			// maybe inside 2-NNFC ring
			cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
			if (iter_Rc_find != Rc_map.end()) // candidate has already been created
				iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			else // need to create a set for the candidate
			{
				Rc_map[iter_cand->second] = std::make_pair(0.0f, std::set<int>()); // 0 as the initial \Delta+(c)' value, as at most inside 2-NNFC ring
				Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			}
		}
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	time_costs[0] = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count(); // Spatial Join time cost
	begin_time = std::chrono::high_resolution_clock::now();

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// construct Max-heap of candidates & iterate candidate Max-heap
	cand_max_fibonacci_heap cand_max_heap; // ordered by \Delta+(c)
	for (cover_hash_map::iterator iter_Rc = Rc_map.begin(); iter_Rc != Rc_map.end(); ++iter_Rc) // for each candidate inside 2-NNFC
		cand_max_heap.push(candidate(iter_Rc->first, iter_Rc->second.first)); // push each candidate into max-heap

	facility_replacement op(-1, -1, 0.0f); // initialize optimal facility-replacement pair
	checked_num = 0; // initialize

	// create a top event visitor for NSJ algorithm
	NSJ_visitor top_visitor;
	offset_c_f_hash_map offset_c_f; // offset for \Delta(f, c)
	delta_c_hash_map delta_c; // \Delta(c)
	top_visitor.set_data_structures(&NNFCs, &offset_c_f, &delta_c); // set data structures

	while (!cand_max_heap.empty())
	{
		// retrieve the top candidate in max-heap
		candidate top_c = cand_max_heap.top();
#ifdef ALGO_TESTING
		algo_out_file << "c" << top_c.id << ", ERD = " << top_c.ERD << '\n';
#endif

#ifndef ALGO_TESTING_NO_EARLY_STOPPING
		if (top_c.ERD < op.ED) // current op <f,c> is the optimal facility-replacement pair
			break;
#endif
		
		// traverse each reference location whose 2-NNFC covers top_c
		for (std::set<int>::iterator iter_refloc = Rc_map[top_c.id].second.begin(); iter_refloc != Rc_map[top_c.id].second.end(); ++iter_refloc) // iterate each r \in Rc'
		{
			int refloc_id = *iter_refloc;
			if (NNFCs[refloc_id].get<COVER_NUM>() == 0) // all candidates that can affect the reference location have been traversed by a certain candidate prior to this top_c
				continue; // traverse from next reference location

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(reflocs[refloc_id].get<VS_ID>(), 'v'), vertex(reflocs[refloc_id].get<VE_ID>(), 'v'), refloc_v,
				reflocs[refloc_id].get<VS_DIST>(), reflocs[refloc_id].get<VE_DIST>(), &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id); // now the ED of top <f,c> is initialized as \Delta(f)
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

		// now, \Delta(top_c) and offset of \Delta(<f,c>) have been obtained
		// calculate max_delta_f w.r.t. top_c
		std::pair<int, float> fac(max_delta_f);
		for (std::map<int, float>::iterator iter_f = offset_c_f[top_c.id].begin(); iter_f != offset_c_f[top_c.id].end(); ++iter_f)
		{
			++checked_num; // need to check more facility-replacement pairs
			float delta_f = iter_f->second + f_r_hash_map[iter_f->first]; // offset + \Delta(f)
#ifdef ALGO_TESTING
			algo_out_file << "<" << iter_f->first << ", " << top_c.id << "> offset: " << iter_f->second << ", Delta(f) = " << f_r_hash_map[iter_f->first].first
				<< ">, Delta(f) wrs c: " << delta_f << '\n';
#endif
			if (delta_f > fac.second)
			{
				fac.first = iter_f->first;
				fac.second = delta_f;
			}
		}
		
		// update optima
#ifdef ALGO_TESTING
		algo_out_file << "c" << top_c.id << ", ED = " << delta_c[top_c.id] << '\n';
#endif
		float delta = delta_c[top_c.id] + fac.second; // \Delta(c) + \Delta(f) w.r.t. c
		if (delta > op.ED)
		{
			op.f_id = fac.first;
			op.c_id = top_c.id;
			op.ED = delta;
		}
#ifdef ALGO_TESTING
		algo_out_file << "current optima: <" << op.f_id << "," << op.c_id << ">, ED = " << op.ED << '\n';
#endif

		cand_max_heap.pop(); // pop the top candidate
	}

	end_time = std::chrono::high_resolution_clock::now();
	time_costs[1] = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count(); // iteration of Max-heap time cost

	return op; // the optimal facility-replacement pair
}
#pragma endregion ALGO_NSJ

//================================================================================
// comment reason: the RIC may be not useful and necessary,
//				   because the interference relationship, which is evaluated by RIC,
//				   cannot avoid the network traversal for evaluating network distance \Delta(c) as that in Euclidean space,
//				   which is the dominant time cost in road network space
//#define RID	// the define is now in "typed_directed_graph.h", as it is also used in "typed_directed_graph.cpp"

#pragma region ALGO_RID
//================================================================================
// RIC hash map for RID algorithm;
// an assistant hash map for sets of f.rid

typedef boost::unordered_map < int, // facility id
	float > // f.rid = max{dnn + d2nn}
rid_hash_map;

typedef boost::unordered_map < int, // facility id
	geo_box > // RIC MBR
ric_hash_map; // RIC hash map

//--------------------------------------------------------------------------------
// a struct for \Delta(f)s or \Delta+(c)s fibonacci max-heap

struct delta_value
{
	int id;	// the id of a facility or a candidate
	float val; // the value of \Delta(f) or \Delta(c) for a facility or candidate respectively

	delta_value(int id_in, float val_in)
		: id(id_in)
		, val(val_in) {}

	delta_value(const delta_value &rhs)
	{
		if (this == &rhs)
			return;
		id = rhs.id;
		val = rhs.val;
	}

	// remark: the same vertex validation is not implemented
	bool operator<(const delta_value &rhs) const
	{
		return val < rhs.val; // only compare \Delta(f) or \Delta(c)
	};
};

typedef boost::heap::fibonacci_heap<delta_value> max_delta_fibonacci_heap;

//--------------------------------------------------------------------------------
// top event visitor for calculate NNFC & RIC
// remark: if a reference location has no nearest or sub-nearest facility (namely isolate sub-grahp), we ignore it and regard its contribute on ED as 0;
//		   then NNFC hash map will never record the reference location

class NNFC_RIC_visitor
	: public NNFC_visitor
{
public:
	virtual bool operator()(const target_vertex &v, const definite_hash_map &hash_map) // 2nd parameter "hash_map" is unused in this visitor implementation
	{
		if (v.type == 'f')
		{
			if (nn == -1) // nearest facility is found
			{
				nn = v.id; // nearest neighbor facility
				dnn = v.dist; // nearset neighbor distance, namely the NNFC radius
			}
			else // nn != -1, sub-nearest facility is found
			{
				// construct <f, \Delta(f)>
				delta_f_hash_map::iterator iter_find = ptr_f_delta_f_hash_map->find(nn); // iter_find (must) != end(), as initialized in "NSJ_construct_NNFC"
				iter_find->second += (dnn - v.dist) * prob; // \Delta(f) += (dnn - d2nn) * Pr(r), where r \in rnn(F,f)

				// compute max{dnn + d2nn} as f.rid
				rid_hash_map::iterator iter_rid = ptr_rid_hash_map->find(nn); // iter_rid (must) != end(), as initialized in "NSJ_construct_NNFC"
				if (iter_rid->second < dnn + v.dist)
					iter_rid->second = dnn + v.dist; // update max{dnn + d2nn}

				// calculate offsets coordinate for NNFC & 2-NNFC MBR
				float min_lon, min_lat, max_lon, max_lat, // for NNFC
					min_lon_2, min_lat_2, max_lon_2, max_lat_2; // for 2-NNFC
				geo_offset(lon, lat, dnn, WESTWARD, min_lon);
				geo_offset(lon, lat, dnn, SOUTHWARD, min_lat);
				geo_offset(lon, lat, dnn, EASTWARD, max_lon);
				geo_offset(lon, lat, dnn, NORTHWARD, max_lat);
				geo_offset(lon, lat, v.dist, WESTWARD, min_lon_2);
				geo_offset(lon, lat, v.dist, SOUTHWARD, min_lat_2);
				geo_offset(lon, lat, v.dist, EASTWARD, max_lon_2);
				geo_offset(lon, lat, v.dist, NORTHWARD, max_lat_2);

				(*ptr_NNFC_hash_map)[refloc_id] = // reference location id
					boost::make_tuple(prob, // present probability
					nn, // nn facility id
					(v.dist - dnn) * prob, // DIFF = (d2nn - dnn) * Pr(r)
					geo_point(lon, lat), // coordinate
					dnn, // dnn
					v.dist, // d2nn
					geo_box(geo_point(min_lon, min_lat), geo_point(max_lon, max_lat)), // NNFC MBR
					geo_box(geo_point(min_lon_2, min_lat_2), geo_point(max_lon_2, max_lat_2)), // 2-NNFC MBR
					0); // candidates being covered by 2-NNFC (HM.r)

				return true; // a flag to terminate graph traversal
			}
		}
		return false; // a flat to continue graph traversal
	}

	void set_data_structure(nnfc_hash_map *ptr_NNFC_hash_map_in, delta_f_hash_map *ptr_f_delta_f_hash_map_in,
		rid_hash_map *ptr_rid_hash_map_in) {
		ptr_NNFC_hash_map = ptr_NNFC_hash_map_in;
		ptr_f_delta_f_hash_map = ptr_f_delta_f_hash_map_in;
		ptr_rid_hash_map = ptr_rid_hash_map_in;
	}

protected:
	rid_hash_map *ptr_rid_hash_map; // a pointer to a hash map for every f.rid
};


//--------------------------------------------------------------------------------
// RID algorithm: compute dr(f), i.e., \Delta(f), construct NNFC MBRs for reference locations and RIC MBRs for facilities
// remark: the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
//		   also, the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve lon lat"
// [in] is_to_bidirectional: indicate whether to convert (true) undirected graph into bidirectional graph or not (false)
// [in] edges_file_path: file path of edges
// [in] facs_file_path: file path of facilities
// [in] reflocs_file_path: file path of reference locations
// [out] f_r_hash_map: a hash map for <f, \Delta(f), rnn(F,f)> tuples
// [out] reflocs: a hash map for the network location information of reference locations, which uses a "fac_hash_map" as alternative data structure
// [out] NNFCs; a hash map for NNFCs, which needs to be constructed in this function
// [out] RICs: a hash map for RICs, which needs to be constructed in this function
// [out] delta_F: a fibonacci heap for max-heap \Delta(f)s, the values are non-positive
// [out] graph: the typed directed graph, which needs to be constructed in this function

void RID_construct_NNFC_and_RIC(bool is_to_bidirectional, const char *edges_file_path, const char *facs_file_path, const char *reflocs_file_path,
	delta_f_hash_map &f_r_hash_map, fac_hash_map &reflocs, nnfc_hash_map &NNFCs, ric_hash_map &RICs, max_delta_fibonacci_heap &delta_F, typed_directed_graph &graph)
{
	// init graph
	graph.set_to_bidirectional(is_to_bidirectional);
	graph.init_graph(edges_file_path);

	// deploy facilities and initialize <f, \Delta(f), rnn(F,f)> tuples
	fac_ex_hash_map facs; // record all facilities with <lon, lat>
	rid_hash_map f_rids; // record f.rid for each facility
	graph.deploy_facilities_ex(facs_file_path, &facs);
	for (fac_ex_hash_map::iterator iter_f = facs.begin(); iter_f != facs.end(); ++iter_f)
	{
		f_r_hash_map[iter_f->first] = 0.0f; // must initialize the hash map, as some facility may affect no reference location
		f_rids[iter_f->first] = 0.0f; // initialize f.rid for every facility
	}

	if (!graph.get_is_facilities_deployed()) // check graph and facilities status
		return;

	// create a top event visitor for NNFCs and f.rid(s)
	NNFC_RIC_visitor top_visitor;
	top_visitor.set_data_structure(&NNFCs, &f_r_hash_map, &f_rids); // set data structure

	// compute NNFC/2NNFC of each reference location
	std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
	if (!ifs_reflocs.fail()) // reference locations file exists
	{
		while (!ifs_reflocs.bad() && ifs_reflocs.good())
		{
			char buf[1024];
			ifs_reflocs.getline(buf, 1024);
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int refloc_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			int vs_id = atoi(str_buf.substr(vs_pos + 1, ve_pos - vs_pos - 1).c_str());
			int ve_id = atoi(str_buf.substr(ve_pos + 1, dist_vs_pos - ve_pos - 1).c_str());
			float dist_vs = static_cast<float>(atof(str_buf.substr(dist_vs_pos + 1, dist_ve_pos - dist_vs_pos - 1).c_str()));
			float dist_ve = static_cast<float>(atof(str_buf.substr(dist_ve_pos + 1, prob_pos - dist_ve_pos - 1).c_str()));
			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// record network location information for each reference location
			reflocs[refloc_id] = boost::make_tuple(vs_id, ve_id, dist_vs, dist_ve);

			// insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(vs_id, 'v'), vertex(ve_id, 'v'), refloc_v, dist_vs, dist_ve, &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

			// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id, lon, lat, prob); // set id, geo-information and present probability of a reference location
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}
	}

	// order \Delta(f)s
	for (delta_f_hash_map::iterator iter_delta_f = f_r_hash_map.begin(); iter_delta_f != f_r_hash_map.end(); ++iter_delta_f)
		delta_F.push(delta_value(iter_delta_f->first, iter_delta_f->second)); // push to max-heap

#ifdef RID
	// calculate RICs
	for (rid_hash_map::iterator iter_f = f_rids.begin(); iter_f != f_rids.end(); ++iter_f)
	{
		fac_ex_hash_map::iterator iter_find_f = facs.find(iter_f->first); // must exist
		float f_lon = iter_find_f->second.get<F_LON>(); // lon of facility
		float f_lat = iter_find_f->second.get<F_LAT>(); // lat of facility
		float f_rid = iter_f->second; // f.rid = max{dnn + d2nn}
		if (f_rid == 0.0f) // the facility influence nobody
			continue; // don't create RIC for the facility
		float min_lon, min_lat, max_lon, max_lat; // for RIC
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, WESTWARD, min_lon);
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, SOUTHWARD, min_lat);
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, EASTWARD, max_lon);
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, NORTHWARD, max_lat);

		RICs[iter_f->first] = // facility id
			geo_box(geo_point(min_lon, min_lat), geo_point(max_lon, max_lat)); // RIC MBR
	}
#endif
}

//--------------------------------------------------------------------------------
// RID algorithm: query optimal candidate based on RICs
// remark: the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
// [re] facility_replacement: the id of optimal facility replacement pair
// [out] checked_num: the number of candidates actually checked
// [in/out] time_cost: the time cost (ms) of querying
// [in] f_r_hash_map: a hash map for <f, \Delta(f), rnn(F,f)> tuples
// [in] reflocs: a hash map for the network location information of reference locations, which uses a "fac_hash_map" as alternative data structure
// [in/out] NNFCs; a hash map for NNFCs, COVER_NUM element will be changed
// [in] RICs: a hash map for RICs
// [in] delta_F: a fibonacci heap for max-heap \Delta(f)s, the values are non-positive
// [in] delta_C: 
// [in/out] graph: the constructed graph, candidates will be further deployed
// [in] cands_file_path: file path of candidates

facility_replacement RID_query(int &checked_num, __int64 &time_cost, delta_f_hash_map &f_r_hash_map, std::pair<int, float> &max_delta_f,
	 fac_hash_map &reflocs, nnfc_hash_map &NNFCs, ric_hash_map &RICs, max_delta_fibonacci_heap &delta_F,
	typed_directed_graph &graph, const char *cands_file_path)
{
	auto begin_time = std::chrono::high_resolution_clock::now();
	boost::unordered_map<int, geo_point> cand_geos;
	geo_cand_rtree cand_rtree;
	cover_hash_map Rc_map;
	// deploy candidates onto the directed graph and construct candidates R*-tree
	graph.deploy_candidates_rtree(cands_file_path, cand_rtree, cand_geos);
	
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join for \Delta+(c)s (alternative by R*-tree range query): construct hash maps for <c, <\Delta+(c), Rc'> >
	// sets <c, <\Delta+(c), Rc'> >
	for (nnfc_hash_map::iterator iter_NNFC = NNFCs.begin(); iter_NNFC != NNFCs.end(); ++iter_NNFC)
	{
		// range query for candidates which are covered by a 2-NNFC
		std::vector<geo_cand> covered_cands; // covered by 2-NNFC
		cand_rtree.query(boost::geometry::index::within(iter_NNFC->second.get<NNFC_2_MBR>()), std::back_inserter(covered_cands));
		iter_NNFC->second.get<COVER_NUM>() = static_cast<int>(covered_cands.size()); // record the maximum number of (potentially) covered (by 2-NNFC) candidates

																					 // record the set of ids of reference locations whose 2-NNFCs cover the candidate, namely Rc'
		for (std::vector<geo_cand>::iterator iter_cand = covered_cands.begin(); iter_cand != covered_cands.end(); ++iter_cand)
		{
			// Rc+': inside NNFC
			if (boost::geometry::within(iter_cand->first, iter_NNFC->second.get<NNFC_MBR>()))
			{
				// compute Ed+ for pair <r,c>
				float distE = boost::geometry::distance(iter_NNFC->second.get<COORDINATE>(), iter_cand->first) * EARTH_RADIUS; // dE(r,c)
				float Ed_u = (iter_NNFC->second.get<DNN>() - distE) * iter_NNFC->second.get<PROB>(); // Ed+r(c) = (dnn - dE(r,c)) * Pr(r)
				if (Ed_u > 0.0f) // possible benefit
								 // >: candidate locates inside NNFC
								 // <: candidate locates inside the 4 corners between NNFC and its MBR
								 // ==: exactly on the NNFC bound, maybe < for road network distance
				{
					// \Delta+(c)'
					cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
					if (iter_Rc_find != Rc_map.end()) // candidate has already been created
					{
						iter_Rc_find->second.first += Ed_u; // \Delta+(c) += (dnn - dE(r,c)) * Pr(r)
						iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					else // need to create a set for the candidate
					{
						Rc_map[iter_cand->second] = std::make_pair(Ed_u, std::set<int>()); // Ed_u as the initial \Delta+(c)' value
						Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					continue; // next candidate
				}
			}

			// maybe inside 2-NNFC ring
			cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
			if (iter_Rc_find != Rc_map.end()) // candidate has already been created
				iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			else // need to create a set for the candidate
			{
				Rc_map[iter_cand->second] = std::make_pair(0.0f, std::set<int>()); // 0 as the initial \Delta+(c)' value, as at most inside 2-NNFC ring
				Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			}
		}
	}

	// order \Delta+(c)s
	max_delta_fibonacci_heap delta_C; // a fibonacci heap for max - heap \Delta + (c)s
	for (cover_hash_map::iterator iter_Rc_c = Rc_map.begin(); iter_Rc_c != Rc_map.end(); ++iter_Rc_c)
		delta_C.push(delta_value(iter_Rc_c->first, iter_Rc_c->second.first));

#ifdef RID
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join for RICs (alternative by R*-tree range query): record <f,c> pairs which need validation
	boost::unordered_set<facility_replacement> f_c_pairs_need_valid; // <f,c> pairs requiring validation
	for (ric_hash_map::iterator iter_RIC = RICs.begin(); iter_RIC != RICs.end(); ++iter_RIC)
	{
		// range query for candidates which are covered by an RIC
		std::vector<geo_cand> covered_cands_by_RIC; // covered by RIC
		cand_rtree.query(boost::geometry::index::within(iter_RIC->second), std::back_inserter(covered_cands_by_RIC));

		// record the <f,c> pairs which need validation
		for (std::vector<geo_cand>::iterator iter_cand_by_RIC = covered_cands_by_RIC.begin(); iter_cand_by_RIC != covered_cands_by_RIC.end(); ++iter_cand_by_RIC)
			f_c_pairs_need_valid.insert(facility_replacement(iter_RIC->first, iter_cand_by_RIC->second, 0.0f)); // a pair needs validation
	}
#endif

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// iterate Max-heaps of \Delta(f) and \Delta+(c) as <f,c> pairs, and validate if necessary
	facility_replacement op(-1, -1, 0.0f); // initialize optimal facility-replacement pair
	checked_num = 0; // initialize

					 // create a top event visitor for RID algorithm
	NSJ_visitor top_visitor;
	offset_c_f_hash_map offset_c_f; // offset for \Delta(f, c)
	delta_c_hash_map delta_c_road; // \Delta(c), not \Delta+(c) (the corresponding variable is delta_C)
	top_visitor.set_data_structures(&NNFCs, &offset_c_f, &delta_c_road); // set data structures

	while (!delta_C.empty())
	{
		// retrieve the top candidate in the \Delta+(c) max-heap
		delta_value top_c = delta_C.top();

		if (top_c.val < op.ED) // current op <f,c> is the optimal facility-replacement pair
			break;

		// traverse each reference location whose 2-NNFC covers top_c
		for (std::set<int>::iterator iter_refloc = Rc_map[top_c.id].second.begin(); iter_refloc != Rc_map[top_c.id].second.end(); ++iter_refloc) // iterate each r \in Rc'
		{
			int refloc_id = *iter_refloc;
			if (NNFCs[refloc_id].get<COVER_NUM>() == 0) // all candidates that can affect the reference location have been traversed by a certain candidate prior to this top_c
				continue; // traverse from next reference location

						  // insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(reflocs[refloc_id].get<VS_ID>(), 'v'), vertex(reflocs[refloc_id].get<VE_ID>(), 'v'), refloc_v,
				reflocs[refloc_id].get<VS_DIST>(), reflocs[refloc_id].get<VE_DIST>(), &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

						// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id); // now the ED of top <f,c> is initialized as \Delta(f)
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

		for (max_delta_fibonacci_heap::ordered_iterator iter_ordered_f = delta_F.ordered_begin(); iter_ordered_f != delta_F.ordered_end(); ++iter_ordered_f)
		{
			float delta = delta_c_road[top_c.id] + iter_ordered_f->val; // \Delta(c) + \Delta(f) (un-offset value, which may be less than the actual ED due to the interference)
			++checked_num; // need to check more facility-replacement pairs
#ifdef RID
						   // check if validation is needed
			facility_replacement f_c_checking(iter_ordered_f->id, top_c.id, delta); // the current checking <f,c> pair and its un-offset ED
			if (f_c_pairs_need_valid.find(f_c_checking) == f_c_pairs_need_valid.end()) // no interference, then no need to offset
			{
				if (delta > op.ED) // here delta is the actual ED
				{
					op.f_id = iter_ordered_f->id;
					op.c_id = top_c.id;
					op.ED = delta;
				}
			}
			else // need offset
			{
				float delta_with_offset = delta  // \Delta(c) + \Delta(f)
					+ offset_c_f[top_c.id][iter_ordered_f->id]; // + offset of \Delta(<f,c>)

				if (delta_with_offset > op.ED)
				{
					op.f_id = iter_ordered_f->id;
					op.c_id = top_c.id;
					op.ED = delta_with_offset;
				}
			}
#endif		
		}

		delta_C.pop(); // pop the top candidate
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	time_cost = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count(); // iteration of Max-heap time cost

	return op; // the optimal facility-replacement pair
}

//Greedy RIQ using
/**
facility_replacement RID_query(int &checked_num, __int64 &time_cost, delta_f_hash_map &f_r_hash_map, std::pair<int, float> &max_delta_f, cover_hash_map &Rc_map, geo_cand_rtree &cand_rtree,
	boost::unordered_map<int, geo_point> &cand_geos,fac_hash_map &reflocs, nnfc_hash_map &NNFCs, ric_hash_map &RICs, max_delta_fibonacci_heap &delta_F,
	typed_directed_graph &graph, const char *cands_file_path)
{
	auto begin_time = std::chrono::high_resolution_clock::now();
	
	
	
	// deploy candidates onto the directed graph and construct candidates R*-tree
	graph.deploy_candidates_rtree(cands_file_path, cand_rtree, cand_geos);

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join for \Delta+(c)s (alternative by R*-tree range query): construct hash maps for <c, <\Delta+(c), Rc'> >
	// sets <c, <\Delta+(c), Rc'> >
	for (nnfc_hash_map::iterator iter_NNFC = NNFCs.begin(); iter_NNFC != NNFCs.end(); ++iter_NNFC)
	{
		// range query for candidates which are covered by a 2-NNFC
		std::vector<geo_cand> covered_cands; // covered by 2-NNFC
		cand_rtree.query(boost::geometry::index::within(iter_NNFC->second.get<NNFC_2_MBR>()), std::back_inserter(covered_cands));
		iter_NNFC->second.get<COVER_NUM>() = static_cast<int>(covered_cands.size()); // record the maximum number of (potentially) covered (by 2-NNFC) candidates

																					 // record the set of ids of reference locations whose 2-NNFCs cover the candidate, namely Rc'
		for (std::vector<geo_cand>::iterator iter_cand = covered_cands.begin(); iter_cand != covered_cands.end(); ++iter_cand)
		{
			// Rc+': inside NNFC
			if (boost::geometry::within(iter_cand->first, iter_NNFC->second.get<NNFC_MBR>()))
			{
				// compute Ed+ for pair <r,c>
				float distE = boost::geometry::distance(iter_NNFC->second.get<COORDINATE>(), iter_cand->first) * EARTH_RADIUS; // dE(r,c)
				float Ed_u = (iter_NNFC->second.get<DNN>() - distE) * iter_NNFC->second.get<PROB>(); // Ed+r(c) = (dnn - dE(r,c)) * Pr(r)
				if (Ed_u > 0.0f) // possible benefit
								 // >: candidate locates inside NNFC
								 // <: candidate locates inside the 4 corners between NNFC and its MBR
								 // ==: exactly on the NNFC bound, maybe < for road network distance
				{
					// \Delta+(c)'
					cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
					if (iter_Rc_find != Rc_map.end()) // candidate has already been created
					{
						iter_Rc_find->second.first += Ed_u; // \Delta+(c) += (dnn - dE(r,c)) * Pr(r)
						iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					else // need to create a set for the candidate
					{
						Rc_map[iter_cand->second] = std::make_pair(Ed_u, std::set<int>()); // Ed_u as the initial \Delta+(c)' value
						Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					continue; // next candidate
				}
			}

			// maybe inside 2-NNFC ring
			cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
			if (iter_Rc_find != Rc_map.end()) // candidate has already been created
				iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			else // need to create a set for the candidate
			{
				Rc_map[iter_cand->second] = std::make_pair(0.0f, std::set<int>()); // 0 as the initial \Delta+(c)' value, as at most inside 2-NNFC ring
				Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			}
		}
	}

	// order \Delta+(c)s
	max_delta_fibonacci_heap delta_C; // a fibonacci heap for max - heap \Delta + (c)s
	for (cover_hash_map::iterator iter_Rc_c = Rc_map.begin(); iter_Rc_c != Rc_map.end(); ++iter_Rc_c)
		delta_C.push(delta_value(iter_Rc_c->first, iter_Rc_c->second.first));

#ifdef RID
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join for RICs (alternative by R*-tree range query): record <f,c> pairs which need validation
	boost::unordered_set<facility_replacement> f_c_pairs_need_valid; // <f,c> pairs requiring validation
	for (ric_hash_map::iterator iter_RIC = RICs.begin(); iter_RIC != RICs.end(); ++iter_RIC)
	{
		// range query for candidates which are covered by an RIC
		std::vector<geo_cand> covered_cands_by_RIC; // covered by RIC
		cand_rtree.query(boost::geometry::index::within(iter_RIC->second), std::back_inserter(covered_cands_by_RIC));

		// record the <f,c> pairs which need validation
		for (std::vector<geo_cand>::iterator iter_cand_by_RIC = covered_cands_by_RIC.begin(); iter_cand_by_RIC != covered_cands_by_RIC.end(); ++iter_cand_by_RIC)
			f_c_pairs_need_valid.insert(facility_replacement(iter_RIC->first, iter_cand_by_RIC->second, 0.0f)); // a pair needs validation
	}
#endif

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// iterate Max-heaps of \Delta(f) and \Delta+(c) as <f,c> pairs, and validate if necessary
	facility_replacement op(-1, -1, 0.0f); // initialize optimal facility-replacement pair
	checked_num = 0; // initialize

					 // create a top event visitor for RID algorithm
	NSJ_visitor top_visitor;
	offset_c_f_hash_map offset_c_f; // offset for \Delta(f, c)
	delta_c_hash_map delta_c_road; // \Delta(c), not \Delta+(c) (the corresponding variable is delta_C)
	top_visitor.set_data_structures(&NNFCs, &offset_c_f, &delta_c_road); // set data structures

	while (!delta_C.empty())
	{
		// retrieve the top candidate in the \Delta+(c) max-heap
		delta_value top_c = delta_C.top();

		if (top_c.val < op.ED) // current op <f,c> is the optimal facility-replacement pair
			break;

		// traverse each reference location whose 2-NNFC covers top_c
		for (std::set<int>::iterator iter_refloc = Rc_map[top_c.id].second.begin(); iter_refloc != Rc_map[top_c.id].second.end(); ++iter_refloc) // iterate each r \in Rc'
		{
			int refloc_id = *iter_refloc;
			if (NNFCs[refloc_id].get<COVER_NUM>() == 0) // all candidates that can affect the reference location have been traversed by a certain candidate prior to this top_c
				continue; // traverse from next reference location

						  // insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(reflocs[refloc_id].get<VS_ID>(), 'v'), vertex(reflocs[refloc_id].get<VE_ID>(), 'v'), refloc_v,
				reflocs[refloc_id].get<VS_DIST>(), reflocs[refloc_id].get<VE_DIST>(), &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

						// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id); // now the ED of top <f,c> is initialized as \Delta(f)
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

		for (max_delta_fibonacci_heap::ordered_iterator iter_ordered_f = delta_F.ordered_begin(); iter_ordered_f != delta_F.ordered_end(); ++iter_ordered_f)
		{
			float delta = delta_c_road[top_c.id] + iter_ordered_f->val; // \Delta(c) + \Delta(f) (un-offset value, which may be less than the actual ED due to the interference)
			++checked_num; // need to check more facility-replacement pairs
#ifdef RID
						   // check if validation is needed
			facility_replacement f_c_checking(iter_ordered_f->id, top_c.id, delta); // the current checking <f,c> pair and its un-offset ED
			if (f_c_pairs_need_valid.find(f_c_checking) == f_c_pairs_need_valid.end()) // no interference, then no need to offset
			{
				if (delta > op.ED) // here delta is the actual ED
				{
					op.f_id = iter_ordered_f->id;
					op.c_id = top_c.id;
					op.ED = delta;
				}
			}
			else // need offset
			{
				float delta_with_offset = delta  // \Delta(c) + \Delta(f)
					+ offset_c_f[top_c.id][iter_ordered_f->id]; // + offset of \Delta(<f,c>)

				if (delta_with_offset > op.ED)
				{
					op.f_id = iter_ordered_f->id;
					op.c_id = top_c.id;
					op.ED = delta_with_offset;
				}
			}
#endif		
		}

		delta_C.pop(); // pop the top candidate
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	time_cost = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count(); // iteration of Max-heap time cost

	return op; // the optimal facility-replacement pair
}
*/
void RID_update_NNFC_and_RIC(typed_directed_graph &graph, delta_f_hash_map &f_r_hash_map, nnfc_hash_map &NNFCs, ric_hash_map &RICs, max_delta_fibonacci_heap &delta_F, facility_replacement &op_RID, cover_hash_map &Rc_map, int &i)
{


	vertex f_v(op_RID.f_id, 'f');
	vertex v_v(21047 + i, 'v');
	graph.replace_vertex(f_v, v_v, NULL, NULL);
	vertex c_v(op_RID.c_id, 'c');
	vertex ff_v(op_RID.f_id, 'f');
	graph.replace_vertex(c_v, ff_v, NULL, NULL);

	fac_ex_hash_map facs;
	rid_hash_map f_rids;
	NNFC_RIC_visitor top_visitor;
	top_visitor.set_data_structure(&NNFCs, &f_r_hash_map, &f_rids); // set data structure
	for (nnfc_hash_map::iterator iter_NNFC = NNFCs.begin(); iter_NNFC != NNFCs.end(); ++iter_NNFC)
	{
		if (boost::geometry::within(iter_NNFC->second.get<3>(), RICs[op_RID.f_id]))
		{
			vertex refloc_v(iter_NNFC->first, 'r');
			graph.dijkstra(refloc_v, &top_visitor);
		}
	}
	for (std::set<int>::iterator ite1 = Rc_map[op_RID.c_id].second.begin(); ite1 != Rc_map[op_RID.c_id].second.end(); ite1++)
	{

		vertex refloc_v(*ite1, 'r');
		graph.dijkstra(refloc_v, &top_visitor);
	}
	//-----------------------------------------------------------------------

	for (delta_f_hash_map::iterator iter_delta_f = f_r_hash_map.begin(); iter_delta_f != f_r_hash_map.end(); ++iter_delta_f)
		delta_F.push(delta_value(iter_delta_f->first, iter_delta_f->second));
#ifdef RID
	// calculate RICs
	for (rid_hash_map::iterator iter_f = f_rids.begin(); iter_f != f_rids.end(); ++iter_f)
	{
		fac_ex_hash_map::iterator iter_find_f = facs.find(iter_f->first); // must exist
		float f_lon = iter_find_f->second.get<F_LON>(); // lon of facility
		float f_lat = iter_find_f->second.get<F_LAT>(); // lat of facility
		float f_rid = iter_f->second; // f.rid = max{dnn + d2nn}
		if (f_rid == 0.0f) // the facility influence nobody
			continue; // don't create RIC for the facility
		float min_lon, min_lat, max_lon, max_lat; // for RIC
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, WESTWARD, min_lon);
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, SOUTHWARD, min_lat);
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, EASTWARD, max_lon);
		NNFC_visitor::geo_offset(f_lon, f_lat, f_rid, NORTHWARD, max_lat);

		RICs[iter_f->first] = // facility id
			geo_box(geo_point(min_lon, min_lat), geo_point(max_lon, max_lat)); // RIC MBR
	}
#endif
}

facility_replacement RID_new_query(int &checked_num, __int64 &time_cost, delta_f_hash_map &f_r_hash_map, std::pair<int, float> &max_delta_f, cover_hash_map &Rc_map,
	fac_hash_map &reflocs, nnfc_hash_map &NNFCs, ric_hash_map &RICs, max_delta_fibonacci_heap &delta_F,typed_directed_graph &graph, facility_replacement &op_RID, geo_cand_rtree &cand_rtree, boost::unordered_map<int, geo_point> &cand_geos)
{
	auto begin_time = std::chrono::high_resolution_clock::now();
	for (boost::unordered_map<int, geo_point>::iterator iter_c = cand_geos.begin(); iter_c != cand_geos.end(); ++iter_c)
	{
		if (iter_c->first == op_RID.c_id)
		{
			cand_rtree.remove(geo_cand(iter_c->second, op_RID.c_id));
			break;
		}
	}

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join for \Delta+(c)s (alternative by R*-tree range query): construct hash maps for <c, <\Delta+(c), Rc'> >

	for (nnfc_hash_map::iterator iter_NNFC = NNFCs.begin(); iter_NNFC != NNFCs.end(); ++iter_NNFC)
	{
		// range query for candidates which are covered by a 2-NNFC
		std::vector<geo_cand> covered_cands; // covered by 2-NNFC
		cand_rtree.query(boost::geometry::index::within(iter_NNFC->second.get<NNFC_2_MBR>()), std::back_inserter(covered_cands));
		iter_NNFC->second.get<COVER_NUM>() = static_cast<int>(covered_cands.size()); // record the maximum number of (potentially) covered (by 2-NNFC) candidates

																					 // record the set of ids of reference locations whose 2-NNFCs cover the candidate, namely Rc'
		for (std::vector<geo_cand>::iterator iter_cand = covered_cands.begin(); iter_cand != covered_cands.end(); ++iter_cand)
		{
			// Rc+': inside NNFC
			if (boost::geometry::within(iter_cand->first, iter_NNFC->second.get<NNFC_MBR>()))
			{
				// compute Ed+ for pair <r,c>
				float distE = boost::geometry::distance(iter_NNFC->second.get<COORDINATE>(), iter_cand->first) * EARTH_RADIUS; // dE(r,c)
				float Ed_u = (iter_NNFC->second.get<DNN>() - distE) * iter_NNFC->second.get<PROB>(); // Ed+r(c) = (dnn - dE(r,c)) * Pr(r)
				if (Ed_u > 0.0f) // possible benefit
								 // >: candidate locates inside NNFC
								 // <: candidate locates inside the 4 corners between NNFC and its MBR
								 // ==: exactly on the NNFC bound, maybe < for road network distance
				{
					// \Delta+(c)'
					cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
					if (iter_Rc_find != Rc_map.end()) // candidate has already been created
					{
						iter_Rc_find->second.first += Ed_u; // \Delta+(c) += (dnn - dE(r,c)) * Pr(r)
						iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					else // need to create a set for the candidate
					{
						Rc_map[iter_cand->second] = std::make_pair(Ed_u, std::set<int>()); // Ed_u as the initial \Delta+(c)' value
						Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose NNFC (i.e., 2-NNFC) covers the candidate
					}
					continue; // next candidate
				}
			}

			// maybe inside 2-NNFC ring
			cover_hash_map::iterator iter_Rc_find = Rc_map.find(iter_cand->second); // find candidate
			if (iter_Rc_find != Rc_map.end()) // candidate has already been created
				iter_Rc_find->second.second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			else // need to create a set for the candidate
			{
				Rc_map[iter_cand->second] = std::make_pair(0.0f, std::set<int>()); // 0 as the initial \Delta+(c)' value, as at most inside 2-NNFC ring
				Rc_map[iter_cand->second].second.insert(iter_NNFC->first); // record id of reference location whose 2-NNFC covers the candidate
			}
		}
	}

	// order \Delta+(c)s
	max_delta_fibonacci_heap delta_C; // a fibonacci heap for max - heap \Delta + (c)s
	for (cover_hash_map::iterator iter_Rc_c = Rc_map.begin(); iter_Rc_c != Rc_map.end(); ++iter_Rc_c)
		delta_C.push(delta_value(iter_Rc_c->first, iter_Rc_c->second.first));

#ifdef RID
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// Spatial Join for RICs (alternative by R*-tree range query): record <f,c> pairs which need validation
	boost::unordered_set<facility_replacement> f_c_pairs_need_valid; // <f,c> pairs requiring validation
	for (ric_hash_map::iterator iter_RIC = RICs.begin(); iter_RIC != RICs.end(); ++iter_RIC)
	{
		// range query for candidates which are covered by an RIC
		std::vector<geo_cand> covered_cands_by_RIC; // covered by RIC
		cand_rtree.query(boost::geometry::index::within(iter_RIC->second), std::back_inserter(covered_cands_by_RIC));

		// record the <f,c> pairs which need validation
		for (std::vector<geo_cand>::iterator iter_cand_by_RIC = covered_cands_by_RIC.begin(); iter_cand_by_RIC != covered_cands_by_RIC.end(); ++iter_cand_by_RIC)
			f_c_pairs_need_valid.insert(facility_replacement(iter_RIC->first, iter_cand_by_RIC->second, 0.0f)); // a pair needs validation
	}
#endif

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// iterate Max-heaps of \Delta(f) and \Delta+(c) as <f,c> pairs, and validate if necessary
	facility_replacement op(-1, -1, 0.0f); // initialize optimal facility-replacement pair
	checked_num = 0; // initialize

					 // create a top event visitor for RID algorithm
	NSJ_visitor top_visitor;
	offset_c_f_hash_map offset_c_f; // offset for \Delta(f, c)
	delta_c_hash_map delta_c_road; // \Delta(c), not \Delta+(c) (the corresponding variable is delta_C)
	top_visitor.set_data_structures(&NNFCs, &offset_c_f, &delta_c_road); // set data structures

	while (!delta_C.empty())
	{
		// retrieve the top candidate in the \Delta+(c) max-heap
		delta_value top_c = delta_C.top();

		if (top_c.val < op.ED) // current op <f,c> is the optimal facility-replacement pair
			break;

		// traverse each reference location whose 2-NNFC covers top_c
		for (std::set<int>::iterator iter_refloc = Rc_map[top_c.id].second.begin(); iter_refloc != Rc_map[top_c.id].second.end(); ++iter_refloc) // iterate each r \in Rc'
		{
			int refloc_id = *iter_refloc;
			if (NNFCs[refloc_id].get<COVER_NUM>() == 0) // all candidates that can affect the reference location have been traversed by a certain candidate prior to this top_c
				continue; // traverse from next reference location

						  // insert a reference location vertex to <vs, ve>
			vertex refloc_v(refloc_id, 'r');
			std::vector<directed_edge> removed_edges, inserted_edges;
			graph.insert_vertex(vertex(reflocs[refloc_id].get<VS_ID>(), 'v'), vertex(reflocs[refloc_id].get<VE_ID>(), 'v'), refloc_v,
				reflocs[refloc_id].get<VS_DIST>(), reflocs[refloc_id].get<VE_DIST>(), &removed_edges, &inserted_edges,
				false); // "false" means not to replace overlapped vertex for reference location, we will create a virtual vertex for it

						// traverse graph from the reference location by dijkstra algorithm
			top_visitor.set_refloc(refloc_id); // now the ED of top <f,c> is initialized as \Delta(f)
			graph.dijkstra(refloc_v, &top_visitor);

			// restore the graph by removing the inserted reference location vertex
			graph.restore_graph(refloc_v, removed_edges, inserted_edges);
		}

		for (max_delta_fibonacci_heap::ordered_iterator iter_ordered_f = delta_F.ordered_begin(); iter_ordered_f != delta_F.ordered_end(); ++iter_ordered_f)
		{
			float delta = delta_c_road[top_c.id] + iter_ordered_f->val; // \Delta(c) + \Delta(f) (un-offset value, which may be less than the actual ED due to the interference)
			++checked_num; // need to check more facility-replacement pairs
#ifdef RID
						   // check if validation is needed
			facility_replacement f_c_checking(iter_ordered_f->id, top_c.id, delta); // the current checking <f,c> pair and its un-offset ED
			if (f_c_pairs_need_valid.find(f_c_checking) == f_c_pairs_need_valid.end()) // no interference, then no need to offset
			{
				if (delta > op.ED) // here delta is the actual ED
				{
					op.f_id = iter_ordered_f->id;
					op.c_id = top_c.id;
					op.ED = delta;
				}
			}
			else // need offset
			{
				float delta_with_offset = delta  // \Delta(c) + \Delta(f)
					+ offset_c_f[top_c.id][iter_ordered_f->id]; // + offset of \Delta(<f,c>)

				if (delta_with_offset > op.ED)
				{
					op.f_id = iter_ordered_f->id;
					op.c_id = top_c.id;
					op.ED = delta_with_offset;
				}
			}
#endif		
		}

		delta_C.pop(); // pop the top candidate
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	time_cost = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count(); // iteration of Max-heap time cost

	return op; // the optimal facility-replacement pair
}

#pragma endregion ALGO_RID