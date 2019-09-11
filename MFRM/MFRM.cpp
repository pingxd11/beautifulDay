#include "stdafx.h"

#include <chrono> // std::chrono::high_resolution_clock
#include <string> // std::string
#include <fstream> // std::ifstream, std::ofstream
#include <iostream> // std::cout

#include "ir_utility.h"
//#include "greedy_k.h" // need to un-comment greedy_k part and comment out all solutions part
#include "solutions.h" // needed for all solutions and compare_EN_with_MinDist
#include "mindist.h"
#include "gen_datasets.h"
#include "convert_datasets.h"
#include "real_datasets.h"

// remark: use "#define RID" for facilities with <lon, lat> or comment it in "typed_directed_graph.h", which is included by "solutions.h"

// argv[2]: config file for datasets files
bool gen_datasets(_TCHAR* argv); // forward declaration
bool convert_datasets(_TCHAR* argv); // forward declaration
bool real_datasets(_TCHAR* argv); // forward declaration
int compare_EN_with_MinDist(_TCHAR* argv[])
{
	// read config file for datasets files
	bool is_directed_graph = false; // default is bidirected graph
	std::string edges_file_path, facs_file_path, facs_xy_file_path, reflocs_avg_file_path, reflocs_file_path, cands_file_path;
	std::string results_EN_file_path, results_MD_file_path;
	std::ifstream ifs_config(argv[2]);
	if (!ifs_config.fail())	// Config file exists.
	{
		while (!ifs_config.bad() && ifs_config.good())
		{
			char buf[1024];
			ifs_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "is_directed_graph")
				is_directed_graph = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "edges")
				edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "facs")
				facs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "facs_xy")
				facs_xy_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "reflocs_avg")
				reflocs_avg_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "reflocs")
				reflocs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "cands")
				cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "results_EN")
				results_EN_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "results_MD")
				results_MD_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
		}
	}

	// EN
	std::ofstream ofs_re_EN(results_EN_file_path, std::ofstream::out | std::ofstream::trunc); // create results file for EN
	if (ofs_re_EN.fail())	// creating results file for EN fails
		return -1;
	__int64 traversal_time; // here, this variable is just to complete param list
	std::vector<facility_replacement> topk;
	EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
		edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str(), traversal_time,
		0, // 0 means output all FR pairs
		&topk); // for all facility-replacement pairs
	int size = topk.size();
	for (int i = 0; i < size; ++i)
		ofs_re_EN << "<" << topk[i].f_id << ", " << topk[i].c_id << ">\t" << topk[i].ED << '\n';
	ofs_re_EN.flush();

	// Min-dist
	std::ofstream ofs_re_MD(results_MD_file_path, std::ofstream::out | std::ofstream::trunc); // create results file for Min-Dist
	if (ofs_re_MD.fail())	// creating results file for Min-dist fails
		return -1;
	std::vector<geo_point> users;
	std::map<int, geo_point> facs;
	geo_cand_rtree facs_tree;
	std::vector<geo_point> cands;
	load_users_facs_cands(reflocs_avg_file_path.c_str(), reflocs_file_path.c_str(), facs_file_path.c_str(), cands_file_path.c_str(),
		users, facs, facs_tree, cands);
	std::vector<facility_replacement> FRs;
	float total_dist = 0.0f;
	MINDIST(FRs, total_dist, users, facs, facs_tree, cands);
	ofs_re_MD << "total: " << total_dist << '\n';
	for (int i = 0; i < FRs.size(); ++i)
		ofs_re_MD << "<" << FRs[i].f_id << ", " << FRs[i].c_id << ">\t" << FRs[i].ED << '\n';
	ofs_re_MD.flush();

	return 0;
}

int compare_EN_with_MinDist(_TCHAR* argv[])
{
	std::string CA_file_path, BJ_file_path;
	std::ifstream ifs_config(argv[2]);
	if (!ifs_config.fail())	// Config file exists.
	{
		while (!ifs_config.bad() && ifs_config.good())
		{
			char buf[1024];
			ifs_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "CA")
				CA_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "BJ")
				BJ_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
		}
	}
	std::cout << CA_file_path << BJ_file_path;
	std::string edges_file_path("gen_edges.txt");
	std::string reflocs_avg_file_path[] = { "gen_avg_reflocs_20000.txt",
											"gen_avg_reflocs_60000.txt",
											"gen_avg_reflocs_100000.txt" };
	std::string reflocs_file_path[] = { "gen_rand_reflocs_20000.txt",
										"gen_rand_reflocs_60000.txt",
										"gen_rand_reflocs_100000.txt" };
	std::string reflocs[] = { "_20k", "_60k", "_100k" };
	std::string cands_file_path[] = {   "gen_cands_100a.txt",
										"gen_cands_100b.txt",
										"gen_cands_100c.txt",
										"gen_cands_100d.txt",
										"gen_cands_100e.txt",
										"gen_cands_100f.txt",
										"gen_cands_100g.txt",
										"gen_cands_100h.txt",
										"gen_cands_100i.txt",
										"gen_cands_100j.txt" ,
		                                "gen_cands_100k.txt",
		                                "gen_cands_100l.txt",
		                                "gen_cands_100m.txt",
		                                "gen_cands_100n.txt",
		                                "gen_cands_100o.txt",
		                                "gen_cands_100p.txt",
		                                "gen_cands_100q.txt",
		                                "gen_cands_100r.txt",
		                                "gen_cands_100s.txt",
		                                "gen_cands_100t.txt" };
	std::string cands[] = { "_0", "_1", "_2", "_3", "_4", "_5", "_6", "_7", "_8", "_9","_10","_11","_12","_13","_14","_15","_16","_17","_18","_19","_20" };
	
	// CA
	bool is_directed_graph = false; // CA is bidirected graph
									
	std::string CA_facs_file_path[] = { "800_Gas_xy.txt",
									    "800_pos_xy.txt",
									    "800_ATM_xy.txt",
									    "800_Sub_xy.txt" }; 
	std::string CA_facs[] = { "_GA", "_PO", "_AT", "_SU" };


	for (int iref = 0; iref < 1; ++iref) // only 20k
	{
		std::vector<geo_point> users;
		md_load_realusers((CA_file_path + reflocs_file_path[iref]).c_str(), users);
		//md_load_users((CA_file_path + reflocs_avg_file_path[iref]).c_str(), (CA_file_path + reflocs_file_path[iref]).c_str(), users);
		for (int ifac = 0; ifac < 4; ++ifac)
		{
			std::map<int, geo_point> facs;
			geo_cand_rtree facs_tree;
			boost::unordered_map<int, float> fac_neg_effect;
			md_load_facs((CA_file_path + CA_facs_file_path[ifac]).c_str(), facs, facs_tree, fac_neg_effect);

			for (int icand = 0; icand < 20; ++icand)
			{
				std::string results_EN_file_path(CA_file_path);
				results_EN_file_path += CA_facs[ifac] + reflocs[iref] + cands[icand];
				std::string results_MD_file_path(results_EN_file_path);
				results_EN_file_path += "_EN.txt";
				results_MD_file_path += "_MD.txt";

				// EN
				std::cout << results_EN_file_path << std::endl;
				std::ofstream ofs_re_EN(results_EN_file_path, std::ofstream::out | std::ofstream::trunc); // create results file for EN
				if (ofs_re_EN.fail())	// creating results file for EN fails
					return -1;
				__int64 traversal_time; // here, this variable is just to complete param list
				std::vector<facility_replacement> topk;
				EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
					(CA_file_path + edges_file_path).c_str(),
					(CA_file_path + CA_facs_file_path[ifac]).c_str(),
					(CA_file_path + reflocs_file_path[iref]).c_str(),
					(CA_file_path + cands_file_path[icand]).c_str(),
					traversal_time,
					0, // 0 means output all FR pairs
					&topk); // for all facility-replacement pairs
				int size = topk.size();
				for (int i = 0; i < size; ++i)
					ofs_re_EN << "<" << topk[i].f_id << ", " << topk[i].c_id << ">\t" << topk[i].ED << '\n';
				ofs_re_EN.flush();

				// Min-dist
				std::cout << results_MD_file_path << std::endl;
				std::ofstream ofs_re_MD(results_MD_file_path, std::ofstream::out | std::ofstream::trunc); // create results file for Min-Dist
				if (ofs_re_MD.fail())	// creating results file for Min-dist fails
					return -1;
				std::vector<geo_point> cands;
				geo_cand_rtree cands_tree;
				boost::unordered_map<int, float> cand_benefit;
				md_load_cands((CA_file_path + cands_file_path[icand]).c_str(), cands, cands_tree, cand_benefit);
				std::vector<facility_replacement> FRs;
				//float total_dist = 0.0f;
				//MINDIST(FRs, total_dist, users, facs, facs_tree, cands);
				//ofs_re_MD << "total: " << total_dist << '\n';
				MINDIST_FAST(FRs, users, facs, facs_tree, fac_neg_effect, cands, cands_tree, cand_benefit);
				for (int i = 0; i < FRs.size(); ++i)
					ofs_re_MD << "<" << FRs[i].f_id << ", " << FRs[i].c_id << ">\t" << FRs[i].ED << '\n';
				ofs_re_MD.flush();
			}
		}
	}
	
	
	is_directed_graph = true; // BJ is directed graph
	std::string BJ_facs_file_path[] = { "1500_stations_xy.txt",
										"1500_logistics_xy.txt",
										"1500_parking_xy.txt", 
										"1500_cafes_xy.txt" };
	std::string BJ_facs[] = { "_st", "_lo", "_pa", "_ca" };

	for (int iref = 0; iref < 1; ++iref) // only 20k
	{
		std::vector<geo_point> users;
		md_load_realusers((BJ_file_path + reflocs_file_path[iref]).c_str(), users);
		//md_load_users((BJ_file_path + reflocs_avg_file_path[iref]).c_str(),	(BJ_file_path + reflocs_file_path[iref]).c_str(), users);
		
		for (int ifac = 0; ifac < 4; ++ifac)
		{
			std::map<int, geo_point> facs;
			geo_cand_rtree facs_tree;
			boost::unordered_map<int, float> fac_neg_effect;
			md_load_facs((BJ_file_path + BJ_facs_file_path[ifac]).c_str(), facs, facs_tree, fac_neg_effect);

			for (int icand = 0; icand < 20; ++icand)
			{
				std::string results_EN_file_path(BJ_file_path);
				results_EN_file_path += BJ_facs[ifac] + reflocs[iref] + cands[icand];
				std::string results_MD_file_path(results_EN_file_path);
				results_EN_file_path += "_EN.txt";
				results_MD_file_path += "_MD.txt";

				// EN
				std::cout << results_EN_file_path << std::endl;
				std::ofstream ofs_re_EN(results_EN_file_path, std::ofstream::out | std::ofstream::trunc); // create results file for EN
				if (ofs_re_EN.fail())	// creating results file for EN fails
					return -1;
				__int64 traversal_time; // here, this variable is just to complete param list
				std::vector<facility_replacement> topk;
				EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
					(BJ_file_path + edges_file_path).c_str(),
					(BJ_file_path + BJ_facs_file_path[ifac]).c_str(),
					(BJ_file_path + reflocs_file_path[iref]).c_str(),
					(BJ_file_path + cands_file_path[icand]).c_str(),
					traversal_time,
					0, // 0 means output all FR pairs
					&topk); // for all facility-replacement pairs
				int size = topk.size();
				for (int i = 0; i < size; ++i)
					ofs_re_EN << "<" << topk[i].f_id << ", " << topk[i].c_id << ">\t" << topk[i].ED << '\n';
				ofs_re_EN.flush();

				// Min-dist
				std::cout << results_MD_file_path << std::endl;
				std::ofstream ofs_re_MD(results_MD_file_path, std::ofstream::out | std::ofstream::trunc); // create results file for Min-Dist
				if (ofs_re_MD.fail())	// creating results file for Min-dist fails
					return -1;				
				std::vector<geo_point> cands;
				geo_cand_rtree cands_tree;
				boost::unordered_map<int, float> cand_benefit;
				md_load_cands((BJ_file_path + cands_file_path[icand]).c_str(), cands, cands_tree, cand_benefit);
				std::vector<facility_replacement> FRs;
				//float total_dist = 0.0f;
				//MINDIST(FRs, total_dist, users, facs, facs_tree, cands);
				//ofs_re_MD << "total: " << total_dist << '\n';
				MINDIST_FAST(FRs, users, facs, facs_tree, fac_neg_effect, cands, cands_tree, cand_benefit);
				for (int i = 0; i < FRs.size(); ++i)
					ofs_re_MD << "<" << FRs[i].f_id << ", " << FRs[i].c_id << ">\t" << FRs[i].ED << '\n';
				ofs_re_MD.flush();
			}
		}
	}
	
	return 0;
}

int results_for_EN_MinDist(_TCHAR* argv[])
{
	typedef boost::unordered_map<std::pair<int, int>, float> en_results_map;
	float CA_RANGE = 1380.76f, BJ_RANGE = 258.859f;
	std::string CA_file_path, CA_results_file_path, BJ_file_path, BJ_results_file_path;
	std::ifstream ifs_config(argv[2]);
	if (!ifs_config.fail())	// Config file exists.
	{
		while (!ifs_config.bad() && ifs_config.good())
		{
			char buf[1024];
			ifs_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "CA")
			{
				CA_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				CA_results_file_path = CA_file_path + "results_CA\\CA_20k\\";
			}
			else if (str_param == "BJ")
			{

				BJ_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				BJ_results_file_path = BJ_file_path + "results_BJ\\BJ_20k\\";

			}
		}
	}
	std::cout << CA_results_file_path;
	std::ofstream ofs_re(CA_results_file_path + "__results.txt", std::ofstream::out | std::ofstream::trunc); // create results file
	if (ofs_re.fail())
		return -1;
	ofs_re << "\t\tf dist\t%\tc dist\t%\tED EN\tMD\tDiff\t%\n";
	ofs_re << "--------------------------------------------\n";

	std::string reflocs_avg_file_path[] = {
		"gen_avg_reflocs_20000.txt",
		"gen_avg_reflocs_60000.txt",
		"gen_avg_reflocs_100000.txt" };
	std::string reflocs_file_path[] = {
		"gen_rand_reflocs_20000.txt",
		"gen_rand_reflocs_60000.txt",
		"gen_rand_reflocs_100000.txt" };
	std::string reflocs[] = { "_20k", "_60k", "_100k" }; // only reflocs[0] has results
	std::string cands_file_path[] = {
		"gen_cands_100a.txt",
		"gen_cands_100b.txt",
		"gen_cands_100c.txt",
		"gen_cands_100d.txt",
		"gen_cands_100e.txt",
		"gen_cands_100f.txt",
		"gen_cands_100g.txt",
		"gen_cands_100h.txt",
		"gen_cands_100i.txt",
		"gen_cands_100j.txt",
		"gen_cands_100k.txt",
		"gen_cands_100l.txt",
		"gen_cands_100m.txt",
		"gen_cands_100n.txt",
		"gen_cands_100o.txt",
		"gen_cands_100p.txt",
		"gen_cands_100q.txt",
		"gen_cands_100r.txt",
		"gen_cands_100s.txt",
		"gen_cands_100t.txt" };
	std::string cands[] = { "_0", "_1", "_2", "_3", "_4", "_5", "_6", "_7", "_8", "_9","_10","_11","_12","_13","_14","_15","_16","_17","_18","_19","_20" };

	// CA
	std::string CA_facs_file_path[] = {
		"800_Gas_xy.txt",
		"800_pos_xy.txt",
		"800_ATM_xy.txt",
		"800_Sub_xy.txt" };
	std::string CA_facs[] = { "_GA", "_PO", "_AT", "_SU" };

	for (int iref = 0; iref < 1; ++iref) // only 20k
	{
		for (int ifac = 0; ifac < 4; ++ifac)
		{
			std::map<int, geo_point> facs;
			md_load_facs((CA_file_path + CA_facs_file_path[ifac]).c_str(), facs);

			for (int icand = 0; icand < 20; ++icand)
			{
				std::map<int, geo_point> candidates;
				md_load_cands((CA_file_path + cands_file_path[icand]).c_str(), candidates);

				std::string results_EN_file_path(CA_results_file_path);
				results_EN_file_path += CA_facs[ifac] + reflocs[iref] + cands[icand];
				std::string results_MD_file_path(results_EN_file_path);
				results_EN_file_path += "_EN.txt";
				results_MD_file_path += "_MD.txt";

				en_results_map en_results;
				std::pair<int, int> en_top_result;
				md_load_en_results(results_EN_file_path.c_str(), en_results, en_top_result);
				std::pair<int, int> md_result;
				md_load_md_result(results_MD_file_path.c_str(), md_result);

				float dist_f = boost::geometry::distance(facs[en_top_result.first], facs[md_result.first]) * EARTH_RADIUS;
				float dist_f_p = dist_f / CA_RANGE;
				float dist_c = boost::geometry::distance(candidates[en_top_result.second], candidates[md_result.second]) * EARTH_RADIUS;
				float dist_c_p = dist_c / CA_RANGE;
				float ED_EN = en_results[en_top_result];
				float ED_MD = en_results[md_result];
				float diff = ED_EN - ED_MD;
				float diff_p = diff / ED_EN;
				ofs_re << CA_facs[ifac] + reflocs[iref] + cands[icand] << "\t" // left column
					<< dist_f << "\t"
					<< dist_f_p << "\t"
					<< dist_c << "\t"
					<< dist_c_p << "\t"
					<< ED_EN << "\t"
					<< ED_MD << "\t"
					<< diff << "\t"
					<< diff_p << std::endl;
			}
		}

	}


	ofs_re << "--------------------------------------------\n";

	// BJ
	std::string BJ_facs_file_path[] = {
		"1500_stations_xy.txt",
		"1500_logistics_xy.txt",
		"1500_parking_xy.txt",
		"1500_cafes_xy.txt" };
	std::string BJ_facs[] = { "_st", "_lo", "_pa", "_ca" };

	for (int iref = 0; iref < 1; ++iref) // only 20k
	{
		for (int ifac = 0; ifac < 4; ++ifac)
		{
			std::map<int, geo_point> facs;
			md_load_facs((BJ_file_path + BJ_facs_file_path[ifac]).c_str(), facs);

			for (int icand = 0; icand < 20; ++icand)
			{
				std::map<int, geo_point> candidates;
				md_load_cands((BJ_file_path + cands_file_path[icand]).c_str(), candidates);

				std::string results_EN_file_path(BJ_results_file_path);
				results_EN_file_path += BJ_facs[ifac] + reflocs[iref] + cands[icand];
				std::string results_MD_file_path(results_EN_file_path);
				results_EN_file_path += "_EN.txt";
				results_MD_file_path += "_MD.txt";

				en_results_map en_results;
				std::pair<int, int> en_top_result;
				md_load_en_results(results_EN_file_path.c_str(), en_results, en_top_result);
				std::pair<int, int> md_result;
				md_load_md_result(results_MD_file_path.c_str(), md_result);

				float dist_f = boost::geometry::distance(facs[en_top_result.first], facs[md_result.first]) * EARTH_RADIUS;
				float dist_f_p = dist_f / BJ_RANGE;
				float dist_c = boost::geometry::distance(candidates[en_top_result.second], candidates[md_result.second]) * EARTH_RADIUS;
				float dist_c_p = dist_c / BJ_RANGE;
				float ED_EN = en_results[en_top_result];
				float ED_MD = en_results[md_result];
				float diff = ED_EN - ED_MD;
				float diff_p = diff / ED_EN;
				std::cout << BJ_facs[ifac] + reflocs[iref] + cands[icand] << std::endl;
				ofs_re << BJ_facs[ifac] + reflocs[iref] + cands[icand] << "\t" // left column
					<< dist_f << "\t"
					<< dist_f_p << "\t"
					<< dist_c << "\t"
					<< dist_c_p << "\t"
					<< ED_EN << "\t"
					<< ED_MD << "\t"
					<< diff << "\t"
					<< diff_p << std::endl;
			}
		}
	}

	return 0;
}
int _tmain(int argc, _TCHAR* argv[])
{
	// parameter "-md" means calculate for comparing with min-dist FR problems

	if (std::wstring(argv[1]) == _T("-md"))
	    return compare_EN_with_MinDist(argv);
	if (std::wstring(argv[1]) == _T("-mdr"))
		return results_for_EN_MinDist(argv);

	// parameter "-ap" means calculate AP for two ranked sequences
	// argv[2]: real data file path
	// argv[3]: query data file path
	// argv[4]: AP result file path
	if (std::wstring(argv[1]) == _T("-ap"))
	{
		std::ifstream in_real_file(argv[2]);
		std::map<unsigned, std::string> real_map;
		load_ranked_data(in_real_file, real_map, 0, '\t');

		std::ifstream in_query_file(argv[3]);
		std::map<unsigned, std::string> query_map;
		load_ranked_data(in_query_file, query_map, 0, '\t');

		double ap = 0.0;
		APatK(real_map, query_map, ap);

		std::ofstream ofs_ap(argv[4], std::ofstream::out | std::ofstream::trunc); // create AP result file
		if (!ofs_ap.fail())	// creating AP result file successfully
			ofs_ap << "AP@" << real_map.size() << ": " << ap << std::endl;
		else
			return -1;

		return 0;
	}

	// read config file for datasets files
	bool need_gen, need_cov, need_real;
	bool is_directed_graph = false; // default is bidirected graph
	std::string edges_file_path, facs_file_path, reflocs_file_path, cands_file_path;
	std::string results_file_path;
	bool is_EN_50 = false; // EN for top-50 candidates
	bool is_execute_EN = false, is_execute_LNB = false, is_execute_NSJ = false;
	bool is_Greedy_K = false;
	int K = 0; // default is 5, however, here is an init value to verify config file format (K cannot be 0)
	bool is_execute_RID = false;
	std::ifstream ifs_config(argv[1]);
	if (!ifs_config.fail())	// Config file exists.
	{
		while (!ifs_config.bad() && ifs_config.good())
		{
			char buf[1024];
			ifs_config.getline(buf, 1024);
			std::string str_buf(buf);

			std::string::size_type begin_pos = 0, mid_pos;
			mid_pos = str_buf.find(' ', begin_pos);

			std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
			if (str_param == "need_gen")
				need_gen = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "need_cov")
				need_cov = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "need_real")
				need_real = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			if (str_param == "is_directed_graph")
				is_directed_graph = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "edges")
				edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "facs")
				facs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "reflocs")
				reflocs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "cands")
				cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "results")
				results_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			else if (str_param == "is_EN_50")
				is_EN_50 = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_EN")
				is_execute_EN = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_LNB")
				is_execute_LNB = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_execute_NSJ")
				is_execute_NSJ = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "is_Greedy_K")
				is_Greedy_K = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
			else if (str_param == "K")
				K = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
			else if (str_param == "is_execute_RID")
				is_execute_RID = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
		}
	}
	if (need_gen)
		if (!gen_datasets(argv[2])) // generate candidates and reference locations
			return -1;

	if (need_cov)
		if (!convert_datasets(argv[3])) // convert facilities, then generate candidates and reference locations
			return -1;

	if (need_real)
		if (!real_datasets(argv[4])) // snap real reference locations to road network
			return -1;

	// output results
	std::ofstream ofs_re(results_file_path, std::ofstream::out | std::ofstream::trunc); // create results file
	if (ofs_re.fail())	// creating results file fails
		return -1;

	// time recorders
	auto begin_time = std::chrono::high_resolution_clock::now();
	auto end_time = std::chrono::high_resolution_clock::now();
	/**
	// Greedy-K LNB

	std::cout << is_Greedy_K << ", " << K << std::endl;
	std::vector<facility_replacement> FRs; // record the optimal facility-replacement pairs result
	if (is_Greedy_K)
	{
		// 1st time LNB
		std::cout << "1st time LNB" << std::endl;
		facility_replacement op_LNB;
		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_LNB;
		fac_hash_map facs_map;
		LNT_hash_map LNT;
		L2NT_hash_map L2NT;
		vertex_dnn_hash_map dnn_map;
		delta_f_hash_map delta_f_map;
		std::pair<int, float> max_delta_f;
		refloc_hash_map reflocs_map;
		refloc_LN_or_L2N_hash_map refloc_LN_map;
		refloc_LN_or_L2N_hash_map refloc_L2N_map;
		LNB_construct_LNT(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), LNT, L2NT, dnn_map, delta_f_map, max_delta_f, graph_for_LNB,
			facs_map, reflocs_map, refloc_LN_map, refloc_L2N_map);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "LNB construct LNT & L2NT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		begin_time = std::chrono::high_resolution_clock::now();
		int checked_num;
		cand_info_hash_map cands_map;
		op_LNB = LNB_query(checked_num, LNT, L2NT, dnn_map, delta_f_map, max_delta_f, graph_for_LNB, facs_map, reflocs_map, cands_file_path.c_str(), cands_map);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "LNB took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			<< "    check num = " << checked_num << '\n'
			<< "    result is <" << op_LNB.f_id << ", " << op_LNB.c_id << ">, ED = " << op_LNB.ED << std::endl;

		//for (cand_info_hash_map::iterator iter = cands_map.begin(); iter != cands_map.end(); ++iter)
		//	std::cout << iter->first << ", " << iter->second.get<C_VS>() << ", " << iter->second.get<C_VE>() << std::endl;

		FRs.push_back(op_LNB);
		//std::cout << op_LNB.f_id << ", " << op_LNB.c_id << std::endl;

		// 2nd - Kth times LNB
		for (int i = 1; i < K; ++i)
		{
			std::cout << i + 1 << "st time LNB" << std::endl;
			begin_time = std::chrono::high_resolution_clock::now();
			std::set<int> obsolete_facs;

			LNB_update_LNT_L2NT_Delta(op_LNB, LNT, L2NT, delta_f_map, max_delta_f, graph_for_LNB, facs_map, reflocs_map, refloc_LN_map, refloc_L2N_map,
				cands_map, obsolete_facs, reflocs_file_path.c_str());
			end_time = std::chrono::high_resolution_clock::now();
			ofs_re << "LNB re-construct LNT & L2NT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;
			//std::cout << delta_f_map[op_LNB.c_id + 10000].first << std::endl;
			//std::cout << max_delta_f.first << ", " << max_delta_f.second << std::endl;

			begin_time = std::chrono::high_resolution_clock::now();
			int checked_num;
			op_LNB = LNB_query(checked_num, LNT, L2NT, dnn_map, delta_f_map, max_delta_f, graph_for_LNB, facs_map, reflocs_map, NULL, cands_map);
			end_time = std::chrono::high_resolution_clock::now();
			ofs_re << "Greedy-K LNB took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
				<< "    check num = " << checked_num << '\n'
				<< "    result is <" << op_LNB.f_id << ", " << op_LNB.c_id << ">, ED = " << op_LNB.ED << std::endl;

			FRs.push_back(op_LNB);
			//std::cout << op_LNB.f_id << ", " << op_LNB.c_id << std::endl;
		}
	}
	*/
	// EN for top-100 FR pairs
	/*if (is_EN_50)
	{
		int top_num = 100;
		__int64 traversal_time;
		begin_time = std::chrono::high_resolution_clock::now();
		std::vector<facility_replacement> topk;
		EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str(), traversal_time,
			top_num, &topk); // for top-50 facility-replacement pairs
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "EN took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n";
		for (int i = 0; i < top_num; ++i)
			ofs_re << "<" << topk[i].f_id << ", " << topk[i].c_id << ">\t" << topk[i].ED << '\n';
		ofs_re << "ED difference: " << topk[0].ED - topk[top_num - 1].ED << '\n';
		ofs_re.flush();

		return 0;
	}*/

	// the optimal facility-replacement pair result
	facility_replacement op_EN, op_LNB, op_NSJ, op_RID;
	
	if (is_execute_EN)
	{
		__int64 traversal_time;
		begin_time = std::chrono::high_resolution_clock::now();
		op_EN = EN(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), cands_file_path.c_str(), traversal_time);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "Network traversals of all reference locations took " << traversal_time << " ms\n"
			<< "EN took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			<< "   result is <" << op_EN.f_id << ", " << op_EN.c_id << ">, ED = " << op_EN.ED << std::endl;
	}

	if (is_execute_LNB)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_LNB;
		fac_hash_map facs_map;
		LNT_hash_map LNT;
		L2NT_hash_map L2NT;
		vertex_dnn_hash_map dnn_map;
		delta_f_hash_map delta_f_map;
		std::pair<int, float> max_delta_f;
		refloc_hash_map reflocs_map;
		LNB_construct_LNT(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), LNT, L2NT, dnn_map, delta_f_map, max_delta_f, graph_for_LNB, facs_map, reflocs_map);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "LNB construct LNT & L2NT took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		begin_time = std::chrono::high_resolution_clock::now();
		int checked_num;
		op_LNB = LNB_query(checked_num, LNT, L2NT, dnn_map, delta_f_map, max_delta_f, graph_for_LNB, facs_map, reflocs_map, cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "LNB took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			<< "    check num = " << checked_num << '\n'
			<< "    result is <" << op_LNB.f_id << ", " << op_LNB.c_id << ">, ED = " << op_LNB.ED << std::endl;
	}

	if (is_execute_NSJ)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_NSJ;
		cand_max_fibonacci_heap fac_ordered_by_delta_f;
		delta_f_hash_map f_r_hash_map;
		std::pair<int, float> max_delta_f;
		fac_hash_map reflocs;
		nnfc_hash_map NNFCs;
		NSJ_construct_NNFC(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), f_r_hash_map, max_delta_f, reflocs, NNFCs, graph_for_NSJ);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "NSJ construct NNFCs took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		begin_time = std::chrono::high_resolution_clock::now();
		int checked_num;
		__int64 time_costs[2]; // [0]R*-tree time, [1]iterate max-heap
		op_NSJ = NSJ_query(checked_num, time_costs, f_r_hash_map, max_delta_f, reflocs, NNFCs, graph_for_NSJ, cands_file_path.c_str());
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "NSJ took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << " ms\n"
			<< "    R*-tree query took " << time_costs[0] << "ms, \n"
			<< "    iteration of candidates Max-heap took " << time_costs[1] << "ms, \n"
			<< "    check num = " << checked_num << '\n'
			<< "    result is <" << op_NSJ.f_id << ", " << op_NSJ.c_id << ">, ED = " << op_NSJ.ED << std::endl;
	}

	if (is_execute_RID)
	{
		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_RID;
		//cand_max_fibonacci_heap fac_ordered_by_delta_f;
		delta_f_hash_map f_r_hash_map;
		std::pair<int, float> max_delta_f; // not used in functions, however, is still reserved
		fac_hash_map reflocs;
		nnfc_hash_map NNFCs;
		ric_hash_map RICs;
		max_delta_fibonacci_heap delta_F;
		RID_construct_NNFC_and_RIC(!is_directed_graph, // this param indicates to convert (true) undirected graph into bidirectional graph or not (false)
			edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), f_r_hash_map, reflocs, NNFCs, RICs, delta_F, graph_for_RID);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "RID construct RICs and NNFCs took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		int checked_num;
		__int64 time_cost;
		op_RID = RID_query(checked_num, time_cost, f_r_hash_map, max_delta_f, reflocs, NNFCs, RICs, delta_F, graph_for_RID, cands_file_path.c_str());
		ofs_re << "RID took " << time_cost << "ms, \n"
			<< "    check num = " << checked_num << '\n'
			<< "    result is <" << op_RID.f_id << ", " << op_RID.c_id << ">, ED = " << op_RID.ED << std::endl;
	}
	/**
	//Greedy-k RID
	std::vector<facility_replacement> FRs;
	if (is_execute_RID)
	{
		// 1st time RID
		std::cout << "1st time RID" << std::endl;
		facility_replacement op_RID;

		begin_time = std::chrono::high_resolution_clock::now();
		typed_directed_graph graph_for_RID;
		delta_f_hash_map f_r_hash_map;//cand_max_fibonacci_heap fac_ordered_by_delta_f;
		std::pair<int, float> max_delta_f; // not used in functions, however, is still reserved
		fac_hash_map reflocs;
		nnfc_hash_map NNFCs;
		ric_hash_map RICs;
		max_delta_fibonacci_heap delta_F;
		cover_hash_map Rc_map;
		geo_cand_rtree cand_rtree;
		boost::unordered_map<int, geo_point> cand_geos;

		RID_construct_NNFC_and_RIC(!is_directed_graph, edges_file_path.c_str(), facs_file_path.c_str(), reflocs_file_path.c_str(), f_r_hash_map, reflocs, NNFCs, RICs, delta_F, graph_for_RID);
		end_time = std::chrono::high_resolution_clock::now();
		ofs_re << "RID construct RICs and NNFCs took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

		int checked_num;
		__int64 time_cost;
		op_RID = RID_query(checked_num, time_cost, f_r_hash_map, max_delta_f, Rc_map, cand_rtree, cand_geos, reflocs, NNFCs,RICs, delta_F, graph_for_RID, cands_file_path.c_str());
		ofs_re << "RID took " << time_cost << "ms, \n"
			<< "    check num = " << checked_num << '\n'
			<< "    result is <" << op_RID.f_id << ", " << op_RID.c_id << ">, ED = " << op_RID.ED << std::endl;


		FRs.push_back(op_RID);

		for (int i = 1; i < K; ++i)
		{
			std::cout << i + 1 << "st time RICs and NNFCs " << std::endl;

			begin_time = std::chrono::high_resolution_clock::now();
			RID_update_NNFC_and_RIC(graph_for_RID, f_r_hash_map,
				NNFCs, RICs, delta_F, op_RID, Rc_map, i);
			end_time = std::chrono::high_resolution_clock::now();

			ofs_re << "RID re-construct RISs & NNFCs took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time).count() << std::endl;

			op_RID = RID_new_query(checked_num, time_cost, f_r_hash_map, max_delta_f, Rc_map,
				reflocs, NNFCs, RICs, delta_F, graph_for_RID, op_RID, cand_rtree, cand_geos);
			ofs_re << "RID took " << time_cost << " ms\n"
				<< "    check num = " << checked_num << '\n'
				<< "    result is <" << op_RID.f_id << ", " << op_RID.c_id << ">, ED = " << op_RID.ED << std::endl;

			FRs.push_back(op_RID);
		}
	}
	*/
	return 0;
}
	bool gen_datasets(_TCHAR* argv)
	{
		// read gen_config file for generating datasets files
		std::string geo_range_file_path;
		bool need_integrate_facs = false;
		bool fac_with_coordinates = false;
		int bar, hospital, po, park, school;
		std::string bars_file_path, hospitals_file_path, pos_file_path, parks_file_path, schools_file_path;
		std::string bars_200_file_path, hospitals_800_file_path, pos_800_file_path, parks_3200_file_path, schools_3200_file_path;
		bool need_gen_cands = false;
		std::string str_cands_count_1, str_cands_count_2, str_cands_count_3, str_cands_count_4, str_cands_count_5;
		int cands_count_1 = 0, cands_count_2 = 0, cands_count_3 = 0, cands_count_4 = 0, cands_count_5 = 0;
		std::string gen_cands_file_path;
		bool need_gen_ref_loc = false;
		std::string str_clients_count_1, str_clients_count_2, str_clients_count_3, str_clients_count_4, str_clients_count_5;
		int clients_count_1 = 0, clients_count_2 = 0, clients_count_3 = 0, clients_count_4 = 0, clients_count_5 = 0;
		std::string ref_loc_points, gen_random_ref_locs_file_path, gen_avg_ref_locs_file_path;
		bool need_check_directionality = false;
		bool need_gen_edges = false;
		std::string gen_edges_file_path;
		std::string vertices_file_path, edges_file_path, pois_file_path;

		std::ifstream ifs_gen_config(argv);
		if (!ifs_gen_config.fail())	// gen_config file exists.
		{
			while (!ifs_gen_config.bad() && ifs_gen_config.good())
			{
				char buf[1024];
				ifs_gen_config.getline(buf, 1024);
				std::string str_buf(buf);

				std::string::size_type begin_pos = 0, mid_pos;
				mid_pos = str_buf.find(' ', begin_pos);

				std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
				if (str_param == "geo_range")
					geo_range_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for facilities
				else if (str_param == "need_integrate_facs")
					need_integrate_facs = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "fac_with_coordinates")
					fac_with_coordinates = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "bar" && need_integrate_facs)
					bar = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "bars" && need_integrate_facs)
					bars_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "bars_200" && need_integrate_facs)
					bars_200_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "hospital" && need_integrate_facs)
					hospital = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "hospitals" && need_integrate_facs)
					hospitals_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "hospitals_800" && need_integrate_facs)
					hospitals_800_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "po" && need_integrate_facs)
					po = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "pos" && need_integrate_facs)
					pos_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "pos_800" && need_integrate_facs)
					pos_800_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "park" && need_integrate_facs)
					park = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "parks" && need_integrate_facs)
					parks_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "parks_3200" && need_integrate_facs)
					parks_3200_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "school" && need_integrate_facs)
					school = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "schools" && need_integrate_facs)
					schools_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "schools_3200" && need_integrate_facs)
					schools_3200_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for candidates
				else if (str_param == "need_gen_cands")
					need_gen_cands = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "cands_count_1" && need_gen_cands)
				{
					str_cands_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_1 = atoi(str_cands_count_1.c_str());
				}
				else if (str_param == "cands_count_2" && need_gen_cands)
				{
					str_cands_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_2 = atoi(str_cands_count_2.c_str());
				}
				else if (str_param == "cands_count_3" && need_gen_cands)
				{
					str_cands_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_3 = atoi(str_cands_count_3.c_str());
				}
				else if (str_param == "cands_count_4" && need_gen_cands)
				{
					str_cands_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_4 = atoi(str_cands_count_4.c_str());
				}
				else if (str_param == "cands_count_5" && need_gen_cands)
				{
					str_cands_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_5 = atoi(str_cands_count_5.c_str());
				}
				else if (str_param == "gen_cands" && need_gen_cands)
					gen_cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for reference locations (clients)
				else if (str_param == "need_gen_ref_loc")
					need_gen_ref_loc = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "clients_count_1" && need_gen_ref_loc)
				{
					str_clients_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_1 = atoi(str_clients_count_1.c_str());
				}
				else if (str_param == "clients_count_2" && need_gen_ref_loc)
				{
					str_clients_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_2 = atoi(str_clients_count_2.c_str());
				}
				else if (str_param == "clients_count_3" && need_gen_ref_loc)
				{
					str_clients_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_3 = atoi(str_clients_count_3.c_str());
				}
				else if (str_param == "clients_count_4" && need_gen_ref_loc)
				{
					str_clients_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_4 = atoi(str_clients_count_4.c_str());
				}
				else if (str_param == "clients_count_5" && need_gen_ref_loc)
				{
					str_clients_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_5 = atoi(str_clients_count_5.c_str());
				}
				else if (str_param == "ref_loc_points" && need_gen_ref_loc)
					ref_loc_points = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "gen_random_ref_locs" && need_gen_ref_loc)
					gen_random_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "gen_avg_ref_locs" && need_gen_ref_loc)
					gen_avg_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for edges
				else if (str_param == "need_gen_edges")
					need_gen_edges = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "gen_edges")
					gen_edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for background datasets
				else if (str_param == "need_check_directionality")
					need_check_directionality = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "vertices")
					vertices_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "edges")
					edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "pois")
					pois_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			}
		}
		else
			return false;

		if (need_integrate_facs || need_gen_cands || need_gen_ref_loc || need_gen_edges)
		{
			raw_vertices vertices;
			raw_edges edges;
			geo_edge_rtree rtree;
			if (!load_vertices_and_edges(vertices, edges, rtree, geo_range_file_path.c_str(), vertices_file_path.c_str(), edges_file_path.c_str(), need_check_directionality))
				return false;

			boost::unordered_set<fac_or_cand_loc> gen_facs; // locations of generated facilities
			if (need_integrate_facs)
			{
				// initialize counts array
				int pois_count[CATEGORY];
				for (int i = 0; i < CATEGORY; ++i)
					pois_count[i] = 0;

				// load raw pois
				std::vector<geo_point> bars, hospitals, pos, parks, schools; // raw poi containers
				if (!load_raw_pois(pois_file_path.c_str(), pois_count, bars, hospitals, pos, parks, schools))
					return false;

				// integrate pois
				/*if (!integrate_pois(edges, rtree, bars, bars_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, hospitals, hospitals_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, pos, pos_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, parks, parks_file_path.c_str(), gen_facs)
				|| !integrate_pois(edges, rtree, schools, schools_file_path.c_str(), gen_facs))*/
				if (!integrate_pois(edges, rtree, hospitals, hospitals_file_path.c_str(), gen_facs, fac_with_coordinates)
					|| !integrate_pois(edges, rtree, parks, parks_file_path.c_str(), gen_facs, fac_with_coordinates)
					|| !integrate_pois(edges, rtree, schools, schools_file_path.c_str(), gen_facs, fac_with_coordinates))
					return false;

				// pick pois
				/*if (!pick_pois(bars_file_path.c_str(), bars_200_file_path.c_str(), bar)
				|| !pick_pois(hospitals_file_path.c_str(), hospitals_800_file_path.c_str(), hospital)
				|| !pick_pois(pos_file_path.c_str(), pos_800_file_path.c_str(), po)
				|| !pick_pois(parks_file_path.c_str(), parks_3200_file_path.c_str(), park)
				|| !pick_pois(schools_file_path.c_str(), schools_3200_file_path.c_str(), school))*/
				if (!pick_pois(hospitals_file_path.c_str(), hospitals_800_file_path.c_str(), hospital)
					|| !pick_pois(parks_file_path.c_str(), parks_3200_file_path.c_str(), park)
					|| !pick_pois(schools_file_path.c_str(), schools_3200_file_path.c_str(), school))
					return false;
			}

			// generate candidates
			if (need_gen_cands)
			{
				if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_1 + ".txt").c_str(), cands_count_1))
					return false;
				if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_2 + ".txt").c_str(), cands_count_2))
					return false;
				if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_3 + ".txt").c_str(), cands_count_3))
					return false;
				if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_4 + ".txt").c_str(), cands_count_4))
					return false;
				if (!gen_cands(vertices, edges, gen_facs, std::string(gen_cands_file_path + str_cands_count_5 + ".txt").c_str(), cands_count_5))
					return false;
			}

			// generate reference locations
			if (need_gen_ref_loc)
			{
				if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(), clients_count_1, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(), clients_count_2, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(), clients_count_3, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(), clients_count_4, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(vertices, edges, std::string(gen_random_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(), clients_count_5, ref_loc_points.c_str()))
					return false;
			}

			// just output edges, whose distances are already replaced with Euclidean distance between vertices
			if (need_gen_edges
				&& !gen_edges(edges, gen_edges_file_path.c_str()))
				return false;
		}

		return true;
	}

	// convert and generate datasets based on Beijing datasets
	bool convert_datasets(_TCHAR* argv)
	{
		// read cov_config file for generating datasets files
		std::string cov_range_file_path;
		bool need_integrate_facs = false;
		bool fac_with_coordinates = false;
		int station, cafe, logistic, bank, school;
		std::string stations_file_path, cafes_file_path, logistics_file_path, banks_file_path, schools_file_path;
		std::string stations_1500_file_path, cafes_1500_file_path, logistics_2500_file_path, banks_2500_file_path, schools_2500_file_path;
		bool need_gen_cands = false;
		std::string str_cands_count_1, str_cands_count_2, str_cands_count_3, str_cands_count_4, str_cands_count_5;
		int cands_count_1 = 0, cands_count_2 = 0, cands_count_3 = 0, cands_count_4 = 0, cands_count_5 = 0;
		std::string gen_cands_file_path;
		bool need_gen_ref_loc = false;
		std::string str_clients_count_1, str_clients_count_2, str_clients_count_3, str_clients_count_4, str_clients_count_5;
		int clients_count_1 = 0, clients_count_2 = 0, clients_count_3 = 0, clients_count_4 = 0, clients_count_5 = 0;
		std::string ref_loc_points, gen_random_ref_locs_file_path, gen_avg_ref_locs_file_path;
		bool need_gen_edges = false;
		std::string gen_edges_file_path;
		std::string edges_file_path, geos_file_path;

		std::ifstream ifs_cov_config(argv);
		if (!ifs_cov_config.fail())	// cov_config file exists.
		{
			while (!ifs_cov_config.bad() && ifs_cov_config.good())
			{
				char buf[1024];
				ifs_cov_config.getline(buf, 1024);
				std::string str_buf(buf);

				std::string::size_type begin_pos = 0, mid_pos;
				mid_pos = str_buf.find(' ', begin_pos);

				std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
				if (str_param == "cov_range")
					cov_range_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for facilities
				else if (str_param == "need_integrate_facs")
					need_integrate_facs = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "fac_with_coordinates")
					fac_with_coordinates = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "station" && need_integrate_facs)
					station = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "stations" && need_integrate_facs)
					stations_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "stations_1500" && need_integrate_facs)
					stations_1500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "cafe" && need_integrate_facs)
					cafe = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "cafes" && need_integrate_facs)
					cafes_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "cafes_1500" && need_integrate_facs)
					cafes_1500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "logistic" && need_integrate_facs)
					logistic = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "logistics" && need_integrate_facs)
					logistics_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "logistics_2500" && need_integrate_facs)
					logistics_2500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "bank" && need_integrate_facs)
					bank = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "banks" && need_integrate_facs)
					banks_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "banks_2500" && need_integrate_facs)
					banks_2500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "school" && need_integrate_facs)
					school = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "schools" && need_integrate_facs)
					schools_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "schools_2500" && need_integrate_facs)
					schools_2500_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for candidates
				else if (str_param == "need_gen_cands")
					need_gen_cands = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "cands_count_1" && need_gen_cands)
				{
					str_cands_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_1 = atoi(str_cands_count_1.c_str());
				}
				else if (str_param == "cands_count_2" && need_gen_cands)
				{
					str_cands_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_2 = atoi(str_cands_count_2.c_str());
				}
				else if (str_param == "cands_count_3" && need_gen_cands)
				{
					str_cands_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_3 = atoi(str_cands_count_3.c_str());
				}
				else if (str_param == "cands_count_4" && need_gen_cands)
				{
					str_cands_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_4 = atoi(str_cands_count_4.c_str());
				}
				else if (str_param == "cands_count_5" && need_gen_cands)
				{
					str_cands_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					cands_count_5 = atoi(str_cands_count_5.c_str());
				}
				else if (str_param == "gen_cands" && need_gen_cands)
					gen_cands_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for reference locations (clients)
				else if (str_param == "need_gen_ref_loc")
					need_gen_ref_loc = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "clients_count_1" && need_gen_ref_loc)
				{
					str_clients_count_1 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_1 = atoi(str_clients_count_1.c_str());
				}
				else if (str_param == "clients_count_2" && need_gen_ref_loc)
				{
					str_clients_count_2 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_2 = atoi(str_clients_count_2.c_str());
				}
				else if (str_param == "clients_count_3" && need_gen_ref_loc)
				{
					str_clients_count_3 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_3 = atoi(str_clients_count_3.c_str());
				}
				else if (str_param == "clients_count_4" && need_gen_ref_loc)
				{
					str_clients_count_4 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_4 = atoi(str_clients_count_4.c_str());
				}
				else if (str_param == "clients_count_5" && need_gen_ref_loc)
				{
					str_clients_count_5 = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
					clients_count_5 = atoi(str_clients_count_5.c_str());
				}
				else if (str_param == "ref_loc_points" && need_gen_ref_loc)
					ref_loc_points = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "gen_random_ref_locs" && need_gen_ref_loc)
					gen_random_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "gen_avg_ref_locs" && need_gen_ref_loc)
					gen_avg_ref_locs_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for edges
				else if (str_param == "need_gen_edges")
					need_gen_edges = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str()) == 1 ? true : false;
				else if (str_param == "gen_edges")
					gen_edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for background datasets
				else if (str_param == "edges")
					edges_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "geos")
					geos_file_path = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			}
		}
		else
			return false;

		if (need_integrate_facs || need_gen_cands || need_gen_ref_loc || need_gen_edges)
		{
			raw_edges_ext edges;
			seg_of_edge subsegments; // sub-segments of edges
			geo_subsegment_rtree rtree; // sub-segments with edge id and sub-segment id
			if (!load_edges_and_geos(edges, subsegments, rtree, cov_range_file_path.c_str(), edges_file_path.c_str(), geos_file_path.c_str()))
				return false;

			// just output edges, whose distances are already accumulated by sub-segments
			if (need_gen_edges
				&& !gen_edges(edges, gen_edges_file_path.c_str()))
				return false;

			// load and pick pois
			boost::unordered_set<fac_or_cand_loc> picked_facs; // locations of picked facilities
			if (need_integrate_facs)
			{
				//first line zhushidiao
				boost::unordered_set<fac_or_cand_loc> picked_stations; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, stations_file_path.c_str(), stations_1500_file_path.c_str(), station, picked_stations, fac_with_coordinates))
					return false;
				boost::unordered_set<fac_or_cand_loc> picked_cafes; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, cafes_file_path.c_str(), cafes_1500_file_path.c_str(), cafe, picked_cafes, fac_with_coordinates))
					return false;
				boost::unordered_set<fac_or_cand_loc> picked_logistics; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, logistics_file_path.c_str(), logistics_2500_file_path.c_str(), logistic, picked_logistics, fac_with_coordinates))
					return false;
				boost::unordered_set<fac_or_cand_loc> picked_banks; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, banks_file_path.c_str(), banks_2500_file_path.c_str(), bank, picked_banks, fac_with_coordinates))
					return false;


				/*if (!load_and_pick_pois(edges, subsegments, rtree, stations_file_path.c_str(), stations_1500_file_path.c_str(), station, picked_stations))
				return false;

				boost::unordered_set<fac_or_cand_loc> picked_cafes; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, cafes_file_path.c_str(), cafes_1500_file_path.c_str(), cafe, picked_cafes))
				return false;

				boost::unordered_set<fac_or_cand_loc> picked_logistics; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, logistics_file_path.c_str(), logistics_2500_file_path.c_str(), logistic, picked_logistics))
				return false;

				boost::unordered_set<fac_or_cand_loc> picked_banks; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, banks_file_path.c_str(), banks_2500_file_path.c_str(), bank, picked_banks))
				return false;

				boost::unordered_set<fac_or_cand_loc> picked_schools; // locations of picked stations
				if (!load_and_pick_pois(edges, subsegments, rtree, schools_file_path.c_str(), schools_2500_file_path.c_str(), school, picked_schools))
				return false;

				// construct for all types of facilities
				picked_facs.insert(picked_stations.begin(), picked_stations.end());
				picked_facs.insert(picked_cafes.begin(), picked_cafes.end());
				picked_facs.insert(picked_logistics.begin(), picked_logistics.end());
				picked_facs.insert(picked_banks.begin(), picked_banks.end());
				picked_facs.insert(picked_schools.begin(), picked_schools.end()); */

			}

			// generate candidates
			if (need_gen_cands)
			{
				if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_1 + ".txt").c_str(), cands_count_1))
					return false;
				if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_2 + ".txt").c_str(), cands_count_2))
					return false;
				if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_3 + ".txt").c_str(), cands_count_3))
					return false;
				if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_4 + ".txt").c_str(), cands_count_4))
					return false;
				if (!gen_cands(edges, subsegments, picked_facs, std::string(gen_cands_file_path + str_cands_count_5 + ".txt").c_str(), cands_count_5))
					return false;
			}

			// generate reference locations
			if (need_gen_ref_loc)
			{
				if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_1 + ".txt").c_str(), clients_count_1, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_2 + ".txt").c_str(), clients_count_2, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_3 + ".txt").c_str(), clients_count_3, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_4 + ".txt").c_str(), clients_count_4, ref_loc_points.c_str()))
					return false;
				if (!gen_ref_locs(edges, subsegments, std::string(gen_random_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(),
					std::string(gen_avg_ref_locs_file_path + str_clients_count_5 + ".txt").c_str(), clients_count_5, ref_loc_points.c_str()))
					return false;
			}
		}

		return true;
	}

	// integrate real BJ and CA datasets
	bool real_datasets(_TCHAR* argv)
	{
		// read real_config file for integrate datasets files
		int BJ_user_count, CA_user_count;
		std::string BJ_edges, BJ_geos, BJ_refloc_pre, BJ_reflocs;
		std::string CA_edges, CA_vertices, CA_refloc_pre, CA_reflocs;

		std::ifstream ifs_real_config(argv);
		if (!ifs_real_config.fail())	// real_config file exists.
		{
			while (!ifs_real_config.bad() && ifs_real_config.good())
			{
				char buf[1024];
				ifs_real_config.getline(buf, 1024);
				std::string str_buf(buf);

				std::string::size_type begin_pos = 0, mid_pos;
				mid_pos = str_buf.find(' ', begin_pos);

				std::string str_param = str_buf.substr(begin_pos, mid_pos - begin_pos);
				//--------------------------------------------------------------------------------
				// for BJ
				if (str_param == "BJ_user_count")
					BJ_user_count = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "BJ_edges")
					BJ_edges = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "BJ_geos")
					BJ_geos = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "BJ_refloc_pre")
					BJ_refloc_pre = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "BJ_reflocs")
					BJ_reflocs = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				//--------------------------------------------------------------------------------
				// for CA
				else if (str_param == "CA_user_count")
					CA_user_count = atoi(str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1).c_str());
				else if (str_param == "CA_edges")
					CA_edges = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "CA_vertices")
					CA_vertices = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "CA_refloc_pre")
					CA_refloc_pre = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
				else if (str_param == "CA_reflocs")
					CA_reflocs = str_buf.substr(mid_pos + 1, str_buf.size() - mid_pos - 1);
			}
		}
		else
			return false;

		//--------------------------------------------------------------------------------
		// for CA real dataset
		raw_vertices vertices_CA;
		raw_edges edges_CA;
		geo_edge_rtree rtree_CA;
		std::cout << "load_vertices_and_edges is running ..." << std::endl;
		if (!load_vertices_and_edges(vertices_CA, edges_CA, rtree_CA, NULL, CA_vertices.c_str(), CA_edges.c_str(), false)) // NULL means geo-range file is ignored
			return false;
		if (!load_and_snap_reflocs_CA(edges_CA, rtree_CA, CA_user_count, CA_refloc_pre.c_str(), CA_reflocs.c_str()))
			return false;

		//--------------------------------------------------------------------------------
		// for BJ real dataset
		raw_edges_ext edges_BJ;
		seg_of_edge subsegments_BJ; // sub-segments of edges
		geo_subsegment_rtree rtree_BJ; // sub-segments with edge id and sub-segment id
		std::cout << "load_vertices_and_edges is running ..." << std::endl;
		if (!load_edges_and_geos(edges_BJ, subsegments_BJ, rtree_BJ, NULL, BJ_edges.c_str(), BJ_geos.c_str())) // NULL means geo-range file is ignored
			return false;
		if (!load_and_snap_reflocs_BJ(edges_BJ, subsegments_BJ, rtree_BJ, BJ_user_count, BJ_refloc_pre.c_str(), BJ_reflocs.c_str()))
			return false;

		return true;
	}


	void append_candidates_with_reflocs(const char *reflocs_file_path, const char *cands_file_path)
	{
		// create a new filename for candidates
		std::string of_path(cands_file_path);
		size_t pos = of_path.rfind('_');
		of_path.erase(pos, of_path.length() - pos);
		of_path.append("_10000.txt");

		std::ofstream ofs_cand(of_path, std::ofstream::out | std::ofstream::trunc); // open a new file to output
		if (!ofs_cand.fail())	// creating file successfully
		{
			// deal with each candidate
			std::ifstream ifs_cand(cands_file_path); // read candidates file
			if (!ifs_cand.fail()) // candidates file exists
			{
				while (!ifs_cand.bad() && ifs_cand.good())
				{
					char buf[1024];
					ifs_cand.getline(buf, 1024);
					std::string str_buf(buf);
					ofs_cand << str_buf << '\n';
				}
			}
		}

		// insert reference locations as vertices into graph
		std::ifstream ifs_reflocs(reflocs_file_path); // read reference locations file
		if (!ifs_reflocs.fail())	// reference locations file exists
		{
			while (!ifs_reflocs.bad() && ifs_reflocs.good())
			{
				char buf[1024];
				ifs_reflocs.getline(buf, 1024);
				std::string str_buf(buf);

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
				//float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
				float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
				float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

				// construct data structure for reference locations
				if (refloc_id > 9999)
					break;
				if (refloc_id > 799)
					ofs_cand << refloc_id << ' ' << vs_id << ' ' << ve_id << ' ' << dist_vs << ' ' << dist_ve << ' ' << lon << ' ' << lat << '\n';
			}
		}
	}

