#pragma once

#include <iostream>

#include "solutions.h"


void load_users_facs_cands(const char *reflocs_avg_file_path, const char *reflocs_file_path, const char *facs_file_path, const char *cands_file_path,
						   std::vector<geo_point>& users, std::map<int, geo_point>& facs, geo_cand_rtree& facs_tree, std::vector<geo_point>& cands)
{
	// load refloc_avg
	std::vector<int> reflocs_count_per_user;
	std::ifstream ifs_reflocs_avg(reflocs_avg_file_path); // read reference locations file
	if (!ifs_reflocs_avg.fail()) // reference locations file exists
	{
		int count = 0;
		while (!ifs_reflocs_avg.bad() && ifs_reflocs_avg.good())
		{
			char buf[1024];
			ifs_reflocs_avg.getline(buf, 1024);
			if (count > 0)
			{
				--count;
				continue;
			}
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);

			std::string prob = str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1);
			if (prob == "1")
			{
				reflocs_count_per_user.push_back(1);
				count = 1;
			}
			else if (prob == "0.5")
			{
				reflocs_count_per_user.push_back(2);
				count = 2;
			}
			else if (prob == "0.333333")
			{
				reflocs_count_per_user.push_back(3);
				count = 3;
			}
			else if (prob == "0.25")
			{
				reflocs_count_per_user.push_back(4);
				count = 4;
			}
			else if (prob == "0.2")
			{
				reflocs_count_per_user.push_back(5);
				count = 5;
			}
			else if (prob == "0.166667")
			{
				reflocs_count_per_user.push_back(6);
				count = 6;
			}
			--count;
		}
	}

	// load reflos
	std::vector<boost::tuple<float, float, float>> reflocs;
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

			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			reflocs.push_back(boost::make_tuple(prob, lon, lat));
		}
	}

	// identify users
	users.clear(); // reset
	size_t i_of_tuples = 0; // global index of reflocs
	for (size_t i = 0; i < reflocs_count_per_user.size(); ++i)
	{
		boost::tuple<float, float, float> max_tuple(FLT_MIN, 0.0f, 0.0f);
		for (int j = 0; j < reflocs_count_per_user[i]; ++j)
		{
			boost::tuple<float, float, float> the_tuple = reflocs[i_of_tuples];
			++i_of_tuples;

			if (the_tuple.get<0>() > max_tuple.get<0>())
				max_tuple = the_tuple;
		}
		users.push_back(geo_point(max_tuple.get<1>(), max_tuple.get<2>()));
	}

	// load facs
	facs.clear(); // reset
	facs_tree.clear(); // reset fac tree
	std::ifstream ifs_facs(facs_file_path); // read facilities file
	if (!ifs_facs.fail())	// facilities file exists
	{
		while (!ifs_facs.bad() && ifs_facs.good())
		{
			char buf[1024];
			ifs_facs.getline(buf, 1024);
			std::string str_buf(buf);

			// the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int fac_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct facility R*-tree and map
			facs[fac_id] = geo_point(lon, lat);
			facs_tree.insert(std::make_pair(geo_point(lon, lat), fac_id));
		}
	}

	// load cands
	cands.clear(); // reset
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail())	// candidates file exists
	{
		while (!ifs_cands.bad() && ifs_cands.good())
		{
			char buf[1024];
			ifs_cands.getline(buf, 1024);
			std::string str_buf(buf);

			// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct candidates vector
			cands.push_back(geo_point(lon, lat));
		}
	}
}


void md_load_users(const char *reflocs_avg_file_path, const char *reflocs_file_path,	std::vector<geo_point>& users)
{
	// load refloc_avg
	std::vector<int> reflocs_count_per_user;
	std::ifstream ifs_reflocs_avg(reflocs_avg_file_path); // read reference locations file
	if (!ifs_reflocs_avg.fail()) // reference locations file exists
	{
		int count = 0;
		while (!ifs_reflocs_avg.bad() && ifs_reflocs_avg.good())
		{
			char buf[1024];
			ifs_reflocs_avg.getline(buf, 1024);
			if (count > 0)
			{
				--count;
				continue;
			}
			std::string str_buf(buf);

			// the reference locations file format must be "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, prob_pos, lon_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			prob_pos = str_buf.find(' ', dist_ve_pos + 1);
			lon_pos = str_buf.find(' ', prob_pos + 1);

			std::string prob = str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1);
			if (prob == "1")
			{
				reflocs_count_per_user.push_back(1);
				count = 1;
			}
			else if (prob == "0.5")
			{
				reflocs_count_per_user.push_back(2);
				count = 2;
			}
			else if (prob == "0.333333")
			{
				reflocs_count_per_user.push_back(3);
				count = 3;
			}
			else if (prob == "0.25")
			{
				reflocs_count_per_user.push_back(4);
				count = 4;
			}
			else if (prob == "0.2")
			{
				reflocs_count_per_user.push_back(5);
				count = 5;
			}
			else if (prob == "0.166667")
			{
				reflocs_count_per_user.push_back(6);
				count = 6;
			}
			--count;
		}
	}

	// load reflos
	std::vector<boost::tuple<float, float, float>> reflocs;
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

			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			reflocs.push_back(boost::make_tuple(prob, lon, lat));
		}
	}

	// identify users
	users.clear(); // reset
	size_t i_of_tuples = 0; // global index of reflocs
	for (size_t i = 0; i < reflocs_count_per_user.size(); ++i)
	{
		boost::tuple<float, float, float> max_tuple(FLT_MIN, 0.0f, 0.0f);
		for (int j = 0; j < reflocs_count_per_user[i]; ++j)
		{
			boost::tuple<float, float, float> the_tuple = reflocs[i_of_tuples];
			++i_of_tuples;

			if (the_tuple.get<0>() > max_tuple.get<0>())
				max_tuple = the_tuple;
		}
		users.push_back(geo_point(max_tuple.get<1>(), max_tuple.get<2>()));
	}
}

void load_realusers_facs_cands(const char *reflocs_file_path, const char *facs_file_path, const char *cands_file_path,
	std::vector<geo_point>& users, std::map<int, geo_point>& facs, geo_cand_rtree& facs_tree, std::vector<geo_point>& cands)
{

	std::vector<int> reflocs_count_per_user;
	std::vector<boost::tuple<float, float, float>> reflocs;
	std::ifstream ifs_reflocs(reflocs_file_path);
	if (!ifs_reflocs.fail()) // reference locations file exists
	{

		int count = 0;
		float prob_sum = 0;
		float p_threshold = 0.001f;
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

			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			reflocs.push_back(boost::make_tuple(prob, lon, lat));
			count++;
			prob_sum = prob_sum + prob;
			if (abs(prob_sum - 1.0f)< p_threshold)
			{
				reflocs_count_per_user.push_back(count);
				count = 0;
				prob_sum = 0;
			}


		}
	}


	// identify users
	users.clear(); // reset
	size_t i_of_tuples = 0; // global index of reflocs
	for (size_t i = 0; i < reflocs_count_per_user.size(); ++i)
	{
		boost::tuple<float, float, float> max_tuple(FLT_MIN, 0.0f, 0.0f);
		for (int j = 0; j < reflocs_count_per_user[i]; ++j)
		{
			boost::tuple<float, float, float> the_tuple = reflocs[i_of_tuples];
			++i_of_tuples;

			if (the_tuple.get<0>() > max_tuple.get<0>())
				max_tuple = the_tuple;
		}
		users.push_back(geo_point(max_tuple.get<1>(), max_tuple.get<2>()));
	}


	facs.clear(); // reset
	facs_tree.clear(); // reset fac tree
	std::ifstream ifs_facs(facs_file_path); // read facilities file
	if (!ifs_facs.fail())	// facilities file exists
	{
		while (!ifs_facs.bad() && ifs_facs.good())
		{
			char buf[1024];
			ifs_facs.getline(buf, 1024);
			std::string str_buf(buf);

			// the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int fac_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct facility R*-tree and map
			facs[fac_id] = geo_point(lon, lat);
			facs_tree.insert(std::make_pair(geo_point(lon, lat), fac_id));
		}
	}

	// load cands
	cands.clear(); // reset
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail())	// candidates file exists
	{
		while (!ifs_cands.bad() && ifs_cands.good())
		{
			char buf[1024];
			ifs_cands.getline(buf, 1024);
			std::string str_buf(buf);

			// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct candidates vector
			cands.push_back(geo_point(lon, lat));
		}
	}
}

void md_load_realusers(const char *reflocs_file_path, std::vector<geo_point>& users)
{
	// load refloc_avg
	std::vector<int> reflocs_count_per_user;
	std::vector<boost::tuple<float, float, float>> reflocs;
	std::ifstream ifs_reflocs(reflocs_file_path);
	if (!ifs_reflocs.fail()) // reference locations file exists
	{

		int count = 0;
		float prob_sum = 0;
		float p_threshold = 0.001f;
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

			float prob = static_cast<float>(atof(str_buf.substr(prob_pos + 1, lon_pos - prob_pos - 1).c_str()));
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			reflocs.push_back(boost::make_tuple(prob, lon, lat));
			count++;
			prob_sum = prob_sum + prob;
			if (abs(prob_sum - 1.0f)< p_threshold)
			{
				reflocs_count_per_user.push_back(count);
				count = 0;
				prob_sum = 0;
			}


		}
	}


	// identify users
	users.clear(); // reset
	size_t i_of_tuples = 0; // global index of reflocs
							//std::ofstream out_file;
							//out_file.open("pro_max.txt", std::ofstream::app);
	for (size_t i = 0; i < reflocs_count_per_user.size(); ++i)
	{
		boost::tuple<float, float, float> max_tuple(FLT_MIN, 0.0f, 0.0f);
		for (int j = 0; j < reflocs_count_per_user[i]; ++j)
		{
			boost::tuple<float, float, float> the_tuple = reflocs[i_of_tuples];
			++i_of_tuples;

			if (the_tuple.get<0>() > max_tuple.get<0>())
				max_tuple = the_tuple;
		}
		users.push_back(geo_point(max_tuple.get<1>(), max_tuple.get<2>()));
		//out_file << max_tuple.get<0>() << std::endl;
	}
	std::cout << "success" << std::endl;

}

void md_load_facs(const char *facs_file_path, std::map<int, geo_point>& facs, geo_cand_rtree& facs_tree, boost::unordered_map<int, float>& fac_neg_effect)
{
	facs.clear(); // reset
	facs_tree.clear(); // reset fac tree
	fac_neg_effect.clear();
	std::ifstream ifs_facs(facs_file_path); // read facilities file
	if (!ifs_facs.fail())	// facilities file exists
	{
		while (!ifs_facs.bad() && ifs_facs.good())
		{
			char buf[1024];
			ifs_facs.getline(buf, 1024);
			std::string str_buf(buf);

			// the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int fac_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct facility R*-tree and map
			facs[fac_id] = geo_point(lon, lat);
			facs_tree.insert(std::make_pair(geo_point(lon, lat), fac_id));
			fac_neg_effect[fac_id] = 0.0f;
		}
	}
}


void md_load_facs(const char *facs_file_path, std::map<int, geo_point>& facs)
{
	facs.clear(); // reset
	std::ifstream ifs_facs(facs_file_path); // read facilities file
	if (!ifs_facs.fail())	// facilities file exists
	{
		while (!ifs_facs.bad() && ifs_facs.good())
		{
			char buf[1024];
			ifs_facs.getline(buf, 1024);
			std::string str_buf(buf);

			// the facilities file format must be "fac_id vs_id ve_id dist_vs_fac dist_fac_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int fac_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			facs[fac_id] = geo_point(lon, lat);
		}
	}
}


void md_load_cands(const char *cands_file_path, std::vector<geo_point>& cands, geo_cand_rtree& cands_tree, boost::unordered_map<int, float>& cand_benefit)
{
	cands.clear(); // reset
	cands_tree.clear(); // reset cand tree
	cand_benefit.clear();
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail())	// candidates file exists
	{
		while (!ifs_cands.bad() && ifs_cands.good())
		{
			char buf[1024];
			ifs_cands.getline(buf, 1024);
			std::string str_buf(buf);

			// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			// construct candidates vector
			cands_tree.insert(std::make_pair(geo_point(lon, lat), int(cands.size())));
			cand_benefit[int(cands.size())] = 0.0f;
			cands.push_back(geo_point(lon, lat));			
		}
	}
}


void md_load_cands(const char *cands_file_path, std::map<int, geo_point>& cands)
{
	cands.clear(); // reset
	std::ifstream ifs_cands(cands_file_path); // read candidates file
	if (!ifs_cands.fail())	// candidates file exists
	{
		while (!ifs_cands.bad() && ifs_cands.good())
		{
			char buf[1024];
			ifs_cands.getline(buf, 1024);
			std::string str_buf(buf);

			// the candidates file format must be "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat"
			std::string::size_type begin_pos = 0, vs_pos, ve_pos, dist_vs_pos, dist_ve_pos, lon_pos, lat_pos;
			vs_pos = str_buf.find(' ', begin_pos);
			ve_pos = str_buf.find(' ', vs_pos + 1);
			dist_vs_pos = str_buf.find(' ', ve_pos + 1);
			dist_ve_pos = str_buf.find(' ', dist_vs_pos + 1);
			lon_pos = str_buf.find(' ', dist_ve_pos + 1);
			lat_pos = str_buf.find(' ', lon_pos + 1);

			int cand_id = atoi(str_buf.substr(begin_pos, vs_pos - begin_pos).c_str());
			float lon = static_cast<float>(atof(str_buf.substr(lon_pos + 1, lat_pos - lon_pos - 1).c_str()));
			float lat = static_cast<float>(atof(str_buf.substr(lat_pos + 1, str_buf.size() - lat_pos - 1).c_str()));

			cands[cand_id] = geo_point(lon, lat);
		}
	}
}


void md_load_en_results(const char *results_file_path, boost::unordered_map<std::pair<int, int>, float>& en_results, std::pair<int, int>& md_result)
{
	en_results.clear(); // reset
	std::ifstream ifs_results(results_file_path); // read file
	if (!ifs_results.fail())	// file exists
	{
		bool bTop = true;
		while (!ifs_results.bad() && ifs_results.good())
		{
			char buf[1024];
			ifs_results.getline(buf, 1024);
			std::string str_buf(buf);
			if (str_buf == "")
				break;

			// the file format must be "<f, c>\tED"
			std::string::size_type f_pos = 1, comma_pos, c_pos, r_pos, ED_pos;
			comma_pos = str_buf.find(',', f_pos);
			c_pos = comma_pos + 2;
			r_pos = str_buf.find('>', c_pos);
			ED_pos = r_pos + 2;

			int f_id = atoi(str_buf.substr(f_pos, comma_pos - f_pos).c_str());
			int c_id = atoi(str_buf.substr(c_pos, r_pos - c_pos).c_str());
			float ED = static_cast<float>(atof(str_buf.substr(ED_pos, str_buf.size() - ED_pos).c_str()));

			en_results[std::make_pair(f_id, c_id)] = ED;
			if (bTop)
			{
				md_result = std::make_pair(f_id, c_id);
				bTop = false;
			}
		}
	}
}


void md_load_md_result(const char *results_file_path, std::pair<int, int>& md_result)
{
	std::ifstream ifs_results(results_file_path); // read file
	if (!ifs_results.fail())	// file exists
	{
		if (!ifs_results.bad() && ifs_results.good())
		{
			char buf[1024];
			ifs_results.getline(buf, 1024);
			std::string str_buf(buf);

			// the file format must be "<f, c>\tED"
			std::string::size_type f_pos = 1, comma_pos, c_pos, r_pos;
			comma_pos = str_buf.find(',', f_pos);
			c_pos = comma_pos + 2;
			r_pos = str_buf.find('>', c_pos);

			int f_id = atoi(str_buf.substr(f_pos, comma_pos - f_pos).c_str());
			int c_id = atoi(str_buf.substr(c_pos, r_pos - c_pos).c_str());

			md_result = std::make_pair(f_id, c_id);
		}
	}
}


void MINDIST(std::vector<facility_replacement>& FRs, float& total_dist, std::vector<geo_point>& users, std::map<int, geo_point>& facs, geo_cand_rtree& facs_tree, std::vector<geo_point>& cands)
{
	fac_rep_max_fibonacci_heap max_heap;

	// calculate total distance
	total_dist = 0.0f;
	for (size_t i = 0; i < users.size(); ++i)
	{
		std::vector<geo_cand> nearest_facility;
		facs_tree.query(boost::geometry::index::nearest(users[i], 1), std::back_inserter(nearest_facility));
		float distE = boost::geometry::distance(users[i], nearest_facility[0].first) * EARTH_RADIUS; // dE(r,c)
		total_dist += distE;
	}

	for (std::map<int, geo_point>::iterator iter_facs = facs.begin(); iter_facs != facs.end(); ++iter_facs)
	{
		std::cout << "fac: " << iter_facs->first << std::endl;

		geo_cand_rtree the_facs_tree(facs_tree);
		the_facs_tree.remove(geo_cand(iter_facs->second, iter_facs->first));
		
		for (int i = 0; i < cands.size(); ++i)
		{
			the_facs_tree.insert(geo_cand(cands[i], 10000 + i));

			float the_total_dist = 0.0f;
			for (size_t j = 0; j < users.size(); ++j)
			{
				std::vector<geo_cand> nearest_facility;
				the_facs_tree.query(boost::geometry::index::nearest(users[j], 1), std::back_inserter(nearest_facility));
				float distE = boost::geometry::distance(users[j], nearest_facility[0].first) * EARTH_RADIUS; // dE(r,c)
				the_total_dist += distE;
			}
			float diff = total_dist - the_total_dist;
			max_heap.push(facility_replacement(iter_facs->first, i, diff));

			the_facs_tree.remove(geo_cand(cands[i], 10000 + i));
		}
	}

	FRs.clear();
	while (!max_heap.empty())
	{
		FRs.push_back(max_heap.top());
		max_heap.pop(); // pop the top vertex
	}
}


void geo_offset(float lon_c, float lat_c, float distance, direction dir, float &lon_or_lat)
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


typedef boost::unordered_map < int, // user id
	boost::tuple < int, // nn facility id
	float, // dnn
	float, // d2nn
	geo_box > > // MBR of 2nfc
	nfc_hash_map; // NFC hash map
nfc_hash_map nfc_map;
enum nfc { _nn = 0, _dnn, _d2nn, _mbr };

typedef std::map < int, // cand id
	std::map < int, // fac id
	float > > // offset
	offset_map;


void MINDIST_FAST(std::vector<facility_replacement>& FRs, std::vector<geo_point>& users, std::map<int, geo_point>& facs, geo_cand_rtree& facs_tree,
	boost::unordered_map<int, float>& fac_neg_effect, std::vector<geo_point>& cands, geo_cand_rtree& cands_tree, boost::unordered_map<int, float>& cand_benefit)
{
	// nfc for each user
	for (size_t i = 0; i < users.size(); ++i)
	{
		std::vector<geo_cand> nearest_facility; // .first is geo_point, .second is fac id
		facs_tree.query(boost::geometry::index::nearest(users[i], 2), std::back_inserter(nearest_facility));
		float dnn = boost::geometry::distance(users[i], nearest_facility[0].first) * EARTH_RADIUS; // dE(r,c)
		float d2nn = boost::geometry::distance(users[i], nearest_facility[1].first) * EARTH_RADIUS; // dE(r,c)
		fac_neg_effect[nearest_facility[0].second] += dnn - d2nn;

		// offsets for nfc geo_box
		float lon = users[i].get<LON>(), lat = users[i].get<LAT>();
		float min_lon, min_lat, max_lon, max_lat; // 2nfc
		//geo_offset(lon, lat, dnn, WESTWARD, min_lon);
		//geo_offset(lon, lat, dnn, SOUTHWARD, min_lat);
		//geo_offset(lon, lat, dnn, EASTWARD, max_lon);
		//geo_offset(lon, lat, dnn, NORTHWARD, max_lat);
		geo_offset(lon, lat, d2nn, WESTWARD, min_lon);
		geo_offset(lon, lat, d2nn, SOUTHWARD, min_lat);
		geo_offset(lon, lat, d2nn, EASTWARD, max_lon);
		geo_offset(lon, lat, d2nn, NORTHWARD, max_lat);

		nfc_map[i] = boost::make_tuple(nearest_facility[0].second, // nn f id
			dnn, d2nn, geo_box(geo_point(min_lon, min_lat), geo_point(max_lon, max_lat)));
	}

	int max_neg_f = -1;
	float max_neg = -FLT_MAX;
	for (boost::unordered_map<int, float>::iterator iter = fac_neg_effect.begin(); iter != fac_neg_effect.end(); ++iter)
	{
		if (iter->second > max_neg)
		{
			max_neg_f = iter->first;
			max_neg = iter->second;
		}
	}

	// for each user, which cands can affect her
	offset_map offsets;
	for (nfc_hash_map::iterator iter_nfc = nfc_map.begin(); iter_nfc != nfc_map.end(); ++iter_nfc)
	{
		// range query for candidates which are covered by a 2nfc
		std::vector<geo_cand> covered_cands; // covered by 2nfc
		cands_tree.query(boost::geometry::index::within(iter_nfc->second.get<_mbr>()), std::back_inserter(covered_cands));
		for (std::vector<geo_cand>::iterator iter_cand = covered_cands.begin(); iter_cand != covered_cands.end(); ++iter_cand)
		{
			float distE = boost::geometry::distance(users[iter_nfc->first], cands[iter_cand->second]) * EARTH_RADIUS; // dE(r,c)
			if (distE < iter_nfc->second.get<_dnn>())
			{
				cand_benefit[iter_cand->second] += iter_nfc->second.get<_dnn>() - distE; // c positive effect

				if (offsets.find(iter_cand->second) != offsets.end()) // cand exists
				{
					if (offsets[iter_cand->second].find(iter_nfc->second.get<_nn>()) != offsets[iter_cand->second].end()) // fac exists
						offsets[iter_cand->second][iter_nfc->second.get<_nn>()] += iter_nfc->second.get<_d2nn>() - iter_nfc->second.get<_dnn>();
					else
						offsets[iter_cand->second][iter_nfc->second.get<_nn>()] = iter_nfc->second.get<_d2nn>() - iter_nfc->second.get<_dnn>();
				}
				else
				{
					offsets[iter_cand->second] = std::map<int, float>();
					offsets[iter_cand->second][iter_nfc->second.get<_nn>()] = iter_nfc->second.get<_d2nn>() - iter_nfc->second.get<_dnn>();
				}
			}
			else if (distE < iter_nfc->second.get<_d2nn>())
			{
				if (offsets.find(iter_cand->second) != offsets.end()) // cand exists
				{
					if (offsets[iter_cand->second].find(iter_nfc->second.get<_nn>()) != offsets[iter_cand->second].end()) // fac exists
						offsets[iter_cand->second][iter_nfc->second.get<_nn>()] += iter_nfc->second.get<_d2nn>() - distE;
					else
						offsets[iter_cand->second][iter_nfc->second.get<_nn>()] = iter_nfc->second.get<_d2nn>() - distE;
				}
				else
				{
					offsets[iter_cand->second] = std::map<int, float>();
					offsets[iter_cand->second][iter_nfc->second.get<_nn>()] = iter_nfc->second.get<_d2nn>() - distE;
				}
			}
		}
	}

	cand_max_fibonacci_heap cand_max_heap;
	for (boost::unordered_map<int, float>::iterator iter = cand_benefit.begin(); iter != cand_benefit.end(); ++iter)
	{
		cand_max_heap.push(candidate(iter->first, iter->second));
	}

	fac_rep_max_fibonacci_heap max_heap;
	facility_replacement ol;
	while (!cand_max_heap.empty())
	{
		candidate c = cand_max_heap.top();

		if (c.ERD <= ol.ED)
			break;
		std::cout << "cand: " << c.id << std::endl;

		float max_neg_effect = max_neg;
		int f_max_neg_effect = max_neg_f;
		for (std::map<int, float>::iterator iter = offsets[c.id].begin(); iter != offsets[c.id].end(); ++iter)
		{
			float neg = iter->second + fac_neg_effect[iter->first];
			if (neg > max_neg_effect)
			{
				f_max_neg_effect = iter->first;
				max_neg_effect = neg;
			}
		}

		ol.f_id = f_max_neg_effect;
		ol.c_id = c.id;
		ol.ED = c.ERD + max_neg_effect;

		max_heap.push(ol);

		cand_max_heap.pop();
	}

	FRs.clear();
	while (!max_heap.empty())
	{
		FRs.push_back(max_heap.top());
		max_heap.pop(); // pop the top vertex
	}
}