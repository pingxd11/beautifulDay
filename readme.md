Frost-TIST
===
Environment
---------
1. The project is implemented in C++.<br>
2. Microsoft Visual Studio Express 2015 for Windows Desktop IDE environment.<br>
3. Boost libraries need to be introduced,whose version is 1.68.0. compile x64.

Supplement
-----

`targetver.h` :show Windows platform adaption<br>
`stdafx.h`:include common libraries<br>
`stdafx.cpp`<br>

`typed_directed_graph.h`:Base class and function for initializing graph and deploying facilities<br>
`typed_directed_graph.cpp`:candidates and reference locations data<br>

`geo_defines.h`:Geo data types<br>
`greedy_k.h`:Base class and function for greedy_k LNB<br>
`gen_datasets.h`:Generate CA related data	                                                                                                                 
`calculate_vertical_point`:given a point,calculate its vertical point on a straight line,which is represented by two points<br>
`gen_cands`:output generated candidates<br>
`gen_edges`:output generated edges based on raw edges dataset,where only distances are replaced with Euclidean distance between vertices<br>
`gen_ref_locs`:Reference location data generation with distance limitation<br>
`integrate_pois`:integrate raw pois dataset into the data we need,and output results to files<br>
`load_raw_pois`:load raw pois dataset, and pick pois on demand<br>
`load_vertices_and_edges`:load vertices and edges datasets also construct an R*-tree to organize edges<br>
`pick_pois`:pick pois from integrated pois dataset,and output results to files<br>
                                  
`convert_datasets.h`: Generate CA related data<br>
`gen_cands`:output generated candidates<br>
`gen_edges`：output generated edges based on raw edges dataset,where only distances are replaced with Euclidean distance between vertices<br>
`gen_ref_locs`:output generated reference locations<br>
`load_and_pick_pois`:load and pick pois, also including integration facilities<br>
`load_edges_and_geos`:load edges and sub-segments (geos) datasets;also construct an R*-tree to organize sub-segments<br>
		
`real_datasets.h`:integrate real BJ and CA datasets	<br>
`load_and_snap_reflocs_BJ`:load and snap reference locations of BJ dataset<br>
`load_and_snap_reflocs_CA`:load and snap reference locations of CA dataset<br>

`mindist.h`: calculate for comparing with min-dist FR problems<br>
`load_users_facs_cands`:EN Load synthetic users, facilities,candidates data<br>
`load_realusers_facs_cands`:EN:Load real users,facilities,candidates data<br>
`md_load_users:min-dist`:Load synthetic users data<br>
`md_load_realusers`:min-dist:Load real users data <br>
`md_load_facs:min-dist`:Load facilities data<br>
`md_load_cands:min-dist`:Load candidates data md_load_en_results<br>                              
`EN`:Load the resulting data of EN<br>
`md_load_md_result`:mid-dist:Load the resulting data of mid-dist<br>
`MINDIST`:mid-dist:Calculate the resulting data<br>
`solutions.h`:Base class and function for EN,LNB,RID,RID-Q,greedy_k RID<br>
`EN`:EN algorithm<br>
`LNB_top_visito`r:top event visitor for LNB algorithm<br>
`LNB_extend_visitor`:extend event visitor for LNB algorithm<br>
`LNB_construct_LNT`: LNB algorithm:construct LNT and L2NT<br>
`LNB_query`:LNB algorithm:query optimal facility-replacement pair based on LNT & L2NT<br>
`NNFC_visitor`:top event visitor for calculate NNFC<br>
`NSJ_visitor`:top event visitor for NSJ algorithm<br>
`NSJ_construct_NNFC`:NSJ algorithm:construct NNFCs R-tree<br>
`NSJ_query`:NSJ algorithm:query optimal candidate based on NNFCs<br>
`NNFC_RIC_visitor`:top event visitor for calculate NNFC & RIC<br>
`RID_construct_NNFC_and_RIC`: RID algorithm:compute dr(f), i.e., \Delta(f), construct NNFC MBRs for reference locations and RIC MBRs for facilities<br>                                            		
`RID_query`:RID algorithm:query optimal candidate based on RICs<br>
`RID_new_query`:greedy_k RID:query optimal candidate based on RICs<br>
`RID_update_NNFC_and_RIC`:greedy_k RID:compute dr(f), i.e., \Delta(f), construct NNFC MBRs for references<br>       
`load_ranked_data`:Load real or query data<br>
`APatK`:Calculate Precision (P),Recall (R) and Average Precision (AP)<br>
		
`MFRM.cpp`:Main function<br>
`compare_EN_with_MinDist`:Calculate for comparing with min-dist FR problems<br>
`results_for_EN_MinDist `<br>                                                        
		                 
`load_ranked_data`:Calculate AP for two ranked sequences<br>
`APatK`:Implementation of all algorithoms (eg., EN,LNB,RID) in our paper<br>        

Data
-----
1.There are two datasets which are recorded by two folder:`BJ_datasets `and` CA_datasets`.<br>
2.The specific data file includes the following contents:<br>
`gen_edges`:road network data.<br>
`gen_cands_`:candidates location data.<br>
`gen_avg_reflocs_`:uniform distribution of reference locations data.<br>
`gen_rand_reflocs_`:random distribution of reference locations data.<br>
`real_reflocs`:real users'reference locations data.<br>
`800_fac_xy`:facility location data with coordinate.<br>
`800_fac`:facility location data without coordinate.

data format
-------
1.The edges file format is "edge_id vs_id ve_id dist".<br>
2.The facilities file format is "fac_id vs_id ve_id dist_vs_fac dist_fac_ve" or "fac_id vs_id ve_id dist_vs_fac dist_fac_ve lon lat".<br> 
3.The candidates file format is "cand_id vs_id ve_id dist_vs_cand dist_cand_ve lon lat".<br>
4.The reference locations file format is "refloc_id vs_id ve_id dist_vs_refloc dist_refloc_ve prob lon lat".<br>

Usage
-----
1.The program must be run in the specified directory:`'D:\Experiment\MFRM'`.<br>
2.In the configuration file,1 means execute and 0 means not execute,it can be modified according to the experimental requirements.
3. parameter setting:<br>

    Implements of all algorithoms:'...\config_CA.txt',' ...\gen_config.txt', '...\cov_config.txt', '...\real_config.txt'
    Calculate for comparing with min-dist FR problems:first step:'-md' '...\config_MD.txt' second step:'-mdr' '...\config_MD.txt'
