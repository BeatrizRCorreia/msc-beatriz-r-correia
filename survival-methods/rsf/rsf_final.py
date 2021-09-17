import time
import random
import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from sklearn.model_selection import train_test_split
from pysurvival.models.simulations import SimulationModel
from pysurvival.models.survival_forest import RandomSurvivalForestModel
from pysurvival.utils.metrics import concordance_index
# from pysurvival.utils.display import integrated_brier_score
from pysurvival.utils.metrics import integrated_brier_score
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale
from pandas import DataFrame
from collections import Counter
from operator import itemgetter

seed_value = random.seed(5)

def prepareData():
	# To run on my pc
	# barcode_mapping = pd.read_csv("/Users/beatrizcorreia/Documents/msc-beatriz-r-correia/data/BRCA_primary_solid_tumor/data2_sampleMap.csv", index_col = 0)
	# transposed_data = pd.read_csv("/Users/beatrizcorreia/Documents/msc-beatriz-r-correia/data/BRCA_primary_solid_tumor/data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv")
	# col_data = pd.read_csv("/Users/beatrizcorreia/Documents/msc-beatriz-r-correia/data/BRCA_primary_solid_tumor/data2_colData.csv", dtype = {"vital_status": int})

	# To run on server
	barcode_mapping = pd.read_csv("/home/draconian_mule/msc-beatriz-correia/data/BRCA_primary_solid_tumor/data2_sampleMap.csv", index_col = 0)
	transposed_data = pd.read_csv("/home/draconian_mule/msc-beatriz-correia/data/BRCA_primary_solid_tumor/data2_01_BRCA_RNASeq2GeneNorm-20160128_transposed.csv")
	col_data = pd.read_csv("/home/draconian_mule/msc-beatriz-correia/data/BRCA_primary_solid_tumor/data2_colData.csv", dtype = {"vital_status": int})

	# print(col_data)
	
	# print("female: ", col_data[col_data.gender == 'female'].shape[0])
	# print("male: ", col_data[col_data.gender == 'male'].shape[0])

	# print("vital_status 0: ", col_data[col_data.vital_status == 1].shape[0])
	# print("vital_status 1: ", col_data[col_data.vital_status == 0].shape[0])
	
	# print(col_data["pathologic_stage"].value_counts())
	# print(col_data["gender"].value_counts())
	# print(col_data["vital_status"].value_counts())
	# print(col_data["date_of_initial_pathologic_diagnosis"].value_counts().sort_index())
	# print(col_data["patient.age_at_initial_pathologic_diagnosis"].value_counts().sort_index().to_dict())

	# print(col_data["date_of_initial_pathologic_diagnosis"].isna().sum())
	# print(col_data["patient.age_at_initial_pathologic_diagnosis"].isna().sum())

	col_data_filtered = col_data[["patientID", "vital_status", "Days.to.Date.of.Last.Contact", "Days.to.date.of.Death", "days_to_death", "days_to_last_followup"]]
	col_data_filtered = col_data_filtered.set_index("patientID")
	# print(col_data_filtered["Days.to.Date.of.Last.Contact"].tolist())
	new_short_barcodes = barcode_mapping.loc[:,"primary"]
	new_transposed_data = transposed_data.assign(patientID = new_short_barcodes.tolist())
	new_transposed_data = new_transposed_data.drop("Unnamed: 0", axis = 1)
	new_transposed_data = new_transposed_data.set_index("patientID")
	genes = new_transposed_data.columns.tolist()
	# Scaling before train and test splits
	# new_transposed_data = DataFrame(scale(new_transposed_data), index=new_transposed_data.index, columns=new_transposed_data.columns)
	# new_transposed_data = new_transposed_data.loc[:, new_transposed_data.std() > 0]
	# print(new_transposed_data)
	# print(col_data_filtered.loc["TCGA-E9-A245",:])
	col_data_filtered["Days.to.Date.of.Last.Contact"] = col_data_filtered["Days.to.Date.of.Last.Contact"].fillna(-1)
	col_data_filtered["Days.to.date.of.Death"] = col_data_filtered["Days.to.date.of.Death"].fillna(-1)
	col_data_filtered["days_to_death"] = col_data_filtered["days_to_death"].fillna(-1)
	col_data_filtered["days_to_last_followup"] = col_data_filtered["days_to_last_followup"].fillna(-1)
	# print(col_data_filtered.loc["TCGA-E9-A245",:])
	# print(col_data_filtered["Days.to.Date.of.Last.Contact"])
	col_data_filtered["Days.to.Date.of.Last.Contact"] = col_data_filtered["Days.to.Date.of.Last.Contact"].replace("[Not Available]", -1)
	col_data_filtered["Days.to.Date.of.Last.Contact"] = col_data_filtered["Days.to.Date.of.Last.Contact"].replace("[Completed]", -1)

	col_data_filtered["Days.to.Date.of.Last.Contact"] = col_data_filtered["Days.to.Date.of.Last.Contact"].astype(str).astype(float)
	col_data_filtered["Days.to.date.of.Death"] = col_data_filtered["Days.to.date.of.Death"].astype(str).astype(float)
	col_data_filtered["days_to_death"] = col_data_filtered["days_to_death"].astype(str).astype(float)
	col_data_filtered["days_to_last_followup"] = col_data_filtered["days_to_last_followup"].astype(str).astype(float)
	
	col_data_filtered["Days.to.Date.of.Last.Contact"] = col_data_filtered["Days.to.Date.of.Last.Contact"].replace(-1, np.nan)
	col_data_filtered["Days.to.date.of.Death"] = col_data_filtered["Days.to.date.of.Death"].replace(-1, np.nan)
	col_data_filtered["days_to_death"] = col_data_filtered["days_to_death"].replace(-1, np.nan)
	col_data_filtered["days_to_last_followup"] = col_data_filtered["days_to_last_followup"].replace(-1, np.nan)

	# print(col_data_filtered[["Days.to.Date.of.Last.Contact", "Days.to.date.of.Death", "days_to_death", "days_to_last_followup"]].max(axis = 1).astype(int))
	col_data_filtered["time"] = col_data_filtered[["Days.to.Date.of.Last.Contact", "Days.to.date.of.Death", "days_to_death", "days_to_last_followup"]].max(axis = 1).astype(int)
	# print(col_data_filtered.loc["TCGA-E9-A245",:])
	col_data_filtered = col_data_filtered.drop("Days.to.Date.of.Last.Contact", axis = 1)
	col_data_filtered = col_data_filtered.drop("Days.to.date.of.Death", axis = 1)
	col_data_filtered = col_data_filtered.drop("days_to_death", axis = 1)
	col_data_filtered = col_data_filtered.drop("days_to_last_followup", axis = 1)
	# print("COLDATA BEFORE ELIMINATING PATIENTS:")
	# print(col_data_filtered)
	# print("TIME NA PATIENT:")
	# print(col_data_filtered[col_data_filtered['time'].isna()])
	# print("TIME <= 0:")
	# print(col_data_filtered[col_data_filtered['time'] <= 0])
	col_data_filtered = col_data_filtered[col_data_filtered['time'].notna()]
	# print("COLDATA AFTER DROPPING NA PATIENTS:")
	# print(col_data_filtered)
	col_data_filtered = col_data_filtered[col_data_filtered['time'] > 0]
	# print("COLDATA AFTER DROPPING TIME <= 0 PATIENTS")
	# print(col_data_filtered)
	all_data = new_transposed_data.merge(col_data_filtered, how = "inner", on = "patientID")
	# print(all_data)
	N = all_data.shape[0]
	return all_data, genes, N
	# return new_transposed_data, col_data_filtered

def get_variable_importance(rsf, results_file_name):
	#print("\nVARIABLE IMPORTANCE TABLE ---------------------------", file=open(results_file_name, "a"))
	variable_importance_table = getattr(rsf, "variable_importance_table")
	#print(variable_importance_table, file=open(results_file_name, "a"))
	#print("\nFEATURES FROM TABLE:", file=open(results_file_name, "a"))
	features_from_table = variable_importance_table['feature'].tolist()
	#print(features_from_table, file=open(results_file_name, "a"))
	#print("\nIMPORTANCE FROM TABLE:", file=open(results_file_name, "a"))
	importance_from_table = variable_importance_table['importance'].tolist()
	#print(importance_from_table, file=open(results_file_name, "a"))
	#print("\nPCT IMPORTANCE FROM TABLE:", file=open(results_file_name, "a"))
	#pctimportance_from_table = variable_importance_table['pct_importance'].tolist()
	#print(pctimportance_from_table, file=open(results_file_name, "a"))

	zip_iterator = zip(features_from_table, importance_from_table)
	dictionary_of_vimps = dict(zip_iterator)
	return dictionary_of_vimps

def std_calc(list_of_values, average):
	N = len(list_of_values)
	numerator = 0
	for i in list_of_values:
		numerator = numerator + (i - average)**2
	all_std_expression = math.sqrt(numerator/N)
	return all_std_expression

def RSF_test_number_of_trees(number_of_trees, iterations_per_nr_of_trees, max_features, max_depth, min_node_size, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test):
	start_time1 = time.time()
	all_vimps_dict = {}
	all_cindexes_test = []
	all_cindexes_train = []
	all_ibss = []
	all_std_cindexes_test = []
	all_std_cindexes_train = []
	all_std_ibss = []
	all_results_file_name = "RSF_all_results.txt"
	for i in number_of_trees:
		iterations_vimps_dict = {}
		iterations_cindexes_test = []
		iterations_cindexes_train = []
		iterations_ibss = []
		iterations_file_name = "RSF_nr_of_trees_" + str(i)
		iterations_file_name = iterations_file_name + ".txt"
		print("\nNUMBER OF TREES: ", i, "******************************************************")
		print("\nNUMBER OF TREES: ", i, "******************************************************", file=open(all_results_file_name, "a"))
		start_time3 = time.time()
		for x in range(0, iterations_per_nr_of_trees):
			print(random.randint(1, 999999))
			# iterations_file_name = "RSF_nr_of_trees_testtt.txt"
			print("\nITERATION ", x+1, "-------------------------------------------------------")
			print("\nITERATION ", x+1, "-------------------------------------------------------", file=open(iterations_file_name, "a"))
			start_time2 = time.time()
			rsf = RandomSurvivalForestModel(num_trees=i)
			rsf.fit(scaled_X_train, T_train, E_train, max_features=max_features, max_depth=max_depth, min_node_size=min_node_size, seed=random.randint(1, 999999))

			vimp_dict = get_variable_importance(rsf, iterations_file_name)
			# print("\nvimp_dict: ", vimp_dict, file=open(iterations_file_name, "a"))
			iterations_vimps_dict = Counter(iterations_vimps_dict)
			iterations_vimps_dict.update(Counter(vimp_dict))

			# MODEL PERFORMANCE METRICS
			# c_index test
			c_index_test = concordance_index(rsf, scaled_X_test, T_test, E_test) #0.81
			iterations_cindexes_test.append(c_index_test)
			print('C-index test: {:.4f}'.format(c_index_test), file=open(iterations_file_name, "a"))

			# c_index_train
			c_index_train = concordance_index(rsf, scaled_X_train, T_train, E_train) #0.81
			iterations_cindexes_train.append(c_index_train)
			print('C-index train: {:.4f}'.format(c_index_train), file=open(iterations_file_name, "a"))

			# ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test, t_max=30, figure_size=(20, 6.5))
			ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test)
			iterations_ibss.append(ibs)
			print('IBS: {:.4f}'.format(ibs), file=open(iterations_file_name, "a"))

			end_time2 = time.time()
			hours2, rem2 = divmod(end_time2-start_time2, 3600)
			minutes2, seconds2 = divmod(rem2, 60)
			print("\nFinished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2))
			print("Finished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2), file=open(iterations_file_name, "a"))

		iterations_vimps_dict = dict(iterations_vimps_dict)
		iterations_vimps_dict = {k: v / iterations_per_nr_of_trees for k, v in iterations_vimps_dict.items() if v != 0}
		# print("\niterations_vimps_dict: ", iterations_vimps_dict, file=open(iterations_file_name, "a"))
		list_from_dict1 = [(k, v) for k, v in iterations_vimps_dict.items()]
		vimps_list_ordered1 = sorted(list_from_dict1, key=itemgetter(1))
		descending_vimps_list1 = vimps_list_ordered1[::-1]
		print("\ndescending_vimps_list: ", descending_vimps_list1, file=open(iterations_file_name, "a"))

		average_iterations_cindexes_test = sum(iterations_cindexes_test) / len(iterations_cindexes_test)
		all_cindexes_test.append(average_iterations_cindexes_test)
		print('\nAverage c-index test: {:.4f}'.format(average_iterations_cindexes_test), file=open(all_results_file_name, "a"))
		std_cindexes_test = std_calc(iterations_cindexes_test, average_iterations_cindexes_test)
		all_std_cindexes_test.append(std_cindexes_test)
		print('Standard deviation: {:.4f}'.format(std_cindexes_test), file=open(all_results_file_name, "a"))

		average_iterations_cindexes_train = sum(iterations_cindexes_train) / len(iterations_cindexes_train)
		all_cindexes_train.append(average_iterations_cindexes_train)
		print('\nAverage c-index train: {:.4f}'.format(average_iterations_cindexes_train), file=open(all_results_file_name, "a"))
		std_cindexes_train = std_calc(iterations_cindexes_train, average_iterations_cindexes_train)
		all_std_cindexes_train.append(std_cindexes_train)
		print('Standard deviation: {:.4f}'.format(std_cindexes_train), file=open(all_results_file_name, "a"))

		average_iterations_ibs = sum(iterations_ibss) / len(iterations_ibss)
		all_ibss.append(average_iterations_ibs)
		print('\nAverage IBS: {:.4f}'.format(average_iterations_ibs), file=open(all_results_file_name, "a"))
		std_ibss = std_calc(iterations_ibss, average_iterations_ibs)
		all_std_ibss.append(std_ibss)
		print('Standard deviation: {:.4f}'.format(std_ibss), file=open(all_results_file_name, "a"))

		all_vimps_dict = Counter(all_vimps_dict)
		all_vimps_dict.update(Counter(iterations_vimps_dict))

		end_time3 = time.time()
		hours3, rem3 = divmod(end_time3-start_time3, 3600)
		minutes3, seconds3 = divmod(rem3, 60)
		print("\nFinished the iterations for", i, "trees in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3))
		print("\nFinished the iterations for", i, "trees in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3), file=open(all_results_file_name, "a"))

	all_vimps_dict = dict(all_vimps_dict)
	all_vimps_dict = {k: v / len(number_of_trees) for k, v in all_vimps_dict.items() if v != 0}
	# print("\nall_vimps_dict: ", all_vimps_dict, file=open(all_results_file_name, "a"))
	list_from_dict2 = [(k, v) for k, v in all_vimps_dict.items()]
	vimps_list_ordered2 = sorted(list_from_dict2, key=itemgetter(1))
	descending_vimps_list2 = vimps_list_ordered2[::-1]
	print("\ndescending_vimps_list: ", descending_vimps_list2, file=open(all_results_file_name, "a"))
	positive_vimps_genes = [item for item in descending_vimps_list2 if item[1] > 0]
	negative_vimps_genes = [item for item in descending_vimps_list2 if item[1] < 0]
	print("\nGenes and their positive vimps: ", positive_vimps_genes, file=open(all_results_file_name, "a"))
	list_of_genes_with_positive_vimps = [tup[0] for tup in positive_vimps_genes]
	list_of_genes_with_negative_vimps = [tup[0] for tup in negative_vimps_genes]
	print("\nList of genes with positive vimps: ", list_of_genes_with_positive_vimps, file=open(all_results_file_name, "a"))
	print("\nNumber of genes with positive vimps: ", len(list_of_genes_with_positive_vimps), file=open(all_results_file_name, "a"))
	print("\nNumber of genes with negative vimps: ", len(list_of_genes_with_negative_vimps), file=open(all_results_file_name, "a"))

	print("\n***********\n**SUMARY**\n***********", file=open(all_results_file_name, "a"))
	print("\nNumber of trees: ", number_of_trees[::-1], file=open(all_results_file_name, "a"))
	print("\nAll c-indexes test: ", all_cindexes_test[::-1], file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_test, file=open(all_results_file_name, "a"))
	print("\nAll c-indexes train: ", all_cindexes_train[::-1], file=open(all_results_file_name, "a"))
	print("\nAll ibs: ", all_ibss[::-1], file=open(all_results_file_name, "a"))

	end_time1 = time.time()
	hours1, rem1 = divmod(end_time1-start_time1, 3600)
	minutes1, seconds1 = divmod(rem1, 60)
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1))
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1), file=open(all_results_file_name, "a"))

	plt.figure(figsize=(12, 6))
	line1, = plt.plot(number_of_trees[::-1], all_cindexes_train[::-1], color = (1, 0, 0.8), marker='D', label="Train set c-index")
	line2, = plt.plot(number_of_trees[::-1], all_cindexes_test[::-1], color = (0.2, 0.5, 0.9), marker='D', label="Test set c-index")
	plt.grid(linestyle = '-.')
	plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
	plt.ylabel("C-index")
	plt.xlabel("Number of trees in the forest")
	# plt.show()
	plt.savefig("RSF_all_results_1.png")
	return

def RSF_test_tree_depth(max_depths, iterations_per_nr_of_trees, max_features, number_of_trees, min_node_size, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test):
	start_time1 = time.time()
	all_vimps_dict = {}
	all_cindexes_test = []
	all_cindexes_train = []
	all_ibss = []
	all_std_cindexes_test = []
	all_std_cindexes_train = []
	all_std_ibss = []
	all_results_file_name = "RSF_all_results_tree_depth.txt"
	for i in max_depths:
		iterations_vimps_dict = {}
		iterations_cindexes_test = []
		iterations_cindexes_train = []
		iterations_ibss = []
		iterations_file_name = "RSF_max_tree_depth_" + str(i)
		iterations_file_name = iterations_file_name + ".txt"
		print("\nMAX TREE DEPTH: ", i, "******************************************************")
		print("\nMAX TREE DEPTH: ", i, "******************************************************", file=open(all_results_file_name, "a"))
		start_time3 = time.time()
		for x in range(0, iterations_per_nr_of_trees):
			print(random.randint(1, 999999))
			# iterations_file_name = "RSF_nr_of_trees_testtt.txt"
			print("\nITERATION ", x+1, "-------------------------------------------------------")
			print("\nITERATION ", x+1, "-------------------------------------------------------", file=open(iterations_file_name, "a"))
			start_time2 = time.time()
			rsf = RandomSurvivalForestModel(num_trees=number_of_trees)
			rsf.fit(scaled_X_train, T_train, E_train, max_features=max_features, max_depth=i, min_node_size=min_node_size, seed=random.randint(1, 999999))

			vimp_dict = get_variable_importance(rsf, iterations_file_name)
			# print("\nvimp_dict: ", vimp_dict, file=open(iterations_file_name, "a"))
			iterations_vimps_dict = Counter(iterations_vimps_dict)
			iterations_vimps_dict.update(Counter(vimp_dict))

			# MODEL PERFORMANCE METRICS
			# c_index test
			c_index_test = concordance_index(rsf, scaled_X_test, T_test, E_test) #0.81
			iterations_cindexes_test.append(c_index_test)
			print('C-index test: {:.4f}'.format(c_index_test), file=open(iterations_file_name, "a"))

			# c_index_train
			c_index_train = concordance_index(rsf, scaled_X_train, T_train, E_train) #0.81
			iterations_cindexes_train.append(c_index_train)
			print('C-index train: {:.4f}'.format(c_index_train), file=open(iterations_file_name, "a"))

			# ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test, t_max=30, figure_size=(20, 6.5))
			ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test)
			iterations_ibss.append(ibs)
			print('IBS: {:.4f}'.format(ibs), file=open(iterations_file_name, "a"))

			end_time2 = time.time()
			hours2, rem2 = divmod(end_time2-start_time2, 3600)
			minutes2, seconds2 = divmod(rem2, 60)
			print("\nFinished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2))
			print("Finished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2), file=open(iterations_file_name, "a"))

		iterations_vimps_dict = dict(iterations_vimps_dict)
		iterations_vimps_dict = {k: v / iterations_per_nr_of_trees for k, v in iterations_vimps_dict.items() if v != 0}
		# print("\niterations_vimps_dict: ", iterations_vimps_dict, file=open(iterations_file_name, "a"))
		list_from_dict1 = [(k, v) for k, v in iterations_vimps_dict.items()]
		vimps_list_ordered1 = sorted(list_from_dict1, key=itemgetter(1))
		descending_vimps_list1 = vimps_list_ordered1[::-1]
		print("\ndescending_vimps_list: ", descending_vimps_list1, file=open(iterations_file_name, "a"))

		average_iterations_cindexes_test = sum(iterations_cindexes_test) / len(iterations_cindexes_test)
		all_cindexes_test.append(average_iterations_cindexes_test)
		print('\nAverage c-index test: {:.4f}'.format(average_iterations_cindexes_test), file=open(all_results_file_name, "a"))
		std_cindexes_test = std_calc(iterations_cindexes_test, average_iterations_cindexes_test)
		all_std_cindexes_test.append(std_cindexes_test)
		print('Standard deviation: {:.4f}'.format(std_cindexes_test), file=open(all_results_file_name, "a"))

		average_iterations_cindexes_train = sum(iterations_cindexes_train) / len(iterations_cindexes_train)
		all_cindexes_train.append(average_iterations_cindexes_train)
		print('\nAverage c-index train: {:.4f}'.format(average_iterations_cindexes_train), file=open(all_results_file_name, "a"))
		std_cindexes_train = std_calc(iterations_cindexes_train, average_iterations_cindexes_train)
		all_std_cindexes_train.append(std_cindexes_train)
		print('Standard deviation: {:.4f}'.format(std_cindexes_train), file=open(all_results_file_name, "a"))

		average_iterations_ibs = sum(iterations_ibss) / len(iterations_ibss)
		all_ibss.append(average_iterations_ibs)
		print('\nAverage IBS: {:.4f}'.format(average_iterations_ibs), file=open(all_results_file_name, "a"))
		std_ibss = std_calc(iterations_ibss, average_iterations_ibs)
		all_std_ibss.append(std_ibss)
		print('Standard deviation: {:.4f}'.format(std_ibss), file=open(all_results_file_name, "a"))

		all_vimps_dict = Counter(all_vimps_dict)
		all_vimps_dict.update(Counter(iterations_vimps_dict))

		end_time3 = time.time()
		hours3, rem3 = divmod(end_time3-start_time3, 3600)
		minutes3, seconds3 = divmod(rem3, 60)
		print("\nFinished the iterations for max depth =", i, "in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3))
		print("\nFinished the iterations for max depth =", i, "in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3), file=open(all_results_file_name, "a"))

	all_vimps_dict = dict(all_vimps_dict)
	all_vimps_dict = {k: v / len(max_depths) for k, v in all_vimps_dict.items() if v != 0}
	# print("\nall_vimps_dict: ", all_vimps_dict, file=open(all_results_file_name, "a"))
	list_from_dict2 = [(k, v) for k, v in all_vimps_dict.items()]
	vimps_list_ordered2 = sorted(list_from_dict2, key=itemgetter(1))
	descending_vimps_list2 = vimps_list_ordered2[::-1]
	print("\ndescending_vimps_list: ", descending_vimps_list2, file=open(all_results_file_name, "a"))
	positive_vimps_genes = [item for item in descending_vimps_list2 if item[1] > 0]
	negative_vimps_genes = [item for item in descending_vimps_list2 if item[1] < 0]
	print("\nGenes and their positive vimps: ", positive_vimps_genes, file=open(all_results_file_name, "a"))
	list_of_genes_with_positive_vimps = [tup[0] for tup in positive_vimps_genes]
	list_of_genes_with_negative_vimps = [tup[0] for tup in negative_vimps_genes]
	print("\nList of genes with positive vimps: ", list_of_genes_with_positive_vimps, file=open(all_results_file_name, "a"))
	print("\nNumber of genes with positive vimps: ", len(list_of_genes_with_positive_vimps), file=open(all_results_file_name, "a"))
	print("\nNumber of genes with negative vimps: ", len(list_of_genes_with_negative_vimps), file=open(all_results_file_name, "a"))

	print("\n***********\n**SUMARY**\n***********", file=open(all_results_file_name, "a"))
	print("\nMax depths: ", max_depths[::-1], file=open(all_results_file_name, "a"))
	print("\nAll c-indexes test: ", all_cindexes_test[::-1], file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_test, file=open(all_results_file_name, "a"))
	print("\nAll c-indexes train: ", all_cindexes_train[::-1], file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_train, file=open(all_results_file_name, "a"))
	print("\nAll ibs: ", all_ibss[::-1], file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_ibss, file=open(all_results_file_name, "a"))

	end_time1 = time.time()
	hours1, rem1 = divmod(end_time1-start_time1, 3600)
	minutes1, seconds1 = divmod(rem1, 60)
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1))
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1), file=open(all_results_file_name, "a"))

	plt.figure(figsize=(12, 6))
	line1, = plt.plot(max_depths[::-1], all_cindexes_train[::-1], color = (1, 0, 0.8), marker='D', label="Train set c-index")
	line2, = plt.plot(max_depths[::-1], all_cindexes_test[::-1], color = (0.2, 0.5, 0.9), marker='D', label="Test set c-index")
	plt.grid(linestyle = '-.')
	plt.xticks(max_depths[::-1])
	plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
	plt.ylabel("C-index")
	plt.xlabel("Max tree depths")
	# plt.show()
	plt.savefig("RSF_all_results_2.png")

def RSF_test_min_node_size(min_node_sizes, iterations_per_nr_of_trees, max_features, number_of_trees, max_depth, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test):
	start_time1 = time.time()
	all_vimps_dict = {}
	all_cindexes_test = []
	all_cindexes_train = []
	all_ibss = []
	all_std_cindexes_test = []
	all_std_cindexes_train = []
	all_std_ibss = []
	all_results_file_name = "RSF_all_results_min_node_size.txt"
	for i in min_node_sizes:
		iterations_vimps_dict = {}
		iterations_cindexes_test = []
		iterations_cindexes_train = []
		iterations_ibss = []
		iterations_file_name = "RSF_min_node_size_" + str(i)
		iterations_file_name = iterations_file_name + ".txt"
		print("\nMIN NODE SIZE: ", i, "******************************************************")
		print("\nMIN NODE SIZE: ", i, "******************************************************", file=open(all_results_file_name, "a"))
		start_time3 = time.time()
		for x in range(0, iterations_per_nr_of_trees):
			print(random.randint(1, 999999))
			# iterations_file_name = "RSF_nr_of_trees_testtt.txt"
			print("\nITERATION ", x+1, "-------------------------------------------------------")
			print("\nITERATION ", x+1, "-------------------------------------------------------", file=open(iterations_file_name, "a"))
			start_time2 = time.time()
			rsf = RandomSurvivalForestModel(num_trees=number_of_trees)
			rsf.fit(scaled_X_train, T_train, E_train, max_features=max_features, max_depth=max_depth, min_node_size=i, seed=random.randint(1, 999999))

			vimp_dict = get_variable_importance(rsf, iterations_file_name)
			# print("\nvimp_dict: ", vimp_dict, file=open(iterations_file_name, "a"))
			iterations_vimps_dict = Counter(iterations_vimps_dict)
			iterations_vimps_dict.update(Counter(vimp_dict))

			# MODEL PERFORMANCE METRICS
			# c_index test
			c_index_test = concordance_index(rsf, scaled_X_test, T_test, E_test) #0.81
			iterations_cindexes_test.append(c_index_test)
			print('C-index test: {:.4f}'.format(c_index_test), file=open(iterations_file_name, "a"))

			# c_index_train
			c_index_train = concordance_index(rsf, scaled_X_train, T_train, E_train) #0.81
			iterations_cindexes_train.append(c_index_train)
			print('C-index train: {:.4f}'.format(c_index_train), file=open(iterations_file_name, "a"))

			# ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test, t_max=30, figure_size=(20, 6.5))
			ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test)
			iterations_ibss.append(ibs)
			print('IBS: {:.4f}'.format(ibs), file=open(iterations_file_name, "a"))

			end_time2 = time.time()
			hours2, rem2 = divmod(end_time2-start_time2, 3600)
			minutes2, seconds2 = divmod(rem2, 60)
			print("\nFinished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2))
			print("Finished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2), file=open(iterations_file_name, "a"))

		iterations_vimps_dict = dict(iterations_vimps_dict)
		iterations_vimps_dict = {k: v / iterations_per_nr_of_trees for k, v in iterations_vimps_dict.items() if v != 0}
		# print("\niterations_vimps_dict: ", iterations_vimps_dict, file=open(iterations_file_name, "a"))
		list_from_dict1 = [(k, v) for k, v in iterations_vimps_dict.items()]
		vimps_list_ordered1 = sorted(list_from_dict1, key=itemgetter(1))
		descending_vimps_list1 = vimps_list_ordered1[::-1]
		print("\ndescending_vimps_list: ", descending_vimps_list1, file=open(iterations_file_name, "a"))

		average_iterations_cindexes_test = sum(iterations_cindexes_test) / len(iterations_cindexes_test)
		all_cindexes_test.append(average_iterations_cindexes_test)
		print('\nAverage c-index test: {:.4f}'.format(average_iterations_cindexes_test), file=open(all_results_file_name, "a"))
		std_cindexes_test = std_calc(iterations_cindexes_test, average_iterations_cindexes_test)
		all_std_cindexes_test.append(std_cindexes_test)
		print('Standard deviation: {:.4f}'.format(std_cindexes_test), file=open(all_results_file_name, "a"))

		average_iterations_cindexes_train = sum(iterations_cindexes_train) / len(iterations_cindexes_train)
		all_cindexes_train.append(average_iterations_cindexes_train)
		print('\nAverage c-index train: {:.4f}'.format(average_iterations_cindexes_train), file=open(all_results_file_name, "a"))
		std_cindexes_train = std_calc(iterations_cindexes_train, average_iterations_cindexes_train)
		all_std_cindexes_train.append(std_cindexes_train)
		print('Standard deviation: {:.4f}'.format(std_cindexes_train), file=open(all_results_file_name, "a"))

		average_iterations_ibs = sum(iterations_ibss) / len(iterations_ibss)
		all_ibss.append(average_iterations_ibs)
		print('\nAverage IBS: {:.4f}'.format(average_iterations_ibs), file=open(all_results_file_name, "a"))
		std_ibss = std_calc(iterations_ibss, average_iterations_ibs)
		all_std_ibss.append(std_ibss)
		print('Standard deviation: {:.4f}'.format(std_ibss), file=open(all_results_file_name, "a"))

		all_vimps_dict = Counter(all_vimps_dict)
		all_vimps_dict.update(Counter(iterations_vimps_dict))

		end_time3 = time.time()
		hours3, rem3 = divmod(end_time3-start_time3, 3600)
		minutes3, seconds3 = divmod(rem3, 60)
		print("\nFinished the iterations for min node size =", i, "in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3))
		print("\nFinished the iterations for min node size =", i, "in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3), file=open(all_results_file_name, "a"))

	all_vimps_dict = dict(all_vimps_dict)
	all_vimps_dict = {k: v / len(min_node_sizes) for k, v in all_vimps_dict.items() if v != 0}
	# print("\nall_vimps_dict: ", all_vimps_dict, file=open(all_results_file_name, "a"))
	list_from_dict2 = [(k, v) for k, v in all_vimps_dict.items()]
	vimps_list_ordered2 = sorted(list_from_dict2, key=itemgetter(1))
	descending_vimps_list2 = vimps_list_ordered2[::-1]
	print("\ndescending_vimps_list: ", descending_vimps_list2, file=open(all_results_file_name, "a"))
	positive_vimps_genes = [item for item in descending_vimps_list2 if item[1] > 0]
	negative_vimps_genes = [item for item in descending_vimps_list2 if item[1] < 0]
	print("\nGenes and their positive vimps: ", positive_vimps_genes, file=open(all_results_file_name, "a"))
	list_of_genes_with_positive_vimps = [tup[0] for tup in positive_vimps_genes]
	list_of_genes_with_negative_vimps = [tup[0] for tup in negative_vimps_genes]
	print("\nList of genes with positive vimps: ", list_of_genes_with_positive_vimps, file=open(all_results_file_name, "a"))
	print("\nNumber of genes with positive vimps: ", len(list_of_genes_with_positive_vimps), file=open(all_results_file_name, "a"))
	print("\nNumber of genes with negative vimps: ", len(list_of_genes_with_negative_vimps), file=open(all_results_file_name, "a"))

	print("\n***********\n**SUMARY**\n***********", file=open(all_results_file_name, "a"))
	print("\nMin node sizes: ", min_node_sizes, file=open(all_results_file_name, "a"))
	print("\nAll c-indexes test: ", all_cindexes_test, file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_test, file=open(all_results_file_name, "a"))
	print("\nAll c-indexes train: ", all_cindexes_train, file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_train, file=open(all_results_file_name, "a"))
	print("\nAll ibs: ", all_ibss, file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_ibss, file=open(all_results_file_name, "a"))

	end_time1 = time.time()
	hours1, rem1 = divmod(end_time1-start_time1, 3600)
	minutes1, seconds1 = divmod(rem1, 60)
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1))
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1), file=open(all_results_file_name, "a"))

	plt.figure(figsize=(12, 6))
	line1, = plt.plot(min_node_sizes, all_cindexes_train, color = (1, 0, 0.8), marker='D', label="Train set c-index")
	line2, = plt.plot(min_node_sizes, all_cindexes_test, color = (0.2, 0.5, 0.9), marker='D', label="Test set c-index")
	plt.grid(linestyle = '-.')
	plt.xticks(min_node_sizes)
	plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
	plt.ylabel("C-index")
	plt.xlabel("Min node sizes")
	# plt.show()
	plt.savefig("RSF_all_results_3.png")

def RSF_test_best_model(min_node_sizes, iterations_per_nr_of_trees, max_features, number_of_trees, max_depth, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test):
	start_time1 = time.time()
	all_vimps_dict = {}
	all_cindexes_test = []
	all_cindexes_train = []
	all_ibss = []
	all_std_cindexes_test = []
	all_std_cindexes_train = []
	all_std_ibss = []
	all_results_file_name = "RSF_all_results_test_best_model.txt"
	for i in min_node_sizes:
		iterations_vimps_dict = {}
		iterations_cindexes_test = []
		iterations_cindexes_train = []
		iterations_ibss = []
		iterations_file_name = "RSF_min_node_size_" + str(i)
		iterations_file_name = iterations_file_name + ".txt"
		print("\nMIN NODE SIZE: ", i, "******************************************************")
		print("\nMIN NODE SIZE: ", i, "******************************************************", file=open(all_results_file_name, "a"))
		start_time3 = time.time()
		for x in range(0, iterations_per_nr_of_trees):
			print(random.randint(1, 999999))
			# iterations_file_name = "RSF_nr_of_trees_testtt.txt"
			print("\nITERATION ", x+1, "-------------------------------------------------------")
			print("\nITERATION ", x+1, "-------------------------------------------------------", file=open(iterations_file_name, "a"))
			start_time2 = time.time()
			rsf = RandomSurvivalForestModel(num_trees=number_of_trees)
			rsf.fit(scaled_X_train, T_train, E_train, max_features=max_features, max_depth=max_depth, min_node_size=i, seed=random.randint(1, 999999))

			vimp_dict = get_variable_importance(rsf, iterations_file_name)
			# print("\nvimp_dict: ", vimp_dict, file=open(iterations_file_name, "a"))
			iterations_vimps_dict = Counter(iterations_vimps_dict)
			iterations_vimps_dict.update(Counter(vimp_dict))

			# MODEL PERFORMANCE METRICS
			# c_index test
			c_index_test = concordance_index(rsf, scaled_X_test, T_test, E_test) #0.81
			iterations_cindexes_test.append(c_index_test)
			print('C-index test: {:.4f}'.format(c_index_test), file=open(iterations_file_name, "a"))

			# c_index_train
			c_index_train = concordance_index(rsf, scaled_X_train, T_train, E_train) #0.81
			iterations_cindexes_train.append(c_index_train)
			print('C-index train: {:.4f}'.format(c_index_train), file=open(iterations_file_name, "a"))

			# ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test, t_max=30, figure_size=(20, 6.5))
			ibs = integrated_brier_score(rsf, scaled_X_test, T_test, E_test)
			iterations_ibss.append(ibs)
			print('IBS: {:.4f}'.format(ibs), file=open(iterations_file_name, "a"))

			end_time2 = time.time()
			hours2, rem2 = divmod(end_time2-start_time2, 3600)
			minutes2, seconds2 = divmod(rem2, 60)
			print("\nFinished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2))
			print("Finished the iteration in {:0>2}:{:0>2}:{:05.2f}".format(int(hours2),int(minutes2),seconds2), file=open(iterations_file_name, "a"))

		iterations_vimps_dict = dict(iterations_vimps_dict)
		iterations_vimps_dict = {k: v / iterations_per_nr_of_trees for k, v in iterations_vimps_dict.items() if v != 0}
		# print("\niterations_vimps_dict: ", iterations_vimps_dict, file=open(iterations_file_name, "a"))
		list_from_dict1 = [(k, v) for k, v in iterations_vimps_dict.items()]
		vimps_list_ordered1 = sorted(list_from_dict1, key=itemgetter(1))
		descending_vimps_list1 = vimps_list_ordered1[::-1]
		print("\ndescending_vimps_list: ", descending_vimps_list1, file=open(iterations_file_name, "a"))

		average_iterations_cindexes_test = sum(iterations_cindexes_test) / len(iterations_cindexes_test)
		all_cindexes_test.append(average_iterations_cindexes_test)
		print('\nAverage c-index test: {:.4f}'.format(average_iterations_cindexes_test), file=open(all_results_file_name, "a"))
		std_cindexes_test = std_calc(iterations_cindexes_test, average_iterations_cindexes_test)
		all_std_cindexes_test.append(std_cindexes_test)
		print('Standard deviation: {:.4f}'.format(std_cindexes_test), file=open(all_results_file_name, "a"))

		average_iterations_cindexes_train = sum(iterations_cindexes_train) / len(iterations_cindexes_train)
		all_cindexes_train.append(average_iterations_cindexes_train)
		print('\nAverage c-index train: {:.4f}'.format(average_iterations_cindexes_train), file=open(all_results_file_name, "a"))
		std_cindexes_train = std_calc(iterations_cindexes_train, average_iterations_cindexes_train)
		all_std_cindexes_train.append(std_cindexes_train)
		print('Standard deviation: {:.4f}'.format(std_cindexes_train), file=open(all_results_file_name, "a"))

		average_iterations_ibs = sum(iterations_ibss) / len(iterations_ibss)
		all_ibss.append(average_iterations_ibs)
		print('\nAverage IBS: {:.4f}'.format(average_iterations_ibs), file=open(all_results_file_name, "a"))
		std_ibss = std_calc(iterations_ibss, average_iterations_ibs)
		all_std_ibss.append(std_ibss)
		print('Standard deviation: {:.4f}'.format(std_ibss), file=open(all_results_file_name, "a"))

		all_vimps_dict = Counter(all_vimps_dict)
		all_vimps_dict.update(Counter(iterations_vimps_dict))

		end_time3 = time.time()
		hours3, rem3 = divmod(end_time3-start_time3, 3600)
		minutes3, seconds3 = divmod(rem3, 60)
		print("\nFinished the iterations for min node size =", i, "in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3))
		print("\nFinished the iterations for min node size =", i, "in {:0>2}:{:0>2}:{:05.2f}".format(int(hours3),int(minutes3),seconds3), file=open(all_results_file_name, "a"))

	all_vimps_dict = dict(all_vimps_dict)
	all_vimps_dict = {k: v / len(min_node_sizes) for k, v in all_vimps_dict.items() if v != 0}
	# print("\nall_vimps_dict: ", all_vimps_dict, file=open(all_results_file_name, "a"))
	list_from_dict2 = [(k, v) for k, v in all_vimps_dict.items()]
	vimps_list_ordered2 = sorted(list_from_dict2, key=itemgetter(1))
	descending_vimps_list2 = vimps_list_ordered2[::-1]
	print("\ndescending_vimps_list: ", descending_vimps_list2, file=open(all_results_file_name, "a"))
	positive_vimps_genes = [item for item in descending_vimps_list2 if item[1] > 0]
	negative_vimps_genes = [item for item in descending_vimps_list2 if item[1] < 0]
	print("\nGenes and their positive vimps: ", positive_vimps_genes, file=open(all_results_file_name, "a"))
	list_of_genes_with_positive_vimps = [tup[0] for tup in positive_vimps_genes]
	list_of_genes_with_negative_vimps = [tup[0] for tup in negative_vimps_genes]
	print("\nList of genes with positive vimps: ", list_of_genes_with_positive_vimps, file=open(all_results_file_name, "a"))
	print("\nNumber of genes with positive vimps: ", len(list_of_genes_with_positive_vimps), file=open(all_results_file_name, "a"))
	print("\nNumber of genes with negative vimps: ", len(list_of_genes_with_negative_vimps), file=open(all_results_file_name, "a"))

	print("\n***********\n**SUMARY**\n***********", file=open(all_results_file_name, "a"))
	print("\nMin node sizes: ", min_node_sizes, file=open(all_results_file_name, "a"))
	print("\nAll c-indexes test: ", all_cindexes_test, file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_test, file=open(all_results_file_name, "a"))
	print("\nAll c-indexes train: ", all_cindexes_train, file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_cindexes_train, file=open(all_results_file_name, "a"))
	print("\nAll ibs: ", all_ibss, file=open(all_results_file_name, "a"))
	print("Standard deviation:", all_std_ibss, file=open(all_results_file_name, "a"))

	end_time1 = time.time()
	hours1, rem1 = divmod(end_time1-start_time1, 3600)
	minutes1, seconds1 = divmod(rem1, 60)
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1))
	print("\nFinished all in {:0>2}:{:0>2}:{:05.2f}".format(int(hours1),int(minutes1),seconds1), file=open(all_results_file_name, "a"))

	iterations_for_plot = list(range(1, iterations_per_nr_of_trees+1))
	plt.figure(figsize=(12, 6))
	line1, = plt.plot(iterations_for_plot, iterations_cindexes_train, color = (1, 0, 0.8), marker='D', label="Train set c-index")
	line2, = plt.plot(iterations_for_plot, iterations_cindexes_test, color = (0.2, 0.5, 0.9), marker='D', label="Test set c-index")
	plt.grid(linestyle = '-.')
	plt.xticks(iterations_for_plot)
	plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
	plt.ylabel("C-index")
	plt.xlabel("iterations executed")
	# plt.show()
	plt.savefig("RSF_best_model.png")

if __name__ == '__main__':

	all_data, genes, N = prepareData()

	print(N)

	all_data = all_data.reset_index()
	print(all_data)

	# Building training and testing sets
	index_train, index_test = train_test_split(range(N), test_size = 0.2, random_state=5)
	data_train = all_data.loc[index_train].reset_index(drop = True)
	data_test  = all_data.loc[index_test].reset_index(drop = True)

	# Creating the X, T and E input
	X_train, X_test = data_train[genes], data_test[genes]
	print(X_train.shape)
	print(X_test.shape)
	X_train = X_train.drop(columns=["AAA1", "AADACL4", "ACCSL", "ACTL7A", "ACTL9", "ACTRT1", "ACTRT2", "ADAD1", "ADAM30", "ADAM3A", "ADAM5P", "ADAM7", "ADIG", "AFM", "AKAP4", "AMELX", "AMELY", "APCS", "APOA4", "APOH", "AQP12A", "ARGFX", "ASB10", "ASB11", "ASB17", "ASB18", "ATOH1", "ATXN3L", "ATXN8OS", "AVP", "BANF2", "BCORL2", "BCYRN1", "BEYLA", "BHLHE23", "BIRC8", "BLID", "BMP10", "BMP15", "BPESC1", "BPY2", "BSPH1", "BSX", "C10orf120", "C10orf40", "C10orf53", "C10orf71", "C10orf96", "C11orf36", "C11orf40", "C11orf64", "C12orf12", "C13orf28", "C14orf165", "C14orf177", "C14orf23", "C14orf48", "C14orf70", "C15orf43", "C16orf78", "C17orf105", "C18orf20", "C18orf62", "C19orf41", "C19orf75", "C1orf100", "C1orf141", "C1orf146", "C1orf185", "C1orf68", "C20orf166", "C20orf185", "C20orf71", "C20orf79", "C21orf131", "C21orf54", "C21orf94", "C22orf33", "C22orf42", "C2orf27B", "C2orf53", "C2orf80", "C2orf83", "C3P1", "C3orf22", "C3orf51", "C3orf65", "C3orf77", "C3orf79", "C4orf11", "C4orf17", "C4orf35", "C4orf45", "C5orf48", "C5orf52", "C6orf146", "C7orf66", "C7orf72", "C8A", "C8orf22", "C8orf71", "C8orf74", "C9orf27", "C9orf57", "C9orf79", "CABP2", "CABP5", "CACNG2", "CCDC105", "CCL14-CCL15", "CD300LD", "CDH9", "CDX4", "CDY1B", "CDY1", "CDY2B", "CEACAM18", "CELA2A", "CELA2B", "CELA3A", "CELA3B", "CETN1", "CFHR2", "CGB1", "CHAT", "CHRND", "CLCA1", "CLDN17", "CLDN22", "CLRN2", "CNBD1", "CPN1", "CPXCR1", "CRYGA", "CRYGB", "CRYGC", "CRYGD", "CSH1", "CSH2", "CSHL1", "CSN1S2A", "CSPG4PY2", "CSRP3", "CST8", "CT45A2", "CT45A3", "CT45A4", "CT45A6", "CT47A10", "CT47A11", "CT47A1", "CT47A2", "CT47A6", "CT47A7", "CT47A9", "CT47B1", "CTRB1", "CTSL3", "CXorf51", "CXorf66", "CYLC1", "CYLC2", "CYMP", "CYP11B2", "CYorf15A", "CYorf15B", "DAOA", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DBX1", "DCAF8L1", "DEFA4",   "DEFA5", "DEFA6", "DEFB103B", "DEFB104A", "DEFB105A", "DEFB106A", "DEFB107A", "DEFB108B", "DEFB109P1B", "DEFB109P1", "DEFB110", "DEFB112", "DEFB113", "DEFB114", "DEFB115", "DEFB116", "DEFB118", "DEFB119", "DEFB121", "DEFB122", "DEFB123", "DEFB125", "DEFB127", "DEFB128", "DEFB129", "DEFB130", "DEFB131", "DEFB134", "DEFB135", "DEFB136", "DEFB4A", "DHRS7C", "DKFZP434H168", "DKFZP434K028", "DKFZp434L192", "DMRTB1", "DNAJB3", "DNAJB8", "DPPA2", "DPPA5", "DRGX", "DSCR10", "DSCR4", "DUSP21", "DUX4", "DUXA", "DYTN", "EIF1AY", "ELSPBP1", "EVX2", "FABP12", "FABP9", "FAM138B", "FAM138D", "FAM138E", "FAM194B", "FAM197Y2", "FAM24A", "FAM27B", "FAM27L", "FAM41AY1", "FAM47B", "FAM71B", "FAM74A1", "FAM74A3", "FAM74A4", "FAM75A6", "FAM92A3", "FAM99A", "FAM99B", "FAM9B", "FGF4", "FGF6",  "FKSG73", "FKSG83", "FLJ25328", "FLJ25363", "FLJ25758", "FLJ36000", "FLJ40434", "FLJ43859", "FLJ44082", "FLJ46321", "FLJ46361", "FOXB2", "FOXR2", "FSHB", "FSHR", "FTHL17", "FTLP10", "GABRA6",  "GABRR3", "GAGE12F", "GAGE12J", "GAGE13", "GAGE1", "GAGE2A", "GAGE2B", "GAGE2C",  "GAGE2D",  "GAGE2E",  "GAGE8", "GALNTL5", "GALP", "GAST", "GCG", "GC",  "GDEP",  "GDF2",  "GFRA4", "GJA10",   "GJA8", "GJD2", "GK2", "GLT6D1",  "GLYCAM1", "GML", "GOLGA2P3", "GOLGA6C", "GOLGA6D", "GOLGA8E", "GOLGA8F", "GOT1L1", "GPHB5", "GPR101", "GPR119", "GPR148",  "GPR149",  "GPR32", "GPRC6A", "GPX5", "GPX6", "GRXCR1",  "GSC2", "GSTTP1",  "GSTTP2",  "GSX1", "GSX2", "GUCA1C",  "GUCA2B",  "GUCY2GP", "GYPA", "GYPB", "H1FOO", "H2BFWT", "HBBP1", "HBII-52-24", "HBII-52-27", "HBII-52-28", "HBII-52-45", "HBII-52-46", "HBM", "HELT", "HES3", "HIGD1C", "HIST1H2AA", "HIST1H2BA", "HIST1H4G", "HMGB4", "HNRNPCL1", "HPVC1", "HSFX1", "HSFY2", "HSFYL1", "HTN1", "HTN3", "HTR3B",  "HTR3D", "HTR5A", "HYALP1", "IAPP", "IFNA10", "IFNA14", "IFNA16", "IFNA17",  "IFNA2", "IFNA4", "IFNA6", "IFNA7",   "IFNA8",   "IL17F",  "IL1F10", "IL1F6", "IL22", "IL31", "IL3", "IL9", "INS", "IQCF2", "IQCF3", "IQCF5",   "KCNK18",  "KIF2B",   "KIR3DP1", "KRT25", "KRT26", "KRT28", "KRT76", "KRTAP1-3", "KRTAP10-10", "KRTAP10-11", "KRTAP10-12", "KRTAP10-1", "KRTAP10-2", "KRTAP10-3", "KRTAP10-5", "KRTAP10-6", "KRTAP10-7", "KRTAP10-8", "KRTAP10-9", "KRTAP11-1", "KRTAP12-1", "KRTAP12-2", "KRTAP12-3", "KRTAP12-4", "KRTAP13-1", "KRTAP13-2", "KRTAP13-3", "KRTAP13-4", "KRTAP15-1", "KRTAP19-1", "KRTAP19-2", "KRTAP19-3", "KRTAP19-4", "KRTAP19-5", "KRTAP19-6", "KRTAP19-7", "KRTAP19-8", "KRTAP2-4", "KRTAP20-1", "KRTAP20-2", "KRTAP20-3", "KRTAP20-4", "KRTAP21-1", "KRTAP21-2", "KRTAP22-1", "KRTAP23-1", "KRTAP24-1", "KRTAP25-1", "KRTAP26-1", "KRTAP27-1", "KRTAP4-11", "KRTAP4-12", "KRTAP4-2",  "KRTAP4-3",  "KRTAP4-4", "KRTAP4-5", "KRTAP4-7", "KRTAP4-8", "KRTAP4-9", "KRTAP5-11", "KRTAP5-3",  "KRTAP5-5",  "KRTAP6-1", "KRTAP6-2",  "KRTAP6-3",  "KRTAP7-1",  "KRTAP8-1", "KRTAP9-2",  "KRTAP9-4",  "KRTAP9-8",  "KRTAP9-9", "LCE1A", "LCE1B", "LCE1D", "LCE1E", "LCE1F", "LCE2A", "LCE2B", "LCE2C", "LCE2D",   "LCE3A", "LCE3B",   "LCE3C",  "LCE3D",   "LCE3E",   "LCE4A",   "LCE6A",  "LCN15",   "LCN9", "LELP1",   "LEUTX",  "LGALS13", "LIM2", "LIPM", "LOC100128554", "LOC100129935", "LOC100169752", "LOC100287704", "LOC144742", "LOC144776", "LOC150185", "LOC151300", "LOC151658", "LOC154449", "LOC157627", "LOC200726", "LOC254312", "LOC255025", "LOC283914", "LOC284788", "LOC285194", "LOC285370", "LOC285375", "LOC285627", "LOC285735", "LOC286094", "LOC286135", "LOC286238", "LOC29034", "LOC338588", "LOC339568", "LOC340094", "LOC340357", "LOC348021", "LOC360030", "LOC388946", "LOC390858", "LOC402644", "LOC642929", "LOC643486", "LOC643923", "LOC643955", "LOC644145", "LOC646498", "LOC646813", "LOC650293", "LOC653544", "LOC653545", "LOC727677", "LOC727924", "LOC728410", "LOC729121", "LOC730811", "LOC731779", "LOC732275", "LRIT1",   "LRRC10", "LRRC30",  "LYZL6",   "MAGEB10", "MAGEB3", "MAS1", "MBD3L1", "MBD3L2", "MBD3L5", "MC2R", "MGC34034", "MLN", "MMD2", "MMP26", "MOGAT3", "MOS", "MRGPRG", "MRGPRX1", "MS4A12",  "MS4A13",  "MS4A5", "MT1B", "MT1IP",   "MT4",  "MYADML", "MYF5", "MYL10", "MYL1", "NAA11", "NCRNA00029", "NCRNA00051", "NCRNA00099", "NCRNA00111", "NCRNA00112", "NCRNA00114", "NCRNA00157", "NCRNA00159", "NCRNA00185", "NCRNA00207", "NCRNA00230B", "NEU2", "NEUROD4", "NEUROD6", "NEUROG1", "NF1P1", "NKX2-4", "NKX2-6", "NKX6-2", "NMS", "NOBOX", "NOX3", "NPBWR2", "NPSR1", "NPS",  "NPVF", "NXF2B", "NXF2", "OCM2", "ONECUT3", "OPALIN", "OPN1LW", "OPN1MW", "OPRM1", "OR10A2", "OR10A4", "OR10A5", "OR10A6", "OR10A7", "OR10AG1", "OR10C1", "OR10G3", "OR10G4", "OR10G7", "OR10G8",  "OR10G9", "OR10H2", "OR10H3", "OR10H4", "OR10J1", "OR10J3", "OR10J5", "OR10K1",  "OR10K2", "OR10P1", "OR10R2", "OR10S1", "OR10T2", "OR10W1", "OR10X1", "OR10Z1", "OR11A1", "OR11G2", "OR11H1", "OR11H4",  "OR11H6",  "OR11L1",  "OR12D2", "OR12D3",  "OR13C2",  "OR13C3",  "OR13C4", "OR13C5",  "OR13C8",  "OR13C9",  "OR13D1", "OR13F1", "OR13H1", "OR14A16", "OR14C36",  "OR14I1",  "OR14J1",  "OR1A1",   "OR1A2", "OR1D2", "OR1D4", "OR1E2", "OR1I1", "OR1L1", "OR1L4",  "OR1L6", "OR1M1", "OR1N2", "OR1S1", "OR1S2", "OR2A12", "OR2A2", "OR2A5", "OR2AG1", "OR2AK2", "OR2AT4",  "OR2B3",   "OR2D2",   "OR2D3", "OR2F1", "OR2F2", "OR2G2", "OR2G3", "OR2G6", "OR2H1",  "OR2J2", "OR2J3", "OR2L8", "OR2M1P", "OR2M2", "OR2M3", "OR2M5",  "OR2M7",  "OR2T10", "OR2T11", "OR2T12", "OR2T1", "OR2T27", "OR2T29", "OR2T2",   "OR2T34", "OR2T35", "OR2T3", "OR2T5",  "OR2T6",   "OR2V2",   "OR2W1", "OR2W5",   "OR2Y1",   "OR2Z1",   "OR3A3", "OR3A4", "OR4A15",  "OR4A16",  "OR4A47", "OR4A5",  "OR4B1",   "OR4C11",  "OR4C12", "OR4C13",  "OR4C15",  "OR4C16",  "OR4C3", "OR4C45",  "OR4C46",  "OR4D10",  "OR4D11", "OR4D1",   "OR4D2",   "OR4D5",   "OR4D6", "OR4D9",   "OR4E2",   "OR4F15",  "OR4F17", "OR4F21",  "OR4F4",   "OR4F5",   "OR4F6",  "OR4K13",  "OR4K14",  "OR4K17",  "OR4K1",  "OR4K2",   "OR4K5",   "OR4L1",   "OR4M1", "OR4M2",   "OR4N2",   "OR4N3P",  "OR4N5", "OR4P4",   "OR4Q3",   "OR4S1",   "OR4S2", "OR4X1",   "OR4X2",   "OR51A2",  "OR51A4", "OR51A7",  "OR51B2",  "OR51B6",  "OR51D1", "OR51F1", "OR51F2", "OR51G1", "OR51G2", "OR51I2",  "OR51L1",  "OR51M1",  "OR51S1", "OR51T1",  "OR51V1",  "OR52A1",  "OR52A4", "OR52A5",  "OR52B2",  "OR52B4",  "OR52E2", "OR52E8",  "OR52I2",  "OR52J3",  "OR52K1", "OR52L1",  "OR52M1",  "OR52R1", "OR56A1", "OR56A4",  "OR5A1",   "OR5A2",   "OR5AC2", "OR5AK2",  "OR5AN1",  "OR5AP2",  "OR5AR1", "OR5AS1",  "OR5B12",  "OR5B17",  "OR5B21", "OR5B2",   "OR5B3",   "OR5D13",  "OR5D14", "OR5D16",  "OR5D18",  "OR5F1",   "OR5H14", "OR5H15",  "OR5H1",   "OR5H2",   "OR5H6", "OR5I1",   "OR5J2",   "OR5K3",   "OR5K4", "OR5L1", "OR5L2", "OR5M10", "OR5M11", "OR5M1", "OR5M3", "OR5M8", "OR5M9", "OR5R1", "OR5T1", "OR5T2", "OR5T3", "OR5V1", "OR5W2",   "OR6A2",   "OR6B1", "OR6C1",   "OR6C2",   "OR6C3",   "OR6C4", "OR6C65", "OR6C68",  "OR6C6",   "OR6C70", "OR6C74",  "OR6C75",  "OR6C76",  "OR6K2",  "OR6K6",   "OR6M1",   "OR6N1",   "OR6N2", "OR6P1",   "OR6Q1",   "OR6S1",   "OR6T1",  "OR6V1",   "OR6X1",   "OR6Y1",   "OR7A10", "OR7A17",  "OR7C2",   "OR7D4",   "OR7E156P", "OR7E24",  "OR7E5P",  "OR7G1",   "OR7G2", "OR7G3",   "OR8B12",  "OR8B2",   "OR8B3", "OR8B4",   "OR8B8",   "OR8D1",   "OR8D2",  "OR8D4",   "OR8G2",   "OR8G5",   "OR8H1",  "OR8H2",   "OR8H3",   "OR8I2",   "OR8J1", "OR8J3",   "OR8K1",   "OR8K3",   "OR8K5", "OR8U1",   "OR9G4",   "OR9G9",   "OR9I1", "OR9K2",   "OR9Q1",   "OR9Q2",   "OSTN", "OTOL1",   "OTOP1",   "OTUD6A",  "P2RX6P", "PAGE1",   "PAGE3",   "PAGE4",   "PAR4", "PASD1",   "PATE1",   "PATE3",   "PAX4", "PCGEM1",  "PCNAP1",  "PDHA2",   "PDILT", "PER4", "PFN3", "PHOX2B", "PISRT1", "PLA2G2E", "PLAC1L", "PLGLA", "PLG", "PMCHL1", "PNLIPRP1",  "POM121L12", "POTEA", "POU3F4",  "PPBPL2",  "PPIAL4B", "PPIAL4E", "PPP1R3A", "PRAC", "PRAMEF10",  "PRAMEF11", "PRAMEF12",  "PRAMEF13",  "PRAMEF14",  "PRAMEF16", "PRAMEF17",  "PRAMEF18",  "PRAMEF1", "PRAMEF20", "PRAMEF22",  "PRAMEF2", "PRAMEF3", "PRAMEF4", "PRAMEF5", "PRAMEF6", "PRAMEF9", "PRB4", "PRDM14",  "PRDM9", "PRG3", "PRHOXNB", "PRM1", "PRM2", "PRM3", "PRNT", "PRO1768", "PRODH2",  "PRR20A",  "PRR20B", "PRR20C",  "PRR20D",  "PRR21",   "PRR23A", "PRR23B",  "PRR23C",  "PRSS38",  "PRSS48", "PRSS54",  "PRY2", "PSG10",   "PSG11", "PSG2", "PSG6", "PSG7", "PTH2", "PTH", "PWRN1", "PWRN2", "PYDC2", "RAB9BP1", "RAG2", "RBMXL3",  "RBMY1A1", "RBMY1A3P",  "RBMY1B",  "RBMY1E",  "RBMY1F", "RBMY1J",  "RBMY2EP", "RBMY2FP", "RBMY3AP", "REG1A",   "REG1B",   "REG1P",   "REG3A",  "REG3G",   "RESP18",  "RETNLB",  "RGS21", "RNASE11", "RNASE12", "RNASE8",  "RNASE9", "RNU4ATAC",  "RNU5E",   "RNU6ATAC",  "RNY4", "RNY5", "RPS4Y2",  "RXFP2", "S100A7L2", "SAA3P",   "SCARNA11",  "SCARNA14",  "SCARNA15", "SCARNA18",  "SCARNA1", "SCARNA20",  "SCARNA21", "SCARNA22",  "SCARNA23",  "SCARNA27",  "SCARNA3", "SCARNA4", "SCARNA8", "SCGB1D1", "SCGB1D4", "SDC4P",   "SEMG1", "SERPINA13", "SERPINA7", "SERPINB11", "SFTA3",   "SI",   "SLC10A2", "SLC17A1", "SLC17A2", "SLC17A6", "SLC18A3", "SLC22A6", "SLC2A2",  "SLC6A18", "SLC6A5", "SMCP", "SMEK3P", "SMR3A", "SNAR-A13", "SNAR-A2", "SNAR-A3", "SNAR-A4", "SNAR-B2", "SNAR-C2", "SNAR-C3", "SNAR-C4", "SNAR-D", "SNAR-E", "SNAR-F", "SNAR-G1", "SNAR-G2", "SNAR-H", "SNAR-I", "SNORA11B", "SNORA11C", "SNORA11D", "SNORA11E", "SNORA13", "SNORA14A", "SNORA15", "SNORA16A",  "SNORA16B",  "SNORA19", "SNORA1",  "SNORA25", "SNORA26", "SNORA28",  "SNORA29", "SNORA2A", "SNORA2B", "SNORA30", "SNORA32", "SNORA35", "SNORA36A",  "SNORA36B", "SNORA36C", "SNORA37", "SNORA38B",  "SNORA38", "SNORA41", "SNORA42", "SNORA44", "SNORA46", "SNORA47", "SNORA4",  "SNORA50", "SNORA51",  "SNORA54", "SNORA55", "SNORA56", "SNORA58", "SNORA5B", "SNORA5C", "SNORA66", "SNORA69", "SNORA70B", "SNORA70C", "SNORA71C", "SNORA71D", "SNORA75", "SNORA77", "SNORA78", "SNORA79", "SNORA80", "SNORA84", "SNORD102",  "SNORD103A", "SNORD104", "SNORD105B", "SNORD105",  "SNORD109B", "SNORD110",  "SNORD111B", "SNORD111",  "SNORD113-1", "SNORD113-2", "SNORD113-4", "SNORD113-5", "SNORD113-6", "SNORD113-7", "SNORD113-9", "SNORD114-10", "SNORD114-11", "SNORD114-12", "SNORD114-13", "SNORD114-14", "SNORD114-15", "SNORD114-16", "SNORD114-17", "SNORD114-18", "SNORD114-19", "SNORD114-1", "SNORD114-20", "SNORD114-21", "SNORD114-22", "SNORD114-23",  "SNORD114-24",  "SNORD114-25",  "SNORD114-26", "SNORD114-27",  "SNORD114-28",  "SNORD114-29",  "SNORD114-2", "SNORD114-30", "SNORD114-31", "SNORD114-3", "SNORD114-4", "SNORD114-5", "SNORD114-6", "SNORD114-7", "SNORD114-8", "SNORD114-9", "SNORD115-10",  "SNORD115-11",  "SNORD115-13", "SNORD115-14",  "SNORD115-16",  "SNORD115-17",  "SNORD115-1", "SNORD115-20", "SNORD115-22", "SNORD115-25", "SNORD115-2", "SNORD115-30", "SNORD115-31", "SNORD115-32", "SNORD115-33", "SNORD115-35",  "SNORD115-37",  "SNORD115-38",  "SNORD115-39", "SNORD115-3", "SNORD115-40",  "SNORD115-41", "SNORD115-44", "SNORD115-48", "SNORD115-4", "SNORD115-5", "SNORD115-6", "SNORD115-7", "SNORD115-8", "SNORD115-9", "SNORD116-10", "SNORD116-11", "SNORD116-12", "SNORD116-13", "SNORD116-14", "SNORD116-15", "SNORD116-16", "SNORD116-18",  "SNORD116-19", "SNORD116-1", "SNORD116-22", "SNORD116-23", "SNORD116-24", "SNORD116-25", "SNORD116-26", "SNORD116-27", "SNORD116-29", "SNORD116-2", "SNORD116-3", "SNORD116-5", "SNORD116-8", "SNORD117",  "SNORD119",  "SNORD11B",  "SNORD11", "SNORD121A", "SNORD121B", "SNORD123",  "SNORD124", "SNORD125",  "SNORD126",  "SNORD127",  "SNORD12B", "SNORD12C",  "SNORD12", "SNORD16", "SNORD18A", "SNORD18B",  "SNORD18C",  "SNORD19B",  "SNORD19", "SNORD1A", "SNORD1B", "SNORD20", "SNORD21", "SNORD23", "SNORD24", "SNORD25", "SNORD26", "SNORD27", "SNORD28", "SNORD29", "SNORD2", "SNORD30", "SNORD31", "SNORD32A", "SNORD32B", "SNORD33", "SNORD34", "SNORD35A",  "SNORD35B", "SNORD36A", "SNORD36B",  "SNORD36C",  "SNORD37", "SNORD38A",  "SNORD38B",  "SNORD41", "SNORD42A", "SNORD42B", "SNORD43", "SNORD44", "SNORD45A", "SNORD45B",  "SNORD45C",  "SNORD46", "SNORD47", "SNORD48", "SNORD49A",  "SNORD49B",  "SNORD4A", "SNORD4B", "SNORD50A",  "SNORD50B",  "SNORD51",  "SNORD52", "SNORD53", "SNORD54", "SNORD55", "SNORD56B", "SNORD56", "SNORD57", "SNORD58A", "SNORD58C", "SNORD59A", "SNORD59B", "SNORD5", "SNORD60", "SNORD61", "SNORD62A",  "SNORD63", "SNORD65", "SNORD66", "SNORD67", "SNORD68", "SNORD69", "SNORD6",  "SNORD70", "SNORD71", "SNORD72", "SNORD73A",  "SNORD74", "SNORD75", "SNORD76", "SNORD77", "SNORD78", "SNORD79", "SNORD7",  "SNORD80", "SNORD81", "SNORD82", "SNORD84", "SNORD85", "SNORD86", "SNORD87", "SNORD88A",  "SNORD88B",  "SNORD88C",  "SNORD89", "SNORD8",  "SNORD90", "SNORD91A",  "SNORD91B", "SNORD92", "SNORD93", "SNORD95", "SNORD96A", "SNORD96B",  "SNORD98", "SNORD99", "SNORD9", "SNTG1",   "SOX14",   "SOX1", "SPACA1", "SPAG11A", "SPAG11B", "SPAM1",  "SPANXE", "SPANXN1", "SPANXN2", "SPANXN4", "SPANXN5", "SPATA16", "SPATA19", "SPINK14", "SPINK7", "SPINK9",  "SPINT3",  "SPINT4",  "SPO11", "SPP2", "SPRR2B",  "SPRR2C",  "SPRR4",  "SPZ1", "SRY",  "SSTR4",   "SSX3", "SSX7", "SSX8", "ST8SIA3", "STATH", "SULT1C3", "SULT2A1", "SULT6B1", "SUN5", "SYCN", "TAAR1",  "TAAR2",  "TAAR5", "TAAR6",   "TAAR8",   "TAAR9",   "TARM1", "TAS2R16", "TAS2R39", "TAS2R41", "TAS2R7", "TBC1D21", "TBC1D28", "TBL1Y", "TBPL2", "TCP10L2", "TCP10",   "TECRL",   "TECTB",  "TEX13A",  "TEX28",   "TFAP2D",  "TGIF2LX", "TGIF2LY", "TGM6", "TINAG", "TMCO5A", "TMEM174", "TMEM202", "TMEM207", "TMEM225", "TMEM30C", "TMEM8C",  "TMEM95",  "TMIGD1", "TMPRSS11BNL", "TMPRSS11B", "TMSB4Y",  "TNP1", "TNP2", "TPD52L3", "TPRX1",  "TREML2P1", "TRIM42",  "TRIM48",  "TRIM49L", "TRIM53", "TRIM60",  "TRIM64",  "TRIM77",  "TRPC7",  "TRYX3",   "TSPY1",   "TSPY2",   "TSPY3", "TSPY4",   "TSSK2",   "TTLL8",   "TTTY10", "TTTY11",  "TTTY12",  "TTTY13",  "TTTY14", "TTTY15",  "TTTY16",  "TTTY17A", "TTTY17B", "TTTY18",  "TTTY19",  "TTTY1B",  "TTTY20", "TTTY21",  "TTTY22",  "TTTY23",  "TTTY2", "TTTY3B",  "TTTY4C",  "TTTY5",   "TTTY6B", "TTTY6",   "TTTY7",   "TTTY8",   "TTTY9B", "TXNDC8",  "UBE2DNL", "UBE2U",   "UBTFL1", "UCMA", "UCN3", "UGT1A3",  "UGT1A4", "UGT1A5", "UGT1A8",  "UNCX", "USP17L6P", "USP17", "USP29", "USP9Y", "UTY", "VCX2", "VCY", "VENTXP1", "VENTXP7", "VN1R2", "VN1R4", "VPREB1", "VSTM2B", "VTRNA1-1",  "VTRNA1-2",  "VTRNA1-3",  "VWC2L",  "WFDC11",  "WFDC8",   "WFDC9",   "XAGE5", "XKRY2", "XKRY", "ZCCHC13", "ZCCHC16", "ZIM3", "ZNF479",  "ZNF645",  "ZNF679", "ZNF705D", "ZNF735",  "ZNRF4", "ZSWIM2"])
	X_test = X_test.drop(columns=["AAA1", "AADACL4", "ACCSL", "ACTL7A", "ACTL9", "ACTRT1", "ACTRT2", "ADAD1", "ADAM30", "ADAM3A", "ADAM5P", "ADAM7", "ADIG", "AFM", "AKAP4", "AMELX", "AMELY", "APCS", "APOA4", "APOH", "AQP12A", "ARGFX", "ASB10", "ASB11", "ASB17", "ASB18", "ATOH1", "ATXN3L", "ATXN8OS", "AVP", "BANF2", "BCORL2", "BCYRN1", "BEYLA", "BHLHE23", "BIRC8", "BLID", "BMP10", "BMP15", "BPESC1", "BPY2", "BSPH1", "BSX", "C10orf120", "C10orf40", "C10orf53", "C10orf71", "C10orf96", "C11orf36", "C11orf40", "C11orf64", "C12orf12", "C13orf28", "C14orf165", "C14orf177", "C14orf23", "C14orf48", "C14orf70", "C15orf43", "C16orf78", "C17orf105", "C18orf20", "C18orf62", "C19orf41", "C19orf75", "C1orf100", "C1orf141", "C1orf146", "C1orf185", "C1orf68", "C20orf166", "C20orf185", "C20orf71", "C20orf79", "C21orf131", "C21orf54", "C21orf94", "C22orf33", "C22orf42", "C2orf27B", "C2orf53", "C2orf80", "C2orf83", "C3P1", "C3orf22", "C3orf51", "C3orf65", "C3orf77", "C3orf79", "C4orf11", "C4orf17", "C4orf35", "C4orf45", "C5orf48", "C5orf52", "C6orf146", "C7orf66", "C7orf72", "C8A", "C8orf22", "C8orf71", "C8orf74", "C9orf27", "C9orf57", "C9orf79", "CABP2", "CABP5", "CACNG2", "CCDC105", "CCL14-CCL15", "CD300LD", "CDH9", "CDX4", "CDY1B", "CDY1", "CDY2B", "CEACAM18", "CELA2A", "CELA2B", "CELA3A", "CELA3B", "CETN1", "CFHR2", "CGB1", "CHAT", "CHRND", "CLCA1", "CLDN17", "CLDN22", "CLRN2", "CNBD1", "CPN1", "CPXCR1", "CRYGA", "CRYGB", "CRYGC", "CRYGD", "CSH1", "CSH2", "CSHL1", "CSN1S2A", "CSPG4PY2", "CSRP3", "CST8", "CT45A2", "CT45A3", "CT45A4", "CT45A6", "CT47A10", "CT47A11", "CT47A1", "CT47A2", "CT47A6", "CT47A7", "CT47A9", "CT47B1", "CTRB1", "CTSL3", "CXorf51", "CXorf66", "CYLC1", "CYLC2", "CYMP", "CYP11B2", "CYorf15A", "CYorf15B", "DAOA", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DBX1", "DCAF8L1", "DEFA4",   "DEFA5", "DEFA6", "DEFB103B", "DEFB104A", "DEFB105A", "DEFB106A", "DEFB107A", "DEFB108B", "DEFB109P1B", "DEFB109P1", "DEFB110", "DEFB112", "DEFB113", "DEFB114", "DEFB115", "DEFB116", "DEFB118", "DEFB119", "DEFB121", "DEFB122", "DEFB123", "DEFB125", "DEFB127", "DEFB128", "DEFB129", "DEFB130", "DEFB131", "DEFB134", "DEFB135", "DEFB136", "DEFB4A", "DHRS7C", "DKFZP434H168", "DKFZP434K028", "DKFZp434L192", "DMRTB1", "DNAJB3", "DNAJB8", "DPPA2", "DPPA5", "DRGX", "DSCR10", "DSCR4", "DUSP21", "DUX4", "DUXA", "DYTN", "EIF1AY", "ELSPBP1", "EVX2", "FABP12", "FABP9", "FAM138B", "FAM138D", "FAM138E", "FAM194B", "FAM197Y2", "FAM24A", "FAM27B", "FAM27L", "FAM41AY1", "FAM47B", "FAM71B", "FAM74A1", "FAM74A3", "FAM74A4", "FAM75A6", "FAM92A3", "FAM99A", "FAM99B", "FAM9B", "FGF4", "FGF6",  "FKSG73", "FKSG83", "FLJ25328", "FLJ25363", "FLJ25758", "FLJ36000", "FLJ40434", "FLJ43859", "FLJ44082", "FLJ46321", "FLJ46361", "FOXB2", "FOXR2", "FSHB", "FSHR", "FTHL17", "FTLP10", "GABRA6",  "GABRR3", "GAGE12F", "GAGE12J", "GAGE13", "GAGE1", "GAGE2A", "GAGE2B", "GAGE2C",  "GAGE2D",  "GAGE2E",  "GAGE8", "GALNTL5", "GALP", "GAST", "GCG", "GC",  "GDEP",  "GDF2",  "GFRA4", "GJA10",   "GJA8", "GJD2", "GK2", "GLT6D1",  "GLYCAM1", "GML", "GOLGA2P3", "GOLGA6C", "GOLGA6D", "GOLGA8E", "GOLGA8F", "GOT1L1", "GPHB5", "GPR101", "GPR119", "GPR148",  "GPR149",  "GPR32", "GPRC6A", "GPX5", "GPX6", "GRXCR1",  "GSC2", "GSTTP1",  "GSTTP2",  "GSX1", "GSX2", "GUCA1C",  "GUCA2B",  "GUCY2GP", "GYPA", "GYPB", "H1FOO", "H2BFWT", "HBBP1", "HBII-52-24", "HBII-52-27", "HBII-52-28", "HBII-52-45", "HBII-52-46", "HBM", "HELT", "HES3", "HIGD1C", "HIST1H2AA", "HIST1H2BA", "HIST1H4G", "HMGB4", "HNRNPCL1", "HPVC1", "HSFX1", "HSFY2", "HSFYL1", "HTN1", "HTN3", "HTR3B",  "HTR3D", "HTR5A", "HYALP1", "IAPP", "IFNA10", "IFNA14", "IFNA16", "IFNA17",  "IFNA2", "IFNA4", "IFNA6", "IFNA7",   "IFNA8",   "IL17F",  "IL1F10", "IL1F6", "IL22", "IL31", "IL3", "IL9", "INS", "IQCF2", "IQCF3", "IQCF5",   "KCNK18",  "KIF2B",   "KIR3DP1", "KRT25", "KRT26", "KRT28", "KRT76", "KRTAP1-3", "KRTAP10-10", "KRTAP10-11", "KRTAP10-12", "KRTAP10-1", "KRTAP10-2", "KRTAP10-3", "KRTAP10-5", "KRTAP10-6", "KRTAP10-7", "KRTAP10-8", "KRTAP10-9", "KRTAP11-1", "KRTAP12-1", "KRTAP12-2", "KRTAP12-3", "KRTAP12-4", "KRTAP13-1", "KRTAP13-2", "KRTAP13-3", "KRTAP13-4", "KRTAP15-1", "KRTAP19-1", "KRTAP19-2", "KRTAP19-3", "KRTAP19-4", "KRTAP19-5", "KRTAP19-6", "KRTAP19-7", "KRTAP19-8", "KRTAP2-4", "KRTAP20-1", "KRTAP20-2", "KRTAP20-3", "KRTAP20-4", "KRTAP21-1", "KRTAP21-2", "KRTAP22-1", "KRTAP23-1", "KRTAP24-1", "KRTAP25-1", "KRTAP26-1", "KRTAP27-1", "KRTAP4-11", "KRTAP4-12", "KRTAP4-2",  "KRTAP4-3",  "KRTAP4-4", "KRTAP4-5", "KRTAP4-7", "KRTAP4-8", "KRTAP4-9", "KRTAP5-11", "KRTAP5-3",  "KRTAP5-5",  "KRTAP6-1", "KRTAP6-2",  "KRTAP6-3",  "KRTAP7-1",  "KRTAP8-1", "KRTAP9-2",  "KRTAP9-4",  "KRTAP9-8",  "KRTAP9-9", "LCE1A", "LCE1B", "LCE1D", "LCE1E", "LCE1F", "LCE2A", "LCE2B", "LCE2C", "LCE2D",   "LCE3A", "LCE3B",   "LCE3C",  "LCE3D",   "LCE3E",   "LCE4A",   "LCE6A",  "LCN15",   "LCN9", "LELP1",   "LEUTX",  "LGALS13", "LIM2", "LIPM", "LOC100128554", "LOC100129935", "LOC100169752", "LOC100287704", "LOC144742", "LOC144776", "LOC150185", "LOC151300", "LOC151658", "LOC154449", "LOC157627", "LOC200726", "LOC254312", "LOC255025", "LOC283914", "LOC284788", "LOC285194", "LOC285370", "LOC285375", "LOC285627", "LOC285735", "LOC286094", "LOC286135", "LOC286238", "LOC29034", "LOC338588", "LOC339568", "LOC340094", "LOC340357", "LOC348021", "LOC360030", "LOC388946", "LOC390858", "LOC402644", "LOC642929", "LOC643486", "LOC643923", "LOC643955", "LOC644145", "LOC646498", "LOC646813", "LOC650293", "LOC653544", "LOC653545", "LOC727677", "LOC727924", "LOC728410", "LOC729121", "LOC730811", "LOC731779", "LOC732275", "LRIT1",   "LRRC10", "LRRC30",  "LYZL6",   "MAGEB10", "MAGEB3", "MAS1", "MBD3L1", "MBD3L2", "MBD3L5", "MC2R", "MGC34034", "MLN", "MMD2", "MMP26", "MOGAT3", "MOS", "MRGPRG", "MRGPRX1", "MS4A12",  "MS4A13",  "MS4A5", "MT1B", "MT1IP",   "MT4",  "MYADML", "MYF5", "MYL10", "MYL1", "NAA11", "NCRNA00029", "NCRNA00051", "NCRNA00099", "NCRNA00111", "NCRNA00112", "NCRNA00114", "NCRNA00157", "NCRNA00159", "NCRNA00185", "NCRNA00207", "NCRNA00230B", "NEU2", "NEUROD4", "NEUROD6", "NEUROG1", "NF1P1", "NKX2-4", "NKX2-6", "NKX6-2", "NMS", "NOBOX", "NOX3", "NPBWR2", "NPSR1", "NPS",  "NPVF", "NXF2B", "NXF2", "OCM2", "ONECUT3", "OPALIN", "OPN1LW", "OPN1MW", "OPRM1", "OR10A2", "OR10A4", "OR10A5", "OR10A6", "OR10A7", "OR10AG1", "OR10C1", "OR10G3", "OR10G4", "OR10G7", "OR10G8",  "OR10G9", "OR10H2", "OR10H3", "OR10H4", "OR10J1", "OR10J3", "OR10J5", "OR10K1",  "OR10K2", "OR10P1", "OR10R2", "OR10S1", "OR10T2", "OR10W1", "OR10X1", "OR10Z1", "OR11A1", "OR11G2", "OR11H1", "OR11H4",  "OR11H6",  "OR11L1",  "OR12D2", "OR12D3",  "OR13C2",  "OR13C3",  "OR13C4", "OR13C5",  "OR13C8",  "OR13C9",  "OR13D1", "OR13F1", "OR13H1", "OR14A16", "OR14C36",  "OR14I1",  "OR14J1",  "OR1A1",   "OR1A2", "OR1D2", "OR1D4", "OR1E2", "OR1I1", "OR1L1", "OR1L4",  "OR1L6", "OR1M1", "OR1N2", "OR1S1", "OR1S2", "OR2A12", "OR2A2", "OR2A5", "OR2AG1", "OR2AK2", "OR2AT4",  "OR2B3",   "OR2D2",   "OR2D3", "OR2F1", "OR2F2", "OR2G2", "OR2G3", "OR2G6", "OR2H1",  "OR2J2", "OR2J3", "OR2L8", "OR2M1P", "OR2M2", "OR2M3", "OR2M5",  "OR2M7",  "OR2T10", "OR2T11", "OR2T12", "OR2T1", "OR2T27", "OR2T29", "OR2T2",   "OR2T34", "OR2T35", "OR2T3", "OR2T5",  "OR2T6",   "OR2V2",   "OR2W1", "OR2W5",   "OR2Y1",   "OR2Z1",   "OR3A3", "OR3A4", "OR4A15",  "OR4A16",  "OR4A47", "OR4A5",  "OR4B1",   "OR4C11",  "OR4C12", "OR4C13",  "OR4C15",  "OR4C16",  "OR4C3", "OR4C45",  "OR4C46",  "OR4D10",  "OR4D11", "OR4D1",   "OR4D2",   "OR4D5",   "OR4D6", "OR4D9",   "OR4E2",   "OR4F15",  "OR4F17", "OR4F21",  "OR4F4",   "OR4F5",   "OR4F6",  "OR4K13",  "OR4K14",  "OR4K17",  "OR4K1",  "OR4K2",   "OR4K5",   "OR4L1",   "OR4M1", "OR4M2",   "OR4N2",   "OR4N3P",  "OR4N5", "OR4P4",   "OR4Q3",   "OR4S1",   "OR4S2", "OR4X1",   "OR4X2",   "OR51A2",  "OR51A4", "OR51A7",  "OR51B2",  "OR51B6",  "OR51D1", "OR51F1", "OR51F2", "OR51G1", "OR51G2", "OR51I2",  "OR51L1",  "OR51M1",  "OR51S1", "OR51T1",  "OR51V1",  "OR52A1",  "OR52A4", "OR52A5",  "OR52B2",  "OR52B4",  "OR52E2", "OR52E8",  "OR52I2",  "OR52J3",  "OR52K1", "OR52L1",  "OR52M1",  "OR52R1", "OR56A1", "OR56A4",  "OR5A1",   "OR5A2",   "OR5AC2", "OR5AK2",  "OR5AN1",  "OR5AP2",  "OR5AR1", "OR5AS1",  "OR5B12",  "OR5B17",  "OR5B21", "OR5B2",   "OR5B3",   "OR5D13",  "OR5D14", "OR5D16",  "OR5D18",  "OR5F1",   "OR5H14", "OR5H15",  "OR5H1",   "OR5H2",   "OR5H6", "OR5I1",   "OR5J2",   "OR5K3",   "OR5K4", "OR5L1", "OR5L2", "OR5M10", "OR5M11", "OR5M1", "OR5M3", "OR5M8", "OR5M9", "OR5R1", "OR5T1", "OR5T2", "OR5T3", "OR5V1", "OR5W2",   "OR6A2",   "OR6B1", "OR6C1",   "OR6C2",   "OR6C3",   "OR6C4", "OR6C65", "OR6C68",  "OR6C6",   "OR6C70", "OR6C74",  "OR6C75",  "OR6C76",  "OR6K2",  "OR6K6",   "OR6M1",   "OR6N1",   "OR6N2", "OR6P1",   "OR6Q1",   "OR6S1",   "OR6T1",  "OR6V1",   "OR6X1",   "OR6Y1",   "OR7A10", "OR7A17",  "OR7C2",   "OR7D4",   "OR7E156P", "OR7E24",  "OR7E5P",  "OR7G1",   "OR7G2", "OR7G3",   "OR8B12",  "OR8B2",   "OR8B3", "OR8B4",   "OR8B8",   "OR8D1",   "OR8D2",  "OR8D4",   "OR8G2",   "OR8G5",   "OR8H1",  "OR8H2",   "OR8H3",   "OR8I2",   "OR8J1", "OR8J3",   "OR8K1",   "OR8K3",   "OR8K5", "OR8U1",   "OR9G4",   "OR9G9",   "OR9I1", "OR9K2",   "OR9Q1",   "OR9Q2",   "OSTN", "OTOL1",   "OTOP1",   "OTUD6A",  "P2RX6P", "PAGE1",   "PAGE3",   "PAGE4",   "PAR4", "PASD1",   "PATE1",   "PATE3",   "PAX4", "PCGEM1",  "PCNAP1",  "PDHA2",   "PDILT", "PER4", "PFN3", "PHOX2B", "PISRT1", "PLA2G2E", "PLAC1L", "PLGLA", "PLG", "PMCHL1", "PNLIPRP1",  "POM121L12", "POTEA", "POU3F4",  "PPBPL2",  "PPIAL4B", "PPIAL4E", "PPP1R3A", "PRAC", "PRAMEF10",  "PRAMEF11", "PRAMEF12",  "PRAMEF13",  "PRAMEF14",  "PRAMEF16", "PRAMEF17",  "PRAMEF18",  "PRAMEF1", "PRAMEF20", "PRAMEF22",  "PRAMEF2", "PRAMEF3", "PRAMEF4", "PRAMEF5", "PRAMEF6", "PRAMEF9", "PRB4", "PRDM14",  "PRDM9", "PRG3", "PRHOXNB", "PRM1", "PRM2", "PRM3", "PRNT", "PRO1768", "PRODH2",  "PRR20A",  "PRR20B", "PRR20C",  "PRR20D",  "PRR21",   "PRR23A", "PRR23B",  "PRR23C",  "PRSS38",  "PRSS48", "PRSS54",  "PRY2", "PSG10",   "PSG11", "PSG2", "PSG6", "PSG7", "PTH2", "PTH", "PWRN1", "PWRN2", "PYDC2", "RAB9BP1", "RAG2", "RBMXL3",  "RBMY1A1", "RBMY1A3P",  "RBMY1B",  "RBMY1E",  "RBMY1F", "RBMY1J",  "RBMY2EP", "RBMY2FP", "RBMY3AP", "REG1A",   "REG1B",   "REG1P",   "REG3A",  "REG3G",   "RESP18",  "RETNLB",  "RGS21", "RNASE11", "RNASE12", "RNASE8",  "RNASE9", "RNU4ATAC",  "RNU5E",   "RNU6ATAC",  "RNY4", "RNY5", "RPS4Y2",  "RXFP2", "S100A7L2", "SAA3P",   "SCARNA11",  "SCARNA14",  "SCARNA15", "SCARNA18",  "SCARNA1", "SCARNA20",  "SCARNA21", "SCARNA22",  "SCARNA23",  "SCARNA27",  "SCARNA3", "SCARNA4", "SCARNA8", "SCGB1D1", "SCGB1D4", "SDC4P",   "SEMG1", "SERPINA13", "SERPINA7", "SERPINB11", "SFTA3",   "SI",   "SLC10A2", "SLC17A1", "SLC17A2", "SLC17A6", "SLC18A3", "SLC22A6", "SLC2A2",  "SLC6A18", "SLC6A5", "SMCP", "SMEK3P", "SMR3A", "SNAR-A13", "SNAR-A2", "SNAR-A3", "SNAR-A4", "SNAR-B2", "SNAR-C2", "SNAR-C3", "SNAR-C4", "SNAR-D", "SNAR-E", "SNAR-F", "SNAR-G1", "SNAR-G2", "SNAR-H", "SNAR-I", "SNORA11B", "SNORA11C", "SNORA11D", "SNORA11E", "SNORA13", "SNORA14A", "SNORA15", "SNORA16A",  "SNORA16B",  "SNORA19", "SNORA1",  "SNORA25", "SNORA26", "SNORA28",  "SNORA29", "SNORA2A", "SNORA2B", "SNORA30", "SNORA32", "SNORA35", "SNORA36A",  "SNORA36B", "SNORA36C", "SNORA37", "SNORA38B",  "SNORA38", "SNORA41", "SNORA42", "SNORA44", "SNORA46", "SNORA47", "SNORA4",  "SNORA50", "SNORA51",  "SNORA54", "SNORA55", "SNORA56", "SNORA58", "SNORA5B", "SNORA5C", "SNORA66", "SNORA69", "SNORA70B", "SNORA70C", "SNORA71C", "SNORA71D", "SNORA75", "SNORA77", "SNORA78", "SNORA79", "SNORA80", "SNORA84", "SNORD102",  "SNORD103A", "SNORD104", "SNORD105B", "SNORD105",  "SNORD109B", "SNORD110",  "SNORD111B", "SNORD111",  "SNORD113-1", "SNORD113-2", "SNORD113-4", "SNORD113-5", "SNORD113-6", "SNORD113-7", "SNORD113-9", "SNORD114-10", "SNORD114-11", "SNORD114-12", "SNORD114-13", "SNORD114-14", "SNORD114-15", "SNORD114-16", "SNORD114-17", "SNORD114-18", "SNORD114-19", "SNORD114-1", "SNORD114-20", "SNORD114-21", "SNORD114-22", "SNORD114-23",  "SNORD114-24",  "SNORD114-25",  "SNORD114-26", "SNORD114-27",  "SNORD114-28",  "SNORD114-29",  "SNORD114-2", "SNORD114-30", "SNORD114-31", "SNORD114-3", "SNORD114-4", "SNORD114-5", "SNORD114-6", "SNORD114-7", "SNORD114-8", "SNORD114-9", "SNORD115-10",  "SNORD115-11",  "SNORD115-13", "SNORD115-14",  "SNORD115-16",  "SNORD115-17",  "SNORD115-1", "SNORD115-20", "SNORD115-22", "SNORD115-25", "SNORD115-2", "SNORD115-30", "SNORD115-31", "SNORD115-32", "SNORD115-33", "SNORD115-35",  "SNORD115-37",  "SNORD115-38",  "SNORD115-39", "SNORD115-3", "SNORD115-40",  "SNORD115-41", "SNORD115-44", "SNORD115-48", "SNORD115-4", "SNORD115-5", "SNORD115-6", "SNORD115-7", "SNORD115-8", "SNORD115-9", "SNORD116-10", "SNORD116-11", "SNORD116-12", "SNORD116-13", "SNORD116-14", "SNORD116-15", "SNORD116-16", "SNORD116-18",  "SNORD116-19", "SNORD116-1", "SNORD116-22", "SNORD116-23", "SNORD116-24", "SNORD116-25", "SNORD116-26", "SNORD116-27", "SNORD116-29", "SNORD116-2", "SNORD116-3", "SNORD116-5", "SNORD116-8", "SNORD117",  "SNORD119",  "SNORD11B",  "SNORD11", "SNORD121A", "SNORD121B", "SNORD123",  "SNORD124", "SNORD125",  "SNORD126",  "SNORD127",  "SNORD12B", "SNORD12C",  "SNORD12", "SNORD16", "SNORD18A", "SNORD18B",  "SNORD18C",  "SNORD19B",  "SNORD19", "SNORD1A", "SNORD1B", "SNORD20", "SNORD21", "SNORD23", "SNORD24", "SNORD25", "SNORD26", "SNORD27", "SNORD28", "SNORD29", "SNORD2", "SNORD30", "SNORD31", "SNORD32A", "SNORD32B", "SNORD33", "SNORD34", "SNORD35A",  "SNORD35B", "SNORD36A", "SNORD36B",  "SNORD36C",  "SNORD37", "SNORD38A",  "SNORD38B",  "SNORD41", "SNORD42A", "SNORD42B", "SNORD43", "SNORD44", "SNORD45A", "SNORD45B",  "SNORD45C",  "SNORD46", "SNORD47", "SNORD48", "SNORD49A",  "SNORD49B",  "SNORD4A", "SNORD4B", "SNORD50A",  "SNORD50B",  "SNORD51",  "SNORD52", "SNORD53", "SNORD54", "SNORD55", "SNORD56B", "SNORD56", "SNORD57", "SNORD58A", "SNORD58C", "SNORD59A", "SNORD59B", "SNORD5", "SNORD60", "SNORD61", "SNORD62A",  "SNORD63", "SNORD65", "SNORD66", "SNORD67", "SNORD68", "SNORD69", "SNORD6",  "SNORD70", "SNORD71", "SNORD72", "SNORD73A",  "SNORD74", "SNORD75", "SNORD76", "SNORD77", "SNORD78", "SNORD79", "SNORD7",  "SNORD80", "SNORD81", "SNORD82", "SNORD84", "SNORD85", "SNORD86", "SNORD87", "SNORD88A",  "SNORD88B",  "SNORD88C",  "SNORD89", "SNORD8",  "SNORD90", "SNORD91A",  "SNORD91B", "SNORD92", "SNORD93", "SNORD95", "SNORD96A", "SNORD96B",  "SNORD98", "SNORD99", "SNORD9", "SNTG1",   "SOX14",   "SOX1", "SPACA1", "SPAG11A", "SPAG11B", "SPAM1",  "SPANXE", "SPANXN1", "SPANXN2", "SPANXN4", "SPANXN5", "SPATA16", "SPATA19", "SPINK14", "SPINK7", "SPINK9",  "SPINT3",  "SPINT4",  "SPO11", "SPP2", "SPRR2B",  "SPRR2C",  "SPRR4",  "SPZ1", "SRY",  "SSTR4",   "SSX3", "SSX7", "SSX8", "ST8SIA3", "STATH", "SULT1C3", "SULT2A1", "SULT6B1", "SUN5", "SYCN", "TAAR1",  "TAAR2",  "TAAR5", "TAAR6",   "TAAR8",   "TAAR9",   "TARM1", "TAS2R16", "TAS2R39", "TAS2R41", "TAS2R7", "TBC1D21", "TBC1D28", "TBL1Y", "TBPL2", "TCP10L2", "TCP10",   "TECRL",   "TECTB",  "TEX13A",  "TEX28",   "TFAP2D",  "TGIF2LX", "TGIF2LY", "TGM6", "TINAG", "TMCO5A", "TMEM174", "TMEM202", "TMEM207", "TMEM225", "TMEM30C", "TMEM8C",  "TMEM95",  "TMIGD1", "TMPRSS11BNL", "TMPRSS11B", "TMSB4Y",  "TNP1", "TNP2", "TPD52L3", "TPRX1",  "TREML2P1", "TRIM42",  "TRIM48",  "TRIM49L", "TRIM53", "TRIM60",  "TRIM64",  "TRIM77",  "TRPC7",  "TRYX3",   "TSPY1",   "TSPY2",   "TSPY3", "TSPY4",   "TSSK2",   "TTLL8",   "TTTY10", "TTTY11",  "TTTY12",  "TTTY13",  "TTTY14", "TTTY15",  "TTTY16",  "TTTY17A", "TTTY17B", "TTTY18",  "TTTY19",  "TTTY1B",  "TTTY20", "TTTY21",  "TTTY22",  "TTTY23",  "TTTY2", "TTTY3B",  "TTTY4C",  "TTTY5",   "TTTY6B", "TTTY6",   "TTTY7",   "TTTY8",   "TTTY9B", "TXNDC8",  "UBE2DNL", "UBE2U",   "UBTFL1", "UCMA", "UCN3", "UGT1A3",  "UGT1A4", "UGT1A5", "UGT1A8",  "UNCX", "USP17L6P", "USP17", "USP29", "USP9Y", "UTY", "VCX2", "VCY", "VENTXP1", "VENTXP7", "VN1R2", "VN1R4", "VPREB1", "VSTM2B", "VTRNA1-1",  "VTRNA1-2",  "VTRNA1-3",  "VWC2L",  "WFDC11",  "WFDC8",   "WFDC9",   "XAGE5", "XKRY2", "XKRY", "ZCCHC13", "ZCCHC16", "ZIM3", "ZNF479",  "ZNF645",  "ZNF679", "ZNF705D", "ZNF735",  "ZNRF4", "ZSWIM2"])
	print(X_train.shape)
	print(X_test.shape)
	T_train, T_test = data_train['time'].values, data_test['time'].values
	E_train, E_test = data_train['vital_status'].values, data_test['vital_status'].values

	scaler = StandardScaler()

	# Get scaling parameters with the train sample exclusively, using the Scaler.fit() function
	scaler.fit(X_train)

	# Scale data using Scaler.transform()
	scaled_X_train = pd.DataFrame(scaler.transform(X_train), columns = X_train.columns)
	scaled_X_test = pd.DataFrame(scaler.transform(X_test), columns = X_test.columns)

	# print("\nGENES -----------------")
	# print("scaled_X_train", scaled_X_train)
	# print("\nscaled_X_test", scaled_X_test)

	# print("\nTIMES -----------------")
	# print("T_train", T_train)
	# print("\nT_test", T_test)

	print("\nVITAL STATUS -----------------")
	print("E_train")
	print(E_train)
	unique, counts = np.unique(E_train, return_counts=True)
	print(np.asarray((unique, counts)).T)
	print("\nE_test")
	print(E_test)
	unique, counts = np.unique(E_test, return_counts=True)
	print(np.asarray((unique, counts)).T)

	# TEST NUMBER OF TREES ---------------------------------------------------------------------------------------------------------------------------------------------------------------
	# RSF_test_number_of_trees(number_of_trees, iterations_per_nr_of_trees, max_features, max_depth, min_node_size, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)
	# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# number_of_trees = [1600, 1400, 1200, 1000, 800, 600, 400, 200]
	# RSF_test_number_of_trees(number_of_trees, 2, "sqrt", 953, 10, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)

	# TEST TREE DEPTH --------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# RSF_test_tree_depth(max_depths, iterations_per_nr_of_trees, max_features, number_of_trees, min_node_size, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)
	# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# max_depths = [1905, 1715, 1524, 1334, 1143, 953, 762, 572, 381]
	max_depths = [7618, 6666, 5714, 4762, 3809, 2857, 1905, 953]
	RSF_test_tree_depth(max_depths, 2, "sqrt", 1600, 10, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)
	# RSF_test_tree_depth(max_depths, 5, "sqrt", 1600, 10, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)

	# TEST MIN NODE SIZE -----------------------------------------------------------------------------------------------------------------------------------------------------------------
	# RSF_test_min_node_size(min_node_sizes, iterations_per_nr_of_trees, max_features, number_of_trees, max_depth, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)
	# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# min_node_sizes = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
	# min_node_sizes = [2, 10, 18]
	# RSF_test_min_node_size(min_node_sizes, 10, "sqrt", 1600, 7618, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)

	# TEST BEST MODEL --------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# RSF_test_best_model(min_node_sizes, iterations_per_nr_of_trees, max_features, number_of_trees, max_depth, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)
	# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# min_node_sizes = [10]
	# RSF_test_best_model(min_node_sizes, 30, "sqrt", 1600, 7618, scaled_X_train, scaled_X_test, T_train, T_test, E_train, E_test)
