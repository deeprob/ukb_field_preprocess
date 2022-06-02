#!/usr/bin/env python

FILE_OBJECTIVE = """Prepare meta table of binarized phenotypes"""

import argparse
import os
import pandas as pd
import utils as ut
import multiprocessing as mp


def format_pheno_table(
	exome_index, 
	pheno_info_root, pheno_storage_root, pheno_type, pheno_cat, pheno_id, strategy):

	if pheno_type in {"categorical_single", "categorical_multiple"}:
		strategy = ""

	pheno_df_path = ut.get_binarized_table_path(pheno_storage_root, pheno_type, pheno_cat, pheno_id, strategy)
	pheno_df = ut.read_binarized_table(pheno_df_path)

	if pheno_type in {"integer", "continuous", "categorical_single"}:
		pheno_df =  ut.reindex_binarized_table1(pheno_df, pheno_id, exome_index)

	elif pheno_type in {"categorical_multiple"}:
		pheno_encodings = ut.get_field_encodings(pheno_info_root, pheno_type, pheno_cat, pheno_id)
		if type(pheno_encodings) == dict:
			pheno_df =  ut.reindex_binarized_table2(pheno_df, pheno_id, pheno_encodings, exome_index)
	return pheno_df


def main(phenos_of_interest_file, exome_file, pheno_info_root, pheno_storage_root, strategy):
	
	# read phenos of interest df
	phenos_of_interest_df = ut.read_phenos_of_interest_data(phenos_of_interest_file)
	# get exome index
	exome_index = ut.get_exome_index(exome_file)

	pool_iter = [(exome_index, pheno_info_root, pheno_storage_root, t, c, i, strategy) for t,c,i in zip(
 		phenos_of_interest_df.Type, phenos_of_interest_df.Phenotype_group, 
        phenos_of_interest_df.Phenotype_ID)]

	pool = mp.Pool(64)
	pheno_dfs = [pdf for pdf in list(pool.starmap(format_pheno_table, pool_iter)) if not pdf.empty]
	pool.close()
	pool.join()

	meta_df = pd.concat(pheno_dfs, axis=1, join="inner")
	meta_df_path = os.path.join(pheno_storage_root, "meta_pheno_table.csv")
	meta_df.to_csv(meta_df_path)

	return


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description=FILE_OBJECTIVE)
	parser.add_argument(
		"phenos_of_interest_file", 
		type=str, 
        help="""The file path of the manually prepared excel file that contains information about 
        the phenotypes' type (column name: Type), category (column name: Phenotype_group) and 
        field id: (column name: Phenotype_ID)"""
        )
	parser.add_argument(
        "id2exome_file", 
        type=str, 
        help="""The file path of the csv file downloaded from UKB that maps sample ids to their 
        exome vcf values"""
        )
	parser.add_argument("pheno_info_root", type=str, help="The folder where previously downloaded fields and their encodings are stored")
	parser.add_argument("pheno_storage_root", type=str, help="The folder where binarized phenotype tables are stored and the meta table will be stored")
	parser.add_argument("-s", "--strategy", type=str, help="The binarizing strategy for integer and continuous type; can be either quantile or median", default="quantile")
	args = parser.parse_args()

	main(
		args.phenos_of_interest_file,
		args.id2exome_file, 
		args.pheno_info_root,
		args.pheno_storage_root,
		args.strategy
		)
