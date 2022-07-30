#!/usr/bin/env python

FILE_OBJECTIVE = """Binarize phenotypes based on their field values"""

import argparse
import multiprocessing as mp
import utils as ut


def create_binarized_table(
    exome_df, 
    pheno_info_root, pheno_type, pheno_cat, pheno_id, pheno_ordinal, pheno_storage_root, 
    strategy
    ):
    # get the pheno table path where sample to pheno info is stored
    pheno_table_path = ut.get_pheno_table_filepath(pheno_info_root, pheno_type, pheno_cat, pheno_id)
    pheno_all_df = ut.read_pheno_table(pheno_table_path)
    # get rid of all negative pheno values except for categorical multiples
    pheno_no_negative_vals_df = ut.filter_pheno_table_no_negs(pheno_all_df, pheno_type)
    
    if pheno_type in {"categorical_single", "categorical_multiple"}:
        # merge pheno info values depending on the type of phenotype
        pheno_merged_fields_df = ut.merge_values_categorical(pheno_no_negative_vals_df, pheno_type)
        # binarize categoricals
        pheno_encoding_path = ut.get_pheno_encoding_filepath(pheno_info_root, pheno_type, pheno_cat)
        pheno_encodings = ut.read_pheno_encodings(pheno_encoding_path, pheno_id)
        ohe_encodings,  ordinal_encodings = None, None
        if pheno_ordinal == "B":
            # this type of phenos have both ordinal and ohe type encodings example: field 4537: Work/job satisfaction
            # a separately prepared modified field encodings json file is required for these
            ohe_encoding_path = ut.get_modified_pheno_encoding_filepath(pheno_storage_root, "ohe")
            ordinal_encoding_path = ut.get_modified_pheno_encoding_filepath(pheno_storage_root, "ordinal")
            ohe_encodings = ut.read_pheno_encodings(ohe_encoding_path, pheno_id)
            ordinal_encodings = ut.read_pheno_encodings(ordinal_encoding_path, pheno_id)
        pheno_binarized_df = ut.binarize_categoricals(pheno_merged_fields_df, pheno_type, pheno_encodings, pheno_ordinal, ohe_encodings, ordinal_encodings)
    elif pheno_type in {"integer", "continuous"}:
        # merge pheno info values depending on the type of phenotype
        pheno_merged_fields_df = ut.merge_values_numerical(pheno_no_negative_vals_df)
        # binarize numerical pheno values based on the strategy selected 
        pheno_binarized_df = ut.binarize_numericals(pheno_merged_fields_df, strategy=strategy, quantile_low=0.05, quantile_high=0.95)
    # keep pheno values only for the samples with exome data
    pheno_binarized_df = ut.filter_pheno_with_exomes(pheno_binarized_df, exome_df)
    # save the binarized table of the field in root -> type -> category dir 
    ut.save_pheno_table(pheno_binarized_df, pheno_storage_root, pheno_type, pheno_cat, pheno_id, strategy)
    return

def main(
    phenos_of_interest_file, exome_file, pheno_info_root, pheno_storage_root,
    pheno_type, strategy, threads
    ):
    
    # read phenos of interest df
    phenos_of_interest_df = ut.read_phenos_of_interest_data(phenos_of_interest_file)
    # read exome df
    exome_df = ut.read_exome_data(exome_file)
    # from the phenos of interest, select the ones that fall under an user defined type
    phenos_of_interest_df = phenos_of_interest_df.loc[phenos_of_interest_df.Type.isin([pheno_type])]
    
    pool_iter = [(exome_df, pheno_info_root, t, c, i, o, pheno_storage_root, strategy) for t,c,i,o in zip(
        phenos_of_interest_df.Type, phenos_of_interest_df.Phenotype_group, 
        phenos_of_interest_df.Phenotype_ID, phenos_of_interest_df.not_ordinal)]

    pool = mp.Pool(threads)
    pool.starmap(create_binarized_table, pool_iter)
    pool.close()
    pool.join()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=FILE_OBJECTIVE)
    parser.add_argument(
        "phenos_of_interest_file", 
        type=str, 
        help="""The file path of the manually prepared excel file that contains information about 
        the phenotypes' type (column name: Type), category (column name: Phenotype_group) and 
        field id: (column name: Phenotype_ID), field ordinality (column name: not_ordinal)"""
        )
    parser.add_argument(
        "id2exome_file", 
        type=str, 
        help="""The file path of the csv file downloaded from UKB that maps sample ids to their 
        exome vcf values"""
        )
    parser.add_argument("pheno_info_root", type=str, help="The folder where previously downloaded fields and their encodings are stored")
    parser.add_argument("pheno_storage_root", type=str, help="The folder where binarized phenotype tables will be stored")
    parser.add_argument("pheno_type", type=str, help="The phenotype type which will be binarized eg: categorical_single/integer/continuous/categorical_multiple")
    parser.add_argument("-s", "--strategy", type=str, help="The binarizing strategy for integer and continuous type; can be either quantile or median", default="")
    parser.add_argument("-n", "--n_threads", type=int, help="number of cores to use", default=64)

    args = parser.parse_args()

    main(
        args.phenos_of_interest_file,
        args.id2exome_file, 
        args.pheno_info_root, 
        args.pheno_storage_root, 
        args.pheno_type,
        args.strategy,
        args.n_threads
        )
