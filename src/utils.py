import os
import json
import pandas as pd
import numpy as np


#######################
# reading input files #
#######################

def read_phenos_of_interest_data(file, min_samples=2000):
    """
    This function reads interesting phenotype file of type xlsx prepared manually, 
    selects the shortlisted phenotypes with at least 2000 exome samples and
    returns a filtered dataframe 
    """
    df = pd.read_excel(file)
    df = df.loc[df.shortlist=="X"]
    df = df.loc[df.Num_exome_samples_with_phenotype>=min_samples]
    return df

def read_exome_data(sample_to_exome_file):
    """
    This function reads the sample mapped to their exome vcf files,
    selects only those sample ids, which have both vcf and vcf index,
    return a filtered dataframe
    """
    df = pd.read_csv(sample_to_exome_file, index_col=0)
    df = df.dropna(how="all").reset_index()
    return df

def get_pheno_table_filepath(root_dir, pheno_type, pheno_cat, pheno_id):
    """
    This function accepts
    1) root_dir: where the sample to phenotype value table is stored under type -> category -> id hierarchy
    2) pheno_type: the type which the field belongs to
    3) pheno_cat: the category which the field belongs to 
    4) pheno_id: the UKBiobank id of the field
    It returns
    the concatenated filepath of the sample to field value table
    """
    pheno_table_path = os.path.join(
        root_dir, pheno_type, pheno_cat, "tables", f"{pheno_id}.csv"
        )
    assert os.path.exists(pheno_table_path)
    return pheno_table_path


def read_pheno_table(pheno_table_path):
    df =  pd.read_csv(pheno_table_path, index_col=0)
    return df


##############################
# filtering phenotype values #
##############################

def filter_pheno_table_no_negs(pheno_df_all, pheno_type):
    """
    Filter all samples which have negative field values
    Note: negative field values are associated with irrelavant information
    except in case of categorical multiple type of fields where -7 
    which denotes None of the above has a meaning 
    """
    if pheno_type in {"categorical_multiple"}:
        df = pheno_df_all.loc[(pheno_df_all==-7).any(axis=1)|(pheno_df_all>=0).any(axis=1)]
    else:
        df = pheno_df_all.loc[(pheno_df_all>=0).any(axis=1)]
    return df


########################
# merging field values #
########################

def merge_values_categorical(pheno_df_no_negative_vals, pheno_type):
    """
    For categorical_single:
    Merges all the field values for UKBiobank into a single field value 
    by taking the consensus of all field values. If there is not consensus, 
    it is marked as NaN.

    For categorical_multiple:
    Since these types of fields can take multiple values, there is not need to merge
    into a consensus value. Therefore it is returned as is for downstream processing
    """
    if pheno_type == "categorical_single":
        pheno_df_no_negative_vals_nona = pheno_df_no_negative_vals.fillna(method='bfill', axis=1).iloc[:,0]
        pheno_df_no_negative_vals_consistent = pheno_df_no_negative_vals.nunique(axis=1)==1
        pheno_df_no_negative_vals_consensus = pd.Series([i if j==True else np.nan for i,j in zip(pheno_df_no_negative_vals_nona, pheno_df_no_negative_vals_consistent)])
        assert len(pheno_df_no_negative_vals) ==  len(pheno_df_no_negative_vals_consensus)
        pheno_df_no_negative_vals["merged"] = pheno_df_no_negative_vals_consensus.values
        # drop rows that do not reach a consensus
        pheno_df_no_negative_vals = pheno_df_no_negative_vals.dropna(subset=["merged"])
    return pheno_df_no_negative_vals

def merge_values_numerical(pheno_df_no_negative_vals):
    """
    Merges all the field values for UKBiobank into a single field value 
    by taking the mean of all field values.
    """
    pheno_df_no_negative_vals["merged"] = pheno_df_no_negative_vals.mean(axis=1, numeric_only=True)
    return pheno_df_no_negative_vals


#####################
# binarizing values #
#####################

def binarize_numericals(df, strategy="median", quantile_low=0.25, quantile_high=0.75):

    ser = df["merged"]

    if strategy == "median":
        thresh = ser.median()
        binarized_low = (ser<=thresh).astype(int)
        binarized_high = (ser>thresh).astype(int)
        df[f"binarized_{thresh}_low"] = binarized_low
        df[f"binarized_{thresh}_high"] = binarized_high

    elif strategy == "quantile":
        qlow = ser.quantile(quantile_low)
        qhigh = ser.quantile(quantile_high)

        binarized_low = (ser <= qlow).astype(int)
        binarized_high = (ser >= qhigh).astype(int)
        df[f"binarized_{qlow}_low"] = binarized_low
        df[f"binarized_{qhigh}_high"] = binarized_high

    return df


def decide_categorical_bins(ser, field_encodings):

    field_encodings_relevant = sorted([int(fe) for fe in field_encodings.keys() if int(fe)>=0])
    # case 1, if there are two types of relevant categories for this field,
    # just return low and high
    if len(field_encodings_relevant) == 2:
        # low category
        fe_low = field_encodings_relevant[0]
        fe_low_val = field_encodings[str(fe_low)]
        fe_low_val = "-".join(fe_low_val.split())
        # high category
        fe_high = field_encodings_relevant[-1]
        fe_high_val = field_encodings[str(fe_high)]
        fe_high_val = "-".join(fe_high_val.split())
        # binarize
        binarized_low = (ser == fe_low).astype(int)
        binarized_high = (ser == fe_high).astype(int)

    # case 2, if there are more than two types of relevant categories for this field
    if len(field_encodings_relevant) > 2:
        # check how many low bins will give us at least 10% samples 
        # excluding the highest bin
        for lower_bin_range in range(1, len(field_encodings_relevant)):
            fe_lows = field_encodings_relevant[:lower_bin_range]
            fe_low_val = "|".join(["-".join(field_encodings[str(fe_low)].replace(",", "").split()) for fe_low in fe_lows])
            binarized_low = ser.isin(fe_lows).astype(int)
            if (sum(binarized_low)/len(binarized_low)) > 0.1:
                break
        
        # check how many high bins will give us at least 10% samples 
        # excluding all the bins selected by the lower bin range
        for higher_bin_range in range(len(field_encodings_relevant) - 1, lower_bin_range - 1, -1):
            fe_highs = field_encodings_relevant[higher_bin_range:]
            fe_high_val = "|".join(["-".join(field_encodings[str(fe_high)].replace(",", "").split()) for fe_high in fe_highs])            
            binarized_high = ser.isin(fe_highs).astype(int)
            if (sum(binarized_high)/len(binarized_high)) > 0.1:
                break
 
    return binarized_low, binarized_high, fe_low_val, fe_high_val


def binarize_categoricals(df, field_type, field_encodings):


    if field_type == "categorical_single":
        
        ser = df["merged"]
        binarized_low, binarized_high, fe_low_val, fe_high_val = decide_categorical_bins(ser, field_encodings)

        # check if the number of samples in the lowest and highest bins 
        # are less than 10% of the original number of samples
        if (sum(binarized_low)/len(binarized_low)) < 0.1:
            print(f"Warning:: lowest category has less than 10% samples for field id {df.columns[0]}")
            print(f"Warning:: lowest category has {sum(binarized_low)} samples")

        if (sum(binarized_high)/len(binarized_high)) < 0.1:
            print(f"Warning:: highest category has less than 10% samples for field id {df.columns[0]}")
            print(f"Warning:: highest category has {sum(binarized_high)} samples w/ or w/o exomes")

        df[f"binarized_{fe_low_val}_low"] = binarized_low
        df[f"binarized_{fe_high_val}_high"] = binarized_high

    if field_type == "categorical_multiple":
        if type(field_encodings) == dict:
            # only get the relevant field encoding values
            field_encodings_relevant = sorted([int(fe) for fe in field_encodings.keys() if int(fe)>=0])
            
            # add "None of the above", encoding value "-7", here if it is present
            if "-7" in field_encodings.keys():
                field_encodings_relevant.append(-7)

            for fe in field_encodings_relevant:
                fe_ser = df.isin([fe]).any(axis=1).astype(int)
                fe_val = field_encodings[str(fe)]
                fe_val = "-".join(fe_val.replace(",", "").split())
                df[f"binarized_{fe_val}"] = fe_ser
        else:
            print(f"Warning :: Field id {df.columns[0]} has incorrect field encoding type: {field_encodings}.")
            print(f"Warning :: It will not be binarized.")    

    return df


########################
# filtering for exomes #
########################

def filter_pheno_with_exomes(binarized_df, exome_df):
    binarized_df = binarized_df.loc[binarized_df.index.isin(exome_df.eid)]
    return binarized_df


##########################
# saving binarized table #
##########################

def save_pheno_table(
    binarized_pheno_df, storage_root_dir, 
    pheno_type, pheno_cat, pheno_id, 
    method):

    table_basename = f"{pheno_id}_{method}" if method else str(pheno_id)
    binarized_pheno_path = os.path.join(
        storage_root_dir, pheno_type, pheno_cat, "tables", f"{table_basename}.csv"
        )
    os.makedirs(os.path.dirname(binarized_pheno_path), exist_ok=True)
    binarized_pheno_df.to_csv(binarized_pheno_path)
    return


##########################
# field encodings parser #
##########################

def get_pheno_encoding_filepath(root_dir, pheno_type, pheno_cat):
    """
    This function accepts
    1) root_dir: where the sample to phenotype value table is stored under type -> category -> id hierarchy
    2) pheno_type: the type which the field belongs to
    3) pheno_cat: the category which the field belongs to 
    It returns
    the filepath of the json file that contains fields encodings of all fields 
    that fall under the type and category specified
    """
    pheno_json_path = os.path.join(
        root_dir, pheno_type, pheno_cat, f"fields_data_coding.json"
        )
    assert os.path.exists(pheno_json_path)
    return pheno_json_path


def read_pheno_encodings(pheno_json_path, pheno_field_id):
    pheno_field_id = str(pheno_field_id)
    with open(pheno_json_path, "r") as f:
        field_encoding_dict = json.load(f)
    return field_encoding_dict[pheno_field_id]
