# Prepare one hot encoded phenotypes for all samples with exome data
UKBiobank denotes a field as any specific information about an individual, such as *cooked vegetable intake*. A field is divided into a category, for example cooked vegetable intake falls under the category of *diet*. Each field also has a type, for instance the cooked vegetable intake field values are of type *integer*. Some of the field values also have special data encodings.

These fields along with their encodings and annotations provided by UKBiobank were previously downloaded and described here: https://github.com/deeprob/ukb_field_download.

The current repository aims to "one-hot-encode" these fields based on their field values. The final result will be a table that maps all UKB samples with exomes to their one hot encoded fields. The format of the table will be:

```
eid, {field-id-1}_{one-hot-encoded-field-val-1}, ..., {field-id-x}_{one-hot-encoded-field-val-x}, ...
```

# Download preselected fields of interest file
The fields of interest were manually chosen. The file with marked fields of interest is present in:

Dropbox -> UK_BioBank -> BMI_project -> lifestyle.xlsx

**Note**: This excel file must contain the following 5 columns:

1. *Type*: Denotes a UKB field's type, for e.g. categorical_single
2. *Phenotype_group*: Denotes a UKB field's category, for e.g. Mental Health
3. *Phenotype_ID*: Denotes a UKB field's unique id, for e.g. 6156
4. *shortlist*: Denotes the shortlisted fields. Should be marked as "X"
5. *Num_exome_samples_with_phenotype*: Denotes the number of exome samples present for the fields.

# Download the bulk file that contains sample id to exome vcf mappings
This file was downloaded as part of the download exomes from UKB pipeline as described here: https://github.com/deeprob/ukb_bulk_download.

Location -> /data5/deepro/ukbiobank/download/download_bulk/data/ukb48799.csv

# Steps to one hot encode UKB fields
## Filtering negative field values
At first we filtered out all the negative field values of field types a) *integer*, b) *continuous* and c) *categorical single* for each individual since negative field values for these types denote irrelevant information such as "Do not know" or "Prefer not to answer". 

For *categorical multiple* field types, along with the non-negative field values, we also included "None of the above" field value encoded as -7.

## Merging multiple field values
Each individual may have multiple field values for the same field. 

For field types a) *integer* and b) *continuous*, we took the mean of these multiple values as the final merged value. 

For field type, *continuous single*, we took the consensus of these multiple values as the final merged value if there is a consensus, else we ignored that particular individual. 

For field type, *continuous multiple*, since an individual might have multiple field values all of which are non-redundant information about that individual, for example for field id 20544 which denotes "Mental Health problem diagnosed by a professional", an individual can have "Panic attacks", "Depression" and "Anxiety, nerves or generalized anxiety disorder", we did not merge these field values.

## One hot encoding field values
### Integer and Continuous field types
For the numerical field types a) *integer* and b) *continuous*, we first generated the 25th and 75th quantiles from the field values using equal frequency binning. Then, we marked the 25th quantile as low and 75th quantile as high thresholds. Finally, for each of these fields, we created two one-hot encoded columns, 1) for the low threshold category where all individuals with field value lower than the 25th quantile were marked as 1 and higher were marked as 0 and 2) for the high threshold category where all individuals with field value higher than the 75th quantile were marked as 1 and lower were marked as 0. Thus, in the final meta table, each field of type *integer* or *continuous* will be represented by two columns, low and high.

### Categorical single field types
For the ordinal field type, *categorical single*, individuals were characterized into two levels, low and high. The low level was defined as the field value(s) with the lowest rank(s) assigned to it by UKB such that this level has at least 10% of total number of samples in that particular field and excluding the field value with the highest rank. For example field 20407 has the following ranked field values:

*Never, Less-than-monthly, Monthly, Weekly, Daily-or-almost-daily*

Here, the lowest ranked field value is *Never*. We first define the low level as Never and check if that field value has greater than 10% of the total samples present in this field. If it does, then we stop and keep the low level for this field as *Never*. If it does not, then we include the next ranked field value *Less-than-monthly*, check again for number of samples, and we keep on doing it until we reach the field value ranked just lower than the highest ranked field value, *Weekly* thus excluding the inclusion of the highest ranked field value *Daily-or-almost-daily* in the low level. Once the field values for the low level are selected, we one hot encode each individual on the basis of whether their field value is present among the low level ranked field values.

For the high level, we start with the highest ranked field value, check for the number of samples and only include those field values which were not previously included in the low level. Once the high level field values are defined, we again one hot encode each individual on the basis of whether their field value is present among the high level ranked field values. Thus, in the final meta table, each field of type *categorical single* will be represented by two columns, low and high.

**Note**: There might be cases where the lowest (highest) ranked field value has greater than 90% samples. In that case the high (low) level will contain less than 10% of the total samples irrespective of the number of ranked field values included in that level. These cases can be filtered later based on end user requirements.

### Categorical multiple field types
The *categorical multiple* fields were encoded in the true one-hot-encoding sense where each field value were separately converted to zeros or ones based on their absence or presence in an individual. For example field 6155, Vitamin and mineral supplements has field values:

*Vitamin A, Vitamin B, Vitamin C, Vitamin D, Vitamin E, Folic acid or Folate (Vit B9), Multivitamins +/- minerals, None of the above*

We started by one-hot-encoding *Vitamin A* at first by checking if an individual has reported taking *Vitamin A*, then we encoded *Vitamin B* and so on. In the final meta table, each field of type *categorical multiple* will be represented by N columns, N denoting the number of unique field values for that field.

# Script descriptions
The script *0_binarize_phenos.py* has 5 required arguments

1. *phenos_of_interest_file*: preselected field of interest file described in section 2
2. *id2exome_file*: sample id to exome vcf file described in section 3
3. *pheno_info_root*: The path to the root folder where previously downloaded fields and their encodings are stored as mentioned in section 1
4. *pheno_storage_root*: The path to the root folder where one-hot-encoded field tables created by this script will be stored
5. *pheno_type*: The type of the field to be one-hot-encoded. Can be one of integer/continuous/categorical_single/categorical_multiple

This script will convert all fields present in the phenos of interest file that are of the specific user mentioned type into their one hot encoded values and stored them under the field type and category. The field type -> field category -> field id hierarchy is represented as a directory structure within the root dir, *pheno_storage_root*.
