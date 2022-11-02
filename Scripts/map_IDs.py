import pandas
import sys


epigenomics_values_path = sys.argv[1]
mtblmcs_values_path = sys.argv[2]
ids_path = sys.argv[3]


output_dir = sys.argv[4]
output_dir2 = sys.argv[5]
output_dir3 = sys.argv[6]

# Read in dataframes as Pandas dataframes
epigenetics = pandas.read_csv(epigenomics_values_path, index_col=0)
metabolomics = pandas.read_csv(mtblmcs_values_path, index_col=0)
IDs = pandas.read_csv(ids_path, sep=',', index_col=False)


# Write rows and columns vice-versa
epigenetics = epigenetics.transpose()


# Limit IDs df to MethylID and metaboID, drop NA values, remove deplicate rows (samples IDs)
IDs = IDs.loc[:,['XOmicsPhenoID', 'XOmicsMethylID', 'XOmicsmetaboID']]
IDs = IDs.dropna()


# Save duplicate IDs that are discarded on other df
Duplicates_discarded = IDs.duplicated(subset=['XOmicsPhenoID'], keep = 'first')
Duplicates_discarded_df = IDs[Duplicates_discarded]


# Drop duplicates (last) and keep first
IDs = IDs.drop_duplicates(subset=['XOmicsMethylID'])
IDs = IDs.drop_duplicates(subset=['XOmicsmetaboID'])


# Subset IDs that are present in the respective dataframes
IDs = IDs[IDs.XOmicsMethylID.isin(epigenetics.index)]
IDs = IDs[IDs.XOmicsmetaboID.isin(metabolomics.index)]


# Subset samples that are in both dataframes
epigenetics = epigenetics[epigenetics.index.isin(IDs.XOmicsMethylID)]
metabolomics = metabolomics[metabolomics.index.isin(IDs.XOmicsmetaboID)]


# Set PhenoID as rownames 
epigenetics = epigenetics.set_index(IDs.XOmicsPhenoID)
metabolomics = metabolomics.set_index(IDs.XOmicsPhenoID)


# Write dataframes with same IDs to csv files
epigenetics.to_csv(output_dir)
metabolomics.to_csv(output_dir2)

# Write duplicate IDs that are discarded to csv file
Duplicates_discarded_df.to_csv(output_dir3)
