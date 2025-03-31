import tiledbsoma    
import cellxgene_census    
from cellxgene_census.experimental.pp import mean_variance    
import pandas as pd  
import argparse
from pathlib import Path
import sys
import json
from datetime import datetime
import pyarrow.compute as pc
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
import tempfile,os

from constants import CL_PINNED_CONFIG_URL, CL_BASIC_OWL_NAME
from cellXgene_backend_utils import get_pinned_ontology_url
import owlready2
from constants import INCLUDED_ASSAYS
from constants import CENSUS_VERSION

import pyarrow as pa

# ontology object
ontology = owlready2.get_ontology(get_pinned_ontology_url(CL_BASIC_OWL_NAME))
ontology.load()
#path=Path(__file__).parent.absolute()
# read in cell_type_ontology.json
json_path='/opt/cell_type_ontology.json'
#json_path = '/workspaces/indapa-CellXGene/cell_type_ontology.json'
with open(json_path) as f:
    cell_type_ontology = json.load(f)

#reverse the cell_type_ontology dictionary
ontology_term_celltype = {v: k for k, v in cell_type_ontology.items()}



"""
Functions to help with processing data from cellXgene and tiledb-soma API
"""

def get_mean_expression_normal(curr_tissue: str, curr_cell_type:str, census_version:str= CENSUS_VERSION) -> pd.DataFrame:
    """
    Calculate the mean expression from raw counts using mean_variance function from cellxgene_census.experimental.pp
   
    """
    with cellxgene_census.open_soma(census_version=census_version) as census:
        human_experiment = census["census_data"]["homo_sapiens"]
        with human_experiment.axis_query(
            measurement_name = "RNA",
            obs_query = tiledbsoma.AxisQuery(                                                                   
               value_filter = f"tissue_general == '{curr_tissue}' and cell_type == '{curr_cell_type}' and is_primary_data == True and disease == 'normal'"), 
        
        ) as query:
            if (query.n_obs > 0 and query.n_vars>0):
                mean_genes = mean_variance(query, axis=0, calculate_mean = True, calculate_variance = False, layer="normalized")
                var_df = query.var().concat().to_pandas()
                var_df_mean = pd.concat([var_df.set_index("soma_joinid"), mean_genes], axis=1)
                var_df_mean["tissue_general"] = curr_tissue
                var_df_mean['cell_type'] = curr_cell_type
            else:
                # return empty dataframe
                var_df_mean =  pd.DataFrame()
        return var_df_mean





def cppt(df: pd.DataFrame, scaling_value:int=1e4) -> pd.DataFrame:
    """
    This function computes the counts per ten thousand (CPTT) of all genes in a given dataframe 
    @param df: dataframe with columns feature_id, feature_name, feature_length, soma_dim_0, soma_data
    @param scaling_value: scaling value for counts per ten thousand (CPTT)

    soma_dim_0 -- corresponding to the cell's soma_joinids
    soma_dim_1 -- corresponding to the gene's soma_joinids
    soma_data -- corresponding to expression value for this gene and cell

    """
    
    # Convert feature_length to base pairs if needed
    df['feature_length'] = df['feature_length'] / 1000

    # Calculate count_per_length using vectorized operations
    df['count_per_length'] = df['soma_data'] / df['feature_length']

    # Calculate the total count_per_length for each soma_dim_0 using groupby and transform
    df['cptt_scaling'] = df.groupby('soma_dim_0')['count_per_length'].transform('sum')

    # Calculate cptt_scaling_x without merging
    df['cptt_scaling_x'] = scaling_value / df['feature_length']

    # Calculate cptt using vectorized operations
    df['cptt'] = df['soma_data'] / df['cptt_scaling'] * df['cptt_scaling_x']

    # Drop intermediate columns if not needed for further processing
    df.drop(['count_per_length', 'cptt_scaling', 'cptt_scaling_x'], axis=1, inplace=True)

    return df


def calc_filtered_average_expression_log_cptt_slack(df:pd.DataFrame, cptt_cutoff:int=3):
    df = df.dropna(subset=['soma_data'])
    df['cptt_1'] = df['cptt'] + 1  # add pseudo count of 1
    df['log_cptt_1'] = np.log(df['cptt_1'])  # calculate log of cptt+1

     # if column cptt is <= cutoff set it to zero
    df['log_cptt_1_masked'] = np.where(df['log_cptt_1'] <= cptt_cutoff, 0, df['log_cptt_1'])

    # get the non-zero cptt_1 values
    non_zero_log_cptt_1 = df.loc[df['cptt_1'] > cptt_cutoff]

    #calculate mean expression of non-zero log_cptt_1 values by feature name
    mean_res = non_zero_log_cptt_1.groupby('feature_name')['log_cptt_1'].mean().reset_index()

    # count the number of unique soma_dim_0 values
    total_cells = df['soma_dim_0'].nunique()

    # count the number of rows with cptt > 0 by feature_id
    cell_expressed_df = (
        df.loc[df['cptt_1'] > cptt_cutoff
        ]
        .groupby('feature_name')['soma_dim_0']
        .count()
        .reset_index(name='cell_expressed_count')
    )

    cell_expressed_df['total_cells'] = total_cells
    cell_expressed_df['percentage_expressed'] = (
        cell_expressed_df['cell_expressed_count'] / total_cells * 100
    )

    mean_exp_pct_exp_res = mean_res.merge(cell_expressed_df, on=['feature_name'])
    return mean_exp_pct_exp_res


def calc_filtered_average_expression_log_cptt_gpt(df:pd.DataFrame, cptt_cutoff:int=3):
    """
    Calculate the mean expression and percentage of cells expressed for each gene in a given dataframe
    """
    df = df.dropna(subset=['soma_data'])
    # if column cptt is <= cutoff set it to zero
    df['cptt'] = np.where(df['cptt'] < cptt_cutoff, 0, df['cptt'])


    df['cptt_1'] = df['cptt'] + 1  # add pseudo count of 1
    df['log_cptt_1'] = np.log(df['cptt_1'])  # calculate log of cptt+1

    # count the number of unique soma_dim_0 values
    total_cells = df['soma_dim_0'].nunique()

    # calculate the mean expression for each gene
    mean_res = df.groupby('feature_name')['log_cptt_1'].mean().reset_index()

    # count the number of rows with cptt > 0 by feature_id
    cell_expressed_df = (
        df.loc[df['cptt'] > 0]
        .groupby('feature_name')['soma_dim_0']
        .count()
        .reset_index(name='cell_expressed_count')
    )

    cell_expressed_df['total_cells'] = total_cells
    cell_expressed_df['percentage_expressed'] = (
        cell_expressed_df['cell_expressed_count'] / total_cells * 100
    )

    mean_exp_pct_exp_res = mean_res.merge(cell_expressed_df, on=['feature_name'])

    return mean_exp_pct_exp_res




def filter_low_expressing_cells(df:pd.DataFrame, threshold:int=500) -> pd.DataFrame:

    """ 
    given a dataframe with columns soma_dim_0, soma_1, and  soma_data, filter out cells with soma_data with <=500 genes expressed
    @param df - pandas dataframe with columns soma_dim_0, soma_1, and  soma_data
    @param threshold - threshold for filtering out cells with soma_data with <=threshold genes expressed

    soma_dim_0 -- corresponding to the cell's soma_joinids
    soma_dim_1 -- corresponding to the gene's soma_joinids
    soma_data -- corresponding to expression value for this gene and cell
    """

    # Group the DataFrame by 'cell_id' and count the number of genes for each cell
    gene_counts = df.groupby('soma_dim_0')['soma_dim_1'].count()
    # Filter the cell IDs that have at least 500 genes
    valid_cell_ids = gene_counts[gene_counts >= threshold].index

    
    output_string=(f"total cells after removing cells with less than {threshold} genes expressed: {len(valid_cell_ids)}\n")
    sys.stderr.write(output_string)
    # Filter the original DataFrame to keep only the rows with cell IDs in valid_cell_ids
    res = df[df['soma_dim_0'].isin(valid_cell_ids)]
    return res

def load_tissue_parquet_file(parquet_path: Path,  tissue_list: list) -> pd.DataFrame:
    """
    This function loads a parquet file into a pandas dataframe and returns a dataframe with cells from the tissue_list
    
    @param path: Path object to parquet directory
    @param tissue_list: list of tissues

    """

    filters = [
            ('organ', 'in', tissue_list)
    ]

    res=pd.read_parquet(parquet_path, filters=filters)
    
    return res

def load_parquet_file(parquet_path: Path, feature_list: list, tissue_list: list, cell_list: list ) -> pd.DataFrame:
    """
    This function loads a parquet file into a pandas dataframe and returns a dataframe with  genes from the feature_list, and cell types from the cell_list and tissues from the tissue_list
   
    """

    filters = [
            ('organ', 'in', tissue_list),
            ('gene_symbol', 'in', feature_list),
            ('cell_type', 'in', cell_list)
    ]

    res=pd.read_parquet(parquet_path, filters=filters)
    return res

def get_census_versions():

    """
    return a json/dictionary with all available census versions
    """
    
    return cellxgene_census.get_census_version_directory

def get_cell_tissue_info(census_version=CENSUS_VERSION) -> pd.DataFrame:

    """
    return dataframe with cell and tissue type information from census
    
    """

    census = cellxgene_census.open_soma(census_version=census_version)

    human = census["census_data"]["homo_sapiens"]

    # Read entire _obs_ into a pandas dataframe.
    obs_df = human.obs.read(column_names=["cell_type","tissue_general"]).concat().to_pandas()
    obs_df.drop_duplicates(inplace=True)
    return obs_df


def get_cell_type_info(census_version:str=CENSUS_VERSION) -> pd.DataFrame:
    """
    return a dataframe with cell type and ontology term information from census
    """
    census = cellxgene_census.open_soma(census_version=census_version)

    human = census["census_data"]["homo_sapiens"]

    # Read entire _obs_ into a pandas dataframe.
    obs_df = human.obs.read(column_names=["cell_type","cell_type_ontology_term_id"]).concat().to_pandas()
    return obs_df
    

def descendants(cell_type:str) -> list:
    """
    @param cell_type string describing the cell type of interest

    This function performs the ontology rollup  using owl2ready API. It is a recusive function to collect all the descendant cell type identifies 


    """
    global ontology
    cell_type_iri = cell_type.replace(":", "_")
    entity = ontology.search_one(iri=f"http://purl.obolibrary.org/obo/{cell_type_iri}")
    descendants = [i.name.replace("_", ":") for i in entity.descendants()] if entity else [cell_type]
    return descendants

def get_ontology_label(term_id:str) -> str:
    """
    @param term_id CL ontology term id

    This function returns the string description for a CL id
    """
    global ontology
    term_id_iri = term_id.replace(":", "_")
    entity = ontology.search_one(iri=f"http://purl.obolibrary.org/obo/{term_id_iri}")
    return entity.label[0] if entity else None

def cell_typeTo_ontology_term(cell_type):
    global cell_type_ontology
    try:
        return cell_type_ontology[cell_type]
    except KeyError:
        sys.stderr.write(f"Cell type {cell_type} not found in cell type ontology\n")
        return None

def ontology_termTo_cell_type(ontology_term:str) -> str:
    global ontology_term_celltype
    try:
        return ontology_term_celltype[ontology_term]
    except KeyError:
        sys.stderr.write(f"Cell type {ontology_term} not found in cell type ontology\n")
        return None
    

def filter_cell_descendants_in_tissue(descendant_list: list, tissue: str) -> pd.DataFrame:

    """
    filter out cell types in descendant_list that are not in census
    """

    curr_path=Path(__file__).parent 
    cell_tissue_file=curr_path / 'cell_tissue_pairs.csv' # file with cell type and tissue information from census

    #open the cell_tissue_file and read each line in a for loop
    #each line is a cell_type and tissue pair

    df=pd.read_csv(cell_tissue_file)

    df_filtered=df.query('tissue_general == @tissue')

    census_cell_types_list=df_filtered['cell_type'].tolist() # list of all cell types in census

    filtered_descendants_list=[cell_type for cell_type in descendant_list if cell_type in census_cell_types_list]

    df=pd.DataFrame(filtered_descendants_list, columns=['descendant_cell_type'])
    df['tissue']=tissue

    #make column indicating if cell type is root or descendant
    # the first cell type in the list is the root cell type
    df['is_root']=False
    df.loc[0,'is_root']=True
    

    return df




def get_cell_ontology_descendants(cell_type: str) -> list:

    """
    return list of cell ontology descendants for a given cell type

    @param cell_type: cell type name (e.g. "T cell")

    """
    

    ontology_term_id = cell_typeTo_ontology_term(cell_type)
    descendants_list=descendants(ontology_term_id)
    cell_types_terms_descendants = [get_ontology_label(ontology_term_id) for ontology_term_id in descendants_list]
    #remove None valuues from cell_types list
    cell_type_descendants = [cell_descendant for cell_descendant in cell_types_terms_descendants if cell_descendant]
    
    #remove the current cell type from the list
    cell_type_descendants=[cell_descendant for cell_descendant in cell_types_terms_descendants if cell_descendant != cell_type]

    #insert cell_type as first element in list
    cell_type_descendants.insert(0, cell_type)

    return cell_type_descendants


def YYYYMMMDD_date()->str:
    """
    This function returns the date in YYYYMMDD format
    """
    current_date = datetime.now()
    formatted_date = current_date.strftime("%Y%m%d")
    return formatted_date

def collect_census_queries_gene(tissue:str, cell_type:str, feature_name:str, census_version:str=CENSUS_VERSION ) -> pd.DataFrame:
    """
    return concatenated dataframe of all cells returned for a given census query for tissue and cell_type and feature_name
    
    The only difference between this function and collect_census_queries 
    is that here we are only considering a single gene instead of all genes
    
    1. get ontology_term_id for cell_type
    2. get all descendants for ontology_term_id
    3. get all cell_types for descendants
    4. query census for cell_types in tissue for feature_name
    
    5. concatenate all cells returned for each cell_type

    """


    ontology_term_id = cell_typeTo_ontology_term(cell_type)
    descendants_list=descendants(ontology_term_id)
    cell_types = [get_ontology_label(ontology_term_id) for ontology_term_id in descendants_list]
    #remove None valuues from cell_types list
    cell_types = [cell_type for cell_type in cell_types if cell_type]

    cfg = {
    "tiledb_config": {
        "soma.init_buffer_bytes": .25 * 1024**3,
        "vfs.s3.no_sign_request": True,
        },
    }
    ctx = tiledbsoma.SOMATileDBContext().replace(**cfg)
    
    sys.stderr.write(f"querying census for {cell_type} in {tissue} for {feature_name}\n")
    #list of dataframes retured from query.X("raw") after filtering low expressing cells
    census_query_results=[]
    
    obs_value_filter_string = f"tissue_general == '{tissue}' and  cell_type in {cell_types} and is_primary_data == True and disease == 'normal' and  assay_ontology_term_id in {INCLUDED_ASSAYS}"
    var_value_filter_string = f"feature_name == '{feature_name}'"
    with cellxgene_census.open_soma(census_version=census_version, context=ctx) as census:
        human_experiment = census["census_data"]["homo_sapiens"]
        with human_experiment.axis_query(
            measurement_name = "RNA",
            obs_query = tiledbsoma.AxisQuery( value_filter=obs_value_filter_string),
                                                           
        
        ) as query:
            var_df = query.var().concat().to_pandas()
            #obs_df = query.obs().concat().to_pandas()
            #n_vars = query.n_vars
            n_obs = query.n_obs
        
            
            if n_obs == 0:
                sys.stderr.write(f"no cells returned for {cell_type} in {tissue}")
                sys.exit(0)

            # query.X() returns an iterator of pyarrow.Table, with X data in COO format (only non-zero values).
            # each arrow_tbl will have three columns:
            #   1. soma_dim_0 -- corresponding to the cell's soma_joinids
            #   2. soma_dim_1 -- corresponding to the gene's soma_joinids
            #   3. soma_data -- corresponding to expression value for this gene and cell
           
            iteration=0
           
            
            for arrow_tbl in query.X("raw").tables():
                sys.stderr.write(f"iterating expression count matrix for {cell_type} in {tissue} for iteration {iteration}\n")      
                # Convert the arrow table to a pandas dataframe, this is the slice of the expression count matrix
               
                arrow_tbl_df=arrow_tbl.to_pandas()
                # filter cells with gene count <=500
                arrow_tbl_filtered_df=filter_low_expressing_cells(arrow_tbl_df, threshold=500)
                #join dataframes with gene info(var_df) with expression count data (arrow_tbl_df)
                joined_df=pd.merge(var_df, arrow_tbl_filtered_df, left_on='soma_joinid', right_on='soma_dim_1', how='inner' ) # join var df to arrow_tbl_df by ids representing the genes
                
                joined_df.reset_index(inplace=True)
                # Append the dataframe to list 
                
                census_query_results.append(joined_df)
                iteration+=1
                
            #concat all the query slices into single df  
            census_query_results_df = pd.concat(census_query_results, ignore_index=True)
            #calculate cppt, which normalizes count to counts per ten thousand
            census_query_results_cppt=cppt(census_query_results_df)
            #return result
            return census_query_results_cppt


def collect_census_total_cells(tissue:str, cell_type:str, census_version:str= CENSUS_VERSION ) -> pd.DataFrame:

    """
    Return the list of unique soma_dim_0 ids for a given cell type and tissue

    """

    ontology_term_id = cell_typeTo_ontology_term(cell_type)
    descendants_list=descendants(ontology_term_id)
    cell_types = [get_ontology_label(ontology_term_id) for ontology_term_id in descendants_list]
    #remove None valuues from cell_types list
    cell_types = [cell_type for cell_type in cell_types if cell_type]
    #remove locstr from cell_types list
    #https://owlready2.readthedocs.io/en/latest/annotations.html?highlight=locstr#language-specific-annotations
    #locstr causes TileDB to throw an error
    cell_types = [i for i in cell_types if not isinstance(i, owlready2.util.locstr)]

    cfg = {
    "tiledb_config": {
        "soma.init_buffer_bytes": .25 * 1024**3,
        "vfs.s3.no_sign_request": True,
        },
    }
    ctx = tiledbsoma.SOMATileDBContext().replace(**cfg)
    
    sys.stderr.write(f"querying census for {cell_type} in {tissue}\n")
    #list of dataframes retured from query.X("raw") after filtering low expressing cells
    #census_query_results=[]
    
    obs_value_filter_string = f"tissue_general == '{tissue}' and  cell_type in {cell_types} and is_primary_data == True and disease == 'normal' and  assay_ontology_term_id in {INCLUDED_ASSAYS}"
   
    with cellxgene_census.open_soma(census_version=census_version, context=ctx) as census:
        human_experiment = census["census_data"]["homo_sapiens"]
        with human_experiment.axis_query(
            measurement_name = "RNA",
            obs_query = tiledbsoma.AxisQuery( value_filter=obs_value_filter_string),
                                                                      
            #value_filter = f"tissue_general == '{tissue}' and  cell_type in {cell_types} and is_primary_data == True and disease == 'normal' and  assay_ontology_term_id in {INCLUDED_ASSAYS}"), 
        
        ) as query:
            var_df = query.var().concat().to_pandas()
            #obs_df = query.obs().concat().to_pandas()
            #n_vars = query.n_vars
            n_obs = query.n_obs
        
            
            if n_obs == 0:
                sys.stderr.write(f"no cells returned for {cell_type}")
                sys.exit(0)

            # query.X() returns an iterator of pyarrow.Table, with X data in COO format (only non-zero values).
            # each arrow_tbl will have three columns:
            #   1. soma_dim_0 -- corresponding to the cell's soma_joinids
            #   2. soma_dim_1 -- corresponding to the gene's soma_joinids
            #   3. soma_data -- corresponding to expression value for this gene and cell
           
            iteration=0
            soma_dim_0_ids_set=()
            
            for arrow_tbl in query.X("raw").tables():
                sys.stderr.write(f"iterating expression count matrix for {cell_type} in {tissue} for iteration {iteration}\n")      
                # Convert the arrow table to a pandas dataframe, this is the slice of the expression count matrix
               
                arrow_tbl_df=arrow_tbl.to_pandas()
                # filter cells with gene count <=500
                arrow_tbl_filtered_df=filter_low_expressing_cells(arrow_tbl_df, threshold=500)
                #get unique number of soma_dim_0 values

                soma_dim_0_ids=arrow_tbl_filtered_df['soma_dim_0'].unique().tolist()
                soma_dim_0_ids_set.update(soma_dim_0_ids)
                iteration+=1

                

            #count unique soma dim0 ids
            unique_count=len(soma_dim_0_ids_set)
            

            #make a dataframe with curr_tissue, curr_cell_type, and unique_count
            df = pd.DataFrame({'tissue': [tissue], 'cell_type': [cell_type], 'unique_count': [unique_count]})
            return df
           


        
        





def collect_census_queries(tissue:str, cell_type:str,  census_version:str= "latest" ) -> pd.DataFrame:
    """
    return concatenated dataframe of all cells returned for a given census query for tissus and cell_type
    

    1. get ontology_term_id for cell_type
    2. get all descendants for ontology_term_id
    3. get all cell_types for descendants
    4. query census for cell_types in tissue
    5. filter out low expressing cells
    6. compute CPPT for each cell
    7. contatenate cppt results
    8. return concatenated cppt results
    
    """

    
    ontology_term_id = cell_typeTo_ontology_term(cell_type)
    descendants_list=descendants(ontology_term_id)
    cell_types = [get_ontology_label(ontology_term_id) for ontology_term_id in descendants_list]
    #remove None valuues from cell_types list
    cell_types = [cell_type for cell_type in cell_types if cell_type]
    #remove locstr from cell_types list
    #https://owlready2.readthedocs.io/en/latest/annotations.html?highlight=locstr#language-specific-annotations
    #locstr causes TileDB to throw an error
    cell_types = [i for i in cell_types if not isinstance(i, owlready2.util.locstr)]

    # convert the locstr to english string
    cell_types_locstr = [i for i in cell_types if  isinstance(i, owlready2.util.locstr)]
    cell_types_locstr = [x.__str__() for x in cell_types_locstr]

    #add the english string to the cell_types list
    cell_types.extend(cell_types_locstr)

    cfg = {
    "tiledb_config": {
        "soma.init_buffer_bytes": .25 * 1024**3,
        "vfs.s3.no_sign_request": True,
        },
    }
    ctx = tiledbsoma.SOMATileDBContext().replace(**cfg)
    
    sys.stderr.write(f"querying census for {cell_type} in {tissue}\n")
    #list of dataframes retured from query.X("raw") after filtering low expressing cells
    #census_query_results=[]
    
    obs_value_filter_string = f"tissue_general == '{tissue}' and  cell_type in {cell_types} and is_primary_data == True and disease == 'normal' and  assay_ontology_term_id in {INCLUDED_ASSAYS}"
   
    with cellxgene_census.open_soma(census_version=census_version, context=ctx) as census:
        human_experiment = census["census_data"]["homo_sapiens"]
        with human_experiment.axis_query(
            measurement_name = "RNA",
            obs_query = tiledbsoma.AxisQuery( value_filter=obs_value_filter_string),
                                                                      
            #value_filter = f"tissue_general == '{tissue}' and  cell_type in {cell_types} and is_primary_data == True and disease == 'normal' and  assay_ontology_term_id in {INCLUDED_ASSAYS}"), 
        
        ) as query:
            var_df = query.var().concat().to_pandas()
            obs_df = query.obs().concat().to_pandas()
            #obs_df.to_csv("/workspaces/indapa-CellXGene/Results/obs_df.csv", index=False)
            #n_vars = query.n_vars
            n_obs = query.n_obs
        
            
            #if n_obs == 0:
            #    sys.stderr.write(f"no cells returned for {cell_type}")
            #    sys.exit(0)

            # query.X() returns an iterator of pyarrow.Table, with X data in COO format (only non-zero values).
            # each arrow_tbl will have three columns:
            #   1. soma_dim_0 -- corresponding to the cell's soma_joinids
            #   2. soma_dim_1 -- corresponding to the gene's soma_joinids
            #   3. soma_data -- corresponding to expression value for this gene and cell
           
            iteration=0
            census_query_results_df=pd.DataFrame()
            
            for arrow_tbl in query.X("raw").tables():
                sys.stderr.write(f"iterating expression count matrix for {cell_type} in {tissue} for iteration {iteration}\n")      
                # Convert the arrow table to a pandas dataframe, this is the slice of the expression count matrix
               
                arrow_tbl_df=arrow_tbl.to_pandas()
                # filter cells with gene count <=500
                arrow_tbl_filtered_df=filter_low_expressing_cells(arrow_tbl_df, threshold=500)
                #join dataframes with gene info(var_df) with expression count data (arrow_tbl_df)
                joined_df=pd.merge(var_df, arrow_tbl_filtered_df, left_on='soma_joinid', right_on='soma_dim_1', how='inner' ) # join var df to arrow_tbl_df by ids representing the genes
                
                # Append the joined DataFrame to the results DataFrame
                census_query_results_df = pd.concat([census_query_results_df, joined_df], ignore_index=True)

                iteration+=1
                
            #calculate cppt, which normalizes count to counts per ten thousand
            census_query_results_cppt=cppt(census_query_results_df)
            #return result
            return census_query_results_cppt


def collect_census_queries_compute_mean_exp(curr_tissue, curr_cell_type):
    """
    wrapper function for collect_census_queries and calc_filtered_average_expression_log_cptt
    """
    cptt_df=collect_census_queries(curr_tissue, curr_cell_type)
    #res=calc_filtered_average_expression_log_cptt_gpt(cptt_df)
    res=calc_filtered_average_expression_log_cptt_slack(cptt_df)
    res['tissue']=curr_tissue
    res['cell_type']=curr_cell_type

    return res









def mean_exp_pct_toParquet(df:pd.DataFrame, output_path:Path,  partition:bool=True) -> None:
    """ 
        @param: df - dataframe of results from cellXgene_census_mean_exp.py
        @param: output_path - Path object to parquet directory where results will be written
        @param: partition - boolean value to determine if parquet should be partitioned by organ and cell_type
        This function writes the dataframe to parquet file. If partition is True, the parquet file is partitioned by organ and cell_type

        This is the function to write the results of the mean expression and percentage of cells expressed to parquet file
        
    
    """
    
    if 'organ' not in df.columns:
        sys.stderr.write("organ not in column names of dataframe. No parquet generated")
        
    elif 'cell_type' not in df.columns:
        sys.stderr.write("cell_type not in column of dataframe. No parquet generated.")
    else:
        date_foldername  = YYYYMMMDD_date()
    
        my_schema = pa.Schema.from_pandas(df)
        #print (my_schema)
        
        res_parquet = output_path / date_foldername
        if partition:
            df.to_parquet(res_parquet, schema=my_schema, partition_cols=['organ', 'cell_type'])
        else:
            
            res_parquet = output_path / date_foldername
            df.to_parquet(res_parquet, schema=my_schema)