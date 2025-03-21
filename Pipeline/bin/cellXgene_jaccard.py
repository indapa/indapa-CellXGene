#!/bin/env python
#create main function
import tiledbsoma    
import cellxgene_census    
import pandas as pd  
import argparse
from pathlib import Path
import sys
from datetime import datetime
import numpy as np
from cellXgene_utils import filter_low_expressing_cells, tpm, cppt


def jaccard_index(curr_tissue, curr_cell_type, gene_one, gene_two, census_version='2023-10-18'):
    """
    This computes the co-expression as defined by jaccard indexof two genes in a 
    given tissue and cell type from Cell x Gene Census from normal samples
    
    Jaccard index: https://en.wikipedia.org/wiki/Jaccard_index
    
    @param curr_tissue: tissue name
    @param curr_cell_type: cell type name
    @param census_version: census version
    @param gene_one: gene you want to find co-expression with
    @param gene_two: gene you want to find co-expression with
    
    """
    with cellxgene_census.open_soma(census_version=census_version) as census:
        human_experiment = census["census_data"]["homo_sapiens"]
        with human_experiment.axis_query(
            measurement_name = "RNA",
            
            obs_query = tiledbsoma.AxisQuery(
                value_filter = f"tissue_general == '{curr_tissue}' and cell_type == '{curr_cell_type}' and is_primary_data == True and disease == 'normal'"),
            #var_query = tiledbsoma.AxisQuery(value_filter = f"feature_name in ['{gene_one}', '{gene_two}']"), 
            
        ) as query:
            var_df = query.var().concat().to_pandas()
            obs_df = query.obs().concat().to_pandas().set_index("soma_joinid")
            n_vars = query.n_vars
            n_obs = query.n_obs
            # the total number of genes returned form query is less than 2 or there are zero cells
            #return an empty dataframe 
            if n_vars < 2 or n_obs == 0:
                sys.stderr.write(f"{curr_tissue} and {curr_cell_type} for {gene_one} and {gene_two} doesn't have sufficient data")
                #make an empy dataframe
                return pd.DataFrame({
                    'Gene_one':gene_one,
                    'Gene_two':gene_two,
                    'Cell_type':curr_cell_type,
                    'Tissue':curr_tissue,
                    'Jaccard_index':0,
                    'Jaccard_intersecton':0,
                    'Jaccard_union':0}, index=[0])


            obs_df['jaccard_count']=0

           
            # query.X() returns an iterator of pyarrow.Table, with X data in COO format (only non-zero values).
            # each arrow_tbl will have three columns:
            #   1. soma_dim_1 -- corresponding to the cell's soma_joinids
            #   2. soma_dim_2 -- corresponding to the gene's soma_joinids
            #   3. soma_data -- corresponding to expression value for this gene and cell

            for arrow_tbl in query.X("raw").tables():            
                # pull the slice out of the arrow table that corresponds to the two genes we are interested in
               
                
                arrow_tbl_df=arrow_tbl.to_pandas()
                #filter out cells with soma_data with <=2 genes expressed
                arrow_tbl_df=filter_low_expressing_cells(arrow_tbl_df, threshold=500)
                #count the number of unique soma_dim_0 (cell ids)
                valid_cell_ids, counts = np.unique(arrow_tbl_df["soma_dim_0"], return_counts=True)
                #sys.stderr.write(f"total cells after removing cells with less than 500 genes expressed: {len(valid_cell_ids)}\n")
                joined_df=pd.merge(var_df, arrow_tbl_df, left_on='soma_joinid', right_on='soma_dim_1', how='left' ) # join var df to arrow_tbl_df by ids representing the genes
                joined_df.reset_index(inplace=True)
                #calculate tpm
                res2=tpm(joined_df)
                res2['tpm_1'] = res2['tpm'] + 1
                #take the log of tpm_1
                res2['log_tpm_1'] = np.log(res2['tpm_1'])
                #remove any genes with log_tpm_1 <=3
                res2=res2[res2['log_tpm_1'] > 3]
                # filter out genes that are not gene_one or gene_two
                res2=res2[res2['feature_name'].isin([gene_one, gene_two])]
                if res2.shape[0] == 0:
                    sys.stderr.write(f"{curr_tissue} and {curr_cell_type} for {gene_one} and {gene_two} doesn't have sufficient data")
                    return pd.DataFrame({
                    'Gene_one':gene_one,
                    'Gene_two':gene_two,
                    'Cell_type':curr_cell_type,
                    'Tissue':curr_tissue,
                    'Jaccard_index':0,
                    'Jaccard_intersecton':0,
                    'Jaccard_union':0}, index=[0])
                
                cell_ids, counts = np.unique(res2["soma_dim_0"], return_counts=True)
                
                obs_df.loc[cell_ids, 'jaccard_count'] += counts
            
            # count the occurrence of unique jaccard counts
            value_counts = obs_df['jaccard_count'].value_counts().reset_index()
            value_counts['cell_type']=curr_cell_type
            value_counts['tissue_general']=curr_tissue
            value_counts['gene_1']=gene_one
            value_counts['gene_2']=gene_two
            
            jaccard_dict=value_counts[['jaccard_count','count']].set_index('jaccard_count')['count'].to_dict()
            if 2 not in jaccard_dict.keys():
                jaccard_dict[2]=0
            if 1 not in jaccard_dict.keys():
                jaccard_dict[1]=0
            #jaccard index is the number of cells expressing both genes over the number cells expressing either gene
            jaccard_index=jaccard_dict[2] / jaccard_dict[1]
            # create dataframe gene_one, gene_two, cell type, tissue, jaccard_index
            # return dataframe
            dict={'Gene_one':gene_one, 
                  'Gene_two':gene_two, 
                  'Cell_type':curr_cell_type,
                  'Tissue':curr_tissue,
                  'Jaccard_index':jaccard_index,
                  'Jaccard_intersecton':jaccard_dict[2],
                  'Jaccard_union':jaccard_dict[1]} 
            res=pd.DataFrame(dict, index=[0])
            
            return res




def YYYYMMMDD_date():
    """
    This function returns the date in YYYYMMDD format
    """
    current_date = datetime.now()
    formatted_date = current_date.strftime("%Y%m%d")
    return formatted_date


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--tissue", help="tissue name")
    argparser.add_argument("--cell_type", help="cell type name")
    argparser.add_argument("--gene_one", help="gene one")
    argparser.add_argument("--gene_two", help="gene two")

    
    args = argparser.parse_args()
    curr_tissue = args.tissue
    curr_cell_type = args.cell_type
    gene_one = args.gene_one
    gene_two = args.gene_two

    # get the date
    date = YYYYMMMDD_date()
    # get the jaccard index
    res = jaccard_index(curr_tissue, curr_cell_type, gene_one, gene_two)
    
    #replace spaces in curr_cell_type name
    curr_cell_type = curr_cell_type.replace(" ", "_")
    curr_cell_type = curr_cell_type.replace("/", "_")
    curr_cell_type = curr_cell_type.replace("-", "_")
    curr_cell_type = curr_cell_type.replace("(", "_")
    curr_cell_type = curr_cell_type.replace(")", "_")
    curr_cell_type = curr_cell_type.replace(",", "_")
    curr_cell_type = curr_cell_type.replace("+", "positive")
    
    #replace spaces in curr_tissue name
    curr_tissue = curr_tissue.replace(" ", "_")
    #replace - in curr_tissue name
    curr_tissue = curr_tissue.replace("/", "_")
    curr_tissue = curr_tissue.replace("-", "_")
    curr_tissue = curr_tissue.replace("(", "_")
    curr_tissue = curr_tissue.replace(")", "_")
    curr_tissue = curr_tissue.replace(",", "_")
    curr_tissue = curr_tissue.replace("+", "positive")

    parent_path= Path("/home/aindap/streamlitapps/CellXGeneData/Data/Jaccard")

    output_path = parent_path / date

    if output_path.exists() is False:
        output_path.mkdir(parents=True)
    fname = f'jaccard_index_{curr_tissue}_{curr_cell_type}_{gene_one}_{gene_two}'
    fout = output_path / fname
    # write to file
    res.to_csv(fout, index=False)
    

    
    




if __name__ == '__main__':
    main()