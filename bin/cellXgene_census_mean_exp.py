#!/bin/env python

# calculates the mean expression and percentage of expressed cells for a given cell type and tissue
  

import argparse
import pandas as pd
from pathlib import Path
from cellXgene_utils import collect_census_queries_compute_mean_exp
from cellXgene_utils import YYYYMMMDD_date


def _rename_reformat_res_columns(res:pd.DataFrame) -> pd.DataFrame:
    """
    Rename and reformat columns of res dataframe to match the desired output format
    """
    new_res = res.rename(columns={"feature_name": "gene_symbol",
                          'cell_type': 'cell_type',
                          'tissue':'organ',
                          'total_cells':'total_num_of_cells',
                          'percentage_expressed':'pct_of_expressed',
                          'total_cells':'total_num_of_cells',
                          'log_cptt_1': 'average_expression_of_expressed'})
    new_res['disease']='normal'
    new_res=new_res[['gene_symbol', 'cell_type', 'organ', 'disease','cell_expressed_count', 'total_num_of_cells', 'pct_of_expressed', 'average_expression_of_expressed']]
    return new_res

def main():
    argparser = argparse.ArgumentParser(description="Calculates the mean expression and percentage of expressed cells for a given cell type and tissue from CellXGene Census")
    argparser.add_argument("--tissue", help="tissue name")
    argparser.add_argument("--cell_type", help="cell type name")
    
   
    
    
    args = argparser.parse_args()
    curr_tissue = args.tissue
    curr_cell_type = args.cell_type

   

    

    res=collect_census_queries_compute_mean_exp(curr_tissue, curr_cell_type)
    res=_rename_reformat_res_columns(res)

    
    
    #replace any spaces in curr_cell_type
    curr_cell_type = curr_cell_type.replace(" ", "_")
    
    #replace any spaces in curr_tissue
    curr_tissue = curr_tissue.replace(" ", "_")
    
    fname=curr_tissue + "_" + curr_cell_type + ".csv"
    
    res.to_csv(fname, index=False)
   




if __name__ == '__main__':
    main()