#!/usr/bin/env python3
import cellxgene_census
import argparse


"""
return  a Pandas dataframe of with cell and tissue type information from census
"""

CENSUS_VERSION='2023-12-15'

def _get_census_version(version: str):
    census = cellxgene_census.open_soma(census_version=version)    
    human = census["census_data"]["homo_sapiens"]

    # Read entire _obs_ into a pandas dataframe.
    obs_df = human.obs.read(column_names=["cell_type","tissue_general"]).concat().to_pandas()
    obs_df.drop_duplicates(inplace=True)
    return obs_df


def main ():

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--version", help="Census version", default ='2023-12-15')
    
    args = argparser.parse_args()
    census_version = args.version


    res = _get_census_version(census_version)
    
    outfile = ".".join(['cellXgene_census', census_version, 'csv'])
    res.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()