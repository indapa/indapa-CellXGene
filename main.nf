#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { get_census_versions; get_mean_expression } from './modules/CellXGene'

workflow {

    tissue = "skin of body"
    cell_type = "macrophage"

    get_mean_expression(tissue, cell_type)
    //get_census_versions(params.census_version)
  
}