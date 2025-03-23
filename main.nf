#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { get_census_versions; get_mean_expression } from './modules/CellXGene'

workflow {



    get_mean_expression(params.tissue, params.cell_type)
    
}