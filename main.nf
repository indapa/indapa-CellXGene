#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { get_census_versions; get_mean_expression; count_unique_cells } from './modules/CellXGene'

assert params.samplesheet : "Parameter 'samplesheet' is required!"
assert params.census_version : "Parameter 'census_version' is required!"
assert params.output_dir : "Parameter 'output_dir' is required!"

def checkSamplesheet(samplesheet_file) {
    assert file(samplesheet_file).exists() : "Samplesheet file not found: ${samplesheet_file}"
    return file(samplesheet_file)
}

def samplesheet = checkSamplesheet(params.samplesheet)

Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        assert row.tissue : "Missing 'tissue' column in samplesheet"
        assert row.cell_type : "Missing 'cell_type' column in samplesheet"
        
        def tissue = row.tissue
        def cell_type = row.cell_type
        
        return tuple(tissue, cell_type)
    }
    .set { input_cellxgene_ch }

workflow {
    get_mean_expression(input_cellxgene_ch)
    count_unique_cells(input_cellxgene_ch)
    
}