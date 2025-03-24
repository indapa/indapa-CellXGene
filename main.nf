#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { get_census_versions; get_mean_expression } from './modules/CellXGene'

def required_params = [
    'samplesheet', 
    'census_version',
    'output_dir'
]

for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}

def checkSamplesheet(samplesheet_file) {
    if (!file(samplesheet_file).exists()) {
        exit 1, "Samplesheet file not found: ${samplesheet_file}"
    }
    return file(samplesheet_file)
}

ss_status = checkSamplesheet(params.samplesheet)

Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        
        def tissue = file(row.tissue)
        def cell_type = file(row.cell_type)
        
        return tuple(tissue, cell_type)
    }
    .set { input_cellxgene_ch }

workflow {



    get_mean_expression(input_cellxgene_ch)
    //get_census_versions(params.census_version)
    
}