
process get_mean_expression {
   
    publishDir "${params.output_dir}", mode: 'copy'
    
    container 'indapa/indapa-cellxgene:latest'
    
    input:
    val(tissue)
    val(cell_type)

    output:
    path('*.csv'), emit: tissue_cell_type_expression

    script:
    """
    python ${workflow.projectDir}/bin/cellXgene_census_mean_exp.py --tissue "${tissue}" --cell_type "${cell_type}"
    """
}

process get_census_versions {
   
    publishDir "${params.output_dir}", mode: 'copy'
    
    //container 'community.wave.seqera.io/library/cellxgene-census_pip_tiledbsoma:602e534d4ed2a75c'
    container 'indapa/indapa-cellxgene:latest'

    
    input:
    val(census_version)

    output:
    path('*.csv'), emit: census_version_info

    script:
    """
    python ${workflow.projectDir}/bin/get_census_versions.py --version ${census_version} 
    """
}