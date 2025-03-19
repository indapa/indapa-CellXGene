INCLUDED_ASSAYS  = ['EFO:0010550',
                          'EFO:0009901',
                          'EFO:0011025',
                          'EFO:0009899',
                          'EFO:0009900',
                          'EFO:0009922',
                          'EFO:0030003',
                          'EFO:0030004',
                          'EFO:0008995',
                          'EFO:0008919',
                          'EFO:0008722',
                          'EFO:0010010']
CENSUS_VERSION='2023-12-15'

CL_PINNED_CONFIG_URL = "https://raw.githubusercontent.com/chanzuckerberg/single-cell-curation/v3.1.3/cellxgene_schema_cli/cellxgene_schema/ontology_files/owl_info.yml"
CL_BASIC_OBO_NAME = "cl-basic.obo"
CL_BASIC_OWL_NAME = "cl-basic.owl"

WMG_PINNED_SCHEMA_VERSION = "3.1.0"

DEPLOYMENT_STAGE_TO_API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
}

TOTAL_STREAMLIT_TISSUES=3
TOTAL_STREAMLIT_CELL_TYPES=193
TOTAL_STREAMLIT_UNIQ_CELLS=4341697