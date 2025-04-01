This repo contains code query data from [CZ CellXGene Discover](https://cellxgene.cziscience.com/) by using [CZ CellXGene Census API](https://chanzuckerberg.github.io/cellxgene-census//) and [TileDB-SOMA API](https://github.com/single-cell-data/TileDB-SOMA). A Nextflow workflow is used to generate results from the CZ CellXGene Discover API and TileDB-SOMA API. 

## Running the pipeline

The easiest way to run this pipeline is on Sequera Cloud. I personally recommend setting up [Batch Forge](https://docs.seqera.io/platform/25.1/compute-envs/aws-batch#tower-forge). Note, you need to have your own AWS account to set this up and you will be charged for the resources you use.

Another option is is to run the pipieline locally. You can create a Github Codespace and once it's created you can run the pipeline using the following command:

```
nextflow run indapa-CellXGene/main.nf -c nextflow.config --output_dir <output_dir> --samplesheet <samplesheet>
``` 

## Parameters required for pipeline

Required parameters are output_dir and samplesheet. An example [samplesheet](https://github.com/indapa/indapa-CellXGene/blob/master/Samplesheets/samplesheet-test.csv) is provided in the repo. 

See the [schema](https://github.com/indapa/indapa-CellXGene/blob/master/nextflow_schema.json) for the required parameters.







