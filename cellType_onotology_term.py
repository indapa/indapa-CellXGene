from cellXgene_utils import get_cell_type_info
import pprint


res= get_cell_type_info()
# convert res to dict with cell type as key and ontology term as value
result_dict = res.set_index('cell_type')['cell_type_ontology_term_id'].to_dict()

# save result_dict as json file
import json
with open('cell_type_ontology.json', 'w') as fp:
    json.dump(result_dict, fp)
