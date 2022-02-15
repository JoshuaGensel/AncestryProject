import json
import tskit

basic_schema = tskit.MetadataSchema({'codec': 'json'})
complex_schema = tskit.MetadataSchema({
    'codec': 'json',
    'additionalProperties': False,
    'properties': {'accession': {'description': 'ENA accession number',
                                 'type': 'string'},
                   'pcr': {'description': 'Was PCR used on this sample',
                           'name': 'PCR Used',
                           'type': 'boolean'}},
    'required': ['accession', 'pcr'],
    'type': 'object',
})

# Make an example ts with struct metadata
tables = tskit.TableCollection(sequence_length=1)
tables.individuals.metadata_schema = complex_schema
row_id = tables.individuals.add_row(0, metadata={"accession": "Bob1234", "pcr": True})
ts = tables.tree_sequence()

def convert_to_json_metadata(ts):
    # make a new ts with json metadata
    new_tables = ts.dump_tables()
    # iterate through (nearly) all the tables
    for table_name, table in new_tables.name_map.items():
        # provenance table doesn't have metadata
        if table_name not in ["provenances"]:
            # packset_metadata doesn't validate, so dump json in here and switch schema after
            table.packset_metadata([json.dumps(row.metadata).encode() for row in table])
            table.metadata_schema = basic_schema
    # May also need to convert top level metadata?
    return new_tables.tree_sequence()

# test
print(convert_to_json_metadata(ts).individual(0))