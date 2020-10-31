##############################
## Map data sources
##############################
datasource_files = []
source_files = []
if datasources is not None:
    datasource_files = datasources["data"].tolist()
    source_files = datasources["source"].tolist()
src_map = dict(zip(datasource_files, source_files))

rule manticore_get_external:
    """Get external resources with appropriate method"""
    wildcard_constraints:
        resource_file = "({})".format("|".join(x for x in datasource_files))
    output:
        resource_file = "{resource_file}"
    input: uri = lambda wildcards: manticore_get_external_input(src_map[wildcards.resource_file])
    params:
        uri = lambda wildcards: parse_uri(src_map[wildcards.resource_file]),
        scheme = lambda wildcards: get_uri_scheme(src_map[wildcards.resource_file])
    log: "logs/manticore_get_external/{resource_file}.log"
    wrapper:
        f"{WRAPPER_PREFIX}/utils/manticore_get_external"
