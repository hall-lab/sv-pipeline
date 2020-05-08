#!/usr/bin/env python

import sys, textwrap, os

import click
import jinja2

def input_json_template():
    text = '''\
    {
      "##_COMMENT1": "DATA INPUTS",
      "Pre_Merge_SV.aligned_crams": [
          {% for cram in cramList -%}
          "{{ cram }}" {{ "," if not loop.last }}
          {% endfor -%}
      ],
      "Pre_Merge_SV.aligned_cram_suffix": ".cram",
      "##_COMMENT2": "REFERENCE FILES",
      "Pre_Merge_SV.ref_fasta": "gs://human-b38/GRCh38DH/all_sequences.fa",
      "Pre_Merge_SV.ref_fasta_index": "gs://human-b38/GRCh38DH/all_sequences.fa.fai",
      "Pre_Merge_SV.call_regions_bed": "gs://human-b38/GRCh38DH/annotations/canonical_chromosome.bed.gz",
      "Pre_Merge_SV.call_regions_bed_index": "gs://human-b38/GRCh38DH/annotations/canonical_chromosome.bed.gz.tbi",
      "Pre_Merge_SV.ref_cache": "gs://human-b38/cache.tar.gz",
      "Pre_Merge_SV.exclude_regions": "gs://human-b38/GRCh38DH/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed",
      "Pre_Merge_SV.center":  "washu",
      "Pre_Merge_SV.cohort":  "{{ cohort }}",
      "##_COMMENT3": "SYSTEM PARAMETERS",
      "Pre_Merge_SV.preemptible_tries": 5
    }
    '''
    return textwrap.dedent(text)

def option_json_template():
    text = '''\
    {
      "workflow_failure_mode": "ContinueWhilePossible",
      "enable_ssh_access": true,
      "jes_gcs_root": "gs://{{ bucket }}/{{ cohort }}/pre-merge/cromwell/cromwell-executions",
      "read_from_cache": true,
      "monitoring_script": "gs://wustl-ccdg-cromwell/mem_monitor.sh",
      "final_workflow_outputs_dir": "gs://{{ bucket }}/{{ cohort }}/pre-merge/cromwell/final-outputs",
      "default_runtime_attributes": {
        "zones": "us-central1-a us-central1-b us-central1-c us-central1-f"
      },
      "google_labels":{
        "user": "idas",
        "project": "sv-pre-merge-2020-05-08",
        "cohort": "{{ cohort }}",
        "pipeline": "sv-pipeline"
      }
    }
    '''
    return textwrap.dedent(text)

def derive_cram_list(sample_map, cohort):
    cram_list = []
    columns = [0, 2]
    with open(sample_map, 'r') as f:
        for line in f:
            items = line.rstrip().split("\t")
            (cohort_label, cram_path) = [ items[c] for c in columns ]
            if cohort_label == cohort:
                cram_list.append(cram_path)
    return cram_list

def derive_cohort_input_path(cohort):
    # get the directory path of currently running script
    script_path = os.path.dirname(os.path.realpath(__file__))
    input_path = os.path.join(script_path, '..', 'inputs')
    cohort_path = os.path.normpath(os.path.join(input_path, cohort))
    return cohort_path

def create_input_json_file_path(cohort):
    cohort_path = derive_cohort_input_path(cohort)
    if not os.path.isdir(cohort_path):
        os.makedirs(cohort_path)
    input_json = os.path.join(cohort_path, "{}.input.json".format(cohort))
    return input_json

def create_option_json_file_path(cohort):
    cohort_path = derive_cohort_input_path(cohort)
    if not os.path.isdir(cohort_path):
        os.makedirs(cohort_path)
    option_json = os.path.join(cohort_path, "{}.options.json".format(cohort))
    return option_json

def create_input_json(cohort, sample_map):
    cram_list = derive_cram_list(sample_map, cohort)
    if not cram_list:
        sys.exit("[err] Found no crams for cohort: {}".format(cohort))

    json_template = input_json_template()
    template = jinja2.Template(json_template)
    rendered_input_json = template.render(cramList=cram_list, cohort=cohort)
    input_json_file_path = create_input_json_file_path(cohort)
    with open(input_json_file_path, 'w') as f:
        print(rendered_input_json, file=f)

def create_option_json(cohort, bucket):
    json_template = option_json_template()
    template = jinja2.Template(json_template)
    rendered_option_json = template.render(bucket=bucket, cohort=cohort)
    option_json_file_path = create_option_json_file_path(cohort)
    with open(option_json_file_path, 'w') as f:
        print(rendered_option_json, file=f)

@click.command()
@click.option('--json-type', required=True, type=click.Choice(['input', 'option']),
              help="The type of json file to create")
@click.option('--sample-map', default=None, type=click.Path(exists=True),
              help="The input sample map containing samples/crams to process")
@click.option('--bucket', default="eocad-sv-052020", type=click.STRING,
              help="The bucket name to place cromwell inputs/outputs to")
@click.option('--cohort', required=True, type=click.STRING,
              help="The cohort to create the input/option json file for")
def main(json_type, sample_map, bucket, cohort):
    if json_type == 'input' and sample_map == None:
        msg = ("[err] To create an input json you need to "
               "supply a sample-map list via --sample-map")
        sys.exit(msg)

    if json_type == 'input':
        create_input_json(cohort, sample_map)
    else:
        create_option_json(cohort, bucket)

if __name__ == "__main__":
    main()
