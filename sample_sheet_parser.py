import re
import sys

def sampleSheetParser (handle):
    """
    Parse the contents of SampleSheet.csv, convert contents into a
    Python dictionary object.
    Samples are tracked by sample name and sample number (e.g., S9).
    This is to distinguish replicates of the same sample.
    """

    # FIXME: Conan is going to start annotating samples as either amplicon or nextera
    # FIXME: We may need to change this code to handle this

    tag = None
    get_header = False
    header = [] # store Data block column labels
    sample_number = 1 # 1-indexing
    run_info = {} # return object

    for line in handle:
        # parse tags
        if line.startswith('['):
            tag = line.strip('\n').rstrip(',').strip('[]')
            if tag == 'Data':
                get_header = True
            continue

        # else not a tag line, parse contents
        tokens = line.strip('\n').split(',')

        # process tokens according to current tag
        if tag == 'Header':
            # inside [Header] block
            key, value = tokens[:2]
            run_info.update({key: value})

        elif tag == 'Data':
            # inside [Data] block
            if not run_info.has_key('Data'):
                run_info.update({'Data': {}})

            if get_header:
                # parse the first line as the header row
                header = tokens
                if not 'Sample_Name' in header:
                    sys.stderr.write("ERROR: SampleSheet.csv Data header does not include Sample_Name")
                    sys.exit()
                get_header = False
                continue

            index1 = tokens[header.index('index')]
            index2 = tokens[header.index('index2')]

            # parse Sample_Name field
            filename = tokens[header.index('Sample_Name')]
            clean_filename = re.sub('[_.;]', '-', filename)
            clean_filename += '_S%d' % sample_number # should correspond to FASTQ filename

            run_info['Data'].update({clean_filename: {'index1': index1,
                                                      'index2': index2}})

            # parse Description field
            # FIXME: currently this is partitioned by sub-samples (semi-colon delimited)
            # FIXME: but the pipeline does not support differential handling of sub-samples

            desc = tokens[header.index('Description')]
            desc_fields = desc.split() # whitespace-delimited
            for desc_field in desc_fields:
                desc_field_label = desc_field.split(':')[0]
                desc_subfields = desc_field.replace(desc_field_label+':', '').split(';')

                for desc_subfield in desc_subfields:
                    desc_subfield_tokens = desc_subfield.split(':')
                    #sample_name = desc_subfield_tokens[0]
                    #clean_sample_name = re.sub('[_.]', '-', sample_name)

                    if desc_field_label == 'Research':
                        is_research = (desc_subfield_tokens[1] == 'TRUE')
                        # would have been keying by clean_sample_name here and below
                        run_info['Data'][clean_filename].update({
                            'research': is_research,
                            'is_T_primer': 'TPRIMER' in desc_subfield_tokens})

                    if desc_field_label == 'Disablecontamcheck':
                        run_info['Data'][clean_filename].update({'disable_contamination_check': (desc_subfield_tokens[-1]=='TRUE')})

                    if desc_field_label == 'Comments':
                        run_info['Data'][clean_filename].update({'comments': desc_subfield_tokens[-1]})

                if not run_info['Data'][clean_filename].has_key('disable_contamination_check'):
                    # 'Disablecontamcheck' subfield not always under Description
                    # default to contamination check active
                    run_info['Data'][clean_filename].update({'disable_contamination_check': False})

                # FIXME: for the time being, apply ONLY first sub-sample description to entire sample
                break

            sample_number += 1

        else:
            # ignore other tags
            pass
    return run_info