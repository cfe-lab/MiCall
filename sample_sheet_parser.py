#! /usr/bin/env python

import re

def sample_sheet_parser (handle):
    """
    Parse the contents of SampleSheet.csv, convert contents into a
    Python dictionary object.
    Samples are tracked by sample name and sample number (e.g., S9).
    This is to distinguish replicates of the same sample.
    """

    tag = None
    get_header = False
    header = [] # store Data block column labels
    sample_number = 1 # 1-indexing
    run_info = {} # return object
    sample_sheet_version = None #Version 1 is there we used underscores and semicolons, version 2 we use tilde's and hashes
    sample_delimiter_v1 = ';' # Between multiple samples in one well
    project_delimiter_v1 = '_' # Between the sample and project names
    sample_delimiter_v2 = '---'
    project_delimiter_v2 = '--'
    sample_delimiter = sample_delimiter_v1
    project_delimiter = project_delimiter_v1

    for line in handle:
        stripped_line = line.rstrip('\n,')
        # parse tags
        if stripped_line.startswith('['):
            tag = stripped_line.strip('[]')
            if tag == 'Data':
                get_header = True
            continue # pragma: no cover (because of optimizer)

        # else not a tag line, parse contents
        tokens = stripped_line.split(',')

        # process tokens according to current tag
        if tag == 'Header':
            # inside [Header] block
            key, value = tokens[:2]
            run_info.update({key: value})

        elif tag == 'Data':
            # inside [Data] block
            if not run_info.has_key('Data'):
                run_info.update({'Data': {}})
            if not run_info.has_key('DataSplit'):
                run_info.update({'DataSplit': []})

            if get_header:
                # parse the first line as the header row
                header = tokens
                if not 'Sample_Name' in header:
                    raise ValueError("sample sheet data header does not include Sample_Name")
                get_header = False
                continue

            index1 = tokens[header.index('index')]
            index2 = tokens[header.index('index2')]

            # parse Sample_Name field
            filename = tokens[header.index('Sample_Name')]
            if(sample_sheet_version == None and 
               (sample_delimiter_v2 in filename or 
                project_delimiter_v2 in filename)):
                
                sample_sheet_version = 2
                sample_delimiter = sample_delimiter_v2
                project_delimiter = project_delimiter_v2
            elif(sample_sheet_version == None and
                 (sample_delimiter_v1 in filename or
                  project_delimiter_v1 in filename)):
                
                sample_sheet_version = 1
                
            clean_filename = re.sub('[_.;]', '-', filename)
            clean_filename += '_S%d' % sample_number # should correspond to FASTQ filename

            # July 9, 2014: also want to keep track of the *original*
            # Sample_Name as it appeared in the sample sheet.
            # July 10, 2014: clean up to fix some naming problems.
            run_info['Data'].update(
                {clean_filename: {
                    "index1": index1,
                    "index2": index2,
                    "comments": "",
                    "disable_contamination_check": False,
                    "research": True,
                    "chemistry": run_info["Assay"],
                    "orig_sample_name": filename
                }})
            run_info['sample_sheet_version'] = sample_sheet_version

            # Parse Description field.  This uses some version-specific
            # code to handle version 1 (where semicolons and underscores were used)
            # to version 2 (where tildes and hashes are used).

            desc = tokens[header.index('Description')]
            desc_fields = desc.split() # whitespace-delimited
            for desc_field in desc_fields:
                desc_field_label = desc_field.split(':')[0] #research/chemistry/comments/disable_contam_check
                desc_subfields = desc_field.replace(desc_field_label+':', '').split(sample_delimiter)
                
                for desc_subfield in desc_subfields:
                    desc_subfield_tokens = desc_subfield.split(':')
                    if(sample_sheet_version == 2): #for compatibility
                        sample, project, value = desc_subfield.split(project_delimiter)
                        desc_subfield_tokens = [
                            sample + project_delimiter + project, 
                            value]
                    
                    if desc_field_label == 'Research':
                        is_research = (desc_subfield_tokens[1] == 'TRUE')
                        # would have been keying by clean_sample_name here and below
                        run_info['Data'][clean_filename].update({
                            'research': is_research,
                            'is_T_primer': 'TPRIMER' in desc_subfield_tokens})

                    if desc_field_label == 'Disablecontamcheck':
                        run_info['Data'][clean_filename].update({'disable_contamination_check': (desc_subfield_tokens[-1]=='TRUE')})

                    if desc_field_label == 'Chemistry':
                        run_info['Data'][clean_filename].update({'chemistry': desc_subfield_tokens[-1]})

                    if desc_field_label == 'Comments':
                        run_info['Data'][clean_filename].update({'comments': desc_subfield_tokens[-1]})

                    # FIXME: for the time being, apply ONLY first sub-sample description to entire sample
                    break

            
            #Okay, lets populate the datasplit based on the data entry. 
            for sampproj in filename.split(sample_delimiter):
                # July 11, 2014: we don't want to change the above too much in fear
                # of breaking backwards-compatibility, but here we can insert some
                # more stuff into the dictionary.

                # Start with a copy of the dictionary we defined above and punch it up with
                # some more details.
                entry = dict(run_info["Data"][clean_filename])
                
                tmp = sampproj.split(project_delimiter)
                # We need to update some parts of this.  In particular note that
                # the data is being split up by sample when there are multiple
                # samples in the same row.
                entry.update({'sample' : tmp[0]})
                entry.update({'project' : tmp[1]})
                entry.update({'filename' : clean_filename})
                entry.update({'sample_number' : 'S%d' % sample_number })

                # Now, because we've split by sample, we must also split the stuff in
                # the description field up and change the values in entry
                # appropriately.
                for desc_field in desc_fields:
                    name,value,tmp = None,None,None
                    if(sample_sheet_version == 1):
                        name = desc_field.split(':')[0] #slice #actually this is wrong...
                        value = desc_field.replace(name + ':','')
                        tmp = value.split(sample_delimiter_v1)
                    elif(sample_sheet_version == 2):
                        name,value = desc_field.split(':')
                        tmp = value.split(sample_delimiter_v2)
                    
                    for elem in tmp:
                        samp, proj, val = None, None, None
                        if(sample_sheet_version == 1):
                            sj, val = elem.split(':')
                            samp, proj = sj.split(project_delimiter_v1)
                        elif(sample_sheet_version == 2):
                            samp, proj, val = elem.split(project_delimiter_v2)
                        
                        if(samp == entry['sample'] and proj == entry['project']):
                            if(name == 'Research'):
                                entry['research'] = (val == 'TRUE')
                            elif(name == 'Comments' and val != None):
                                entry['comments'] = val 
                            elif(name == 'Chemistry'):
                                entry['chemistry'] = val
                            elif(name == 'Disablecontamcheck'):
                                entry['disable_contamination_check'] = (val == 'TRUE')
                
                run_info['DataSplit'].append(entry)
            
            sample_number += 1
    return run_info

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Read a sample sheet")
    parser.add_argument("samplesheet")
    args = parser.parse_args()

    with open(args.samplesheet, "rb") as f:
        ss = sample_sheet_parser(f)
        print(ss)

if __name__ == "__main__":
    main()