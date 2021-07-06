import os
import re
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
print(sys.path)

import constants

def parse_domain(domain):
    if domain != "gene":
        if domain in constants.domains:
            return constants.domains[domain]
        elif re.search("^\d+$", domain):
            print(f"Attempting to parse {domain} as CDD PSSM-Id")
        else:
            print(f"Unexpected domain input: {domain}")
            print(f"Please provide a PSSM-Id (numeric) or one of the following supported domain names: {','.join(constants.domains)}")
            sys.exit(1) ## 1: invalid args
    return domain

def get_homologue(genes, out_dir, out_pref):
    pass

def get_reference():
    print("booya")
    return
    # pass
