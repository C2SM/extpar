#! /usr/bin/env python3

import re
import os

files = [ "mo_soil_tg_fields.f90",
          "mo_glcc_tg_fields.f90", 
          "mo_landuse_tg_fields.f90",       
          "mo_isa_tg_fields.f90",      
          "mo_ahf_tg_fields.f90",      
          "mo_era_tg_fields.f90",      
          "mo_topo_tg_fields.f90",     
          "mo_aot_target_fields.f90",  
          "mo_python_target_fields.f90",  
          "mo_flake_tg_fields.f90"
        ]

leading_pattern = re.compile(r'^allocate *\(', re.I)
trailing_pattern = re.compile(r', *stat *= *errorcode\) *$', re.I)
             
for file in files:
    backup = file + ".bak"
    os.rename(file, backup)
    print("[INFO] {}".format(backup))
    with open(backup, "r") as f:
        lines = [line.rstrip() for line in f]
    with open(file, "w") as of:
        for line in lines:
            stripped_line = line.strip()
            if re.match(r'[\s]*allocate', line, re.I):
                print("-------------------------------------------------------------------------------")
                print("[INFO]", stripped_line)
                # strip leading 'allocate *(' and trailing ', *STAT=errorcode)'
                stripped_line = re.sub(leading_pattern, '', stripped_line)
                stripped_line = re.sub(trailing_pattern, '', stripped_line)                
                arrays = re.findall(r'([\w_%]+\([\w,:%]+\))', stripped_line)
                for array in arrays:
                    name = re.search(r'([\w_%]+)', array)[1]
                    dims = re.search(r'[\w_%]+\(([\w,:%]+)\)', array)[1]
                    #print("[INFO]", array, " -> ", name, re.sub(r'1:', '', dims))
                    
                    print("if (l_use_array_cache) then", file=of)
                    print("   call allocate_cached('{0}', {0}, [{1}])".format(name, re.sub(r'1:', '', dims)), file=of)
                    print("else", file=of)
                    print("   allocate({}({}), stat=errorcode)".format(name, re.sub(r'1:', '', dims)), file=of)
                    print("endif", file=of)                                        
            else:
                print(line, file=of)
