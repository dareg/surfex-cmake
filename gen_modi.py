#!/usr/bin/env python3

import sys
import re
from pathlib import Path

pattern_sub_name = re.compile(r"^\s*SUBROUTINE (\w+)", re.IGNORECASE)
pattern_end_of_decl = re.compile(r"^\s*if", re.IGNORECASE)


def gen_modi(f90):
    src = open(f90, "r")
    modi_path = f90.parent / Path("modi_" + str(f90.name))
    modi = open(modi_path, "w")
    sub_name = ""
    for line in src:
        if not sub_name:
            matches = pattern_sub_name.match(line)
            if matches:
                sub_name = matches.group(1)
                modi.write(f"MODULE MODI_{sub_name}\nCONTAINS\n")

        matches = pattern_end_of_decl.match(line)
        if matches:
            modi.write(f"END MODULE MODI_{sub_name}\n")
            return
        else:
            modi.write(line)


for file in sys.argv[1:]:
    gen_modi(Path(file))
