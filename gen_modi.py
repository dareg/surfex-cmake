#!/usr/bin/env python3

import sys
import re
from pathlib import Path

reg_module = re.compile(r"^\s*MODULE (\w+)", re.IGNORECASE)
reg_continuation_line = re.compile(r"&\s*$")
reg_sub_func_name = re.compile(r"^\s*(SUBROUTINE|FUNCTION)\s+(\w+)", re.IGNORECASE)
reg_starts_with_ifdef = re.compile(r"^#ifdef|#if\s+defined|#ifndef", re.IGNORECASE)

all_decl_regex = [
    r"^\s*$",
    r"^\s*!",
    r"^\s*SUBROUTINE ",
    r"^\s*FUNCTION ",
    r"^\s* ",
    r"^\s*INTEGER\s*[\(*,:\s]",
    r"^\s*REAL\s*[\(*,:\s]",
    r"^\s*DOUBLE PRECISION\s*[\(*,:\s]",
    r"^\s*CHARACTER\s*[\(*,:\s]",
    r"^\s*LOGICAL\s*[\(*,:\s]",
    r"^\s*TYPE\s*\(",
    r"^\s*USE ",
    r"^\s*#",
    r"^\s*IMPLICIT NONE",
    r"^\s*include ",
]

reg_is_still_decl = re.compile("|".join(x for x in all_decl_regex), re.IGNORECASE)


def is_modi_needed(f90):
    if f90.name.startswith("modi_"):
        return False
    src = open(f90, "r")
    for line in src:
        matches = reg_module.match(line)
        if matches:
            return False
    return True


def simplify_code(src):
    """Remove continuation lines and comments"""
    lines = []
    buf = ""
    for line in src:
        line = line.strip()
        if line.startswith("!"):
            continue
        if '!' in line:
            line=line.split('!')[0].strip()
        if line.endswith("&"):
            buf += line.replace("&", "")
        else:
            lines.append(buf + line.replace("&", "") + "\n")
            buf = ""
    return lines


def remove_dangling_ifdef(lines):
    ifdef = []
    for idx, line in enumerate(lines):
        if reg_starts_with_ifdef.match(line):
            ifdef.append(idx)
        if line.startswith("#endif"):
            ifdef.pop()
    for idx in reversed(ifdef):
        lines.pop(idx)


def gen_modi(f90):
    src = open(f90, "r")
    modi_path = f90.parent / Path("modi_" + str(f90.name))
    sub_name = ""
    func_or_sub = ""
    lines = simplify_code(src)
    ending = ""
    modi_lines = []
    for line in lines:
        if not sub_name:
            matches = reg_sub_func_name.match(line)
            if matches:
                sub_name = matches.group(2)
                ending = f"END {matches.group(1)} {sub_name}"
                modi_lines.append(f"MODULE MODI_{sub_name}\nCONTAINS\n")

        matches = reg_is_still_decl.match(line)
        if not matches:
            # print("end of decl", line)
            modi_lines.append(f"{ending}\n")
            modi_lines.append(f"END MODULE MODI_{sub_name}\n")
            break
        else:
            modi_lines.append(line)

    remove_dangling_ifdef(modi_lines)
    modi = open(modi_path, "w")
    modi.write("".join(modi_lines))


for arg in sys.argv[1:]:
    file = Path(arg)
    if is_modi_needed(file):
        gen_modi(file)
