#!/usr/bin/env python
#
# convert cython wrapping code to nanobind wrapping code
#
# usage:  ./cython2nanobind.py simplex_tree.pxd > simplex_tree.cc
#
#


import re
import sys


def generate_nanobind_bindings(cython_code):
    # Extract all method declarations
    method_pattern = re.compile(
        r'\s*([^\(]+)\(([^)]*)\)(?:.*nogil)?(?:.*except\s\+)?'
    )

    # List to hold the generated nanobind code
    nanobind_code = []

    for line in cython_code.split('\n'):
        # Skip empty lines or comments
        if not line.strip() or line.strip().startswith('//') or 'cdef cppclass' in line:
            if 'cdef cppclass' in line:
                print(line)
                match = re.match(r'^\s*cdef cppclass .*?"(.+?)".*$', line)
                print(match)
                if match:
                    class_name = match.group(1)
            continue

        # Extract method information
        match = method_pattern.match(line)
        if not match:
            print(f"Could not parse line: {line}")
            continue

        return_type, params = match.groups()
        method_name = return_type.split()[-1]

        # Process parameters for nb::arg()
        param_ingredients = []
        for param in re.findall(r'([^\s]+)\s+([^\s,)]+)', params):
            type_, name = param
            if 'vector[' in type_:
                type_ = type_.replace('vector[', 'std::vector<').replace(']', '>')
            param_ingredients.append(f'nb::arg("{name}"),')

        # Create parameter string for nb::arg()
        param_str = " ".join(param_ingredients) if param_ingredients else ""
        if (param_str): param_str = "          " + param_str + "\n"

        # Generate documentation string
        docstring = 'R"pbdoc(TODO)pbdoc"'

        # Generate nanobind binding code
        binding_code = f'    .def("{method_name}",\n\
          &{class_name}::{method_name},\n{param_str}\
          {docstring})'
        nanobind_code.append(binding_code)

    return '\n'.join(nanobind_code)

if __name__ == "__main__":
    # filename = sys.argv[1:] # Get all arguments as files
    if len(sys.argv) != 2:
        print("Error: must provide a single filename")
        sys.exit(1)

    cppclass_regexp = re.compile(
        r'^\s*cdef\scppclass'
    )
    ignore_regexp = re.compile(
        r'^\s*pass#.*$'
    )
    empty_regexp = re.compile(
        r'^\s*pass#?.*$'
    )

    cython_code = None
    cpt = 0
    with open(sys.argv[1], "r") as file:
        for line in file:
            cpt += 1
            line.strip()

            # new code bloc
            if re.match(cppclass_regexp, line):
                #print(f"MATCH {cpt}: {line}")

                # old bloc must be translated
                if cython_code:
                    bindings = generate_nanobind_bindings(cython_code)
                    print(f"{bindings} \n")
                    #print(f"TRANSLATE {cpt}")
                    cyton_code = None

                # and a new bloc must be recorded
                #print(f"START {cpt}: {line}")
                cython_code = line
                continue

            # ignore empty_lines and lines contening 'pass'
            if re.match(ignore_regexp, line) or re.match(empty_regexp, line):
                continue

            else:
                # append the new line to the current bloc
                if cython_code:
                    #print(f"APPEND {cpt}: {line}")
                    cython_code += line

        # must translate the final bloc if necessary
        if cython_code:
            bindings = generate_nanobind_bindings(cython_code)
            print(f"{bindings} \n")
            #print(f"FINAL TRANSLATE {cpt}")
            cyton_code = None

    sys.exit(0)
