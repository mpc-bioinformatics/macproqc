#!/usr/bin/env python

import platform
import shlex
import argparse

def parse_ext_args(args_string):
    """Parse arguments supplied via nf-core's ext.args."""
    if args_string in (None, "", "null"):
        args_string = ""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-peptide_mass_tolerance_upper",
        default=5.0,
        type=float,
        required=False,
        help="Comet parameter: peptide_mass_tolerance_upper",
    )
    parser.add_argument(
        "-peptide_mass_tolerance_lower",
        default=-5.0,
        type=float,
        required=False,
        help="Comet parameter: peptide_mass_tolerance_lower",
    )
    parser.add_argument(
        "-peptide_mass_units",
        default=2,
        type=int,
        required=False,
        help="Comet parameter: peptide_mass_units",
    )
    parser.add_argument(
        "-isotope_error",
        default=2,
        type=int,
        required=False,
        help="Comet parameter: isotope_error",
    )
    parser.add_argument(
        "-fragment_bin_tol",
        default=0.02,
        type=float,
        required=False,
        help="Comet parameter: fragment_bin_tol",
    )
    parser.add_argument(
        "-fragment_bin_offset",
        default=0.0,
        type=float,
        required=False,
        help="Comet parameter: fragment_bin_offset",
    )
    parser.add_argument(
        "-theoretical_fragment_ions",
        default=0,
        type=int,
        required=False,
        help="Comet parameter: theoretical_fragment_ions",
    )
    parser.add_argument(
        "-psms_per_spectrum",
        default=5,
        type=int,
        required=False,
        help="Comet parameter: num_output_lines",
    )

    parser.add_argument(
        "-variable_modifications",
        nargs="?",
        const="",
        default="",
        type=str,
        required=False,
        help='Adjust variable  modifications, given in the Comet annotation and separated by a semicolon, e.g. "variable_mod01 = 15.9949 M 0 3 -1 0 0 0.0;variable_mod02 = 79.966331 STY 0 3 -1 0 0 97.976896"',
    )
    parser.add_argument(
        "-static_modifications",
        nargs="?",
        const="",
        default="",
        type=str,
        required=False,
        help='Adding static modifications, given in the Comet annotation and separated by a semicolon, e.g. "add_K_lysine = 8.014199;add_R_arginine = 10.008269"',
    )

    return parser.parse_args(shlex.split(args_string))


if __name__ == "__main__":
    args = parse_ext_args("${task.ext.args}")

    comet_config_template = "${comet_config_template}"
    params_out_file = "${params_out_file}"

    # split the variable modifications into a dictionary
    variable_modifications = {}
    if args.variable_modifications and args.variable_modifications.strip() not in {"''", '""'}:
        for line in args.variable_modifications.split(";"):
            line = line.strip()
            if not line:
                continue
            vals = line.split("=")
            if len(vals) == 2:
                key = vals[0].strip()
                value = vals[1].strip()
                variable_modifications[key] = value

    # decode the static modifications
    static_modifications = {}
    if args.static_modifications and args.static_modifications.strip() not in {"''", '""'}:
        for line in args.static_modifications.split(";"):
            line = line.strip()
            if not line:
                continue
            vals = line.split("=")
            if len(vals) == 2:
                key = vals[0].strip()
                value = float(vals[1].strip())
                static_modifications[key] = value

    with open(comet_config_template, 'r') as params_in, open(params_out_file, 'w') as params_out:
        # adjust and write the parameter file line-wise
        for line in params_in:
            line = line.strip()
            # the first few parameters are static / always used like this in macproqc
            if line.startswith("decoy_search"):                 # disable decoy search
                line = "decoy_search = 0"
            elif line.startswith("output_sqtfile"):             # disable SQT output
                line = "output_sqtfile = 0"
            elif line.startswith("output_txtfile"):             # disable TXT output
                line = "output_txtfile = 0"
            elif line.startswith("output_pepxmlfile"):          # disable pepXML output
                line = "output_pepxmlfile = 0"
            elif line.startswith("output_mzidentmlfile"):       # enable mzIdentML output
                line = "output_mzidentmlfile = 1"
            elif line.startswith("output_percolatorfile"):      # disable percolator output
                line = "output_percolatorfile = 0"
            elif line.startswith("max_duplicate_proteins"):     # export all proteins per PSM
                line = "max_duplicate_proteins = -1"
            elif line.startswith("equal_I_and_L"):              # I and L are not equal for us here
                line = "equal_I_and_L = 0"

            # now set the variable search parameters
            elif line.startswith("peptide_mass_tolerance_upper"):
                line = f"peptide_mass_tolerance_upper = {args.peptide_mass_tolerance_upper}"
            elif line.startswith("peptide_mass_tolerance_lower"):
                line = f"peptide_mass_tolerance_lower = {args.peptide_mass_tolerance_lower}"
            elif line.startswith("peptide_mass_units"):
                line = f"peptide_mass_units = {args.peptide_mass_units}"
            elif line.startswith("isotope_error"):
                line = f"isotope_error = {args.isotope_error}"
            elif line.startswith("fragment_bin_tol"):
                line = f"fragment_bin_tol = {args.fragment_bin_tol}"
            elif line.startswith("fragment_bin_offset"):
                line = f"fragment_bin_offset = {args.fragment_bin_offset}"
            elif line.startswith("theoretical_fragment_ions"):
                line = f"theoretical_fragment_ions = {args.theoretical_fragment_ions}"
            elif line.startswith("num_output_lines"):
                line = f"num_output_lines = {args.psms_per_spectrum}"

            # adjust variable  modifications
            for key, value in variable_modifications.items():
                if line.startswith(key):
                    line = f"{key} = {str(value)}"

            # set static modifications
            for key, value in static_modifications.items():
                if line.startswith(key):
                    line = f"{key} = {str(value)}"

            params_out.write(line + "\\n")

    # finally write out the versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {platform.python_version()}\\n")
