#! /bin/env python3

"""Script to read version of the CMIP7 data request spreadsheet, check for
any field requests which are not availble from the CAM CMIP7 model
configurations, and produce the requested diagnostic sections of CAM's
runtime namelist.
Note that the data request spreadsheet must be in CSV format
"""

import argparse
import contextlib
import csv
import os
import re
import sys

_FREQUENCY_COLNAME = "CMIP7 Freq."
_MODELTYPE_COLNAME = "Modelling Realm - Primary"
_REGION_COLNAME = "Region"
_CAM_DIAG_COLNAME = "NorESM3 name (dependency)"
_REQUIRED_HEADERS = [_FREQUENCY_COLNAME, _MODELTYPE_COLNAME, _REGION_COLNAME, _CAM_DIAG_COLNAME]
_HIST_FILEORDER = ['mon', 'day', '6hr', '3hr', '1hr', 'subhr']
_HIST_FRQCODES = {'mon':'0', 'day':'-24', '6hr':'-6', '3hr':'-3', '1hr':'-1', 'subhr':'1'}
_HIST_MFILT = {'mon':'1', 'day':'30', '6hr':'30', '3hr':'30', '1hr':'30', 'subhr':'30'}
_HIST_TITLES =  {'mon':'! monthly output', 'day':'! daily output', '6hr':'! 6-hourly output',
                 '3hr':'! 3-hourly output', '1hr':'! 1-hourly output',
                 'subhr':'! timestep output'}

def is_number(text):
    """Return True if <text> represents a literal numeric constant.
    Return False otherwise."""
    val = False
    try:
        flt = float(text)
        val = True
    except ValueError as verr:
        val = False
    # end try
    return val

def command_line(args):
    """Read the command line arguments (args) to retrieve the path to the
    CMIP7 data request spreadsheet.
    Return the path to the spreadsheet file."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("csv_file", metavar='<path to CMIP7 data request file>', type=str)
    parser.add_argument("--namelist-file", type=str, default="-",
                        help="Path to write namelist file entries (Default: stdout)")
    parser.add_argument("--overwrite", action='store_true', default=False,
                        help="Overwrite namelist file if it exists")
    parser.add_argument("--include-cosp", action='store_true', default=False,
                        help="Include COSP diagnostic fields in output")
    pargs = parser.parse_args(args)
    return pargs.csv_file, pargs.namelist_file, pargs.overwrite, pargs.include_cosp

@contextlib.contextmanager
def flex_open(filename=None, mode='w'):
    if filename and filename != '-':
        fh = open(filename, mode)
    else:
        fh = sys.stdout
    # end if

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def read_diagnostic_fieldnames(include_cosp=False):
    """Read the master list of CAM diagnostic (history) fieldnames from the
       saved master list.
    If <include_cosp> is True, include the COSP diagnostic fieldnames in the
       master list. Otherwise, keep a separate list of COSP fieldnames
    Return the list of fieldnames and the list of COSP fieldnames."""

    cam_diag_fieldname_file = os.path.join(os.path.dirname(__file__),
                                           "master_fieldlist.txt")
    all_fieldnames = []
    cosp_fieldnames = []
    with open(cam_diag_fieldname_file, mode='r') as infile:
        for line in infile:
            fieldnames = [x.strip() for x in line.split()]
            all_fieldnames.extend(fieldnames)
        # end for
    # end with
    # Handle COSP fieldnames separately
    cosp_diag_fieldname_file = os.path.join(os.path.dirname(__file__),
                                            "cosp_fieldlist.txt")
    with open(cosp_diag_fieldname_file, mode='r') as infile:
        for line in infile:
            fieldnames = [x.strip() for x in line.split()]
            if include_cosp:
                all_fieldnames.extend(fieldnames)
            else:
                cosp_fieldnames.extend(fieldnames)
            # end for
        # end with
    # end if
    return all_fieldnames, cosp_fieldnames

def parse_spreadsheet(csvfile, model_name="atmos"):
    """Parse <csvfile> and return a dictionary of the requested CAM fields at
    different output frequencies.
    The dictionary keys are the frequency and the value is a list of
    fieldnames."""
    cmip_dict = {}
    with open(csvfile, mode='r', newline="") as infile:
        reader = csv.reader(infile)
        headers = next(reader)
        # Create a dictionary with the column number for each required column
        col_dirs = {}
        for colnum, col in enumerate(headers):
            if col in _REQUIRED_HEADERS:
                if col in col_dirs:
                    emsg = (f"Duplicate column entry, '{col}', in columns "
                            f"{col_dirs[col]} and {colnum}")
                    raise ValueError(emsg)
                # end if
                col_dirs[col] = colnum
            # end if
        # end for
        if len(col_dirs) != len(_REQUIRED_HEADERS):
            missing = ', '.join(set(_REQUIRED_HEADERS) - set(col_dirs.keys()))
            raise ValueError(f"Missing headers: {missing}")
        # end if
        rownum = 1
        freq_col = col_dirs[_FREQUENCY_COLNAME]
        model_col = col_dirs[_MODELTYPE_COLNAME]
        region_col = col_dirs[_REGION_COLNAME]
        name_col = col_dirs[_CAM_DIAG_COLNAME]
        for row in reader:
            rownum += 1
            if row[model_col] != model_name:
                continue
            # end if
            if (row[region_col] != "GLB") and row[name_col].strip():
                print(f"Field {row[name_col]} on row {rownum} has region, "
                      "{row[region_col]},  skipping")
            else:
                # First, make sure there is a dictionary entry for this frequency
                if row[freq_col] not in cmip_dict:
                    cmip_dict[row[freq_col]] = []
                # end if
                names = [x.strip() for x in re.split(r'[+/,*-]', row[name_col])
                         if x.strip() and (not is_number(x.strip()))]
                cmip_dict[row[freq_col]].extend(names)
            # end if
        # end for
    # end with
    # Cleanup each request to remove duplicates and sort
    for freq in cmip_dict:
        cmip_dict[freq] = sorted(set(cmip_dict[freq]))
    # end for
    return cmip_dict

def check_for_missing_fieldnames(masterlist, data_request):
    """Given a data request dictionary (<data_request>),
    check to see if any are not in <masterlist>.
    Return a list of missing fields names (an empty list means none).
    Clean <data_request> to remove missing field entries (side effect)."""
    # Gather the set of all fields (combine different frequencies)
    all_reqfields = set()
    for fields in data_request.values():
        all_reqfields |= set(fields)
    # end for
    # Any fields not in <masterlist> are missing
    missing = all_reqfields - set(masterlist)
    # Remove missing fields from data_request
    for key, value in data_request.items():
        data_request[key] = [x for x in value if x not in missing]
    # end for
    return sorted(missing)

def generate_namelist_entries(data_request, nl_filename, maxline=125, hist_files=_HIST_FILEORDER):
    """Write the set of namelist entries represented by <data_request> to
    <nl_filename> (which may be standard output)."""
    lbreak = ''
    with flex_open(nl_filename, mode="w") as outfile:
        for index, freq in enumerate(hist_files):
            if freq in data_request:
                if freq == 'subhr':
                    avgflag = 'I'
                else:
                    avgflag = 'A'
                # end if
                # Write history file config info
                outfile.write(f"{lbreak}{_HIST_TITLES[freq]}\n")
                outfile.write(f"nhtfrq({index + 1}) = {_HIST_FRQCODES[freq]}\n")
                outfile.write(f"mfilt({index + 1}) = {_HIST_MFILT[freq]}\n")
                outfile.write(f"empty_htapes({index + 1}) = .true.\n")
                fields = data_request[freq]
                fldstring = ', '.join([f"{x}:{avgflag}" for x in fields])
                nlstr = f"fincl{index + 1} = {fldstring}"
                # Write the fincl string with appropriate line breaks
                begpos = 0
                strlen = len(nlstr)
                while begpos < strlen:
                    endpos = strlen
                    if endpos - begpos > maxline:
                        endpos = nlstr[0:begpos + maxline].rfind(' ')
                        if endpos < begpos:
                            endpos = strlen
                        # end if
                    # end if
                    outfile.write(f"{nlstr[begpos:endpos]}\n")
                    begpos = endpos
                # end while
            # end if
            lbreak = '\n'
        # end for
        # Now, write out combined namelist items
    # end with


###############################################################################

if __name__ == "__main__":
    csvfile, nl_filename, overwrite, include_cosp = command_line(sys.argv[1:])
    if not overwrite and os.path.exists(nl_filename):
        raise ValueError(f"namelist file, '{nl_filename}', exists, aborting")
    # end if
    all_fieldnames, cosp_fieldnames = read_diagnostic_fieldnames(include_cosp)
    data_request = parse_spreadsheet(csvfile)
    missing = check_for_missing_fieldnames(all_fieldnames, data_request)
    # Separate the COSP fields from the others
    cosp_missing = set(missing) & set(cosp_fieldnames)
    missing = set(missing) - cosp_missing
    if cosp_missing:
        print(f"The following {len(cosp_missing)} fields will not be included "
              "because they are from the COSP module:")
        for field in sorted(cosp_missing):
            print(f"  {field}")
        # end for
    # end if
    if missing:
        print(f"The following {len(missing)} fields are not output from CAM:")
        for field in sorted(missing):
            print(f"  {field}")
        # end for
    # end if
    generate_namelist_entries(data_request, nl_filename, maxline=80)
    sys.exit(0)
