This directory contains scripts which help create specialized
diagnostic (history) output lists for CMIP7 experiments.

The master list of CAM fields is in the file,
**master_fieldlist.txt**. This file can be regenerated from a CAM
logfile using the script in **extract_master_list.py**:

```
./extract_master_list.py <path to logfile> --write-names > master_fieldlist.txt
```

The script, **cmip_diagnostic_namelists.py** is responsible for
reading the CMIP7 data request spreadsheet and creating CAM history
namelist files. This script collects all the global atmosphere
diagnostic field requests in the spreadsheet and groups them by
requested frequency of output. A check is then performed to look for
and omit any fields which are not available in the CAM master field
list. The rest of the fields are output as various `fincl<n>`
lists. Finally, the script outputs some general history namelist
variables, e.g., `nhtfrq`.

To run the script, first download a current version of the [CMIP7 data
request
spreadsheet](https://docs.google.com/spreadsheets/d/1XUdCTl1zKsWi_yTvMZqsMtBnnUPIhIeTPVSLHgL0Hdo)
(File `==>` Download `==>` .csv). The script is then run as:

```
./cmip_diagnostic_namelists.py <path-to-downloaded-spreadsheet-file>
```

The full interface is:
```
usage: cmip_diagnostic_namelists.py [-h] [--namelist-file NAMELIST_FILE]
                                    [--overwrite] [--include-cosp]
                                    <path to CMIP7 data request file>

Script to read version of the CMIP7 data request spreadsheet, check for
any field requests which are not availble from the CAM CMIP7 model
configurations, and produce the requested diagnostic sections of CAM's
runtime namelist.
Note that the data request spreadsheet must be in CSV format

positional arguments:
  <path to CMIP7 data request file>

options:
  -h, --help            show this help message and exit
  --namelist-file NAMELIST_FILE
                        Path to write namelist file entries (Default: stdout)
  --overwrite           Overwrite namelist file if it exists
  --include-cosp        Include COSP diagnostic fields in output
```
