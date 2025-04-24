# Glider_AZFP_Prototyping.py
# Converting the raw data from the glider AZFP to netCDF format using echopype

# Created by Jace Marquardt
# Last updated: 2025-04-23

import echopype as ep
import echopype.visualize as epviz
import glob
import numpy as np
import os

import warnings
warnings.filterwarnings('ignore')

def echopype_processing():
    """
    This function processes the raw echogram data from the glider AZFP and converts it to netCDF format.
    """
    
    data_directory = r"C:\Users\marqjace\OneDrive - Oregon State University\Desktop\Python\azfp\data\may_2023\to_process"
    output_directory = r"C:\Users\marqjace\OneDrive - Oregon State University\Desktop\Python\azfp\data\may_2023\processed"

    xml_file = 'tweaked.xml'

    file_list = glob.glob(os.path.join(data_directory, '*.01?'))
    file_list.sort()

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Convert the list of .01X files using echopype and save the output as netCDF files
    ed_list = []
    for raw_file in file_list:
        # Load the raw data into an xarray dataset
        ed = ep.open_raw(raw_file, sonar_model='AZFP', xml_path=xml_file)

        # Update some of the platform metadata
        ed.platform.attrs['platform_name'] = 'Glider 592'
        ed.platform.attrs['platform_type'] = 'Sub-Surface Glider'
        ed.platform.attrs['platform_code_ICES'] = '27'  # ICES code: sub-surface gliders
                                                        # (buoyancy-based propulsion capable
                                                        # at variable depths not constrained 
                                                        # to be near sea surface)

        # change the Sv_offset values from NaN to 0's
        ed.vendor['Sv_offset'] = ed.vendor['DS'] * 0.0

        # save the converted data to a SONAR-netCDF4 data file
        ed.to_netcdf(os.path.join(output_directory, os.path.split(raw_file)[1] + '.raw.nc'))
        ed_list.append(ed)

    # combine the data into a single dataset
    ds = ep.combine_echodata(ed_list)

    # setup values that can be used to spoof a temperature and salinity values. would want to pull 
    # this information from the glider CTD record if we were doing this processing for real.
    temperature = np.array([15, 7])
    salinity = np.array([32, 34])

    # calculate Sv for each file and save the results per file to disk for further analysis
    for ds in ed_list:
        time_record = np.array([ds.environment.time1.values[0], ds.environment.time1.values[-1]]).astype(float)
        t = np.interp(ds.environment.time1.values.astype(float), time_record, temperature)
        s = np.interp(ds.environment.time1.values.astype(float), time_record, salinity)
        env_params = {
            'pressure': 0.0,  # will always set to 0, replace after the fact
            'temperature': np.mean(t),
            'salinity': np.mean(s)
        }
        ds_sv = ep.calibrate.compute_Sv(ds, env_params=env_params)
        path, file = os.path.split(ds.provenance.source_filenames.values[0])
        ds_sv.to_netcdf(os.path.join(output_directory, file + '.proc.nc'))