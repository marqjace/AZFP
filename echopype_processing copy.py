# echopype_processing.py
# Converting the raw data from the glider AZFP to netCDF format using echopype

# Created by Jace Marquardt
# Last updated: 2025-04-23

import echopype as ep
import echopype.visualize as epviz
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

import warnings
warnings.filterwarnings('ignore')

def convert_raw(data_directory, output_directory, xml_file):
    """
    This function takes the raw echogram data from the glider AZFP, calibrates it using the XML file, and converts it to netCDF format.

    :param data_directory: str, path to the directory containing the raw echogram data files
    :param output_directory: str, path to the directory where the converted netCDF files will be saved
    :param xml_file: str, path to the XML file containing the sonar configuration

    :return: Converted and calibrated bioacoustic sonar data in a xarray dataset object
    """

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

    return ed_list


def glider_process(data_directory, output_directory, xml_file, glider_temp, glider_salinity, glider_depth):
    """
    This function takes the converted and calibrated echogram data from the glider AZFP, processes it, and generates echograms.
    It also saves the echograms as PNG files in the specified figures directory within '/processed'.

    :param data_directory: str, path to the directory containing the raw echogram data files
    :param output_directory: str, path to the directory where the converted netCDF files will be saved
    :param xml_file: str, path to the XML file containing the sonar configuration
    :param glider_temp: numpy array, temperature values for the glider
    :param glider_salinity: numpy array, salinity values for the glider
    
    :return: Saved echograms as PNG files in the figures directory
    """

    figures_directory = os.path.join(output_directory, 'figures/')
    if not os.path.isdir(figures_directory):
        os.mkdir(figures_directory)

    ed_list = convert_raw(data_directory, output_directory, xml_file)

    # setup values that can be used to spoof a temperature and salinity values. would want to pull 
    # this information from the glider CTD record if we were doing this processing for real.
    temperature = glider_temp
    salinity = glider_salinity
    depth = glider_depth

    # calculate Sv for each file and save the results per file to disk for further analysis
    for ds in ed_list:
        # Create a time_record that matches the length of glider_depth
        time_record = np.linspace(
            ds.environment.time1.values[0].astype(float), 
            ds.environment.time1.values[-1].astype(float), 
            len(glider_depth)
        )

        # Interpolate temperature and salinity to match the length of glider_depth
        temperature_interp = np.linspace(temperature[0], temperature[-1], len(glider_depth))
        salinity_interp = np.linspace(salinity[0], salinity[-1], len(glider_depth))

        # Debugging: Print time_record and ds.environment.time1.values
        print("time_record:", time_record)
        print("ds.environment.time1.values:", ds.environment.time1.values)

        # Interpolate depth over the full range of time_record
        t = np.interp(ds.environment.time1.values.astype(float), time_record, temperature_interp)
        s = np.interp(ds.environment.time1.values.astype(float), time_record, salinity_interp)
        d = np.interp(ds.environment.time1.values.astype(float), time_record, glider_depth)

        # Debugging: Print glider_depth and interpolated depth
        print("glider_depth:", glider_depth)
        print("Interpolated depth (d):", d)

        env_params = {
            'pressure': 0.0,  # will always set to 0, replace after the fact
            'temperature': np.mean(t),
            'salinity': np.mean(s)
        }
        ds_sv = ep.calibrate.compute_Sv(ds, env_params=env_params)
        path, file = os.path.split(ds.provenance.source_filenames.values[0])
        ds_sv.to_netcdf(os.path.join(output_directory, file + '.proc.nc'))  # Save the processed data as .proc.nc

        if 'echo_range' in ds_sv:
            echo_range = ds_sv['echo_range']
            # Ensure depth is reshaped to match echo_range dimensions
            depth_reshaped = np.broadcast_to(d, echo_range.shape)

            print("time_record:", time_record)
            print("ds.environment.time1.values:", ds.environment.time1.values)
            print("glider_depth:", glider_depth)
            print("Interpolated depth (d):", d)
            print("Length of time_record:", len(time_record))
            print("Length of glider_depth:", len(glider_depth))

            ds_sv['echo_range'] = echo_range + depth_reshaped

        ds_sv_clean = ep.preprocess.remove_noise(ds_sv, range_sample_num=30, ping_num=5) # remove noise from the Sv data

        time_min = ds.platform.time2.values[0]
        time_max = ds.platform.time2.values[-1]

        # Create a new directory for each file
        file_directory = os.path.join(figures_directory, file)
        if not os.path.exists(file_directory):
            os.mkdir(file_directory)

        echograms = epviz.create_echogram(ds_sv_clean, get_range=True, robust=True, vmin=-90, vmax=-50)
        frequency_labels = ["67_kHz", "120_kHz", "200_kHz"]  

        for i, echogram in enumerate(echograms):
            frequency_label = frequency_labels[i]
            echogram.fig.suptitle(f'{time_min} - {time_max} ({frequency_label})')
            echogram.fig.savefig(os.path.join(file_directory, f'{file}_{frequency_label}.png'), dpi=300, bbox_inches='tight')
            plt.close(echogram.fig)


data_directory = r"C:\Users\marqjace\OneDrive - Oregon State University\Desktop\Python\azfp\data\may_2023\to_process"
output_directory = r"C:\Users\marqjace\OneDrive - Oregon State University\Desktop\Python\azfp\data\may_2023\processed"
xml_file = 'tweaked.xml'

glider_temp = np.array([15, 7])
glider_salinity = np.array([32, 34])
glider_depth = np.arange(0, 251, 1)

# gdata = scipy.io.loadmat(r"C:\Users\marqjace\OneDrive - Oregon State University\Desktop\Python\azfp\data\may_2023\WA_202305241820-deployment_osu592_pass3.mat")
# glider_temp = gdata['Temp']
# glider_salinity = gdata['Salt']

glider_process(data_directory, output_directory, xml_file, glider_temp, glider_salinity, glider_depth)

