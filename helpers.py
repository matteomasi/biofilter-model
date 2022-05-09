""" 
BIOFILTER-MODEL - Helper functions
==================================

Collection of classes:
1) Unit conversions: convert ppb to g/m3 and vice-versa
2) FileIO: import and filter the time series
3) Parameters: import parameters from .yml file

(c) Matteo M. 2022

"""

import pandas as pd 
import numpy as np
import os, fnmatch
import yaml

class Conversion():
    # Conversion  
    @staticmethod
    def ppbtog(mw = 100.0):
        """
        Convert ppb to g/m3. Returns the conversion factor.

        Parameters
        ----------
        mw : float, default 100.0
            Molar weight
        """
        m = mw/24.45*1e-6
        return m

    @staticmethod
    def gtoppb(mw = 100.0): 
        """
        Convert g/m3 to ppb. Returns the conversion factor.

        Parameters
        ----------
        mw : float, default 100.0
            Molar weight
        """
        m = 24.45*1e6/mw
        return m


class FileIO():
    def __init__(self, csv_file):
        self.csv_file = csv_file        # Path to the csv file
        self.c_name = 'tVOC'            # Default column name in the header (e.g., tVOC)
        self.sm_name = 'Smoothed'       # Column name of the smoothed data
        self.smooth = False             # Internal variable for storing smmothing preference
        self.dateformat = True          # Is True when the time column is a date; False if time is in hours
        self.mw = 100.0


    
    def readfile(self, smooth = False, n_samples = 10):
        self.data = pd.read_csv(self.csv_file)
        if self.dateformat: # If time is formatted 
            self.data.Time = pd.to_datetime(self.data.Time, dayfirst=True)
            self.data = self.data.set_index('Time')
            self.data['t_diff'] = self.data.index - self.data.index[0]     # Calculate time difference
            self.data['t_s'] = self.data['t_diff'].dt.total_seconds()      # Time in seconds
        else:
            self.data['t_s'] = self.data.Time*3600 # Convert hours to seconds


        # Smooth data (moving window)
        if smooth:
            self.smooth = True
            self.data[self.sm_name] = self.data[self.c_name].rolling(n_samples, center=True, min_periods=1).mean()   # min_periods=1 avoid nan at the beginning and at the end of the ts

        return self.data

    def interp(self,ts):
        # Interpolate according to new time vector 
        if self.smooth:
            self.data_int = np.interp(ts,self.data['t_s'],self.data[self.sm_name])
        else:
            self.data_int = np.interp(ts,self.data['t_s'],self.data[self.c_name])
        return self.data_int


class Parameters():
    def __init__(self):
        self.filename = 'parameters/parameters.yml'
        self.directory = 'parameters/'
    
    # List .yml files in the specified directory
    def listfiles(self,directory):
        self.directory = directory
        listOfFiles = os.listdir(self.directory)
        pattern = "*.yml"
        l = []
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, pattern):
               l.append(entry)
        if len(l) < 1:
            print('Error. No .yml file found.')
        return l
    
    # Read .yml file and return the content as a dictionary
    def yamlread(self,f_name):
        self.filename = f_name
        with open(self.filename, 'r') as f:
            params = yaml.load(f, Loader=yaml.FullLoader)
            return params
    
    # Read .yml file and return raw content
    def yamlraw(self,f_name):
        self.filename = f_name
        with open(self.filename,'r') as f:
            raw = f.read()
            return raw

