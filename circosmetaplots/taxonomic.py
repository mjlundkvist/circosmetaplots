# -*- coding: utf-8 -*-
"""
Created on Mon May 29 08:21:08 2017

@author: Moa
"""
from .weatherdata import WeatherData
from .otus import OTUreads
from .circosconf import CircosConf
import subprocess, sys, os, shutil

class TaxonomyPlot(object):
    def __init__(self, data_files, weather_file, place, year, level, plot_level):
        self.place, self.year, self.start_level, self.plot_level = place, str(year), level, plot_level
        self.OTU=OTUreads(data_files)               
        self.weather=WeatherData(weather_file, self.year, False)
        self.conf=CircosConf()
        self.perl=shutil.which('perl')
        self.folder_name=self.create_folder()        

    def create_folder(self):
        """
        Creates new (and uniquely named) folder which is also made the work directory,
        so that all text files will be collected at one place.
        """
        cur_dir=os.getcwd()
        unique=False
        dirlist= [item for item in os.listdir(cur_dir) if os.path.isdir(os.path.join(cur_dir,item))]
        folder_name='taxonomy_{}_{}'.format(self.place,self.year)
        j=1
        while not unique:
            if folder_name in dirlist:
                folder_name='taxonomy_{}_{}({})'.format(self.place,self.year,str(j))
                j+=1
            else:
                unique=True
        new_folder=os.path.join(cur_dir,folder_name)
        os.mkdir(new_folder)
        os.chdir(new_folder)
        return folder_name
    
    def taxonomy_files(self):
        """
        Writes the text files with the otu information needed for the phylogenetic plot.
        Returns the number of destinct groups at the highest level of the plot (which will be
        equivalent to the ideograms in the circos plot.
        """
        location=self.place.capitalize()+'-'+str(self.year)+'-'
        no_of_ideograms=self.OTU.make_tree(location,self.start_level,self.plot_level)
        return no_of_ideograms
    
    def taxonomy_plot(self,seasons):
        """Runs all functions for the configuration and data files for the circos plot."""
        print('Formatting data.')
        no_of_ideograms=self.taxonomy_files()
        location=self.place.capitalize()+'-'+str(self.year)
        if seasons==True:
            seasons=self.weather.seasons(self.place)
        print('Done')
        self.conf.taxo_conf(no_of_ideograms, location, self.start_level, self.plot_level, seasons)
    
    def run(self,circos_location, seasons=True):
        """
        Runs all function needed to get the plot.
        """
        self.taxonomy_plot(seasons)
        print(seasons)
        filename=self.place+'_'+self.year+'.png'
        self.run_circos(filename, circos_location)
        print('In folder ', self.folder_name)

    def run_circos(self, output_image, circos_location):
        """
        Calls perl and circos to generate an image by name of the argument.
        """
        #perl = r"C:\Strawberry\perl\bin\perl.exe"
        #perl_script=r'C:\Users\Moa\Downloads\circos-0.69-5\bin\circos'
        params=' -outputfile', ' {}'.format(output_image)
        pl_script = subprocess.Popen([self.perl, circos_location, params])#,stdout=sys.stdout,
        pl_script.communicate()
        if pl_script.returncode!=0:
            print('Something went wrong with running circos from python')
        else:
            print('Circos image generated.')

