# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:43:57 2017

@author: Moa
"""
from .reads import YearReads
from .circosconf import CircosConf
from .weatherdata import WeatherData
import subprocess, sys, os
from .node import Node_y
import shutil

class YearPlot(object):
    def __init__(self, data_files, location, year, weather_file, wind_file):
        self.location_id=location+'-'+str(year)+'-'
        self.conf=CircosConf()
        self.reads=YearReads(data_files, self.location_id)
        self.weather=WeatherData(weather_file, year, wind_file)
        self.root=Node_y(None, 'root', None)
        self.perl=shutil.which('perl')
        self.circos_location=None
        self.place=location
        self.read_totals=None
        self.stop_rank=None
        self.folder,self.path=self.create_folder()
        
    def run_single_mix(self,taxa, circos_location,seasons=True):
        """
        Makes a single plot year plot of a mix of taxa.
        Doesn't handle unspecifieds or incertae sedis.
        Argument: a list of tuples/lists with rank and taxa for each to be plotted.
        """
        self.circos_location=circos_location
        self.weathers(seasons)
        read_totals=self.reads.total_by_week()
        self.conf.year_base(read_totals)
        filenames,maximums=self.reads.reads_single(taxa)
        level,name='',''
        self.conf.year_conf_s(filenames, maximums, self.location_id, level, name, seasons)
        self.run_circos('circos.conf','{}.png'.format(self.location_id))
        
    def run_single(self,level,name, circos_location,seasons=True):
        """
        Makes a single plot, of all direct descendants of the given taxa. 
        Takes the taxa and its rank as arguments
        """
        self.circos_location=circos_location
        print('Formatting data.')
        self.weathers(seasons)
        filenames,maximums=self.reads.reads_single2(level,name)
        read_totals=self.reads.total_by_week()
        self.conf.year_base(read_totals)
        self.conf.year_conf_s(filenames, maximums, self.location_id, level, name,seasons)
        print('Calling Circos.')
        self.run_circos('circos.conf','{}_{}.png'.format(level,name))
        print('Plot in folder ', self.folder)
    
    def run(self,stop_rank, circos_location,seasons=True):
        """
        The hierarchal run of the year plot. Function that start the whole folder tree where all the ranks (until the stop) gets plotted.
        """
        self.circos_location=circos_location
        print('start folder:',self.folder)
        print('Starting plotting.')
        start_rank='k'
        if len(stop_rank)>1:
            self.stop_rank=stop_rank[0]
        else:
            self.stop_rank=stop_rank
        self.starting(start_rank,seasons)
        path=''
        self.layer(self.root,path)
        print('Done')
        print('Plots in folder ', self.folder)
    
    def starting(self,start_rank,seasons):
        """
        Writes the weather files that don't change for all the plots, plus the read files for the very first.
        """
        self.weathers(seasons)
        self.read_totals=self.reads.total_by_week()
        self.conf.year_base(self.read_totals)
        self.conf.year_conf_background(seasons)
        self.run_circos('background.conf','background.png')
        labels,maximums=self.reads.uppermost_layer(start_rank,self.root)
        path=''
        name=' '
        self.conf.year_conf(labels, maximums, path, self.location_id, start_rank, name)
        self.run_circos('circos.conf','{}_rank.png'.format(self.root.children[0].level))
        #self.remove_files(path,name,labels)
        
    def layer(self,node,path):
        """
        Recursive function that continues down the clades of the highest level, 
        making folders, extracting data and making plots for all descendants until the stop level.
        Toggled comment of run_circos.
        """
        for child in node.children:
            name=child.name
            folder_name=child.level+'_'+name.replace(' ','_')
            new_path=os.path.join(path,folder_name)
            os.mkdir(new_path)
            labels,maximums,next_letter=self.reads.lower_layer(child, new_path)
            self.conf.year_conf(labels, maximums, new_path, self.location_id, child.level, child.name)
            self.run_circos(os.path.join(new_path,'circos.conf'),os.path.join(new_path,'{}_{}.png'.format(child.level, child.name.replace(' ','_'))))
            self.remove_files(new_path,name,labels)
            if next_letter!=self.stop_rank:
                self.layer(child,new_path)
                
    def weathers(self,seasons):
        """
        Extracts data and writes the files neccessary for the weather plots in the year circos plot.
        """
        self.weather.find_data('temperature.txt', self.place, 'Temperature', 'DD')
        self.weather.find_data('precipitation.txt', self.place, 'Precipitation', 'DD')
        self.weather.wind_track(self.place)
        if seasons==True:
            self.weather.seasons(self.place)
    
    def create_folder(self):
        """
        Creates new (and uniquely named) folder which is also made the work directory,
        so that all text files will be collected at one place.
        """
        cur_dir=os.getcwd()
        unique=False
        dirlist= [item for item in os.listdir(cur_dir) if os.path.isdir(os.path.join(cur_dir,item))]
        folder_name=self.location_id[:-1].replace('-','_')
        j=1
        while not unique:
            if folder_name in dirlist:
                folder_name=self.location_id[:-1].replace('-','_')+'('+str(j)+')'
                j+=1
            else:
                unique=True
        new_folder=os.path.join(cur_dir,folder_name)
        os.mkdir(new_folder)
        os.chdir(new_folder)
        return folder_name,new_folder  
    
    def run_circos(self,input_file, output_image):
        """
        Calls Circos and generates the plot.
        """
        params='-conf', ' {}'.format(input_file), ' -outputfile', ' {}'.format(output_image), ' -silent'
        pl_script = subprocess.Popen([self.perl, self.circos_location, params])
        pl_script.communicate()
        if pl_script.returncode!=0:
            print('Something went wrong. Image {} not generated.'.format(output_image))
            #sys.exit()
        else:
            print('{} generated.'.format(output_image))
            
    def remove_files(self, path, name, names):
        """
        Removes the unneccessary files from the hierarchal folders.
        """
        os.remove(os.path.join(path,'legends.{}.0.txt'.format(name)))
        os.remove(os.path.join(path,'legends.{}.1.txt'.format(name)))
        os.remove(os.path.join(path,'circos.conf'))
        for name in names:
            os.remove(os.path.join(path,'{}.txt'.format(name)))
        
            

taxonomy_ranks={'kingdom':'k', 'phylum':'p','order':'o','class':'c','family':'f','genus':'g', 'species':'s'} 
       
