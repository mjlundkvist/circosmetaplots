# -*- coding: utf-8 -*-
"""
Created on Wed May 24 11:14:12 2017

@author: Moa
"""
import pandas as pd
import numpy as np

class WeatherData(object):
    """
    Deals with the weather data for the circos plots.
    """
    def __init__(self,file,year, wind):
        self.locations=self.read_places(file)        
        self.df=self.read_data(file)
        self.last_year,self.next_year=self.set_year(year)
        if wind!=False:
            print('Reading in wind direction data.')
            self.wind_locations=self.read_wind_places(wind)
            self.wind_df=self.read_wind(wind)
            self.last_wind_df,self.next_wind_df=self.set_wind_year(year)
            print('Done')
        
    def read_places(self,file):
        """
        Reads the weather data excel file one preliminary time to see what places are in it.
        The information is used because all the weather parameters are named the same for the
        different locaitons.
        """
        df = pd.read_csv(file,sep='\t',encoding = "ISO-8859-1",nrows=1)
        places=[]
        for column in df.columns:
            if column[0:7]!='Unnamed':
                column=column.split('.')[0].lower()
                if column not in places:
                    places.append(column)
                
        return places
    
    def read_wind_places(self,file):
        """
        Reads the weather data excel file one preliminary time to see what places are in it.
        The information is used because all the weather parameters are named the same for the
        different locaitons.
        """
        df = pd.read_csv(file,sep='\t',encoding = "ISO-8859-1",nrows=1)
        places=[]
        for column in df.columns:
            if column[0:7]!='Unnamed':
                places.append(column.lower())
        return places
            
    def read_data(self,file):
        """
        Reads the data and reurns it as a pandas dataframe.
        Also renames the unnamed columns to something more descriptive.
        """
        print('Reading in general weather data.')
        df = pd.read_csv(file, sep='\t',header=1, converters={'YYYY':int,'MM':int,'DD':int,'WW':int},encoding = "ISO-8859-1")#index_col=[4] to use the dates as index
        #df=df.rename(columns={'Unnamed: 3':'Date', 'Unnamed: 4':'WW'})
        print('Done.\n')
        return df

    def read_wind(self,file):
        
        df2 = pd.read_csv(file,header=1,sep='\t',encoding = "ISO-8859-1",converters={'YYYY':int,'MM':int,'DD':int,'WW':int}) #för sidan med vindriktningen     
        #df2=df2.rename(columns={'Unnamed: 3':'Date', 'Unnamed: 4':'WW'})
        
        return df2    
    
    def set_year(self,year):
        """Discards information from other years than the one requested."""
        year=int(year)
        last_year=self.df.loc[self.df['YYYY']==year-1]
        next_year=self.df.loc[self.df['YYYY']==year+1]
        self.df=self.df.loc[self.df['YYYY']==year]
        return last_year, next_year
    
    def set_wind_year(self, year):
        year=int(year)
        next_wind_df=self.wind_df.loc[self.wind_df['YYYY']==year+1]
        last_wind_df=self.wind_df.loc[self.wind_df['YYYY']==year-1]
        self.wind_df=self.wind_df.loc[self.wind_df['YYYY']==year]
        return last_wind_df, next_wind_df
    
    def karyotype_week(self):
        """OTU reads writes the other kind of karyotype file. Considering reads are used to color, use the one in conf instead."""
        f=open('karyotype.weeks.txt','w')
        for i in range(1,53):
            f.write('chr\t-\tweek{0}\t{0}\t0\t7\tlgrey\n'.format(i))
        f.close()
        
    def find_data(self, filename, place, parameter, resolution, statistic=None):
        #parameter=self.add_units(parameter)
        parameter=self.find_parameter(parameter,place)
        if resolution=='DD':
            days=self.day_res(parameter)
            self.days_file(days, filename)
        elif resolution=='WW':
            weeks=self.weeks(statistic.lower(),parameter)
            self.weeks_file(weeks, filename)
        elif resolution=='MM':
            months=self.month(statistic.lower(),parameter)
            self.months_file(months, filename)
        else:
            print('Resolution wrong format')
            
    def weeks(self,statistic,parameter):
        if statistic=='sum':
            weeks=self.sum_week(parameter)
        elif statistic=='max':
            weeks=self.max_week(parameter)
        elif statistic=='min':
            weeks=self.min_week(parameter)
        elif statistic=='avg' or statistic=='mean':
            weeks=self.mean_week(parameter)
        else:
            raise ValueError('Weather data did not get an accepted statistic to use for weeks')
        return weeks
           
    def month(self,statistic,parameter):
        if statistic=='sum':
            months=self.sum_month(parameter)
        elif statistic=='max':
            months=self.max_month(parameter)
        elif statistic=='min':
            months=self.min_month(parameter)
        elif statistic=='avg' or statistic=='mean':
            months=self.mean_month(parameter)
        else:
            raise ValueError('Weather data did not get an accepted statistic to use for weeks')
        return months

    def find_parameter_wind(self, param, place):
        """
        Gives the number corresponding to the location to identify the correct column.
        """
        place=place.lower()
        if place in self.wind_locations:
            lst_ind=self.wind_locations.index(place)
        else:
            print('Whoops, location not in data.')
        if lst_ind==0:
            return param
        else:
            param+='.'+str(lst_ind)
        return param
    
    def find_parameter(self,param,place):
        """
        Gives the number corresponding to the location to identify the correct column.
        """
        place=place.lower()
        if place in self.locations:
            lst_ind=self.locations.index(place)
        else:
            print('Whoops, location not in data.')
        if lst_ind==0:
            return param
        else:
            param+='.'+str(lst_ind)
        return param
    
    def days_file(self,data, filename):
        """
        Writes the text file for the weather data at daily resolution.
        Range starts with 1 because the first position is empty.
        """
        f=open(filename,'w')
        for i in range(1,53):
            week=data[i]
            week_l=len(week)
            if week_l!=7:
                print('why is the length of the week wrong?:',week_l)
                """This should not be needed, if it prints something is wrong."""
                if i==1:
                    week=week[0:7]
                elif i==52:
                    week=week[week_l-7:week_l]
            for j in range(len(week)):
                line='week{}\t{}\t{}\t{}\n'.format(i,j,j+1,week[j])
                f.write(line)
        f.close()                
    
    def weeks_file(self, data, filename):
        """
        Writes the text file with the weekly resolution.
        Range starts with 1 because the first position is empty.
        """
        f=open(filename,'w')
        for i in range(1,53):#range(1,53) for all weeks
            line='week{}\t0\t7\t{}\n'.format(i,data[i])
            f.write(line)
        f.close()
    
    def months_file(self,data, filename):
        """
        Finds the weeks over which all the months span and sends that information
        to the function that writes in the file. 
        Range starts with 1 because the first position is empty. 
        """
        f=open(filename,'w')
        for i in range(1,13):
            val=data[i]
            days=self.df.loc[self.df['MM']==i, 'WW'].value_counts().tolist()
            weeks=self.df.loc[self.df['MM']==i, 'WW'].value_counts().index.tolist()
            self.weeks_in_months(f,days,weeks,val)
        f.close()
            
    def weeks_in_months(self, f,days,weeks, val):
        """
        Writes in the text file for use in a circos plot with the data at the monthly resolution.
        """
        for j in range(len(weeks)):
            cur_week=weeks[j]
            if days[j]==7:
                start=0
                stop=7
            else:
                if cur_week==min(weeks):
                    start=7-days[j]
                    stop=7
                else:
                    start=0
                    stop=days[j]
            line='week{}\t{}\t{}\t{}\n'.format(weeks[j],start,stop,val)
            f.write(line)
            
            
    """Functions that calculate the monthly data."""           
    def sum_month(self,parameter):
        months=[0,]
        for i in range(1,13):
            val=self.df.loc[self.df['MM']==i,parameter].sum()
            months.append(val)
        return months
            
    def min_month(self,parameter):
        months=[0,]
        for i in range(1,13):
            val=self.df.loc[self.df['MM']==i,parameter].min()
            months.append(val)
        return months
            
    def max_month(self,parameter):
        months=[0,]
        for i in range(1,13):
            val=self.df.loc[self.df['MM']==i,parameter].max()
            months.append(val)
        return months
            
    def mean_month(self,parameter):
        months=[0,]
        for i in range(1,13):
            val=self.df.loc[self.df['MM']==i,parameter].mean()
            months.append(val)
        return months
    
    
    """Functions that calculate the weekly data."""
    def sum_week(self,parameter):
        weeks=[0,]
        val=self.df.loc[(self.df['WW']==1)&(self.df['MM']==1)][parameter].sum()
        val2=self.last_year.loc[(self.last_year['WW']==1)&(self.last_year['MM']==12)][parameter].sum()
        weeks.append(val+val2)
        for i in range(2,52):
            val=self.df.loc[self.df['WW']==i,parameter].sum()
            weeks.append(val)
        val=self.df.loc[(self.df['WW']==52)&(self.df['MM']==12)][parameter].sum()
        val2=self.next_year.loc[(self.next_year['WW']==52)&(self.next_year['MM']==1)][parameter].sum()
        weeks.append(val+val2)
        return weeks
            
    def min_week(self,parameter):
        weeks=[0,]
        first=self.df.loc[(self.df['WW']==1)&(self.df['MM']==1)][parameter]
        val=first.min()
        if len(first)<7:
            val2=self.last_year.loc[(self.last_year['WW']==1)&(self.last_year['MM']==12)][parameter].min()
            weeks.append(min(val,val2))
        else:
            weeks.append(val)
        for i in range(2,52):
            val=self.df.loc[self.df['WW']==i,parameter].min()
            weeks.append(val)
        last=self.df.loc[(self.df['WW']==52)&(self.df['MM']==12)][parameter]
        val=last.min()
        if len(last)<7:
            val2=self.next_year.loc[(self.next_year['WW']==52)&(self.next_year['MM']==1)][parameter].min()
            weeks.append(min(val,val2))
        else:
            weeks.append(val)
        return weeks
            
    def max_week(self,parameter):
        weeks=[0,]
        first=self.df.loc[(self.df['WW']==1)&(self.df['MM']==1)][parameter]
        val=first.max()
        if len(first)<7:
            val2=self.last_year.loc[(self.last_year['WW']==1)&(self.last_year['MM']==12)][parameter].max()
            weeks.append(max(val,val2))
        else:
            weeks.append(val)
        for i in range(2,52):
            val=self.df.loc[self.df['WW']==i,parameter].max()
            weeks.append(val)
        last=self.df.loc[(self.df['WW']==52)&(self.df['MM']==12)][parameter]
        val=last.max()
        if len(last)<7:
            val2=self.next_year.loc[(self.next_year['WW']==52)&(self.next_year['MM']==1)][parameter].max()
            weeks.append(max(val,val2))
        else:
            weeks.append(val)
        return weeks
            
    def mean_week(self,parameter):
        weeks=[0,]
        first=self.df.loc[(self.df['WW']==1)&(self.df['MM']==1)][parameter]
        if len(first)<7:
            val2=self.last_year.loc[(self.last_year['WW']==1)&(self.last_year['MM']==12)][parameter].tolist()
            first=first.tolist()
            first.extend(val2)
            weeks.append(np.mean(first))
        else:
            weeks.append(first.mean())
        for i in range(2,52):
            val=self.df.loc[self.df['WW']==i,parameter].mean()
            weeks.append(val)
        last=self.df.loc[(self.df['WW']==52)&(self.df['MM']==12)][parameter]
        if len(last)<7:
            val2=self.next_year.loc[(self.next_year['WW']==52)&(self.next_year['MM']==1)][parameter].tolist()
            last=last.tolist()
            last.extend(val2)
            weeks.append(np.mean(last))
        else:
            weeks.append(last.mean())
        return weeks
            
    def day_res(self,parameter):
        """
        Finds data for each day of each week over the year.
        """
        weeks=[[],]
        days=self.days_first_week(parameter)
        weeks.append(days)
        for i in range(2,52):
            days=self.df.loc[self.df['WW']==i,parameter].tolist()
            #print(days)
            weeks.append(days)
        days=self.days_last_week(parameter)
        weeks.append(days)
        return weeks
    
    def days_first_week(self,parameter):
        """
        The first week gets its own function since it may overlap with the last year.
        """
        days=self.df.loc[(self.df['WW']==1)&(self.df['MM']==1)][parameter].tolist()
        if len(days)==7:
            return days
        elif len(days)<7:
            days2=self.last_year.loc[(self.last_year['WW']==1) & (self.last_year['MM']==12)][parameter].tolist()
            return days2+days
        else:
            #Should never happen.
            print('week 1 in trouble')
    
    def days_last_week(self,parameter):
        """
        The last week gets its own function since it may overlap with the next year.
        """
        week=[]
        days=self.df.loc[(self.df['WW']==52)&(self.df['MM']==12)][parameter].tolist()
        if len(days)==7:
            return days
        elif len(days)<7:
            days2=self.next_year.loc[(self.next_year['WW']==52) & (self.next_year['MM']==1)][parameter].tolist()
            return days+days2
        else:
            print('week 52 in trouble')
        return week
    
    def seasons(self, location):
        """
        Finds the seasons over the year.
        """
        winter_start, winter_end,summer,dates=self.season_dates(location)
        colors=self.season_colors(winter_start,winter_end, summer,len(dates))
        seasons=self.seasons_file(colors,dates)
        return seasons
    
    def season_dates(self,location):
        """
        Finds the start days for each season.
        """
        temp=self.find_parameter('Temperature',location)
        dates=[]
        winter_start=self.last_winter(temp)
        index_lst=self.df.index.tolist()
        try:
            start=self.df[0:7]['WW'].value_counts()[52]
        except KeyError:
            start=0
        dates.append(index_lst[start])        
        if winter_start==False:
           w1_ind=self.late_winter(temp,index_lst)
           dates.append(w1_ind)
        spring_ind = self.spring(temp,index_lst)
        if spring_ind:
            dates.append(spring_ind)
        sum_ind=self.summer(temp,index_lst,spring_ind)
        if sum_ind:
            dates.append(sum_ind)
            summer=True
        else:
            summer=False
            sum_ind=spring_ind
        aut_ind=self.autumn(temp,index_lst,sum_ind)
        if aut_ind:
            dates.append(aut_ind)
            w_ind=self.winter(temp,index_lst,aut_ind)
            if w_ind:
                dates.append(w_ind)
                winter_end=True
            else:
                winter_end=False
        else:
            winter_end=False  
        try:
            end=index_lst[-1]-self.df[-7:]['WW'].value_counts()[1]+1
        except KeyError:
            end=index_lst[-1]+1
        dates.append(end)
        return winter_start, winter_end,summer,dates
        

    def spring(self,temp,index_lst):
        """
        Finds the index of the day spring starts. February 15 is the first possible day according to SMHI, so it starts from there.
        """
        start=31+14#because spring can't start before feb 15
        last=self.df.loc[self.df['MM']==7,'Date'].index.tolist()[-1]#latest: 31 juli
        stop=index_lst.index(last)
        for i in index_lst[start:stop+1]:
            if self.df.loc[i:i+6][temp].min()>0.0 and self.df.loc[i:i+6][temp].max()<10.0:# and <10.0? #the seventh will be included in pandas
                return i
        return None, start
            
    def summer(self,temp,index_lst,start):
        """
        Finds the index of the day summer starts.
        """
        index_start=index_lst.index(start)
        for i in index_lst[index_start::]:
            if self.df.loc[i:i+4][temp].min()>10.0:
                return i
        return None
            
    def autumn(self,temp,index_lst,summer_start):
        """
        Finds the index of the day autumn starts. August 1 is the first possible day.
        """
        start=index_lst.index(summer_start)
        first=self.df.loc[self.df['MM']==8,'Date'].index.tolist()[0]#earliest: 1aug
        start2=index_lst.index(first)
        if start<start2:
            start=start2
        #latest:15 feb
        for i in index_lst[start::]:
            if self.df.loc[i:i+4][temp].max()<10.0:
                return i
        raise ValueError('Endless summer... in Sweden?!')
        return None
        
    def winter(self,temp,index_lst,autumn_start):
        """
        Finds the index of the day winter starts.
        """
        start=index_lst.index(autumn_start)      
        for i in index_lst[start::]:
            if self.df.loc[i:i+4][temp].max()<0.0:
                return i       
        #if winter never starts... sadly possible
        return None
    
    def late_winter(self,temp,index_lst):
        """
        Returns the index for the first day of winter in the beginning of the year. Is used in case the last year didn't have a winter.
        """
        end=31+14
        for i in index_lst[0:end]:
            if self.df.loc[i:i+4][temp].max()<0.0:
                return i
        return None, None  
    
    def last_winter(self,temp):
        """
        Checks if winter started during the last year or not.
        """
        index_lst=self.last_year.index.tolist()
        first=self.last_year.loc[self.last_year['MM']==8,'Date'].index.tolist()[0]
        #earliest: 1aug
        start=index_lst.index(first)
        for i in index_lst[start::]:
            if self.last_year.loc[i:i+4][temp].max()<0.0:
                return True
        return False
    
    def season_colors(self,last_winter,end_winter,summer, length):
        if last_winter==True:
            if end_winter==True:
                if summer==True:
                    ids=[1,2,3,4,1]
                else:
                    ids=[1,2,4,1]
            else:
                if summer==True:
                    ids=[1,2,3,4]
                else:
                    ids=[1,2,4]
        elif last_winter==False:
            if end_winter==True:
                if summer==True:
                    if length==7:
                        ids=[4,1,2,3,4,1]
                    elif length==6:
                        ids=[4,2,3,4,1]
                else:
                    if length==6:
                        ids=[4,1,2,4,1]
                    elif length==5:
                        ids=[4,2,4,1]
            else:
                if summer==True:
                    if length==6:
                        ids=[4,1,2,3,4]
                    elif length==5:
                        ids=[4,2,3,4,1]
                else:
                    if length==5:
                        ids=[4,1,2,4]
                    elif length==4:
                        ids=[4,2,4]
        else:
            print('Something is not quite right with the seasons...')
        return ids
    
    def seasons_file(self,colors,dates):
        """
        Writes a file for the season highlight plot.
        """
        f=open('highlight.seasons.txt','w')
        j=0
        seasons=[]
        for i in range(0, len(dates)-1): 
            pair=(dates[i],dates[i+1]-1)
            days,weeks=self.days_weeks(pair)
            season=self.phylo_seasons(days,weeks,colors[j])
            seasons.append(season)
            self.weeks_in_seasons(f,days,weeks,colors[j])
            j+=1
        f.close()
        return seasons
    
    def phylo_seasons(self,days,weeks,color):
        """
        Calculates how many weeks the seasons cover in percentages, for use in the taxonomy plot.
        """
        colors={1:'pastel1-9-qual-2',2:'pastel1-9-qual-6',3:'pastel1-9-qual-3',4:'pastel1-9-qual-5'}
        mini=min(weeks)
        maxi=max(weeks)
        a=weeks.index(mini)
        b=weeks.index(maxi)
        first_part=days[a]/7.0
        last_part=days[b]/7.0
        return [first_part,mini,last_part, maxi, colors[color]]
        
    def days_weeks(self,pair):
        """
        Finds the weeks and days for the indexes corresponding to the start and end of a season.
        """
        start,stop=pair[0],pair[1]
        days=self.df.loc[start:stop]['WW'].value_counts().tolist()
        weeks=self.df.loc[start:stop]['WW'].value_counts().index.tolist()
        return days,weeks

    def weeks_in_seasons(self, f, days, weeks, s_id):
        """
        Writes a season into the highlight file
        """
        for j in range(len(weeks)):
            cur_week=weeks[j]
            if days[j]==7:
                start=0
                stop=7
            else:
                if cur_week==min(weeks):
                    start=7-days[j]
                    stop=7
                else:
                    start=0
                    stop=days[j]
            line='week{}\t{}\t{}\tid={}\n'.format(weeks[j],start,stop,s_id)
            f.write(line)
            
    def wind_track(self,place):
        """
        Collects the data needed for the wind plot and writes the neccessary data files.
        """
        id_weeks=self.wind_direction(place)
        data_weeks=self.wind(place)
        self.wind_track_file(data_weeks, id_weeks)
        
    def wind_direction(self,place):
        """
        Takes mean wind direction for all days for all weeks during the year.
        """
        parameter=self.find_parameter_wind('Wind Direction',place)
        weeks=[[] for i in range(53)]
        weeks[1]=self.first_wind(parameter)
        for week in range(2,52):
            days=self.wind_df.loc[self.wind_df['WW']==week,'DD'].unique().tolist()
            for day in days:
                wind_dir=self.wind_df.loc[(self.wind_df['WW']==week) & (self.wind_df['DD']==day)][parameter].mean()
                weeks[week].append(self.direction(wind_dir))
        weeks[52]=self.last_wind(parameter)
        return weeks
    
    def first_wind(self,parameter):
        """
        Finds the wind direction for the first week of the year.(Own function since the week may overlap the last year.)
        """
        week=[]
        days=self.wind_df.loc[(self.wind_df['WW']==1) & (self.wind_df['MM']==1)]['DD'].unique().tolist()
        if len(days)==7:
            for day in days:
                wind_dir=self.wind_df.loc[(self.wind_df['WW']==1) & (self.wind_df['DD']==day)][parameter].mean()
                week.append(self.direction(wind_dir))
        elif len(days)<7:
            days2=self.last_wind_df.loc[(self.last_wind_df['WW']==1) & (self.last_wind_df['MM']==12)]['DD'].unique().tolist()
            for day in days2:
                wind_dir=self.last_wind_df.loc[(self.last_wind_df['WW']==1) & (self.last_wind_df['DD']==day)][parameter].mean()
                week.append(self.direction(wind_dir))
            for day in days:
                wind_dir=self.wind_df.loc[(self.wind_df['WW']==1) & (self.wind_df['DD']==day)][parameter].mean()
                week.append(self.direction(wind_dir))
        else:
            print('week 1 in trouble')
        return week
    
    def last_wind(self,parameter):
        """
        Finds the wind direction for the last week of the year.(Own function since hte week may overlap the next year.)
        """
        week=[]
        days=self.wind_df.loc[(self.wind_df['WW']==52) & (self.wind_df['MM']==12)]['DD'].unique().tolist()
        if len(days)==7:
            for day in days:
                wind_dir=self.wind_df.loc[(self.wind_df['WW']==52) & (self.wind_df['DD']==day)][parameter].mean()
                week.append(self.direction(wind_dir))
        elif len(days)<7:
            days2=self.next_wind_df.loc[(self.next_wind_df['WW']==52) & (self.next_wind_df['MM']==1)]['DD'].unique().tolist()
            for day in days:
                wind_dir=self.wind_df.loc[(self.wind_df['WW']==52) & (self.wind_df['DD']==day)][parameter].mean()
                week.append(self.direction(wind_dir))
            for day in days2:
                #print(day)
                wind_dir=self.next_wind_df.loc[(self.next_wind_df['WW']==52) & (self.next_wind_df['DD']==day)][parameter].mean()
                week.append(self.direction(wind_dir))
        else:
            print('week 52 in trouble')
        return week
            
    def wind_track_file(self,data,id_weeks):
        """
        Writes the text file for the wind data at daily resolution. IDs corresponding to the wind direction is included.
        Range starts with 1 because the first position is empty.
        """
        f=open('wind.txt','w')
        for i in range(1,53):
            id_days=id_weeks[i]
            week=data[i]
            week_l=len(week)
            if week_l!=7:
                if i==1:
                    week=week[0:7]
                elif i==52:
                    week=week[week_l-7:week_l]
            for j in range(len(week)):
                line='week{}\t{}\t{}\t{}\tid={}\n'.format(i,j,j+1,week[j],id_days[j])
                f.write(line)
        f.close()
        
    def wind(self,place):
        """
        Finds data for each day of each week over the year.
        """
        weeks=[[],]
        #param='Medelvind (m/s)'
        param = 'Mean Wind'
        param=self.find_parameter(param, place)
        days=self.df.loc[(self.df['WW']==1)&(self.df['MM']==1)][param].tolist()
        weeks.append(days)
        for i in range(2,52):
            days=self.df.loc[self.df['WW']==i,param].tolist()
            weeks.append(days)
        days=self.df.loc[(self.df['WW']==52)&(self.df['MM']==12)][param].tolist()
        weeks.append(days)
        return weeks                 
        
    def direction(self,value):
        """
        Gives the direction corresponding to the degree of the wind direction. 
        From: https://www.campbellsci.com/blog/convert-wind-directions
        Jacob Davis from Campbell Scientific
        """

        directions=['north','northeast','east','southeast','south','southwest','west','northwest','north']
        return directions[round(value/45)]
            
if __name__ == "__main__":
    pass

