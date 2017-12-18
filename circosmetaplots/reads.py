# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:22:23 2017

@author: Moa
"""

import pandas as pd
import math, os
from node import Node_y

class YearReads(object):
    """
    Reads and formats OTU data for the circos plot depicting the year.
    """
    def __init__(self, data_files, location_id):
        self.taxonomic_levels=['k','p','c','o','f','g','s']
        self.location_id=location_id
        print('Reading OTU data.')
        self.df=self.read_data(data_files)
        self.year_subset()
        self.taxonomy_columns()
        self.df_format()
        self.start_rank=None
        print('Done')

    def read_data(self,data_files):
        """Reads the data files(s) as a pandas dataframe and decodes them. 
        If more than one file, they're cmbined to one dataframe."""
        i=0
        df= [0 for i in range(len(data_files))]
        for data_file in data_files:
            df[i] = pd.read_csv(data_file,sep='\t',header=1,dtype={'#OTU ID': str},index_col='#OTU ID')
            i+=1
        total_df = pd.concat(df)   
        return total_df
    
    def year_subset(self):
        """
        Filters the dataframe so that it only contains one place, one year 
        (every other week), only species with at least 100
        reads in one filter.
        """
        index_list=['{}{}'.format(self.location_id,str(i)) for i in range(1,53,2)]
        index_list.append('taxonomy')
        self.df=self.df.loc[:,index_list]
        self.df=self.df.loc[self.df.max(axis=1)>100]

    def complexity_norm(self,df):
        """
        A normalization for complexity. 
        """
        for i in range(1,53,2):
            df.loc[df[self.location_id+str(i)]<100.0,self.location_id+str(i)]=0.0
            norm_sum=df.loc[df[self.location_id+str(i)]>0.0,self.location_id+str(i)].apply(math.log10).sum()
            normalize=lambda x:math.pow(10,(math.floor(math.log10(x*norm_sum))-1))
            df[self.location_id+str(i)]=df.loc[df[self.location_id+str(i)]>0.0,self.location_id+str(i)].apply(normalize)
        df.fillna(0.0,inplace=True)        
        return df
    
    def taxonomy_columns(self):
        """
        Adds columns for phylums, class, order etc to the dataframe.
        (Only named after the first letter of the order. That letter is still in front of the name, so that
        a species named Rosa blanda will be called s__Rosa blanda).
        """
        df2=pd.DataFrame(self.df.taxonomy.str.split('; ').tolist(),columns=self.taxonomic_levels, index=self.df.index)
        self.df=pd.concat([self.df,df2],axis=1,join_axes=[self.df.index])
    
    def df_format(self):
        """
        Renames cathegories so that the three different files will match.
        Also gives 'No blast hit' a name for ech level (unclassified).
        Those that are unclassified for only a few ranks are called unspecified to differentiate them.
        """
        for i in self.taxonomic_levels[1:]:
            self.df.loc[self.df[i]==i+'__',i]=i+'__unspecified'
            self.df.loc[self.df[i]==i+'__unidentified',i]=i+'__unspecified'
            self.df.loc[self.df[i].isnull(),i]=i+'__unspecified'
            self.df.loc[self.df.k=='No blast hit',i]=i+'__unclassified'
            self.df.loc[self.df.k=='Unassigned',i]=i+'__unclassified'
            self.df.loc[self.df[i]==i+'__Incertae sedis',i]=i+'__Incertae_sedis'
        self.df.loc[self.df.k=='No blast hit','k']='k__unclassified'
        self.df.loc[self.df.k=='Unassigned','k']='k__unclassified'
        self.df.loc[self.df.g=='g__unspecified','s']='s__unspecified'

    def reads_single(self,organisms):
        """
        Sums reads for a single plot with given taxa at given ranks. Doesn't use a tree, a specific function was needed.
        """
        maximums=[]
        names=[]
        filenames=[]
        change_name=False
        for pair in organisms:
            rank,taxon=pair[0],pair[1]
            filename='{}.txt'.format(taxon)
            if filename in filenames:
                filename='{}_{}'.format(rank,filename)
                change_name=True
            filenames.append(filename)
            maximum=self.sum_single(taxon,rank,filename)
            maximums.append(maximum)
            if change_name==True:
                taxon='{}_{}'.format(rank,taxon)
            names.append(taxon)
        return names, maximums
    
    def reads_single2(self,rank,taxon):
        """
        Sums reads for the rank below the given taxon and rank. Doesn't use a tree, a specific function was needed.
        """
        maximums=[]
        names=[]
        if len(rank)>1:
            rank=rank[0].lower()
        else:
             rank=rank.lower()
        new_rank=self.next_level(rank)
        indexes=self.df.loc[self.df[rank] == '{}__{}'.format(rank,taxon), new_rank].unique().tolist()
        for group in indexes:
            group=group[3::]
            filename='{}.txt'.format(group)
            if group=='unspecified' or group=='Incertae_sedis':
                maximum=self.sum_single_unknown(group, new_rank, filename, taxon, rank)
            else:
                maximum=self.sum_single(group,new_rank,filename)
            maximums.append(maximum)
            names.append(group)
        return names, maximums
                    
    def uppermost_layer(self,rank,root):  
        """
        For the hierarchal run of plots.
        Makes the first layer of nodes, that are children to the root.
        Sums reads for all OTU assigned the specific rank for the year.
        """
        indexes=self.df[rank].unique().tolist()
        self.start_rank=rank
        j=1
        labels=[]
        maximums=[]
        for group in indexes:
            group=group[3::]
            child=Node_y(root,group,rank)
            root.children.append(child)
            maximum=self.sum_by_week(child,'{}.txt'.format(group))
            labels.append(group)
            maximums.append(maximum)
            j+=1
        return labels, maximums
    
    def lower_layer(self,node,path):
        """
        For the hierarchal run of plots. Makes the lower layers of nodes.
        Sums reads for all OTU assigned the specific rank for the year.
        """
        new_rank=self.next_level(node.level)
        if node.name!='Incertae_sedis':
            indexes=self.df.loc[self.df[node.level] == '{}__{}'.format(node.level,node.name), new_rank].unique().tolist()
        else:
            indexes=self.incertae_children(node,new_rank)
        labels=[]
        maximums=[]
        for group in indexes:
            taxon=group[3::]
            child=Node_y(node,taxon,new_rank)
            node.children.append(child)
            filename=os.path.join(path,'{}.txt'.format(taxon))
            maximum=self.sum_by_week(child, filename)
            labels.append(taxon)
            maximums.append(maximum)
        return labels, maximums, new_rank
    
    def incertae_children(self,node,cur_level):
        """
        Finds lower ranks for OTUs that are incertae sedis for the current rank.
        """
        parent=node.parent
        while parent.name=='Incertae_sedis'and parent.level!=self.start_rank:
            parent=parent.parent
        extra=parent.children[0].level
        names=self.df[(self.df[node.level]=='{}__Incertae_sedis'.format(node.level)) & (self.df[parent.level]=='{}__{}'.format(parent.level,parent.name))& 
                        (self.df[extra]=='{}__Incertae_sedis'.format(extra))][cur_level].unique().tolist()
        return names
    
    def next_level(self,cur_letter):
        """
        Finds the next organizational level for the tree.
        """
        list_index=self.taxonomic_levels.index(cur_letter)+1
        next_letter=self.taxonomic_levels[list_index]
        return next_letter
    
    def sum_by_week(self,node,filename):
        """
        For the hierarchal run.
        Sums reads for one group at some rank for one year (every other week).
        """
        group=node.name
        rank=node.level
        if group=='unspecified':
            sums=self.sum_unspecified(node)
        elif group=='Incertae_sedis':
            sums=self.sum_incertae_sedis(node)        
        else:
            sums=[]
            for i in range(1,53,2):
                index_line=self.location_id+str(i)
                filter_sum=self.df.loc[self.df[rank] == rank+'__'+group, index_line].sum()
                sums.append(filter_sum)            
        max_reads=max(sums)
        self.by_week_file(sums,filename)
        return max_reads
    
    def sum_single(self,taxon, rank, filename):
        """
        For a single plot.
        Sums reads for one group at some level for one year (every other week).
        """
        sums=[]
        for i in range(1,53,2):
            index_line=self.location_id+str(i)
            filter_sum=self.df.loc[self.df[rank] == rank+'__'+taxon, index_line].sum()
            sums.append(filter_sum)            
        max_reads=max(sums)
        self.by_week_file(sums,filename)
        return max_reads
    
    def sum_single_unknown(self,taxon, rank, filename, old_taxon, old_rank):
        """
        For a single plot.
        Sums reads for one group at some level for one year (every other week).
        """
        sums=[]
        for i in range(1,53,2):
            index_line=self.location_id+str(i)
            #filter_sum=self.df.loc[self.df[rank] == rank+'__'+taxon, index_line].sum()
            filter_sum=self.df[(self.df[rank]=='{}__{}'.format(rank,taxon)) & (self.df[old_rank]=='{}__{}'.format(old_rank,old_taxon))][index_line].sum()
            sums.append(filter_sum)            
        max_reads=max(sums)
        self.by_week_file(sums,filename)
        return max_reads
    
    def sum_unspecified(self,node):
        """
        For a hierarchal run.
        Sums reads for OTUs that are unspecified for the current rank.
        """
        parent=node.parent
        while parent.name=='unspecified'and parent.level!=self.start_rank:
            parent=parent.parent
        sums=[]
        extra=parent.children[0].level
        for i in range(1,53,2):
            index_line=self.location_id+str(i)
            filter_sum=self.df[(self.df[node.level]=='{}__unspecified'.format(node.level)) & (self.df[parent.level]=='{}__{}'.format(parent.level,parent.name))& 
                        (self.df[extra]=='{}__unspecified'.format(extra))][index_line].sum()
            sums.append(filter_sum)
        return sums
    
    def sum_incertae_sedis(self,node):
        """
        For a hierarchal run.
        Sums reads for OTUs that are incertae sedis for the current rank.
        """
        parent=node.parent
        while parent.name=='Incertae_sedis'and parent.level!=self.start_rank:
            parent=parent.parent
        extra=parent.children[0].level
        sums=[]
        for i in range(1,53,2):
            index_line=self.location_id+str(i)
            #print('extra',extra)
            filter_sum=self.df[(self.df[node.level]=='{}__Incertae_sedis'.format(node.level)) & (self.df[parent.level]=='{}__{}'.format(parent.level,parent.name))& 
                        (self.df[extra]=='{}__Incertae_sedis'.format(extra))][index_line].sum()
            sums.append(filter_sum)
        return sums

    def total_by_week(self):
        """
        Returns the logarithm of ten of the sums of the total number of reads for the weeks.
        """
        sums=[]
        for i in range(1,53,2):
            index_line=self.location_id+str(i)
            week=self.df[index_line].sum(axis=0)
            sums.append(math.log10(week))
        return sums 
    
    def by_week_file(self,data,filename):
        """Makes a file for circos to use to plot data for every other week of a year."""
        f=open(filename,'w')
        j=0
        for i in range(1,52,2):
            line1='week{}\t3\t7\t{}\n'.format(i-1,str(data[j]))
            line2='week{}\t0\t7\t{}\n'.format(i,str(data[j]))
            line3='week{}\t0\t3\t{}\n'.format(i+1,str(data[j]))
            f.write(line1+line2+line3)                
            j+=1
        f.write('week52\t3\t7\t{}\n'.format(str(data[0])))
        f.close()
        
if __name__ == "__main__":
    b=YearReads([r'C:\Users\Moa\Documents\Data\16S_otu_table_sorted2.txt'], 'Kiruna-2006-')
    n, m=b.reads_single2('k','Bacteria')
    print(n,m)
