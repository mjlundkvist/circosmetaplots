# -*- coding: utf-8 -*-
"""
Created on Thu May 11 13:46:33 2017

@author: Moa
"""

from node import Node
import pandas as pd
import math 

class OTUreads(object):
    """
    Reads and formats OTU data for the taxonomic circos plot.
    """
    def __init__(self, data_files):
        self.taxonomic_levels=['k','p','c','o','f','g','s']
        print('Reading otu data')
        self.df=self.read_data(data_files)
        print('Done')
        self.taxonomy_columns()
        self.df_format()
        self.root=Node(None, 'root', None)
        self.leaves=[]
        self.k_index=0
        self.start_rank=None
        self.stop_rank=None
        self.index_list=None
        
    def read_data(self,data_files):
        """Reads the data as a pandas dataframe."""
        i=0
        df= [0 for i in range(len(data_files))]
        for data_file in data_files:
            df[i] = pd.read_csv(data_file,sep='\t',index_col='#OTU ID', header=1,dtype={'#OTU ID': str})
            i+=1
        total_df = pd.concat(df)    
        return total_df
        
    def complexity_norm(self,df):
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
        for i in self.taxonomic_levels[1:]:
            self.df.loc[self.df.k=='No blast hit',i]=i+'__unclassified'
        self.df.loc[self.df.k=='No blast hit','k']='k__unclassified'

    def df_format(self):
        for i in self.taxonomic_levels[1:]:
            self.df.loc[self.df[i]==i+'__',i]=i+'__unspecified'
            self.df.loc[self.df[i]==i+'__unidentified',i]=i+'__unspecified'
            self.df.loc[self.df[i].isnull(),i]=i+'__unspecified'
            self.df.loc[self.df.k=='No blast hit',i]=i+'__unclassified'
            self.df.loc[self.df.k=='Unassigned',i]=i+'__unclassified'
        self.df.loc[self.df.k=='Unassigned','k']='k__unclassified'
        self.df.loc[self.df.k=='No blast hit','k']='k__unclassified'
        self.df.loc[self.df.g=='g__unspecified','s']='s__unspecified' #Remember this!
    
    def year_subset(self,year_pl):
        """
        Filters the dataframe so that it only contains one place, one year 
        (every other week), only species with at least 100
        reads in one filter.
        """
        index_list=[year_pl+str(i) for i in range(1,53,2)]
        index_list.extend(self.taxonomic_levels)
        df=self.df.loc[:,index_list]
        self.df=df.loc[df.max(axis=1)>100]
        
    def make_tree(self,year_pl,start_rank,plot_level):
        """
        Builds a tree to organize the data of the phylogenetic tree 
        and writes files for circos to use.
        arguments: year_pl -- year and place in format like Kiruna-2006-
        """
        self.year_subset(year_pl)
        self.index_list=[year_pl+str(i) for i in range(1,53,2)]
        if len(start_rank)==1:
            self.start_rank=start_rank.lower()
        else:
            self.start_rank=start_rank[0].lower()
        if len(plot_level)==1:
            self.stop_rank=plot_level.lower()
        else:
            self.stop_rank=plot_level[0].lower()
        self.uppermost_layer(start_rank)
        for child in self.root.children:
            self.tree_layer(child)

        self.set_count(start_rank)
        
        for node in self.root.children:
            self.find_positions(node)
            self.write_karyotype(node)
            
        self.write_file()
        return len(self.root.children)

    
    def uppermost_layer(self,letter):
        """
        Makes the first layer of the tree, and writes out a file to 
        use as karyotype in circos.
        """       
        indexes=self.df[letter].unique().tolist()
        for i in range(len(indexes)):
            label='k'+str(i)
            cur_node=Node(self.root,indexes[i],letter)
            cur_node.type=label          
            self.root.children.append(cur_node)

    def tree_layer(self, node):
        """
        Recursive function that continues to make new nodes until the layer 's' for species is reached.
        unspecifieds not a problem since they continue to be unspecified, and the there's a function that makes 
            certain unspecifieds from other branches aren't included when species are reached
        """ 
        cur_level=self.next_rank(node)        
        if node.name[3::]!='Incertae sedis':
            children=self.df.loc[self.df[node.level] == node.name, cur_level].unique().tolist()
        else:
            children=self.incertae_children(node,cur_level)
        for name in children:
            cur_node=Node(node,name,cur_level)
            node.children.append(cur_node)
            if cur_level!=self.stop_rank:
                self.tree_layer(cur_node)
            else:
                self.leaves.append(cur_node)
                self.find_reads(cur_node)
                
    def incertae_children(self,node,cur_level):
        """
        Finds lower ranks of taxons that are classified as Incertae Sedis.
        """
        parent=node.parent
        while parent.name[3::]=='Incertae sedis' and parent.level!=self.start_rank:
            parent=parent.parent
        extra=parent.children[0].level
        indexes=self.df[(self.df[node.level]=='{}__Incertae sedis'.format(node.level)) & (self.df[parent.level]==parent.name)& 
                        (self.df[extra]=='{}__Incertae sedis'.format(extra))][cur_level].unique().tolist()
        #print(parent.name, indexes)
        return indexes
    
    def find_reads(self, node):
        """
        Finds the index codes that each OTU has and saves it in the node, 
        Note that the same species may have several different OTU IDs
        #species get a name in that file, but since they are going to be messing around with it you can't know for sure
        elif name =='s__unidentified':
            indexes=self.find_unidentified(df,node)
        """
        
        if node.name[3::]=='unspecified': 
            serie=self.find_unspecifieds(node)
        elif node.name[3::]=='Incertae sedis':
            serie=self.find_incertae(node)
        else:
            serie=self.df.loc[self.df[node.level]==node.name,self.index_list].sum()
        node.serie, node.max=self.normalization(serie)        
        
    def find_unspecifieds(self, node):
        """
        Finds indexes for unspecifieds (since several layers above species can also be unspecified).
        """
        parent=node.parent
        while parent.name[3::]=='unspecified' and parent.level!=self.start_rank:
            #'unidentified'
            parent=parent.parent
        extra=parent.children[0].level
        reads=self.df.loc[(self.df.s=='s__unspecified') & (self.df[parent.level]==parent.name) & (self.df[extra]=='{}__unspecified'.format(extra)),self.index_list].sum()
        return reads
    
    def find_incertae(self, node):
        """
        Finds indexes for unspecifieds (since several layers above species can also be unspecified).
        """
        parent=node.parent
        while parent.name[3::]=='Incertae sedis' and parent.level!=self.start_rank:
            parent=parent.parent
        extra=parent.children[0].level
        reads=self.df.loc[(self.df[node.level]=='{}__Incertae sedis'.format(node.level)) & (self.df[parent.level]==parent.name)& 
                        (self.df[extra]=='{}__Incertae sedis'.format(extra)),self.index_list].sum()
        return reads

    def normalization(self, serie):
        """Normalizes a pandas serie so that the values are between 0 and 1."""
        x_max=serie.max()
        x_min=serie.min()        
        normalize=lambda x:(x-x_min)/(x_max-x_min) 
        norm=serie.map(normalize)

        return norm, x_max
    
    def next_rank(self,node):
        """
        Finds the next and organizational level for the tree, given a node.
        """
        last_level=node.level
        list_index=self.taxonomic_levels.index(last_level)+1
        cur_level=self.taxonomic_levels[list_index]
        return cur_level

    def next_level(self,last_level):
        """
        Finds all the nodes on the level above the current one.
        """
        level=set()     
        for node in last_level:
            level.add(node.parent)
        return list(level)

   
    def write_file(self):
        """
        Writes files for the heatmaps.
        """
        j=0
        f=open('heatmap.weeks.txt','w')
        for i in range(1,52,2):
            self.write_week(j)
            j+=1
        f.close()
    
    def write_week(self,j):
        """Writes one heatmap for one week for all the nodes at the level."""
        f=open('week.{}.txt'.format(j),'w')
        for node in self.leaves:
            line='{}\t{}\t{}\t{}\n'.format(node.type, node.start, node.stop, node.serie[j])
            f.write(line)  
        f.close()
        
    def find_positions(self,node):
        """Sets positions in the nodes, given how many counts the node has."""
        start=node.start
        for child in node.children:
            child.start=start
            start+=child.counts
            child.stop=start
            if child.level!=self.stop_rank:
                self.find_positions(child)
                
    def set_count(self,start_rank):
        """
        Sets the width for the lowest level to be plotted to be 1.
        Sets the counts for the second lowest level to be its numbers of children.
        """
        for node in self.leaves:
            node.counts=1
        self.set_higher_counts(self.leaves, start_rank)               
            
    def set_higher_counts(self,level, start_rank):
        """
        Sets the width of the higher levels to be the sum of its 
        childrens widths. Recursive.
        """
        new=self.next_level(level)
        for node in new:
            counts=0
            for child in node.children:
                counts+=child.counts
            node.counts=counts
        if new[0].level!=start_rank:
            self.set_higher_counts(new, start_rank)
        
    def write_karyotype(self, node):
        """
        Writes the karyotype file which is the basis for the ideograms for the tree file. 
        Also writes a highlight file which covers all the ideograms.
        """
        end=node.counts#-1
        color='chr'+str(self.k_index)
        label='k'+str(self.k_index)
        name=node.name[3::]
        line='chr\t-\t{}\t{}\t0\t{}\t{}\n'.format(label,name,str(end),color)
        with open(r'karyotype.otu.txt', "a") as myfile:
            myfile.write(line)
        self.season_highlight(label,end)
        self.k_index+=1
        for child in node.children:
            self.write_levels(child)
            
    def season_highlight(self,label,end):
        """
        Writes the size of one ideogram into the season highlight file.
        """
        line='{}\t0\t{}\n'.format(label,str(end))
        with open(r'heat.highlight.txt', "a") as myfile:
            myfile.write(line)

    def write_levels(self, node):
        """
        Recursive function that writes highlight and text files for levels up to and including the stop level.
        """
        cur_level=node.level
        filename=r'text.{}.txt'.format(cur_level)
        filename2=r'highlight.{}.txt'.format(cur_level)
        label=node.type
        end=str(node.stop)
        name=node.name[3::]
        if name=='' or name=='unspecified' or name=='unidentified':
            name='unclassified'
        start=node.start
        line_hl='{}\t{}\t{}\n'.format(label, start, end)
        line_label='{}\t{}\t{}\t{}\n'.format(label, start, end, name)
        with open(filename, "a") as myfile:
            myfile.write(line_label)
        with open(filename2, "a") as myfile:
            myfile.write(line_hl)
        if cur_level==self.stop_rank:
            line_max='{}\t{}\t{}\t{}\n'.format(label, start, end, math.log10(node.max))
            with open('max.txt', "a") as myfile:
                myfile.write(line_max)
        else:
            for child in node.children:
                self.write_levels(child)
    
    def rank_legend(self,rank):
        ind=self.taxonomic_levels.index(rank)
        if ind>0:
            legend_rank=self.taxonomic_levels[ind-1]
            return self.df[legend_rank].value_counts().index[0][3::]
        else:
            return None
                
if __name__ == "__main__":
    #data_files=[r'C:\Users\Moa\Documents\Data\OTU_tables\16S_otu_table_mc2_w_tax_sorted.txt',	r'C:\Users\Moa\Documents\Data\OTU_tables\ITS_otu_table_mc2_w_tax_sorted.txt', r'C:\Users\Moa\Documents\Data\OTU_tables\rbcL_otu_table_mc2_w_tax_sorted.txt']
    data_files=[r'C:\Users\Moa\Documents\Data\OTU_tables\ITS_otu_table_mc2_w_tax_sorted.txt']
    ID=r'..\Data\IDs.txt'
    year='Ljungbyhed-2006-'
    a=OTUreads(data_files,ID)
    a.make_uneven_tree(year,'p','f')    
    #a.clear_files('c','g')
    #a.make_tree(year,'p','f')
    #




