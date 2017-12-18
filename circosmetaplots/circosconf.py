# -*- coding: utf-8 -*-
"""
Created on Mon May 29 08:55:34 2017

@author: Moa
"""
import os
import math

class CircosConf(object):
    """
    Is used to write the configuration files to two kinds of circos plots: the year plot, and the taxonomic plot.
    """
    def __init__(self):
        self.colors=['lblue', 'lyellow','lpurple','lgreen', 'lorange','lred']
        self.ranks=['k','p','c','o','f','g','s']
        self.file=None        
        
    def taxo_conf(self,no_of_ideograms, location, start, stop, seasons):
        """
        Writes the entire circos configuration file for the taxonomy plot.
        """
        self.file=open('circos.conf','w')
        self.karyotype('karyotype.otu.txt')
        self.taxo_colors(no_of_ideograms)
        self.taxo_rules(no_of_ideograms)
        edge=self.overall_taxo_plots(start,stop, location, seasons)
        label_radius=edge+0.3
        ideogram_rule='<rules>\n<rule>\ncondition= var(chr)=~ /[13579]/\nlabel_radius={0:.2f}r+20p\n</rule>\n</rules>\n'.format(label_radius)
        self.ideogram(0.005, 0.13, 50, 25, '{}r'.format(label_radius), ideogram_rule)
        self.file_end(True)
        self.file.close()
        
    def uneven_conf(self,no_of_ideograms, legends, start, stop, seasons):
        """
        Writes the entire circos configuration file for a variation on taxonomy plot.
        """
        self.file=open('circos.conf','w')
        self.karyotype('karyotype.otu.txt')
        self.taxo_colors(no_of_ideograms)
        self.taxo_rules(no_of_ideograms)
        edge=self.uneven_plots(start,stop, legends, seasons)
        label_radius='dims(ideogram,radius_inner)+20p'
        label_radius2=edge+0.05
        ideogram_rule='<rules>\n<rule>\ncondition= var(label) eq "unclassified"\nlabel_radius={}r\n</rule>\n</rules>\n'.format(label_radius2)
        self.ideogram(0.005, 0.13, 75, 25, label_radius, ideogram_rule)
        self.file_end(True)
        self.file.close()
        
    def uneven_plots(self,start,stop, legends, seasons):
        """
        Writes the plot blocks and legend files for a variation of the taxonomy plot.
        """
        self.file.write('<plots>\n')
        edge=self.taxo_plots2(start,stop)
        self.auto_heatmaps(edge,0.15)
        self.taxo_seasons(seasons)
        #'{}r'.format(edge+0.8),'10.0r'
        self.text_block('legends.0.txt','0.1r','3.0r','bold', 30,'-75p','0p','yes','')
        self.uneven_legends_file(legends)
        self.file.write('</plots>\n\n')
        return edge
        
    def auto_heatmaps(self,start,gap):
        """
        Writes the plot blocks for the auto heatmap tracks.
        """
        track_step=(7.55-(start+gap))/26.0
        auto_tracks='track_width   = {0}\ntrack_start   = {1}\ntrack_step    = {0}\n\n'.format(round(track_step,2),round(start+gap,2))
        self.file.write(auto_tracks)
        self.auto_heatmap_block()    
        self.auto_heatmap_file()
        self.auto_include()
        
    def overall_taxo_plots(self,start,stop,place,seasons):
        """
        Writes all the subplot blocks and legend files for the taxonomy plot.
        """
        self.file.write('<plots>\n')
        if seasons!=False:             
            self.taxo_seasons(seasons)
        edge=self.taxo_plots(start,stop)
        self.auto_heatmaps(edge,0.6)
        self.text_block('legends.0.txt','0.15r','1.5r','bold', 22,'-75p','0p','yes','')
        self.taxo_legends_file(place)
        self.file.write('</plots>\n\n')
        return edge
        
    def taxo_plots(self,start,stop):
        """
        Writes the hierarchal donut plots in the configuration file. For the taxonomy plot.
        """
        r0=1
        step=1.3
        start_i=self.ranks.index(start)+1
        stop_i=self.ranks.index(stop)
        size_step=1
        padding=-5
        rpadding=10
        if stop_i-start_i>2:
            size=10
            step=1.05
            size_step=0
            padding=-3
            rpadding=5
        elif stop_i-start_i>1:
            size=13
            step=1.1
        else:
            size=19
        for level in self.ranks[start_i:stop_i+1]:
            if level==stop and size==10:
                step=0.8
            r1=r0+step
            file_text='text.{}.txt'.format(level)
            file_hl='highlight.{}.txt'.format(level)
            self.text_block(file_text,'{}r'.format(round(r0,2)), '{}r+100p'.format(round(r1,2)), 'light', size,'{}p'.format(rpadding),'{}p'.format(padding), 'no','')
            self.plot_block('highlight', file_hl, '{}r'.format(round(r0,2)), '{}r'.format(round(r1,2)), 'dgrey', 0, '<<include rules.txt>>\n')
            size+=size_step
            r0=r1
        r1=r0+0.5    
        self.plot_block('histogram', 'max.txt', '{0:.2f}r'.format(r0), '{0:.2f}r'.format(r1), 'black', 1, '')
        return r0
    
    def taxo_plots2(self,start,stop):
        """
        Writes the hierarchal donut plots in the configuration file. For a variation on the taxonomy plot.
        """
        r0=1
        step=1.2
        start_i=self.ranks.index(start)+1
        stop_i=self.ranks.index(stop)
        padding=-5
        rpadding=10
        size=22
        for level in self.ranks[start_i:stop_i+1]:
            if level=='f':
                step=1.4
            elif level=='p':
                step=1.1
            else:
                step=1.2
            r1=r0+step
            file_text='text.{}.txt'.format(level)
            file_hl='highlight.{}.txt'.format(level)
            self.text_block(file_text,'{0:.2f}r'.format(r0,2), '{0:.2f}r'.format(r1,2),'light', size,'{}p'.format(rpadding),'{}p'.format(padding),'no', '')
            self.plot_block('highlight', file_hl, '{0:.2f}r'.format(r0,2), '{0:.2f}r'.format(r1,2), 'dgrey', 0, '<<include rules.txt>>\n')
            r0=r1
        r1=r0+0.5    
        return r0
    
    def auto_heatmap_file(self):
        """
        Writes a heatmap text file for the taxonomy plot.
        """
        f=open('heatmap.conf','w')

        plot='<plot>\ntype=heatmap\ncolor=black_a7,black_a6,black_a5, black_a4,black_a3,black_a2,black_a1,black\n'
        auto_r0='eval(sprintf("%fr",conf(plots,track_start)+counter(heatcounter)*conf(plots,track_step)))\n'
        auto_r1='eval(sprintf("%fr",conf(plots,track_start)+conf(plots,track_width)+counter(heatcounter)*conf(plots,track_step)))\n'
        heatmap_block='file=week.counter(heatcounter).txt\nr1={}\nr0={}\n'.format(auto_r1, auto_r0)
        minmax='min=0\nmax=1\n'
        counters='post_increment_counter = heatcounter:1\n'
        end='</plot>\n\n'

        f.write(plot+heatmap_block+minmax+counters+end)
        f.close()
        
    def auto_heatmap_block(self):
        """
        Writes the first auto heatmap into the main circos configuration file. Exists to initialize the counters.
        """        
        auto_r0='eval(sprintf("%fr",conf(plots,track_start)+counter(heatcounter)*conf(plots,track_step)))'
        auto_r1='eval(sprintf("%fr",conf(plots,track_start)+conf(plots,track_width)+counter(heatcounter)*conf(plots,track_step)))'
        color='black_a7,black_a6,black_a5, black_a4,black_a3,black_a2,black_a1,black'
        options='init_counter = heatcounter:0\npost_increment_counter = heatcounter:1\nmin=0\nmax=1\n'
        self.plot_block('heatmap', 'week.counter(heatcounter).txt', auto_r0, auto_r1, color, 1, options)
    
    def auto_include(self):
        """
        25 lines to include the heatmaps around the taxonomy donut. One block is kept inside the conf-file, to initiate the counter used.
        """
        for i in range(25):
            self.file.write('<<include heatmap.conf>>\n')
        self.file.write('\n')
            
    def taxo_rules(self,j):
        """
        Makes a rule file that colors the plot after the ideogram it is to. 
        """
        f=open('rules.txt','w')
        rules='<rules>\n'
        for i in range(j):
            rule_hl='<rule>\ncondition  = on(k{})\nfill_color= {}\nflow = continue\n</rule>\n'.format(i,'chr{}'.format(i))
            rules+=rule_hl
        rules+='</rules>\n\n'
        f.write(rules)
        f.close()
        
    def taxo_seasons(self,seasons):
        for season in seasons:
            perc1,first,perc2,last,color=season[0],season[1],season[2],season[3],season[4]
            if first%2:
                if perc1==1.0:
                    start=math.floor(first/2.0)
                else:
                    start=math.ceil(first/2.0)-perc1
                if last%2:
                    if perc2==1.0:
                        stop=math.ceil(last/2.0)
                    else:
                        stop=math.floor(last/2.0)+perc2
                else:
                    stop=math.floor(last/2.0)
            else:
                start=math.floor(first/2.0)
                if last%2:
                    if perc2==1.0:
                        stop=math.ceil(last/2.0)
                    else:
                        stop=math.floor(last/2.0)+perc2
                else:
                    stop=math.floor(last/2.0)
            self.plot_block('highlight', 'heat.highlight.txt', 'eval(sprintf("%fr",conf(plots,track_start)+{0}*conf(plots,track_step)))'.format(start), 
                            'eval(sprintf("%fr",conf(plots,track_start)+{}*conf(plots,track_step)))'.format(stop), color, -1, '')
    
    def taxo_colors(self, no_of_ideograms):
        """
        Writes a block that defines the colors for the ideograms.
        """
        color_block='<colors>\n'
        for i in range(no_of_ideograms):
            color=self.ideogram_color(i)
            line='chr{}* = {}\n'.format(i,color)
            color_block+=line
        color_block+='</colors>\n\n'
        self.file.write(color_block)
        
    def taxo_legends_file(self,location):
        """
        Writes the legend text files for the taxonomy plot.
        """
        f=open('legends.0.txt','w')
        place,year=location.split('-')[0],location.split('-')[1]
        if len(place)<10:
            option='rpadding=20p'
        else:
            option=''
        option2='rpadding=-150p'
        f.write('k0\t0\t0\t{}\t{}\n'.format(year,option2))
        f.write('k0\t0\t0\t{}\t{}\n'.format(place,option))
        f.close()  
        
    def uneven_legends_file(self, legend):
        """
        Writes the legend text files for the taxonomy plot variant.
        legend contains the location and year.
        """
        f=open('legends.0.txt','w')
        legends=legend.split('-')
        if len(legends[0])>6:
            rpadding='rpadding=-100p'
        else:
            rpadding='rpadding=-50p'
        options=['rpadding=-50p',rpadding]
        if legends[2]!='':
            options.append('rpadding=-30p')
        i=0
        for option in options:
            f.write('k0\t0\t0\t{}\t{}\n'.format(legends[i],option))
            i+=1
        f.close()
        
    def text_block(self,filename, r0, r1, font, label_size, rpadding, padding, parallel, option):
        """Writes a text plot configuration block for the configuration file."""
        text1='<plot>\ntype  = text\nfile  = {}\nr1 = {}\nr0 = {}\nlabel_font = {}\nlabel_parallel={}\n'.format(filename, r1, r0, font, parallel)
        text2='label_size = {}p\nz=3\nrpadding={}\npadding={}\n{}\n</plot>\n\n'.format(label_size,rpadding,padding,option)  
        self.file.write(text1+text2)
            
    def ideogram_color(self, i):
        """
        Gives a color.
        """
        length=len(self.colors)
        if i<length:
            return self.colors[i]
        elif i<(length*2):
            return self.colors[i-length][1:]
        elif i<(length*3):
            return 'd'+self.colors[i-length*2][1:]
        else:
            return self.colors[i%2]
    
    def ideogram(self, spacing, radius, thickness, label_size, l_radius, include_rules):
        """
        Writes the ideogram block for the circos configuration file.
        """
        i_start='<ideogram>\n<spacing>\ndefault = {}r\n</spacing>\nradius={}r\nthickness={}p\n'.format(spacing,radius,thickness)
        i_labels='show_label = yes\nlabel_font = default\nlabel_radius={}\nlabel_size= {}\nlabel_parallel= yes\n'.format(l_radius,label_size)
        i_end='fill = yes\nstroke_color= dgrey\nstroke_thickness = 2p\n</ideogram>\n\n'
        ideogram=i_start+i_labels+include_rules+i_end
        self.file.write(ideogram)
        
    def file_end(self,svg=True):
        """Writes some neccessary things to the circos configuration file. If argument svg == False, no .svg image will be created, only a .png
        auto_alpha_steps is for transparency used in heatmaps for taxonomy plots.
        """
        if svg==True:
            comment='#'
        else:
            comment=''
        image_block='\n<image>\n<<include etc/image.conf>>\n{}svg* = no\nauto_alpha_steps* = 7\n</image>\n\n'.format(comment)
        last_blocks='<<include etc/colors_fonts_patterns.conf>>\n<<include etc/housekeeping.conf>>\n' 
        extra=r'file_delim* = \t'
        self.file.write(image_block+last_blocks+extra+'\n')
        
    def file_end_season(self,depth):
        """
        Writes some neccessary things to the circos configuration file. 
        This one includes a background image which is for the folder tree plots. 
        The image contains the plots that are the same for all plots, so that they only need be plotted once.
        """
        image_block='\n<image>\n<<include etc/image.conf>>\nsvg* = no\nbackground* = ./{}background.png\n</image>\n\n'.format('../'*depth)
        last_blocks='<<include etc/colors_fonts_patterns.conf>>\n<<include etc/housekeeping.conf>>\n' 
        extra=r'file_delim* = \t'
        self.file.write(image_block+last_blocks+extra+'\n')
        
    def karyotype(self,filename):
        karyotype='karyotype = {}\n\n'.format(filename)
        self.file.write(karyotype)
            
    def ticks_week(self):
        """
        Writes the tick file for the year circos plot.
        """
        f=open('ticks.conf','w')
        ticks='show_ticks = yes\n<ticks>\nradius = 1r\ncolor = dgrey\nthickness = 2p\n<tick>\nspacing = 7\nsize = 10p\n</tick>\n</ticks>\n'
        f.write(ticks)
        f.close()
        
    def season_block(self,depth,seasons):
        """
        Writes the weather plots for the year circos plot:
        Season track, temperature histogram, precipitation histogram, wind(+wind direction) histogram, compass rose highlight plot and the axis for the read heatmaps.
        """
        r0=0.8        
        r1=0.65
        r=0.45
        if seasons==True:
            s_rule=self.season_rules()
            self.plot_block('highlight', './{}highlight.seasons.txt'.format('../'*depth), 'dims(ideogram,radius_inner)', '{}r'.format(r0), 'red', 1, s_rule+'\nstroke_thickness=0\n')

        temp_rules='<rules>\n<rule>\ncondition  = var(value) <0\nfill_color       = vdblue\n</rule>\n</rules>\n'    
        axes='<<include ./{}axes.conf>>\n'.format('../'*depth)
        self.plot_block('histogram', './{}temperature.txt'.format('../'*depth), '{}r'.format(r0), 'dims(ideogram,radius_inner)', 'dred', 2, axes+temp_rules+'min=-20\nmax=30\nthickness=0p\n')
        self.plot_block('histogram', './{}precipitation.txt'.format('../'*depth), '{}r'.format(r1), '{}r'.format(r0), 'pblue', 2, 'orientation=in\n')
    
        w_rule=self.wind_rules()
        self.plot_block('histogram', './{}wind.txt'.format('../'*depth), '{}r'.format(r), '{}r'.format(r1), 'pblue', 2, w_rule)
        self.compass_rose(depth)
        
    def season_rules(self):
        """
        Returns a string with the rules used to color the season highlight.
        """
        s_rule='<rules>\n'
        seasons=[('1','pastel1-9-qual-2'),('2','pastel1-9-qual-6'),('3','pastel1-9-qual-3'),('4','pastel1-9-qual-5')]
        for season in seasons:
            s_rule+='<rule>\ncondition = var(id) == {0}\nfill_color={1}\ncolor={1}\n</rule>\n'.format(season[0],season[1])
        s_rule+='</rules>\n'
        return s_rule
    
    def wind_rules(self):
        """
        Returns a string with the rules used to color the wind histogram.
        """
        rules=[('north','my_blue'),('northwest','greenblue'),('west','my_green'),('southwest','greenyellow'),('south','my_yellow'),('southeast','my_orange'),('east','my_red'),('northeast','my_purple')]
        w_rule='<rules>\n'
        for rule in rules:
            w_rule+='<rule>\ncondition = var(id) eq "{0}"\nfill_color	= {1}\ncolor	= {1}\n</rule>\n\n'.format(rule[0],rule[1])
        w_rule+='</rules>\n\n'        
        return w_rule
        
    def temp_axes_file(self,):
        """
        Writes the axes file for the temperature histogram.
        """
        f=open('axes.conf','w')
        a1='<axes>\nthickness = 1\ncolor     = dgrey\n<axis>\nposition  = -10\n</axis>\n<axis>\nposition  = 0\n'
        a2='</axis>\n<axis>\nposition  = 10\n</axis>\n<axis>\nposition  = 20\n</axis>\n</axes>\n'
        f.write(a1+a2)
        f.close()
    
    def year_legend_plots(self,names,maximums,depth,name):
        """
        Calls the functions that write the plot blocks for the year circos plot.
        """
        
        placement=self.reads_heatmaps(names,maximums,depth)
        if name:
            filename1='./legends.{}.0.txt'.format(name)
            filename2='./legends.{}.1.txt'.format(name)
        else:
            filename1='./legends.0.txt'
            filename2='./legends.1.txt'
        if len(maximums)<14:
            label_size=20
        else:
            label_size=15

        self.text_block(filename1, '0.25r', '1.0r', 'light', 35, '0p', '0p', 'yes', '')
        self.text_block(filename2, '1.0r', '2.5r', 'light', label_size, '0p', '0p', 'yes', '')
        return placement

    def reads_heatmaps(self,names,maximums,depth):
        """
        Writes the plot blocks for the read heatmaps.
        """
        #print(names)
        placements,upper,width=self.read_tracks(maximums)
        self.read_axis(upper,width,depth)
        text_place=self.heatmaps_taxa(names,placements)
        self.white_highlight(depth)
        return text_place
    
    def heatmaps_taxa(self,names,placements):
        """
        Writes the plot blocks for the read heatmaps for the year plot into the main configuration file.
        """
        colors=['blues','reds']
        text_place=[]
        for i in range(len(placements)):
            p0=placements[i][0]
            p1=placements[i][1]
            self.plot_block('heatmap', './{}.txt'.format(names[i]), '1.1r+{}p'.format(p0), '1.1r+{}p'.format(p1), '{}-8-seq'.format(colors[i%2]), 1, '')
            text_place.append(p1)
        return text_place
        
    def labels_reads_file(self, path, organisms,maximums,placements,name):
        """
        Writes the text file with the names of the taxa displayed at the right distance from the center.
        """
        if name:
            filename=os.path.join(path,'legends.{}.1.txt'.format(name))
        else:
            filename=os.path.join(path,'legends.1.txt')
        f=open(filename,'w')
        length=len(organisms)
        pos=[0,]
        for i in range(1,int(length/2)+1):
            pos.extend([i,-i])
        for j in range(length):
            f.write('week{}\t0\t0\t{}\tr0=1.1r+{}p\n'.format(27-pos[j],organisms[j],placements[j]))
        f.close()
        
    def read_axis(self,upper,width,depth):
        spacing=1.0/float(upper)
        self.write_axis(width,depth)
        self.axis_labels(upper,width,spacing,depth)
            
    def axis_labels(self,upper,width,spacing,depth):
        """
        Writes the axis block for the heatmaps.
        """
        distance=round(width*spacing)
        option='label_snuggle=yes\nmax_snuggle_distance  = 3r'
        self.text_block('./{}text_axis.{}.txt'.format('../'*depth,upper),'0.15r+{}p'.format(distance*upper),'0.15r+{}p'.format(distance*upper+70), 'light', 25,'0p','0p', 'yes', option)
            
    def reads_axis_files(self):
        """Writes the axis file for the heatmap axis."""
        for i in range(3,9):
            filename='text_axis.{}.txt'.format(i)
            with open(filename, "w") as myfile:
                myfile.write('week10\t0\t0\t{}\n'.format('10^{}'.format(i)))        

    def write_axis(self,width,depth):
        """
        Writes the empty plot block used to give the axes for the heatmap widths.
        """
        axes='<axes>\nthickness = 1\ncolor     = dgrey\n<axis>\nspacing   = 1r\n</axis>\n</axes>\n'
        self.plot_block('histogram', './{}empty.txt'.format('../'*depth), '0.15r', '0.15r+{}p'.format(width), 'dgrey', -1, axes)             
    
    def read_tracks(self,maximums):
        """
        Calculates the width of the read heatmaps based on the maximum number of reads. 
        The maximums are logged so that the differences won't be ridiculus.
        """
        label_space=30
        if len(maximums)>14:
            total_width=675
            label_space=25
        elif len(maximums)>3:
            total_width=675
        elif len(maximums)==3:
            total_width=400
        elif len(maximums)==2:
            total_width=300
        else:
            total_width=150
        
        length=len(maximums)
        maximum=max(maximums)
        minimum=min(maximums)
        space=[]
        plot_space=total_width-(label_space*length)
        maximums=[math.log10(x) for x in maximums]
        total=sum(maximums)
        upper,width=self.axis_values(maximum,minimum,total,plot_space)
        maximums=[int(plot_space*(maximum/total)) for maximum in maximums]
        tot_pixels=0
        for maximum in maximums:
            space.append((tot_pixels,tot_pixels+maximum))
            tot_pixels+=maximum+label_space
        return space,upper,width
    
    def axis_values(self,maximum,minimum,total,plot_space):
        """
        Calculates the width neccessary for the heatmap axis.
        """
        upper=math.floor(math.log10(maximum))
        width=int(plot_space*(upper/total))
        return upper,width
    
    def season_file_legends(self, location, level, name, path):
        """
        Writes the text files for the legends for the year circos plot.
        """

        filename=os.path.join(path,'legends.{}.0.txt'.format(name))         
        f=open(filename,'w')
        split=location.split('-')
        place_year=split[0]+' '+split[1]
        f.write('week1\t0\t0\t{}\trpadding=10p\n'.format(place_year))
        level=self.taxonomy_names(level)
        option='rpadding=-50p'
        if len(name)<11:
            option2='rpadding=-50p'
        else:
            option2='rpadding=10p,label_snuggle=yes,max_snuggle_distance=3r'
        f.write('week27\t0\t0\t{}\t{}\n'.format(level,option))
        if name:
            f.write('week27\t0\t0\t{}\t{}\n'.format(name,option2))
        f.close()
        
    def season_file_legends_s(self, location, level, name, path):
        """
        Writes the text files for the legends for the year circos plot.
        """
        if name:
            filename=os.path.join(path,'legends.{}.0.txt'.format(name))
        else:
            filename=os.path.join(path,'legends.0.txt')           
        f=open(filename,'w')
        if level:
            level=self.taxonomy_names(level)
        option='rpadding=-50p'
        if len(name)<11:
            option2='rpadding=-75p'
        else:
            option2='rpadding=10p,label_snuggle=yes,max_snuggle_distance=3r'	
        split=location.split('-')
        place_year=split[0]+' '+split[1]
        if level:
            f.write('week27\t0\t0\t{}\t{}\n'.format(level,option))
        if name:
            f.write('week27\t0\t0\t{}\t{}\n'.format(name,option2))
        f.write('week1\t0\t0\t{}\trpadding=10p\n'.format(place_year))
        f.close()
        
    def year_background_legend(self, location):
        """
        Writes the text files for the legends for the year circos plot.
        """
        f=open('place.txt','w')
        split=location.split('-')
        place_year=split[0]+' '+split[1]
        f.write('week1\t0\t0\t{}\trpadding=10p\n'.format(place_year))
        f.close()
        
    def taxonomy_names(self,letter):
        taxonomy_names={'k':'Kingdom', 'p':'Phylum','o':'Order','c':'Class','f':'Family','g':'Genus', 's':'Species'}
        return taxonomy_names[letter]
    
    def find_depth(self,level):
        """Finds the depth in the folder system for the year circos plot."""
        depth=self.ranks.index(level)+1
        return depth
    
    def year_base(self,read_totals):
        """
        Used to write the configuration files that don't need to be remade for each folder.
        """ 
        colors=self.greys(read_totals)
        self.karyotype_week(colors)
        self.own_colors()
        self.compass_file()
        self.ticks_week()
        self.temp_axes_file()
        self.reads_axis_files()
        f=open('empty.txt','w')
        f.close()

    def year_conf(self, organisms, maximums, path, location, level, name):
        """
        Writes the configuration file for the year plot.
        """
        self.file=open(os.path.join(path,'circos.conf'),'w')
        depth=self.find_depth(level)
        self.karyotype('./{}karyotype.weeks.txt'.format('../'*depth))
        self.file.write('<<include ./{}colors.conf>>\n'.format('../'*depth))
        self.file.write('<plots>\n')
        placements=self.year_legend_plots(organisms, maximums, depth,name)
        self.file.write('</plots>\n\n')
        self.labels_reads_file(path, organisms, maximums, placements,name)
        self.season_file_legends(location, level, name, path) 
        self.ideogram(0, 0.5, 20, 30, '1.02r', 'show=no\n')
        self.file_end_season(depth)
        self.file.close()
        
    def year_conf_background(self,seasons):
        """
        Writes the configuration file for a year plot without any reads.
        Used for the hierarchal run as background instead of generating the weather plots over and over again.
        """
        self.file=open('background.conf','w')
        depth=0
        self.karyotype('karyotype.weeks.txt')
        self.file.write('<<include colors.conf>>\n')
        self.file.write('<plots>\n')
        self.season_block(depth,seasons)
        self.file.write('</plots>\n')
        self.ideogram(0, 0.5, 20, 30, '1.02r', '')
        self.file.write('<<include ticks.conf>>\n')
        self.file_end(False)
        self.file.close()
        
    def year_conf_s(self, organisms, maximums,location, level, name, seasons):
        """
        Writes the configuration file for a single year.
        """
        depth=0
        path=''
        self.file=open('circos.conf','w')
        self.karyotype('karyotype.weeks.txt')
        #self.own_colors()
        self.file.write('<<include colors.conf>>\n')
        self.file.write('<plots>\n')

        self.season_block(depth,seasons)
        placements=self.year_legend_plots(organisms, maximums, depth,name)
        self.file.write('</plots>\n\n')
        self.labels_reads_file(path, organisms, maximums, placements,name)
        self.season_file_legends_s(location, level, name, path) 
        self.ideogram(0, 0.5, 20, 30, '1.02r', '')
        self.file.write('<<include ticks.conf>>\n')
        self.file_end()
        self.file.close()
        
    def own_colors(self):
        """
        Writes a file with the definition of the colors used for the directions in the wind plot.
        """
        c1='<colors>\ngreenblue= 0, 140, 110\nmy_purple= 160, 0, 160\nmy_green=10,215,10\nmy_blue=45,45,255\nmy_red=255,30,30\nmy_yellow=255,255,0\n'
        c2='greenyellow=165,225,0\nmy_orange=255,110,0\ngrey1=247,247,247\ngrey2=217,217,217\ngrey3=189,189,189\ngrey4=150,150,150\ngrey5=99,99,99\n'
        c3='grey6=37,37,37\n</colors>\n'
        color_block=c1+c2+c3
        f=open('colors.conf','w')
        f.write(color_block)
        f.close()
        
    def white_highlight(self,depth):
        """
        Writes a white highlight plot block. It is used to cover the heatmap axis so that it only shows up on a part of the circle and not all the way round.
        (Not possible to only display part of axis in Circos.)
        """
        
        rules='<rules>\n<rule>\ncondition = on(week10)||on(week11)||on(week12)||on(week13)||on(week14)||on(week15)||on(week16)||on(week17)\nshow=no\n</rule>\n</rules>\n'
        self.plot_block('highlight','./{}highlight.seasons.txt'.format('../'*depth),'0.1r', '0.45r', 'white', 0, rules)
        
    def compass_file(self):
        """
        Writes the compass rose highlight file.
        """
        f=open('compass_rose.txt','w')
        weeks=[(1,0,7),(6,4,7),(7,0,7),(8,0,3),(13,0,7),(14,0,7),(19,4,7),(20,0,7),(21,0,3),
               (26,0,7),(27,0,7),(32,4,7),(33,0,7),(34,0,3),(39,0,7),(40,0,7),(45,4,7),(46,0,7),(47,0,3),(52,0,7)]
        for week in weeks:
            f.write('week{}\t{}\t{}\n'.format(week[0],week[1],week[2]))
        f.close()
        
    def greys(self, totals):
        """
        Calculates what color the weeks that have corresponding airfilters should have in the heatmap that displays total reads.
        Not elegant, but the result looks the same as the heatmap of the totals.
        """
        number=6
        
        colors=[]
        maximum=max(totals)
        minimum=min(totals)
        tot=maximum-minimum
        part=tot/number
        stops=[]
        for i in range(number-1):
            minimum+=part
            stops.append(minimum)
        minimum=min(totals)
        for item in totals:
            if item<=stops[0]:
                colors.append('grey1')
            elif item>stops[0]and item<=stops[1]:
                colors.append('grey2')
            elif item>stops[1]and item<=stops[2]:
                colors.append('grey3')
            elif item>stops[2]and item<=stops[3]:
                colors.append('grey4')
            elif item>stops[3]and item<=stops[4]:
                colors.append('grey5')
            elif item>stops[4]:
                colors.append('grey6')
            else:
                print('greys: why??',item)
        return colors
    
    def karyotype_week(self, colors):
        """
        Writes the year 'karyotype' file.
        """
        f=open('karyotype.weeks.txt','w')
        j=0
        for i in range(1,53):
            if i%2:
                f.write('chr\t-\tweek{0}\t{0}\t0\t7\t{1}\n'.format(i,colors[j]))
                j+=1
            else:
                f.write('chr\t-\tweek{0}\t{0}\t0\t7\twhite\n'.format(i))
        f.close()
        
    def compass_rose(self,depth):
        """
        Writes the compass rose block for the year plot into the configuration file.
        """
        first=[("40","39",'my_green'),("1","52",'my_blue'),("26","27",'my_yellow'),('13','14','my_red')]
        second=[(33,34,32,'greenyellow'),(46,47,45,'greenblue'),(20,19,21,'my_orange'),(7,6,8,'my_purple')]
        rules='<rules>\n'
        for direction in first:
            rules+='<rule>\ncondition = var(chr) eq "week{0}"| var(chr) eq "week{1}"\nfill_color	= {2}\ncolor	= {2}\n</rule>\n'.format(direction[0],direction[1],direction[2])
        for color in second:
            rules+='<rule>\ncondition = var(chr) eq "week{0}" | var(chr) eq "week{1}"| var(chr) eq "week{2}"\nfill_color	= {3}\ncolor	= {3}\n</rule>\n'.format(color[0],color[1],color[2],color[3])
        rules+='</rules>\n'
        self.plot_block('highlight', './{}compass_rose.txt'.format('../'*depth), '0.01r', '0.10r', 'grey', 0, rules)
        
    def plot_block(self, plot_type, file_name, r0, r1, color, z, rules):
        plot_block='<plot>\ntype={0}\nfile={1}\nr0 = {2}\nr1 = {3}\ncolor={4}\nfill_color={4}\n{5}\nz  = {6}\n</plot>\n\n'.format(plot_type, file_name, r0, r1, color, rules, z)
        self.file.write(plot_block)
        
if __name__ == "__main__":
    a=CircosConf()





