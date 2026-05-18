#sampledescription_Replicatenumber_EarlyorLate_Illumina-Lane_Illumina-read_filenumber.fastq.gz.
import os
import seaborn as sns
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
import warnings
warnings.filterwarnings('ignore')
configfile: "config.yaml"
samples_folder = config['samples_folder']
samples_folder = samples_folder+'/'
output_folder = config['output_folder']
output_folder = output_folder+'/'
target = config['target']
target_list = []


for value in target:
 
    target_list.append(value.split(".fq")[0].split(".fastq")[0].split('_'))
    #print(target_list[-1])
    
#print(target)
#print(expand("{targets}.clip", targets = target))

output_list = []
stats_list = []
figure_list = []
acf_list = []
rpkm_list = []
raw_kde = []
genome = config['genome'].split('.fa')[0]
genomeid = config['genome'].split('/')[-1].split('.fa')[0]

barcode1= str(config['barcode1']) #AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
barcode2= str(config['barcode2']) #AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

for file in target_list:
    temp = '_'.join(file[0:1])
    
    temp = '_'.join(file[0:2])
    raw_kde = raw_kde + expand(output_folder+'figures/raw/'+temp+'_{sizes}_'+genomeid+'_log2RT_rawkdeplot.png',  sizes = config['window_sizes'])
    output_list = output_list + expand(output_folder+'bg/norm/'+temp+'_{sizes}_'+genomeid+'_RTnorm.bg', sizes = config['window_sizes'])
    acf_list = acf_list + expand(output_folder+'figures/acf/'+temp+'_{sizes}_'+genomeid+'_log2RT.svg',  sizes = config['window_sizes'])
    
    temp = '_'.join(file[0:3])
    stats_list.append(output_folder+'bam/stats/rmdup/'+temp+".rdbamstats")
    stats_list.append(output_folder+'bam/stats/fbstats/'+temp+".fbstats")
    figure_list = figure_list + expand(output_folder+'figures/read_count/all/raw/'+temp+'.csv',  sizes = config['window_sizes'])
    rpkm_list = rpkm_list + expand(output_folder+'bg/rpkm/'+temp+'_{sizes}_'+genomeid+'_rpkm.bg',  sizes = config['window_sizes'])
    

    #output_list.append('samples/'+"_".join(file))
output_list = list(set(output_list))







# for i in range(0,len(output_list)):
#     if len(output_list[i].split("_E_")) > 1:
#         output_list[i] = output_list[i].split('_E_')[0]+output_list[i].split('_E_')[1]
#     elif len(output_list[i].split("_L_")) > 1:
#         output_list[i] = output_list[i].split('_L_')[0]+output_list[i].split('_L_')[1]
#     else:
#         output_list
#         #print("error in input")
# output_list = list(set(output_list))
#print(output_list)
rule all:
    input:        
        output_list,
        stats_list,
        figure_list,
        acf_list,
        raw_kde,
        output_folder+'figures/log2/log2RT_kdeplot.png',
        
        output_folder+'figures/read_count/reads.png',
        output_folder+'figures/windows/window_counts.png',
        
        
rule raw_and_acfplot:
    resources:        
        threads=1,
        mem_mb=4000,
        runtime=60
    input:
        output_folder+'bg/log2/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT.bg',

    output:
        svg = output_folder+'figures/acf/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT.svg',
        png = output_folder+'figures/acf/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT.png',
        svg_raw = output_folder+'figures/raw/10Mb/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT_10Mb_raw.svg',
        png_raw = output_folder+'figures/raw/10Mb/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT_10Mb_raw.png',
        svg_raw_2Mb = output_folder+'figures/raw/2Mb/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT_2Mb_raw.svg',
        png_raw_2Mb = output_folder+'figures/raw/2Mb/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT_2Mb_raw.png',
    
    
    run:
        df = pd.read_csv(str(input), sep = '\t', header = None)
        df.columns = ['chr', 'start', 'end', 'RT']
        ax = pd.plotting.autocorrelation_plot(df['RT'])
        lag2 = ax.lines[5].get_ydata()[1]
        plt.axhline(y = lag2, color = 'r', linestyle = '-')
        ax.set_xlim([0, 1000])
        ax.set_ylim([0, 1])
        plt.text(650, 0.81, "ACF = " + str(round(lag2, 2)), size='x-large')
        plt.savefig(output.png)
        plt.savefig(output.svg)
        plt.close()
        chromo=list(df['chr'])[0]
        df = df.loc[df['chr'] == chromo]
        start = 5*1000000
        end = 105*1000000
        df = df.loc[df['start'] >= start]
        df = df.loc[df['end'] <= end]
        end = 25*1000000
        df2 = df.loc[df['end'] <= end]
        df['start'] = df['start']/1000000
        df['end'] = df['end']/1000000
        df2['start'] = df2['start']/1000000
        df2['end'] = df2['end']/1000000
        ax = sns.scatterplot(data = df, x = "start", y = 'RT', s=3, linewidth=0)
        ax.set(xlabel='NC_000001.11' + " (Mb)")
        plt.ylim(-7.5, 7.5)
        plt.axhline(y = 0, color = 'b', linestyle = '-')
        plt.savefig(output.png_raw)
        plt.savefig(output.svg_raw)
        plt.close()
        
        ax = sns.scatterplot(data = df2, x = "start", y = 'RT', s=3, linewidth=0)
        ax.set(xlabel=str(chromo) + " (Mb)")
        plt.ylim(-7.5, 7.5)
        plt.axhline(y = 0, color = 'b', linestyle = '-')
        plt.savefig(output.png_raw_2Mb)
        plt.savefig(output.svg_raw_2Mb)
    
        
rule combine_figures:
    input:
        figure_list
    output:
        svg = output_folder+'figures/read_count/reads.svg',
        png = output_folder+'figures/read_count/reads.png',
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
        
    run:
        combineddf = pd.DataFrame()
        for figure in input:
            fig_file = pd.read_csv(figure, header=0)
            #fig_file.columns = ['name', 'Raw']
            combineddf = pd.concat([combineddf, fig_file])
        #combineddf.columns = ['name', 'Raw']
        names = list(combineddf['name'])
        desc = []
        rep = []
        earlylate = []
        for value in names:
            temp = value.split(r'/')[-1].split('_')
            desc.append(temp[0])
            rep.append(temp[1])
            earlylate.append(temp[2].split('.')[0])
        combineddf['desc'] = desc
        combineddf['rep'] = rep
        combineddf['el'] = earlylate
        plotnumber = len(set(desc))*2
        descs = list(set(desc))
        fig,ax=plt.subplots(nrows=plotnumber, figsize=(5,plotnumber*5))
        #subfigs = fig.subfigures(1, len(descs))
        i = 0
        while i <len(set(desc))*2:
            try:
                plotdf = combineddf[combineddf['desc']==descs[int(i/2)]]
                plotdf = plotdf.drop(labels=['Unnamed: 0', 'desc','index','name'], axis = 1)
                melted = plotdf.melt(id_vars = ['el', 'rep'])
                melted1 = melted[melted['el'] == 'E']
                barplot = sns.barplot(data = melted1, ax = fig.axes[i], y = 'variable', x = 'value', hue = 'rep')
                barplot.legend_.remove()
                barplot.set(xlabel=None)
                #barplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)
                #barplot.xticks(rotation=45)
                fig.axes[i].set_title(descs[int(i/2)]+" - Early")
                fig.axes[i].set_yticklabels(['Raw', 'Mapped', 'Filtered', 'Deduplicated'], rotation=0)
                melted2 = melted[melted['el'] == 'L']
                barplot = sns.barplot(data = melted2, ax = fig.axes[i+1], y = 'variable', x = 'value', hue = 'rep')
                barplot.set(yticklabels=[])
                barplot.set(ylabel=None)
                barplot.set(xlabel=None)
                barplot.legend_.remove()
                #barplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)
                fig.axes[i+1].set_yticklabels(['Raw', 'Mapped', 'Filtered', 'Deduplicated'], rotation=0)
                fig.axes[i+1].set_title(descs[int(i/2)]+" - Late")
                #barplot.xticks(rotation=45)
            except:
                i=i
    
    
            i+=2

        
        #for i in range(0,len(subfigs)):
        #    subfigs[i].suptitle(descs[i])
        plt.tight_layout()
        
        plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
        plt.savefig(output.svg, bbox_inches='tight',pad_inches=0.1)
        plt.close()

        
rule rawKDE:
    input:
        output_folder+'bg/log2/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT.bg',
    output:
        png = output_folder+'figures/raw/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT_rawkdeplot.png',

    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    run:
        df = pd.DataFrame()
        
        
        temp_df = pd.DataFrame()
        path = str(input)
        temp_df = pd.read_csv(path, sep = '\t', header = None)
        temp_df.columns = ['chr', 'start', 'end', 'RT']
        temp_df['name'] = value.split(r'/')[-1].split('.bg')[0]
        ax = sns.kdeplot(data=temp_df, x="RT")
        #sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.4))
        plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
        plt.savefig(str(output.png).split('.png')[0]+'.svg', bbox_inches='tight', pad_inches=0.1)
        plt.close()
      

rule KDElog2Plot:
    input:
        output_list
    output:
        svg = output_folder+'figures/log2/log2RT_kdeplot.svg',
        png = output_folder+'figures/log2/log2RT_kdeplot.png',
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    run:
        df = pd.DataFrame()
        sizes = config['window_sizes']
        for value in output_list:
            temp_df = pd.DataFrame()
            
            temp_df = pd.read_csv(value, sep = '\t', header = None)
            temp_df.columns = ['chr', 'start', 'end', 'RT']
            temp_df['name'] = value.split(r'/')[-1].split('.bg')[0]
            ax = sns.kdeplot(data=temp_df, x="RT")
            #sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.4))
            plt.savefig(output_folder+'figures/log2/'+value.split(r'/')[-1].split('RTnorm.bg')[0]+'log2RT.png', bbox_inches='tight', pad_inches=0.1)
            plt.savefig(output_folder+'figures/log2/'+value.split(r'/')[-1].split('RTnorm.bg')[0]+'log2RT.svg', bbox_inches='tight', pad_inches=0.1)
            plt.close()
            if len(df) == 0:
                df = temp_df
            else:
                df = pd.concat([df, temp_df], axis=0)
        df.columns=['chr', 'start', 'end', 'RT', 'name']
        ax = sns.kdeplot(data=df, x="RT", hue = 'name')
        sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.4))
        plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
        plt.savefig(output.svg, bbox_inches='tight', pad_inches=0.1)
        plt.close()
        
rule plot_read_counts:
    input:
        clip = output_folder+'clip/{desc}_{rep}_{EL}_R1.clip',
        bam = output_folder+'bam/raw/{desc}_{rep}_{EL}.bam',
        filtered = output_folder+'bam/filtered/{desc}_{rep}_{EL}.filtered',
        rmdup = output_folder+'bam/rmdup/{desc}_{rep}_{EL}_rmdup.bam',
    output:
        png = output_folder+'figures/read_count/all/{desc}_{rep}_{EL}.png',
        svg = output_folder+'figures/read_count/all/{desc}_{rep}_{EL}.svg',
        csv = output_folder+'figures/read_count/all/raw/{desc}_{rep}_{EL}.csv'
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    
        
    run:
        name = input.rmdup.split('rmdup/')[1].split('.bam')[0]
        raw_count = int(os.popen("echo $(cat '"+input.clip+"'|wc -l)/4|bc").read())
        bam_count = int(os.popen("samtools view -c '"+input.bam+"'").read())
        filtered_count = int(os.popen("samtools view -c '"+input.filtered+"'").read())        
        rmdup_count = int(os.popen("samtools view -c '"+input.rmdup+"'").read())        
        if len(config['reads']) == 2:
            bam_count = bam_count*2
            filtered_count = filtered_count*2
            rmdup_count = rmdup_count*2
        if raw_count != 0:
            ylist = [raw_count/raw_count, bam_count/raw_count, filtered_count/raw_count, rmdup_count/raw_count]
            xlist = ['Raw:\n'+str(raw_count), 'Mapped:\n'+str(bam_count), 'Filtered:\n'+str(filtered_count), 'Deduplicated:\n'+str(rmdup_count)]
            xlistnumb=[raw_count, bam_count, filtered_count, rmdup_count]
        else:
            ylist = [0,0,0,0]
            xlist = [0,0,0,0]
            xlistnumb = [0,0,0,0]
        fig = sns.barplot(x = xlist, y = ylist)
        plt.savefig(output.png)
        plt.savefig(output.svg)
        plt.close()
        df = pd.DataFrame()
        df['ratio'] = ylist
        
        #df['count'] = xlistnumb
        df = df.T
        df.columns = ['Raw Count', 'Bam Count', 'Filtered Count', 'Deduplicated Count']
        df['name']=name
        df=df.reset_index()
        df.to_csv(output.csv)
    
rule plot_window_counts:
    input:
        rpkm_list
    output:
        png = output_folder+'figures/windows/window_counts.png',
        svg = output_folder+'figures/windows/window_counts.svg',
    threads: 1
    resources:
        mem_mb=4000,
        runtime=60
    
    run:
        plotdf = pd.DataFrame()
        covered = []
        none = []
        plotdf_index = []
        file_list = list(input)
        for file in file_list:
            fileplot = pd.DataFrame()
            file_cov = []
            file_none = []
            file_idx = []
            df = pd.read_csv(str(file),  sep = '\t', header = None)
            df.columns = ['chr', 'start', 'end', 'RT']
            df['label'] = '1'
            df['covered'] = [0 if x <= 0 else 1 for x in df['RT']]
            non = df['covered'].value_counts()[0]
            cov = df['covered'].value_counts()[1]
            tot = non + cov
            non = non/(tot)
            cov = cov/tot        
            covered.append(cov)
            none.append(non)
            plotdf_index.append('_'.join(str(file).split(r'/')[-1].split('_')[0:3]))
            file_cov.append(cov)
            file_none.append(non)
            file_idx.append('_'.join(str(file).split(r'/')[-1].split('_')[0:3]))
            fileplot['covered']  = file_cov
            fileplot['none'] = file_none
            fileplot.index = file_idx
            ax = fileplot.plot.barh(rot=0, stacked = True)
            plt.axvline(x = 0.8, color = 'r', linestyle = '-')
            #sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.25))
            plt.savefig(output_folder+'figures/windows/'+'_'.join(str(file).split(r'/')[-1].split('_')[0:3])+'.png', bbox_inches='tight', pad_inches=0.1)
            plt.savefig(output_folder+'figures/windows/'+'_'.join(str(file).split(r'/')[-1].split('_')[0:3])+'.svg', bbox_inches='tight', pad_inches=0.1)
            plt.close()
        plotdf['covered'] = covered
        plotdf['none'] = none
        plotdf.index = plotdf_index
        ax = plotdf.plot.barh(rot=0, stacked = True)
        plt.axvline(x = 0.8, color = 'r', linestyle = '-')
        sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.25))
        plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
        plt.savefig(output.svg, bbox_inches='tight', pad_inches=0.1)

#def get_bwa_map_input_fastqs(wildcards):
#    return config["samples"][wildcards.sample]

#rule bcftools_call:
#    input:
#        fa="data/genome.fa",
#        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
#        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
#    output:
#        "calls/all.vcf"
#    shell:
#        "bcftools mpileup -f {input.fa} {input.bam} | "
#        "bcftools call -mv - > {output}"

rule cutadapt:
    input:
        samples_folder+'{desc}_{rep}_{EL}_{reads}.fastq.gz',
        #fq = expand('samples/{sample}.fastq.gz', sample = config["samples"]),
        #get_bwa_map_input_fastqs
    output:
        clip=output_folder+'clip/{desc}_{rep}_{EL}_{reads}.clip',
        txt=output_folder+'clip/{desc}_{rep}_{EL}_{reads}.txt',
        #clip = expand('clip/{sample}.clip', sample = config["samples"]),
        #'clip/{sample}.clip'
    threads: 4
    resources:
        mem_mb=64000,
        runtime=480
    run:
        #print(len(config['reads']))
        #print(input)
        output_string = output.clip
        input_string = input
        if len(config['reads']) > 1:
            if '_R1.fastq' in str(input):
                #print("r1")                
                shell("cutadapt -a "+barcode1+" -q 0 -O 1 -m 0 -j 0 -o '"+output_string+"' '"+input_string+"'")
                with open(str(output.txt), 'w') as f:
                     f.write("clip success")
                f.close()

            if '_R2.fastq' in str(input):
                #print("r2")
                shell("cutadapt -a "+barcode2+" -q 0 -O 1 -m 0 -j 0 -o '"+output_string+"' '"+input_string+"'")
                with open(str(output.txt), 'w') as f:
                     f.write("clip success")
                f.close()
        
        else:
            shell("cutadapt -a "+barcode1+" -q 0 -O 1 -m 0 -j 0 -o '"+output_string+"' '"+input_string+"'")
            with open(str(output.txt), 'w') as f:
                f.write("clip success")
            f.close()


        
rule genome_unzip:
    input:
        expand(genome+'.fa.gz'),
    output:
        expand(genome+'.fa'),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
        
    shell:
        "gunzip -c '{input}' > '{output}'"


        
rule genome_index:
    input:
        expand(genome+'.fa'),
    output:
        expand(genome+'.fa.amb'),
    threads: 1
    resources:
        mem_mb=64000    
    run:
        shell("bwa index '{input}'"),

        
        
rule faidx:
    input:
        fa=expand(genome+'.fa'),
    output:        
        expand(genome+'.fa.fai'),
        
    threads: 4
    resources:
        mem_mb=16000    
    run:
        
        fa = str(input.fa)
        shell("samtools faidx '"+fa+"'")
        
rule sizes:
    input:
        expand(genome+'.fa.fai'),
    output:
        sizes = expand(genome+'.fa.fai.sizes'),
    threads: 1
    resources:
        mem_mb=8000    
    run:
        shell("cut -f1,2 '{input}' > '{output.sizes}'")
        

        
rule makewindows:
    input:
        expand(genome+'.fa.fai.sizes'),
        #sizes = expand('{sizes}', sizes = config['window_sizes'])
    output:
        file = expand(output_folder+'bed/{sizes}_windows_'+genomeid+'.bed', sizes = config['window_sizes']),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
        
    run:
        commands = expand("bedtools makewindows -w {sizes} -s {sizes} -g {{input}} > '"+output_folder+"bed/{sizes}_windows_"+genomeid+".bed'", sizes = config['window_sizes'])
        for i in range(0,len(commands)):
            shell(commands[i])          
  
        
        
        
        
        
        
        
        
        
        
# rule makewindows:
#     input:
#         expand('{genome}.fa.fai.sizes', genome = config['genome']),
        
#     output:
#         fivekb = expand('{genome}_5kb_windows.bed', genome = config['genome']),
#         twentykb = expand('{genome}_20kb_windows.bed', genome = config['genome']),
#         hundredkb = expand('{genome}_100kb_windows.bed', genome = config['genome']),
        
#     run:
#         shell("bedtools makewindows -w 5000 -s 5000 -g {input} > {output.fivekb}"),
#         shell("bedtools makewindows -w 20000 -s 20000 -g {input} > {output.twentykb}"),
#         shell("bedtools makewindows -w 100000 -s 100000 -g {input} > {output.hundredkb}")
        
        
        
            
#

#bedtools makewindows -w 20000 -s 20000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 10kb_windows_hg38.bed

#bedtools makewindows -w 100000 -s 100000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 50kb_windows_hg38.bed
        
rule align:
    input:
        clip = expand(output_folder+'clip/{{desc}}_{{rep}}_{{EL}}_{reads}.clip', reads = config['reads']),
        txt = expand(output_folder+'clip/{{desc}}_{{rep}}_{{EL}}_{reads}.txt', reads = config['reads']),
        fa = expand(genome+'.fa'),
        fai = expand(genome+'.fa.fai'),
        amb = expand(genome+'.fa.amb'),
    output:
        output_folder+'bam/raw/{desc}_{rep}_{EL}.bam',
        #'samples/{desc}_{rep}_{EL}_{reads}.fastq.gz',
    threads: 8
    resources:
        tasks=1,
        cpus_per_task=8,
        mem_mb=128000,
        runtime=960,
        disk_mb = 128000
    run:
        
        
        command = "bwa mem -v 2 -t 8 "+genome+".fa '{input.clip}' > '{output}'"
        command = command[0]
        input_str = "'"
        for value in input.clip:
            input_str = input_str+value+"' '"
        input_str = input_str[0:-2]
        
        if len(config['reads']) > 1:
            #shell("bwa mem -v 2 -t 8 "+genome+" '{input.clip}' > '{output}.bam'")
            shell("bwa mem -v 2 -t 8 "+genome+".fa "+input_str+" > '{output}.bam'")
            shell("samtools fixmate -m '{output}.bam' '{output}'")
            shell("rm '{output}.bam'")
        else:
            shell("bwa mem -v 2 -t 8 "+genome+".fa '{input.clip}' > '{output}'")
            


#find *.clip | parallel --jobs 24 "bwa mem -v 2 -t 2 /panfs/roc/groups/0/riveramj/shared/bwaIndex_hg38/genome {} > {}_hg38.sam"

rule bamstats:
    input:
        output_folder+'bam/raw/{desc}_{rep}_{EL}.bam',
    output:
        output_folder+'bam/stats/raw/{desc}_{rep}_{EL}.bamstats',
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "samtools stats '{input}' > '{output}'"

        
rule filterbam:
    input:
        output_folder+'bam/raw/{desc}_{rep}_{EL}.bam',
    output:
        output_folder+'bam/filtered/{desc}_{rep}_{EL}.filtered',
        
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "samtools view -bhq 20 '{input}' -o '{output}'"

        
# rule sortbam:
#     input:
#         'bam/filtered/{desc}_{rep}_{EL}.filtered',
#     output:
#         allcontigs = 'bam/sorted/{desc}_{rep}_{EL}.sortedallcontigs.bam',
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     shell:
#         "samtools sort -m 10G {input} -o {output.allcontigs}"

        

# rule contigbam:
#     input:
#         bam='bam/sorted/{desc}_{rep}_{EL}.sortedallcontigs.bam',
#         index='bam/sorted/{desc}_{rep}_{EL}.sortedallcontigs.bam.bai'
#     output:        
#         'bam/sorted/{desc}_{rep}_{EL}.sorted',
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     run:        
#         shell('samtools view -bh {input.bam} $(cat contigs.txt | tr "\n" " ") > {output}.tmp'),
#         shell("samtools sort -m 10G {output}.tmp -o {output}"),
#         shell('rm {output}.tmp')


# rule bamindex:
#     input:
#         'bam/sorted/{desc}_{rep}_{EL}.sortedallcontigs.bam' 
#     output:
#         'bam/sorted/{desc}_{rep}_{EL}.sortedallcontigs.bam.bai'
        
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     shell:
#         'samtools index {input}'


rule contigbam:
    input:
        output_folder+'bam/filtered/{desc}_{rep}_{EL}.filtered',
    output:        
        output_folder+'bam/sorted/{desc}_{rep}_{EL}.sorted',
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    run: 
        contigs_list = config['contigs']
        with open('contigs.temp.txt', 'w') as f:
            f.write('\n'.join(contigs_list))
        shell("samtools sort -m 10G '{input}' -o '{output}.allcontigs.tmp'"),
        shell("samtools index '{output}.allcontigs.tmp'"),
        #shell('samtools view -bh {output}.allcontigs.tmp $(cat contigs.temp.txt | tr "\n" " ") > {output}.tmp'),
        shell("samtools view -bh '{output}.allcontigs.tmp' "+" ".join(contigs_list)+"> '{output}.tmp'"),
        shell("samtools sort -m 10G '{output}.tmp' -o '{output}'"),
        shell("samtools index '{output}'"),
        shell("rm '{output}.allcontigs.tmp'"),



#find *.bam | parallel --jobs 24 "samtools stats {} > {}.bamstats"

#find *.bam | parallel --jobs 24 "samtools view -bhq 20 {} -o {}.filtered"

#find *.filtered | parallel --jobs 24 "samtools sort -m 10G {} -o {}.sorted"

rule fbstats:
    input:
        output_folder+'bam/sorted/{desc}_{rep}_{EL}.sorted',
    output:
        output_folder+'bam/stats/fbstats/{desc}_{rep}_{EL}.fbstats',
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60

    shell:
        "samtools stats '{input}' > '{output}'"



rule rmdup:
    input:
        output_folder+'bam/sorted/{desc}_{rep}_{EL}.sorted',        
    output:
        output_folder+'bam/rmdup/{desc}_{rep}_{EL}_rmdup.bam',
        
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    run:
        shell("mkdir -p '"+output_folder+"bam/rmdup'"),
        shell("samtools markdup -r '{input}' '{output}'")

        
rule rmdup_stats:
    input:
        output_folder+'bam/rmdup/{desc}_{rep}_{EL}_rmdup.bam',        
    output:
        output_folder+'bam/stats/rmdup/{desc}_{rep}_{EL}.rdbamstats',
        
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    shell:
        "samtools stats '{input}' > '{output}'"

        
        
        
rule rpkm:
    input:
        bam = output_folder+'bam/rmdup/{desc}_{rep}_{EL}_rmdup.bam',
        windows = output_folder+'bed/{sizes}_windows_'+genomeid+'.bed',
    output:
        output_path=output_folder+'bg/rpkm/{desc}_{rep}_{EL}_{sizes}_'+genomeid+'_rpkm.bg',
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    run:
        import subprocess

        bam_path = input.bam
        #bam_path = "'"+str(bam_path)+"'"
        windows_path = input.windows
        #windows_path = "'"+str(windows_path)+"'"
        output_path = output.output_path
        #output_path = "'"+str(output_path)+"'"

        total_reads = int(
            subprocess.check_output(
                ["samtools", "view", "-c", bam_path],
                text=True
            ).strip()
        )

        if total_reads == 0:
            raise ValueError(f"No reads found in BAM: {bam_path}")

        scale = 1_000_000 / total_reads

        cmd = ["bedtools", "intersect", "-c", "-sorted", "-a", windows_path, "-b", bam_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

        with open(output_path, "w") as out:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue

                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                count = int(fields[3])

                window_len = end - start
                if window_len <= 0:
                    continue

                rpkm = count * 1000 * (scale / window_len)
                out.write(f"{chrom}\t{start}\t{end}\t{rpkm}\n")

        proc.stdout.close()
        return_code = proc.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(return_code, cmd)
        
rule normalize_smooth_rt:
    input:
        rt=output_folder+'bg/log2/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT.bg'
    output:
        rt=output_folder+'bg/norm/{desc}_{rep}_{sizes}_'+genomeid+'_RTnorm.bg'
    params:
        gap_file=lambda wc: config["rt_processing"].get("gap_file", ""),
        target_file=lambda wc: config["rt_processing"].get("target_file", ""),
        norm_method=lambda wc: config["rt_processing"].get("normalization", "quantile"),
        median_center=lambda wc: config["rt_processing"].get("median_center", True),
        exclude_chrY=lambda wc: config["rt_processing"].get("exclude_chrY", True),
        loess_span_bp=lambda wc: config["rt_processing"].get("loess_span_bp", 500000),
        mask_gaps=lambda wc: config["rt_processing"].get("mask_gaps", True),
        clip_min=lambda wc: config["rt_processing"].get("clip_min", -8),
        clip_max=lambda wc: config["rt_processing"].get("clip_max", 8),
        bin_size_map=lambda wc: config["rt_processing"].get("bin_size_map", {}),
        smooth=lambda wc: config["rt_processing"].get("smoothing", {}),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=120
    script:
        "scripts/normalize_smooth_rt.py"
        
rule log:
    input:
        early = output_folder+'bg/rpkm/{desc}_{rep}_E_{sizes}_'+genomeid+'_rpkm.bg',
        late = output_folder+'bg/rpkm/{desc}_{rep}_L_{sizes}_'+genomeid+'_rpkm.bg',
    output:
        output_folder+'bg/log2/{desc}_{rep}_{sizes}_'+genomeid+'_log2RT.bg',
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    run:
        shell(r"paste '{input.early}' '{input.late}' | awk '{{if($8 != 0 && $4 != 0){{print$1,$2,$3,log($4/$8)/log(2)}}}}' OFS='\t' > '{output}.tmp'")
        shell("cp '{output}.tmp' '{output}'")
        
        

        
#for file in *_E.bg; do
#	paste $file ${file%_E.bg}_L.bg | awk '{if($8 != 0 && $4 != 0){print$1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > ${file%_E.bg}_hg38_w50kb_RT.bg

#	done
        
#bedtools makewindows -w 5000 -s 5000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 10kb_windows_hg38.bed

#bedtools makewindows -w 20000 -s 20000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 10kb_windows_hg38.bed

#bedtools makewindows -w 100000 -s 100000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 50kb_windows_hg38.bed



#find *.sorted | parallel --jobs 24 "samtools stats {} > {}.fbstats"

#find *.sorted | parallel --jobs 24 "samtools markdup -r {} {}_rmdup.bam"

#find *_rmdup.bam | parallel --jobs 24 "samtools stats {} > {}.rdbamstats"

##module load bedtools

#bedtools makewindows -w 5000 -s 5000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 10kb_windows_hg38.bed

#bedtools makewindows -w 20000 -s 20000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 10kb_windows_hg38.bed

#bedtools makewindows -w 100000 -s 100000 -g /panfs/roc/groups/0/riveramj/shared/hg38.chrom.sizes > 50kb_windows_hg38.bed
