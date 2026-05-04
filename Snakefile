#sampledescription_Replicatenumber_EarlyorLate_Illumina-Lane_Illumina-read_filenumber.fastq.gz.
import os
import seaborn as sns
import pandas as pd
import matplotlib
from functools import reduce
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import pysam
configfile: "config.yaml"
target = config['target']
plot_scale = config['plot_scale']
target_list = []
import pyBigWig
for value in target:
    target_list.append(value)
import numpy as np
from scipy.stats import ttest_ind
from sklearn.preprocessing import StandardScaler
#print(target)
#print(expand("{targets}.clip", targets = target))
transcripts = pd.read_csv('transcripts.csv')

names = []
# for i in range(0,len(transcripts)):
    
#     try:
#         names.append(transcripts['info'][i].split('gene "')[1].split('";')[0])
#     except IndexError:
#         names.append(transcripts['info'][i].split('gene_id "')[1].split('";')[0])

#transcripts['name'] = names
bias_list = []
csv_bias_list = []
figure_bias_list = []
sliding_list = []
nucpos_list = []
tss_list = []
qc_list = []

#'ngmerge/{desc}_{strain}_{rep}_merged.fastq'


for file in target_list:
    temp = file.split('_R1.fastq.gz')[0].split('_R2.fastq.gz')[0]

    bias_list = bias_list + expand('bigwig/filtered/'+temp+'_rev_flip.bigwig')

    csv_bias_list = csv_bias_list + expand('csv/espan_bias/bias/'+temp+'.csv')
    figure_bias_list = figure_bias_list + expand('figures/espan_bias/png/'+temp+'.png')
    sliding_list = sliding_list + expand('figures/sliding/png/'+temp+'.png')
    nucpos_list = nucpos_list + expand('figures/nucpos/png/'+temp+'.png')
    tss_list = tss_list+expand('figures/tss/png/'+temp+'.png')
    qc_list = qc_list + expand('fastq/fastqc/'+temp+'_R1_fastqc.html')
    
#     acf_list = acf_list + expand('figures/acf/'+temp+'_{sizes}_{genomeid}_log2RT.svg', genomeid = config['genomeid'], sizes = config['window_sizes'])
#     temp = '_'.join(file[0:3])
#     stats_list.append('bam/stats/rmdup/'+temp+".rdbamstats")
#     stats_list.append('bam/stats/fbstats/'+temp+".fbstats")
#     figure_list = figure_list + expand('figures/read_count/all/raw/'+temp+'.csv', genomeid = config['genomeid'], sizes = config['window_sizes'])
#     rpkm_list = rpkm_list + expand('bg/rpkm/'+temp+'_5000_{genomeid}_rpkm.bg', genomeid = config['genomeid'], sizes = config['window_sizes'])
#     print(expand('bg/rpkm/'+temp+'_5000_{genomeid}_rpkm.bg', genomeid = config['genomeid'], sizes = config['window_sizes']))
    #output_list.append('samples/'+"_".join(file))
#output_list = list(set(output_list))


rule all:
    input:        
        bias_list,
        csv_bias_list,
        sliding_list,
        nucpos_list,
        tss_list,
        qc_list,
        figure_bias_list,
        
#         acf_list,
#         'figures/log2/log2RT_kdeplot.png',
#         'figures/read_count/reads.png',
#         'figures/windows/window_counts.png',
 
    
    

rule genome_index:
    input:
        expand('{genome}', genome = config['genome']),
    output:
        expand('{genome}.1.bt2', genome = config['genome']),
        fai=expand('{genome}.fai', genome = config['genome']),
        bwa_index = expand('{genome}.amb', genome = config['genome']),
        sizes = expand('{genome}.fai.sizes', genome = config['genome']),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=120
    run:
        #os.system('cd '+outputfolder+'bowtie/;'+'bowtie2-build '+genome_path+ ' ' +outputfolder+'bowtie/'+'genome.index')
        shell('bowtie2-build {input} {input}'),        
        shell("bwa index {input}"),
        shell("samtools faidx {input}"),
        shell("cut -f1,2 {input}.fai > {output.sizes}")    
    
    
    
rule fastqc:
    input:
        'samples/{number}_{desc}_R1.fastq.gz',
    output:
        html='fastq/fastqc/{number}_{desc}_R1_fastqc.html',
    threads: 1
    resources:
        mem_mb=16000,
        runtime=480
    run:
        if len(config['reads']) > 1:
            shell("fastqc "+str(input)+" "+str(input).split("R1")[0]+"R2.fastq.gz -t 2 --outdir=fastq/fastqc/")
        else:
            shell("fastqc "+str(input)+" -t 2 --outdir=fastq/fastqc/")

rule read_trimming:
    input:
        sample = 'samples/{number}_{desc}_{reads}.fastq.gz',
        
        #fq = expand('samples/{sample}.fastq.gz', sample = config["samples"]),
        #get_bwa_map_input_fastqs
    output:
        clip='fastq/trim/{number}_{desc}_{reads}.trim',
        txt='fastq/trim/{number}_{desc}_{reads}.txt',
        #clip = expand('clip/{sample}.clip', sample = config["samples"]),
        #'clip/{sample}.clip'
    threads: 4
    resources:
        mem_mb=64000,
        runtime=480
    run:
        r1 = expand('{r1_adapter}', r1_adapter = config['r1_adapter']),
        r2 = expand('{r2_adapter}', r2_adapter = config['r2_adapter']),
        print(len(config['reads']))
        print(input)
        if len(config['reads']) > 1:
            if '_R1.fastq' in str(input):
                print("r1")
                shell("cutadapt -a "+str(r1[0][0])+" -q 0 -O 1 -m 0 -j 0 -o {output.clip} {input.sample}")
                with open(str(output.txt), 'w') as f:
                     f.write("clip success")
                f.close()

            if '_R2.fastq' in str(input):
                print(r2)
                shell("cutadapt -a "+str(r2[0][0])+" -q 0 -O 1 -m 0 -j 0 -o {output.clip} {input.sample}")
                with open(str(output.txt), 'w') as f:
                     f.write("clip success")
                f.close()
        
        else:
            shell("cutadapt -a {input.r1} -q 0 -O 1 -m 0 -j 0 -o {output.clip} {input.sample}")
            with open(str(output.txt), 'w') as f:
                f.write("clip success")
            f.close()               
    
#rule read_trimming:
    
rule read_merge:
    input:
        #trim='clip/{desc}_{strain}_{rep}_{reads}.trim'
        ngmerge = expand('{ngmerge}', ngmerge = config['ngmerge']),
        trim = expand('fastq/trim/{{number}}_{{desc}}_{reads}.trim', reads = config['reads']),
        txt = expand('fastq/trim/{{number}}_{{desc}}_{reads}.txt', reads = config['reads']),
    output:
        expand('fastq/ngmerge/{{number}}_{{desc}}_merged.fastq', reads = config['reads']),
        
        
    threads: 8
    resources:
        mem_mb=16000,
        runtime=960,
        disk_mb = 16000,
           
    run:
        wd = str(config['working_dir'])
        if len(config['reads']) > 1:
            r1=str(input.trim[0])
            r2=str(input.trim[1])
            print(str(input.trim))
            out_r1=str(output[0])
            
            
            
            #shell('cd {input.ngmerge};')
            shell('{input.ngmerge}/NGmerge -1 '+r1+' -2 '+r2+' -o {output}')
            shell('touch {output}')
            
        elif config['centers']==False:
            shell('ln -s1 '+r1+' {output}')
            output_str= output.split('_R1_')[0]+'_R2_merged.fastq'
            shell('ln -s1 '+r2+' {output}')
        elif len(config['reads']) <= 1:
            shell('ln -s1 '+r1+' '+ '{output}')
        else:
            print("Error in Config File")
        
                

rule read_centers:
    input:
        'fastq/ngmerge/{number}_{desc}_merged.fastq',
    output:
        fastq='fastq/centers/{number}_{desc}_centers.fastq',
    threads: 8
    resources:
        mem_mb=128000,
        runtime=960,
        disk_mb = 128000
    run: 
        import time
        if config['centers'] == True:
            


           
            middle_seq = []
            try:
                with open(str(input), "rt") as handle:
                    for record in SeqIO.parse(handle, "fastq"):
                        if len(record.seq) > 149 and len(record.seq) < 170:
                            middle_seq.append(record[int(len(record.seq)/2)-37:int(len(record.seq)/2)+37])
            except UnicodeDecodeError:
                print('error')
                with gzip.open(str(input), "rt") as handle:
                    for record in SeqIO.parse(handle, "fastq"):
                        if len(record.seq) > 149 and len(record.seq) < 170:
                            middle_seq.append(record[int(len(record.seq)/2)-37:int(len(record.seq)/2)+37])

            print(middle_seq[0])
            with open('{output.fastq}.tmp', "w") as output_handle:
                SeqIO.write(middle_seq, output_handle, "fastq")
                
            print('true')
            
            while(os.path.isfile(output.fastq+".tmp")==False):
                start = time.time()
                time.sleep(60)
            end = time.time()
            print(end - start)
            time.sleep(60)
            shell('cp {output.fastq}.tmp {output}')
            shell('rm {output.fastq}.tmp')
        else:
            print('false')
            shell("cp  {input} {output.fastq}") 
        
        
            


            
 
            
                
rule align:
    input:
        reads = 'fastq/centers/{number}_{desc}_centers.fastq',
        
        #txt = expand('clip/{{desc}}_{{strain}}_{{rep}}_{reads}.txt', reads = config['reads']),
        fa = expand('{genome}', genome = config['genome']),
        #index = expand('{genome}', genome = config['genome']),
    output:
        fix = 'bam/raw/fixmate/{number}_{desc}.fixme',
        bam = 'bam/raw/bam/{number}_{desc}.bam',
        #'samples/{desc}_{rep}_{EL}_{reads}.fastq.gz',
    threads: 8
    resources:
        mem_mb=128000,
        runtime=960,
        disk_mb = 128000
    run:
        
        shell('bowtie2 -x {input.fa} -U {input.reads} -S {output.fix}')
        shell("samtools fixmate -m {output.fix} {output.bam}")
       
       
            
            
            
            
            
            
            
            
            
            
            
rule filterbam:
    input:
        'bam/raw/bam/{number}_{desc}.bam',
    output:
        'bam/filtered/{number}_{desc}.filtered.bam',
        
    threads: 1
    resources:
        mem_mb=8000,
        runtime=60
    shell:
        "samtools view -bhq 20 {input} -o {output}"            
            
            
            
            
rule sortbam:
    input:
        'bam/filtered/{number}_{desc}.filtered.bam',
    output:
        'bam/sorted/{number}_{desc}.sorted.bam',
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    shell:
        "samtools sort -m 10G {input} -o {output}"

        


rule bamindex:
    input:
        'bam/sorted/{number}_{desc}.sorted.bam' 
    output:
        'bam/sorted/{number}_{desc}.sorted.bam.bai'
        
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    shell:
        'samtools index {input}'     
        

rule coverage:
    input:
        bam = 'bam/sorted/{number}_{desc}.sorted.bam',
        index = 'bam/sorted/{number}_{desc}.sorted.bam.bai',
    output:
        'bigwig/raw/{number}_{desc}.bigwig'
    shell:   
        'bamCoverage -b {input.bam} -bs 1 -o {output}'
        
                

    

rule for_rev:
    input:
        bam = 'bam/sorted/{number}_{desc}.sorted.bam',
        index = 'bam/sorted/{number}_{desc}.sorted.bam.bai',
    output:
        forward='bam/filtered/{number}_{desc}_for.bam',
        rev='bam/filtered/{number}_{desc}_rev.bam',
        finished = 'bam/filtered/{number}_{desc}_rev.finished'
        
    run:
        bam=input.bam
        samfile = pysam.AlignmentFile(bam, 'rb')
       
        reads_for = pysam.Samfile(output.forward, "wb", template=samfile)
        reads_rev = pysam.Samfile(output.rev, "wb", template=samfile)
        for read in samfile.fetch():            
            if not read.is_reverse:
                reads_for.write(read)
            else:
                reads_rev.write(read)

        reads_for.close()
        reads_rev.close()
        samfile.close()
        pd.DataFrame().to_csv(output.finished)
        
rule filteredbamindex:
    input:
        forward='bam/filtered/{number}_{desc}_for.bam',
        rev='bam/filtered/{number}_{desc}_rev.bam',
        finished = 'bam/filtered/{number}_{desc}_rev.finished',
    output:
        'bam/filtered/{number}_{desc}_for.bam.bai',
        'bam/filtered/{number}_{desc}_rev.bam.bai',
        
    #threads: 1
    #resources:
    #    mem_mb=16000,
    #    runtime=60
    run:
        shell('samtools index {input.forward}')
        shell('samtools index {input.rev}')
    
rule filtered_coverage:
    input:
        bam1='bam/filtered/{number}_{desc}_for.bam',
        bam2='bam/filtered/{number}_{desc}_rev.bam',
        index1='bam/filtered/{number}_{desc}_for.bam.bai',
        index2='bam/filtered/{number}_{desc}_rev.bam.bai',
       
    output:
        bigwig1='bigwig/filtered/{number}_{desc}_for.bigwig',
        bigwig2='bigwig/filtered/{number}_{desc}_rev.bigwig',
    run:   
        shell('bamCoverage -b {input.bam1} -bs 1 -o {output.bigwig1}')
        shell('bamCoverage -b {input.bam2} -bs 1 -o {output.bigwig2}')
rule sizes:
    input:
        fa = expand('{genome}', genome = config['genome']),
        fai = expand('{genome}.fai', genome = config['genome'])
    output:
        "{input.fa}.sizes"
    shell:
        "cut -f1,2 {input.fai} > {input.fa}.sizes"
rule flip:
    input:
        bigwig='bigwig/filtered/{number}_{desc}_rev.bigwig',
        sizes = expand('{genome}.fai.sizes', genome = config['genome'])
    output:
        'bigwig/filtered/{number}_{desc}_rev_flip.bigwig'
    run:
        
        shell('bigWigToWig {input.bigwig} {input.bigwig}.wig')
        newlines=[]
        with open(input.bigwig+".wig") as f:
            lines = f.readlines()
        for line in lines:
            if len(line.split('\t')) == 4:
                splitline = line.split('\t')
                newline = splitline[0]+'\t'+splitline[1]+'\t'+splitline[2]+'\t'+str(int(splitline[3].split('\n')[0])*-1)+'\n'
                newlines.append(newline)
            else:
                newlines.append(line)
        f.close()        
        with open(input.bigwig+".wig", 'w') as f:
            for line in newlines:
                f.write(line)
        f.close()
        shell('wigToBigWig {input.bigwig}.wig {input.sizes} '+ "{output}")
        shell('rm {input.bigwig}.wig')
# rule sicer:
    
# rule nuchunter:



    
rule gene_cov:
    input:
        bigwig1='bigwig/filtered/{number}_{desc}_for.bigwig',
        bigwig2='bigwig/filtered/{number}_{desc}_rev.bigwig',
        bam = 'bam/sorted/{number}_{desc}.sorted.bam',
        index = 'bam/sorted/{number}_{desc}.sorted.bam.bai',
        transcripts = config['transcripts']
    output:
        png = 'figures/tss/png/{number}_{desc}.png',
        svg = 'figures/tss/svg/{number}_{desc}.svg',
        csv = 'figures/tss/csv/{number}_{desc}.csv',
        
    run:
        key= (input.index).split('/')[-1].split('.')[0]
        gene_coverage_dict = {}
        transcripts = pd.read_csv(input.transcripts)
        names = []
        for i in range(0,len(transcripts)):

            try:
                names.append(transcripts['info'][i].split('gene "')[1].split('";')[0])
            except IndexError:
                names.append(transcripts['info'][i].split('gene_id "')[1].split('";')[0])
        transcripts['name'] = names
        
        bw1 = pyBigWig.open(input.bigwig1)
        bw2 = pyBigWig.open(input.bigwig2)
        #total_cov = reduce(lambda x, y: x + y, [ int(l.rstrip('\n').split('\t')[0]) for l in pysam.idxstats(input.bam) ])
        total_cov = os.popen('samtools view -c -F 260 '+input.bam).read()
        total_cov = int(total_cov.split('\n')[0])
        upcover = []
        oncover = []
        temp_chr_coverage = pd.DataFrame()
        for j in range(0,len(transcripts)):

            if transcripts.iloc[j]['strand'] == '-':
                tss = transcripts.iloc[j]['end']
                try:
                    
                    
                    temp_chr_coverage[transcripts.iloc[j]['name']]=np.flip(np.array(bw1.values(transcripts.iloc[j]['chr'], int(tss-2000), int(tss+2000))+np.array(bw2.values(transcripts.iloc[j]['chr'], int(tss-2000), int(tss+2000)))))
                    
#                     upcover.append(np.mean(np.array(bw1.values(transcripts.iloc[j]['chr'], int(tss), int(tss+1000))+np.array(bw2.values(transcripts.iloc[j]['chr'], int(tss), int(tss+1000))))))
#                     oncover.append(np.mean(np.array(bw1.values(transcripts.iloc[j]['chr'], int(transcripts.iloc[j]['start']), int(transcripts.iloc[j]['end']))+np.array(bw2.values(transcripts.iloc[j]['chr'], int(transcripts.iloc[j]['start']), int(transcripts.iloc[j]['end']))))))
                except RuntimeError:
                    j=j
                    upcover.append(np.nan)
                    oncover.append(np.nan)
            else:
                tss = transcripts.iloc[j]['start']
                
                try:
                    temp_chr_coverage[transcripts.iloc[j]['name']]=(np.array(bw1.values(transcripts.iloc[j]['chr'], int(tss-2000), int(tss+2000)))+np.array(bw2.values(transcripts.iloc[j]['chr'], int(tss-2000), int(tss+2000))))

#                     upcover.append(np.mean(np.array(bw1.values(transcripts.iloc[j]['chr'], int(tss-1000), int(tss))+np.array(bw2.values(transcripts.iloc[j]['chr'], int(tss-1000), int(tss))))))
#                     oncover.append(np.mean(np.array(bw1.values(transcripts.iloc[j]['chr'], int(transcripts.iloc[j]['start']), int(transcripts.iloc[j]['end']))+np.array(bw2.values(transcripts.iloc[j]['chr'], int(transcripts.iloc[j]['start']), int(transcripts.iloc[j]['end']))))))

                except RuntimeError:
                    upcover.append(np.nan)
                    oncover.append(np.nan)
                    j=j
            
            #transcripts['up_'+key] = np.array(upcover)/total_cov*1000000
            #transcripts_orig['up_'+key] = np.array(upcover)/total_cov*1000000

            #transcripts['on_'+key] = np.array(oncover)/total_cov*1000000
            #transcripts_orig['on_'+key] = np.array(oncover)/total_cov*1000000

            #transcripts_orig.to_csv(output.csv)
        avg_cov = np.mean(temp_chr_coverage)
        temp_chr_coverage = temp_chr_coverage.T
        temp_chr_coverage['avg_cov'] = avg_cov
        temp_chr_coverage = temp_chr_coverage.sort_values(by='avg_cov', ascending = False)
        temp_chr_coverage = temp_chr_coverage.T
        scaler = StandardScaler()
        temp_chr_coverage = pd.DataFrame(scaler.fit_transform(temp_chr_coverage), index = temp_chr_coverage.index, columns = temp_chr_coverage.columns)


            



        plotdf = pd.DataFrame()
        
        
        plotdf[str(key)] = np.mean(temp_chr_coverage.iloc[0:4000].T)
        print(str(key))
        print(str(list(plotdf.index)))
        #x = plotdf.index,
        sns.lineplot(plotdf).set(title=key)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        plt.legend([],[], frameon=False)

        plt.savefig(output.png, bbox_inches="tight")
        plt.savefig(output.svg, bbox_inches="tight")
        plotdf.to_csv(output.csv)
 
 



rule nucpos:
    input:
        bigwig1='bigwig/filtered/{number}_{desc}_for.bigwig',
        bigwig2='bigwig/filtered/{number}_{desc}_rev.bigwig',
        bam = 'bam/sorted/{number}_{desc}.sorted.bam',
        index = 'bam/sorted/{number}_{desc}.sorted.bam.bai',
        originfile = config['origins'],
        #expand('{genome}.index', genome = config['genome'])
    output:
        png = 'figures/nucpos/png/{number}_{desc}.png',
        svg = 'figures/nucpos/svg/{number}_{desc}.svg',
        csv = 'figures/nucpos/csv/{number}_{desc}.csv',
    run:
        key= (input.index).split('/')[-1].split('.')[0]
        import matplotlib.pyplot as plt
        ori_g1 = pd.read_csv(input.originfile, delim_whitespace=True, header = None)
        position_plotdf=pd.DataFrame()
        origins_coverage_dict = {}
        
        forfile = input.bigwig1
        revfile = input.bigwig2
        bw1 = pyBigWig.open(forfile)
        bw2 = pyBigWig.open(revfile)
        temp_chr_coverage = pd.DataFrame()
        for j in range(0,len(ori_g1)):           
            temp_chr_coverage[ori_g1.iloc[j][3]]=(np.array(bw1.values(ori_g1.iloc[j][0], ori_g1.iloc[j][1]-1000, ori_g1.iloc[j][1]+1000))+np.array(bw2.values(ori_g1.iloc[j][0], ori_g1.iloc[j][1]-1000, ori_g1.iloc[j][1]+1000)))

        avg_cov = np.mean(temp_chr_coverage)
        temp_chr_coverage = temp_chr_coverage.T
        temp_chr_coverage['avg_cov'] = avg_cov
        temp_chr_coverage = temp_chr_coverage.sort_values(by='avg_cov', ascending = False)
        temp_chr_coverage = temp_chr_coverage.T
        scaler = StandardScaler()
        temp_chr_coverage = pd.DataFrame(scaler.fit_transform(temp_chr_coverage), index = temp_chr_coverage.index, columns = temp_chr_coverage.columns)


        origins_coverage_dict[key] = temp_chr_coverage
        


        plotdf = pd.DataFrame()
        positions = []
        for key in list(origins_coverage_dict.keys()):
            plotdf[str(key)] = np.mean(origins_coverage_dict[key].iloc[0:2000].T.iloc[0:50])

            sns.lineplot(plotdf).set(title=key)
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
            plt.legend([],[], frameon=False)

            plt.savefig(output.png, bbox_inches="tight")
            plt.savefig(output.svg, bbox_inches="tight")
            plotdf.to_csv(output.csv)
            positions = positions + list(range(0,2000))
            plt.show()

            
            
            
rule sliding:
    input:
        bigwig1='bigwig/filtered/{number}_{desc}_for.bigwig',
        bigwig2='bigwig/filtered/{number}_{desc}_rev.bigwig',
        bam = 'bam/sorted/{number}_{desc}.sorted.bam',
        index = 'bam/sorted/{number}_{desc}.sorted.bam.bai',
        originfile =  config['origins'],
    output:
        png = 'figures/sliding/png/{number}_{desc}.png',
        svg = 'figures/sliding/svg/{number}_{desc}.svg',
        csv = 'figures/sliding/csv/{number}_{desc}.csv',
    run:
        key= (input.index).split('/')[-1].split('.')[0]
        import matplotlib.pyplot as plt
        ori_g1 = pd.read_csv(input.originfile, delim_whitespace=True, header = None)
        position_plotdf=pd.DataFrame()
        origins_coverage_dict = {}

        forfile = input.bigwig1
        revfile = input.bigwig2
        bw1 = pyBigWig.open(forfile)
        bw2 = pyBigWig.open(revfile)
        temp_chr_coverage_for = pd.DataFrame()
        temp_chr_coverage_rev = pd.DataFrame()
        temp_chr_coverage = pd.DataFrame()

        for j in range(0,len(ori_g1)):           
            try:
                temp_chr_coverage_for[ori_g1.iloc[j][3]]=(np.array(bw1.values(ori_g1.iloc[j][0], ori_g1.iloc[j][1]-10100, ori_g1.iloc[j][1]+10100)))


                temp_chr_coverage_rev[ori_g1.iloc[j][3]]=np.array(bw2.values(ori_g1.iloc[j][0], ori_g1.iloc[j][1]-10100, ori_g1.iloc[j][1]+10100))

                temp_chr_coverage[ori_g1.iloc[j][3]]=(np.array(bw1.values(ori_g1.iloc[j][0], ori_g1.iloc[j][1]-10100, ori_g1.iloc[j][1]+10100))+
                                                     np.array(bw2.values(ori_g1.iloc[j][0], ori_g1.iloc[j][1]-10100, ori_g1.iloc[j][1]+10100)))
            except:
                print('end of chr')
        #print(temp_chr_coverage)
        avg_cov = np.mean(temp_chr_coverage)
        temp_chr_coverage_for = temp_chr_coverage_for.T
        temp_chr_coverage_for['avg_cov'] = avg_cov
        temp_chr_coverage_for = temp_chr_coverage_for.sort_values(by='avg_cov', ascending = False)
        temp_chr_coverage_for = temp_chr_coverage_for.T
        #scaler = StandardScaler()
        #temp_chr_coverage_for = pd.DataFrame(scaler.fit_transform(temp_chr_coverage), index = temp_chr_coverage.index, columns = temp_chr_coverage.columns)

        avg_cov = np.mean(temp_chr_coverage)
        temp_chr_coverage_rev = temp_chr_coverage_rev.T
        temp_chr_coverage_rev['avg_cov'] = avg_cov
        temp_chr_coverage_rev = temp_chr_coverage_rev.sort_values(by='avg_cov', ascending = False)
        temp_chr_coverage_rev = temp_chr_coverage_rev.T


        temp_for = pd.DataFrame()
        for j in range(100,len(temp_chr_coverage_for)-100,20):
            temp_for[j]=(np.nanmean(temp_chr_coverage_for[j-100:j+100]+0.1, axis=0))

        temp_rev = pd.DataFrame()
        for j in range(100,len(temp_chr_coverage_rev)-100,20):
            temp_rev[j]=(np.nanmean(temp_chr_coverage_rev[j-100:j+100]+0.1, axis=0))


        temp = pd.DataFrame()

        for j in range(100,len(temp_chr_coverage_rev)-100,20):
            temp[j]=np.log2(temp_for[j]/temp_rev[j]).replace([np.inf, -np.inf], np.nan)
        origins_coverage_dict = {}
        origins_coverage_dict[key] = temp.T
        plotdf = pd.DataFrame()
        #print(key)
        plotdf[str(key)] = pd.DataFrame(np.nanmean(origins_coverage_dict[key].iloc[0:20000].T, axis=0)).rolling(10).mean()
        #plotdf[str(key)] = pd.DataFrame(np.array(origins_coverage_dict[key]))
        #plotdf=plotdf.drop('avg_cov') 
        fig, ax = plt.subplots()
        sns.lineplot(data=plotdf, x = plotdf.index, y = key, linewidth=5).set(title=key)
        plt.vlines(x=500, colors='black', ls='--', ymin=-1000, ymax=1000)
        plt.hlines(y=0, colors='black', ls='--', xmin=-1000, xmax=1000)
        ax.set_ylim(-3,3)
        ax.set_xlim(0,1000)
        ax.set_xticklabels([-10000, -6000, -2000, 2000, 6000, 10000])
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        plt.legend([],[], frameon=False)
        #plt.figure(figsize=(8,16))
        plt.savefig(output.png)#, bbox_inches="tight")
        plt.savefig(output.svg)#, bbox_inches="tight")
        plotdf.to_csv(output.csv)
    
    
    
rule nucleosome_coverage:
    input:
        bigwig1='bigwig/filtered/{number}_{desc}_for.bigwig',
        bigwig2='bigwig/filtered/{number}_{desc}_rev.bigwig',
        
        nucleosomefile = config['nucleosomes']
    output:
        bed = 'bed/nucleosome_bias/{number}_{desc}.bed',
    
    
    
    run:
        #test new nucleosomes method

        nucleosome_dict = {}

        forfile = input.bigwig1
        revfile = input.bigwig2
        #nucbed=str(inputdf.iloc[i][0])


        #nucleosome coverage
        nucleosomes_bed= input.nucleosomefile
        nucleosomes_bed = pd.read_csv(nucleosomes_bed, sep = '\t', header = None, names = ['chr', 'pos1', 'pos2', 'name', 'pos3', 'strand'])
        #nucleosomes_bed = nucleosomes_bed.sort_values(by = 'chr')
        #nucleosomes_bed.header=
        #file = "/home/cyu/Jupyter/Espan_pipeline/output/bigwig/K4_s30_2nd_999_CKDL220015030-1A_HY7MKDSX3_L2_1_merged_middle_sort_filtered_for.bigwig"
        bw = pyBigWig.open(forfile)
        for_coverage = []
        for chrom in bw.chroms().keys():
            nucleosomes_chr = nucleosomes_bed[nucleosomes_bed['chr'] == chrom]
            for ind in nucleosomes_chr.index:
                start = nucleosomes_chr['pos1'][ind]-73
                end = nucleosomes_chr['pos2'][ind]+73
                try:
                    for_coverage.append(np.mean(bw.values(chrom, start, end)))#+1)
                except RuntimeError:
                    try:
                        for_coverage.append(np.mean(bw.values(chrom, 0, end)))#+1)
                    except RuntimeError:
                        for_coverage.append(np.mean(bw.values(chrom, start, bw.chroms()[chrom])))#+1)

        bw = pyBigWig.open(revfile)
        rev_coverage = []                
        for chrom in bw.chroms().keys():
            nucleosomes_chr = nucleosomes_bed[nucleosomes_bed['chr'] == chrom]
            for ind in nucleosomes_chr.index:
                start = nucleosomes_chr['pos1'][ind]-73
                end = nucleosomes_chr['pos2'][ind]+73
                try:
                    rev_coverage.append(np.mean(bw.values(chrom, start, end)))#+1)
                except RuntimeError:
                    try:
                        rev_coverage.append(np.mean(bw.values(chrom, 0, end)))#+1)
                    except RuntimeError:
                        rev_coverage.append(np.mean(bw.values(chrom, start, bw.chroms()[chrom])))#+1)

        nucleosomes_bed['for_coverage'] = for_coverage
        nucleosomes_bed['rev_coverage'] = rev_coverage

        for_coverage = np.array(for_coverage)
        for_coverage = for_coverage.astype('float')
        for_coverage[for_coverage == 0] = 'nan'            

        rev_coverage = np.array(rev_coverage)
        rev_coverage = rev_coverage.astype('float')
        rev_coverage[rev_coverage == 0] = 'nan'


        nucleosomes_bed['log2ratio'] = np.log2(np.array(for_coverage)/np.array(rev_coverage))
        averagecoverage= []
        for cov_index in range(0,len(for_coverage)):
            averagecoverage.append(np.mean([for_coverage[cov_index], rev_coverage[cov_index]]))

        nucleosomes_bed['average'] = averagecoverage
        nucleosomes_bed['for_coverage'] = for_coverage
        nucleosomes_bed['rev_coverage'] = rev_coverage
        #nucleosome_dict[inputdf[sampletype][i]] = nucleosomes_bed

        nucleosomes_bed.to_csv(output.bed, sep = '\t')
 
    
rule bias_dict:
    input:
        bigwig1='bigwig/filtered/{number}_{desc}_for.bigwig',
        bigwig2='bigwig/filtered/{number}_{desc}_rev.bigwig',
        bam = 'bam/sorted/{number}_{desc}.sorted.bam',
        index = 'bam/sorted/{number}_{desc}.sorted.bam.bai',
        originfile =  config['origins'],
        bed = 'bed/nucleosome_bias/{number}_{desc}.bed'
    output:
        bias = 'csv/espan_bias/bias/{number}_{desc}.csv',
        nuc = 'csv/espan_bias/nuc_coverage/all/{number}_{desc}.csv',
        forw = 'csv/espan_bias/nuc_coverage/for/{number}_{desc}.csv',
        rev = 'csv/espan_bias/nuc_coverage/rev/{number}_{desc}.csv',
    run:  
        #use nucleosome bias dict to get bias of each origin of replication
        do_brdu = True
        bw = pyBigWig.open(input.bigwig1)
        bias_dict = {}
        nuc_cover_dict = {}
        for_cover_dict = {}
        rev_cover_dict = {}
        defined_origins = True
        nucleosome_bed = pd.read_csv(input.bed, sep='\t', header=0)
        bias_df = pd.DataFrame()
        nuc_cover_df = pd.DataFrame()
        for_cover_df = pd.DataFrame()
        rev_cover_df = pd.DataFrame()
        #key = inputdf.iloc[input_index][1]

        ori_g1 = pd.read_csv(input.originfile, delim_whitespace=True, header = None)
        origin_added = 0
        origin_failed = 0

        for chromkey in bw.chroms().keys():

            for ori_index in range(0,len(ori_g1[ori_g1[0] == chromkey])):
                ori_name = ori_g1[ori_g1[0] == chromkey].iloc[ori_index][3]
                ori_pos = ori_g1[ori_g1[0] == chromkey].iloc[ori_index][4]

                upstream_nucleosomes = (nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos2']<ori_pos]['log2ratio']) 
                upstream_nucleosomes_coverage = nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos2']<ori_pos]['average']
                upstream_for_coverage = nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos2']<ori_pos]['for_coverage']
                upstream_rev_coverage = nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos2']<ori_pos]['rev_coverage']


                temp_nuc = list(upstream_nucleosomes.fillna(value=1).iloc[-11:-1])
                temp_nuc_cov = list(upstream_nucleosomes_coverage.fillna(value=1).iloc[-11:-1])
                temp_nuc_for = list(upstream_for_coverage.fillna(value=1).iloc[-11:-1])
                temp_nuc_rev = list(upstream_rev_coverage.fillna(value=1).iloc[-11:-1])




                downstream_nucleosomes = (nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos1']>ori_pos]['log2ratio'])
                downstream_nucleosomes_coverage = (nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos1']>ori_pos]['average']) 
                downstream_for_coverage = (nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos1']>ori_pos]['for_coverage']) 
                downstream_rev_coverage = (nucleosome_bed[nucleosome_bed['chr'] == chromkey][nucleosome_bed['pos1']>ori_pos]['rev_coverage']) 




                temp_nuc = temp_nuc+ list(downstream_nucleosomes.fillna(value=1).iloc[0:10])
                temp_nuc_cov = temp_nuc_cov+ list(downstream_nucleosomes_coverage.fillna(value=1).iloc[0:10])
                temp_nuc_for = temp_nuc_for+ list(downstream_for_coverage.fillna(value=1).iloc[0:10])
                temp_nuc_rev = temp_nuc_rev+ list(downstream_rev_coverage.fillna(value=1).iloc[0:10])

                

                if len(temp_nuc) == 20 and len(temp_nuc_cov) == 20:
                    print("origin_added")
                    origin_added+=1
                    temp_nuc.append(np.mean(temp_nuc_cov))
                    bias_df[ori_name] = temp_nuc
                    nuc_cover_df[ori_name] = temp_nuc_cov
                    for_cover_df[ori_name] = temp_nuc_for
                    rev_cover_df[ori_name] = temp_nuc_rev
                                   
                else:
                    origin_failed+=1
                    print("len="+str(len(temp_nuc))+" "+str(len(temp_nuc_cov)))
                    
                
                   


            print(origin_added)
            print(origin_failed)
            bias_df.T.to_csv(output.bias)
            nuc_cover_df.T.to_csv(output.nuc)
            for_cover_df.T.to_csv(output.forw)
            rev_cover_df.T.to_csv(output.rev)

rule plot_bias:
    input:
        bias = 'csv/espan_bias/bias/{number}_{desc}.csv',
        nuc = 'csv/espan_bias/nuc_coverage/all/{number}_{desc}.csv',
        forw = 'csv/espan_bias/nuc_coverage/for/{number}_{desc}.csv',
        rev = 'csv/espan_bias/nuc_coverage/rev/{number}_{desc}.csv',
    output:
        png = 'figures/espan_bias/png/{number}_{desc}.png',
        svg = 'figures/espan_bias/svg/{number}_{desc}.svg',
        #csv = 'figures/espan_bias/csv/{number}-{desc}.csv',
        
    run:
        count = 0
        key= (input.bias).split('/')[-1].split('.')[0]
        #plt=reload(plt)
        
            
    #         print(count)

    #         fig, ((ax3, ax1), (ax4, ax2), (ax6,ax5)) = plt.subplots(nrows=3, ncols=2,  
    #                                                 gridspec_kw={'height_ratios': [5, 3, 3], 
    #                                                               'width_ratios': [1, 20]})
            
        fig, ((ax3, ax1), (ax4, ax2)) = plt.subplots(nrows=2, ncols=2,  
                                                        gridspec_kw={'height_ratios': [5, 3], 
                                                                      'width_ratios': [1, 20]})


        bias_df = pd.read_csv(input.bias, header=0, index_col=0)
        for_df = pd.read_csv(input.forw, header=0, index_col=0)
        rev_df = pd.read_csv(input.rev, header=0, index_col=0)

        combined_df = pd.concat([bias_df, for_df, rev_df])
        combined_df = combined_df[combined_df['20']>1]
        combined_df = combined_df.sort_values(by = '20', ascending=False)
        print(len(combined_df.columns))
        #biasdf = combined_df.sort_values(by = 20, ascending = False).T[0:20].T.iloc[0:50]
        #fordf = combined_df.sort_values(by = 20, ascending = False).T[21:41].T.iloc[0:50]
        #revdf = combined_df.sort_values(by = 20, ascending = False).T[41:61].T.iloc[0:50]
        #avg_cover = np.mean(bias_dict[key].sort_values(by = 20, ascending = False).T[20:21].T).iloc[0:50]
        biasdf = combined_df.T[0:20].T
        fordf = combined_df.T[21:41].T
        revdf = combined_df.T[41:61].T
        avg_cover = np.mean(bias_df.T[20:21].T, axis=0)
        lead_cover_avg = np.mean([np.mean(fordf.T.iloc[10:21].T.melt()['value']),
                                  np.mean(revdf.T.iloc[0:10].T.melt()['value'])])
        lag_cover_avg = np.mean([np.mean(fordf.T.iloc[0:10].T.melt()['value']), 
                                 np.mean(revdf.T.iloc[10:21].T.melt()['value'])])
        plotdf = biasdf.T[0:20].T.melt()
        color = []
        means = []
        for i in range(0,len(biasdf.columns)):
            means.append(np.mean(biasdf))


        up = np.array(list(biasdf.iloc[:,0:10].melt()['value'])).astype(float)
        
        print(len(up))
        down =  np.array(list(biasdf.iloc[:,10:20].melt()['value'])).astype(float)
        print(len(down))
        tstat, p = ttest_ind(up, down)

        p='{:.2e}'.format(p)

        avg_bias = np.mean([-1*np.mean(up),np.mean(down)])
        avg_bias='{:.2}'.format(avg_bias)

        from matplotlib.colors import LinearSegmentedColormap
        colormap = LinearSegmentedColormap.from_list('rg', ['r', 'k', 'lime'], N=256)
        print(list(biasdf.columns))
        sns.heatmap(biasdf, ax=ax2, cmap = 'RdBu',center=0, cbar_ax=ax4, vmin = -1, vmax=1)
        ax3.axis('off')
        #ax6.axis('off')
        ax2.set_yticks([])
        ax2.set_xticklabels(['-10', '-9', '-8', '-7' ,'-6', '-5', '-4', '-3', '-2', '-1', 
                             '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
        # Get the images on an axis
        im = ax4.collections

        sns.barplot(data = plotdf, x = 'variable', y = 'value', ax = ax1, color = 'grey', ci=None)#, palette = colormeans)
        ax1.set_xticks([])
#         print(lead_cover_avg)
#         print(lag_cover_avg)
        ax1.xaxis.set_visible(False)
        ax4.yaxis.set_ticks_position('left')


        ax1.set_ylim(-1*float(plot_scale), float(plot_scale))

        #ax4.get_children()[-13].set_visible(False)


        #plt.title(key)
        labels = list(range(-10,11,1))
        labels.pop(10)
        #ax5.set_xticklabels(labels)
        plt.tight_layout()
        plt.xlabel("nucleosome position")
        ax1.set(title=key+'\n' +'log2(lead/lag): ' + str(avg_bias)+'\n'+"p-value:" +p + '\n'+"avg cov:" + str(np.round(float(avg_cover), 2)))
        ax1.set(ylabel='log2(Watson/Crick)')
        plt.subplots_adjust(wspace=0, hspace=0.05)
        plt.savefig(output.png, bbox_inches="tight")
        plt.savefig(output.svg, bbox_inches="tight")
        #plt.show()
        #biasdf.T[0:20].T.to_csv(output.csv)
                                         
                                          
            

# rule origin_cov:
    
# rule bias:
    

        
# rule raw_and_acfplot:
#     resources:        
#         threads=1,
#         mem_mb=4000,
#         runtime=60
#     input:
#         'bg/log2/{desc}_{rep}_{sizes}_{genomeid}_log2RT.bg',
#     output:
#         svg = 'figures/acf/{desc}_{rep}_{sizes}_{genomeid}_log2RT.svg',
#         png = 'figures/acf/{desc}_{rep}_{sizes}_{genomeid}_log2RT.png',
#         svg_raw = 'figures/raw/10Mb/{desc}_{rep}_{sizes}_{genomeid}_log2RT_10Mb_raw.svg',
#         png_raw = 'figures/raw/10Mb/{desc}_{rep}_{sizes}_{genomeid}_log2RT_10Mb_raw.png',
#         svg_raw_2Mb = 'figures/raw/2Mb/{desc}_{rep}_{sizes}_{genomeid}_log2RT_2Mb_raw.svg',
#         png_raw_2Mb = 'figures/raw/2Mb/{desc}_{rep}_{sizes}_{genomeid}_log2RT_2Mb_raw.png',
    
    
#     run:
#         df = pd.read_csv(str(input), sep = '\t', header = None)
#         df.columns = ['chr', 'start', 'end', 'RT']
#         ax = pd.plotting.autocorrelation_plot(df['RT'])
#         lag2 = ax.lines[5].get_ydata()[1]
#         plt.axhline(y = lag2, color = 'r', linestyle = '-')
#         ax.set_xlim([0, 1000])
#         ax.set_ylim([0, 1])
#         plt.text(650, 0.81, "ACF = " + str(round(lag2, 2)), size='x-large')
#         plt.savefig(output.png)
#         plt.savefig(output.svg)
#         plt.close()
#         chromo="NC_000001.11"
#         df = df.loc[df['chr'] == chromo]
#         start = 5*1000000
#         end = 105*1000000
#         df = df.loc[df['start'] >= start]
#         df = df.loc[df['end'] <= end]
#         end = 25*1000000
#         df2 = df.loc[df['end'] <= end]
#         df['start'] = df['start']/1000000
#         df['end'] = df['end']/1000000
#         df2['start'] = df2['start']/1000000
#         df2['end'] = df2['end']/1000000
#         ax = sns.scatterplot(data = df, x = "start", y = 'RT', s=3, linewidth=0)
#         ax.set(xlabel='NC_000001.11' + " (Mb)")
#         plt.ylim(-7.5, 7.5)
#         plt.axhline(y = 0, color = 'b', linestyle = '-')
#         plt.savefig(output.png_raw)
#         plt.savefig(output.svg_raw)
#         plt.close()
        
#         ax = sns.scatterplot(data = df2, x = "start", y = 'RT', s=3, linewidth=0)
#         ax.set(xlabel='NC_000001.11' + " (Mb)")
#         plt.ylim(-7.5, 7.5)
#         plt.axhline(y = 0, color = 'b', linestyle = '-')
#         plt.savefig(output.png_raw_2Mb)
#         plt.savefig(output.svg_raw_2Mb)
    
        
# rule combine_figures:
#     input:
#         figure_list
#     output:
#         svg = 'figures/read_count/reads.svg',
#         png = 'figures/read_count/reads.png',
#     threads: 1
#     resources:
#         mem_mb=4000,
#         runtime=60
        
#     run:
#         combineddf = pd.DataFrame()
#         for figure in input:
#             combineddf = pd.concat([combineddf, pd.read_csv(figure)])

#         names = list(combineddf['Unnamed: 0'])
#         desc = []
#         rep = []
#         earlylate = []
#         for value in names:
#             temp = value.split(r'/')[-1].split('_')
#             desc.append(temp[0])
#             rep.append(temp[1])
#             earlylate.append(temp[2].split('.')[0])
#         combineddf['desc'] = desc
#         combineddf['rep'] = rep
#         combineddf['el'] = earlylate
#         plotnumber = len(set(desc))*2
#         descs = list(set(desc))
#         fig,ax=plt.subplots(nrows=plotnumber, figsize=(5,plotnumber*5))
#         #subfigs = fig.subfigures(1, len(descs))
#         i = 0
#         while i <len(set(desc))*2:
#             plotdf = combineddf[combineddf['desc']==descs[int(i/2)]]
#             plotdf = plotdf.drop(labels=['Unnamed: 0', 'desc'], axis = 1)
#             melted = plotdf.melt(id_vars = ['el', 'rep'])
#             melted1 = melted[melted['el'] == 'E']
#             barplot = sns.barplot(data = melted1, ax = fig.axes[i], y = 'variable', x = 'value', hue = 'rep')
#             barplot.legend_.remove()
#             barplot.set(xlabel=None)
#             #barplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)
#             #barplot.xticks(rotation=45)
#             fig.axes[i].set_title(descs[int(i/2)]+" - Early")
#             fig.axes[i].set_yticklabels(['Raw', 'Mapped', 'Filtered', 'Deduplicated'], rotation=0)
#             melted2 = melted[melted['el'] == 'L']
#             barplot = sns.barplot(data = melted2, ax = fig.axes[i+1], y = 'variable', x = 'value', hue = 'rep')
#             barplot.set(yticklabels=[])
#             barplot.set(ylabel=None)
#             barplot.set(xlabel=None)
#             barplot.legend_.remove()
#             #barplot.legend(loc='center left', bbox_to_anchor=(1.25, 0.5), ncol=1)
#             fig.axes[i+1].set_yticklabels(['Raw', 'Mapped', 'Filtered', 'Deduplicated'], rotation=0)
#             fig.axes[i+1].set_title(descs[int(i/2)]+" - Late")
#             #barplot.xticks(rotation=45)
    
    
#             i+=2

        
#         #for i in range(0,len(subfigs)):
#         #    subfigs[i].suptitle(descs[i])
#         plt.tight_layout()
        
#         plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
#         plt.savefig(output.svg, bbox_inches='tight',pad_inches=0.1)
#         plt.close()

        
# rule plotRT:
#     input:
#         output_list
#     output:
#         svg = 'figures/log2/log2RT_kdeplot.svg',
#         png = 'figures/log2/log2RT_kdeplot.png',
#     threads: 1
#     resources:
#         mem_mb=4000,
#         runtime=60
#     run:
#         df = pd.DataFrame()
#         for value in output_list:
#             temp_df = pd.DataFrame()
#             if "5000" in value:
#                 temp_df = pd.read_csv(value, sep = '\t', header = None)
#                 temp_df.columns = ['chr', 'start', 'end', 'RT']
#                 temp_df['name'] = value.split(r'/')[-1].split('.bg')[0]
#             if len(df) == 0:
#                 df = temp_df
#             else:
#                 df = pd.concat([df, temp_df], axis=0)
#         ax = sns.kdeplot(data=df, x="RT", hue = 'name')
#         sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.4))
#         plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
#         plt.savefig(output.svg, bbox_inches='tight', pad_inches=0.1)
#         plt.close()

        
# rule plot_read_counts:
#     input:
#         clip = 'clip/{desc}_{rep}_{EL}_R1.clip',
#         bam = 'bam/raw/{desc}_{rep}_{EL}.bam',
#         filtered = 'bam/filtered/{desc}_{rep}_{EL}.filtered',
#         rmdup = 'bam/rmdup/{desc}_{rep}_{EL}_rmdup.bam',
#     output:
#         png = 'figures/read_count/all/{desc}_{rep}_{EL}.png',
#         svg = 'figures/read_count/all/{desc}_{rep}_{EL}.svg',
#         csv = 'figures/read_count/all/raw/{desc}_{rep}_{EL}.csv'
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
    
        
#     run:
        
#         raw_count = int(os.popen('echo $(cat '+input.clip+'|wc -l)/4|bc').read())
#         bam_count = int(os.popen('samtools view -c '+input.bam).read())
#         filtered_count = int(os.popen('samtools view -c '+input.filtered).read())/2        
#         rmdup_count = int(os.popen('samtools view -c '+input.rmdup).read())/2        
#         if len(config['reads']) > 1:
#             bam_count = bam_count/2
#             filtered_count = filtered_count/2
#             rmdup_count = rmdup_count/2
#         ylist = [raw_count/raw_count, bam_count/raw_count, filtered_count/raw_count, rmdup_count/raw_count]
#         xlist = ['Raw:\n'+str(raw_count), 'Mapped:\n'+str(bam_count), 'Filtered:\n'+str(filtered_count), 'Deduplicated:\n'+str(rmdup_count)]
#         fig = sns.barplot(x = xlist, y = ylist)
#         plt.savefig(output.png)
#         plt.savefig(output.svg)
#         plt.close()
#         df = pd.DataFrame()
#         df[input.clip.split('.clip')[0]] = ylist
#         df = df.T
#         df.columns = xlist
#         df.to_csv(output.csv)
    
# rule plot_window_counts:
#     input:
#         rpkm_list
#     output:
#         png = 'figures/windows/window_counts.png',
#         svg = 'figures/windows/window_counts.svg',
#     threads: 1
#     resources:
#         mem_mb=4000,
#         runtime=60
    
#     run:
#         plotdf = pd.DataFrame()
#         covered = []
#         none = []
#         plotdf_index = []
#         file_list = list(input)
#         for file in file_list:
#             df = pd.read_csv(str(file),  sep = '\t', header = None)
#             df.columns = ['chr', 'start', 'end', 'RT']
#             df['label'] = '1'
#             df['covered'] = [0 if x <= 0 else 1 for x in df['RT']]
#             non = df['covered'].value_counts()[0]
#             cov = df['covered'].value_counts()[1]
#             tot = non + cov
#             non = non/(tot)
#             cov = cov/tot        
#             covered.append(cov)
#             none.append(non)
#             plotdf_index.append('_'.join(str(file).split(r'/')[-1].split('_')[0:3]))
#         plotdf['covered'] = covered
#         plotdf['none'] = none
#         plotdf.index = plotdf_index
#         ax = plotdf.plot.barh(rot=0, stacked = True)
#         plt.axvline(x = 0.8, color = 'r', linestyle = '-')
#         sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, -0.25))
#         plt.savefig(output.png, bbox_inches='tight', pad_inches=0.1)
#         plt.savefig(output.svg, bbox_inches='tight', pad_inches=0.1)

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

# rule cutadapt:
#     input:
#         'samples/{desc}_{strain}_{rep}_{reads}.fastq.gz',
#         #fq = expand('samples/{sample}.fastq.gz', sample = config["samples"]),
#         #get_bwa_map_input_fastqs
#     output:
#         clip='clip/{desc}_{strain}_{rep}_{reads}.clip',
#         txt='clip/{desc}_{strain}_{rep}_{reads}.txt',
#         #clip = expand('clip/{sample}.clip', sample = config["samples"]),
#         #'clip/{sample}.clip'
#     threads: 4
#     resources:
#         mem_mb=64000,
#         runtime=480
#     run:
#         print(len(config['reads']))
#         print(input)
#         if len(config['reads']) > 1:
#             if '_R1.fastq' in str(input):
#                 print("r1")
#                 shell("cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -q 0 -O 1 -m 0 -j 0 -o {output.clip} {input}")
#                 with open(str(output.txt), 'w') as f:
#                      f.write("clip success")
#                 f.close()

#             if '_R2.fastq' in str(input):
#                 print("r2")
#                 shell("cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 0 -O 1 -m 0 -j 0 -o {output.clip} {input}")
#                 with open(str(output.txt), 'w') as f:
#                      f.write("clip success")
#                 f.close()
        
#         else:
#             shell("cutadapt -a AGATCGGAAGAGCACACGTCTG -q 0 -O 1 -m 0 -j 0 -o {output.clip} {input}")
#             with open(str(output.txt), 'w') as f:
#                 f.write("clip success")
#             f.close()

            
#             else: #3
#             if go == True and not (os.path.exists(output_trim_path+file1+"_forward_paired.fq.gz")): 
#                 os.system("cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o "+
#                           output_trim_path+file1+"_forward_paired.fq.gz"+
#                           ' -p '+output_trim_path+file2+"_reverse_paired.fq.gz"+
#                           ' '+inputfolder+r1_combined[i] + " "+
#                           inputfolder+r2_combined[i])


        
rule genome_unzip:
    input:
        expand('{genome}.gz', genome = config['genome']),
    output:
        expand('{genome}', genome = config['genome']),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120
        
    shell:
        "gunzip -c {input} > {output}"


        
# rule genome_index:
#     input:
#         expand('{genome}', genome = config['genome']),
#     output:
#         fai=expand('{genome}.fai', genome = config['genome']),
#         bwa_index = expand('{genome}.amb', genome = config['genome']),
#         sizes = expand('{genome}.fai.sizes', genome = config['genome']),
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=120
#     run:
#         shell("bwa index {input}"),
#         shell("samtools faidx {input}"),
#         shell("cut -f1,2 {input}.fai > {output.sizes}")

        

        
# rule makewindows:
#     input:
#         expand('{genome}.fai.sizes', genome = config['genome']),
#         #sizes = expand('{sizes}', sizes = config['window_sizes'])
#     output:
#         file = expand('bed/{sizes}_windows_{genomeid}.bed', genomeid = config['genomeid'], sizes = config['window_sizes']),
#     threads: 1
#     resources:
#         mem_mb=8000,
#         runtime=60
        
#     run:
#         commands = expand("bedtools makewindows -w {sizes} -s {sizes} -g {{input}} > bed/{sizes}_windows_{genomeid}.bed", sizes = config['window_sizes'], genomeid = config['genomeid'])
#         for i in range(0,len(commands)):
#             shell(commands[i])          
  
        
        
        
        
        
        
        
        
        
        
# rule makewindows:
#     input:
#         expand('{genome}.fai.sizes', genome = config['genome']),
        
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
        
# rule align:
#     input:
#         clip = expand('clip/{{desc}}_{{strain}}_{{rep}}_{reads}.clip', reads = config['reads']),
#         txt = expand('clip/{{desc}}_{{strain}}_{{rep}}_{reads}.txt', reads = config['reads']),
#         fa = expand('{genome}', genome = config['genome']),
#         fai = expand('{genome}.fai', genome = config['genome']),
#     output:
#         'bam/raw/{desc}_{strain}_{rep}.bam',
#         #'samples/{desc}_{rep}_{EL}_{reads}.fastq.gz',
#     threads: 8
#     resources:
#         mem_mb=128000,
#         runtime=960,
#         disk_mb = 128000
#     run:
        
#         if len(config['reads']) > 1:
#             shell("bwa mem -v 2 -t 8 {input.fa} {input.clip} > {output}.bam")
#             shell("samtools fixmate -m {output}.bam {output}")
#         else:
#             shell("bwa mem -v 2 -t 8 {input.fa} {input.clip} > {output}")
            


#find *.clip | parallel --jobs 24 "bwa mem -v 2 -t 2 /panfs/roc/groups/0/riveramj/shared/bwaIndex_hg38/genome {} > {}_hg38.sam"

# rule bamstats:
#     input:
#         'bam/raw/{desc}_{strain}_{rep}.bam',
#     output:
#         'bam/stats/raw/{desc}_{strain}_{rep}.bamstats',
#     threads: 1
#     resources:
#         mem_mb=8000,
#         runtime=60
#     shell:
#         "samtools stats {input} > {output}"

        
# rule filterbam:
#     input:
#         'bam/raw/{desc}_{strain}_{rep}.bam',
#     output:
#         'bam/filtered/{desc}_{strain}_{rep}.filtered',
        
#     threads: 1
#     resources:
#         mem_mb=8000,
#         runtime=60
#     shell:
#         "samtools view -bhq 20 {input} -o {output}"

        
# rule sortbam:
#     input:
#         'bam/filtered/{desc}_{strain}_{rep}.filtered',
#     output:
#         allcontigs = 'bam/sorted/{desc}_{strain}_{rep}.sortedallcontigs.bam',
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
#     shell:        
#         'samtools view -bh {input.bam} $(cat contigs.txt | tr "\n" " ") > {output}'


# rule bamindex:
#     input:
#         'bam/sorted/{desc}_{strain}_{rep}.sortedallcontigs.bam' 
#     output:
#         'bam/sorted/{desc}_{strain}_{rep}.sortedallcontigs.bam.bai'
        
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     shell:
#         'samtools index {input}'


#find *.bam | parallel --jobs 24 "samtools stats {} > {}.bamstats"

#find *.bam | parallel --jobs 24 "samtools view -bhq 20 {} -o {}.filtered"

#find *.filtered | parallel --jobs 24 "samtools sort -m 10G {} -o {}.sorted"

# rule fbstats:
#     input:
#         'bam/sorted/{desc}_{strain}_{rep}.sorted',
#     output:
#         'bam/stats/fbstats/{desc}_{strain}_{rep}.fbstats',
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60

#     shell:
#         "samtools stats {input} > {output}"



# rule rmdup:
#     input:
#         'bam/sorted/{desc}_{strain}_{rep}.sorted',        
#     output:
#         'bam/rmdup/{desc}_{strain}_{rep}_rmdup.bam',
        
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     shell:
#         "samtools markdup -r {input} {output}"

        
# rule rmdup_stats:
#     input:
#         'bam/rmdup/{desc}_{strain}_{rep}_rmdup.bam',        
#     output:
#         'bam/stats/rmdup/{desc}_{strain}_{rep}.rdbamstats',
        
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     shell:
#         "samtools stats {input} > {output}"

        
        
        
# rule rpkm:
#     input:
#         bam = 'bam/rmdup/{desc}_{rep}_{EL}_rmdup.bam',
#         windows = 'bed/{sizes}_windows_{genomeid}.bed',
#     output:
#         'bg/rpkm/{desc}_{rep}_{EL}_{sizes}_{genomeid}_rpkm.bg',
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     run:
#         shell("bash scripts/rpkm.sh {input.bam} {input.windows} {output}")

        
        
# rule log:
#     input:
#         early = 'bg/rpkm/{desc}_{rep}_E_{sizes}_{genomeid}_rpkm.bg',
#         late = 'bg/rpkm/{desc}_{rep}_L_{sizes}_{genomeid}_rpkm.bg',
#     output:
#         'bg/log2/{desc}_{rep}_{sizes}_{genomeid}_log2RT.bg',
#     threads: 1
#     resources:
#         mem_mb=16000,
#         runtime=60
#     run:
#         shell(r"paste {input.early} {input.late} | awk '{{if($8 != 0 && $4 != 0){{print$1,$2,$3,log($4/$8)/log(2)}}}}' OFS='\t' > {output}")

        
        

        
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