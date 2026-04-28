<img src="resources/logo.png" style="width:3.91313in;height:1.48088in" />

**REPLAY: A standalone, user-friendly application for DNA replication
timing analysis from Repli-seq data**

REPLAY is a fast, reproducible, and fully automated application for DNA
replication timing analysis from Repli-seq data.

REPLAY is a standalone, end-to-end executable application designed to democratize Repli-seq analysis. By bridging a high-performance Snakemake/Apptainer backend with an intuitive PySide6 graphical interface, REPLAY allows researchers to transform raw FASTQ files into publication-quality RT profiles and diagnostic reports with a single click.

Key features:

- Zero-Installation: Standalone executable for Linux (Windows via WSL2 and macOS via virtualization). No manual environment setup required.

- End-to-End: Handles everything from raw zipped FASTQ processing (adapter trimming, alignment, filtering) to final smoothing and normalization.

- Scientific Rigor: Integrated QC metrics to validate biological data fidelity.

- Reproducible: Encapsulated in Apptainer containers to ensure bit-wise identical results across different systems.

- Data normalization: Integrated strategies for RT profiles normalization  (quantile, IQR, median) and smoothing (LOESS, Gaussian) strategies.


\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

**Quick start guide**

REQUIREMENTS:

Processor (CPU): x86-64 architecture. Minimum 2 cores. 6 cores or higher recommended. The Snakemake backend automatically parallelizes tasks based on available threads.

Memory (RAM): 10 GB minimum. 32 GB or higher is recommended for processing large mammalian genomes (e.g., human or mouse) at high resolutions ($<10\text{ kb}$ bins).

Storage: 
- Containers: ~1 GB for the executable and Apptainer/Singularity image.
- Genome + Index: ~8.3 GB for hg38
- Input Data: Variable (depends on FASTQ size).
- Output Data: Assume at least ~13x input for intermediate files. e.g. 1.21GB fastq.gz -> 14.8GB trimmed reads/bam/bedgraphs
  


1\. Download REPLAY

[<u>https://github.com/Rivera-Mulia-Lab/REPLAY/</u>](https://github.com/Rivera-Mulia-Lab/REPLAY/)

The files REPLAY (the executable) and snakemake_repliseq.sif must be downloaded separately and placed into
the folder due to their size.

2\. Launch the application

> Linux:
>
> Double-click REPLAY
>
> Windows:
>
> Install Windows Subsystem for Linux (WSL)
> 
> Can be installed through the Microsoft Store
> or Open Command Prompt and type wsl --install
> 
> Run replay.bat
> 
> macOS:
>
> Download VirtualBox, Lima, or Parallels Desktop from the Mac App Store and install linux distributions.
> 
> For Parallels:
> Open Parallels Desktop and Click "+", Free Systems, Select Ubuntu, and Download
> 
> In the VM, navigate to the REPLAY folder and double-click REPLAY
> 
> 
> 

3\. Run analysis 

Using the GUI:

1.  Set file path for "Raw Reads Folder" (fastq.gz files)

        This is where your original files are kept.
    
2.  Set file path for "Pipeline Reads Folder"
   
        This is where reads are referenced in the pipeline.
    
3.  Select "Output Folder"
   
        This is where all outputs will be generated

4.  Choose Read Options
   
        Default barcodes are standard Illumina sequences.
    
        Only Single-Read or Paired-End can be performed at once.

5.  Set Configuration: Genome and Analysis
    
        Use buttons to download the human and mouse genomes. This will take some minutes depending on the internet connection
        speed and processor.
    
        Use Presets for each genome after they are downloaded.
    
        Custom genomes can be defined under the Configuration tab.
    
        Set Smoothing and Normalization options
    
7.  Set Configuration: System Resources
    
        By default, 12GB of memory and 2 cores are used
    
        Max resources can be determined and applied with the "Max" button
    
        At least 10GB of memory is required for human genome alignment.
    
8.  Apply Advanced Configuration and Settings
    
        Configuration tab allows for the definition of custom genomes, contigs, masking, normalization options, and window
        sizes
     
        Advanced Settings tab allows for additional smoothing options, defining snakemake latency, and to unlock the snakemake
        directory if a crash occurs.
    
9.  Define Files
    
        Click the Scan + Autopopulate Button
    
        This will generate filenames for each fastq.gz file in your "Raw Reads Folder"
    
        Edit the fields to match your samples. For each unique description+replicate, there must be
        an E and an L. For paired-end, there must be an R1 and R2 for each.
     
        The app will automatically check to see if the file naming is correct.
    
10. Copy or Link the Files

        Pressing the "Copy" button will copy your files from the "Raw Reads Folder" to the "Pipeline Reads Folder"
     
        Pressing the "Symlink" button will create a link to your files in the "Pipeline Reads Folder"
     
        If your "Raw Reads Folder" and "Pipeline Reads Folder" are the same, you can use the rename button to rename your
        original files.
     
12.  Run
    
         Press "Check Required Files" to ensure everything is correct.
     
         Press "Run Local" to start the analysis.
     
         Check on the status with the progress bar and by expanding the "Application Log"
     
13. QC Analysis
    
        After running is complete, the "Results" tab will contain QC images for each sample and filepaths for the data.


\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

**Outputs and Diagnostics**

REPLAY automatically generates:

1.  Replication Timing Profiles ready for visualization using the UCSC Genome Browser or IGV.

    1.  Raw RT (log2 Early/Late) in BedGraph format.

    2.  Normalized and smoothed RT profiles in BedGraph format.

2.  Quality Control Metrics

    1.  Mapping and filtering statistics

    2.  Coverage analysis

    3.  Autocorrelation (ACF)

    4.  RT signal on distinct exemplary genomic regions

    5.  RT signal distributions (Raw and Normalized)

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

**Citation**

If you use REPLAY in your research, please cite our preprint:
Dickinson, Q., Yu, C., Rivera-Mulia, J.C. REPLAY: A reproducible and user-friendly application for DNA replication timing analysis from Repli-seq data.
