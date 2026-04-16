<img src="resources/logo.png" style="width:3.91313in;height:1.48088in" />

**REPLAY: A standalone, user-friendly application for DNA replication
timing analysis from Repli-seq data**

REPLAY is a fast, reproducible, and fully automated application for DNA
replication timing analysis from Repli-seq data.

REPLAY is distributed as a standalone executable application, allowing
users to perform complete end-to-end analysis—from raw FASTQ files to
genome-wide replication timing profiles—without coding, installing dependencies
or using the command line.

Through an intuitive graphical interface, users can:

- Select input and output directories

- Choose the reference genome

- Configure normalization (quantile, median, IQR)

- Adjust smoothing parameters

- Generate RT profiles and quality metrics

REPLAY integrates all processing steps, including quality control,
trimming, alignment, binning, RT log2 calculation, normalization,
smoothing, and generation of RT profiles ready visualization and
analysis, while ensuring full reproducibility.

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

**Quick start guide**

REQUIREMENTS:
10GB Memory
x86-64 architecture

1\. Download REPLAY

[<u>https://github.com/Rivera-Mulia-Lab/REPLAY/</u>](https://github.com/Rivera-Mulia-Lab/REPLAY/)

2\. Launch the application

> Linux:
>
> Double-click REPLAY
>
> Windows:
>
> Install Windows Subsystem for Linux (WSL)
> 
>     Can be installed through the Microsoft Store
>     or Open Command Prompt and type wsl --install
> 
> Run replay.bat
> 
> macOS:
>
> Download Parallels Desktop from the Mac App Store
> 
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
    
        Use buttons to download human and mouse genomes. This will take a moment.
    
        Use Presets for each genome after they are downloaded.
    
        Custom genomes can be defined under the Configuration tab.
    
        Set Smoothing and Normalization options
    
6.  Set Configuration: System Resources
    
        By default, 12GB of memory and 2 cores are used
    
        Max resources can be determined and applied with the "Max" button
    
        At least 10GB of memory is required for human genome alignment.
    
7.  Apply Advanced Configuration and Settings
    
        Configuration tab allows for the definition of custom genomes, contigs, masking, normalization options, and window sizes
     
        Advanced Settings tab allows for additional smoothing options, defining snakemake latency, and to unlock the snakemake directory if a crash occurs.
    
8.  Define Files
    
        Click the Scan + Autopopulate Button
    
        This will generate filenames for each fastq.gz file in your "Raw Reads Folder"
    
        Edit the fields to match your samples. For each unique description+replicate, there must be
        an E and an L. For paired-end, there must be an R1 and R2 for each.
     
        The app will automatically check to see if the file naming is correct.
    
9.  Copy or Link the Files

        Pressing the "Copy" button will copy your files from the "Raw Reads Folder" to the "Pipeline Reads Folder"
     
        Pressing the "Symlink" button will create a link to your files in the "Pipeline Reads Folder"
     
        If your "Raw Reads Folder" and "Pipeline Reads Folder" are the same, you can use the rename button to rename your original files.
     
10.  Run
    
         Press "Check Required Files" to ensure everything is correct.
     
         Press "Run Local" to start the analysis.
     
         Check on the status with the progress bar and by expanding the "Application Log"
     
11. QC Analysis
    
        After running is complete, the "Results" tab will contain QC images for each sample and filepaths for the data.


\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

**Outputs**

REPLAY automatically generates:

1.  Replication Timing Profiles

    1.  Raw RT (log2 Early/Late)

    2.  Normalized and smoothed RT profiles

2.  Quality Control Metrics

    1.  Mapping and filtering statistics

    2.  Coverage analysis

    3.  Autocorrelation (ACF)

    4.  RT signal on distinct exemplary genomic regions

    5.  RT signal distributions (Raw and Normalized)

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_
