<img src="resources/logo.png" style="width:3.91313in;height:1.48088in" />

**REPLAY: A standalone, user-friendly application for DNA replication
timing analysis from Repli-seq data**

REPLAY is a fast, reproducible, and fully automated application for DNA
replication timing analysis from Repli-seq data.

REPLAY is distributed as a standalone executable application, allowing
users to perform complete end-to-end analysis—from raw FASTQ files to
genome-wide replication timing profiles—without installing dependencies
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

1\. Download REPLAY

[<u>https://github.com/Rivera-Mulia-Lab/REPLAY/</u>](https://github.com/Rivera-Mulia-Lab/REPLAY/)

2\. Launch the application

> Linux:
>
> Double-click REPLAY.exe
>
> Windows:
>
> Install Windows Subsystem for Linux (WSL)
>
> XXX…..
>
> macOS:
>
> Install a virtual machine (Parallels or Oracle VM VirtualBox)
>
> XXX…..

3\. Run analysis (no coding required)

Using the GUI:

1.  Select input directory (FASTQ files)

2.  Select output directory

3.  Choose reference genome

4.  Select normalization method:

    1.  None

    2.  Quantile

    3.  Median

    4.  IQR

5.  Adjust smoothing levels (span)

6.  Click Run

The pipeline will execute automatically.

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

**Graphical User Interface**

REPLAY is designed for accessibility and ease of use.

The GUI allows users to:

- Configure all required parameters

- Launch the full pipeline with a single click

- Monitor progress in real time

- Access output files and QC metrics

No programming or command-line experience is required.

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

    5.  RT signal distributions

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_
