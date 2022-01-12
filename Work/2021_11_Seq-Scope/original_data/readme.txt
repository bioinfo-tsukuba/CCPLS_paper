Date: Apr 19, 2021

Dataset Title: Seq-Scope processed datasets for liver and colon results (RDS), H&E images and method videos

Dataset Creators: Chun-Seok Cho, Jingyue Xi, Yichen Si, Sung-Rye Park, Jer-En Hsu, Myungjin Kim, Goo Jun, Hyun Min Kang & Jun Hee Lee

Dataset Contact: Jun Hee Lee leeju@umich.edu

Funding: The work was supported by the NIH (T32AG000114 to C.S.C., K01AG061236 to M.K, U01HL137182 to H.M.K., J.X. and Y.S, R01DK114131 and R01DK102850 to J.H.L., P30AG024824, P30AG013283, P30DK034933, P30DK089503, P30CA046592, P30AR069620, and U2CDK110768), the Chan Zuckerberg Initiative (to H.M.K.), the Frankel Cardiovascular Center Inaugural grant (to J.H.L. and M.K.), an American Association for the Study of Liver Diseases pilot research award (to J.H.L. and H.M.K.), Mcubed initiatives (to M.K., H.M.K. and J.H.L.), Glenn Foundation grants (to J.H.L.), Taiwan Government Scholarship for Overseas Study (to J.E.H.), Frankel Cardiovascular Center Inaugural Grant Award (to J.H.L. and M.K.), and the ADVANCE Program and the Michigan Translational Research and Commercialization for Life Sciences (MTRAC) Program (to J.H.L.), funded by the Michigan Economic Development Corporation (MEDC).

Overview: There are three experimental outputs from Seq-Scope. (1) High definition map coordinate identifier (HDMI) sequence, tile and spatial coordinate information from 1st-Seq, (2) HDMI sequence, coupled with cDNA sequence from 2nd-Seq, and (3) Histological image obtained from Hematoxylin and Eosin (H&E) staining of the tissue slice. (1) and (2) were uploaded to GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169706). (3) will be included in the current deposit. In addition, the current deposit will also include the processed RDS (single R object) data files.

Methodology: Read alignment was performed using STAR/STARsolo 2.7.5c (https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md), from which the digital gene expression (DGE) matrix was generated. Data binning was performed by dividing the imaging space into 100 um2 (10 um-sided) square grids and collapsing all HDMI-UMI information into one barcode per grid. The binned and processed DGE matrix was analyzed in the Seurat v4 package (https://github.com/satijalab/seurat/). Multiscale analysis was employed to fine tune the annotation. In addition, we employed sliding windows analysis; after the initial 10 um grid sampling, the grid was shifted both horizontally and vertically with 5 um, 2 um or 1 um intervals, producing 4, 25 and 100 times more data, respectively, sampled from the same Seq-Scope results. Then the original 10 um grid dataset was used to guide these sliding windows datasets to perform high-resolution cell type annotation. Details of the analysis could be found in our paper (see "related publications" below).

Files contained here:

Colon_10um_annotated.rds
>>> Colon dataset binned with 10um-sided square grid.
Colon_5um_2112_anchored.rds
>>> Colon dataset binned with 5um-sided square grid. Annotation was guided by 10um grid dataset.
Colon_2112_SW4X_anchored.rds
>>> Colon dataset processed through sliding windows with 5um intervals. Tile 2112. Annotation was guided by 10um grid dataset.
Colon_2112_topleft_SW25X_anchored.RDS
>>> Colon dataset processed through sliding windows with 2um intervals. Part of Tile 2112. Annotation was guided by 10um grid dataset.
Colon_2112_inset_SW100X_anchored.RDS
>>> Colon dataset processed through sliding windows with 1um intervals. Part of Tile 2112. Annotation was guided by 10um grid dataset.
Liver_normal_10um_annotated.rds
>>> Normal liver dataset binned with 10um-sided square grid
Liver_normal_Segmentation_Hepatocytes.rds
>>> Normal liver dataset binned according to H&E-based single hepatocyte segmentation.
Liver_TD_10um_annotated.rds
>>> TD liver dataset binned with 10um-sided square grid
Liver_TD2117_SW4X_anchored.RDS
>>> TD liver dataset processed through sliding windows with 5um intervals. Tile 2117. Annotation was guided by 10um grid dataset.
Liver_TD2118_SW4X_anchored.RDS
>>> TD liver dataset processed through sliding windows with 5um intervals. Tile 2118. Annotation was guided by 10um grid dataset.
Liver_TD2117_bottom_SW25X_anchored.RDS
>>> TD liver dataset processed through sliding windows with 2um intervals. Part of Tile 2117. Annotation was guided by 10um grid dataset.
Liver_TD2118_middle_SW25X_anchored.RDS
>>> TD liver dataset processed through sliding windows with 2um intervals. Part of Tile 2118. Annotation was guided by 10um grid dataset.
DKO*.jpg
>>> Raw H&E images for TD liver.
wt*.jpg
>>> Raw H&E images for Normal liver.
Colon*.jpg
>>> Raw H&E images for Normal colon.
MiSeq-DraI-100pM-mbcore-RD2-revHDMIs-pos-uniq.txt
>>> 1st-Seq coordinate information for liver dataset
MiSeq-DraI-100pM-mbcore-RD4-revHDMIs-pos-uniq.txt
>>> 1st-Seq coordinate information for colon dataset
DraI-100pM-mbcore-RD2.fastq.gz
>>> 1st-Seq raw MiSeq FASTQ output file for liver dataset
DraI_100pM_RD4.fastq.gz
>>> 1st-Seq raw MiSeq FASTQ output file for colon dataset

Related publication(s):
Cho CS, Xi J, Si Y, Park SR, Hsu JE, Kim M, Jun G, Kang HM, Lee JH. Microscopic examination of spatial transcriptome using Seq-Scope. Cell. 2021 Jun 24;184(13):3559-3572.e22. doi: 10.1016/j.cell.2021.05.010. Epub 2021 Jun 10. PubMed PMID: 34115981; PubMed Central PMCID: PMC8238917. 

Chun-Seok Cho, Jingyue Xi, Sung-Rye Park, Jer-En Hsu, Myungjin Kim, Goo Jun, Hyun-Min Kang, Jun Hee Lee "Seq-Scope: Submicrometer-resolution spatial transcriptomics for single cell and subcellular studies" (preprint) bioRxiv https://doi.org/10.1101/2021.01.25.427807 

Related data sets in NCBI Gene Expression Omnibus (GEO) repository: Cho C, Xi J, Si Y, Lee JH, Kang HM, Park S, Hsu J, Kim M, Jun G "Seq-Scope: Submicrometer-resolution spatial barcoding technology that enables microscopic examination of tissue transcriptome at single cell and subcellular levels" https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169706 

Use and Access:
This data set is made available under a Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)