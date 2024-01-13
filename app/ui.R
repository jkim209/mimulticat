library(shiny)
library(BiocManager)
library(seqinr)
library(shinydashboard)
# library(tidyverse)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(phangorn)
library(bios2mds)
library(zip)
library(ape)
library(zCompositions)
library(compositions)
library(stringr)
library(caret)
library(ggplot2)
library(randomForest)
library(data.table)
library(xgboost)
library(SHAPforxgboost)
library(fontawesome)
library(grid)
library(ggplotify)
library(remotes)
library(FSA)
library(VGAM)

source("MiDataProc.Data.Upload.R")
source("MiDataProc.Data.Input.R")
source("MiDataProc.Alpha.Diversity.R")
source("MiDataProc.Beta.Diversity.R")
source("MiDataProc.Taxa.R")
source("MiDataProc.ML.Models.R")
source("MiDataProc.ML.RF.R")
source("MiDataProc.ML.XGB.R")
source("setSliderColor.R")
source("shinyDashboardThemeDIY.R")
options(scipen=999)

# COMMENTS ------
{
  ## HOME COMMENTS -----
  
  TITLE = p("MiMultiCat: A Unified Cloud Platform for the Analysis of Microbiome Data with Multi-Categorical Responses", style = "font-size:17pt")
  HOME_COMMENT = p(strong("MiMultiCat", style = "font-size:12pt"), "is a unified cloud platform for the analysis of microbiome data with multi-categorical responses. 
                   The two key features of MiMultiCat are as follows. First, MiMultiCat streamlines a long sequence of microbiome data preprocessing and analytic procedures on user-friendly web interfaces; 
                   as such, it is easy to use for many people in various disciplines (e.g., biology, medicine, public health). 
                   Second, MiMultiCat performs both association testing and prediction modeling extensively. 
                   For association testing, MiMultiCat handles both ecological (e.g., alpha- and beta-diversity) and taxonomical (e.g., phylum, class, order, family, genus, species) contexts through covariate-adjusted or unadjusted analysis. 
                   For prediction modeling, MiMultiCat employs random forest and gradient boosting algorithms that are well-suited to microbiome data with nice visual interpretations.", style = "font-size:12pt")
  HOME_COMMENT2 = p(strong("URLs:"), "Web server (online implementation):", tags$a(href = "http://mimulticat.micloud.kr", "http://mimulticat.micloud.kr"), 
                    "; GitHub repository (local implementation):", tags$a(href = "https://github.com/jkim209/MiMultiCatGit", "https://github.com/jkim209/MiMultiCatGit"), style = "font-size:12pt")
  HOME_COMMENT3 = p(strong("Maintainers:"), "Jihun Kim (", tags$a(href = "toujours209@gmail.com", "toujours209@gmail.com"), ")", style = "font-size:12pt")
  HOME_COMMENT4 = p(strong("Reference:"), "Kim J, Jang H, Koh H. MiMultiCat: A Unified Cloud Platform for the Analysis of Microbiome Data with Multi-Categorical Responses. Bioengineering. 2024, 11(1), 60.", style = "font-size:12pt")
  
  # INPUT COMMENTS -----
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), 
                              "This should be an '.Rdata' or '.rds' file, and the data should be in the 'phyloseq' format (see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, meta/sample information, and phylogenetic tree",br(), br(),
                              
                              "Details:",br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) 
                              and columns are subjects (row names are feature IDs and column names are subject IDs).", br(), br(),
                              
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                              
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, 
                              where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                              
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiMultiCat automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                              
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiMultiCat will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in the  'phyloseq' format. For more details about 'phyloseq', see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)",  br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), 
                              " > tree <- phy_tree(biom)",br(), br(),
                              
                              "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                              and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(colnames(otu.tab), tree$tip.label)", br(),
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                              
                              strong("Reference:"), "Goodrich JK, Waters JL, Poole AC, Sutter JL, Koren O, Blekhman R, et al. Human genetics shape the gut microbiome. Cell. 2014:159(4):789-799.")
  
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects (row names are feature IDs and column names are subject IDs). 
                                    Alternatively, you can upload .biom file processed by QIIME.", br(), br(),
                                    
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). 
                                    Alternatively, you can upload .tsv file processed by QIIME.", br(), br(),
                                    
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, 
                                    where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                                    
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiMultiCat automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                                    The tip labels of the phylogenetic tree are feature IDs." ,br(), br(),
                                    
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiMultiCat will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, 
                                     feature table (otu.tab.txt), taxonomic table (tax.tab.txt), and metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE)", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     "> sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE)", br(),
                                     "> tree <- read.table(file = 'tree.tre')", br(), br(),
                                     
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                                     
                                     strong("Reference:"), "Goodrich JK, Waters JL, Poole AC, Sutter JL, Koren O, Blekhman R, et al. Human genetics shape the gut microbiome. Cell. 2014:159(4):789-799.")
  
  # QC COMMENTS -----
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. 
                         Alternatively, you can type 'Fungi' for ITS data or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove units that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per unit.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.02% (0.0002).",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p("Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. 
                            Default is \"\", \"metagenome\", \"gut metagenome\", \"mouse gut metagenome\".",
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p("Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings 
                            (i.e., taxonomic names that contain the specified character strings). Multiple character strings should be separated by a comma. 
                            Default is \"uncultured\", \"incertae\", \"Incertae\", \"unidentified\", \"unclassified\", \"unknown\".",
                            style = "font-size:11pt")
  
  # MiMultiCat Reference -----
  
  MiMultiCat_REFERENCE <- p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)")
  
  # ALPHA COMMENTS -----
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  ALPHA_ANOVAF_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                             "2. Tukey JW. Commparing Individual Means in the Analysis of Variance. Biometrics. 1949;5(2):99-114", br())
  ALPHA_KRUSKAL_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                              "2. Kruskal WH, Wallis WA. Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association. 1952;47(260):583-621", br(), 
                              "3. Dunn OH. Multiple Comparisons Using Rank Sums. Technometrics. 1964;6(3):241-252")
  ALPHA_PROPODDS_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                               "2. McCullaph P. Regression models for ordinal data. J R Stat Soc Series B. 1980;42(2):109-142.")
  ALPHA_MULTINOM_REFERENCE <- p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)")
  
  # BETA COMMENTS -----
  
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative β-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  BETA_DA_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                        "2.Jiang Z, He M, Chen J, Zhao N, Zhan X. MiRKAT-MC: A Distance-Based Microbiome Kernel Association Test With Multi-Categorical Outcomes. Front Genet. 2022;13:841764.")
  BETA_PERMANOVA_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                               "2. Anderson MA. new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46.", br(),
                               "3. McArdle BH, Anderson MJ. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology. 2001;82(1):290-297.")
  
  # DT COMMENTS -----
  
  DATA_TRANSFORM_COMMENT = p("Transform the data into four different formats (1) CLR (centered log ratio) (Aitchison, 1982), (2) Count (Rarefied) (Sanders, 1968), (3) Proportion, (4) Arcsine-root 
                             for each taxonomic rank (phylum, class, order, familiy, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77", br(),
                               "2. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")

  TAXA_ANOVAF_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                             "2. Tukey JW. Commparing Individual Means in the Analysis of Variance. Biometrics. 1949;5(2):99-114", br(),
                            "3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289-300.")
  TAXA_KRUSKAL_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                              "2. Kruskal WH, Wallis WA. Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association. 1952;47(260):583-621", br(), 
                              "3. Dunn OH. Multiple Comparisons Using Rank Sums. Technometrics. 1964;6(3):241-252", br(),
                             "4. Benjamini Y, Hochberg Y. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289-300.")
  TAXA_PROPODDS_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                               "2. McCullaph P. Regression models for ordinal data. J R Stat Soc Series B. 1980;42(2):109-142.", br(),
                              "3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289-300.")
  TAXA_MULTINOM_REFERENCE <- p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                               "2. Benjamini Y, Hochberg Y. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289-300.")
  
  # RF COMMENTS -----
  
  RF_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                   "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32", br())
  RF_REFERENCE_CLR = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                       "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32", br(),
                       "3. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77")
  RF_REFERENCE_RC = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                      "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32", br(),
                      "3. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")
  
  # XGB COMMENTS -----
  
  XGB_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                    "2. Friedman JH. Greedy function approximation: A gradient boosting machine. Ann Stat. 2001;29(5):1189-1232",br(),
                    "3. Chen T, Guestrin C. XGBoost: A scalable tree boosting system. in Proc the 22nd ACM SIGKDD Int Conf KDD. ACM. 2016;785-794", br(),
                    "4. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. in Proc Adv Neural Inf Process Syst. 2017;4765-4774.")
  XGB_REFERENCE_CLR = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                        "2. Friedman JH. Greedy function approximation: A gradient boosting machine. Ann Stat. 2001;29(5):1189-1232",br(),
                        "3. Chen T, Guestrin C. XGBoost: A scalable tree boosting system. in Proc the 22nd ACM SIGKDD Int Conf KDD. ACM. 2016;785-794", br(),
                        "4. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. in Proc Adv Neural Inf Process Syst. 2017;4765-4774.",br(),
                        "5. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77")
  XGB_REFERENCE_RC = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                       "2. Friedman JH. Greedy function approximation: A gradient boosting machine. Ann Stat. 2001;29(5):1189-1232",br(),
                       "3. Chen T, Guestrin C. XGBoost: A scalable tree boosting system. in Proc the 22nd ACM SIGKDD Int Conf KDD. ACM. 2016;785-794", br(),
                       "4. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. in Proc Adv Neural Inf Process Syst. 2017;4765-4774.",br(),
                       "5. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")
}

# UI ---------------------------------------------------------------------------
{
  ui = dashboardPage(
    title = "MiMultiCat",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      sidebarMenu(id = "side_menu",
                  menuItem("Home", tabName = "home"),
                  menuItem("Data Processing",
                           menuSubItem("Data Input", tabName = "step1", 
                                       icon = fontawesome::fa("upload", margin_left = "0.3em", margin_right = "0.1em")),
                           menuSubItem("Quality Control", tabName = "step2", 
                                       icon = fontawesome::fa("chart-bar", margin_left = "0.3em")),
                           menuSubItem(span(span("Diversity Calculation /", br()), 
                                            span("Data Transformation", style = "margin-left: 23px")), tabName = "divdtCalculation", 
                                       icon = fontawesome::fa("calculator", margin_left = "0.3em", margin_right = "0.25em"))),
                  menuItem("Association",
                           menuSubItem("Alpha Diversity", tabName = "alphaDiv", 
                                       icon = fontawesome::fa("font", margin_left = "0.2em", margin_right = "0.1em")),
                           menuSubItem("Beta Diversity", tabName = "betaDiv", 
                                       icon = fontawesome::fa("bold", margin_left = "0.3em", margin_right = "0.15em")),
                           menuSubItem("Taxonomic Analysis", tabName = "taxa", 
                                       icon = fontawesome::fa("diagram-project"))),
                  menuItem("Prediction",
                           menuSubItem("Random Forest", tabName = "rf", 
                                       icon = fontawesome::fa("network-wired")),
                           menuSubItem("Gradient Boosting", tabName = "xgb", 
                                       icon = fontawesome::fa("diagram-project")))
                  )
      ),
    
    dashboardBody(
      
      # THEME -----
      shinyDashboardThemeDIY(
        
        # general
        appFontFamily = "Arial"
        ,appFontColor = "rgb(0,0,0)"
        ,primaryFontColor = "rgb(0,0,0)"
        ,infoFontColor = "rgb(0,0,0)"
        ,successFontColor = "rgb(0,0,0)"
        ,warningFontColor = "rgb(0,0,0)"
        ,dangerFontColor = "rgb(0,0,0)"
        ,bodyBackColor = "rgb(255,255,255)"
        
        # header
        ,logoBackColor = "rgb(65,186,81)"
        
        ,headerButtonBackColor = "rgb(65,186,81)"
        ,headerButtonIconColor = "rgb(65,186,81)"
        ,headerButtonBackColorHover = "rgb(65,186,81)"
        ,headerButtonIconColorHover = "rgb(0,0,0)"
        
        ,headerBackColor = "rgb(65,186,81)"
        ,headerBoxShadowColor = "#aaaaaa"
        ,headerBoxShadowSize = "0px 0px 0px"
        
        # sidebar
        ,sidebarBackColor = "rgb(24,31,41)"
        ,sidebarPadding = 0
        
        ,sidebarMenuBackColor = "transparent"
        ,sidebarMenuPadding = 0
        ,sidebarMenuBorderRadius = 0
        
        ,sidebarShadowRadius = ""
        ,sidebarShadowColor = "0px 0px 0px"
        
        ,sidebarUserTextColor = "rgb(24,31,41)"
        
        ,sidebarSearchBackColor = "rgb(255, 255, 255)"
        ,sidebarSearchIconColor = "rgb(24,31,41)"
        ,sidebarSearchBorderColor = "rgb(24,31,41)"
        
        ,sidebarTabTextColor = "rgb(210,210,210)"
        ,sidebarTabTextSize = 14
        ,sidebarTabBorderStyle = "none"
        ,sidebarTabBorderColor = "none"
        ,sidebarTabBorderWidth = 0
        
        ,sidebarTabBackColorSelected = "rgb(45,52,63)"
        ,sidebarTabTextColorSelected = "rgb(252,255,255)"
        ,sidebarTabRadiusSelected = "0px"
        
        ,sidebarTabBackColorHover = "rgb(67,75,86)"
        ,sidebarTabTextColorHover = "rgb(252,255,255)"
        ,sidebarTabBorderStyleHover = "none"
        ,sidebarTabBorderColorHover = "none"
        ,sidebarTabBorderWidthHover = 0
        ,sidebarTabRadiusHover = "0px"
        
        # boxes
        ,boxBackColor = "rgb(245,245,245)"
        ,boxBorderRadius = 3
        ,boxShadowSize = "0px 0px 0px"
        ,boxShadowColor = "rgba(0,0,0,0)"
        ,boxTitleSize = 16
        ,boxDefaultColor = "rgb(210,214,220)"
        ,boxPrimaryColor = "rgb(35, 49, 64)"
        ,boxInfoColor = "rgb(65,186,81)"
        ,boxSuccessColor = "rgb(102,199,115)"
        ,boxWarningColor = "rgb(244,156,104)"
        ,boxDangerColor = "rgb(255,88,55)"
        
        ,tabBoxTabColor = "rgb(255,255,255)"
        ,tabBoxTabTextSize = 14
        ,tabBoxTabTextColor = "rgb(0,0,0)"
        ,tabBoxTabTextColorSelected = "rgb(35, 49, 64)"
        ,tabBoxBackColor = "rgb(255,255,255)"
        ,tabBoxHighlightColor = "rgb(65,186,81)"
        ,tabBoxBorderRadius = 0
        
        # inputs
        ,buttonBackColor = "rgb(245,245,245)"
        ,buttonTextColor = "rgb(0,0,0)"
        ,buttonBorderColor = "rgb(24,31,41)"
        ,buttonBorderRadius = 3
        
        ,buttonBackColorHover = "rgb(227,227,227)"
        ,buttonTextColorHover = "rgb(100,100,100)"
        ,buttonBorderColorHover = "rgb(200,200,200)"
        
        ,textboxBackColor = "rgb(255,255,255)"
        ,textboxBorderColor = "rgb(200,200,200)"
        ,textboxBorderRadius = 0
        ,textboxBackColorSelect = "rgb(245,245,245)"
        ,textboxBorderColorSelect = "rgb(200,200,200)"
        
        # tables
        ,tableBackColor = "rgb(255, 255, 255)"
        ,tableBorderColor = "rgb(245, 245, 245)"
        ,tableBorderTopSize = 1
        ,tableBorderRowSize = 1
        
      ),
      
      # STYLE -----
      
      ## CONTENTS -----
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      
      ## PROGRESS BAR -----
      tags$head(tags$style(HTML('.progress-bar {background-color: (102,199,115);}'))),
      
      ## PRETTY RADIO BUTTON -----
      tags$head(tags$style(HTML('
      .pretty input:checked~.state.p-primary label:after, .pretty.p-toggle .state.p-primary label:after {
    background-color: #66C773!important;}'))),
      
      ## SLIDER -----
      setSliderColor(rep("#66C773", 100), seq(1, 100)),
      chooseSliderSkin("Flat"),
      
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      tabItems(
        
        # Home -----
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT, 
                    div(tags$img(src="MiMultiCat_Home_Img.png", height = 720, width = 600), style = "text-align: center"),
                    HOME_COMMENT2, HOME_COMMENT3, HOME_COMMENT4)),
        
        ## 0. DATA INPUT -----
        
        tabItem(tabName = "step1", br(),
                fluidRow(column(width = 6,
                                box(
                                  width = NULL, status = "info", solidHeader = TRUE,
                                  title = strong("Data Input", style = "color:white"),
                                  selectInput("inputOption", h4(strong("Data type")), 
                                              c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                                  div(id = "optionsInfo", 
                                      tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), 
                                      tags$p("", style = "margin-bottom:-8px"), style = "margin-top: -15px"),
                                  uiOutput("moreOptions")
                                  )
                                ),
                         column(width = 6, style='padding-left:0px', 
                                uiOutput("addDownloadinfo")
                                )
                         )
                ),
        
        ## 1-1. QC ----
        
        tabItem(tabName = "step2", br(),
            fluidRow(column(width = 3,  style = "padding-left:+15px",
                            # Quality Control
                            box(
                              width = NULL, status = "info", solidHeader = TRUE,
                              title = strong("Quality Control", style = "color:white"),
                              textInput("kingdom", h4(strong("Kingdom")), value = "Bacteria"),
                              QC_KINGDOM_COMMENT,
                              
                              tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                              tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'),
                              
                              sliderInput("slider1", h4(strong("Library size")), 
                                          min=0, max=10000, value = 3000, step = 1000),
                              QC_LIBRARY_SIZE_COMMENT1,
                              QC_LIBRARY_SIZE_COMMENT2,
                              
                              sliderInput("slider2", h4(strong("Mean proportion")), 
                                          min = 0, max = 0.1, value = 0.02, step = 0.001,  post  = " %"),
                              QC_MEAN_PROP_COMMENT1,
                              QC_MEAN_PROP_COMMENT2,
                              
                              br(),
                              p(" ", style = "margin-bottom: -20px;"),
                              
                              h4(strong("Errors in taxonomic names")),
                              textInput("rem.str", label = "Complete match", value = ""),
                              QC_TAXA_NAME_COMMENT1,
                              
                              textInput("part.rem.str", label = "Partial match", value = ""),
                              QC_TAXA_NAME_COMMENT2,
                              
                              actionButton("run", (strong("Run!")), class = "btn-info"), 
                              p(" ", style = "margin-bottom: +10px;"), 
                              p(strong("Attention:"),"You have to click this Run button to perform diversity calculation, 
                                data transformation and further analyses.", style = "margin-bottom:-10px"), br()),
                            
                            uiOutput("moreControls"),
                            uiOutput("qc_reference")),
                     
                     column(width = 9, style = "padding-left:+10px",
                            box(
                              width = NULL, status = "info", solidHeader = TRUE,
                              fluidRow(width = 12,
                                       status = "info", solidHeader = TRUE, 
                                       valueBoxOutput("sample_Size", width = 3),
                                       valueBoxOutput("OTUs_Size", width = 3),
                                       valueBoxOutput("phyla", width = 3),
                                       valueBoxOutput("classes", width = 3)),
                              
                              fluidRow(width = 12, 
                                       status = "info", solidHeader = TRUE,
                                       valueBoxOutput("orders", width = 3),
                                       valueBoxOutput("families", width = 3),
                                       valueBoxOutput("genera", width = 3),
                                       valueBoxOutput("species", width = 3)),
                              
                              fluidRow(style = "position:relative",
                                       tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                              tabPanel("Histogram",
                                                       plotlyOutput("hist"),
                                                       sliderInput("binwidth", "# of Bins:",
                                                                   min = 0, max = 100, value = 50, width = "100%")),
                                              tabPanel("Box Plot", 
                                                       plotlyOutput("boxplot"))),
                                       tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                              tabPanel("Histogram",
                                                       plotlyOutput("hist2"),
                                                       sliderInput("binwidth2", "# of Bins:",
                                                                   min = 0, max = 100, value = 50, width = "100%")),
                                              tabPanel("Box Plot",
                                                       plotlyOutput("boxplot2")))))))
            ),
        
        ## 1-2. DC & DT -----
        
        tabItem(tabName = "divdtCalculation", br(),
                fluidRow(column(width = 6, style = "padding-left:+15px",
                                box(title = strong("Diversity Calculation & Data Transformation", style = "color:white"), 
                                    width = NULL, status = "info", solidHeader = TRUE, 
                                    
                                    strong(p("Diversity Calculation", style = "font-size:12pt")),
                                    ALPHA_COMMENT, 
                                    BETA_COMMENT, 
                                    
                                    strong(p("Data Transformation", style = "font-size:12pt")),
                                    DATA_TRANSFORM_COMMENT, 
                                    actionButton("divdtCalcRun", (strong("Run!")), class = "btn-info")),
                                uiOutput("divdtDownload")),
                         
                         column(width = 6, style='padding-left:0px',
                                box(title = strong("References", style = "color:white"), 
                                    width = NULL, status = "info", solidHeader = TRUE,
                                    strong(p("Alpha Diversity", style = "font-size:12pt")), ALPHA_REFERENCES,
                                    strong(p("Beta Diversity", style = "font-size:12pt")), BETA_REFERENCES,
                                    strong(p("Data Transformation", style = "font-size:12pt")), DATA_TRANSFORM_REFERENCE)))
              ),
        
        ## 2. ALPHA ------
        
        tabItem(tabName = "alphaDiv", br(),
                tabPanel(
                  title = NULL,
                  sidebarLayout( 
                    position = "left",
                    sidebarPanel(width = 3,
                                 shinyjs::hidden(
                                   uiOutput("alpha_primvars"),
                                   uiOutput("alpha_var_type"),
                                   uiOutput("alpha_primvars_ref"),
                                   uiOutput("alpha_primvars_rename"),
                                   uiOutput("alpha_covariate"),
                                   uiOutput("alpha_method"),
                                   uiOutput("alpha_run"),
                                   uiOutput("alpha_downloadTable"),
                                   uiOutput("alpha_references"))),
                    
                    mainPanel(width = 9,
                              fluidPage(width = NULL, uiOutput("alpha_display_results"), uiOutput("alpha_display_pair")),
                              uiOutput("barPanel"), br(), br())))),
        
        ## 3. BETA ------
        
        tabItem(tabName = "betaDiv", br(),
                tabPanel(
                  title = NULL,
                  sidebarLayout( 
                    position = "left",
                    sidebarPanel(width = 3,
                                 shinyjs::hidden(
                                   uiOutput("beta_primvars"),
                                   uiOutput("beta_var_type"),
                                   uiOutput("beta_primvars_ref"),
                                   uiOutput("beta_primvars_rename"),
                                   uiOutput("beta_covariate"),
                                   uiOutput("beta_method"),
                                   uiOutput("beta_run"),
                                   uiOutput("beta_downloadTabUI"),
                                   uiOutput("beta_reference"))),
                    
                    mainPanel(width = 9,
                              fluidPage(width = 12, uiOutput("beta_nom_results")),
                              uiOutput("beta_barPanel"), br(), br())))),
        
        ## 4. TAXA ------
        
        tabItem(tabName = "taxa", br(),
                tabPanel(
                  title = NULL,
                  sidebarLayout( 
                    position = "left",
                    sidebarPanel(width = 3,
                                 shinyjs::hidden(
                                   uiOutput("taxa_primvars"),
                                   uiOutput("taxa_var_type"),
                                   uiOutput("taxa_primvars_ref"),
                                   uiOutput("taxa_primvars_rename"),
                                   uiOutput("taxa_covariate"),
                                   uiOutput("taxa_method"),
                                   uiOutput("taxa_run"),
                                   uiOutput("taxa_downloadTable"),
                                   uiOutput("taxa_references"))),
                    
                    mainPanel(width = 9,
                              fluidPage(width = NULL, 
                                        uiOutput("taxa_display_global"),
                                        uiOutput("taxa_display_forest"),
                                        div(id = "taxa_display_area", style='height:550px;overflow-y: scroll;', uiOutput("taxa_display")), 
                                        uiOutput("taxa_display_dend")),
                              br(), 
                              br(), 
                              br(), 
                              br(), 
                              uiOutput("taxa_barPanel"))))),
        
        ## 5-1. RF ------
        
        tabItem(tabName = "rf", br(),
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               shinyjs::hidden(
                                 uiOutput("rf_nom_data_input"),
                                 uiOutput("rf_nom_data_input_opt"),
                                 uiOutput("rf_nom_train_setting"),
                                 uiOutput("rf_nom_downloadTabUI"),
                                 uiOutput("rf_nom_reference")
                                 )),
                  mainPanel(width = 9,
                            fluidRow(width = 12, uiOutput("rf_nom_results"))))
                ),
        
        ## 5-2. XGB ------
        
        tabItem(tabName = "xgb", br(),
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               shinyjs::hidden(
                                 uiOutput("xgb_nom_data_input"),
                                 uiOutput("xgb_nom_data_input_opt"),
                                 uiOutput("xgb_nom_train_setting"),
                                 uiOutput("xgb_nom_downloadTabUI"),
                                 uiOutput("xgb_nom_reference")
                                 )),
                  mainPanel(width = 9,
                            fluidRow(width = 12, uiOutput("xgb_nom_results"))))
                )
      )
    )
  )
}
