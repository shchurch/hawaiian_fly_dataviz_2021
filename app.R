library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(ggcorrplot)
library(png)
library(viridis)
theme_set(theme_classic())
set.seed(84095)
load(file="data/genetree_all_dmel_id.RData")
load(file="data/ovary_changes.RData")
load(file="data/full_correlation_matrix.RData")
load(file="data/no_outliers_correlation_matrix.RData")
load(file="data/ovary_average_ratio.RData")
source("data/qgraph_as_ggraph.R")

species_names <- c('Drosophila_primaeva','Drosophila_sproati','Drosophila_picticornis','Drosophila_macrothrix','Drosophila_mimica','Drosophila_cfdives','Drosophila_nanella','Scaptomyza_varia','Scaptomyza_cyrtandrae','Scaptomyza_varipicta','Drosophila_atroscutellata','Drosophila_tanythrix')
species_id <- c('008D','106A','025A','055A','040C','16_1','002D','CFB','088B','020A','029A','043D');names(species_id) <- species_names
species_codes <- c('Dpri','Dspr','Dpic','Dmac','Dmim','Dcfd','Dnan','Svar','Scyr','Svpt','Datr','Dtan');names(species_codes) <- species_names
species_order <-  c('Drosophila_primaeva','Drosophila_mimica','Drosophila_nanella','Drosophila_atroscutellata','Drosophila_tanythrix','Drosophila_cfdives','Drosophila_sproati','Drosophila_macrothrix','Drosophila_picticornis','Scaptomyza_cyrtandrae','Scaptomyza_varipicta','Scaptomyza_varia')

selected_gene_families <- c("Yp1;Yp2;Yp3","Octbeta1R;Octbeta2R","Doa","Haspin","nos")
gene_family_names <- genetree_all_dmel_id %>% filter(genetree %in% rownames(full_cormat)) %>% filter(!name %in% selected_gene_families) %>% pull(name) %>% unique %>% sort

gene_family_names <- c(selected_gene_families,gene_family_names)

ui <- dashboardPage(skin="black",
  dashboardHeader(title="Evolutionary changes in ovary-biased expression across 12 species of Hawaiian Drosophilidae flies",titleWidth = 1000),
  dashboardSidebar(
    selectizeInput(inputId="variable",choices=gene_family_names,label="Select a gene or gene group",options = list(create = TRUE)),
    div(style="padding-left:10px;width:90%;text-align:left","Genes have been grouped together by homology, as inferred with the software agalma (bitbucket.org/caseywdunn), which uses sequence similarity to cluster genes.",br(),br(),"Note: not all genes will be represented. For the purposes of this app we have restricted to only genes represented in all the datasets of twelve species."),
    checkboxInput(inputId="checkbox",label = "filter outliers before calculating correlations",value = FALSE)
  ),
  dashboardBody(
    tags$script("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';"),
    tags$head(tags$style(
        type="text/css",
        "#dissection img {max-width: 100%; width: 100%; padding: 10px; height: auto: middle}",
        "#plot_tree img {max-width: 100%; width: 100%; padding: 10px; height: auto}",
        "#cor_legend img {max-width: 75%; width: 100%; padding: 10px; height: auto}"
    )),
    fixedRow(
      column(width=4,
        box(
          #height=485,
          width = NULL,
          HTML('<h2>Overview</h1>This R::Shiny app allows you to explore the correlation of evolutionary changes in the ovary-biased expression using RNA sequence data from twelve species of Hawaiian Drosophilidae flies. You can select any gene  for which to display expression data. To generate this dataset, we generated RNA sequence datasets from 12 species of Hawaiian flies, and compared them using phylogenetic methods. Check out our paper describing these data and analyses <a href="https://doi.org/10.1101/2021.11.30.470652">here</a>.')
        ),
        box(
          solidHeader=TRUE,
          width=NULL,
          height=330,
          imageOutput(outputId="dissection")
        ),
      ),
      column(width=3,
        box(
          solidHeader=TRUE,
          width=NULL,
          imageOutput(outputId="plot_tree")
        ),
        box(
          width = NULL,
          HTML("<h4>Evolutionary tree</h4>The phylogeny of the twelve species studied here, numbered by branch.")
        )
      ),
      column(width=5,
        box(
          width=NULL,
          solidHeader=TRUE,
          plotOutput(outputId="plot_exp")
        ),
        box(
          width = NULL,
          HTML('<p><h4>Figure 1</h4>The expression ratio between the ovary and carcass for each transcript.<br/><span style="color: red">Red</span> = expressed more in the ovary, <span style="color: blue">Blue</span>  = expressed more in the carcass.</p>')
        )
      )
    ),
    fixedRow(
      column(width=7,
        box(
          #title="Expression changes",
          width = NULL,
          solidHeader=TRUE,
          #status="primary",
          plotOutput(outputId="plot_changes")
        ),
        box(
          width = NULL,
          HTML("<p><h4>Figure 2</h4>The distribution of evolutionary changes in expression bias along each branch in the phylogeny. The black point indicates the selected gene. Other points are changes in other genes, colored by strength and direction of correlation.")
        )
      ),
      column(width=2,
        box(
          solidHeader=TRUE,
          width=NULL,
          imageOutput(outputId="cor_legend")
        ),
        box(
          width = NULL,
          "Legend for both correlation plots."
        )
      ),
      column(width=3,
        box(
          #title="Expression correlation",
          width = NULL,
          solidHeader=TRUE,
          #status="primary",
          plotOutput(outputId="plot_cor")
        ),
        box(
          width = NULL,
          HTML("<h4>Figure 3</h4>Correlation heatmap for genes with a strong evolutionary correlation (abs. value > 0.825) of expression to the selected gene.")
        )
      ),

    )
  )
)

server <- function(input, output) { 

  output$dissection <- renderImage({
    return(list(
      src="data/dissection-01.png",
      contentType="image/png",
      Alt="dissection diagram"))
  },deleteFile=FALSE)

  output$cor_legend <- renderImage({
    return(list(
      src="data/correlation_legend-01.png",
      contentType="image/png",
      Alt="correlation legend"))
  },deleteFile=FALSE)

  output$plot_tree <- renderImage({
    return(list(
      src="data/labeled_ultrametric_tree-01.png",
      contentType="image/png",
      Alt="species tree"))
  },deleteFile=FALSE)
 
  output$plot_exp <- renderPlot({
    check <- input$checkbox   
    if(check == TRUE){
      cormat <- cormat_no_outliers
    } else {
      cormat <- full_cormat
    }

    target <- input$variable
    if(target %in% genetree_all_dmel_id$name){
      target_genetree <- genetree_all_dmel_id %>% 
        filter(name == target) %>% 
        pull(genetree) %>% unique
    }
  
    if(any(target_genetree %in% rownames(cormat))){
  
      ovary_changes_cor <- left_join(ovary_changes,
        data.frame(key = names(cormat[target_genetree,]),
        correlation = cormat[target_genetree,]),by="key") %>% 
        na.omit
      bias_colors <- setNames(c("red","blue"),c("ovary","carcass"))
    
  
      plot_exp <- ggplot(ave_ratio_summary %>% filter(genetree == target_genetree),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + geom_point(size=3) +
        scale_color_manual(values = bias_colors) + 
        geom_vline(xintercept=0,linetype="dashed",size=0.5) + 
        scale_x_continuous(limits = c(-5,5)) +
        ggtitle(paste("expression bias of ",target," transcripts",sep="")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(text = element_text(size=15)) +
        theme(axis.title.y = element_blank()) +
        theme(axis.title.x = element_blank()) +
        theme(legend.position="none") 
      plot_exp

    } else {
      text = paste("\n   gene family not present in dataset")
      ggplot() + annotate("text", x = 0, y = 0, size=6, label = text) + 
        theme_void()      
    }
  })

  output$plot_changes <- renderPlot({
    check <- input$checkbox   
    if(check == TRUE){
      cormat <- cormat_no_outliers
    } else {
      cormat <- full_cormat
    }

    target <- input$variable
    if(target %in% genetree_all_dmel_id$name){
      target_genetree <- genetree_all_dmel_id %>% 
        filter(name == target) %>% 
        pull(genetree) %>% unique
    }
  
    if(any(target_genetree %in% rownames(cormat))){
  
      ovary_changes_cor <- left_join(ovary_changes,
        data.frame(key = names(cormat[target_genetree,]),
          correlation = cormat[target_genetree,]),by="key") %>% 
        na.omit
      correlation_threshold <- 0.825
      target_changes <- ovary_changes_cor %>% filter(abs(correlation) > correlation_threshold)
      target_cormat <- target_changes %>% select(name,child_node,scaled_change) %>% 
        spread(.,name,scaled_change)  %>% 
        select(-child_node) %>% 
        cor(.,method="pearson",use="pairwise.complete.obs")

      node_order <- c("21","13","22","19","20","14","15","16","18","17","3","2","1","6","5","4","12","11","10","9","8","7")

      plot_changes <- ggplot(ovary_changes_cor,aes(x=factor(branch,levels=node_order),y=scaled_change,group=name,color=correlation)) +
        viridis::scale_color_viridis(option = "B") + 
        geom_jitter(size=1.5,pch=16,aes(alpha=abs(correlation))) + 
        geom_point(data=ovary_changes %>% filter(name == target),color="white",fill="black",size=3,stroke=1,pch=21) + 
        theme(legend.position="none") + 
        theme(text = element_text(size=10)) +
        ylab("scaled evolutionary\nchange") + 
        xlab("phylogenetic branch")
      plot_changes

    } else {
      text = paste("\n   gene family not present in dataset")
      ggplot() + annotate("text", x = 0, y = 0, size=6, label = text) + 
        theme_void()      
    }
  })

  output$plot_cor <- renderPlot({
    check <- input$checkbox   
    if(check == TRUE){
      cormat <- cormat_no_outliers
    } else {
      cormat <- full_cormat
    }

    target <- input$variable
    if(target %in% genetree_all_dmel_id$name){
      target_genetree <- genetree_all_dmel_id %>% 
        filter(name == target) %>% 
        pull(genetree) %>% unique
    }
  
    if(any(target_genetree %in% rownames(cormat))){
  
      ovary_changes_cor <- left_join(ovary_changes,
        data.frame(key = names(cormat[target_genetree,]),
          correlation = cormat[target_genetree,]),by="key") %>% 
        na.omit
      correlation_threshold <- 0.825
      target_changes <- ovary_changes_cor %>% filter(abs(correlation) > correlation_threshold)
      target_cormat_head <- target_changes %>% filter(name == target) %>% 
        select(name,child_node,scaled_change) %>% 
        spread(.,name,scaled_change)  %>% 
        select(-child_node)
      target_cormat_tail <- target_changes %>% filter(name != target) %>% 
        select(name,child_node,scaled_change) %>% 
        spread(.,name,scaled_change)  %>% 
        select(-child_node)

      if(length(target_cormat_tail) > 0){
        target_cormat <- cbind(target_cormat_tail,target_cormat_head) %>% 
          cor(.,method="pearson",use="pairwise.complete.obs")
        plot_cor <- ggcorrplot(target_cormat, hc.order = FALSE,type="upper",show.legend=F) + scale_fill_viridis(option="B",limits=c(-1,1),name="correlation")
      } else {
        text = paste("\n   No strong correlations found")
        plot_cor <- ggplot() + 
          annotate("text", x = 0, y = 0, size=8, label = text) + 
          theme_void()
      }
      plot_cor

    } else {
      text = paste("\n   gene family not present in dataset")
      ggplot() + annotate("text", x = 0, y = 0, size=6, label = text) + 
        theme_void()      
    }
  })
}

shinyApp(ui, server)