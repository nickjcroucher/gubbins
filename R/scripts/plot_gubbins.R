#!/usr/bin/env Rscript

#############
# Libraries #
#############

suppressPackageStartupMessages({
  require(argparser)
  require(magrittr)
  require(tidyverse)
  require(treeio)
  require(ggtree)
  require(aplot)
  require(patchwork)
  require(cowplot)
  require(RColorBrewer)
}
)

##############################
# Functions: Processing tree #
##############################

process_gubbins_tree <- function(tree_fn) {
  
  # Parse file
  gubbins_tree <- treeio::read.tree(tree_fn)

  return(gubbins_tree)
  
}

generate_gubbins_ggtree <- function(gubbins_tree, show_taxa = FALSE, bl_threshold = Inf, truncated_col = "red") {
  
  # Colour for truncated branches
  branch_palette <-
    c("Truncated" = truncated_col,
      "Normal" = "black")
  
  
  # Truncate long branches
  gubbins_tree$edge.length[gubbins_tree$edge.length > bl_threshold] <- bl_threshold
  
  # Read in ggtree
  gubbins_ggtree <- 
    ggtree(gubbins_tree, size = 0.25)
  
  # Label truncated branches
  gubbins_ggtree$data %<>%
    dplyr::mutate(truncated = dplyr::if_else(branch.length>=bl_threshold,
                                             "Truncated",
                                             "Normal"))
  # Scale tree
  gubbins_ggtree <-
    gubbins_ggtree +
    aes(colour = truncated) +
    scale_colour_manual(values = branch_palette) +
    theme_tree2() +
    theme(legend.position = "none",
          axis.line.x = element_line(linewidth = 0.25))
  
  # Add taxon labels
  if (show_taxa) {
    gubbins_ggtree <-
      gubbins_ggtree + 
      geom_tiplab()
  }
  
  return(gubbins_ggtree)
}

##################################
# Functions: Processing metadata #
##################################

make_palette <- function(i,df=meta.df) {
  
  unique_values <- unique(df[,i])
  unlimited_palettes <- c("Blues", "Greens", "Oranges", "Purples", "Greys", "Reds")
  continuous_palettes <-c("RdYlGn", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu")
  discrete_palettes <- c("Accent", "Dark2", "Paired", "Set1", "Set2", "Set3", "Spectral", "Dark1")
  if (class(df[,i]) == "numeric") {
    if (i <= length(continuous_palettes)) {
      palette_object <- colorRampPalette(brewer.pal(5,continuous_palettes[i]))
    } else {
      palette_object <- colorRampPalette(brewer.pal(5,unlimited_palettes[i-length(continuous_palettes)]))
    }
  } else {
    if (i <= length(discrete_palettes)) {
      palette_object <- RColorBrewer::brewer.pal(max(3,length(unique_values)),discrete_palettes[i])
    } else {
      palette_object <- RColorBrewer::brewer.pal(max(3,length(unique_values)),unlimited_palettes[i-length(discrete_palettes)])
    }
    palette_object <- palette_object[1:length(unique_values)]
    names(palette_object) <- unique_values
  }
  return(palette_object)
  
}

make_individual_heatmap <- function(i, df = meta.df, tree_plot = gubbins_tree, legend_direction = "horizontal", ncol = 1, label_size = NULL) {
  
  # Get legend
  individual_palette <- make_palette(i, df = df)
  
  # Make title grob
  title_text <- as.expression(parse(text=colnames(df)[i]))

  # Make plot
  individual_heatmap <-
    ggplot(data = data.frame(
              x = 1,
              y = factor(rownames(df),
                         levels = rev(ggtree::get_taxa_name(tree_plot))),
              vals = df[,i]
          ),
          aes(
            x = x,
            y = y,
            colour = vals,
            fill= vals
          )
    ) +
    geom_tile(height = 1)
  
  label_y_position <- ggplot_build(individual_heatmap)$layout$panel_params[[1]]$y.range[2]
  individual_heatmap <-
    individual_heatmap +
    geom_text(aes(x = 1, y = label_y_position),
              hjust = 0.0,
              angle = 90,
              nudge_y = 0.0,
              inherit.aes = FALSE,
              label = colnames(df)[i],
              parse = TRUE,
              size = label_size
    )
    
  individual_heatmap <-
    individual_heatmap +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) + 
    coord_cartesian(clip = "off")
  
    # Add palette
    if (class(df[,i]) == "numeric") {
      individual_heatmap <-
        individual_heatmap +
        scale_colour_gradientn(colours = individual_palette(100),
                               aesthetics = c("fill","colour"),
                               name = title_text) +
        guides(colour = guide_colorbar(title.hjust=0.5, ncol = ncol, direction = legend_direction, title.position = "top"),
               fill = guide_colorbar(title.hjust=0.5, ncol = ncol, direction = legend_direction, title.position = "top"))
    } else {
      individual_heatmap <-
        individual_heatmap +
        scale_colour_manual(values = individual_palette,
                              aesthetics = c("fill","colour"),
                              name = title_text) +
        guides(colour = guide_legend(title.hjust=0.5, ncol = ncol, direction = legend_direction, title.position = "top"),
                 fill = guide_legend(title.hjust=0.5, ncol = ncol, direction = legend_direction, title.position = "top"))
    }
  
  # Return plot
  return(individual_heatmap)
  
}

return_metadata_plot <- function(metadata_fn, gubbins_tree, legend_direction = NULL, meta_label_size = NULL) {
  
  # Read file
  meta.df <-
    read.csv(metadata_fn, row.names = 1, check.names = FALSE)
  
  # Optimise direction
  if (is.null(legend_direction)) {
    if (ncol(meta.df) > 2) {
      legend_direction = "vertical"
    } else {
      legend_direction = "horizontal"
    }
  }
  
  # Generate individual plots
  individual_heatmaps_with_legends <- lapply(1:ncol(meta.df),make_individual_heatmap, df = meta.df, tree_plot = gubbins_tree, label_size = meta_label_size)
  individual_heatmaps <- lapply(individual_heatmaps_with_legends, function(x) {x + theme(legend.position = "none",plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm"))})
  individual_legends <- lapply(individual_heatmaps_with_legends, function(x) {suppressWarnings(cowplot::get_legend(x))})
  
  # Plot tree and heatmaps
  combined_heatmap <- patchwork::wrap_plots(individual_heatmaps,nrow = 1)
  
  # Return
  return(list(combined_heatmap,individual_legends))

}

####################################
# Functions: Processing annotation #
####################################

trim_start <- function(df,start_pos) {
  
  if (!is.null(start_pos)) {
    df %<>%
      dplyr::filter(end > start_pos) %>%
      dplyr::mutate(start = dplyr::case_when(
        start < start_pos & end > start_pos ~ start_pos,
        TRUE ~ start
      )
      )
  }
  
  return(df)
}

trim_end <- function(df,end_pos) {
  
  if (!is.null(end_pos)) {
    df %<>%
      dplyr::filter(start < end_pos) %>%
      dplyr::mutate(end = dplyr::case_when(
        start < end_pos & end > end_pos ~ end_pos,
        TRUE ~ end
      )
      )
  }
  
  return(df)
  
}

# Read Gubbins recombination GFF file - taken from RCandy to enable conda installation
load.gubbins.GFF<-function(gubbins.gff.file,recom.input.type="Gubbins"){
  # Check if correct input recombination data type is specified
  if( !recom.input.type %in% c("Gubbins","BRATNextGen") ){
    stop("Invalid recombination data specified. Choose from 'Gubbins' or 'BRATNextGen'")
  }
  # Check if the input data was generated by Gubbins (GFF file) or BRATNextGen (tabular file)
  if( recom.input.type=="Gubbins" ){
    # Check if a valid GFF Gubbins recombination file is specified
    if( !is.null(gubbins.gff.file) ){
      # Check if the Gubbins recombination file name is a string or character
      if( is.character(gubbins.gff.file) ){
        # Check if the Gubbins recombination file exists in the file path
        if( file.exists(gubbins.gff.file) ){
          
          # Read the GFF Gubbins recombination file, skips the first two lines which contain comments
          tree.rec.data1<-dplyr::as_tibble(read.table(gubbins.gff.file,header=FALSE,sep="\t",comment.char="?",skip=2,fill=T,row.names=NULL))
          colnames(tree.rec.data1)<-c("SEQ","PROG","TYPE","START","END","XX","YY","ZZ","REC")
          
          # Check if at least one recombination event was found in the Gubbins GFF file
          if( length(tree.rec.data1$SEQ)>1 ){
            # Extract the taxon names in the Gubbins GFF file
            tree.rec.data.tmp<-tree.rec.data1 %>% dplyr::mutate(REC1=.data$REC) %>%
              dplyr::group_by(.data$SEQ,.data$PROG,.data$TYPE,.data$START,.data$END,.data$XX,.data$YY,.data$ZZ,.data$REC) %>%
              tidyr::nest() %>% dplyr::rowwise() %>%
              dplyr::mutate(gene=list(setdiff(stringr::str_split(stringr::str_trim(gsub("taxa=","",gsub("taxa= ","",stringr::str_trim(stringr::str_split(regmatches(data[[1]],regexpr("(taxa=).*;",data[[1]])),";")[[1]],side="both"))),side="both")," ")[[1]],c(""," "))) )
            
            # Return a data frame containg recombination events in appropriate format
            return(tree.rec.data.tmp)
          }else{
            # Exit the program when no valid GFF Gubbins recombination file is found
            stop("It appears that there are no recombination events in the file '",gubbins.gff.file,"'")
          }
        }else{
          # Exit the program when no valid GFF Gubbins recombination file is found
          stop("Cannot find the Gubbins file '",gubbins.gff.file,"' containing recombination events")
        }
      }else{
        if( length(setdiff(class(gubbins.gff.file),c("tbl_df","tbl","data.frame","rowwise_df","grouped_df")))==0 ){
          return(gubbins.gff.file)
        }else{
          # Exit the program when an invalid GFF Gubbins recombination file is specified
          stop("Gubbins recombination events file name does not appear to be a character or string")
        }
      }
    }else{
      # Return nothing if no GFF Gubbins recombination file is provided
      return(NULL)
    }
  }else{
    tree.rec.data.tmp<-as_tibble(read.table(gubbins.gff.file,
                                            fill=TRUE,sep="\t",comment.char="?",skip=2,header=FALSE)) %>%
      dplyr::mutate(V2=.data$V1) %>% dplyr::group_by(.data$V1) %>% tidyr::nest() %>%
      dplyr::ungroup() %>% dplyr::select(-.data$V1) %>% dplyr::rowwise() %>%
      dplyr::mutate(rec.events=list(stringr::str_split(.data$data[[1]]," ")[[1]][!stringr::str_split(.data$data[[1]]," ")[[1]] %in% c("")])) %>%
      dplyr::mutate( Start=.data$rec.events[1],
                     End=.data$rec.events[2],
                     Origin=.data$rec.events[3],
                     HomeCluster=.data$rec.events[4],
                     BAPSIndex=.data$rec.events[5],
                     StrainName=.data$rec.events[6]) %>%
      dplyr::mutate(SEQ="SEQUENCE",
                    PROG="GUBBINS",
                    TYPE="CDS",
                    START=.data$Start,
                    END=.data$End,
                    XX=0.000,
                    YY=".",
                    ZZ=0) %>% dplyr::select(-.data$data,-.data$rec.events,-.data$Start,-.data$End,-.data$Origin,-.data$HomeCluster,-.data$BAPSIndex) %>%
      dplyr::group_by(.data$SEQ,.data$PROG,.data$TYPE,.data$START,.data$END,.data$XX,.data$YY,.data$ZZ) %>% tidyr::nest() %>%
      dplyr::mutate( REC=paste("taxa='",stringr::str_trim(paste0(.data$data[[1]]$StrainName,collapse="..."),side="both"),"';",sep=" ",collapse=" ")) %>%
      dplyr::mutate(REC=gsub("\\.\\.\\."," ",gsub(" ","",.data$REC))) %>% dplyr::mutate(REC=list(.data$REC)) %>%
      dplyr::mutate(gene=list(.data$data[[1]]$StrainName),data=.data$REC)
    return(tree.rec.data.tmp)
  }
}

# Check if a valid reference genome name is provided - taken from RCandy to enable conda installation
load.genome.GFF<-function(reference.genome){
  V1<-seqname<-feature<-score<-strand<-REC<-SEQ<-PROG<-TYPE<-START<-END<-XX<-YY<-ZZ<-NULL
  if( !is.null(reference.genome) & is.character(reference.genome) ){
    # Same coordinates for the genome region to show, default whole genome
    # Read the reference genome GFF annotation file
    tmp.ref.df<-dplyr::as_tibble(read.csv(reference.genome,comment.char="#",header=FALSE,sep="\t",fill=TRUE,row.names=NULL)) %>%
      dplyr::filter((!grepl("#",V1)) | V1!="seqname" )
    colnames(tmp.ref.df)<-c("seqname","source","feature","start","end","score","strand","frame","attributes")
    tmp.ref.df<-tmp.ref.df[!tmp.ref.df$source %in% c("source"),]
    tmp.ref.df[1,3]<-"source"
    reference.genome.obj1<-tmp.ref.df %>%
      dplyr::mutate(seqname=gsub("# ","",.data$seqname)) %>% mutate(seqname=gsub("^#","",.data$seqname)) %>%
      dplyr::filter(!grepl("gff-version",.data$seqname)) %>%
      dplyr::filter(!.data$feature %in% c("ORIGIN","NA","","##") & .data$seqname!="") %>%
      dplyr::mutate(start=as.integer(.data$start),end=as.integer(.data$end)) %>%
      dplyr::filter(.data$feature %in% c("source","locus_tag","gene","CDS"))
    
    # Filter out lines not containing information about the genetic features
    colnames(reference.genome.obj1)<-c("seqname","source","feature","start","end","score","strand","frame","attributes")
    reference.genome.obj<-reference.genome.obj1 %>%
      dplyr::group_by(.data$seqname,.data$source,.data$feature,.data$start,.data$end,.data$score,.data$strand,.data$frame) %>%
      dplyr::filter(!.data$feature %in% c("ORIGIN","NA","","##")) %>%
      dplyr::mutate(attributes=gsub("="," ",.data$attributes)) %>%
      dplyr::mutate(attributes=gsub(":"," ",.data$attributes)) %>% #dplyr::filter(feature %in% c("CDS","gene","source")) %>%
      tidyr::nest() %>% dplyr::mutate(gene=stringr::str_split(stringr::str_trim(stringr::str_split(regmatches(data[[1]],regexpr("(gene|locus_tag|Parent|db_xref|mol_type|organism|ID).*",data[[1]])),";")[[1]][1],side="both")," ")[[1]][-1][1] )
  }else{
    if( length(setdiff(class(reference.genome),c("tbl_df","tbl","data.frame","rowwise_df","grouped_df")))==0 ){
      reference.genome.obj<-reference.genome
      colnames(reference.genome.obj)<-c("seqname","source","feature","start","end","score","strand","frame","attributes")
    }else{
      # Exit the program when valid genome length is found
      stop("Could not find a feature labelled 'source' in the genome annotation file")
    }
  }
  if( !"source" %in% reference.genome.obj$feature ){
    # Exit the program when valid genome length is found
    stop("Could not find a feature labelled 'source' in the genome annotation file")
  }
  return(reference.genome.obj)
}

return_annotation_df <- function(anno_gff_fn, start_pos = NULL, end_pos = NULL, features = c("CDS")) {
  
  gff.df <-
    load.genome.GFF(anno_gff_fn) %>%
      dplyr::ungroup() %>%
      dplyr::filter(feature %in% features) %>%
      dplyr::mutate(strand = dplyr::if_else(strand == "+",1,-1)) %>%
      dplyr::select(start,end,strand,gene)
  
  # Trim
  gff.df %<>%
    trim_start(start_pos) %>%
    trim_end(end_pos)
  
  return(gff.df)
  
}

generate_annotation_plot <- function(gubbins_anno,anno_features = c("CDS")) {
  
  # Get end coordinate
  anno_max <-
    max(gubbins_anno$end)
  
  # Generate plot
  anno_plot <-
    ggplot(gubbins_anno,
           aes(x = start,
               xend = end,
               y = strand,
               yend = strand)) +
    geom_segment(linewidth = 5) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(limits = c(1, anno_max), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1)) +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()
    )
  
  # Return
  return(anno_plot)
}

markup_df_from_annotation <- function(anno_df) {
  
  # Read data
  locus.df <-
    anno_df %>%
      rename(label = gene) %>%
      dplyr::mutate(midpoint = start + 0.5*(end - start + 1))
  
}

markup_df_from_csv <- function(markup_fn, start_pos = NULL, end_pos = NULL) {

  # Read data
  locus.df <-
    read.csv(markup_fn) %>%
    {if("Start" %in% names(.)) dplyr::rename(.,start = Start) else .} %>%
    {if("End" %in% names(.)) dplyr::rename(.,end = End) else .} %>%
    {if("Label" %in% names(.)) dplyr::rename(.,label = Label) else .}
  
  # Trim
  locus.df %<>%
    trim_start(start_pos) %>%
    trim_end(end_pos)
  
  # Add midpoints
  locus.df %<>%
    dplyr::mutate(midpoint = start + 0.5*(end - start + 1),
                  strand = 1)
  
  return(locus.df)
}

generate_markup_plot <- function(markup_df) {

  # Generate plot
  labels_plot <-
    ggplot(markup_df,
           aes(x = start,
               xend = end,
               y = strand,
               yend = strand)) +
    geom_segment() +
    geom_text(aes(x = midpoint,
                  y = strand*1.01,
                  label = label),
              angle = 45,
              hjust = 0,
              vjust = 0,
              size = 3,
              parse=TRUE) +
    theme_bw() +
    #scale_x_continuous(limits=c(1, anno_max), expand = c(0, 0)) +
    scale_y_continuous(limits=c(1,2)) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()
    )
  
  # Return
  return(labels_plot)
  
}

#######################################
# Functions: Processing recombination #
#######################################

process_gubbins_recombination_df <- function(rec_gff_fn) {
  
  # Load and process recombination
  gubbins_rec <-
    load.gubbins.GFF(rec_gff_fn) %>%
      ungroup() %>%
      select(START,END,gene) %>%
      dplyr::rename(start = START) %>%
      dplyr::rename(end = END)

}

plot_gubbins_recombination <- function(gubbins_rec,gubbins_tree, start_pos = NULL, end_pos = NULL) {
  
  # Process data frame
  gubbins_rec %<>%
    dplyr::mutate(Colour = dplyr::if_else(length(gene) == 1,
                                          "blue",
                                          "red")) %>%
    tidyr::unnest_longer(gene,
                         values_to = "Taxa") %>%
    dplyr::mutate(Taxa = factor(Taxa,
                                levels = rev(get_taxa_name(gubbins_tree)))) %>%
    dplyr::mutate(length = end - start + 1) %>%
    dplyr::arrange(rev(length)) %>%
    dplyr::select(-length)
  
  # Trim to visualised region
  gubbins_rec %<>%
    trim_start(start_pos) %>%
    trim_end(end_pos)
  
  # Plot recombination
  rec_plot <-
    ggplot(gubbins_rec,
           aes(x = start,
               xend = end,
               y = Taxa,
               yend = Taxa,
               colour = Colour)) +
    geom_segment(alpha = 0.25) +
    scale_colour_manual(values = c("red" = "red",
                                   "blue" = "blue")) +
    scale_y_discrete(drop = FALSE) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title = element_blank())
  
  return(rec_plot)
  
}

generate_heatmap <- function(rec_df,start_coordinate,end_coordinate,max_val=10) {
  
  # Filter dataframe
  rec_df %<>%
    dplyr::select(start,end) %>%
    dplyr::mutate(bases = purrr::map2(start,end,seq)) %>%
    dplyr::select(bases) %>%
    tidyr::unnest(bases) %>%
    dplyr::filter(bases > start_coordinate & bases < end_coordinate) %>%
    dplyr::group_by(bases) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(count = dplyr::if_else(count > max_val, max_val, count)) %>%
    dplyr::mutate(y = 1)
  
  rec_heatmap_plot <-
    ggplot(rec_df,
           aes(x = bases,
               y = 1,
               colour = count,
               fill = count)) +
      geom_tile() +
      scale_x_continuous(limits = c(start_coordinate-0.5,end_coordinate+0.5),expand = c(0,0)) +
      scale_y_continuous(limits = c(0.5,1.5),expand = c(0,0)) +
      scale_colour_gradient2(low = "navy",
                             mid = "orange",
                             high = "red",
                             limits = c(0,max_val),
                             midpoint = max_val/2,
                             name = "Recombination\nevents",
                             aesthetics = c("fill","colour")) +
      guides(colour = guide_colorbar(title.position = "left",
                                   title.hjust = 0.5,
                                   direction = "horizontal"),
             fill = guide_colorbar(title.position = "left",
                                 title.hjust = 0.5,
                                 direction = "horizontal")) +
      theme(panel.background = element_rect(fill = 'navy'),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            legend.title = element_text(size = 9))
    
  return(rec_heatmap_plot)
}

##############################
# Functions: Combining plots #
##############################

plot_gubbins <- function(tree = NULL,
                         rec = NULL,
                         markup = NULL,
                         anno = NULL,
                         meta = NULL,
                         plot_heatmap = TRUE,
                         start_coordinate = NULL,
                         end_coordinate = NULL,
                         show_taxa = FALSE,
                         max_branch_length = Inf,
                         annotation_labels = FALSE,
                         tree_width = 0.4,
                         meta_width = 0.25,
                         annotation_height = 0.05,
                         markup_height = 0.075,
                         heatmap_height = 0.025,
                         legend_height = 0.2,
                         meta_label_size = 6,
                         heatmap_y_nudge = 0.0,
                         heatmap_x_nudge = 0.0,
                         legend_direction = NULL) {
  
  # Generate individual components
  gubbins_tree_obj <- process_gubbins_tree(tree)
  gubbins_tree <- generate_gubbins_ggtree(gubbins_tree_obj, bl_threshold = max_branch_length, show_taxa = show_taxa)
  gubbins_rec_df <- process_gubbins_recombination_df(rec)
  genome_length <- max(gubbins_rec_df$end)
  
  gubbins_meta <- NULL
  gubbins_anno <- NULL
  if (!is.null(meta)) {
    meta_plot_info <- return_metadata_plot(meta,gubbins_tree,legend_direction = legend_direction, meta_label_size = meta_label_size)
    gubbins_meta <- meta_plot_info[[1]]
    gubbins_legends <- meta_plot_info[[2]]
  }
  if (!is.null(anno)) {
    anno_features <- c("CDS")
    anno_df <- return_annotation_df(anno,
                                    start_pos = start_coordinate,
                                    end_pos = end_coordinate,
                                    features = anno_features)
    gubbins_anno <- generate_annotation_plot(anno_df)
    genome_length <- max(max(anno_df$end),genome_length)
  }
  if (!is.null(markup)) {
    markup_df <- markup_df_from_csv(markup,
                                    start_pos = start_coordinate,
                                    end_pos = end_coordinate)
    gubbins_markup <- generate_markup_plot(markup_df)
    genome_length <- max(max(gubbins_markup$data$end,genome_length))
  } else if (annotation_labels) {
    markup_df <- markup_df_from_annotation(anno_df)
    gubbins_markup <- generate_markup_plot(markup_df)
  }
  
  # Set region to be plotted
  if (is.null(start_coordinate) | is.null(end_coordinate)) {
    start_coordinate <- 1
    end_coordinate <- genome_length
  }
  
  gubbins_rec <- 
    plot_gubbins_recombination(gubbins_rec_df,
                               gubbins_tree,
                               start_pos = start_coordinate,
                               end_pos = end_coordinate)
  gubbins_rec <-
    gubbins_rec + scale_x_continuous(limits=c(start_coordinate, end_coordinate), expand = c(0, 0))
  
  # Calculate heatmap
  if (plot_heatmap) {
    heatmap_plot <- generate_heatmap(gubbins_rec_df,start_coordinate, end_coordinate)
  }
  
  # Combine components
  options("aplot_guides" = "keep")
  if (!is.null(meta)) {
    combined_plot <-
      gubbins_rec %>%
        aplot::insert_left(gubbins_meta, width = meta_width) %>%
        aplot::insert_left(gubbins_tree,width = tree_width)
  } else {
    combined_plot <-
      gubbins_rec %>%
        aplot::insert_left(gubbins_meta, width = tree_width)
  }
  
  # Add heatmap
  if (plot_heatmap) {
    suppressWarnings( # Suppress warnings here because of the missing first and last tiles when subsetting a region
      combined_plot <-
        combined_plot %>%
          aplot::insert_top(heatmap_plot + theme(legend.position = "none"), height = heatmap_height)
    )
  }
  
  # Add annotation to plot
  if (!is.null(anno)) {
    combined_plot <-
      combined_plot %>%
        aplot::insert_top(gubbins_anno, height = annotation_height)
  }

  # Add markup to plot
  if (!is.null(markup) | annotation_labels) {
    combined_plot <-
      combined_plot %>%
        aplot::insert_top(gubbins_markup, height = markup_height)
  }
  
  # Now add legends
  gubbins_legends_plot <- cowplot::plot_grid(plotlist = gubbins_legends,
                                              nrow = 1)
  gubbins_plot_with_legends <-
    cowplot::plot_grid(plotlist = list(aplot::as.patchwork(combined_plot),
                                       gubbins_legends_plot),
                       nrow = 2,
                       rel_heights = c(1,legend_height)) +
    annotation_custom(cowplot::get_legend(heatmap_plot), xmin  = -0.8 + heatmap_x_nudge, ymin = 0.925 + heatmap_y_nudge)
  
  return(gubbins_plot_with_legends)
}

#######################
# Parse commmand line #
#######################

parse_command_line <- function() {
  
  # Define parser
  p <- arg_parser("plot_gubbins.R: Producing publication-ready figures of Gubbins analyses")
  
  # Input files
  p <- add_argument(p,
                    "--tree",
                    "Gubbins tree (Newick file)",
                    default = NULL,
                    type = "character"
                    )
  p <- add_argument(p,
                    "--rec",
                    "Gubbins recombination inference (GFF file)",
                    default = NULL,
                    type = "character"
  )
  p <- add_argument(p,
                    "--annotation",
                    "Reference genome annotation (GFF file)",
                    default = NULL,
                    type = "character"
  )
  p <- add_argument(p,
                    "--markup",
                    "Genome loci to mark (CSV file)",
                    default = NULL,
                    type = "character"
  )
  p <- add_argument(p,
                    "--meta",
                    "Metadata for each sequence (CSV file)",
                    default = NULL,
                    type = "character"
  )
  p <- add_argument(p,
                    "--output",
                    "Output file name (PNG or PDF suffix)",
                    default = NULL,
                    type = "character"
  )
  
  # Relative sizes
  p <- add_argument(p,
                    "--tree-width",
                    "Width of tree relative to recombination panel",
                    default = 0.4,
                    type = "numeric"
  )
  p <- add_argument(p,
                    "--meta-width",
                    "Width of metadata panel relative to recombination panel",
                    default = 0.25,
                    type = "numeric"
  )
  p <- add_argument(p,
                    "--anno-height",
                    "Height of annotation panel relative to recombination panel",
                    default = 0.05,
                    type = "numeric"
  )
  p <- add_argument(p,
                    "--markup-height",
                    "Height of markup panel relative to recombination panel",
                    default = 0.075,
                    type = "numeric"
  )
  p <- add_argument(p,
                    "--heatmap-height",
                    "Height of heatmap relative to recombination panel",
                    default = 0.025,
                    type = "numeric"
  )
  p <- add_argument(p,
                    "--legend-height",
                    "Height of legends relative to recombination panel",
                    default = 0.25,
                    type = "numeric"
  )
  
  # Graph format options
  p <- add_argument(p,
                    "--start-coordinate",
                    "Left boundary of genomic region to plot",
                    default = NULL,
                    type = "integer"
  )
  p <- add_argument(p,
                    "--end-coordinate",
                    "Right boundary of genomic region to plot",
                    default = NULL,
                    type = "integer"
  )
  p <- add_argument(p,
                    "--no-heatmap",
                    "Do not plot recombination heatmap",
                    flag=TRUE
  )
  p <- add_argument(p,
                    "--legend-direction",
                    "Orientation of legends (horizontal or vertical)",
                    default = NULL,
                    type = "character"
  )
  p <- add_argument(p,
                    "--output-height",
                    "Height of output file (inches)",
                    default = 8,
                    type = "numeric"
  )
  p <- add_argument(p,
                    "--output-width",
                    "Width of output file (inches)",
                    default = 11,
                    type = "numeric"
  )
  
  # Return parser
  return(p)
}

evaluate_args <- function(args) {

  args_error <- FALSE
  
  if (is.na(args[["tree"]])) {
    message("Tree file is required")
    args_error <- TRUE
  } else if (!file.exists(args[["tree"]])) {
    message("Tree file does not exist")
    args_error <- TRUE
  }
  
  if (is.na(args[["rec"]])) {
    message("Recombination GFF file is required")
    args_error <- TRUE
  } else if (!file.exists(args[["rec"]])) {
    message("Recombination GFF file does not exist")
    args_error <- TRUE
  }
  
  return(args_error)
}

###############
# Parse input #
###############

parser <- parse_command_line()
args <- parse_args(parser)
args_error <- evaluate_args(args)

if (args_error) {
  print(parser)
  quit()
}

gubbins_plot <-
  plot_gubbins(tree = args[["tree"]],
               rec = args[["rec"]],
               markup = args[["markup"]],
               anno = args[["annotation"]],
               meta = args[["meta"]])

ggsave(file=args[["output"]],
       height = args[["output_height"]],
       width = args[["output_width"]])
