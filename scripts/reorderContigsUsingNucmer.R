#!/usr/bin/Rscript --vanilla
# Author: Sean McIlwain
date_tag="2013_08_28"


#This script attempts to reorder contigs that have been aligned against a reference using nucmer/mummer
#To use, pass in the nucmer.coords and the fasta of the original contigs.
#unplaced contigs to be placed in a separate file and the decisions made will be outputted in a csv file.
#placed contigs will be in separate files with their implied ordering.  Contigs will also be reverse complemented if
#more than 50% of the aligned contig is in reverse order.

#Will need Biostrings library.
## try http if https is not available
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")


if (FALSE) { # Let's make this standalone
getDir<-function(name, dir_in="") {

  dir=dir_in;
  while (!(name %in% list.dirs(dir, full.names=FALSE, recursive = FALSE))) {
    if ("git" %in% list.dirs(dir, full.names=FALSE, recursive = FALSE)) {
       stop("src directory not found!");
    }
    dir = paste0("../",dir);
    
  }
  dir = paste0(dir, "/", name, "/");
  cat(dir, "\n");

  return(dir);
  

}

#############
# Gets the executable directory
############
getExecDir <-function() {
    temp = commandArgs(trailingOnly=FALSE);
    file_string = temp[grep("--file", temp)];
    
    if (length(temp) == 0) {
        stop("Cannot find executable path\n");
    }
    
    temp2 = strsplit(file_string, "=")[[1]][2];
    return(temp2);
    
}

getCommonDir<-function(dir_in="") {
  return(getDir("common_functions", dir_in));
}
getCommonSrcDir<-function(dir_in="") {
  return(paste0(getCommonDir(dir_in), "src/"));
}

getSrcDir <- function(dir_in="") {
  return(getDir("src", dir_in));
}

getDataDir <- function(dir_in="") {
  return(getDir("data", dir_in));
}

src_dir = getSrcDir();
data_dir = getDataDir();
common_src_dir = getCommonSrcDir();


# My custom functions
source(paste0(common_src_dir, "miscFunctions.R"));

source(paste0(common_src_dir, "nucmer.R"));



} else {

  
nucmer=list();
nucmer$readCoords <-function(alignment_path, ref_path=NA, contig_path=NA) {

    library(Biostrings)

    #Read reference sequences
    if (!is.na(ref_path)) {
      ref_seqs = readDNAStringSet(ref_path);
    }
    #Read contig sequences
    if (!is.na(contig_path)) {
      contig_seqs = readDNAStringSet(contig_path);
    }

    #Read alignments from mumer
    temp = read.table(alignment_path, skip=5, sep = "|",stringsAsFactors=FALSE);


    alignments = data.frame(
        ref.name = rep("", nrow(temp)),
        ref.start = rep(-1, nrow(temp)),
        ref.stop = rep(-1, nrow(temp)),
        ref.length = rep(-1, nrow(temp)),
        contig.name = rep("", nrow(temp)),
        contig.start = rep(-1, nrow(temp)),
        contig.stop = rep(-1, nrow(temp)),
        contig.length = rep(NA, nrow(temp)),
        percent.id = rep(0, nrow(temp)),
        contig.full.length = rep(0, nrow(temp)),
        contig.coverage = rep(0, nrow(temp)),
        combined = rep(0, nrow(temp)),
        stringsAsFactors = FALSE
        );


    for (row_idx in 1:dim(temp)[1]) {
        tempsplit = strsplit(temp[row_idx,1], " ")[[1]];
        start_stop = tempsplit[tempsplit != ""];
        alignments$ref.start[row_idx] = as.integer(start_stop[1]);
        alignments$ref.stop[row_idx] = as.integer(start_stop[2]);
        
        tempsplit = strsplit(temp[row_idx,2], " ")[[1]];
        start_stop = tempsplit[tempsplit != ""];
        alignments$contig.start[row_idx] = as.integer(start_stop[1]);
        alignments$contig.stop[row_idx] = as.integer(start_stop[2]);
    

    
        tempsplit = strsplit(temp[row_idx,3], " ")[[1]];
        lengths = tempsplit[tempsplit != ""];
        alignments$ref.length[row_idx] = as.integer(lengths[1]);
        alignments$contig.length[row_idx] = as.integer(lengths[2]);

        alignments$percent.id[row_idx] = as.numeric(temp[row_idx,4]);
    
        tags = strsplit(temp[row_idx,5], "\t")[[1]];
        alignments$ref.name[row_idx] = substr(tags[1],2,nchar(tags[1])); 
        alignments$contig.name[row_idx] = tags[2];        
        if (!is.na(contig_path)) {
	    #print(names(contig_seqs));
            contig_full_length = nchar(as.character(contig_seqs[names(contig_seqs)==alignments$contig.name[row_idx]]));
	    alignments$contig.full.length[row_idx] = contig_full_length;
            alignments$contig.coverage[row_idx] = as.numeric(alignments$contig.length[row_idx]) / contig_full_length * 100;
        }
    }
    if (is.na(contig_path)) {
       #estimate the contig length by the alignments
       contigs = unique(alignments$contig.name);
       for (contig in contigs) {
         idx = alignments$contig.name == contig;
         max.length = max(c(alignments$contig.stop[idx], alignments$contig.start[idx]));
         alignments$contig.full.length[idx] = max.length;
       }
       alignments$contig.coverage = alignments$contig.length / alignments$contig.full.length * 100;
    
    }

    return(alignments);
}

getArguments<-function(arguments_list, option_list=c(), option_defaults=c()) {

    ans = list();
    cargs = c();
    for (idx in 1:length(arguments_list)) {
        ans[arguments_list[idx]] = NA;
    }

    for (idx in 1:length(option_defaults)) {
    	ans[option_list[idx]] = option_defaults[idx];
    }

    args = commandArgs(trailingOnly=TRUE);    
    idx = 1;
    
    while (idx <= length(args)) {
        found_option = FALSE;
        if (length(option_list) > 0) {
            for (option in option_list) {
                if (paste0("--", option) == args[idx]) {
                    ans[option] = args[idx+1];
                    idx = idx + 1;
                    found_option = TRUE;
                    break;
                }
            }
        }
        
        if (!found_option) {
            cargs = c(cargs, args[idx]);
        }
        idx = idx + 1;
    }
    
    if (length(cargs) != length(arguments_list)) {
        cat("Command expects ",length(arguments_list), " " , length(cargs), " were provided\n");
        for (idx in 1:length(cargs)) {
          cat(idx,":",cargs[idx], "\n");
        }
        stop(usage());
    }
    if (length(arguments_list) > 0) {}
    for (idx in 1:length(arguments_list)) {
        ans[arguments_list[idx]] = cargs[idx];
    }
    
    return(ans);
    
}


}

########################################
# FUNCTIONS AND PROCEDURES
########################################



##############################
# MAIN SCRIPT
##############################
#Clean result directories
#cleanup();

args = getArguments(c("nucmer","contigs"));


contig_fasta_path = args$contigs;


coords = nucmer$readCoords(args$nucmer, contig_path=args$contigs);

contig.names = unique(coords$contig.name);
coords$rc = FALSE;
coords$rc[coords$contig.start > coords$contig.stop] = TRUE;
coords$max.coverage = 0;

coords$status = "unplaced";

ref.names = unique(coords$ref.name);
#Step 1, find contigs that uniquely map

unsorted_contig_assignments = list();


for (contig.name in contig.names) {
  ref.maps = table(coords$ref.name[coords$contig.name == contig.name]);
  contig.indices = coords$contig.name == contig.name;
  if(length(ref.maps) == 1) {
    cat("Found unique contig\n");
    ref.name = names(ref.maps)[1];
    coords$status[contig.indices] = ref.name;
    coords$max.coverage[contig.indices] = sum(coords$contig.coverage[contig.indices]);
    if (ref.name %in% names(unsorted_contig_assignments)) {
      unsorted_contig_assignments[[ref.name]] = rbind(unsorted_contig_assignments[[ref.name]], coords[coords$contig.name == contig.name,]);
    } else {
      unsorted_contig_assignments[[ref.name]] = coords[coords$contig.name == contig.name,];
    }  
  } else {
    #For each reference, what is the combined contig coverage?
    combined_refs = list();
    max_coverage = vector();
    for (ref.name in names(ref.maps)) {
      combined_refs[[ref.name]] = coords[contig.indices & coords$ref.name == ref.name,];
      max_coverage[[ref.name]] = sum(combined_refs[[ref.name]]$contig.coverage);
    }
    norm_max_coverage = max_coverage / sum(max_coverage);
    
    best.assignment = names(max_coverage)[which(max_coverage == max(max_coverage))];
    if (length(best.assignment) == 1) {
      coords$status[contig.indices] = best.assignment;
      coords$max.coverage[contig.indices] = max(max_coverage);
      if (ref.name %in% names(unsorted_contig_assignments)) {
        unsorted_contig_assignments[[best.assignment]] = rbind(unsorted_contig_assignments[[best.assignment]], coords[contig.indices & coords$ref.name == best.assignment,]);
      } else {
        unsorted_contig_assignments[[best.assignment]] = coords[contig.indices & coords$ref.name == best.assignment,];
      }
      #stop();
    }
    
  }
}

sorted_contig_assignments = unsorted_contig_assignments;
coords$ref.order = "";
coords$contig.rc = "";
contig_order = list();
contig_rc = c();

for (ref.name in names(sorted_contig_assignments)) {

  sorted_contig_assignments[[ref.name]] = sorted_contig_assignments[[ref.name]][order(sorted_contig_assignments[[ref.name]]$ref.start),];
  contig_order[[ref.name]] = sorted_contig_assignments[[ref.name]]$contig.name[1];
  
  coords$ref.order[coords$contig.name == sorted_contig_assignments[[ref.name]]$contig.name[1]] = 1;
  contigs = sorted_contig_assignments[[ref.name]]$contig.name;
  if (length(contigs) > 1) {
    for (idx in 2:length(contigs)) {
      if (!(contigs[idx] %in% contig_order[[ref.name]])) {
        contig_order[[ref.name]] = c(contig_order[[ref.name]], contigs[idx]);
	coords$ref.order[coords$contig.name == contigs[idx]] = idx;
      }
    }
  }
  
  all_entries =sorted_contig_assignments[[ref.name]] 
  for (contig in contig_order[[ref.name]]) {
    entries = all_entries[all_entries$contig.name == contig,];
    prc = sum(entries$contig.coverage[entries$rc]); #If more than 50% of the mapped contig say to rc, then rc
    if (prc > 50) {
      contig_rc = c(contig_rc, contig);
      coords$contig.rc[coords$contig.name == contig] = "rc";
    }
  }

}
cat("========================\n");
print(sorted_contig_assignments)
cat("========================\n");
print(contig_order);

old.contigs = readDNAStringSet(args$contigs);


for (ref.name in names(contig_order)) {
  
  current_contig_order = contig_order[[ref.name]];
  cat("ref.name:", ref.name,"\n");
  new.bss = DNAStringSet();
  for (idx in 1:length(current_contig_order)) {
    current.contig = current_contig_order[idx];
    current.rc = current.contig %in% contig_rc;
    cat("  current.contig:", current.contig, "\n");
    cat("  rc:", current.rc, "\n");
    current.seq = old.contigs[names(old.contigs)==current.contig];
    if (current.rc) {
      current.seq = reverseComplement(current.seq);
    }
    new.bss = c(new.bss, current.seq);
  }
  writeXStringSet(new.bss, paste0(ref.name, ".ordered.contigs.fasta"));
  cat("write\n");


}
unplaced.contig.names = unique(coords$contig.name[coords$status=="unplaced"]);

if (length(unplaced.contig.names)) {

  unplaced.bss = old.contigs[names(old.contigs) %in% unplaced.contig.names];
  writeXStringSet(unplaced.bss, paste0("unplaced.contigs.fasta"));
}

write.csv(coords, "remap.coords.result.csv", row.names=FALSE);





