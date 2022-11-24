library(optparse)
library(tidyverse)
library(data.table)
library(dplyr)
library(stringr)
library(ShortRead)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


source("RIO.R")

input <- function(infile) {
         pfix = prefix()
   parameters <<- readParameters(infile)
   sample_name <<- parameters['sample_name', 2]
   fa1 <<- paste(pfix, parameters['fa1', 2], sep="/")
   fa2 <<- paste(pfix, parameters['fa2', 2], sep="/")
   microbiome_output_file <<- paste(pfix, parameters['microbiome_output_file', 2], sep="/")
   cb_len <<- as.numeric(parameters['cb_len', 2])
   umi_len <<- as.numeric(parameters['umi_len', 2])
   kraken_report <<- paste(pfix, parameters['kraken_report', 2], sep="/")
   mpa_report <<- paste(pfix, parameters['mpa_report', 2], sep="/")
   host <<- as.numeric(parameters['host', 2])
   min_frac <<- as.numeric(parameters['min_frac', 2])
   kmer_len <<- as.numeric(parameters['kmer_len', 2])
   nsample <<- as.numeric(parameters['nsample', 2])

   ranks <<- c('G', 'S')
}

run <- function() {}

output <- function(out_path) {


reads1 = readFasta(fa1)
sequences1 = sread(reads1) 
reads2 = readFasta(fa2)
sequences2 = sread(reads2) 

headers = ShortRead::id(reads1)
barcode = substr(sequences1, 1, cb_len)
umi = substr(sequences1, cb_len+1, cb_len + umi_len)
taxid = gsub('.*taxid\\|', '', headers)
id = gsub('\\s.*', '', headers)

headers2 = ShortRead::id(reads2)
id2 = gsub('\\s.*', '', headers2)

i = intersect(id, id2)
ii = which(id %in% i)
barcode = barcode[ii]
umi = umi[ii]
taxid = taxid[ii]
id = id[ii]
sequences1 = sequences1[ii]

ii = which(id2 %in% i)
sequences2 = sequences2[ii]

kr = read.delim(kraken_report, header = F)
kr = kr[-c(1:2), ] %>% mutate(V8 = trimws(V8)) %>% mutate(V8 = str_replace_all(V8, '[^[:alnum:]]', '_'))
mpa = read.delim(mpa_report, header = F)
mpa$taxid = NA

for(i in 2:nrow(mpa)){
  t_names = mpa[i,1] %>% as.character() %>% 
    strsplit('\\|') %>% 
    unlist() %>% 
    str_remove('.*__') %>% 
    str_replace_all('[^[:alnum:]]', '_') 
  mpa$taxid[i] = paste0('*', paste(kr$V7[match(t_names, kr$V8)], collapse = '*'), '*')
}

microbiome_output_file = read.delim(microbiome_output_file, header = F)

microbiome_output_file = microbiome_output_file %>% 
  select(-V1) %>% 
  separate(V3, into = c('name', 'taxid'), sep = '\\(taxid') %>% 
  mutate(taxid = str_remove(taxid, '\\)') %>% trimws(),
         name = trimws(name))

tx = kr$V7[kr$V6 %in% ranks] %>% setdiff(host)
tx = microbiome_output_file$taxid[microbiome_output_file$taxid %in% tx] %>% unique()


barcode_kmer = list()
counter = 0 
for(taxa in tx){
  counter = counter + 1
  cat(paste('\r', 'taxa processed:', round(counter/length(tx)*100, 3), '%'))
  
  lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
    str_extract(paste0('\\*', taxa, '\\*.*')) %>%
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() %>% 
    unique()
  
  full.lin = str_subset(mpa$taxid, paste0('\\*', taxa, '\\*')) %>% 
    str_remove('^\\*') %>% 
    str_remove('\\*$') %>% 
    str_split('\\*') %>% 
    unlist() %>% 
    as.numeric() %>% 
    unique()
  
  out = subset(microbiome_output_file, taxid %in% lin) %>% separate(V5, into = c('r1', 'r2'), sep = '\\|\\:\\|') 
  out$r1[str_which(out$r1, paste0(' ', host, ':'))] = NA
  out$r2[str_which(out$r2, paste0(' ', host, ':'))] = NA
  out$r1[out$r1 == ''] = NA; out$r2[out$r2 == ''] = NA
  out = subset(out, !is.na(r1) | !is.na(r2))
  out$r1 = trimws(out$r1)
  out$r2 = trimws(out$r2)
  
  
  if(nrow(out) == 0){next}
  
  i = which(id %in% out$V2)
  seq = data.frame(r1 = sequences1[i] %>% as.character(), r2 = sequences2[i] %>% as.character())
  barcode.x = barcode[i]
  umi.x = umi[i]
  
  tax.df = list()
  if(nrow(out) > nsample){n = sample(nrow(out), nsample)} else {n=1:nrow(out)}
  counter2 = 0
  
  for(i in n){
    # cat(paste('\r', 'barcodes processed:', round(counter2/length(n)*100, 3), '%   '))
    for(mate in c('r1', 'r2')){
      
      r = data.frame(pos = out[[mate]][i] %>% strsplit('\\s') %>% unlist()) %>% 
        separate(pos, into = c('taxid', 'nkmer'), sep = ':', convert = T) 
      # if(any(r$taxid %in% c(0, full.lin) == F)){next}
      
      r = r %>% 
        mutate(fkmer = nkmer/sum(nkmer),
               nt_start = cumsum(nkmer) - nkmer + 1,
               nt_end = cumsum(nkmer) + kmer_len - 1) %>% 
        mutate(nt_len = nt_end - nt_start + 1) %>% 
        ungroup() %>% 
        # subset(taxid == taxa)
        subset(taxid %in% c(0,full.lin))
      
      if(sum(r$fkmer) < min_frac){next}
      
      counter2 = counter2 + 1
      
      if(nrow(r) > 0){
        kmer = c()
        for(k in 1:nrow(r)){
          for(m in 1:r$nkmer[k]){
            kmer = c(kmer, substr(seq[[mate]][i], r$nt_start[k] + m - 1, r$nt_start[k] + m + kmer_len - 2))
          }
        }
        tax.df[[counter2]] = data.frame(barcode = barcode.x[i], taxid = taxa, k = kmer, n = sum(r$nt_len[r$taxid %in% lin]))
      } else {tax.df[[counter2]] = data.frame(barcode = barcode.x[i], taxid = taxa, k = NA, n = NA)}
    }
  }
  
  if(length(tax.df) == 0){next}
  
  tax.df = bind_rows(tax.df) %>%
    tibble() %>% 
    subset(!is.na(k)) %>%
    group_by(barcode, taxid) %>%
    summarize(kmer = length(k),
              uniq = length(unique(k)),
              .groups = 'keep')
  
  barcode_kmer[[counter]] = tax.df
  # cat('\n')
}

barcode_kmer = rbindlist(barcode_kmer)

# c = barcode_kmer %>%
#   subset(kmer > 1) %>%
#   group_by(taxid) %>%
#   mutate(nn = n()) %>%
#   subset(nn > 2) %>%
#   group_by(taxid) %>%
#   summarize(r = cor.test(kmer, uniq, method = 'spearman', use = 'pairwise.complete')$estimate,
#             p = cor.test(kmer, uniq, method = 'spearman', use = 'pairwise.complete')$p.value) %>%
#   mutate(p = p.adjust(p)) 

write.table(barcode_kmer, file = paste0(out_path, sample_name, '.sckmer.txt'), quote = F, row.names = F)

paste('Done')


}
