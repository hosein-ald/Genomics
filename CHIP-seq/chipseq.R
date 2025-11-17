library(GenomicRanges)
library(rtracklayer)
library(here)
library(Vennerable)
library(readr)
library(ggpmisc)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)

TAZ_peaks=import(here("/Users/hossein.allahdadi/Downloads/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks=import(here("/Users/hossein.allahdadi/Downloads/data/fastq/YAP_peak/YAP_peaks.narrowPeak"))


#overlaps
TAZ_overlap_YAP_peaks = subsetByOverlaps(TAZ_peaks, YAP_peaks)
YAP_overlap_TAZ_peaks = subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(TAZ_overlap_YAP_peaks)
n_overlap=length(YAP_overlap_TAZ_peaks)


n_YAP <- length(YAP_peaks)  # Total peaks 
n_TAZ <- length(TAZ_peaks)  # Total peaks 


venn_data <- Venn(SetNames = c("YAP", "TAZ"),
                  Weight = c(
                    "10" = n_YAP, # Unique to A
                    "01" = n_TAZ, # Unique to B
                    "11" = n_overlap         # Intersection
                  ))
plot(venn_data)




#TEAD4 vs YAP-TAZ overlaps
TEAD4_peak<- import(here("/Users/hossein.allahdadi/Downloads/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)

n_YAP_TAZ <- length(YAP_overlap_TAZ_peaks)  # Total peaks 
n_TEAD4 <- length(TEAD4_peak)  # Total peaks 
n_overlap2<- length(YAP_overlap_TAZ_peaks_overlap_TEAD4)

venn_data2 <- Venn(SetNames = c("YAP/TAZ", "TEAD4"),
                   Weight = c(
                     "10" = n_YAP_TAZ, # Unique to A
                     "01" = n_TEAD4, # Unique to B
                     "11" = n_overlap2        # Intersection
                   ))

plot(venn_data2)




#Summit plot
TAZ_summit<- import(here("/Users/hossein.allahdadi/Downloads/data/fastq/TAZ_peak/TAZ_summits.bed"))
TAZ_summit<- TAZ_summit[TAZ_summit$name %in% TAZ_overlap_YAP_peaks$name]
TEAD4_summit<- import(here("/Users/hossein.allahdadi/Downloads/data/fastq/TEAD4_peak/TEAD4_summits.bed"))
TEAD4_summit

#expand the TAZ summit to a 500bp window
TAZ_500bp_window<- resize(TAZ_summit, width = 500, fix="center")
hits<- findOverlaps(TEAD4_summit, TAZ_500bp_window)
# a hits object with the indices of the overlapping query and subject
hits



summit_distance<- distance(TEAD4_summit[queryHits(hits)], TAZ_summit[subjectHits(hits)])
table(summit_distance)

TEAD4_summit[queryHits(hits)][summit_distance ==0]
TAZ_summit[subjectHits(hits)][summit_distance ==0]


# Compute signed distances
signed_distance <- function(A, B) {
  # Compute unsigned distance
  dist <- distance(A, B)

  # Determine signs based on whether A precedes or follows B
  sign <- ifelse(start(A) < start(B), -1, 1)
  
  # Apply sign to distance
  dist * sign
}

library(dplyr)
library(ggplot2)
summit_distance<- signed_distance(TEAD4_summit[queryHits(hits)],
                                  TAZ_summit[subjectHits(hits)])

distance_df<- table(summit_distance) %>%
  tibble::as_tibble() 

distance_df


distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y = n)) +
  geom_point()

#line graph
distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y = n)) +
  geom_line()



#FIGURE 1f
export(YAP_overlap_TAZ_peaks_overlap_TEAD4,con = here("/Users/hossein.allahdadi/Downloads/data/fastq/YAP_TAZ_TEAD4_common.bed"))

counts<- read_tsv(here("/Users/hossein.allahdadi/Downloads/data/fastq/YAP_TAZ_TEAD4_counts.tsv"), col_names = FALSE)
colnames(counts)<- c("chr", "start", "end", "name", "score", "value", "YAP1", "TAZ", "TEAD4")
View(counts)[1:100,1:9]
head(counts)

#normalize the counts to CPM (counts per million).
counts<- counts %>%
  mutate(YAP1 = YAP1/23653961 * 10^6,
         TAZ = TAZ/26789648 * 10^6,
         TEAD4 = TEAD4/34332907 * 10^6)
head(counts)

ggplot(counts, aes(x=TEAD4, y= YAP1)) +
  geom_point()

#see the outlier
counts %>%
  filter(TEAD4 > 60)

#We can remove that outlier, or use log2 scale
ggplot(counts, aes(x=TEAD4, y= YAP1)) +
  geom_point(color = "#ff4000") +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("YAP1 signal")


#adding R2 to the previous plot, YAP1_TEAD4
ggplot(counts, aes(x=TEAD4, y= YAP1)) +
  geom_point(color = "#ff4000") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
  stat_poly_eq(
    aes(label = ..rr.label..),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("YAP1 signal")



correlation_coefficent<- cor(log2(counts$TEAD4), log2(counts$YAP1))
correlation_coefficent

R_squared<- correlation_coefficent^2
R_squared


#adding R2 to the previous plot, TAZ_TEAD4
ggplot(counts, aes(x=TEAD4, y= TAZ)) +
  geom_point(color = "#ff4000") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
  stat_poly_eq(
    aes(label = ..rr.label..),
    formula = y ~ x,
    parse = TRUE,
    color = "black"
  ) +
  scale_x_continuous(trans = 'log2') +
  scale_y_continuous(trans = 'log2') +
  theme_classic(base_size = 14) +
  xlab("TEAD4 signal") +
  ylab("TAZ signal")


#########
#Drawing the  stacked bar plot. It shows the proportion of the peaks grouped by 
#their distance to the closest TSS 
#########

# Get the TSS
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
head(hg38_transcripts)
# get the TSS.
tss_gr <- promoters(hg38_transcripts, upstream=0, downstream=1)
head(tss_gr)

# Calculate the distance to the nearest TSS for each
distance_to_tss= distanceToNearest(YAP_peaks, tss_gr)
mcols(distance_to_tss)
YAP_dist<- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist<- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist<- mcols(distanceToNearest(TEAD4_peak, tss_gr))$distance

#put them in a single dataframe
tss_distance_df<- bind_rows(data.frame(factor = "YAP", distance = YAP_dist),
                            data.frame(factor = "TAZ", distance = TAZ_dist),
                            data.frame(factor = "TEAD4", distance = TEAD4_dist))

head(tss_distance_df)

#categorise the based on the distances
counts_per_category<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  group_by(factor, category) %>%
  count()

counts_per_category




total_counts<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  count(factor, name = "total")

total_counts


merged_df<- left_join(counts_per_category, total_counts)
merged_df %>%
  mutate(Percentage = n/total * 100) %>%
  ggplot(aes(x= factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distance to TSS",
    x = "Group",
    y = "Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_classic(base_size = 14)



library(rtracklayer) # for reading in bed file

TAZ_peaks<- import(here("data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks<- import(here("data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peak<- import(here("data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

YAP_overlap_TAZ_peaks<- subsetByOverlaps(YAP_peaks, TAZ_peaks)

YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)
YAP_overlap_TAZ_peaks_overlap_TEAD4

H3K4me1<- import(here("/Users/hossein.allahdadi/Downloads/data/GSE49561/H3K4me1.bed"))
H3K4me3<- import(here("/Users/hossein.allahdadi/Downloads/data/GSE49561/H3K4me3.bed"))
H3K27ac<- import(here("/Users/hossein.allahdadi/Downloads/data/GSE49561/H3K27Ac.bed"))

active_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac)
inactive_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)



n_active_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                      active_enhancers) %>% length()

n_inactive_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                        inactive_enhancers) %>% length()

n_promoters<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
                               promoters) %>% length()

n_unclassified<- length(YAP_overlap_TAZ_peaks_overlap_TEAD4) - n_active_enhancers -
  n_inactive_enhancers - n_promoters


#Put the numbers in a dataframe:
  
annotation_df<- data.frame(category = c("active_enhancers", "inactive_enhancers",
                                          "promoters", "unclassified"),
                             peak_number = c(n_active_enhancers, n_inactive_enhancers, 
                                             n_promoters, n_unclassified))


View(annotation_df)

#plotting
ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_brewer(palette = "Set3") 


#new plot
annotation_df$category<- factor(annotation_df$category, 
                                levels = c("promoters", "active_enhancers",
                                           "inactive_enhancers", "unclassified"))

colors<- c("#8D1E0F", "#F57D2B", "#FADAC4", "#D4DADA")

ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors)




#Add the percentage to the pie chart
# Calculate percentages and cumulative positions for labeling
annotation_df <- annotation_df %>%
  dplyr::mutate(
    percentage = peak_number / sum(peak_number) * 100,
    label = paste0(round(percentage, 1), "%")
  )

annotation_df

# Create the pie chart
ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) # Add percentage labels




