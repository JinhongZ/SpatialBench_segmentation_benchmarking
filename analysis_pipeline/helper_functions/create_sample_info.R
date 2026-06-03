# Purpose:  Create a default sample_info data.table of full file path to 
#           transcript metadata for each segmentation and platform if not 
#           provided

# list file directories for MERSCOPE's Default and Cellpose2 samples
merscope_dir <- "/vast/projects/SpatialBench/segmentation/Cellpose/resegmentation_cellpose_merscope/MERSCOPE_cellpose_reseg"
merscope_sub_dir <- list.files(merscope_dir, pattern = "^Batch|^Sample")
merscope_name <- merscope_sub_dir %>% 
  gsub(pattern = "^Batch\\d+_region\\d_|^Sample-", replacement = "") %>% 
  tolower() %>% 
  gsub(pattern = "^7", replacement = "wt7")
merscope_vizgen <- file.path(merscope_sub_dir, "Vizgen")
merscope_cellpose <- file.path(merscope_sub_dir, "Cellpose")

merscope_vizgen_info <- data.frame(
  sample = merscope_name,
  segmentation_dir = merscope_dir,
  sub_dir = file.path(merscope_vizgen, "detected_transcripts.csv")
)
merscope_cellpose_info <- data.frame(
  sample = merscope_name,
  segmentation_dir = merscope_dir
) %>% 
  mutate(sub_dir = case_when(
    sample %in% c("wt709", "wt713") ~ file.path(merscope_cellpose, "detected_transcripts_resegmented.csv"),
    .default = file.path(merscope_cellpose, "detected_transcripts.csv")
  ))

# list file directories for MERSCOPE's Proseg samples
merscope_dir <- "/vast/projects/SpatialBench/segmentation/Proseg"
merscope_sub_dir <- list.files(merscope_dir, pattern = "^proseg_merscope")
merscope_sample <- c(
  "ctrl172c", "ctrl174", "ko166", "ctrl173", "wt710b", "ko168", 
  "wt713b", "ctrl172", "ko167", "wt709", "wt710", "wt713"
)
merscope_proseg_info <- data.frame(
  sample = merscope_sample,
  segmentation_dir = merscope_dir,
  sub_dir = file.path(merscope_sub_dir, "output/transcript-metadata.csv.gz")
) %>% filter(sample %in% merscope_vizgen_info$sample)

# merge multiple MERSCOPE samples
merscope_info <- data.table::rbindlist(
  Map(function(df, seg) {
    df %>% mutate(platform = "MERSCOPE", model = NA_character_, segmentation = seg)
  }, list(merscope_vizgen_info, merscope_proseg_info, merscope_cellpose_info), 
  c("Default", "Proseg", "Cellpose2"))
)

# list file directories for Xenium unimodal default segmentation
xenium_dir <- "/vast/projects/SpatialBench/data/xenium/G000218_Benchmarking/resegmented_xeniumranger_v2_0_1"
xenium_uni_default <- tribble(
  ~sub_dir,                                 ~sample,     
  "resegmented_v2_0_1__0017329__Region_2", "ctrl172",
  "resegmented_v2_0_1__0011456__Region_4", "ctrl173",
  "resegmented_v2_0_1__0017329__Region_4", "ctrl174",
  "resegmented_v2_0_1__0011456__Region_1", "ko166",
  "resegmented_v2_0_1__0011456__Region_2", "ko167",
  "resegmented_v2_0_1__0017329__Region_3", "ko168",
  "resegmented_v2_0_1__0017329__Region_5", "wt709",
  "resegmented_v2_0_1__0017329__Region_1", "wt710",
  "resegmented_v2_0_1__0011456__Region_3", "wt713"
) %>% 
  mutate(
    segmentation_dir = xenium_dir, 
    sub_dir = file.path(sub_dir, "outs/transcripts.csv.gz")
  )

# list file directories for Xenium unimodal Proseg segmentation
xenium_dir <- "/vast/projects/SpatialBench/segmentation/Proseg"
xenium_uni_proseg <- tribble(
  ~sub_dir,                                        ~sample,
  "proseg_xenium__Batch27__0017329__Region_2",     "ctrl172",
  "proseg_xenium__Batch24__0011456__Region_4",     "ctrl173",
  "proseg_xenium__Batch27__0017329__Region_4",     "ctrl174",
  "proseg_xenium__Batch24__0011456__Region_1",     "ko166",
  "proseg_xenium__Batch24__0011456__Region_2",     "ko167",
  "proseg_xenium__Batch27__0017329__Region_3",     "ko168",
  "proseg_xenium__Batch27__0017329__Region_5",     "wt709",
  "proseg_xenium__Batch27__0017329__Region_1",     "wt710",
  "proseg_xenium__Batch24__0011456__Region_3",     "wt713"
) %>% mutate(
  segmentation_dir = xenium_dir,
  sub_dir = file.path(sub_dir, "output/transcript-metadata.csv.gz")
)

# list file directories for Xenium multimodal default segmentation
xenium_dir <- "/vast/projects/SpatialBench/data/xenium/G000218_Benchmarking/20241031__014144__20241031_G000218_Batch34_Xen"
xenium_multi_default <- tribble(
  ~sub_dir,                                                  ~sample,
  "output-XETG00068__0032118__Region_1__20241031__014253",   "ctrl173",
  "output-XETG00068__0032118__Region_2__20241031__014253",   "wt710",
  "output-XETG00068__0032118__Region_3__20241031__014253",   "wt709",
  "output-XETG00068__0032118__Region_4__20241031__014253",   "ctrl174",
  "output-XETG00068__0032118__Region_5__20241031__014253",   "ctrl172"
) %>% mutate(
  segmentation_dir = xenium_dir,
  sub_dir = file.path(sub_dir, "transcripts.parquet")
)

# list file directories for Xenium multimodal Proseg segmentation
xenium_dir <- "/vast/projects/SpatialBench/segmentation/Proseg"
xenium_multi_proseg <- tribble(
  ~sub_dir,                                      ~sample,
  "proseg_xenium__Batch34__0032118__Region_1",   "ctrl173",
  "proseg_xenium__Batch34__0032118__Region_2",   "wt710",
  "proseg_xenium__Batch34__0032118__Region_3",   "wt709",
  "proseg_xenium__Batch34__0032118__Region_4",   "ctrl174",
  "proseg_xenium__Batch34__0032118__Region_5",   "ctrl172"
) %>% mutate(
  segmentation_dir = xenium_dir,
  sub_dir = file.path(sub_dir, "output/transcript-metadata.csv.gz")
)

# list file directories for Xenium's Cellpose2 segmentation
xenium_dir <- "/vast/projects/SpatialBench/segmentation/Cellpose/resegmentation_cellpose_xenium"
xenium_cellpose <- tribble(
  ~sub_dir,                                                                 ~sample,
  "output-XETG00068__0032118__Region_1__20241031__014253_CP_resegmented",   "ctrl173",
  "output-XETG00068__0032118__Region_2__20241031__014253_CP_resegmented",   "wt710",
  "output-XETG00068__0032118__Region_3__20241031__014253_CP_resegmented",   "wt709",
  "output-XETG00068__0032118__Region_4__20241031__014253_CP_resegmented",   "ctrl174",
  "output-XETG00068__0032118__Region_5__20241031__014253_CP_resegmented",   "ctrl172"
) %>% mutate(
  segmentation_dir = xenium_dir,
  sub_dir = file.path(sub_dir, "outs/transcripts.parquet")
)

# merge multiple Xenium samples
xenium_info <- data.table::rbindlist(
  Map(function(df, model, seg) {
    df %>% mutate(platform = "Xenium", model = model, segmentation = seg)
  }, 
  list(xenium_uni_default, xenium_uni_proseg, xenium_multi_default, xenium_multi_proseg, xenium_cellpose), 
  c("unimodal", "unimodal", "multimodal", "multimodal", "multimodal"),
  c("Default", "Proseg", "Default", "Proseg", "Cellpose2"))
)

# merge MERSCOPE and Xenium samples
sample_info <- data.table::rbindlist(
  list(merscope_info, xenium_info),
  use.names = TRUE
) %>% 
  mutate(file_path = file.path(segmentation_dir, sub_dir))
