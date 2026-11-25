
# ============================================================
# Batch SRA Processing Pipeline
# Windows + SRA Toolkit + DADA2 + SILVA v138.2
#
# Workflow:
# SRA
#  ↓
# FASTQ
#  ↓
# Quality filtering
#  ↓
# Error-rate learning
#  ↓
# FASTQ dereplication
#  ↓
# Forward and reverse DADA2 denoising
#  ↓
# Paired-end merging
#  ↓
# V3-V4 ASV table construction
#  ↓
# Chimera removal
#  ↓
# V4 region extraction
#  ↓
# Collapsing identical V4 sequences
#  ↓
# Taxonomic assignment with SILVA v138.2
#  ↓
# Genus-level abundance and relative abundance
#  ↓
# PICRUSt2 input files
#

# Output structure:
# One output directory per SRA sample.
# All outputs are written directly to the sample directory.
# ============================================================


# ============================================================
# 0. Configuration
# ============================================================

# Root directory containing SRA files.
input_root <- "C:/Users/Ruihong/Desktop/gut-microbiota/SRA_data"

# SRA Toolkit
sra_tool_dir <- "C:/sratoolkit.3.4.1/bin"

fastq_dump_path <- file.path(
  sra_tool_dir,
  "fastq-dump-orig.exe"
)

# SILVA v138.2
silva_genus <- file.path(
  input_root,
  "SILVA",
  "silva_nr99_v138.2_toGenus_trainset.fa.gz"
)


# ------------------------------------------------------------
# CPU configuration
# ------------------------------------------------------------

available_cores <- parallel::detectCores(
  logical = TRUE
)

if (is.na(available_cores)) {
  available_cores <- 4L
}

# Use up to 8 threads while reserving one logical core.
DADA_THREADS <- as.integer(
  max(
    1,
    min(
      8,
      available_cores - 1
    )
  )
)

# Windows filtering configuration
FILTER_MULTITHREAD <- FALSE

# OMP is used only by filterAndTrim.
FILTER_OMP <- TRUE

# FASTQ filtering chunk size.
FILTER_CHUNK_SIZE <- 100000L

# FASTQ dereplication chunk size.
DEREP_CHUNK_SIZE <- 1000000L

# Maximum number of bases used for error-rate learning.
ERROR_NBASES <- 1e8


# ------------------------------------------------------------
# DADA2 parameters
# ------------------------------------------------------------

MAX_EE <- c(2, 2)

MAX_N <- 0

TRUNC_Q <- 2

# Retain full read lengths (no truncation).
TRUNC_LEN <- c(0, 0)


# ------------------------------------------------------------
# V4 extraction parameters
# ------------------------------------------------------------

# Primer-based V4 extraction using the 515F/806R pair.
# Primer sequences should match the amplicon-specific primers.

PRIMER_515F <- "GTGYCAGCMGCCGCGGTAA"

PRIMER_806R <- "GGACTACNVGGGTWTCTAAT"

PRIMER_MAX_MISMATCH <- 3

V4_MIN_LENGTH <- 180

V4_MAX_LENGTH <- 300

V4_WARNING_RETENTION <- 90


# ------------------------------------------------------------
# Batch processing options
# ------------------------------------------------------------

# TEST_ONLY controls single-sample versus full-batch processing.

TEST_ONLY <- FALSE

TEST_SRA_ID <- "11_Ctrl_1_ERR11788629"

# Skip samples with complete output files.
SKIP_COMPLETED <- TRUE

# FALSE: remove FASTQ and filtered files after successful processing.
# TRUE: retain FASTQ and filtered files.
KEEP_INTERMEDIATE_FILES <- FALSE


# ============================================================
# 1. Install required R packages
# ============================================================

cat("\n检查必要的 R 包...\n")

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
  
}

required_packages <- c(
  "dada2",
  "Biostrings"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  cat(
    "需要安装：",
    paste(missing_packages, collapse = ", "),
    "\n"
  )
  
  BiocManager::install(
    missing_packages,
    update = FALSE,
    ask = FALSE
  )
  
} else {
  
  cat("必要的 R 包已经安装。\n")
  
}

suppressPackageStartupMessages({
  
  library(dada2)
  
  library(Biostrings)
  
})

# Configure the OpenMP environment.
Sys.setenv(
  OMP_NUM_THREADS = DADA_THREADS
)

cat("\nCPU 设置：\n")

cat(
  "逻辑核心数：",
  available_cores,
  "\n"
)

cat(
  "DADA2 线程数：",
  DADA_THREADS,
  "\n"
)


# ============================================================
# 2. Validate input paths
# ============================================================

if (!dir.exists(input_root)) {
  stop("SRA 总目录不存在！")
}

if (!file.exists(fastq_dump_path)) {
  stop("找不到 fastq-dump-orig.exe！")
}

if (!file.exists(silva_genus)) {
  stop("找不到 SILVA v138.2 数据库！")
}

cat("\nSRA 转换程序：\n")
cat(fastq_dump_path, "\n")

cat("\nSILVA 数据库：\n")
cat(silva_genus, "\n")


# ============================================================
# 3. Define the FASTQ discovery function
# ============================================================

find_fastq <- function(sample_dir, mate) {
  
  pattern <- paste0(
    "_",
    mate,
    "\\.fastq(\\.gz)?$"
  )
  
  files <- list.files(
    
    sample_dir,
    
    pattern = pattern,
    
    full.names = TRUE,
    
    recursive = TRUE,
    
    ignore.case = TRUE
    
  )
  
  if (length(files) == 0) {
    return(NA_character_)
  }
  
  root_path <- normalizePath(
    sample_dir,
    winslash = "/",
    mustWork = FALSE
  )
  
  file_dirs <- normalizePath(
    dirname(files),
    winslash = "/",
    mustWork = FALSE
  )
  
  # Prefer FASTQ files in the sample directory root.
  root_files <- files[
    file_dirs == root_path
  ]
  
  if (length(root_files) == 1) {
    return(root_files)
  }
  
  if (length(root_files) > 1) {
    stop(
      paste(
        "发现多个正向或反向 FASTQ 文件：",
        mate
      )
    )
  }
  
  # Support FASTQ files in legacy subdirectories.
  if (length(files) == 1) {
    return(files)
  }
  
  stop(
    paste(
      "发现多个无法区分的 FASTQ 文件：",
      mate
    )
  )
  
}


# ============================================================
# 4. Validate paired-end FASTQ files
# ============================================================

valid_fastq_pair <- function(fnF, fnR) {
  
  if (
    length(fnF) != 1 ||
    length(fnR) != 1 ||
    is.na(fnF) ||
    is.na(fnR)
  ) {
    return(FALSE)
  }
  
  if (
    !file.exists(fnF) ||
    !file.exists(fnR)
  ) {
    return(FALSE)
  }
  
  sizes <- file.info(
    c(fnF, fnR)
  )$size
  
  if (
    length(sizes) != 2 ||
    any(is.na(sizes)) ||
    any(sizes <= 0)
  ) {
    return(FALSE)
  }
  
  TRUE
  
}


# ============================================================
# 5. Define V4 primers
# ============================================================

primer_515F <- DNAString(
  PRIMER_515F
)

primer_806R <- DNAString(
  PRIMER_806R
)

primer_806R_rc <- reverseComplement(
  primer_806R
)


# ============================================================
# 6. Define the V4 extraction function
# ============================================================

extract_V4 <- function(sequence) {
  
  subject <- DNAString(
    toupper(sequence)
  )
  
  # Locate the 515F primer.
  hits_F <- matchPattern(
    
    primer_515F,
    
    subject,
    
    max.mismatch = PRIMER_MAX_MISMATCH,
    
    with.indels = FALSE,
    
    fixed = "subject"
    
  )
  
  # Locate the reverse complement of the 806R primer.
  hits_R <- matchPattern(
    
    primer_806R_rc,
    
    subject,
    
    max.mismatch = PRIMER_MAX_MISMATCH,
    
    with.indels = FALSE,
    
    fixed = "subject"
    
  )
  
  if (
    length(hits_F) == 0 ||
    length(hits_R) == 0
  ) {
    
    return(NA_character_)
    
  }
  
  f_end <- end(hits_F)
  
  r_start <- start(hits_R)
  
  # Generate candidate forward/reverse primer pairs.
  pairs <- expand.grid(
    
    fi = seq_along(f_end),
    
    ri = seq_along(r_start)
    
  )
  
  # Retain pairs with the forward primer upstream of the reverse primer.
  pairs <- pairs[
    
    f_end[pairs$fi] < r_start[pairs$ri],
    
    ,
    
    drop = FALSE
    
  ]
  
  if (nrow(pairs) == 0) {
    return(NA_character_)
  }
  
  insert_start <- f_end[pairs$fi] + 1
  
  insert_end <- r_start[pairs$ri] - 1
  
  insert_length <- (
    insert_end - insert_start + 1
  )
  
  # Retain sequences within the expected V4 length range.
  valid <- which(
    
    insert_length >= V4_MIN_LENGTH &
      
      insert_length <= V4_MAX_LENGTH
    
  )
  
  if (length(valid) == 0) {
    return(NA_character_)
  }
  
  # Select the region closest to 250 bp.
  best <- valid[
    which.min(
      abs(
        insert_length[valid] - 250
      )
    )
  ]
  
  result <- subseq(
    
    subject,
    
    start = insert_start[best],
    
    end = insert_end[best]
    
  )
  
  as.character(result)
  
}


# ============================================================
# 7. Identify completed samples
# ============================================================

is_sample_completed <- function(sample_dir) {
  
  required_files <- c(
    
    "V4_ASV_taxonomy_table.csv",
    
    "V4_genus_abundance.csv",
    
    "V4_genus_proportion.csv",
    
    "V4_DADA2_read_tracking.csv",
    
    "V4_ASV_sequences.fasta",
    
    "V4_ASV_abundance.tsv",
    
    "V4_extraction_quality_report.csv",
    
    "V4_complete_analysis_results.rds"
    
  )
  
  required_paths <- file.path(
    sample_dir,
    required_files
  )
  
  if (!all(file.exists(required_paths))) {
    return(FALSE)
  }
  
  sizes <- file.info(
    required_paths
  )$size
  
  if (
    any(is.na(sizes)) ||
    any(sizes <= 0)
  ) {
    return(FALSE)
  }
  
  TRUE
  
}


# ============================================================
# 8. Process an individual SRA sample
# ============================================================

process_one_sra <- function(sra_file) {
  
  sra_id <- tools::file_path_sans_ext(
    basename(sra_file)
  )
  
  # Create one output directory per SRA sample.
  sample_dir <- file.path(
    dirname(sra_file),
    sra_id
  )
  
  success_marker <- file.path(
    sample_dir,
    "PROCESSING_SUCCESS.txt"
  )
  
  
  # ----------------------------------------------------------
  # 8.1 Skip completed samples
  # ----------------------------------------------------------
  
  if (
    SKIP_COMPLETED &&
    dir.exists(sample_dir) &&
    is_sample_completed(sample_dir)
  ) {
    
    cat(
      "\n完整结果已经存在，跳过：",
      sra_id,
      "\n"
    )
    
    return(
      
      list(
        
        status = "SKIPPED_COMPLETE",
        
        sra_id = sra_id,
        
        original_asv = NA,
        
        v4_asv = NA,
        
        v4_retention = NA,
        
        merged_reads = NA,
        
        nonchim_reads = NA,
        
        genus_count = NA,
        
        output_dir = sample_dir,
        
        message = "完整结果已经存在"
        
      )
      
    )
    
  }
  
  
  # ----------------------------------------------------------
  # 8.2 Create the sample output directory
  # ----------------------------------------------------------
  
  dir.create(
    
    sample_dir,
    
    recursive = TRUE,
    
    showWarnings = FALSE
    
  )
  
  # Remove any previous completion marker before processing.
  if (file.exists(success_marker)) {
    file.remove(success_marker)
  }
  
  cat("\n")
  cat("====================================================\n")
  
  cat("开始处理：", sra_id, "\n")
  
  cat("SRA 文件：", sra_file, "\n")
  
  cat("结果目录：", sample_dir, "\n")
  
  cat("====================================================\n")
  
  
  # ----------------------------------------------------------
  # 8.3 Convert SRA to FASTQ
  # ----------------------------------------------------------
  
  cat("\n检查 FASTQ 文件...\n")
  
  fnF <- find_fastq(
    sample_dir,
    "1"
  )
  
  fnR <- find_fastq(
    sample_dir,
    "2"
  )
  
  if (!valid_fastq_pair(fnF, fnR)) {
    
    cat("\n需要转换 SRA 为 FASTQ...\n")
    
    # Remove incomplete FASTQ files from previous runs.
    old_fastq <- list.files(
      
      sample_dir,
      
      pattern = "_[12]\\.fastq(\\.gz)?$",
      
      full.names = TRUE,
      
      recursive = TRUE,
      
      ignore.case = TRUE
      
    )
    
    if (length(old_fastq) > 0) {
      file.remove(old_fastq)
    }
    
    sra_path <- normalizePath(
      
      sra_file,
      
      winslash = "\\",
      
      mustWork = TRUE
      
    )
    
    output_path <- normalizePath(
      
      sample_dir,
      
      winslash = "\\",
      
      mustWork = TRUE
      
    )
    
    # Use fastq-dump-orig.exe for SRA conversion.
    cmd <- paste(
      
      shQuote(
        fastq_dump_path,
        type = "cmd"
      ),
      
      "--split-files",
      
      "--outdir",
      
      shQuote(
        output_path,
        type = "cmd"
      ),
      
      shQuote(
        sra_path,
        type = "cmd"
      )
      
    )
    
    cat("\n执行命令：\n")
    cat(cmd, "\n")
    
    dump_status <- system(
      
      cmd,
      
      wait = TRUE,
      
      intern = FALSE,
      
      show.output.on.console = TRUE
      
    )
    
    if (
      is.null(dump_status) ||
      length(dump_status) != 1 ||
      is.na(dump_status) ||
      dump_status != 0
    ) {
      
      stop("fastq-dump-orig.exe 转换失败！")
      
    }
    
  } else {
    
    cat("\n发现已有完整 FASTQ，跳过转换。\n")
    
  }
  
  
  # ----------------------------------------------------------
  # 8.4 Revalidate FASTQ files
  # ----------------------------------------------------------
  
  fnF <- find_fastq(
    sample_dir,
    "1"
  )
  
  fnR <- find_fastq(
    sample_dir,
    "2"
  )
  
  if (!valid_fastq_pair(fnF, fnR)) {
    
    stop(
      "FASTQ 转换结束，但没有找到完整的双端文件！"
    )
    
  }
  
  # Move legacy FASTQ files to the sample directory root.
  for (f in c(fnF, fnR)) {
    
    root_normalized <- normalizePath(
      sample_dir,
      winslash = "/"
    )
    
    current_normalized <- normalizePath(
      dirname(f),
      winslash = "/"
    )
    
    if (current_normalized != root_normalized) {
      
      target <- file.path(
        sample_dir,
        basename(f)
      )
      
      if (!file.copy(
        f,
        target,
        overwrite = TRUE
      )) {
        
        stop("无法将 FASTQ 移动到样本文件夹根目录！")
        
      }
      
      file.remove(f)
      
    }
    
  }
  
  fnF <- file.path(
    sample_dir,
    basename(fnF)
  )
  
  fnR <- file.path(
    sample_dir,
    basename(fnR)
  )
  
  
  # ----------------------------------------------------------
  # 8.5 Quality filtering
  # ----------------------------------------------------------
  
  cat("\n开始质量过滤...\n")
  
  filtF <- file.path(
    
    sample_dir,
    
    paste0(
      sra_id,
      "_1_filtered.fastq.gz"
    )
    
  )
  
  filtR <- file.path(
    
    sample_dir,
    
    paste0(
      sra_id,
      "_2_filtered.fastq.gz"
    )
    
  )
  
  # Remove previous filtered files before processing.
  unlink(
    c(filtF, filtR)
  )
  
  ft_out <- dada2::filterAndTrim(
    
    fwd = fnF,
    
    filt = filtF,
    
    rev = fnR,
    
    filt.rev = filtR,
    
    truncLen = TRUNC_LEN,
    
    maxN = MAX_N,
    
    maxEE = MAX_EE,
    
    truncQ = TRUNC_Q,
    
    rm.phix = TRUE,
    
    compress = TRUE,
    
    # Disable fork-based parallelization on Windows.
    multithread = FILTER_MULTITHREAD,
    
    # Enable OMP only for filterAndTrim.
    OMP = FILTER_OMP,
    
    n = FILTER_CHUNK_SIZE,
    
    verbose = TRUE
    
  )
  
  cat("\n质量过滤结果：\n")
  
  print(ft_out)
  
  if (
    is.null(ft_out) ||
    nrow(ft_out) != 1 ||
    is.na(ft_out[1, "reads.out"]) ||
    ft_out[1, "reads.out"] <= 0
  ) {
    
    stop("质量过滤失败或没有剩余 reads！")
    
  }
  
  if (
    !file.exists(filtF) ||
    !file.exists(filtR)
  ) {
    
    stop("没有生成完整的过滤后 FASTQ 文件！")
    
  }
  
  
  # ----------------------------------------------------------
  # 8.6 Learn forward error rates
  # ----------------------------------------------------------
  
  cat("\n开始学习正向错误率...\n")
  
  errF <- dada2::learnErrors(
    
    filtF,
    
    nbases = ERROR_NBASES,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # 8.7 Learn reverse error rates
  # ----------------------------------------------------------
  
  cat("\n开始学习反向错误率...\n")
  
  errR <- dada2::learnErrors(
    
    filtR,
    
    nbases = ERROR_NBASES,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # 8.8 Dereplicate FASTQ reads
  # ----------------------------------------------------------
  
  cat("\n开始正向 FASTQ 去重复...\n")
  
  # derepFastq does not accept an OMP argument.
  derepF <- dada2::derepFastq(
    
    filtF,
    
    n = DEREP_CHUNK_SIZE,
    
    verbose = TRUE
    
  )
  
  cat("\n开始反向 FASTQ 去重复...\n")
  
  derepR <- dada2::derepFastq(
    
    filtR,
    
    n = DEREP_CHUNK_SIZE,
    
    verbose = TRUE
    
  )
  
  # Validate dereplication object classes.
  if (!inherits(derepF, "derep")) {
    
    stop("正向 derepFastq 没有返回 derep 对象！")
    
  }
  
  if (!inherits(derepR, "derep")) {
    
    stop("反向 derepFastq 没有返回 derep 对象！")
    
  }
  
  cat("\nFASTQ 去重复成功！\n")
  
  
  # ----------------------------------------------------------
  # 8.9 Perform forward DADA2 denoising
  # ----------------------------------------------------------
  
  cat("\n开始正向 DADA2 去噪...\n")
  
  dadaF <- dada2::dada(
    
    derepF,
    
    err = errF,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  if (!inherits(dadaF, "dada")) {
    
    stop("正向 DADA2 没有返回 dada 对象！")
    
  }
  
  
  # ----------------------------------------------------------
  # 8.10 Perform reverse DADA2 denoising
  # ----------------------------------------------------------
  
  cat("\n开始反向 DADA2 去噪...\n")
  
  dadaR <- dada2::dada(
    
    derepR,
    
    err = errR,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  if (!inherits(dadaR, "dada")) {
    
    stop("反向 DADA2 没有返回 dada 对象！")
    
  }
  
  cat("\n正反向 DADA2 去噪成功！\n")
  
  
  # ----------------------------------------------------------
  # 8.11 Merge paired-end reads
  # ----------------------------------------------------------
  
  cat("\n开始双端合并...\n")
  
  # Use the individual sample objects directly.
  # mergePairs returns a data.frame for a single sample.
  
  merger <- dada2::mergePairs(
    
    dadaF,
    
    derepF,
    
    dadaR,
    
    derepR,
    
    minOverlap = 12,
    
    maxMismatch = 0,
    
    verbose = TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # Validate the mergePairs output.
  # ----------------------------------------------------------
  
  if (!is.data.frame(merger)) {
    
    stop(
      paste(
        "mergePairs 返回的数据类型不正确：",
        paste(class(merger), collapse = ", ")
      )
    )
    
  }
  
  if (!"abundance" %in% colnames(merger)) {
    
    stop("mergePairs 结果中不存在 abundance 列！")
    
  }
  
  # Sum the abundances of successfully merged reads.
  
  total_merged <- sum(
    
    merger$abundance,
    
    na.rm = TRUE
    
  )
  
  cat("\n双端合并统计：\n")
  
  cat(
    "合并成功的 unique pairings：",
    nrow(merger),
    "\n"
  )
  
  cat(
    "合并成功 reads：",
    total_merged,
    "\n"
  )
  
  merge_rate <- (
    total_merged /
      as.numeric(ft_out[1, "reads.out"])
  ) * 100
  
  cat(
    "双端合并成功率：",
    round(merge_rate, 2),
    "%\n"
  )
  
  if (
    length(total_merged) != 1 ||
    is.na(total_merged) ||
    !is.finite(total_merged) ||
    total_merged <= 0
  ) {
    
    stop("双端合并没有产生有效 reads！")
    
  }
  
  # Save the mergePairs result.
  saveRDS(
    
    merger,
    
    file.path(
      sample_dir,
      "mergePairs_result.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.12 Construct the V3-V4 ASV table
  # ----------------------------------------------------------
  
  cat("\n开始构建 V3-V4 ASV 表...\n")
  
  # makeSequenceTable requires a named list of merge results.
  # Wrap the single-sample merge result in a named list.
  
  merger_list <- setNames(
    
    list(merger),
    
    sra_id
    
  )
  
  seqtab <- dada2::makeSequenceTable(
    merger_list
  )
  
  if (
    is.null(seqtab) ||
    nrow(seqtab) != 1 ||
    ncol(seqtab) == 0
  ) {
    
    stop("V3-V4 ASV 表构建失败！")
    
  }
  
  rownames(seqtab) <- sra_id
  
  cat("\nV3-V4 ASV 表构建成功！\n")
  
  cat(
    "V3-V4 ASV 数量：",
    ncol(seqtab),
    "\n"
  )
  
  saveRDS(
    
    seqtab,
    
    file.path(
      sample_dir,
      "V3V4_ASV_table_before_chimera.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.13 Remove chimeric sequences
  # ----------------------------------------------------------
  
  cat("\n开始去除嵌合体...\n")
  
  seqtab.nochim <- dada2::removeBimeraDenovo(
    
    seqtab,
    
    method = "consensus",
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  if (
    is.null(seqtab.nochim) ||
    nrow(seqtab.nochim) != 1 ||
    ncol(seqtab.nochim) == 0
  ) {
    
    stop("去除嵌合体后没有剩余 ASV！")
    
  }
  
  rownames(seqtab.nochim) <- sra_id
  
  saveRDS(
    
    seqtab.nochim,
    
    file.path(
      sample_dir,
      "V3V4_ASV_table_nochim.rds"
    )
    
  )
  
  original_asv_count <- ncol(
    seqtab.nochim
  )
  
  nonchim_reads <- sum(
    seqtab.nochim
  )
  
  cat("\n去嵌合体成功！\n")
  
  cat(
    "非嵌合 ASV 数量：",
    original_asv_count,
    "\n"
  )
  
  cat(
    "非嵌合 reads：",
    nonchim_reads,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 8.14 Extract the V4 region
  # ----------------------------------------------------------
  
  cat("\n开始提取统一 V4 区域...\n")
  
  v3v4_sequences <- colnames(
    seqtab.nochim
  )
  
  v4_sequences <- vapply(
    
    v3v4_sequences,
    
    extract_V4,
    
    FUN.VALUE = character(1)
    
  )
  
  v4_success <- (
    
    !is.na(v4_sequences) &
      
      nzchar(v4_sequences)
    
  )
  
  if (!any(v4_success)) {
    
    stop(
      paste(
        "所有 ASV 均未能提取 V4。",
        "请检查实际扩增引物和 V4 定位引物。"
      )
    )
    
  }
  
  # Calculate the ASV-level extraction rate.
  v4_asv_success_rate <- mean(
    v4_success
  ) * 100
  
  # Calculate the read-weighted V4 retention rate.
  v4_retained_reads <- sum(
    
    seqtab.nochim[
      ,
      v4_success,
      drop = FALSE
    ]
    
  )
  
  v4_retention <- (
    v4_retained_reads /
      nonchim_reads
  ) * 100
  
  cat("\nV4 提取统计：\n")
  
  cat(
    "原始 V3-V4 ASV：",
    original_asv_count,
    "\n"
  )
  
  cat(
    "成功提取 V4 的 ASV：",
    sum(v4_success),
    "\n"
  )
  
  cat(
    "ASV 提取成功率：",
    round(v4_asv_success_rate, 2),
    "%\n"
  )
  
  cat(
    "V4 reads 保留率：",
    round(v4_retention, 2),
    "%\n"
  )
  
  if (v4_retention < V4_WARNING_RETENTION) {
    
    warning(
      paste0(
        sra_id,
        " V4 reads 保留率低于 ",
        V4_WARNING_RETENTION,
        "%，建议检查引物定位。"
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # 8.15 Export the V4 extraction quality report
  # ----------------------------------------------------------
  
  v4_report <- data.frame(
    
    Original_ASV = paste0(
      "V3V4_ASV_",
      seq_along(v3v4_sequences)
    ),
    
    V3V4_Sequence = v3v4_sequences,
    
    V3V4_Length = nchar(v3v4_sequences),
    
    V4_Sequence = v4_sequences,
    
    V4_Length = ifelse(
      v4_success,
      nchar(v4_sequences),
      NA
    ),
    
    Extracted = v4_success,
    
    Total_Reads = as.numeric(
      seqtab.nochim[1, ]
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  write.csv(
    
    v4_report,
    
    file.path(
      sample_dir,
      "V4_extraction_quality_report.csv"
    ),
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  # ----------------------------------------------------------
  # 8.16 Collapse identical V4 sequences
  # ----------------------------------------------------------
  
  cat("\n开始合并相同 V4 序列...\n")
  
  asv_counts_kept <- t(
    
    seqtab.nochim[
      ,
      v4_success,
      drop = FALSE
    ]
    
  )
  
  # Sum counts for ASVs sharing an identical V4 sequence.
  v4_collapsed <- rowsum(
    
    asv_counts_kept,
    
    group = as.character(
      v4_sequences[v4_success]
    ),
    
    reorder = FALSE
    
  )
  
  v4_seqtab <- t(
    as.matrix(v4_collapsed)
  )
  
  rownames(v4_seqtab) <- sra_id
  
  colnames(v4_seqtab) <- rownames(
    v4_collapsed
  )
  
  v4_unique_sequences <- colnames(
    v4_seqtab
  )
  
  if (ncol(v4_seqtab) == 0) {
    stop("合并相同 V4 序列后没有 ASV！")
  }
  
  cat(
    "\n最终 V4 ASV 数量：",
    ncol(v4_seqtab),
    "\n"
  )
  
  saveRDS(
    
    v4_seqtab,
    
    file.path(
      sample_dir,
      "V4_ASV_table.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.17 Assign standardized V4 ASV identifiers
  # ----------------------------------------------------------
  
  v4_ids <- paste0(
    
    "V4_ASV_",
    
    sprintf(
      "%06d",
      seq_along(v4_unique_sequences)
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.18 Construct the V4 ASV abundance matrix
  # ----------------------------------------------------------
  
  v4_count_matrix <- t(
    v4_seqtab
  )
  
  rownames(v4_count_matrix) <- v4_unique_sequences
  
  colnames(v4_count_matrix) <- sra_id
  
  
  # ----------------------------------------------------------
  # 8.19 Export the PICRUSt2 FASTA file
  # ----------------------------------------------------------
  
  cat("\n生成 PICRUSt2 FASTA...\n")
  
  fasta_set <- DNAStringSet(
    v4_unique_sequences
  )
  
  names(fasta_set) <- v4_ids
  
  Biostrings::writeXStringSet(
    
    fasta_set,
    
    filepath = file.path(
      sample_dir,
      "V4_ASV_sequences.fasta"
    ),
    
    format = "fasta"
    
  )
  
  
  # ----------------------------------------------------------
  # 8.20 Export the PICRUSt2 abundance table
  # ----------------------------------------------------------
  
  cat("\n生成 PICRUSt2 abundance TSV...\n")
  
  picrust_abundance <- data.frame(
    
    "#OTU ID" = v4_ids,
    
    as.data.frame(
      v4_count_matrix,
      check.names = FALSE
    ),
    
    check.names = FALSE
    
  )
  
  write.table(
    
    picrust_abundance,
    
    file.path(
      sample_dir,
      "V4_ASV_abundance.tsv"
    ),
    
    sep = "\t",
    
    quote = FALSE,
    
    row.names = FALSE,
    
    col.names = TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # 8.21 Assign taxonomy using SILVA v138.2
  # ----------------------------------------------------------
  
  cat("\n开始 SILVA v138.2 分类注释...\n")
  
  v4_taxa <- dada2::assignTaxonomy(
    
    v4_unique_sequences,
    
    refFasta = silva_genus,
    
    minBoot = 50,
    
    tryRC = FALSE,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  if (
    is.null(v4_taxa) ||
    nrow(v4_taxa) != length(v4_unique_sequences)
  ) {
    
    stop("SILVA 分类注释失败！")
    
  }
  
  saveRDS(
    
    v4_taxa,
    
    file.path(
      sample_dir,
      "V4_SILVA_taxonomy.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.22 Construct the V4 ASV taxonomy table
  # ----------------------------------------------------------
  
  cat("\n构建 V4 ASV taxonomy table...\n")
  
  if (is.null(rownames(v4_taxa))) {
    
    stop("分类注释结果缺少序列行名！")
    
  }
  
  taxa_index <- match(
    
    v4_unique_sequences,
    
    rownames(v4_taxa)
    
  )
  
  if (anyNA(taxa_index)) {
    
    stop("分类注释结果与 V4 ASV 序列无法匹配！")
    
  }
  
  taxa_match <- v4_taxa[
    taxa_index,
    ,
    drop = FALSE
  ]
  
  taxa_df <- as.data.frame(
    
    taxa_match,
    
    stringsAsFactors = FALSE,
    
    check.names = FALSE
    
  )
  
  rownames(taxa_df) <- NULL
  
  # Store the representative V4 sequence in the first column.
  # Retain ASV_ID values for matching PICRUSt2 inputs.
  
  final_table <- data.frame(
    
    Sequence = v4_unique_sequences,
    
    ASV_ID = v4_ids,
    
    as.data.frame(
      v4_count_matrix,
      check.names = FALSE
    ),
    
    taxa_df,
    
    check.names = FALSE
    
  )
  
  write.csv(
    
    final_table,
    
    file.path(
      sample_dir,
      "V4_ASV_taxonomy_table.csv"
    ),
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  # ----------------------------------------------------------
  # 8.23 Aggregate genus-level abundance
  # ----------------------------------------------------------
  
  cat("\n开始属水平丰度汇总...\n")
  
  if (!"Genus" %in% colnames(taxa_match)) {
    
    stop("SILVA 注释结果没有 Genus 列！")
    
  }
  
  genus_label <- as.character(
    taxa_match[, "Genus"]
  )
  
  missing_genus <- (
    
    is.na(genus_label) |
      
      !nzchar(trimws(genus_label))
    
  )
  
  if ("Family" %in% colnames(taxa_match)) {
    
    family_label <- as.character(
      taxa_match[, "Family"]
    )
    
    family_label[
      
      is.na(family_label) |
        
        !nzchar(trimws(family_label))
      
    ] <- "Family_NA"
    
    genus_label[missing_genus] <- paste0(
      
      "Unassigned_",
      
      family_label[missing_genus]
      
    )
    
  } else {
    
    genus_label[missing_genus] <- "Unassigned"
    
  }
  
  # Aggregate ASV read counts by genus.
  genus_matrix <- rowsum(
    
    v4_count_matrix,
    
    group = genus_label,
    
    reorder = TRUE
    
  )
  
  genus_order <- order(
    
    rowSums(genus_matrix),
    
    decreasing = TRUE
    
  )
  
  genus_matrix <- genus_matrix[
    genus_order,
    ,
    drop = FALSE
  ]
  
  genus_abundance <- data.frame(
    
    Genus = rownames(genus_matrix),
    
    as.data.frame(
      genus_matrix,
      check.names = FALSE
    ),
    
    check.names = FALSE
    
  )
  
  write.csv(
    
    genus_abundance,
    
    file.path(
      sample_dir,
      "V4_genus_abundance.csv"
    ),
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  # ----------------------------------------------------------
  # 8.24 Calculate genus-level relative abundance
  # ----------------------------------------------------------
  
  cat("\n计算属水平相对丰度...\n")
  
  genus_totals <- colSums(
    genus_matrix
  )
  
  if (any(genus_totals <= 0)) {
    
    stop("属水平丰度总数为0，无法计算相对丰度！")
    
  }
  
  genus_prop_matrix <- sweep(
    
    genus_matrix,
    
    2,
    
    genus_totals,
    
    "/"
    
  ) * 100
  
  genus_proportion <- data.frame(
    
    Genus = rownames(genus_prop_matrix),
    
    as.data.frame(
      genus_prop_matrix,
      check.names = FALSE
    ),
    
    check.names = FALSE
    
  )
  
  write.csv(
    
    genus_proportion,
    
    file.path(
      sample_dir,
      "V4_genus_proportion.csv"
    ),
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  # ----------------------------------------------------------
  # 8.25 Summarize DADA2 read counts
  # ----------------------------------------------------------
  
  cat("\n生成 reads 统计表...\n")
  
  denoised_F <- sum(
    dada2::getUniques(dadaF)
  )
  
  denoised_R <- sum(
    dada2::getUniques(dadaR)
  )
  
  read_tracking <- data.frame(
    
    Sample = sra_id,
    
    Raw_pairs = as.numeric(
      ft_out[1, "reads.in"]
    ),
    
    Filtered_pairs = as.numeric(
      ft_out[1, "reads.out"]
    ),
    
    Denoised_F = denoised_F,
    
    Denoised_R = denoised_R,
    
    Merged = total_merged,
    
    Merge_rate_percent = round(
      merge_rate,
      2
    ),
    
    Nonchim = nonchim_reads,
    
    V4_retained = v4_retained_reads,
    
    V4_retention_percent = round(
      v4_retention,
      2
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  write.csv(
    
    read_tracking,
    
    file.path(
      sample_dir,
      "V4_DADA2_read_tracking.csv"
    ),
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  # ----------------------------------------------------------
  # 8.26 Save the complete analysis results
  # ----------------------------------------------------------
  
  cat("\n保存完整分析结果...\n")
  
  complete_results <- list(
    
    SRA_ID = sra_id,
    
    DADA2_forward = dadaF,
    
    DADA2_reverse = dadaR,
    
    Merged_reads = merger,
    
    V3V4_ASV_before_chimera = seqtab,
    
    V3V4_ASV = seqtab.nochim,
    
    V4_ASV = v4_seqtab,
    
    V4_taxonomy = v4_taxa,
    
    V4_final_table = final_table,
    
    V4_genus_abundance = genus_abundance,
    
    V4_genus_proportion = genus_proportion,
    
    PICRUSt2_abundance = picrust_abundance,
    
    read_tracking = read_tracking,
    
    V4_extraction_report = v4_report
    
  )
  
  saveRDS(
    
    complete_results,
    
    file.path(
      sample_dir,
      "V4_complete_analysis_results.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.27 Remove intermediate files
  # ----------------------------------------------------------
  
  if (!KEEP_INTERMEDIATE_FILES) {
    
    cat("\n清理中间 FASTQ 文件...\n")
    
    intermediate_files <- c(
      
      fnF,
      
      fnR,
      
      filtF,
      
      filtR
      
    )
    
    intermediate_files <- intermediate_files[
      file.exists(intermediate_files)
    ]
    
    if (length(intermediate_files) > 0) {
      
      file.remove(
        intermediate_files
      )
      
    }
    
  }
  
  
  # ----------------------------------------------------------
  # 8.28 Write the completion marker
  # ----------------------------------------------------------
  
  writeLines(
    
    c(
      
      paste("SRA:", sra_id),
      
      paste("Completed:", Sys.time()),
      
      paste("Merged reads:", total_merged),
      
      paste("V3-V4 ASV:", original_asv_count),
      
      paste("V4 ASV:", ncol(v4_seqtab)),
      
      paste(
        "V4 reads retention:",
        round(v4_retention, 2),
        "%"
      )
      
    ),
    
    success_marker
    
  )
  
  
  # ----------------------------------------------------------
  # 8.29 Return sample-level results
  # ----------------------------------------------------------
  
  cat("\n")
  
  cat("====================================================\n")
  
  cat("样本处理成功：", sra_id, "\n")
  
  cat("双端合并 reads：", total_merged, "\n")
  
  cat("V3-V4 ASV：", original_asv_count, "\n")
  
  cat("V4 ASV：", ncol(v4_seqtab), "\n")
  
  cat(
    "V4 reads 保留率：",
    round(v4_retention, 2),
    "%\n"
  )
  
  cat(
    "属数量：",
    nrow(genus_matrix),
    "\n"
  )
  
  cat("结果目录：", sample_dir, "\n")
  
  cat("====================================================\n")
  
  return(
    
    list(
      
      status = "SUCCESS",
      
      sra_id = sra_id,
      
      original_asv = original_asv_count,
      
      v4_asv = ncol(v4_seqtab),
      
      v4_retention = round(
        v4_retention,
        2
      ),
      
      merged_reads = total_merged,
      
      nonchim_reads = nonchim_reads,
      
      genus_count = nrow(genus_matrix),
      
      output_dir = sample_dir,
      
      message = "处理成功"
      
    )
    
  )
  
}


# ============================================================
# 9. Discover all SRA files
# ============================================================

cat("\n开始搜索 SRA 文件...\n")

sra_files <- list.files(
  
  path = input_root,
  
  pattern = "\\.sra$",
  
  full.names = TRUE,
  
  recursive = TRUE,
  
  ignore.case = TRUE
  
)

if (length(sra_files) == 0) {
  
  stop("没有找到 SRA 文件！")
  
}

sra_files <- normalizePath(
  
  sra_files,
  
  winslash = "/",
  
  mustWork = TRUE
  
)

cat(
  "\n找到 SRA 文件：",
  length(sra_files),
  "\n"
)


# ============================================================
# 10. Configure single-sample mode
# ============================================================

if (TEST_ONLY) {
  
  sra_ids <- tools::file_path_sans_ext(
    basename(sra_files)
  )
  
  selected <- which(
    sra_ids == TEST_SRA_ID
  )
  
  if (length(selected) != 1) {
    
    stop(
      paste(
        "测试样本没有找到或存在重复：",
        TEST_SRA_ID
      )
    )
    
  }
  
  sra_files <- sra_files[selected]
  
  cat("\n当前为单样本测试模式！\n")
  
  cat(
    "测试文件：",
    basename(sra_files),
    "\n"
  )
  
} else {
  
  cat("\n当前为全部 SRA 批处理模式！\n")
  
}


# ============================================================
# 11. Initialize the batch summary table
# ============================================================

summary_file <- file.path(
  
  input_root,
  
  "batch_processing_summary.csv"
  
)

summary_columns <- c(
  
  "SRA_ID",
  
  "SRA_File",
  
  "Status",
  
  "Original_V3V4_ASV",
  
  "Final_V4_ASV",
  
  "V4_Read_Retention_Percent",
  
  "Merged_Reads",
  
  "Nonchim_Reads",
  
  "Genus_Count",
  
  "Elapsed_Seconds",
  
  "Output_Directory",
  
  "Message"
  
)

empty_summary <- data.frame(
  
  SRA_ID = character(),
  
  SRA_File = character(),
  
  Status = character(),
  
  Original_V3V4_ASV = numeric(),
  
  Final_V4_ASV = numeric(),
  
  V4_Read_Retention_Percent = numeric(),
  
  Merged_Reads = numeric(),
  
  Nonchim_Reads = numeric(),
  
  Genus_Count = numeric(),
  
  Elapsed_Seconds = numeric(),
  
  Output_Directory = character(),
  
  Message = character(),
  
  stringsAsFactors = FALSE
  
)

if (file.exists(summary_file)) {
  
  old_summary <- tryCatch(
    
    read.csv(
      
      summary_file,
      
      stringsAsFactors = FALSE,
      
      check.names = FALSE
      
    ),
    
    error = function(e) NULL
    
  )
  
  if (
    is.null(old_summary) ||
    !all(summary_columns %in% colnames(old_summary))
  ) {
    
    old_summary <- empty_summary
    
  } else {
    
    old_summary <- old_summary[
      ,
      summary_columns,
      drop = FALSE
    ]
    
  }
  
} else {
  
  old_summary <- empty_summary
  
}


# ============================================================
# 12. Process SRA files sequentially
# ============================================================

for (i in seq_along(sra_files)) {
  
  current_sra <- sra_files[i]
  
  current_id <- tools::file_path_sans_ext(
    basename(current_sra)
  )
  
  start_time <- Sys.time()
  
  cat("\n\n")
  
  cat("####################################################\n")
  
  cat(
    "当前进度：",
    i,
    "/",
    length(sra_files),
    "\n"
  )
  
  cat(
    "当前文件：",
    basename(current_sra),
    "\n"
  )
  
  cat("####################################################\n")
  
  
  # ----------------------------------------------------------
  # Isolate failures at the sample level.
  # ----------------------------------------------------------
  
  result <- tryCatch(
    
    {
      
      process_one_sra(current_sra)
      
    },
    
    error = function(e) {
      
      error_message <- conditionMessage(e)
      
      error_dir <- file.path(
        dirname(current_sra),
        current_id
      )
      
      dir.create(
        
        error_dir,
        
        recursive = TRUE,
        
        showWarnings = FALSE
        
      )
      
      cat("\n处理失败：", current_id, "\n")
      
      cat("错误信息：", error_message, "\n")
      
      # Record the error call stack.
      calls <- sys.calls()
      
      trace_text <- vapply(
        
        calls,
        
        function(x) {
          
          paste(
            deparse(x),
            collapse = " "
          )
          
        },
        
        character(1)
        
      )
      
      writeLines(
        
        c(
          
          paste("SRA:", current_sra),
          
          paste("Error:", error_message),
          
          paste("Time:", Sys.time()),
          
          "",
          
          "Error traceback:",
          
          trace_text
          
        ),
        
        file.path(
          error_dir,
          "error_traceback.txt"
        )
        
      )
      
      list(
        
        status = "FAILED",
        
        sra_id = current_id,
        
        original_asv = NA,
        
        v4_asv = NA,
        
        v4_retention = NA,
        
        merged_reads = NA,
        
        nonchim_reads = NA,
        
        genus_count = NA,
        
        output_dir = error_dir,
        
        message = error_message
        
      )
      
    }
    
  )
  
  
  # ----------------------------------------------------------
  # Calculate elapsed processing time.
  # ----------------------------------------------------------
  
  elapsed_seconds <- as.numeric(
    
    difftime(
      
      Sys.time(),
      
      start_time,
      
      units = "secs"
      
    )
    
  )
  
  
  # ----------------------------------------------------------
  # Assemble the current sample summary row.
  # ----------------------------------------------------------
  
  result_row <- data.frame(
    
    SRA_ID = result$sra_id,
    
    SRA_File = current_sra,
    
    Status = result$status,
    
    Original_V3V4_ASV = result$original_asv,
    
    Final_V4_ASV = result$v4_asv,
    
    V4_Read_Retention_Percent = result$v4_retention,
    
    Merged_Reads = result$merged_reads,
    
    Nonchim_Reads = result$nonchim_reads,
    
    Genus_Count = result$genus_count,
    
    Elapsed_Seconds = round(
      elapsed_seconds,
      2
    ),
    
    Output_Directory = result$output_dir,
    
    Message = result$message,
    
    stringsAsFactors = FALSE
    
  )
  
  
  # ----------------------------------------------------------
  # Update the batch summary table.
  # ----------------------------------------------------------
  
  existing <- which(
    
    !is.na(old_summary$SRA_File) &
      
      old_summary$SRA_File == current_sra
    
  )
  
  # Preserve existing summary records for completed samples.
  
  if (
    result$status != "SKIPPED_COMPLETE" ||
    length(existing) == 0
  ) {
    
    old_summary <- old_summary[
      
      is.na(old_summary$SRA_File) |
        
        old_summary$SRA_File != current_sra,
      
      ,
      
      drop = FALSE
      
    ]
    
    old_summary <- rbind(
      
      old_summary,
      
      result_row
      
    )
    
  }
  
  # Save the summary after each sample.
  write.csv(
    
    old_summary,
    
    summary_file,
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  cat("\n当前样本状态：", result$status, "\n")
  
  cat(
    "处理耗时：",
    round(elapsed_seconds / 60, 2),
    "分钟\n"
  )
  
  cat(
    "已完成进度：",
    i,
    "/",
    length(sra_files),
    "\n"
  )
  
}


# ============================================================
# 13. Generate the final batch summary
# ============================================================

write.csv(
  
  old_summary,
  
  summary_file,
  
  row.names = FALSE,
  
  fileEncoding = "UTF-8"
  
)

run_summary <- old_summary[
  
  old_summary$SRA_File %in% sra_files,
  
  ,
  
  drop = FALSE
  
]

success_count <- sum(
  run_summary$Status == "SUCCESS",
  na.rm = TRUE
)

skipped_count <- sum(
  run_summary$Status == "SKIPPED_COMPLETE",
  na.rm = TRUE
)

failed_count <- sum(
  run_summary$Status == "FAILED",
  na.rm = TRUE
)

cat("\n")

cat("####################################################\n")

cat("             SRA 批处理结束！\n")

cat("####################################################\n")

cat("\n本次处理文件数量：", length(sra_files), "\n")

cat("成功：", success_count, "\n")

cat("跳过：", skipped_count, "\n")

cat("失败：", failed_count, "\n")

cat("\n汇总文件：\n")

cat(summary_file, "\n")

cat("\n全部任务结束！\n")
