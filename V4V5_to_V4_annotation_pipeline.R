
# ============================================================
# Batch SRA Analysis Pipeline
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
# V4-V5 ASV table construction
#  ↓
# Chimera removal
#  ↓
# Unified V4 region extraction
#  ↓
# Collapse identical V4 sequences
#  ↓
# SILVA v138.2 taxonomic assignment
#  ↓
# Genus-level abundance and relative abundance
#  ↓
# PICRUSt2 input generation
#

# Output structure:
# One output directory per SRA file
# All results are saved directly in the sample directory
# ============================================================


# ============================================================
# 0. Configuration
# ============================================================

# Root directory containing SRA files
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

# Use up to 8 threads while reserving one logical core
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

# Applied only to filterAndTrim
FILTER_OMP <- TRUE

# FASTQ filtering chunk size
FILTER_CHUNK_SIZE <- 100000L

# FASTQ dereplication chunk size
DEREP_CHUNK_SIZE <- 1000000L

# Maximum number of bases used for error-rate learning
ERROR_NBASES <- 1e8


# ------------------------------------------------------------
# DADA2 parameters
# ------------------------------------------------------------

MAX_EE <- c(2, 2)

MAX_N <- 0

TRUNC_Q <- 2

# Disable sequence truncation
TRUNC_LEN <- c(0, 0)


# ------------------------------------------------------------
# V4 extraction parameters
# ------------------------------------------------------------

# V4 extraction from V4-V5 amplicons uses the 515F/806R primer pair
# Primer sequences must match the amplicon primers

PRIMER_515F <- "GTGYCAGCMGCCGCGGTAA"

PRIMER_806R <- "GGACTACNVGGGTWTCTAAT"

PRIMER_MAX_MISMATCH <- 3

V4_MIN_LENGTH <- 180

V4_MAX_LENGTH <- 300

V4_WARNING_RETENTION <- 90


# ------------------------------------------------------------
# Batch processing parameters
# ------------------------------------------------------------

# Enable single-sample test mode for an initial run
# Disable after validation to process all SRA files

TEST_ONLY <- FALSE

TEST_SRA_ID <- "11_Ctrl_1_ERR11788629"

# Skip samples with complete results
SKIP_COMPLETED <- TRUE

# FALSE: remove FASTQ and filtered files after successful processing
# TRUE: retain FASTQ and filtered files
KEEP_INTERMEDIATE_FILES <- FALSE


# ============================================================
# 1. Check and install required R packages
# ============================================================

cat("\nChecking required R packages...\n")

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
    "Packages to install: ",
    paste(missing_packages, collapse = ", "),
    "\n"
  )
  
  BiocManager::install(
    missing_packages,
    update = FALSE,
    ask = FALSE
  )
  
} else {
  
  cat("All required R packages are already installed.\n")
  
}

suppressPackageStartupMessages({
  
  library(dada2)
  
  library(Biostrings)
  
})

# Configure the OpenMP environment
Sys.setenv(
  OMP_NUM_THREADS = DADA_THREADS
)

cat("\nCPU configuration:\n")

cat(
  "Logical cores: ",
  available_cores,
  "\n"
)

cat(
  "DADA2 threads: ",
  DADA_THREADS,
  "\n"
)


# ============================================================
# 2. Validate input paths
# ============================================================

if (!dir.exists(input_root)) {
  stop("The root SRA directory does not exist.")
}

if (!file.exists(fastq_dump_path)) {
  stop("fastq-dump-orig.exe was not found.")
}

if (!file.exists(silva_genus)) {
  stop("The SILVA v138.2 reference database was not found.")
}

cat("\nSRA conversion executable:\n")
cat(fastq_dump_path, "\n")

cat("\nSILVA database:\n")
cat(silva_genus, "\n")


# ============================================================
# 3. Define the FASTQ search function
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
  
  # Prefer FASTQ files in the sample directory root
  root_files <- files[
    file_dirs == root_path
  ]
  
  if (length(root_files) == 1) {
    return(root_files)
  }
  
  if (length(root_files) > 1) {
    stop(
      paste(
        "Multiple FASTQ files were found for mate:",
        mate
      )
    )
  }
  
  # Support legacy FASTQ subdirectories
  if (length(files) == 1) {
    return(files)
  }
  
  stop(
    paste(
      "Multiple ambiguous FASTQ files were found for mate:",
      mate
    )
  )
  
}


# ============================================================
# 4. Validate FASTQ input files
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
# 5. Configure V4 primers
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
  
  # Locate the 515F primer
  hits_F <- matchPattern(
    
    primer_515F,
    
    subject,
    
    max.mismatch = PRIMER_MAX_MISMATCH,
    
    with.indels = FALSE,
    
    fixed = "subject"
    
  )
  
  # Locate the reverse complement of the 806R primer
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
  
  # Pair forward- and reverse-primer positions
  pairs <- expand.grid(
    
    fi = seq_along(f_end),
    
    ri = seq_along(r_start)
    
  )
  
  # Require the forward primer to precede the reverse primer
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
  
  # Retain sequences within the expected V4 length range
  valid <- which(
    
    insert_length >= V4_MIN_LENGTH &
      
      insert_length <= V4_MAX_LENGTH
    
  )
  
  if (length(valid) == 0) {
    return(NA_character_)
  }
  
  # Select the region closest to 250 bp
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
# 7. Check sample completion status
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
# 8. Define the per-SRA processing function
# ============================================================

process_one_sra <- function(sra_file) {
  
  sra_id <- tools::file_path_sans_ext(
    basename(sra_file)
  )
  
  # Create one output directory per SRA file
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
      "\nComplete results already exist; skipping: ",
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
        
        message = "Complete results already exist."
        
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
  
  # Remove the previous success marker before reprocessing
  if (file.exists(success_marker)) {
    file.remove(success_marker)
  }
  
  cat("\n")
  cat("====================================================\n")
  
  cat("Processing: ", sra_id, "\n")
  
  cat("SRA file: ", sra_file, "\n")
  
  cat("Output directory: ", sample_dir, "\n")
  
  cat("====================================================\n")
  
  
  # ----------------------------------------------------------
  # 8.3 Convert SRA to FASTQ
  # ----------------------------------------------------------
  
  cat("\nChecking FASTQ files...\n")
  
  fnF <- find_fastq(
    sample_dir,
    "1"
  )
  
  fnR <- find_fastq(
    sample_dir,
    "2"
  )
  
  if (!valid_fastq_pair(fnF, fnR)) {
    
    cat("\nConverting SRA to FASTQ...\n")
    
    # Remove incomplete FASTQ files from previous attempts
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
    
    # Use fastq-dump-orig.exe explicitly
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
    
    cat("\nExecuting command:\n")
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
      
      stop("fastq-dump-orig.exe conversion failed.")
      
    }
    
  } else {
    
    cat("\nComplete FASTQ files already exist; skipping conversion.\n")
    
  }
  
  
  # ----------------------------------------------------------
  # 8.4 Validate FASTQ files
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
      "FASTQ conversion completed, but a complete paired-end file set was not found."
    )
    
  }
  
  # Move legacy subdirectory FASTQ files to the sample root
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
        
        stop("Failed to move FASTQ files to the sample output directory root.")
        
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
  
  cat("\nStarting quality filtering...\n")
  
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
  
  # Remove previous filtered files to ensure a clean run
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
    
    # Disable fork-based parallelism on Windows
    multithread = FILTER_MULTITHREAD,
    
    # Apply OMP only to filterAndTrim
    OMP = FILTER_OMP,
    
    n = FILTER_CHUNK_SIZE,
    
    verbose = TRUE
    
  )
  
  cat("\nQuality-filtering results:\n")
  
  print(ft_out)
  
  if (
    is.null(ft_out) ||
    nrow(ft_out) != 1 ||
    is.na(ft_out[1, "reads.out"]) ||
    ft_out[1, "reads.out"] <= 0
  ) {
    
    stop("Quality filtering failed or no reads remain.")
    
  }
  
  if (
    !file.exists(filtF) ||
    !file.exists(filtR)
  ) {
    
    stop("Complete filtered FASTQ files were not generated.")
    
  }
  
  
  # ----------------------------------------------------------
  # 8.6 Learn forward error rates
  # ----------------------------------------------------------
  
  cat("\nLearning forward error rates...\n")
  
  errF <- dada2::learnErrors(
    
    filtF,
    
    nbases = ERROR_NBASES,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # 8.7 Learn reverse error rates
  # ----------------------------------------------------------
  
  cat("\nLearning reverse error rates...\n")
  
  errR <- dada2::learnErrors(
    
    filtR,
    
    nbases = ERROR_NBASES,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  
  # ----------------------------------------------------------
  # 8.8 FASTQ dereplication
  # ----------------------------------------------------------
  
  cat("\nStarting forward FASTQ dereplication...\n")
  
  # derepFastq does not support the OMP argument
  derepF <- dada2::derepFastq(
    
    filtF,
    
    n = DEREP_CHUNK_SIZE,
    
    verbose = TRUE
    
  )
  
  cat("\nStarting reverse FASTQ dereplication...\n")
  
  derepR <- dada2::derepFastq(
    
    filtR,
    
    n = DEREP_CHUNK_SIZE,
    
    verbose = TRUE
    
  )
  
  # Validate returned object classes
  if (!inherits(derepF, "derep")) {
    
    stop("Forward derepFastq did not return a derep object.")
    
  }
  
  if (!inherits(derepR, "derep")) {
    
    stop("Reverse derepFastq did not return a derep object.")
    
  }
  
  cat("\nFASTQ dereplication completed successfully.\n")
  
  
  # ----------------------------------------------------------
  # 8.9 Forward-read DADA2 denoising
  # ----------------------------------------------------------
  
  cat("\nStarting forward DADA2 denoising...\n")
  
  dadaF <- dada2::dada(
    
    derepF,
    
    err = errF,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  if (!inherits(dadaF, "dada")) {
    
    stop("Forward DADA2 did not return a dada object.")
    
  }
  
  
  # ----------------------------------------------------------
  # 8.10 Reverse-read DADA2 denoising
  # ----------------------------------------------------------
  
  cat("\nStarting reverse DADA2 denoising...\n")
  
  dadaR <- dada2::dada(
    
    derepR,
    
    err = errR,
    
    multithread = DADA_THREADS,
    
    verbose = TRUE
    
  )
  
  if (!inherits(dadaR, "dada")) {
    
    stop("Reverse DADA2 did not return a dada object.")
    
  }
  
  cat("\nForward and reverse DADA2 denoising completed successfully.\n")
  
  
  # ----------------------------------------------------------
  # 8.11 Merge paired-end reads
  # ----------------------------------------------------------
  
  cat("\nStarting paired-end merging...\n")
  
  # Use a single-sample object
  # mergePairs returns a data.frame
  
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
  # Validate the mergePairs return type
  # ----------------------------------------------------------
  
  if (!is.data.frame(merger)) {
    
    stop(
      paste(
        "mergePairs returned an unexpected data type:",
        paste(class(merger), collapse = ", ")
      )
    )
    
  }
  
  if (!"abundance" %in% colnames(merger)) {
    
    stop("The mergePairs result does not contain an abundance column.")
    
  }
  
  # Sum reads from successfully merged pairs
  # Avoid applying vapply to the data.frame
  
  total_merged <- sum(
    
    merger$abundance,
    
    na.rm = TRUE
    
  )
  
  cat("\nPaired-end merging statistics:\n")
  
  cat(
    "Successfully merged unique pairings: ",
    nrow(merger),
    "\n"
  )
  
  cat(
    "Successfully merged reads: ",
    total_merged,
    "\n"
  )
  
  merge_rate <- (
    total_merged /
      as.numeric(ft_out[1, "reads.out"])
  ) * 100
  
  cat(
    "Paired-end merge rate: ",
    round(merge_rate, 2),
    "%\n"
  )
  
  if (
    length(total_merged) != 1 ||
    is.na(total_merged) ||
    !is.finite(total_merged) ||
    total_merged <= 0
  ) {
    
    stop("Paired-end merging produced no valid reads.")
    
  }
  
  # Save the mergePairs output for downstream inspection
  saveRDS(
    
    merger,
    
    file.path(
      sample_dir,
      "mergePairs_result.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.12 Build the V4-V5 ASV table
  # ----------------------------------------------------------
  
  cat("\nConstructing the V4-V5 ASV table...\n")
  
  # makeSequenceTable expects a list of sample-level merge results
  # Wrap the single-sample merger in a list here
  
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
    
    stop("Failed to construct the V4-V5 ASV table.")
    
  }
  
  rownames(seqtab) <- sra_id
  
  cat("\nV4-V5 ASV table constructed successfully.\n")
  
  cat(
    "Number of V4-V5 ASVs: ",
    ncol(seqtab),
    "\n"
  )
  
  saveRDS(
    
    seqtab,
    
    file.path(
      sample_dir,
      "V4V5_ASV_table_before_chimera.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.13 Remove chimeric sequences
  # ----------------------------------------------------------
  
  cat("\nRemoving chimeras...\n")
  
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
    
    stop("No ASVs remained after chimera removal.")
    
  }
  
  rownames(seqtab.nochim) <- sra_id
  
  saveRDS(
    
    seqtab.nochim,
    
    file.path(
      sample_dir,
      "V4V5_ASV_table_nochim.rds"
    )
    
  )
  
  original_asv_count <- ncol(
    seqtab.nochim
  )
  
  nonchim_reads <- sum(
    seqtab.nochim
  )
  
  cat("\nChimera removal completed successfully.\n")
  
  cat(
    "Number of non-chimeric ASVs: ",
    original_asv_count,
    "\n"
  )
  
  cat(
    "Non-chimeric reads: ",
    nonchim_reads,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 8.14 Extract the V4 region
  # ----------------------------------------------------------
  
  cat("\nExtracting the standardized V4 region...\n")
  
  v4v5_sequences <- colnames(
    seqtab.nochim
  )
  
  v4_sequences <- vapply(
    
    v4v5_sequences,
    
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
        "V4 extraction failed for all ASVs.",
        "Check the amplification primers and V4 extraction primers."
      )
    )
    
  }
  
  # ASV-level extraction success rate
  v4_asv_success_rate <- mean(
    v4_success
  ) * 100
  
  # Read-weighted V4 retention rate
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
  
  cat("\nV4 extraction statistics:\n")
  
  cat(
    "Original V4-V5 ASVs: ",
    original_asv_count,
    "\n"
  )
  
  cat(
    "ASVs with successfully extracted V4 regions: ",
    sum(v4_success),
    "\n"
  )
  
  cat(
    "ASV extraction success rate: ",
    round(v4_asv_success_rate, 2),
    "%\n"
  )
  
  cat(
    "V4 read retention rate: ",
    round(v4_retention, 2),
    "%\n"
  )
  
  if (v4_retention < V4_WARNING_RETENTION) {
    
    warning(
      paste0(
        sra_id,
        " V4 read retention is below ",
        V4_WARNING_RETENTION,
        "%; verify primer matching."
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # 8.15 Export the V4 extraction quality report
  # ----------------------------------------------------------
  
  v4_report <- data.frame(
    
    Original_ASV = paste0(
      "V4V5_ASV_",
      seq_along(v4v5_sequences)
    ),
    
    V4V5_Sequence = v4v5_sequences,
    
    V4V5_Length = nchar(v4v5_sequences),
    
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
  
  cat("\nCollapsing identical V4 sequences...\n")
  
  asv_counts_kept <- t(
    
    seqtab.nochim[
      ,
      v4_success,
      drop = FALSE
    ]
    
  )
  
  # Sum read counts across ASVs sharing the same V4 sequence
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
    stop("No ASVs remained after collapsing identical V4 sequences.")
  }
  
  cat(
    "\nFinal V4 ASV count: ",
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
  # 8.17 Generate standardized V4 ASV identifiers
  # ----------------------------------------------------------
  
  v4_ids <- paste0(
    
    "V4_ASV_",
    
    sprintf(
      "%06d",
      seq_along(v4_unique_sequences)
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.18 Build the V4 ASV abundance matrix
  # ----------------------------------------------------------
  
  v4_count_matrix <- t(
    v4_seqtab
  )
  
  rownames(v4_count_matrix) <- v4_unique_sequences
  
  colnames(v4_count_matrix) <- sra_id
  
  
  # ----------------------------------------------------------
  # 8.19 Export the PICRUSt2 FASTA file
  # ----------------------------------------------------------
  
  cat("\nGenerating the PICRUSt2 FASTA file...\n")
  
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
  # 8.20 Export the PICRUSt2 abundance TSV
  # ----------------------------------------------------------
  
  cat("\nGenerating the PICRUSt2 abundance TSV...\n")
  
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
  
  cat("\nAssigning taxonomy with SILVA v138.2...\n")
  
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
    
    stop("SILVA taxonomic assignment failed.")
    
  }
  
  saveRDS(
    
    v4_taxa,
    
    file.path(
      sample_dir,
      "V4_SILVA_taxonomy.rds"
    )
    
  )
  
  
  # ----------------------------------------------------------
  # 8.22 Build the V4 ASV taxonomy table
  # ----------------------------------------------------------
  
  cat("\nConstructing the V4 ASV taxonomy table...\n")
  
  if (is.null(rownames(v4_taxa))) {
    
    stop("The taxonomy results are missing sequence row names.")
    
  }
  
  taxa_index <- match(
    
    v4_unique_sequences,
    
    rownames(v4_taxa)
    
  )
  
  if (anyNA(taxa_index)) {
    
    stop("The taxonomy results could not be matched to the V4 ASV sequences.")
    
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
  
  # Store the representative V4 sequence in the first column
  # Retain ASV_ID for correspondence with PICRUSt2 files
  
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
  # 8.23 Summarize genus-level abundance
  # ----------------------------------------------------------
  
  cat("\nSummarizing genus-level abundance...\n")
  
  if (!"Genus" %in% colnames(taxa_match)) {
    
    stop("The SILVA taxonomy results do not contain a Genus column.")
    
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
  
  # Aggregate read counts for ASVs assigned to the same genus
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
  
  cat("\nCalculating genus-level relative abundance...\n")
  
  genus_totals <- colSums(
    genus_matrix
  )
  
  if (any(genus_totals <= 0)) {
    
    stop("Total genus-level abundance is zero; relative abundance cannot be calculated.")
    
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
  
  cat("\nGenerating the read-tracking table...\n")
  
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
  # 8.26 Save complete analysis results
  # ----------------------------------------------------------
  
  cat("\nSaving the complete analysis results...\n")
  
  complete_results <- list(
    
    SRA_ID = sra_id,
    
    DADA2_forward = dadaF,
    
    DADA2_reverse = dadaR,
    
    Merged_reads = merger,
    
    V4V5_ASV_before_chimera = seqtab,
    
    V4V5_ASV = seqtab.nochim,
    
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
    
    cat("\nRemoving intermediate FASTQ files...\n")
    
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
  # 8.28 Write the success marker
  # ----------------------------------------------------------
  
  writeLines(
    
    c(
      
      paste("SRA:", sra_id),
      
      paste("Completed:", Sys.time()),
      
      paste("Merged reads:", total_merged),
      
      paste("V4-V5 ASV:", original_asv_count),
      
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
  
  cat("Sample processed successfully: ", sra_id, "\n")
  
  cat("Paired-end merged reads: ", total_merged, "\n")
  
  cat("V4-V5 ASV count: ", original_asv_count, "\n")
  
  cat("V4 ASV count: ", ncol(v4_seqtab), "\n")
  
  cat(
    "V4 read retention rate: ",
    round(v4_retention, 2),
    "%\n"
  )
  
  cat(
    "Number of genera: ",
    nrow(genus_matrix),
    "\n"
  )
  
  cat("Output directory: ", sample_dir, "\n")
  
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
      
      message = "Processing completed successfully."
      
    )
    
  )
  
}


# ============================================================
# 9. Discover all SRA files
# ============================================================

cat("\nSearching for SRA files...\n")

sra_files <- list.files(
  
  path = input_root,
  
  pattern = "\\.sra$",
  
  full.names = TRUE,
  
  recursive = TRUE,
  
  ignore.case = TRUE
  
)

if (length(sra_files) == 0) {
  
  stop("No SRA files were found.")
  
}

sra_files <- normalizePath(
  
  sra_files,
  
  winslash = "/",
  
  mustWork = TRUE
  
)

cat(
  "\nNumber of SRA files found: ",
  length(sra_files),
  "\n"
)


# ============================================================
# 10. Single-sample test mode
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
        "The test sample was not found or its ID is not unique:",
        TEST_SRA_ID
      )
    )
    
  }
  
  sra_files <- sra_files[selected]
  
  cat("\nSingle-sample test mode is enabled.\n")
  
  cat(
    "Test file: ",
    basename(sra_files),
    "\n"
  )
  
} else {
  
  cat("\nBatch processing mode is enabled for all SRA files.\n")
  
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
  
  "Original_V4V5_ASV",
  
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
  
  Original_V4V5_ASV = numeric(),
  
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
    "Current progress: ",
    i,
    "/",
    length(sra_files),
    "\n"
  )
  
  cat(
    "Current file: ",
    basename(current_sra),
    "\n"
  )
  
  cat("####################################################\n")
  
  
  # ----------------------------------------------------------
  # Isolate errors at the sample level
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
      
      cat("\nProcessing failed: ", current_id, "\n")
      
      cat("Error details: ", error_message, "\n")
      
      # Save the error call stack
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
  # Calculate elapsed time
  # ----------------------------------------------------------
  
  elapsed_seconds <- as.numeric(
    
    difftime(
      
      Sys.time(),
      
      start_time,
      
      units = "secs"
      
    )
    
  )
  
  
  # ----------------------------------------------------------
  # Build the current sample summary row
  # ----------------------------------------------------------
  
  result_row <- data.frame(
    
    SRA_ID = result$sra_id,
    
    SRA_File = current_sra,
    
    Status = result$status,
    
    Original_V4V5_ASV = result$original_asv,
    
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
  # Update the batch summary table
  # ----------------------------------------------------------
  
  existing <- which(
    
    !is.na(old_summary$SRA_File) &
      
      old_summary$SRA_File == current_sra
    
  )
  
  # Preserve existing records for completed samples that are skipped
  # Retain the prior summary entry when one is already available
  
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
  
  # Save the summary after each sample
  write.csv(
    
    old_summary,
    
    summary_file,
    
    row.names = FALSE,
    
    fileEncoding = "UTF-8"
    
  )
  
  
  cat("\nCurrent sample status: ", result$status, "\n")
  
  cat(
    "Elapsed time: ",
    round(elapsed_seconds / 60, 2),
    " minutes\n"
  )
  
  cat(
    "Completed progress: ",
    i,
    "/",
    length(sra_files),
    "\n"
  )
  
}


# ============================================================
# 13. Final batch summary
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

cat("             SRA batch processing completed!\n")

cat("####################################################\n")

cat("\nNumber of files processed in this run: ", length(sra_files), "\n")

cat("Successful: ", success_count, "\n")

cat("Skipped: ", skipped_count, "\n")

cat("Failed: ", failed_count, "\n")

cat("\nSummary file:\n")

cat(summary_file, "\n")

cat("\nAll tasks completed!\n")
