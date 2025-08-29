process DADA2_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.svg"), emit: svg
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "versions.yml"               , emit: versions
    path "*.args.txt"                 , emit: args


    script:
    def prefix = task.ext.prefix ?: "prefix"
    def args = task.ext.args ?: ''
    def seed = task.ext.seed ?: '100'
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        safe_learn <- function(files, fname) {
            total_reads <- sum(sapply(files, function(f) length(readLines(gzfile(f)))/4))
            if(total_reads < 50) {   # threshold for DADA2 error learning
                cat("WARNING:", fname, "has too few reads (", total_reads, "). Skipping learnErrors.\n")
                saveRDS(NULL, paste0(fname, ".rds"))
                pdf(paste0(fname, ".pdf")); plot.new(); dev.off()
                svg(paste0(fname, ".svg")); plot.new(); dev.off()
                writeLines("Too few reads to estimate errors.", paste0(fname, ".convergence.txt"))
                return(NULL)
            } else {
                err <- learnErrors(files, nbases = 1e8, nreads = NULL, randomize = TRUE,
                                MAX_CONSIST = 10, OMEGA_C = 0, qualityType = "Auto",
                                errorEstimationFunction = loessErrfun, multithread = 6, verbose = TRUE)
                saveRDS(err, paste0(fname, ".rds"))
                pdf(paste0(fname, ".pdf")); plotErrors(err, nominalQ = TRUE); dev.off()
                svg(paste0(fname, ".svg")); plotErrors(err, nominalQ = TRUE); dev.off()
                sink(paste0(fname, ".convergence.txt"))
                dada2:::checkConvergence(err)
                sink()
                return(err)
            }
        }

        safe_learn <- function(files, fname) {
            total_reads <- sum(sapply(files, function(f) length(readLines(gzfile(f)))/4))
            if(total_reads < 50) {   # threshold for DADA2 error learning
                cat("WARNING:", fname, "has too few reads (", total_reads, "). Skipping learnErrors.\n")
                saveRDS(NULL, paste0(fname, ".rds"))
                pdf(paste0(fname, ".pdf")); plot.new(); dev.off()
                svg(paste0(fname, ".svg")); plot.new(); dev.off()
                writeLines("Too few reads to estimate errors.", paste0(fname, ".convergence.txt"))
                return(NULL)
            } else {
                err <- learnErrors(files, nbases = 1e8, nreads = NULL, randomize = TRUE,
                                MAX_CONSIST = 10, OMEGA_C = 0, qualityType = "Auto",
                                errorEstimationFunction = loessErrfun, multithread = 6, verbose = TRUE)
                saveRDS(err, paste0(fname, ".rds"))
                pdf(paste0(fname, ".pdf")); plotErrors(err, nominalQ = TRUE); dev.off()
                svg(paste0(fname, ".svg")); plotErrors(err, nominalQ = TRUE); dev.off()
                sink(paste0(fname, ".convergence.txt"))
                dada2:::checkConvergence(err)
                sink()
                return(err)
            }
        }

        # Create empty placeholder files first
        files_to_create <- c(
            "${prefix}_1.err.rds",
            "${prefix}_2.err.rds",
            "${prefix}_1.err.pdf",
            "${prefix}_2.err.pdf",
            "${prefix}_1.err.svg",
            "${prefix}_2.err.svg",
            "${prefix}_1.err.log",
            "${prefix}_2.err.log",
            "${prefix}_1.err.convergence.txt",
            "${prefix}_2.err.convergence.txt"
        )

        sapply(files_to_create, function(f) {
            if(!file.exists(f)) file.create(f)
        })

        errF <- safe_learn(fnFs, "${prefix}_1")
        errR <- safe_learn(fnRs, "${prefix}_2")
        sink(file = NULL)

        write.table('learnErrors        nbases = 1e8, nreads = NULL, randomize = TRUE, MAX_CONSIST = 10, OMEGA_C = 0,qualityType = "Auto",errorEstimationFunction = loessErrfun',
                    file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

        writeLines(c("\"NFCORE_AMPLISEQ:AMPLISEQ:DADA2_ERR\":",
                    paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                    paste0("    dada2: ", packageVersion("dada2")) ),
                "versions.yml")


        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${prefix}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${prefix}.err.rds")
        sink(file = NULL)

        pdf("${prefix}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}


        // fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE), method = "radix")
        // fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE), method = "radix")

        // sink(file = "${prefix}.err.log")
        // errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        // saveRDS(errF, "${prefix}_1.err.rds")
        // errR <- learnErrors(fnRs, $args, multithread = $task.cpus, verbose = TRUE)
        // saveRDS(errR, "${prefix}_2.err.rds")
        // sink(file = NULL)

        // pdf("${prefix}_1.err.pdf")
        // plotErrors(errF, nominalQ = TRUE)
        // dev.off()
        // svg("${prefix}_1.err.svg")
        // plotErrors(errF, nominalQ = TRUE)
        // dev.off()

        // pdf("${prefix}_2.err.pdf")
        // plotErrors(errR, nominalQ = TRUE)
        // dev.off()
        // svg("${prefix}_2.err.svg")
        // plotErrors(errR, nominalQ = TRUE)
        // dev.off()

        // sink(file = "${prefix}_1.err.convergence.txt")
        // dada2:::checkConvergence(errF)
        // sink(file = NULL)

        // sink(file = "${prefix}_2.err.convergence.txt")
        // dada2:::checkConvergence(errR)
        // sink(file = NULL)

        // write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        // writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
