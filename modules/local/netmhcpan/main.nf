process NETMHCPAN {
    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/de9c5fbcc5583f3c096617ef2c8f84c5e69b479cc5a5944f10d0e1d226779662/data' :
        'community.wave.seqera.io/library/bash_gawk_perl_tcsh:a941b4e9bd4b8805' }"

    input:
    tuple val(meta), path(tsv), path(software)

    output:
    tuple val(meta), path("*.xls"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "I") {
        error "NETMHCPAN only supports MHC class I. Use NETMHCIIPAN for MHC class II."
    }
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def alleles = meta.alleles_supported.tokenize(';').collect { it.replace('*', '').replace('H2','H-2') }.join(',')

    """
    # Create SHORT symlinks to inputs to reduce path lengths
    ln -s $tsv input.tsv
    
    # Determine platform
    UNIX=\$(uname -s)
    AR=\$(uname -m)
    PLATFORM="\${UNIX}_\${AR}"
    
    # Create unique short symlink to netmhcpan using process-specific identifier
    UNIQUE_ID=\$(mktemp -u XXXXXX)
    NM_LINK=/tmp/nm.\${UNIQUE_ID}
    ln -s \${PWD}/netmhcpan \${NM_LINK}
    
    # Set up NetMHCpan environment with SHORT paths
    export NMHOME=\${NM_LINK}
    export NETMHCpan=\${NMHOME}/\${PLATFORM}
    export TMPDIR=/tmp/t.\${UNIQUE_ID}
    mkdir -p \${TMPDIR}
    
    export DTUIBSWWW=www
    export NetMHCpanWWWPATH=/services/NetMHCpan/tmp/
    export NetMHCpanWWWDIR=/usr/opt/www/pub/CBS/services/NetMHCpan/tmp
    
    # Call the binary with short paths - mode 0 (antigen presentation)
    \${NETMHCpan}/bin/netMHCpan-4.2 \
        -p input.tsv \
        -a $alleles \
        -xls \
        -mode 0 \
        -xlsfile mode0.xls \
        $args
    
    # Mode 1 (pathogen epitope)
    \${NETMHCpan}/bin/netMHCpan-4.2 \
        -p input.tsv \
        -a $alleles \
        -xls \
        -mode 1 \
        -xlsfile mode1.xls \
        $args
    
    # TODO: This isn't quite working at the moment
    # Merge horizontally: mode0 columns + mode1 columns (excluding Pos, Peptide, ID)
    paste mode0.xls <(cut -f4- mode1.xls) | \\
        sed '1s/core\\ticore\\tScore\\tRank\\tBA_score\\tBA_Rank\\tAve\\tNB/core_AP\\ticore_AP\\tScore_AP\\tRank_AP\\tBA_score_AP\\tBA_Rank_AP\\tAve_AP\\tNB_AP\\tcore_PE\\ticore_PE\\tScore_PE\\tRank_PE\\tAve_PE\\tNB_PE/' \\
        > ${prefix}_predicted_netmhcpan.xls
    
    # Clean up
    rm -f mode0.xls mode1.xls
    rm -rf \${TMPDIR}
    rm -f \${NM_LINK}
    
    # Verify output was created
    if [ ! -f "${prefix}_predicted_netmhcpan.xls" ]; then
        echo "ERROR: NetMHCpan failed to create output file" >&2
        ls -la
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhcpan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_netmhcpan.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhcpan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}


// **What it does:**
// 1. Uses `paste` to concatenate the files horizontally
// 2. Uses `cut -f4-` to skip the first 3 columns (Pos, Peptide, ID) from mode1 to avoid duplication
// 3. Uses `sed` to rename the headers: `_AP` for Antigen Presentation (mode 0) and `_PE` for Pathogen Epitope (mode 1)

// **Output will look like:**

// Pos	Peptide	ID	core_AP	icore_AP	Score_AP	Rank_AP	BA_score_AP	BA_Rank_AP	Ave_AP	NB_AP	core_PE	icore_PE	Score_PE	Rank_PE	Ave_PE	NB_PE
// 1	QLDIINKN	PEPLIST	QLD-IINKN	QLDIINKN	0.000007	67.313835	0.007122	94.232468	0.000007	0	QLD-IINKN	QLDIINKN	0.375298	64.453041	0.375298	0