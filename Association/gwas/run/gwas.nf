nextflow.enable.dsl = 2

// Utilities
def toCsvList(list) { list ? list.join(',') : '' }


process COMPUTE_PCA {
    tag "plink pca"
    publishDir params.outdir + '/run', mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam)

    output:
    path 'geno.eigenvec', emit: eigen
    path 'geno.eigenval'

    script:
    """
    set -euo pipefail
    plink --bfile geno --pca 10 --out geno
    """
}

process PREPARE_PHENO_AND_COVAR {
    tag "prep pheno/covar"
    publishDir params.outdir + '/run', mode: 'copy'

    input:
    path pheno
    val covar_path
    path fam
    path eigen optional true

    output:
    tuple path('pheno.txt'), path('covar.txt'), emit: pcv

        script:
        def idcol = params.pheno_id_col
        def trait = params.trait
        def covnames = toCsvList(params.covar_names)
        """
        set -euo pipefail
        # Build PLINK pheno.txt with header: FID IID <trait>
        awk -F'[\t,]' 'BEGIN{OFS="\t"}
            NR==1{
                for(i=1;i<=NF;i++){h[tolower($i)]=i}
                print "FID","IID","${trait}"; next
            }
            {
                id=$h[tolower("${idcol}")]; tr=$h[tolower("${trait}")];
                if(id>0 && tr>0){print $id,$id,$tr}
            }
        ' ${pheno} > pheno.txt

        # Build covar.txt: either user covar (FID IID covs...) or PCs from eigenvec
        if [ -n "${covar_path}" ] && [ -s "${covar_path}" ]; then
                    awk -F'[\t,]' 'BEGIN{OFS="\t"}
                NR==1{
                    for(i=1;i<=NF;i++){h[tolower($i)]=i; names[i]=$i}
                    # select covariate columns
                            selc = split(tolower("${covnames}"), sel, /,/); n=0
                            if(selc>0 && sel[1]!=""){
                        for(i=1;i<=NF;i++){ for(j in sel){ if(tolower(names[i])==sel[j]){ pick[++n]=i } } }
                    } else {
                        for(i=1;i<=NF;i++){ if(tolower(names[i])!=tolower("${idcol}")) pick[++n]=i }
                    }
                    # header
                    printf "FID\tIID"; for(i=1;i<=n;i++){ printf "\t%s", names[pick[i]] } printf "\n"; next
                }
                {
                    id=h[tolower("${idcol}")]; if(id>0){ printf "%s\t%s", $id,$id; for(i=1;i<=n;i++){ printf "\t%s", $(pick[i]) } printf "\n" }
                }
            ' ${covar_path} > covar.txt
        else
            if [ -s "${eigen}" ]; then
                # eigenvec columns: FID IID PC1..PC10
                awk 'BEGIN{OFS="\t"} {printf "%s\t%s", $1,$2; for(i=3;i<=12;i++){printf "\t%s", $i} printf "\n"}' geno.eigenvec > covar.txt
                sed -i '1iFID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10' covar.txt
            else
                # minimal covar with only IDs
                awk 'BEGIN{OFS="\t"} {print $1,$2}' geno.fam | sed '1iFID\tIID' > covar.txt
            fi
        fi
        """
}

process RUN_PLINK_GLM {
    tag "plink glm ${params.trait}"
    publishDir params.outdir + "/results/${params.trait}/plink_glm", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam), path(pheno), path(covar)

    output:
    path 'plink.glm.linear', optional: true
    path 'plink.glm.logistic', optional: true

    script:
    """
    set -euo pipefail
    # Try linear; if fails, try logistic
    if plink --bfile geno --pheno pheno.txt --pheno-name ${params.trait} --covar covar.txt --glm --allow-no-sex --out plink ; then
      true
    else
      plink --bfile geno --pheno pheno.txt --pheno-name ${params.trait} --covar covar.txt --logistic --allow-no-sex --out plink || true
    fi
    """
}

process RUN_RMVP_OR_GAPIT {
    tag "${params.gwas_engine} ${params.trait}"
    publishDir params.outdir + "/results/${params.trait}/${params.gwas_engine}", mode: 'copy'

    input:
    tuple path(bed), path(bim), path(fam), path(pheno), path(covar)

    output:
    path 'gwas_results', type: 'dir'

    script:
    def models = toCsvList(params.models*.toUpperCase())
    def engine = params.gwas_engine.toString().toLowerCase()
    def runner = engine == 'gapit' ? "${projectDir}/src/r/gwas/run_gapit.R" : "${projectDir}/src/r/gwas/run_rmvp.R"
    """
    set -euo pipefail
    Rscript ${runner} \
      --bed geno.bed \
      --bim geno.bim \
      --fam geno.fam \
      --pheno pheno.txt \
      --idcol ${params.pheno_id_col} \
      --trait ${params.trait} \
      --covar covar.txt \
      --models ${models} \
      --outdir gwas_results
    """
}

process RUN_GWAS_GAPIT {
    tag "GWAS GAPIT: ${meta.id}"
    publishDir "${params.outdir}/run/gapit", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(pca), path(pheno)

        output:
        path "GAPIT.*", emit: gapit_outputs
        path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Running GWAS with GAPIT for ${meta.id}" > ${prefix}.log

    # 1. Convert VCF to a simple numeric format for R
    # Here we use PLINK's recode A format (0,1,2 for genotypes)
    plink2 \\
        --vcf ${vcf} \\
        --recode A \\
        --out ${prefix}_geno

    # 2. Run the external R script
        Rscript ${projectDir}/src/r/run_gapit.R ${prefix}_geno.raw ${pheno} ${pca} ${prefix} >> ${prefix}.log 2>&1
    """
}

process RUN_GWAS_PLINK {
    tag "GWAS PLINK: ${meta.id}"
    publishDir "${params.outdir}/run/plink", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(pca), path(pheno)

        output:
        path "*.glm.linear", emit: plink_linear
        path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Running GWAS with PLINK for ${meta.id}" > ${prefix}.log

    # Convert VCF to plink format first
    plink2 \\
        --vcf ${vcf} \\
        --make-pgen \\
        --out ${prefix}

    # Run GLM association
    plink2 \\
        --pfile ${prefix} \\
        --pheno ${pheno} \\
        --covar ${pca} \\
        --glm allow-no-covars \\
        --out ${prefix}.gwas

    # The output is ${prefix}.gwas.<pheno_name>.glm.linear
    # We will just publish all glm.linear files
    find . -name "*.glm.linear" -exec mv {} . \\;

    echo "PLINK GWAS complete." >> ${prefix}.log
    """
}

workflow GWAS {
    take:
    ch_vcf
    ch_pheno
    ch_covar

    main:
    ch_plink = VCF_TO_PLINK(ch_vcf).plink
    ch_eigen = COMPUTE_PCA(ch_plink).eigen

    // fam path for pheno/covar prep
    ch_fam = ch_plink.map { bed,bim,fam -> fam }
    ch_pcv = PREPARE_PHENO_AND_COVAR(ch_pheno, ch_covar, ch_fam, ch_eigen).pcv

    // Pair genotype tuple with pheno/covar tuple (single combination)
    ch_pair = ch_plink.combine(ch_pcv).map { pl, pcv ->
        // pl is [bed,bim,fam], pcv is [pheno,covar]
        tuple(pl[0], pl[1], pl[2], pcv[0], pcv[1])
    }

    ch_run1  = RUN_PLINK_GLM(ch_pair)
    ch_run2  = RUN_RMVP_OR_GAPIT(ch_pair)

    emit:
    ch_run1.mix(ch_run2)
}
