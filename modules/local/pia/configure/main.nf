process PIA_CONFIGURE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pia:1.5.9--hdfd78af_0':
        'quay.io/biocontainers/pia:1.5.9--hdfd78af_0' }"

    input:
    val meta

    output:
    tuple val(meta), path("*.pia.json"), emit: pia_json
    tuple val("${task.process}"), val('pia'), eval('pia --version 1>&1 | head -1  | sed "s/.*version //"'), topic: versions, emit: versions_pia

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def psm_export = meta.psm_export != null ? meta.psm_export : true
    def peptide_export = meta.peptide_export != null ? meta.peptide_export : true
    def protein_export = meta.protein_export != null ? meta.protein_export : true
    def fdr_filter = meta.fdr_filter != null ? meta.fdr_filter : true
    def remove_decoys = meta.remove_decoys != null ? meta.remove_decoys : true
    def fdr_threshold = meta.fdr_threshold ?: 0.01
    def psm_export_format = meta.psm_export_format ?: 'mzTab'

    """
    pia --example > pia_analysis.json

    # delete the first row of the file, which contains an explanatory comment and does not start with the json
    sed -i 1d pia_analysis.json

    sed -i 's;"createPSMsets": .*,;"createPSMsets": false,;g' pia_analysis.json
    sed -i 's;"psmLevelFileID": .*,;"psmLevelFileID": 1,;g' pia_analysis.json
    sed -i 's;"topIdentifications": .*,;"topIdentifications": 1,;g' pia_analysis.json
    sed -i 's;"calculateCombinedFDRScore": .*,;"calculateCombinedFDRScore": false,;g' pia_analysis.json
    if [ ${psm_export} = true ];
    then
      sed -i 's;"psmExportFile":.*,;"psmExportFile": "piaExport-PSMs.${psm_export_format}",;g' pia_analysis.json
    else
      sed -i 's;"psmExportFile":.*;;g' pia_analysis.json
    fi

    if [ ${peptide_export} = true ];
    then
      sed -i 's;"inferePeptides":.*,;"inferePeptides": true,;g' pia_analysis.json
      sed -i 's;"peptideExportFile":.*,;"peptideExportFile": "piaExport-peptides.csv",;g' pia_analysis.json
    else
      sed -i 's;"inferePeptides":.*,;"inferePeptides": false,;g' pia_analysis.json
      sed -i 's;"peptideExportFile":.*,;;g' pia_analysis.json
    fi
    sed -i 's;"peptideLevelFileID":.*,;"peptideLevelFileID": 1,;g' pia_analysis.json

    if [ ${protein_export} = true ];
    then
      sed -i 's;"infereProteins":.*,;"infereProteins": true,;g' pia_analysis.json
      sed -i 's;"inferenceMethod":.*,;"inferenceMethod": "inference_spectrum_extractor",;g' pia_analysis.json
      sed -i '/inferenceFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
      sed -i 's;"scoringBaseScore":.*,;"scoringBaseScore": "psm_fdr_score",;g' pia_analysis.json
      sed -i 's;"scoringPSMs":.*,;"scoringPSMs": "best",;g' pia_analysis.json
      sed -i 's;"proteinExportFile":.*,;"proteinExportFile": "piaExport-proteins.mzTab",;g' pia_analysis.json
      sed -i 's;"proteinExportWithPSMs":.*,;"proteinExportWithPSMs": true,;g' pia_analysis.json
    else
      sed -i 's;"infereProteins":.*,;"infereProteins": false,;g' pia_analysis.json
      sed -i 's;"proteinExportFile":.*,;;g' pia_analysis.json
    fi

    if [ ${fdr_filter} = true ];
    then
      if [ ${remove_decoys} = true ];
      then
        sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}",\\n    "psm_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
        sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}",\\n    "peptide_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
        sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= ${fdr_threshold}",\\n    "protein_accessions_filter !regex_only DECOY_.*"/g;}' pia_analysis.json
      else
        sed -i '/psmFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
        sed -i '/peptideFilters/{n;s/.*/    "psm_score_filter_psm_fdr_score <= ${fdr_threshold}"/g;}' pia_analysis.json
        sed -i '/proteinFilters/{n;s/.*/    "protein_q_value_filter <= ${fdr_threshold}"/g;}' pia_analysis.json
      fi
    else
      sed -i '/psmFilters/{n;s/.*//g;}' pia_analysis.json
      sed -i '/peptideFilters/{n;s/.*//g;}' pia_analysis.json
      sed -i '/proteinFilters/{n;s/.*//g;}' pia_analysis.json
    fi

    mv pia_analysis.json ${prefix}.pia.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pia.json
    """
}
