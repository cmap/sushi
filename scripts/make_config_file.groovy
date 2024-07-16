import hudson.model.*
import jenkins.model.*
import groovy.json.JsonSlurper

pipeline {
    agent any

    // Define parameters that can be edited via the Jenkins UI
    parameters {
        booleanParam(name: 'TRIGGER_BUILD', defaultValue: true, description: 'Check this to trigger the build. If unchecked, the build will not be triggered and only the config.json will be generated.')
        string(name: 'GIT_BRANCH', defaultValue: 'podman_dev', description: 'Pipeline branch to use')
        string(name: 'COMMIT_HASH', defaultValue: '', description: 'Specific commit hash to use (leave empty to use the latest commit in the branch)')
        string(name: 'BUILD_DIR', defaultValue: '', description: 'Output path to deposit build. Format should be /directory/PROJECT_CODE/BUILD_NAME')
        string(name: 'BUILD_NAME', defaultValue: '', description: 'Build name')
        choice(name: 'SEQ_TYPE', choices: ['DRAGEN', 'MiSeq', 'HiSeq', 'NovaSeq'], description: 'MiSeq and HiSeq/NovaSeq return files named differently. This setting sets the INDEX_1, INDEX_2, and BARCODE_SUFFIX parameters in fastq2readcount. Select DRAGEN if fastq files are from the DRAGEN pipeline from GP. Choosing NovaSeq reverses index 2.')
        string(name: 'CTL_TYPES', defaultValue: 'negcon', description: 'Type to mark as control in compute_LFC')
        string(name: 'DAYS', defaultValue: '', description: 'If running the sushi_to_mts module, provide any days/timepoints (separated by commas) that should be dropped from output data. No quotes needed (ie, 2,8).')
        string(name: 'CELL_SET_META', defaultValue: '', description: 'Cell set metadata')
        string(name: 'ID_COLS', defaultValue: 'cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep', description: 'Columns to concat to create unique ID for each sample-replicate')
        string(name: 'SAMPLE_COLS', defaultValue: 'cell_set,treatment,dose,dose_unit,day,bio_rep', description: 'Sample columns')
        string(name: 'SIG_COLS', defaultValue: 'cell_set,treatment,dose,dose_unit,day', description: 'Signature columns')
        string(name: 'SEQUENCING_INDEX_COLS', defaultValue: 'index_1,index_2,flowcell_names', description: 'Sequencing index columns')
        string(name: 'CONTROL_BARCODE_META', defaultValue: '/data/CB_meta.csv', description: 'Metadata for control barcodes.')
        string(name: 'CONTROL_COLS', defaultValue: 'cell_set,day', description: 'Set of columns that define individual controls')
        booleanParam(name: 'REMOVE_DATA', defaultValue: false, description: 'Select if there is experimental data that needs to be removed before normalization. TODO: expand on this.')
        string(name: 'COUNT_COL_NAME', defaultValue: 'normalized_n', description: 'Field used to calculate L2FC')
        booleanParam(name: 'RUN_NORM', defaultValue: true, description: 'Run normalization module on data.')
        booleanParam(name: 'PULL_POOL_ID', defaultValue: false, description: 'Flag indicating whether to pull pool IDs from CellDB - only applicable to cell sets (i.e. EXT.PR500.CS01.1.A, EXT.PR500.CS01.1.B, etc).')
        booleanParam(name: 'CONVERT_SUSHI', defaultValue: false, description: 'Convert output column headers to format for MTS pipeline.')
        booleanParam(name: 'RUN_EPS_QC', defaultValue: false, description: 'Run EPS QC')
        string(name: 'CELL_SET_META', defaultValue: 'cell_set_meta.csv', description: 'Cell Set Metadata. Static cell_line_meta location: /data/vdb/prismSeq/cell_set_meta.csv')
        string(name: 'SAMPLE_META', defaultValue: 'sample_meta.csv', description: 'File name of sample metadata within the BUILD_DIR directory.')
        string(name: 'COUNT_THRESHOLD', defaultValue: '40', description: 'Minimum threshold to filter cell line counts by.')
        string(name: 'PSEUDOCOUNT', defaultValue: '20', description: 'Pseudocount for normalization.')
        string(name: 'CELL_LINE_META', defaultValue: 'cell_line_meta.csv', description: 'File in BUILD_DIR containing cell line metadata')
        string(name: 'RAW_COUNTS', defaultValue: 'raw_counts.csv', description: 'Filename in BUILD_DIR containing raw counts')
        string(name: 'FILTERED_COUNTS', defaultValue: 'filtered_counts.csv', description: 'File in BUILD_DIR containing filtered counts')
        string(name: 'LFC', defaultValue: 'l2fc.csv', description: 'File containing log2 fold change values')
        string(name: 'ANNOTATED_COUNTS', defaultValue: 'annotated_counts.csv', description: 'File in BUILD_DIR containing annotated counts')
        string(name: 'NORMALIZED_COUNTS', defaultValue: 'normalized_counts.csv', description: 'File in BUILD_DIR containing normalized counts')
        string(name: 'COLLAPSED_VALUES', defaultValue: 'collapsed_values.csv', description: 'File in BUILD_DIR containing replicate collapsed l2fc values')
        string(name: 'API_URL', defaultValue: 'https://api.clue.io/api/', description: 'API URL')
    }

    environment {
        CONFIG_FILE_PATH = "${params.BUILD_DIR}/config.json"
    }

    stages {
        stage('Checkout') {
            steps {
                script {
                    if (params.COMMIT_HASH?.trim()) {
                        // Checkout the specific commit
                        checkout([$class: 'GitSCM',
                                  branches: [[name: params.COMMIT_HASH]],
                                  doGenerateSubmoduleConfigurations: false,
                                  extensions: [[$class: 'RelativeTargetDirectory', relativeTargetDir: '']],
                                  submoduleCfg: [],
                                  userRemoteConfigs: scm.userRemoteConfigs
                        ])
                    } else {
                        // Checkout the branch
                        checkout scm: [$class: 'GitSCM',
                                       branches: [[name: "*/${params.GIT_BRANCH}"]],
                                       doGenerateSubmoduleConfigurations: false,
                                       extensions: [],
                                       userRemoteConfigs: scm.userRemoteConfigs
                        ]
                    }
                }
            }
        }

        stage('Generate JSON Config') {
            steps {
                script {
                    def config = [:]
                    if (fileExists(env.CONFIG_FILE_PATH)) {
                        def configText = readFile(file: env.CONFIG_FILE_PATH)
                        config = new HashMap(new JsonSlurper().parseText(configText))
                    }

                    def paramsMap = [
                        GIT_BRANCH: params.GIT_BRANCH,
                        COMMIT_HASH: params.COMMIT_HASH,
                        BUILD_DIR: params.BUILD_DIR,
                        BUILD_NAME: params.BUILD_NAME,
                        SEQ_TYPE: params.SEQ_TYPE,
                        CTL_TYPES: params.CTL_TYPES,
                        DAYS: params.DAYS,
                        CELL_SET_META: params.CELL_SET_META,
                        ID_COLS: params.ID_COLS,
                        SAMPLE_COLS: params.SAMPLE_COLS,
                        SIG_COLS: params.SIG_COLS,
                        SEQUENCING_INDEX_COLS: params.SEQUENCING_INDEX_COLS,
                        CONTROL_BARCODE_META: params.CONTROL_BARCODE_META,
                        CONTROL_COLS: params.CONTROL_COLS,
                        REMOVE_DATA: params.REMOVE_DATA,
                        COUNT_COL_NAME: params.COUNT_COL_NAME,
                        RUN_NORM: params.RUN_NORM,
                        PULL_POOL_ID: params.PULL_POOL_ID,
                        CONVERT_SUSHI: params.CONVERT_SUSHI,
                        RUN_EPS_QC: params.RUN_EPS_QC,
                        SAMPLE_META: params.SAMPLE_META,
                        COUNT_THRESHOLD: params.COUNT_THRESHOLD,
                        PSEUDOCOUNT: params.PSEUDOCOUNT,
                        CELL_LINE_META: params.CELL_LINE_META,
                        RAW_COUNTS: params.RAW_COUNTS,
                        FILTERED_COUNTS: params.FILTERED_COUNTS,
                        LFC: params.LFC,
                        ANNOTATED_COUNTS: params.ANNOTATED_COUNTS,
                        NORMALIZED_COUNTS: params.NORMALIZED_COUNTS,
                        COLLAPSED_VALUES: params.COLLAPSED_VALUES,
                        API_URL: params.API_URL
                    ]

                    paramsMap.each { key, value ->
                        if (value) {
                            config[key] = value
                        }
                    }

                    if (!config.containsKey('API_KEY')) {
                        config.API_KEY = sh(script: 'cat /local/jenkins/.clue_api_key', returnStdout: true).trim()
                    }

                    writeFile file: env.CONFIG_FILE_PATH, text: groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(config))
                    echo "Generated config.json: ${config}"
                }
            }
        }

        stage('Add Commit Hash to Config') {
            steps {
                script {
                    def commitHash = sh(script: 'git rev-parse HEAD', returnStdout: true).trim()
                    def configFilePath = env.CONFIG_FILE_PATH
                    sh """
                        jq --arg commit "$commitHash
